#!/usr/bin/env python3
"""Run a reproducible low-pT MadGraph campaign one sample at a time.

The runner deliberately owns no physics defaults beyond the process registry.
All user-facing production values live in the TOML configuration. Eta is not
patched: MadGraph's process-specific default eta cuts are preserved.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import tomllib
from pathlib import Path
from typing import Any, Iterator


class CampaignError(RuntimeError):
    """Raised when a campaign cannot safely continue."""


@dataclasses.dataclass(frozen=True)
class Sample:
    name: str
    process: str
    directory: str
    target_label: str | None
    heavy_b: bool = False


SAMPLES: dict[str, Sample] = {
    "bbbar": Sample("bbbar", "p p > b b~", "QCD_bbbar_lowpt", "b", heavy_b=True),
    "ccbar": Sample("ccbar", "p p > c c~", "QCD_ccbar_lowpt", "c"),
    "gg": Sample("gg", "p p > g g", "QCD_gg_lowpt", "g"),
    "uubar": Sample("uubar", "p p > u u~", "QCD_uds_u_lowpt", "uds"),
    "ddbar": Sample("ddbar", "p p > d d~", "QCD_uds_d_lowpt", "uds"),
    "ssbar": Sample("ssbar", "p p > s s~", "QCD_uds_s_lowpt", "uds"),
    "inclusive": Sample("inclusive", "p p > j j", "QCD_inclusive_lowpt", None),
}

REQUIRED_BRANCHES = (
    "Jet.PT",
    "Jet.Flavor",
    "Jet.Eta",
    "Jet.Phi",
    "Jet.Mass",
    "Jet.NCharged",
    "Jet.NNeutrals",
    "Track.PT",
    "Track.Eta",
    "Track.Phi",
    "Track.Mass",
    "Track.Charge",
    "Track.D0",
    "Track.DZ",
)

TRUTH_LABELS = {5: "b", 4: "c", 21: "g", 1: "uds", 2: "uds", 3: "uds"}


@dataclasses.dataclass(frozen=True)
class CampaignConfig:
    name: str
    mg5_root: Path
    output_root: Path
    nevents: int
    runs_per_sample: int
    cores: int
    ebeam1: float
    ebeam2: float
    parton_pt_min: float
    light_parton_dr_min: float
    bb_dr_min: float
    jet_radius: float
    reco_jet_pt_min: float
    keep_hepmc: bool
    enabled_samples: tuple[str, ...]
    seeds: dict[str, int]
    config_hash: str


def load_config(path: Path) -> CampaignConfig:
    """Load and validate the intentionally small editable campaign surface."""
    try:
        raw = tomllib.loads(path.read_text())
    except FileNotFoundError as error:
        raise CampaignError(f"Campaign config does not exist: {path}") from error
    except tomllib.TOMLDecodeError as error:
        raise CampaignError(f"Campaign config is not valid TOML: {error}") from error

    campaign = raw.get("campaign")
    samples = raw.get("samples")
    seeds = raw.get("seeds")
    if not isinstance(campaign, dict) or not isinstance(samples, dict) or not isinstance(seeds, dict):
        raise CampaignError("Config requires [campaign], [samples], and [seeds] sections.")
    if "parton_eta_max" in campaign or "eta" in campaign:
        raise CampaignError("Eta is intentionally not configurable; the runner preserves MadGraph defaults.")

    def required(name: str) -> Any:
        if name not in campaign:
            raise CampaignError(f"Missing campaign setting: {name}")
        return campaign[name]

    enabled = samples.get("enabled")
    if not isinstance(enabled, list) or not enabled:
        raise CampaignError("[samples].enabled must be a non-empty list.")
    unknown = set(enabled) - set(SAMPLES)
    if unknown:
        raise CampaignError(f"Unknown samples in config: {', '.join(sorted(unknown))}")
    if len(set(enabled)) != len(enabled):
        raise CampaignError("[samples].enabled cannot contain duplicates.")

    parsed_seeds: dict[str, int] = {}
    for sample in enabled:
        value = seeds.get(sample)
        if not isinstance(value, int) or value < 1:
            raise CampaignError(f"[seeds].{sample} must be a positive integer.")
        parsed_seeds[sample] = value
    if len(set(parsed_seeds.values())) != len(parsed_seeds):
        raise CampaignError("Seeds must be unique across enabled samples.")

    values = {
        "nevents": required("nevents"),
        "runs_per_sample": required("runs_per_sample"),
        "cores": required("cores"),
        "ebeam1": required("ebeam1"),
        "ebeam2": required("ebeam2"),
        "parton_pt_min": required("parton_pt_min"),
        "light_parton_dr_min": required("light_parton_dr_min"),
        "bb_dr_min": required("bb_dr_min"),
        "jet_radius": required("jet_radius"),
        "reco_jet_pt_min": required("reco_jet_pt_min"),
    }
    if not all(isinstance(values[name], int) and values[name] > 0 for name in ("nevents", "runs_per_sample", "cores")):
        raise CampaignError("nevents, runs_per_sample, and cores must be positive integers.")
    if not all(isinstance(values[name], (int, float)) and float(values[name]) > 0 for name in ("ebeam1", "ebeam2", "parton_pt_min", "jet_radius", "reco_jet_pt_min")):
        raise CampaignError("Beam energies, pT thresholds, and jet radius must be positive.")
    if not all(isinstance(values[name], (int, float)) and float(values[name]) >= 0 for name in ("light_parton_dr_min", "bb_dr_min")):
        raise CampaignError("Delta-R thresholds cannot be negative.")
    if not isinstance(required("keep_hepmc"), bool):
        raise CampaignError("keep_hepmc must be true or false.")

    serialized = json.dumps(raw, sort_keys=True, separators=(",", ":")).encode()
    return CampaignConfig(
        name=str(required("name")),
        mg5_root=Path(str(required("mg5_root"))).expanduser(),
        output_root=Path(str(required("output_root"))).expanduser(),
        nevents=int(values["nevents"]),
        runs_per_sample=int(values["runs_per_sample"]),
        cores=int(values["cores"]),
        ebeam1=float(values["ebeam1"]),
        ebeam2=float(values["ebeam2"]),
        parton_pt_min=float(values["parton_pt_min"]),
        light_parton_dr_min=float(values["light_parton_dr_min"]),
        bb_dr_min=float(values["bb_dr_min"]),
        jet_radius=float(values["jet_radius"]),
        reco_jet_pt_min=float(values["reco_jet_pt_min"]),
        keep_hepmc=required("keep_hepmc"),
        enabled_samples=tuple(enabled),
        seeds=parsed_seeds,
        config_hash=hashlib.sha256(serialized).hexdigest(),
    )


def run_name(index: int) -> str:
    return f"run_{index:02d}"


def process_directory(config: CampaignConfig, sample: Sample) -> Path:
    return config.output_root / "processes" / sample.directory


def render_build_commands(sample: Sample, process_dir: Path) -> str:
    """Build only the matrix-element directory; card execution happens later."""
    lines = [
        "set automatic_html_opening False --no_save",
        "set notification_center False --no_save",
        "import model sm",
        "define p = g u d s u~ d~ s~",
    ]
    if sample.name == "inclusive":
        lines.append("define j = g u d s c b u~ d~ s~ c~ b~")
    lines.extend((f"generate {sample.process}", f"output {process_dir} -f"))
    return "\n".join(lines) + "\n"


def render_launch_commands(config: CampaignConfig, process_dir: Path, current_run: str) -> str:
    """Launch one process with internal multicore execution, never parallel samples."""
    return "\n".join(
        (
            "set automatic_html_opening False --no_save",
            "set notification_center False --no_save",
            "set run_mode 2 --no_save",
            f"set nb_core {config.cores} --no_save",
            f"launch {process_dir} --name={current_run} -f",
            "",
        )
    )


def _replace_assignment(text: str, name: str, value: float | int) -> str:
    pattern = re.compile(rf"^(?P<indent>\s*)(?P<old>[^!\n=]+?)\s*=\s*{re.escape(name)}\b(?P<tail>.*)$", re.MULTILINE)
    replacement, count = pattern.subn(lambda match: f"{match.group('indent')}{value} = {name}{match.group('tail')}", text, count=1)
    if count != 1:
        raise CampaignError(f"Could not locate exactly one '{name}' setting in run_card.dat.")
    return replacement


def _replace_delphes_module(text: str, module_name: str, radius: float, pt_min: float) -> str:
    pattern = re.compile(rf"(module FastJetFinder {re.escape(module_name)} \{{.*?\n\}})", re.DOTALL)
    match = pattern.search(text)
    if not match:
        raise CampaignError(f"Could not locate Delphes {module_name} module.")
    block = match.group(1)
    block, radius_count = re.subn(r"(?m)^(\s*set ParameterR\s+)\S+", rf"\g<1>{radius}", block, count=1)
    block, pt_count = re.subn(r"(?m)^(\s*set JetPTMin\s+)\S+", rf"\g<1>{pt_min}", block, count=1)
    if radius_count != 1 or pt_count != 1:
        raise CampaignError(f"Could not configure Jet parameters in Delphes {module_name} module.")
    return text[: match.start()] + block + text[match.end() :]


def patch_delphes_card(text: str, radius: float, pt_min: float) -> str:
    """Patch the only two jet collections consumed by the analysis, not FatJets."""
    text = _replace_delphes_module(text, "GenJetFinder", radius, pt_min)
    return _replace_delphes_module(text, "FastJetFinder", radius, pt_min)


def prepare_cards(config: CampaignConfig, sample: Sample, process_dir: Path, seed: int) -> None:
    """Reset generated cards to their defaults and apply the controlled settings."""
    cards = process_dir / "Cards"
    defaults = {
        "run_card.dat": "run_card_default.dat",
        "delphes_card.dat": "delphes_card_default.dat",
        "pythia8_card.dat": "pythia8_card_default.dat",
    }
    for target, default in defaults.items():
        source = cards / default
        if not source.exists():
            raise CampaignError(f"Missing MadGraph card template: {source}")
        shutil.copy2(source, cards / target)

    run_card = cards / "run_card.dat"
    text = run_card.read_text()
    for name, value in (
        ("nevents", config.nevents),
        ("iseed", seed),
        ("ebeam1", config.ebeam1),
        ("ebeam2", config.ebeam2),
    ):
        text = _replace_assignment(text, name, value)
    if sample.heavy_b:
        text = _replace_assignment(text, "ptb", config.parton_pt_min)
        text = _replace_assignment(text, "drbb", config.bb_dr_min)
        text = _replace_assignment(text, "maxjetflavor", 4)
    else:
        text = _replace_assignment(text, "ptj", config.parton_pt_min)
        text = _replace_assignment(text, "drjj", config.light_parton_dr_min)
        text = _replace_assignment(text, "maxjetflavor", 5)
    run_card.write_text(text)

    delphes_card = cards / "delphes_card.dat"
    delphes_card.write_text(patch_delphes_card(delphes_card.read_text(), config.jet_radius, config.reco_jet_pt_min))


def _setting_from_banner(banner: str, name: str) -> str:
    pattern = re.compile(rf"(?m)^\s*(\S+)\s*=\s*{re.escape(name)}\b")
    match = pattern.search(banner)
    if not match:
        raise CampaignError(f"Banner is missing '{name}'.")
    return match.group(1)


def _as_float(value: str, name: str) -> float:
    try:
        return float(value)
    except ValueError as error:
        raise CampaignError(f"Banner value for {name} is not numeric: {value}") from error


def validate_banner(config: CampaignConfig, sample: Sample, banner_path: Path, seed: int) -> dict[str, Any]:
    """Confirm the executed banner, not mutable Cards files, matches the campaign."""
    if not banner_path.exists():
        raise CampaignError(f"Missing MadGraph banner: {banner_path}")
    banner = banner_path.read_text(errors="replace")
    expected_scalars = {
        "nevents": config.nevents,
        "iseed": seed,
        "ebeam1": config.ebeam1,
        "ebeam2": config.ebeam2,
    }
    for name, expected in expected_scalars.items():
        actual = _as_float(_setting_from_banner(banner, name), name)
        if actual != float(expected):
            raise CampaignError(f"Banner {name}={actual} does not match expected {expected}.")
    if f"generate {sample.process}" not in banner:
        raise CampaignError(f"Banner does not contain expected process: {sample.process}")

    cut_name, cut_value = ("ptb", config.parton_pt_min) if sample.heavy_b else ("ptj", config.parton_pt_min)
    actual_cut = _as_float(_setting_from_banner(banner, cut_name), cut_name)
    if actual_cut != cut_value:
        raise CampaignError(f"Banner {cut_name}={actual_cut} does not match expected {cut_value}.")
    dr_name, dr_value = ("drbb", config.bb_dr_min) if sample.heavy_b else ("drjj", config.light_parton_dr_min)
    actual_dr = _as_float(_setting_from_banner(banner, dr_name), dr_name)
    if actual_dr != dr_value:
        raise CampaignError(f"Banner {dr_name}={actual_dr} does not match expected {dr_value}.")

    radius_matches = [float(value) for value in re.findall(r"(?m)^\s*set ParameterR\s+([0-9.eE+-]+)", banner)]
    jet_pt_matches = [float(value) for value in re.findall(r"(?m)^\s*set JetPTMin\s+([0-9.eE+-]+)", banner)]
    if radius_matches.count(config.jet_radius) < 2:
        raise CampaignError(f"Banner does not contain two R={config.jet_radius} jet collections.")
    if jet_pt_matches.count(config.reco_jet_pt_min) < 2:
        raise CampaignError(f"Banner does not contain two JetPTMin={config.reco_jet_pt_min} jet collections.")
    return {"nevents": config.nevents, "seed": seed, "parton_cut": {cut_name: cut_value, dr_name: dr_value}}


def truth_label(flavor: int) -> str | None:
    return TRUTH_LABELS.get(abs(flavor) if abs(flavor) != 21 else 21)


def validate_root(root_path: Path, sample: Sample, pt_min: float, pt_max: float = 20.0) -> dict[str, Any]:
    """Validate all required branches and count usable truth-labelled low-pT jets."""
    try:
        import awkward as ak
        import numpy as np
        import uproot
    except ImportError as error:
        raise CampaignError("ROOT validation needs awkward, numpy, and uproot in the Python environment.") from error
    if not root_path.exists() or root_path.stat().st_size == 0:
        raise CampaignError(f"Missing or empty Delphes ROOT file: {root_path}")
    try:
        tree = uproot.open(root_path)["Delphes"]
        # Uproot resolves Delphes dot aliases even though physical branch names use slashes.
        tree.arrays(REQUIRED_BRANCHES, entry_stop=1, library="ak")
    except Exception as error:
        raise CampaignError(f"ROOT does not satisfy the required Delphes branch contract: {error}") from error

    counts = {label: 0 for label in ("b", "c", "g", "uds", "unknown")}
    jet_total = 0
    low_pt_total = 0
    for arrays in tree.iterate(["Jet.PT", "Jet.Flavor"], step_size="100 MB", library="ak"):
        pt = ak.to_numpy(ak.flatten(arrays["Jet.PT"], axis=None))
        flavor = ak.to_numpy(ak.flatten(arrays["Jet.Flavor"], axis=None))
        jet_total += int(len(pt))
        selected = (pt >= pt_min) & (pt < pt_max)
        low_pt_total += int(np.count_nonzero(selected))
        selected_flavor = np.abs(flavor[selected].astype(np.int64, copy=False))
        counts["b"] += int(np.count_nonzero(selected_flavor == 5))
        counts["c"] += int(np.count_nonzero(selected_flavor == 4))
        counts["g"] += int(np.count_nonzero(selected_flavor == 21))
        counts["uds"] += int(np.count_nonzero(np.isin(selected_flavor, (1, 2, 3))))
        known = counts["b"] + counts["c"] + counts["g"] + counts["uds"]
        counts["unknown"] = low_pt_total - known

    d0_finite = True
    dz_finite = True
    for arrays in tree.iterate(["Track.D0", "Track.DZ"], step_size="100 MB", library="ak"):
        d0 = ak.to_numpy(ak.flatten(arrays["Track.D0"], axis=None))
        dz = ak.to_numpy(ak.flatten(arrays["Track.DZ"], axis=None))
        d0_finite = d0_finite and bool(np.isfinite(d0).all())
        dz_finite = dz_finite and bool(np.isfinite(dz).all())
    if low_pt_total == 0:
        raise CampaignError(f"ROOT has no jets in {pt_min} <= jet_pt < {pt_max} GeV.")
    if sample.target_label and counts[sample.target_label] == 0:
        raise CampaignError(f"ROOT has no low-pT truth-labelled {sample.target_label} jets for {sample.name}.")
    if not d0_finite or not dz_finite:
        raise CampaignError("ROOT contains non-finite Track.D0 or Track.DZ values.")
    return {
        "events": int(tree.num_entries),
        "jets": jet_total,
        "jets_in_fiducial_region": low_pt_total,
        "truth_counts": counts,
        "track_d0_finite": d0_finite,
        "track_dz_finite": dz_finite,
    }


def locate_run_artifacts(run_directory: Path, current_run: str) -> tuple[Path, Path]:
    """Locate the matching tag, including a new tag assigned after a retry."""
    banners = sorted(run_directory.glob(f"{current_run}_tag_*_banner.txt"))
    roots = sorted(run_directory.glob("tag_*_delphes_events.root"))
    pairs: list[tuple[Path, Path]] = []
    for banner in banners:
        match = re.fullmatch(rf"{re.escape(current_run)}_(tag_[^_]+)_banner\.txt", banner.name)
        if not match:
            continue
        root = run_directory / f"{match.group(1)}_delphes_events.root"
        if root in roots:
            pairs.append((banner, root))
    if len(pairs) != 1:
        raise CampaignError(f"Expected exactly one matching banner/ROOT pair in {run_directory}, found {len(pairs)}.")
    return pairs[0]


def _atomic_json_write(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, prefix=f".{path.name}.", delete=False) as temporary:
        json.dump(data, temporary, indent=2, sort_keys=True)
        temporary.write("\n")
        temporary_path = Path(temporary.name)
    temporary_path.replace(path)


def load_manifest(path: Path, config: CampaignConfig, resume: bool) -> dict[str, Any]:
    if not path.exists():
        return {
            "schema_version": 1,
            "campaign": config.name,
            "config_hash": config.config_hash,
            "created_at": int(time.time()),
            "samples": {},
        }
    if not resume:
        raise CampaignError(f"Campaign manifest already exists: {path}. Use --resume after reviewing it.")
    manifest = json.loads(path.read_text())
    if manifest.get("config_hash") != config.config_hash:
        raise CampaignError("Config differs from the existing campaign manifest; create a new campaign directory instead.")
    return manifest


@contextlib.contextmanager
def campaign_lock(path: Path) -> Iterator[None]:
    """Acquire a non-blocking filesystem lock shared by all campaign invocations."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise CampaignError(f"Another campaign runner owns lock {path}.") from error
        handle.seek(0)
        handle.truncate()
        handle.write(f"pid={os.getpid()}\nstarted_at={int(time.time())}\n")
        handle.flush()
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def running_madgraph_processes(mg5_root: Path) -> list[str]:
    """Best-effort protection against a manually started MadGraph generation."""
    running: list[str] = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit() or int(proc.name) == os.getpid():
            continue
        try:
            command = (proc / "cmdline").read_bytes().replace(b"\0", b" ").decode(errors="replace")
        except OSError:
            continue
        if str(mg5_root) in command and any(token in command for token in ("mg5_aMC", "generate_events", "madevent")):
            running.append(f"pid {proc.name}: {command}")
    return running


def ensure_environment(config: CampaignConfig) -> None:
    executable = config.mg5_root / "bin" / "mg5_aMC"
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise CampaignError(f"MadGraph executable is not available: {executable}")
    if any(" " in str(path) for path in (config.mg5_root, config.output_root)):
        raise CampaignError("MadGraph command paths cannot contain spaces.")
    active = running_madgraph_processes(config.mg5_root)
    if active:
        raise CampaignError("Refusing to overlap another MadGraph process:\n" + "\n".join(active))


def _run_mg5(config: CampaignConfig, command_path: Path, log_path: Path) -> None:
    executable = config.mg5_root / "bin" / "mg5_aMC"
    active = running_madgraph_processes(config.mg5_root)
    if active:
        raise CampaignError("Refusing to overlap another MadGraph process:\n" + "\n".join(active))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("wb") as log:
        completed = subprocess.run(
            [str(executable), str(command_path)],
            cwd=config.mg5_root,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        raise CampaignError(f"MadGraph failed; inspect {log_path}")


def _command_path(config: CampaignConfig, sample: Sample, action: str, current_run: str | None = None) -> Path:
    suffix = f"_{current_run}" if current_run else ""
    return config.output_root / "commands" / f"{sample.name}_{action}{suffix}.mg5"


def _write_command(path: Path, command: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(command)


def ensure_process(config: CampaignConfig, sample: Sample, resume: bool) -> Path:
    directory = process_directory(config, sample)
    if directory.exists():
        if not resume:
            raise CampaignError(f"Process directory already exists: {directory}. Use --resume or a new output_root.")
        process_card = directory / "Cards" / "proc_card_mg5.dat"
        if not process_card.exists() or sample.process not in process_card.read_text(errors="replace"):
            raise CampaignError(f"Existing process directory does not match {sample.process}: {directory}")
        return directory
    command_path = _command_path(config, sample, "build")
    _write_command(command_path, render_build_commands(sample, directory))
    _run_mg5(config, command_path, config.output_root / "logs" / f"{sample.name}_build.log")
    if not directory.exists():
        raise CampaignError(f"MadGraph did not create process directory: {directory}")
    return directory


def clean_failed_run(process_dir: Path, current_run: str) -> None:
    target = process_dir / "Events" / current_run
    resolved_events = (process_dir / "Events").resolve()
    if not target.exists():
        return
    if resolved_events not in target.resolve().parents:
        raise CampaignError(f"Refusing to remove unexpected run path: {target}")
    shutil.rmtree(target)


def remove_hepmc(run_directory: Path) -> None:
    for candidate in run_directory.glob("*_pythia8_events.hepmc*"):
        candidate.unlink()


def validate_completed_run(config: CampaignConfig, sample: Sample, state: dict[str, Any]) -> None:
    """Revalidate artifacts before trusting a completed resume entry."""
    validate_banner(config, sample, Path(state["banner"]), state["seed"])
    validate_root(Path(state["root"]), sample, config.reco_jet_pt_min)


def run_campaign(config: CampaignConfig, args: argparse.Namespace) -> dict[str, Any]:
    selected = tuple(args.only) if args.only else config.enabled_samples
    unknown = set(selected) - set(config.enabled_samples)
    if unknown:
        raise CampaignError(f"--only can select only enabled samples: {', '.join(sorted(unknown))}")
    if len(set(selected)) != len(selected):
        raise CampaignError("--only cannot contain duplicate samples.")
    if args.cores is not None:
        if args.cores < 1:
            raise CampaignError("--cores must be positive.")
        config = dataclasses.replace(config, cores=args.cores)

    if args.dry_run:
        return {
            "mode": "dry-run",
            "campaign": config.name,
            "selected_samples": list(selected),
            "runs_per_sample": config.runs_per_sample,
            "cores_per_madgraph_instance": config.cores,
            "sequential_samples": True,
            "eta_behavior": "preserve MadGraph defaults",
            "output_root": str(config.output_root),
        }

    ensure_environment(config)
    config.output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = config.output_root / "campaign_manifest.json"
    with campaign_lock(config.mg5_root / ".btag_campaign.lock"):
        manifest = load_manifest(manifest_path, config, args.resume)
        for sample_name in selected:
            sample = SAMPLES[sample_name]
            process_dir = ensure_process(config, sample, args.resume)
            sample_state = manifest["samples"].setdefault(sample.name, {"process": sample.process, "runs": {}})
            for index in range(1, config.runs_per_sample + 1):
                current_run = run_name(index)
                seed = config.seeds[sample.name] + index - 1
                state = sample_state["runs"].get(current_run, {})
                if state.get("status") == "completed":
                    validate_completed_run(config, sample, state)
                    continue
                if state.get("status") == "failed" and not args.retry_failed:
                    raise CampaignError(f"{sample.name}/{current_run} previously failed. Inspect its log or pass --retry-failed.")
                run_directory = process_dir / "Events" / current_run
                if run_directory.exists() and state.get("status") != "completed":
                    if not args.retry_failed:
                        raise CampaignError(f"Partial directory exists: {run_directory}. Inspect it or use --retry-failed.")
                    clean_failed_run(process_dir, current_run)

                prepare_cards(config, sample, process_dir, seed)
                command_path = _command_path(config, sample, "launch", current_run)
                _write_command(command_path, render_launch_commands(config, process_dir, current_run))
                sample_state["runs"][current_run] = {"status": "running", "seed": seed, "started_at": int(time.time())}
                _atomic_json_write(manifest_path, manifest)
                log_path = config.output_root / "logs" / f"{sample.name}_{current_run}.log"
                try:
                    _run_mg5(config, command_path, log_path)
                    banner_path, root_path = locate_run_artifacts(run_directory, current_run)
                    banner_result = validate_banner(config, sample, banner_path, seed)
                    root_result = validate_root(root_path, sample, config.reco_jet_pt_min)
                    if not config.keep_hepmc:
                        remove_hepmc(run_directory)
                except Exception as error:
                    sample_state["runs"][current_run] = {
                        "status": "failed",
                        "seed": seed,
                        "failed_at": int(time.time()),
                        "error": str(error),
                        "log": str(log_path),
                    }
                    _atomic_json_write(manifest_path, manifest)
                    raise
                sample_state["runs"][current_run] = {
                    "status": "completed",
                    "seed": seed,
                    "completed_at": int(time.time()),
                    "banner": str(banner_path),
                    "root": str(root_path),
                    "log": str(log_path),
                    "banner_validation": banner_result,
                    "root_validation": root_result,
                }
                _atomic_json_write(manifest_path, manifest)
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("lowpt_v1.toml"), help="Campaign TOML configuration.")
    parser.add_argument("--execute", action="store_true", help="Execute MadGraph. Without this flag the command performs a dry-run.")
    parser.add_argument("--dry-run", action="store_true", help="Print the plan without writing or launching MadGraph.")
    parser.add_argument("--resume", action="store_true", help="Resume only runs marked completed in the existing manifest.")
    parser.add_argument("--retry-failed", action="store_true", help="Remove a failed partial run directory and retry it.")
    parser.add_argument("--only", nargs="+", choices=tuple(SAMPLES), help="Limit execution to enabled sample names.")
    parser.add_argument("--cores", type=int, help="Override the config's internal MadGraph core count.")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.execute and args.dry_run:
        parser.error("Use either --execute or --dry-run, not both.")
    if not args.execute:
        args.dry_run = True
    try:
        config = load_config(args.config)
        result = run_campaign(config, args)
    except CampaignError as error:
        print(f"Campaign failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
