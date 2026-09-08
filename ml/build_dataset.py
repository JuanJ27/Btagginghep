from __future__ import annotations

import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

try:
    import awkward as ak
    import uproot
except ModuleNotFoundError:
    ak = None
    uproot = None


# ============================================================
# Configuration
# ============================================================

DEFAULT_DATA_DIR = Path("../data")
DEFAULT_OUTPUT_DIR = Path("../outputs")
SOURCE_PROFILES_PATH = Path(__file__).with_name("source_profiles.json")

# If you want CSV instead of parquet, change to "csv"
SAVE_FORMAT = "parquet"   # "parquet" or "csv"

SAMPLE_MAP = {
    "QCD_bbbar": "b",
    "QCD_ccbar": "c",
    "QCD_gg": "g",
    "QCD_uds": "uds",
}

BRANCHES = [
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
]

JET_CONE = 0.4
MAX_JETS_PER_EVENT = None
SEED = 42
SPLIT_RATIOS = {"train": 0.70, "val": 0.15, "test": 0.15}
PT_REGIONS = {
    "all": ("all jet pT values", lambda jet_pt: True),
    "lt20": ("jet_pt < 20", lambda jet_pt: jet_pt < 20),
    "gt20": ("jet_pt > 20", lambda jet_pt: jet_pt > 20),
}
LABEL_SOURCES = ("sample", "jet_flavor")
JET_FLAVOR_MAPPING = {
    "abs(flavor) == 5": "b",
    "abs(flavor) == 4": "c",
    "abs(flavor) in {1, 2, 3}": "uds",
    "flavor == 21": "g",
}
COMPANION_METADATA_COLUMNS = {
    "source_sample_label", "truth_label", "truth_known", "is_b", "jet_flavor",
    "root_file", "event_in_file", "jet_rank", "global_event_id",
}


# ============================================================
# Small helpers
# ============================================================

class SourceProfileError(RuntimeError):
    """Raised when a reproducible source profile cannot be used."""


def load_source_profiles(path: Path = SOURCE_PROFILES_PATH) -> tuple[str, dict[str, dict]]:
    with path.open() as profile_file:
        config = json.load(profile_file)
    if config.get("version") != 1 or not isinstance(config.get("profiles"), dict):
        raise SourceProfileError(f"Unsupported source profile configuration: {path}")
    return config["default_profile"], config["profiles"]

def delta_phi(phi1: np.ndarray, phi2: float) -> np.ndarray:
    dphi = phi1 - phi2
    return (dphi + np.pi) % (2 * np.pi) - np.pi


def delta_r(jet_eta: float, jet_phi: float, trk_eta: np.ndarray, trk_phi: np.ndarray) -> np.ndarray:
    deta = trk_eta - jet_eta
    dphi = delta_phi(trk_phi, jet_phi)
    return np.sqrt(deta**2 + dphi**2)


def safe_mean(x: np.ndarray) -> float:
    return float(np.mean(x)) if len(x) > 0 else 0.0


def safe_std(x: np.ndarray) -> float:
    return float(np.std(x)) if len(x) > 0 else 0.0


def weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    if len(x) == 0 or np.sum(w) <= 0:
        return 0.0
    return float(np.average(x, weights=w))


def first_radius(sorted_dr: np.ndarray, cum_frac: np.ndarray, target: float) -> float:
    mask = cum_frac >= target
    if np.any(mask):
        return float(sorted_dr[np.argmax(mask)])
    return 0.0


def save_dataframe(df: pd.DataFrame, path_stem: Path) -> Path:
    if SAVE_FORMAT == "parquet":
        out_path = path_stem.with_suffix(".parquet")
        df.to_parquet(out_path, index=False)
    elif SAVE_FORMAT == "csv":
        out_path = path_stem.with_suffix(".csv")
        df.to_csv(out_path, index=False)
    else:
        raise ValueError("SAVE_FORMAT must be 'parquet' or 'csv'")
    return out_path


# ============================================================
# Feature extraction
# ============================================================

def compute_jet_features(
    jet_pt: float,
    jet_eta: float,
    jet_phi: float,
    jet_mass: float,
    jet_ncharged: float,
    jet_nneutral: float,
    trk_pt: np.ndarray,
    trk_eta: np.ndarray,
    trk_phi: np.ndarray,
    trk_d0: np.ndarray,
    trk_dz: np.ndarray,
) -> dict:
    dr = delta_r(jet_eta, jet_phi, trk_eta, trk_phi)
    mask = dr < JET_CONE

    trk_pt = trk_pt[mask]
    trk_eta = trk_eta[mask]
    trk_phi = trk_phi[mask]
    trk_d0 = trk_d0[mask]
    trk_dz = trk_dz[mask]
    dr = dr[mask]

    n_tracks = len(trk_pt)

    if n_tracks == 0:
        return {
            "jet_pt": jet_pt,
            "jet_eta": jet_eta,
            "jet_abs_eta": abs(jet_eta),
            "jet_mass": jet_mass,
            "jet_ncharged": jet_ncharged,
            "jet_nneutral": jet_nneutral,
            "n_tracks_cone": 0,
            "track_pt_sum": 0.0,
            "track_pt_sum_over_jet_pt": 0.0,
            "avg_track_pt": 0.0,
            "std_track_pt": 0.0,
            "n_tracks_below_avg_pt": 0,
            "n_tracks_above_avg_pt": 0,
            "max_pt_ratio": 0.0,
            "min_pt_ratio": 0.0,
            "max_dr_ratio": 0.0,
            "min_dr_ratio": 0.0,
            "dr_max_pt": 0.0,
            "dr_min_pt": 0.0,
            "dr_max_dr": 0.0,
            "dr_min_dr": 0.0,
            "pt_difference": 0.0,
            "r50": 0.0,
            "r95": 0.0,
            "mean_abs_d0": 0.0,
            "max_abs_d0": 0.0,
            "std_abs_d0": 0.0,
            "mean_abs_dz": 0.0,
            "max_abs_dz": 0.0,
            "std_abs_dz": 0.0,
            "pt_weighted_abs_d0": 0.0,
            "pt_weighted_abs_dz": 0.0,
        }

    track_pt_sum = float(np.sum(trk_pt))
    avg_track_pt = float(np.mean(trk_pt))
    std_track_pt = float(np.std(trk_pt))

    n_tracks_below_avg_pt = int(np.sum(trk_pt < avg_track_pt))
    n_tracks_above_avg_pt = int(np.sum(trk_pt >= avg_track_pt))

    idx_max_pt = int(np.argmax(trk_pt))
    idx_min_pt = int(np.argmin(trk_pt))
    idx_max_dr = int(np.argmax(dr))
    idx_min_dr = int(np.argmin(dr))

    max_pt_ratio = float(trk_pt[idx_max_pt] / jet_pt) if jet_pt > 0 else 0.0
    min_pt_ratio = float(trk_pt[idx_min_pt] / jet_pt) if jet_pt > 0 else 0.0
    max_dr_ratio = float(trk_pt[idx_max_dr] / jet_pt) if jet_pt > 0 else 0.0
    min_dr_ratio = float(trk_pt[idx_min_dr] / jet_pt) if jet_pt > 0 else 0.0

    dr_max_pt = float(dr[idx_max_pt])
    dr_min_pt = float(dr[idx_min_pt])
    dr_max_dr = float(dr[idx_max_dr])
    dr_min_dr = float(dr[idx_min_dr])

    pt_difference = float(trk_pt[idx_max_pt] - trk_pt[idx_min_pt])

    order = np.argsort(dr)
    sorted_dr = dr[order]
    sorted_pt = trk_pt[order]

    cum_pt = np.cumsum(sorted_pt)
    cum_frac = cum_pt / cum_pt[-1]

    r50 = first_radius(sorted_dr, cum_frac, 0.50)
    r95 = first_radius(sorted_dr, cum_frac, 0.95)

    abs_d0 = np.abs(trk_d0)
    abs_dz = np.abs(trk_dz)

    return {
        "jet_pt": jet_pt,
        "jet_eta": jet_eta,
        "jet_abs_eta": abs(jet_eta),
        "jet_mass": jet_mass,
        "jet_ncharged": jet_ncharged,
        "jet_nneutral": jet_nneutral,
        "n_tracks_cone": n_tracks,
        "track_pt_sum": track_pt_sum,
        "track_pt_sum_over_jet_pt": (track_pt_sum / jet_pt) if jet_pt > 0 else 0.0,
        "avg_track_pt": avg_track_pt,
        "std_track_pt": std_track_pt,
        "n_tracks_below_avg_pt": n_tracks_below_avg_pt,
        "n_tracks_above_avg_pt": n_tracks_above_avg_pt,
        "max_pt_ratio": max_pt_ratio,
        "min_pt_ratio": min_pt_ratio,
        "max_dr_ratio": max_dr_ratio,
        "min_dr_ratio": min_dr_ratio,
        "dr_max_pt": dr_max_pt,
        "dr_min_pt": dr_min_pt,
        "dr_max_dr": dr_max_dr,
        "dr_min_dr": dr_min_dr,
        "pt_difference": pt_difference,
        "r50": r50,
        "r95": r95,
        "mean_abs_d0": safe_mean(abs_d0),
        "max_abs_d0": float(np.max(abs_d0)),
        "std_abs_d0": safe_std(abs_d0),
        "mean_abs_dz": safe_mean(abs_dz),
        "max_abs_dz": float(np.max(abs_dz)),
        "std_abs_dz": safe_std(abs_dz),
        "pt_weighted_abs_d0": weighted_mean(abs_d0, trk_pt),
        "pt_weighted_abs_dz": weighted_mean(abs_dz, trk_pt),
    }


# ============================================================
# File processing
# ============================================================

def collect_all_files(data_dir: Path) -> list[tuple[Path, str]]:
    items = []

    for folder, label in SAMPLE_MAP.items():
        root_files = sorted((data_dir / folder).glob("*.root"))
        if not root_files:
            print(f"Warning: no ROOT files found in {data_dir / folder}")
            continue

        for rf in root_files:
            items.append((rf, label))

    return items


def select_source_files(data_dir: Path, profile_name: str, profiles: dict[str, dict]) -> tuple[dict, list[tuple[Path, str]]]:
    if profile_name not in profiles:
        raise SourceProfileError(f"Unknown source profile '{profile_name}'. Available profiles: {', '.join(sorted(profiles))}.")
    profile = profiles[profile_name]
    if profile.get("discover_all"):
        return profile, collect_all_files(data_dir)

    items = []
    for folder, filenames in profile.get("files", {}).items():
        if folder not in SAMPLE_MAP:
            raise SourceProfileError(f"Source profile '{profile_name}' names unknown sample folder '{folder}'.")
        for filename in filenames:
            root_file = data_dir / folder / filename
            if not root_file.is_file():
                raise SourceProfileError(f"Source profile '{profile_name}' expected file is missing: {root_file}")
            items.append((root_file, SAMPLE_MAP[folder]))
    return profile, items


def preflight_inputs(items: list[tuple[Path, str]], profile_name: str, profile: dict) -> tuple[list[tuple[Path, str]], dict]:
    if uproot is None:
        raise RuntimeError("ROOT preflight requires the uproot package.")

    valid_items = []
    audit = {"source_profile": profile_name, "source_profile_description": profile["description"], "valid_files": [], "invalid_files": []}
    for root_file, label in items:
        record = {"file": str(root_file.resolve()), "label": label}
        try:
            with uproot.open(root_file) as source:
                if "Delphes" not in source:
                    raise ValueError("missing Delphes tree")
                tree = source["Delphes"]
                if tree.num_entries == 0:
                    raise ValueError("Delphes tree has no entries")
                missing = [branch for branch in BRANCHES if branch not in tree]
                if missing:
                    raise ValueError(f"missing required branches: {', '.join(missing)}")
        except Exception as error:
            record["reason"] = str(error)
            audit["invalid_files"].append(record)
            if profile.get("discover_all"):
                print(f"WARNING: skipping invalid ROOT input {root_file}: {record['reason']}")
        else:
            valid_items.append((root_file, label))
            audit["valid_files"].append(record)
    audit["selected_file_count"] = len(items)
    audit["valid_file_count"] = len(valid_items)
    audit["invalid_file_count"] = len(audit["invalid_files"])
    return valid_items, audit


def write_audit_summary(audit: dict, output_dir: Path) -> Path:
    audit_path = output_dir / "input_audit.json"
    with audit_path.open("w") as audit_file:
        json.dump(audit, audit_file, indent=2)
    print(f"Input preflight: {audit['valid_file_count']} valid, {audit['invalid_file_count']} invalid")
    print(" -", audit_path)
    return audit_path


def passes_pt_region(jet_pt: float, pt_region: str) -> bool:
    return PT_REGIONS[pt_region][1](jet_pt)


def map_jet_flavor(flavor: float) -> tuple[str | None, str]:
    """Map the Delphes Jet.Flavor PDG-like parton code without coercion."""
    if not np.isfinite(flavor):
        return None, "non_finite"
    if flavor == 0:
        return None, "zero"
    if abs(flavor) == 5:
        return "b", "matched"
    if abs(flavor) == 4:
        return "c", "matched"
    if abs(flavor) in {1, 2, 3}:
        return "uds", "matched"
    if flavor == 21:
        return "g", "matched"
    return None, "unsupported"


def raw_flavor_key(flavor: float) -> str:
    if not np.isfinite(flavor):
        return str(flavor).lower()
    return str(int(flavor)) if flavor.is_integer() else str(flavor)


def derive_jet_labels(source_sample_label: str, jet_flavor: float, label_source: str) -> dict | None:
    truth_label, _ = map_jet_flavor(jet_flavor)
    if label_source == "jet_flavor" and truth_label is None:
        return None
    sample_label = source_sample_label if label_source == "sample" else truth_label
    return {
        "source_sample_label": source_sample_label,
        "sample_label": sample_label,
        "is_b": int(sample_label == "b"),
        "jet_flavor": jet_flavor,
    }


def derive_companion_labels(source_sample_label: str, jet_flavor: float) -> dict:
    """Return truth metadata for every selected reconstructed jet."""
    truth_label, _ = map_jet_flavor(jet_flavor)
    return {
        "source_sample_label": source_sample_label,
        "truth_label": truth_label if truth_label is not None else pd.NA,
        "truth_known": truth_label is not None,
        "is_b": int(truth_label == "b") if truth_label is not None else pd.NA,
        "jet_flavor": jet_flavor,
    }


def empty_truth_summary() -> dict:
    return {"matched_jets": 0, "unmatched_jets": 0, "unmatched_by_reason_and_raw_flavor": {}}


def record_truth_label(summary: dict, jet_flavor: float) -> None:
    _, reason = map_jet_flavor(jet_flavor)
    if reason == "matched":
        summary["matched_jets"] += 1
        return
    summary["unmatched_jets"] += 1
    counts = summary["unmatched_by_reason_and_raw_flavor"].setdefault(reason, {})
    raw_value = raw_flavor_key(jet_flavor)
    counts[raw_value] = counts.get(raw_value, 0) + 1


def process_one_file(
    root_file: Path,
    source_sample_label: str,
    pt_region: str,
    label_source: str,
    max_jets_per_event: int | None = MAX_JETS_PER_EVENT,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    if ak is None or uproot is None:
        raise RuntimeError("ROOT processing requires the awkward and uproot packages.")
    supervised_rows = []
    companion_rows = []
    truth_summary = empty_truth_summary()
    event_counter = 0

    for batch in uproot.iterate(f"{root_file}:Delphes", BRANCHES, step_size="100 MB", library="ak"):
        n_events = len(batch["Jet.PT"])

        for ievt in range(n_events):
            jet_pt = np.asarray(ak.to_numpy(batch["Jet.PT"][ievt]), dtype=np.float32)
            jet_flavor = np.asarray(ak.to_numpy(batch["Jet.Flavor"][ievt]), dtype=np.float64)
            jet_eta = np.asarray(ak.to_numpy(batch["Jet.Eta"][ievt]), dtype=np.float32)
            jet_phi = np.asarray(ak.to_numpy(batch["Jet.Phi"][ievt]), dtype=np.float32)
            jet_mass = np.asarray(ak.to_numpy(batch["Jet.Mass"][ievt]), dtype=np.float32)
            jet_ncharged = np.asarray(ak.to_numpy(batch["Jet.NCharged"][ievt]), dtype=np.float32)
            jet_nneutral = np.asarray(ak.to_numpy(batch["Jet.NNeutrals"][ievt]), dtype=np.float32)

            trk_pt = np.asarray(ak.to_numpy(batch["Track.PT"][ievt]), dtype=np.float32)
            trk_eta = np.asarray(ak.to_numpy(batch["Track.Eta"][ievt]), dtype=np.float32)
            trk_phi = np.asarray(ak.to_numpy(batch["Track.Phi"][ievt]), dtype=np.float32)
            trk_d0 = np.asarray(ak.to_numpy(batch["Track.D0"][ievt]), dtype=np.float32)
            trk_dz = np.asarray(ak.to_numpy(batch["Track.DZ"][ievt]), dtype=np.float32)

            n_jets = len(jet_pt) if max_jets_per_event is None else min(len(jet_pt), max_jets_per_event)

            for j in range(n_jets):
                if not passes_pt_region(float(jet_pt[j]), pt_region):
                    continue
                flavor = float(jet_flavor[j])
                record_truth_label(truth_summary, flavor)
                feats = compute_jet_features(
                    jet_pt=float(jet_pt[j]),
                    jet_eta=float(jet_eta[j]),
                    jet_phi=float(jet_phi[j]),
                    jet_mass=float(jet_mass[j]),
                    jet_ncharged=float(jet_ncharged[j]),
                    jet_nneutral=float(jet_nneutral[j]),
                    trk_pt=trk_pt,
                    trk_eta=trk_eta,
                    trk_phi=trk_phi,
                    trk_d0=trk_d0,
                    trk_dz=trk_dz,
                )

                feats["root_file"] = str(root_file)
                feats["event_in_file"] = event_counter
                feats["jet_rank"] = j
                feats["global_event_id"] = f"{root_file}::evt::{event_counter}"
                companion_feats = feats.copy()
                companion_feats.update(derive_companion_labels(source_sample_label, flavor))
                companion_rows.append(companion_feats)

                labels = derive_jet_labels(source_sample_label, flavor, label_source)
                if labels is not None:
                    supervised_feats = feats.copy()
                    supervised_feats.update(labels)
                    supervised_rows.append(supervised_feats)

            event_counter += 1

    companion_df = pd.DataFrame(companion_rows)
    if not companion_df.empty:
        companion_df["truth_label"] = companion_df["truth_label"].astype("string")
        companion_df["is_b"] = companion_df["is_b"].astype("Int64")
        companion_df["truth_known"] = companion_df["truth_known"].astype(bool)
    return pd.DataFrame(supervised_rows), companion_df, truth_summary


# ============================================================
# Splitting
# ============================================================

class InsufficientEventsForSplitError(ValueError):
    """Raised when an event-grouped stratified split cannot retain every label."""


class TruthSourceCoverageError(RuntimeError):
    """Raised when a selected source stratum has no matched truth jets."""


def validate_split_ratios(ratios: dict[str, float]) -> dict[str, float]:
    required = {"train", "val", "test"}
    if set(ratios) != required:
        raise ValueError(f"Split ratios must define exactly {sorted(required)}.")
    if any(ratio <= 0 for ratio in ratios.values()):
        raise ValueError("All split ratios must be positive.")
    if not np.isclose(sum(ratios.values()), 1.0):
        raise ValueError("Split ratios must sum to 1.")
    return {name: float(ratios[name]) for name in ("train", "val", "test")}


def validate_event_split_preflight(event_df: pd.DataFrame, ratios: dict[str, float], stratum_column: str, stratum_name: str) -> None:
    """Ensure each label can appear in every stratified split after grouping."""
    minimum_events = max(int(np.ceil(2 / ratio)) for ratio in ratios.values())
    counts = event_df[stratum_column].value_counts().sort_index()
    deficient = {label: int(count) for label, count in counts.items() if count < minimum_events}
    if deficient:
        details = ", ".join(f"{label}={count}" for label, count in deficient.items())
        raise InsufficientEventsForSplitError(
            f"Insufficient unique events per {stratum_name} for the configured "
            f"event-grouped train/validation/test split (need at least {minimum_events} each): {details}."
        )


def split_by_event(full_df: pd.DataFrame, ratios: dict[str, float] = SPLIT_RATIOS, stratum_column: str = "sample_label"):
    ratios = validate_split_ratios(ratios)
    event_df = full_df[["global_event_id", stratum_column]].drop_duplicates().reset_index(drop=True)
    stratum_name = "sample label" if stratum_column == "sample_label" else stratum_column.replace("_", " ")
    validate_event_split_preflight(event_df, ratios, stratum_column, stratum_name)
    temporary_ratio = ratios["val"] + ratios["test"]

    # First split: train vs temp
    train_events, temp_events = train_test_split(
        event_df,
        test_size=temporary_ratio,
        random_state=SEED,
        stratify=event_df[stratum_column],
    )

    # Second split: val vs test
    val_events, test_events = train_test_split(
        temp_events,
        test_size=ratios["test"] / temporary_ratio,
        random_state=SEED,
        stratify=temp_events[stratum_column],
    )

    train_df = full_df[full_df["global_event_id"].isin(train_events["global_event_id"])].copy()
    val_df = full_df[full_df["global_event_id"].isin(val_events["global_event_id"])].copy()
    test_df = full_df[full_df["global_event_id"].isin(test_events["global_event_id"])].copy()

    return train_df, val_df, test_df, train_events, val_events, test_events


def split_supervised_and_companions(
    supervised_df: pd.DataFrame,
    companion_df: pd.DataFrame,
    ratios: dict[str, float] = SPLIT_RATIOS,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Assign event splits once from all selected reconstructed jets."""
    companion_train, companion_val, companion_test, train_events, val_events, test_events = split_by_event(
        companion_df, ratios, "source_sample_label"
    )
    split_ids = (train_events["global_event_id"], val_events["global_event_id"], test_events["global_event_id"])
    train_df, val_df, test_df = (
        supervised_df[supervised_df["global_event_id"].isin(event_ids)].copy()
        for event_ids in split_ids
    )
    return (
        train_df, val_df, test_df, companion_train, companion_val, companion_test,
        train_events, val_events, test_events,
    )


# ============================================================
# Main
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build event-grouped b-tagging datasets from Delphes ROOT files.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR, help="Directory containing QCD sample folders.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Directory for dataset files and manifest.")
    default_profile, profiles = load_source_profiles()
    parser.add_argument("--source-profile", choices=sorted(profiles), default=default_profile, help="Input provenance profile; independent from --pt-region. Default: discover_all.")
    parser.add_argument("--pt-region", choices=PT_REGIONS, default="all", help="Jet pT selection: all, lt20 (< 20), or gt20 (> 20).")
    parser.add_argument("--label-source", choices=LABEL_SOURCES, default="sample", help="Target label source: sample (folder provenance) or jet_flavor (Delphes Jet.Flavor truth).")
    parser.add_argument("--max-jets-per-event", type=int, default=MAX_JETS_PER_EVENT, help="Optional maximum selected jets per event; default: no cap.")
    parser.add_argument("--audit-only", action="store_true", help="Validate selected ROOT inputs, write input_audit.json, and exit before jet processing.")
    return parser.parse_args()


def class_counts(df: pd.DataFrame) -> dict[str, int]:
    return {label: int(count) for label, count in df["sample_label"].value_counts().sort_index().items()}


def source_target_contingency(df: pd.DataFrame) -> dict[str, dict[str, int]]:
    table = pd.crosstab(df["source_sample_label"], df["sample_label"])
    return {source: {target: int(count) for target, count in row.items() if count} for source, row in table.to_dict(orient="index").items()}


def companion_counts(df: pd.DataFrame) -> dict[str, int]:
    known = int(df["truth_known"].sum())
    return {"total_jets": int(len(df)), "truth_known_jets": known, "unknown_truth_jets": int(len(df) - known)}


def validate_binary_targets(df: pd.DataFrame, label_source: str) -> None:
    present_labels = set(df["sample_label"]) if not df.empty else set()
    missing = []
    if "b" not in present_labels:
        missing.append("b")
    if not (present_labels - {"b"}):
        missing.append("non-b")
    if missing:
        prefix = "Truth filtering left no" if label_source == "jet_flavor" else "Dataset contains no"
        raise RuntimeError(f"{prefix} required target jets: {', '.join(missing)}.")


def validate_truth_source_coverage(audit: dict, df: pd.DataFrame) -> None:
    expected_sources = {record["label"] for record in audit["valid_files"]}
    matched_sources = set(df["source_sample_label"])
    missing_sources = sorted(expected_sources - matched_sources)
    if missing_sources:
        raise TruthSourceCoverageError(
            "Selected source strata have zero matched truth jets: "
            f"{', '.join(missing_sources)}."
        )


def run(args: argparse.Namespace) -> None:
    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    label_source = getattr(args, "label_source", "sample")
    max_jets_per_event = getattr(args, "max_jets_per_event", MAX_JETS_PER_EVENT)
    if max_jets_per_event is not None and max_jets_per_event <= 0:
        raise ValueError("max_jets_per_event must be positive when set.")
    output_dir.mkdir(parents=True, exist_ok=True)
    _, profiles = load_source_profiles()
    profile = profiles.get(args.source_profile)
    selected_items = []
    try:
        profile, selected_items = select_source_files(data_dir, args.source_profile, profiles)
        items, audit = preflight_inputs(selected_items, args.source_profile, profile)
    except Exception as error:
        selected_files = [{"file": str(root_file.resolve()), "label": label} for root_file, label in selected_items]
        missing_files = []
        if profile and not profile.get("discover_all"):
            for folder, filenames in profile.get("files", {}).items():
                for filename in filenames:
                    root_file = data_dir / folder / filename
                    if not root_file.is_file():
                        missing_files.append({"file": str(root_file.resolve()), "label": SAMPLE_MAP.get(folder), "reason": str(error)})
        audit = {
            "source_profile": args.source_profile,
            "source_profile_description": profile.get("description") if profile else None,
            "selected_files": selected_files,
            "missing_files": missing_files,
            "valid_files": [],
            "invalid_files": missing_files,
            "selected_file_count": len(selected_files) + len(missing_files),
            "valid_file_count": 0,
            "invalid_file_count": len(missing_files),
            "status": "failed",
            "error": {"type": type(error).__name__, "message": str(error)},
        }
        write_audit_summary(audit, output_dir)
        raise
    strict_error = None
    if audit["invalid_file_count"] and not profile.get("discover_all"):
        first_invalid = audit["invalid_files"][0]
        strict_error = (
            f"Source profile '{args.source_profile}' has invalid file "
            f"{first_invalid['file']}: {first_invalid['reason']}"
        )
        audit["status"] = "failed"
        audit["error"] = {"type": "SourceProfileError", "message": strict_error}
    else:
        audit["status"] = "success"
        audit["error"] = None
    write_audit_summary(audit, output_dir)

    if strict_error:
        raise SourceProfileError(strict_error)
    if args.audit_only:
        return
    if len(items) == 0:
        raise RuntimeError(f"No valid ROOT files found under {data_dir} for source profile '{args.source_profile}'.")

    all_dfs = []
    all_companion_dfs = []
    truth_summary = empty_truth_summary()

    for root_file, label in items:
        print(f"Processing {root_file.name} as {label}")
        df, companion_df, file_truth_summary = process_one_file(
            root_file, label, args.pt_region, label_source, max_jets_per_event
        )
        truth_summary["matched_jets"] += file_truth_summary["matched_jets"]
        truth_summary["unmatched_jets"] += file_truth_summary["unmatched_jets"]
        for reason, raw_counts in file_truth_summary["unmatched_by_reason_and_raw_flavor"].items():
            counts = truth_summary["unmatched_by_reason_and_raw_flavor"].setdefault(reason, {})
            for raw_value, count in raw_counts.items():
                counts[raw_value] = counts.get(raw_value, 0) + count
        print(f"  -> {len(df)} supervised jets, {len(companion_df)} companion jets")
        all_dfs.append(df)
        all_companion_dfs.append(companion_df)

    full_df = pd.concat(all_dfs, ignore_index=True)
    full_companion_df = pd.concat(all_companion_dfs, ignore_index=True)

    if full_df.empty and label_source == "sample":
        raise RuntimeError(f"No jets passed pT region '{args.pt_region}' ({PT_REGIONS[args.pt_region][0]}).")
    if label_source == "jet_flavor":
        validate_truth_source_coverage(audit, full_df)
    validate_binary_targets(full_df, label_source)

    print("\nTotal jets in full dataset:", len(full_df))
    print("\nJets per class:")
    print(full_df["sample_label"].value_counts())

    split_ratios = validate_split_ratios(SPLIT_RATIOS)
    split_stratum = "source_sample_label"
    (
        train_df, val_df, test_df, companion_train_df, companion_val_df, companion_test_df,
        train_events, val_events, test_events,
    ) = split_supervised_and_companions(full_df, full_companion_df, split_ratios)

    print("\nSplit summary:")
    print(f"Train jets: {len(train_df)}")
    print(f"Val jets:   {len(val_df)}")
    print(f"Test jets:  {len(test_df)}")

    print("\nTrain class counts:")
    print(train_df["sample_label"].value_counts())

    print("\nVal class counts:")
    print(val_df["sample_label"].value_counts())

    print("\nTest class counts:")
    print(test_df["sample_label"].value_counts())

    train_path = save_dataframe(train_df, output_dir / "train")
    val_path = save_dataframe(val_df, output_dir / "val")
    test_path = save_dataframe(test_df, output_dir / "test")
    companion_train_path = save_dataframe(companion_train_df, output_dir / "train_companion")
    companion_val_path = save_dataframe(companion_val_df, output_dir / "val_companion")
    companion_test_path = save_dataframe(companion_test_df, output_dir / "test_companion")

    source_inventory = [str(root_file.resolve()) for root_file, _ in items]
    metadata_columns = ["source_sample_label", "sample_label", "is_b", "jet_flavor", "root_file", "event_in_file", "jet_rank", "global_event_id"]
    feature_columns = [column for column in full_df.columns if column not in metadata_columns]
    companion_metadata_columns = sorted(COMPANION_METADATA_COLUMNS)
    companion_feature_columns = [column for column in full_companion_df.columns if column not in COMPANION_METADATA_COLUMNS]

    manifest = {
        "data_dir": str(data_dir),
        "output_dir": str(output_dir),
        "save_format": SAVE_FORMAT,
        "source_profile": args.source_profile,
        "source_profile_description": profile["description"],
        "pt_region": args.pt_region,
        "pt_expression": PT_REGIONS[args.pt_region][0],
        "label_source": {
            "version": 1,
            "selected": label_source,
            "jet_flavor_mapping": JET_FLAVOR_MAPPING,
            "truth_summary": truth_summary,
        },
        "source_root_inventory": source_inventory,
        "input_audit": audit,
        "files": [{"file": str(f.resolve()), "label": label} for f, label in items],
        "generator_settings": {
            "seed": SEED,
            "jet_cone": JET_CONE,
            "max_jets_per_event": max_jets_per_event,
            "max_jets_policy": "unlimited" if max_jets_per_event is None else "explicit_cap",
            "branches": BRANCHES,
        },
        "split": {
            "method": "event_grouped_stratified",
            "stratification_basis": split_stratum,
            "assignment_scope": "all selected reconstructed-jet events; reused by supervised and companion splits",
            "seed": SEED,
            "ratios": split_ratios,
            "applied_event_ratios": {
                "train": len(train_events) / len(full_companion_df[["global_event_id"]].drop_duplicates()),
                "val": len(val_events) / len(full_companion_df[["global_event_id"]].drop_duplicates()),
                "test": len(test_events) / len(full_companion_df[["global_event_id"]].drop_duplicates()),
            },
        },
        "schema": {"columns": full_df.columns.tolist(), "feature_columns": feature_columns, "metadata_columns": metadata_columns},
        "companion": {
            "definition": "All reconstructed jets that pass the configured pT selection; no per-event jet cap is applied." if max_jets_per_event is None else "All reconstructed jets that pass the configured pT selection, limited by the explicit per-event jet cap.",
            "label_contract": "truth_label and is_b are NA when Jet.Flavor is unmapped; truth_known identifies recognised b/c/uds/g truth.",
            "schema": {
                "columns": full_companion_df.columns.tolist(),
                "scoring_feature_columns": companion_feature_columns,
                "metadata_and_label_columns": companion_metadata_columns,
            },
            "counts": {
                "full": companion_counts(full_companion_df),
                "train": companion_counts(companion_train_df),
                "val": companion_counts(companion_val_df),
                "test": companion_counts(companion_test_df),
            },
            "outputs": {
                "train": str(companion_train_path), "val": str(companion_val_path), "test": str(companion_test_path),
            },
        },
        "sample_counts": {"full": class_counts(full_df), "train": class_counts(train_df), "val": class_counts(val_df), "test": class_counts(test_df)},
        "target_class_counts_by_split": {"full": class_counts(full_df), "train": class_counts(train_df), "val": class_counts(val_df), "test": class_counts(test_df)},
        "source_target_contingency": source_target_contingency(full_df),
        "binary_class_counts": {
            "full": {str(key): int(value) for key, value in full_df["is_b"].value_counts().sort_index().items()},
            "train": {str(key): int(value) for key, value in train_df["is_b"].value_counts().sort_index().items()},
            "val": {str(key): int(value) for key, value in val_df["is_b"].value_counts().sort_index().items()},
            "test": {str(key): int(value) for key, value in test_df["is_b"].value_counts().sort_index().items()},
        },
        "n_full_jets": int(len(full_df)),
        "n_train_jets": int(len(train_df)),
        "n_val_jets": int(len(val_df)),
        "n_test_jets": int(len(test_df)),
        "n_train_events": int(len(train_events)),
        "n_val_events": int(len(val_events)),
        "n_test_events": int(len(test_events)),
        "train_output": str(train_path),
        "val_output": str(val_path),
        "test_output": str(test_path),
    }

    with open(output_dir / "dataset_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    print("\nSaved:")
    print(" -", train_path)
    print(" -", val_path)
    print(" -", test_path)
    print(" -", companion_train_path)
    print(" -", companion_val_path)
    print(" -", companion_test_path)
    print(" -", output_dir / "dataset_manifest.json")
    print("\nDone.")


def main():
    run(parse_args())


if __name__ == "__main__":
    main()
