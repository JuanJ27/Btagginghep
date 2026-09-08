"""Dataset and semantic guards for the conditional b-versus-c study."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


SMOKE_TEST = "SMOKE_TEST"
THESIS_RUN_AUTHORITATIVE_LOWPT = "THESIS_RUN_AUTHORITATIVE_LOWPT"
KEY_COLUMNS = ("global_event_id", "root_file", "event_in_file", "jet_rank")
TRUTH_COLUMNS = {"jet_flavor", "sample_label", "source_sample_label", "is_b"}


class HybridContractError(ValueError):
    """Raised when a dataset cannot support the requested hybrid study."""


def load_manifest(dataset_dir: Path) -> dict:
    path = dataset_dir / "dataset_manifest.json"
    if not path.is_file():
        raise HybridContractError(f"Missing dataset manifest: {path}")
    return json.loads(path.read_text())


def validate_mode(manifest: dict, splits: dict[str, pd.DataFrame], mode: str) -> None:
    if mode not in {SMOKE_TEST, THESIS_RUN_AUTHORITATIVE_LOWPT}:
        raise HybridContractError(f"Unknown mode '{mode}'. Use {SMOKE_TEST} or {THESIS_RUN_AUTHORITATIVE_LOWPT}.")
    selected = manifest.get("label_source", {}).get("selected")
    if selected != "jet_flavor":
        raise HybridContractError("Hybrid b-vs-c requires truth-labelled data: manifest label_source.selected must be 'jet_flavor'.")
    for split_name, frame in splits.items():
        missing = (set(KEY_COLUMNS) | {"jet_flavor", "sample_label", "is_b", "jet_pt"}) - set(frame.columns)
        if missing:
            raise HybridContractError(f"{split_name} is missing required columns: {sorted(missing)}")
        if frame.loc[:, list(KEY_COLUMNS)].duplicated().any():
            raise HybridContractError(f"{split_name} has duplicate jet keys; keyed model comparisons would be ambiguous.")
    if mode == THESIS_RUN_AUTHORITATIVE_LOWPT:
        if manifest.get("pt_region") != "lt20":
            raise HybridContractError(
                "THESIS_RUN_AUTHORITATIVE_LOWPT requires manifest pt_region='lt20'; "
                "all-pT or high-pT data are only valid for SMOKE_TEST."
            )
        offending = {name: int((frame["jet_pt"] >= 20).sum()) for name, frame in splits.items()}
        if any(offending.values()):
            raise HybridContractError(f"Authoritative low-pT data contain jet_pt >= 20 entries: {offending}.")
