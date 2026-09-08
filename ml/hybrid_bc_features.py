"""Fixed, low-dimensional physical inputs for the b-versus-c experiment."""
from __future__ import annotations

import numpy as np
import pandas as pd

from ml.hybrid_bc_contract import KEY_COLUMNS, TRUTH_COLUMNS, HybridContractError


# These are reconstructed jet/track observables, not labels, provenance, or identifiers.
BC_FEATURES = ("jet_pt", "n_tracks_cone", "mean_abs_d0", "max_abs_d0")


def validate_feature_frame(frame: pd.DataFrame) -> None:
    forbidden = set(BC_FEATURES) & (set(KEY_COLUMNS) | TRUTH_COLUMNS)
    if forbidden:  # Defensive: this should remain impossible if the constant is edited.
        raise HybridContractError(f"Feature contract includes metadata or truth columns: {sorted(forbidden)}")
    missing = set(BC_FEATURES) - set(frame.columns)
    if missing:
        raise HybridContractError(f"Dataset is missing b-vs-c physical features: {sorted(missing)}")
    matrix = frame.loc[:, BC_FEATURES].to_numpy(dtype=float)
    if not np.isfinite(matrix).all():
        raise HybridContractError("b-vs-c physical features must be finite.")


def bc_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Return truth-matched b/c rows only; g/uds are never mapped to c."""
    result = frame.loc[frame["sample_label"].isin(("b", "c"))].copy()
    if set(result["sample_label"].unique()) != {"b", "c"}:
        raise HybridContractError("The conditional b-vs-c study requires both truth b and truth c rows.")
    return result


def feature_matrix(frame: pd.DataFrame) -> np.ndarray:
    validate_feature_frame(frame)
    return frame.loc[:, BC_FEATURES].to_numpy(dtype=float)
