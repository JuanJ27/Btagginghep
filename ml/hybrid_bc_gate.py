"""Observable RF gate based on grouped out-of-fold training scores."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold

from ml.hybrid_bc_contract import HybridContractError
from ml.hybrid_bc_features import BC_FEATURES, feature_matrix


def make_gate_rf(seed: int, n_estimators: int) -> RandomForestClassifier:
    return RandomForestClassifier(
        n_estimators=n_estimators, max_depth=8, min_samples_leaf=10,
        class_weight="balanced", random_state=seed, n_jobs=-1,
    )


def grouped_oof_rf_scores(train: pd.DataFrame, seed: int, n_estimators: int, folds: int) -> np.ndarray:
    groups = train["global_event_id"].to_numpy()
    if len(np.unique(groups)) < folds:
        raise HybridContractError(f"RF OOF gate requires at least {folds} unique training events.")
    x, y = feature_matrix(train), train["is_b"].to_numpy(dtype=int)
    if set(y) != {0, 1}:
        raise HybridContractError("RF gate training requires both b and non-b truth classes.")
    scores = np.empty(len(train), dtype=float)
    splitter = GroupKFold(n_splits=folds)
    for fold, (fit_index, holdout_index) in enumerate(splitter.split(x, y, groups)):
        model = make_gate_rf(seed + fold, n_estimators)
        model.fit(x[fit_index], y[fit_index])
        scores[holdout_index] = model.predict_proba(x[holdout_index])[:, 1]
    return scores


def select_rf_threshold(validation_scores: np.ndarray, validation_is_b: np.ndarray, target_b_efficiency: float) -> float:
    if not 0 < target_b_efficiency < 1:
        raise HybridContractError("target_b_efficiency must be between zero and one.")
    b_scores = validation_scores[validation_is_b == 1]
    if not len(b_scores):
        raise HybridContractError("Validation threshold selection requires truth b rows.")
    return float(np.quantile(b_scores, 1 - target_b_efficiency))


def select_band_half_width(validation_scores: np.ndarray, threshold: float, target_fraction: float = 0.20) -> float:
    """Select an observable score-band width using validation score coverage only."""
    candidates = np.array([0.02, 0.05, 0.10, 0.15])
    coverage = np.array([np.mean(np.abs(validation_scores - threshold) <= width) for width in candidates])
    return float(candidates[np.argmin(np.abs(coverage - target_fraction))])


def gate_mask(scores: np.ndarray, threshold: float, half_width: float) -> np.ndarray:
    return np.abs(scores - threshold) <= half_width


def contamination_report(frame: pd.DataFrame, selected: np.ndarray) -> dict:
    labels = frame["sample_label"].to_numpy()
    report = {"gate_rows": int(selected.sum()), "gate_fraction": float(selected.mean())}
    for label in ("b", "c", "g", "uds"):
        label_mask = labels == label
        count = int(label_mask.sum())
        report[label] = {"rows": count, "selected": int((selected & label_mask).sum()), "selected_fraction": None if not count else float(selected[label_mask].mean())}
    report["warning"] = "g/uds rates are truth-only diagnostics; the b-vs-c reranker is never evaluated as c for these jets."
    return report


def fit_final_gate(train: pd.DataFrame, seed: int, n_estimators: int) -> RandomForestClassifier:
    model = make_gate_rf(seed, n_estimators)
    model.fit(feature_matrix(train), train["is_b"].to_numpy(dtype=int))
    return model
