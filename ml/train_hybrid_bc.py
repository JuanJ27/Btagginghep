"""Run the bounded, truth-labelled conditional b-versus-c hybrid experiment."""
from __future__ import annotations

import argparse
import hashlib
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

from ml.hybrid_bc_contract import KEY_COLUMNS, SMOKE_TEST, THESIS_RUN_AUTHORITATIVE_LOWPT, HybridContractError, load_manifest, validate_mode
from ml.hybrid_bc_features import BC_FEATURES, bc_rows, feature_matrix, validate_feature_frame
from ml.hybrid_bc_gate import contamination_report, fit_final_gate, gate_mask, grouped_oof_rf_scores, select_band_half_width, select_rf_threshold
from ml.hybrid_bc_models import ANGLE_VQC_PRESETS, angle_vqc_preset, fit_angle_vqc, fit_controls, reject_amplitude_encoding


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Conditional b-vs-c reranker study; it never emits b-vs-all probabilities.")
    parser.add_argument("--dataset-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mode", choices=(SMOKE_TEST, THESIS_RUN_AUTHORITATIVE_LOWPT), default=SMOKE_TEST)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gate-folds", type=int, default=5)
    parser.add_argument("--gate-estimators", type=int, default=100)
    parser.add_argument("--target-b-efficiency", type=float, default=0.70)
    parser.add_argument("--max-gate-train", type=int, default=5000, help="Deterministic bounded smoke-study training rows; use 0 for all rows.")
    parser.add_argument("--max-bc-eval", type=int, default=2000, help="Common b/c validation/test rows per control; use 0 for all rows.")
    parser.add_argument("--quantum", choices=("disabled", "angle", "amplitude"), default="disabled")
    parser.add_argument("--quantum-maxiter", type=int, default=50, help="Bounded local COBYLA iterations for --quantum angle.")
    parser.add_argument("--quantum-preset", choices=tuple(sorted(ANGLE_VQC_PRESETS)), default="zz1_real1_linear", help="Controlled four-feature angle-VQC architecture.")
    parser.add_argument("--quantum-shots", type=int, default=1024, help="Seeded local-Aer shots per circuit evaluation.")
    return parser.parse_args()


def load_split(dataset_dir: Path, name: str) -> pd.DataFrame:
    path = dataset_dir / f"{name}.parquet"
    if not path.is_file():
        raise FileNotFoundError(f"Missing Parquet split: {path}")
    return pd.read_parquet(path)


def bounded_stratified(frame: pd.DataFrame, maximum: int, seed: int) -> pd.DataFrame:
    if maximum <= 0 or len(frame) <= maximum:
        return frame.copy()
    pieces = []
    for label, part in frame.groupby("sample_label", sort=True):
        n = max(1, round(maximum * len(part) / len(frame)))
        pieces.append(part.sample(n=min(n, len(part)), random_state=seed + (0 if label == "b" else 1)))
    return pd.concat(pieces).sort_index().head(maximum).copy()


def keyed_output(frame: pd.DataFrame, scores: dict[str, np.ndarray]) -> pd.DataFrame:
    validate_bc_rows(frame, "test")
    output = frame.loc[:, [*KEY_COLUMNS, "sample_label", "jet_flavor"]].copy()
    for name, values in scores.items():
        validate_bc_scores(frame, values, name, "test")
        output[f"score_bc_b_{name}"] = values
    output["score_semantics"] = "P(b | truth-selected b/c study); never a b-vs-all score and never a g/uds classification."
    return output


def validate_bc_rows(frame: pd.DataFrame, split_name: str) -> None:
    """Require the exact, unambiguous truth-selected b/c rows used for each comparison."""
    missing = (set(KEY_COLUMNS) | {"sample_label"}) - set(frame.columns)
    if missing:
        raise HybridContractError(f"{split_name} b/c rows are missing required columns: {sorted(missing)}")
    if frame.loc[:, list(KEY_COLUMNS)].duplicated().any():
        raise HybridContractError(f"{split_name} b/c rows have duplicate jet keys.")
    labels = set(frame["sample_label"])
    if labels != {"b", "c"}:
        raise HybridContractError(f"{split_name} comparison rows must contain only b and c labels, with both classes present; got {sorted(labels)}.")


def validate_bc_scores(frame: pd.DataFrame, scores: np.ndarray, model_name: str, split_name: str) -> np.ndarray:
    """Validate conditional P(b | b/c) scores against their keyed b/c row set."""
    values = np.asarray(scores, dtype=float)
    if values.ndim != 1 or len(values) != len(frame):
        raise HybridContractError(f"{model_name} {split_name} scores must have exactly one value per keyed b/c row.")
    if not np.isfinite(values).all() or (values < 0).any() or (values > 1).any():
        raise HybridContractError(f"{model_name} {split_name} scores must be finite conditional probabilities in [0, 1].")
    return values


def keyed_row_fingerprint(frame: pd.DataFrame, split_name: str) -> str:
    """Fingerprint the ordered key set so independent preset reports are auditable."""
    validate_bc_rows(frame, split_name)
    payload = frame.loc[:, list(KEY_COLUMNS)].to_csv(index=False, lineterminator="\n").encode()
    return hashlib.sha256(payload).hexdigest()


def select_bc_threshold(validation_scores: np.ndarray, validation_labels: np.ndarray, target_b_efficiency: float) -> float:
    """Choose the highest validation cutoff reaching target b efficiency without splitting score ties."""
    if not 0 < target_b_efficiency <= 1:
        raise HybridContractError("target b efficiency must be in (0, 1].")
    labels = np.asarray(validation_labels, dtype=int)
    scores = np.asarray(validation_scores, dtype=float)
    if labels.ndim != 1 or len(labels) != len(scores) or set(labels) != {0, 1}:
        raise HybridContractError("Threshold selection requires aligned validation b/c labels containing both b (1) and c (0).")
    if not np.isfinite(scores).all() or (scores < 0).any() or (scores > 1).any():
        raise HybridContractError("Threshold selection requires finite conditional probabilities in [0, 1].")
    required_b = target_b_efficiency * labels.sum()
    selected_b = 0
    for score in np.unique(scores)[::-1]:
        selected_b += int(labels[scores == score].sum())
        if selected_b >= required_b:
            return float(score)
    raise HybridContractError("Validation threshold selection could not reach the requested b efficiency.")


def evaluate_bc_model(
    validation_frame: pd.DataFrame,
    validation_scores: np.ndarray,
    test_frame: pd.DataFrame,
    test_scores: np.ndarray,
    target_b_efficiency: float,
    model_name: str,
) -> dict:
    validate_bc_rows(validation_frame, "validation")
    validate_bc_rows(test_frame, "test")
    validation_values = validate_bc_scores(validation_frame, validation_scores, model_name, "validation")
    test_values = validate_bc_scores(test_frame, test_scores, model_name, "test")
    validation_labels = (validation_frame["sample_label"] == "b").to_numpy(int)
    test_labels = (test_frame["sample_label"] == "b").to_numpy(int)
    threshold = select_bc_threshold(validation_values, validation_labels, target_b_efficiency)
    selected = test_values >= threshold
    true_positive = int(np.count_nonzero(selected & (test_labels == 1)))
    false_negative = int(np.count_nonzero(~selected & (test_labels == 1)))
    false_positive = int(np.count_nonzero(selected & (test_labels == 0)))
    true_negative = int(np.count_nonzero(~selected & (test_labels == 0)))
    selected_count = int(selected.sum())
    return {
        "score_semantics": "P(b | truth-selected b/c study); never a b-vs-all score.",
        "validation_roc_auc": float(roc_auc_score(validation_labels, validation_values)),
        "test_roc_auc": float(roc_auc_score(test_labels, test_values)),
        "frozen_validation_threshold": threshold,
        "target_validation_b_efficiency": target_b_efficiency,
        "test_b_efficiency": true_positive / int(test_labels.sum()),
        "test_c_mistag": false_positive / int((test_labels == 0).sum()),
        "test_bc_purity": true_positive / selected_count if selected_count else 0.0,
        "test_selected_count": selected_count,
        "test_confusion_counts": {"b_selected": true_positive, "b_rejected": false_negative, "c_selected": false_positive, "c_rejected": true_negative},
    }


def write_comparison_figure(output_path: Path, test_frame: pd.DataFrame, test_predictions: dict[str, np.ndarray], comparison: dict[str, dict]) -> None:
    """Write a conditional b/c ROC comparison using the same keyed test rows for every model."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = (test_frame["sample_label"] == "b").to_numpy(int)
    figure, axis = plt.subplots(figsize=(7, 5))
    for name, scores in test_predictions.items():
        fpr, tpr, _ = roc_curve(labels, scores)
        axis.plot(fpr, tpr, label=f"{name} (AUC={comparison[name]['test_roc_auc']:.3f})")
    axis.plot([0, 1], [0, 1], "k--", linewidth=1, label="chance")
    axis.set(xlabel="c mistag rate", ylabel="b efficiency", title="Conditional b-vs-c comparison")
    axis.legend(loc="lower right")
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def run(args: argparse.Namespace) -> dict:
    dataset_dir, output_dir = args.dataset_dir.resolve(), args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    splits = {name: load_split(dataset_dir, name) for name in ("train", "val", "test")}
    manifest = load_manifest(dataset_dir)
    validate_mode(manifest, splits, args.mode)
    for frame in splits.values():
        validate_feature_frame(frame)
    gate_train = bounded_stratified(splits["train"], args.max_gate_train, args.seed)
    oof_scores = grouped_oof_rf_scores(gate_train, args.seed, args.gate_estimators, args.gate_folds)
    final_gate = fit_final_gate(gate_train, args.seed, args.gate_estimators)
    val_scores = final_gate.predict_proba(feature_matrix(splits["val"]))[:, 1]
    test_scores = final_gate.predict_proba(feature_matrix(splits["test"]))[:, 1]
    threshold = select_rf_threshold(val_scores, splits["val"]["is_b"].to_numpy(int), args.target_b_efficiency)
    half_width = select_band_half_width(val_scores, threshold)
    train_selected = gate_train.loc[gate_mask(oof_scores, threshold, half_width)].copy()
    val_selected = splits["val"].loc[gate_mask(val_scores, threshold, half_width)].copy()
    test_selected = splits["test"].loc[gate_mask(test_scores, threshold, half_width)].copy()
    bc_train = bounded_stratified(bc_rows(train_selected), args.max_gate_train, args.seed)
    bc_val = bounded_stratified(bc_rows(val_selected), args.max_bc_eval, args.seed)
    bc_test = bounded_stratified(bc_rows(test_selected), args.max_bc_eval, args.seed)
    validate_bc_rows(bc_train, "train")
    validate_bc_rows(bc_val, "validation")
    validate_bc_rows(bc_test, "test")
    models = fit_controls(feature_matrix(bc_train), (bc_train["sample_label"] == "b").to_numpy(int), args.seed)
    validation_predictions = {name: model.predict_proba(feature_matrix(bc_val))[:, 1] for name, model in models.items()}
    test_predictions = {name: model.predict_proba(feature_matrix(bc_test))[:, 1] for name, model in models.items()}
    if args.quantum == "amplitude":
        reject_amplitude_encoding()
    if args.quantum == "angle":
        vqc_name = f"vqc_angle_{args.quantum_preset}"
        _, _, validation_predictions[vqc_name], test_predictions[vqc_name] = fit_angle_vqc(
            feature_matrix(bc_train),
            (bc_train["sample_label"] == "b").to_numpy(int),
            feature_matrix(bc_val),
            feature_matrix(bc_test),
            maxiter=args.quantum_maxiter,
            seed=args.seed,
            preset=args.quantum_preset,
            shots=args.quantum_shots,
        )
    keyed_output(bc_test, test_predictions).to_parquet(output_dir / "bc_test_scores.parquet", index=False)
    comparison = {
        name: evaluate_bc_model(bc_val, validation_predictions[name], bc_test, test_predictions[name], args.target_b_efficiency, name)
        for name in test_predictions
    }
    write_comparison_figure(output_dir / "hybrid_bc_model_comparison.png", bc_test, test_predictions, comparison)
    report = {
        "mode": args.mode,
        "feature_columns": list(BC_FEATURES),
        "semantic_contract": "Primary RF is b-vs-all. Conditional outputs are P(b | b/c study), not b-vs-all probabilities; g/uds are excluded from reranker outputs.",
        "gate": {"score_source_train": "grouped_out_of_fold_rf", "group_key": "global_event_id", "threshold_source": "validation", "threshold": threshold, "half_width": half_width, "train": contamination_report(gate_train, gate_mask(oof_scores, threshold, half_width)), "validation": contamination_report(splits["val"], gate_mask(val_scores, threshold, half_width)), "test": contamination_report(splits["test"], gate_mask(test_scores, threshold, half_width))},
        "same_keyed_bc_rows": {
            "train": {"count": int(len(bc_train)), "fingerprint": keyed_row_fingerprint(bc_train, "train")},
            "validation": {"count": int(len(bc_val)), "fingerprint": keyed_row_fingerprint(bc_val, "validation")},
            "test": {"count": int(len(bc_test)), "fingerprint": keyed_row_fingerprint(bc_test, "test")},
        },
        "model_comparison": {"row_set": "Each model is evaluated on the same keyed, truth-selected b/c validation and test gate rows.", "models": comparison},
        "classical_controls": list(models),
        "quantum": {
            "requested": args.quantum,
            "encoding": "angle (4 physical features, 4 qubits)" if args.quantum == "angle" else "disabled",
            "preset": {"name": args.quantum_preset, **angle_vqc_preset(args.quantum_preset)} if args.quantum == "angle" else None,
            "optimizer": {"name": "COBYLA", "maxiter": args.quantum_maxiter} if args.quantum == "angle" else None,
            "shots": args.quantum_shots if args.quantum == "angle" else None,
            "seed": args.seed if args.quantum == "angle" else None,
            "library_versions": {"qiskit": version("qiskit"), "qiskit-aer": version("qiskit-aer"), "qiskit-machine-learning": version("qiskit-machine-learning")} if args.quantum == "angle" else None,
            "claim": "Optional local-simulator control only; no quantum-advantage claim.",
        },
    }
    (output_dir / "hybrid_bc_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    print(json.dumps(run(parse_args()), indent=2))


if __name__ == "__main__":
    main()
