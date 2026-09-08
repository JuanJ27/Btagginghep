"""Run the bounded, truth-labelled conditional b-versus-c hybrid experiment."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ml.hybrid_bc_contract import KEY_COLUMNS, SMOKE_TEST, THESIS_RUN_AUTHORITATIVE_LOWPT, HybridContractError, load_manifest, validate_mode
from ml.hybrid_bc_features import BC_FEATURES, bc_rows, feature_matrix, validate_feature_frame
from ml.hybrid_bc_gate import contamination_report, fit_final_gate, gate_mask, grouped_oof_rf_scores, select_band_half_width, select_rf_threshold
from ml.hybrid_bc_models import fit_angle_vqc, fit_controls, reject_amplitude_encoding


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
    output = frame.loc[:, [*KEY_COLUMNS, "sample_label", "jet_flavor"]].copy()
    for name, values in scores.items():
        output[f"score_bc_b_{name}"] = values
    output["score_semantics"] = "P(b | truth-selected b/c study); never a b-vs-all score and never a g/uds classification."
    return output


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
    models = fit_controls(feature_matrix(bc_train), (bc_train["sample_label"] == "b").to_numpy(int), args.seed)
    test_predictions = {name: model.predict_proba(feature_matrix(bc_test))[:, 1] for name, model in models.items()}
    if args.quantum == "amplitude":
        reject_amplitude_encoding()
    if args.quantum == "angle":
        _, _, _, test_predictions["vqc_angle"] = fit_angle_vqc(
            feature_matrix(bc_train),
            (bc_train["sample_label"] == "b").to_numpy(int),
            feature_matrix(bc_val),
            feature_matrix(bc_test),
            maxiter=args.quantum_maxiter,
        )
    keyed_output(bc_test, test_predictions).to_parquet(output_dir / "bc_test_scores.parquet", index=False)
    report = {
        "mode": args.mode,
        "feature_columns": list(BC_FEATURES),
        "semantic_contract": "Primary RF is b-vs-all. Conditional outputs are P(b | b/c study), not b-vs-all probabilities; g/uds are excluded from reranker outputs.",
        "gate": {"score_source_train": "grouped_out_of_fold_rf", "group_key": "global_event_id", "threshold_source": "validation", "threshold": threshold, "half_width": half_width, "train": contamination_report(gate_train, gate_mask(oof_scores, threshold, half_width)), "validation": contamination_report(splits["val"], gate_mask(val_scores, threshold, half_width)), "test": contamination_report(splits["test"], gate_mask(test_scores, threshold, half_width))},
        "same_keyed_bc_rows": {"train": int(len(bc_train)), "validation": int(len(bc_val)), "test": int(len(bc_test))},
        "classical_controls": list(models),
        "quantum": {"requested": args.quantum, "encoding": "angle (4 qubits, shallow circuit)" if args.quantum == "angle" else "disabled", "optimizer": {"name": "COBYLA", "maxiter": args.quantum_maxiter} if args.quantum == "angle" else None, "claim": "Optional local-simulator control only; no quantum-advantage claim."},
    }
    (output_dir / "hybrid_bc_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    print(json.dumps(run(parse_args()), indent=2))


if __name__ == "__main__":
    main()
