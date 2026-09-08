from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve

from ml.evaluation import evaluate_working_points, plot_working_point_report, select_working_points


DEFAULT_DATASET_DIR = Path(__file__).resolve().parents[1] / "outputs"
METADATA_COLUMNS = {"global_event_id", "source_sample_label", "sample_label", "is_b", "jet_flavor", "root_file", "event_in_file", "jet_rank"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train a reproducible Random Forest b-tagging baseline.")
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR, help="Directory containing train, val, and test datasets.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for model artifacts and plots.")
    parser.add_argument("--seed", type=int, default=42, help="Random Forest random seed.")
    parser.add_argument("--n-estimators", type=int, default=500, help="Number of trees.")
    parser.add_argument("--max-depth", type=int, default=None, help="Maximum tree depth; omit for unrestricted depth.")
    parser.add_argument("--min-samples-leaf", type=int, default=1, help="Minimum samples in each leaf.")
    parser.add_argument("--target-efficiency", type=float, default=0.70, help="b-efficiency operating point for the report.")
    return parser.parse_args()


def load_split(dataset_dir: Path, name: str) -> pd.DataFrame:
    parquet_path = dataset_dir / f"{name}.parquet"
    csv_path = dataset_dir / f"{name}.csv"
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    if csv_path.exists():
        return pd.read_csv(csv_path)
    raise FileNotFoundError(f"Missing {name} dataset: expected {parquet_path} or {csv_path}")


def load_optional_split(dataset_dir: Path, name: str) -> pd.DataFrame | None:
    parquet_path = dataset_dir / f"{name}.parquet"
    csv_path = dataset_dir / f"{name}.csv"
    if parquet_path.exists() or csv_path.exists():
        return load_split(dataset_dir, name)
    return None


def select_features(df: pd.DataFrame) -> list[str]:
    if "is_b" not in df:
        raise ValueError("Dataset is missing binary target column 'is_b'.")
    features = [column for column in df.select_dtypes(include=[np.number]).columns if column not in METADATA_COLUMNS]
    if not features:
        raise ValueError("No numeric model features remain after excluding metadata columns.")
    return features


def select_score_columns(df: pd.DataFrame) -> list[str]:
    columns = [column for column in ("global_event_id", "root_file", "event_in_file", "jet_rank") if column in df]
    if "source_sample_label" in df:
        columns.append("source_sample_label")
    columns.extend(column for column in ("jet_flavor", "truth_label", "sample_label", "truth_known", "is_b") if column in df)
    return columns


def score_companion_split(model: RandomForestClassifier, companion: pd.DataFrame, features: list[str], split_name: str) -> np.ndarray:
    missing = sorted(set(features) - set(companion.columns))
    if missing:
        raise ValueError(f"{split_name} companion dataset is missing feature columns: {missing}")
    matrix = companion[features].to_numpy(dtype=np.float64)
    if not np.isfinite(matrix).all():
        bad_columns = [feature for feature in features if not np.isfinite(companion[feature].to_numpy(dtype=np.float64)).all()]
        raise ValueError(f"{split_name} companion dataset has non-finite values in: {bad_columns}")
    return model.predict_proba(matrix)[:, 1]


def validate_inputs(df: pd.DataFrame, features: list[str], split_name: str) -> tuple[np.ndarray, np.ndarray]:
    missing = sorted(set(features) - set(df.columns))
    if missing:
        raise ValueError(f"{split_name} dataset is missing feature columns: {missing}")
    if not {"is_b", "sample_label", "global_event_id"}.issubset(df.columns):
        raise ValueError(f"{split_name} dataset must include is_b, sample_label, and global_event_id.")
    matrix = df[features].to_numpy(dtype=np.float64)
    if not np.isfinite(matrix).all():
        bad_columns = [feature for feature in features if not np.isfinite(df[feature].to_numpy(dtype=np.float64)).all()]
        raise ValueError(f"{split_name} dataset has non-finite values in: {bad_columns}")
    target = df["is_b"].to_numpy(dtype=np.int8)
    if not np.isin(target, [0, 1]).all() or len(np.unique(target)) != 2:
        raise ValueError(f"{split_name} target 'is_b' must contain both binary classes 0 and 1.")
    return matrix, target


def manifest_hash(dataset_dir: Path) -> str | None:
    manifest = dataset_dir / "dataset_manifest.json"
    if not manifest.exists():
        return None
    return hashlib.sha256(manifest.read_bytes()).hexdigest()


def select_operating_threshold(scores: np.ndarray, is_b: np.ndarray, target_efficiency: float) -> float:
    if not 0 < target_efficiency < 1:
        raise ValueError("target_efficiency must be between 0 and 1.")
    if not np.any(is_b == 1):
        raise ValueError("Threshold selection requires at least one b jet.")
    candidates = np.unique(scores)
    eligible = candidates[(scores[:, None] >= candidates)[is_b == 1].mean(axis=0) >= target_efficiency]
    if not len(eligible):
        raise ValueError("No observed score threshold meets the requested b efficiency.")
    return float(eligible.max())


def operating_point(test: pd.DataFrame, scores: np.ndarray, threshold: float) -> dict[str, float | None]:
    is_b = test["is_b"].to_numpy()
    selected = scores >= threshold
    result: dict[str, float | None] = {"threshold": threshold, "b_efficiency": float(selected[is_b == 1].mean())}
    for label in sorted(test["sample_label"].unique()):
        if label == "b":
            continue
        rate = float(selected[test["sample_label"].to_numpy() == label].mean())
        result[f"{label}_mistag"] = rate
        result[f"{label}_rejection"] = None if rate == 0 else 1 / rate
    return result


def plot_score_distribution(test: pd.DataFrame, scores: np.ndarray, output_dir: Path) -> None:
    plt.figure(figsize=(8, 5))
    for label in sorted(test["sample_label"].unique()):
        plt.hist(scores[test["sample_label"].to_numpy() == label], bins=40, density=True, histtype="step", linewidth=1.8, label=label)
    plt.xlabel("Random Forest b score")
    plt.ylabel("Density")
    plt.legend(title="Sample")
    plt.tight_layout()
    plt.savefig(output_dir / "rf_score_by_sample.png", dpi=150)
    plt.close()


def plot_feature_importance(importances: pd.DataFrame, output_dir: Path) -> None:
    top = importances.head(20).sort_values("importance")
    plt.figure(figsize=(8, 6))
    plt.barh(top["feature"], top["importance"])
    plt.xlabel("Random Forest feature importance")
    plt.tight_layout()
    plt.savefig(output_dir / "rf_feature_importance.png", dpi=150)
    plt.close()


def plot_performance(y_test: np.ndarray, scores: np.ndarray, output_dir: Path) -> None:
    false_positive_rate, true_positive_rate, _ = roc_curve(y_test, scores)
    background_rejection = np.divide(1, false_positive_rate, out=np.full_like(false_positive_rate, np.nan), where=false_positive_rate > 0)
    plt.figure(figsize=(6, 5))
    plt.semilogy(true_positive_rate, background_rejection)
    plt.xlabel("b efficiency")
    plt.ylabel("Background rejection")
    plt.ylim(bottom=1)
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "rf_b_efficiency_background_rejection.png", dpi=150)
    plt.close()


def main() -> None:
    args = parse_args()
    dataset_dir = args.dataset_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    train, val, test = (load_split(dataset_dir, name) for name in ("train", "val", "test"))
    val_companion = load_optional_split(dataset_dir, "val_companion")
    test_companion = load_optional_split(dataset_dir, "test_companion")
    features = select_features(train)
    x_train, y_train = validate_inputs(train, features, "train")
    x_val, y_val = validate_inputs(val, features, "validation")
    x_test, y_test = validate_inputs(test, features, "test")
    model = RandomForestClassifier(n_estimators=args.n_estimators, max_depth=args.max_depth, min_samples_leaf=args.min_samples_leaf, class_weight="balanced", random_state=args.seed, n_jobs=-1)
    model.fit(x_train, y_train)
    val_scores = model.predict_proba(x_val)[:, 1]
    test_scores = model.predict_proba(x_test)[:, 1]
    threshold = select_operating_threshold(val_scores, y_val, args.target_efficiency)
    importances = pd.DataFrame({"feature": features, "importance": model.feature_importances_}).sort_values("importance", ascending=False)
    joblib.dump(model, output_dir / "rf_model.joblib")
    joblib.dump(features, output_dir / "rf_features.joblib")
    importances.to_csv(output_dir / "rf_feature_importance.csv", index=False)
    keyed_scores = test[select_score_columns(test)].copy()
    keyed_scores["score_rf"] = test_scores
    keyed_scores.to_csv(output_dir / "rf_test_scores.csv", index=False)
    validation_keyed_scores = val[select_score_columns(val)].copy()
    validation_keyed_scores["score_rf"] = val_scores
    validation_keyed_scores.to_csv(output_dir / "rf_validation_scores.csv", index=False)
    test_companion_scores = None
    if val_companion is not None:
        validation_companion_scores = val_companion[select_score_columns(val_companion)].copy()
        validation_companion_scores["score_rf"] = score_companion_split(model, val_companion, features, "validation")
        validation_companion_scores.to_csv(output_dir / "rf_validation_companion_scores.csv", index=False)
    if test_companion is not None:
        test_companion_scores = test_companion[select_score_columns(test_companion)].copy()
        test_companion_scores["score_rf"] = score_companion_split(model, test_companion, features, "test")
        test_companion_scores.to_csv(output_dir / "rf_test_companion_scores.csv", index=False)
    working_points = select_working_points(validation_keyed_scores, "score_rf")
    working_point_report = evaluate_working_points(keyed_scores, working_points, "score_rf", event_scores=test_companion_scores)
    event_score_source = "test_companion" if test_companion_scores is not None else "test"
    (output_dir / "rf_working_points.json").write_text(json.dumps({"threshold_source": "validation", "event_score_source": event_score_source, "working_points": working_point_report}, indent=2) + "\n")
    plot_working_point_report(working_point_report, str(output_dir / "rf_working_points.png"))
    plot_score_distribution(test, test_scores, output_dir)
    plot_feature_importance(importances, output_dir)
    plot_performance(y_test, test_scores, output_dir)
    metrics = {
        "dataset_dir": str(dataset_dir), "dataset_manifest_sha256": manifest_hash(dataset_dir),
        "dataset_manifest": str(dataset_dir / "dataset_manifest.json"),
        "feature_count": len(features), "validation_auc": float(roc_auc_score(y_val, val_scores)), "test_auc": float(roc_auc_score(y_test, test_scores)),
        "operating_point": {
            "threshold_source": "validation",
            "target_b_efficiency": args.target_efficiency,
            "validation_b_efficiency": float((val_scores[y_val == 1] >= threshold).mean()),
            "test_metrics": operating_point(test, test_scores, threshold),
        },
        "working_points": {"threshold_source": "validation", "event_score_source": event_score_source, "report": "rf_working_points.json"},
        "rf_configuration": {"seed": args.seed, "n_estimators": args.n_estimators, "max_depth": args.max_depth, "min_samples_leaf": args.min_samples_leaf, "class_weight": "balanced", "n_jobs": -1},
    }
    (output_dir / "rf_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
