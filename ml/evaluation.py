"""Frozen-threshold object and event evaluation for b-tagging scores."""
from __future__ import annotations

from math import comb
from typing import Any

import numpy as np
import pandas as pd


TRUTH_LABELS = ("b", "c", "g", "uds")
WORKING_POINTS = {
    "purity_80": ("purity", 0.80),
    "purity_90": ("purity", 0.90),
    "efficiency_80": ("efficiency", 0.80),
    "efficiency_90": ("efficiency", 0.90),
}


def _truth_labels(scores: pd.DataFrame) -> pd.Series:
    """Return truth labels without ever consulting source provenance."""
    if "jet_flavor" in scores:
        flavor = pd.to_numeric(scores["jet_flavor"], errors="coerce")
        labels = pd.Series(index=scores.index, dtype="object")
        labels.loc[flavor.abs() == 5] = "b"
        labels.loc[flavor.abs() == 4] = "c"
        labels.loc[flavor.abs().isin([1, 2, 3])] = "uds"
        labels.loc[flavor == 21] = "g"
        return labels
    if "sample_label" not in scores:
        raise ValueError("Score table needs jet_flavor or truth-derived sample_label.")
    return scores["sample_label"].where(scores["sample_label"].isin(TRUTH_LABELS))


def normalized_scores(scores: pd.DataFrame, score_column: str = "score") -> pd.DataFrame:
    """Normalize every scored candidate and explicitly mark recognized truth."""
    required = {"global_event_id", score_column}
    missing = sorted(required - set(scores.columns))
    if missing:
        raise ValueError(f"Score table is missing required columns: {missing}")
    result = scores.copy()
    result["_truth_label"] = _truth_labels(result)
    explicit_known = result["truth_known"].eq(True).fillna(False) if "truth_known" in result else True
    result["_truth_known"] = result["_truth_label"].isin(TRUTH_LABELS) & explicit_known
    result.loc[~result["_truth_known"], "_truth_label"] = None
    result["_score"] = pd.to_numeric(result[score_column], errors="coerce")
    if result["_score"].isna().any() or not np.isfinite(result["_score"]).all():
        raise ValueError("Scores must be finite for every candidate jet.")
    return result


def known_truth_scores(scores: pd.DataFrame, score_column: str = "score") -> pd.DataFrame:
    """Normalize a score table and retain only jets with recognized b/c/g/uds truth."""
    result = normalized_scores(scores, score_column)
    result = result.loc[result["_truth_known"]].copy()
    return result


def select_working_points(validation_scores: pd.DataFrame, score_column: str = "score") -> dict[str, dict[str, Any]]:
    """Select four WPs from distinct validation score groups, never splitting ties."""
    frame = known_truth_scores(validation_scores, score_column)
    if not (frame["_truth_label"] == "b").any():
        raise ValueError("Working-point selection requires at least one known-truth b jet.")
    candidates = np.sort(frame["_score"].unique())
    rows = []
    for threshold in candidates:
        selected = frame["_score"] >= threshold
        tp = int((selected & (frame["_truth_label"] == "b")).sum())
        n_selected = int(selected.sum())
        rows.append((float(threshold), tp / int((frame["_truth_label"] == "b").sum()), tp / n_selected))

    report: dict[str, dict[str, Any]] = {}
    for name, (kind, target) in WORKING_POINTS.items():
        eligible = [row for row in rows if row[1 if kind == "efficiency" else 2] >= target]
        if not eligible:
            report[name] = {"criterion": kind, "target": target, "feasible": False, "threshold": None}
            continue
        # Strictest means highest threshold; loosest means lowest threshold.
        chosen = max(eligible, key=lambda row: row[0]) if kind == "efficiency" else min(eligible, key=lambda row: row[0])
        report[name] = {
            "criterion": kind, "target": target, "feasible": True, "threshold": chosen[0],
            "validation_b_efficiency": chosen[1], "validation_b_purity": chosen[2],
        }
    return report


def object_metrics(scores: pd.DataFrame, threshold: float, score_column: str = "score") -> dict[str, Any]:
    frame = known_truth_scores(scores, score_column)
    selected = frame["_score"] >= threshold
    b = frame["_truth_label"] == "b"
    non_b = ~b
    tp, fp = int((selected & b).sum()), int((selected & non_b).sum())
    fn, tn = int((~selected & b).sum()), int((~selected & non_b).sum())
    metrics: dict[str, Any] = {
        "threshold": float(threshold), "known_truth_jets": int(len(frame)), "selected_yield": int(selected.sum()),
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "b_efficiency": None if tp + fn == 0 else tp / (tp + fn),
        "b_purity": None if tp + fp == 0 else tp / (tp + fp),
        "global_non_b_fpr": None if fp + tn == 0 else fp / (fp + tn),
    }
    metrics["global_non_b_rejection"] = None if not metrics["global_non_b_fpr"] else 1 / metrics["global_non_b_fpr"]
    for label in ("c", "g", "uds"):
        members = frame["_truth_label"] == label
        rate = None if not members.any() else float(selected[members].mean())
        metrics[f"{label}_mistag"] = rate
        metrics[f"{label}_rejection"] = None if not rate else 1 / rate
    return metrics


def _top_k_hit(group: pd.DataFrame, k: int) -> float:
    """Expected b hit under uniform random resolution of a boundary score tie."""
    ordered = group.sort_values("_score", ascending=False)
    before = ordered.iloc[:k]
    if len(ordered) <= k or before.iloc[-1]["_score"] > ordered.iloc[k]["_score"]:
        return float((before["_truth_label"] == "b").any())
    boundary = before.iloc[-1]["_score"]
    higher = ordered[ordered["_score"] > boundary]
    tied = ordered[ordered["_score"] == boundary]
    if (higher["_truth_label"] == "b").any():
        return 1.0
    slots, total, b_count = k - len(higher), len(tied), int((tied["_truth_label"] == "b").sum())
    return 1.0 - comb(total - b_count, slots) / comb(total, slots)


def _top_k_boundary_unknown_count(group: pd.DataFrame, k: int) -> int:
    """Count unknown candidates at or above the score boundary for top-k."""
    if len(group) < k:
        return 0
    boundary = group.nlargest(k, "_score")["_score"].iloc[-1]
    return int(((group["_score"] >= boundary) & ~group["_truth_known"]).sum())


def event_metrics(scores: pd.DataFrame, threshold: float, score_column: str = "score") -> dict[str, Any]:
    frame = normalized_scores(scores, score_column)
    events = list(frame.groupby("global_event_id", sort=False))
    any_b_eff, macro_eff, macro_purity, pairwise = [], [], [], []
    top1, top2 = [], []
    top1_unknown_boundary, top2_unknown_boundary = 0, 0
    zero_selected = 0
    pair_count = 0
    unknown_selected_events = 0
    unknown_selected_jets = 0
    any_b_purity_known = []
    for _, event in events:
        known = event.loc[event["_truth_known"]]
        selected_candidates = event["_score"] >= threshold
        selected_known = known["_score"] >= threshold
        selected_unknown = selected_candidates & ~event["_truth_known"]
        b = known["_truth_label"] == "b"
        unknown_selected_jets += int(selected_unknown.sum())
        unknown_selected_events += int(selected_unknown.any())
        if not selected_candidates.any():
            zero_selected += 1
        if b.any():
            any_b_eff.append(float(selected_candidates.any()))
            macro_eff.append(float((selected_known & b).sum() / b.sum()))
            top1.append(_top_k_hit(event, 1))
            top1_unknown_boundary += _top_k_boundary_unknown_count(event, 1)
            if len(event) >= 2:
                top2.append(_top_k_hit(event, 2))
                top2_unknown_boundary += _top_k_boundary_unknown_count(event, 2)
        if selected_known.any():
            any_b_purity_known.append(float(b.any()))
        if selected_known.any():
            macro_purity.append(float((selected_known & b).sum() / selected_known.sum()))
        positives, negatives = known.loc[b, "_score"].to_numpy(), known.loc[~b, "_score"].to_numpy()
        if len(positives) and len(negatives):
            comparisons = (positives[:, None] > negatives).mean() + 0.5 * (positives[:, None] == negatives).mean()
            pairwise.append(float(comparisons))
            pair_count += len(positives) * len(negatives)
    has_unknown_selected = unknown_selected_jets > 0
    return {
        "threshold": float(threshold), "candidate_events": len(events),
        "known_truth_events": int(sum(event["_truth_known"].any() for _, event in events)),
        "candidate_jets": int(len(frame)), "known_truth_jets": int(frame["_truth_known"].sum()),
        "unknown_candidate_jets": int((~frame["_truth_known"]).sum()),
        "selected_candidate_jets": int((frame["_score"] >= threshold).sum()),
        "selected_known_truth_jets": int(((frame["_score"] >= threshold) & frame["_truth_known"]).sum()),
        "unknown_selected_candidate_jets": unknown_selected_jets, "unknown_selected_candidate_events": unknown_selected_events,
        "any_b_event_efficiency": _mean_or_none(any_b_eff), "any_b_event_efficiency_denominator": len(any_b_eff),
        "any_b_event_purity": None if has_unknown_selected else _mean_or_none(any_b_purity_known),
        "any_b_event_purity_denominator": None if has_unknown_selected else len(any_b_purity_known),
        "any_b_event_purity_known_truth_only": _mean_or_none(any_b_purity_known),
        "any_b_event_purity_known_truth_only_denominator": len(any_b_purity_known),
        "macro_per_event_jet_efficiency": _mean_or_none(macro_eff), "macro_per_event_jet_efficiency_denominator": len(macro_eff),
        "macro_per_event_jet_purity": _mean_or_none(macro_purity), "macro_per_event_jet_purity_denominator": len(macro_purity),
        "macro_per_event_jet_purity_known_truth_only": _mean_or_none(macro_purity),
        "zero_selection_count": zero_selected, "zero_selection_rate": None if not events else zero_selected / len(events),
        "pairwise_b_over_non_b_macro": _mean_or_none(pairwise), "pairwise_b_over_non_b_denominator_events": len(pairwise), "pairwise_b_over_non_b_denominator_pairs": pair_count,
        "top1_b_hit": None if top1_unknown_boundary else _mean_or_none(top1),
        "top1_confirmed_b_hit": _mean_or_none(top1), "top1_b_hit_denominator_events": len(top1),
        "top1_unknown_candidates_at_or_above_boundary": top1_unknown_boundary,
        "top2_b_hit": None if top2_unknown_boundary else _mean_or_none(top2),
        "top2_confirmed_b_hit": _mean_or_none(top2), "top2_b_hit_denominator_events": len(top2),
        "top2_unknown_candidates_at_or_above_boundary": top2_unknown_boundary,
    }


def _mean_or_none(values: list[float]) -> float | None:
    return None if not values else float(np.mean(values))


def evaluate_working_points(
    scores: pd.DataFrame, working_points: dict[str, dict[str, Any]], score_column: str = "score", event_scores: pd.DataFrame | None = None,
) -> dict[str, dict[str, Any]]:
    """Apply already-selected thresholds; this function never selects a threshold."""
    result = {}
    for name, definition in working_points.items():
        threshold = definition.get("threshold")
        result[name] = {"definition": definition, "object": None, "event": None} if threshold is None else {
            "definition": definition, "object": object_metrics(scores, threshold, score_column),
            "event": event_metrics(event_scores if event_scores is not None else scores, threshold, score_column),
        }
    return result


def bootstrap_event_metrics(
    scores: pd.DataFrame, threshold: float, score_column: str = "score", n_resamples: int = 200, seed: int = 42,
) -> dict[str, Any]:
    """Event-bootstrap fixed-threshold metrics; thresholds are never reselected."""
    if n_resamples <= 0:
        raise ValueError("n_resamples must be positive.")
    frame = normalized_scores(scores, score_column)
    grouped = [group for _, group in frame.groupby("global_event_id", sort=False)]
    if not grouped:
        raise ValueError("Event bootstrap requires at least one candidate event.")
    rng = np.random.default_rng(seed)
    metric_names = ("any_b_event_efficiency", "any_b_event_purity", "macro_per_event_jet_efficiency", "macro_per_event_jet_purity", "zero_selection_rate", "pairwise_b_over_non_b_macro", "top1_b_hit", "top2_b_hit")
    samples: dict[str, list[float]] = {name: [] for name in metric_names}
    for _ in range(n_resamples):
        picks = rng.integers(0, len(grouped), len(grouped))
        resampled = []
        for copy_id, picked in enumerate(picks):
            event = grouped[picked].copy()
            event["global_event_id"] = f"bootstrap-{copy_id}"
            resampled.append(event)
        values = event_metrics(pd.concat(resampled, ignore_index=True), threshold, "_score")
        for name in metric_names:
            if values[name] is not None:
                samples[name].append(float(values[name]))
    return {
        "threshold": float(threshold), "n_resamples": n_resamples, "seed": seed,
        "metrics": {
            name: None if not values else {"mean": float(np.mean(values)), "lower_95": float(np.quantile(values, 0.025)), "upper_95": float(np.quantile(values, 0.975))}
            for name, values in samples.items()
        },
    }


def plot_working_point_report(report: dict[str, dict[str, Any]], output_path: str) -> None:
    """Plot object b efficiency/purity and event any-b efficiency for feasible WPs."""
    import matplotlib.pyplot as plt
    names = [name for name, values in report.items() if values["object"] is not None]
    if not names:
        return
    x = np.arange(len(names))
    plt.figure(figsize=(8, 4.5))
    plt.bar(x - 0.25, [report[n]["object"]["b_efficiency"] for n in names], 0.25, label="Object b efficiency")
    plt.bar(x, [report[n]["object"]["b_purity"] for n in names], 0.25, label="Object b purity")
    plt.bar(x + 0.25, [report[n]["event"]["any_b_event_efficiency"] for n in names], 0.25, label="Any-b event efficiency")
    plt.xticks(x, names, rotation=15)
    plt.ylim(0, 1.05)
    plt.ylabel("Metric")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
