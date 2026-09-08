import numpy as np
import pandas as pd

from ml.evaluation import bootstrap_event_metrics, event_metrics, evaluate_working_points, object_metrics, select_working_points
from ml.train_rf import score_companion_split, select_score_columns


def scores() -> pd.DataFrame:
    return pd.DataFrame({
        "global_event_id": ["one", "one", "two", "two", "three", "three"],
        "jet_flavor": [5, 4, 5, 1, 5, 21],
        "sample_label": ["b", "c", "b", "uds", "b", "g"],
        "score": [0.9, 0.9, 0.8, 0.3, 0.6, 0.6],
    })


def test_wp_selection_uses_distinct_score_groups_and_validation_only():
    result = select_working_points(scores())
    assert result["efficiency_80"]["threshold"] == 0.6
    assert result["purity_80"]["feasible"] is False
    assert result["purity_90"]["feasible"] is False


def test_object_metrics_use_jet_flavor_not_source_or_wrong_sample_label():
    frame = scores().assign(sample_label="b", source_sample_label="b")
    metrics = object_metrics(frame, 0.8)
    assert (metrics["tp"], metrics["fp"], metrics["fn"], metrics["tn"]) == (2, 1, 1, 2)
    assert metrics["c_mistag"] == 1.0


def test_event_metrics_have_explicit_denominators_and_fractional_ties():
    metrics = event_metrics(scores(), 0.8)
    assert metrics["known_truth_events"] == 3
    assert metrics["any_b_event_efficiency_denominator"] == 3
    assert metrics["top1_b_hit"] == 2 / 3  # events one and three are b/non-b score ties
    assert metrics["top2_b_hit"] == 1.0
    assert metrics["zero_selection_count"] == 1
    assert metrics["pairwise_b_over_non_b_denominator_pairs"] == 3


def test_unknown_truth_is_excluded_not_counted_as_non_b():
    frame = scores().copy()
    frame.loc[1, "jet_flavor"] = 0
    frame["truth_known"] = [True, False, True, True, True, True]
    metrics = object_metrics(frame, 0.8)
    assert metrics["fp"] == 0
    assert metrics["known_truth_jets"] == 5


def test_unknown_high_score_companion_changes_event_candidates_not_truth_denominators():
    frame = pd.DataFrame({
        "global_event_id": ["one", "one", "two", "two"],
        "jet_flavor": [5, 0, 1, 0],
        "truth_known": [True, False, True, False],
        "score": [0.2, 0.99, 0.1, 0.98],
    })
    metrics = event_metrics(frame, 0.8)
    assert metrics["candidate_events"] == 2
    assert metrics["known_truth_events"] == 2
    assert metrics["candidate_jets"] == 4
    assert metrics["unknown_candidate_jets"] == 2
    assert metrics["unknown_selected_candidate_jets"] == 2
    assert metrics["any_b_event_efficiency"] == 1.0
    assert metrics["any_b_event_purity"] is None
    assert metrics["any_b_event_purity_known_truth_only"] is None
    assert metrics["top1_b_hit"] is None
    assert metrics["top1_confirmed_b_hit"] == 0.0
    assert metrics["top1_unknown_candidates_at_or_above_boundary"] == 1
    assert metrics["zero_selection_count"] == 0
    assert object_metrics(frame, 0.8)["fp"] == 0


def test_companion_scoring_keeps_unknown_truth_in_keyed_output():
    class Model:
        def predict_proba(self, matrix):
            return np.array([[0.1, 0.9] for _ in matrix])

    companion = pd.DataFrame({
        "global_event_id": ["one"], "root_file": ["input.root"], "event_in_file": [0], "jet_rank": [1],
        "jet_flavor": [0], "truth_label": [None], "truth_known": [False], "is_b": [None], "feature": [2.0],
    })
    assert score_companion_split(Model(), companion, ["feature"], "test").tolist() == [0.9]
    assert select_score_columns(companion) == [
        "global_event_id", "root_file", "event_in_file", "jet_rank", "jet_flavor", "truth_label", "truth_known", "is_b",
    ]


def test_evaluation_applies_frozen_threshold_without_reselection():
    frozen = {"fixed": {"threshold": 0.8, "feasible": True}}
    report = evaluate_working_points(scores(), frozen)
    assert report["fixed"]["object"]["threshold"] == 0.8


def test_event_bootstrap_is_seeded_and_keeps_the_frozen_threshold():
    first = bootstrap_event_metrics(scores(), 0.8, n_resamples=10, seed=7)
    second = bootstrap_event_metrics(scores(), 0.8, n_resamples=10, seed=7)
    assert first == second
    assert first["threshold"] == 0.8
