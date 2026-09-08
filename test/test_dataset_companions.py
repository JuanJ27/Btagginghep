from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from ml import build_dataset


def test_companion_rows_retain_unknown_truth_without_capping_jets(tmp_path, monkeypatch):
    jet_count = 5
    batch = {
        "Jet.PT": [list(range(10, 10 + jet_count))],
        "Jet.Flavor": [[5.0, 0.0, 4.0, 21.0, 1.0]],
        "Jet.Eta": [[0.0] * jet_count],
        "Jet.Phi": [[0.0] * jet_count],
        "Jet.Mass": [[1.0] * jet_count],
        "Jet.NCharged": [[1.0] * jet_count],
        "Jet.NNeutrals": [[0.0] * jet_count],
        "Track.PT": [[]],
        "Track.Eta": [[]],
        "Track.Phi": [[]],
        "Track.D0": [[]],
        "Track.DZ": [[]],
    }
    monkeypatch.setattr(build_dataset, "ak", SimpleNamespace(to_numpy=lambda values: values))
    monkeypatch.setattr(build_dataset, "uproot", SimpleNamespace(iterate=lambda *args, **kwargs: [batch]))

    supervised, companions, summary = build_dataset.process_one_file(
        Path(tmp_path / "sample.root"), "b", "all", "jet_flavor"
    )

    assert len(companions) == jet_count
    assert len(supervised) == jet_count - 1
    assert companions["jet_rank"].tolist() == list(range(jet_count))
    assert companions["truth_known"].tolist() == [True, False, True, True, True]
    assert pd.isna(companions.loc[1, "truth_label"])
    assert summary["unmatched_jets"] == 1


def test_companion_and_supervised_splits_share_event_assignment():
    companions = pd.DataFrame(
        {
            "global_event_id": [f"event-{event}" for event in range(24)],
            "source_sample_label": ["b" if event % 2 else "uds" for event in range(24)],
        }
    )
    supervised = companions.iloc[::2].assign(sample_label="b", is_b=1).copy()

    split = build_dataset.split_supervised_and_companions(
        supervised, companions, {"train": 0.5, "val": 0.25, "test": 0.25}
    )
    supervised_splits, companion_splits = split[:3], split[3:6]

    for supervised_frame, companion_frame in zip(supervised_splits, companion_splits, strict=True):
        assert set(supervised_frame["global_event_id"]).issubset(set(companion_frame["global_event_id"]))
