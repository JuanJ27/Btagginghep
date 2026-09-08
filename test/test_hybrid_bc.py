import numpy as np
import pandas as pd
import pytest
import sys
from types import ModuleType

from ml.hybrid_bc_contract import HybridContractError, SMOKE_TEST, THESIS_RUN_AUTHORITATIVE_LOWPT, validate_mode
from ml.hybrid_bc_features import BC_FEATURES, bc_rows, feature_matrix
from ml.hybrid_bc_gate import gate_mask, select_band_half_width
from ml.hybrid_bc_models import build_angle_encoded_vqc, fit_angle_vqc, reject_amplitude_encoding, require_qiskit_vqc


def frame(labels=("b", "c", "g", "uds")):
    rows = []
    for index, label in enumerate(labels):
        rows.append({"global_event_id": f"event-{index}", "root_file": "x.root", "event_in_file": index, "jet_rank": 0, "jet_flavor": {"b": 5, "c": 4, "g": 21, "uds": 1}[label], "sample_label": label, "source_sample_label": label, "is_b": int(label == "b"), "jet_pt": 10.0, "n_tracks_cone": 2.0, "mean_abs_d0": .1, "max_abs_d0": .2, "pt_weighted_abs_d0": .15})
    return pd.DataFrame(rows)


def test_bc_features_are_physical_and_bc_rows_do_not_relabel_other_flavours():
    data = frame()
    assert len(BC_FEATURES) == 4
    assert list(bc_rows(data)["sample_label"]) == ["b", "c"]
    assert feature_matrix(data).shape == (4, 4)


def test_authoritative_mode_rejects_all_pt_manifest():
    splits = {name: frame() for name in ("train", "val", "test")}
    manifest = {"pt_region": "all", "label_source": {"selected": "jet_flavor"}}
    with pytest.raises(HybridContractError, match="pt_region='lt20'"):
        validate_mode(manifest, splits, THESIS_RUN_AUTHORITATIVE_LOWPT)
    validate_mode(manifest, splits, SMOKE_TEST)


def test_authoritative_mode_checks_actual_pt_values():
    splits = {name: frame() for name in ("train", "val", "test")}
    splits["test"].loc[0, "jet_pt"] = 20.0
    manifest = {"pt_region": "lt20", "label_source": {"selected": "jet_flavor"}}
    with pytest.raises(HybridContractError, match="jet_pt >= 20"):
        validate_mode(manifest, splits, THESIS_RUN_AUTHORITATIVE_LOWPT)


def test_gate_band_is_score_observable():
    scores = np.array([.1, .49, .51, .9])
    width = select_band_half_width(scores, .5)
    selected = gate_mask(scores, .5, width)
    assert width in {.02, .05, .10, .15}
    assert selected.tolist() == [False, True, True, False]


def test_quantum_paths_are_explicit_when_dependencies_are_absent():
    with pytest.raises(HybridContractError, match="Amplitude encoding"):
        reject_amplitude_encoding()
    with pytest.raises(HybridContractError, match="Quantum backend requested"):
        require_qiskit_vqc()


def test_angle_vqc_uses_qiskit_2_circuit_factories_lazily(monkeypatch):
    calls = {}
    library = ModuleType("qiskit.circuit.library")
    library.zz_feature_map = lambda **kwargs: calls.setdefault("feature_map", kwargs) or "feature-map"
    library.real_amplitudes = lambda **kwargs: calls.setdefault("ansatz", kwargs) or "ansatz"
    aer = ModuleType("qiskit_aer")
    aer.AerSimulator = lambda: calls.setdefault("simulator", "local") or "simulator"
    primitives = ModuleType("qiskit_aer.primitives")

    class SamplerV2:
        @classmethod
        def from_backend(cls, backend):
            calls["sampler_backend"] = backend
            return "local-aer-sampler"

    primitives.SamplerV2 = SamplerV2
    algorithms = ModuleType("qiskit_machine_learning.algorithms")
    algorithms.VQC = lambda **kwargs: calls.setdefault("vqc", kwargs) or "vqc"
    optimizers = ModuleType("qiskit_machine_learning.optimizers")

    def cobyla(**kwargs):
        calls["optimizer"] = kwargs
        return "cobyla"

    optimizers.COBYLA = cobyla
    monkeypatch.setitem(sys.modules, "qiskit", ModuleType("qiskit"))
    monkeypatch.setitem(sys.modules, "qiskit.circuit", ModuleType("qiskit.circuit"))
    monkeypatch.setitem(sys.modules, "qiskit.circuit.library", library)
    monkeypatch.setitem(sys.modules, "qiskit_aer", aer)
    monkeypatch.setitem(sys.modules, "qiskit_aer.primitives", primitives)
    monkeypatch.setitem(sys.modules, "qiskit_machine_learning", ModuleType("qiskit_machine_learning"))
    monkeypatch.setitem(sys.modules, "qiskit_machine_learning.algorithms", algorithms)
    monkeypatch.setitem(sys.modules, "qiskit_machine_learning.optimizers", optimizers)

    build_angle_encoded_vqc()

    assert calls["feature_map"] == {"feature_dimension": 4, "reps": 1}
    assert calls["ansatz"] == {"num_qubits": 4, "reps": 1, "entanglement": "linear"}
    assert calls["simulator"] == "local"
    assert calls["sampler_backend"] == "local"
    assert calls["vqc"]["sampler"] == "local-aer-sampler"
    assert calls["optimizer"] == {"maxiter": 50}
    assert calls["vqc"]["optimizer"] == "cobyla"


def test_angle_vqc_fits_train_only_angular_preprocessing_and_scores_test(monkeypatch):
    from ml import hybrid_bc_models
    calls = {}

    class FakeVQC:
        def fit(self, x, y):
            calls["fit"] = (x.copy(), y.copy())
            return self

        def predict_proba(self, x):
            calls["predict"] = x.copy()
            return np.tile([0.25, 0.75], (len(x), 1))

    def build_fake_vqc(maxiter):
        calls["maxiter"] = maxiter
        return FakeVQC()

    monkeypatch.setattr(hybrid_bc_models, "build_angle_encoded_vqc", build_fake_vqc)
    x_train = np.array([[10., 10., 10., 10.], [20., 20., 20., 20.]])
    x_val = np.array([[0., 0., 0., 0.]])
    x_test = np.array([[30., 30., 30., 30.]])

    _, preprocessor, angular_val, test_b_probability = fit_angle_vqc(x_train, np.array([0, 1]), x_val, x_test, maxiter=7)

    assert preprocessor.data_min_.tolist() == [10.] * 4
    assert preprocessor.data_max_.tolist() == [20.] * 4
    assert np.all(calls["fit"][0] >= -np.pi) and np.all(calls["fit"][0] <= np.pi)
    assert angular_val.tolist() == [[-np.pi] * 4]
    assert calls["predict"].tolist() == [[np.pi] * 4]
    assert calls["maxiter"] == 7
    assert test_b_probability.tolist() == [0.75]
