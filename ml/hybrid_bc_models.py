"""Train-only-preprocessed classical b/c controls and an optional simulator-only VQC hook."""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.svm import SVC

from ml.hybrid_bc_contract import HybridContractError


def classical_controls(seed: int) -> dict:
    return {
        "logistic_regression": make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, class_weight="balanced", random_state=seed)),
        "rbf_svm": make_pipeline(StandardScaler(), CalibratedClassifierCV(SVC(kernel="rbf", class_weight="balanced", random_state=seed), ensemble=False)),
        "restricted_rf": RandomForestClassifier(n_estimators=150, max_depth=6, min_samples_leaf=10, class_weight="balanced", random_state=seed, n_jobs=1),
        "compact_mlp": make_pipeline(StandardScaler(), MLPClassifier(hidden_layer_sizes=(12,), max_iter=300, early_stopping=True, random_state=seed)),
    }


def fit_controls(x_train: np.ndarray, y_train: np.ndarray, seed: int) -> dict:
    if set(y_train) != {0, 1}:
        raise HybridContractError("b-vs-c controls require both b and c training rows.")
    models = classical_controls(seed)
    for model in models.values():
        model.fit(x_train, y_train)
    return models


def require_qiskit_vqc():
    """Lazy simulator-only import; this module contains no provider or hardware path."""
    try:
        from qiskit.circuit.library import real_amplitudes, zz_feature_map
        from qiskit_aer import AerSimulator
        from qiskit_aer.primitives import SamplerV2
        from qiskit_machine_learning.algorithms import VQC
        from qiskit_machine_learning.optimizers import COBYLA
    except ModuleNotFoundError as error:
        raise HybridContractError(
            "Quantum backend requested but unavailable. Install compatible qiskit, qiskit-machine-learning, and qiskit-aer packages; no hardware provider is used."
        ) from error
    return zz_feature_map, real_amplitudes, AerSimulator, SamplerV2, COBYLA, VQC


def build_angle_encoded_vqc(maxiter: int = 50, seed: int = 42):
    """Build the planned four-qubit, shallow angle-encoding VQC configuration."""
    if maxiter < 1:
        raise HybridContractError("VQC COBYLA maxiter must be positive.")
    zz_feature_map, real_amplitudes, AerSimulator, SamplerV2, COBYLA, VQC = require_qiskit_vqc()
    simulator = AerSimulator()  # Explicitly local simulator; never a provider/backend service.
    sampler = SamplerV2.from_backend(simulator, default_shots=1024, seed=seed)
    feature_map = zz_feature_map(feature_dimension=4, reps=1)
    ansatz = real_amplitudes(num_qubits=4, reps=1, entanglement="linear")
    initial_point = np.random.default_rng(seed).uniform(-0.1, 0.1, ansatz.num_parameters)
    return VQC(feature_map=feature_map, ansatz=ansatz, sampler=sampler, optimizer=COBYLA(maxiter=maxiter), initial_point=initial_point)


def fit_angle_vqc(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_validation: np.ndarray,
    x_test: np.ndarray,
    maxiter: int = 50,
    seed: int = 42,
) -> tuple[object, MinMaxScaler, np.ndarray, np.ndarray]:
    """Fit angular scaling on b/c training rows only, then score validation and test rows."""
    if x_train.shape[1] != 4:
        raise HybridContractError("The four-qubit angle VQC requires exactly four physical features.")
    if set(y_train) != {0, 1}:
        raise HybridContractError("Angle VQC requires both b and c training rows.")
    preprocessor = MinMaxScaler(feature_range=(-np.pi, np.pi), clip=True).fit(x_train)
    angular_train = preprocessor.transform(x_train)
    angular_validation = preprocessor.transform(x_validation)
    angular_test = preprocessor.transform(x_test)
    model = build_angle_encoded_vqc(maxiter=maxiter, seed=seed)
    model.fit(angular_train, y_train)
    validation_probabilities = np.asarray(model.predict_proba(angular_validation), dtype=float)
    test_probabilities = np.asarray(model.predict_proba(angular_test), dtype=float)
    for split_name, probabilities in (("validation", validation_probabilities), ("test", test_probabilities)):
        if probabilities.ndim != 2 or probabilities.shape != (len(x_validation) if split_name == "validation" else len(x_test), 2):
            raise HybridContractError(f"VQC {split_name} predict_proba must return one b/c probability pair per row.")
        if not np.isfinite(probabilities).all() or (probabilities < 0).any() or (probabilities > 1).any():
            raise HybridContractError(f"VQC {split_name} probabilities must be finite values in [0, 1].")
    return model, preprocessor, validation_probabilities[:, 1], test_probabilities[:, 1]


def reject_amplitude_encoding() -> None:
    raise HybridContractError(
        "Amplitude encoding is intentionally unsupported/config-disabled in this scaffold. "
        "Roadmap: add a normalized state-preparation ablation with explicit feature dimension and noise controls."
    )
