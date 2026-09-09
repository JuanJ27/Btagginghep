# Btagginghep

Truth-labelled b-tagging workflow for Delphes reconstructed jets. The primary
tagger is a Random Forest (RF) b-vs-all classifier; the optional hybrid study
is a bounded, conditional b-vs-c reranker and never replaces the RF score.

## Dataset and truth contract

`ml/build_dataset.py` reads per-jet `Jet.Flavor` and maps its PDG-like values
to `b` (`abs == 5`), `c` (`abs == 4`), `uds` (`abs in {1,2,3}`), and `g`
(`== 21`). This truth label is distinct from `source_sample_label`, which
records the input-folder provenance. With `--label-source jet_flavor`, only
mapped truth jets enter the supervised b-vs-all splits.

The builder also writes `train_companion`, `val_companion`, and
`test_companion` files containing every reconstructed jet passing the selected
pT region, keyed by `global_event_id`. Companion truth is nullable:
unmapped `Jet.Flavor` values have `truth_known=False` and `truth_label`/`is_b`
set to `NA`, so they remain scoreable but are excluded from truth metrics.
There is no per-event jet limit by default. Use
`--max-jets-per-event N` only to apply an explicit positive cap.

## Evaluation

`ml/evaluation.py` freezes four working points on validation scores and applies
them unchanged to test scores:

- `purity_80`, `purity_90`
- `efficiency_80`, `efficiency_90`

Reports include per-jet b efficiency, b purity, global non-b rejection, and
c/g/uds mistag and rejection rates. Event reports group candidates by
`global_event_id` and include any-b efficiency/purity, macro per-event jet
metrics, zero-selection rate, pairwise b-over-non-b ranking, and top-1/top-2
b-hit diagnostics. Unknown selected candidates are reported explicitly and
prevent a claimed all-truth event purity.

## Reproducible smoke commands

The checked-in `high20` dataset is an **all-pT smoke test**, not low-pT thesis
evidence: its manifest declares `pt_region: all`, even though its source
profile uses high-threshold inputs.

```bash
# RF all-pT smoke study
python -m ml.train_rf \
  --dataset-dir outputs/high20_jet_flavor \
  --output-dir outputs/rf_high20_jet_flavor

# Conditional b-vs-c hybrid smoke study (local classical path)
python -m ml.train_hybrid_bc \
  --dataset-dir outputs/high20_jet_flavor \
  --output-dir outputs/hybrid_bc_high20_smoke \
  --mode SMOKE_TEST --quantum disabled

# Execute the RF event-smoke presentation notebook
jupyter nbconvert --to notebook --execute --inplace \
  notebooks/rf_high20_event_smoke_presentation.ipynb
```

The historical `high20_jet_flavor` files do not contain companion splits, so
the presentation notebook honestly evaluates its events from supervised test
rows only. Newly built datasets include companions for all-reconstructed-event
evaluation.

## Local-Aer VQC presets

The conditional b/c VQC has three four-feature, four-qubit presets for bounded
ablation studies. They retain the same train-only scaling, keyed b/c rows,
seed, shots, and optimizer, so each run changes one circuit-depth dimension:

- `zz1_real1_linear`: baseline, one ZZ feature-map and one RealAmplitudes layer.
- `zz2_real1_linear`: two ZZ feature-map layers; tests additional data encoding depth.
- `zz1_real2_linear`: two RealAmplitudes layers; tests additional variational capacity.

For example, run a local simulator-only smoke comparison with:

```bash
python -m ml.train_hybrid_bc \
  --dataset-dir outputs/high20_jet_flavor \
  --output-dir outputs/hybrid_bc_vqc_zz1_real2 \
  --mode SMOKE_TEST --quantum angle \
  --quantum-preset zz1_real2_linear --quantum-shots 1024 \
  --quantum-maxiter 50
```

Each report persists the preset, package versions, seed, shots, and a hash of
the b/c rows used in each split. These are local-Aer controls only and do not
establish quantum advantage.

## Authoritative low-pT gate

`THESIS_RUN_AUTHORITATIVE_LOWPT` in `ml/hybrid_bc_contract.py` is
non-negotiable: the manifest must select `label_source: jet_flavor` and
`pt_region: lt20`; every train/validation/test row must have `jet_pt < 20`; and
each split must retain unique jet keys (`global_event_id`, `root_file`,
`event_in_file`, `jet_rank`) plus `jet_flavor`, `sample_label`, `is_b`, and
`jet_pt`. All-pT and high-pT inputs are accepted only in `SMOKE_TEST` mode.

Qiskit is optional and used only for the local `AerSimulator` angle-encoded
VQC hook in `ml/hybrid_bc_models.py`. Install the compatible `qiskit`,
`qiskit-aer`, and `qiskit-machine-learning` dependencies only when using that
path; there is no provider or hardware execution path.

Generated datasets, models, reports, plots, ROOT files, and PDFs are
intentionally ignored by `.gitignore` (including `/outputs/`).
