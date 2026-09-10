from __future__ import annotations

from pathlib import Path

import awkward as ak
import pytest
import uproot

from production.madgraph.run_campaign import CampaignError, SAMPLES, locate_run_artifacts, load_config, patch_delphes_card, prepare_cards, render_build_commands, render_launch_commands, validate_root


def write_config(path: Path, *, include_eta: bool = False) -> None:
    eta = "parton_eta_max = 5.0\n" if include_eta else ""
    path.write_text(
        """[campaign]
name = "test"
mg5_root = "/tmp/mg5"
output_root = "/tmp/campaign"
nevents = 50
runs_per_sample = 1
cores = 2
ebeam1 = 6800.0
ebeam2 = 6800.0
parton_pt_min = 5.0
light_parton_dr_min = 0.4
bb_dr_min = 0.0
jet_radius = 0.4
reco_jet_pt_min = 5.0
keep_hepmc = false
"""
        + eta
        + """[samples]
enabled = ["bbbar", "ccbar"]
[seeds]
bbbar = 1001
ccbar = 2001
"""
    )


def test_config_preserves_eta_as_an_uneditable_madgraph_default(tmp_path):
    config_path = tmp_path / "campaign.toml"
    write_config(config_path)
    config = load_config(config_path)
    assert config.cores == 2
    assert config.enabled_samples == ("bbbar", "ccbar")
    write_config(config_path, include_eta=True)
    with pytest.raises(CampaignError, match="Eta is intentionally not configurable"):
        load_config(config_path)


def test_build_commands_define_a_clean_initial_proton_and_only_inclusive_defines_jets(tmp_path):
    bbbar = render_build_commands(SAMPLES["bbbar"], tmp_path / "bbbar")
    inclusive = render_build_commands(SAMPLES["inclusive"], tmp_path / "inclusive")
    assert "define p = g u d s u~ d~ s~" in bbbar
    assert "define p = p b b~" not in bbbar
    assert "define j =" not in bbbar
    assert "define j = g u d s c b u~ d~ s~ c~ b~" in inclusive


def test_delphes_patch_changes_only_analysis_jet_collections():
    card = """module FastJetFinder GenJetFinder {
  set ParameterR 0.5
  set JetPTMin 20.0
}
module FastJetFinder FastJetFinder {
  set ParameterR 0.5
  set JetPTMin 20.0
}
module FastJetFinder FatJetFinder {
  set ParameterR 0.8
  set JetPTMin 200.0
}
"""
    patched = patch_delphes_card(card, 0.4, 5.0)
    assert patched.count("set ParameterR 0.4") == 2
    assert patched.count("set JetPTMin 5.0") == 2
    assert "set ParameterR 0.8" in patched
    assert "set JetPTMin 200.0" in patched


def test_prepare_cards_uses_b_specific_cuts_without_touching_eta(tmp_path):
    config_path = tmp_path / "campaign.toml"
    write_config(config_path)
    config = load_config(config_path)
    cards = tmp_path / "process" / "Cards"
    cards.mkdir(parents=True)
    (cards / "run_card_default.dat").write_text(
        "50 = nevents ! events\n0 = iseed ! seed\n6500 = ebeam1 ! beam\n6500 = ebeam2 ! beam\n0 = ptb ! b cut\n-1 = etab ! untouched eta\n0 = drbb ! bb cut\n4 = maxjetflavor ! flavor\n"
    )
    (cards / "delphes_card_default.dat").write_text(
        "module FastJetFinder GenJetFinder {\n set ParameterR 0.5\n set JetPTMin 20.0\n}\nmodule FastJetFinder FastJetFinder {\n set ParameterR 0.5\n set JetPTMin 20.0\n}\n"
    )
    (cards / "pythia8_card_default.dat").write_text("HEPMCoutput:file = hepmc.gz\n")
    prepare_cards(config, SAMPLES["bbbar"], tmp_path / "process", 1001)
    run_card = (cards / "run_card.dat").read_text()
    assert "5.0 = ptb" in run_card
    assert "-1 = etab" in run_card
    assert "ptj" not in run_card
    assert (cards / "pythia8_card.dat").read_text() == "HEPMCoutput:file = hepmc.gz\n"


def test_launch_commands_use_one_multicore_madgraph_instance(tmp_path):
    config_path = tmp_path / "campaign.toml"
    write_config(config_path)
    config = load_config(config_path)
    command = render_launch_commands(config, tmp_path / "process", "run_01")
    assert "set run_mode 2 --no_save" in command
    assert "set nb_core 2 --no_save" in command
    assert "launch " in command
    assert " -m" not in command
    assert "run_01" in command


def test_root_validation_accepts_the_delphes_branch_contract(tmp_path):
    root_path = tmp_path / "delphes.root"
    jet_values = ak.Array([[6.0, 9.0], [12.0]])
    track_values = ak.Array([[0.1, -0.1], [0.2]])
    with uproot.recreate(root_path) as root_file:
        root_file["Delphes"] = {
            "Jet.PT": jet_values,
            "Jet.Flavor": ak.Array([[5, 21], [4]]),
            "Jet.Eta": jet_values,
            "Jet.Phi": jet_values,
            "Jet.Mass": jet_values,
            "Jet.NCharged": ak.Array([[2, 1], [1]]),
            "Jet.NNeutrals": ak.Array([[1, 1], [1]]),
            "Track.PT": track_values,
            "Track.Eta": track_values,
            "Track.Phi": track_values,
            "Track.Mass": track_values,
            "Track.Charge": ak.Array([[1, -1], [1]]),
            "Track.D0": track_values,
            "Track.DZ": track_values,
        }
    result = validate_root(root_path, SAMPLES["bbbar"], pt_min=5.0)
    assert result["jets_in_fiducial_region"] == 3
    assert result["truth_counts"] == {"b": 1, "c": 1, "g": 1, "uds": 0, "unknown": 0}


def test_artifact_lookup_handles_a_retry_with_a_new_madgraph_tag(tmp_path):
    banner = tmp_path / "run_01_tag_2_banner.txt"
    root = tmp_path / "tag_2_delphes_events.root"
    banner.touch()
    root.touch()
    assert locate_run_artifacts(tmp_path, "run_01") == (banner, root)
