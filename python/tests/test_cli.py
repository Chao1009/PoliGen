"""`lipolgen-run` -- option resolution and the outputs it writes."""

import json

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import cli


def test_defaults_and_switch_precedence():
    o = cli.resolve(["--isotope", "7Li", "--events", "1234"])
    assert o["isotope"] == "7Li"
    assert o["events"] == 1234
    assert o["channel"] == "inclusive"          # default survived
    assert o["seed"] == cli.DEFAULTS["seed"]


def test_config_file_then_switches(tmp_path):
    path = tmp_path / "run.json"
    path.write_text(json.dumps({"isotope": "7Li", "channel": "tagged-alpha",
                                "events": 500, "optics": "tagging"}))
    o = cli.resolve(["--config-file", str(path), "--events", "77"])
    assert o["isotope"] == "7Li"
    assert o["channel"] == "tagged-alpha"
    assert o["optics"] == "tagging"
    assert o["events"] == 77                    # the switch wins over the file


def test_lumi_alone_selects_luminosity_mode(tmp_path):
    """`--lumi X` alone must not collide with the `events` DEFAULT.

    Luminosity mode used to be reachable only as `--events 0 --lumi X`: the
    option defaults carry events = 100000, so `--lumi X` alone arrived at
    `main()` with both set and was refused with "--events and --lumi are
    exclusive" -- a message about a switch the user never typed, and a mode no
    help text or document described (fixed 2026-09-05).
    """
    o = cli.resolve(["--lumi", "0.02"])
    assert o["lumi"] == 0.02
    assert o["events"] == 0                 # luminosity mode, not a collision

    # an --events typed anywhere still collides, and is still refused there
    o = cli.resolve(["--lumi", "0.02", "--events", "500"])
    assert (o["events"], o["lumi"]) == (500, 0.02)
    path = tmp_path / "run.json"
    path.write_text(json.dumps({"events": 500}))
    o = cli.resolve(["--config-file", str(path), "--lumi", "0.02"])
    assert (o["events"], o["lumi"]) == (500, 0.02)

    # ... and neither flag leaves the defaults alone
    assert cli.resolve([])["events"] == cli.DEFAULTS["events"]
    assert cli.resolve(["--events", "77"])["lumi"] == cli.DEFAULTS["lumi"]


def test_config_file_rejects_typos(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps({"isotop": "6Li"}))
    with pytest.raises(SystemExit):
        cli.resolve(["--config-file", str(path)])


def test_coherent_scenario_switches():
    o = cli.resolve(["--channel", "coherent", "--coherent-f0", "0.06",
                     "--coherent-slope-b", "42"])
    assert o["coherent"] == {"f0": 0.06, "slope_b": 42.0}


def test_yaml_config_file(tmp_path):
    yaml = pytest.importorskip("yaml")
    path = tmp_path / "run.yaml"
    path.write_text(yaml.safe_dump({"channel": "tagged-alpha", "events": 400}))
    o = cli.resolve(["--config-file", str(path)])
    assert o["channel"] == "tagged-alpha" and o["events"] == 400


def test_events_and_lumi_are_exclusive():
    with pytest.raises(SystemExit):
        cli.main(["--events", "10", "--lumi", "1.0"])


def test_full_run_writes_npz_and_hepmc(tmp_path, capsys):
    npz = tmp_path / "out.npz"
    argv = ["--isotope", "6Li", "--config", "1", "--channel", "tagged-alpha",
            "--plan", "tensor-thirds", "--events", "3000", "--seed", "1",
            "--npz", str(npz), "--quiet"]
    if lg._lipolgen.HAVE_HEPMC3:
        argv += ["--hepmc", str(tmp_path / "out.hepmc")]
    assert cli.main(argv) == 0
    with np.load(str(npz), allow_pickle=False) as f:
        assert f["x"].size == 3000
        meta = json.loads(str(f["meta"]))
    assert meta["channel"] == "tagged-6Li-alpha"
    assert meta["seed"] == 1
    if lg._lipolgen.HAVE_HEPMC3:
        assert (tmp_path / "out.hepmc").stat().st_size > 0


def test_hfs_npz_needs_hadronize(tmp_path):
    with pytest.raises(SystemExit):
        cli.main(["--events", "10", "--hfs-npz", str(tmp_path / "h.npz"),
                  "--quiet"])


@pytest.mark.skipif(not lg._lipolgen.HAVE_PYTHIA8,
                    reason="built without PYTHIA 8")
def test_hadronize_writes_an_hfs_sample(tmp_path):
    path = tmp_path / "hfs.npz"
    assert cli.main(["--channel", "inclusive", "--events", "100", "--seed", "8",
                     "--hadronize", "--hfs-npz", str(path), "--quiet"]) == 0
    s = lg.export.load_hfs_npz(str(path))
    assert s["offsets"].size == 101
    assert s["pid"].size > 100                 # hadrons really were written
    assert s["meta"]["sigma_gen_mb"] > 0


def test_7li_rejects_the_spin1_thirds_plan():
    with pytest.raises(ValueError):
        lg.make_plan("tensor-thirds", j=1.5)


def test_plan_names_all_build():
    for name in ("tensor-thirds", "helicity-flip", "transverse-tensor",
                 "tensor-flip"):
        plan = lg.make_plan(name, j=1.0, pz=0.7, pzz=0.6, pe=0.7)
        assert len(plan) >= 1
        assert abs(sum(c.lumi_fraction for c in plan.categories) - 1.0) < 1e-9


def test_cluster_wave_flag_reaches_the_config():
    opts = cli.resolve(["--channel", "tagged-alpha", "--cluster-wave", "vmc"])
    assert opts["cluster_wave"] == "vmc"
    assert cli.resolve(["--channel", "tagged-alpha"])["cluster_wave"] == "hulthen"
