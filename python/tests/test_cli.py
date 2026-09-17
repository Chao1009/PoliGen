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


# --------------------------------------------------------------------------
# `--isotope`: the spelling, and who typed it.
#
# `choices=` arrived on `--isotope` on 2026-09-16 so that an unknown species
# is one argparse line instead of a KeyError traceback.  On its own it also
# NARROWED what runs: MEASURED the same day against b647c16 (that commit's
# `python/lipolgen/` against this build's extension module),
# `--isotope 6li --channel tagged-d-p` exits 0 there -- a channel that implies
# its own species drops the typed one, so any casing reached a run -- and
# exited 2 with the bare `choices=`.  `cli._isotope_name` normalises the case
# before the check; these two tests are the pin.  What does NOT come back is a
# spelling `_ION_SPIN` does not carry at all: `zzz` now stops at argparse.

ISOTOPE_SPELLINGS = [("6li", "6Li"), ("6LI", "6Li"), ("6Li", "6Li"),
                     ("7li", "7Li"), ("7LI", "7Li"), ("d", "d"), ("D", "d"),
                     ("3he", "3He"), ("3He", "3He"), ("p", "p"), ("P", "p")]


@pytest.mark.parametrize("typed,canonical", ISOTOPE_SPELLINGS)
def test_every_isotope_casing_resolves_to_its_canonical_name(typed, canonical):
    assert cli.resolve(["--isotope", typed])["isotope"] == canonical
    # and it is the same isotope the run is built on
    assert lg.make_config(isotope=canonical, channel="inclusive",
                          events=2).isotope == canonical


@pytest.mark.parametrize("typed,resolved", [("6li", "6Li"), ("6LI", "6Li"),
                                            ("D", "d"), ("3he", "3He")])
def test_a_lower_case_isotope_still_runs(typed, resolved, capsys):
    """Exit status AND resolved isotope, on a channel that implies no species
    (the typed one is the run's) and on one that implies its own (the typed
    one is dropped, which is how these spellings ran before `choices=`)."""
    # `helicity-flip` is the one plan that takes any j, so the same command
    # line covers the spin-1/2 beam species as well as 6Li and d.
    assert cli.main(["--isotope", typed, "--channel", "inclusive",
                     "--plan", "helicity-flip",
                     "--events", "2", "--quiet"]) == 0
    assert cli.resolve(["--isotope", typed])["isotope"] == resolved
    capsys.readouterr()
    assert cli.main(["--isotope", typed, "--channel", "tagged-d-p",
                     "--events", "2", "--quiet"]) == 0


def test_an_unknown_isotope_is_argparses_one_line_refusal():
    for argv in (["--isotope", "zzz", "--channel", "tagged-d-p"],
                 ["--isotope", "zzz", "--channel", "inclusive"]):
        with pytest.raises(SystemExit) as e:
            cli.main(argv + ["--events", "2", "--quiet"])
        assert e.value.code == 2


def test_the_isotope_override_notice_needs_a_TYPED_isotope(capsys, tmp_path):
    """`DEFAULTS` fills 6Li, so `opts["isotope"] is not None` was true on every
    run and `lipolgen-run --channel tagged-d-p` announced an override of a
    species nobody typed (2026-09-16)."""
    assert cli.main(["--isotope", "7Li", "--channel", "tagged-d-p",
                     "--events", "2"]) == 0
    assert "--isotope 7Li overridden: --channel tagged-d-p implies d" \
        in capsys.readouterr().out
    assert cli.main(["--channel", "tagged-d-p", "--events", "2"]) == 0
    assert "overridden" not in capsys.readouterr().out
    # the config file types it too, under the same option name
    path = tmp_path / "iso.json"
    path.write_text(json.dumps({"isotope": "7Li", "channel": "tagged-d-p",
                                "events": 2}))
    assert cli.main(["--config-file", str(path)]) == 0
    assert "--isotope 7Li overridden" in capsys.readouterr().out
    assert cli.resolve([])["isotope"] == "6Li"          # the fill still fills
    assert cli.resolve([])["isotope_typed"] is False


# --------------------------------------------------------------------------
# The flags nothing typed.
#
# Seven `--` switches had ZERO occurrences anywhere under python/tests until
# 2026-09-16: --rc-c0-shape, --rc-a-transfer-frac, --rc-delta-high-x,
# --rel-lumi-offset, --run, --nthreads and --coherent-amp.  The knob matrix
# reaches the same KNOBS, but it reaches them through `make_config(**kw)` and
# never through `cli.resolve` / `cli.main`, so the argparse-to-config wiring
# of these seven was untested: severing all seven in an isolated checkout left
# the whole pytest suite green -- IDENTICAL to the control run in the same
# checkout, which is the whole of the finding.  These two tests type them.
#
# THE COUNT THAT STOOD HERE IS GONE (2026-09-16).  It read "1081 passed / 154
# skipped", and no reproducible state of this repository gives that pair: the
# tree it was measured on is not identified, and the experiment itself -- the
# seven severed, these two tests not yet written -- cannot be re-run on a tree
# that carries them.  What IS reproducible is the CONTROL tally, so that is
# what is quoted, with its recipe.  Copy the tracked tree with its
# working-tree content and no `.git`:
#     git ls-files -z | xargs -0 tar cf - | (cd <dir> && tar xf -)
# then `cmake -S <dir> -B <dir>/build -DLIPOLGEN_DEPS_PREFIX=$LIPOLGEN_DEPS`,
# `cmake --build <dir>/build -j`, and, with PYTHONPATH = $LIPOLGEN_DEPS/lib
# plus <dir>/build/python, `python -m pytest python/tests -q` from <dir>.
# MEASURED 2026-09-16 on this tree: **1103 passed / 155 skipped**, against
# 1107 / 151 in the tree itself -- 1258 collected either way.  The four that
# skip in the copy and run in the tree are the two git-checkout tests of
# `test_doc_link_gate.py`, `test_release_metadata.py`'s 'origin' remote test,
# and `test_hfs.py`'s polligen consumer, which looks for a sibling
# PolarizedLithiumSim checkout that a copy in a scratch directory has no
# sibling to find.

UNTYPED_FLAGS = [
    ("--rc-c0-shape", "vmc-ft", "rc_c0_shape", "vmc-ft", "ho"),
    ("--rc-a-transfer-frac", "0.5", "rc_a_transfer_frac", 0.5, 0.0),
    ("--rc-delta-high-x", "0.05", "rc_delta_high_x", 0.05, 0.015),
    ("--rel-lumi-offset", "0.02", "rel_lumi_offset", 0.02, 0.0),
    ("--run", "2", "run", 2, 1),
    ("--nthreads", "2", "nthreads", 2, 1),
]


def test_the_seven_untyped_flags_resolve():
    """Each of the seven reaches `resolve()` under its own key, and each falls
    back to the documented default when it is not typed."""
    argv = []
    for flag, text, _key, _want, _default in UNTYPED_FLAGS:
        argv += [flag, text]
    o = cli.resolve(argv + ["--coherent-amp", "0.02"])
    for flag, _text, key, want, _default in UNTYPED_FLAGS:
        assert key in o, "%s resolves to no option key" % flag
        assert o[key] == want, (flag, key, o[key], want)
    # --coherent-amp lands inside the `coherent` sub-dict, not at the top
    # level, which is itself part of the wiring nothing typed.
    assert o["coherent"] == {"amp": 0.02}

    bare = cli.resolve([])
    for flag, _text, key, _want, default in UNTYPED_FLAGS:
        assert bare[key] == default, (flag, key, bare[key], default)
    assert bare["coherent"] is None


def test_the_seven_untyped_flags_reach_the_run(tmp_path):
    """And each reaches the FILE: the meta records them and the npz moves.

    `resolve()` alone would not catch a switch that is parsed and then not
    passed to `make_config`, which is exactly the break the experiment above
    simulated.  Two runs, because the tensor-RC band knobs are REFUSED on the
    coherent channel (the band prices nothing there) and `--coherent-amp` is
    refused off it.
    """
    # (a) the five RC / bookkeeping knobs, on the channel the band applies to
    base, moved = tmp_path / "base.npz", tmp_path / "moved.npz"
    common = ["--channel", "inclusive", "--isotope", "6Li",
              "--plan", "tensor-thirds", "--events", "2000", "--seed", "4",
              "--rc", "tensor-band", "--quiet"]
    assert cli.main(common + ["--npz", str(base)]) == 0
    assert cli.main(common + [
        "--rc-c0-shape", "vmc-ft",
        "--rc-a-transfer-frac", "0.5",
        "--rc-delta-high-x", "0.05",
        "--rel-lumi-offset", "0.02",
        "--run", "2",
        "--npz", str(moved)]) == 0

    with np.load(str(base), allow_pickle=False) as f:
        m0 = json.loads(str(f["meta"]))
        x0 = np.asarray(f["x"])
    with np.load(str(moved), allow_pickle=False) as f:
        m1 = json.loads(str(f["meta"]))
        x1 = np.asarray(f["x"])

    assert m0["rc_c0_shape"] == "ho" and m1["rc_c0_shape"] == "vmc-ft"
    assert m0["run"] == 1 and m1["run"] == 2
    rows0 = m0["knob_provenance"]
    rows1 = m1["knob_provenance"]
    for name in ("rc_c0_shape", "rc_a_transfer_frac", "rc_delta_high_x",
                 "rel_lumi_offset"):
        assert name in rows1, name
        assert rows0[name]["value"] != rows1[name]["value"], name
    # ... and the EVENTS moved, which is the half `resolve()` cannot see:
    # --run changes the counter-based stream, so the kinematics differ.
    assert not np.array_equal(x0, x1)

    # (b) --coherent-amp, on the channel that reads it
    cbase, cmoved = tmp_path / "cbase.npz", tmp_path / "cmoved.npz"
    ccommon = ["--channel", "coherent", "--isotope", "6Li",
               "--plan", "tensor-thirds", "--events", "2000", "--seed", "4",
               "--quiet"]
    assert cli.main(ccommon + ["--npz", str(cbase)]) == 0
    assert cli.main(ccommon + ["--coherent-amp", "0.02",
                               "--npz", str(cmoved)]) == 0
    with np.load(str(cbase), allow_pickle=False) as f:
        c0 = json.loads(str(f["meta"]))
    with np.load(str(cmoved), allow_pickle=False) as f:
        c1 = json.loads(str(f["meta"]))
    r0 = c0["knob_provenance"]
    r1 = c1["knob_provenance"]
    assert "coherent_amp" in r1
    assert r0["coherent_amp"]["value"] != r1["coherent_amp"]["value"]


def test_nthreads_changes_no_number(tmp_path):
    """--nthreads is the one of the seven that must NOT move the file.

    `Pipeline` maps event index -> event as a pure function, so
    `generate(nthreads=8)` is bit-identical to `nthreads=1`; that is the
    promise `python/README.md` makes and the reason the flag is safe.
    """
    one, many = tmp_path / "t1.npz", tmp_path / "t8.npz"
    common = ["--channel", "tagged-6Li-alpha", "--events", "2000",
              "--seed", "6", "--quiet"]
    assert cli.main(common + ["--nthreads", "1", "--npz", str(one)]) == 0
    assert cli.main(common + ["--nthreads", "8", "--npz", str(many)]) == 0
    with np.load(str(one), allow_pickle=False) as a, \
            np.load(str(many), allow_pickle=False) as b:
        for key in ("x", "q2", "phi", "weight", "k", "cos_theta_k", "route"):
            assert np.array_equal(np.asarray(a[key]), np.asarray(b[key])), key
