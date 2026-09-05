"""D4 -- the `PDF:PomSet` systematic of the coherent T2 tier, pinned.

The full scan is 15 sets x 20 000 events (17 s wall, seeds 4242 / 777 /
31337); it lives in `docs/open_items/run_2026-09-03/phase_D_numbers.md` sec.
D4 and the band it produced is quoted in `pythia_bridge.hpp`, `--pom-set`'s
help, the run banner, `docs/USAGE.md` sec. 4 (Coherent 6Li -- it said "sec. 8"
until 2026-09-05; sec. 8 is "Command-line generators") and
`docs/OPEN_ITEMS_SOLUTIONS.md` sec. 5.  THE BAND'S WINDOW IS THE TWELVE DPDF
FITS (3-10, 12-15) ABOUT SET 6, not all fifteen: sets 1 and 2 are not Pomeron
fits and set 11 is refused.  THIS file re-runs a reduced version so that the
band cannot rot silently: the same four sets, 4 000 events each.

Two of the assertions are EXACT rather than statistical and are the ones that
matter most:

  * the T0 columns are bit-identical across sets -- so a `PomSet` band on |t|,
    x_P, M_X or the cross section is identically zero and quoting one would be
    meaningless (the observable had to move to the hadronic final state);
  * the light-only e_q^2 flavour fallback produces a bit-identical final
    state -- so the ~20 % fallback share of a default run is a bookkeeping
    artefact and NOT a systematic in its own right.
"""

import hashlib

import numpy as np
import pytest

import lipolgen as lg

pytestmark = pytest.mark.skipif(not lg._lipolgen.HAVE_PYTHIA8,
                                reason="built without PYTHIA 8")

N_EVENTS = 4000
SEED = 4242
#: default, the <n_ch> maximum, the <n_ch> minimum, the kaon maximum.
SETS = (6, 5, 9, 10)


def _run(pom_set, q2_pdf_min=None, include_charm=None, events=N_EVENTS):
    cfg = lg.make_config(isotope="6Li", config=1, channel="coherent",
                         events=events, seed=SEED)
    beams = lg.default_configs(cfg.isotope)[cfg.beam_config]
    po = lg.PythiaBridgeOptions()
    po.coherent_t2 = lg.CoherentT2.Pomeron
    po.pom_set = pom_set
    if q2_pdf_min is not None:
        po.q2_pdf_min = q2_pdf_min
    if include_charm is not None:
        po.include_charm = include_charm
    bridge = lg.PythiaBridge(beams, po)
    lg.set_pythia_hadronizer(cfg, bridge)
    cols = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6)).generate(
        0, True, 1)
    a = lg.hfs_arrays(cols["events"], False, False)
    off = np.asarray(a["offsets"])
    ch = np.asarray(a["charge"])
    pid = np.asarray(a["pid"])
    p4 = np.asarray(a["p4"])
    n_ev = len(off) - 1
    n_ch = np.array([np.count_nonzero(ch[off[i]:off[i + 1]] != 0)
                     for i in range(n_ev)], dtype=float)
    apid = np.abs(pid)
    kaon = np.isin(apid, (321, 311, 130, 310))
    return dict(
        meta=dict(cols["meta"]),
        t0=hashlib.md5(b"".join(
            np.ascontiguousarray(cols[c]).tobytes()
            for c in ("t", "x_pom", "q2", "x", "weight"))).hexdigest(),
        hfs=hashlib.md5(pid.tobytes()
                        + np.ascontiguousarray(p4).tobytes()).hexdigest(),
        mean_n_ch=float(n_ch.mean()),
        sem_n_ch=float(n_ch.std(ddof=1) / np.sqrt(n_ev)),
        mean_n_had=float(np.diff(off).mean()),
        mean_pt=float(np.hypot(p4[:, 1], p4[:, 2]).mean()),
        f_kaon=float(kaon.mean()),
        n_particles=int(pid.size),
    )


@pytest.fixture(scope="module")
def scan():
    return {s: _run(s) for s in SETS}


def test_the_pomeron_set_moves_nothing_at_t0(scan):
    """EXACT.  t, x_pom, q2, x and the event weight are one md5 over every
    set, so there is no PomSet band on |t|, x_P, M_X or sigma."""
    assert len({r["t0"] for r in scan.values()}) == 1


def test_the_pomeron_set_moves_the_whole_hadronic_final_state(scan):
    assert len({r["hfs"] for r in scan.values()}) == len(SETS)


def test_the_charged_multiplicity_band(scan):
    """<n_charged>: measured -2.6 % / +9.9 % about set 6 over the twelve DPDF
    fits, with set 5 the maximum and set 9 the minimum.  Reduced statistics,
    so only the ordering and a conservative floor on the spread are pinned."""
    hi, lo, dflt = scan[5], scan[9], scan[6]
    assert hi["mean_n_ch"] > dflt["mean_n_ch"] > 0.0
    assert hi["mean_n_ch"] - lo["mean_n_ch"] > 0.25      # measured 0.495
    # ... which is many sigma, not a fluctuation
    sem = np.hypot(hi["sem_n_ch"], lo["sem_n_ch"])
    assert hi["mean_n_ch"] - lo["mean_n_ch"] > 4.0 * sem
    assert hi["mean_n_ch"] / dflt["mean_n_ch"] - 1.0 > 0.05   # measured +9.9 %


def test_the_kaon_fraction_is_the_largest_mover(scan):
    """The secondary observable: a FACTOR 2.8 end to end (set 9 to set 10),
    the one most likely to matter for a PID-based analysis."""
    assert scan[10]["f_kaon"] / scan[9]["f_kaon"] > 2.2
    assert scan[10]["f_kaon"] > scan[6]["f_kaon"] > scan[9]["f_kaon"]


def test_meta_says_which_set_made_the_file(scan):
    for s in SETS:
        assert scan[s]["meta"]["pom_set"] == s
        assert scan[s]["meta"]["coherent_t2"] == "pomeron"
        assert scan[s]["meta"]["t2_bridge"] == "pythia8"
    metas = [scan[s]["meta"] for s in SETS]
    for i in range(len(metas)):
        for j in range(i + 1, len(metas)):
            assert metas[i] != metas[j]


def test_the_fallback_counter_is_on_the_record(scan):
    m = scan[6]["meta"]
    assert 0.15 < m["pom_flavour_fallback_frac"] < 0.25       # measured 0.207
    assert m["n_pom_flavour_fallback"] == pytest.approx(
        m["pom_flavour_fallback_frac"] * m["t2_n_pomeron"], abs=1.0)
    assert scan[5]["meta"]["n_pom_flavour_fallback"] == 0


def test_the_light_only_fallback_is_not_a_systematic():
    """EXACT.  Raising `q2_pdf_min` from 1.0 to 1.75 takes set 6's fallback
    share from ~20 % to 0 % and leaves the final state BIT-IDENTICAL: every
    Pomeron DPDF carries a single light-quark singlet, so e_q^2 x f_q is
    proportional to e_q^2 exactly and the 'fallback' IS the true draw."""
    base = _run(6, events=2000)
    lifted = _run(6, q2_pdf_min=1.75, events=2000)
    assert base["meta"]["n_pom_flavour_fallback"] > 0
    assert lifted["meta"]["n_pom_flavour_fallback"] == 0
    assert base["hfs"] == lifted["hfs"]
    assert base["n_particles"] == lifted["n_particles"]


def test_pom_set_11_is_refused():
    """It needs PYTHIA's setXPom, which the bridge never calls: measured,
    100.00 % of events took the e_q^2 fallback and the Pomeron PDF was never
    consulted, so `meta["pom_set"] = 11` would record a knob that did not
    run."""
    beams = lg.default_configs("6Li")[1]
    po = lg.PythiaBridgeOptions()
    po.coherent_t2 = lg.CoherentT2.Pomeron
    po.pom_set = 11
    with pytest.raises(RuntimeError, match="PomSet = 11"):
        lg.PythiaBridge(beams, po)
    # ... but only when the Pomeron instance is actually built.
    po.coherent_t2 = lg.CoherentT2.Off
    lg.PythiaBridge(beams, po)


def test_the_cli_reports_a_bridge_refusal_as_a_command_line_error():
    """A `--pom-set` the bridge refuses is a configuration error typed at the
    command line, so the CLI shows the bridge's own message rather than a
    traceback (and does not restate it -- one definition)."""
    from lipolgen import cli
    with pytest.raises(SystemExit) as ei:
        cli.main(["--channel", "coherent", "--events", "20", "--hadronize",
                  "--pom-set", "11", "--quiet"])
    assert "PomSet = 11" in str(ei.value)
    assert "did not run" in str(ei.value)
