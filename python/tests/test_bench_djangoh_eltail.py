# SPDX-License-Identifier: GPL-3.0-or-later
"""T4 row 7 of docs/benchmarking/BENCHMARK_PLAN.md, re-wired 2026-10-26 (run
2026-10-26, phase B3): DJANGOH 4.6.21's ISOLATED elastic radiative tail
(IEL31 + IEL32 + IEL33, the matched window) against LiPolGen's `rc_tail` on
a proton at e 18 x p 275 GeV -- `validation/benchmarks/t4_djangoh_rad_noRad.py`,
reference vendored in `validation/benchmarks/data/djangoh_eltail_ep18x275.json`
(the ePIC Rad/noRad table, `djangoh_rad_norad_ep18x275.json`, stays as the
recorded history; its own guards are in test_bench_chain_rc.py).

These tests pin what the harness MEASURES and that its verdict follows from
the numbers by the rule its header fixed BEFORE the DJANGOH runs (the design
block's sha256 is pinned too).  No test asserts a tree default: the t-peak
option's (the shipped default until D5, 2026-10-27) 39 % / 10 % shortfall at
Q^2 = 1-10 / 10-100 is recorded as measured; the table's +2 % is pinned.

The full row costs 22-26 s (2026-10-26, a shared machine), under the 30 s line,
so it runs always; the table-vs-quadrature diagnostic (~110 s) is opt-in
(`bench_expensive`, LIPOLGEN_BENCH=1).  Nothing here needs the network.
"""
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys

import pytest

lg = pytest.importorskip("lipolgen")

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(
    __file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")
_NAME = "t4_djangoh_rad_noRad"
_ELTAIL = "djangoh_eltail_ep18x275.json"
_EPIC = "djangoh_rad_norad_ep18x275.json"

#: measured 2026-10-26 (phase_B3_djangoh.md); the harness's own numbers
R_POLRAD_FULL = (1.01989075, 1.02613455, 1.02191885, 1.00978744)
R_T_PEAK = (0.61181928, 0.89680803, 1.0111558, 1.00674599)
DJ_FRACTION = (0.04072928, 0.0164353, 0.01709873, 0.02188814)
DJ_TAIL_NB = (23.3808849182754, 1.111466620112, 0.05516243362021,
              0.00151825766947284)
DJ_BORN_NB = (574.056, 67.6268, 3.226113, 0.069364414)
#: the 2026-09-23 history numbers (phase_B3_chain_rc.md sec. 1.3), which the
#: re-wired harness must still reproduce on the W_l^2 >= 9 window
HIST_T_PEAK = (1.01616918, 1.01317215, 1.01497026, 1.02068610)
HIST_POLRAD_FULL = (1.02709574, 1.01507182, 1.01513763, 1.02075698)
EPIC_RATIO = (1.07245359, 1.08940691, 1.16062561, 1.2279651)


def _load():
    spec = importlib.util.spec_from_file_location(
        "bench_" + _NAME, os.path.join(_BENCH, _NAME + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_harness(script, *args, cwd=None):
    """The harness as a script, with the imported lipolgen on its path."""
    env = dict(os.environ)
    pkg_parent = os.path.dirname(os.path.dirname(os.path.abspath(lg.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(
        [pkg_parent] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep)
                        if p])
    return subprocess.run([sys.executable, script] + list(args),
                          capture_output=True, text=True, env=env, cwd=cwd,
                          timeout=600)


def _copy_bench(tmp_path):
    """The harness and both vendored files, copied into a fresh tree."""
    d = tmp_path / "bench"
    (d / "data").mkdir(parents=True)
    shutil.copy(os.path.join(_BENCH, _NAME + ".py"), d / (_NAME + ".py"))
    for f in (_ELTAIL, _EPIC):
        shutil.copy(os.path.join(_BENCH, "data", f), d / "data" / f)
    return d


@pytest.fixture(scope="module")
def dj():
    return _load()


@pytest.fixture(scope="module")
def row(dj):
    return dj.run(verbose=False)


# ------------------------------------------------------------ file guards

def test_eltail_file_guards(dj, tmp_path):
    """Point counts, first and last value, and ONE flipped byte refused."""
    doc, bins, control = dj.load_eltail(m_p=lg.PROTON_MASS)
    assert len(bins) == 4 and len(control) == 8
    assert len(doc["runs"]) == 40 and len(doc["cards"]) == 40
    assert [(b["q2_lo"], b["q2_hi"]) for b in bins] == [
        (1.0, 10.0), (10.0, 100.0), (100.0, 1000.0), (1000.0, 10000.0)]
    # first and last value (recomputed from the raw run records)
    assert bins[0]["born_nb"] == pytest.approx(DJ_BORN_NB[0], rel=1e-12)
    assert bins[-1]["tail_novp_nb"] == pytest.approx(DJ_TAIL_NB[-1], rel=1e-12)
    for b, t, bo in zip(bins, DJ_TAIL_NB, DJ_BORN_NB):
        assert b["tail_novp_nb"] == pytest.approx(t, rel=1e-12)
        assert b["born_nb"] == pytest.approx(bo, rel=1e-12)
    raw = open(dj.ELTAIL_DATA, "rb").read()
    assert hashlib.sha256(raw).hexdigest() == dj.ELTAIL_SHA256
    old = b'"born_nb": 574.056,'
    assert raw.count(old) == 1
    bad = tmp_path / _ELTAIL
    bad.write_bytes(raw.replace(old, b'"born_nb": 574.057,'))
    with pytest.raises(RuntimeError, match="sha256"):
        dj.load_eltail(bad)
    good = tmp_path / "good.json"
    good.write_bytes(raw)
    dj.load_eltail(good)


def test_eltail_provenance_and_every_card_verbatim(dj):
    doc, _bins, _control = dj.load_eltail()
    prov = doc["provenance"]
    for key in ("generator", "source_repository", "source_commit",
                "source_licence", "source_files_sha256", "manual", "build",
                "pdf", "epic_templates", "design", "timeline_utc", "papers",
                "publishing_convention"):
        assert prov[key], key
    assert prov["source_commit"].startswith(
        "e356b847bcacf2f1596ac4ea42c677cef0614176")
    assert "10150" in prov["pdf"]["set"]
    assert prov["epic_templates"]["commit"] == (
        "9869d9a0d21b9191a1bbec7294b4851feb5ee20a")
    # every run's card is in the file, byte for byte what was run
    assert set(doc["cards"]) == set(doc["runs"])
    for name, card in doc["cards"].items():
        r = doc["runs"][name]
        assert hashlib.sha256(card.encode("utf-8")).hexdigest() == \
            r["card_sha256"], name
        assert r["exit"] == 0 and r["stdout_sha256"] and r["out_sha256"], name
        assert "RNDM-SEEDS\n            0 0" in card, name
        assert "START\n            0" in card, name


def test_isolating_cards_are_the_declared_design(dj):
    """Design block (a)-(d), read off the cards that ran: IEL2 = 0, INC off,
    IEL31..33 on, WMIN 0, YMAX 0.985, Q2MAX inactive, the x ceilings, the
    fixed-alpha GSW-PARAM -- and the ePIC settings in the control."""
    doc, _bins, _control = dj.load_eltail()
    xmax = {"1_10": "0.1096D0", "10_100": "0.5518D0", "100_1000": "0.9249D0",
            "1000_10000": "0.95D0"}
    n = 0
    for name, card in doc["cards"].items():
        lines = [l.strip() for l in card.split("\n")]

        def val(code):
            return lines[lines.index(code) + 1]
        if name.startswith("control_"):
            rad = name.split("_")[1]
            assert val("INT-OPT-NC").split() == (
                ["1", "18", "18", "18", "18", "0", "0", "0", "0"]
                if rad == "Rad1" else ["1"] + ["0"] * 8)
            assert val("KINEM-CUTS").split()[-1] == "3.0D0"     # WMIN = 3
            assert val("INT-ONLY") == "1" and val("INT-POINTS") == "100000"
            continue
        kind, tag = name.split("_")[0], [t for t in xmax if "_%s_" % t in name][0]
        cuts = val("KINEM-CUTS").split()
        assert cuts[0] == "3" and cuts[1] == "1D-5" and cuts[2] == xmax[tag]
        assert cuts[3] == "1D-4" and cuts[4] == "0.985D0"
        assert cuts[7] == "0.0D0"                              # WMIN
        assert val("INT-ONLY") == "-1" and val("INT-POINTS") == "100000"
        opt = val("INT-OPT-NC").split()
        if kind == "born":
            assert opt == ["1"] + ["0"] * 8
            assert val("GSW-PARAM") == "2 0 3 1 0 0 2 1 1 1 1"
        else:
            assert opt == ["0"] * 6 + ["20", "20", "20"]      # IEL2 = 0
            assert val("GSW-PARAM") == ("2 1 2 1 1 1 2 1 1 1 1" if "_novp_" in
                                        name else "2 1 3 1 1 1 2 1 1 1 1")
        if "_direct_q2max" not in name:
            assert cuts[6] == "100000000D0"                   # Q2MAX off
        n += 1
    assert n == 32


def test_design_block_was_fixed_before_the_runs(dj):
    """The design block hashes to the recorded digest, the vendored file
    records the same digest, and every DJANGOH run started after the block
    was written; one edited character in the block is refused."""
    assert dj.design_block_sha256() == dj.DESIGN_BLOCK_SHA256
    doc, _bins, _control = dj.load_eltail()
    design = doc["provenance"]["design"]
    assert design["design_block_sha256"] == dj.DESIGN_BLOCK_SHA256
    t0 = design["design_written_utc"]
    assert t0 == "2026-10-26T19:11:53Z"
    for name, r in doc["runs"].items():
        assert r["started_utc"].replace("+00:00", "Z") > t0, name
    src = open(os.path.join(_BENCH, _NAME + ".py"), encoding="utf-8").read()
    old = "PASS(polrad-full)  iff  |R_b - 1| <= 0.03 + 2 delta_b"
    assert src.count(old) == 1
    with pytest.raises(RuntimeError, match="design block"):
        dj.design_block_sha256(src.replace(old, old.replace("0.03", "0.05")))


def test_window_rule_and_matched_cuts(dj):
    """XMAX_b is the declared rule, and the whole rectangle has W_l^2 >= 9
    (its corner x = XMAX_b, Q^2 = Q2_lo is the minimum)."""
    m_p = lg.PROTON_MASS
    got = [dj.xmax_rule(lo, m_p) for lo in (1.0, 10.0, 100.0, 1000.0)]
    assert got == [0.1096, 0.5518, 0.9249, 0.95]
    for lo, x in zip((1.0, 10.0, 100.0, 1000.0), got):
        assert lo * (1.0 / x - 1.0) + m_p * m_p >= 9.0


def test_control_reproduces_the_epic_samples(dj, row):
    _doc, _bins, control = dj.load_eltail()
    assert len(control) == 8
    assert all(c["reproduced"] for c in control)
    assert max(abs(c["pull"]) for c in control) == pytest.approx(1.7417, abs=1e-3)
    readme = {(c["q2_lo"], c["rad"]): c["epic_readme_nb"] for c in control}
    assert readme[(10.0, 1)] == pytest.approx(74.05000239, rel=1e-12)
    assert readme[(1000.0, 0)] == pytest.approx(0.07015936779, rel=1e-12)
    ctl = [s for s in row["subrows"] if s["name"].startswith("control[")]
    assert len(ctl) == 8 and all(s["status"] == "pass" for s in ctl)
    assert not any(s["in_verdict"] for s in ctl)


def test_q2max_is_ignored_by_the_tail_and_honoured_by_the_born(dj):
    """Design block (c)'s reason for the subtraction, as measured."""
    _doc, bins, _control = dj.load_eltail()
    for b in bins:
        assert b["tail_novp_direct_q2max_nb"] == b["tail_novp_q2min_lo_nb"]
        assert b["born_direct_q2max_nb"] / b["born_nb"] == pytest.approx(
            1.0, abs=3e-5)


# ------------------------------------------------------------ the row

def test_row_passes_by_the_declared_rule(dj, row):
    assert row["name"] == _NAME
    assert row["status"] == "pass"
    for d in row["detail"]:
        allow = dj.TOL_NUMERICAL + dj.N_DELTA * d["djangoh_tail_rel_err"]
        assert d["allowance"] == pytest.approx(allow, rel=1e-12)
        r_full, r_tpk = d["polrad-full"]["R"], d["t-peak"]["R"]
        assert d["polrad-full"]["pass"] == (abs(r_full - 1.0) <= allow)
        assert d["t-peak"]["pass"] == (r_tpk <= 1.0 + allow)
        # R is the same whether read as tails or as fractions on ONE Born
        assert d["polrad-full"]["fraction_on_djangoh_born"] / \
            d["djangoh_fraction"] == pytest.approx(r_full, rel=1e-12)
    ver = [s for s in row["subrows"] if s["in_verdict"]]
    assert len(ver) == 8 and all(s["status"] == "pass" for s in ver)


def test_tolerance_constants_match_the_hashed_block(dj, row):
    """Verification of 2026-10-26: TOL_NUMERICAL and N_DELTA sat outside the
    hashed block.  Now the block's two "PASS(...) iff" sentences are parsed
    on every run and must equal the constants; the allowances are pinned."""
    src = open(os.path.join(_BENCH, _NAME + ".py"), encoding="utf-8").read()
    i0 = src.index(dj.DESIGN_BEGIN)
    i1 = src.index(dj.DESIGN_END) + len(dj.DESIGN_END)
    assert dj.tolerance_from_design_block(src[i0:i1]) == (0.03, 2.0)
    assert (dj.TOL_NUMERICAL, dj.N_DELTA) == (0.03, 2.0)
    assert [d["allowance"] for d in row["detail"]] == pytest.approx(
        [0.0307467, 0.0309598, 0.0308831, 0.0307745], abs=5e-8)
    assert "0.03 + 2 delta_b" in row["tolerance"]
    assert "one-sided" in row["tolerance"] and "CEILING" in row["reason"]


def test_row_numbers_are_the_measured_ones(row):
    got = row["generator_value"]
    for g, w in zip(got["polrad-full"], R_POLRAD_FULL):
        assert g == pytest.approx(w, rel=2e-7)
    for g, w in zip(got["t-peak"], R_T_PEAK):
        assert g == pytest.approx(w, rel=2e-7)
    assert row["reference_value"] == pytest.approx(list(DJ_FRACTION), rel=1e-7)
    # the t-peak option's (shipped until D5) shortfall at low Q^2, recorded
    # (not asserted as a defect): 39 % and 10 % missing in bins 1 and 2
    det = row["detail"]
    assert 1.0 - det[0]["t-peak"]["R"] == pytest.approx(0.3882, abs=1e-4)
    assert 1.0 - det[1]["t-peak"]["R"] == pytest.approx(0.1032, abs=1e-4)
    # DJANGOH's running alpha (recorded sub-row): +1.3 ... +2.1 %
    vp = [d["running_alpha_over_fixed"] for d in det]
    assert vp == pytest.approx([1.0151737, 1.0132801, 1.0151167, 1.0206271],
                               rel=1e-6)
    # nothing outside the tail table; the clip bites only in bin 4, 5.5e-4
    for d in det:
        for tm in ("t-peak", "polrad-full"):
            assert d[tm]["uncovered_born_fraction"] == 0.0
    assert det[3]["polrad-full"]["clipped_tail_fraction"] == pytest.approx(
        5.46e-4, rel=1e-2)
    cfg = row["config"]
    assert cfg["xmax_b"] == [0.1096, 0.5518, 0.9249, 0.95]
    assert cfg["y_max"] == pytest.approx(0.985)
    assert cfg["x_min"] == 1e-5 and cfg["y_min"] == 1e-4
    # D5 (2026-10-27): the LABELS follow the shipped default; the hashed
    # rules name the models, so the numbers and the verdict did not move.
    assert "polrad-full (shipped)" in row["generator"]
    assert "t-peak (option)" in row["generator"]


def test_history_reproduces_2026_09_23(row):
    hist = [s for s in row["subrows"] if s["name"].startswith("history_epic[")]
    assert len(hist) == 4
    for s, t, f, e in zip(hist, HIST_T_PEAK, HIST_POLRAD_FULL, EPIC_RATIO):
        assert s["lipolgen_t-peak_ratio"] == pytest.approx(t, rel=2e-7)
        assert s["lipolgen_polrad-full_ratio"] == pytest.approx(f, rel=2e-7)
        assert s["epic_rad_over_norad"] == pytest.approx(e, rel=1e-8)
        assert s["status"] == "recorded" and not s["in_verdict"]


def test_report_line_and_normalised_row(dj, row, capsys):
    dj.print_row(row)
    out = capsys.readouterr().out
    lines = [l for l in out.splitlines() if l.startswith("REPORT |")]
    assert len(lines) == 1 and lines[0].startswith("REPORT | %s |" % _NAME)
    assert "| PASS --" in lines[0]
    assert "np.float64" not in out
    assert out.count("SUBROW | %s[" % _NAME) == len(row["subrows"]) == 24
    sys.path.insert(0, _BENCH)
    try:
        import _report
    finally:
        sys.path.remove(_BENCH)
    norm = _report.normalise_row(row, _NAME)
    assert norm["status"] == "pass" and norm["notes"] == []
    assert len(norm["subrows"]) == 24
    json.dumps(_report.jsonable(norm))                  # strict JSON


def test_plot_writes_a_png(dj, row, tmp_path):
    pytest.importorskip("matplotlib")
    path = dj.plot(row, str(tmp_path))
    assert os.path.dirname(path) == str(tmp_path)
    raw = open(path, "rb").read()
    assert len(raw) > 1000 and raw[:8] == b"\x89PNG\r\n\x1a\n"


# ------------------------------------------------------------ exit codes

def test_exit_0_on_a_run():
    p = _run_harness(os.path.join(_BENCH, _NAME + ".py"))
    assert p.returncode == 0, p.stderr
    assert "REPORT | %s |" % _NAME in p.stdout


def test_exit_2_on_an_altered_file(tmp_path):
    d = _copy_bench(tmp_path)
    f = d / "data" / _ELTAIL
    raw = f.read_bytes()
    old = b'"born_nb": 574.056,'
    assert raw.count(old) == 1
    f.write_bytes(raw.replace(old, b'"born_nb": 574.057,'))
    p = _run_harness(str(d / (_NAME + ".py")))
    assert p.returncode == 2 and "sha256" in p.stderr


def test_exit_2_on_an_edited_design_block(tmp_path):
    d = _copy_bench(tmp_path)
    h = d / (_NAME + ".py")
    src = h.read_text(encoding="utf-8")
    old = "PASS(polrad-full)  iff  |R_b - 1| <= 0.03 + 2 delta_b"
    assert src.count(old) == 1
    h.write_text(src.replace(old, old.replace("0.03", "0.05")),
                 encoding="utf-8")
    p = _run_harness(str(h))
    assert p.returncode == 2 and "design block" in p.stderr


def test_exit_2_on_an_edited_tolerance_constant(tmp_path):
    """The constants the verdict uses are read back from the hashed block on
    every run (they sat outside it until the verification of 2026-10-26)."""
    d = _copy_bench(tmp_path)
    h = d / (_NAME + ".py")
    src = h.read_text(encoding="utf-8")
    old = "\nTOL_NUMERICAL = 0.03\n"
    assert src.count(old) == 1
    h.write_text(src.replace(old, "\nTOL_NUMERICAL = 0.05\n"), encoding="utf-8")
    p = _run_harness(str(h))
    assert p.returncode == 2 and "tolerance constants" in p.stderr


def test_exit_2_without_numpy(tmp_path):
    """A missing dependency is a broken harness (exit 2), numpy included --
    its import is the first thing the module does."""
    fake = tmp_path / "numpy"
    fake.mkdir()
    (fake / "__init__.py").write_text(
        'raise ImportError("simulated: numpy not installed")\n')
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(tmp_path)] + [p for p in [env.get("PYTHONPATH")] if p])
    p = subprocess.run([sys.executable, os.path.join(_BENCH, _NAME + ".py")],
                       env=env, capture_output=True, text=True, timeout=120)
    assert p.returncode == 2, p.stderr
    assert "simulated: numpy not installed" in p.stderr


# ------------------------------------------------------------ opt-in

@pytest.mark.bench_expensive
def test_diagnostics_table_vs_quadrature(dj):
    """`--diagnostics` (~110 s): the exact polrad-full tail straight from the
    quadrature agrees with DJANGOH to <= 1.4 %; the shipped table lifts the
    window integral by +1.9 ... +2.5 % (both models), through its TAIL
    column's interpolation (x1.020 ... x1.037), not its Born's; in bin 4 the
    tail_max clip takes 1.3 % of polrad-full's back.  Recorded, never in the
    verdict."""
    out = dj.diagnostics()
    full = out["polrad-full"]["bins"]
    tpk = out["t-peak"]["bins"]
    assert [b["table_over_direct"] for b in full] == pytest.approx(
        [1.01949, 1.02508, 1.02242, 1.02349], abs=2e-5)
    assert [b["R_direct"] for b in full] == pytest.approx(
        [1.00039, 1.00103, 0.99951, 0.98661], abs=2e-5)
    assert [b["tail_interp_factor"] for b in full] == pytest.approx(
        [1.02029, 1.02517, 1.02259, 1.03652], abs=2e-5)
    assert [b["born_interp_and_clip_factor"] for b in full] == pytest.approx(
        [0.99922, 0.99991, 0.99984, 0.98743], abs=2e-5)
    assert [b["table_over_direct"] for b in tpk] == pytest.approx(
        [1.02096, 1.02171, 1.02190, 1.02244], abs=2e-5)
    assert [b["R_direct"] for b in tpk] == pytest.approx(
        [0.59926, 0.87775, 0.98948, 0.98465], abs=2e-5)
    assert [b["tail_interp_factor"] for b in tpk] == pytest.approx(
        [1.02169, 1.02181, 1.02206, 1.02242], abs=2e-5)
    assert out["polrad-full"]["quadrature_max_rel_change"] == pytest.approx(
        1.80e-4, abs=5e-6)
    assert out["t-peak"]["quadrature_max_rel_change"] == pytest.approx(
        9.23e-5, abs=5e-7)
