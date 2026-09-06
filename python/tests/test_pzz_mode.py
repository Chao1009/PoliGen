# SPDX-License-Identifier: GPL-3.0-or-later
"""`--pzz-mode` -- honouring the typed `--pzz` on the `helicity-flip` plan.

REGISTRY ROW 20 / `OPEN_ITEMS_SOLUTIONS.md` sec. 15.5 **D13**, priced in
`docs/open_items/run_2026-09-06/phase_A_numbers.md` sec. A3.  Until 2026-09-06
`--plan helicity-flip --pzz 0.6` at `--pz 0.7` filled the max-entropy ladder's
own alignment (T = 0.4 at J = 3/2) and the typed 0.6 was NOT READ; the run
banner said so, but no switch honoured it.  `--pzz-mode typed` is that switch,
and `--pzz-mode ladder` -- the DEFAULT -- is the branch that was always there.

WHAT THIS FILE GATES, in the order the decision needs it:

  P1  the default is the ladder, and it is the SAME OBJECT as the pre-flag
      call -- populations double for double, not to a tolerance;
  P2  `typed` honours what was typed, at both spins, and leaves the VECTOR
      moment alone (so A_par's own divisor does not move);
  P3  outside the plan's domain it is REFUSED with the edge named, never
      clamped -- at `--pz 0.7` the J = 3/2 domain ends at T = 0.58 and the
      CLI's own default `--pzz 0.6` is outside it by 0.02, while the same
      0.6 is INSIDE the spin-1 domain and runs;
  P4  the measured price at the standard configuration, pinned;
  P5  `knob_provenance` -- the mode's own row, and the `--pzz` row that
      flips from not-read to read with it;
  P6  the CLI reaches all of it, and its refusal is a SystemExit carrying
      the edge rather than a traceback.
"""
import numpy as np
import pytest

import lipolgen as lg
from lipolgen import _lipolgen as _l
from lipolgen import cli

PZ, PE, SEED = 0.7, 0.7, 20260713


def _plan(j, pzz=0.5, mode="ladder", **kw):
    return lg.make_plan("helicity-flip", j=j, pz=PZ, pzz=pzz, pe=PE,
                        pzz_mode=mode, **kw)


def _pipeline(iso, j, pzz=0.5, mode="ladder", events=200):
    cfg = lg.make_config(isotope=iso, channel="inclusive", config=1,
                         events=events, seed=SEED)
    return _l.Pipeline(cfg, _plan(j, pzz, mode))


# ------------------------------------------------------------------- P1

@pytest.mark.parametrize("j", [1.0, 1.5])
def test_the_default_is_the_ladder_double_for_double(j):
    """P1.  Not "agrees to 1e-12": the SAME doubles.

    The default takes the branch that was always there -- `make_plan` does
    not touch `HelicityFlipOptions` unless the mode says `typed` -- so this
    is bit for bit by construction, and the test says which construction.
    """
    default = _plan(j)
    named = _plan(j, mode="ladder")
    raw = _l.helicity_flip_plan(j, PZ, PE)           # the pre-flag call
    for other in (named, raw):
        assert (list(default.categories[0].populations)
                == list(other.categories[0].populations))
        assert default.pzz_true == other.pzz_true
        assert default.pz_true == other.pz_true
    # ... and that fill IS the max-entropy ladder, not merely close to it
    assert (list(default.categories[0].populations)
            == list(_l.populations_maxent(j, PZ)))
    assert default.pzz_true == _l.spin_temperature_pzz(j, PZ)


@pytest.mark.parametrize("j,iso", [(1.0, "6Li"), (1.5, "7Li")])
def test_the_default_run_is_bit_identical_to_the_pre_flag_run(j, iso):
    """P1, one level up: a whole run, not just a fill."""
    cfg = lg.make_config(isotope=iso, channel="inclusive", config=1,
                         events=400, seed=SEED)
    a = _l.Pipeline(cfg, _plan(j)).generate(0, False, 1)
    b = _l.Pipeline(cfg, _l.helicity_flip_plan(j, PZ, PE)).generate(0, False, 1)
    # BYTES, not `array_equal`: the `R` column is all-NaN on this channel and
    # NaN != NaN, which would make a value-wise comparison silently vacuous
    # in one direction and false in the other.  This is the comparison the
    # knob-provenance matrix hashes.
    for k in sorted(a):
        if isinstance(a[k], np.ndarray) and a[k].dtype != object:
            assert a[k].dtype == b[k].dtype and a[k].shape == b[k].shape, k
            assert (np.ascontiguousarray(a[k]).tobytes()
                    == np.ascontiguousarray(b[k]).tobytes()), k


# ------------------------------------------------------------------- P2

def test_typed_honours_the_typed_alignment_at_both_spins():
    """P2.  `pzz_true` is what every estimator divides by; at `typed` it is
    the number that was typed, and the populations are the closed-form
    solution rather than the ladder's."""
    t1 = _plan(1.0, 0.5, "typed")
    assert t1.pzz_true == pytest.approx(0.5, abs=1e-15)
    assert (list(t1.categories[0].populations)
            == list(_l.spin1_populations(PZ, 0.5)))
    t32 = _plan(1.5, 0.5, "typed")
    assert t32.pzz_true == pytest.approx(0.5, abs=1e-15)
    assert (list(t32.categories[0].populations)
            == list(_l.spin32_populations(PZ, 0.5, 0.0)))
    # ... and both differ from the ladder's fill at the same --pz
    for j, t in ((1.0, t1), (1.5, t32)):
        assert (list(t.categories[0].populations)
                != list(_plan(j).categories[0].populations))


@pytest.mark.parametrize("j", [1.0, 1.5])
def test_the_mode_is_a_tensor_choice_and_never_a_vector_one(j):
    """P2.  `--pz` is honoured by BOTH fills, so A_par's divisor P_e P_z is
    the same number whichever mode is named -- which is why the A_par shift
    measured in sec. A3 is a RATE effect and not a normalisation one."""
    for mode in ("ladder", "typed"):
        p = _plan(j, 0.5, mode)
        assert p.pz_true == PZ
        assert p.categories[0].vector_moment() == pytest.approx(PZ, abs=1e-12)
        assert p.pe_true == PE
        # the luminosity pattern is the plan's, not the fill's
        assert [c.name for c in p.categories] == ["apar+", "apar-"]
        assert [c.lumi_fraction for c in p.categories] == [0.5, 0.5]


# ------------------------------------------------------------------- P3

def test_outside_the_domain_it_is_refused_with_the_edge_named():
    """P3.  The registry's own command, `--plan helicity-flip --pzz 0.6` at
    `--pz 0.7`: REFUSED at J = 3/2 (outside by 0.02, the edge being T = 0.58
    where p(-1/2) = 0) and ACCEPTED at J = 1, where 0.6 is inside the wider
    spin-1 domain 0.1 <= P_zz <= 1.  Nothing is clamped: a clamped fill would
    publish an alignment nobody typed."""
    with pytest.raises(RuntimeError) as e:
        _plan(1.5, 0.6, "typed")
    msg = str(e.value)
    assert "spin32_populations" in msg
    assert "0.26 <= t <= 0.58" in msg          # both edges at this --pz
    assert "p(m = -1/2) = -0.005" in msg       # which m, and by how much
    # the edge itself is REACHED, not approached
    assert _plan(1.5, 1.0 - 0.6 * PZ, "typed").pzz_true == pytest.approx(
        0.58, abs=1e-12)
    # ... and at J = 1 the same 0.6 runs
    assert _plan(1.0, 0.6, "typed").pzz_true == pytest.approx(0.6, abs=1e-15)


def test_the_spin_1_refusal_names_its_own_edges():
    """P3.  `spin1_populations` used to throw "unphysical (pz, pzz): negative
    population", which named neither the offending m nor the domain.  It is
    reachable from the command line since `--pzz-mode typed` existed."""
    with pytest.raises(RuntimeError) as e:
        _plan(1.0, -0.5, "typed")
    msg = str(e.value)
    assert "spin1_populations" in msg
    assert "p(m = -1) = -0.1 < 0" in msg
    assert "0.1 <= pzz <= 1" in msg
    assert "0.26 <= T <= 0.58" in msg          # the smaller J = 3/2 domain


def test_an_unknown_mode_is_refused_by_name():
    with pytest.raises(ValueError) as e:
        _plan(1.0, 0.5, "maxent")
    assert "unknown pzz_mode" in str(e.value)
    assert "ladder" in str(e.value) and "typed" in str(e.value)


# ------------------------------------------------------------------- P4

#: THE PRICE, measured 2026-09-06 at the standard configuration -- inclusive,
#: beam config 1, seed 20260713, `--plan helicity-flip --pz 0.7 --pe 0.7`,
#: the shipped `--b1-model miller`.  `sigma_pb` is the cell-grid integral and
#: does not depend on the event count, so these pin at 200 events exactly what
#: sec. A3 measured at 100 000.
PRICE = {
    # (isotope, j, typed pzz): (sigma ladder, sigma typed, pzz ladder)
    ("6Li", 1.0, 0.5): (591783.2520093301, 591769.3290332475,
                        0.4094026279413209),
    ("6Li", 1.0, 0.6): (591783.2520093301, 591753.9610642146,
                        0.4094026279413209),
    ("7Li", 1.5, 0.5): (590952.426415097, 590952.4264150971,
                        0.40000000000000147),
}


@pytest.mark.parametrize("key", sorted(PRICE))
def test_the_measured_price_reproduces(key):
    """P4.  The numbers sec. A3 quotes, from committed state."""
    iso, j, pzz = key
    s_lad, s_typ, pzz_lad = PRICE[key]
    a = _pipeline(iso, j, pzz, "ladder")
    b = _pipeline(iso, j, pzz, "typed")
    assert a.sigma_pb() == pytest.approx(s_lad, rel=1e-12)
    assert b.sigma_pb() == pytest.approx(s_typ, rel=1e-12)
    assert a.plan.pzz_true == pytest.approx(pzz_lad, rel=1e-12)
    assert b.plan.pzz_true == pytest.approx(pzz, abs=1e-15)


def test_on_7Li_the_typed_fill_moves_no_physics_only_the_last_bit():
    """P4, and it is the finding the decision turns on.

    ⁷Li's rank-2 sector is identically zero (`rank2_input_report`), so the
    only fill moments the kernel sums are `sum p_m` = 1 and `<J_z>/J` = P_z --
    both fixed by `--pz`, which BOTH modes honour.  The typed fill therefore
    changes no observable the tree computes: the cross sections agree to
    within one ulp and the difference is arithmetic, not physics.  What DOES
    move is the number the file records as the fill's alignment, T = 0.4 ->
    0.5, which is the divisor of any tensor estimator.
    """
    lad, typ = _pipeline("7Li", 1.5, 0.5, "ladder"), _pipeline("7Li", 1.5, 0.5,
                                                               "typed")
    assert abs(typ.sigma_pb() / lad.sigma_pb() - 1.0) < 1e-15
    for x, y in zip(lad.sigma_per_category_pb(), typ.sigma_per_category_pb()):
        assert abs(y / x - 1.0) < 1e-15
    # the two moments the kernel actually sums, from the populations
    ms = [1.5, 0.5, -0.5, -1.5]
    for p in (lad.plan, typ.plan):
        pops = list(p.categories[0].populations)
        assert sum(pops) == pytest.approx(1.0, abs=1e-15)
        assert (sum(a * m for a, m in zip(pops, ms)) / 1.5
                == pytest.approx(PZ, abs=1e-14))
    # ... and the recorded alignment is the thing that moved
    assert lad.plan.pzz_true == pytest.approx(0.4, abs=1e-12)
    assert typ.plan.pzz_true == pytest.approx(0.5, abs=1e-15)


def test_the_tensor_estimator_divisor_moves_with_the_fill():
    """P4.  Every tensor estimator divides by the fill's own alignment, so at
    fixed N the statistical error goes as 1/P_zz: -18.1195 % on ⁶Li
    (0.409403 -> 0.5) and -20 % on ⁷Li (0.4 -> 0.5)."""
    n = 1.0e5
    for j, want in ((1.0, -18.119474), (1.5, -20.0)):
        lad = _l.spin_temperature_pzz(j, PZ)
        for f in (_l.err_azz, _l.err_cos2phi_amplitude):
            got = 100.0 * (f(n, 0.5) / f(n, lad) - 1.0)
            assert got == pytest.approx(want, abs=1e-4)


# ------------------------------------------------------------------- P5

def _rows(plan_name, iso, j, mode, pzz=0.5):
    cfg = lg.make_config(isotope=iso, channel="inclusive", config=1,
                         events=60, seed=7)
    plan = lg.make_plan(plan_name, j=j, pz=PZ, pzz=pzz, pe=PE, pzz_mode=mode)
    ctx = _l.KnobRunContext()
    ctx.plan_name, ctx.pz, ctx.pzz, ctx.pe = plan_name, PZ, pzz, PE
    ctx.pzz_mode, ctx.rel_lumi_offset = mode, 0.0
    return {r.name: r for r in _l.Pipeline(cfg, plan).knob_provenance(ctx)}


@pytest.mark.parametrize("mode", ["ladder", "typed"])
def test_the_mode_is_read_under_helicity_flip(mode):
    """P5.  Two different fills is what `read` means, and the row says so by
    MEASUREMENT -- it rebuilds the other mode's fill and counts the
    populations that differ."""
    r = _rows("helicity-flip", "6Li", 1.0, mode)["pzz_mode"]
    assert r.status == _l.KnobStatus.Read
    assert r.value == mode
    assert r.at_default == (mode == "ladder")
    assert "populations are different doubles" in r.reason


@pytest.mark.parametrize("plan_name", ["tensor-thirds", "transverse-tensor",
                                       "tensor-flip"])
def test_the_mode_is_labelled_not_read_under_the_tensor_plans(plan_name):
    """P5.  The three tensor factories have no ladder branch, so the mode
    leaves their fill exactly as it is.  LABELLED, not refused: it names a
    member of a family a plan scan sets uniformly (the criterion in
    `KnobProvenance`'s header)."""
    r = _rows(plan_name, "6Li", 1.0, "typed")["pzz_mode"]
    assert r.status == _l.KnobStatus.NotRead
    assert r.label == "not read by plan %s" % plan_name
    assert r.meta_value == r.label


def test_the_pzz_row_follows_the_mode():
    """P5, the other half: `--pzz` itself was not-read under helicity-flip and
    IS read at `--pzz-mode typed`."""
    lad = _rows("helicity-flip", "6Li", 1.0, "ladder")["pzz"]
    assert lad.status == _l.KnobStatus.NotRead
    assert lad.label == "not read by plan helicity-flip at --pzz-mode ladder"
    assert "populations ARE populations_maxent" in lad.reason
    typ = _rows("helicity-flip", "6Li", 1.0, "typed")["pzz"]
    assert typ.status == _l.KnobStatus.Read
    assert "NOT the max-entropy ladder" in typ.reason
    # unchanged where it always was read
    assert (_rows("tensor-thirds", "6Li", 1.0, "ladder")["pzz"].status
            == _l.KnobStatus.Read)


def test_a_context_naming_a_mode_that_does_not_exist_is_refused():
    """P5.  The table writes claims about what the run did, so a third mode
    string reported as the ladder would be the very defect it exists to
    prevent.  Empty stays legal and omits the row."""
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=60, seed=7)
    p = _l.Pipeline(cfg, _plan(1.0))
    ctx = _l.KnobRunContext()
    ctx.plan_name, ctx.pz, ctx.pzz, ctx.pe = "helicity-flip", PZ, 0.5, PE
    ctx.pzz_mode = "maxent"
    with pytest.raises(RuntimeError) as e:
        p.knob_provenance(ctx)
    assert "ladder" in str(e.value) and "typed" in str(e.value)
    # empty = not supplied: the row is omitted, not guessed
    ctx.pzz_mode = ""
    assert "pzz_mode" not in {r.name for r in p.knob_provenance(ctx)}
    # ... and a context with no plan name omits the whole fill block
    assert "pzz_mode" not in {r.name for r in p.knob_provenance()}


def test_the_meta_records_the_mode():
    """P5.  A knob that ran and is recorded nowhere is the defect the
    provenance table exists to prevent."""
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=60, seed=7)
    plan = _plan(1.0, 0.5, "typed")
    ctx = _l.KnobRunContext()
    ctx.plan_name, ctx.pz, ctx.pzz, ctx.pe = "helicity-flip", PZ, 0.5, PE
    ctx.pzz_mode, ctx.rel_lumi_offset = "typed", 0.0
    meta = _l.Pipeline(cfg, plan).generate(1, False, 1, context=ctx)["meta"]
    row = meta["knob_provenance"]["pzz_mode"]
    assert row["value"] == "typed"
    assert row["status"] == "read"
    assert row["flag"] == "--pzz-mode"
    assert row["at_default"] is False


# ------------------------------------------------------------------- P6

def test_the_cli_default_is_ladder_and_typed_is_reachable(capsys):
    """P6.  Both banners name the mode and the counterfactual fill, so a
    reader of either run knows which of the two fills it used."""
    assert cli.main(["--isotope", "7Li", "--plan", "helicity-flip",
                     "--events", "40", "--seed", "3"]) == 0
    out = capsys.readouterr().out
    assert "fill helicity-flip: J = 1.5, P_z = 0.7, T = 0.4" in out
    assert "--pzz-mode ladder (the default)" in out
    assert "--pzz = 0.6 is NOT read" in out

    assert cli.main(["--isotope", "7Li", "--plan", "helicity-flip",
                     "--pzz", "0.5", "--pzz-mode", "typed",
                     "--events", "40", "--seed", "3"]) == 0
    out = capsys.readouterr().out
    assert "fill helicity-flip: J = 1.5, P_z = 0.7, T = 0.5" in out
    assert "--pzz-mode typed: the fill is built at the TYPED T = 0.5" in out
    assert "ladder at this --pz would give T = 0.4" in out


def test_the_cli_refusal_is_a_systemexit_carrying_the_edge():
    """P6.  `--plan helicity-flip --pzz-mode typed` at the CLI's own defaults
    (--pz 0.7, --pzz 0.6) on ⁷Li: refused, with the edge, and not as a
    traceback."""
    with pytest.raises(SystemExit) as e:
        cli.main(["--isotope", "7Li", "--plan", "helicity-flip",
                  "--pzz-mode", "typed", "--events", "10", "--quiet"])
    msg = str(e.value)
    assert "0.26 <= t <= 0.58" in msg
    assert "outside here by 0.02" in msg
