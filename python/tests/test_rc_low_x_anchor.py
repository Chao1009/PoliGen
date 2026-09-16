# SPDX-License-Identifier: GPL-3.0-or-later
"""The RC band's LOW-x ANCHOR, `--rc-delta-low-x` -- the three priced values.

REGISTRY ROW 17 / `AUTHOR_DECISIONS.md` sec. B13, priced in
`docs/open_items/run_2026-09-06/phase_A_numbers.md` sec. A2.  The shipped band
is `w = 1 +- delta(x) tau` with `RC_DELTA_LOW_X = 0.30` hung at
`RC_X_LOW = 0.01`, and phase B established what that 0.30 rests on: it is
Gakh-Shekhovtsova's panel (a) LOWEST-x value CARRIED UPWARD IN x.  The panel
spans x = 0.00226 - 0.00966 at Q^2 = 0.1 and reads |delta| = 0.266 at its
bottom and **0.113 at x = 0.00966, the x NEAREST the anchor** -- 0.01 lies
ABOVE the panel's top.  So the shipped 0.30 errs WIDE: it is 2.65x the value
the panel actually reads next to `RC_X_LOW`.  The two alternatives are those
two panel readings, and nothing here decides between them.

The three numeric literals below (0.30, 0.266, 0.113) have ONE home each and
this file defines none of them: 0.30 is `lipolgen.RC_DELTA_LOW_X` and is read
from the module; 0.266 and 0.113 are the PAPER's two panel readings, stated in
`include/lipolgen/rc.hpp` at `RC_DELTA_LOW_X` / `RC_X_LOW` and measured off the
arXiv figure arrays in `run_2026-09-03/phase_B_numbers.md` sec. B6.3 -- they are
inputs to a price here, not shipped constants.

WHAT THIS FILE GATES, in the order the decision needs it:

  A1  the shipped default has NOT moved -- 0.30 at 0.01, 0.015 at 0.16;
  A2  the three-row table: the band half-width on A_zz at the four standard
      points, pinned to rtol 1e-12 (a REDUCED re-run of sec. A2's table --
      same anchors, same x, same Q^2, same configuration, 2000 events instead
      of the document's own construction-only evaluation);
  A3  how the anchor propagates in x -- EXACTLY proportional at x <= 0.01,
      compressed by the log-linear interpolation between the anchors, and
      IDENTICALLY ZERO at x >= 0.16, where the E12-13-011 anchor pins it;
  A4  the anchor moves exactly TWO columns and ONE meta key, and no other;
  A5  on the tagged band, where the clamp bites, the anchor does NOT move the
      clipped-event count at all -- `RcModel::clamp_tau` clips |tau| against
      `band_tau_max` and delta never enters it -- and what it DOES move is the
      width of the band on those clipped events, exactly `1 -+ delta_low`.

Run: `python -m pytest python/tests/test_rc_low_x_anchor.py -q`
"""

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import _lipolgen as _l

# The two ALTERNATIVES the registry names -- the paper's own panel (a)
# readings.  0.113 is the one at x = 0.00966, the x NEAREST `RC_X_LOW` = 0.01.
GS_PANEL_A_AT_TOP_X = 0.113        # x = 0.00966 -- the x nearest the anchor
GS_PANEL_A_AT_BOTTOM_X = 0.266     # x = 0.00226 -- what rounds to the "30 %"
ANCHORS = (0.30, GS_PANEL_A_AT_BOTTOM_X, GS_PANEL_A_AT_TOP_X)

XS = (0.01, 0.063, 0.10, 0.16)
Q2 = 5.0

# phase_A_numbers.md sec. A2, table A2.2: delta(x) * |A_zz| at P_zz = +1,
# Q^2 = 5 GeV^2, 6Li, --config 1, --plan tensor-thirds --pzz 0.6, this
# generator's own 6Li b1 (NOT phase_C sec. 8.2's x3.253983 deuteron-b1 column).
# The document prints these to seven digits (1.358472e-04 ...); pinned here
# to the full double, because rtol 1e-12 is what stops a silent drift.
HALF_WIDTH = {
    0.30:  (0.0001358472494654602, 0.00014590673286701302,
            9.680015652051017e-05, 1.920690621629717e-05),
    0.266: (0.00012045122785937472, 0.00013085664856817475,
            8.798803601270058e-05, 1.920690621629717e-05),
    0.113: (5.116913063199001e-05, 6.313126922340236e-05,
            4.83334937275573e-05, 1.920690621629717e-05),
}
# A_zz(Born) at those four x, same configuration (phase_B_numbers.md sec. B3.2
# corrected the x = 0.01 / 0.10 entries; 0.063 and 0.16 are new here).
AZZ = (-4.528241648849e-04, -1.316774298646e-03,
       -1.528923484258e-03, -1.280460414420e-03)


def _pipeline(anchor=None, channel="inclusive", events=2000, seed=99):
    kw = {} if anchor is None else {"rc_delta_low_x": anchor}
    cfg = lg.make_config(isotope="6Li", channel=channel, config=1,
                         events=events, seed=seed, rc="tensor-band", **kw)
    return lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))


# ---------------------------------------------------------------- A1

def test_the_shipped_anchor_has_not_moved():
    """A1.  Nothing this pricing did touches the default."""
    assert _l.RC_DELTA_LOW_X == 0.30
    assert _l.RC_X_LOW == 0.01
    assert _l.RC_DELTA_HIGH_X == 0.015
    assert _l.RC_X_HIGH == 0.16
    # ... and the shipped 0.30 is WIDER than both alternatives, which is the
    # only direction the reader has to be told about: `w = 1 -+ delta tau`
    # uses the MAGNITUDE alone, so the shipped choice errs conservative.
    assert GS_PANEL_A_AT_TOP_X < GS_PANEL_A_AT_BOTTOM_X < _l.RC_DELTA_LOW_X
    assert _pipeline().rc_model.delta(0.01) == _l.RC_DELTA_LOW_X


# ---------------------------------------------------------------- A2

def test_the_three_row_table_of_band_half_widths_on_azz():
    """A2.  The priced table, re-run reduced, at rtol 1e-12."""
    p0 = _pipeline()
    kern, s = p0.dis_sampler.kernel, p0.dis_sampler.s
    azz = []
    for x in XS:
        t = kern.tables(x, Q2)
        azz.append(lg.azz(t.b1, t.f1, t.f2, x, Q2 / (x * s), t.b2, 0.0))
    assert azz == pytest.approx(list(AZZ), rel=1e-12)

    for anchor in ANCHORS:
        m = (p0 if anchor == 0.30 else _pipeline(anchor)).rc_model
        got = [m.delta(x) * abs(a) for x, a in zip(XS, azz)]
        assert got == pytest.approx(list(HALF_WIDTH[anchor]), rel=1e-12), anchor
        # The band PEAKS at x = 0.063 on every anchor -- lowering the anchor
        # sharpens that peak, it does not move it to the edge.
        assert got[1] == max(got), anchor


# ---------------------------------------------------------------- A3

def test_how_the_anchor_propagates_in_x():
    """A3.  Proportional at the anchor, compressed between, zero above."""
    base = _pipeline().rc_model
    for anchor in ANCHORS[1:]:
        m = _pipeline(anchor).rc_model
        # at and below RC_X_LOW: EXACTLY the ratio of the anchors
        assert m.delta(0.01) == anchor
        assert m.delta(0.005) == anchor
        # at and above RC_X_HIGH: the E12-13-011 anchor pins it, no motion
        assert m.delta(0.16) == base.delta(0.16) == _l.RC_DELTA_HIGH_X
        assert m.delta(0.30) == base.delta(0.30) == _l.RC_DELTA_HIGH_X
        # in between: strictly compressed towards 1 -- the log-linear
        # interpolation carries only part of the anchor's move, so scaling the
        # x = 0.01 reduction across the whole low-x range OVERSTATES it.
        r_low = anchor / _l.RC_DELTA_LOW_X
        for x in (0.02, 0.063, 0.10, 0.12):
            r = m.delta(x) / base.delta(x)
            assert r_low < r < 1.0, (anchor, x, r)


# ---------------------------------------------------------------- A4

def test_the_anchor_moves_two_columns_and_two_meta_keys():
    """A4.  It is a band knob and nothing else: no kinematics, no weight.

    TWO `meta` keys, of 60: `rc_delta_low_x` and the `knob_provenance` row
    that RECORDS it -- which is what the assertion below has always said.
    (Renamed 2026-09-15; it read `..._and_one_meta_key`.)"""
    base = _pipeline().generate(0)
    for anchor in ANCHORS[1:]:
        cols = _pipeline(anchor).generate(0)
        moved = [k for k in base if k != "meta"
                 and np.asarray(base[k]).tobytes()
                 != np.asarray(cols[k]).tobytes()]
        assert moved == ["rc_tensor_lo", "rc_tensor_hi"], (anchor, moved)
        mmoved = [k for k in base["meta"]
                  if repr(base["meta"][k]) != repr(cols["meta"][k])]
        # `knob_provenance` moves because it RECORDS the value -- that is the
        # point of it; nothing else in `meta` does.
        assert sorted(mmoved) == ["knob_provenance", "rc_delta_low_x"], anchor
        assert cols["meta"]["rc_delta_low_x"] == anchor
        assert np.all(np.asarray(cols["rc_tensor_lo"])
                      + np.asarray(cols["rc_tensor_hi"]) == 2.0)


# ---------------------------------------------------------------- A5

def test_the_anchor_does_not_move_the_tagged_clipped_fraction():
    """A5.  The clip is a `tau` property; the anchor sets the CLIPPED WIDTH.

    Reduced re-run of sec. A2's 20 k tagged-alpha rows: 4000 events, seed 1,
    the same channel and plan.  The full-size numbers live in the document;
    what is gated here is the STRUCTURE, which is what the decision turns on.
    """
    runs = {a: _pipeline(None if a == 0.30 else a, channel="tagged-alpha",
                         events=4000, seed=1).generate(0) for a in ANCHORS}
    base = runs[0.30]
    n_clipped = base["meta"]["rc_clipped_band_events"]
    assert n_clipped > 0                      # the clamp bites on this channel
    for anchor, cols in runs.items():
        m = cols["meta"]
        # IDENTICAL, not merely close: `RcModel::clamp_tau` clips |tau| against
        # `band_tau_max` and `delta` never enters it.
        assert m["rc_clipped_band_events"] == n_clipped, anchor
        assert (m["rc_clipped_band_event_fraction"]
                == base["meta"]["rc_clipped_band_event_fraction"]), anchor
        assert m["rc_clipped_tail_events"] == 0        # no tail on a tag
        moved = [k for k in base if k != "meta"
                 and np.asarray(base[k]).tobytes()
                 != np.asarray(cols[k]).tobytes()]
        assert moved == ([] if anchor == 0.30
                         else ["rc_tensor_lo", "rc_tensor_hi"]), anchor
        # What the anchor DOES move: the width of the band on the clipped
        # events.  |tau| is clamped to exactly 1 there, so the edge is exactly
        # 1 -+ delta_low -- 0.700 / 0.734 / 0.887 on the three anchors.
        hi = np.asarray(cols["rc_tensor_hi"])
        lo = np.asarray(cols["rc_tensor_lo"])
        assert hi.min() == 1.0 - anchor, (anchor, hi.min())
        assert lo.max() == 1.0 + anchor, (anchor, lo.max())
        assert hi.min() > 0.0             # the whole point of `band_tau_max`
