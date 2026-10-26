#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""One figure per benchmark harness, in one style, for docs/benchmarking/REPORT.md.

Not a harness (its name does not match `t[1-6]_*.py`).  A harness that can
draw its comparison implements

    def plot(report, out_dir=None) -> str      # the written file's path

using the helpers here, and `python3 <harness>.py --plot` writes it;
`make_report.py --plots` calls `plot()` after `run()` in each child and links
the file in the row's Detail.  Figures are NOT part of `make_report.py
--check` (PNG bytes depend on the matplotlib build); every number a figure
draws is in REPORT.json, so the figure is a view, never the record.

STYLE, kept deliberately plain and identical across rows (the dataviz
reference palette, light surface, validated slots in FIXED order -- a
configuration keeps its colour on every figure):
  * the REFERENCE (data, or the published calculation) is primary-ink points
    with error bars; each GENERATOR configuration a 2 px line or >= 8 px
    marker in the slot colour below; a model band a light grey fill;
  * one panel, or a comparison panel over a pull / ratio panel (3 : 1);
  * title = `<harness>  --  STATUS` (the status word in its status colour),
    axes labelled with units, a legend whenever there are >= 2 series, no
    number printed on every point;
  * a one-line footer: the date, the git revision, 'LiPolGen benchmark';
  * 6.4 in wide, dpi 110, saved as PNG under docs/benchmarking/figures/.

Usage in a harness:

    from _plots import new_figure, save, style_series, status_title
    fig, (ax, axp) = new_figure(two_panel=True)
    ax.errorbar(x, y, yerr=dy, **style_series("reference"))
    ax.plot(x, model, **style_series("default", label="ToyF2 (shipped)"))
    status_title(fig, NAME, report["status"])
    return save(fig, NAME)                      # docs/benchmarking/figures/<NAME>.png
"""

import datetime
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
FIGURES = os.path.join(ROOT, "docs", "benchmarking", "figures")

#: the reference palette's categorical slots, in their fixed order; a
#: configuration name maps to ONE slot everywhere (never cycled)
PALETTE = {
    "default": "#2a78d6",    # slot 1 blue   -- the shipped default
    "mstw": "#eb6834",       # slot 2 orange -- MSTW2008 LO
    "ct18nlo": "#1baf7a",    # slot 3 aqua   -- CT18NLO
    "alt": "#eda100",        # slot 4 yellow -- the alternative wave / model (AV18, instant form, ...)
    "r1998": "#e87ba4",      # slot 5 magenta
    "nnpdfpol": "#6250d6",   # slot 7 violet
    "alt2": "#008300",       # slot 6 green  -- a third configuration (CD-Bonn, polrad-full, ...)
    "alt3": "#e34948",       # slot 8 red    -- a fourth, last resort
}
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BAND = "#e1e0d9"
SURFACE = "#fcfcfb"
STATUS_COLOUR = {"pass": "#0ca30c", "fail": "#d03b3b", "blocked": MUTED,
                 "blocked (environment)": MUTED, "error": "#d03b3b"}


def mpl():
    """matplotlib with the non-interactive backend; imported lazily so a
    harness without `--plot` never needs it."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK_SECONDARY,
        "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_SECONDARY,
        "ytick.labelcolor": INK_SECONDARY, "text.color": INK,
        "legend.frameon": False, "legend.fontsize": 8, "font.size": 9,
        "lines.linewidth": 2.0, "lines.markersize": 5,
    })
    return plt


def new_figure(two_panel=False, height=4.0):
    """A 6.4-in-wide figure; with `two_panel`, a comparison axis over a
    pull / ratio axis sharing x (height ratios 3 : 1)."""
    plt = mpl()
    if two_panel:
        fig, axes = plt.subplots(2, 1, figsize=(6.4, height + 1.2), sharex=True,
                                 layout="constrained",
                                 gridspec_kw={"height_ratios": [3, 1]})
        axes[1].axhline(0.0, color=AXIS, linewidth=1.0, zorder=0)
    else:
        fig, axes = plt.subplots(figsize=(6.4, height), layout="constrained")
    # room for the title row (status_title) above and the footer (save) below
    fig.get_layout_engine().set(rect=(0.0, 0.03, 1.0, 0.905), hspace=0.02)
    return fig, (tuple(axes) if two_panel else axes)


def style_series(kind, label=None, marker=None):
    """Keyword arguments for `ax.plot` / `ax.errorbar`: the reference as
    primary-ink points with error bars, a configuration in its slot colour."""
    if kind == "reference":
        kw = dict(color=INK, fmt=marker or "o", markersize=4, capsize=2,
                  linewidth=1.0, elinewidth=1.0, zorder=5, label=label or "reference")
        return kw
    colour = PALETTE.get(kind, PALETTE["alt3"])
    kw = dict(color=colour, zorder=3, label=label or kind)
    if marker:
        kw.update(marker=marker, linestyle="none", markersize=6)
    return kw


def band(ax, x, lo, hi, label=None):
    """A model band as a light grey fill (never a saturated colour)."""
    ax.fill_between(x, lo, hi, color=BAND, alpha=0.8, linewidth=0, zorder=1, label=label)


def status_title(fig, name, status, extra=""):
    """`<name>  --  STATUS` with the status word in its colour; `extra` is a
    short second line (the window, the tolerance)."""
    s = str(status).lower()
    colour = STATUS_COLOUR.get(s, INK)
    fig.suptitle(name, x=0.01, y=0.995, ha="left", fontsize=10, color=INK)
    fig.text(0.99, 0.995, s.upper(), ha="right", va="top", fontsize=10,
             color=colour, fontweight="bold")
    if extra:
        fig.text(0.01, 0.955, extra, ha="left", va="top", fontsize=8, color=INK_SECONDARY)


def revision():
    try:
        out = subprocess.run(["git", "--no-optional-locks", "rev-parse", "--short", "HEAD"],
                             cwd=ROOT, capture_output=True, text=True, timeout=10)
        rev = out.stdout.strip() or "unknown"
        dirty = subprocess.run(["git", "--no-optional-locks", "status", "--porcelain"],
                               cwd=ROOT, capture_output=True, text=True, timeout=10).stdout.strip()
        return rev + ("+" if dirty else "")
    except Exception:
        return "unknown"


def save(fig, name, out_dir=None, date=None):
    """Footer, then PNG at dpi 110 under `out_dir` (default FIGURES);
    returns the path.  `date` defaults to today (UTC)."""
    out_dir = out_dir or FIGURES
    os.makedirs(out_dir, exist_ok=True)
    date = date or datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d")
    fig.text(0.01, 0.005, "LiPolGen benchmark  %s  rev %s" % (date, revision()),
             ha="left", va="bottom", fontsize=7, color=MUTED)
    path = os.path.join(out_dir, name + ".png")
    fig.savefig(path, dpi=110)
    mpl().close(fig)
    return path


def figure_rel(path):
    """The path as REPORT.md links it (relative to docs/benchmarking/)."""
    return os.path.relpath(path, os.path.dirname(FIGURES))
