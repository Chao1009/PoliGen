#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Freeze `default_inclusive_kernel(LI6())`'s tensor slots as a 1e-12 gate.

Design D section 8, Agent B "Step 0".  Nothing else in the tree stores what
`pipeline.cpp::default_inclusive_kernel` returns:

  * `validation/reference/*.json` come from `dump_polligen_reference.py`,
    which dumps the PYTHON generator and builds its own kernel options with
    `toy_b1(mode="toy")` -- a different b1 from the default kernel's
    `Li6B1(MillerB1)`;
  * `tests/test_reference.cpp::build_kernel` likewise constructs its own
    `InclusiveKernel::Options` and never calls `default_inclusive_kernel`.

So the existing rtol-1e-12 gates do NOT cover that function at all, and a
regression in it -- for instance while adding the `B1Model` overload of
design D section 4.2 -- would be invisible.  This script dumps
`InclusiveKernel::tables(x, q2)`'s f1 / b1 / b2 / delta on a 200-point x grid
at three Q2 into `validation/reference/b1_default_li6.json`, which
`tests/test_b1_nuclear.cpp` (T9) reads back at rtol 1e-12.

It runs through the INSTALLED pybind11 module (`import lipolgen`), i.e. the
same C++ `default_inclusive_kernel` the test calls, so the file records what
the library does TODAY.  Regenerate it only when a default is deliberately
changed, and say so in the commit message: that is the whole point of the pin.

Floats are written with json's repr() formatter, i.e. the shortest decimal
that round-trips a double, so the C++ side (`jsonmin`, strtod) reads back the
identical bits.

THE PROVENANCE LABEL (review 2026-09-26).  The file's `provenance` string
describes the numbers it sits next to, so this script never writes a label
that could describe other numbers:

  * if every numeric block (x, q2, tables, kernel) that build() returns is
    byte-identical to the file's, the file's own `provenance` is KEPT
    verbatim -- it described exactly these numbers and still does;
  * if any of them moved, nothing is written unless
    `--accept-changed-numbers` is given (the deliberate-change case above),
    and then the label is DERIVED at that moment: the fixed SELF-PIN text
    (LABEL_FIXED) plus the dump date (UTC), `git describe --always --dirty`
    of the tree and the blocks that moved.  The same derived label is used
    when there is no earlier file.

Until 2026-09-26 the whole label was a constant here, so the first
legitimate regeneration would have re-stamped "last committed 2026-09-03,
b1071b1 ... byte-identical to that dump" onto numbers that were neither.

Writes are atomic (build first, a temporary file in the same directory,
then os.replace), for the pin and for `_manifest.json`: until 2026-09-26
the pin was opened for writing before build() ran, so a failing build()
left it at 0 bytes.  A file whose new text equals the old byte for byte is
not rewritten at all.

Run:
    source env.sh
    python3 validation/dump_b1_default_li6.py --check   # compare, write nothing
    python3 validation/dump_b1_default_li6.py           # regenerate
    python3 validation/dump_b1_default_li6.py --accept-changed-numbers

Exit status: 0 when the file matches (`--check`) or was written / left
unchanged; 1 when `--check` finds a difference, or a regeneration is refused
(numbers moved without `--accept-changed-numbers`, or the existing file
cannot be read); 2 on a usage error.
"""

import argparse
import datetime
import json
import math
import os
import pathlib
import subprocess
import sys
import tempfile

import lipolgen as lg

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "reference"
NAME = "b1_default_li6.json"

# 200 log-spaced x from 1e-3 to 0.95 -- wide enough to cross MillerB1's own
# structure (its sign change sits near x = 0.4) and to reach the large-x tail
# where the digitized tables end.  The grid is WRITTEN INTO the file and read
# back by the test, so it never has to be reproduced by a C++ linspace.
N_X = 200
X_LO, X_HI = 1.0e-3, 0.95
Q2_GRID = [1.0, 2.5, 10.0]

# The blocks the label vouches for; the label is kept only while all four
# are byte-identical to the file's.
NUMERIC_BLOCKS = ("x", "q2", "tables", "kernel")

# The self-pin label (BENCHMARK_PLAN.md sec. 5.5), the part that holds for
# ANY dump of this script.  The file committed on 2026-09-23 carries a longer
# hand-written label (dated to the b1071b1 dump of 2026-09-03), which is kept
# as long as its numbers are; a derived label is this text plus
# dated_label_suffix().
LABEL_FIXED = (
    'SELF-PIN, NOT A REFERENCE OF ANY KIND.  Every number in this file '
    "was dumped from LiPolGen's OWN C++ "
    '(default_inclusive_kernel(LI6()) through the pybind11 module) by '
    'validation/dump_b1_default_li6.py; no external source, no polligen '
    'output and no measurement enters it.  It is a REGRESSION GUARD: T9 '
    '(tests/test_b1_nuclear.cpp) and python/tests/test_b1_model.py '
    'check that the library still returns what it returned when these '
    "numbers were dumped, at rtol 1e-12.  Agreement with it proves "
    "'unchanged', never 'correct'.  Unlike every other file in this "
    'directory it does not come from dump_polligen_reference.py.'
)


def x_grid():
    lo, hi = math.log(X_LO), math.log(X_HI)
    return [math.exp(lo + (hi - lo) * i / (N_X - 1)) for i in range(N_X)]


def build():
    ion = lg.li6()
    kernel = lg.default_inclusive_kernel(ion)
    xs = x_grid()
    blocks = []
    for q2 in Q2_GRID:
        f1, b1, b2, delta = [], [], [], []
        for x in xs:
            t = kernel.tables(x, q2)
            f1.append(t.f1)
            b1.append(t.b1)
            b2.append(t.b2)
            delta.append(t.delta)
        blocks.append({"q2": q2, "f1": f1, "b1": b1, "b2": b2, "delta": delta})
    return {
        "generator": "LiPolGen/validation/dump_b1_default_li6.py",
        "description": (
            "InclusiveKernel::tables(x, q2) of default_inclusive_kernel(LI6()) "
            "with the LIBRARY DEFAULT b1 backend, B1Model::Miller = "
            "Li6B1(MillerB1) through LI6_B1_RANK2_TRANSFER, and delta_func = "
            "toy_delta_gluon(scale = 1e-2).  Design D section 8 step 0; read "
            "by tests/test_b1_nuclear.cpp T9 at rtol 1e-12.  This is the pin "
            "that makes '--b1-model miller is bit-for-bit unchanged' a "
            "testable statement."
        ),
        "kernel": {
            "function": "default_inclusive_kernel(const Ion&)",
            "ion": "6Li",
            "spin": ion.spin,
            "b1_model": "miller",
            "b1_func": "Li6B1(MillerB1)::b1_func()",
            "delta_func": "toy_delta_gluon(x, q2, f1, 1e-2)",
            "f2_source": "default (ToyF2)",
            "target_mass": kernel.target_mass,
            "tensor_gamma": kernel.tensor_gamma,
        },
        "x": x_grid(),
        "q2": Q2_GRID,
        "tables": blocks,
    }


def render(doc):
    """The file's exact text (what json.dump with these arguments plus a
    final newline always wrote)."""
    return json.dumps(doc, indent=1, sort_keys=True) + "\n"


def _canon(block):
    return json.dumps(block, indent=1, sort_keys=True)


def changed_blocks(new, old):
    """The NUMERIC_BLOCKS whose serialized text differs (or is missing)."""
    return [k for k in NUMERIC_BLOCKS
            if k not in old or k not in new or _canon(new[k]) != _canon(old[k])]


def worst_table_rel(new, old):
    """Worst |new - old| / |old| over the f1/b1/b2/delta tables (a message
    aid for a refused regeneration; inf when the shapes differ)."""
    try:
        worst = 0.0
        if len(new["tables"]) != len(old["tables"]):
            return math.inf
        for bn, bo in zip(new["tables"], old["tables"]):
            for key in ("f1", "b1", "b2", "delta"):
                if len(bn[key]) != len(bo[key]):
                    return math.inf
                for a, b in zip(bn[key], bo[key]):
                    d = abs(a - b)
                    r = d / abs(b) if b else d
                    if not math.isfinite(r):
                        return math.inf
                    worst = max(worst, r)
        return worst
    except (KeyError, TypeError, ValueError):
        return math.inf


def git_describe():
    """`git describe --always --dirty` of the tree this script sits in, or
    "" when git cannot say (no git, not a checkout)."""
    try:
        res = subprocess.run(
            ["git", "-C", str(HERE), "describe", "--always", "--dirty",
             "--abbrev=7"],
            capture_output=True, text=True, timeout=30, check=True)
        return res.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def _tree(rev):
    if rev:
        return ("the LiPolGen tree at %s (git describe --always --dirty; "
                "-dirty = uncommitted changes somewhere in the tree)" % rev)
    return "a LiPolGen tree whose revision git describe could not give"


def dated_label_suffix(today, rev, moved, unreadable=None):
    """The part of a derived label that only holds for THIS dump."""
    text = "  Numbers dumped %s (UTC) from %s." % (today, _tree(rev))
    if unreadable is not None:
        text += ("  The file they replaced could not be read (%s); written "
                 "with --accept-changed-numbers." % unreadable)
    elif moved is None:
        text += "  No earlier pin was on disk."
    else:
        text += ("  They replaced an earlier pin whose %s block(s) differed; "
                 "written with --accept-changed-numbers, i.e. as a deliberate "
                 "change of a default." % ", ".join(moved))
    return text


def load_existing(path):
    """(text, doc) of the pin on disk, (None, None) when there is none.
    Raises ValueError when it exists but is not a JSON object."""
    if not path.is_file():
        return None, None
    text = path.read_text()
    doc = json.loads(text)
    if not isinstance(doc, dict):
        raise ValueError("%s is not a JSON object" % path)
    return text, doc


def assemble(data, old, accept_changed, today=None, rev=None,
             unreadable=None):
    """build()'s `data` plus its `provenance`, by the rule in the docstring.

    Returns (doc, moved): `moved` is None when there was no (readable)
    earlier pin -- `unreadable` then names why, if a file was there -- else
    the list of NUMERIC_BLOCKS that differ (empty when none).  Raises
    SystemExit when numbers moved and `accept_changed` is false."""
    doc = dict(data)
    if old is None:
        moved = None
    else:
        moved = changed_blocks(data, old)
        if not moved and isinstance(old.get("provenance"), str):
            doc["provenance"] = old["provenance"]
            return doc, moved
        if moved and not accept_changed:
            raise SystemExit(
                "REFUSED: the library's default_inclusive_kernel(LI6()) no "
                "longer reproduces %s -- block(s) %s differ (worst relative "
                "difference in the tables %.3e).  Nothing was written.  If a "
                "default was changed DELIBERATELY, rerun with "
                "--accept-changed-numbers (the label is then derived afresh) "
                "and say so in the commit message; otherwise this is the "
                "regression the pin exists to catch."
                % (NAME, ", ".join(moved), worst_table_rel(data, old)))
    if today is None:
        today = datetime.datetime.now(datetime.timezone.utc).date().isoformat()
    if rev is None:
        rev = git_describe()
    if moved == []:
        # Same numbers, but the earlier file had no label to keep.
        suffix = ("  Label written %s (UTC) from %s; the numbers are "
                  "byte-identical to the unlabelled pin it replaced."
                  % (today, _tree(rev)))
    else:
        suffix = dated_label_suffix(today, rev, moved, unreadable)
    doc["provenance"] = LABEL_FIXED + suffix
    return doc, moved


def _new_file_mode(path):
    try:
        return path.stat().st_mode & 0o7777
    except FileNotFoundError:
        umask = os.umask(0)
        os.umask(umask)
        return 0o666 & ~umask


def write_atomically(path, text):
    """Replace `path` by `text` in one os.replace; False (and no write) when
    the file already holds exactly `text`."""
    try:
        if path.read_text() == text:
            return False
    except FileNotFoundError:
        pass
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix="." + path.name + ".",
                               suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.chmod(tmp, _new_file_mode(path))
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise
    return True


def update_manifest(path):
    """Add this file to `_manifest.json`'s `files` list, idempotently.

    `dump_polligen_reference.py` writes the same manifest, and since the
    review of 2026-09-03 it MERGES rather than rewriting it wholesale: entries
    it did not produce (this one) and the whole `generators` map survive a
    polligen re-dump, in either order.  The `generators` map says which script
    owns which file, because the manifest's top-level `generator` key names
    only the polligen dump.  Written atomically, and not at all when nothing
    changes.
    """
    if not path.is_file():
        return False
    doc = json.loads(path.read_text())
    files = doc.setdefault("files", [])
    if NAME not in files:
        files.append(NAME)
    gens = doc.setdefault("generators", {})
    gens[NAME] = "LiPolGen/validation/dump_b1_default_li6.py"
    for f in files:
        gens.setdefault(f, doc.get("generator", ""))
    return write_atomically(path, json.dumps(doc, indent=1, sort_keys=True)
                            + "\n")


def check(out):
    """--check: compare, write nothing.  0 when the file is exactly what a
    regeneration would leave, 1 otherwise."""
    data = build()
    try:
        text, old = load_existing(out)
    except ValueError as err:
        print("%s: cannot read it (%s)" % (out, err))
        return 1
    if old is None:
        print("%s: missing" % out)
        return 1
    moved = changed_blocks(data, old)
    if moved:
        print("%s: numeric block(s) %s differ from build() (worst relative "
              "difference in the tables %.3e); a regeneration would be "
              "REFUSED without --accept-changed-numbers"
              % (out, ", ".join(moved), worst_table_rel(data, old)))
        return 1
    doc, _ = assemble(data, old, accept_changed=False)
    if render(doc) != text:
        print("%s: every numeric block (%s) is byte-identical to build(), "
              "but the file's other text differs from what a regeneration "
              "writes" % (out, ", ".join(NUMERIC_BLOCKS)))
        return 1
    print("%s: byte-identical to what a regeneration writes (numeric blocks "
          "%s, %d x %d table points; provenance label kept); nothing written"
          % (out, ", ".join(NUMERIC_BLOCKS), len(data["q2"]), len(data["x"])))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--check", action="store_true",
                    help="compare build() with the file; write nothing")
    ap.add_argument("--accept-changed-numbers", action="store_true",
                    help="write even though numeric blocks moved (a "
                         "deliberate default change); the label is derived "
                         "afresh")
    args = ap.parse_args(argv)
    out = OUT / NAME
    if args.check:
        return check(out)
    data = build()                      # before anything is opened
    unreadable = None
    try:
        _, old = load_existing(out)
    except ValueError as err:
        if not args.accept_changed_numbers:
            raise SystemExit(
                "REFUSED: cannot read the existing %s (%s); restore it with "
                "`git checkout -- validation/reference/%s`, or pass "
                "--accept-changed-numbers to write a fresh pin"
                % (out, err, NAME))
        old, unreadable = None, "%s: %s" % (type(err).__name__, err)
    doc, moved = assemble(data, old, args.accept_changed_numbers,
                          unreadable=unreadable)
    wrote = write_atomically(out, render(doc))
    manifest = update_manifest(OUT / "_manifest.json")
    if not wrote:
        print("%s: unchanged (numbers and label byte-identical); nothing "
              "written" % out)
    elif moved is None:
        print("wrote %s (no readable earlier pin; label derived)" % out)
    elif moved:
        print("wrote %s: block(s) %s moved (--accept-changed-numbers); label "
              "derived afresh" % (out, ", ".join(moved)))
    else:
        print("wrote %s (numbers unchanged)" % out)
    if manifest:
        print("updated %s" % (OUT / "_manifest.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
