#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 7 -- the change -> tier map, as data -- and which
harnesses a diff must re-run (milestone M5: "the tiers the change touches").

Not a harness (its name does not match `t[1-6]_*.py`).  `SECTION7` below is
the plan's sec. 7 list transcribed entry by entry; `python/tests/
test_bench_report.py` re-reads sec. 7 from the plan and fails if the two
drift apart, so the plan stays the spec and this file its machine form.

A tier maps to harness files by the file-name prefix alone: `T3`, `T3-EMC`
and `T3-cluster` all select every `t3_*.py` (the sub-tier is printed, but no
harness carries it in its name, so none is guessed).

Two rules are this file's own, NOT sec. 7's, and are labelled as such in the
output: a changed harness re-runs itself, and a changed vendored file under
`validation/benchmarks/data/` re-runs every harness whose source names it.
Two sec. 7 items are not files and cannot be read off a diff: "any PDF-set
version" (an LHAPDF set upgrade in the environment) is printed as a standing
note; "the Pomeron tier" is the gamma*-Pomeron code in `pythia_bridge.*`
(docs/PYTHIA_BRIDGE.md sec. 12), which that entry's own glob already covers.
A changed code file that matches no sec. 7 pattern (e.g. `src/core/sf.cpp`:
sec. 7 names `sf.hpp` only) is listed as UNMAPPED -- decide by hand; nothing
is guessed.

Run:  python3 validation/benchmarks/tiers.py --since <rev>   # vs the work tree
      python3 validation/benchmarks/tiers.py --files src/core/rc.cpp ...
      python3 validation/benchmarks/tiers.py --since <rev> --names  # for
          make_report.py --only
      python3 validation/benchmarks/tiers.py --map                  # sec. 7

`--since` runs `git --no-optional-locks diff --name-only <rev>` (the work
tree against <rev>) plus the untracked files (`git ls-files --others
--exclude-standard`); it never writes to the repository.  Exit status: 0
when it printed its answer; 2 when broken (git failed, a bad argument).
"""

import argparse
import fnmatch
import glob
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
BENCH_REL = "validation/benchmarks"
DATA_REL = BENCH_REL + "/data"

#: BENCHMARK_PLAN.md sec. 7, one entry per bullet:
#: (file patterns, non-file items, tiers, note).  A pattern without '/'
#: matches a file's basename; one with '/' matches its repository path.
SECTION7 = (
    (("sf.hpp", "lhapdf_sf.cpp", "mstw_sf.cpp"), ("any PDF-set version",),
     ("T1", "T2", "T3-EMC"), None),
    (("b1_nuclear.*", "cluster.*", "cluster_config.*", "tagged.*",
      "data/vmc/**"), (), ("T2", "T3"), None),
    (("breakup.*", "spectator.*", "fsi.*", "triton_sf.*"), (),
     ("T3-cluster", "T4-tagged"), None),
    (("coherent.*", "pythia_bridge.*"), ("the Pomeron tier",),
     ("T4-inclusive", "T4-VM", "T6"),
     "the Pomeron tier is the gamma*-Pomeron code in pythia_bridge.* "
     "(docs/PYTHIA_BRIDGE.md sec. 12), covered by this entry's glob"),
    (("rc.*",), (), ("T4-RC", "T2"),
     "HERMES b1 was RC-corrected -- the comparison states what RC it assumes"),
    (("spin.*", "bookkeeping.*", "sampler.*", "rng.*"), (), ("T5",),
     "and T1 #1: plan row 1 (t5_epios_source_modes, tier T5/T1), already a "
     "t5_ harness"),
    (("hepmc_writer.*", "cli.py", "export.py"), (), ("T6",), None),
)

#: where sec. 7's patterns apply: the code, not the docs or the tests
CODE_PREFIXES = ("include/", "src/", "python/lipolgen/", "data/",
                 "third_party/")
CODE_FILES = ("python/bindings.cpp", "CMakeLists.txt", "pyproject.toml",
              "env.sh")
#: never re-run anything by themselves (docs, tests, the report itself)
QUIET_PREFIXES = ("docs/", "tests/", "python/tests/", "examples/")

STANDING_NOTE = ("not readable off a diff (sec. 7): 'any PDF-set version' -- "
                 "an LHAPDF set upgraded in the environment re-runs T1, T2, "
                 "T3-EMC by hand")


def tier_digit(tier):
    m = re.match(r"T(\d)", tier)
    if not m:
        raise ValueError("not a tier: %r" % tier)
    return int(m.group(1))


def harness_files(bench_dir=None):
    """Every harness stem, sorted: `t[1-6]_*.py` in the benchmarks dir."""
    bench_dir = bench_dir or HERE
    return sorted(os.path.splitext(os.path.basename(p))[0]
                  for p in glob.glob(os.path.join(bench_dir, "t[1-6]_*.py")))


def harnesses_for_tiers(tiers, bench_dir=None):
    digits = {tier_digit(t) for t in tiers}
    return [h for h in harness_files(bench_dir) if int(h[1]) in digits]


def _is_code(path):
    return path.startswith(CODE_PREFIXES) or path in CODE_FILES


def match_section7(path):
    """[(pattern, tiers, note)] of every sec. 7 entry `path` matches."""
    base = path.rsplit("/", 1)[-1]
    hits = []
    for patterns, _items, tiers, note in SECTION7:
        for pat in patterns:
            ok = (fnmatch.fnmatchcase(path, pat.replace("**", "*"))
                  if "/" in pat else fnmatch.fnmatchcase(base, pat))
            if ok:
                hits.append((pat, tiers, note))
                break
    return hits


def _harnesses_naming(basename, bench_dir):
    out = []
    for h in harness_files(bench_dir):
        try:
            with open(os.path.join(bench_dir, h + ".py"), encoding="utf-8") as f:
                if basename in f.read():
                    out.append(h)
        except OSError:                              # pragma: no cover
            continue
    return out


def plan_for(files, bench_dir=None):
    """What a set of changed repository paths must re-run.

    Returns a dict: `matches` [(file, pattern, tiers)] by sec. 7;
    `unmapped` code files no sec. 7 pattern matches; `harness_changes`
    (re-run themselves), `data_changes` {data file: [harnesses]} and
    `helper_changes` {_report.py: [harnesses naming it]} (this file's own
    rules); `tiers`; `harnesses` (sorted stems); `notes`."""
    bench_dir = bench_dir or HERE
    known = set(harness_files(bench_dir))
    matches, unmapped, harness_changes, data_changes = [], [], [], {}
    helper_changes = {}
    tiers, harnesses, notes = set(), set(), []
    for f in sorted(set(files)):
        f = f.replace(os.sep, "/")
        if f.startswith(DATA_REL + "/"):
            names = _harnesses_naming(f.rsplit("/", 1)[-1], bench_dir)
            data_changes[f] = names
            harnesses.update(names)
            continue
        if f.startswith(BENCH_REL + "/"):
            stem = os.path.splitext(f.rsplit("/", 1)[-1])[0]
            if stem in known:
                harness_changes.append(stem)
                harnesses.add(stem)
            elif stem == "_report":
                names = _harnesses_naming("_report", bench_dir)
                helper_changes[f] = names
                harnesses.update(names)
                # make_report.py imports _report, whose normalise_row /
                # jsonable / extract_subrows shape every REPORT.json row: a
                # change there can stale the report even when no harness
                # names it (so far none does)
                note = ("%s changed: no measurement moves, but regenerate "
                        "the report (make_report.py --all) and read the "
                        "--check diff" % f)
                if note not in notes:
                    notes.append(note)
            elif stem in ("make_report", "tiers"):
                note = ("%s changed: no measurement moves, but regenerate "
                        "the report (make_report.py --all) and read the "
                        "--check diff" % f)
                if note not in notes:
                    notes.append(note)
            continue
        if f.startswith(QUIET_PREFIXES) or not _is_code(f):
            continue
        hits = match_section7(f)
        if not hits:
            unmapped.append(f)
        for pat, t, note in hits:
            matches.append((f, pat, t))
            tiers.update(t)
            if note and note not in notes:
                notes.append(note)
    for t in tiers:
        harnesses.update(harnesses_for_tiers([t], bench_dir))
    notes.append(STANDING_NOTE)
    return dict(matches=matches, unmapped=unmapped,
                harness_changes=sorted(harness_changes),
                data_changes=data_changes, helper_changes=helper_changes,
                tiers=sorted(tiers, key=lambda t: (tier_digit(t), t.lower())),
                harnesses=sorted(h for h in harnesses if h in known),
                notes=notes)


def _git(args, root):
    out = subprocess.run(["git", "--no-optional-locks"] + args, cwd=root,
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("git %s failed: %s" % (" ".join(args),
                                                   out.stderr.strip()))
    return [l for l in out.stdout.splitlines() if l.strip()]


def changed_since(rev, root=None, untracked=True):
    """Repository paths changed in the work tree against `rev` (plus the
    untracked, non-ignored files unless `untracked=False`).  Read-only."""
    root = root or ROOT
    files = _git(["diff", "--name-only", rev, "--"], root)
    if untracked:
        files += _git(["ls-files", "--others", "--exclude-standard"], root)
    return sorted(set(files))


def print_plan(plan, files, out=None):
    out = out or sys.stdout
    w = lambda s="": print(s, file=out)             # noqa: E731
    w("changed files: %d" % len(files))
    for f, pat, t in plan["matches"]:
        w("  sec.7  %-44s %-16s -> %s" % (f, pat, ", ".join(t)))
    for f in plan["unmapped"]:
        w("  UNMAPPED %s  (code, but no sec. 7 pattern: decide by hand)" % f)
    for h in plan["harness_changes"]:
        w("  harness %s.py changed -> re-runs itself (this file's rule)" % h)
    for f, hs in sorted(plan["data_changes"].items()):
        w("  vendored %s -> %s (this file's rule)" % (
            f, ", ".join(hs) if hs else "no harness names it"))
    for f, hs in sorted(plan["helper_changes"].items()):
        w("  helper %s -> %s (this file's rule)" % (
            f, ", ".join(hs) if hs else "no harness imports it"))
    w("tiers: %s" % (", ".join(plan["tiers"]) or "none"))
    w("harnesses to re-run: %s" % (", ".join(plan["harnesses"]) or "none"))
    for n in plan["notes"]:
        w("note: %s" % n)


def print_map(out=None):
    out = out or sys.stdout
    for patterns, items, tiers, note in SECTION7:
        print("%s -> %s%s" % (", ".join(list(patterns) + list(items)),
                              ", ".join(tiers),
                              "  (%s)" % note if note else ""), file=out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--since", metavar="REV",
                   help="changed files: the work tree against REV")
    g.add_argument("--files", nargs="+", metavar="PATH",
                   help="changed files, as repository paths")
    g.add_argument("--map", action="store_true", help="print sec. 7 and exit")
    ap.add_argument("--names", action="store_true",
                    help="print only the harness names, comma-separated")
    ap.add_argument("--no-untracked", action="store_true",
                    help="with --since: ignore untracked files")
    try:
        a = ap.parse_args(argv)
    except SystemExit as exc:
        return 0 if exc.code == 0 else 2
    try:
        if a.map:
            print_map()
            return 0
        files = (changed_since(a.since, untracked=not a.no_untracked)
                 if a.since else list(a.files))
        plan = plan_for(files)
    except Exception as exc:                         # broken: exit 2
        print("tiers.py: %s" % exc, file=sys.stderr)
        return 2
    if a.names:
        print(",".join(plan["harnesses"]))
    else:
        print_plan(plan, files)
    return 0


if __name__ == "__main__":
    sys.exit(main())
