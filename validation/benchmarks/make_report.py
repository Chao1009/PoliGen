#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regenerate docs/benchmarking/REPORT.md and REPORT.json from the harnesses.

BENCHMARK_PLAN.md sec. 3: each harness "writes one row of
docs/benchmarking/REPORT.md -- value, reference, tolerance, pass/fail, and the
configuration that made it ... The report is regenerated, never
hand-edited."  This is the generator.  It is not a harness (its name does
not match `t[1-6]_*.py`) and it changes no number: it runs every harness's
own `run()` and writes down what came back.

WHAT IT DOES.  It discovers `validation/benchmarks/t[1-6]_*.py` (sorted),
imports each with importlib (as the pytests' `_load` does) and times
`run(verbose=False)` with the harness's output captured, each harness in its
OWN child process by default, so that a crash -- an exception, a
`sys.exit`, even a segfault in the C++ library -- costs that harness's row
and nothing else:
  * an exception whose message says an LHAPDF set is missing ("Info file not
    found for PDF set ...") -> status `blocked (environment)`, reason
    "environment: LHAPDF set <name> not installed";
  * any other exception, or a child that dies or times out -> status
    `error`, reason = the exception's type and FIRST line (the traceback
    goes to this script's stderr, not into the report);
  * an EXPENSIVE harness -- a row that costs more than 30 s, or one the
    plan counts as expensive by kind (a generator-vs-generator sample, as
    t4_pythia_ep_closure, ~10 s); its module sets `EXPENSIVE = True`, or it
    is in this file's `EXPENSIVE` set -- runs only with
    `LIPOLGEN_BENCH=1` or `--all`; otherwise its row is `skipped (opt-in)`
    (validation/benchmarks/README.md, "The opt-in variable").
Each row is normalised (`_report.normalise_row`) to one schema: name, tier
(from the t<N>_ prefix), plan_row (the BENCHMARK_PLAN.md sec. 4 row number
when the harness is one of the ten, read from the plan's table, else null),
status, generator, reference, tolerance, reason, config, data_sha256 (sha256
of every byte of each vendored file the module names in an upper-case
constant; null for a named file that is not on disk), runtime_s (wall time of
`run()` alone, time.perf_counter), expensive, subrows, notes.  REPORT.json is that, plus a
`generated` header (date, git revision and dirty flag, the environment:
Python / numpy / PYTHIA / LHAPDF / HepMC3 versions and the LHAPDF sets
present); REPORT.md renders the same data -- a GENERATED header, one table
row per harness, the sub-row detail, and a "planned, not wired" section
built from the blocked rows of this run and, verbatim, BENCHMARK_PLAN.md
sec. 6 (the U-rows) and the sec. 4 rows that have no harness.

Run (the committed report is made with --all, i.e. every row run):
    cd LiPolGen && source env.sh
    python3 validation/benchmarks/make_report.py --all        # regenerate
    python3 validation/benchmarks/make_report.py --check      # still true?
    python3 validation/benchmarks/make_report.py --all \\
        --compare docs/benchmarking/REPORT.json               # M2 rule
Options:
    --out-dir DIR   write REPORT.md / REPORT.json there (default
                    docs/benchmarking; a partial report -- --only / --since --
                    is never written there, it needs an explicit --out-dir)
    --only A,B      run only these harnesses (stems, comma-separated)
    --since REV     run only what tiers.py says a diff against REV re-runs
                    (BENCHMARK_PLAN.md sec. 7)
    --all           also run the opt-in (expensive) rows; = LIPOLGEN_BENCH=1
    --compare OLD   the M2 rule: print every row (and sub-row) whose status
                    changed or whose headline numbers moved against OLD
                    (a REPORT.json); exit 1 if any.  Writes nothing unless
                    --out-dir is given.  The headline numbers of a row are the
                    numbers printed in its generator field, at the precision
                    the harness prints; of a sub-row, its numeric values.  A
                    row that ran in OLD but not here (skipped, error) is one
                    change, printed with its reason.
    --check         regenerate into a temporary directory -- in the mode the
                    committed REPORT.json records (its `opt_in`) -- and diff
                    against the committed REPORT.md / REPORT.json, ignoring
                    the date, the git revision, the runtimes and the
                    environment header; exit 1 on any other difference
                    (a hand edit, a moved number, a new or vanished row).
    --in-process    run the harnesses in this process (no isolation)
    --timeout S     per-harness child timeout (default 1800 s)

Exit status: 0 when the report was generated (whatever the verdicts it
carries, like a harness: a FAIL row or an error row is a recorded fact);
1 when --compare or --check found a difference; 2 when this script itself is
broken or misused (bad option, unreadable or unwritable file).

The generated files carry no SPDX gate obligation (validation/
check_spdx_headers.py covers `validation/benchmarks/*.py`, i.e. this file,
`_report.py` and `tiers.py`, not docs/); REPORT.md still opens with the
same SPDX comment as validation/benchmarks/README.md, and REPORT.json names
the licence in its `generated` block.
"""

import argparse
import contextlib
import copy
import datetime
import difflib
import glob
import hashlib
import importlib.util
import io
import json
import math
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if HERE not in sys.path:                  # `import _report` / `import tiers`
    sys.path.insert(0, HERE)

import _report as R                       # noqa: E402
import tiers as T                         # noqa: E402

DEFAULT_OUT = os.path.join(ROOT, "docs", "benchmarking")
PLAN = os.path.join(ROOT, "docs", "benchmarking", "BENCHMARK_PLAN.md")
ENV_SH = os.path.join(ROOT, "env.sh")
MD_NAME, JSON_NAME = "REPORT.md", "REPORT.json"
REPORT_RELS = ("docs/benchmarking/" + MD_NAME, "docs/benchmarking/" + JSON_NAME)

#: opt-in harnesses whatever their module says (the module's own
#: `EXPENSIVE = True` is honoured as well): a row that costs more than 30 s,
#: or one the plan counts as expensive by kind -- t4_pythia_ep_closure is a
#: generator-vs-generator sample and takes ~10 s (README.md, "The opt-in
#: variable")
EXPENSIVE = frozenset({"t4_pythia_ep_closure"})
OPT_IN_ENV = "LIPOLGEN_BENCH"
OPT_IN_REASON = ("expensive row (> 30 s, or a generator-vs-generator sample): "
                 "runs only with %s=1 or make_report.py --all" % OPT_IN_ENV)
SCHEMA = "lipolgen-benchmark-report/1"
SPDX = "GPL-3.0-or-later"
CHILD_TAG = "MAKE_REPORT_ROW "
SYS_PATH_ENV = "MAKE_REPORT_SYS_PATH"
DEFAULT_TIMEOUT = 1800.0
CMD = "python3 validation/benchmarks/make_report.py"
CMD_PREFIX = "cd LiPolGen && source env.sh && "

#: the header fields --check and --compare ignore
VOLATILE_GENERATED = ("date", "git_rev", "git_dirty", "environment")


# ================================================================ discovery

def discover(bench_dir=HERE):
    """Every harness file, sorted by name."""
    return sorted(glob.glob(os.path.join(bench_dir, "t[1-6]_*.py")))


def stem_of(path):
    return os.path.splitext(os.path.basename(path))[0]


def tier_of(stem):
    return "T%s" % stem[1]


def load_module(path):
    """Import a harness file by path (importlib, like the pytests' `_load`),
    registered in sys.modules so that dataclasses/pickle inside it work."""
    name = "lipolgen_bench_" + stem_of(path)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    try:
        spec.loader.exec_module(mod)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    return mod


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def data_files(mod, bench_dir=HERE):
    """{`data/<file>`: sha256 of every byte, or None if not on disk} for each
    vendored file the module names in an upper-case (DATA-like) constant --
    an absolute path into data/, or a bare file name that exists there --
    looking through tuples, lists and dicts up to three levels deep."""
    data_dir = os.path.normpath(os.path.join(bench_dir, "data"))
    try:
        present = set(os.listdir(data_dir))
    except OSError:
        present = set()
    found = set()

    def visit(v, depth):
        if isinstance(v, (str, os.PathLike)):
            s = os.fspath(v)
            if not s or len(s) > 4096 or "\n" in s:
                return
            if os.path.isabs(s):
                p = os.path.normpath(s)
                if p.startswith(data_dir + os.sep):
                    found.add(p)
            elif os.sep not in s and "/" not in s and s in present:
                found.add(os.path.join(data_dir, s))
        elif depth > 0 and isinstance(v, dict):
            for x in v.values():
                visit(x, depth - 1)
        elif depth > 0 and isinstance(v, (list, tuple)):
            for x in v:
                visit(x, depth - 1)

    for k, v in sorted(vars(mod).items()):
        if k.isupper() and not k.startswith("_"):
            visit(v, 3)
    out = {}
    for p in sorted(found):
        rel = "data/" + os.path.relpath(p, data_dir).replace(os.sep, "/")
        out[rel] = sha256_file(p) if os.path.isfile(p) else None
    return out


# ================================================================ one harness

def _base_row(stem):
    return dict(name=stem, tier=tier_of(stem), plan_row=None, status=None,
                generator=None, reference=None, tolerance=None, reason=None,
                config={}, data_sha256={}, runtime_s=None, expensive=False,
                subrows=[], notes=[])


def _fail_row(row, exc):
    row["status"], row["reason"] = R.classify_exception(exc)
    traceback.print_exception(type(exc), exc, exc.__traceback__,
                              file=sys.stderr)
    return row


def run_harness(path, opt_in, bench_dir=HERE):
    """The normalised row of one harness, run in THIS process.  Never
    raises for anything the harness does (KeyboardInterrupt excepted)."""
    stem = stem_of(path)
    row = _base_row(stem)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            mod = load_module(path)
    except (Exception, SystemExit) as exc:
        return _fail_row(row, exc)
    try:
        row["data_sha256"] = data_files(mod, bench_dir)
    except Exception as exc:                        # pragma: no cover
        row["notes"].append("data files not read: %s" % R.first_line(exc))
    row["expensive"] = bool(getattr(mod, "EXPENSIVE", False)) or stem in EXPENSIVE
    if row["expensive"] and not opt_in:
        row["status"], row["reason"] = R.SKIPPED_OPT_IN, OPT_IN_REASON
        return row
    if not callable(getattr(mod, "run", None)):
        row["status"], row["reason"] = R.ERROR, "the module defines no run()"
        return row
    t0 = time.perf_counter()
    try:
        with contextlib.redirect_stdout(buf):
            raw = mod.run(verbose=False)
    except (Exception, SystemExit) as exc:
        row["runtime_s"] = round(time.perf_counter() - t0, 3)
        return _fail_row(row, exc)
    row["runtime_s"] = round(time.perf_counter() - t0, 3)
    try:
        norm = R.normalise_row(raw, stem)
    except Exception as exc:                        # a bug here, not there
        row["status"] = R.ERROR
        row["reason"] = "make_report could not normalise run()'s return: %s" % (
            R.classify_exception(exc)[1])
        return row
    row["notes"].extend(norm.pop("notes"))
    row.update(norm)
    return row


def _signal_name(num):
    try:
        return signal.Signals(num).name
    except ValueError:                              # pragma: no cover
        return "signal %d" % num


def run_harness_isolated(path, opt_in, bench_dir=HERE, timeout=DEFAULT_TIMEOUT):
    """`run_harness` in a child process: a harness that segfaults, exits or
    hangs costs its own row only."""
    stem = stem_of(path)
    env = dict(os.environ)
    env[SYS_PATH_ENV] = json.dumps([p for p in sys.path if p])
    cmd = [sys.executable, os.path.abspath(__file__), "--child", path,
           "--bench-dir", bench_dir, "--child-opt-in", "1" if opt_in else "0"]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, env=env,
                           timeout=timeout)
    except subprocess.TimeoutExpired:
        row = _base_row(stem)
        row["status"] = R.ERROR
        row["reason"] = "harness process timed out after %g s" % timeout
        return row
    line = None
    for l in p.stdout.splitlines():
        if l.startswith(CHILD_TAG):
            line = l[len(CHILD_TAG):]
    if line is not None:
        row = json.loads(line)
        if row.get("status") == R.ERROR and p.stderr.strip():
            sys.stderr.write("---- %s (child stderr, tail) ----\n%s\n" % (
                stem, "\n".join(p.stderr.strip().splitlines()[-25:])))
        return row
    row = _base_row(stem)
    row["status"] = R.ERROR
    if p.returncode < 0:
        row["reason"] = "harness process killed by signal %d (%s)" % (
            -p.returncode, _signal_name(-p.returncode))
    else:
        row["reason"] = ("harness process exited with status %d before "
                         "returning its row" % p.returncode)
    if p.stderr.strip():
        sys.stderr.write("---- %s (child stderr, tail) ----\n%s\n" % (
            stem, "\n".join(p.stderr.strip().splitlines()[-25:])))
    return row


def child_main(path, opt_in, bench_dir):
    """The --child entry: everything the harness prints (Python or C++) is
    sent to stderr; the one row goes to the ORIGINAL stdout, tagged."""
    extra = json.loads(os.environ.get(SYS_PATH_ENV, "[]"))
    sys.path[:0] = [p for p in extra if p not in sys.path]
    sys.stdout.flush()
    saved = os.dup(1)
    os.dup2(2, 1)
    try:
        row = run_harness(path, opt_in, bench_dir)
        text = json.dumps(R.jsonable(row), sort_keys=False)
    finally:
        sys.stdout.flush()
    with os.fdopen(saved, "w") as out:
        out.write(CHILD_TAG + text + "\n")
    return 0


# ================================================================ the plan

def _cells(line):
    parts = re.split(r"(?<!\\)\|", line.strip())
    return [c.strip() for c in parts[1:-1]]


def _section(text, number):
    m = re.search(r"^## %d\..*$" % number, text, re.M)
    if not m:
        return ""
    nxt = re.search(r"^## ", text[m.end():], re.M)
    return text[m.end(): m.end() + nxt.start() if nxt else len(text)]


def read_plan(path=PLAN):
    """What the report takes from BENCHMARK_PLAN.md, READ (never invented):
    the sec. 4 rows (number, benchmark, tier, status column, harnesses named
    in the row) and the sec. 6 U-rows.  None if the plan cannot be read."""
    try:
        with open(path, encoding="utf-8") as f:
            text = f.read()
    except OSError:
        return None
    rows4 = []
    for line in _section(text, 4).splitlines():
        if not re.match(r"^\|\s*\d+\s*\|", line):
            continue
        c = _cells(line)
        names = []
        for n in re.findall(r"\bt[1-6]_[A-Za-z0-9_]+", line):
            n = n[:-3] if n.endswith(".py") else n
            if n not in names:
                names.append(n)
        rows4.append(dict(row=int(c[0]), benchmark=c[1], tier=c[2],
                          status=c[-1], harnesses=names))
    u_rows = []
    for line in _section(text, 6).splitlines():
        if re.match(r"^\|\s*U-\d+\s*\|", line):
            c = _cells(line)
            u = dict(id=c[0], observable=c[1], closest_proxy=c[2],
                     distance=c[-1])
            if len(c) > 4:
                # An unescaped '|' in the plan's own text (U-3's "a2(|t|)"
                # until the plan escaped it, 2026-09-27) splits the row into
                # more cells than the table has columns; the surplus is
                # re-joined, pipe escaped, into the closest-proxy column.
                u["closest_proxy"] = "\\|".join(c[2:-1])
                u["parse_note"] = (
                    "the plan's row has %d cells for 4 columns (an unescaped "
                    "'|' in its text); the surplus was re-joined into the "
                    "closest-proxy column" % len(c))
            u_rows.append(u)
    return dict(rows4=rows4, u_rows=u_rows)


def plan_row_map(plan):
    out = {}
    for r in (plan or {}).get("rows4", []):
        for h in r["harnesses"]:
            out.setdefault(h, r["row"])
    return out


# ================================================================ environment

def _cmd(args):
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=30)
    except Exception:
        return None
    return p.stdout.strip() if p.returncode == 0 and p.stdout.strip() else None


def _module_version(name):
    try:
        mod = importlib.import_module(name)
    except Exception as exc:
        return "not importable (%s)" % R.first_line(exc)
    return str(getattr(mod, "__version__", "unknown"))


def lhapdf_search_path():
    dirs = [d for d in os.environ.get("LHAPDF_DATA_PATH", "").split(os.pathsep)
            if d]
    d = _cmd(["lhapdf-config", "--datadir"])
    if d:
        dirs.append(d)
    out = []
    for d in dirs:
        d = os.path.normpath(d)
        if d not in out:
            out.append(d)
    return out


def lhapdf_sets(dirs):
    """Every set in the search path: a directory holding `<name>.info`."""
    sets = set()
    for d in dirs:
        try:
            for n in os.listdir(d):
                if os.path.isfile(os.path.join(d, n, n + ".info")):
                    sets.add(n)
        except OSError:
            continue
    return sorted(sets)


def pythia_version():
    v = _cmd(["pythia8-config", "--version"])
    if v:
        return v, "pythia8-config --version"
    xml = os.path.join(os.environ.get("PYTHIA8DATA", ""), "Version.xml")
    try:
        with open(xml, encoding="utf-8") as f:
            m = re.search(r'Pythia:versionNumber"\s+default="([0-9.]+)"', f.read())
        if m:
            return m.group(1), "$PYTHIA8DATA/Version.xml"
    except OSError:
        pass
    return None, "not found"


def pinned_pythia():
    """The PYTHIA version env.sh's header says the project pins, if any."""
    try:
        with open(ENV_SH, encoding="utf-8") as f:
            m = re.search(r"PYTHIA\s+(\d+\.\d+)", f.readline())
        return m.group(1) if m else None
    except OSError:
        return None


def environment():
    """The environment the rows were measured in (REPORT.json `generated.
    environment`; ignored by --check, compared by the bench test)."""
    lg = {}
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            import lipolgen as _lg
        lg = dict(version=str(getattr(_lg, "__version__", "unknown")),
                  have_lhapdf=bool(getattr(_lg, "HAVE_LHAPDF", False)),
                  have_pythia8=bool(getattr(_lg, "HAVE_PYTHIA8", False)))
    except Exception as exc:
        lg = dict(version="not importable (%s)" % R.first_line(exc))
    pv, psrc = pythia_version()
    dirs = lhapdf_search_path()
    return dict(
        python=platform.python_version(),
        numpy=_module_version("numpy"),
        scipy=_module_version("scipy"),
        pyhepmc=_module_version("pyhepmc"),
        lipolgen=lg,
        pythia=pv, pythia_source=psrc, pythia_pinned=pinned_pythia(),
        lhapdf=_cmd(["lhapdf-config", "--version"]),
        hepmc3=_cmd(["HepMC3-config", "--version"]),
        lhapdf_search_path=dirs,
        lhapdf_sets_present=lhapdf_sets(dirs),
    )


def git_info(root=ROOT):
    """(short revision, dirty?) -- read-only (`--no-optional-locks`); the
    report files themselves do not count as dirt."""
    def git(*args):
        try:
            p = subprocess.run(["git", "--no-optional-locks"] + list(args),
                               cwd=root, capture_output=True, text=True,
                               timeout=60)
        except Exception:
            return None
        return p.stdout if p.returncode == 0 else None
    rev = (git("rev-parse", "--short", "HEAD") or "").strip() or "unknown"
    status = git("status", "--porcelain")
    if status is None:
        return rev, None
    dirty = False
    for l in status.splitlines():
        path = l[3:].split(" -> ")[-1].strip().strip('"')
        if path not in REPORT_RELS:
            dirty = True
            break
    return rev, dirty


# ================================================================ generate

def summarise(rows):
    out = {s: 0 for s in R.REPORT_STATUSES}
    for r in rows:
        out[r["status"]] = out.get(r["status"], 0) + 1
    return out


def blocked_list(rows):
    out = []
    for r in rows:
        if r["status"] in ("blocked", R.ENV_BLOCKED):
            out.append(dict(row=r["name"], status=r["status"],
                            reason=r["reason"]))
        for s in r.get("subrows", []):
            if s["status"] == "blocked":
                out.append(dict(row="%s[%s]" % (r["name"], s["name"]),
                                status="blocked", reason=s["reason"]))
    return out


def command_line(opt_in, only):
    c = CMD + (" --all" if opt_in else "")
    if only:
        c += " --only " + ",".join(sorted(set(only)))
    return c


def generate(bench_dir=HERE, opt_in=False, only=None, isolate=True,
             timeout=DEFAULT_TIMEOUT, plan_path=PLAN, now=None, git=None,
             env=None, echo=True):
    """The whole report as a dict (REPORT.json's content).  `now`, `git`
    ((rev, dirty)) and `env` may be injected (tests, --check)."""
    paths = discover(bench_dir)
    if only:
        by = {stem_of(p): p for p in paths}
        unknown = [n for n in only if n not in by]
        if unknown:
            raise ValueError("no such harness: %s (have: %s)" % (
                ", ".join(unknown), ", ".join(sorted(by))))
        paths = [by[n] for n in sorted(set(only))]
    plan = read_plan(plan_path)
    prow = plan_row_map(plan)
    rows = []
    for p in paths:
        row = (run_harness_isolated(p, opt_in, bench_dir, timeout) if isolate
               else run_harness(p, opt_in, bench_dir))
        row = R.jsonable(row)
        row["plan_row"] = prow.get(row["name"])
        rows.append(row)
        if echo:
            print(R.report_line(row["name"], row["generator"] or "-",
                                row["reference"] or "-",
                                row["tolerance"] or "-", row["status"])
                  + ("  (%s)" % row["reason"]
                     if row["status"] in (R.ERROR, R.ENV_BLOCKED,
                                          R.SKIPPED_OPT_IN) else ""))
            sys.stdout.flush()
    rev, dirty = git if git is not None else git_info()
    if now is None:
        now = datetime.datetime.now(datetime.timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ")
    generated = dict(
        notice=("GENERATED by validation/benchmarks/make_report.py -- do not "
                "edit; regenerate (BENCHMARK_PLAN.md sec. 3)"),
        spdx_license_identifier=SPDX,
        date=now, git_rev=rev, git_dirty=dirty,
        command=command_line(opt_in, only), opt_in=bool(opt_in),
        only=sorted(set(only)) if only else None,
        plan=None if plan is None else "docs/benchmarking/BENCHMARK_PLAN.md",
        environment=env if env is not None else environment(),
    )
    pnw = dict(
        blocked=blocked_list(rows),
        plan_rows_without_harness=[
            dict(row=r["row"], benchmark=r["benchmark"], tier=r["tier"],
                 status=r["status"])
            for r in (plan or {}).get("rows4", []) if not r["harnesses"]],
        u_rows=(plan or {}).get("u_rows", []),
    )
    return R.jsonable(dict(schema=SCHEMA, generated=generated,
                           summary=summarise(rows), rows=rows,
                           planned_not_wired=pnw))


# ================================================================ rendering

def md_text(s):
    """Plain text made safe for a markdown table cell / list item."""
    if s is None:
        return "—"
    s = R.one_line(s)
    for a, b in (("\\", "\\\\"), ("|", "\\|"), ("<", "&lt;"), (">", "&gt;"),
                 ("*", "\\*")):
        s = s.replace(a, b)
    return s


def md_status(status, in_verdict=None):
    s = "**%s**" % str(status).upper()
    if in_verdict is False:
        s += " (recorded, not in the verdict)"
    return s


def _val(v, limit=120):
    if isinstance(v, str):
        return v if len(v) <= limit else v[:limit - 1] + "…"
    return R.fmt(v)


def md_values(values):
    if not values:
        return "—"
    return md_text("; ".join("%s = %s" % (k, _val(v)) for k, v in values.items()))


def fmt_runtime(r):
    return "—" if r.get("runtime_s") is None else "%.2f s" % r["runtime_s"]


def render_md(rep):
    g, env = rep["generated"], rep["generated"]["environment"]
    out = []
    w = out.append
    w("<!-- SPDX-License-Identifier: %s -->" % SPDX)
    w("<!-- GENERATED FILE -- written by validation/benchmarks/make_report.py; "
      "DO NOT EDIT. Regenerate it instead. -->")
    w("# Benchmark report")
    w("")
    w("> **Generated — do not edit.** Written by "
      "`validation/benchmarks/make_report.py` from the harnesses "
      "`validation/benchmarks/t[1-6]_*.py`, together with `REPORT.json` "
      "beside it (the same rows, machine-readable: the M2 baseline). "
      "`BENCHMARK_PLAN.md` §3: *the report is regenerated, never "
      "hand-edited* — `make_report.py --check` fails on any difference "
      "other than the date, the revision, a runtime or the environment.")
    w("")
    if g.get("only"):
        w("**PARTIAL REPORT:** only %s." % ", ".join("`%s`" % n for n in g["only"]))
        w("")
    w("- **generated:** %s" % g["date"])
    w("- **git revision:** `%s`%s" % (
        g["git_rev"], {True: " + uncommitted changes (dirty)", False: " (clean)",
                       None: " (dirty state unknown)"}[g["git_dirty"]]))
    lg = env.get("lipolgen") or {}
    pin = env.get("pythia_pinned")
    w("- **environment:** Python %s; numpy %s; scipy %s; PYTHIA %s (%s%s); "
      "LHAPDF %s; HepMC3 %s (pyhepmc %s); lipolgen %s (LHAPDF tier %s, "
      "PYTHIA tier %s)" % (
          env.get("python"), env.get("numpy"), env.get("scipy"),
          env.get("pythia") or "not found", env.get("pythia_source"),
          "; env.sh pins %s" % pin if pin and pin != env.get("pythia") else "",
          env.get("lhapdf") or "not found", env.get("hepmc3") or "not found",
          env.get("pyhepmc"), lg.get("version"),
          "on" if lg.get("have_lhapdf") else "off",
          "on" if lg.get("have_pythia8") else "off"))
    sets = env.get("lhapdf_sets_present") or []
    w("- **LHAPDF sets present:** %s (searched: %s)" % (
        ", ".join("`%s`" % s for s in sets) if sets else "**none**",
        ", ".join("`%s`" % d for d in env.get("lhapdf_search_path") or [])
        or "nothing"))
    w("- **opt-in (expensive) rows:** %s" % (
        "run (`--all` / `%s=1`)" % OPT_IN_ENV if g["opt_in"] else
        "**not run** — listed as skipped (opt-in); `--all` or `%s=1` runs "
        "them" % OPT_IN_ENV))
    w("- **this file was made by:** `%s%s`" % (CMD_PREFIX, g["command"]))
    w("- **check it:** `%s%s --check`; **M2 compare:** `%s%s --all --compare "
      "docs/benchmarking/REPORT.json`" % (CMD_PREFIX, CMD, CMD_PREFIX, CMD))
    w("")
    s = rep["summary"]
    w("**%d harness%s:** %s." % (
        len(rep["rows"]), "" if len(rep["rows"]) == 1 else "es",
        ", ".join("%d %s" % (n, k) for k, n in s.items())))
    w("")
    w("## The rows")
    w("")
    w("One row per harness. *plan row* is the `BENCHMARK_PLAN.md` §4 row "
      "number when the harness is one of the first ten. The generator column "
      "is the harness's own value text; the sub-rows and the configuration "
      "are under *Detail*.")
    w("")
    w("| harness | tier | plan row | status | generator (the number) | "
      "reference | tolerance | runtime |")
    w("|---|---|---|---|---|---|---|---|")
    for r in rep["rows"]:
        gen = md_text(r["generator"])
        if r["generator"] is None and r.get("reason"):
            # no value came back (error / environment / opt-in): the reason
            # is the row's content
            gen = "— (*%s*)" % md_text(r["reason"])
        w("| [`%s`](#%s) | %s | %s | %s | %s | %s | %s | %s |" % (
            r["name"], r["name"].lower(), r["tier"],
            "—" if r["plan_row"] is None else r["plan_row"],
            md_status(r["status"]), gen,
            md_text(r["reference"]), md_text(r["tolerance"]),
            fmt_runtime(r)))
    w("")
    w("## Detail")
    for r in rep["rows"]:
        w("")
        w("### %s" % r["name"])
        w("")
        w("- **status:** %s%s" % (md_status(r["status"]),
                                  "; expensive (opt-in)" if r.get("expensive")
                                  else ""))
        if r.get("reason"):
            w("- **reason:** %s" % md_text(r["reason"]))
        for n in r.get("notes") or []:
            w("- **note:** %s" % md_text(n))
        data = r.get("data_sha256") or {}
        if data:
            w("- **vendored data (sha256 of every byte):** %s" % "; ".join(
                "`%s` %s" % (k, "`%s`" % v if v else "**not on disk**")
                for k, v in data.items()))
        if r.get("config"):
            w("- **config:** %s" % md_values(r["config"]))
        subs = r.get("subrows") or []
        if subs:
            w("")
            w("| sub-row | status | values / reason |")
            w("|---|---|---|")
            for sr in subs:
                body = md_values(sr.get("values"))
                if sr.get("reason"):
                    body = "**reason:** %s%s" % (
                        md_text(sr["reason"]),
                        "" if body == "—" else "; " + body)
                w("| `%s` | %s | %s |" % (sr["name"].replace("|", "\\|"),
                                          md_status(sr["status"],
                                                    sr.get("in_verdict")),
                                          body))
    pnw = rep["planned_not_wired"]
    w("")
    w("## Planned, not wired")
    w("")
    w("Nothing below is a measurement. The first table is this run's BLOCKED "
      "rows and sub-rows, with the reason each one gives; the other two are "
      "read verbatim from `BENCHMARK_PLAN.md` (§4's rows that have no "
      "harness, and §6, what cannot be benchmarked until data exists).")
    w("")
    w("### Blocked in this run")
    w("")
    if pnw["blocked"]:
        w("| row | status | reason |")
        w("|---|---|---|")
        for b in pnw["blocked"]:
            w("| `%s` | %s | %s |" % (b["row"], md_status(b["status"]),
                                     md_text(b["reason"])))
    else:
        w("None.")
    w("")
    w("### `BENCHMARK_PLAN.md` §4 rows with no harness (status column, verbatim)")
    w("")
    if pnw["plan_rows_without_harness"]:
        w("| # | benchmark | tier | status (the plan's words) |")
        w("|---|---|---|---|")
        for p in pnw["plan_rows_without_harness"]:
            w("| %s | %s | %s | %s |" % (p["row"], p["benchmark"], p["tier"],
                                        p["status"]))
    else:
        w("None (or the plan could not be read).")
    w("")
    w("### `BENCHMARK_PLAN.md` §6 — what cannot be benchmarked until data "
      "exists (verbatim)")
    w("")
    if pnw["u_rows"]:
        w("| # | observable | closest proxy | distance |")
        w("|---|---|---|---|")
        for u in pnw["u_rows"]:
            w("| %s | %s | %s | %s |" % (u["id"], u["observable"],
                                        u["closest_proxy"], u["distance"]))
        for u in pnw["u_rows"]:
            if u.get("parse_note"):
                w("")
                w("*%s:* %s." % (u["id"], md_text(u["parse_note"])))
    else:
        w("The plan could not be read.")
    w("")
    return "\n".join(out)


def dump_json(rep):
    return json.dumps(rep, indent=1, ensure_ascii=False) + "\n"


def write(rep, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    md, js = os.path.join(out_dir, MD_NAME), os.path.join(out_dir, JSON_NAME)
    with open(js, "w", encoding="utf-8") as f:
        f.write(dump_json(rep))
    with open(md, "w", encoding="utf-8") as f:
        f.write(render_md(rep))
    return md, js


# ================================================================ compare

_NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def headline_numbers(text):
    """The numbers printed in a generator field, in order."""
    return [float(x) for x in _NUM.findall(text or "")]


def _num_moved(x, y, rel):
    """A purely RELATIVE move (no absolute floor, so a 1e-14 headline that
    doubles is a move; until 2026-09-27 the scale was max(1, |x|, |y|), which
    missed any change below ~1e-9).  0 == 0 and nan, nan are not moves; a
    change to or from nan / inf is."""
    x, y = float(x), float(y)
    if x == y or (math.isnan(x) and math.isnan(y)):
        return False
    if not (math.isfinite(x) and math.isfinite(y)):
        return True
    return abs(x - y) > rel * max(abs(x), abs(y))


def _moved(a, b, rel=1e-9):
    if len(a) != len(b):
        return True
    return any(_num_moved(x, y, rel) for x, y in zip(a, b))


def _numeric(values):
    return {k: v for k, v in (values or {}).items()
            if isinstance(v, (int, float)) and not isinstance(v, bool)}


#: statuses under which a row carries no value of its own
_NO_VALUE = (R.SKIPPED_OPT_IN, R.ERROR, R.ENV_BLOCKED)


def compare(old, new, out=None):
    """Print every change of status or headline number between two reports
    (dicts); return the number of changes that count (the M2 rule)."""
    out = out or sys.stdout
    w = lambda s: print(s, file=out)                # noqa: E731
    changes = 0
    oldrows = {r["name"]: r for r in old.get("rows", [])}
    newrows = {r["name"]: r for r in new.get("rows", [])}
    partial = bool(new.get("generated", {}).get("only"))
    for name in sorted(set(oldrows) | set(newrows)):
        o, n = oldrows.get(name), newrows.get(name)
        if n is None:
            if not partial:
                w("GONE     %s (in the baseline as %s; not in this report)"
                  % (name, str(o["status"]).upper()))
                changes += 1
            continue
        if o is None:
            w("NEW      %s: %s (not in the baseline; informational)"
              % (name, str(n["status"]).upper()))
            continue
        if o["status"] != n["status"]:
            hint = ("  -- not run here; rerun with --all"
                    if n["status"] == R.SKIPPED_OPT_IN else
                    "  -- %s" % n.get("reason") if n["status"] in _NO_VALUE
                    else "")
            w("STATUS   %s: %s -> %s%s" % (name, str(o["status"]).upper(),
                                           str(n["status"]).upper(), hint))
            changes += 1
            if o["status"] in _NO_VALUE or n["status"] in _NO_VALUE:
                continue        # one side has no numbers: nothing to move
        a, b = headline_numbers(o.get("generator")), headline_numbers(n.get("generator"))
        if _moved(a, b):
            w("MOVED    %s:\n    was: %s\n    now: %s" % (
                name, o.get("generator"), n.get("generator")))
            changes += 1
        osub = {s["name"]: s for s in o.get("subrows", [])}
        nsub = {s["name"]: s for s in n.get("subrows", [])}
        for sn in sorted(set(osub) | set(nsub)):
            os_, ns_ = osub.get(sn), nsub.get(sn)
            label = "%s[%s]" % (name, sn)
            if ns_ is None or os_ is None:
                w("SUBROW   %s: %s" % (label, "gone" if ns_ is None else
                                       "new (informational)"))
                changes += ns_ is None
                continue
            if os_["status"] != ns_["status"]:
                w("SUBROW   %s: %s -> %s" % (label, str(os_["status"]).upper(),
                                             str(ns_["status"]).upper()))
                changes += 1
            ov, nv = _numeric(os_.get("values")), _numeric(ns_.get("values"))
            moved = [k for k in sorted(set(ov) | set(nv))
                     if k not in ov or k not in nv or _moved([ov[k]], [nv[k]])]
            if moved:
                w("SUBMOVED %s: %s" % (label, "; ".join(
                    "%s %s -> %s" % (k, R.fmt(ov.get(k)), R.fmt(nv.get(k)))
                    for k in moved)))
                changes += 1
    w("compare: %d change%s (the M2 rule: any change exits 1)" % (
        changes, "" if changes == 1 else "s"))
    return changes


# ================================================================ check

def overlay_volatile(new, old):
    """`new` with the fields --check ignores (the date, the revision, the
    environment, every runtime) taken from `old`."""
    ov = copy.deepcopy(new)
    for k in VOLATILE_GENERATED:
        if k in old.get("generated", {}):
            ov["generated"][k] = copy.deepcopy(old["generated"][k])
    oldrt = {r["name"]: r.get("runtime_s") for r in old.get("rows", [])}
    for r in ov["rows"]:
        if r["name"] in oldrt:
            r["runtime_s"] = oldrt[r["name"]]
    return ov


def _diff(a, b, fa, fb, limit=200):
    d = list(difflib.unified_diff(a.splitlines(), b.splitlines(), fa, fb,
                                  lineterm="", n=1))
    return d[:limit] + (["... (%d more diff lines)" % (len(d) - limit)]
                        if len(d) > limit else [])


def check(out_dir=DEFAULT_OUT, bench_dir=HERE, isolate=True,
          timeout=DEFAULT_TIMEOUT, plan_path=PLAN, out=None):
    """0 when the committed report is what regeneration gives (modulo the
    volatile fields), 1 when not, 2 when there is no committed report."""
    out = out or sys.stdout
    md_path, js_path = os.path.join(out_dir, MD_NAME), os.path.join(out_dir, JSON_NAME)
    try:
        with open(js_path, encoding="utf-8") as f:
            old = json.load(f)
        with open(md_path, encoding="utf-8") as f:
            old_md = f.read()
    except (OSError, ValueError) as exc:
        print("check: cannot read the committed report: %s" % exc, file=out)
        return 2
    og = old.get("generated", {})
    if og.get("only"):
        print("check: the committed report is PARTIAL (only %s); regenerate "
              "it in full" % og["only"], file=out)
        return 1
    new = generate(bench_dir, opt_in=bool(og.get("opt_in")), isolate=isolate,
                   timeout=timeout, plan_path=plan_path, echo=False)
    tmp = tempfile.mkdtemp(prefix="make_report_check_")
    write(new, tmp)
    print("check: regenerated into %s (opt-in rows %s, as committed)" % (
        tmp, "run" if og.get("opt_in") else "skipped"), file=out)
    if new["generated"]["environment"] != og.get("environment"):
        print("check: note -- the environment differs from the committed "
              "report's (ignored here; rows may differ because of it)",
              file=out)
    ov = overlay_volatile(new, old)
    dj = _diff(dump_json(old), dump_json(ov), js_path, "regenerated/" + JSON_NAME)
    dm = _diff(old_md, render_md(ov), md_path, "regenerated/" + MD_NAME)
    for line in dj + dm:
        print(line, file=out)
    if dj or dm:
        print("check: DIFFERS -- REPORT.json %s, REPORT.md %s; the "
              "regenerated files are kept in %s" % (
                  "differs" if dj else "same", "differs" if dm else "same",
                  tmp), file=out)
        return 1
    shutil.rmtree(tmp, ignore_errors=True)
    print("check: OK -- the committed REPORT.md / REPORT.json are what "
          "regeneration gives (date, revision, runtimes, environment ignored)",
          file=out)
    return 0


# ================================================================ CLI

def _names(values):
    out = []
    for v in values or []:
        for n in v.split(","):
            n = n.strip()
            if n:
                out.append(n[:-3] if n.endswith(".py") else n)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Regenerate docs/benchmarking/REPORT.md and REPORT.json "
                    "from validation/benchmarks/t[1-6]_*.py (see the module "
                    "docstring).")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--only", action="append", metavar="NAMES")
    ap.add_argument("--since", metavar="REV")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--compare", metavar="OLD_JSON")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--in-process", action="store_true")
    ap.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT)
    ap.add_argument("--bench-dir", default=HERE, help=argparse.SUPPRESS)
    ap.add_argument("--plan", default=PLAN, help=argparse.SUPPRESS)
    ap.add_argument("--child", metavar="PATH", help=argparse.SUPPRESS)
    ap.add_argument("--child-opt-in", default="0", help=argparse.SUPPRESS)
    try:
        a = ap.parse_args(argv)
    except SystemExit as exc:
        return 0 if exc.code == 0 else 2
    if a.child:
        return child_main(a.child, a.child_opt_in == "1", a.bench_dir)
    opt_in = a.all or os.environ.get(OPT_IN_ENV) == "1"
    isolate = not a.in_process
    if a.check:
        if a.only or a.since or a.compare:
            print("make_report: --check takes no --only/--since/--compare",
                  file=sys.stderr)
            return 2
        return check(a.out_dir or DEFAULT_OUT, a.bench_dir, isolate,
                     a.timeout, a.plan)
    only = _names(a.only)
    if a.since:
        try:
            files = T.changed_since(a.since)
            plan = T.plan_for(files, a.bench_dir)
        except Exception as exc:
            print("make_report: --since: %s" % exc, file=sys.stderr)
            return 2
        T.print_plan(plan, files, out=sys.stderr)
        only = sorted(set(only) | set(plan["harnesses"]))
        if not only:
            print("make_report: nothing to re-run since %s" % a.since)
            return 0
    if only and a.out_dir is None and not a.compare:
        print("make_report: a partial report (--only/--since) is never "
              "written over docs/benchmarking/REPORT.*; pass --out-dir",
              file=sys.stderr)
        return 2
    old = None
    if a.compare:
        try:
            with open(a.compare, encoding="utf-8") as f:
                old = json.load(f)
        except (OSError, ValueError) as exc:
            print("make_report: --compare: %s" % exc, file=sys.stderr)
            return 2
    try:
        rep = generate(a.bench_dir, opt_in=opt_in, only=only or None,
                       isolate=isolate, timeout=a.timeout, plan_path=a.plan)
    except ValueError as exc:
        print("make_report: %s" % exc, file=sys.stderr)
        return 2
    if a.out_dir is not None or old is None:
        try:
            md, js = write(rep, a.out_dir or DEFAULT_OUT)
        except OSError as exc:
            print("make_report: cannot write: %s" % exc, file=sys.stderr)
            return 2
        print("wrote %s and %s" % (md, js))
    print("summary: %s" % ", ".join("%d %s" % (n, k)
                                    for k, n in rep["summary"].items()))
    if old is not None:
        return 1 if compare(old, rep) else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
