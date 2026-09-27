#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Shared helpers for the benchmark harnesses and `make_report.py`.

Not a harness (its name does not match `t[1-6]_*.py`, so `make_report.py`
does not discover it).  Three things live here so that they are written once:

1. **numpy-scalar coercion** -- `to_builtin()` turns numpy scalars and arrays
   (and tuples, sets, paths) into plain Python, recursively.  Under numpy >= 2
   `repr(np.float64(1.0))` is `np.float64(1.0)`, which is how the DJANGOH
   harness's `ROW |` line came to print `np.float64(...)` (review 2026-09-26);
   `fmt()` formats a value for a REPORT line after coercing it.
2. **the REPORT line** -- `report_line()` / `emit()` print the one line of
   validation/benchmarks/README.md convention 4:
   `REPORT | name | generator value | reference | tolerance | STATUS`.
3. **row normalisation** -- `normalise_row()` maps whatever a harness's
   `run()` returned onto the one schema `make_report.py` writes to
   docs/benchmarking/REPORT.json.  The harnesses were written by different
   hands and disagree (DJANGOH returns `generator_value` / `reference_value`;
   HERMES and the T1 rows carry per-config sub-rows under `configs`; the
   deuteron statics under `models`; the de Vries rms under `de_vries`; T5 sum
   rules and T6 under `subrows`); the adaptation is done HERE, never by
   editing a harness.

New harnesses MAY `import _report` (put `validation/benchmarks` on
`sys.path` first, as `make_report.py` does) and call `emit(row)`; existing
harnesses are not required to.  Nothing here changes a number: the helpers
format and re-shape, and the only rounding (`round_sig`, 10 significant
digits) is applied to the copy written into REPORT.json, never to a verdict.

Run:  nothing to run; `python -m pytest -p no:cacheprovider
python/tests/test_bench_report.py` exercises it.
"""

import math
import os
import re

#: what a harness's `run()` may return as its status
HARNESS_STATUSES = ("pass", "fail", "blocked")
#: what make_report.py may write for a row: the harness statuses plus the
#: three the generator itself assigns
ENV_BLOCKED = "blocked (environment)"
ERROR = "error"
SKIPPED_OPT_IN = "skipped (opt-in)"
REPORT_STATUSES = HARNESS_STATUSES + (ENV_BLOCKED, ERROR, SKIPPED_OPT_IN)

#: LHAPDF's message when a set is not installed; an exception carrying it is
#: an environment fact, not a broken harness
LHAPDF_MISSING = "Info file not found for PDF set"
_LHAPDF_SET = re.compile(r"Info file not found for PDF set '?([^'\s]+)'?")

#: significant digits kept for numbers written into REPORT.json
SIG_DIGITS = 10

#: keys never copied into a row's `config` or a sub-row's `values`: wall-clock
#: times differ from run to run, and `--check` must reproduce the file
VOLATILE_KEY = re.compile(r"(runtime|seconds|elapsed|wall_?time|build_s|time_s)",
                          re.IGNORECASE)

#: containers of per-config / per-model sub-rows, in the order they are read
#: (the harness each one serves: HERMES, t1_* -> configs; t2_deuteron_static ->
#: models; t3_li6_rms_devries -> de_vries; t5_sum_rules, t6_*, t4_pythia_* ->
#: subrows)
SUBROW_CONTAINERS = ("subrows", "configs", "models", "de_vries")
#: a sub-row may nest one more level (t1_g1d_world: config -> experiments)
SUBROW_NESTED = ("experiments",)
#: keys that name a sub-row held in a list
SUBROW_NAME_KEYS = ("name", "clause", "source_key", "key", "id", "label")
#: two keys that name a sub-row together when none of the above is present
#: (t4_pythia_ep_closure: window + observable)
SUBROW_COMPOSITE_KEYS = ("window", "observable")


# ------------------------------------------------------------ coercion

def _numpy():
    try:
        import numpy
        return numpy
    except Exception:                               # pragma: no cover
        return None


def to_builtin(obj):
    """`obj` with every numpy scalar/array, tuple, set and path turned into
    plain Python (float, int, bool, str, list, dict with str keys),
    recursively.  Sets are sorted when their items allow it, so the result
    is deterministic."""
    np = _numpy()
    if np is not None:
        if isinstance(obj, np.generic):
            return to_builtin(obj.item())
        if isinstance(obj, np.ndarray):
            return [to_builtin(x) for x in obj.tolist()]
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, os.PathLike):
        return os.fspath(obj)
    if isinstance(obj, dict):
        return {str(k) if not isinstance(k, str) else k: to_builtin(v)
                for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_builtin(x) for x in obj]
    if isinstance(obj, (set, frozenset)):
        items = [to_builtin(x) for x in obj]
        try:
            return sorted(items)
        except TypeError:
            return sorted(items, key=repr)
    if isinstance(obj, bytes):
        return obj.decode("utf-8", "replace")
    return str(obj)


def is_scalar(v):
    return v is None or isinstance(v, (bool, int, float, str))


def round_sig(v, digits=SIG_DIGITS):
    """A float kept to `digits` significant digits; a non-finite float as the
    string 'nan' / 'inf' / '-inf' (strict JSON has no NaN).  Anything else is
    returned unchanged."""
    if isinstance(v, bool) or not isinstance(v, float):
        return v
    if not math.isfinite(v):
        return "nan" if math.isnan(v) else ("inf" if v > 0 else "-inf")
    return float("%.*g" % (digits, v))


def jsonable(obj):
    """`to_builtin(obj)` with every float passed through `round_sig`."""
    obj = to_builtin(obj)
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [jsonable(x) for x in obj]
    return round_sig(obj)


def fmt(v):
    """One value as REPORT-line text: numpy coerced first; floats `%.6g`;
    lists joined with ' / '; dicts as 'k: v; k: v'; None as 'n/a'."""
    v = to_builtin(v)
    if v is None:
        return "n/a"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        if not math.isfinite(v):
            return "nan" if math.isnan(v) else ("inf" if v > 0 else "-inf")
        return "%.6g" % v
    if isinstance(v, (int, str)):
        return str(v)
    if isinstance(v, list):
        return " / ".join(fmt(x) for x in v)
    if isinstance(v, dict):
        return "; ".join("%s: %s" % (k, fmt(x)) for k, x in v.items())
    return str(v)                                   # pragma: no cover


def one_line(text):
    """`text` with every run of whitespace (newlines included) folded to one
    space, so a field can never break a `|`-separated line."""
    return " ".join(str(text).split())


# ------------------------------------------------------------ the REPORT line

def report_line(name, generator, reference, tolerance, status):
    """validation/benchmarks/README.md convention 4, numpy-safe."""
    return "REPORT | %s | %s | %s | %s | %s" % tuple(
        one_line(fmt(x)) for x in (name, generator, reference, tolerance,
                                   str(status).upper()))


def emit(row, stream=None):
    """Print `row`'s REPORT line (to `stream`, default stdout) and return it.
    Accepts DJANGOH-style `generator_value` / `reference_value` as well."""
    line = report_line(
        row.get("name"),
        row["generator"] if row.get("generator") is not None
        else row.get("generator_value"),
        row["reference"] if row.get("reference") is not None
        else row.get("reference_value"),
        row.get("tolerance"), row.get("status"))
    print(line, file=stream)
    return line


# ------------------------------------------------------------ exceptions

def lhapdf_missing_set(message):
    """The LHAPDF set a message says is not installed, else None."""
    if LHAPDF_MISSING not in str(message):
        return None
    m = _LHAPDF_SET.search(str(message))
    return m.group(1) if m else "(name not in the message)"


def first_line(exc):
    text = str(exc).strip()
    return text.splitlines()[0] if text else ""


def classify_exception(exc):
    """(status, reason) for an exception raised by importing or running a
    harness.  A missing LHAPDF set is `blocked (environment)` with the
    harness convention's reason text; anything else is `error` with the
    exception's type and FIRST line."""
    missing = lhapdf_missing_set(exc)
    if missing is not None:
        return ENV_BLOCKED, "environment: LHAPDF set %s not installed" % missing
    fl = first_line(exc)
    return ERROR, ("%s: %s" % (type(exc).__name__, fl) if fl
                   else type(exc).__name__)


# ------------------------------------------------------------ normalisation

def flatten_scalars(d, depth=1, _prefix=""):
    """The scalar fields of dict `d` (and of its nested dicts, `depth` levels
    down, as 'outer.inner'), plus short lists (<= 8) of scalars; longer lists
    and deeper structures are left to the harness's own output.  Volatile
    (timing) keys are dropped.  Values are JSON-ready (`jsonable`)."""
    out = {}
    if not isinstance(d, dict):
        return out
    for k, v in d.items():
        k = str(k)
        if VOLATILE_KEY.search(k):
            continue
        v = to_builtin(v)
        if is_scalar(v):
            out[_prefix + k] = round_sig(v)
        elif isinstance(v, list) and len(v) <= 8 and all(is_scalar(x) for x in v):
            out[_prefix + k] = [round_sig(x) for x in v]
        elif isinstance(v, dict) and depth > 0:
            out.update(flatten_scalars(v, depth - 1, _prefix + k + "."))
    return out


def _subrow_status(d):
    """The status of a sub-row dict: `status`, else the first `*_status`
    field, else a boolean `within` / `ok` read as pass/fail."""
    s = d.get("status")
    if s is None:
        for k, v in d.items():
            if str(k).endswith("_status") and isinstance(v, str):
                s = v
                break
    if s is None:
        for k in ("within", "ok"):
            if isinstance(d.get(k), bool):
                s = "pass" if d[k] else "fail"
                break
    return None if s is None else str(to_builtin(s)).strip().lower()


def _subrow_name(d, index):
    for k in SUBROW_NAME_KEYS:
        if isinstance(d.get(k), str) and d[k]:
            return d[k]
    parts = [d[k] for k in SUBROW_COMPOSITE_KEYS if isinstance(d.get(k), str)]
    if parts:
        return "/".join(parts)
    return "#%d" % index


def _items(container):
    if isinstance(container, dict):
        return [(str(k), v) for k, v in container.items()]
    if isinstance(container, (list, tuple)):
        return [(_subrow_name(v, i) if isinstance(v, dict) else "#%d" % i, v)
                for i, v in enumerate(container)]
    return []


def _subrow(name, d, source):
    status = _subrow_status(d)
    reason = d.get("reason")
    row = dict(name=name, source=source,
               status=status if status is not None else "n/a",
               reason=None if reason is None else one_line(fmt(reason)))
    if isinstance(d.get("in_verdict"), bool):
        row["in_verdict"] = d["in_verdict"]
    skip = {"status", "reason", "in_verdict"} | set(SUBROW_NESTED)
    body = {k: v for k, v in d.items() if k not in skip}
    row["values"] = flatten_scalars(body)
    if not row["values"]:
        # nothing at the first two levels (t5_sum_rules' Bjorken sub-row
        # keeps its numbers under per_q2 -> Q2 -> reference/verdict): look
        # three levels further before recording an empty sub-row
        row["values"] = flatten_scalars(body, depth=3)
    return row


def extract_subrows(raw):
    """Every sub-row a harness's return dict carries, in the harness's own
    order, as {name, source, status, reason, [in_verdict], values}."""
    out = []
    for key in SUBROW_CONTAINERS:
        for name, d in _items(raw.get(key)):
            if not isinstance(d, dict):
                continue
            out.append(_subrow(name, d, key))
            for nk in SUBROW_NESTED:
                for ename, e in _items(d.get(nk)):
                    if isinstance(e, dict):
                        out.append(_subrow("%s/%s" % (name, ename), e,
                                           "%s.%s" % (key, nk)))
    return out


def _text(raw, key, fallback_key=None):
    v = raw.get(key)
    if v is None and fallback_key is not None:
        v = raw.get(fallback_key)
    return None if v is None else one_line(fmt(v))


def normalise_row(raw, name):
    """The generator-independent part of a REPORT.json row, from whatever
    `run()` returned.  `name` is the harness file's stem, which is the row's
    name whatever the harness calls itself (a disagreement is noted)."""
    notes = []
    if not isinstance(raw, dict):
        return dict(name=name, status=ERROR,
                    reason="run() returned %s, not a dict" % type(raw).__name__,
                    generator=None, reference=None, tolerance=None,
                    config={}, subrows=[], notes=notes)
    raw_status = raw.get("status")
    status = None if raw_status is None else str(to_builtin(raw_status)).strip().lower()
    reason = _text(raw, "reason")
    if status not in HARNESS_STATUSES:
        notes.append("run() returned status %r; expected one of %s"
                     % (raw_status, ", ".join(HARNESS_STATUSES)))
        reason = "harness returned status %r" % (raw_status,)
        status = ERROR
    if raw.get("name") not in (None, name):
        notes.append("run() calls itself %r" % (to_builtin(raw.get("name")),))
    return dict(
        name=name, status=status, reason=reason,
        generator=_text(raw, "generator", "generator_value"),
        reference=_text(raw, "reference", "reference_value"),
        tolerance=_text(raw, "tolerance"),
        config=flatten_scalars(raw.get("config")),
        subrows=extract_subrows(raw),
        notes=notes,
    )
