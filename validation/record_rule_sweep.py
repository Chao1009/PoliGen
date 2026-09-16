#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Measure R6 -- the relaxed rule `check_physics_channels_links.py --records`
applies to the dated run records -- against the alternatives to it.

This is the script behind the tables in
`docs/open_items/run_2026-09-06/phase_D_numbers.md` sec. D2.4 (why a shared
token is four characters) and sec. D3.3 (why the two-token rule replaced the
one-token rule).  It was a scratch file when D2.4 was written and the numbers
were published without it, which is the defect D3.4 records; it lives here now
because a published measurement whose recipe is a deleted file is not a
measurement.

    python3 validation/record_rule_sweep.py tokens     # the D2.4 table
    python3 validation/record_rule_sweep.py drift      # the D3.3 table
    python3 validation/record_rule_sweep.py both

`--tree <dir>` measures another checkout (`git archive <commit> | tar -x -C
<dir>`) with THIS gate: that is how the "pre-fix" columns are taken.  The gate
module is loaded by path, so the tree under measurement supplies the records
and the target files and this file supplies the rule.

WHAT IS MEASURED, and why each number means something.

  PASS/FAIL -- how many citations the rule accepts where they actually point.
  A rule that accepts nothing discriminates perfectly and is useless, so this
  is the cost column of every widening.

  DECOY -- the same citation pointed at a block of the same length in the same
  file, shifted clear of the true one: four decoys per citation, at
  +(L+5), -(L+5), +(L+40) and -(L+40) lines, where L is the block length.  A
  decoy that falls outside the file is dropped, not counted as a refusal.  A
  decoy must satisfy R6 (i) -- at least ONE non-blank line -- to be scored at
  all, exactly as a true block must.  This is the discrimination column: a rule
  that accepts a decoy as readily as the truth is measuring nothing.

  DRIFT -- the true block moved by +-k lines, which is what actually happens to
  a record when someone edits the file above it.  This is the failure mode R6
  exists for, and the one the one-token rule was weakest against.

  BREAKAGE -- of the citations the SHIPPING rule accepts today (the `ok` class,
  every one of which was read by a human during the 2026-09-15 repair), how
  many a candidate rule would refuse.  Every one of those is an honest citation
  that would have to be re-pointed, exempted or annotated, so it is the price
  of the candidate.
"""
import argparse
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
GATE = HERE.parent / "check_physics_channels_links.py"


def load_gate(path: Path):
    spec = importlib.util.spec_from_file_location("gate_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------- candidates
#
# Each candidate is a predicate over the SHARED tokens (the intersection of the
# citing sentence's tokens and the target's, both already filtered by the
# gate's own `relax_tokens`: >= `n` characters, stoplisted, case-folded).
def rule_any(shared, n=1):
    return len(shared) >= n


def rule_two_or_long(shared, long=7):
    return len(shared) >= 2 or any(len(t) >= long for t in shared)


def rule_distinctive(shared):
    """One shared token that is not a plain lower-case English word -- it
    carries a digit, an underscore or a dot, or was written with a capital --
    or two shared tokens of any shape."""
    return len(shared) >= 2 or any(not t.isalpha() for t in shared)


def rule_quoted(shared, sw, cite, a, b):
    """R5's own evidence -- a phrase the citing sentence quotes occurs in the
    cited block -- with the one-token rule as the fallback it already is."""
    return sw.quoted(cite, a, b)


def rule_quoted_or_two(shared, sw, cite, a, b):
    return len(shared) >= 2 or sw.quoted(cite, a, b)


def rule_quoted_or_long(shared, sw, cite, a, b):
    return any(len(t) >= 7 for t in shared) or sw.quoted(cite, a, b)


def rule_quoted_or_two_or_long(shared, sw, cite, a, b):
    return (len(shared) >= 2 or any(len(t) >= 7 for t in shared)
            or sw.quoted(cite, a, b))


def rule_rare(shared, sw, cite, a, b, k=3):
    """One shared token that names WHERE in the file the target is: it occurs
    on at most `k` lines of the whole file.  A token that is on forty lines of
    the file cannot tell the cited line from the thirty-nine others, which is
    exactly the drift the record's citations suffer."""
    return any(sw.line_count(cite[2], t) <= k for t in shared)


def rule_offneighbour(shared, sw, cite, a, b):
    """One shared token that is NOT on either line adjacent to the block: the
    minimum a rule must have to notice a one-line drift at all."""
    lines = sw.source_lines(cite[2])
    near = " ".join(lines[max(0, a - 2):a - 1] + lines[b:b + 1]).lower()
    return any(t not in near for t in shared)


def rule_quoted_or_rare(shared, sw, cite, a, b):
    return sw.quoted(cite, a, b) or rule_rare(shared, sw, cite, a, b)


def rule_quoted_or_rare_or_two(shared, sw, cite, a, b):
    return (sw.quoted(cite, a, b) or rule_rare(shared, sw, cite, a, b)
            or len(shared) >= 2)


CANDIDATES = [
    ("R6 as shipped in D2: >= 1 token of 4+", lambda s: rule_any(s, 1), 4),
    (">= 1 token of 5+", lambda s: rule_any(s, 1), 5),
    (">= 1 token of 6+", lambda s: rule_any(s, 1), 6),
    (">= 1 token of 7+", lambda s: rule_any(s, 1), 7),
    (">= 2 tokens of 4+", lambda s: rule_any(s, 2), 4),
    (">= 2 tokens, or one of 7+", lambda s: rule_two_or_long(s, 7), 4),
    (">= 2 tokens, or one of 9+", lambda s: rule_two_or_long(s, 9), 4),
    (">= 2 tokens, or one non-word", rule_distinctive, 4),
    (">= 3 tokens of 4+", lambda s: rule_any(s, 3), 4),
    ("a QUOTED phrase in the target (R5's rule)", rule_quoted, 4),
    ("quoted phrase, or >= 2 tokens", rule_quoted_or_two, 4),
    ("quoted phrase, or one token of 7+", rule_quoted_or_long, 4),
    ("quoted phrase, or >= 2 tokens, or one of 7+", rule_quoted_or_two_or_long, 4),
    ("one token on <= 3 lines of the file", rule_rare, 4),
    ("one token on <= 10 lines of the file",
     lambda sh, sw, c, a, b: rule_rare(sh, sw, c, a, b, 10), 4),
    ("one token off both adjacent lines", rule_offneighbour, 4),
    ("quoted phrase, or a token on <= 3 lines", rule_quoted_or_rare, 4),
    ("quoted, or <= 3 lines, or >= 2 tokens", rule_quoted_or_rare_or_two, 4),
    ("off-adjacent AND (>= 2 tokens or non-word)",
     lambda sh, sw, c, a, b: rule_offneighbour(sh, sw, c, a, b)
     and rule_distinctive(sh), 4),
    ("quoted, or off-adjacent AND >= 2 tokens",
     lambda sh, sw, c, a, b: sw.quoted(c, a, b)
     or (rule_offneighbour(sh, sw, c, a, b) and len(sh) >= 2), 4),
    ("quoted, or off-adjacent AND (2 or non-word)",
     lambda sh, sw, c, a, b: sw.quoted(c, a, b)
     or (rule_offneighbour(sh, sw, c, a, b) and rule_distinctive(sh)), 4),
    ("quoted, or off-adjacent AND one of 6+",
     lambda sh, sw, c, a, b: sw.quoted(c, a, b)
     or (rule_offneighbour(sh, sw, c, a, b)
         and any(len(t) >= 6 for t in sh)), 4),
    ("one token of 5+ off both adjacent lines",
     lambda sh, sw, c, a, b: sw.off_adjacent(c, a, b, sh, 5), 4),
    ("one token of 6+ off both adjacent lines",
     lambda sh, sw, c, a, b: sw.off_adjacent(c, a, b, sh, 6), 4),
    ("one token of 7+ off both adjacent lines",
     lambda sh, sw, c, a, b: sw.off_adjacent(c, a, b, sh, 7), 4),
    ("quoted, or one token of 6+ off adjacent",
     lambda sh, sw, c, a, b: sw.quoted(c, a, b)
     or sw.off_adjacent(c, a, b, sh, 6), 4),
    ("quoted, or one of 5+ off adjacent, or >= 2 off",
     lambda sh, sw, c, a, b: sw.quoted(c, a, b)
     or sw.off_adjacent(c, a, b, sh, 5)
     or len(sw.off_adjacent_set(c, a, b, sh, 4)) >= 2, 4),
    ("one token of 6+ or a NUMBER of 4+, off adjacent",
     lambda sh, sw, c, a, b: bool(sw.off_adjacent_set(c, a, b, sh, 6)
                                  or {t for t in sw.off_adjacent_set(c, a, b, sh, 4)
                                      if any(ch.isdigit() for ch in t)}), 4),
    ("one token of 6+ or a NUMBER of 5+, off adjacent",
     lambda sh, sw, c, a, b: bool(sw.off_adjacent_set(c, a, b, sh, 6)
                                  or {t for t in sw.off_adjacent_set(c, a, b, sh, 5)
                                      if any(ch.isdigit() for ch in t)}), 4),
    ("off-adjacent AND one of 6+",
     lambda sh, sw, c, a, b: rule_offneighbour(sh, sw, c, a, b)
     and any(len(t) >= 6 for t in sh), 4),
]


class Sweep:
    """Every R6 citation of a tree, with the sentence and target text the rule
    reads, collected once so that a candidate is a pure function over them."""

    def __init__(self, gate, stop=True):
        self.g = gate
        self.stop = stop
        self.cites = []          # Cite records, below
        self.lines = {}
        self.srcs = {}
        self._counts = {}
        self._collect()

    def source(self, rel):
        """The gate's own `Source` for a target file -- what the R5 quoted
        -phrase branches read (`declares`, `in_code`)."""
        if rel not in self.srcs:
            f = self.g.ROOT / rel
            self.srcs[rel] = self.g.Source(f, rel) if f.is_file() else None
        return self.srcs[rel]

    def source_lines(self, rel):
        if rel not in self.lines:
            f = self.g.ROOT / rel
            self.lines[rel] = (f.read_text(errors="replace").splitlines()
                               if f.is_file() else None)
        return self.lines[rel]

    def tokens(self, s, n):
        g = self.g
        out = set()
        for m in g.RELAX_TOK.finditer(s):
            t = m.group(0).lower()
            parts = [t] + (t.split(".") if "." in t else [])
            for p in parts:
                if len(p) >= n and (not self.stop or p not in g.RELAX_STOP):
                    out.add(p)
        return out

    def _collect(self):
        g = self.g
        for doc in g.record_files():
            rel = str(doc.relative_to(g.ROOT))
            text = doc.read_text()
            marks = g.marked_sections(text, g.ASOF_SECTION)
            readm = g.marked_sections(text, g.READ_SECTION)
            cites = [(m, int(m.group(2)), int(m.group(2)))
                     for m in g.RECORD_REF.finditer(text)]
            cites += [(m, int(m.group(2)), int(m.group(3)))
                      for m in g.RANGE.finditer(text)]
            cites.sort(key=lambda c: c[0].start())
            for m, a, b in cites:
                path = m.group(1)
                if path is None:
                    path, online = g.record_inherited_path(text, m.start())
                    if path is None and not online:
                        path = g.column_path(text, m.start())
                if path is None:
                    continue
                lines = self.source_lines(path)
                if lines is None or not 1 <= a <= b <= len(lines):
                    continue
                covered = [k for k in
                           [(lo, hi) for lo, hi, _ in marks + readm]
                           if k[0] <= m.start() < k[1]]
                historical = bool(g.ASOF_INLINE.match(text, m.end())) or any(
                    (lo, hi, n) in marks for lo, hi in covered
                    for n in [next((x[2] for x in marks + readm
                                    if (x[0], x[1]) == (lo, hi)), None)])
                if historical:
                    klass = "historical"
                elif not any(x.strip() for x in lines[a - 1:b]):
                    klass = "blank"
                elif path.startswith(g.RECORDS_REL + "/"):
                    klass = "record target"
                elif g.READ_INLINE.match(text, m.end()) or covered:
                    klass = "read"
                else:
                    klass = "live"
                self.cites.append(
                    (rel, text.count("\n", 0, m.start()) + 1, path, a, b,
                     g.relax_sentence(text, m.start(), m.end()), klass,
                     text, m.start(), m.end()))

    # ---------------------------------------------------------------- scoring
    def live(self, only="all"):
        out = [c for c in self.cites if c[6] == "live"]
        if only == "points":
            return [c for c in out if c[3] == c[4]]
        if only == "ranges":
            return [c for c in out if c[3] != c[4]]
        return out

    def line_count(self, rel, token):
        """How many lines of a file carry a token -- the token's power to name
        a place in it."""
        key = (rel, token)
        if key not in self._counts:
            lines = self.source_lines(rel) or []
            self._counts[key] = sum(1 for x in lines if token in x.lower())
        return self._counts[key]

    def off_adjacent_set(self, cite, a, b, shared, n):
        """The shared tokens of at least `n` characters that do NOT occur on
        either line adjacent to the block: the evidence that can tell this
        block from the one a single-line drift would have hit."""
        lines = self.source_lines(cite[2])
        near = " ".join(lines[max(0, a - 2):a - 1] + lines[b:b + 1]).lower()
        return {t for t in shared if len(t) >= n and t not in near}

    def off_adjacent(self, cite, a, b, shared, n):
        return bool(self.off_adjacent_set(cite, a, b, shared, n))

    def quoted(self, cite, a, b):
        """R5's evidence, applied here: does a phrase the citing sentence
        QUOTES occur in the cited block, as text or as a symbol?"""
        _, _, path, _, _, _, _, text, start, end = cite
        src = self.source(path)
        if src is None:
            return False
        for q in self.g.quoted_phrases(text, start, end):
            if self.g.phrase_in_block(src, a, b, q):
                return True
        return False

    def accepts(self, cite, rule, n, a=None, b=None, decoy=False):
        rel, dl, path, ta, tb, sentence, _, _, _, _ = cite
        a = ta if a is None else a
        b = tb if b is None else b
        lines = self.source_lines(path)
        if not 1 <= a <= b <= len(lines):
            return None                                  # out of bounds: drop
        if not any(x.strip() for x in lines[a - 1:b]):
            # R6 (i).  A cited block that is blank is REFUSED; a DECOY block
            # that is blank is not a block anyone could have cited, so it is
            # dropped from the decoy population rather than scored as a
            # refusal the rule can take credit for.
            return None if decoy else False
        shared = self.tokens(sentence, n) & self.tokens(" ".join(lines[a - 1:b]), n)
        try:
            return bool(rule(shared))
        except TypeError:                # a rule that also reads the citation
            return bool(rule(shared, self, cite, a, b))

    def decoys(self, cite):
        _, _, _, a, b, _, _, _, _, _ = cite
        L = b - a + 1
        for d in (L + 5, -(L + 5), L + 40, -(L + 40)):
            yield a + d, b + d


def table_tokens(sw, args):
    print(f"R6 token-length sweep over {len(sw.live())} live-file citations "
          f"in {sw.g.ROOT}  (stoplist {'on' if sw.stop else 'OFF'})")
    print(f"{'n':>3} {'pass':>6} {'fail':>6} {'decoys accepted':>20} {'true/decoy':>11}")
    for n in args.lengths:
        ok = bad = 0
        dok = dn = 0
        for c in sw.live():
            v = sw.accepts(c, lambda s: rule_any(s, 1), n)
            ok, bad = (ok + 1, bad) if v else (ok, bad + 1)
            for a, b in sw.decoys(c):
                v = sw.accepts(c, lambda s: rule_any(s, 1), n, a, b, decoy=True)
                if v is None:
                    continue
                dn += 1
                dok += 1 if v else 0
        rate = 100.0 * dok / dn if dn else 0.0
        ratio = (100.0 * ok / len(sw.live())) / rate if rate else float("inf")
        print(f"{n:>3} {ok:>6} {bad:>6} {dok:>8} / {dn:<6} = {rate:5.1f} % "
              f"{ratio:>10.2f}")


def table_drift(sw, args):
    live = sw.live(args.only)
    print(f"R6 candidate rules over {len(live)} live-file citations "
          f"({args.only}) in {sw.g.ROOT}")
    print(f"{'rule':<38} {'pass':>5} {'breaks':>7} " +
          " ".join(f"{d:+d}".rjust(7) for d in args.drifts) +
          f" {'decoy':>7}")
    base = [sw.accepts(c, lambda s: rule_any(s, 1), 4) for c in live]
    for label, rule, n in CANDIDATES:
        ok = sum(1 for c in live if sw.accepts(c, rule, n))
        breaks = sum(1 for c, b in zip(live, base)
                     if b and not sw.accepts(c, rule, n))
        cells = []
        for d in args.drifts:
            num = den = 0
            for c in live:
                v = sw.accepts(c, rule, n, c[3] + d, c[4] + d)
                if v is None:
                    continue
                den += 1
                num += 1 if v else 0
            cells.append(f"{100.0 * num / den:6.1f}%" if den else "     -")
        dok = dn = 0
        for c in live:
            for a, b in sw.decoys(c):
                v = sw.accepts(c, rule, n, a, b, decoy=True)
                if v is None:
                    continue
                dn += 1
                dok += 1 if v else 0
        print(f"{label:<38} {ok:>5} {breaks:>7} " + " ".join(x.rjust(7) for x in cells)
              + f" {100.0 * dok / dn:6.1f}%")


# ------------------------------------------------------------------- R5
#
# The same shape of measurement for the STRICT documents' recording rule: how
# many range citations a candidate evidences, and how often it evidences a
# DECOY block (the same citation pointed at a same-length block shifted clear
# of the true one).  This is the table behind D3.5.
def r5_cites(g):
    cache = {}
    def src(p):
        if p not in cache:
            f = g.ROOT / p
            cache[p] = g.Source(f, p) if f.exists() else None
        return cache[p]
    out = []
    docs = [g.DOC] + [g.ROOT / p for p, _ in g.EXTRA_DOCS]
    for doc in docs:
        if not doc.exists():
            continue
        text = doc.read_text()
        for m in g.RANGE.finditer(text):
            path = m.group(1) or g.inherited_path(text, m.start())
            if path is None:
                continue
            a, b, name = int(m.group(2)), int(m.group(3)), m.group(4)
            s = src(path)
            if s is None or not 1 <= a < b <= len(s):
                continue
            out.append((doc.name, text.count("\n", 0, m.start()) + 1,
                        path, a, b, name, text, m, s))
    return out


def r5_variants(g, **kw):
    """The gate's OWN R5, with its published constants set to `kw`.  Nothing is
    re-implemented here: the table below measures the shipping rule and the rule
    it replaced, by setting `LIST_ITEM_PARAGRAPH`, `EVIDENCE_WORD_LINES` and
    `EVIDENCE_MAX_BLOCK` and calling `range_evidence`."""
    def evidence(cite, a=None, b=None):
        docn, dl, path, ta, tb, name, text, m, src = cite
        a = ta if a is None else a
        b = tb if b is None else b
        if not 1 <= a < b <= len(src):
            return None
        saved = {k: getattr(g, k) for k in kw}
        for k, v in kw.items():
            setattr(g, k, v)
        try:
            return g.range_evidence(text, m, src, a, b, name)[0]
        finally:
            for k, v in saved.items():
                setattr(g, k, v)
    return evidence


R5_VARIANTS = (
    ("as landed in D1 (no minimum evidence)",
     dict(LIST_ITEM_PARAGRAPH=False, EVIDENCE_WORD_LINES=None,
          EVIDENCE_MAX_BLOCK=None)),
    ("+ a list item is its own paragraph",
     dict(LIST_ITEM_PARAGRAPH=True, EVIDENCE_WORD_LINES=None,
          EVIDENCE_MAX_BLOCK=None)),
    ("+ a one-word text match on <= 5 lines of the file",
     dict(LIST_ITEM_PARAGRAPH=True, EVIDENCE_WORD_LINES=5,
          EVIDENCE_MAX_BLOCK=None)),
    ("+ on <= 3 lines instead",
     dict(LIST_ITEM_PARAGRAPH=True, EVIDENCE_WORD_LINES=3,
          EVIDENCE_MAX_BLOCK=None)),
    ("+ no text match for a block over 150 lines (AS SHIPPED)",
     dict(LIST_ITEM_PARAGRAPH=True, EVIDENCE_WORD_LINES=5,
          EVIDENCE_MAX_BLOCK=150)),
    ("+ over 60 lines instead",
     dict(LIST_ITEM_PARAGRAPH=True, EVIDENCE_WORD_LINES=5,
          EVIDENCE_MAX_BLOCK=60)),
)


def decoy_rate(cites, ev, lo=0, hi=10 ** 6):
    dok = dn = 0
    for c in cites:
        L = c[4] - c[3] + 1
        if not lo <= L < hi:
            continue
        for d in (L + 5, -(L + 5), L + 40, -(L + 40)):
            a, b = c[3] + d, c[4] + d
            if not 1 <= a < b <= len(c[8]):
                continue
            if (not c[8].lines[a - 1].strip()
                    and not c[8].lines[b - 1].strip()):
                continue
            dn += 1
            dok += 1 if ev(c, a, b) else 0
    return dok, dn


def table_r5(gate, args):
    cites = r5_cites(gate)
    print(f"R5 candidate rules over {len(cites)} range citations in {gate.ROOT}")
    print(f"{'rule':<54} {'evid':>5} {'breaks':>7} {'decoys':>18}")
    base = None
    for label, kw in R5_VARIANTS:
        ev = r5_variants(gate, **kw)
        verdicts = [ev(c) for c in cites]
        if base is None:
            base = verdicts
        ok = sum(1 for v in verdicts if v)
        breaks = sum(1 for v, b in zip(verdicts, base) if b and not v)
        dok, dn = decoy_rate(cites, ev)
        print(f"{label:<54} {ok:>5} {breaks:>7} {dok:>6} / {dn:<5} "
              f"= {100.0 * dok / dn:4.1f} %")
    ev = r5_variants(gate, LIST_ITEM_PARAGRAPH=True, EVIDENCE_WORD_LINES=None,
                     EVIDENCE_MAX_BLOCK=None)
    print("\nwhat a TEXT match is worth, by block length "
          "(no minimum evidence):")
    for lo, hi, label in ((0, 20, "1-20 lines"), (20, 60, "20-60"),
                          (60, 150, "60-150"), (150, 10 ** 6, "150+")):
        dok, dn = decoy_rate(cites, ev, lo, hi)
        print(f"  {label:<12} decoys accepted {dok:>4} / {dn:<5} = "
              f"{(100.0 * dok / dn) if dn else 0:4.1f} %")


def unnamed_cites(g):
    """Every UNNAMED point citation of the three strictly gated documents --
    rule B's population."""
    cache = {}
    def src(path):
        if path not in cache:
            f = g.ROOT / path
            cache[path] = g.Source(f, path) if f.exists() else None
        return cache[path]
    out = []
    for doc in [g.DOC] + [g.ROOT / p for p, _ in g.EXTRA_DOCS]:
        if not doc.exists():
            continue
        text = doc.read_text()
        for m in g.REF.finditer(text):
            path, line, name = m.group(1), int(m.group(2)), m.group(3)
            if name and g.CITATION.match(name):
                name = None
            if name:
                continue
            s = src(path)
            if s is None or not 1 <= line <= len(s):
                continue
            out.append((doc.name, text.count("\n", 0, m.start()) + 1,
                        path, line, text, m, s))
    return out


def table_unnamed(g, args):
    """Rule B2: what a quoted-phrase test and a shared-token test cost and buy
    over the unnamed point citations of the strictly gated documents."""
    cites = unnamed_cites(g)
    def quoted(c, line):
        return g.range_evidence(c[4], c[5], c[6], line, line, None)[0] is not None
    def token(c, line):
        sent = g.relax_sentence(c[4], c[5].start(), c[5].end())
        return bool(g.relax_tokens(sent) & g.relax_tokens(c[6].lines[line - 1]))
    rules = (("non-blank only (rule B before 2026-09-15)", lambda c, l: True),
             ("a QUOTED phrase on the line (AS SHIPPED)", quoted),
             ("a shared token of 4+ on the line", token),
             ("either", lambda c, l: quoted(c, l) or token(c, l)))
    print(f"rule B over {len(cites)} unnamed point citations in {g.ROOT}")
    print(f"{'rule':<42} {'flags':>6} " +
          " ".join(f"{d:+d}".rjust(7) for d in args.drifts))
    for label, rule in rules:
        flags = 0
        for c in cites:
            if not (c[6].substantive(c[3]) and rule(c, c[3])):
                flags += 1
        cells = []
        for d in args.drifts:
            num = den = 0
            for c in cites:
                t = c[3] + d
                if not 1 <= t <= len(c[6]):
                    continue
                den += 1
                if c[6].substantive(t) and rule(c, t):
                    num += 1
            cells.append(f"{100.0 * num / den:6.1f}%" if den else "     -")
        print(f"{label:<42} {flags:>6} " + " ".join(x.rjust(7) for x in cells))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("what",
                   choices=["tokens", "drift", "r5", "unnamed", "both"],
                   default="both", nargs="?")
    p.add_argument("--tree", default=None,
                   help="a checkout to measure (default: this repository)")
    p.add_argument("--no-stoplist", action="store_true",
                   help="measure with RELAX_STOP emptied")
    p.add_argument("--lengths", type=int, nargs="+", default=[3, 4, 5])
    p.add_argument("--drifts", type=int, nargs="+", default=[1, -1, 2, -2, 5, -5])
    p.add_argument("--only", choices=["all", "points", "ranges"], default="all",
                   help="a POINT citation names one line and a RANGE a block; "
                        "a one-line drift is a different event for the two")
    args = p.parse_args(argv)
    gate = load_gate(Path(args.tree).resolve() / "validation" /
                     GATE.name if args.tree else GATE)
    if args.tree:                       # measure THAT tree's records and files
        gate.ROOT = Path(args.tree).resolve()
    if args.what == "r5":
        table_r5(gate, args)
        return 0
    if args.what == "unnamed":
        table_unnamed(gate, args)
        return 0
    sw = Sweep(gate, stop=not args.no_stoplist)
    if args.what in ("tokens", "both"):
        table_tokens(sw, args)
    if args.what in ("drift", "both"):
        table_drift(sw, args)
    if args.what in ("r5", "both"):
        table_r5(gate, args)
    if args.what in ("unnamed", "both"):
        table_unnamed(gate, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
