#!/usr/bin/env python3
"""Gate for docs/PHYSICS_CHANNELS.md: every line citation must still point at
what the row says it points at.  Run from the repository root; exit status 1
lists every broken citation.

FOUR CITATION SHAPES ARE CHECKED, and each has its own rule.

A. `` `path:line` `name` `` -- a NAMED point citation (1075 of them).  The
   name must be ADJACENT: one space between the two backticked spans, which is
   how all 1075 are written.
   STRICT (the default).  Three conditions, the first two on the referenced
   line ITSELF -- there is no window:

     S1  the name occurs on that exact line as CODE: not inside a `//`, `///`
         or `/* */` comment, not inside a `#` comment, and not inside a string
         that merely mentions it.  A comment that names the symbol is
         documentation ABOUT it, not the place it lives, and a citation that
         lands on one reads as correct while pointing one or two lines off.

     S2  if that file declares or defines the name ANYWHERE, the referenced
         line must be one of those declaration lines -- a call is not a
         declaration.

     S3  if it declares the name NOWHERE, the citation is a use site, which is
         legitimate (`src/core/xsec.cpp:88` `g1_nucleus` is the one place a
         kernel calls it) -- but S1 then accepts EVERY line of that file that
         carries the name as code.  42 citations are of that shape and 29 of
         them accepted more than one line: 14 lines for
         `python/lipolgen/__init__.py:249` `isotope`, 11 for
         `src/core/pipeline.cpp:531` `B1Model`, then 9, 6, 6, 6.  So a use site
         with more than one accepting line must be PINNED: the cited line is
         fingerprinted in `validation/physics_channels_ranges.json` exactly as
         a range's block is (R4 below), and a drift that lands the citation on
         a different use of the same name is reported -- MOVED, and `--fix`
         re-points it, or EDITED, and a human re-reads the row.  A use site
         with exactly one accepting line needs no entry: the rule already pins
         it, since any drift moves that one line out from under the citation.
         The 13 lone ones carry no sidecar cost.

   S2 deliberately bites only in the declaring file.  A use site cited in the
   file that also defines the symbol is the ambiguous case, and those are named
   in ALLOW below with the reason they are use sites.

   Recognised declaration forms: `#define`; `class`/`struct`/`union`/
   `namespace`/`enum [class] N`; `using N =`; `typedef ... N;`; an enumerator
   inside an enum body; a constructor/destructor; a function or variable
   declaration or definition (including a qualified out-of-line definition
   `A::f(`), a function parameter, an all-caps registration macro whose first
   argument is the name (`LIPOLGEN_CHANNEL(LI6_ALPHA_TAG, ...)`), a Python
   `def`/`class`/module-level assignment, and a name that exists only as a
   string -- a CLI flag, or a binding/metadata key in first-argument or
   subscript position (`m.def("f", ...)`, `meta["k"] = ...`).  A name in a BASE
   CLASS LIST is NOT a declaration of that name: `class ToyF2 : public UnpolSF`
   declares `ToyF2` only (see RESIDUALS 3).

B. `` `path:line` `` with NO name (19 of them).  Nothing names the target, so
   the only thing that can be asserted is that the line is a line somebody
   could have meant: in bounds, NOT BLANK, and not a bare delimiter (`}`, `};`,
   `{`, `)`, `,`, `*/` ... -- a line with no alphanumeric character on it).
   That is the minimum, and it is the whole of it; see RESIDUALS 1.

C. `` `path:first-last` `` -- a RANGE, a block citation: 97 citations of 89
   distinct blocks, 5 of them written `` `:first-last` `` and inheriting the
   path of the citation before them on the same document line.  Checked for:

     R1  in bounds and well ordered: 1 <= first < last <= len(file);
     R2  non-blank first and last line -- a block whose edge has slid onto the
         blank line above or below it is the exact shape this gate exists to
         catch;
     R3  a name, when one is given in the adjacent `` `path:first-last` `name` ``
         form, must be DECLARED inside the range (rule A's declaration test,
         applied to every line of the block).  0 of the 97 currently supply
         one: the document's ranges are block citations whose subject is named
         in the PROSE BEFORE them, and attaching that prose name was measured
         on 2026-09-04 to be wrong 9 times in 11 by backward proximity and 3
         times in 4 by forward proximity, so neither is done;
     R4  THE BLOCK IS STILL THE BLOCK.  R1-R3 cannot see a range whose file
         grew above it: 20 of the 89 blocks had gone stale that way by
         2026-09-04 and only 4 of the 20 had a blank edge to give them away.
         So the content of every cited block is fingerprinted in
         `validation/physics_channels_ranges.json` (sha256 of the block, plus
         its first and last line in clear so the sidecar's own diff is
         readable).  A mismatch is reported one of two ways:
           * the recorded block is found ELSEWHERE in the file -- it MOVED;
             the message names the new range and `--fix` rewrites it;
           * the recorded block is found NOWHERE -- it was EDITED in place;
             re-read the row, confirm the citation still says something true,
             then re-record with `--record-ranges`.
         An entry with no citation left in the document is reported too: a dead
         fingerprint is dead weight that hides the next drift.

   R4 has a running cost, and it is the point: editing a block that
   PHYSICS_CHANNELS.md cites breaks this gate until a human confirms the row
   and re-records.  `--record-ranges` prints every entry it adds, changes or
   drops.  The same sidecar holds the S3 use-site pins and the D fingerprints
   below: 123 entries in all -- 89 blocks, 27 use-site lines (29 citations, two
   pairs of which cite one line) and 7 upstream lines.

D. `` `File.cc:...` `` -- an EXTERNAL citation, a PYTHIA upstream source
   pointer that resolves outside this repository: 5 point citations naming 6
   lines (`BeamRemnants.cc:662,935` names two) and 1 range.  Every C++ file in
   this tree is `.cpp`/`.hpp` and every in-repo citation carries a directory,
   so a bare `.cc` name is always upstream.  They are resolved against the
   dependency tree `env.sh` sets up -- `$LIPOLGEN_PYTHIA_SRC` if it is set,
   otherwise `$LIPOLGEN_DEPS` or `$LIPOLGEN_DEPS_PREFIX` walked back to
   `<prefix>/../src/pythia8*/src` -- and checked for bounds, a
   non-blank edge, and (strict) the same sha256 fingerprint as a range, so a
   PYTHIA upgrade that moves the cited line says so instead of leaving the row
   quietly wrong.  When no PYTHIA source tree is on the machine NOTHING is
   asserted, and each citation is listed as `skipped` by name with the exit
   status unchanged: an unchecked citation you cannot see is a hole.  For the
   same reason `--record-ranges` on such a machine carries the recorded
   external entries across instead of dropping them.  `--fix` does not rewrite
   an external citation: the rows quote the upstream text, so a human re-reads
   it.

MODES.  `--loose` restores the historical behaviour for named point citations
(any textual occurrence within +-2 lines satisfies the citation) and skips R4,
the S3 pins and the D fingerprints.  It is kept for editing a file whose line
numbers are in flux; it is NOT the gate.  The +-2 window is what let 29 anchors
drift onto comments, blank lines and neighbouring declarations while the gate
still reported "0 broken" (run 2026-09-04, phase C).  `--fix` rewrites a broken
citation: a named point onto the nearest DECLARATION of the name (refusing an
exact tie between two declarations, and never touching an ALLOW entry), a range
or a pinned use site onto its new home when the fingerprint locates it, and in
`--loose` mode a named point onto the nearest textual match.
`--record-ranges` (re)writes the sidecar.

RESIDUALS -- what this gate still cannot see, measured 2026-09-04.

  A DRIFT'S SIGN, stated once because these figures are meaningless without it.
  A "+k drift" is k lines LOST above the citation, so the document's number now
  points at what used to be k lines LATER: the citation survives if line + k
  was an accepting line.  A "-k drift" is k lines INSERTED above it.  Both are
  measured end to end -- insert or delete k lines at the top of all 77 cited
  files, run the gate, count the citations it does NOT report.

  1. An UNNAMED point citation (rule B) is pinned only to "a non-blank,
     non-trivial line".  Of the 19, 14 survive a +1 drift and 12 a -1 (10 both);
     at 2 lines it is 17 and 14.  The neighbouring line carries text too.
     Naming the symbol in the document is the only real fix, and is preferred
     whenever the target has a name.

  2. A NAMED point citation is pinned to the set of lines that declare its
     name, which is usually but not always one line.  120 of the 1075 accept
     more than one line (overloads, and a name declared twice in one file);
     2 survive a +1 drift and 34 a +2 -- almost all of them `class X {` on line
     N with `virtual ~X() = default;` on N+2 -- and 0 survive a drift in the
     other direction.  Those 34 are invisible to this gate by construction.
     All 120 and all 34 are declaration-pinned citations: since S3 the 42
     use-site citations accept exactly one line each, and 0 of them survive any
     drift of 1 or 2 lines either way.
     (The base-clause guard of RESIDUALS 3 changes NEITHER figure -- see there.
     A previous note here claimed it took them from 156 and 36; that is
     withdrawn.)

  3. The declaration test is textual, not a parse.  It reads `int f(int x)` as
     declaring both `f` and `x`, which is deliberate (parameters are cited).
     Until 2026-09-04 it also read a BASE CLASS as a declaration -- `class
     ToyF2 : public UnpolSF {` counted as declaring `UnpolSF`, so
     `sf.hpp:359` `TensorSF` would have been accepted on any of four derived
     classes' headers (`:421`, `:432`, `:444`, `:494`).  No citation sat on
     one, but `--fix` could have landed there; the base-clause guard in
     `cxx_declares` closes it, and `TensorSF` now declares at :359 and :361
     only.  MEASURED REACH over the cited set: toggling the guard changes the
     acceptance set of exactly 6 citations -- `sf.hpp:77` `UnpolSF`, `:161`
     `PolSF`, `:359` `TensorSF`, `triton_sf.hpp:226` `TritonSpectralFunction`,
     `fsi.hpp:214` `FsiWeight`, `rc.hpp:356` `Spin1ElasticFF` -- every one of
     which already accepted more than one line and already survived a +2 drift
     WITH the guard on, because each is a `class X {` whose destructor sits on
     X+2.  So the guard closes a `--fix` hazard and changes RESIDUALS 2's two
     counts by 0.  What remains: a name inside a template argument list or a
     trailing return type can still read as a parameter declaration.

  4. A range or use site or upstream line with a valid fingerprint that a
     maintainer re-records without re-reading the row is blessed.  The gate can
     force the question; it cannot answer it.

  5. Citations to files outside the recognised path set are not checked at all.
     The set is include/ src/ python/ data/ validation/ docs/ tests/ examples/
     plus README.md, and -- since 2026-09-04, by rule D -- bare `File.cc`
     PYTHIA upstream pointers.  `pyproject.toml` and `CMakeLists.txt` are cited
     by line range in the packaging row and are NOT in the set; that row says
     so.
"""
import hashlib
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "PHYSICS_CHANNELS.md"
RANGES_REL = "validation/physics_channels_ranges.json"
PATH = r"(?:(?:include|src|python|data|validation|docs|tests|examples)/[^`:\s]+|README\.md)"
# A point citation `path:line`, optionally followed by `name`.  The name must be
# ADJACENT -- one space, nothing else.  All 1075 named citations in the document
# are written that way; the twelve-character window this replaces bought nothing
# and attached one citation's neighbouring prose symbol to it.  And the name is
# matched by LOOKAHEAD, never consumed: the consuming group ate the citation
# that followed it, which hid ten citations (one point, nine ranges) from the
# gate entirely until 2026-09-04.
REF = re.compile(r"`(" + PATH + r"):(\d+)`(?:(?= `([^`\n]+)`))?")
# A range citation `path:first-last`, or a bare `:first-last` continuing the
# citation before it on the same document line.  A range's name must be
# adjacent -- one space, nothing else; the prose around a range names its
# subject too often for a wider window to mean anything (see R3).
RANGE = re.compile(r"`(" + PATH + r")?:(\d+)-(\d+)`(?:(?= `([^`\n]+)`))?")
# `docs/USAGE.md:11-50`, `:343-347`: a citation, never a symbol name.
CITATION = re.compile(r"^(?:" + PATH + r")?:\d+(?:-\d+)?$")
# An EXTERNAL citation: a bare `File.cc:...` with no directory.  Every C++ file
# in this repository is `.cpp`/`.hpp` and every in-repo citation carries a
# directory, so a bare `.cc` name is always a PYTHIA upstream source pointer.
# Six of them (`BeamRemnants.cc:662,935` names two lines) resolve outside the
# repository, so until 2026-09-04 they matched no rule and were checked by
# nothing.  The spec is one or more comma-separated lines or ranges.
EXT = re.compile(r"`([A-Za-z][A-Za-z0-9_]*\.cc):(\d+(?:[-,]\d+)*)`")
PYTHIA_KEY = "pythia8:"
WIN = 2

# Point citations that cannot land on a declaration.  An exemption is keyed on
# the EXACT reference -- (path, line, name) -- so it covers that one citation
# and no other, and it carries a PIN: a distinctive substring of the line it
# names.  The pin is what makes the exemption drift-proof.  `mentions()` on the
# name alone was not: it matched comments and strings anywhere on whatever line
# the number happened to land on, and being keyed on (path, name) it silently
# exempted the OTHER citations of the same name in the same file as well.
# Each value is (line-pin, reason).
#
# ONE OF THE SIX WAS NOT PER-CITATION IN EFFECT, and the leak was not in this
# table.  Five of these names are DECLARED in the file that holds their
# exemption, so S2 governs every other citation of them there and says so by
# name ("an exemption covers one citation, not a name").  `coherent` is not:
# `src/core/pipeline.cpp` declares no such symbol, so before S3 the use-site
# carve-out accepted `:996` and `:1009` -- both `cfg_.coherent` -- with nothing
# printed at all.  The exemption's effective reach over that file was 3 lines,
# not 1, and 2 of them were invisible.  S3 closes it: those two lines are an
# ambiguous use site and now need a fingerprint of their own.  Measured after
# the fix, each of the six covers exactly one citation and nothing else in its
# file passes unpinned.
ALLOW = {
    ("include/lipolgen/xsec.hpp", 220, "g1_rank3"): (
        "no rank-3 (octupole) slot",
        "no such symbol exists in the tree, deliberately: the row cites the "
        "comment that states why there is no rank-3 slot"),
    ("src/core/pipeline.cpp", 773, "coherent"): (
        "on the coherent channel the tensor signal is in the recoil",
        "the cited line is inside PipelineConfig::validate()'s refusal "
        "message; `coherent` there is the English word in that message, which "
        "is the text the row quotes"),
    ("src/core/pipeline.cpp", 1692, "optics_lumi_factor"): (
        "cfg_.lumi_pb * optics_lumi_factor()",
        "the row's claim is that luminosity mode multiplies the optics factor "
        "in, so it cites the multiplication, not the accessor (declared at "
        "include/lipolgen/pipeline.hpp:779, cited there too)"),
    ("src/core/rc.cpp", 1032, "is_tagged_channel"): (
        "is_tagged_channel(channel_)",
        "one of three consecutive branch lines of RcModel::exclusion_reason "
        "the row cites together; the neighbours are use sites in files that "
        "do not declare their symbol"),
    ("src/core/breakup.cpp", 332, "proton_fraction"): (
        "proton_fraction(in.x, in.q2, 1, 3)",
        "the triton branch's own proton/neutron draw, which is what the row "
        "describes; the definition is generic over (Z, A)"),
    ("src/pythia/pythia_bridge.cpp", 700, "dis_parton_fraction"): (
        "dis_parton_fraction(q, p_n, &xi_tmp, 1.0, mq)",
        "the re-solve with the chosen quark's mass, cited beside the "
        "definition at src/pythia/pythia_bridge.cpp:170 in the same row"),
}

CTRL = {"if", "for", "while", "switch", "return", "else", "do", "case", "new",
        "delete", "throw", "goto", "sizeof", "static_cast", "const_cast",
        "dynamic_cast", "reinterpret_cast", "co_return", "co_await", "catch",
        "try", "break", "continue"}
BAD_PREFIX = set("()=;.?!+/%|\"'[]")
ANGLE = re.compile(r"<[^<>]*>")
IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
FLAG = re.compile(r"^--?[A-Za-z0-9][\w-]*$")
ALNUM = re.compile(r"[A-Za-z0-9]")
# `class D : public B` / `struct D final : B` -- everything after that colon is
# a base-specifier list, which declares D and mentions B.
CLASS_HEAD = re.compile(r"\b(?:class|struct)\s+[A-Za-z_][A-Za-z0-9_]*\s*(?:final\s*)?:")
LIT = r"\"((?:[^\"\\]|\\.)*)\"|'((?:[^'\\]|\\.)*)'"


def strip_angles(s: str) -> str:
    prev = None
    while prev != s:
        prev, s = s, ANGLE.sub(" ", s)
    return s


def word(base: str) -> re.Pattern:
    return re.compile(r"(?<![A-Za-z0-9_])" + re.escape(base) + r"(?![A-Za-z0-9_])")


def mask_cxx(lines: list[str], keep_strings: bool = False) -> list[str]:
    """Blank out comment and string/char-literal characters, keeping columns.

    `keep_strings` blanks the COMMENTS ONLY and leaves string and character
    literals standing.  That is the line `literal_defines` has to read: a name
    that lives only as a string is introduced by the line that spells it, and
    `mask_cxx`'s ordinary output has already erased it -- but the RAW line,
    which is what this function was called with until 2026-09-05, also carries
    the comments, and a `///` line is DOCUMENTATION and not a declaration.
    That is exactly how `docs/PHYSICS_CHANNELS.md` came to cite
    `include/lipolgen/pythia_bridge.hpp:338` for `n_pom_flavour_fallback`:
    :338 reads `/// Recorded per run since 2026-09-04 as
    meta["n_pom_flavour_fallback"]`, whose bracketed literal satisfied
    `literal_defines`, so the doc comment one line ABOVE the declaration
    counted as a declaration and `--fix` re-anchored onto it -- invisibly,
    because the gate then agreed with itself.  With the comments gone that
    line offers nothing and :340, the declaration, is the only target left.
    """
    out, in_block = [], False
    for ln in lines:
        chars, i, n = list(ln), 0, len(ln)
        while i < n:
            c = chars[i]
            if in_block:
                if c == "*" and i + 1 < n and chars[i + 1] == "/":
                    chars[i] = chars[i + 1] = " "
                    i += 2
                    in_block = False
                else:
                    chars[i] = " "
                    i += 1
                continue
            if c == "/" and i + 1 < n and chars[i + 1] == "/":
                for j in range(i, n):
                    chars[j] = " "
                break
            if keep_strings and c in "\"'":
                q, j = c, i + 1
                while j < n:
                    if chars[j] == "\\":
                        j += 2
                        continue
                    if chars[j] == q:
                        j += 1
                        break
                    j += 1
                i = j
                continue
            if c == "/" and i + 1 < n and chars[i + 1] == "*":
                chars[i] = chars[i + 1] = " "
                i += 2
                in_block = True
                continue
            if c in "\"'":
                q, j = c, i + 1
                chars[i] = " "
                while j < n:
                    if chars[j] == "\\":
                        chars[j] = " "
                        if j + 1 < n:
                            chars[j + 1] = " "
                        j += 2
                        continue
                    if chars[j] == q:
                        chars[j] = " "
                        j += 1
                        break
                    chars[j] = " "
                    j += 1
                i = j
                continue
            i += 1
        out.append("".join(chars))
    return out


def literal_defines(line: str, name: str) -> bool:
    """The site that introduces a name living only as a string: a CLI flag
    anywhere on the line, or a key in first-argument / subscript position."""
    if name.startswith("-"):
        return any((m.group(1) or m.group(2) or "") == name
                   for m in re.finditer(LIT, line))
    return any((m.group(1) or m.group(2) or "") == name
               for m in re.finditer(r"[(\[]\s*(?:" + LIT + r")", line))


def macro_registers(masked_line: str, name: str) -> bool:
    """`REGISTER_THING(NAME, ...)`: the X-macro / registration idiom declares
    through an all-caps wrapper."""
    return bool(re.match(r"^\s*[A-Z][A-Z0-9_]*\s*\(\s*" + re.escape(name)
                         + r"\s*[,)]", masked_line))


def enum_bodies(masked: list[str]) -> tuple[set[int], set[int]]:
    """(lines inside a multi-line `enum { ... }`, one-line `enum E { A, B };`)."""
    inside, depth, body, single = False, 0, set(), set()
    for i, s in enumerate(masked, 1):
        if not inside:
            if re.search(r"\benum\b", s) and "{" in s:
                if "}" in s:
                    single.add(i)
                else:
                    inside, depth = True, s.count("{") - s.count("}")
            continue
        depth += s.count("{") - s.count("}")
        if depth <= 0:
            inside = False
        else:
            body.add(i)
    return body, single


def cxx_tags(masked: list[str]) -> set[str]:
    return set(re.findall(r"\b(?:class|struct)\s+([A-Za-z_][A-Za-z0-9_]*)",
                          " \n".join(masked)))


def cxx_declares(masked, i, base, enums, tags, code):
    """`code` is the line with its COMMENTS stripped and its string literals
    kept (`mask_cxx(..., keep_strings=True)`), never the raw line: a `//` or
    `///` line is documentation, and a name it merely spells inside a string
    is not declared there."""
    s = masked[i - 1]
    body, single = enums
    if literal_defines(code, base):
        return "literal"
    if macro_registers(s, base):
        return "macro-call"
    if re.search(r"^\s*#\s*define\s+" + re.escape(base) + r"(?![A-Za-z0-9_])", s):
        return "macro"
    # `struct A::B {` declares B, not A -- hence the (?!\s*::) guard
    if re.search(r"\b(?:struct|class|union|namespace)\s+(?:alignas\s*\([^)]*\)\s*)?"
                 + re.escape(base) + r"(?![A-Za-z0-9_])(?!\s*::)", s):
        return "tag"
    if re.search(r"\benum\s+(?:class\s+|struct\s+)?" + re.escape(base)
                 + r"(?![A-Za-z0-9_])(?!\s*::)", s):
        return "tag"
    if re.search(r"\busing\s+" + re.escape(base) + r"\s*=", s):
        return "alias"
    if re.search(r"\btypedef\b.*(?<![A-Za-z0-9_])" + re.escape(base) + r"\s*[;\[]", s):
        return "alias"
    if i in body or i in single:
        seg = s.split("{", 1)[1] if i in single else s
        for m in word(base).finditer(seg):
            pre, post = seg[:m.start()], seg[m.end():].lstrip()
            if pre.rstrip().endswith((".", "->", "::")):
                continue
            if post[:1] in ("", ",", "=", "}") and not post.startswith("=="):
                return "enumerator"
    if base in tags and re.match(r"^\s*(?:explicit\s+|constexpr\s+|inline\s+)*~?\s*"
                                 + re.escape(base) + r"\s*\(", s):
        return "ctor"
    for m in word(base).finditer(s):
        pre, post = s[:m.start()], s[m.end():]
        if post.startswith("::"):
            if re.match(r"::\s*~?\s*" + re.escape(base) + r"\s*\(", post):
                return "ctor"
            continue
        if pre.rstrip().endswith((".", "->")):
            continue
        # `class D : public B` names B in a base-specifier list; the line
        # declares D.  Without this, four derived classes each counted as a
        # declaration of `UnpolSF` / `TensorSF`.  The list ends at the opening
        # brace, so a one-line `class D : B { int x; };` still declares `x`.
        head = CLASS_HEAD.search(pre)
        if head and "{" not in pre[head.end():]:
            continue
        p = post.lstrip()
        if p[:1] not in ("(", ";", "=", "{", "[", ",", ")", ":") and p != "":
            continue
        if p.startswith("==") or p.startswith("=>"):
            continue
        qual = re.sub(r"(?:[A-Za-z_][A-Za-z0-9_]*\s*::\s*)+$", "", pre)
        for seg, kind in ((qual, "decl"),
                          (qual.rsplit("(", 1)[-1].rsplit(",", 1)[-1], "param")):
            q = strip_angles(seg).strip()
            if not q or any(c in q for c in BAD_PREFIX) or "->" in q:
                continue
            if any(t in CTRL for t in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", q)):
                continue
            return kind
        if not strip_angles(qual).strip() and qual != pre:
            return "def"
    return None


def py_depth0(lines: list[str]) -> set[int]:
    """1-based lines that begin at bracket depth 0, i.e. start a statement."""
    res, depth = set(), 0
    for i, ln in enumerate(lines, 1):
        if depth == 0:
            res.add(i)
        s = re.sub(r"'''.*?'''|\"\"\".*?\"\"\"", "", ln)
        s = re.sub(LIT, "", s)
        s = re.sub(r"#.*$", "", s)
        depth += s.count("(") + s.count("[") + s.count("{")
        depth -= s.count(")") + s.count("]") + s.count("}")
        depth = max(depth, 0)
    return res


def py_declares(lines, i, base, depth0):
    raw = lines[i - 1]
    if literal_defines(raw, base):
        return "literal"
    if base.startswith("-"):
        return None
    s = re.sub(r"#.*$", "", raw)
    if re.match(r"^\s*(?:async\s+)?def\s+" + re.escape(base) + r"\s*\(", s):
        return "def"
    if re.match(r"^\s*class\s+" + re.escape(base) + r"\b", s):
        return "class"
    if i in depth0 and re.match(r"^\s*" + re.escape(base) + r"\s*(?::[^=]+)?=(?!=)", s):
        return "assign"
    return None


class Source:
    """One referenced file, with the per-language tables the rules need."""

    def __init__(self, path: Path, rel: str):
        self.lines = path.read_text(errors="replace").splitlines()
        if rel.endswith((".hpp", ".cpp", ".h", ".cc")):
            self.kind = "cxx"
            self.masked = mask_cxx(self.lines)
            # ... and the same line with the COMMENTS gone and the string
            # literals kept, which is the only line `literal_defines` may read.
            self.code = mask_cxx(self.lines, keep_strings=True)
            self.enums = enum_bodies(self.masked)
            self.tags = cxx_tags(self.masked)
        elif rel.endswith(".py"):
            self.kind = "py"
            self.depth0 = py_depth0(self.lines)
        else:
            self.kind = "txt"
        self._decls: dict[str, list[int]] = {}
        self._uses: dict[str, list[int]] = {}
        self._windows: dict[int, dict[str, list[int]]] = {}

    def __len__(self):
        return len(self.lines)

    def declares(self, i: int, base: str):
        if not 1 <= i <= len(self.lines):
            return None
        if not (IDENT.match(base) or FLAG.match(base)):
            # not an identifier at all (`pe=0.7`, `Mo-Tsai`, `c.wgt`): the
            # reference is to the text, and the exact line must carry it
            return "text" if base in self.lines[i - 1] else None
        if self.kind == "cxx":
            return cxx_declares(self.masked, i, base, self.enums, self.tags,
                                self.code[i - 1])
        if self.kind == "py":
            return py_declares(self.lines, i, base, self.depth0)
        return "text" if word(base).search(self.lines[i - 1]) else None

    def in_code(self, i: int, base: str) -> bool:
        if not 1 <= i <= len(self.lines):
            return False
        if not (IDENT.match(base) or FLAG.match(base)):
            return base in self.lines[i - 1]
        if literal_defines(self.code[i - 1] if self.kind == "cxx"
                           else self.lines[i - 1], base):
            return True
        if self.kind == "cxx":
            s = self.masked[i - 1]
        elif self.kind == "py":
            s = re.sub(r"#.*$", "", self.lines[i - 1])
        else:
            s = self.lines[i - 1]
        return bool(word(base).search(s))

    def decl_lines(self, base: str) -> list[int]:
        if base not in self._decls:
            self._decls[base] = [i for i in range(1, len(self.lines) + 1)
                                 if self.declares(i, base)]
        return self._decls[base]

    def use_lines(self, base: str) -> list[int]:
        """Every line carrying the name as CODE.  In a file that declares the
        name nowhere this is the whole acceptance set of a use-site citation,
        which is why S3 pins it when it holds more than one line."""
        if base not in self._uses:
            self._uses[base] = [i for i in range(1, len(self.lines) + 1)
                                if self.in_code(i, base)]
        return self._uses[base]

    def substantive(self, i: int) -> bool:
        """A line somebody could have meant: not blank, and not a bare
        delimiter (`}`, `};`, `)`, `*/`, a lone comma)."""
        if not 1 <= i <= len(self.lines):
            return False
        return bool(ALNUM.search(self.lines[i - 1]))

    def block_sha(self, a: int, b: int) -> str:
        return hashlib.sha256("\n".join(self.lines[a - 1:b]).encode()).hexdigest()

    def find_block(self, sha: str, length: int) -> list[int]:
        """1-based starts of every window of `length` lines hashing to `sha`."""
        table = self._windows.get(length)
        if table is None:
            table = {}
            for i in range(0, len(self.lines) - length + 1):
                h = hashlib.sha256("\n".join(self.lines[i:i + length]).encode()).hexdigest()
                table.setdefault(h, []).append(i + 1)
            self._windows[length] = table
        return table.get(sha, [])


def nearest_decl(src: Source, base: str, old: int):
    """Declaration line nearest to `old`; None when there is none or the
    nearest two are equidistant."""
    cands = src.decl_lines(base)
    if not cands:
        return None
    cands = sorted(cands, key=lambda i: abs(i - old))
    if len(cands) > 1 and abs(cands[0] - old) == abs(cands[1] - old):
        return None
    return cands[0]


def nearest_text(src: Source, base: str, old: int):
    pat = word(base)
    cands = [i for i, ln in enumerate(src.lines, 1) if pat.search(ln)]
    if not cands:
        return None
    cands = sorted(cands, key=lambda i: abs(i - old))
    if len(cands) > 1 and abs(cands[0] - old) == abs(cands[1] - old):
        return None
    return cands[0]


PATH_AT = re.compile(r"`(" + PATH + r"):\d")


def inherited_path(text: str, pos: int):
    """A bare `` `:first-last` `` continues the citation before it on the same
    document line; this returns that citation's path."""
    start = text.rfind("\n", 0, pos) + 1
    last = None
    for m in PATH_AT.finditer(text, start, pos):
        last = m.group(1)
    return last


def pythia_source_dir():
    """Where PYTHIA's own `.cc` files live, or None when they are not on this
    machine.  `$LIPOLGEN_PYTHIA_SRC` wins if it is set.  Otherwise the
    dependency prefix is taken from `$LIPOLGEN_DEPS` (what `env.sh` exports) or
    `$LIPOLGEN_DEPS_PREFIX` (the name CMakeLists.txt and the pip instructions
    in README.md / docs/USAGE.md use for the same directory) -- both point at
    `<...>/deps/install`, and the unpacked sources sit beside it at
    `<...>/deps/src/pythia8*/src`.  Returning None is not a failure: the
    citations are reported as SKIPPED, by name, so an unchecked citation is
    still visible."""
    env = os.environ.get("LIPOLGEN_PYTHIA_SRC")
    if env:
        d = Path(env)
        return d if d.is_dir() else None
    for var in ("LIPOLGEN_DEPS", "LIPOLGEN_DEPS_PREFIX"):
        deps = os.environ.get(var)
        if not deps:
            continue
        base = Path(deps)
        for cand in sorted(base.parent.glob("src/pythia8*/src")) + \
                sorted(base.glob("src/pythia8*/src")):
            if cand.is_dir():
                return cand
    return None


def ranges_file() -> Path:
    return ROOT / RANGES_REL


def load_ranges() -> dict:
    f = ranges_file()
    if not f.exists():
        return {}
    return json.loads(f.read_text())


def write_ranges(data: dict) -> None:
    f = ranges_file()
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(dict(sorted(data.items())), indent=1,
                            ensure_ascii=False) + "\n")


def clip(s: str, n: int = 100) -> str:
    s = s.rstrip()
    return s if len(s) <= n else s[:n - 1] + "…"


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    fix = "--fix" in argv
    record = "--record-ranges" in argv
    strict = "--loose" not in argv
    if not DOC.exists():
        print(f"missing {DOC}")
        return 1
    text = DOC.read_text()
    cache: dict[str, Source | None] = {}
    bad: list[str] = []
    fixed: list[str] = []
    allowed: list[str] = []
    skipped: list[str] = []
    n = nr = nx = 0
    recorded = load_ranges()
    fresh: dict[str, dict] = {}
    seen: set[str] = set()
    remap: dict[str, str] = {}
    ext_dir = pythia_source_dir()

    def source(path: str):
        if path not in cache:
            f = ROOT / path
            cache[path] = Source(f, path) if f.exists() else None
        return cache[path]

    def range_repl(m: re.Match) -> str:
        nonlocal nr
        path = m.group(1) or inherited_path(text, m.start())
        a, b, name = int(m.group(2)), int(m.group(3)), m.group(4)
        nr += 1
        shown = f"{m.group(1) or ''}:{a}-{b}"
        if path is None:
            bad.append(f"{shown}  (bare range with no citation before it on "
                       f"the same line to take a path from)")
            return m.group(0)
        src = source(path)
        if src is None:
            bad.append(f"{path}:{a}-{b}  (file missing)")
            return m.group(0)
        if not 1 <= a < b <= len(src):                                    # R1
            bad.append(f"{path}:{a}-{b}  (not a range inside a "
                       f"{len(src)}-line file)")
            return m.group(0)
        if not src.lines[a - 1].strip() or not src.lines[b - 1].strip():  # R2
            which = "first" if not src.lines[a - 1].strip() else "last"
            bad.append(f"{path}:{a}-{b}  the {which} line of the range is blank")
            return m.group(0)
        if name and not CITATION.match(name):                             # R3
            base = re.split(r"[(<]", name)[0].split("::")[-1].strip()
            if IDENT.match(base) or FLAG.match(base):
                if not any(src.declares(i, base) for i in range(a, b + 1)):
                    bad.append(f"{path}:{a}-{b}  `{name}` is not declared "
                               f"anywhere inside the range")
                    return m.group(0)
        key = f"{path}:{a}-{b}"
        seen.add(key)
        if not strict:
            return m.group(0)
        sha = src.block_sha(a, b)
        fresh[key] = {"sha256": sha,
                      "first": clip(src.lines[a - 1]),
                      "last": clip(src.lines[b - 1])}
        if record:
            return m.group(0)
        entry = recorded.get(key)
        if entry is None:                                                 # R4
            bad.append(f"{path}:{a}-{b}  block is not fingerprinted; check the "
                       f"row, then re-record with --record-ranges")
            return m.group(0)
        if entry.get("sha256") == sha:
            return m.group(0)
        hits = src.find_block(entry.get("sha256", ""), b - a + 1)
        if len(hits) == 1:
            na, nb = hits[0], hits[0] + (b - a)
            if fix:
                fixed.append(f"{path}:{a}-{b} -> :{na}-{nb}  (block moved)")
                remap[key] = f"{path}:{na}-{nb}"
                return m.group(0).replace(f":{a}-{b}`", f":{na}-{nb}`", 1)
            bad.append(f"{path}:{a}-{b}  the cited block MOVED to "
                       f":{na}-{nb} (\"{entry.get('first', '')}\" …)")
            return m.group(0)
        if len(hits) > 1:
            bad.append(f"{path}:{a}-{b}  the cited block now occurs at "
                       + ", ".join(f":{h}-{h + b - a}" for h in hits))
            return m.group(0)
        bad.append(f"{path}:{a}-{b}  the cited block was EDITED in place "
                   f"(recorded \"{entry.get('first', '')}\" … "
                   f"\"{entry.get('last', '')}\"); confirm the row still says "
                   f"something true, then re-record with --record-ranges")
        return m.group(0)

    def repl(m: re.Match) -> str:
        nonlocal n
        path, line, name = m.group(1), int(m.group(2)), m.group(3)
        if name and CITATION.match(name):
            name = None                # the next citation, not this one's name
        n += 1
        src = source(path)
        if src is None:
            bad.append(f"{path}:{line}  (file missing)")
            return m.group(0)
        if not 1 <= line <= len(src):
            bad.append(f"{path}:{line}  (only {len(src)} lines)")
            return m.group(0)
        if not name:
            # Rule B: nothing names the target, so the only assertion left is
            # that the line is one somebody could have meant.
            if not src.substantive(line):
                what = "blank" if not src.lines[line - 1].strip() else \
                       f"only punctuation ({src.lines[line - 1].strip()!r})"
                bad.append(f"{path}:{line}  unnamed reference to a line that "
                           f"is {what}")
            return m.group(0)
        base = re.split(r"[(<]", name)[0].split("::")[-1].strip()
        if not base:
            return m.group(0)

        if not strict:
            lo, hi = max(1, line - WIN), min(len(src), line + WIN)
            if any(base in src.lines[i - 1] for i in range(lo, hi + 1)):
                return m.group(0)
            if fix:
                new = nearest_text(src, base, line)
                if new is not None:
                    fixed.append(f"{path}:{line} -> :{new}  `{name}`")
                    return m.group(0).replace(f"{path}:{line}`", f"{path}:{new}`", 1)
            bad.append(f"{path}:{line}  `{name}` not within +-{WIN} lines")
            return m.group(0)

        if src.declares(line, base):
            return m.group(0)
        decls = src.decl_lines(base)
        if not decls and src.in_code(line, base):
            # S3.  A use site in a file that declares the name NOWHERE is a
            # legitimate citation (`src/core/xsec.cpp:88` `g1_nucleus`), but
            # the rule that admits it accepts EVERY line of that file carrying
            # the name as code -- 14 lines for `python/lipolgen/__init__.py`
            # `isotope`, 11 for `src/core/pipeline.cpp` `B1Model`.  So when
            # there is more than one such line the citation is pinned by a
            # fingerprint of the cited LINE, exactly as a range is pinned by a
            # fingerprint of its block; when there is only one, the rule
            # already pins it and no sidecar entry is needed.
            uses = src.use_lines(base)
            if len(uses) < 2 or not strict:
                return m.group(0)
            key = f"{path}:{line}"
            seen.add(key)
            sha = src.block_sha(line, line)
            fresh[key] = {"sha256": sha, "line": clip(src.lines[line - 1]),
                          "name": base}
            if record:
                return m.group(0)
            entry = recorded.get(key)
            if entry is None:
                bad.append(f"{path}:{line}  `{name}` is a use site in a file "
                           f"that declares it nowhere and names it as code on "
                           f"{len(uses)} lines, so the line must be "
                           f"fingerprinted; check the row, then re-record with "
                           f"--record-ranges")
                return m.group(0)
            if entry.get("sha256") == sha:
                return m.group(0)
            hits = src.find_block(entry.get("sha256", ""), 1)
            was = entry.get("line", "")
            if len(hits) == 1:
                if fix:
                    fixed.append(f"{path}:{line} -> :{hits[0]}  `{name}` "
                                 f"(use site moved)")
                    remap[key] = f"{path}:{hits[0]}"
                    return m.group(0).replace(f"{path}:{line}`",
                                              f"{path}:{hits[0]}`", 1)
                bad.append(f"{path}:{line}  `{name}` -- the cited use site "
                           f"MOVED to :{hits[0]} (\"{was}\")")
                return m.group(0)
            if len(hits) > 1:
                bad.append(f"{path}:{line}  `{name}` -- the cited use site now "
                           f"occurs at " + ", ".join(f":{h}" for h in hits))
                return m.group(0)
            bad.append(f"{path}:{line}  `{name}` -- the cited use site was "
                       f"EDITED in place (recorded \"{was}\"); confirm the row "
                       f"still says something true, then re-record with "
                       f"--record-ranges")
            return m.group(0)

        # The exemption is reached only by a citation the rule already
        # rejected, and it is keyed on (path, LINE, name): it covers exactly
        # this citation, leaves every other citation of the same name in the
        # same file fully checked, and holds only while the line it names
        # still carries the pinned text.
        key = (path, line, base) if (path, line, base) in ALLOW else \
              (path, line, name) if (path, line, name) in ALLOW else None
        if key is not None:
            pin, why = ALLOW[key]
            if pin in src.lines[line - 1]:
                allowed.append(f"{path}:{line}  `{name}` -- {why}")
                return m.group(0)
            bad.append(f"{path}:{line}  `{name}` is allow-listed on this line "
                       f"but the line no longer contains its pin {pin!r}; the "
                       f"block has moved, or the code changed")
            return m.group(0)
        if any(p == path and nm in (base, name) for p, _, nm in ALLOW):
            bad.append(f"{path}:{line}  `{name}` is allow-listed in this file "
                       f"at another line only; an exemption covers one "
                       f"citation, not a name")
            return m.group(0)
        if fix:
            new = nearest_decl(src, base, line)
            if new is not None:
                fixed.append(f"{path}:{line} -> :{new}  `{name}`")
                return m.group(0).replace(f"{path}:{line}`", f"{path}:{new}`", 1)
        why = (f"`{name}` is declared at {', '.join(':%d' % d for d in decls)} "
               f"but the reference points at a "
               + ("use site" if src.in_code(line, base) else
                  "line that does not carry it as code")) if decls else \
              f"`{name}` is not on that line and is declared nowhere in the file"
        bad.append(f"{path}:{line}  {why}")
        return m.group(0)

    def ext_check(m: re.Match) -> None:
        """A PYTHIA upstream pointer.  Nothing in this repository can pin it,
        so it is checked against the dependency tree `env.sh` points at: in
        bounds, not blank, and -- in strict mode -- still hashing to what was
        recorded, which is what makes a PYTHIA upgrade that moves the line say
        so instead of leaving the row quietly wrong.  `--fix` does NOT rewrite
        these: the row quotes the upstream text, so a human re-reads it."""
        nonlocal nx
        fname, spec = m.group(1), m.group(2)
        nx += 1
        if ext_dir is None:
            skipped.append(f"{fname}:{spec}  external (PYTHIA upstream): no "
                           f"PYTHIA source tree on this machine -- set "
                           f"$LIPOLGEN_PYTHIA_SRC, or source env.sh so "
                           f"$LIPOLGEN_DEPS names the dependency prefix")
            return
        f = ext_dir / fname
        if not f.exists():
            bad.append(f"{fname}:{spec}  (external: no such file under "
                       f"{ext_dir})")
            return
        src = Source(f, fname)
        for part in spec.split(","):
            if "-" in part:
                a, b = (int(v) for v in part.split("-", 1))
            else:
                a = b = int(part)
            if not 1 <= a <= b <= len(src):
                bad.append(f"{fname}:{part}  (external: not inside a "
                           f"{len(src)}-line file)")
                continue
            if not src.lines[a - 1].strip() or not src.lines[b - 1].strip():
                bad.append(f"{fname}:{part}  (external: an edge of the "
                           f"citation is a blank line)")
                continue
            key = f"{PYTHIA_KEY}{fname}:{part}"
            seen.add(key)
            if not strict:
                continue
            sha = src.block_sha(a, b)
            fresh[key] = {"sha256": sha, "first": clip(src.lines[a - 1]),
                          "last": clip(src.lines[b - 1]),
                          "src": ext_dir.parent.name}
            if record:
                continue
            entry = recorded.get(key)
            if entry is None:
                bad.append(f"{fname}:{part}  (external: not fingerprinted; "
                           f"read the upstream line, then re-record with "
                           f"--record-ranges)")
                continue
            if entry.get("sha256") == sha:
                continue
            hits = src.find_block(entry.get("sha256", ""), b - a + 1)
            where = (f"; it now reads {clip(src.lines[a - 1])!r}" if a == b
                     else "")
            moved = (f" -- it MOVED to :{hits[0]}"
                     + ("" if a == b else f"-{hits[0] + b - a}")
                     if len(hits) == 1 else "")
            bad.append(f"{fname}:{part}  (external: the text recorded "
                       f"against {entry.get('src', 'the PYTHIA tree')}, "
                       f"{entry.get('first', '')!r}, is no longer "
                       f"there{moved}{where}; re-read the row against the "
                       f"PYTHIA version in use, then re-record with "
                       f"--record-ranges)")

    for m in EXT.finditer(text):
        ext_check(m)
    new_text = RANGE.sub(range_repl, text)
    new_text = REF.sub(repl, new_text)
    if fix and fixed:
        DOC.write_text(new_text)
        if remap:
            # the block is byte-identical, only its line numbers moved, so the
            # fingerprint is carried over rather than re-taken
            moved = {remap[k]: v for k, v in recorded.items() if k in remap}
            kept = {k: v for k, v in recorded.items() if k not in remap}
            write_ranges({**kept, **moved})
    if ext_dir is None:
        # nothing was read, so nothing can be said: carry the external
        # fingerprints across untouched rather than dropping them.
        for k, v in recorded.items():
            if k.startswith(PYTHIA_KEY):
                fresh[k] = v
                seen.add(k)
    if strict and not record:
        for key in sorted(set(recorded) - seen - set(remap)):
            bad.append(f"{key}  fingerprinted in {RANGES_REL} but no longer "
                       f"cited by the document; drop it with --record-ranges")
    if record:
        added = sorted(set(fresh) - set(recorded))
        dropped = sorted(set(recorded) - set(fresh))
        changed = sorted(k for k in set(fresh) & set(recorded)
                         if fresh[k]["sha256"] != recorded[k].get("sha256"))
        write_ranges(fresh)
        for k in added:
            print("  recorded  ", k)
        for k in changed:
            print("  re-recorded", k)
        for k in dropped:
            print("  dropped   ", k)
        print(f"{len(fresh)} fingerprints written to {RANGES_REL} "
              f"({len(added)} new, {len(changed)} changed, {len(dropped)} "
              f"dropped): cited BLOCKS, pinned USE SITES and EXTERNAL lines")
    mode = "strict" if strict else "loose"
    print(f"{n} references checked ({mode}), {nr} ranges, {nx} external, "
          f"{len(bad)} broken"
          + (f", {len(allowed)} allow-listed" if strict else "")
          + (f", {len(skipped)} skipped" if skipped else "")
          + (f", {len(fixed)} fixed" if fix else ""))
    for x in fixed:
        print("  fixed", x)
    for a in allowed:                     # an exemption you cannot see is a hole
        print("  allowed", a)
    for sk in skipped:                    # nor one you cannot see
        print("  skipped", sk)
    for b in bad:
        print("  ", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
