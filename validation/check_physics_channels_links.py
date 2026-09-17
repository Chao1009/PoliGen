#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Gate for docs/PHYSICS_CHANNELS.md: every line citation must still point at
what the row says it points at.  Run from the repository root; exit status 1
lists every broken citation.

Since 2026-09-05 (Phase E item E1) the same checks also cover
docs/theory/SPIN32_FINITE_GAMMA.md and docs/PYTHIA_BRIDGE.md -- see EXTRA_DOCS
below for why those two and not others, and how each gets its own allow-list.
Every document's report is printed in turn (the primary document's line is
unlabelled, exactly as before this change; an extra document's line carries
`[filename]` so a broken citation says which document it is in).

Since 2026-09-15 (Phase D item D2) `--records` adds a fifth pass, over the
DATED RUN RECORDS in `docs/open_items/` -- history, gated by nothing until that
date -- under one RELAXED rule, R6, described in full beside `RECORDS_REL`
below.  It is opt-in: `--records` runs it after the strict documents,
`--records-only` runs it alone (the reproducible measurement recipe), and
neither runs in a `--record-ranges`/`--audit-ranges` pass, so every count D1
published with `| tail -1` still names what it named.

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

B. `` `path:line` `` with NO name (43 of them).  Nothing names the target, so
   the assertions are about the line and about the sentence that cites it.

     B1  in bounds, NOT BLANK, and not a bare delimiter (`}`, `};`, `{`, `)`,
         `,`, `*/` ... -- a line with no alphanumeric character on it).  That
         was the whole of rule B until 2026-09-15, and it asserts almost
         nothing: 76.2 % of these citations survive a +1 drift under it, and
         two had been stale for weeks with the gate green
         (`docs/PHYSICS_CHANNELS.md:301` -> `docs/CONVENTIONS.md:357` for
         a_nn, which is now at :383, and `:206` ->
         `include/lipolgen/constants.hpp:86` for the [SS90] bag-model
         coefficient, which is at :136-138).

     B2  THE CITING SENTENCE MUST EVIDENCE THE LINE, by R5's rule (below)
         applied to the cited line and nothing either side of it: a phrase the
         sentence QUOTES must occur on that line as text or as a symbol.  It
         flags 30 of the 42 such citations a3c9ecb carries, against 19 for the
         obvious alternative (R6's shared-token test) -- and lets 2.4 % of +1
         drifts through against that alternative's 16.7 %.  A flagged citation
         is not necessarily wrong; like an R5 refusal it means "this gate
         cannot tell", and the answer is an `UNNAMED_ALLOW` exemption with a
         PIN and a REASON, read and signed one at a time.  28 of the 43 in the
         tree this phase leaves carry one; see RESIDUALS 1.

C. `` `path:first-last` `` -- a RANGE, a block citation: 95 ranges on the primary
   document at the time of writing (the run prints the current count), 5 of them written `` `:first-last` `` and inheriting the
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
     R5  THE BLOCK CONTAINS WHAT THE CITING SENTENCE SAYS IT DOES.  R4 pins a
         block to the content it had WHEN IT WAS RECORDED, and nothing checked
         that the content was the right content: a range recorded onto the
         wrong block is "0 broken" for ever after (three `docs/CONVENTIONS.md`
         ranges on 2026-09-05, two `docs/USAGE.md` ranges on 2026-09-06, every
         one of them found by a human reading the row, none by this gate).  So
         since 2026-09-06 `--record-ranges` REFUSES to record a block that the
         citing sentence does not evidence, and refusing means the entry is not
         written, so the next ordinary run reports the range as "not
         fingerprinted" and the gate stays red until a human acts.  The
         evidence is taken from THE CITING SENTENCE -- the `RECORD_WIN`
         characters on either side of the citation, clipped to the enclosing
         paragraph, where a markdown TABLE ROW is its own paragraph (the claim
         cell and the provenance cell of one row are one citing sentence; the
         row above is not).  Three branches, in order:

           (a) an ADJACENT NAME (`` `path:first-last` `name` ``, the same
               one-space adjacency rule the point citations use) must be
               DECLARED OR USED inside the block.  Declared-or-used, not R3's
               declared: a range is a block citation and citing the three lines
               that CALL a thing is legitimate.  R3 is the stricter rule and
               still runs first, so an identifier name that fails R3 never
               reaches this;
           (b) failing that, a QUOTED PHRASE -- backticks, "straight" or
               "curly" quotation marks -- inside the window, excluding the
               citations themselves.  At least one such phrase must occur in
               the block, either as TEXT (whitespace-normalised substring, 3
               characters or more) or as a SYMBOL (an identifier, qualified
               name, call or CLI flag of 2 characters or more, declared or used
               in the block).  At least one, not all: a sentence quotes several
               things and only one of them need be the block;
           (c) neither: the citation says nothing a machine can check against
               the block, and the entry needs an explicit per-entry exemption
               in `RANGE_ALLOW` below -- a PIN and a REASON.  The pin is a
               substring that occurs inside the block and NOWHERE ELSE in that
               file, so the exemption names the block by its content and cannot
               slide onto another one; the reason is what the human who read
               the row concluded.  The refusal prints what it would have
               needed, quoted phrase by quoted phrase.

         PER CITATION, not per block.  Until 2026-09-15 the verdict was pooled
         on the block key -- "a block cited by two documents is evidenced if
         EITHER sentence evidences it" -- and `RANGE_ALLOW` was keyed on the
         range, so a citation re-pointed onto an already-blessed block was
         accepted whatever its own sentence said (three of five deliberately
         wrong-block ranges passed that way, and 11 citations on 7 keys in the
         committed documents had never been judged against their own sentence
         at all).  The verdict is now keyed on (document, citing line, range),
         the exemption with it, and a block is fingerprinted only when EVERY
         citation of it is evidenced or exempted.  The citing line is part of
         the key on purpose: an exemption records that a human read THAT
         sentence, so a document edit that moves the sentence retires the
         exemption loudly ("no range citation reaches it") instead of quietly
         covering a sentence nobody read.

         MEASURED against the committed tree this rule was written on
         (a3c9ecb), measured 2026-09-15 with `RANGE_ALLOW` suppressed: of the 180
         distinct blocks the three covered documents cite, 91 were evidenced
         (71 by quoted text, 20 by quoted symbol, 0 by an adjacent name -- no
         range in any of the three carried one) and **89 were REFUSED**.
         Twenty of those 89 turned out to be pointed at the wrong block and
         were re-pointed; the other 69 are citations whose sentence says
         nothing a machine can check -- the common shape in these documents is
         prose that names the subject in unquoted words ("provenance at", "the
         inventory table at") -- and they are the exemptions below.  So a
         refusal is not a synonym for a wrong citation; it is "this gate
         cannot tell", and 22 % of the time here it was a wrong citation.
         Recipe: `git archive a3c9ecb | tar -x -C <dir>`, copy this file into
         `<dir>/validation/`, and run `python3
         validation/check_physics_channels_links.py --audit-ranges
         --no-range-allow` there.  `--audit-ranges` runs the rule and writes
         nothing; `--no-range-allow` suppresses the exemption table.

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

   R4 has a running cost, and it is the point: editing a block that a covered
   document cites breaks this gate until a human confirms the row and
   re-records (and since R5, re-recording is itself refused unless the citing
   sentence still evidences the block).  `--record-ranges` prints every entry it adds, changes or
   drops.  The same sidecar holds the S3 use-site pins and the D fingerprints
   below, now SHARED across every document this gate covers (a block or
   upstream line cited by more than one document hashes to one entry, not
   one per document): 220 entries in all -- 180 blocks, 31 use-site lines
   and 9 upstream lines (re-counted 2026-09-15 from the file itself: a range
   entry is the one with a `first`/`last` pair, a use-site entry the one with
   a `line`, an upstream entry the one keyed `pythia8:`; the 181/30 this line
   carried was one out in each), as of 2026-09-05 once docs/theory/
   SPIN32_FINITE_GAMMA.md and docs/PYTHIA_BRIDGE.md joined
   docs/PHYSICS_CHANNELS.md under this gate (EXTRA_DOCS, below).

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
`--record-ranges` (re)writes the sidecar, refusing any range R5 does not
evidence (the refused entry is left unwritten, so the gate reports it as
unfingerprinted until it is re-pointed or exempted).  `--audit-ranges` applies
R5 and writes nothing; adding `--no-range-allow` suppresses `RANGE_ALLOW`, so
`--audit-ranges --no-range-allow` measures how many entries the rule refuses on
their own evidence.

RESIDUALS -- what this gate still cannot see, measured 2026-09-04.

  A DRIFT'S SIGN, stated once because these figures are meaningless without it.
  A "+k drift" is k lines LOST above the citation, so the document's number now
  points at what used to be k lines LATER: the citation survives if line + k
  was an accepting line.  A "-k drift" is k lines INSERTED above it.  Both are
  measured end to end -- insert or delete k lines at the top of all 77 cited
  files, run the gate, count the citations it does NOT report.

  1. An UNNAMED point citation (rule B) was pinned only to "a non-blank,
     non-trivial line" until 2026-09-15.  Of the 19 there were then, 14
     survived a +1 drift and 12 a -1 (10 both); at 2 lines it was 17 and 14.
     B2 closes most of that -- 2.4 % of +1 drifts survive it -- but what
     remains is the 28 citations B2 cannot read, which are exempted in
     `UNNAMED_ALLOW` and are pinned by their exemption's pin and by nothing
     else: a drift that carries the pinned text with it is invisible, exactly
     as it is for a RANGE_ALLOW pin.  Naming the symbol in the document is
     still the real fix, and is preferred whenever the target has a name.

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
from collections import Counter
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
    ("include/lipolgen/xsec.hpp", 221, "g1_rank3"): (
        "no rank-3 (octupole) slot",
        "no such symbol exists in the tree, deliberately: the row cites the "
        "comment that states why there is no rank-3 slot"),
    ("src/core/pipeline.cpp", 847, "coherent"): (
        "on the coherent channel the tensor signal is in the recoil",
        "the cited line is inside PipelineConfig::validate()'s refusal "
        "message; `coherent` there is the English word in that message, which "
        "is the text the row quotes"),
    ("src/core/pipeline.cpp", 1877, "optics_lumi_factor"): (
        "cfg_.lumi_pb * optics_lumi_factor()",
        "the row's claim is that luminosity mode multiplies the optics factor "
        "in, so it cites the multiplication, not the accessor (declared at "
        "include/lipolgen/pipeline.hpp:780, cited there too)"),
    ("src/core/rc.cpp", 1532, "is_tagged_channel"): (
        "is_tagged_channel(channel_)",
        "one of three consecutive branch lines of RcModel::exclusion_reason "
        "the row cites together; the neighbours are use sites in files that "
        "do not declare their symbol"),
    ("src/core/breakup.cpp", 346, "proton_fraction"): (
        "proton_fraction(in.x, in.q2, 1, 3)",
        "the triton branch's own proton/neutron draw, which is what the row "
        "describes; the definition is generic over (Z, A)"),
    ("src/pythia/pythia_bridge.cpp", 701, "dis_parton_fraction"): (
        "dis_parton_fraction(q, p_n, &xi_tmp, 1.0, mq)",
        "the re-solve with the chosen quark's mass, cited beside the "
        "definition at src/pythia/pythia_bridge.cpp:171 in the same row"),
}

# R5, the recording rule.  `RECORD_WIN` is the "within N characters of the
# citation" of branch (b): the citing sentence is taken as this many characters
# on either side of the citation, clipped to the enclosing paragraph (a
# markdown table row being its own paragraph).  200 was chosen by measuring the
# alternatives over the three covered documents (2026-09-06, on the tree this
# rule was landed in): 120 / 200 / 400 characters evidence 94 / 112 / 126 of the
# 211 range citations while accepting a DECOY block -- the same citation pointed
# at a block of the same length shifted clear of the true one -- 10.8 % / 14.3 %
# / 20.4 % of the time, so the true-to-decoy ratio falls 4.5 / 4.1 / 3.3 and
# every widening is paid for.  The whole of the gain from 120 to 200 is a table
# row's claim cell reaching its own provenance cell, which is the same sentence.
# Clipping at the table-cell `|` instead (so that only the cell holding the
# citation counts) evidences 80 of the 211 rather than 112: those 32 are rows
# whose claim cell quotes the very symbol the cited block declares, so the cell
# boundary is not the sentence boundary in these documents and is not used as
# one.  The numbers and their recipes are in
# `docs/open_items/run_2026-09-06/phase_D_numbers.md` sec. D1.6.
RECORD_WIN = 200
# A quoted phrase: a backticked span, or a "straight" or “curly” quotation.
# A phrase may WRAP: these documents are hard-wrapped at ~80 columns, so the
# phrase `docs/CONVENTIONS.md:72-73` is cited for -- "b₁ > 0 means the m = 0
# state has the LARGER cross section" -- runs across a line break, and a regex
# that stopped at the newline refused that citation for quoting nothing while
# the quotation sat in front of it.  A blank line still ends a phrase (nothing
# in these documents quotes across a paragraph), and 240 characters caps it.
QUOTED = re.compile(r"`([^`]{1,240}?)`|\"([^\"]{3,240}?)\"|“([^”]{3,240}?)”")
BLANK_LINE = re.compile(r"\n[ \t]*\n")
MIN_PHRASE = 3                      # characters, for the TEXT branch of (b)
MIN_SYMBOL = 2                      # characters, for the SYMBOL branch of (b)
# The MINIMUM EVIDENCE a text match has to be (see `text_evidence_ok`): a
# one-word phrase on more lines than this of the target file evidences nothing,
# and no text match evidences a block longer than this.  Both were measured in
# D3.5 and both are reproducible with `validation/record_rule_sweep.py r5`.
EVIDENCE_WORD_LINES = 5
EVIDENCE_MAX_BLOCK = 150

# R5 (c): ranges whose citing sentence evidences nothing a machine can read.
# Keyed on the EXACT range -- `path:first-last`, the sidecar's own key -- so an
# exemption covers that one block and no other, and carrying, like ALLOW above,
# a PIN and a REASON.  The pin must occur inside the block and NOWHERE ELSE in
# the file: that is what makes the exemption name a block by its content rather
# than by two line numbers that drift.  Every entry here was read against its
# block on 2026-09-15, one at a time, before it was written:
# `--audit-ranges --no-range-allow` prints each refusal's citing sentence and
# the phrases the rule looked for, and that listing is the worksheet they were
# filled in from.
RANGE_ALLOW: dict[tuple[str, int, str], tuple[str, str]] = {
    ("PHYSICS_CHANNELS.md", 83, "include/lipolgen/spin.hpp:19-21"): (
        "///   octupole J=3/2   O    = <J_z^3 - (41/20) J_z>/(3/10)",
        "the row says the (41/20, 3/10) octupole normalisation is "
        "stated only here; the block is the three moment lines, the "
        "last of which is that normalisation"),
    ("PHYSICS_CHANNELS.md", 90, "tests/test_tagged.cpp:396-413"): (
        "CHECK_CLOSE_AT(m.p2_moment_mixture(pv), -mom.tensor / 5.0, "
        "0.0, tol);",
        "the row calls ⟨P₂(cos θ_k)⟩ = −T/5 a test-pinned identity; the "
        "block is the test case that pins it, including the mixture "
        "form"),
    ("PHYSICS_CHANNELS.md", 110, "src/core/xsec.cpp:302-305"): (
        "const double helicity = state.lam_e * state.pe;",
        "the row's claim is that the vector-L term is evaluated only "
        "when lam_e*P_e != 0 and m != 0, so an m = 0 state has no "
        "vector-L term for any spin; the block is that guard and the "
        "term it protects"),
    ("PHYSICS_CHANNELS.md", 123, "docs/open_items/physics_literature.md:143-151"): (
        "**(d) Independent fourth check, POLRAD 2.0**",
        "the row cites [POLRAD]'s Born Eqs. (9)/(10) as an independent "
        "fourth sign check; the block is the register entry that makes "
        "that check, (b) through (d)"),
    ("PHYSICS_CHANNELS.md", 191, "include/lipolgen/sf.hpp:43-45"): (
        "/// R = sigma_L/sigma_T, simplified R1990-like magnitude "
        "(TOY).  The DEFAULT R",
        "the row says both defaults are explicitly labelled TOY here "
        "and at :109-112; the block is the R default's comment, which "
        "carries that label"),
    ("PHYSICS_CHANNELS.md", 203, "docs/CONVENTIONS.md:207-211"): (
        "the EPPS21 depletion is the single stored number",
        "the row's provenance for EMC_VALENCE_DEPLETION_EPPS21 = "
        "0.031052077003862335 and for the superseded CT18NLO 0.02979; "
        "the block states both and why the denominator is CT18ANLO"),
    ("PHYSICS_CHANNELS.md", 203, "include/lipolgen/sf.hpp:322-336"): (
        "/// PROVENANCE: "
        "`polli_fastsim.polarized.valence_depletion(mode=\"epps21\")` --",
        "the row's provenance for EMC_VALENCE_DEPLETION_EPPS21, its "
        "301-point window and the superseded CT18NLO value; the block "
        "is that provenance comment and the constant"),
    ("PHYSICS_CHANNELS.md", 204, "include/lipolgen/sf.hpp:341-343"): (
        "/// eq = 23 is the R^{3/2 3/2}_{As} of their Eq. (23), eq = 26 "
        "the",
        "the row reads eq = 23 as R^{3/2 3/2}_As and eq = 26 as R^{(3/2 "
        "1)}_As 'per' this block; the block is the comment that assigns "
        "exactly those two"),
    ("PHYSICS_CHANNELS.md", 208, "include/lipolgen/pipeline.hpp:931-934"): (
        "-0.0818(17) fm^2 (`LI6_QUADRUPOLE_FM2`, TUNL A = 6, 1998CE04",
        "the row says the measured Q(⁶Li) and Q_d are quoted here; the "
        "block is the comment carrying −0.0818(17) fm² against Q_d = "
        "+0.2859(3) fm² with their sources"),
    ("PHYSICS_CHANNELS.md", 210, "include/lipolgen/sf.hpp:331-334"): (
        "4.2 % shallower",
        "the row says the only quantified CT18NLO-vs-CT18ANLO "
        "difference in the repo is the 4.2 % depletion shift; the block "
        "is the comment that states it, with both constants"),
    ("PHYSICS_CHANNELS.md", 230, "include/lipolgen/generator.hpp:22-28"): (
        "Bacchetta et al. (JHEP 02 (2007) 093) azimuth",
        "the row says [Bacchetta07]'s phi_S convention is cited only "
        "here and echoed in sampler.hpp, and is registered nowhere; the "
        "block is that citation and the sign convention it fixes"),
    ("PHYSICS_CHANNELS.md", 230, "include/lipolgen/sampler.hpp:26-29"): (
        "phi_S of the alignment axis exactly (massless target; "
        "`reco.py`",
        "two rows say [Bacchetta07]'s φ_S convention is echoed here; "
        "the block is that echo"),
    ("PHYSICS_CHANNELS.md", 236, "include/lipolgen/bookkeeping.hpp:170-172"): (
        "/// Spin-1 A_zz run plan: equal-thirds fills (pz, +pzz), (-pz, "
        "+pzz) and the",
        "the row cites [HERMES05] for the thirds pattern; the block is "
        "the doc comment that names the HERMES-style pattern and the "
        "three fills"),
    ("PHYSICS_CHANNELS.md", 266, "data/vmc/README.md:85-118"): (
        "the cluster MOMENTUM DISTRIBUTIONS (fetched 2026-08-29, second "
        "pass)",
        "the row says the repo carries no journal citation for the "
        "momenta/ AV18+UX files beyond Wiringa's ANL page and the "
        "Wayback URLs in this block; the block is that section of the "
        "inventory"),
    ("PHYSICS_CHANNELS.md", 266, "docs/open_items/vmc_reconciliation.md:120-131"): (
        "| k < 0.678 fm⁻¹ = 0.134 GeV (below the **S** node) | **−1** |",
        "the row says the VMC node values are recorded here; the block "
        "is the sign table with the S node at 0.678 fm⁻¹ = 0.134 GeV "
        "and the D node at 2.250"),
    ("PHYSICS_CHANNELS.md", 267, "data/vmc/README.md:65-103"): (
        "## Files fetched, one row per subdirectory",
        "the row cites the VMC inventory table and its Wayback URLs in "
        "unquoted prose; the block is that inventory (both fetch "
        "passes) with the Wayback-snapshot column"),
    ("PHYSICS_CHANNELS.md", 271, "docs/open_items/vmc_reconciliation.md:197-206"): (
        "| VMC D-sign flipped (= the PRE-FIX physics column) |",
        "the row calls this 'that table, REGENERATED on the fixed "
        "library'; the block is the regenerated k-table whose header "
        "names the pre-fix physics column"),
    ("PHYSICS_CHANNELS.md", 274, "docs/surveys/beagle_survey.md:120-127"): (
        "### 2f. The eD mode with spectator tagging",
        "the row's [BeAGLE22] provenance for the "
        "struck-cluster/spectator split; the block is the survey "
        "section on BeAGLE's eD spectator-tagging mode"),
    ("PHYSICS_CHANNELS.md", 301, "docs/surveys/beagle_survey.md:120-127"): (
        "- **No FSI at all in eD** (paper, §II).",
        "the row's provenance for the impulse approximation and for "
        "BeAGLE's four-momentum-balance-only closure; the block is the "
        "eD section, whose last line is the no-FSI statement the row "
        "leans on"),
    ("PHYSICS_CHANNELS.md", 304, "docs/DEVELOPMENT_PLAN.md:61-65"): (
        "| tier | content | source of truth |",
        "the row points at the in-repo pointers to plans/05 step 5.D; "
        "the block is the T0/T1/T2 tier table whose T1 row names "
        "`plans/05` 5.D"),
    ("PHYSICS_CHANNELS.md", 345, "include/lipolgen/fsi.hpp:213-220"): (
        "a^2 = (r_ch^2(cluster) - r_ch^2(p)) / 3 .",
        "the row says the a² values are quoted from measured charge "
        "radii with no bibliographic source; the block is that formula "
        "and the α number 0.7001 fm² it produces"),
    ("PHYSICS_CHANNELS.md", 346, "docs/open_items/physics_literature.md:81-85"): (
        "**σ_eff = 30–70 mb** rising with W (Cosyn–Sargsian Deeps fit)",
        "two rows cite this block: for the register's only fitted σ(W) "
        "(30–70 mb rising) and for its calling the coherent-α amplitude "
        "a variant; both sentences are inside it"),
    ("PHYSICS_CHANNELS.md", 346, "include/lipolgen/fsi.hpp:176-179"): (
        "THE OPEN PHYSICS INPUT.  sigma_XN(W) is the weakest number in "
        "the model.",
        "the row says the decreasing σ anchors encode a "
        "formation-length argument; the block is the comment that makes "
        "that argument and calls σ_XN(W) the weakest number in the "
        "model"),
    ("PHYSICS_CHANNELS.md", 425, "include/lipolgen/coherent.hpp:631-637"): (
        "ONE thing, the t_min kinematic cut; it never forms M_X",
        "the row says the Python recopseudo.CoherentResponse drew x_P "
        "only for a t_min cut and never formed M_X; the block is the "
        "comment that states exactly that and why LiPolGen cannot"),
    ("PHYSICS_CHANNELS.md", 428, "docs/PYTHIA_BRIDGE.md:574-575"): (
        "veto = 0 at `M_X ≥ 1.4`.",
        "cited as the CHAIN TEST for the M_X floor; the block is the "
        "sentence that says the T2 chain test re-measures the veto "
        "table per run and pins veto = 0 at M_X >= 1.4"),
    ("PHYSICS_CHANNELS.md", 429, "src/core/coherent.cpp:290-291"): (
        "// 6Li: TUNL A=6 (Tilley et al. NPA 708:3).  7Li: recomputed "
        "from",
        "the row cites [TUNL6] for the ⁶Li breakup thresholds; the "
        "block is the comment that names Tilley et al. NPA 708:3 as "
        "their source"),
    ("PHYSICS_CHANNELS.md", 449, "include/lipolgen/cluster_config.hpp:19-71"): (
        "Q_matter(6Li, M) = (3M^2-2)",
        "the row says the (G1)/(G4)/(G5)/(G6) identities are derived in "
        "place here; the block is that derivation, (G1) through (G8)"),
    ("PHYSICS_CHANNELS.md", 451, "docs/USAGE.md:2795-2992"): (
        "## 9. Polarized ⁶Li configurations for coherent-diffraction "
        "codes",
        "the row cites USAGE sec. 9 in full for the configuration "
        "format and its committed example sidecar; the block is that "
        "section (198 lines, so no text match evidences it since "
        "EVIDENCE_MAX_BLOCK)"),
    ("PHYSICS_CHANNELS.md", 489, "docs/PYTHIA_BRIDGE.md:145-146"): (
        "worst |λ − 1| = 1.2 × 10⁻¹³ over",
        "cited as the MEASUREMENT behind the row's worst |λ−1| ≈ "
        "1.2e-13; the block is that measurement over 400 events"),
    ("PHYSICS_CHANNELS.md", 489, "docs/PYTHIA_BRIDGE.md:219-230"): (
        "## 5. The frame map",
        "cited as the DERIVATION of the frame map the row describes; "
        "the block is §5 and its triad construction"),
    ("PHYSICS_CHANNELS.md", 489, "docs/PYTHIA_BRIDGE.md:429-432"): (
        "| worst relative 4-momentum deviation | 8.4 × 10⁻¹⁴ |",
        "cited as a MEASUREMENT for the row's 8.4e-14 and 1.2e-13; the "
        "block is the performance table carrying both"),
    ("PHYSICS_CHANNELS.md", 491, "docs/PYTHIA_BRIDGE.md:137-161"): (
        "P(q) ∝ e_q² · x f_q(ζ, max(Q², q2_pdf_min))",
        "the row's subject is the bridge's flavour choice and its "
        "q2_pdf_min floor (and the c/b table mass it sits beside); the "
        "block is §2's m_q rule and §3 Flavour choice, including the "
        "sampling probability that carries q2_pdf_min"),
    ("PHYSICS_CHANNELS.md", 495, "docs/T2_CHAIN.md:79-96"): (
        "`Role::StruckNucleon` verbatim; its `Role::StruckCluster` "
        "branch is",
        "the row is about which nucleon the bridge is handed; the block "
        "is the T2 record of the struck-nucleon hook, its deprecated "
        "cluster branch and the measured disappearance of the "
        "no-surrogate tail"),
    ("PHYSICS_CHANNELS.md", 495, "include/lipolgen/pythia_bridge.hpp:358-373"): (
        "DEPRECATED (2026-08-30, superseded by the T1 tier)",
        "the row is about which nucleon the bridge is handed; the block "
        "is the deprecated NucleonInCluster hook's own comment, which "
        "is the alternative the row's NucleonChoice replaced"),
    ("PHYSICS_CHANNELS.md", 496, "docs/PYTHIA_BRIDGE.md:473-475"): (
        "**990 is a legal user beam**",
        "the row cites this for the two PYTHIA upstream pointers it "
        "names (BeamSetup.cc:869-875, BeamParticle.cc:178); the block "
        "is the sentence that carries both and says why 990 works as a "
        "beam"),
    ("PHYSICS_CHANNELS.md", 499, "docs/surveys/beagle_survey.md:120-127"): (
        "`ISTHKK == 14`",
        "the row says the intact-recoil roles are never read and never "
        "touched -- the BeAGLE light-nucleus rule; the block is the eD "
        "section that shows BeAGLE finding the spectator and leaving it "
        "alone"),
    ("PHYSICS_CHANNELS.md", 500, "docs/PYTHIA_BRIDGE.md:360-375"): (
        "is the statement that is",
        "the row lists the three HFS identities; the block is the "
        "numbered list of exactly those three, ending with the exact "
        "hfs_sigma_empz_exact one"),
    ("PHYSICS_CHANNELS.md", 500, "docs/T2_CHAIN.md:102-106"): (
        "the closed-form collinear `hfs_sigma_empz_truth` is now "
        "approximate at the",
        "the row's first HFS identity is the approximate one; the block "
        "is the T2 measurement of how approximate (1.4 % ⁶Li α, 0.7 % "
        "⁷Li) and why"),
    ("PHYSICS_CHANNELS.md", 514, "docs/CONVENTIONS.md:325-327"): (
        "X must be TIMELIKE on",
        "the row's claim is that the hadronic X is timelike on every "
        "channel and is checked, never clipped; the block is the "
        "convention that says exactly that"),
    ("PHYSICS_CHANNELS.md", 514, "include/lipolgen/pipeline.hpp:13-39"): (
        "WHAT X IS, PER CHANNEL.",
        "the row states the per-channel hadronic-X balance; the block "
        "is the header comment that writes that balance out for "
        "inclusive, tagged and coherent"),
    ("PHYSICS_CHANNELS.md", 515, "include/lipolgen/pipeline.hpp:316-328"): (
        "Fidelity tier of the final state a run writes",
        "the row says no CLI or config-file key sets the tier, only the "
        "PipelineConfig field; the block is that field's own T0/T1/T2 "
        "documentation"),
    ("PHYSICS_CHANNELS.md", 518, "docs/HEPMC3_CONVENTION.md:6-12"): (
        "**proposed convention for ion-spin states in HepMC3** flagged "
        "as open item",
        "the row says the ion-spin attribute convention is LiPolGen's "
        "own proposal and open item #17; the block is the paragraph "
        "that states it and that nothing upstream exists"),
    ("PHYSICS_CHANNELS.md", 518, "docs/HEPMC3_CONVENTION.md:95-121"): (
        "- **Role → status/PDG map**, complete, with the tier that "
        "writes each row:",
        "the row cites the 'role/status table'; the block is that table "
        "and the PartnerSpectator note under it"),
    ("PHYSICS_CHANNELS.md", 1002, "include/lipolgen/sf.hpp:152-158"): (
        "Wandzura-Wilczek g2(x) = -g1(x)",
        "the [WW] bibliography entry says the relation is carried by "
        "name only, with no journal or arXiv identifier; the block is "
        "the g2_ww comment and declaration that carry the name"),
    ("PHYSICS_CHANNELS.md", 1036, "include/lipolgen/beams.hpp:25-41"): (
        "/// Physical nuclear masses [GeV] -- the same AME2020-derived "
        "values as",
        "cited from the [AME2020] bibliography entry as where the "
        "AME2020 masses are tabulated; the block is that table's header "
        "comment and nucleus_mass"),
    ("PHYSICS_CHANNELS.md", 1041, "include/lipolgen/generator.hpp:22-28"): (
        "`reco.azimuth_wrt_lepton_plane` is written and tested",
        "the [Bacchetta07] bibliography entry says the convention is "
        "cited in this block only; the block is that citation"),
    ("SPIN32_FINITE_GAMMA.md", 20, "include/lipolgen/asymmetries.hpp:57-61"): (
        "M is the FREE",
        "the note's conventions paragraph says `M` is the free nucleon "
        "mass because x is per-nucleon; the block is that comment and "
        "`gamma_squared`"),
    ("SPIN32_FINITE_GAMMA.md", 22, "include/lipolgen/xsec.hpp:121-128"): (
        "struct EventSpinState {",
        "the note's convention line defines λ_e and P_e; the block is "
        "EventSpinState, which carries lam_e and pe (and the axis "
        "angles the same sentence names)"),
    ("SPIN32_FINITE_GAMMA.md", 32, "include/lipolgen/xsec.hpp:23-29"): (
        "RANK-2 GEOMETRY",
        "the note says the rank-2 block reuses the spin-1 machinery "
        "verbatim through one Q_NN geometry; the block is the header's "
        "RANK-2 GEOMETRY paragraph that states it for any J"),
    ("SPIN32_FINITE_GAMMA.md", 33, "src/core/xsec.cpp:175-180"): (
        "InclusiveKernel::tensor_moments",
        "the note cites the code that computes the Q_NN geometry; the "
        "block is `tensor_moments`, which returns (Q_NN, 3 Q_NN)"),
    ("SPIN32_FINITE_GAMMA.md", 34, "include/lipolgen/xsec.hpp:204-217"): (
        "(DEFAULT FALSE) selects the tensor b-sector kernel",
        "the note cites the optional exact finite-gamma Cosyn kernel; "
        "the block is the `tensor_gamma` option that selects it and "
        "states the default"),
    ("SPIN32_FINITE_GAMMA.md", 36, "src/core/xsec.cpp:108-118"): (
        "b1f = &b1_32_func_;",
        "the note says the J = 3/2 slots are dispatched here; the block "
        "is the spin-1.5 branch that assigns them"),
    ("SPIN32_FINITE_GAMMA.md", 38, "include/lipolgen/spin.hpp:94-102"): (
        "/// Normalized rank-3 moment for J=3/2: <Jz^3 - (41/20) "
        "Jz>/(3/10).",
        "the note says the rank-3 (octupole) sector is computed in the "
        "spin bookkeeping; the block is octupole_moment and its comment"),
    ("SPIN32_FINITE_GAMMA.md", 38, "src/core/spin.cpp:246-259"): (
        "octupole moment defined for j = 3/2 only",
        "the note says the rank-3 sector is computed in the spin "
        "bookkeeping but never reaches the cross section; the block is "
        "`octupole_moment`"),
    ("SPIN32_FINITE_GAMMA.md", 107, "include/lipolgen/bookkeeping.hpp:59-68"): (
        "/// Populations p_m ordered m = +J ... -J along the axis "
        "(library-wide).",
        "the note says the library stores an axially symmetric fill: "
        "populations p_m along one axis n̂(θ_S, φ_S), m ordered +J … "
        "−J; the block is SpinCategory, which is that storage"),
    ("SPIN32_FINITE_GAMMA.md", 107, "include/lipolgen/spin.hpp:9-14"): (
        "/// n(theta_S, phi_S): populations p_m with sum(p_m) = 1, "
        "ordered m = +J ... -J",
        "the note says the library stores an axially symmetric fill "
        "with p_m ordered +J … −J; the block is the comment that "
        "defines exactly that and the lab-frame rotation"),
    ("SPIN32_FINITE_GAMMA.md", 110, "include/lipolgen/spin.hpp:18-21"): (
        "///   octupole J=3/2   O    = <J_z^3 - (41/20) J_z>/(3/10)",
        "the note says the three numbers are these; the block is the "
        "four comment lines defining P, P_zz / T and O"),
    ("SPIN32_FINITE_GAMMA.md", 110, "src/core/spin.cpp:261-290"): (
        "AxisMoments moments_along_axis(double j, const "
        "std::vector<double>& populations)",
        "the note says the three axial moments are implemented here; "
        "the block is moments_along_axis, which computes all three"),
    ("SPIN32_FINITE_GAMMA.md", 141, "src/core/spin.cpp:269-275"): (
        "out.tensor = t;",
        "the note says that for J = 1 the code uses P_zz = ⟨3J_z² − 2⟩ "
        "with no division by 3; the block is that J = 1 branch, where t "
        "is assigned undivided"),
    ("SPIN32_FINITE_GAMMA.md", 144, "src/core/xsec.cpp:175-180"): (
        "const double q_nn = (3.0 * m * m - j * (j + 1.0)) / 3.0;",
        "the note says the cross-section kernel divides by 3 itself, "
        "Q_NN = (3m^2 - J(J+1))/3 = c_m/3; the block is that line and "
        "its function"),
    ("SPIN32_FINITE_GAMMA.md", 145, "include/lipolgen/xsec.hpp:23-29"): (
        "is ONE geometry for both spins",
        "the note says the kernel gives one rank-2 geometry for both "
        "spins; the block is the header paragraph that says exactly "
        "that"),
    ("SPIN32_FINITE_GAMMA.md", 145, "tests/test_xsec.cpp:297-316"): (
        "the rank-2 geometry is one formula for both spins",
        "the note says the one-geometry statement is pinned by a test; "
        "the block is that TEST_CASE"),
    ("SPIN32_FINITE_GAMMA.md", 146, "tests/test_xsec.cpp:318-335"): (
        "the J = 1 geometry is the HJM transcription digit for digit",
        "the note says the geometry is pinned against the HJM "
        "transcription digit for digit; the block is that TEST_CASE"),
    ("SPIN32_FINITE_GAMMA.md", 165, "src/core/spin.cpp:337-342"): (
        "a[3][i] = (ms[i] * ms[i] * ms[i] - (41.0 / 20.0) * ms[i]) / "
        "0.3;",
        "the note says the solve's rows are exactly the four weight "
        "vectors (1, m/J, Q_NN, R₃); the block is those four rows of "
        "the 4×4 matrix"),
    ("SPIN32_FINITE_GAMMA.md", 178, "include/lipolgen/xsec.hpp:24-26"): (
        "///   t_geo = Q_NN P_2(cos theta_S)   and   c_eff = 3 Q_NN",
        "the note says line 23 carries the t_ij form and 24-25 the Q_NN "
        "definition with t_geo/c_eff; the block is those three lines"),
    ("SPIN32_FINITE_GAMMA.md", 218, "include/lipolgen/bookkeeping.hpp:214-224"): (
        "/// The geometric (spin-temperature) population ladder p_m ~ "
        "t^m with",
        "the note calls the default branch the spin-temperature "
        "(maximum-entropy) ladder p_m ∝ u^m; the block is "
        "SpinTemperatureLadder and the comment that states the same "
        "ladder"),
    ("SPIN32_FINITE_GAMMA.md", 242, "tests/test_spin.cpp:223-232"): (
        "// the J = 3/2 anchor: pz = 0.7 gives tensor 0.4 and octupole "
        "0.2 exactly,",
        "the note says the u = 3 rational anchor -- populations (27, 9, "
        "3, 1)/40, P_z = 0.7, T = 0.4, R₃ = 0.2 -- is pinned here; the "
        "block is that half of the test case"),
    ("SPIN32_FINITE_GAMMA.md", 244, "include/lipolgen/bookkeeping.hpp:217-218"): (
        "/// giving (P_z, rank-2) = (8/13, 4/13) and (7/10, 2/5).",
        "the note says the u = 3 rational anchor is quoted here; the "
        "block is the two comment lines that carry (27,9,3,1)/40 and "
        "(P_z, rank-2) = (7/10, 2/5)"),
    ("SPIN32_FINITE_GAMMA.md", 277, "tests/test_pipeline.cpp:608-634"): (
        "TEST_CASE(\"pipeline: 7Li generated events give <P2(cos "
        "theta_k)> = -T/5\") {",
        "the note cites this as the pipeline-level ⟨P₂(cos θ_k)⟩ = −T/5 "
        "polarimeter test; the block is that test case"),
    ("SPIN32_FINITE_GAMMA.md", 277, "tests/test_tagged.cpp:396-412"): (
        "// the in-situ alignment polarimeter: <P2(cos theta_k)> = -T/5 "
        "for ANY fill",
        "the note cites this as the ⟨P₂(cos θ_k)⟩ = −T/5 polarimeter at "
        "the model level; the block is that test case"),
    ("SPIN32_FINITE_GAMMA.md", 360, "src/core/xsec.cpp:175-180"): (
        "return std::make_pair(q_nn, 3.0 * q_nn);",
        "the note's mapping result: the code extends Q_NN from J = 1 to "
        "J = 3/2 with T = +-1 on the stretched state; the block is "
        "where it does"),
    ("SPIN32_FINITE_GAMMA.md", 400, "src/core/xsec.cpp:303-305"): (
        "out.w_avg = out.w_avg + helicity * v * ct * a_parallel(t, x, "
        "q2, y);",
        "the note writes the same vector term as `w = λ_e P_e (m/J) cos "
        "θ_S · A_∥`; the block is the three lines that compute it -- "
        "re-pointed in run 2026-09-06 from :267-269"),
    ("SPIN32_FINITE_GAMMA.md", 425, "tests/test_xsec.cpp:223-240"): (
        "Cosyn Eq. 27: A_zz(theta_S=0)(1 + eps R) = -(2/3) b1/F1",
        "the note says the -(2/3) b1/F1 sign is pinned at every y by "
        "two tests; the block is the first of them"),
    ("SPIN32_FINITE_GAMMA.md", 425, "tests/test_xsec.cpp:242-255"): (
        "the program sign IS the literature sign, deliberately",
        "the note's second pin for the published sign; the block is the "
        "TEST_CASE that fixes TENSOR_LL_SIGN = -1"),
    ("SPIN32_FINITE_GAMMA.md", 736, "include/lipolgen/xsec.hpp:263-272"): (
        "/// theta_q between q and the beam, so the alignment tensor of "
        "Eq. (9)",
        "the note derives the spin axis in the photon frame at finite "
        "γ; the block is the comment that sets up exactly that frame "
        "and the θ_q it turns on"),
    ("SPIN32_FINITE_GAMMA.md", 741, "include/lipolgen/xsec.hpp:274-275"): (
        "N = (c s_S cos phi' + s c_S,  s_S sin phi',",
        "the note says its Eq. (35) is verbatim these lines, and the "
        "equation is displayed in the note as an indented block, which "
        "carries no quotation marks for the rule to read; the block is "
        "that vector N"),
    ("SPIN32_FINITE_GAMMA.md", 770, "include/lipolgen/asymmetries.hpp:57-61"): (
        "because x is per-nucleon",
        "the variable note says the generator uses x per nucleon "
        "throughout, which is why its x is [2]'s x_d; the block is the "
        "per-nucleon comment that says so"),
    ("SPIN32_FINITE_GAMMA.md", 790, "include/lipolgen/xsec.hpp:99-103"): (
        "for ANY b2 and not only at the tensor Callan-Gross point.",
        "the note says the header states the γ = 0 collapse 'for any "
        "b2'; the block is that header comment, which writes it in "
        "capitals"),
    ("SPIN32_FINITE_GAMMA.md", 790, "src/core/xsec.cpp:29-32"): (
        "o.f_l = (2.0 * onep * x * b1",
        "the note says the code states the γ = 0 collapse F_TLL_L = (2 "
        "x b1 − b2)/x; the block is cosyn_unpolarized_sfs' f_l, the "
        "expression that collapses to it"),
    ("SPIN32_FINITE_GAMMA.md", 810, "src/core/xsec.cpp:108-118"): (
        "std::fabs(j - 1.5) < 1e-9",
        "the note says `tables` dispatches to the _32 slots on "
        "ion().spin == 1.5; the block is that branch"),
    ("SPIN32_FINITE_GAMMA.md", 811, "src/core/xsec.cpp:119-122"): (
        "// b3, b4 are filled for EVERY spin (xsec.py fills them "
        "outside the rank-2",
        "the note says tables() fills b3, b4 for every spin; the block "
        "is those two fills and the comment that says so"),
    ("SPIN32_FINITE_GAMMA.md", 905, "src/core/xsec.cpp:302-305"): (
        "out.w_avg = out.w_avg + helicity * v * ct * a_parallel(t, x, "
        "q2, y);",
        "the note compares the proposed rank-3 term with the vector "
        "term the code already computes; the block is that term"),
    ("SPIN32_FINITE_GAMMA.md", 916, "include/lipolgen/xsec.hpp:13-21"): (
        "///   W = 1 + w_avg + a_1 cos(phi') + a_2 cos(2 phi'),   phi' "
        "= phi - phi_S",
        "the note gives the J = 3/2 cross section in the generator's "
        "own variables and normalisation; the block is that master "
        "formula as the header states it"),
    ("SPIN32_FINITE_GAMMA.md", 916, "src/core/xsec.cpp:327-335"): (
        "return dsigma_unpol(x, q2, s) / (2.0 * kPi) * std::max(w, "
        "0.0);",
        "the note's master formula (46) is dσ/(dx dQ² dφ) = "
        "[σ_U/2π]·W(φ'); the block is InclusiveKernel::dsigma, which is "
        "that line for line -- re-pointed in run 2026-09-06 from "
        ":283-285, which the P8 split had left on the signature of "
        "amplitudes()"),
    ("SPIN32_FINITE_GAMMA.md", 976, "src/core/xsec.cpp:108-118"): (
        "df = &delta_32_func_;",
        "the row says g1_rank3/g2_rank3 would be filled in the same "
        "branch as the _32 slots; the block is that branch"),
    ("SPIN32_FINITE_GAMMA.md", 981, "src/core/sampler.cpp:588-591"): (
        "out[i * nk + k] += p_m * (base + (st.*a1)[c] * std::cos(phip) +",
        "the note's change list says the per-category φ density "
        "evaluates 1 + w_avg + a1 cos φ′ + a2 cos 2φ′ by hand; the "
        "block is that hand evaluation, in the code's own spelling, "
        "which the 2026-09-16 fold into `mix_weights` moved to "
        "pointer-to-member form and split over two statements: the "
        "1 + w_avg is `1.0 + (st.*w)[c]` on the `base` line, the two "
        "harmonics `(st.*a1)[c]` / `(st.*a2)[c]` on the accumulate "
        "line -- re-pointed in run 2026-09-06 from :480-482, which the "
        "file had moved off, and widened from :590-592 (which had lost "
        "the 1 + w_avg line and ended on a bare brace)"),
    ("SPIN32_FINITE_GAMMA.md", 1009, "include/lipolgen/xsec.hpp:183-186"): (
        "(DEFAULT TRUE since 2026-08-29",
        "the note cites the rationale for `target_mass` defaulting to "
        "true; the block is that rationale"),
    ("SPIN32_FINITE_GAMMA.md", 1011, "src/core/xsec.cpp:124-131"): (
        "t.g2 = (g2_mode_ == G2Mode::kWandzuraWilczek)",
        "the note cites this for the Wandzura-Wilczek g2 table that "
        "a_parallel_exact reads; the block is where tables() fills that "
        "g2"),
    ("SPIN32_FINITE_GAMMA.md", 1024, "include/lipolgen/xsec.hpp:188-191"): (
        "chooses the g2 model",
        "the note asks for a g2_rank3_scale mirroring g2_scale; the "
        "block is the g2_mode/g2_scale documentation it mirrors"),
    ("SPIN32_FINITE_GAMMA.md", 1039, "tests/test_xsec.cpp:223-240"): (
        "terms of numerator and denominator combine into",
        "the test inventory names this test for the A_zz identity at "
        "every y; the block is that TEST_CASE"),
    ("SPIN32_FINITE_GAMMA.md", 1041, "tests/test_xsec.cpp:257-278"): (
        "TEST_CASE(\"the kernel thirds combination carries the same sign "
        "as A_zz\") {",
        "the note's test list describes this as 'the thirds combination "
        "carries the sign of A_zz'; the block is the test case of that "
        "name"),
    ("SPIN32_FINITE_GAMMA.md", 1042, "tests/test_xsec.cpp:170-193"): (
        "the population-averaged cross section is the unpolarized one",
        "the test inventory names the population sum-rule test; the "
        "block is that TEST_CASE"),
    ("SPIN32_FINITE_GAMMA.md", 1043, "tests/test_xsec.cpp:297-316"): (
        "{LI6(), {1.0, 0.0, -1.0}}, {LI7(), {1.5, 0.5, -0.5, -1.5}}",
        "the proposed test list cites this test for the one-geometry "
        "statement; the block is that TEST_CASE, which runs both ions"),
    ("SPIN32_FINITE_GAMMA.md", 1043, "tests/test_xsec.cpp:318-335"): (
        "0.5 * azz(t.b1, t.f1, t.f2, x, q2 / (s * x), &t.b2)",
        "the proposed test list cites this test for the J = 1 HJM "
        "transcription; the block is that TEST_CASE"),
    ("SPIN32_FINITE_GAMMA.md", 1044, "tests/test_xsec.cpp:337-360"): (
        "TEST_CASE(\"the spin-3/2 rate and cos-2phi channels are "
        "mutually consistent\") {",
        "the note's test list cites this for the spin-3/2 rate/cos-2φ "
        "consistency; the block is that test case"),
    ("SPIN32_FINITE_GAMMA.md", 1048, "tests/test_spin.cpp:132-145"): (
        "TEST_CASE(\"spin-3/2 populations round trip through the "
        "moments\") {",
        "the note's test list cites this as the populations/moments "
        "round trip; the block is that test case"),
    ("SPIN32_FINITE_GAMMA.md", 1048, "tests/test_spin.cpp:207-232"): (
        "TEST_CASE(\"max-entropy populations for spin 3/2, anchor (0.7, "
        "0.4)\") {",
        "the note's test list cites this as the max-entropy ladder "
        "test; the block is that test case"),
    ("SPIN32_FINITE_GAMMA.md", 1067, "tests/test_xsec.cpp:297-316"): (
        "kp.tensor_moments(0.5).first == 0.0",
        "proposed test 3 is to mirror this test for the octupole "
        "pure-state weights; the block is the rank-2 TEST_CASE it "
        "mirrors"),
    ("SPIN32_FINITE_GAMMA.md", 1071, "tests/test_xsec.cpp:170-193"): (
        "CHECK_CLOSE(tot * 2 * kPi, kern.dsigma_unpol(p.x, p.q2, s), "
        "1e-9);",
        "proposed test 4 is to extend this test to LI7 with all four "
        "sectors live; the block is the test to extend"),
    ("SPIN32_FINITE_GAMMA.md", 1075, "tests/test_xsec.cpp:388-406"): (
        "TEST_CASE(\"the tensor rate follows P_2(cos theta_S), and the "
        "magic angle kills it\") {",
        "the note says the rank-2 magic angle (P₂ = 0 at θ_S = 54.74°) "
        "is already tested here; the block is that test case, whose "
        "`magic` is acos(1/√3)"),
    ("SPIN32_FINITE_GAMMA.md", 1082, "tests/test_xsec.cpp:649-683"): (
        "amplitudes refuse a spin state of the wrong J",
        "proposed test 6 is a spin-1 rejection in the spirit of this "
        "test; the block is that TEST_CASE"),
    ("SPIN32_FINITE_GAMMA.md", 1100, "tests/test_xsec.cpp:223-240"): (
        "for (const Point& p : kConventionPoints)",
        "the note says the J = 3/2 contrast must mirror this spin-1 "
        "test, which must keep pinning -(2/3) b1/F1 unchanged; the "
        "block is that TEST_CASE"),
    ("SPIN32_FINITE_GAMMA.md", 1168, "include/lipolgen/asymmetries.hpp:9-21"): (
        "/// Spin-1 master formula (unpolarized e, target spin at angle "
        "theta_m,",
        "cited from the [HJM89] bibliography entry as where the spin-1 "
        "b1..b4 basis and c_m = 3m² − 2 are transcribed; the block is "
        "that master formula with its c_m line"),
}



# B2 (rule B, above): unnamed point citations whose sentence quotes nothing the
# cited LINE carries.  Keyed on the CITATION -- (document, citing line,
# `path:line`) -- not on the target: two documents, or two sentences in one
# document, can cite the same line for different claims and each claim is its
# own judgement (`src/core/beams.cpp:89` is cited twice by
# docs/PHYSICS_CHANNELS.md, once for the Schellingerhout product and once from
# the bibliography).  Each value is (pin, reason), with the same contract
# RANGE_ALLOW's pin has: the pin must be ON the cited line and NOWHERE ELSE in
# that file, so the exemption names the line by its content and cannot slide.
# The citing line is part of the key deliberately: an exemption is a record
# that a human read THIS sentence against that line, and a document edit that
# moves the sentence retires it, loudly ("no citation reaches it"), instead of
# quietly covering a sentence nobody read.
UNNAMED_ALLOW: dict[tuple[str, int, str], tuple[str, str]] = {
    ("PHYSICS_CHANNELS.md", 86, "docs/open_items/physics_literature.md:95"): (
        "Bissey et al. Eq. (2) reads, verbatim,",
        "the row contrasts the code's per-nucleon 3He numbers with the "
        "register's whole-nucleus ones; the cited line is the convention "
        "sentence that gives 3He as P_n = 0.86, P_p = -0.028"),
    ("PHYSICS_CHANNELS.md", 86, "docs/open_items/physics_literature.md:110"): (
        "JLab PR12-14-001",
        "the row's claim is that the 0.866/-0.037 in `LI7()` are the "
        "higher-precision Argonne online-table values as quoted by JLab "
        "PR12-14-001; the cited line is the Wiringa 2014 row of the "
        "register's reference table, which is where the register says that"),
    ("PHYSICS_CHANNELS.md", 86, "src/core/beams.cpp:89"): (
        "0.86995 x 0.9325 = 0.81123",
        "the row says [Schell93] is named beside the cluster product; the "
        "cited line IS that product, with the Schellingerhout reference on "
        "it"),
    ("PHYSICS_CHANNELS.md", 88, "include/lipolgen/beams.hpp:13"): (
        "Au (110 GeV/u)",
        "the row says the Au = 110 GeV/u number is only quoted in this "
        "header comment (Au is not a supported species); this is the line "
        "that quotes it"),
    ("PHYSICS_CHANNELS.md", 90, "README.md:128"): (
        "⁷Li ⟨P₂⟩ = −T/5",
        "the row says the identity is stated in README.md; the cited line "
        "is the external-anchors bullet that states it"),
    ("PHYSICS_CHANNELS.md", 90, "docs/DEVELOPMENT_PLAN.md:92"): (
        "⁷Li P₂ = −T/5 polarimeter",
        "the row says the identity is stated in the development plan; the "
        "cited line is the P4 row that states it"),
    ("PHYSICS_CHANNELS.md", 90, "tests/test_pipeline.cpp:608"): (
        "pipeline: 7Li generated events give <P2(cos theta_k)> = -T/5",
        "the row calls the -T/5 identity test-pinned and cites the pipeline "
        "test; the cited line is that TEST_CASE"),
    ("PHYSICS_CHANNELS.md", 125, "docs/open_items/physics_literature.md:28"): (
        "What is genuinely absent is the finite-γ",
        "the row's claim is that the finite-gamma J = 3/2 decomposition "
        "does not exist in the literature; the cited line is the register's "
        "`Recommended path` sentence that says it is genuinely absent"),
    ("PHYSICS_CHANNELS.md", 176, "src/core/pipeline.cpp:1116"): (
        "kernel && is_tagged(channel)",
        "the row says validate() refuses a caller-supplied kernel on any "
        "tagged channel rather than record a backend that never ran; the "
        "cited line is that guard, whose throw runs to :1120 and carries "
        "the \"caller-supplied kernel\" meta label the row quotes at :1118 "
        "(re-pointed 2026-09-15 from :1016, which is the unrelated unpol_sf "
        "= toy refusal)"),
    ("PHYSICS_CHANNELS.md", 204, "docs/open_items/physics_literature.md:113"): (
        "reduced-matrix-element factor",
        "the row notes that the register reads CBT Eqs. (26)/(27) as the K "
        "= 1 reduced-matrix-element factor, against the row's own reading; "
        "the cited line is the CBT06 register row that reads them that way"),
    ("PHYSICS_CHANNELS.md", 209, "include/lipolgen/constants.hpp:139"): (
        "Sather-Schmidt bag-model sum-rule coefficient",
        "the row's [SS90] provenance for C_BAG = -0.012; the cited line is "
        "the comment naming Sather-Schmidt PRD 42:1424, with the constant "
        "itself two lines below (re-pointed 2026-09-15 from :86, which is "
        "the b1-normalisation comment and has been since before bd775bc)"),
    ("PHYSICS_CHANNELS.md", 238, "src/core/rng.cpp:2"): (
        "xoshiro256** seeded via splitmix64",
        "the row says xoshiro256**/splitmix64 are named at this line only "
        "and that no Blackman-Vigna citation exists in the repo; this is "
        "that line"),
    ("PHYSICS_CHANNELS.md", 276, "src/core/spectator.cpp:300"): (
        "Yellow Report Table 10.1, HADRON beam",
        "the row's [YR] provenance for the hadron-beam divergences; the "
        "cited line is the Table 10.1 transcription comment above the "
        "numbers"),
    ("PHYSICS_CHANNELS.md", 304, "docs/CONVENTIONS.md:415"): (
        "a_nn = −18.9 fm",
        "the row says a_nn = -18.9 fm has no primary reference and appears "
        "only as an in-code constant and in CONVENTIONS.md; the cited line "
        "is where CONVENTIONS.md carries it (re-pointed 2026-09-15 from "
        ":357, which the file's growth since bd775bc turned into the "
        "deuteron-control sentence)"),
    ("PHYSICS_CHANNELS.md", 342, "docs/open_items/physics_literature.md:245"): (
        "keep the coherent-cluster amplitude as the systematic variant",
        "the row says the register calls the per-nucleon Glauber product "
        "the citable primary and the coherent-cluster amplitude the "
        "systematic variant, the opposite of the code's default; this is "
        "the register line that says it (re-pointed 2026-09-15 from the "
        "range :81-85, which is the W_FSI formula and its starting "
        "parameters -- the other citation of that range, for the 30-70 mb "
        "Deeps fit, is right and stays)"),
    ("PHYSICS_CHANNELS.md", 420, "tests/test_rc.cpp:518"): (
        "the first C0 zero lies in [2.9, 3.3] fm^-1",
        "the row says T11 gates the whole refit window q0 in [2.9, 3.3] "
        "fm^-1; the cited line is that SUBCASE of T11 (re-pointed "
        "2026-09-15 from :496, a comment inside the <r^2>_point subcase, "
        "which gates a different number)"),
    ("PHYSICS_CHANNELS.md", 996, "include/lipolgen/spin.hpp:24"): (
        "(Varshalovich) conventions",
        "the [Varsh] bibliography entry says the Wigner-d / Clebsch-Gordan "
        "conventions are named at this line and nowhere else in the "
        "repository; this is that line"),
    ("PHYSICS_CHANNELS.md", 1017, "src/core/beams.cpp:89"): (
        "(Schellingerhout PRC 48:2714)",
        "the [Schell93] bibliography entry says the reference is named "
        "beside the 6Li cluster product; the cited line is the product "
        "carrying it"),
    ("SPIN32_FINITE_GAMMA.md", 142, "include/lipolgen/spin.hpp:23"): (
        "Hoodbhoy-Jaffe-Manohar NPB 312:571",
        "the note's conventions paragraph cites the HJM c_m = 3m^2 - 2 "
        "convention; the cited line is the header line that states it"),
    ("SPIN32_FINITE_GAMMA.md", 151, "examples/generate_tagged.cpp:123"): (
        "1.0, j >= 1.0 ? 1.0 : 0.0);",
        "the paragraph says this line RECORDED `pzz_true = 1.0` before "
        "commit 0961c60 and that both CLIs now write `j >= 1.0 ? 1.0 : "
        "0.0`; the cited line is that line as it stands after the fix, "
        "which is what the paragraph's next sentence quotes"),
    ("SPIN32_FINITE_GAMMA.md", 743, "tests/test_tensor_gamma.cpp:382"): (
        "tensor_gamma: Table 1 row 2, a target polarized along q",
        "the note says the frame choice is pinned by the paper's own Table "
        "1 rows; the cited line is that TEST_CASE (its sibling `:406` is "
        "the other row)"),
    ("SPIN32_FINITE_GAMMA.md", 771, "src/core/xsec.cpp:116"): (
        "2.0 * x * t.b1",
        "the note's claim is that the generator's default is b2 = 2*x*b1 in "
        "the per-nucleon normalisation; the cited line is that default "
        "branch (the note writes the product with a middle dot, which is "
        "why the quoted phrase does not match the code)"),
    ("SPIN32_FINITE_GAMMA.md", 791, "tests/test_tensor_gamma.cpp:152"): (
        "reduce to the HJM b-sector at gamma = 0",
        "the note says this test pins the gamma -> 0 collapse explicitly "
        "`for any b2`; the cited line is the TEST_CASE, and the `for any "
        "b2` it quotes is in the test's own comment at :155"),
    ("SPIN32_FINITE_GAMMA.md", 905, "include/lipolgen/xsec.hpp:19"): (
        "lam_e P_e (m/J) cos(theta_S) * A_par(x,y)",
        "the note compares the proposed rank-3 term with the vector term "
        "the code already computes; the cited line is that term in the "
        "header's master formula (the note writes it in Unicode)"),
    ("SPIN32_FINITE_GAMMA.md", 948, "src/core/asymmetries.cpp:85"): (
        "TENSOR_LL_SIGN * 2.0 / 3.0",
        "the note says the spin-1 2/3 is an explicit factor in azz(); the "
        "cited line is where azz() multiplies it in (the note quotes "
        "`2.0/3.0`, the code spaces it)"),
    ("SPIN32_FINITE_GAMMA.md", 983, "src/core/bookkeeping.cpp:89"): (
        "leaves the octupole moment R_3 at its default 0",
        "the row asks for an explicit `o` beside `pzz` because the explicit "
        "branch silently sets R_3 = 0 today; the cited line is the comment "
        "that states exactly that (the call itself is at :94, which is what "
        "sec. 2.4 cites)"),
    ("SPIN32_FINITE_GAMMA.md", 1064, "tests/test_tensor_gamma.cpp:152"): (
        "reduce to the HJM b-sector at gamma = 0",
        "test 2 of the proposed suite is to mirror this test at "
        "`target_mass = true`; the cited line is that TEST_CASE"),
    ("SPIN32_FINITE_GAMMA.md", 1167, "include/lipolgen/spin.hpp:23"): (
        "Hoodbhoy-Jaffe-Manohar NPB 312:571",
        "the [HJM89] bibliography entry says the spin-1 basis and the c_m = "
        "3m^2 - 2 convention are carried at this line; this is that line"),
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


# ------------------------------------------------------------- R5: recording evidence
#
# What a range citation has to SAY about its block before the sidecar will
# fingerprint it.  See R5 in the module docstring for the rule and for what it
# measured over the three covered documents.

# `X:12`, `X:12-20`, `BeamRemnants.cc:662,935` -- a citation, never evidence.
CITELIKE = re.compile(r"^\S*:\d+(?:[-,]\d+)*$")


# A markdown LIST ITEM: `- `, `* `, `+ `, `1. `, `2) `.
LIST_ITEM = re.compile(r"^[ \t]*(?:[-*+]|\d+[.)])[ \t]")
# Whether a list item is its own paragraph.  A constant so that
# `validation/record_rule_sweep.py` can measure the rule with it off; the gate
# never turns it off.
LIST_ITEM_PARAGRAPH = True


def paragraph(text: str, start: int, end: int) -> tuple[int, int]:
    """The enclosing paragraph of a citation.  A markdown TABLE ROW is its own
    paragraph -- one row is one citing sentence (its claim cell and its
    provenance cell belong to each other, and the row above belongs to
    neither) -- and since 2026-09-15 so is a markdown LIST ITEM, for the same
    reason: a bulleted list carries no blank lines, so the whole list was one
    "sentence" and a phrase quoted in the NEXT bullet evidenced this bullet's
    block.  That is not a hypothetical (D3.5): `docs/PHYSICS_CHANNELS.md:974`
    cites `docs/USAGE.md:1501-1509` for the four C++ example generators and
    was evidenced by the word `tier`, which belongs to the bullet BELOW it --
    the block is the T1 tier section and the citation is wrong.  Clipping to
    the item costs exactly ONE citation over the three covered documents of
    a3c9ecb -- and that citation is the wrong one (recipe:
    `validation/record_rule_sweep.py r5 --tree <a3c9ecb checkout with this
    file copied in>`, rows one and two).  Anything else is the
    blank-line-delimited paragraph."""
    ls = text.rfind("\n", 0, start) + 1
    le = text.find("\n", end)
    le = len(text) if le < 0 else le
    if text[ls:le].lstrip().startswith("|"):
        return ls, le
    lo = text.rfind("\n\n", 0, start)
    lo = 0 if lo < 0 else lo + 2
    hi = text.find("\n\n", end)
    hi = len(text) if hi < 0 else hi
    if not LIST_ITEM_PARAGRAPH:       # off only for the D3.5 measurement
        return lo, hi
    item = None                       # the enclosing list item, if any
    i = ls
    while i >= lo:
        ie = text.find("\n", i)
        ie = len(text) if ie < 0 else ie
        if LIST_ITEM.match(text[i:ie]):
            item = i
            break
        if i == lo:
            break
        i = text.rfind("\n", 0, i - 1) + 1
        if i < lo:
            break
    if item is None:
        return lo, hi
    nxt, j = hi, text.find("\n", end)
    while j != -1 and j < hi:
        je = text.find("\n", j + 1)
        je = len(text) if je < 0 else je
        if LIST_ITEM.match(text[j + 1:je]):
            nxt = j
            break
        j = je if je < hi else -1
    return max(lo, item), min(hi, nxt)


def citing_sentence(text: str, start: int, end: int, win: int | None = None) -> str:
    """`win` characters either side of the citation, clipped to its paragraph.
    This is what gets PRINTED; the phrases are read off the paragraph itself
    (`quoted_phrases`), because a window cut through the middle of a code span
    pairs the surviving backtick with the next one and reads a run of prose as
    a quotation."""
    win = RECORD_WIN if win is None else win
    lo, hi = paragraph(text, start, end)
    return text[max(lo, start - win):min(hi, end + win)]


def quoted_phrases(text: str, start: int, end: int,
                   win: int | None = None) -> list[str]:
    """Every quoted phrase within `win` characters of the citation, in order,
    minus the citations themselves.  Scanned over the whole paragraph and then
    filtered by distance, so that backticks always pair as the document wrote
    them and only whole phrases are considered."""
    win = RECORD_WIN if win is None else win
    lo, hi = paragraph(text, start, end)
    near_lo, near_hi = start - win, end + win
    out: list[str] = []
    for m in QUOTED.finditer(text, lo, hi):
        if m.end() < near_lo or m.start() > near_hi:
            continue
        q = (m.group(1) or m.group(2) or m.group(3) or "").strip()
        if not q or BLANK_LINE.search(q) or CITELIKE.match(q):
            continue
        q = " ".join(q.split())
        if q not in out:
            out.append(q)
    return out


def symbol_bases(phrase: str) -> list[str]:
    """The identifiers a phrase names, if it names any: `Optics::lumi_fraction`
    -> [`lumi_fraction`], `spin1_populations(p)` -> [`spin1_populations`],
    `--pzz` -> [`--pzz`].  A one-character name (`M`, `x`) is not taken as a
    symbol: it matches everywhere and evidences nothing.  A slash-separated
    run of names -- `HelicityFlipOptions::theta_s/phi_s`, `pe/pz/pzz`, the
    documents' way of writing "these members" -- is every one of them, but
    only when EVERY part is a name: `data/vmc/README.md` is a path, not three
    symbols."""
    head = re.split(r"[(<\[{ ]", phrase.strip(), 1)[0]
    head = head.split("::")[-1].strip().rstrip(".,;:!?")
    parts = head.split("/")
    if len(parts) > 4:
        return []
    out = [q.strip() for q in parts]
    if all((FLAG.match(q) or IDENT.match(q)) and len(q) >= MIN_SYMBOL
           for q in out) and out:
        return out
    return []


def flat(lines: list[str]) -> str:
    return " ".join(" ".join(lines).split())


def text_evidence_ok(src: "Source", a: int, b: int, want: str) -> bool:
    """Whether a TEXT match is worth anything, MEASURED (D3.5, recipe
    `validation/record_rule_sweep.py r5`).

    Two limits, both of them the same statement -- a phrase evidences a block
    only if finding it there says something about THAT block:

      * a block of more than `EVIDENCE_MAX_BLOCK` lines is too big for a text
        match to mean anything.  Measured over the three covered documents of
        the committed tree a3c9ecb by the decoy test (the same citation
        pointed at a same-length block shifted clear of the true one; recipe
        `validation/record_rule_sweep.py r5 --tree <a3c9ecb checkout>`), text
        evidence accepts a decoy 11.2 % of the time for blocks of 1-20 lines,
        18.6 % for 20-60, 38.5 % for 60-150 and 50.0 % for 150+: in a 371-line
        block the rule is a coin toss.  150 is where it becomes one; the two
        blocks between 60 and 150 were read and each rests on a phrase
        specific to it, so the threshold is not taken lower -- 60 costs two
        more exemptions and buys nothing measurable (the last row of that
        table).
      * a ONE-WORD phrase must occur on at most `EVIDENCE_WORD_LINES` lines of
        the target FILE.  `tier` is on 17 lines of `docs/USAGE.md` and `pzz`
        on 39 of `python/bindings.cpp`: a word that common cannot tell one
        block of that file from another, which is exactly how the 12 blocks
        of D3.5 came to rest on `ratio`, `rates`, `applies`, `seed`, `tier`,
        `bound`, `pzz`.  A multi-word phrase (`--x-max 0.95`, a formula, a
        sentence) is not restricted: it is already specific."""
    if EVIDENCE_MAX_BLOCK is not None and b - a + 1 > EVIDENCE_MAX_BLOCK:
        return False
    if EVIDENCE_WORD_LINES is not None and " " not in want:
        if sum(1 for x in src.lines if want in x) > EVIDENCE_WORD_LINES:
            return False
    return True


def phrase_in_block(src: "Source", a: int, b: int, phrase: str):
    """How a quoted phrase is evidence for the block: as TEXT (a
    whitespace-normalised substring, subject to `text_evidence_ok`) or as a
    SYMBOL (declared or used in it)."""
    want = " ".join(phrase.split())
    body = flat(src.lines[a - 1:b])
    if len(want) >= MIN_PHRASE and text_evidence_ok(src, a, b, want):
        # A bare identifier has to match on WORD boundaries -- `rates` inside
        # "generates" is not the block naming it.  Anything longer (a formula,
        # a sentence) is matched as written.
        if IDENT.match(want):
            if word(want).search(body):
                return "text"
        elif want in body:
            return "text"
    for base in symbol_bases(phrase):
        if any(src.declares(i, base) or src.in_code(i, base)
               for i in range(a, b + 1)):
            return "symbol"
    return None


def range_evidence(text: str, m: re.Match, src: "Source", a: int, b: int,
                   name: str | None):
    """R5 branches (a) and (b) on ONE range citation.

    Returns `(how, detail)`: `how` is "name", "text" or "symbol" when the
    citing sentence evidences the block, and None when it does not, `detail`
    saying either what carried it or what the rule looked for and missed."""
    if name and not CITELIKE.match(name):
        bases = symbol_bases(name)
        if bases:                                                        # (a)
            if any(src.declares(i, base) or src.in_code(i, base)
                   for base in bases for i in range(a, b + 1)):
                return "name", f"the adjacent name `{name}` is in the block"
            return None, (f"the adjacent name `{name}` is neither declared nor "
                          f"used anywhere inside the block")
        # a name that is not an identifier (`pe=0.7`) is a quoted phrase like
        # any other, and is picked up by (b) below
    phrases = quoted_phrases(text, m.start(), m.end())
    for q in phrases:                                                    # (b)
        how = phrase_in_block(src, a, b, q)
        if how:
            return how, f"`{clip(q, 60)}` occurs in the block as {how}"
    if phrases:
        shown = ", ".join("`%s`" % clip(q, 40) for q in phrases[:6])
        more = f" and {len(phrases) - 6} more" if len(phrases) > 6 else ""
        return None, (f"the citing sentence quotes {shown}{more}; none of them "
                      f"is in the block")
    return None, ("the citing sentence carries no adjacent name and quotes "
                  "nothing, so there is nothing to check the block against: "
                  + clip(" ".join(citing_sentence(text, m.start(),
                                                  m.end()).split()), 140))


def pin_state(src: "Source", a: int, b: int, pin: str):
    """None when a RANGE_ALLOW pin still names this block, else why it does
    not.  The pin has to be INSIDE the block and nowhere else in the file: an
    exemption that matches two places in a file names neither."""
    whole = "\n".join(src.lines)
    if pin not in "\n".join(src.lines[a - 1:b]):
        return "is not inside the block"
    n = whole.count(pin)
    if n != 1:
        return f"occurs {n} times in the file, so it does not name this block"
    return None


# --------------------------------------------------------------- extra documents
#
# 2026-09-05 (Phase E, item E1).  `docs/PHYSICS_CHANNELS.md` is the document
# this gate was built for, but it is not the only one carrying `path:line`
# citations that drift the exact same way.  Measured that day, counting with
# this file's own REF/RANGE/EXT regexes: `docs/USAGE.md`, `docs/CONVENTIONS.md`
# and `docs/T2_CHAIN.md` carry NONE; `docs/OPEN_ITEMS_SOLUTIONS.md` carries 3
# (1 range, 2 external) -- not "more than a handful", so it stays uncovered and
# the count is recorded instead (`docs/open_items/run_2026-09-03/
# phase_E_numbers.md` §E1); `docs/theory/SPIN32_FINITE_GAMMA.md` carries 134
# (19 named points, 115 ranges) and `docs/PYTHIA_BRIDGE.md` carries 8 external
# PYTHIA upstream pointers -- both extended here, each with its own
# (possibly empty) allow-list, and both fixed to 0 broken on that date: 15
# stale ranges in the theory note (drift since the Phase E code comments
# shifted several of them -- exactly what that document's own top-matter says
# happened, and until now nothing machine-checked it) and one PYTHIA_BRIDGE.md
# citation pair (`BeamRemnants.cc:662`, `:935`) that had been written as two
# separate citations with no inheritance rule for an external file name, so it
# was checked by nothing at all; rewritten as the single comma-joined
# `BeamRemnants.cc:662,935` docs/PHYSICS_CHANNELS.md already uses.
#
# A document in this list absent under `ROOT` is skipped rather than reported
# missing: only the PRIMARY document (`DOC`) is required to exist, so a test
# fixture tree that creates nothing but `docs/PHYSICS_CHANNELS.md` still
# checks exactly what it always did.  Every range/use-site fingerprint below
# lives in the SAME sidecar as the primary document's (`RANGES_REL`) -- a
# block or upstream line cited by more than one document hashes to one entry
# either way, which is the correct reading: the fact is the source line, not
# which markdown file points at it.
ALLOW_SPIN32: dict = {}
ALLOW_PYTHIA_BRIDGE: dict = {}
EXTRA_DOCS = [
    ("docs/theory/SPIN32_FINITE_GAMMA.md", ALLOW_SPIN32),
    ("docs/PYTHIA_BRIDGE.md", ALLOW_PYTHIA_BRIDGE),
]


def check_document(doc: Path, allow: dict, *, fix: bool, record: bool,
                   strict: bool, cache: dict, recorded: dict, fresh: dict,
                   seen: set, remap: dict, ext_dir, evidence: dict | None = None,
                   use_allow: bool = True, unnamed_seen: set | None = None):
    """Check one document's citations against the live tree.

    `cache` (path -> Source|None), `recorded` (the loaded sidecar, read-only
    here), `fresh` (sha256 entries this run actually reached), `seen` (keys
    reached this run), `remap` (moved-key old -> new, `--fix` only) and
    `evidence` (R5's verdict per range key, when recording or auditing) are
    shared across every document `main()` processes in one invocation, since a
    fingerprint, an R5 verdict or a "no longer cited anywhere" verdict is a
    property of the whole run, not of one document: a block cited by two
    documents is evidenced if EITHER sentence evidences it, and `evidence`
    keeps the first verdict that passes, or the first that fails when none
    does.  Returns `(bad, new_text)`: `bad` is the
    list of broken citations (empty means this document is clean) and
    `new_text` is the possibly-rewritten document text (`None` when `--fix`
    changed nothing), for the caller to write back.  Prints this document's
    own report (mirroring the single-document gate's own tail) before
    returning, so callers only need to aggregate the exit code and the
    cross-document housekeeping (dead fingerprints, `--record-ranges`).
    """
    text = doc.read_text()
    bad: list[str] = []
    fixed: list[str] = []
    allowed: list[str] = []
    skipped: list[str] = []
    n = nr = nx = 0

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
            # R5.  A fingerprint is only worth taking of the RIGHT block, and
            # the citing sentence is the only thing that says which that is --
            # THIS citing sentence.  Until 2026-09-15 the verdict was pooled
            # per BLOCK ("a block cited by two documents is evidenced if
            # EITHER sentence evidences it") and `RANGE_ALLOW` was keyed on
            # the range, so a citation re-pointed onto an already-blessed
            # block was accepted whatever its own sentence said: three of five
            # deliberately wrong-block ranges passed that way, and 11
            # citations on 7 keys in the committed documents had never been
            # judged against their own sentence at all (D3.2).  The verdict is
            # now per CITATION -- (document, citing line, range) -- and a
            # block is fingerprinted only when EVERY citation of it is
            # evidenced or exempted.
            how, detail = range_evidence(text, m, src, a, b, name)
            cite = (doc.name, text.count("\n", 0, m.start()) + 1, key)
            if how is None and cite in RANGE_ALLOW and use_allow:
                pin, why = RANGE_ALLOW[cite]
                bad_pin = pin_state(src, a, b, pin)
                if bad_pin is None:
                    how, detail = "allow", why
                else:
                    detail = (f"exempted in RANGE_ALLOW, but the pin "
                              f"{clip(pin, 60)!r} {bad_pin}")
            if evidence is not None:
                evidence[cite] = (how, detail, doc.name, key)
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
            # Rule B.  B1: the line is one somebody could have meant -- in
            # bounds, not blank, not a bare delimiter.
            if not src.substantive(line):
                what = "blank" if not src.lines[line - 1].strip() else \
                       f"only punctuation ({src.lines[line - 1].strip()!r})"
                bad.append(f"{path}:{line}  unnamed reference to a line that "
                           f"is {what}")
                return m.group(0)
            if not strict:
                return m.group(0)
            # B2 (2026-09-15): THE CITING SENTENCE MUST EVIDENCE THE LINE.
            # B1 was the whole of rule B, and "non-blank" is a test nearly
            # every line of every file passes: 14 of the 19 unnamed citations
            # survived a +1 drift and 12 a -1 (RESIDUALS 1), and two had in
            # fact gone stale weeks apart with the gate green throughout --
            # `docs/PHYSICS_CHANNELS.md:301` pointing at
            # `docs/CONVENTIONS.md:357` for a_nn (now at :383; :357 is the
            # deuteron-control sentence) and `:206` pointing at
            # `include/lipolgen/constants.hpp:86` for the [SS90] bag-model
            # coefficient (`C_BAG` is at :136-138).  So an unnamed citation
            # now carries R5's evidence test, applied to the cited LINE and
            # nothing either side of it: a phrase the citing sentence QUOTES
            # must occur on that line, as text or as a symbol.  MEASURED over
            # the 42 unnamed citations the committed tree a3c9ecb carries
            # (recipe: `validation/record_rule_sweep.py unnamed --tree <a
            # checkout of a3c9ecb with this file copied in>`, D3.1): B1 alone
            # flags 0 of them and lets 76.2 % of +1 and 76.2 % of -1 drifts
            # through; this rule flags 30 and lets 2.4 % and 4.8 %; the R6
            # shared-token rule, the obvious alternative, flags 19 and lets
            # 16.7 % and 21.4 %.  The quoted-phrase rule is the one taken.  A
            # citation it flags is not necessarily wrong -- it is "this gate
            # cannot tell" again -- and the answer is the same as R5's: read
            # the row and sign an exemption with a PIN in UNNAMED_ALLOW.
            how, detail = range_evidence(text, m, src, line, line, None)
            cite = (doc.name, text.count("\n", 0, m.start()) + 1,
                    f"{path}:{line}")
            if unnamed_seen is not None:
                unnamed_seen.add(cite)
            if how is not None:
                return m.group(0)
            if use_allow and cite in UNNAMED_ALLOW:
                pin, why = UNNAMED_ALLOW[cite]
                bad_pin = pin_state(src, line, line, pin)
                if bad_pin is None:
                    allowed.append(f"{path}:{line}  (unnamed) -- {why}")
                    return m.group(0)
                bad.append(f"{path}:{line}  unnamed citation exempted in "
                           f"UNNAMED_ALLOW, but the pin {clip(pin, 60)!r} "
                           f"{bad_pin}")
                return m.group(0)
            bad.append(f"{path}:{line}  unnamed citation, and {detail}")
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
        key = (path, line, base) if (path, line, base) in allow else \
              (path, line, name) if (path, line, name) in allow else None
        if key is not None:
            pin, why = allow[key]
            if pin in src.lines[line - 1]:
                allowed.append(f"{path}:{line}  `{name}` -- {why}")
                return m.group(0)
            bad.append(f"{path}:{line}  `{name}` is allow-listed on this line "
                       f"but the line no longer contains its pin {pin!r}; the "
                       f"block has moved, or the code changed")
            return m.group(0)
        if any(p == path and nm in (base, name) for p, _, nm in allow):
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

    mode = "strict" if strict else "loose"
    label = "" if doc == DOC else f"  [{doc.name}]"
    print(f"{n} references checked ({mode}), {nr} ranges, {nx} external, "
          f"{len(bad)} broken"
          + (f", {len(allowed)} allow-listed" if strict else "")
          + (f", {len(skipped)} skipped" if skipped else "")
          + (f", {len(fixed)} fixed" if fix else "")
          + label)
    for x in fixed:
        print("  fixed", x)
    for a in allowed:                     # an exemption you cannot see is a hole
        print("  allowed", a)
    for sk in skipped:                    # nor one you cannot see
        print("  skipped", sk)
    for b in bad:
        print("  ", b)
    return bad, (new_text if fix and fixed else None)


# ------------------------------------------------------- R6: the dated run records
#
# `docs/open_items/` holds the DATED RUN RECORDS -- run_2026-09-02/,
# run_2026-09-03/, run_2026-09-06/ -- and the older standing notes beside them
# (`engineering.md`, `physics_literature.md`, `vmc_*.md`, `code_designs.md`).
# Until 2026-09-15 (Phase D item D2) nothing gated them at all, and they carry
# 503 `path:line` citations into the live tree that drift with every line shift.
# Fourteen were re-pointed by hand at the 2026-09-06 close-out; nothing said
# whether the rest were right.
#
# They are HISTORY, so the rules above are the wrong rules for them.  A record
# says what was true on its date; R4's fingerprints would freeze a block that
# the record never claimed was frozen, R5's evidence test is a RECORDING rule
# and there is nothing here to record, and S1/S2/S3's declaration tests would
# refuse the legitimate citation of a use site or of a line of prose in another
# document.  So R6 is one RELAXED rule, and it asks only what a reader
# following the anchor would ask:
#
#   R6  the target must be a line somebody could have meant, and it must have
#       something to do with the sentence that cites it.
#
#       (i)   IN BOUNDS and NON-BLANK.  A point citation's line must carry
#             something; a block must have at least one non-blank line in it.
#             (This is deliberately weaker than R2, which demands a non-blank
#             first AND last line: R2 protects a fingerprint, and there is no
#             fingerprint here.)
#       (ii)  a SHARED TOKEN, when the target is a LIVE file -- anything
#             outside `docs/open_items/`.  At least one token of
#             `RELAX_TOKEN` = 4 characters or more (a symbol, a number, a
#             word) must occur BOTH in the citing sentence and in the target
#             line or block.  The comparison is case-insensitive and the
#             citations themselves are blanked out of the sentence first, so
#             `src/core/xsec.cpp` cannot evidence `src/core/xsec.cpp`.
#
#             A POINT citation is held to more than that, since 2026-09-15:
#             it names ONE line, so its shared token must be `RELAX_POINT_TOKEN`
#             = 6 characters or more (or 5 with a digit in it -- `1.848` says
#             more than `which`) AND must be absent from the two lines either
#             side of the cited one.  A rule that the neighbouring lines
#             satisfy too cannot see the drift it exists to catch: on the tree
#             this phase leaves the one-token rule accepts 29.6 % of +1 drifts
#             of these citations and this one 12.2 % (D3.3, thirty candidates,
#             recipe `validation/record_rule_sweep.py drift --only points`).
#             A RANGE keeps the one-token rule: a block shifted by a line
#             still holds what the sentence says, and every candidate measured
#             refused a crowd of honest citations to buy little.
#       (iii) a citation into ANOTHER RECORD (a target under
#             `docs/open_items/`) gets (i) and nothing else.  One record
#             quoting another is a pointer inside the history, and the two
#             move together or not at all.
#
# Nothing here is FINGERPRINTED.  A fingerprint is a promise that a block will
# not change; a record makes no such promise about anyone else's file.
#
# THE HISTORICAL ANNOTATION.  Some of these citations point at a state that no
# longer exists -- the `was` column of a re-point table, a survey of warning
# text that has since been deleted, a design clause that was superseded.  The
# honest repair is not to move the anchor (which would make the record claim
# something it never claimed) but to say WHEN it was true.  Two forms, both
# visible in the rendered document:
#
#     `src/core/xsec.cpp:266-270` (as of a3c9ecb)     -- in the citing sentence
#     **Anchors as of `a3c9ecb`.**                    -- a section marker
#
# and, for an anchor a human FOLLOWED and found right, `(anchor read <date>)`,
# which for a POINT citation must carry the text that reader saw on the line:
#
#     `src/core/xsec.cpp:295` (anchor read 2026-09-15 "double helicity")
#     `python/bindings.cpp:1528` (anchor read 2026-09-15 `arg("emc_ratio`)
#
# in "straight quotes" or `backticks`, whichever the line does not contain.
# Without the pin the annotation exempted the citation from (ii) and left only
# (i), so the anchor was free to drift onto any non-blank line of the file --
# 28 point citations were in exactly that state.  A section marker covers
# RANGES only, for the same reason (ii) is weaker for them.
#
# The section marker covers every citation from the marker to the next
# markdown heading.  An annotated citation is checked for (i) and is exempt
# from (ii): the sentence is describing the past, so there is nothing in the
# present for it to share a token with.  The commit named is not resolved here
# -- this gate shells out to nothing -- but the suite does resolve every one of
# them (`python/tests/test_doc_link_gate.py::test_record_asof_commits_resolve`).
#
# WHAT R6 DOES NOT CHECK, and why the count of citations it checks is smaller
# than the count in the records:
#
#   * a citation whose PATH is carried by PROSE.  `:1285-1289` inside "The
#     design (`:1285-1289`) states the window as ..." inherits its path from
#     the word "design" three paragraphs up.  Two inheritance rules ARE
#     implemented -- the same-line one the main gate already has
#     (`inherited_path`) and a markdown-TABLE one (`column_path`: the last
#     backticked path earlier in the row, else the nearest earlier row of the
#     same table, which is how these records write a file column) -- and
#     between them they resolve 81 of the 136 bare citations.  The other 55
#     are prose-carried and no rule of this kind can read them; they are
#     counted and named, not silently dropped.
#   * a citation whose target does not exist under ROOT.  Nine of them, and
#     none is a LiPolGen path: `src/nucleons.cpp`, `src/inputParameters.cpp`,
#     `src/subnucleon_config.hpp`, `src/main.cpp` are eSTARlight's own sources
#     (the design notes that read them), `src/HardDiffraction.cc` and
#     `examples/main234.cc` are PYTHIA's, and `docs/note_*.md` /
#     `docs/consistency_review_2026-09-02.md` belong to the PREDECESSOR
#     project the physics notes quote.  They are listed by name on every run,
#     exactly as the PYTHIA externals are, so deleting a live file shows up
#     here as a new skipped line rather than as silence.

RECORDS_REL = "docs/open_items"
RELAX_WIN = 200                  # citing-sentence window, as R5's RECORD_WIN
RELAX_TOKEN = 4                  # a shared token is this many characters or more
# ... except for a POINT citation, which names ONE line and must therefore say
# something about that line and not about its neighbourhood: the shared token
# has to be this long AND absent from the lines either side of the cited one.
# MEASURED (D3.3, recipe `validation/record_rule_sweep.py drift --only
# points`, thirty candidate rules): over the 110 live-file POINT citations the
# tree this phase leaves carries, the D2 one-token rule lets 29.1 % of +1
# drifts and 23.6 % of -1 drifts through and accepts a decoy 16.1 % of the
# time; this rule lets 12.7 %, 8.2 % and 9.3 % (98 / 29.6 / 23.5 / 15.7 /
# 12.2 / 8.2 / 8.7 on the pre-section-D3 copy it was first measured on).  Applied to the tree D2 left,
# it flagged 27 point citations -- every one read, then re-pointed, annotated
# or pinned in D3.4, where they are listed one by one.  A NUMBER is allowed to
# be shorter because it says more: `1.848`, `0.001` and `0.05` are five
# characters that name one line of a file, and refusing them cost four honest
# citations for nothing.  The rule is NOT applied to a range: a block shifted
# by one line still contains what the sentence says it does, and no candidate
# discriminated a shifted block at a price worth paying -- the cheapest that
# gets ranges below 50 % at +1 refuses 93 of 156 honest citations (this tree;
# 86 of 145 on the copy the rule was first measured on).
RELAX_POINT_TOKEN = 6
RELAX_POINT_NUMBER = 5           # ... or this long, if it carries a digit
# A token: a word, a symbol, or a number.  Dotted runs (`0.440`, `xsec.cpp`,
# `cfg_.coherent`) are one token AND their parts, so `0.440` can match `0.440`
# and `xsec.cpp` can match `xsec`.
RELAX_TOK = re.compile(r"[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)*")
# ... except that a shared COMMON ENGLISH FUNCTION WORD is not evidence of
# anything.  Both sides of this comparison are English prose as often as they
# are code, so without this list `which`, `that` and `because` carry citations
# that have nothing to do with their target.  Measured with
# `validation/record_rule_sweep.py tokens [--no-stoplist]` (D2.4): against the
# committed a3c9ecb these 115 words turn 17 apparent passes into failures --
# every one of them a citation that shares nothing but a function word -- and
# drop the decoy-acceptance rate from 31.4 % to 26.4 %, a true-to-decoy ratio
# of 1.29 -> 1.37; against the tree this phase leaves they cost NO passes at
# all and drop decoys 39.0 % -> 34.7 % (ratio 2.56 -> 2.88).  (Under the rule
# D2 landed, before a list item became its own citing sentence, the same
# figures were 16, 32.6 -> 27.7 % and 39.1 -> 34.6 %.)  Only words with
# NO meaning in this tree are listed: `line`, `time`, `open`, `make`, `show`,
# `state`, `part`, `left` and `still` are all domain words here and are
# deliberately absent.
RELAX_STOP = frozenset("""
about after again against also although always among another because been
before being below between both cannot could does doing done during each
either else enough even ever every from further gave give given gives goes
gone have having hence here however into itself just less many more most much
must neither never none often once only other others otherwise over quite
rather same seem seems several shall should simply since some something such
than that their them themselves then there therefore these they thing things
think this those though through thus together toward under unless until upon
very were what when where whereas whether which while whose will with within
without would your yours
""".split())
# The historical annotation.  It must follow ITS OWN citation -- immediately,
# or immediately after the citation's adjacent `` `name` `` -- and nothing
# else: keyed on the citing SENTENCE it silently exempted every other citation
# within `RELAX_WIN` characters, which in these documents is most of a table
# row (one annotation was measured covering four neighbours).
ASOF_INLINE = re.compile(r"(?: `[^`\n]+`)? \(as of ([0-9a-f]{7,40})\)")
# ... and as a section marker: a line of its own, covering every citation from
# it to the next markdown heading.
ASOF_SECTION = re.compile(r"^\s*(?:>\s*)?\**Anchors as of\s+`?([0-9a-f]{7,40})`?",
                          re.MULTILINE)
# The second annotation: READ BY HAND.  Measured over these records, HALF of
# what R6 refuses is not a stale anchor at all -- it is an anchor that lands
# exactly where the sentence says while the sentence names its subject in words
# the target never uses ("the joint-sampling order", "the `--b1-unpol` help
# string"), or a citation that is the SUBJECT of the sentence rather than a
# pointer out of it (a re-point table's `was` and `is now` columns, a list of
# the sidecar keys an exemption covers).  D1 hit the same wall one rung up and
# answered it the same way: an R5 refusal is "this gate cannot tell", not "this
# citation is wrong", and the answer is an exemption a human signs.  Here the
# exemption is signed in the document itself, with the date it was read, and it
# is counted and named on every run.
# Since 2026-09-15 a POINT citation's annotation must carry a PIN -- the text
# the reader saw on that line -- written `(anchor read 2026-09-15 "…")`.
# Without one the annotation exempted the citation from (ii) and left only
# (i), so an `(anchor read)` point citation drifted onto ANY non-blank line of
# its target was accepted: 28 of them, pinned by nothing (D3.3).  A RANGE keeps
# the bare form: a block that has been read is 62 of the 90 annotations, a
# block's drift is the one R6 cannot discriminate anyway (D3.3's table), and
# the pin would be 62 more transcriptions for a measurement that does not move.
# The pin is written between "straight quotes" or between `backticks`,
# whichever the line it transcribes does not itself contain.
READ_INLINE = re.compile(
    r"(?: `[^`\n]+`)? \(anchor read (\d{4}-\d{2}-\d{2})"
    r"(?: (?:\"([^\"\n]{3,200})\"|`([^`\n]{3,200})`))?\)")
READ_SECTION = re.compile(r"^\s*(?:>\s*)?\**Anchors read (\d{4}-\d{2}-\d{2})\b",
                          re.MULTILINE)
HEADING = re.compile(r"^#{1,6} ", re.MULTILINE)
# Any backticked path, with or without a `:line` after it -- what a markdown
# table's file column looks like.
PATH_CELL = re.compile(r"`(" + PATH + r")(?:[:`])")
# `REF`, but with the path OPTIONAL.  The main gate's `REF` requires one, so a
# bare `` `:160` `` -- the point form of the inheritance these records use
# everywhere -- matches no rule at all and is checked by nothing: 152 of them
# in 16 records (and 4 in the three strictly gated documents, which is a hole
# this item records but does not close, because widening `REF` would move the
# strict gate's own counts).
RECORD_REF = re.compile(r"`(" + PATH + r")?:(\d+)`(?:(?= `([^`\n]+)`))?")


def relax_tokens(s: str) -> set[str]:
    """The tokens of a string, lower-cased, `RELAX_TOKEN` characters or more."""
    out: set[str] = set()
    for m in RELAX_TOK.finditer(s):
        t = m.group(0).lower()
        if len(t) >= RELAX_TOKEN and t not in RELAX_STOP:
            out.add(t)
        if "." in t:
            out.update(p for p in t.split(".")
                       if len(p) >= RELAX_TOKEN and p not in RELAX_STOP)
    return out


def column_path(text: str, pos: int):
    """A bare `` `:first-last` `` inside a markdown TABLE inherits the file
    named in its row: the path in the row's file cell, else the nearest earlier
    row of the same table that names one (a continuation row leaves the file
    cell empty).  Returns None outside a table, and for a row that names no
    file at all.

    ONE path.  A row that names SEVERAL -- `validation/vmc_reconcile.py` and
    `validation/vmc_tag_fractions.py` and a third in the same cell -- has no
    file column, and guessing the last one is how this rule first resolved
    `phase_CW_numbers.md`'s `:133-136` onto `vmc_tag_fractions.py` when the
    sentence meant the GENERATED `docs/open_items/vmc_reconciliation.md`, and
    reported a sound citation broken.  A gate that invents a failure is worse
    than one that says it cannot tell, so an ambiguous row says it cannot
    tell."""
    ls = text.rfind("\n", 0, pos) + 1
    le = text.find("\n", pos)
    le = len(text) if le < 0 else le
    if not text[ls:le].lstrip().startswith("|"):
        return None
    while True:
        seen = {m.group(1) for m in PATH_CELL.finditer(text, ls, pos)}
        if len(seen) == 1:
            return seen.pop()
        if seen or ls == 0:
            return None                       # ambiguous row, or nothing above
        end = ls - 1
        ls = text.rfind("\n", 0, end) + 1
        if not text[ls:end].lstrip().startswith("|"):
            return None
        pos = end


# A `` `name:123` `` of any shape, path or not, for the record-local
# inheritance rule below.
ANY_AT = re.compile(r"`([^`\n]*?):\d")


def record_inherited_path(text: str, pos: int):
    """`inherited_path`, but it will not READ PAST a citation it cannot
    resolve.  These records write `` (`tagged.hpp:255-257` -> `:364-375`) ``
    with the file spelled WITHOUT its directory, which `PATH` does not match;
    the main gate's rule then walks further left and inherits whatever fully
    spelled path came before it on the line -- `src/core/xsec.cpp`, four
    citations earlier -- and reports a sound citation broken (two of them, in
    `run_2026-09-03/STATUS.md`, on the first run of R6).  So: take the LAST
    `x:123` on the line before the citation whatever it looks like, and use it
    only if it is a path this gate knows.  Returns `(path, saw_one)`, so the
    caller can tell "no citation on this line" from "one this gate cannot
    resolve" and not fall through to the table-column rule for the second."""
    ls = text.rfind("\n", 0, pos) + 1
    last = None
    for m in ANY_AT.finditer(text, ls, pos):
        if m.group(1):          # a BARE `:11-22` is another citation doing the
            last = m            # same inheriting; look past it, not at it
    if last is None:
        return None, False
    ok = re.fullmatch(PATH, last.group(1)) is not None
    return (last.group(1) if ok else None), True


def relax_sentence(text: str, start: int, end: int) -> str:
    """`RELAX_WIN` characters either side of the citation, clipped to its
    paragraph (a markdown table row being its own paragraph, as R5 has it),
    with every citation in the window BLANKED OUT -- `src/core/xsec.cpp:295`
    must not be allowed to evidence `src/core/xsec.cpp`.  Only the
    `` `path:line` `` span is blanked; an adjacent `` `name` `` is the
    sentence talking, and stays."""
    lo, hi = paragraph(text, start, end)
    a0, b0 = max(lo, start - RELAX_WIN), min(hi, end + RELAX_WIN)
    out = list(text[a0:b0])
    for rx in (REF, RANGE, EXT):
        for m in rx.finditer(text, lo, hi):
            close = text.find("`", m.start() + 1) + 1
            for i in range(max(m.start() - a0, 0),
                           min(max(close - a0, 0), len(out))):
                out[i] = " "
    return "".join(out)


def marked_sections(text: str, rx: re.Pattern) -> list[tuple[int, int, str]]:
    """`(start, end, note)` for every section marker `rx` finds in a record.
    A marker reaches from its own line to the next markdown heading."""
    out = []
    for m in rx.finditer(text):
        nxt = HEADING.search(text, m.end())
        out.append((m.start(), nxt.start() if nxt else len(text),
                    m.group(1) if m.groups() else "keys"))
    return out


def asof_sections(text: str) -> list[tuple[int, int, str]]:
    return marked_sections(text, ASOF_SECTION)


def record_files() -> list[Path]:
    """Every dated run record and standing note, in a stable order."""
    base = ROOT / RECORDS_REL
    return sorted(base.rglob("*.md")) if base.is_dir() else []


def check_records(*, cache: dict) -> tuple[list[str], dict]:
    """R6 over every record.  Returns `(bad, tally)` and prints the report."""
    tally = Counter()
    bad: list[str] = []
    skipped_path: list[str] = []
    skipped_file: Counter = Counter()
    # An exemption you cannot see is a hole -- the same rule ALLOW and
    # RANGE_ALLOW follow.  Every section marker is named on every run, with
    # the number of citations it covers.
    markers: list[list] = []

    def source_lines(rel: str):
        if rel not in cache:
            f = ROOT / rel
            cache[rel] = (f.read_text(errors="replace").splitlines()
                          if f.is_file() else None)
        return cache[rel]

    for doc in record_files():
        rel = str(doc.relative_to(ROOT))
        text = doc.read_text()
        marks = marked_sections(text, ASOF_SECTION)
        readm = marked_sections(text, READ_SECTION)
        here = {}
        for lo, hi, note in marks + readm:
            row = [rel, text.count("\n", 0, lo) + 1,
                   "as of " + note if (lo, hi, note) in marks
                   else "read " + note, 0]
            here[(lo, hi)] = row
            markers.append(row)
        cites = [(m, int(m.group(2)), int(m.group(2)))
                 for m in RECORD_REF.finditer(text)]
        cites += [(m, int(m.group(2)), int(m.group(3))) for m in RANGE.finditer(text)]
        cites.sort(key=lambda c: c[0].start())
        for m, a, b in cites:
            tally["cited"] += 1
            line = text.count("\n", 0, m.start()) + 1
            path = m.group(1)
            if path is None:
                # the same line first, and if it names a citation this gate
                # cannot resolve, STOP -- a table column two cells to the left
                # is not what `(`tagged.hpp:255-257` -> `:364-375`)` means
                path, online = record_inherited_path(text, m.start())
                if path is None and not online:
                    path = column_path(text, m.start())
            if path is None:
                tally["no path"] += 1
                skipped_path.append(f"{rel}:{line}  {m.group(0)}  (the path is "
                                    f"carried by the prose, not by the citation)")
                continue
            lines = source_lines(path)
            if lines is None:
                tally["no file"] += 1
                skipped_file[path] += 1
                continue
            where = f"{rel}:{line}  -> {path}:{a}" + ("" if a == b else f"-{b}")
            if not 1 <= a <= b <= len(lines):
                tally["broken"] += 1
                bad.append(f"{where}  (out of bounds: {path} has "
                           f"{len(lines)} lines)")
                continue
            covered = [k for k in here if k[0] <= m.start() < k[1]]
            for k in covered:
                here[k][3] += 1
            # The HISTORICAL annotation is read BEFORE (i)'s non-blank test, and
            # deliberately: `(as of <commit>)` says the anchor is about that
            # commit and not about the present tree, and four of these land on
            # a line that has since gone blank -- which is not a defect in the
            # annotation, it is the thing the annotation is FOR.  In bounds is
            # still required: that catches a typed line number.
            if ASOF_INLINE.match(text, m.end()) or any(
                    (lo, hi, n) in marks for lo, hi in covered
                    for n in [here[(lo, hi)][2][6:]]):
                tally["historical"] += 1
                continue
            if not any(x.strip() for x in lines[a - 1:b]):
                tally["broken"] += 1
                bad.append(f"{where}  (the target is blank)")
                continue
            if path.startswith(RECORDS_REL + "/"):
                tally["record target"] += 1      # (iii): existence and no more
                continue
            # `(anchor read <date>)` is the opposite claim -- a human FOLLOWED
            # this anchor and it lands where the sentence says -- so it is read
            # after the non-blank test, which it has to pass like any other.
            ri = READ_INLINE.match(text, m.end())
            if ri or covered:
                if a == b:
                    # a POINT annotation carries the text the reader saw
                    pin = (ri.group(2) or ri.group(3)) if ri else None
                    if pin is None:
                        tally["broken"] += 1
                        bad.append(
                            f"{where}  (an `(anchor read)` on a POINT citation "
                            f"must carry the text the reader saw on that line: "
                            f"`(anchor read <date> \"…\")`; a section marker "
                            f"covers ranges only)")
                        continue
                    if pin not in lines[a - 1]:
                        tally["broken"] += 1
                        bad.append(
                            f"{where}  (the anchor-read pin {clip(pin, 60)!r} "
                            f"is no longer on that line, which now reads "
                            f"{clip(lines[a - 1].strip(), 60)!r}; read the "
                            f"anchor again and re-date it)")
                        continue
                tally["read"] += 1
                continue
            sentence = relax_sentence(text, m.start(), m.end())
            shared = relax_tokens(sentence) & relax_tokens(" ".join(lines[a - 1:b]))
            if a != b:
                if shared:
                    tally["ok"] += 1
                    continue
                tally["broken"] += 1
                bad.append(f"{where}  (no token of {RELAX_TOKEN}+ characters is "
                           f"in both the citing sentence and the target; "
                           f"re-point it, or say when it was true with "
                           f"`(as of <commit>)`)")
                continue
            # (ii) for a POINT citation: the shared token must be long enough
            # to mean something AND absent from the lines either side, so that
            # a one-line drift changes the answer.
            near = " ".join(lines[max(0, a - 2):a - 1]
                            + lines[b:b + 1]).lower()
            strong = {t for t in shared
                      if t not in near
                      and (len(t) >= RELAX_POINT_TOKEN
                           or (len(t) >= RELAX_POINT_NUMBER
                               and any(c.isdigit() for c in t)))}
            if strong:
                tally["ok"] += 1
                continue
            tally["broken"] += 1
            bad.append(f"{where}  (no token of {RELAX_POINT_TOKEN}+ characters "
                       f"(or {RELAX_POINT_NUMBER}+ with a digit in it) is in "
                       f"both the citing sentence and the cited LINE and "
                       f"absent from the lines either side of it"
                       + (f"; the shared tokens are "
                          + ", ".join(sorted(shared)[:5]) if shared else "")
                       + f"; re-point it, record the reading with `(anchor "
                       f"read <date> \"…\")`, or say when it was true with "
                       f"`(as of <commit>)`)")

    print(f"{tally['cited']} citations in {len(record_files())} dated run "
          f"records (relaxed, R6): {tally['ok']} checked against a live file, "
          f"{tally['record target']} into another record, "
          f"{tally['historical']} annotated historical, "
          f"{tally['read']} read by hand, "
          f"{tally['no path'] + tally['no file']} skipped, "
          f"{tally['broken']} broken  [docs/open_items]")
    for s in skipped_path:
        print("  skipped", s)
    for p in sorted(skipped_file):
        print(f"  skipped {p}  ({skipped_file[p]} citation(s)): no such file "
              f"under ROOT -- an upstream or predecessor-project path this "
              f"repository does not carry")
    for rel, line, note, n in markers:
        print(f"  marker  {rel}:{line}  anchors {note}, covering "
              f"{n} citation(s) to the end of the section")
    for x in bad:
        print("  ", x)
    return bad, tally


# Every flag `main` understands.  The parser below tests membership rather
# than parsing, so an argv element that is in neither this set nor a known
# alias would be SILENTLY IGNORED and the plain strict check would run and
# exit 0 -- a mistyped `--record-range` or `--record` would look like a clean
# gate, and `--help` would print a full passing report.  Reject leftovers
# instead: the R5/R6 modes must not be skippable by a typo.
KNOWN_FLAGS = frozenset({
    "--fix", "--record-ranges", "--audit-ranges", "--no-range-allow",
    "--loose", "--records", "--records-only",
})
USAGE = (
    "usage: check_physics_channels_links.py [--fix] [--loose]\n"
    "         [--record-ranges | --audit-ranges] [--no-range-allow]\n"
    "         [--records | --records-only]\n"
    "\n"
    "  (no flags)        the strict documentation gate\n"
    "  --fix             relocate drifted point citations in place\n"
    "  --loose           skip R4, the S3 pins and the D fingerprints\n"
    "  --record-ranges   R5: re-record the range sidecar\n"
    "  --audit-ranges    R5 without writing the sidecar\n"
    "  --no-range-allow  ignore RANGE_ALLOW exemptions during R5\n"
    "  --records         also run R6 over the dated run records\n"
    "  --records-only    run R6 alone\n"
)


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    unknown = [a for a in argv if a not in KNOWN_FLAGS]
    if unknown:
        print(USAGE, end="")
        print(f"\nunknown argument(s): {' '.join(unknown)}")
        return 2
    fix = "--fix" in argv
    write = "--record-ranges" in argv
    audit = "--audit-ranges" in argv
    # R5 runs for both: `--audit-ranges` is `--record-ranges` with the writing
    # taken out, so the rule can be measured without touching the sidecar (the
    # suite runs it over the real documents, and it is the recipe behind every
    # count published about R5).
    record = write or audit
    use_allow = "--no-range-allow" not in argv
    strict = "--loose" not in argv
    # R6, the dated run records.  Opt-in, and OFF in a recording pass: every
    # count D1 published is measured with `--audit-ranges ... | tail -1`, and a
    # second report printed after it would silently invalidate those recipes.
    records_only = "--records-only" in argv
    do_records = (records_only or "--records" in argv) and not record
    if record and not strict:
        # `--loose` skips R4, the S3 pins and the D fingerprints, so there is
        # nothing for a recording pass to collect: it would write an EMPTY
        # sidecar over the real one.  Refuse instead of destroying it.
        print("--record-ranges/--audit-ranges cannot be combined with --loose: "
              "loose mode collects no fingerprints and would write an empty "
              f"{RANGES_REL}")
        return 1
    if not DOC.exists():
        print(f"missing {DOC}")
        return 1

    cache: dict[str, Source | None] = {}
    recorded = load_ranges()
    fresh: dict[str, dict] = {}
    seen: set[str] = set()
    remap: dict[str, str] = {}
    evidence: dict[tuple, tuple] = {}
    unnamed_seen: set = set()
    ext_dir = pythia_source_dir()

    overall_bad = False
    if records_only:
        bad, _ = check_records(cache={})
        return 1 if bad else 0
    for doc, allow in [(DOC, ALLOW)] + [(ROOT / p, a) for p, a in EXTRA_DOCS]:
        if not doc.exists():
            # Only the primary document is required; an extra document that
            # does not exist under this ROOT (a test fixture tree, most of
            # the time) is simply not checked -- see the note above EXTRA_DOCS.
            continue
        bad, new_text = check_document(doc, allow, fix=fix, record=record,
                                       strict=strict, cache=cache,
                                       recorded=recorded, fresh=fresh,
                                       seen=seen, remap=remap, ext_dir=ext_dir,
                                       evidence=evidence, use_allow=use_allow,
                                       unnamed_seen=unnamed_seen)
        if new_text is not None:
            doc.write_text(new_text)
        if bad:
            overall_bad = True

    if fix and remap:
        # the block is byte-identical, only its line numbers moved, so the
        # fingerprint is carried over rather than re-taken
        moved = {remap[k]: v for k, v in recorded.items() if k in remap}
        kept = {k: v for k, v in recorded.items() if k not in remap}
        write_ranges({**kept, **moved})
        for cite in sorted(c for c in RANGE_ALLOW if c[2] in remap):
            # the exemption is keyed on the citation, and the block moved
            print(f"  note      RANGE_ALLOW[{cite!r}] should become "
                  f"{(cite[0], cite[1], remap[cite[2]])!r}: the block moved "
                  f"and the exemption is keyed on (document, citing line, "
                  f"range)")
    if ext_dir is None:
        # nothing was read, so nothing can be said: carry the external
        # fingerprints across untouched rather than dropping them.
        for k, v in recorded.items():
            if k.startswith(PYTHIA_KEY):
                fresh[k] = v
                seen.add(k)
    if strict and use_allow:
        # A dead B2 exemption hides the next drift exactly as a dead
        # fingerprint does, and it is the shape a document edit leaves behind:
        # the key carries the citing line, so a sentence that moved has an
        # exemption nothing reaches.  Only a tree that HAS the file can say so
        # (a fixture ROOT carries neither the document nor the target).
        dead_unnamed = sorted(c for c in set(UNNAMED_ALLOW) - unnamed_seen
                              if (ROOT / c[2].rsplit(":", 1)[0]).exists())
        if dead_unnamed:
            overall_bad = True
            print(f"{len(dead_unnamed)} unnamed exemption(s) no longer reached "
                  f"by any citation:")
            for c in dead_unnamed:
                print(f"   {c[0]}:{c[1]}  {c[2]}  exempted in UNNAMED_ALLOW "
                      f"but no unnamed citation reaches it; the sentence "
                      f"moved or the citation went away -- re-read the row "
                      f"and re-key the exemption, or drop it")
    if strict and not record:
        dead = sorted(set(recorded) - seen - set(remap))
        if dead:
            overall_bad = True
            print(f"{len(dead)} fingerprint(s) no longer cited by any "
                  f"covered document:")
            for key in dead:
                print(f"   {key}  fingerprinted in {RANGES_REL} but no "
                      f"longer cited; drop it with --record-ranges")
    if record:
        # R5 first: a refused range is NOT written, so the next ordinary run
        # reports it as unfingerprinted and the gate stays red until a human
        # re-points it or exempts it in RANGE_ALLOW.  One refused CITATION is
        # enough to refuse the block: the fingerprint is shared, and a block
        # blessed on one sentence while another sentence points at it wrongly
        # is exactly the hole this closed.
        refused = sorted(k for k, v in evidence.items() if v[0] is None)
        for cite in refused:
            print(f"  REFUSED   {cite[0]}:{cite[1]}  {cite[2]} -- "
                  f"{evidence[cite][1]}")
        refused_keys = {evidence[c][3] for c in refused}
        fresh = {k: v for k, v in fresh.items() if k not in refused_keys}
        by = {}
        for k, v in evidence.items():
            if v[0] is not None:
                by[v[0]] = by.get(v[0], 0) + 1
        # A dead exemption hides the next drift, exactly as a dead fingerprint
        # does -- but only a tree that HAS the file can say the exemption is
        # dead.  Under a fixture ROOT (or any tree that is not this repository)
        # the file is simply absent, and nothing is asserted, which is the same
        # rule EXTRA_DOCS and the external fingerprints already follow.
        dead_allow = sorted(k for k in set(RANGE_ALLOW) - set(evidence)
                            if (ROOT / k[2].rsplit(":", 1)[0]).exists()
                            ) if use_allow else []
        print(f"{len(evidence)} range citation(s) checked against their own "
              f"citing sentences (R5), on {len({v[3] for v in evidence.values()})} "
              f"blocks: "
              + ", ".join(f"{by.get(h, 0)} by {w}" for h, w in
                          (("name", "adjacent name"), ("text", "quoted text"),
                           ("symbol", "quoted symbol"),
                           ("allow", "RANGE_ALLOW exemption")))
              + f", {len(refused)} REFUSED")
        for cite in dead_allow:
            print(f"   {cite[0]}:{cite[1]}  {cite[2]}  exempted in "
                  f"RANGE_ALLOW but no range citation reaches it; the "
                  f"document moved or the citation went away -- re-read the "
                  f"row and re-key the exemption, or drop it")
        if refused or dead_allow:
            overall_bad = True
        if not write:
            return 1 if overall_bad else 0
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
            print("  dropped   ", k,
                  " (REFUSED by R5, so it is no longer fingerprinted)"
                  if k in set(refused) else "")
        print(f"{len(fresh)} fingerprints written to {RANGES_REL} "
              f"({len(added)} new, {len(changed)} changed, {len(dropped)} "
              f"dropped): cited BLOCKS, pinned USE SITES and EXTERNAL lines")
    if do_records:
        # printed LAST, after every strict document, so that nothing above it
        # moves and `| tail -1` still names what it always named
        bad, _ = check_records(cache={})
        if bad:
            overall_bad = True
    return 1 if overall_bad else 0


if __name__ == "__main__":
    sys.exit(main())
