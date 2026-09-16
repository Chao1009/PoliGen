# Phase D numbers — run 2026-09-06

The gate residue class: what `validation/check_physics_channels_links.py` still
could not see. Performed 2026-09-15 against the tree as committed at **a3c9ecb**
(clean at start). No file under `src/` or `include/` is touched by this phase —
only the gate script, its sidecar, its pytests and the two covered documents.

---

## D1. `--record-ranges` now refuses a block the citing sentence does not evidence

`PLAN.md` D1: *"`--record-ranges` blesses content at recording time; a range
recorded onto the wrong block stays "0 broken". Require, at recording, that the
citing sentence's quoted text or named symbol occurs inside the block, and
refuse to record otherwise."*

### D1.1 The one sentence

> **R5 is in the gate: `--record-ranges` refuses to fingerprint a cited block
> unless the citing sentence evidences it — an adjacent name that is declared
> or used inside the block, a quoted phrase (backticks or quotation marks,
> within `RECORD_WIN` = 200 characters, clipped to the citing paragraph, where a
> markdown table row is its own paragraph) that occurs inside it, or an explicit
> per-entry `RANGE_ALLOW` exemption carrying a PIN and a REASON — and a refused
> entry is simply not written, so the next ordinary run reports it as
> unfingerprinted and the gate stays red until a human acts. Applied from
> scratch to every range in the sidecar, the rule **REFUSED 89 of the 180 cited
> blocks**; reading all 89, **20 were pointed at the wrong block** (plus 2 more
> found by eye in the same sentences, which the rule accepted), so **22 range
> keys / 23 citations were re-pointed by content**, and the remaining **69 are
> exempted by name in `RANGE_ALLOW`** with a pin and a reason each (74
> exemptions in all, the extra 5 being re-pointed blocks whose sentence still
> evidences nothing). End state: 180 blocks, **0 refused**, gate 0 broken.**

### D1.2 What the rule is

Three branches, in order, evaluated per citation and pooled per block (a block
cited twice is evidenced if *either* sentence evidences it):

| branch | what it takes | what it requires |
|---|---|---|
| (a) adjacent name | `` `path:a-b` `name` `` — one space, the same adjacency the point citations use | the name is **declared or used** anywhere in the block (R3, stricter, already ran and demands *declared*) |
| (b) quoted phrase | every backtick / `"…"` / `“…”` span within `RECORD_WIN` characters, minus the citations themselves | at least one occurs in the block as **text** (whitespace-normalised, ≥ 3 chars; a bare identifier must match on word boundaries) or as a **symbol** (identifier, qualified name, call or flag ≥ 2 chars, declared or used in the block) |
| (c) neither | — | an explicit `RANGE_ALLOW[key] = (pin, reason)`; the **pin must occur inside the block and nowhere else in that file**, so the exemption names the block by content, not by two line numbers that drift |

A refusal prints what it would have needed, phrase by phrase. Two flags were
added: `--audit-ranges` (run R5, write nothing) and `--no-range-allow`
(suppress the exemption table) — the second is what every count below is
measured with.

Details that matter, each of them a bug found while measuring:

* **quoted phrases are read off the whole paragraph, then filtered by
  distance.** A window cut through the middle of a code span pairs the
  surviving backtick with the next one and reads a run of prose as a quotation
  (it invented phrases like `` `), fed through the dedicated J = 3/2 sl…` ``).
* **a quotation may wrap.** These documents are hard-wrapped at ~80 columns, so
  `docs/CONVENTIONS.md:72-73`'s own quotation — *"b₁ > 0 means the m = 0 state
  has the LARGER cross section"* — runs across a line break. A blank line still
  ends a phrase; 240 characters caps one.
* **a slash-run of names is every name in it** (`HelicityFlipOptions::theta_s/phi_s`,
  `pe/pz/pzz`) but only when every part is a name: `data/vmc/README.md` is a
  path, not three symbols. Worth exactly 2 blocks.
* **a bare identifier matches on word boundaries** in the text branch: `rates`
  inside "generates" is not the block naming it. Worth 4 refusals.
* **case is significant.** `include/lipolgen/xsec.hpp:99-103` (anchor read 2026-09-15) is cited for the
  collapse *"for any b2"* while the header writes *"for ANY b2"*; making the
  text branch case-insensitive recovers that one block and costs +0.5 points of
  decoy acceptance (§D1.6), so it is **not** done and the discrepancy is
  recorded in that block's exemption reason instead.

### D1.3 The measurement, and how to reproduce it

**Against the committed pre-fix tree** (this is the "how many would have been
refused" the phase was asked for):

```
git archive a3c9ecb | tar -x -C /tmp/d1head
cp validation/check_physics_channels_links.py /tmp/d1head/validation/
cd /tmp/d1head && python3 validation/check_physics_channels_links.py \
    --audit-ranges --no-range-allow | tail -1
```

> `180 range block(s) checked against their citing sentences (R5): 0 by adjacent
> name, 71 by quoted text, 20 by quoted symbol, 0 by RANGE_ALLOW exemption,
> **89 REFUSED**`

**That is what the rule D1 landed reports. The same recipe with the gate this
tree now carries reports**

> `210 range citation(s) checked against their own citing sentences (R5), on 180
> blocks: 0 by adjacent name, 59 by quoted text, 31 by quoted symbol, 0 by
> RANGE_ALLOW exemption, **120 REFUSED**`

— 210 CITATIONS rather than 180 blocks, because the verdict is per citation
(§D3.3), and more of them refused, because a common word and a very long block
stopped counting as evidence (§D3.5). Both numbers are this recipe; they differ
in the rule, and the rule is named in the line. Everything below this point in
§§D1.3–D1.6 is the D1 rule's measurement, kept as the record of what was
decided on it.

**Against the tree this phase leaves** (`--no-range-allow`, then without it):

```
python3 validation/check_physics_channels_links.py --audit-ranges --no-range-allow | tail -1
python3 validation/check_physics_channels_links.py --audit-ranges | tail -1
```

> `… 1 by adjacent name, 83 by quoted text, 22 by quoted symbol, 0 by
> RANGE_ALLOW exemption, **74 REFUSED**`
> `… 1 by adjacent name, 83 by quoted text, 22 by quoted symbol, 74 by
> RANGE_ALLOW exemption, **0 REFUSED**`

(again, the D1 rule; under the rule this tree leaves the same two runs report
**103 REFUSED** and **0 REFUSED** over 209 citations on 179 blocks — §D3.8)

| | pre-fix (a3c9ecb) | after D1 |
|---|---|---|
| cited blocks | 180 | 180 |
| evidenced by an adjacent name | 0 | 1 |
| evidenced by quoted text | 71 | 83 |
| evidenced by a quoted symbol | 20 | 22 |
| refused on the sentence alone | **89** | 74 |
| exempted in `RANGE_ALLOW` | — | 74 |
| **refused after exemptions** | — | **0** |

The 89 split, one by one (every one of the 89 was read against its block):

| outcome | count |
|---|---|
| the citation was pointing at the **wrong block** → re-pointed by content | **20** |
| the citation is right, the sentence evidences nothing a machine can read → exempted | 69 |
| *(found while reading a refusal's neighbour, accepted by R5, also wrong)* | *(+2)* |

So **a refusal is not a synonym for a wrong citation** — it is *"this gate
cannot tell"* — and here it was a wrong citation 22 % of the time (20/89).
The 69 are the shape of these documents: a range is usually pointed at by prose
that names its subject in unquoted words ("provenance at", "the inventory table
at", "the joint-sampling order"), which no rule of this kind can read.

### D1.4 The 22 re-pointed range citations

All by CONTENT: each new block was read and is the thing the sentence says it
is. Old ranges are what the sidecar had blessed at a3c9ecb.

**Anchors read 2026-09-15** in this section: every range below is the sidecar
key itself, named so that it can be looked up, not a pointer this document
asks you to follow. The `was` column is deliberately a block that is no longer
what the citation means; the `is now` column is gated by R4 and R5 above, in
the documents that actually cite it.

| document | line | was | is now | the old block actually was |
|---|---|---|---|---|
| `PHYSICS_CHANNELS.md` | 107 | `src/core/xsec.cpp:266-270` | `:295-298` | the Δ/RC comment inside `tensor_amplitudes` |
| `PHYSICS_CHANNELS.md` | 274 | `include/lipolgen/tagged.hpp:480-484` | `:496-499` | the `TaggedEvent` struct, not the joint-sampling order |
| `PHYSICS_CHANNELS.md` | 381 | `include/lipolgen/rc.hpp:195-202` | `:305-325` | the HERMES τ formula, not the E12-13-011 anchor |
| `PHYSICS_CHANNELS.md` | 384 | `include/lipolgen/rc.hpp:251-258` | `:414-421` `RC_QE_KF_GEV` | the tagged-τ clamp, not Moniz's k_F(⁶Li) = 0.169 |
| `SPIN32_FINITE_GAMMA.md` | 265 | `python/bindings.cpp:1852-1853` | `:2334-2335` | `LightConeDensities` bindings, not `use_explicit_pzz`/`pzz` |
| `SPIN32_FINITE_GAMMA.md` | 400 | `src/core/xsec.cpp:267-269` | `:296-298` | the Δ/RC comment |
| `SPIN32_FINITE_GAMMA.md` | 415 | `src/core/xsec.cpp:257-260` | `:261-265` | the finite-γ harmonics accumulation, not the `TENSOR_LL_SIGN` term |
| `SPIN32_FINITE_GAMMA.md` | 462 | `include/lipolgen/xsec.hpp:324-326` | `:350-352` | a signature tail, not `positivity_margin` |
| `SPIN32_FINITE_GAMMA.md` | 741 | `include/lipolgen/xsec.hpp:266-267` | `:274-275` | prose above Eq. (35), not Eq. (35) |
| `SPIN32_FINITE_GAMMA.md` | 905 | `src/core/xsec.cpp:266-270` | `:295-298` | *(same key as PHYSICS_CHANNELS.md:107)* |
| `SPIN32_FINITE_GAMMA.md` | 916 | `src/core/xsec.cpp:283-285` | `:320-328` | `amplitudes()`'s signature, not `dsigma` = Eq. (46) |
| `SPIN32_FINITE_GAMMA.md` | 978 | `src/core/sampler.cpp:234-240` | `:276-279` | an A₁-refusal comment, not the `a1n`/`a2n`/`bound`/`margin` fills |
| `SPIN32_FINITE_GAMMA.md` | 979 | `src/core/sampler.cpp:386-396` | `:425-435` | the m draw, not the φ accept–reject |
| `SPIN32_FINITE_GAMMA.md` | 980 | `src/core/sampler.cpp:305-322` | `:335-362` | the `state_tables` cache, not `effective_modulation` |
| `SPIN32_FINITE_GAMMA.md` | 981 | `src/core/sampler.cpp:480-482` | `:519-521` | a thread join, not the per-category φ density |
| `SPIN32_FINITE_GAMMA.md` | 981 | `src/core/sampler.cpp:526-531` | `:595-604` | `tensor_weights_for`, not `phi_histogram_pseudo` |
| `SPIN32_FINITE_GAMMA.md` | 985 | `python/bindings.cpp:1719-1721` | `:2220-2222` | a CD-Bonn sum docstring, not the `_32` `Options` bindings |
| `SPIN32_FINITE_GAMMA.md` | 986 | `python/bindings.cpp:1679-1688` | `:2180-2192` | `cdks_b1_raw_per_nucleon`, not the `Amplitudes` binding |
| `SPIN32_FINITE_GAMMA.md` | 987 | `python/bindings.cpp:2005-2006` | `:2481-2494` | `DeuteronConvolutionB1::densities`, not the `state_tables` dict |
| `SPIN32_FINITE_GAMMA.md` | 995 | `src/core/xsec.cpp:214-276` | `:282-304` | `tensor_amplitudes`, not `amplitudes` — **R5 accepted this one** |
| `SPIN32_FINITE_GAMMA.md` | 996 | `src/core/xsec.cpp:248-264` | `:253-277` | five lines short of the rank-2 branch — **R5 accepted this one** |
| `SPIN32_FINITE_GAMMA.md` | 997 | `src/core/xsec.cpp:266-274` | `:295-302` | the Δ/RC comment, not the λ_e vector block |
| `SPIN32_FINITE_GAMMA.md` | 1031 | `src/core/xsec.cpp:272-274` | `:300-302` | the `with_delta` guard, not the `with_perp` guard |

Two families account for nearly all of it and both are **stale refactors, not
typos**: the P8 split of `InclusiveKernel::amplitudes` into
`tensor_amplitudes` + `amplitudes` moved nine `src/core/xsec.cpp` blocks out
from under `SPIN32_FINITE_GAMMA.md` §§5–6, and `src/core/sampler.cpp` and
`python/bindings.cpp` grew by ~40 and ~500 lines under that note's §6.1 change
list. R4's fingerprints could not see any of it: the blocks the citations
landed on were fingerprinted **as cited**, so the gate agreed with itself.

One citation gained a name in the process: `PHYSICS_CHANNELS.md:384` now reads
`` `include/lipolgen/rc.hpp:414-421` `RC_QE_KF_GEV` ``, which is rule (a)
evidence and R3-checked as well — the only adjacent name on any range in the
three covered documents.

### D1.5 The 74 exemptions

**Anchors read 2026-09-15** in this section, for the same reason as §D1.4: an
exemption is keyed on the range, so the range is spelled out as a key.
**Since §D3.3 the key is (document, citing line, range) and there are 103 of
them**; the 74 below are the blocks D1 exempted, and 64 of their reasons were
carried over unchanged when the table was re-keyed.

`RANGE_ALLOW` (in the gate, beside `ALLOW`, not in the generated sidecar), keyed
on the exact range and carrying `(pin, reason)`. 69 are pre-existing blocks whose
sentence evidences nothing; 5 are blocks created by a re-point in §D1.4 whose
new sentence still evidences nothing (`include/lipolgen/xsec.hpp:274-275`,
`src/core/sampler.cpp:519-521`, `src/core/xsec.cpp:295-298`, `:296-298`,
`:320-328`). Each was read on 2026-09-15 before it was written.

The pin is what makes an exemption drift-proof, and the gate enforces it at
every recording: **the pin must be inside the block and must occur exactly once
in the file.** A pin that slides off its block, or that matches twice, is a
refusal with its own message — tested both ways. `--fix` prints the key that
has to move when it re-points an exempted block.

### D1.6 Why `RECORD_WIN` = 200, measured

Widening the window buys evidence and sells discrimination. Two sweeps over the
three covered documents, both on the tree this phase leaves. **The numbers in
this section are the D1 rule's**; re-taken with the gate this tree now carries
(§§D3.3, D3.5), 120 / 200 / 400 characters evidence **87 / 106 / 119** of the
210 range citations while accepting a decoy **8.7 % / 12.7 % / 18.3 %** of the
time — the same shape, the same conclusion, the same choice of 200, and the
true-to-decoy ratio is uniformly better (4.78 / 3.98 / 3.10 against 4.5 / 4.1 /
3.3). The recipe for the re-take is `validation/record_rule_sweep.py r5` with
`RECORD_WIN` set as the loop below sets it.

```
for W in 120 200 300 400 600; do
  sed -i "s/^RECORD_WIN = .*/RECORD_WIN = $W/" validation/check_physics_channels_links.py
  python3 validation/check_physics_channels_links.py --audit-ranges --no-range-allow | tail -1
done
sed -i "s/^RECORD_WIN = .*/RECORD_WIN = 200/" validation/check_physics_channels_links.py
```

and, for the discriminating power, the decoy test — the same citation pointed at
a block of the same length in the same file, shifted clear of the true one:

```python
# scratch script; paste and run from the repository root
import importlib.util, sys
from pathlib import Path
R = Path(".").resolve()
sp = importlib.util.spec_from_file_location("g", R/"validation/check_physics_channels_links.py")
g = importlib.util.module_from_spec(sp); sp.loader.exec_module(g)
g.RECORD_WIN = int(sys.argv[1])
true, ok, n = {}, 0, 0
for rel in ["docs/PHYSICS_CHANNELS.md"] + [p for p, _ in g.EXTRA_DOCS]:
    text = (R/rel).read_text()
    for m in g.RANGE.finditer(text):
        path = m.group(1) or g.inherited_path(text, m.start())
        a, b, nm = int(m.group(2)), int(m.group(3)), m.group(4)
        if nm and g.CITATION.match(nm): nm = None
        src, key = g.Source(R/path, path), f"{path}:{a}-{b}"
        how, _ = g.range_evidence(text, m, src, a, b, nm)
        if how or key not in true: true[key] = how
        L = b - a + 1
        for sh in (-(L+40), -(L+5), L+5, L+40):
            na, nb = a+sh, b+sh
            if not (1 <= na < nb <= len(src)): continue
            if not src.lines[na-1].strip() or not src.lines[nb-1].strip(): continue
            n += 1
            if g.range_evidence(text, m, src, na, nb, nm)[0]: ok += 1
t = sum(1 for v in true.values() if v)
print(f"RECORD_WIN={g.RECORD_WIN}: true {t}/{len(true)} = {t/len(true):.1%}; "
      f"decoys {ok}/{n} = {ok/n:.1%}")
```

| `RECORD_WIN` | blocks evidenced (of 180) | refused | decoy blocks accepted (of 603) |
|---|---|---|---|
| 120 | 88 (48.9 %) | 92 | 65 (10.8 %) |
| **200** | **106 (58.9 %)** | **74** | **86 (14.3 %)** |
| 300 | 115 (63.9 %) | 65 | 104 (17.2 %) |
| 400 | 120 (66.7 %) | 60 | 123 (20.4 %) |
| 600 | 125 (69.4 %) | 55 | 135 (22.4 %) |

The true/decoy ratio falls monotonically — 4.5 / **4.1** / 3.7 / 3.3 / 3.1 — so
every widening is paid for. 200 is taken because the whole of the 120 → 200 gain
is a table row's claim cell reaching its own provenance cell, which is the same
sentence; beyond that the window reaches material about other things.

**Clipping at the table-cell `|` instead** (so only the cell holding the
citation counts) was measured too and is NOT used: it evidences **80 of the 211
citations (75 of 180 blocks) against 112 (106 blocks)**, and the 32 it loses are
rows whose claim cell quotes the very symbol the cited block declares. In these
documents the cell boundary is not a sentence boundary; the row is. Recipe: the
snippet above with `g.paragraph` replaced by one that, on a line starting `|`,
narrows `(lo, hi)` to the nearest `|` either side of the citation.

### D1.7 What R5 still cannot see

1. **It can only judge what the sentence quotes.** 74 of 180 blocks (41 %) say
   nothing checkable and live on an exemption; R5 protects the other 59 %.
   (Per citation, and with §D3.5's minimum evidence, it is 103 of 209
   citations — 49 % — on the tree §D3 leaves.)
2. **A decoy is accepted 14.3 % of the time** (§D1.6): a quoted phrase that also
   occurs in a *different* block of the same file evidences that block equally.
   R5 raises the cost of a wrong recording; it does not make one impossible.
3. **A quoted identifier that the block carries only as PROSE counts.** 8 of the
   106 evidenced blocks rest on that alone — the block's comment names the
   symbol but no line uses it as code. `src/core/xsec.cpp:214-276` was wrong
   *and* accepted for exactly this reason (the comment inside `tensor_amplitudes`
   discusses `amplitudes()`), which is how two of the 22 escaped the rule and
   were caught by eye instead. Recipe — the §D1.6 loop with the per-phrase test
   replaced by:

   ```python
   def as_code(src, a, b, q):
       return any(src.declares(i, base) or src.in_code(i, base)
                  for base in g.symbol_bases(q) for i in range(a, b + 1))
   # per phrase q: how = g.phrase_in_block(src, a, b, q)
   #   how is None                       -> no evidence
   #   how == "symbol" or " " in q.strip() or not g.symbol_bases(q)
   #     or as_code(src, a, b, q)        -> solid  (score 2)
   #   otherwise                         -> prose only (score 1)
   # a block's score is the best over all its citations' phrases
   ```

   which prints `evidenced: 106 | of those, prose-only: 8` and names them.
4. **R5 runs at RECORDING time only.** A block that was evidenced when recorded
   and is edited afterwards is caught by R4 (the fingerprint), which then forces
   a re-record — and the re-record re-runs R5. A range nobody re-records is
   never re-judged.
5. **The exemptions are as good as the reading behind them.** The pin makes the
   exemption drift-proof, not correct: it says "a human read this block", and if
   that human was wrong the gate will keep agreeing with them.
6. **(Added 2026-09-15.) The verdict was pooled per BLOCK, and this list did
   not say so.** "A block cited by two documents is evidenced if either
   sentence evidences it" is in the gate's own docstring, and it meant a
   citation re-pointed onto a blessed block was never judged at all: three of
   five wrong-block attacks passed, and 11 citations on 7 keys had ridden in
   free since the day R5 landed. That is §D3.3, and it was the largest single
   hole in this item's own subject.

### D1.8 Suites, gates, tallies

Baseline (a3c9ecb, tree clean) and end state, both measured here:

| | baseline | after D1 |
|---|---|---|
| `build/lipolgen_tests` | 408 cases / 17 241 975 assertions / 0 failed | **408 / 17 241 975 / 0 failed** |
| `python -m pytest python/tests -q` | 1023 passed / 151 skipped | **1039 passed / 151 skipped** (+16, all in `test_doc_link_gate.py`) |
| docs gate (`PHYSICS_CHANNELS.md`) | 1236 refs / 96 ranges / 7 external / 0 broken / 6 allow-listed | **1236 / 96 / 7 / 0 broken / 6 allow-listed** |
| `[SPIN32_FINITE_GAMMA.md]` | 19 refs / 115 ranges / 0 broken | **19 / 115 / 0 broken** |
| `[PYTHIA_BRIDGE.md]` | 8 external / 0 broken | **8 external / 0 broken** |
| R5 audit | *(did not exist)* | **180 blocks, 0 REFUSED** |
| `check_spdx_headers.py` | 102/102 | **102/102** |
| sidecar `validation/physics_channels_ranges.json` | 220 entries | **220 entries** (22 new, 0 changed, 22 dropped — the 22 re-pointed keys), composed of **180 blocks / 31 use-site lines / 9 upstream lines** |

The 16 new pytests are the D1 rule on fixtures — a wrong-block range refused and
not written, a right-block one recorded, a sentence that evidences nothing
refused *with what it would have needed*, a refusal dropping a fingerprint that
was already there, a wrapped quotation, rule (a)'s use-as-well-as-declaration,
an exemption accepted, an exemption whose pin left the block, an exemption whose
pin matches twice, a dead exemption, `--audit-ranges` writing nothing,
`--no-range-allow` actually suppressing the table, `--fix` naming the exemption
key that moved, and `--record-ranges --loose` refused — plus two over the real
documents (every cited block evidenced or exempted; every exemption live, pinned
and reasoned).

**One pre-existing footgun closed on the way**: `--record-ranges --loose` used to
write an EMPTY sidecar over the real one (loose mode collects no fingerprints at
all, by design), and is now refused by name. Nothing else about `--loose`
changed: `--loose` alone still reports 1236 references / 96 ranges / 7 external /
0 broken, as it did before.

**Two counts in the gate's own header were one out and are corrected here**: the
sidecar is 180 blocks + 31 use-site lines + 9 upstream lines, not 181 + 30 + 9
(recount: an entry with a `first`/`last` pair is a block, one with a `line` is a
use site, one keyed `pythia8:` is upstream).

**No source file changed.** `git diff --stat` for this item touches
`validation/check_physics_channels_links.py`,
`validation/physics_channels_ranges.json`, `python/tests/test_doc_link_gate.py`,
`docs/PHYSICS_CHANNELS.md` (4 citations) and
`docs/theory/SPIN32_FINITE_GAMMA.md` (19 citations) — and this file.

---

## D2. The dated run records are gated, under a relaxed rule of their own

`PLAN.md` D2: *"the dated run records under `docs/open_items/run_*/` are not
gated at all and their citations into the live documents drift with every line
shift; fourteen were re-pointed by hand at the last close-out and the rest may
be off."*

### D2.1 The one sentence

> **R6 is in the gate: `--records` checks every `path:line` citation in
> `docs/open_items/` — the three dated run-record directories and the older
> standing notes beside them — under one RELAXED rule, because a record is HISTORY and the
> strict rules are the wrong rules for it: the target must be IN BOUNDS and
> NON-BLANK, and, when it is a LIVE file, must share at least one token of four
> characters or more (a symbol, a number, a word) with the citing sentence;
> a citation into another record is checked for existence and non-blankness
> only; nothing is fingerprinted. Against the committed tree a3c9ecb the rule
> found **243 of the 600 citations broken** — 233 on the token test, 10 on a
> blank target — spread over **46 target documents** and **22 of the 41
> records**, and **218 of the 243 are in `run_2026-09-03/` alone**. All 243
> were read: **137 citations were re-pointed** (the thing they name is still in
> the tree and only the line moved), **95 were annotated `(as of <commit>)`**
> because they describe a state that no longer exists, **37 were annotated
> `(anchor read 2026-09-15)`** because they land exactly where the sentence
> says while the sentence names its subject in words the target never uses, and
> **26 needed both** (137 + 132 − 26 = 243). End state: **664 citations, 0
> broken**, with the strict gate's every published count unchanged.**

*(§D3.4 strengthened the token test for POINT citations and §D3.6 (iii)–(iv)
re-derived the splits; the same recipe on the same tree now reports 266 broken
rather than 243, and the end state is 714 citations and 0 broken. The sentence
above is what D2 landed, and the numbers in it are that rule's.)*

### D2.2 What R6 is, and why it is not R1–R5

A record says what was true on its date. R4's fingerprints would freeze a block
the record never claimed was frozen; R5 is a *recording* rule and there is
nothing here to record; S1/S2/S3's declaration tests would refuse the
legitimate citation of a use site, or of a line of prose in another document.
So R6 asks only what a reader following the anchor would ask:

| | test | applies to |
|---|---|---|
| (i) | in bounds, and **at least one non-blank line** in the target | every citation |
| (ii) | at least one shared **token of ≥ 4 characters** between the citing sentence and the target line or block | a target OUTSIDE `docs/open_items/` |
| (iii) | (i) and nothing else | a target INSIDE `docs/open_items/` — one record quoting another is a pointer inside the history, and the two move together or not at all |

(i) is deliberately weaker than R2, which wants a non-blank **first and last**
line: R2 protects a fingerprint and nothing here is fingerprinted. The citing
sentence is R5's — `RELAX_WIN` = 200 characters either side, clipped to the
paragraph, a markdown table row being its own paragraph — with every citation
in the window **blanked out** first, so `src/core/xsec.cpp:295` (anchor read 2026-09-15 "lam_e * state") cannot
evidence `src/core/xsec.cpp`. Tokens are compared case-insensitively; a dotted run is
itself and its parts (`0.440`, `xsec.cpp` → also `xsec`); and 115 common English
function words are not evidence (`which`, `because`, `that` …), because both
sides of this comparison are prose as often as they are code. Measured: against
a3c9ecb that list turns **16 apparent passes into failures** — every one of them
a citation sharing nothing but a function word — and drops decoy acceptance
**32.6 % → 27.7 %** (§D2.4's test; true/decoy 1.30 → 1.38) — these are the D2-rule figures, measured with the D2 gate that no commit carries; the COMMITTED `record_rule_sweep.py tokens --no-stoplist` gives 17 passes, 31.4 % → 26.4 % (1.29 → 1.37) on a3c9ecb; against the tree
this phase leaves it costs **no** passes and drops decoys **39.0 % → 34.7 %** (2.56 → 2.88; the D2 gate read 39.1 → 34.6)
(2.56 → 2.89). Words that mean something in this tree — `line`, `time`,
`state`, `open`, `make`, `show`, `part`, `left`, `still` — are deliberately
**not** on it.

**The two annotations.** A failure is not always a wrong anchor, and R6's
repair vocabulary says which:

* `` `path:line` (as of <commit>) `` — the anchor is about that commit, not
  about the present tree. Read **before** (i)'s non-blank test, deliberately:
  **two** of the 95 land on a line that has since gone blank — the thing the
  annotation is *for*. (This read "four" until 2026-09-15 §D3.6: four have a
  blank FIRST-OR-LAST line, which is R2's criterion and R2 does not run here;
  under R6's own (i) — *any* non-blank line in the block — it is two,
  `phase_A_gate_mechanics.md:751` and `:754`, both into
  `docs/OPEN_ITEMS_SOLUTIONS.md`.) In bounds is still required, so a typed line number
  is still caught. Every commit named is resolved by the suite
  (`test_every_as_of_commit_in_the_records_resolves`).
* `` `path:line` (anchor read <date>) `` — the opposite claim: a human followed
  this anchor and it lands where the sentence says, but the sentence names its
  subject in unquoted words the target never uses. Read **after** the non-blank
  test, which it has to pass like any other.

Both have a section form — `**Anchors as of `<sha>`**` and
`**Anchors read <date>**`, reaching from the marker to the next markdown
heading — and every marker is named on every run with the number of citations
it covers, the same rule `ALLOW` and `RANGE_ALLOW` follow. Two are in use, both
in §D1.4 and §D1.5 above, where the citations are sidecar **keys** quoted as
keys rather than pointers; together they cover **51 citations, 22 of which
would otherwise be broken** (measured by deleting the two marker paragraphs and
re-running).

**What R6 does not check.** 190 citations are skipped, and both classes are
named on every run rather than dropped:

* **181 whose PATH is carried by the prose.** `` The design (`:1285-1289`)
  states the window as … `` inherits its path from the word "design" three
  paragraphs up. Two inheritance rules *are* implemented — the main gate's
  same-line one, and a markdown-TABLE one (`column_path`: the path in the row's
  file cell, else the nearest earlier row of the same table) — and between them
  they resolve **113 of the 294 bare citations**. The rest no rule of this
  kind can read.
* **9 whose target is not in this repository**, and not one is a LiPolGen path:
  `src/nucleons.cpp`, `src/inputParameters.cpp`, `src/subnucleon_config.hpp`
  and `src/main.cpp` are eSTARlight's own sources (the design notes that read
  them), `src/HardDiffraction.cc` and `examples/main234.cc` are PYTHIA's, and
  `docs/note_cos2phi_coherent_6Li.md`, `docs/note_7li_theory_questions.md` and
  `docs/consistency_review_2026-09-02.md` belong to the predecessor project the
  physics notes quote. They are listed by name every run — exactly as the
  PYTHIA externals are — so deleting a live file shows up here as a new skipped
  line rather than as silence.

### D2.3 The measurement, and how to reproduce it

**Against the committed pre-fix tree:**

```
git archive a3c9ecb | tar -x -C /tmp/d2head
cp validation/check_physics_channels_links.py /tmp/d2head/validation/
cd /tmp/d2head && python3 validation/check_physics_channels_links.py --records-only | head -1
```

> `600 citations in 41 dated run records (relaxed, R6): 120 checked against a
> live file, 31 into another record, 0 annotated historical, 0 read by hand,
> 183 skipped, **266 broken**`

**266, not the 243 this section published on 2026-09-15.** Both numbers are
that recipe; they differ because the RULE differs. D2 landed one shared token
of four characters for every citation and found 243; D3 (§D3.3) holds a POINT
citation to a token of six (or five with a digit) that the lines either side do
not carry, and the same pre-fix tree then has 266 broken. The recipe names the
gate it is run with, so the number a reader gets today is 266; 243 is what the
D2 rule gave and is kept here because §D2.5's repair split is against it.

**Against the tree this phase leaves:**

```
python3 validation/check_physics_channels_links.py --records-only | head -1
```

> `714 citations in 42 dated run records (relaxed, R6): 266 checked against a
> live file, 35 into another record, 107 annotated historical, 115 read by
> hand, 191 skipped, **0 broken**`

The 114-citation difference between the two is `phase_D_numbers.md` itself —
the 55 citations D1 wrote, which a3c9ecb does not carry, the 9 D2 added and the
50 §D3 adds. A record that gates itself is the point: four of D2's nine were
reported broken on its first run, and nineteen of §D3's citations (the section holds 43; 114 − 55 − 9 = 50 overall; the first-run state is not reproducible) on its
first run; all of them are re-pointed or annotated here.

| | pre-fix (a3c9ecb) | after D2 | after D3 |
|---|---|---|---|
| citations found | 600 | 664 | 714 |
| checked against a live file | 120 | 257 | 266 |
| checked into another record | 31 | 31 | 35 |
| annotated `(as of <commit>)` | — | 96 | 107 |
| annotated `(anchor read …)` / under a marker | — | 90 (39 + 51) | 115 (64 + 51) |
| skipped (prose path / not in this tree) | 183 | 190 | 191 |
| **broken** | **266** *(243 under the D2 rule)* | **0** *(under the D2 rule)* | **0** |

**The split by failure reason** (same recipe, on the pre-fix tree, D3 rule):

```
python3 validation/check_physics_channels_links.py --records-only \
  | grep '^   docs/open_items' \
  | grep -o '(no token of 4+\|(no token of 6+\|(the target is blank\|(out of bounds' \
  | sort | uniq -c
```

| reason | count |
|---|---|
| a RANGE with no shared token of 4+ characters | 94 |
| a POINT citation with no token of 6+ (or 5+ with a digit) off its neighbours | 162 |
| the target is blank | 10 |
| out of bounds | 0 |

**THE SPLITS BELOW ARE THE D3 RULE'S**, re-taken 2026-09-15 (§D3.6 (iii),
(iv)): the tables that stood here summed to 203 against a headline of 243 —
`include/lipolgen/b1_nuclear.hpp` was written 5 where the run says 6,
`include/lipolgen/pythia_bridge.hpp`'s 6 was missing altogether, and the
remainder row was wrong in both count and cap. They are re-derived rather than
patched, and under the rule the gate now carries, because **the D2 rule exists
in no committed artifact**: the gate at a3c9ecb is the pre-D1 gate, and the
only gate a reader can run is this one. A 243-sum table would therefore be a
number nobody can reproduce, which is the defect being repaired. The recipe is
the run above, on the same pre-fix checkout, piped through `sed`.

**The split by TARGET file** — 48 distinct targets, 266 broken,
`| grep '^   docs/open_items' | sed 's/.*-> //; s/:[0-9-]*  .*//' | sort |
uniq -c | sort -rn`:

| target | broken | | target | broken |
|---|---|---|---|---|
| `src/core/pipeline.cpp` | 28 | | `python/lipolgen/__init__.py` | 8 |
| `python/bindings.cpp` | 23 | | `tests/test_rc.cpp` | 6 |
| `README.md` | 14 | | `include/lipolgen/b1_nuclear.hpp` | 6 |
| `docs/OPEN_ITEMS_SOLUTIONS.md` | 12 | | `src/core/rc.cpp` | 6 |
| `python/lipolgen/cli.py` | 12 | | `include/lipolgen/pythia_bridge.hpp` | 6 |
| `docs/PHYSICS_CHANNELS.md` | 11 | | `validation/o5_a2_reach.py` | 6 |
| `docs/DEVELOPMENT_PLAN.md` | 11 | | `include/lipolgen/rc.hpp` | 5 |
| `docs/USAGE.md` | 10 | | `include/lipolgen/xsec.hpp` | 5 |
| `tests/test_b1_nuclear.cpp` | 10 | | *(29 more, ≤ 4 each)* | 68 |
| `src/core/b1_nuclear.cpp` | 10 | | `include/lipolgen/pipeline.hpp` | 9 |

19 targets carry 5 or more (198 of the 266); the remaining 29 carry 68 between
them, none more than 4. **198 + 68 = 266**, and that is the arithmetic the
previous table failed.

**And the split by RECORD** (`sed 's#^   docs/open_items/##; s/:.*//'`):
`engineering.md` **2**, `run_2026-09-02/` **15**, `run_2026-09-03/` **236**,
`run_2026-09-06/` **13**. The concentration is not an accident: the 2026-09-03
run is the one whose own subject matter — the A = 2 gate, the b₁ backends, the
SF-injection design — was rewritten by the three runs that came after it.

| record | broken | | record | broken |
|---|---|---|---|---|
| `run_2026-09-03/phase_D_sf_injection.md` | 68 | | `run_2026-09-02/estarlight_li6.md` | 5 |
| `run_2026-09-03/phase_A_gate_mechanics.md` | 41 | | `run_2026-09-03/phase_D_numbers.md` | 5 |
| `run_2026-09-03/phase_D_li7_rank2.md` | 28 | | `run_2026-09-06/phase_C_survey.md` | 5 |
| `run_2026-09-03/phase_E_numbers.md` | 23 | | `run_2026-09-06/phase_B_numbers.md` | 3 |
| `run_2026-09-03/phase_D_small_items.md` | 20 | | `engineering.md` | 2 |
| `run_2026-09-03/phase_C_numbers.md` | 17 | | `run_2026-09-03/PLAN.md` | 2 |
| `run_2026-09-03/AUTHOR_DECISIONS.md` | 11 | | `run_2026-09-03/STATUS.md` | 2 |
| `run_2026-09-03/phase_A_miller_normalisation.md` | 11 | | `run_2026-09-03/phase_A_cdbonn.md` | 2 |
| `run_2026-09-02/design_C_tensor_rc.md` | 8 | | `run_2026-09-06/phase_CW_numbers.md` | 2 |
| `run_2026-09-03/phase_B_numbers.md` | 6 | | `run_2026-09-06/phase_C_numbers.md` | 2 |
| | | | *(3 more, 1 each: `phase_A_li6_tables.md`, `run_2026-09-02/PLAN.md`, `design_D_b1_li6.md`)* | 3 |

23 records in all: the ten heaviest carry 233, the other thirteen 33, none of
them more than 5. **233 + 33 = 266.**

### D2.4 Why four characters, measured

The same shape of measurement D1.6 used for `RECORD_WIN`: widening the rule
buys evidence and sells discrimination, so both are measured. The decoy test
points each citation at a block of the same length in the same file, shifted
clear of the true one — four decoys per citation, at ±(L+5) and ±(L+40) lines
— and a decoy that falls outside the file, or that has no non-blank line in it
at all, is dropped from the population rather than scored as a refusal the rule
can take credit for.

**The script is `validation/record_rule_sweep.py`**, committed with this phase
(§D3.6 (vi)): it was a scratch file when this section was written and its
numbers were published with "the script in `/tmp/…/sweep.py` shape" for a
recipe, which is not a recipe. Every number in the two tables below is

```
python3 validation/record_rule_sweep.py tokens [--tree <dir>] [--no-stoplist]
```

with `--tree` naming a `git archive a3c9ecb | tar -x` checkout carrying this
gate — the same arrangement §D1.3 uses. **Run with the gate as D2 landed it,
the script reproduces both tables below to the digit.** Run with the gate this
tree now carries, the citing sentence of a bulleted item is clipped to that
item (§D3.5) and the numbers shift a little — pre-fix **243 / 136 / 121**
passing at n = 3 / 4 / 5 and decoys **56.9 % / 26.4 % / 21.9 %**; after, **266
/ 266 / 250** and **59.8 % / 34.7 % / 29.2 %**. The shape, and the reason four
is the threshold, are the same either way.

**On the pre-fix tree (a3c9ecb), where the pass rate is the thing being
measured:**

| token length | citations passing | failing | decoy blocks accepted | true/decoy |
|---|---|---|---|---|
| 3 | 245 | 131 | 807 / 1403 = 57.5 % | 1.13 |
| **4** | **143** | **233** | **388 / 1403 = 27.7 %** | **1.38** |
| 5 | 128 | 248 | 329 / 1403 = 23.4 % | 1.45 |

**On the tree this phase leaves:**

| token length | citations passing | failing | decoy blocks accepted | true/decoy |
|---|---|---|---|---|
| 3 | 257 | 0 | 568 / 943 = 60.2 % | 1.66 |
| **4** | **257** | **0** | **326 / 943 = 34.6 %** | **2.89** |
| 5 | 233 | 24 | 276 / 943 = 29.3 % | 3.10 |

The stoplist of §D2.2 is measured the same way — the same loop with
`g.RELAX_STOP` replaced by `frozenset()` — and the four numbers it yields are
quoted there.

Three characters is not a rule: it accepts a decoy 60 % of the time and would
have passed 102 of the 243 broken citations unread. Five is very slightly
sharper per citation but refuses 24 anchors that are demonstrably right and
would buy those 24 exemptions for +0.21 of ratio. **Four, as specified, is the
right threshold and the measurement says so** — which is worth recording,
because it was not obvious in advance. Extending the numeric branch so that
`5.60152` matches `5.601518702` (a document rounding what the code writes out)
was measured too and is **not** done: it buys **0** citations on either tree
and costs +0.1 points of decoy acceptance.

### D2.5 What was repaired, and on what evidence

Every one of the 243 was read against its target. The repair split:

| treatment | citations |
|---|---|
| **re-pointed** — the thing the sentence names is still in the tree; only the line moved | **137** |
| **annotated `(as of <commit>)`** — it describes a state that no longer exists | **95** |
| **annotated `(anchor read 2026-09-15)`** — it lands right; the sentence says nothing a machine can match | **37** |
| *(of which both re-pointed and annotated)* | *(26)* |
| distinct citations touched: 137 + 132 − 26 | **243** |

No citation was added or deleted: the ordered citation list of every record is
the same length as at a3c9ecb, checked file by file.

**86 of the re-points are machine-verified relocations**, under this reading
and no other (the "75" that stood here named no reading and reproduced under
none — §D3.6 (v)):

> take every record citation whose target line differs between `a3c9ecb` and
> this tree; find the OLDEST commit that introduced that citation's text into
> that record (`git log -S`); take the block the citation named AT THAT COMMIT;
> and count it verified when that block's text occurs **exactly once** in the
> target file today and the new citation lands **exactly** on it.

```python
# paste into a file and run from the repository root; prints the four counts
import importlib.util, subprocess
from pathlib import Path
ROOT = Path(".").resolve()
spec = importlib.util.spec_from_file_location(
    "g", ROOT / "validation/check_physics_channels_links.py")
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
BASE = "a3c9ecb"

def git(*a):
    r = subprocess.run(["git", *a], capture_output=True, text=True, cwd=ROOT)
    return r.stdout if r.returncode == 0 else None

def cites(text):
    out = []
    cs = [(m, int(m.group(2)), int(m.group(2)))
          for m in g.RECORD_REF.finditer(text)]
    cs += [(m, int(m.group(2)), int(m.group(3))) for m in g.RANGE.finditer(text)]
    cs.sort(key=lambda c: c[0].start())
    for m, a, b in cs:
        p = m.group(1)
        if p is None:
            p, online = g.record_inherited_path(text, m.start())
            if p is None and not online:
                p = g.column_path(text, m.start())
        out.append((p, a, b, m.group(0)))
    return out

moved = []
for doc in g.record_files():
    rel = str(doc.relative_to(ROOT))
    old = git("show", f"{BASE}:{rel}")
    if old is None:
        continue
    a, b = cites(old), cites(doc.read_text())
    assert len(a) == len(b), rel          # no citation added or removed
    for (op, oa, ob, txt), (np, na, nb, _) in zip(a, b):
        if op == np and (oa, ob) != (na, nb):
            moved.append((rel, op, oa, ob, na, nb, txt))

ok = not_unique = elsewhere = 0
for rel, path, oa, ob, na, nb, txt in moved:
    log = git("log", "--format=%h", "-S", txt, "--", rel).split()
    src = git("show", f"{log[-1]}:{path}")
    lines = src.splitlines() if src else []
    if not 1 <= oa <= ob <= len(lines):
        continue
    block = "\n".join(lines[oa - 1:ob])
    today = (ROOT / path).read_text(errors="replace")
    if today.count(block) != 1:
        not_unique += 1
    elif today[:today.index(block)].count("\n") + 1 == na:
        ok += 1
    else:
        elsewhere += 1
print(len(moved), "moved;", ok, "verified;", not_unique, "no longer unique;",
      elsewhere, "landed elsewhere")
```

> `155 moved; 86 verified; 12 no longer unique; 57 landed elsewhere`

That is R4's MOVED logic sourced from git instead of from a sidecar, which is
the right way round for history — a record makes no promise that anyone else's
file will not change, so there is nothing to record. The 57 that "landed
elsewhere" are not errors: they are the re-points made by READING, where the
anchor was already stale when the record was committed (see the paragraph
below) or where the sentence's subject had moved to a different construct
entirely.

**The other 51 of the 137 were re-pointed by reading (69 of the 155 the committed script counts; the earlier "62" was 137 − 75, the withdrawn number)**, and the byte-identical test is
exactly why they had to be: applying it *without* the second condition proposed
188 relocations, and reading a sample showed a large minority of them were
relocating an anchor that was **already stale when the record was committed**
— `src/core/pipeline.cpp:672-674` (anchor read 2026-09-15) was moved to
`:740-742` (anchor read 2026-09-15) and neither block has
anything to do with `default_inclusive_kernel`, which is what its row is about.
A relocation is only evidence when the anchor was demonstrably right to begin
with.

**The historical annotations are machine-CHECKABLE, and 98 of the 107 are
machine-derived.** The candidate commits come from `git rev-list`; the commit
named is one at which the gate's own rule PASSES for that anchor and that
sentence, so `(as of 66dcda2)` is not an opinion — it is "at 66dcda2,
`docs/OPEN_ITEMS_SOLUTIONS.md:428` (as of 66dcda2) read *the gate's default
today | ratio 0.440 (factor 2.27) | outside G3b*, which is what this row says
it says". The bar is two shared tokens, or one of six characters or more; a
single short match is a coincidence and was sent to be read by hand instead
(22 of them were).

**Nine fail that bar at the commit they name, and were established by READING
instead** (§D3.6 (i); this paragraph said "machine-derived … not an opinion" of
all of them, which was overstated). Six are in the older records and are the
table below; the other three are in this record, where the row quotes what
`include/lipolgen/constants.hpp:86` (as of a3c9ecb) held — the b₁-normalisation
comment — which is the point of the row. Each was read at its commit and is
right there by content:

| annotation | at its commit the target reads | verdict |
|---|---|---|
| `run_2026-09-02/design_C_tensor_rc.md:2146` → `src/hepmc/hepmc_writer.cpp:135` (as of 0145885) | `const double gen_mass = (std::abs(p.pdg) == 11 && p.mass == 0.0) ? 0.51099895e-3 : p.mass;` | the one-line change to use `M_ELECTRON` is that line |
| `run_2026-09-03/phase_A_gate_mechanics.md:754` → `docs/OPEN_ITEMS_SOLUTIONS.md:414` (as of 66dcda2) | `**G3b** — peak magnitude within a factor 2 \| **FAIL**: …ratio 0.440…` | the row cites it as "G3b FAIL at :414" |
| `run_2026-09-03/phase_A_miller_normalisation.md:55` → `include/lipolgen/constants.hpp:86` (as of a98f0a0) | `inline constexpr double B1_PER_DEUTERON_TO_PER_NUCLEON = 0.5;` | the constant whose comment §6.3 rewrites |
| `run_2026-09-03/phase_A_miller_normalisation.md:413` → the same line | the same | the §6.3 heading itself |
| `run_2026-09-03/phase_C_numbers.md:2618` → `python/lipolgen/cli.py:604` (as of adec442) | `if p.rc_model is not None:` | the `rc_model` use site `PHYSICS_CHANNELS.md` cited there at adec442 |
| `run_2026-09-06/phase_C_survey.md:278` → `validation/o5_a2_reach.py:466-468` (as of 7f68339) | the `HfsModel(eff_track=0.95)` comment and `EPS_TRACK_SIBLING_STANDIN` | the lines the section quotes verbatim |

The suite now holds that list: `test_every_as_of_anchor_passes_the_bar_at_its_commit`
checks all 107 against the bar and carries exactly these nine as named
exceptions, so a tenth cannot appear unnoticed. The largest single family is
`run_2026-09-03/phase_A_gate_mechanics.md`'s survey of **where the ⁶Li
publication ban was written** — **31 of the 96**, reaching eleven documents
(`docs/OPEN_ITEMS_SOLUTIONS.md` 9, `docs/USAGE.md` 5, `include/lipolgen/
b1_nuclear.hpp` and `tests/test_b1_nuclear.cpp` 3 each, `README.md`,
`docs/PHYSICS_CHANNELS.md`, `python/bindings.cpp` and `python/lipolgen/cli.py`
2 each, three more 1 each), almost all of them describing text that the ban's
lift (`ac22331`) deleted. Re-pointing those would make the record claim something it
never claimed; dating them is the only honest repair, and it is the treatment
`PLAN.md` D2 asks for.

### D2.6 Five defects in the rule itself, each found by measuring

1. **A table row that names several files has no file column.** `column_path`
   took the last backticked path earlier in the row, and on
   `phase_CW_numbers.md`'s row naming `validation/vmc_reconcile.py` *and*
   `validation/vmc_tag_fractions.py` it resolved `:133-136` onto the second
   while the sentence meant the GENERATED
   `docs/open_items/vmc_reconciliation.md` — and reported a sound citation
   broken. A gate that invents a failure is worse than one that says it cannot
   tell, so an ambiguous row now says it cannot tell. Worth 2 false failures.
2. **The same-line rule read past a citation it could not resolve.** These
   records write `` (`tagged.hpp:255-257` → `:364-375`) `` with the file
   spelled without its directory, which `PATH` does not match; the main gate's
   `inherited_path` then walked further left and inherited a fully spelled path
   four citations earlier. Two sound citations in `run_2026-09-03/STATUS.md`
   were reported broken by exactly that. `record_inherited_path` stops at the
   first `x:123` on the line and uses it only if it is a path this gate knows —
   while still looking *past* a bare `:11-22`, which is another citation doing
   the same inheriting.
3. **Bare `` `:160` `` point citations matched no rule at all.** `REF` requires
   a path, so the point form of the inheritance these records use everywhere
   was checked by NOTHING: **152 of them in 16 records**. `RECORD_REF` is `REF`
   with the path made optional and is used by this pass only — widening `REF`
   itself would move the strict gate's own published counts. **The same hole
   exists in the three strictly gated documents, where it is worth 4
   citations; this item records it and does not close it.** (A "citation count
   to 655, and 29 further broken anchors" stood here until 2026-09-15: it was
   a mid-phase state of an uncommitted tree, reproducible from nothing, and it
   is withdrawn — §D3.6 (vii). What the bare form is worth is measurable from
   committed state and is in §D2.3's totals: 152 bare citations in 16 records
   at a3c9ecb.)
4. **An inline annotation keyed on the citing SENTENCE exempts its
   neighbours.** Written that way, one `(as of …)` silently covered every other
   citation within `RELAX_WIN` characters — most of a table row, four
   neighbours in the worst case measured. Both annotations now have to follow
   THEIR OWN citation, immediately or immediately after its adjacent
   `` `name` ``.
5. **`(as of …)` has to be read before the non-blank test.** Four of the 95
   name a line that has since gone blank; refusing them for it would be
   refusing the annotation for doing its job.

### D2.7 What R6 still cannot see

1. **181 of the 664 citations (27 %) carry no path a rule can recover** — the
   path is in the prose ("the design", "the header"), or in a table row that
   names several files. They are counted and named, never silently dropped, but
   they are not checked.
2. **A decoy is accepted 34.6 % of the time** (§D2.4) — much more often than
   R5's 14.3 %, and necessarily so: R5 gets a *quoted phrase* to work with and
   R6 gets one word. R6 catches gross drift, which is the failure mode that
   actually misleads a reader; it does not certify an anchor.
3. **`(anchor read <date>)` is as good as the reading behind it.** 39 inline
   and 51 under two section markers, 90 in all, rest on a human's say-so on one
   date. This item said a pin per citation "was considered and not done: at 88
   entries the sidecar would be larger than the thing it protects" — which
   confused two things, because a pin does not need a sidecar: it goes in the
   annotation, in the document, where the reader can see it. **§D3.4 puts one
   there for every POINT citation** (28 of these 90 were points, and each was
   free to drift onto any non-blank line of its target); the 62 range
   annotations still carry none, for the reason §D3.4 gives.
4. **`(as of <commit>)` says the anchor passed R6 at that commit, not that the
   sentence around it was true.** The suite checks the commit resolves; nothing
   checks that the reading of it was right.
5. **Nothing re-runs.** R6 is checked on every `--records` run, so unlike R5 it
   does not go stale — but a record that is never edited and whose targets keep
   moving will simply accumulate annotations, and the 137 re-points of this
   phase will need doing again the next time `src/core/pipeline.cpp` grows by
   500 lines. The honest reading of that is: **the re-points are worth less than
   the annotations**, and a future phase should consider dating a whole record
   at its own commit rather than chasing its anchors.

### D2.8 Suites, gates, tallies

Baseline is the tree this phase started from (a3c9ecb plus D1, measured in
§D1.8); the end state is measured here.

| | after D1 | after D2 |
|---|---|---|
| `build/lipolgen_tests` | 408 cases / 17 241 975 assertions / 0 failed | **408 / 17 241 975 / 0 failed** |
| `python -m pytest python/tests -q` | 1039 passed / 151 skipped | **1062 passed / 151 skipped** (+23, all in `test_doc_link_gate.py`) |
| docs gate (`PHYSICS_CHANNELS.md`) | 1236 refs / 96 ranges / 7 external / 0 broken / 6 allow-listed | **1236 / 96 / 7 / 0 broken / 6 allow-listed** |
| `[SPIN32_FINITE_GAMMA.md]` | 19 refs / 115 ranges / 0 broken | **19 / 115 / 0 broken** |
| `[PYTHIA_BRIDGE.md]` | 8 external / 0 broken | **8 external / 0 broken** |
| `--loose` | 1236 / 96 / 7 / 0 broken | **1236 / 96 / 7 / 0 broken** |
| R5 audit (`--audit-ranges`) | 180 blocks, 0 REFUSED | **180 blocks, 0 REFUSED** |
| R5 audit (`--audit-ranges --no-range-allow`) | 74 REFUSED | **74 REFUSED** |
| R6 (`--records-only`) | *(did not exist)* | **664 citations / 42 records / 0 broken** |
| `check_spdx_headers.py` | 102/102 | **102/102** |
| sidecar `validation/physics_channels_ranges.json` | 220 entries | **220 entries, untouched** |

The strict gate's output is **byte-for-byte what D1 left**, on purpose: R6 is
opt-in (`--records`, `--records-only`) and is suppressed in a recording pass,
so every recipe §D1.3 publishes with `| tail -1` still names the line it named.

The 23 new pytests are R6 on fixtures — the shared-token rule and its failure,
a token one character too short, a common English word refused as evidence, a
blank and an out-of-bounds target, a block needing only one non-blank line, a
record-to-record citation checked for existence only, a target outside the tree
skipped by name, `(as of …)` covering its own citation and NOT its neighbour,
`(as of …)` after an adjacent name, `(as of …)` over a target that has since
gone blank while `(anchor read …)` over the same target is still broken,
`(anchor read …)`, a section marker stopping at the next heading, table-column
inheritance and its ambiguous row, the bare-filename stop, a bare citation
inheriting past another bare one, a bare point citation reached at all, and the
sidecar untouched by the pass — plus five over the real records: the whole of
`docs/open_items/` passing, the pass reaching ≥ 240 checked citations and ≤ 200
skipped (a rule that quietly stopped resolving paths would pass by checking
nothing), every `(as of <commit>)` resolving with `git cat-file`, and every
section marker covering at least one citation.

**No file under `src/`, `include/`, `python/lipolgen/`, `tests/` or
`validation/` changed, and no citation was added or removed.** `git diff
--stat` for this item touches
`validation/check_physics_channels_links.py`,
`python/tests/test_doc_link_gate.py` and **22 documents under
`docs/open_items/`** (223 lines changed, all of them citation anchors and
annotations) — and this file.

---

## D3. The three structural gaps D1 and D2 left, closed

Performed 2026-09-15, on the tree D1 and D2 left (a3c9ecb plus §§D1–D2). No
file under `src/` or `include/` is touched by this section either — only the
gate script, its new sibling `validation/record_rule_sweep.py`, the pytests,
the two strictly gated documents and fifteen of the dated run records.

Two adversarial reviews of D1 and D2 found three holes, each of them a rule
that was narrower than the sentence describing it, and seven published numbers
that did not reproduce. §§D3.1–D3.5 close the holes; §D3.6 takes the seven
numbers one at a time.

### D3.1 The one sentence, three times

> **An UNNAMED point citation is now evidenced, not merely non-blank (rule
> B2): the citing sentence must quote something the cited LINE carries. Over
> the 42 such citations `a3c9ecb` holds, the old rule flagged 0 and let 76.2 %
> of one-line drifts through; this one flags 30 and lets 2.4 %. Reading all
> 30: 7 were stale and are re-pointed, 23 are right and are exempted by name,
> with a pin, in `UNNAMED_ALLOW`.**

> **R5's evidence is now PER CITATION, not per block: each citing sentence
> must evidence the block it points at, and `RANGE_ALLOW` is keyed on
> (document, citing line, range). Pooling had accepted three of five
> deliberately wrong-block ranges and had never judged 11 citations on 7 keys
> against their own sentence at all; all five attacks are now refused, all 41
> unjudged-or-newly-refused citations were read, and 2 of them were wrong.**

> **R6 holds a POINT citation to a token of six characters (or five with a
> digit) that the lines either side do not carry, and an `(anchor read)` on a
> point must carry the text the reader saw. Over the records this phase
> leaves, a +1 drift of a point citation was accepted 29.6 % of the time and
> is now accepted 12.2 %; 28 anchor-read points that were pinned by nothing
> now carry a pin; 35 citations the strengthened rule flagged were read and
> repaired.**

### D3.2 B2: an unnamed point citation must say something about its line

Rule B asserted that the cited line was "a line somebody could have meant" —
in bounds, non-blank, not a bare delimiter — and nothing else. That is a test
nearly every line of every file passes, and two citations had been stale for
weeks with the gate green:

| citation | points at | says it points at | actually at |
|---|---|---|---|
| `docs/PHYSICS_CHANNELS.md:301` | `docs/CONVENTIONS.md:357` (anchor read 2026-09-15 "always Hulthen") | "a_nn = −18.9 fm … appears only as an in-code constant and in `docs/CONVENTIONS.md`" | `docs/CONVENTIONS.md:383` (anchor read 2026-09-15 "a_nn = −18.9 fm") — `:357` has been the deuteron-control sentence since `cdd8591` |
| `docs/PHYSICS_CHANNELS.md:206` | `include/lipolgen/constants.hpp:86` (as of a3c9ecb) | "[SS90] bag-model sum-rule coefficient" | `include/lipolgen/constants.hpp:136-138` — line 86 has been the b₁-normalisation comment since before `bd775bc` |

B2 is R5's own rule with the window closed to zero: a phrase the citing
sentence QUOTES must occur on the cited line, as text or as a symbol. The
alternative — R6's shared-token test, which is what a first look suggests —
was measured beside it:

```
git archive a3c9ecb | tar -x -C /tmp/d3head
cp validation/check_physics_channels_links.py validation/record_rule_sweep.py \
   /tmp/d3head/validation/
python3 validation/record_rule_sweep.py unnamed --tree /tmp/d3head --drifts 1 -1 2 -2
```

| rule | flags | +1 | −1 | +2 | −2 |
|---|---|---|---|---|---|
| non-blank only (rule B before this phase) | 0 | 76.2 % | 76.2 % | 81.0 % | 92.7 % |
| **a QUOTED phrase on the line (as shipped)** | **30** | **2.4 %** | **4.8 %** | **2.4 %** | **4.9 %** |
| a shared token of 4+ on the line | 19 | 16.7 % | 21.4 % | 7.1 % | 26.8 % |
| either | 17 | 16.7 % | 21.4 % | 7.1 % | 29.3 % |

The token rule flags eleven fewer citations and lets seven times as many
one-line drifts through. The quoted-phrase rule is the one taken, and its cost
is 30 rows to read.

**The 42 → 30 → 7 split.** All 30 were read against their lines, one at a
time. **7 were stale** and are re-pointed by content; **23 are right** — the
sentence quotes something the line does not carry, usually because the document
writes a formula in Unicode where the code writes ASCII — and are exempted in
`UNNAMED_ALLOW` with a pin and a reason each.

| # | citation | was | is now | why |
|---|---|---|---|---|
| 1 | `PHYSICS_CHANNELS.md:173` | `src/core/pipeline.cpp:1016` | `:1110` | `:1016` is the unrelated `unpol_sf = toy` refusal; `:1110` is the guard that refuses a caller-supplied kernel on a tagged channel |
| 2 | `PHYSICS_CHANNELS.md:206` | `include/lipolgen/constants.hpp:86` (as of a3c9ecb) | `include/lipolgen/constants.hpp:136` (anchor read 2026-09-15 "Sather-Schmidt") | the line naming Sather-Schmidt PRD 42:1424, two above `C_BAG` |
| 3 | `PHYSICS_CHANNELS.md:299` | `python/lipolgen/cli.py:101` (anchor read 2026-09-15 "'ladder' (the DEFAULT)") | `python/lipolgen/cli.py:139` (anchor read 2026-09-15 `add_argument("--p-d"`) | `--p-d` is declared there; line 101 is the `--pzz-mode` help |
| 4 | `PHYSICS_CHANNELS.md:301` | `docs/CONVENTIONS.md:357` (anchor read 2026-09-15 "deuteron control channel") | `docs/CONVENTIONS.md:383` (anchor read 2026-09-15 "KAPPA_NN_VIRTUAL") | where `a_nn = −18.9 fm` is |
| 5 | `PHYSICS_CHANNELS.md:417` | `tests/test_rc.cpp:496` (anchor read 2026-09-15 "derived by SUBTRACTING") | `tests/test_rc.cpp:518` (anchor read 2026-09-15 "the first C0 zero lies") | the SUBCASE that gates q₀ ∈ [2.9, 3.3] fm⁻¹; `:496` is a comment in the ⟨r²⟩ subcase |
| 6 | `SPIN32_FINITE_GAMMA.md:260` | `src/core/bookkeeping.cpp:89` (anchor read 2026-09-15 "octupole moment R_3") | `src/core/bookkeeping.cpp:94` | the sentence quotes the CALL `spin32_populations(pz, opt.pzz)`, which is there; `:89` is the comment, and SPIN32_FINITE_GAMMA.md:983 cites it rightly |
| 7 | `SPIN32_FINITE_GAMMA.md:795` | `src/core/xsec.cpp:41` | `src/core/sf.cpp:52` | `r_sigma_lt` is not in `xsec.cpp` and never has been; `sf.cpp:52` is `return 0.18 / (1.0 + q2 / 50.0);`, which is what "is not zero" means |

Three of the seven then evidence themselves (the sentence quotes `--p-d`,
`spin32_populations(pz, opt.pzz)` and `r_sigma_lt`, and the new lines carry
them); the other four are exempted like the 23. The tree this phase leaves
carries **43** unnamed point citations — the 42 plus the one §D3.3 created —
**28 of them exempted**, each printed on every run:

```
python3 validation/check_physics_channels_links.py | grep -c '(unnamed) --'
```

> `28`

### D3.3 R5 is per citation

`--record-ranges` refused a block "the citing sentence" did not evidence, but
the verdict was pooled on the BLOCK: *a block cited by two documents is
evidenced if EITHER sentence evidences it*, and `RANGE_ALLOW` was keyed on the
range. So a citation re-pointed onto any of the 180 already-blessed blocks was
accepted whatever its own sentence said. Five wrong-block ranges were built to
test it, and three passed.

The verdict and the exemption are now keyed on **(document, citing line,
range)**, and a block is fingerprinted only when EVERY citation of it is
evidenced or exempted. The five attacks, re-run on a copy of this tree, one
edit at a time (`--audit-ranges`):

| # | the citation moved | to | before | now |
|---|---|---|---|---|
| W1 | `PHYSICS_CHANNELS.md:302`'s triton-SF sentence | `docs/CONVENTIONS.md:18-21`, never blessed | REFUSED | **REFUSED** |
| W2 | the same | `docs/CONVENTIONS.md:194-198` (anchor read 2026-09-15), an EXEMPTED block | accepted | **REFUSED** |
| W3 | the same | `docs/CONVENTIONS.md:66-74` (anchor read 2026-09-15), evidenced by another citation | accepted | **REFUSED** |
| W4 | `PHYSICS_CHANNELS.md:78`'s θ_S/φ_S sentence | `include/lipolgen/bookkeeping.hpp:98-100` (anchor read 2026-09-15) | REFUSED | **REFUSED** |
| W5 | `SPIN32_FINITE_GAMMA.md:426`'s b₁ > 0 sentence | `docs/CONVENTIONS.md:508-509`, the luminosity-shares block | accepted | **REFUSED** |

and two controls in the same shape — the same citation moved onto a block its
own sentence DOES evidence (`include/lipolgen/bookkeeping.hpp:175-176`, which
declares the quoted `tensor_thirds_plan`; `include/lipolgen/spin.hpp:58-61`,
which declares the quoted `m_values`) — are **ACCEPTED, 0 REFUSED**, so the

> **Control C2 WITHDRAWN (2026-09-15):** SPIN32:426's sentence quotes no `m_values`, and `spin.hpp:58` / `:61` are blank, so R2 refuses the range before R5 ever judges it — the "0 REFUSED" was R5 never seeing it. A valid second control is PHYSICS_CHANNELS.md:302 → `include/lipolgen/constants.hpp:48-52`, quoting `HBARC_GEV_FM`: ACCEPTED, 0 REFUSED.
rule has not simply become a wall.

**The 11 that had never been judged.** Seven range keys were cited by more than
one sentence with only one of them evidencing; the other 11 citations rode in
free. Those 11, plus the 30 citations the new per-citation keying and §D3.5's
minimum evidence newly refuse, are 41 rows that were read one at a time. **39
are right** and are exempted per citation; **2 were wrong**:

| citation | was | is now | why |
|---|---|---|---|
| `PHYSICS_CHANNELS.md:974` | `docs/USAGE.md:1501-1509` ("§12") | `docs/USAGE.md:2717-2724` ("§8") | the sentence is about the four C++ example generators and their flag spellings; `:1501-1509` is the Tier T1 section, and USAGE.md has no §12. This is the citation the previous session's notes flagged as "probably stale by content" and asked this phase to catch |
| `PHYSICS_CHANNELS.md:339` | `docs/open_items/physics_literature.md:81-85` | `:234` (now a point citation) | the claim is that the register calls the per-nucleon Glauber product the citable primary and the coherent-cluster amplitude the systematic variant; that sentence is at `:234`. `:81-85` is the W_FSI formula and its starting parameters — which is what the OTHER citation of that range, the one at PHYSICS_CHANNELS.md:343 for the 30 → 70 mb Deeps fit, is right about |

`RANGE_ALLOW` grows from 74 range keys to **103 citation keys**: 64 carried
over verbatim (one citation, one key, one reason), 39 written after reading.

### D3.4 R6's point rule, measured

A record's point citation names ONE line, and R6 asked only for one shared
token of four characters — a test the neighbouring lines usually pass too, so
the drift it exists to catch was mostly invisible. Thirty candidate rules were
measured over the live-file point citations of this tree:

```
python3 validation/record_rule_sweep.py drift --only points
```

**Measured on the tree BEFORE §D3 was written (98 live-file point citations); §D3 added twelve more (110 on the tree this phase leaves). The current-tree figures, same recipe (`record_rule_sweep.py drift --only points`): D2 rule 29.1 / 23.6 / 23.6 / 25.5 %, decoy 16.1 %; shipped rule 12.7 / 8.2 / 15.5 / 12.7 %, decoy 9.3 %; "≥ 1 token of 6+" 106 passes / 4 breaks. The conclusion — cheapest rule with both ±1 under 20 %, zero honest breaks — holds on both.**

| rule (of the thirty the script carries; pre-§D3 tree, 98 points) | passes | +1 | −1 | +2 | −2 | decoy |
|---|---|---|---|---|---|---|
| R6 as landed in D2: ≥ 1 token of 4+ | 98 | 29.6 % | 23.5 % | 23.5 % | 24.5 % | 15.7 % |
| ≥ 1 token of 6+ | 94 | 19.4 % | 15.3 % | 17.3 % | 19.4 % | 10.5 % |
| ≥ 1 token of 7+ | 83 | 15.3 % | 11.2 % | 14.3 % | 11.2 % | 7.0 % |
| ≥ 2 tokens of 4+ | 50 | 11.2 % | 11.2 % | 2.0 % | 11.2 % | 4.4 % |
| ≥ 2 tokens, or one of 7+ | 89 | 17.3 % | 14.3 % | 15.3 % | 15.3 % | 8.7 % |
| a QUOTED phrase on the line (B2's rule) | 69 | 10.2 % | 6.1 % | 9.2 % | 6.1 % | 3.5 % |
| one token off both adjacent lines | 98 | 15.3 % | 14.3 % | 17.3 % | 16.3 % | 13.4 % |
| one token of 6+ off both adjacent lines | 94 | 11.2 % | 8.2 % | 13.3 % | 12.2 % | 8.5 % |
| **6+, or a NUMBER of 5+, off both adjacent lines (as shipped)** | **98** | **12.2 %** | **8.2 %** | **14.3 %** | **12.2 %** | **8.7 %** |
| 6+, or a number of 4+, off both adjacent lines | 98 | 12.2 % | 9.2 % | 14.3 % | 12.2 % | 9.0 % |

(the full thirty are in the script's output; these are the families). The
shipped rule is the cheapest that gets both ±1 rates under 20 %: on the tree it
leaves it refuses **none** of the 98 honest point citations, and it takes the
+1 drift acceptance from 29.6 % to 12.2 % and the decoy rate from 15.7 % to
8.7 % on the pre-§D3 tree (29.1 % → 12.7 % and 16.1 % → 9.3 % on the tree this phase leaves). Allowing a five-character NUMBER is worth four honest citations — a
record's evidence for a line is often `1.848`, `0.001` or `0.05`, and those are
more specific than most words.

**A RANGE keeps the one-token rule**, and the measurement is why: the same
sweep with `--only ranges` puts a block's +1 drift acceptance at **89.0 %**
under the one-token rule (156 ranges on the tree this phase leaves; 91.0 % over 145 on the pre-§D3 copy), and the cheapest candidate that brings it under 50 %
(three shared tokens, 37.4 %) refuses **93 of the 156** honest range citations (86 of 145 on the pre-§D3 copy).
A block shifted by a line still contains what the sentence says it does; there
is no cheap rule that says otherwise, and pretending there is would buy a
hundred exemptions for nothing.

**`(anchor read <date>)` on a point citation now carries a pin** — the text the
reader saw on that line, in `"straight quotes"` or `` `backticks` ``, whichever
the line does not itself contain. Without it the annotation exempted the
citation from the token rule and left only "non-blank", so the anchor was free
to drift onto any non-blank line of the file: **28 point citations were in
exactly that state**, and all 28 now carry a pin taken from the line as it
stands in this tree — which is the tree the 2026-09-15 reading was done
against. A section marker covers RANGES only; a point under one is refused by
name until it carries its own annotation.

### D3.5 R5's minimum evidence

Three evidence forms were worth less than they looked, and the decoy test says
by how much (`python3 validation/record_rule_sweep.py r5 --tree /tmp/d3head`,
the same a3c9ecb checkout as §D3.2):

| rule | evidenced | breaks | decoys accepted |
|---|---|---|---|
| as landed in D1 (no minimum evidence) | 97 | 0 | 105 / 762 = 13.8 % |
| + a list item is its own paragraph | 96 | 1 | 102 / 762 = 13.4 % |
| + a one-word text match on ≤ 5 lines of the file | 91 | 6 | 93 / 762 = 12.2 % |
| + on ≤ 3 lines instead | 90 | 7 | 89 / 762 = 11.7 % |
| **+ no text match for a block over 150 lines (as shipped)** | **90** | **7** | **93 / 762 = 12.2 %** |
| + over 60 lines instead | 88 | 9 | 93 / 762 = 12.2 % |

1. **A markdown LIST ITEM is now its own citing sentence**, as a table row
   already was. A bulleted list carries no blank lines, so the whole list was
   one paragraph and a phrase quoted in the NEXT bullet evidenced this bullet's
   block. That is not hypothetical: it is exactly how
   `PHYSICS_CHANNELS.md:974` came to be evidenced — by the word `tier`,
   belonging to the bullet below it — while pointing at the wrong block
   (§D3.3). **The rule costs exactly one citation over the three covered
   documents of a3c9ecb, and that citation is the wrong one.**
2. **A one-word text match must occur on at most five lines of the target
   file.** `tier` is on 17 lines of `docs/USAGE.md` and `pzz` on 39 of
   `python/bindings.cpp`: a word that common cannot tell one block of a file
   from another, and the evidence for a dozen blocks was exactly one such word
   — `ratio`, `rates`, `applies`, `seed`, `tier`, `bound`, `pzz`, `radial`,
   `isotope`, `headroom`. **On a3c9ecb the rule costs five citations**
   (`PHYSICS_CHANNELS.md` at 207, 263 and 512, `SPIN32_FINITE_GAMMA.md` at
   1009 and 1024), each read and exempted. A multi-word phrase is not
   restricted: it is already specific. Nor is a SYMBOL match, which is the
   stronger branch — the word has to be declared or used in the block as code.
3. **A text match does not evidence a block of more than 150 lines.** By the
   decoy test, text evidence accepts a decoy 11.2 % of the time for a block of
   1–20 lines, 18.6 % for 20–60, 38.5 % for 60–150 and **50.0 % for 150+**: in
   a 371-line block it is a coin toss. Three blocks are that long
   (`docs/USAGE.md:163-533`, `:1933-2124`, `:2744-2941`); **the rule costs one
   citation on a3c9ecb** — `PHYSICS_CHANNELS.md:448`, exempted with a reason —
   because the other two are carried by a SYMBOL the block actually declares
   (`--x-max` and `qe_suppression`) and the symbol branch is not restricted.
   60 lines would be defensible too and costs two more exemptions for no
   measurable change in the decoy rate, so 150 is where the line is drawn.

### D3.6 What D2 published that did not reproduce

Seven numbers. Each is stated as it stood, then as it is.

| # | published | actually |
|---|---|---|
| (i) | "the 95 historical annotations are machine-derived … not an opinion" | nine of the 107 fail the bar at the commit they name and were established by READING (six of them written by D2, three by this section); they are named in §D2.5, and `test_every_as_of_anchor_passes_the_bar_at_its_commit` now holds the list |
| (ii) | "four of the 95 land on a line that has since gone blank" | **two** under R6's own (i); four have a blank first-or-last line, which is R2's criterion and R2 does not run here (§D2.2) |
| (iii) | a per-target table summing to 203 against a headline of 243 | re-derived under the D3 rule: 48 targets, 266 broken, 198 + 68 (§D2.3) |
| (iv) | "9 more records, ≤ 4 each" | 13 more, none over 5 — the table is re-derived in §D2.3 |
| (v) | "75 of the 137 re-points are machine-verified" | no reading was stated and none reproduces it; under the reading now written out in §D2.5, with the script beside it: **155 moved, 86 verified, 12 no longer unique, 57 landed elsewhere** |
| (vi) | D2.4's decoy table, recipe "the script in `/tmp/…/sweep.py` shape" | the script is committed as `validation/record_rule_sweep.py`; every D2.4 number reproduces from it to the digit under the gate D2 landed, and §D2.4 now carries the numbers the gate this tree leaves gives as well |
| (vii) | "it raised the citation count to 655 … and found 29 further broken anchors" | a mid-phase state of an uncommitted tree; withdrawn (§D2.6) |

### D3.7 What these three rules still cannot see

1. **An exemption is still a signature, not a proof.** 103 `RANGE_ALLOW` and
   28 `UNNAMED_ALLOW` entries say "a human read this row on 2026-09-15". The
   pin keeps the exemption on the block it was written for; nothing keeps the
   reading honest.
2. **The exemption key is the citing LINE**, so editing a covered document
   retires its exemptions — loudly ("no range citation reaches it"), which is
   the intended behaviour, but it is a running cost: a paragraph inserted above
   a table re-keys every exemption below it. The alternative (an ordinal, or a
   hash of the sentence) either slides silently onto a different sentence or
   retires on every wording change, so the line number is the least-bad key.
3. **A range in a record still drifts invisibly.** +1 acceptance for a block is
   89.0 % (156 ranges, this tree) and no measured rule brings it down at a price worth paying (§D3.4).
4. **67 `(anchor read)` RANGE annotations still carry no pin (17 inline plus 51 under the two markers, minus one counted elsewhere; the 48 anchor-read POINTS are pinned; the earlier 62 was D2's 11 + 51)** — the pin was
   made compulsory for points only, where a drift is unambiguous.
5. **B2 reads what the sentence QUOTES**, so a document that names its target
   in unquoted prose gets an exemption rather than a check: 28 of the 43.

### D3.8 Suites, gates, tallies

Baseline is the tree this section started from — a3c9ecb plus §§D1–D2, the
state §D2.8 measured.

| | after D2 | after D3 |
|---|---|---|
| `build/lipolgen_tests` | 408 cases / 17 241 975 assertions / 0 failed | **408 / 17 241 975 / 0 failed** |
| `python -m pytest python/tests -q` | 1062 passed / 151 skipped | **1085 passed / 151 skipped** (+23, all in `test_doc_link_gate.py`, which goes from 83 to 112 tests) |
| docs gate, primary document | 1236 refs / 96 ranges / 7 external / 0 broken / 6 allow-listed | **1237 / 95 / 7 / 0 broken / 24 allow-listed** |
| `[SPIN32_FINITE_GAMMA.md]` | 19 refs / 115 ranges / 0 broken / 0 allow-listed | **19 / 115 / 0 broken / 10 allow-listed** |
| `[PYTHIA_BRIDGE.md]` | 8 external / 0 broken | **8 external / 0 broken** |
| `--loose` | 1236 / 96 / 7 / 0 broken | **1237 / 95 / 7 / 0 broken** |
| R5 audit (`--audit-ranges`) | 180 blocks, 0 REFUSED | **209 citations on 179 blocks, 0 REFUSED** |
| R5 audit (`--audit-ranges --no-range-allow`) | 74 REFUSED | **103 REFUSED** |
| R6 (`--records-only`) | 664 citations / 42 records / 0 broken | **714 citations / 42 records / 0 broken** (266 against a live file, 35 into another record, 107 historical, 115 read, 191 skipped) |
| `check_spdx_headers.py` | 102/102 | **103/103** (the new sweep script) |
| sidecar `validation/physics_channels_ranges.json` | 220 entries | **219 entries** |
| `RANGE_ALLOW` | 74, keyed on the range | **103, keyed on (document, citing line, range)** |
| `UNNAMED_ALLOW` | *(did not exist)* | **28, keyed on (document, citing line, citation)**, and an entry no citation reaches is reported like a dead fingerprint |

**Why the strict gate's numbers moved, line by line.** The reference count goes
1236 → 1237 and the range count 96 → 95 because one range citation became a
point citation (§D3.3's second re-point). The allow-listed count goes 6 → 24
because rule B2's exemptions are printed beside rule A's, as they must be: an
exemption you cannot see is a hole. The sidecar loses one entry, the block
`docs/USAGE.md:1501-1509` (anchor read 2026-09-15), the Tier T1 section, which
nothing cites any more. **No fingerprint
changed**: the one recording pass this section ran reported *0 new, 0 changed,
1 dropped*, and a second pass reports *0 new, 0 changed, 0 dropped* and leaves
the file byte-identical.

The new pytests are the three rules and their edges: rule B2 evidenced by text
and by symbol, refused when the phrase is on the line next door, exempted with
a pin, the exemption keyed on the citation and not on the line, and skipped in
`--loose`; R5 per citation (a second citation of an evidenced block refused; an
exemption that covers one citation and not the next; a moved sentence retiring
its exemption loudly); R5's minimum evidence (a common word, a very long block,
a list item as its own sentence); R6's point rule (a token the neighbours do
not carry, a five-character number that counts, a range keeping the one-token
rule); the anchor-read pin (missing, slid off, and written in backticks); and a
dead `UNNAMED_ALLOW` entry reported — plus four over the real tree: every unnamed exemption live and still pinned,
every unnamed citation evidenced or exempted, every anchor-read pin still on
its line, and every `(as of …)` anchor passing the bar at its commit or named
in the by-hand list.

**No file under `src/`, `include/`, `python/lipolgen/`, `tests/` or
`validation/reference/` changed.** `git diff --stat` for this section touches
the gate, the new `validation/record_rule_sweep.py`, the pytests, the sidecar,
`docs/PHYSICS_CHANNELS.md`, `docs/theory/SPIN32_FINITE_GAMMA.md`, fifteen dated
run records, and this file.

## D4 Residues recorded at the close of the phase (2026-09-15)

1. **Anchor-read POINT pins are checked for presence, not uniqueness** (`check_physics_channels_links.py`, the `READ_INLINE` branch), unlike `RANGE_ALLOW`/`UNNAMED_ALLOW` pins, which go through `pin_state`. Three pins occur on more than one line of their target: `run_2026-09-03/phase_D_li7_rank2.md` line 552 → `validation/dump_polligen_reference.py:491` (anchor read 2026-09-15 "0.05 * f1"), a pin that occurs on 5 lines; `run_2026-09-06/phase_CW_numbers.md` line 699 → `README.md:109` (anchor read 2026-09-15 "rtol 1e-12"), a pin that occurs on 2 lines (drifted onto README.md's line 40 — the other line carrying that pin — the gate would report 0 broken); record line 1011 → `docs/CONVENTIONS.md:357` (anchor read 2026-09-15 "deuteron control channel"), a pin that occurs on 2 lines. None sits on an adjacent line, so ±1 is caught; a longer drift onto the twin line is not.
2. **A neighbour's annotation text is evidence:** `relax_sentence` blanks citations but not the `(anchor read … "pin")` text beside them, so two bare live citations in the §D3.2 table (its record lines 1012 and 1013, which cite `tests/test_rc.cpp:496` (anchor read 2026-09-15 "derived") and `src/core/bookkeeping.cpp:89` (anchor read 2026-09-15 "octupole")) — pass the point rule only on the pin text of the annotated citation next to them. Bare `:N` points are likewise not blanked (29 live windows carry one; no live citation currently rests on such a digit token).
3. Confirmed as designed and still true: an unpinned `(anchor read)` RANGE drifted +1 is accepted; an `(as of)` point onto a blank line is accepted; a range +1 is accepted 89.0 %.

Not fixed in this phase; the close-out's review will re-judge them.
