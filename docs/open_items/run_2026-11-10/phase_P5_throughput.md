<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-10 — phase P5: the filed LO item and the throughput (decisions J8, J9)

Written by the supervising session, 2026-11-17, on base commit `2bd78e2`
(phases P1–P4 of this run committed). It implements **J8 and J9 — the
SUPERVISOR's choices under the author's delegation** ("for the others go
with your recommendation"; `DECISIONS.md` J8, J9), never the author's words.
No source, test, default or reference number changed in this phase.

## 0. What changed, in one paragraph

**J8** (records only): the two places where the leading-order evaluation of
the NLO sets is filed as a development item — docs/USAGE.md ("What it is,
and what it is not", the D7 paragraph) and docs/benchmarking/08_channel_priority.md
(reading note 3) — each gain one dated clause on the same line: the item
stays filed and is not scheduled, and its symptom is bounded since A3.
**J9**: USAGE §8's throughput line gains, on the same line, the new recipe
(2 M events per run, `--out ""`, a quiet machine, a fresh process per run,
the median of 3, the raw log kept) and the numbers measured with it; the
200 k numbers stay as history.

## 1. Every number or label that moved, before → after

| where | before | after |
|---|---|---|
| USAGE §8, the numbers to quote | A16's 200 k medians 0.54 / 0.33 / 1.07 M ev/s (2026-11-06), with 0.61 / 0.33 / 1.05 at 2 M events as a supplement | **inclusive 0.605, tagged ⁶Li α 0.321, coherent 1.01 M ev/s** (2 M events per run, medians of 3, 2026-11-17 at `2bd78e2`) |
| USAGE §8, the recipe | 200 k events per run | 2 M events per run |
| the LO development item (two places) | "the development item D7 names" | the same, plus: stays filed, not scheduled (J8); the symptom, the pair's A₁ beyond 1 from x ≳ 0.83, is bounded since A3 by the clip of g₁ to F₁ on ~2–5 × 10⁻⁷ of the rate (copied from `../run_2026-11-03/phase_P2_positivity.md` §1; this run's phase P4 adds that the share is a cell-centre estimate and gives its bracket) |

### J9 — the measurement

Script (session scratch, `run_1110/p5/run_throughput_2M.sh`): for each of the
three examples, three runs, each a fresh process, `--events 2000000 --out ""`,
the example's own "events/s" line taken, the median of the three is the
number. The three commands are in the log below. Build: `build/` as phase P4
left it (no source changed since); nothing else of this session running.

Raw log, 2026-11-17 (the header, the nine runs with the 1 / 5 / 15-minute
load at each start, the footer; the examples' banners omitted):

```
# J9 throughput, 2026-11-17T18:14:29Z, HEAD 2bd78e2, 4 cores
# load before: 0.51 1.30 1.67 2/717 135587
# busy processes (>20 % CPU): python3(99.4) 
## inclusive run 1: ./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 2000000 --out  (load 0.51 1.30 1.67)
  time         3.286 s -> 6.09e+05 events/s
## inclusive run 2: ./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 2000000 --out  (load 0.55 1.29 1.66)
  time         3.323 s -> 6.02e+05 events/s
## inclusive run 3: ./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 2000000 --out  (load 0.55 1.29 1.66)
  time         3.305 s -> 6.05e+05 events/s
## tagged run 1: ./build/generate_tagged --channel 6Li-alpha --optics tagging --events 2000000 --out  (load 0.58 1.29 1.66)
  time         5.990 s -> 3.339e+05 events/s
## tagged run 2: ./build/generate_tagged --channel 6Li-alpha --optics tagging --events 2000000 --out  (load 0.65 1.28 1.65)
  time         6.388 s -> 3.131e+05 events/s
## tagged run 3: ./build/generate_tagged --channel 6Li-alpha --optics tagging --events 2000000 --out  (load 0.68 1.27 1.65)
  time         6.228 s -> 3.211e+05 events/s
## coherent run 1: ./build/generate_coherent --config 1 --slope 50 --f0 0.04 --events 2000000 --out  (load 0.70 1.27 1.64)
  time         2.001 s -> 9.997e+05 events/s
## coherent run 2: ./build/generate_coherent --config 1 --slope 50 --f0 0.04 --events 2000000 --out  (load 0.73 1.26 1.64)
  time         1.887 s -> 1.06e+06 events/s
## coherent run 3: ./build/generate_coherent --config 1 --slope 50 --f0 0.04 --events 2000000 --out  (load 0.73 1.26 1.64)
  time         1.978 s -> 1.011e+06 events/s
# load after: 0.75 1.26 1.64 1/710 135668
```

Medians: inclusive 6.05e5 (6.09 / 6.02 / 6.05), tagged ⁶Li α 3.211e5
(3.339 / 3.131 / 3.211), coherent 1.011e6 (0.9997 / 1.06 / 1.011) events/s.
Spread of the three runs (max − min over the median): 1.2 %, 6.5 %, 6.0 %.
The machine is shared with other sessions: the 1-minute load was 0.5 at the
start and 0.75 at the end (the single-threaded example itself accounts for
most of that), and the script's start-of-run sample saw one other process
busy (a `python3`, short-lived, not this session's). "Quiet", not provably
idle.

A first reading with the same script on 2026-11-12 at `8b72cd5` (after phase
P2b; a build had finished just before, 1-minute load 2.4 decaying) is kept
for comparison, not quoted in USAGE:

```
# J9 throughput, 2026-11-12T22:23:58Z, HEAD 8b72cd5, 4 cores
# load before: 2.45 0.87 0.34 1/692 74523
# busy processes (>20 % CPU): 
## inclusive run 1: ./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 2000000 --out  (load 2.45 0.87 0.34)
  time         3.211 s -> 6.23e+05 events/s
## inclusive run 2: ./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 2000000 --out  (load 2.50 0.91 0.35)
  time         3.251 s -> 6.15e+05 events/s
## inclusive run 3: ./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 2000000 --out  (load 2.38 0.91 0.36)
  time         3.162 s -> 6.32e+05 events/s
## tagged run 1: ./build/generate_tagged --channel 6Li-alpha --optics tagging --events 2000000 --out  (load 2.27 0.91 0.36)
  time         6.181 s -> 3.236e+05 events/s
## tagged run 2: ./build/generate_tagged --channel 6Li-alpha --optics tagging --events 2000000 --out  (load 2.17 0.91 0.36)
  time         6.187 s -> 3.232e+05 events/s
## tagged run 3: ./build/generate_tagged --channel 6Li-alpha --optics tagging --events 2000000 --out  (load 2.07 0.91 0.37)
  time         6.580 s -> 3.04e+05 events/s
## coherent run 1: ./build/generate_coherent --config 1 --slope 50 --f0 0.04 --events 2000000 --out  (load 2.13 0.97 0.39)
  time         2.011 s -> 9.943e+05 events/s
## coherent run 2: ./build/generate_coherent --config 1 --slope 50 --f0 0.04 --events 2000000 --out  (load 2.13 0.97 0.39)
  time         2.048 s -> 9.767e+05 events/s
## coherent run 3: ./build/generate_coherent --config 1 --slope 50 --f0 0.04 --events 2000000 --out  (load 2.04 0.97 0.39)
  time         2.053 s -> 9.74e+05 events/s
# load after: 1.95 0.97 0.40 1/700 74605
```

Medians then: 6.23e5 / 3.232e5 / 9.767e5. The 2026-11-17 medians against
them: −2.9 %, −0.6 %, +3.5 % — inside the run-to-run spread; phases P3 and P4
changed no code on the event path of these three examples (P3: harnesses and
tests; P4: one `GenRunInfo` attribute written once per file, a virtual read
once per run, comments). Against run 2026-11-03's 2 M-event supplement
(0.61 / 0.33 / 1.05, 2026-11-06 at `888ca3b`): −0.8 %, −2.7 %, −3.7 %. No
change is attributed to this run's code; none was A/B-measured in this phase.

## 2. What was verified and how

* Line neutrality: `git diff --numstat` — docs/USAGE.md 2 added / 2 removed,
  docs/benchmarking/08_channel_priority.md 1 / 1 (USAGE §8 is cited by line
  ranges; the three citation modes after the edits are in §4).
* The J8 clause states only what is recorded: the bound and the 2–5 × 10⁻⁷
  are A3's (`../run_2026-11-03/phase_P2_positivity.md` §1, "So the clip
  binds on ~2–5 × 10⁻⁷ of the rate"); "nothing ranked 1–7 waits on it" is
  J8's own reason, quoted from `DECISIONS.md`.
* The medians were recomputed from the log lines above.

## 3. Files

docs/USAGE.md (two lines, in place); docs/benchmarking/08_channel_priority.md
(one line, in place); this record; STATUS.md (the P5 row, P4's hash).

## 4. Gates

| gate | result |
|---|---|
| build, `lipolgen_tests`, the full pytest, the bench subset, the report | not run: no source, test, harness or report input changed in this phase (their state is phase P4's: 440/440 C++, 1534 passed / 156 skipped, bench 372 / 1, report 15 / 8 / 0) |
| citations | strict 1243 references / 95 ranges / 7 external / 0 broken / 24 allow-listed; `--records-only` 1432 citations in 98 records, 0 broken; `--audit-ranges` 209, 0 REFUSED |
| SPDX; doc-link gate; `git diff --check` | 136 files; 112 passed (30 s); clean |

## 5. Open points

1. The throughput numbers are machine- and load-dependent at the few-percent
   level (the three runs of one example spread by up to 6.5 %); a change
   smaller than that between two dates is not a finding.
2. The machine is no longer exclusively this session's; a strictly idle
   measurement would need the other sessions stopped — not done.
