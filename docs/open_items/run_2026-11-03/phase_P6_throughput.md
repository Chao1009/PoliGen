<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-03 — phase P6: the examples' throughput re-measured on the idle machine (decision A16)

Supervisor, 2026-11-06 (local; 2026-11-07 UTC), on `claude/run-2026-10-26` at
`888ca3b` (phases P1–P5 committed; this phase uncommitted, committed by the
supervising session). Decision: `DECISIONS.md` A16 — the SUPERVISOR's choice
under the author's delegation ("Go with your recommendations for the 17
items, just record the choice and why"), not the author's words. Every number
below was measured in this phase on this machine (4 cores); the script and
the raw logs are in the session scratch `run_1103/p6/` (`run_throughput.sh`,
`throughput.log`, `ab_clip_build_and_run.sh`, `ab_clip.log`); the raw log is
also copied verbatim in §1 so the record stands without the scratch dir.

## 0. What was done, in one paragraph

The three example generators of USAGE §8 were run with the same recipe as the
2026-10-30 measurement (`--out ""`, 200 k events, a fresh process per run,
three runs per channel, the median is the number) — but on the idle machine
(no gate, build or worker running; the load averages logged before every run
and after the last) and with the raw log kept, which is what A16 asked for.
The 2026-10-30 values were taken while F7's gate runs shared the machine and
their log was not kept. Nothing else changed: no source, no default, no
test. USAGE §8's sentence gains a dated A16 clause on its own line (the
section is cited by line ranges, so the edit is line-neutral) and keeps the
2026-10-30 numbers as history.

## 1. The measurement

Commands (the examples' own streaming rate; each prints `time … s -> …
events/s` for the generation loop alone):

```
./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 200000 --out ""
./build/generate_tagged    --channel 6Li-alpha --optics tagging --events 200000 --out ""
./build/generate_coherent  --config 1 --slope 50 --f0 0.04 --events 200000 --out ""
```

Raw log (`throughput.log`, verbatim; the load triple before each run is the
1/5/15-minute average, the 1-minute value still decaying from the bench
subset that finished about four minutes earlier — no process of note was
running):

```
# P6 throughput, 2026-11-07T04:20:27Z, HEAD 888ca3b, 4 cores
# load before: 0.51 0.91 0.97 1/671 409768
# other processes of note: 0 (expect 0)
## inclusive run 1: ./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 200000 --out  (load 0.51 0.91 0.97)
  time         0.373 s -> 5.36e+05 events/s
## inclusive run 2: ./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 200000 --out  (load 0.51 0.91 0.97)
  time         0.371 s -> 5.4e+05 events/s
## inclusive run 3: ./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 200000 --out  (load 0.51 0.91 0.97)
  time         0.327 s -> 6.12e+05 events/s
## tagged run 1: ./build/generate_tagged --channel 6Li-alpha --optics tagging --events 200000 --out  (load 0.51 0.91 0.97)
  time         0.599 s -> 3.336e+05 events/s
## tagged run 2: ./build/generate_tagged --channel 6Li-alpha --optics tagging --events 200000 --out  (load 0.51 0.91 0.97)
  time         0.605 s -> 3.305e+05 events/s
## tagged run 3: ./build/generate_tagged --channel 6Li-alpha --optics tagging --events 200000 --out  (load 0.55 0.91 0.97)
  time         0.645 s -> 3.101e+05 events/s
## coherent run 1: ./build/generate_coherent --config 1 --slope 50 --f0 0.04 --events 200000 --out  (load 0.55 0.91 0.97)
  time         0.187 s -> 1.071e+06 events/s
## coherent run 2: ./build/generate_coherent --config 1 --slope 50 --f0 0.04 --events 200000 --out  (load 0.55 0.91 0.97)
  time         0.210 s -> 9.509e+05 events/s
## coherent run 3: ./build/generate_coherent --config 1 --slope 50 --f0 0.04 --events 200000 --out  (load 0.55 0.91 0.97)
  time         0.185 s -> 1.082e+06 events/s
# load after: 0.55 0.91 0.97 2/674 409833
```

(The full per-run program output — the σ-per-category banner of each
run — is in the scratch log; it is identical run to run except the timing
line, checked by `diff` on the inclusive runs 1 and 2.)

Per run and median (events/s):

| channel | run 1 | run 2 | run 3 | **median** | 2026-10-30 (under load, no log) | ratio |
|---|---|---|---|---|---|---|
| inclusive ⁶Li config 1, A_zz | 5.36e5 | 5.40e5 | 6.12e5 | **0.540 M ev/s** | 0.65 | ×0.83 |
| tagged ⁶Li α, tagging optics | 3.336e5 | 3.305e5 | 3.101e5 | **0.331 M ev/s** | 0.36 | ×0.92 |
| coherent config 1 | 1.071e6 | 9.509e5 | 1.082e6 | **1.071 M ev/s** | 1.08 | ×0.99 |

A 200 k-event run lasts 0.19–0.65 s, so the run-to-run spread is about
±10 % (inclusive run 3 is 13 % above the other two); the recipe was kept
as on 2026-10-30 for comparability and the median is the number quoted.

## 2. Attribution of the inclusive change (A/B, out of tree) and a 2 M-event supplement

**Does this run's code cost the difference?** The three examples were built
OUT OF TREE from the pre-P2 commit `3d059d3` (P1 only: the factors moved,
no clip) and from HEAD `888ca3b`, with identical options (Release,
LHAPDF/HepMC3/PYTHIA on, Python off; `git archive` of each commit into the
scratch dir, the repository untouched), and run alternately, three rounds,
200 k events (`ab_clip.log`, verbatim):

```
# A/B clip attribution, 2026-11-07T04:21:44Z, 4 cores
# pre (3d059d3) built in 23 s; warnings: 0
# post (888ca3b) built in 24 s; warnings: 0
# load before runs: 1.67 1.05 1.01 1/668 410786
pre inclusive round 1: time 0.333 s -> 6.01e+05 events/s (load 1.62)
pre tagged round 1: time 0.568 s -> 3.521e+05 events/s (load 1.62)
pre coherent round 1: time 0.178 s -> 1.121e+06 events/s (load 1.62)
post inclusive round 1: time 0.349 s -> 5.73e+05 events/s (load 1.62)
post tagged round 1: time 0.594 s -> 3.369e+05 events/s (load 1.62)
post coherent round 1: time 0.187 s -> 1.07e+06 events/s (load 1.62)
pre inclusive round 2: time 0.334 s -> 5.98e+05 events/s (load 1.62)
pre tagged round 2: time 0.575 s -> 3.476e+05 events/s (load 1.73)
pre coherent round 2: time 0.189 s -> 1.06e+06 events/s (load 1.73)
post inclusive round 2: time 0.324 s -> 6.17e+05 events/s (load 1.73)
post tagged round 2: time 0.637 s -> 3.142e+05 events/s (load 1.73)
post coherent round 2: time 0.186 s -> 1.075e+06 events/s (load 1.73)
pre inclusive round 3: time 0.326 s -> 6.14e+05 events/s (load 1.73)
pre tagged round 3: time 0.568 s -> 3.524e+05 events/s (load 1.67)
pre coherent round 3: time 0.184 s -> 1.089e+06 events/s (load 1.67)
post inclusive round 3: time 0.339 s -> 5.9e+05 events/s (load 1.67)
post tagged round 3: time 0.578 s -> 3.459e+05 events/s (load 1.67)
post coherent round 3: time 0.200 s -> 1.001e+06 events/s (load 1.67)
# load after: 1.67 1.08 1.02 2/667 410920
```

| channel | pre `3d059d3` median | post `888ca3b` median | post / pre |
|---|---|---|---|
| inclusive | 6.01e5 | 5.90e5 | **×0.98** |
| tagged ⁶Li α | 3.521e5 | 3.369e5 | **×0.96** |
| coherent | 1.089e6 | 1.070e6 | **×0.98** |

So the run's code changes (P1 → HEAD: the A3 clip wrapper's F₁ per g₁ call,
the A17 table, the A10 attributes, …) cost 2–4 % of the examples' rate, not
the 17 % between today's inclusive median and the 2026-10-30 value. (The
A/B rounds ran with the 1-minute load average still at 1.6–1.7 from the
two `-j3` builds that preceded them; both arms saw the same conditions.)

**Is the 200 k recipe itself the limit?** With 2 M events per run (3–6 s a
run) on the idle machine, the in-tree `build/` binary (Python on) and the
out-of-tree HEAD binary (Python off), two rounds each (`throughput_2M.log`,
verbatim):

```
# P6 supplement: 2 M events per run, in-tree build/ (Python ON) vs out-of-tree HEAD build (Python OFF), idle machine, 2026-11-07T04:23:34Z, load 0.82 0.94 0.97
in-tree inclusive round 1: time 3.293 s -> 6.07e+05 events/s (load 0.83)
in-tree tagged round 1: time 6.057 s -> 3.302e+05 events/s (load 0.85)
in-tree coherent round 1: time 1.886 s -> 1.06e+06 events/s (load 0.86)
out-of-tree inclusive round 1: time 3.333 s -> 6e+05 events/s (load 0.86)
out-of-tree tagged round 1: time 6.078 s -> 3.291e+05 events/s (load 0.88)
out-of-tree coherent round 1: time 1.888 s -> 1.059e+06 events/s (load 0.88)
in-tree inclusive round 2: time 3.315 s -> 6.03e+05 events/s (load 0.97)
in-tree tagged round 2: time 6.069 s -> 3.295e+05 events/s (load 0.97)
in-tree coherent round 2: time 1.942 s -> 1.03e+06 events/s (load 0.98)
out-of-tree inclusive round 2: time 3.256 s -> 6.14e+05 events/s (load 0.98)
out-of-tree tagged round 2: time 6.050 s -> 3.306e+05 events/s (load 0.98)
out-of-tree coherent round 2: time 1.886 s -> 1.061e+06 events/s (load 0.98)
# load after: 0.98 0.97 0.98
```

| channel | in-tree, 2 M (2 runs) | out-of-tree, 2 M (2 runs) | the 200 k recipe's median |
|---|---|---|---|
| inclusive | 6.07e5, 6.03e5 | 6.00e5, 6.14e5 | 5.40e5 |
| tagged ⁶Li α | 3.302e5, 3.295e5 | 3.291e5, 3.306e5 | 3.305e5 |
| coherent | 1.06e6, 1.03e6 | 1.059e6, 1.061e6 | 1.071e6 |

The two builds agree within 2 %, so the build configuration is not a
factor. The tagged and coherent rates agree between the recipes; the
inclusive 200 k runs read 10 % low: a 0.33 s run is short enough to be
clock-ramp limited (run 3 of the recipe, 6.12e5, already matched the 2 M
value). The number A16's recipe gives is 0.540 M ev/s and is quoted as
such; the tighter idle-machine value for the inclusive example is
**0.61 M ev/s**, and the honest comparison with 2026-10-30 (0.65, under
load, no log) is then ×0.93, of which the code accounts for ×0.98; the rest
is not attributable, that day's run having no log.

## 3. Files

- `docs/USAGE.md` §8: one dated A16 clause appended to the 2026-10-30
  sentence on its own line (line-neutral: PHYSICS_CHANNELS and the D7 plan
  cite §8 and §9 by line ranges, and the HepMC3-bottleneck line by number).
- `docs/open_items/run_2026-11-03/STATUS.md`: the P5 row's hash back-filled,
  the P6 row DONE.
- this record.

## 4. Gates

No code changed; the gates that read documents were run after the edits:
`check_physics_channels_links.py` strict, `--records-only`,
`--audit-ranges`; `check_spdx_headers.py`; `python/tests/test_doc_link_gate.py`
(tallies in the P6 board row and the commit message).

## 5. Open points

1. The 2026-10-30 numbers were taken under load and without a log; they are
   kept in USAGE as history, labelled, and are not comparable at the 10 %
   level with today's. §2 attributes ×0.98 of the inclusive change to the
   run's code and leaves the rest unattributed.
2. 200 k events give a ±10 % run-to-run spread and read the inclusive rate
   10 % low on this machine (clock ramp); a future re-measurement should use
   2 M events per run (3–6 s) and say so. The recipe was kept here because
   A16 asked for the 2026-10-30 recipe on an idle machine; the 2 M values
   are recorded beside it.
3. The in-tree build (Python on) and an out-of-tree build (Python off) run
   the examples at the same rate (within 2 %); no action.
