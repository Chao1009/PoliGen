<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Decisions of 2026-11-10

Taken on the nine items run 2026-11-03 left for the author
(`../run_2026-11-03/SUMMARY.md` "For the author"; the same list on
`../run_2026-11-03/STATUS.md`). The author's words, in full:

> 1. go with your recommendation. 2. re-dump, for the others go with your
> recommendation

So: **J1 and J2 are the author's decisions** — on item 1 the author took the
recommendation as it stood in the summary, on item 2 the author's own word is
"re-dump". **J3–J9 are the SUPERVISOR's choices under the author's
delegation** ("for the others go with your recommendation") — to be labelled
so wherever they are quoted, never as the author's words. The numbering
follows the summary's list (J1 is its item 1, …). Each entry gives the choice
(bold), why, what changes in the tree and what does not. The phases that
implement them, their records and gate tallies are on the board, `STATUS.md`.

Ground rules carried over from the previous runs: a default moves only under
an entry below; every number that depends on a moved default is re-taken with
a dated comment naming the entry (`J<n> (2026-11-10)`), or the test passes
the old value explicitly where its purpose is to pin the old physics;
tolerances are declared before verdicts (where one is chosen after the
values were seen, the record says so); a benchmark that disagrees is
recorded, never tuned; `validation/reference/*.json` change only as J2's
re-dump writes them, never by hand; `docs/HEPMC3_CONVENTION.md` is edited by
J6 only, with the T6 row's sha256 re-pinned; the hashed design blocks of
`t4_djangoh_rad_noRad.py` and `t2_boeglin_nk.py` are never edited; the four
author decisions of 2026-10-27 (D4–D7) and the seventeen entries of
2026-11-03 (A1–A17) are not re-opened, except A4's sentence on the rank-2
refusals, which item 1's recommendation names for re-opening; the sibling
repository `PolarizedLithiumSim` is read-only; nothing is published or sent.

**J1 — the bare deuteron gets its own b₁ slot: on an inclusive deuteron run
the kernel's rank-2 slot is the deuteron camp itself (`MillerB1`, per
nucleon, unmapped) instead of `Li6B1(MillerB1)`; then the toy b₁'s high-x
behaviour is re-opened on what that exposes.** *The author's decision, on
the recommendation of `../run_2026-11-03/SUMMARY.md` item 1, which read:
"give the deuteron its own b₁ slot (MillerB1 unmapped) and then re-open the
toy b₁'s high-x behaviour, which a ~3× larger deuteron rank-2 term at the
top cell will expose." Why: `default_inclusive_kernel` keys its rank-2 branch
on spin 1, so a deuteron run reads b₁(⁶Li)/nucleon = (2/6) ×
`LI6_B1_RANK2_TRANSFER` × b₁ᵈ = 0.32753 × MillerB1 — "a third of the
nucleons sit in a deuteron cluster whose alignment the α–d motion dilutes",
which is a statement about ⁶Li and false of a deuteron; HERMES χ²(b₁)/6 is
14.56 on it against 5.26 on MillerB1 (`../run_2026-11-03/phase_P3_report_rules.md`
§5 item 1); the tree itself files it as "documented behaviour … a deuteron
run is not a ⁶Li run … a phase-F decision" (`docs/OPEN_ITEMS_SOLUTIONS.md`
§14, "Adjacent defects"). What changes: the spin-1 branch of
`default_inclusive_kernel` for the deuteron; every inclusive deuteron
rank-2 number (b₁, b₂, A_zz, the tensor weight) by the factor 1/0.32753,
to be measured; the deuteron run's `meta` and banner (the ⁶Li transfer is
no longer read there); the HERMES row's wording (its headline is then the
b₁ a default deuteron run carries, so convention 8's default-run clause
applies again and the 14.56 becomes history); the pins and documents that
state the old reading, each re-taken with `J1 (2026-11-10)`. What does not:
⁶Li, ⁷Li, the tagged and coherent channels (to be measured
column-identical); `validate()`'s refusal of `--b1-model cdks` /
`li6-convolution` off ⁶Li; the Δ slot; any reference JSON (polligen's
kernel takes its b₁ as an argument, so no port gate reads this default; a
test whose purpose is the old mapping passes `Li6B1` explicitly). The
re-opening: the phase measures which default deuteron runs build, and where
the toy — Miller's digitized curve up to x = 0.9 and this code's own taper
above it — leaves the exact spin-1 positivity bound; the treatment is then
decided by a dated note under this entry, on those measurements.*

**J2 — the polligen reference tables are re-dumped in place on this machine
and committed exactly as one run of `validation/dump_polligen_reference.py`
writes them.** *The author's decision ("re-dump"). The facts it rests on
(`../run_2026-11-03/phase_P5_port_gate.md`; the supervisor's note under A13):
run on 2026-11-06 the script rewrites 3170 numbers of six reference JSONs at
the rounding level (≤ 2 ULP in four files, ≤ 1.9e−16 absolute in `spin.json`,
up to 9.4e−13 relative among `xsec.json`'s values above 1e−6) and one stale
manifest string, with the sibling unchanged at 4726812; the likely, unproven
cause is this machine's libc/libpython upgrade of 2026-10-26. What changes:
those numbers and that string; `tagged.json` becomes one dump's output (the
`li6_alpha_vmc` block A13 spliced is the script's own); `validation/README.md`
as the dump writes it. What does not: any tolerance — every port gate must
pass at its present rtol against the re-dumped files, and a gate that does
not is recorded and brought back, never widened in the phase; LiPolGen's own
numbers (no source changes); the sibling. Also recorded, so the next drift
can be traced: the interpreter, numpy and libc versions of the dump, and
that a second dump run reproduces the first byte for byte.*

**J3 — BONuS: the declared CHECK thresholds and the shipped 161 × 161 weight
table stay; that the two misses on the AV18 headline do not move a verdict
becomes a pinned test.** *The supervisor's choice. Why: the CHECKs (table
convergence |dv| < 5e−4; bin averaging < 2 %) were declared on the Hulthén
weight and hold there; on the AV18 headline they read 5.5e−4 / 8.3e−4 and
−3.2 % (`validation/benchmarks/t2_bonus_spectator_shape.py`, A5 NOTE). Moving
the thresholds after the values would be tuning; refining the shipped table
fourfold buys 3–8 × 10⁻⁴ on an opt-in weight whose own ceiling is an
unpriced systematic (A12). What a reader needs is that no verdict depends on
either deviation — checked once by hand on 2026-11-06 (0 of 28 blocks) and
true only as long as nobody moves the weight. What changes: a bench-gated
test that recomputes the headline on the 321 × 321 table and with
bin-averaged weights and asserts that no block's verdict flips. What does
not: the thresholds, the table, the headline, the row's status.*

**J4 — HERMES b₁: the A_zz comparison gains a recorded sub-row on the bins
inside the shipped grid; the headline keeps all six bins and its status.**
*The supervisor's choice. Why: the lowest HERMES bins (⟨Q²⟩ = 0.51, 1.06 and
1.65 GeV²) lie below the shipped CT18NLO grid's Q² floor, where D7 freezes
F₂ at the edge (the author's decision, not re-opened); there the b₁ ↔ A_zz
conversion uses an F₁ᵈ that is not the measured one (12.40 against the 7.04
HERMES's own table implies at x = 0.012). Dropping those bins from the
verdict after seeing them would be tuning, so they stay; the reader is given
the in-grid χ² beside the headline, labelled as chosen after the values were
seen, and the bins on the frozen edge are named. What changes: one sub-row
and its sentence in the harness, the README row, the report. What does not:
the headline, its rule, its status, the freeze.*

**J5 — `g1_clipped_rate_fraction` stays a cell-centre estimate; its
definition says so and carries a measured bracket.** *The supervisor's
choice. Why: the share is 2–5 × 10⁻⁷ of the rate; integrating the bound
inside each cell in the library would cost every run for the next digit of a
diagnostic. What was missing is the statement and one number for how far the
estimate can be off: the share counted when the bound binds at ANY corner or
the centre of a cell (an upper bracket) and at ALL of them (a lower one),
measured once on the default runs and written beside the definition. What
changes: the key's documentation (header, USAGE) and the phase record. What
does not: the key, its value, the banner number, any code path of a run.*

**J6 — a HepMC3 file states its provenance or says that it has none, and the
two `*_sf_extrapolation` keys join the one list.** *The supervisor's choice.
Why: A10 made the HepMC3 file carry what the npz `meta` carries, from one
list (`run_provenance`), but left two gaps (`../run_2026-11-03/phase_P4_infrastructure.md`
§5 item 2): a file written without a `Pipeline` has no `lipolgen.*`
attribute at all, indistinguishable from a file written before A10; and the
two keys that say how each structure-function backend is evaluated outside
its grid — the default run has 39 % of its rate on that frozen edge — are
added by the Python module alone, so C++-written files lack them. What
changes: `HepMC3Writer` writes one attribute saying that no run provenance
was set when `set_run_attributes` was never called; the classification
moves to where `run_provenance` can make it (the backends answer for
themselves, the core still not linking LHAPDF), so both paths write the same
keys from one function; `docs/HEPMC3_CONVENTION.md`'s provenance section
says both, its sha256 re-pinned in the T6 harness with the old digest in a
dated comment. What does not: any npz `meta` value (to be measured
byte-identical on the P4 run list); the attributes of a CLI-written file
(T6 clause C24); event content.*

**J7 — a verdict set that is installed but is not the vendored version, or
lacks its member file, is treated as a missing set: the harness raises,
exits 2, and the report row is `blocked (environment)`.** *The supervisor's
choice. Why: A8's reason applies word for word — "a harness whose verdict
cannot be computed is broken in that environment" — and the set this row
was measured with is then not installed; `t5_delta_moment_bound` today
exits 2 for the one and prints a `BLOCKED` row with exit 0 for the other
(`../run_2026-11-03/phase_P4_infrastructure.md` §5 item 3), two outcomes for
one situation. What changes: that harness's two `check_set` branches, the
pytest that pins them, `validation/benchmarks/README.md` convention 6's
sentence. What does not: the row with the right sets installed; the
tier-off and module-missing `BLOCKED` rows of T4, T6 and this harness
(convention 6 as written).*

**J8 — the leading-order evaluation of the NLO sets stays the filed
development item and is not scheduled here; its entry is brought up to date
with what A3 changed. EPPS21 stays the option A2 recorded.** *The
supervisor's choice. Why: the symptom the item was filed for (|A₁| > 1 from
x ≳ 0.83) is bounded since A3, with the clip on 2–5 × 10⁻⁷ of the rate; NLO
coefficient functions are a design of their own (F_L, heavy flavours, the
scheme of a mixed unpolarised/polarised pair) with their own validation,
and by A3's own reasoning would not by themselves restore the positivity of
the mixed pair; nothing ranked 1–7 in `docs/benchmarking/08_channel_priority.md`
waits on it. Records only.*

**J9 — USAGE §8's throughput recipe becomes 2 M events per run (idle
machine, a fresh process each, median of 3, the raw log kept), and the
numbers are re-measured with it at the end of this run.** *The supervisor's
choice. Why: `../run_2026-11-03/phase_P6_throughput.md` §2 measured the
200 k recipe reading the inclusive example 10 % low (0.3 s runs are
clock-ramp limited) while the two recipes agree on the other two examples.
What changes: the recipe sentence and the headline numbers of USAGE §8, on
the same line; the 200 k numbers stay as history.*
