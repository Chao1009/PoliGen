<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-09-23 — one page

**What this run is.** Four pieces of work since `b647c16`: the thorough
review of 2026-09-16 (`239ac88`), the reference search (`d735abc`), the
polligen port gate (phase A) and the M1 benchmark harnesses (phases B1–B3),
closed by a Documents stage. A, B1–B3 and the Documents stage are
**uncommitted** on base `d735abc`. Board: `STATUS.md`. **No physics default,
registry status label, constant or reference-JSON number moved, and no
decision was taken.** Nothing was pushed, published or sent.

## What it measured

* **The review** (`239ac88`, six lenses over `b647c16`): 78 findings, 28
  serious, 16 confirmed, 15 applied, 34 of 40 minors applied, the eight
  residuals closed. Bit-identical at the defaults on **207** CLI
  configurations (hashed against a from-scratch `b647c16` build; threaded =
  serial at 1/4/8 threads). **One behaviour change, on a non-default
  configuration:** a tilted fill on a tagged channel with a spin-dependent
  struck-cluster term is now refused (it computed the wrong spin orientation).
* **The reference search** (`d735abc`): `docs/references/REFERENCES.md`,
  **266** references, **56** PLEASE DOWNLOAD; the sibling corpus at **200**
  PDFs / **244** index entries. Item **59** (the BONuS supplemental
  material) was added to §1 on 2026-09-23, so **57** now.
* **The port gate** (`phase_A_port_gate.md`): `tagged.json` is `polligen`'s
  again in every block (polligen code as of the sibling's `c0f86a8` = its HEAD `4726812`, which differs only under `refs/`, where the corpus is committed); the two spin-1 model blocks
  moved by ≤ 1.33e−13 / 2.28e−13 relative (both in the 1e−9-gated
  quadrature; every 1e−12-gated entry ≤ 6.1e−15); the C++ gate asserts
  `provenance == "polligen"`.
* **The harnesses** (`validation/benchmarks/`, `BENCHMARK_PLAN.md` §4, each
  re-run 2026-09-23 by the Documents stage with the same output; none tuned):

| row | result | the number, with its window |
|---|---|---|
| 1 EPIOS source modes | **PASS** | round-trip residual 0 over 8 modes, tol 1e-12; the EST fill matches the ideal P_zz at 3 of 8 modes |
| 3 HERMES b₁ᵈ | **FAIL, recorded** | χ²(b₁ᵈ)/6 over 6 bins: 22.26 (A = 2 convolution + MSTW) vs 21.80 for b₁ = 0; 5.26 Miller ×0.5, 5.30 Miller ×1 |
| 4 NMC F₂(⁶Li)/F₂(D) | **PASS** | 4.6269/4 on x ≥ 0.30 (p = 0.328); 26.941/15 on the 15 on-grid points (p = 0.029); tol p ≥ 0.01 |
| 5 UVa FB ⁶Li C0 | **FAIL, recorded** | shipped q₀ = 3.0998 fm⁻¹, +0.2713 above the band [2.6944, 2.8284] fm⁻¹ |
| 6 RFY magnetization | **BLOCKED** | on a document; tree side r_mag(⁶Li) = 2.9469 fm, uncompared |
| 7 DJANGOH Rad/noRad | **BLOCKED — cannot be a benchmark as planned** | all 8 logs: IEL2 = IEL31..33 = 0 (elastic tail off); `rc_tail` is the elastic tail; no shared O(α) term |
| 9 EST identity | **PASS** | 1.138e-14 over 19 999 interior points, tol 2.31e-14; the planned 1e-14 does not hold |
| 10 BONuS spectator | **BLOCKED** | on a document; the CLAS database holds BONuS only as F₂ⁿ/F₂ᵖ |

* **Findings beyond the plan.** POLRAD's deuteron row 2 was stale:
  σ_u = −1.0303e−03 (not −1.001e−03), σ_q/σ_u = +0.0623 (not +0.0642),
  rows 1 and 3 reproducing to every digit — now corrected at the sixteen sites a tree-wide grep finds (`phase_B3_chain_rc.md` §6.3, closing note; two of them only by the fix stage), the
  C0-shape suppression factors becoming 67 / ≈ 255; T8(c′) passes at the same
  tolerance (su/POLRAD 0.9548, measured). The UVa FB C0 shape lies
  **outside** the [ho, vmc-ft] RC band at every priced point (x = 0.01–0.30,
  Q² = 5), moving ΔA_zz by ≤ 2.2e−07. The withdrawn "0.31 GeV²" dip has a
  source (Li *et al.* 1971, q = 2.828 fm⁻¹). COMPASS's "11 % tensor
  polarization" is T = ½(1 − 3p₀), half this tree's P_zz. George–Knutson's η
  is a restricted phase-shift analysis, now called a consistency band on the
  quadrupole dial at every live site a grep finds except the unsent letter, and eSTARlight an exclusive-VM rate scale, not a coherent baseline — the last of each (test comments, `PHYSICS_CHANNELS.md` §13, `OPEN_ITEMS_SOLUTIONS.md` §11.2/§11.3b, `configs.py`, six eSTARlight doc sites, one binding docstring) only by the fix stage; dated records carry an adjacent note (George–Knutson) or are left as written (eSTARlight: two), and one TEST_CASE name is kept — the full list is in the §B25 box.

Phase records: `phase_A_port_gate.md`, `phase_B1_spin_deuteron.md`,
`phase_B2_nuclear.md`, `phase_B3_chain_rc.md`.

## What this run did NOT establish

Carried forward from `../run_2026-09-06/SUMMARY.md`, each with what changed.

1. **Decay-lepton reconstruction efficiency — unbounded.** Unchanged.
2. **Whether O5's band top crosses 3 σ** (2.842 … 3.292 σ). Unchanged.
3. **b₁ beyond A = 2 is not validated.** Unchanged, and sharpened: against
   HERMES the A = 2 gate's own configuration is indistinguishable from b₁ = 0
   (22.26 vs 21.80 over 6 bins) — the gate tests the kernel against CDKS, not
   against data — and HERMES cannot decide registry row 2 (Δχ² = +0.035).
4. **b₁(⁷Li) is not implemented.** Unchanged.
5. **No external exact elastic tail.** Unchanged, and the planned second
   engine is gone as planned: DJANGOH's published table has the elastic tail
   off. New: the RC band does not bracket the one data-derived C0 shape.
6. **O3 from GFMC configurations.** Unchanged.
7. **Nothing pushed, published or sent.** Unchanged; the letter carries two
   over-claims the tree's live sites no longer make (row 21, §B25 box).
8. **The sibling was not modified.** Its `polligen` was imported read-only
   (phase A's dump and comparisons; the Documents stage's scratch re-dump).
   Found there and left for its owner: `refs/refs_dict.json` misquotes
   Koivuniemi's ±51.3 % / ±92.4 % as "51 ± 3" / "92 ± 4".
9. **27 decisions pending; this run took none.** Row 2 is now priced against
   data, rows 26 and 27 carry dated facts, row 21 a note; the C0-zero price
   sits under §B15 (no row names `LI6_FF_HO_*`).
10. **Gate residues** (anchor pins checked for presence, not uniqueness; +1
    range drifts). Unchanged; the SPDX gate now also covers
    `validation/benchmarks/*.py`.

Opened by this run: no measurement of the tensor polarization of ⁶Li was
found by any search (an absence, stated with its searches, not a proof); the
UVa FB row's provenance and George–Knutson's scattering system are
unverified (neither source document held); no `REPORT.md` generator exists
and no test reads `LIPOLGEN_BENCH` (no row needs it yet); phase B3 saw one
unreproduced doctest failure (1 of 17 242 712 assertions, not captured)
while three agents shared the tree.

## For the user

**Download** (links as in `../../references/REFERENCES.md` §1; drop the files
into `validation/benchmarks/data/`):
* Tkachenko *et al.*, PRC **89** (2014) 045206, **Supplemental Material** —
  <https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.89.045206>
  (§1 #59; needs an APS login). Unblocks row 10.
* Rand, Frosch, Yearian, Phys. Rev. **144** (1966) 859 —
  <https://doi.org/10.1103/PhysRev.144.859> (§1 #8), **or** de Jager, de
  Vries, de Vries, ADNDT **14** (1974) 479 —
  <https://doi.org/10.1016/S0092-640X(74)80002-1> (§1 #9). Unblocks row 6.

**Decide:** the batch in `../run_2026-09-03/AUTHOR_DECISIONS.md` — 27 rows,
none taken. New since the last run: row 2's HERMES price (§B2(c) box), rows
26 and 27's facts (§B27, §B28 boxes), row 21's two draft corrections (§B25
box), and the C0-zero price (§B15 box).

**Run or request:** a DJANGOH production with `IEL31..IEL33 ≠ 0` (elastic
radiative tail on) at the cuts of `eic/InclusiveDjangohSamples` — the only
thing that turns row 7 into a benchmark. It needs the DJANGOH source (not
public) or an ePIC production request.

## Tallies

Measured at the end of the fix stage, 2026-09-23, on the working tree (base
`d735abc` + phases A, B1–B3, the Documents stage and its fix stage). Baseline
= base commit `d735abc` where it was measured there (phase A §4a, phase B1 §6);
the strict gate and `--audit-ranges` were first measured at phase A's end.

| gate | measured | baseline |
|---|---|---|
| `cmake --build build -j` | exit 0 | — |
| `build/lipolgen_tests` | **416 cases / 17 242 712 assertions / 0 failed / 0 skipped** | 414 / 17 242 629 at `d735abc` (+2 provenance assertions, phase A; +2 cases / +81, `tests/test_bench_spin.cpp`); comment-only edits since change no count |
| `python -m pytest python/tests -q` | **1134 passed / 153 skipped / 0 failed** | 1107 / 151 at `d735abc` (+27 / +2: the eight harnesses' `test_bench_*.py` and phase A's gate; `test_doc_link_gate.py`'s as-of anchors 108 → 110 moved with the Documents stage's two annotations) |
| strict link gate | **1240 / 95 / 7 external / 0 broken / 24 allow-listed** | the same at phase A's end |
| `--fix` / `--record-ranges` | Documents stage: 0 fixed, 2 fingerprints re-recorded (`docs/USAGE.md` §2a and §9, read); fix stage: not needed (every cited file's line count kept) | — |
| `--audit-ranges` | **209 range citations, 179 blocks, 0 REFUSED** | 0 REFUSED at phase A's end |
| `--records` | **799 citations in 49 dated run records, 0 broken** | 747 / 43 at `d735abc` |
| `check_spdx_headers.py` | **112/112** (the gate now covers `validation/benchmarks/*.py`: +8; +1 `tests/test_bench_spin.cpp`) | 103/103 at `d735abc` |
| the eight harnesses | each re-run three times: exit 0, verdicts unchanged (3 PASS, 2 FAIL recorded, 3 BLOCKED); a harness copied away from its `data/` exits 2 | — |
| `dump_polligen_reference.py` re-dump (scratch copy, the tree's `env.sh`) | all nine `validation/reference/` files and `validation/README.md` byte-identical to the tree; a patched pre-fix sign is refused with **no** file written | — |

(Until the fix stage this table read 798 / 49 and 764 / 46 and SPDX 104 in the
baseline column — the shared tree mid-build, not `d735abc`.)

## Verification (2026-09-23)

Four adversarial lenses verified the stage's output; the measured numbers,
the gates and the defaults held under all four. Their findings were
wording and bookkeeping, every one fixed by the fix stage the same day:
POLRAD row 2's stale +0.064 at `PHYSICS_CHANNELS.md`'s polarised-QE row and
in `run_2026-09-03/STATUS.md` (the "every site" sentences now list the
sixteen sites, `phase_B3_chain_rc.md` §6.3); `06_critic.md` §7 item 4 restated
as a finding of absence with its searches, the EPIOS misattribution removed;
George–Knutson "measured" and eSTARlight "coherent baseline" residues
(tests' comments, §13 of `PHYSICS_CHANNELS.md`, `OPEN_ITEMS_SOLUTIONS.md`,
`configs.py`, six docs, one binding docstring; dated notes in four records;
the §B25 box and this page now list what was and was not changed); these
tallies; the `LIPOLGEN_BENCH` opt-in named as a convention no test reads;
the sibling's state (HEAD `4726812`, `refs/` only; the corpus is committed);
two refs_dict entries that do exist (`file: null`); the CLAS-DB cell of
`06_critic.md`; the F-6 deviation; four IEL flags, not three; one uniform
harness exit convention (0 = ran, whatever the verdict; 2 = broken); runtime
ranges instead of single figures; a stale binding comment; the date of
`tagged.json`'s provenance keys (`a7b3d18`, 2026-09-06); and the dump's sign
guard, now a pre-flight that runs before any file is written.
