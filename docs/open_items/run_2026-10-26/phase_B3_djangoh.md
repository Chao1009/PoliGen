<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — phase B3: row 7, DJANGOH's elastic radiative tail isolated

`BENCHMARK_PLAN.md` §4 row 7; `PLAN.md` phase B3 (the design as reviewed);
`../../benchmarking/08_channel_priority.md` §4, ingredient **F** (radiative
corrections). Base `87787dd`. Owned and written by this phase:
`validation/benchmarks/t4_djangoh_rad_noRad.py` (re-wired),
`validation/benchmarks/data/djangoh_eltail_ep18x275.json` (new),
`python/tests/test_bench_djangoh_eltail.py` (new), this file, and
`docs/benchmarking/figures/t4_djangoh_rad_noRad.png` (written once by the
harness's `--plot`). **No physics default, constant, registry label,
tolerance of another row or reference JSON moved; nothing was committed.**
Every DJANGOH file, the PDF set and every log stay outside the repository in
the session scratch directory (`$SCR` below); only numbers and input cards are
vendored.

## 0. Verdict

| item | verdict |
|---|---|
| row 7 | **BLOCKED → PASS.** DJANGOH 4.6.21 was built from the public source and run with its elastic radiative tail isolated (IEL31 + IEL32 + IEL33 on, everything else off, WMIN ≤ M_p) in a window matched to LiPolGen's W_l² ≥ 9 GeV². Against it, LiPolGen's proton tail per Q² bin is R = **1.0199 / 1.0261 / 1.0219 / 1.0098** (`polrad-full`) and **0.6118 / 0.8968 / 1.0112 / 1.0067** (`t-peak`, the shipped default). The tolerance, fixed in the harness header before any run: `polrad-full` \|R − 1\| ≤ 0.03 + 2δ_b, `t-peak` (a lower bound) R ≤ 1.03 + 2δ_b, δ_b = DJANGOH's integration error 3.7 / 4.8 / 4.4 / 3.9 × 10⁻⁴. Both hold in all four bins. |
| control | **8 / 8 reproduced.** The eight ePIC cards, run on this build, give the ePIC README cross sections to −2.6 × 10⁻⁴ … +4.0 × 10⁻⁴ (\|pull\| ≤ 1.74 against the declared 3σ rule). So this build reproduces the eight ePIC integrated cross sections (Born and inelastic channels) to ≤ 4 × 10⁻⁴; the elastic-tail channels IEL31–33, which carry the verdict, have no external control (every ePIC card runs them off). |
| recorded, not in the verdict | Straight from the quadrature (no table) LiPolGen's exact tail is **1.0004 / 1.0010 / 0.9995 / 0.9866** of DJANGOH's. The shipped (x, y) tail table adds **+1.9 … +2.5 %**, all of it from interpolating the tail column. DJANGOH's running α raises its tail by ×1.013 … 1.021. |

## 1. What was compared, and the design fixed before the runs

**The question row 7 now asks.** Does LiPolGen's `rc_tail` (σ_tail/σ_Born,
applied as a weight) carry the right **elastic radiative tail** e p → e p γ
under a DIS bin? The comparison is on a proton at the DJANGOH beams, e 18 × p
275 GeV (LiPolGen config 2), in four Q² bins [1, 10], [10, 100], [100, 10³],
[10³, 10⁴] GeV². The compared number is the **absolute** tail per bin, in µb per
e-p collision. It needs no Born. The ratio R = σ_tail(LiPolGen)/σ_tail(DJANGOH)
carries the verdict.

The displayed "fraction" divides both tails by one Born: DJANGOH's Rad=0 at
the same cuts. So f_LiPolGen/f_DJANGOH = R exactly. LiPolGen's Born uses the
shipped Toy structure functions; DJANGOH's uses CTEQ6.1M. Their ratio is
1.469 / 1.116 / 1.155 / 1.065 in the window. That is a structure-function
difference, not a radiative correction, so it is recorded and never compared.

**The design and the tolerance were written into the harness header first**
(the `==== DESIGN BLOCK ====` of `t4_djangoh_rad_noRad.py`; PLAN.md ground
rule 4). They were fixed before any DJANGOH run, the control included, and
before any LiPolGen number in the new window was computed. The block's sha256,
`fbe06789…c91b`, is recorded in the vendored file. The harness re-checks it on
every run and exits 2 if the block was edited. The tests pin it, and refuse
one edited character.

| UTC 2026-10-26 | event |
|---|---|
| 18:57:41 | `cteq61` fetched (curl) |
| 18:59:43 | DJANGOH built |
| **19:11:53** | **design block written into the harness header** |
| 19:13:21 | first DJANGOH run (a control card) |
| 19:13:35 – 19:14:22 | the eight control runs |
| 19:14:44 – 19:16:37 | the 32 isolating runs |
| 19:18, 19:30 | vendored. The second vendoring re-parsed the same logs after a parser fix (§4); no run was repeated and no number moved |

The design, in short (the header has the full text):

- **(a) IEL2 = 0.** In the source, IEL2 (`INT2(3)`) integrates `HSEL22`:
  elastic e p → e p **at x_l = 1**, with its virtual and soft-photon
  corrections below the infrared cut `XIRDEL` = 2 × 10⁻² (manual §15.13).
  That is the non-radiative elastic peak at W_l = M_p, outside any
  W_l² ≥ 9 window, so it is not a tail.
- **(b) Isolation by option (i).** `INT-OPT-NC 0 0 0 0 0 0 20 20 20` turns
  INC2 and INC31..34 off and IEL31..33 on. The MAIN program integrates
  exactly the flagged channels, so no on/off pair is needed.
  - IEL31, IEL32 and IEL33 are three phase-space mappings (`HSELK1`,
    `HSELK2`, `HSELCO`) of one squared matrix element: lepton-line
    bremsstrahlung with the elastic F₁, F₂ and the lepton mass. Only their
    sum is compared.
  - The Born is the ePIC Rad=0 card with only KINEM-CUTS changed.
- **(c) The window.** WMIN = 0.0 GeV on input, which `HSPRLG` raises to
  M_p = 0.93828 GeV. The elastic channels are zeroed only when WMIN > M_p.
  - **The manual gives no leptonic-W cut.** KINEM-CUTS cuts x_l, y_l, Q²_l
    and the hadronic W_h, and W_h = M_p on every elastic-tail event. So
    integration alone cannot apply W_l² ≥ 9.
  - The window is matched with the cuts the manual does provide:
    x_l ≤ XMAX_b = Q²_lo/(Q²_lo + 9 − M_p²), rounded down to 4 decimals and
    capped at 0.95. That gives **0.1096 / 0.5518 / 0.9249 / 0.95**, with
    1e−5 ≤ x_l and 1e−4 ≤ y_l ≤ 0.985. Every point of this rectangle has
    W_l² ≥ 9 GeV², and LiPolGen integrates the identical rectangle.
  - It differs from the W_l² ≥ 9 window by a high-x, low-y strip, measured in §5.
  - **`Q2MAX` is not applied** by `HSELK1`, `HSELK2` or `HSELCO`. They read
    Q2MIN, XMIN/XMAX and YMIN/YMAX only. So every per-bin number is the
    difference of two runs at Q2MIN = Q²_lo and Q²_hi, both with Q2MAX = 10⁸.
- **(d) Fixed in advance.**
  - PDF: `cteq61`, LHAPDF6 set 10150, the ePIC set.
  - Seeds: `RNDM-SEEDS 0 0`, RANLUX with the fixed seed 1213141567.
  - `INT-POINTS 100000`, against ePIC's 30000.
  - Control: ePIC's 18 iterations and `INT-ONLY 1`. Isolating runs: 20
    iterations per IEL channel and `INT-ONLY −1`. `START 0` everywhere.
  - A bin whose tail error exceeds 1 % would be rerun at 400000 points.
    None did.
  - The tail runs use **LPARIN(3) = 2**, i.e. no running α in the radiative
    cross section. LiPolGen's tail has fixed α and no vacuum polarisation.
    The ePIC value, 3, was run as a recorded variant.
- **(e) The tolerance**, from the header's term-by-term accounting.
  - **DJANGOH** computes the exact O(α) lepton-line tail: no peaking
    approximation, lepton mass kept, dipole form factors. Sources: manual §7
    and its Ref. [38], Bucoveanu & Spiesberger, EPJ A 55 (2019) 57; the
    matrix element in `HSELK1`, `HSELK2` and `HSELCO`.
  - **LiPolGen `polrad-full`** computes the same term: POLRAD Eq. (18) through
    the Eq. (44) proton slot, with the lepton mass and the same dipole.
  - So the leading-log vs exact 5–10 % does not apply to this pair. It is the
    accuracy of `t-peak+ll`, which this row does not wire.
  - The allowance is numerical: 2 % for LiPolGen's table and clip (chosen,
    not measured), plus < 2 × 10⁻⁴ of constants, rounded to 3 %. Twice
    DJANGOH's own integration error is added on top.
  - **`t-peak`** is POLRAD Eq. (38), the t-peak alone. The tree calls it a
    lower bound (`include/lipolgen/rc.hpp`, the `RcTailModel::TPeak`
    comment), so it is tested one-sided.

## 2. The build

Source: `github.com/spiesber/DJANGOH` at `e356b84`, cloned into
`$SCR/src`. Here `SCR=/tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/djangoh`.

- **The file built is `djangoh_h-4.6.21.f`**, copied to `djangoh_h.f`. Its
  banner, "VERSION 4.6.21, 05.12.2022", is the one every ePIC log prints.
- The repository's own `djangoh_h.f` is **4.6.20**. It differs only by the
  version strings, the HSFLAB/HSACPT zero-transverse-momentum guards of
  manual §15.1, and comments.
- Licence: none stated. The GitHub API reports `license: null`, and there is
  no LICENSE file. The program is not vendored.

```
mkdir -p $SCR/build && cd $SCR/build
cp $SCR/src/djangoh_h-4.6.21.f djangoh_h.f
cp $SCR/src/{djangoh_l,djangoh_u,djangoh_p,djangoh_t,sophia-dj,pythia-6.4.24-dj,jetset7409,gmc_random,cernlibs}.f .
cp $SCR/src/Makefile-sample Makefile
sed 's#^pwd    = .*#LHAPDF = /home/cpeng/Projects/polli/deps/install\npwd    = ${PWD}#' $SCR/src/config.mk > config.mk
make -j2      # 18:58:12Z: error, PYINIT REAL(8)/REAL(4) argument mismatch
# CFFLAGS += -fallow-argument-mismatch   -> 18:58:24Z: same error under -pedantic
# -pedantic removed                       -> 18:59:00-18:59:43Z: built
```

- **Build:** 42.4 s wall, 51.6 s user, 778 warnings; binary sha256 `5aea6fa8…98c72`.
- **Flags:** `-O -Wall -fno-automatic -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero -fallow-argument-mismatch`.
- **Toolchain:** GNU Fortran 11.4.0, LHAPDF 6.5.5 from the deps prefix.
- **Why the flags changed.** gfortran 11 refuses two calls to `PYINIT`:
  `LYSTFU` passes a default-REAL `ROOTS` and `HSWPDF` a DOUBLE one.
  `-pedantic` keeps that an error even with `-fallow-argument-mismatch`.
  **No source line was changed.**

**The PDF set.** The ePIC samples used CTEQ6.1M (`STRUCTFUNC 0 2 10150`),
which the prefix does not carry.

- `python3 $(which lhapdf) --pdfdir $SCR/lhapdf_sets install cteq61` failed.
  The prefix's `lhapdf` CLI needs LHAPDF's Python module, which
  `build_deps.sh` does not build, and its shebang calls `python`, which this
  machine lacks.
- Fetched instead with curl into a fresh directory:
  `curl -sS -L -o cteq61.tar.gz https://lhapdfsets.web.cern.ch/current/cteq61.tar.gz`.
  HTTP 200, 4 032 509 B, sha256 `5a6f1393…0b92`.
- Unpacked into `$SCR/lhapdf_sets`. SetIndex is 10150; member 0 has sha256
  `730c2ce0…c5ac`.
- Runs use `LHAPDF_DATA_PATH=$SCR/lhapdf_sets:<prefix>/share/LHAPDF`, so the
  shared prefix was not modified. The log confirms the load: "LHAPDF 6.5.5
  loading …/cteq61_0000.dat, LHAPDF ID = 10150".
- The ePIC logs load LHAPDF **5.9.1**'s `cteq61.LHgrid`. The difference is in
  interpolation, and is recorded in §3.

## 3. The control — the eight ePIC cards

The cards are `runcards/ep.Rad={0,1}.NC.Q2_*.in` from
`eic/InclusiveDjangohSamples` at `9869d9a`. Each card's git blob sha1 was
checked against the pinned tree. They ran verbatim except for three changes:
`RNDM-SEEDS 0 0` (ePIC: −1 −1), `INT-POINTS 100000` (30000), and `START 0`
(500000; integration only).

The ePIC error is the log's "Estimated relative error". The Rad=1 bin-1 log
is truncated, so the largest of the other three, 0.47 × 10⁻³, stands in, as
the design fixes.

| bin [GeV²] | Rad | this build [nb] | ePIC README [nb] | ratio − 1 | combined rel. error | pull | runtime |
|---|---|---|---|---|---|---|---|
| 1–10 | 0 | 628.04 | 627.8853416 | +2.46 × 10⁻⁴ | 1.41 × 10⁻⁴ | +1.74 | 0.5 s |
| 1–10 | 1 | 673.65 ± 0.154 | 673.3778912 | +4.04 × 10⁻⁴ | 5.22 × 10⁻⁴ | +0.77 | 24.9 s |
| 10–100 | 0 | 67.986 | 67.97276736 | +1.95 × 10⁻⁴ | 1.41 × 10⁻⁴ | +1.38 | 0.2 s |
| 10–100 | 1 | 74.031 ± 0.0145 | 74.05000239 | −2.57 × 10⁻⁴ | 4.72 × 10⁻⁴ | −0.54 | 23.0 s |
| 100–10³ | 0 | 3.2401 | 3.239636940 | +1.43 × 10⁻⁴ | 1.41 × 10⁻⁴ | +1.01 | 0.2 s |
| 100–10³ | 1 | 3.7603 ± 7.6 × 10⁻⁴ | 3.760005602 | +7.8 × 10⁻⁵ | 5.12 × 10⁻⁴ | +0.15 | 21.2 s |
| 10³–10⁴ | 0 | 0.070161 | 0.07015936779 | +2.3 × 10⁻⁵ | 1.41 × 10⁻⁴ | +0.16 | 0.2 s |
| 10³–10⁴ | 1 | 0.086138 ± 1.9 × 10⁻⁵ | 0.08615325533 | −1.77 × 10⁻⁴ | 4.99 × 10⁻⁴ | −0.36 | 21.1 s |

**All eight are reproduced by the declared rule**, \|pull\| ≤ 3. The four
Born numbers come from deterministic quadrature and all sit high, by
+0.2 … +2.5 × 10⁻⁴. That is the size of the LHAPDF 5.9.1 → 6.5.5
interpolation change, and it is not resolved further. Our VEGAS errors are
1.95–2.28 × 10⁻⁴, against ePIC's 4.3–4.7 × 10⁻⁴. The per-channel integrals
also agree with the bin 10–100 ePIC log, within errors:

| channel | ePIC log | this build |
|---|---|---|
| INC31 | 21.175 | 21.1718 |
| INC32 | 14.803 | 14.7799 |
| INC33 | 0.4452 | 0.44486 |

## 4. The isolating runs

The 32 runs are listed in the vendored file's `runs`, with each card verbatim
under `cards`. Each run took 0.1–15.1 s, and the set took 1 min 53 s wall on
2 jobs. Every run exited 0. No tail exceeded the 1 % error rule, so nothing
was rerun.

| bin [GeV²] | xmax | Born, Rad=0 [nb] | IEL31 | IEL32 | IEL33 | tail [nb] | rel. err. | tail / Born |
|---|---|---|---|---|---|---|---|---|
| 1–10 | 0.1096 | 574.056 | 8.7736 | 1.1588 | 13.448 | 23.3809 | 3.7 × 10⁻⁴ | 0.040729 |
| 10–100 | 0.5518 | 67.6268 | 0.14430 | 1.28 × 10⁻⁴ | 0.96704 | 1.11147 | 4.8 × 10⁻⁴ | 0.016435 |
| 100–10³ | 0.9249 | 3.22611 | 5.578 × 10⁻⁴ | 0 | 0.054605 | 0.0551624 | 4.4 × 10⁻⁴ | 0.017099 |
| 10³–10⁴ | 0.95 | 0.0693644 | 3.173 × 10⁻⁶ | 0 | 0.0015151 | 0.00151826 | 3.9 × 10⁻⁴ | 0.021888 |

Observations, all recorded in the file:

- **IEL32 (the elastic FSR mapping) integrates to exactly 0 at Q2MIN ≥ 100
  GeV².** VEGAS stops after one iteration. The cause, found by the
  verification of 2026-10-26 in the source: `HSELK2`
  (`djangoh_h-4.6.21.f` lines 17470–17474) returns 0 at every point whose
  hadronic-t range exceeds GS00 = 20 GeV² (`IF(|TSMIN|>GS00 .OR.
  |TSMAX|>GS00) HSELK2=0`), an undocumented gate — not "the dipole is
  dead", as this record first said. A rebuild with GS00 = 10⁵ changes IEL32
  by +0.0146 / +3.0 × 10⁻⁴ / +2.8 × 10⁻⁶ / +2.6 × 10⁻⁸ nb, i.e.
  +6.2 × 10⁻⁴ / +2.7 × 10⁻⁴ / +5 × 10⁻⁵ / +1.7 × 10⁻⁵ of the bin tail: below
  every allowance, so the vendored numbers and the verdict stand.
- **Q2MAX.** A tail run at Q2MAX = Q²_hi returns the Q2MAX = 10⁸ integral
  to every digit, in all four bins. That confirms the source reading behind
  the subtraction. A Born run at Q2MAX = Q²_hi agrees with the subtraction
  to ≤ 3 × 10⁻⁵, the size of its 5-digit print.
- **Running α** (LPARIN(3) = 3, the ePIC value) multiplies the tail by
  ×1.0152 / 1.0133 / 1.0151 / 1.0206. This is vacuum polarisation at the
  hadronic vertex.
- **The VEGAS log parser was fixed after the first vendoring.**
  - LHAPDF's C++ "Thanks for using LHAPDF" line lands inside the Fortran
    output at a buffer flush. In 12 runs it hid one iteration from the first
    regex, so those runs recorded an iteration count one short.
  - The final integrals, standard deviations and χ² were already right. A
    line-by-line re-parse of all 40 logs reproduced every one of them
    exactly.
  - The count field was corrected with `djrun.py reparse` and the file
    re-vendored. Its sha256 is now `8c0bd1ee…31c7`.

## 5. The comparison, and where the residual comes from

Command: `cd /home/cpeng/Projects/polli/PoliGen && source env.sh && python3 validation/benchmarks/t4_djangoh_rad_noRad.py`.
It took 22.2 / 24.0 / 24.7 / 25.3 s, and 25.7 s with `--plot`. The machine
is shared with another worker; the `polrad-full` tail table took 21 s of
that under load.

**The integration.** Gauss–Legendre on 4 × 12 nodes in ln Q² and 2 × 6 × 12
in y. The y nodes are in ln y below 0.5 and in −ln(1 − y) above it. The Born
falling outside the tail table is 0 in every bin and model.

| bin [GeV²] | f DJANGOH | f polrad-full | R polrad-full | f t-peak | R t-peak | allowance | verdict |
|---|---|---|---|---|---|---|---|
| 1–10 | 0.040729 | 0.041539 | 1.019891 | 0.024919 | 0.611819 | 0.030747 | PASS / PASS |
| 10–100 | 0.016435 | 0.016865 | 1.026135 | 0.014739 | 0.896808 | 0.030960 | PASS / PASS |
| 100–10³ | 0.017099 | 0.017474 | 1.021919 | 0.017290 | 1.011156 | 0.030883 | PASS / PASS |
| 10³–10⁴ | 0.021888 | 0.022102 | 1.009787 | 0.022036 | 1.006746 | 0.030775 | PASS / PASS |

**The diagnostic** (`… t4_djangoh_rad_noRad.py --diagnostics`, 93.8 s and
then 109.2 s) recomputes the same window integral two more ways:

- **direct**: the tail taken from the quadrature at every node
  (`RcModel.tail_sigma_at`), with no (x, y) table and no clip;
- **interpolated tail only**: the table's tail column over the exact Born,
  with no clip.

| bin | table / direct, polrad-full | of which tail interpolation | of which Born interpolation and clip | R direct, polrad-full | table / direct, t-peak | R direct, t-peak |
|---|---|---|---|---|---|---|
| 1–10 | 1.01949 | 1.02029 | 0.99922 | **1.00039** | 1.02096 | 0.59926 |
| 10–100 | 1.02508 | 1.02517 | 0.99991 | **1.00103** | 1.02171 | 0.87775 |
| 100–10³ | 1.02242 | 1.02259 | 0.99984 | **0.99951** | 1.02190 | 0.98948 |
| 10³–10⁴ | 1.02349 | 1.03652 | 0.98743 | **0.98661** | 1.02244 | 0.98465 |

Doubling the panels moves the table-route window integrals by at most
9.2 × 10⁻⁵ (t-peak) and 1.8 × 10⁻⁴ (polrad-full); the direct route moves
by 7 × 10⁻¹⁵. What the diagnostic shows:

1. **The two exact calculations agree in bins 1–3, to ≤ 0.1 %; bin 4 is
   −1.34 %.** LiPolGen's POLRAD Eq. (18) proton tail is 1.0004 / 1.0010 /
   0.9995 / 0.9866 of DJANGOH's HERACLES tail. The bin-4 residual is about
   34 of DJANGOH's δ and is not numerical on either side as far as tested
   (verification of 2026-10-26): LiPolGen's η quadrature at n_η 128 → 256
   gives 0.986608 → 0.986608 and the direct window quadrature with its
   panels doubled moves by 7 × 10⁻¹⁵; DJANGOH at INT-POINTS 200 000 moves
   bin 4 by −2 × 10⁻⁴; `HSELK2`'s gate accounts for 1.7 × 10⁻⁵. Its cause
   is unknown; it is recorded, not explained.
2. **The whole +2 % is LiPolGen's tail table.** The tail column's
   log-bilinear interpolation between nodes gives ×1.020 … ×1.037. The
   proton table has 101 x edges and 77 y rows: about 18 log rows per decade
   plus the 0.9…0.985 refinement rows. The Born's interpolation contributes
   ≤ 0.1 %.
3. **The clip matters only in bin 4.** There the `tail_max` = 10 clip takes
   back 1.3 % of polrad-full's tail.
4. **t-peak's R > 1 in bins 3–4 is the table's.** Read directly, t-peak is
   0.989 and 0.985 of the exact tail, below it as a lower bound should be.
5. **The table's shift exceeds the declared 2 % table term in bins 2–4.**
   The 2 % was declared before this was known; the measured table/direct
   ratio is +1.95 / +2.51 / +2.24 / +2.35 %, so bins 2–4 pass through the
   rounding up to 3 %, and bin 2 (1.0261 against 1.0310) nearly used all of
   it. The declared 3 % held; the 2 % sub-term did not.

**The window.** The strip of the W_l² ≥ 9 window that the rectangle leaves
out carries 4.7 % / 0.28 % / 4 × 10⁻⁶ / 5 × 10⁻⁶ of LiPolGen's Born. It
carries ≤ 5.2 × 10⁻³ of the tail (polrad-full, bin 1) and ≤ 10⁻⁵ for
t-peak. The rectangle therefore holds the generator window's tail almost
whole.

**The history** (recorded sub-rows). On the 2026-09-23 W_l² ≥ 9 window,
LiPolGen's ratios reproduce to every pinned digit:

- t-peak: 1.01616918 / 1.01317215 / 1.01497026 / 1.02068610
- polrad-full: 1.02709574 / 1.01507182 / 1.01513763 / 1.02075698

The ePIC Rad/noRad ratios, 1.0725 / 1.0894 / 1.1606 / 1.2280, stay
history. They measure the inelastic correction with the tail off. The figure
draws them grey.

**The figure** is `docs/benchmarking/figures/t4_djangoh_rad_noRad.png`
(60 205 B, sha256 `5bc8fa40…f97d`), written by the `--plot` run at
19:38:23Z. The top panel shows the tail fraction per bin: DJANGOH as ink
points, t-peak as blue diamonds, polrad-full as green squares, the tolerance
band grey, and ePIC Rad/noRad − 1 as grey history. The bottom panel shows
R with the band.

## 6. What it means for ingredient F

The ranking's ingredient-F row (`08_channel_priority.md` §4) called row 7
"the one generator-vs-generator RC check possible". Its previous status was
only POLRAD-internal: "within 0.25×–4× of POLRAD's … the polarized QE tail
is unpriced".

**Checked externally for the first time, with one open residual.** The
unpolarised proton elastic tail of `PolradFull` now has an external exact
O(α) comparison — at one beam setting, with one shared dipole form factor.

- **The exact formula agrees in bins 1–3 to ≤ 0.1 %, and is 1.34 % low in
  bin 4.** `polrad_full_sigma_qe_u`, POLRAD Eq. (18) in the Eq. (44) slot,
  straight from the quadrature, is 1.0004 / 1.0010 / 0.9995 / 0.9866 of
  DJANGOH's IEL31–33 over Q² = 1–10⁴ GeV² at 18 × 275. The bin-4 residual
  is not numerical on either side as far as tested (§5); its cause is
  unknown.
- **This is the core of ⁶Li's quasi-elastic tail.** That tail is the same
  formula, Z × proton + N × neutron, with the Pauli factor. It is 22 % /
  73 % / 99.9 % of the ⁶Li tail at x = 0.01 / 0.10 / 0.30, per
  `include/lipolgen/rc.hpp`'s header.
- **What a reader can now say:** `PolradFull`'s proton-slot normalisation
  has been compared with an external exact tail: within 0.1 % at
  Q² < 10³ GeV², 1.3 % low at 10³–10⁴ GeV². Until now `rc.hpp`'s
  `RcTailModel::PolradFull` comment, point (3), said no such check existed
  (the comment is now partly stale, §10). It does not say the formula is
  right everywhere: one beam, one form factor, no coherent or polarised
  term, and the shipped table adds +2 % on top.

**Two things this prices without moving anything.** Both are filed as
candidates below.

1. **The shipped default `TPeak` under-weights the proton tail at low Q².**
   It captures 60 % / 88 % / 99 % / 98 % of the exact tail (read directly).
   So the shipped `rc_tail` dilution on a proton is 39 % low at Q² = 1–10
   and 10 % low at 10–100. This is the documented "lower bound" behaviour,
   now with a size from outside the tree.
2. **The tail table lifts the tail by +1.9 … +2.5 % once integrated over
   a bin.** The tail column's log-bilinear interpolation alone gives
   +2.0 … +3.7 %; the clip takes some back in bin 4. It applies to both
   models, at this table resolution.

The coherent spin-1 ⁶Li tail is untouched (§7).

**For the letter's channels.** `rc_tail` is an unpolarised dilution of
the inclusive asymmetries, per `rc.hpp`'s per-channel rule (inclusive:
band and tail). This row prices the size of its elastic-tail part on the
proton: the shipped `TPeak` under-weights it most at low Q², and the table
over-weights it by about 2 %. No ⁶Li number is computed or moved here.

## 7. What it does not validate

- **Other tails.**
  - The coherent spin-1 tail: zero on a proton.
  - The Pauli-suppressed quasi-elastic tail (`qe_kf_gev` = 0 here).
  - Any tensor or polarised tail.
- **Other corrections.**
  - The inelastic O(α) correction. LiPolGen applies none; the ePIC history
    is its size, 7–23 %.
  - Proton-line radiation and two-photon exchange: neither side has them.
  - The running of α: recorded only.
- **The common inputs.** Both sides use the same dipole form factor, so it
  is not checked.
- **Outside the window.** Q² < 1 GeV², y > 0.985, x_l > XMAX_b (the strip),
  and bin 4's x_l > 0.95.
- **The generator itself.** Its event sampling and acceptance are not
  exercised: this is a window integral of the library's weight.
- **Other models and targets.** `t-peak+ll`, other beam configurations, and
  ⁶Li or ⁷Li targets.

## 8. Files, commands, runtimes

| file | what |
|---|---|
| `validation/benchmarks/t4_djangoh_rad_noRad.py` | **Re-wired.** The design block is the tolerance. The primary comparison is the isolated tail. The ePIC rows stay as recorded history, under `history` and `history_epic[…]` sub-rows. `plot()` and `--plot`, `--diagnostics`. It prints `REPORT \|` instead of `ROW \|`. 24 sub-rows; `EXPENSIVE = False`. File sha256 `58b306e1…c0a8`. |
| `validation/benchmarks/data/djangoh_eltail_ep18x275.json` | **New.** 134 263 B, sha256 `8c0bd1ee…31c7`. Provenance (§2), window, the 4 bins, the 8 control numbers, 40 run records with each log's sha256, and 40 cards verbatim. |
| `python/tests/test_bench_djangoh_eltail.py` | **New**, 17 tests (§9). |
| `docs/benchmarking/figures/t4_djangoh_rad_noRad.png` | **New**, written once by `--plot`. |
| this file | **New.** |

Outside the repository, in `$SCR`:

- **Kept:** `src/`, the clone at e356b84 plus the manual; `build/`;
  `lhapdf_dl/` and `lhapdf_sets/`; `epic/`, the ePIC runcards, logs and
  README at 9869d9a, already present and hash-checked; `runs/<name>/`, with
  each card, `stdout.log`, `outfiles/*_out.dat` and `done.json`.
- **Tools:** `tools/djrun.py`, the driver (`control`, `isolate`, `reparse`),
  and `tools/vendor_eltail.py`.
- **Design snapshot:** `design_docstring.txt`, the design block as written,
  and `t4_header_at_design_time.py`, a copy of the harness at 19:11:53Z with
  sha256 `f3dc05ba…6255`.

```
cd $SCR/tools && python3 djrun.py control      # 19:13:35-19:14:22Z
cd $SCR/tools && python3 djrun.py isolate      # 19:14:44-19:16:37Z
cd $SCR/tools && python3 djrun.py reparse      # the parser fix of sec. 4, no rerun
cd $SCR/tools && python3 vendor_eltail.py /home/cpeng/Projects/polli/PoliGen
cd /home/cpeng/Projects/polli/PoliGen && source env.sh
python3 validation/benchmarks/t4_djangoh_rad_noRad.py                 # 22-26 s, PASS
python3 validation/benchmarks/t4_djangoh_rad_noRad.py --plot          # 25.7 s
python3 validation/benchmarks/t4_djangoh_rad_noRad.py --diagnostics   # 109.2 s
python3 validation/benchmarks/make_report.py --only t4_djangoh_rad_noRad \
    --out-dir /tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/report_B3
                                    # 24.4 s: 1 pass, 0 fail, 24 sub-rows rendered
python3 validation/check_spdx_headers.py                              # 130 files, all carry the SPDX line
```

## 9. Tests

```
python3 -m pytest -p no:cacheprovider python/tests/test_bench_djangoh_eltail.py -q
    16 passed, 1 skipped (the opt-in diagnostic) in 42.2-48.9 s (three runs)
LIPOLGEN_BENCH=1 python3 -m pytest -p no:cacheprovider python/tests/test_bench_djangoh_eltail.py -q
    17 passed in 148.6 s
```

What the 17 tests cover:

- **File guards:** point counts, the first and last value, and one flipped
  byte refused.
- **Provenance:** every card checked against its sha256.
- **The design as read off the 32 cards.**
- **The design block:** its hash; every run started after 19:11:53Z; one
  edited character refused.
- **The window rule, the control and the Q2MAX facts.**
- **The verdict:** it follows from the numbers.
- **The headline numbers:** R pinned at rel 2 × 10⁻⁷.
- **The history:** the 2026-09-23 pins reproduce.
- **The REPORT line and `_report.normalise_row`.**
- **The PNG.**
- **Exit codes:** 0 on a run; 2 on an altered data file, on an edited
  design block, and without numpy.
- **Opt-in:** the diagnostic numbers.

The shared `python/tests/test_bench_chain_rc.py` was **not edited**. Run
against the re-wired harness, it gives 2 passed and 2 failed:

- `test_djangoh_row_is_blocked_and_prints_both_sides` fails, because the row
  is no longer BLOCKED;
- `test_djangoh_row_numbers_are_the_recorded_ones` fails, because it reads
  `generator_value`, which is now the R values.

A replacement for both, plus the docstring, is in
`$SCR/integrator/test_bench_chain_rc.patch`. It pins the same 2026-09-23
numbers on the new `history` block. Run on a scratch copy of the file, it
gives 4 passed.

## 10. Integration status, the open questions and the priced candidates

*Added by the supervising session on 2026-10-26 after the worker returned.*

**Applied in the B3 commit on the run branch** (so the run's integrator does
not repeat them): `python/tests/test_bench_chain_rc.py` re-pointed with the
worker's replacement (4 passed); `BENCHMARK_PLAN.md` §4 row 7 (both cells);
`validation/benchmarks/README.md` (convention 4, the blocked-rows and opt-in
paragraphs, the row-7 line of "the plan's first ten");
`08_channel_priority.md` §4 row F and §5 item 2; `REFERENCES.md` (the DJANGOH
source and manual, read and kept outside); `STATUS.md` row B3.

**Left for the integrator (phase E):** `README.md`'s "As of" tallies, from
the regenerated `REPORT.json`; the regeneration itself (`--all --plots`),
where row 7 goes blocked → pass and `data_sha256` gains
`data/djangoh_eltail_ep18x275.json` = `8c0bd1ee…31c7`; the two priced
candidates below into `DECISIONS_PRICED.md`; and, optionally, one sentence in
`docs/PHYSICS_CHANNELS.md` after the RC bullet that begins "The exact Mo–Tsai
elastic radiative tail is not implemented": *"Since 2026-10-26
(BENCHMARK_PLAN.md §4 row 7) `PolradFull`'s proton (quasi-elastic) slot is
checked against an external exact O(α) tail, DJANGOH 4.6.21's IEL31–33:
−1.3 … +0.1 % straight from the quadrature, +1.0 … +2.6 % through the
shipped (x, y) table; `TPeak` is 0.60 / 0.88 / 0.99 / 0.98 of it at
Q² = 1–10 / 10–100 / 100–10³ / 10³–10⁴ GeV² (e 18 × p 275)."* — only with
the line-citation gate re-run, since lines shift.

**Left for the author:** two comments in `include/lipolgen/rc.hpp` are now
partly stale for the proton/quasi-elastic slot — `RcTailModel::PolradFull`
point (3), "ITS NORMALISATION IS STILL NOT CHECKED AGAINST MO-TSAI, or
against any external exact tail", and the first bullet of "FIVE THINGS THIS
IS NOT"; a comment-only edit, left alone here because the PHYSICS_CHANNELS
gate cites lines of that file. The licence question below.

### 10.1 Priced candidates (none taken; for `DECISIONS_PRICED.md`)

1. **The `rc_tail` default model, `TPeak` (shipped) vs `PolradFull`
   (opt-in).** On a proton at e 18 × p 275 in row 7's matched window,
   straight from the quadrature (exact = DJANGOH 4.6.21): TPeak/exact =
   0.599 / 0.878 / 0.989 / 0.985 and PolradFull/exact = 1.000 / 1.001 /
   1.000 / 0.987 at Q² = 1–10 / 10–100 / 100–10³ / 10³–10⁴ GeV²; through the
   shipped table TPeak 0.612 / 0.897 / 1.011 / 1.007, PolradFull 1.020 /
   1.026 / 1.022 / 1.010. Cost of switching: the tail table build goes from
   ~2 s to 10–21 s, and every published number moves (`rc.hpp`
   `RcTailModel::PolradFull` point (6)). No registry row is priced directly;
   row 17 (the RC band's low-x anchor) is untouched.
2. **The RC tail table's interpolation.** The log-bilinear (ln x, ln y)
   interpolation of the tail column on the 101 × 77 node grid raises the
   window-integrated tail by +2.0 … +3.7 % (+1.9 … +2.5 % net of the
   `tail_max` clip) on the proton at 18 × 275, for both models. Options:
   keep it, refine the y rows, or interpolate a smoother variable. The cost
   of a fix was not measured. It is likely present on ⁶Li's quasi-elastic
   column at the same resolution; not measured.

### 10.2 Open questions (the worker's, verbatim in substance)

- The tail table's bias (+2 %) is the whole difference between
  `polrad-full` and DJANGOH; fixing it is a library change for the author.
- Whether the default should become `PolradFull` is a registry question
  (candidate 1).
- The four control Born numbers sit +0.2 … +2.5 × 10⁻⁴ above ePIC's,
  attributed to the LHAPDF 5.9.1 → 6.5.5 interpolation of `cteq61` but not
  proven (no LHAPDF 5 here); within the declared 3σ rule.
- *(closed by the verification, §11)* IEL32 = 0 at Q2MIN ≥ 100 GeV² is
  `HSELK2`'s GS00 = 20 GeV² gate on the hadronic-t range (§4); its effect
  on the bin tails is ≤ 6.2 × 10⁻⁴.
- **Licence:** neither `spiesber/DJANGOH` nor `eic/InclusiveDjangohSamples`
  states one. The vendored file carries numbers plus all 40 input cards
  verbatim, derived from the ePIC runcards, reproduced as facts with
  attribution. The author may want to confirm this.
- Runtime margin: 22–26 s on a shared machine against the 30 s opt-in line;
  the `polrad-full` table build (21 s here under load, ~10 s unloaded) could
  push a loaded machine past it. `EXPENSIVE` stays False on the measurements.
- Observation, no action: the tail table has two y nodes 1 × 10⁻⁹ apart at
  0.985 (y_max·(1 + 10⁻⁹) and the 0.985 refinement row); `std::unique`'s
  10⁻¹² tolerance does not merge them, leaving a harmless zero-width cell.

## 11. Independent verification (2026-10-26) and what it changed

An adversarial reader (Opus, four lenses: physics, pre-registration
honesty, code, regression) was run on the row after the worker returned,
before the B3 commit `372b915` and this record's fixes. Its verdict: **the
PASS stands as a measurement under its declared tolerance**, with caveats.
What it confirmed by reading the DJANGOH source and by running: IEL31–33
together are the full lepton-line matrix element (ISR + FSR + interference +
Compton) with m_e kept and F₁, F₂ at −t_s; IEL2 is `HSEL22` at x = 1; the
IEL integrands ignore Q2MAX; WMIN is raised to M_p and the elastic channels
are zeroed only above it; LPARIN(3) = 2 fixes α; the windows are the same
rectangle on both sides; the generator applies `1 + tail_ratio_at/(1 +
w_avg)` so ∫ Born · r is its absolute tail; the design block hashes to the
pinned digest in the harness, the data file and the design-time copy
(ctime = mtime, 19:11:53Z), no DJANGOH output exists between the build and
the first run, no draft carries another tolerance, and every number of §0
and of the REPORT line reproduces to every digit. Findings and
dispositions:

| # | finding | disposition |
|---|---|---|
| 1 | the executable tolerance (`TOL_NUMERICAL`, `N_DELTA`, `verdict()`, the tolerance string) sat outside the hashed block: a copy with 0.05 ran, passed the hash check and still printed "0.03 + 2 delta_b"; no sign of tuning (every draft 0.03 / 2.0; the first harness run at 19:27:45Z printed the final R) | **fixed:** `design_block_sha256()` now parses the block's two "PASS(…) iff" sentences and exits 2 if the constants differ; the tolerance text is built from the constants; the four allowances 0.0307467 / 0.0309598 / 0.0308831 / 0.0307745 and the edited-constant refusal are pinned in tests |
| 2 | the cause given for IEL32 = 0 ("the dipole is dead") was wrong: `HSELK2` has an undocumented GS00 = 20 GeV² gate on the hadronic-t range; a patched rebuild sizes its effect at ≤ 6.2 × 10⁻⁴ of a bin tail | **fixed** in §4, §10.2 and the harness addenda; the verdict unchanged |
| 3 | bin 4's −1.34 % direct residual is real (34 δ) and unexplained; "now validated", "the exact formula is right", "the two exact calculations agree" were overstated | **fixed:** §5 item 1 and §6 now state 0.1 % in bins 1–3, −1.34 % in bin 4, cause unknown, one beam, one form factor |
| 4 | "the 2 % table allowance was nearly all used" was wrong: the table's shift (+1.95 / +2.51 / +2.24 / +2.35 %) exceeds the 2 % term in bins 2–4; bin 2 passes through the rounding to 3 % | **fixed** in §5 item 5 and the harness header |
| 5 | "this DJANGOH is the ePIC DJANGOH" overstated: every control card runs IEL2–IEL33 off, so the elastic-tail code has no external control; LHAPDF 5.9.1 → 6.5.5 differs | **fixed** in §0 and the harness header |
| 6 | the t-peak half of the PASS cannot fail low (pre-registered as a one-sided ceiling, so no rule broken), and through the table t-peak sits above the exact tail in bins 3–4 | **wording fixed** in the REPORT reason, the plan row, the README row and STATUS: PASS = polrad-full two-sided; t-peak ceiling only, 39 % / 10 % low at Q² < 100 GeV² |
| 7 | design block (b)'s "three phase-space mappings of ONE squared matrix element" is loose (1/A1, 1/A2, 1/t_s partial fractions, interference through 1/(A1 + A2)) | addendum under the block (the block is hashed); the conclusion stands |
| 8 | the pre-declared fallback "rerun at INT-POINTS 400000" cannot run on this build (VEGAS's grid array is dimensioned 50; 200 000 runs); never triggered | noted in the harness addenda |
| 9 | unrelated: `test_bench_report.py::test_a_diff_maps_to_tiers_and_harnesses` fails since `82b3237` (`_report.py` now maps to the harnesses that import `import_pythia8`) | **fixed:** the test now expects exactly the harnesses that name `_report` |

The verifier's scripts, reruns and the GS00-patched build are in the
session scratchpad (`verify_b3/`), outside the repository.
