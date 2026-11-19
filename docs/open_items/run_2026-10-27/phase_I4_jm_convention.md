<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-27 — phase I4: the Jaffe–Manohar normalisation of Δ implemented (decision D4)

Worker I4, 2026-10-27, on `claude/run-2026-10-26` at `ca9b23a`; HEAD moved to
`fbc24f2` during the phase (the supervisor's commit of this run's records,
which touches no file this phase edits). Uncommitted: the supervising session
commits. Decision: `DECISIONS.md` D4. Design:
`plan_D4_jm_convention.md`, AMENDED by its adversarial check (session scratch
`checks_1027/D4.json`, reproduced at the end of the plan); literature:
`phase_R1_delta_convention.md`. Neither record was edited here (a
records-hygiene worker owns their citation forms).

Conventions: **measured** = run on the D4 build of this phase; **pre-D4** =
the `ca9b23a` value (the design's measurement, re-confirmed where this record
says so). Line citations are to the working tree as this record leaves it.

> **A duplicate instance.** A second copy of this worker ran in the same tree
> from about 18:20 to 18:30. It says so in a notice it left in session scratch
> (`i4/DUPLICATE_INSTANCE_NOTICE.txt`; it also made this worker's scratch edit
> helper print that notice). By its own account it renamed one local variable
> in `tests/test_reference.cpp`, relinked the test binary, reverted the rename
> and stopped. Nothing it wrote was taken on trust: the working tree was
> rebuilt from `git archive ca9b23a` plus this worker's own edit scripts and
> compared file by file (identical, all 25 files then changed), and the
> binaries were rebuilt from the sources before the final gates.

## 0. What D4 changed, in one paragraph

The cos 2φ coefficient of the double-helicity-flip term now enters through one
named constant, `include/lipolgen/constants.hpp:158` `DELTA_COS2PHI_NORM` = 0.5,
at the two places the library evaluates it, `src/core/xsec.cpp:282`
`DELTA_COS2PHI_NORM` (the kernel's Δ line) and `src/core/asymmetries.cpp:91`
`DELTA_COS2PHI_NORM` (`a_cos2phi`), so a₂ = −½(1−y)/y² c_eff(m) sin²θ_S Δ/D_φ,
Jaffe–Manohar's. Applied as a factor 2⁻¹ it makes every `a_cos2phi`, and every
a₂ on the massless path (where a₂ is the Δ line alone), exactly half its
pre-D4 value, bit for bit; on the opt-in `tensor_gamma` path, where
a₂ = h₂ + the Δ line, the Δ line halves and h₂ does not (§10, V-6). It touches
nothing else (§2). No other physics constant, registry label, reference-JSON
number or vendored data file moved, and no tolerance was loosened (the C++
port gate's, loosened here, was restored by F4: §10, V-4): `toy_delta_gluon`'s
1e-2 scale, `C_BAG` = −0.012 (which now enters unconverted) and
`make_moment_a/b` are as they were. Provenance: `meta["delta_convention"]` in
the `meta` of `Pipeline.generate()` (`python/bindings.cpp:647` `delta_convention`),
so on every .npz the CLI writes, but on no HepMC3 file (§10, V-2); the same string in the
CLI run banner's "Delta (cos 2phi) convention" line (`python/lipolgen/cli.py:1777`, `say(… % _l.DELTA_CONVENTION)`); and
`meta["rank2_input"]` names the Jaffe–Manohar Δ (`src/core/pipeline.cpp:1292`
"Jaffe-Manohar Delta").

## 1. Every number that moved, before → after

| file:line | what | pre-D4 | D4 |
|---|---|---|---|
| `include/lipolgen/constants.hpp:158` | `DELTA_COS2PHI_NORM` (new) | (1, implicit) | 0.5 |
| `tests/test_xsec.cpp:457` | T1, the closed form a₂/(c_eff sin²θ_S), `eps_of(p.y)` × `d_over_f1` | −(ε/2)(Δ/F₁)/(1+εR) | −(ε/4)(Δ/F₁)/(1+εR) |
| `tests/test_xsec.cpp:819` | T2, the oversized Δ (`big_opt.delta_func`) that must drive W negative at m = +1 | 3.0 F₁ | 6.0 F₁ (the same a₂, bit for bit) |
| `tests/test_rc.cpp:207` | T3, the written-out Δ term of the `RcScope` identity (× `tensor_moments(m)`) | −(1−y)/(y y) … | −0.5 (1−y)/(y y) … |
| `tests/test_tensor_gamma.cpp:361` | T4, the leakage as an equivalent Δ at three points | 0.1628248344088873 / 0.16006748610799001 / 0.15477697871375315 | 0.3256496688177746 / 0.32013497221598003 / 0.3095539574275063 (exactly ×2) |
| `tests/test_tensor_gamma.cpp:371-375` | T4, `d_fake` and the band on `coef` | −h₂ y² D_φ/(1−y); 0.14 < coef < 0.17 | −2 h₂ y² D_φ/(1−y) (JM units); 0.28 < coef < 0.34 |
| `tests/test_reference.cpp:386-387` | T5, the factor the port gate applies to the reference's Δ line | none (plain 1e-12) | 0.5/c_ref, c_ref read off `xsec.json` (§2 item 3) |
| `python/tests/test_bench_t5_delta_bound.py:292` | Y1, the tree's coefficient (`pytest.approx`, rel 1e-15) | −(1−y)/y² | −0.5 (1−y)/y² |
| `python/tests/test_bench_t5_delta_bound.py:300` | Y1, DS16 Eq. (5) over the tree at Q² = 5 (`pytest.approx`) | ≈ 0.5 (rel 0.02) | ≈ 1.0 (rel 0.02; 0.99214) |
| `python/tests/test_bench_t5_delta_bound.py:435` | the missed pin, (`ds16_over_tree`, `pre_d4_over_jm`) | (2.0, —) | (1.0, 2.0) |
| `python/tests/test_bench_t5_delta_bound.py:500` | Y4, min W over 1701 rates | 0.9493997515 | 0.9493997990 |
| `python/tests/test_bench_t5_delta_bound.py:502` | Y4, A_T max | 0.05060024846 | 0.05060020096 |
| `python/tests/test_bench_t5_delta_bound.py:512` | Y5, V6 = factor × central | `2.0 * s` | `mod.PRE_D4_OVER_JM * s` (= 2.0, asserted) |
| `python/tests/test_price_defaults.py:44` | Y9, `A2_TOY` S1…S4 per unit P_zz | −1.302913804135417e−03 / −1.0625622769273495e−03 / −1.4290788342523734e−03 / −1.3290877454563259e−03 | −6.514569020677085e−04 / −5.312811384636748e−04 / −7.145394171261867e−04 / −6.645438727281629e−04 (exactly ½) |
| `validation/benchmarks/t5_delta_moment_bound.py:329` | `DS16_OVER_TREE`, re-derived by the pytest | 2.0 | 1.0 |
| `validation/benchmarks/t5_delta_moment_bound.py:1021` | `PRE_D4_OVER_JM` (new): V6's factor, kept as declared | — | 2.0 |

Re-measured on the D4 build, as predicted: the four sweet-spot a₂ of
`default_inclusive_kernel(li6())` (m = +1, θ_S = 90°, s = 3980 GeV²) at
(x, Q²) = (0.0112, 1.14), (0.0282, 1.14), (0.0708, 3.13), (0.141, 14.3) go from
−1.055948214152054e−03 / −1.299961141425631e−03 / −1.440230754221272e−03 /
−1.326442680639959e−03 to **−5.279741070760268e−04 / −6.499805707128153e−04 /
−7.201153771106358e−04 / −6.632213403199796e−04**, equal to `lg.a_cos2phi`
bit for bit, with m = 0 carrying exactly −2 × m = +1 (pinned by N2,
`python/tests/test_delta_convention.py`). Row 13 under its declared rule: M₂,
B and |M₂|/B = 0.0644 / 0.0803 / 0.0941 bit-identical; V6 0.1289 / 0.1605 /
0.1881 and the TIGHT edge 0.2719 / 0.3369 / 0.3937 numerically unchanged (V6
kept in TIGHT, relabelled: amendment 5); PASS → PASS.

Text that moved with the numbers (no number): the t5 harness's V6 label
("the pre-D4 kernel / the letter as printed, x2"), its REPORT line, the letter
rows' "RECORDED: PASS centrally (A7, the kernel's own Delta since D4), FAIL
under V6 (…)", the figure's legend; the comments of `tests/test_tensor_gamma.cpp`
at :307 and :339-354 (JM units).

## 2. What was verified, and how

1. **The exact-half contract.** The design's checker (session scratch
   `d4/d4_before_after.py check pre_d4_ca9b23a.json`; the pre-D4 dump's sha256
   `0334237e…2c8ba3`) on the D4 build: "D4 contract holds: 5310 states exact
   half, 1593 tensor_gamma states half on the Delta line, every Delta-free
   quantity bit-identical".
2. **The finite-γ statement** (`d4/measure_gamma2.py`, re-run): the Δ line is
   an F[U T_TT] insertion of −1.000000000 Δ at γ² = 3.5 × 10⁻¹⁰ (pre-D4
   −2.000000001), −1.000392377 at (x, y, Q²) = (0.01, 0.5, 1.14) (pre-D4
   −2.000784754), −1.003673920 at the spot-1 kinematics (pre-D4 −2.007347840);
   the b-sector's h₂ per unit F_TTT is unchanged to every printed digit.
3. **The port gate cannot pass for the wrong reason.** Both gates read
   polligen's coefficient off `xsec.json`, c_ref = a_cos2phi D_φ y²/(−(1−y)Δ) at
   every Δ ≠ 0 point: 1 in all five Δ-carrying blocks (⁶Li ×3, ⁷Li ×2; max
   |c − 1| = 2.2e−16), absent in the two Δ-free ⁷Li blocks. The tree is
   compared with a LITERAL 0.5 (`python/tests/test_kernel_reference.py:69`
   `JM_COEFFICIENT`; `tests/test_reference.cpp:386-387`, "0.5 / c_ref"), never with
   `DELTA_COS2PHI_NORM`. `python/tests/test_kernel_reference.py:241`
   (`test_the_port_gate_reads_the_coefficient_off_the_file`) pins the
   reader: a copy with `a_cos2phi` ×0.7 is REFUSED ("unknown Delta
   normalisation"); a copy transformed to Jaffe–Manohar's coefficient reads
   c_ref = ½, factor 1, and the tree meets it with no rescaling — the state
   after the sibling changes, with no code change (for the whole Python
   module only since F4, which removed three pins of today's c_ref = 1:
   §10, V-1).
4. **The dump preflight** (`validation/dump_polligen_reference.py:1560`
   `_preflight_delta_convention`, called from `main()` before anything is
   written), exercised alone (no dump, no file written, `python3 -B`): prints
   "polligen's cos 2phi coefficient: 1 (the LOI's Eq. (9); …)"; with
   `asymmetries.a_cos2phi` alone halved (a half-updated sibling) it REFUSES
   ("0.5 (asymmetries.a_cos2phi) and 1.0 (the kernel's a2) …"), and likewise at
   ×0.7. The README it renders: the `xsec.json` section plus the D4 addendum is
   found verbatim in the hand-synced `validation/README.md`.
5. **The coefficient is pinned independently of `a_cos2phi`.** N1
   (`tests/test_xsec.cpp:882`) writes the closed form with a literal ½ and
   c_eff = 3m² − J(J+1) from the spin algebra at m ∈ {+1, 0, −1} (J = 1) and
   m ∈ {±3/2, ±1/2} (J = 3/2; amendment 4 — never ±1, 0 for ⁷Li), θ_S ∈ {0.7,
   π/2}, over the analysis grid at kRtol, plus `DELTA_COS2PHI_NORM == 0.5`. The
   sampler draws φ from that a₂.
6. **A default run is bit-identical** (`tensor-thirds`, θ_S = 0, where the Δ
   term vanishes). The D5 design worker's pre-D4 file (session scratch
   `d5/run5k_rcoff.npz`: inclusive ⁶Li, 5000 events, seed 20261027, P_z 0.7,
   P_zz 0.6, P_e 0.7, made on the `ca9b23a` build) regenerated by the CLI on
   the D4 build: all 48 array columns bit-identical; `meta` differs only by the
   new `delta_convention` key and the `rank2_input` sentence (which the
   `b1_model` knob-provenance row's reason repeats). The banner printed
   "Delta (cos 2phi) convention: jaffe-manohar-1989; cos2phi coefficient
   -(1/2)(1-y)/y^2 (D4, 2026-10-27)".
7. **Adversarial greps** (broadened per the check): no spelled
   `(1.0 - y) / (y * y)` outside tests without the constant beside it (the two
   code sites, and the dump preflight's deliberate pre-D4 unit); no Δ
   coefficient (1−y)/y² without the ½ outside dated records, except where the
   line says it is the sibling's or the pre-D4 one; no bare "0.14–0.16" leakage
   without its pre-D4 label (the one hit, `validation/benchmarks/README.md` row
   9, is a runtime in seconds); "Δ_DS16 = 2 Δ_tree" survives only in the t5
   header's dated first- and second-pass text, now labelled.

## 3. The files

**Library.** `include/lipolgen/constants.hpp` (the constant and
`DELTA_CONVENTION` after `C_BAG`, which is untouched — its line 139 is an
allow pin; the `TENSOR_LL_SIGN` doc at lines 26–30 now says the +1
transcription is the sibling's slide of the JLab LOI's Eq. (9), long called
HJM's — line-neutral); `src/core/xsec.cpp` (the Δ line and its comment, in
place); `src/core/asymmetries.cpp`; `include/lipolgen/xsec.hpp` (master
formula ½, the c_m attribution, the leakage "0.109 % of the letter's printed
amplitude, 0.22 % of the same Δ's since D4" — in place); `include/lipolgen/
asymmetries.hpp` (HJM for the b-sector, Jaffe–Manohar PLB 223 Eq. (16) for the
Δ term, both coefficients ½, in place; line 9, an allow pin, untouched);
`include/lipolgen/sf.hpp` (the double-helicity-flip toy's doc, `include/lipolgen/sf.hpp:428` and `:529`: `scale`
is the prefactor and the peak Δ/F₁ is 0.337 scale at x = 0.07 — the check's
correction, wrong before D4 too); `src/core/pipeline.cpp` (the `rank2_input`
sentence).
**Bindings and CLI.** `DELTA_COS2PHI_NORM` and `DELTA_CONVENTION` are bound at
the END of the module (`python/bindings.cpp:5161` `DELTA_COS2PHI_NORM`), after
every cited line (amendment 1: the design's insertion after `C_BAG` would have
broken SPIN32's citation of the `Amplitudes` class binding (bindings.cpp lines 2187-2199) and five R6 citations);
`meta["delta_convention"]` is written line-neutrally (the comment above
`meta["rank2_input"]` reflowed by one line); the banner line follows the RC
block, after the last line of `cli.py` any document cites.
**Tests.** As §1; new: N1 in `tests/test_xsec.cpp` (417 C++ cases, was 416);
N2 `python/tests/test_delta_convention.py` (6 tests: the constant and its
label, the closed form, every m and angle, the four amplitudes and the m
pattern, `meta["delta_convention"]` of a generated `transverse-tensor` sample,
the CLI banner); the reader test in `test_kernel_reference.py`. The optional
bound-name check (P4) lives in N2, not in `test_module.py`, whose line 265 a
2026-09-02 record cites.
**Harnesses and the dump script** — every edit is line-neutral up to the last
line this run's design and literature records cite, so none of their R6
citations moved: `t5_delta_moment_bound.py` (in place: the MODELS paragraph's
coefficient, S2's "DS16_OVER_TREE = 2" and "the kernel" labelled pre-D4, the
V6 variation relabelled, "it does not decide which of A7 and V6 … — the author
did, D4: A7", `DS16_OVER_TREE` = 1 at line 329, V6 / TIGHT / `f1ll_critical_v6`
on `PRE_D4_OVER_JM`, the REPORT, SUBROW and figure labels; the D4 note itself,
`PRE_D4_OVER_JM` and `V6_LABEL` above `main()`, `validation/benchmarks/
t5_delta_moment_bound.py:993-1023`, with a pointer from the header's TWO PASSES
paragraph); `t1_r_world.py` and `t5_weighted_unweighted.py` docstrings (½);
`dump_polligen_reference.py` (the preflight and the README's D4 addendum
`README_ADDENDA` defined after `main()`, the call placed after its first
statements, `render_readme` appending the addendum in place) and the generated
`validation/README.md` kept in sync by hand.
**Documents.** `docs/CONVENTIONS.md` (the convention bullet at the end, §4; the
leakage 0.28–0.32; the `TENSOR_LL_SIGN` bullet's attribution at lines 82–83);
`docs/PHYSICS_CHANNELS.md` (§2 prose, the rank-2 and `tensor_gamma` rows, the
Δ row with the convention paragraph and its symbols cell citing the Δ line
instead of the `h2` add, the Δ-models row, §13 — all in place — and references
59–62 [JM89g], [DS16], [LOI18], [MWZ13] appended after 58); `README.md`;
SPIN32 Eq. (47)'s a₂ line (in place, no line added); `docs/USAGE.md` (the §2
code comment and the sentence before "Record:", both outside every
fingerprinted range); `docs/benchmarking/04_theory.md` (T-10, T-16–T-18);
`08_channel_priority.md` (rank 1's significances, §3.1, §4 rows A and C, a
revision-log row: no rank moved); `BENCHMARK_PLAN.md` (row 13);
`validation/benchmarks/README.md` (row 13); `docs/references/REFERENCES.md`
(#3, #4, the foundations note); `docs/references/01_tensor-b1.md` (three HJM
sites); `docs/surveys/needs_survey.md` (lines 52, 59, 70, 198, 443, 452, all
annotated — amendment 3, "choose all"); `docs/DEVELOPMENT_PLAN.md` (the box at
line 28, width kept; line 250).

## 4. The convention, as written into `docs/CONVENTIONS.md`

The bullet "Δ normalisation — the cos 2φ coefficient" (the last of "Physics
defaults that are a CHOICE") carries, corrected per the check: the formula and
its m pattern; the sources that derive it (Jaffe–Manohar Eq. (16), read as the
preprint; Sather–Schmidt (A7) quoted correctly — "− ½xλΔ₄ cos 2φ" beside
"+ ⅙(xy²Δ₁ + λΔ₂)"; Ma–Wang–Zhang (23); Detmold–Shanahan (5) by contraction;
the LOI's (6)); Pascalutsa–Vanderhaeghen kept INFERRED; the LOI's Eq. (9) and
the sibling's slide as the origin of the doubled coefficient; HJM for the
b-sector only; **the one-line rule** — a Δ(x, Q²) or a moment quoted in
Jaffe–Manohar, Sather–Schmidt, Detmold–Shanahan or Ma–Wang–Zhang units enters
LiPolGen as printed; a cos 2φ amplitude computed with the pre-D4 kernel or
with the sibling's `a_cos2phi` is twice LiPolGen's for the same Δ — with the
warning not to apply the factor twice; the parton side as NOT one
normalisation (the JM preprint's Eq. (29) α_s/8π against its own Eqs. (20),
(21), (25) and SS90 (3.1)'s α_s/16π; its Eq. (31) against SS90 (2.11); the
sign of Δ); the b-sector reading (SS90 (A2)'s Δ₁ is the helicity-diagonal b₁
the slot holds, so the massless kernel reproduces (A7)/(A8) term by term:
nothing double-counted, nothing to subtract; only the opt-in `tensor_gamma`
path has an O(γ²Δ) basis question); what moves in a generated file (φ, e′, X,
every later draw of the per-event stream, `rc_tensor_lo/hi` — `spin_weight_k`
is filled by no generator) and what does not; provenance; the reversal.
`docs/PHYSICS_CHANNELS.md`'s Δ row carries the same paragraph, condensed, in
place.

## 5. Gates

**Baseline** (`ca9b23a`, before any edit): strict citation gate 1240 refs /
95 ranges / 7 external / 0 broken / 24 allow-listed; SPIN32 19 / 115 / 0 / 0 /
10; PYTHIA_BRIDGE 0 / 0 / 8 / 0 / 0. `--records-only`: 1389 citations in 74
records, 0 broken. `--audit-ranges`: 209 range citations on 179 blocks,
0 REFUSED. SPDX: 133 files. Targeted pytest with `LIPOLGEN_BENCH=1`
(`test_b1_model` 51, `test_bench_t1_r_world` 51, `test_bench_t5_delta_bound`
23, `test_kernel_reference` 10, `test_meta_provenance` 8, `test_module` 89,
`test_price_defaults` 1): 233 passed, 0 skipped.

**Step-1 sanity gate (library only, before any test was touched).**
`./build/lipolgen_tests`: 416 cases, 412 passed, **4 failed — exactly T1**
("the cos 2phi amplitude carries R and NOT F2", 270 checks), **T2** ("production
scenarios keep a healthy positivity margin", 1), **T3** (`tests/test_rc.cpp`,
"xsec: amplitudes() == tensor_amplitudes() + the vector terms", 162) **and T5**
("reference tables: xsec.json": 30 `a_cos2phi`, 816 a₂ = ⁶Li 3 × 24 × 6 + ⁷Li
2 × 32 × 6, 360 dsigma); **T4 passed, as the check predicted** (its table has
Δ = 0). Full pytest: 16 failed / 1440 passed / 155 skipped = exactly Y1, Y4,
the five Y8 port-gate tests, Y9, and the eight `test_doc_link_gate.py` tests
that run the strict gate, red on the predicted MOVED named citation of
`POLEMC_VALENCE_WINDOW_LO` (constants.hpp line 145 at ca9b23a, 168 now) and the predicted EDITED
blocks. The missed pin (the pytest's line 435) is not a library failure: it fails only when the
harness constant moves, as the check said. Records gate at step 1: 0 broken.

**Final gates of I4** (the tree as I4 left it; F4's changes and the gates
after them are in §10):

- `cmake --build build -j2`: clean, no warning (the step-1 near-full rebuild
  took 2 min 21 s). `./build/lipolgen_tests`: 417 cases, 417 passed, 0 skipped,
  17 244 191 assertions, 7 min 33 s (416 before; N1 added).
- `python3 -m pytest -p no:cacheprovider python/tests -q`: 1463 passed,
  155 skipped, 0 failed, 473 s (step 1: 1440 passed, 16 failed, 155 skipped;
  the seven new tests are N2's six and the reader test; the skips are the
  same 155).
  `LIPOLGEN_BENCH=1 … -k bench`: 351 passed, 1 skipped (the nnpdf-data wheel
  is not on this machine: environment, as before), 0 failed, 306 s.
- Targeted files, `LIPOLGEN_BENCH=1`, per file before → after (no test skips
  on either side): `test_b1_model` 51 → 51, `test_bench_t1_r_world` 51 → 51,
  `test_bench_t5_delta_bound` 23 → 23, `test_kernel_reference` 10 → 11 (the
  reader test), `test_meta_provenance` 8 → 8, `test_module` 89 → 89,
  `test_price_defaults` 1 → 1, `test_delta_convention` — → 6 (new),
  `test_doc_link_gate` 112 passed.
- The strict citation gate: 1241 refs / 95 ranges / 7 external / 0 broken /
  24 allow-listed (1240 before: the Δ row's `h2`-add citation replaced by two
  `DELTA_COS2PHI_NORM` ones); SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE
  0 / 0 / 8 / 0 / 0. On the way: `--fix` re-pointed the one MOVED named
  citation (`POLEMC_VALENCE_WINDOW_LO`, `include/lipolgen/constants.hpp:168`);
  the eight EDITED blocks (`include/lipolgen/xsec.hpp` 13-21, 23-29, 197-233,
  204-217; `include/lipolgen/asymmetries.hpp` 9-21; `src/core/xsec.cpp`
  260-284; `docs/CONVENTIONS.md` 75-83 and 81-83) were each re-read against
  every row that cites them (SPIN32 32, 34, 145, 426, 916, 975, 996, 1007, 1168;
  PHYSICS_CHANNELS 82) — each row still says something true — and re-recorded
  by `--record-ranges` (0 refused), which also dropped the fingerprint of the
  Δ row's former citation of src/core/xsec.cpp line 266 (the `tensor_gamma` add
  `out.a2 = out.a2 + h.h2`, not the Δ line).
- `--records-only`: 0 broken (1416 citations in 75 records, this one
  included; 917 checked against a live file, 110 annotated historical,
  199 skipped). `--audit-ranges`: 209 range citations on 179 blocks, 0 REFUSED.
  `check_spdx_headers.py`: 133 files.
- The report: `make_report.py --all --plots --compare <the saved pre-change
  REPORT.json> --out-dir docs/benchmarking` reports 13 pass, 10 fail, 0 blocked
  (as before) and exactly three changes, all row 13:

      MOVED    t5_delta_moment_bound: the generator string, "(V6, Delta_DS16 = 2 Delta_tree: …)"
               -> "(Jaffe-Manohar Delta, the kernel's since D4; V6, the pre-D4 kernel / the
               letter as printed, x2: …)", every number unchanged
      SUBMOVED t5_delta_moment_bound[bound_side[central]]: ds16_over_tree 2 -> 1; pre_d4_over_jm n/a -> 2
      SUBMOVED t5_delta_moment_bound[helicity_positivity[shipped kernel]]: a_t_max 0.0506002 -> 0.0506002;
               w_min 0.9494 -> 0.9494   (0.05060024846 -> 0.05060020096; 0.9493997515 -> 0.949399799)

  A structured diff of the two `REPORT.json` adds only the generation stamp
  (date, revision `fbc24f2`, the supervisor's records commit, which moved HEAD
  during this phase) and 22 wall-time `runtime_s` values. `make_report.py
  --check`: OK. Figures: `t5_delta_moment_bound.png` changed in content (the
  labels and legend); every other regenerated figure differed from the
  committed one only in its footer stamp ("LiPolGen benchmark <date> rev
  <hash>", measured pixel by pixel: the differences confined to the footer's
  rows), so those 18 were restored from `HEAD` and not re-committed.
  `validation/price_defaults.py --plot`, run once (164 s, N = 200 000, six
  families): the default row's `A2_S1…S4` are exactly the halved `A2_TOY`
  above, σ and A_zz unchanged; `price_defaults.png` (ratios) differs from the
  committed one only in its footer stamp and is kept regenerated, as the task
  asked. The optional U-5 edit of `BENCHMARK_PLAN.md` was withdrawn because
  `make_report.py` renders it into `REPORT.json`'s `planned_not_wired`.

**R6 on the way.** The first complete edit set left 14 broken R6 citations:
one in a 2026-09-02 record (`test_module.py:265`, moved by a two-line
insertion) and 13 in this run's own design and literature records, which cite
the t5 harness, its pytest, the dump script, `t1_r_world.py` and
`t5_weighted_unweighted.py` at `ca9b23a` lines (the run's records-hygiene
worker had removed their `(as of ca9b23a)` marks). Since those records are not
this worker's to edit, every one was answered from this side: the insertion
in `test_module.py` withdrawn, and every edit in the five files made
line-neutral up to the last cited line (the D4 note of the t5 harness and the
dump preflight moved below the cited code). `--records-only` is back at
0 broken, and `test_doc_link_gate.py`'s three records tests pass. This
record's own first draft then tripped two census tests of that file (one
citation over the 200 the records may skip, and three `(as of …)` anchors
beyond the pinned 110); its citations were rewritten, not the census.
R6 does not see every move: six record citations into files this phase
re-numbered still share a token with what now sits at their lines (§10,
V-5).

## 6. The check's items: applied, and not

Applied — every missing site, every "wrong" statement and every underrated
risk, as follows. Missing sites: the pytest's line 435 re-taken and `pre_d4_over_jm` pinned;
the t5 header's MODELS coefficient; `test_tensor_gamma.cpp` :307;
`REFERENCES.md` :424; `constants.hpp` :26-28 (and its siblings `CONVENTIONS.md`
:82-83 and `PHYSICS_CHANNELS.md` :123); `needs_survey.md` 52, 59, 70, 198 (with
443, 452); `sf.hpp` :529 (and :428's same error); `DEVELOPMENT_PLAN.md` :28;
SPIN32 :986 and the five R6 breaks avoided by binding at the module's end; the
Δ row's symbols cell. Wrong statements: P1 placed line-neutrally; T4 out of
the step-1 failure list; SS90 (A7) quoted with ⅙; the dump script's guard
stated truthfully and a refusing preflight added; `spin_weight_k` not claimed;
the b-basis reading documented, nothing subtracted; the parton side stated as
open; the halved significances per bin (21.5 / 15.3 / 17.3 / 10.7σ toy, 13.9 /
9.2 / 13.6 / 9.7σ grid; ranges 10.7–21.5σ and 9.2–13.9σ) in
`08_channel_priority.md`; N1's J = 3/2 projections ±3/2, ±1/2;
Pascalutsa–Vanderhaeghen kept "inferred". Risks: provenance made mandatory
and pinned (N2); the TIGHT edge kept as declared and relabelled; SPIN32 got no
new line; line shifts avoided (§5, "R6 on the way"); the one-line rule
written down; D7's order recorded (§8).

Not applied, with the reason:
- the design's "TIGHT loses V6" and its Y2/Y3 `saturation_tight` re-takes —
  superseded by amendment 5 (the author has not answered plan §12 item 3);
- L7 (an "F[U T_TT] = −Δ as γ → 0" sentence in `include/lipolgen/xsec.hpp`'s
  caveat (ii)) — optional, and the identification is inferred and holds with
  this tree's sign only; it is in `docs/CONVENTIONS.md`, marked inferred;
- P4 in `test_module.py` — moved to N2 (a 2026-09-02 record cites that file's
  line 265);
- the plan's §14 corrections to the D4 entry of `DECISIONS.md` — the author's
  decision record, which this task did not name; proposed for the supervisor
  in §8 item 8, with item 5 reworded (the b-basis reading needs no filing);
- the vendored `delta_moment_bound_sources.json` LOI18 note (V3) — a dated,
  sha256-pinned source record, not re-hashed; the t5 D4 note says it describes
  the pre-D4 kernel;
- a re-dump of `validation/reference/xsec.json` — never from this tree; after
  the sibling changes (§7);
- O8, the optional U-5 note in `BENCHMARK_PLAN.md` ("the same Δ
  normalisation; the distance is the hadron") — it surfaced in `REPORT.json`'s
  `planned_not_wired`, outside row 13, so it was withdrawn;
- the plan's §10 item 4 for a TILTED plan (x, Q², m, cell identical, φ
  different, between a pre-D4 and a D4 file) — no pre-D4 build exists in this
  phase; the exact-half contract (§2 item 1) and the default-run identity
  (§2 item 6) are what was measured.

## 7. Sibling sites for the author (PolarizedLithiumSim, read-only here)

From plan §13, with the check's corrections (the check re-confirmed the line
numbers at the sibling's `4726812`):

| file:line (sibling) | what | change for Jaffe–Manohar |
|---|---|---|
| `fastsim/polli_fastsim/asymmetries.py` lines 4–9 | master-formula docstring: "Hoodbhoy-Jaffe-Manohar NPB 312:571", Δ line without ½ | ½; HJM for the b-sector, JM PLB 223 Eq. (16) for Δ |
| `fastsim/polli_fastsim/asymmetries.py` lines 235–237 | `a_cos2phi` | ½ |
| `evgen/polligen/xsec.py` lines 17, 19–21 | docstring a₂; "Spin-1 (Hoodbhoy-Jaffe-Manohar NPB 312:571 …)" | ½; attribution |
| `evgen/polligen/xsec.py` lines 732–733 | the a₂ line | ½ |
| `evgen/scripts/tensor_gamma_leakage.py` lines 14, 34, 43, 77, 202 | Δ_fake = −a₂ y² D_φ/(1−y); 0.14–0.16; 0.23–0.26; leakage fractions of the amplitude | −2a₂…; 0.28–0.32; 0.46–0.52; fractions ×2 |
| `evgen/tests/test_xsec_identity.py` lines 306, 328 | −(ε/2) closed form | −(ε/4) |
| `evgen/tests/test_recopseudo.py` lines 561, 576 | "HJM NPB 312:571 (1989) Eq. (30)"; −(1−y)/(y y) | attribution; ½ |
| `evgen/tests/test_tensor_gamma.py` | the leakage pins T4 ports | ×2 (JM units) |
| `fastsim/polli_fastsim/fom.py` line 286, `fastsim/scripts/money_delta*.py` (e.g. `money_delta.py` line 89), `evgen/scripts/money_cos2phi.py` | consumers of `a_cos2phi` / the kernel a₂ | amplitudes halve for the same Δ; significances halve |
| `paper/main.tex` lines 199–210 | Eq. (1), Δ term without ½, cited to `Jaffe:1989xy,Hoodbhoy:1988am` | ½; the Δ term is Jaffe–Manohar's |
| `paper/main.tex` lines 232–237 | Eq. (3): A = −(1−y)/y² Δ/D_φ, Δ̂ = −Â y² D_φ/(1−y) | −½…; Δ̂ = −2Â y² D_φ/(1−y) |
| `paper/main.tex` line 116; `paper/table1.tex`; `paper/figs/fig2_caption.tex` | 21–43σ; injected amplitudes; Â = 7.31 × 10⁻³, 8.88 × 10⁻³ | halve for the same Δ (21.5 / 15.3 / 17.3 / 10.7σ toy, 13.9 / 9.2 / 13.6 / 9.7σ grid) — or Δ doubles for the same amplitude: which the letter means is the author's |
| `plans/07_plb_letter_gluonometry.md` lines 52, 105, 108 | sweet-spot amplitudes (0.44–0.95) × 10⁻²; 21–43σ / 32–64σ | as above |
| `docs/Discussions.pptx` p. 5 | LOI Eq. (9) transcribed and later called HJM's | annotate |
| `docs/code_review_2026-08-25.md` §3 | normalisation left open "for want of Jaffe–Manohar 1989" | settled |

Never apply the factor twice: once the sibling's coefficient is halved, its Δ
inputs (`C_BAG`, 𝒩 = −0.292) stay as they are. After the update the
`xsec.json` re-dump is `python3 validation/dump_polligen_reference.py` (that
script has no `--accept-changed-numbers`; its preflight prints the coefficient
it finds and refuses an unknown or a half-updated one; it rewrites all seven
tables, so read the diff). The port gates then read c_ref = ½ and compare
plainly with no code change, and t5's V6 label "the letter as printed" can be
retired by a dated edit.

## 8. Open points

1. **The sign of Δ for a given gluon linear polarization** — JM (2) +, DS16
   (13) − ("differs by a sign from that in Ref. [7]"), SS90 (2.11) −2, KS20
   (68). Not the coefficient; it decides whether the letter's `moment_A` sign
   (opposite the shipped toy's) is SS90's. Settles it: a box-graph calculation
   in one fixed state convention (R1 §6 item 3).
2. **The exact Jaffe–Manohar Eq. (14) against the Bjorken-limit coefficient** —
   the kernel carries the Bjorken limit; the exact contraction differs by
   −0.7 % at x = 0.1, y = 0.5, Q² = 5 GeV² (0.4964 of the pre-D4 coefficient)
   and by 0.04–0.6 % at the four sweet spots (0.49981 / 0.49879 / 0.49718 /
   0.49756). A modelling choice, not a convention (R1 §6 item 6).
3. **The spin-3/2 slot's normalisation** — JM state that their analysis
   applies without modification to spin 3/2, so the ½ applies to
   `delta_32_func` through the same line (N1 pins it at m = ±3/2, ±1/2); the
   c_eff = 3Q_NN geometry's normalisation relative to J = 1 is the open
   plans/04 #14 item, unchecked against JM's spin-J appendix (R1 §6 item 7). No
   shipped ⁷Li number moves (the slot is empty).
4. **The parton side, α_s/8π against α_s/16π** — the JM preprint's Eq. (29)
   prints α_s/(8π) for ∫x Δ dx where its own Eqs. (20), (21), (25) and SS90
   (3.1) give α_s/(16π) (on the printed Eq. (29) the bag moment would be
   −0.023 α_s); its Eq. (31) has a = g_x̂ − g_ŷ against SS90 (2.11)'s
   −2[g_x̂ − g_ŷ]. `C_BAG` is SS90's (α_s/16π). Settles it: the journal text or
   a fresh n = 2 matrix-element calculation (R1 §6 item 2).
5. **The t5 TIGHT edge** — kept with V6 as declared (amendment 5); whether it
   should lose V6 now that V6 is the pre-D4 reading is the author's question
   (plan §12 item 3).
6. **The coherent `amp`** (0.01 at P_zz = 1) — an amplitude, untouched by D4;
   if it was chosen to mirror the inclusive toy's amplitude it now mirrors
   twice it (plan §12 item 4).
7. **Order with D7.** D7 re-points the t5 harness's shipped kernel and
   re-takes `A2_TOY` and the t5 positivity pins; the D4 values here are on the
   ToyF2 default and must be re-measured on whichever kernel ships when D7
   lands (N2's four absolute amplitudes too: they are on
   `default_inclusive_kernel(li6())`).
8. **For the supervisor: the plan's §14 corrections to the D4 entry of
   `DECISIONS.md`**, not applied here: (1) the general formula
   a₂ = −½(1−y)/y² c_eff(m) sin²θ_S Δ/D_φ; (2) "read as its MIT preprint
   CTP#1712, Eq. (16)"; (3) "every **Δ-induced** cos 2φ amplitude" (θ_S = 0
   output unchanged; the leakage's Δ-equivalent doubles); (4) `moment_A/B` are
   the bag model at face value in this tree from D4 on, while the letter's
   printed amplitudes encode 2 × C_BAG until the sibling changes; (5) the
   attribution of the Δ term to HJM is corrected and the b-sector reading is
   documented (no double counting on the default path — nothing to file);
   (6) the port gate reads polligen's coefficient off `xsec.json`.
9. **Δ inside the finite-γ path** — a separate decision, not D4 (plan §12
   item 9).
10. **A HepMC3 record of the Δ convention** (§10, V-2) — HepMC3 output
    (`--hepmc`, `Pipeline.write_hepmc`, the C++ `HepMC3Writer`, every example
    executable) carries none, so a HepMC3 file of a tilted plan made after D4
    cannot be told from one made before. A run-level `GenRunInfo` attribute
    would do it, but `docs/HEPMC3_CONVENTION.md`, which would have to document
    it, is sha256-pinned by the `t6_hepmc3_convention` row: a follow-up with
    its own re-pin, not part of D4.
11. **Six drifted record citations** (§10, V-5) — for the supervisor or the
    records-hygiene worker to re-point or annotate; this worker may not edit
    those records.

## 9. Provenance

Session scratch `i4/`: `edits/step*.py` and `edits2/r_*.py` (every edit as
exact old/new strings, in order; `edits2/` the line-neutral rework of three
files, applied to their `ca9b23a` text) with `edits/apply.py`; `recon/` (the
reconstruction from `git archive`); `pre_REPORT.json` (sha256
`3ce51a5f…b780afa`, the committed REPORT.json); `cpp_step1.txt`,
`cpp_step2.txt`, `cpp_final.txt`; `pytest_full_step1.txt`, `pytest_full_final*.txt`,
`pytest_bench_final.txt`; `pytest_targeted_*.txt`; `gate_*.txt`;
`d4_contract_check.txt`; `measure_gamma2_post.out`; `report_compare2.txt`,
`report_check.txt`, `report_json_diff.txt`; `price_defaults_d4.json`; `run5k_d4.npz`.
Nothing was published or sent; the sibling was read only through `python3 -B`
imports and shell reads.

## 10. Verification

The adversarial verification of the tree I4 left (session scratch `v4/`).
Its verdict: "SOUND, WITH FIXES NEEDED". The D4 physics and numbers are right,
the factor is applied once at the two sites, and every gate tally and suite
count in this record reproduces. It found two medium defects and six low ones.
Worker F4 fixed them on 2026-10-27 (session scratch `f4/`). Each finding was
re-checked against the tree before anything was touched, and none was refuted.
F4 is line-neutral in every file a record or a gated document cites by line,
with two exceptions that nothing cites by line: the end of
`docs/CONVENTIONS.md`, and N2's docstring, which grew by one line. The
verifier labelled no finding with a lens; the lens column is F4's assignment,
using the groups of the verifier's list of what it found clean.

| # | lens | finding | F4's re-check | what was done |
|---|---|---|---|---|
| V-1 (medium) | code | The Python port gate pinned today's c_ref = 1 in three places, so it would fail a correctly re-dumped Jaffe–Manohar `xsec.json` that the C++ gate accepts. The three: `test_amplitudes_every_state` asserted the factor was ½, the reader test asserted c_ref == 1, and the reader's Jaffe–Manohar copy halved the file once more. This contradicts plan §5 item 3 and this record's "with no code change". | Confirmed on F4's own transformed copies of `xsec.json`, bit-identical to the verifier's. Each copy rescales the 1206 Δ-carrying leaves: `a_cos2phi`, the a₂ Δ line, and the θ_S ≠ 0 dsigma rows. On the Jaffe–Manohar re-dump the Python module failed 4 tests; a scratch build of the C++ gate passed 7033 of 7033. | Fixed in place, line-neutral. The `assert` at `python/tests/test_kernel_reference.py:159` now accepts a factor of ½ or 1, and the reader test accepts c_ref = 1 or ½ in a Δ-carrying block. Its Jaffe–Manohar copy is now scaled by `_jm_factor` (0.5/c_file) (`python/tests/test_kernel_reference.py:241-268`). Measured after the fix: the real file and the re-dump each pass all 9 of the module's `xsec.json` tests. The half-updated copies (`a_cos2phi` only, or the kernel's Δ line only) and the ×0.7 copy are refused, with 4, 4 and 6 tests failing. Against the verifier's scratch build at 1.0, all six Δ-carrying tests fail on both files, so the literal 0.5 still anchors the gate. |
| V-2 (medium) | honesty | The documents said every generated file carries `meta["delta_convention"]` and that a file without the key predates D4. HepMC3 output carries no record of the convention. | Confirmed with a 50-event `transverse-tensor` CLI run writing `--npz` and `--hepmc` (scratch `f4/out/`). The npz `meta` and the banner carry the string. In the HepMC3 file, the only "convention" token is the writer's tool string. Two more cases lack the key, both measured on 50 events: an .npz that `write_columns_npz` writes from `columns_from_events`, and one that `write_hfs_npz` writes without that `meta`. | Reworded in place at all six sites: the Provenance part of the `docs/CONVENTIONS.md` bullet, `docs/USAGE.md` §2, the `PHYSICS_CHANNELS.md` Δ row, the label's doc in `constants.hpp`, N2's docstring, and §0 here. The key is in the `meta` of `Pipeline.generate()`, and so on every .npz the CLI writes (`--npz`, and `--hfs-npz`, which is handed that `meta`). A HepMC3 file carries no run-level record, so a missing key dates only such an .npz. A HepMC3 attribute is filed as §8 item 10. |
| V-3 (low) | code | `DELTA_CONVENTION` was a free-standing string, so "setting it back is the whole of the reversal" would leave the label naming the ½ at 1.0. | Confirmed by reading. The verifier ran it: at 1.0 the label, the npz `meta` and the banner all still read "-(1/2)(1-y)/y^2", and N2's two provenance tests pass. | Added `static_assert(DELTA_COS2PHI_NORM == 0.5, …)` beside the label (`include/lipolgen/constants.hpp:164`), line-neutral: the label's doc was reflowed to three lines. "The whole of the reversal" now names both, in `constants.hpp` and `docs/CONVENTIONS.md`. Checked in scratch: a copy of the header with the constant at 1.0 fails to compile with "DELTA_CONVENTION names the 1/2", and the tree's header compiles. |
| V-4 (low) | honesty | The C++ port gate's a₂ and dsigma comparisons gained absolute floors, kRtol·(\|a2_ref\| + \|dline\|) and kRtol·\|ds_ref\|. Against the pre-D4 `kRtol, 1e-300` that is 5× looser on the massless a₂ rows (the tolerance there is 2.5 kRtol·\|Δ line\| against 0.5 kRtol·\|Δ line\|) and about 2× on dsigma; the verifier gave 2.5–5× on the Δ-carrying rows. It was not needed and not recorded. | Confirmed with a scratch build of the gate with `kRtol, 1e-300` restored at both sites. The real file and the re-dump pass 7033 of 7033. Of the copies, the `a_cos2phi`-only and kernel-only ones fail 1176 assertions each, and the ×0.7 one fails 5. Linked against the verifier's library at 1.0, the gate fails 1206 assertions on both files, every Δ-carrying field. | Restored at both sites, line-neutral, each marked "the pre-D4 tolerance, unchanged": the a₂ check at `tests/test_reference.cpp:430-431` and the dsigma check at `tests/test_reference.cpp:454-457`. The Python gate was not loosened and stays as it is. Its old `pytest.approx` carried an absolute floor of 1e-12; the D4 form is rtol 1e-12 with a floor of 1e-12 × the Δ line, which is tighter. |
| V-5 (low) | line-neutrality | Three record citations into files I4 did not keep line-neutral now land on other code. R6 passes them on one shared token, and this record did not name them. | Confirmed. F4's own scan maps every `path:line` citation under `docs/` into a file changed since `ca9b23a` through difflib's equal blocks. It finds the three, plus three that moved only in part: the start is right and the end falls 4 lines short. None is in a record this worker may edit. The scan's list is the same before and after F4. | Named in the table below for the supervisor or the records-hygiene worker (§8 item 11). |
| V-6 (low) | honesty | §0 said every Δ-carrying a₂ is exactly half its pre-D4 value. On the `tensor_gamma` path a₂ is h₂ plus the Δ line, and only the Δ line halves, as §2 item 1 and the contract checker say. | Confirmed by reading. That path adds h₂ to a₂ (`src/core/xsec.cpp:262-266`, `out.a2 = out.a2 + h.h2`), while on the massless path a₂ is the Δ line alone. | §0 reworded. |
| V-7 (low) | docs | Four code comments gave "Jaffe-Manohar PLB 223 (1989) 218 Eq. (16)" as the journal's equation number. Eq. (16) was read in the MIT preprint CTP#1712, and the journal text was never seen. | Confirmed by reading. `docs/CONVENTIONS.md`, the `[JM89g]` entry and 04_theory's T-16 carry the qualifier; the four comments did not. | Added "(preprint CTP#1712 numbering)" in `constants.hpp`, N1's comment in `tests/test_xsec.cpp` and N2's docstring. In `include/lipolgen/asymmetries.hpp:9-21` the change is to lines 10–11 ("preprint CTP#1712 Eq. (16)"), inside a fingerprinted block. Its one citing row, SPIN32's [HJM89] entry at line 1168, was re-read against the edited block and still says something true: the master formula, the b-sector's HJM attribution and the c_m line are all there. `--record-ranges` re-recorded that one fingerprint (0 refused). Line 9, the allow entry's pinned first line, is untouched. |
| V-8 (low) | docs | The new [LOI18] entry gave the authors as "J. Maxwell, D. Crabb, D. Day et al.". | Confirmed on the LOI's text (session scratch `r1a/loi18.txt`, lines 4 and 38). The printed list begins "M. Jones, C. Keith, J. Maxwell∗", where ∗ marks the corresponding author. | Corrected in place to "M. Jones, C. Keith, J. Maxwell (corresponding author) et al.". The older `docs/references/01_tensor-b1.md` entry, which predates D4 and lists the names in another order, was left as it is. |

The six drifted citations (V-5). Every line number is written as plain text, so none of them is a gated citation:

| citing record and line | cites | lines at `ca9b23a` | where that content is now |
|---|---|---|---|
| `plan_D4_jm_convention.md`, line 205 (T5) | `tests/test_reference.cpp`, the port gate | 346–402 | 393–465 |
| `plan_D4_jm_convention.md`, line 224 (Y8) | `python/tests/test_kernel_reference.py`, the port gate | 98–174 | 153–284 |
| `plan_D7_backend_defaults.md`, line 1311 | `python/tests/test_kernel_reference.py`, the class-level kernels | 29–54 | 39–64 |
| `plan_D4_jm_convention.md`, line 225 (Y9) | `python/tests/test_price_defaults.py`, `A2_TOY` | 38–44 | 38–48 (partial: the cited range now holds only `A2_TOY`'s first line) |
| `plan_D7_backend_defaults.md`, lines 637 and 754 | `test_price_defaults.py`, the default row | 39–50 | 39–54 (partial: `RATIOS`' second line now falls outside) |

**Gates after F4** (the tree F4 leaves; uncommitted):

- `cmake --build build -j2`: clean, no warning, 2 min 14 s. The rebuild was
  near-full because `constants.hpp` changed.
- `./build/lipolgen_tests`: 417 cases, 417 passed, 0 skipped, 17 244 191 assertions, in 7 min 30 s. These are I4's counts: restoring the tolerance moved no assertion.
- `python3 -m pytest -p no:cacheprovider python/tests -q`: 1463 passed, 155 skipped, 0 failed, 469 s (I4's counts).
  `LIPOLGEN_BENCH=1 … -k bench`: 351 passed, 1 skipped, 0 failed, 305 s. The skip is the nnpdf-data wheel, which is not on this machine (environment, as before).
- Targeted files with `LIPOLGEN_BENCH=1`: `test_b1_model` 51, `test_bench_t1_r_world` 51, `test_bench_t5_delta_bound` 23, `test_kernel_reference` 11, `test_meta_provenance` 8, `test_module` 89, `test_price_defaults` 1, `test_delta_convention` 6 and `test_doc_link_gate` 112. All passed with no skip, I4's counts.
- The strict citation gate: 1241 refs / 95 ranges / 7 external / 0 broken / 24 allow-listed (as I4); SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0 / 0. SPIN32 was at 1 broken, the
  EDITED `asymmetries.hpp` block (V-7), until it was re-read and
  re-recorded. `--records-only`: 0 broken. There are 1423 citations in 75 records, this one included: 924 checked against a live file, 110 annotated historical and 199 skipped. So the seven citations §10 adds are all checked live, and the census of skips stays under its 200. `--audit-ranges`:
  209 range citations on 179 blocks, 0 REFUSED. `check_spdx_headers.py`: 133 files.
- `make_report.py --check`: OK, regenerated in 1 min 3 s. No harness, `REPORT.json` or
  figure changed in F4.
- The design's exact-half checker (`d4/d4_before_after.py check`) was re-run
  against the pre-D4 dump (sha256 `0334237e…`) on the module F4 rebuilt. It
  reports: "D4 contract holds: 5310 states exact half, 1593 tensor_gamma
  states half on the Delta line, every Delta-free quantity bit-identical".
  So F4 moved no number.

Scratch `f4/` holds the following. The step scripts `edits/step_a_code.py`,
`step_b_docs.py` and `step_c_record.py` (every edit as an exact old/new pair)
and `tkr_edits.py`. `pre_f4/`, the edited files as I4 left them.
`make_variants.py` and `variants/` (the transformed `xsec.json` copies).
`run_pygate.py`, `build_ref.sh`, `run_cpp_variant.sh` and the `reftest_*`
binaries. `drift_scan.py`. `out/` (the provenance sample). `sa/` (the
`static_assert` check). `build1.log`, `cpp_full.txt`, `pytest_*.txt`,
`gate_*.txt`, `report_check.txt` and `d4_contract_after_f4.txt`. Nothing was committed, published or sent, and the sibling was not
opened.
