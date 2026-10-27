<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-27 — plan D4: the Jaffe–Manohar normalisation of Δ in the kernel's cos 2φ term (design only)

Worker D4, 2026-10-27, on `claude/run-2026-10-26` at `ca9b23a` (the build in
`build/` is that commit's: `_lipolgen*.so` 2026-10-26 13:04, every C++ source
older, no C++ commit since 2026-09-27). **Design only: nothing in the tree was
edited except this file.** Decision: `DECISIONS.md` D4, the author's words
"Then we go with the Jaffe-Manohar convention (your recommendation)". Inputs
read: `DECISIONS.md`, `STATUS.md`, `../run_2026-10-26/phase_C2_delta_bound.md`
§§1–10 (the first derivation of the factor 2, §4.2; the priced options, §9),
`../run_2026-10-26/DECISIONS_PRICED.md` item 4, the readers' text files of
phase R1 (session scratch `r1a/`, `r1b/`; `phase_R1_delta_convention.md` was
not on disk when this record was written, so the readers' findings are taken
from the task brief and re-checked where this record says so), every site
listed in §3, and the sibling `PolarizedLithiumSim` read-only (shell reads;
its HEAD is `4726812`, the commit the C2 record names).

Conventions of this record: **read** = I read the text or the page; **ran** =
I measured it on the `ca9b23a` build through the Python API; **inferred** =
reasoning, flagged as such. Line citations are to the tree at `ca9b23a`. The
citations of lines D4 itself will change first carried `(as of ca9b23a)`; that
annotation was taken out on 2026-10-27, because the records gate's as-of census
(`test_every_as_of_anchor_passes_the_bar_at_its_commit`) admits no new
historical anchor, so every citation here is now checked against the live
tree, and one whose line the implementation moves fails that gate until it is
re-pointed (or annotated, which needs that census raised).

## 0. The change in one paragraph

The kernel's double-helicity-flip term is today
a₂ += −(1−y)/y² · c_eff(m) sin²θ_S · Δ/D_φ (`src/core/xsec.cpp:280-283`, the `with_delta` branch),
and the standalone identity `a_cos2phi` spells the same coefficient
(`src/core/asymmetries.cpp:90-92`). These are the only two places the
library evaluates the coefficient (a grep for the spelled factor finds exactly
these two lines; the sampler, the RC band, the Mode-W weights and every
harness read the kernel's `a2`). D4 multiplies both by ½, so that
**a₂ = −½ (1−y)/y² · c_eff(m) sin²θ_S · Δ/D_φ** — for a pure m = ±1 state with
its axis transverse to the beam a₂ = −½(1−y)/y² Δ/D_φ — which is Jaffe–Manohar's
Eq. (16), Sather–Schmidt's (A7), Detmold–Shanahan's Eqs. (3)–(5) contracted,
and the LOI's own Eq. (6). Δ itself does not move: every Δ model, `C_BAG` and
the shipped `toy_delta_gluon(…, 1e−2)` keep their values, so **every
Δ-induced cos 2φ amplitude the generator produces halves**, and only those
(fills with sin θ_S ≠ 0 — `transverse-tensor`, `tensor-flip`, a tilted
`helicity-flip` — on a spin-1 inclusive kernel; the default `tensor-thirds`
plan sits at θ_S = 0 and its output does not change). The b-sector, A_zz, the
finite-γ path's b₁…b₄ pieces, every R or F₂ *ratio* of the cos 2φ amplitude,
and the coherent channel's c₂ are untouched. Recommended form: one named
constant, `DELTA_COS2PHI_NORM = 0.5` in `constants.hpp` (the `TENSOR_LL_SIGN`
precedent — one copy, value pinned, "setting it back is the whole of the
reversal"), applied as a factor at the two sites so that the new values are
**exactly** half the old ones bit for bit (multiplying by 2⁻¹ commutes with
IEEE rounding); no restructuring of either expression.

## 1. The evidence, re-checked here

| source | what it prints for a pure m = ±1 state, axis ⊥ q | coefficient of Δ relative to D_φ = F₁ + (1−y)F₂/(xy²) | where read (session scratch) | status |
|---|---|---|---|---|
| Jaffe–Manohar, *Nuclear gluonometry*, PLB 223 (1989) 218, as the MIT preprint CTP#1712 (March 1989) | Eq. (16): lim dσ/dx dy dφ = (e⁴ME/4π²Q⁴){xy²F₁ + (1−y)F₂ − (x(1−y)/2) Δ cos 2φ}, target E_x̂ (m_x = +1) | **−½ (1−y)/y²** | `r1a/jm89.txt` lines 170–188 (OCR of the preprint) | read |
| Sather–Schmidt, PRD 42 (1990) 1424 | App. (A7), m_J = 1 along x: […+ ½(xy²Δ₁ + λΔ₂) − ½ xλΔ₄ cos 2φ …], λ → 1−y; "Δ = Δ₄ in the Bjorken limit" | **−½ (1−y)/y²** | `r1a/ss90.txt` lines 345, 387–392; p. 1428 page image | read |
| Detmold–Shanahan, PRD 94 (2016) 014507 | Eqs. (3)–(5), Δ = A₊₋,₋₊ with F₁ = ⅓Σ_H A₊H,₊H, contracted with the lepton tensor | ratio to the tree's coefficient 0.4999998 / 0.4999984 / 0.4999999 at Q² = 10⁵ GeV² (x, y) = (0.1, 0.5), (0.3, 0.3), (0.01, 0.7) | the t5 pytest's own `_ds16_cos2phi_over_f1`, re-run (`d4/ds16_ratio.out`) | ran |
| JLab LOI, arXiv:1803.11206 | Eq. (6): the JM bracket with −x(1−y)/2 Δ cos 2φ; Eq. (9): −(1−y)/y² b_m sin²θ_m Δ cos 2φ, underived | (6) **−½**; (9) **−1** (twice) | `r1a/loi18.txt` lines 136–143 and 186–202 | read |
| Pascalutsa–Vanderhaeghen, arXiv:2610.09862 | Eq. (26) ε T_TT cos 2φ_TT F_UTTT; (27b) F_UTTT = −(M/3ν)(b₂ − 3b₃); (29b) T_TT cos 2φ_TT = (P_zz/2) sin²θ cos 2φ; (A8a) Δ = (1/ν)(⅓b₂ − b₃) at M = 1 | F_UTTT = −Δ with T_TT = ½ sin²θ ⇒ **−½** | `r1b/arxiv_2610.09862.txt` lines 551, 561, 578, 850 | read |
| the sibling's `docs/Discussions.pptx` p. 5 — what `include/lipolgen/constants.hpp:26-28` calls "the transcription of Hoodbhoy-Jaffe-Manohar" | the LOI's Eq. (9), verbatim | **−1** | `r1a/pptx_slide5/image8.png` | read |
| Hoodbhoy–Jaffe–Manohar, NPB 312 (1989) 571 | no Δ term: b₃, b₄ "do not receive contributions at leading twist … appear at twist-4"; "the double helicity flip amplitude A₊₋,₋₊ falls off like 1/Q²" | — | `r1a/hjm89.txt` lines 90–91, 209 | read |

Also read: SS90 p. 1424, "Later it was realized that there is another tower of
gluonic operators" (`r1a/ss90.txt` line 43); SS90 footnote 8 on the p. 1428
page image, "In the Bjorken limit they are related by b₁ = Δ₁ + ½Δ₄ and
b₂ = Δ₂ + xΔ₄" (§12 item 1); JM's statement that the analysis "applies without
modification to a spin-3/2 target" (`r1a/jm89.txt` lines 93–94). **Not
verified by me:** the journal text of PLB 223 (only the preprint), the 2019
proceedings reprint of LOI Eq. (9), the Kumano–Song / Kumano 2024 parton-side
agreement (reader statements, taken as given), and Cosyn *et al.* 2025 as a Δ
source — its F[U T_TT] is built from b₂, b₃ only (Eq. 17e) and it has no Δ;
the "F_UTTT = −Δ" statement is Pascalutsa–Vanderhaeghen's.

**A correction to the brief's reader-B inference.** The finite-γ path
(`tensor_harmonics_gamma`, `src/core/xsec.cpp:182-219`) carries **no** Δ: its
F[U T_TT] is `cosyn_tensor_sfs`' f_tt = −(γ²/x)(b₂/6 − b₃/2)
(`src/core/xsec.cpp:20-37`), and the Δ line is added outside it on both paths
(`include/lipolgen/xsec.hpp:300-304`, "keeps its massless kinematic factor in
both paths"). What is true (ran, §2.3) is that the Δ line is *equivalent*, as
γ → 0, to an F[U T_TT] insertion of **−2 Δ_tree** in that path's own
normalisation — and to −Δ after D4, Pascalutsa–Vanderhaeghen's value.

## 2. Measured on the `ca9b23a` build (Python API, no build)

### 2.1 a₂ at the four sweet-spot cells (pin material)

`default_inclusive_kernel(li6())`, m = +1, θ_S = 90°, unpolarized electrons,
y = Q²/(s x) at `default_configs('6Li')[1]` (s = 3980 GeV², e 10 × ⁶Li
99.5 GeV/u); the kernel's a₂ equals `lg.a_cos2phi` there bit for bit (difference
0.0), and a₂(m = 0) = −2 a₂(m = ±1) exactly (`d4/measure_a2.out`):

| cell | x | Q² | y | Δ/F₁ | a₂ today | **a₂ after D4 (exactly ½)** |
|---|---|---|---|---|---|---|
| spot 2 | 0.0112 | 1.14 | 0.0255743 | 2.484273e−3 | −1.055948214152054e−03 | **−5.279741070760268e−04** |
| spot 1 | 0.0282 | 1.14 | 0.01015717 | 3.057612e−3 | −1.299961141425631e−03 | **−6.499805707128153e−04** |
| spot 3 | 0.0708 | 3.13 | 0.0111078 | 3.368579e−3 | −1.440230754221272e−03 | **−7.201153771106358e−04** |
| spot 4 | 0.141 | 14.3 | 0.02548202 | 3.025091e−3 | −1.326442680639959e−03 | **−6.632213403199796e−04** |

The recipe's super-bin amplitudes per unit P_zz (`python/tests/test_price_defaults.py:40-43`)
halve exactly too: −6.514569020677085e−04 / −5.312811384636748e−04 /
−7.145394171261867e−04 / −6.645438727281629e−04 for S1 … S4 (a sum of halved
terms is exactly half the sum).

### 2.2 The t5 row's helicity-state positivity sub-row

On the harness's own grid and function (567 points, 1701 rates): today
**min W = 0.9493997515410681** at x = 0.9, Q² = 2, y = 0.95, θ_S = 90°,
φ′ = 90°, m = 0; A_T ∈ [−0.10120030690436055, 0.05060024845893199]. With the Δ
term halved (on the default path a₂ *is* that term): **min W =
0.9493997990444439** (same point), A_T ∈ [−0.10120030690436055 (θ_S = 0,
unmoved), **0.05060020095555617**] (`d4/measure_pos.out`). PASS either way; the
margin is b₁'s, as C2 §6 says.

### 2.3 The finite-γ path against the Δ line

A pure F[U T_TT] insertion read off the kernel itself — b₁…b₄ solved through
`cosyn_tensor_sfs` so that F_TLL_T = F_TLL_L = F_TLT = 0, F_TTT = 1, then
`tensor_harmonics_gamma` at m = +1, θ_S = 90° — against the Δ line per unit Δ
on the same table (`d4/measure_gamma2.out`): the Δ line equals an F_TTT of
**−2.000784754** Δ at (x, y, Q²) = (0.01, 0.5, 1.14), −2.000008946 at Q² = 100,
**−2.000000001** at Q² = 10⁶ (γ² = 3.5 × 10⁻¹⁰); −2.0073478 at the spot-1
kinematics (x = 0.0282, y = 0.01016, Q² = 1.14). So today the Δ line is
F_UTTT = −2Δ_tree as γ → 0 and after D4 it is −Δ; the O(γ²) remainder
(0.04 % at x = 0.01, 4 % at x = 0.1, Q² = 1.14) is the "massless kinematic
factor" the header keeps on purpose. Also measured: the Δ line is
**bit-identical** with `tensor_gamma` on and off (`==`, 14 points), and
m = 0 carries −2.000000 times m = +1 on the γ path too (`d4/measure_gamma.out`).
**What `tensor_harmonics_gamma` must become: nothing.** Its b-pieces do not
contain Δ; on the γ path a₂ = h₂ + the Δ line, h₂ → 0 as γ → 0 (the b-sector
collapses identically onto the massless one, `tests/test_tensor_gamma.cpp:152`),
and the Δ line is the same line on both paths — so halving that line gives the
new a₂ on both paths, the γ → 0 identity survives, and h₀, h₁, h₂ stay
bit-identical (§10 checker). Moving Δ *into*
F[U T_TT] is not part of D4 (§12 item 9).

### 2.4 The two "oversized Δ" guard tests

Emulated on today's build (Δ = 1.5 F₁ under the old coefficient gives exactly
the a₂ of Δ = 3 F₁ under the new one; `d4/measure_oversized.out`):
the m = +1 margin over 0.01 ≤ y ≤ 0.95 goes from −0.4593 (old, Δ = 3F₁) to
**+0.2703** — so `tests/test_xsec.cpp:816-826`, which builds the
oversized kernel at m = +1 and asserts `negative`, **fails after D4**; at
c_eff = 1 the bound is |a₂| ≤ ¾ε/(1+εR) < 1 for every y. The sampler's
version (`tests/test_sampler.cpp:240-248`, `transverse_tensor_plan(0.6)`)
**still throws**, now through the m = 0 state (a₂/den = +1.449, margin
−0.449 at x = 0.6813, Q² = 242.1, "for m=0" instead of "for m=1").

### 2.5 The b₁ against Δ/2 (sizing open question 1)

Shipped ⁶Li kernel (`d4/measure_b1.out`): (Δ/2)/b₁ = **1.142 / 1.018 / 0.664 /
0.676** at the four cells, **9.84** at (0.3, 5), 0.046 at (0.6, 5), 3 × 10⁻⁶ at
(0.9, 2). The toy's Δ/2 is the size of the shipped b₁ wherever claim 1 lives.

## 3. Inventory — every site the ½ touches, and every site checked and left

Actions: **EDIT** (code), **RETAKE** (a pin, with a dated comment naming D4),
**REWORD** (text), **ADD**, **NONE** (checked; reason given). "In place" =
without changing the file's line count above the cited lines (§9 says why).

### 3.1 The library (C++)

| # | file:line | what | action |
|---|---|---|---|
| L1 | `src/core/xsec.cpp:280-283` | `tensor_amplitudes`' Δ line, `-(1.0 - y) / (y * y) * qc.second * st * st * t.delta / std::max(den, 1e-30)`, used by `amplitudes()` (always) and `tensor_amplitudes(..., with_delta)` (`RcScope::TensorAll`) | **EDIT in place**: line 282 becomes `* DELTA_COS2PHI_NORM / std::max(den, 1e-30));` (the factor anywhere in the product chain gives exactly half); the comment `:274-279` reworded within its six lines to name D4 |
| L2 | `src/core/asymmetries.cpp:90-92` | `a_cos2phi` — the standalone identity, bound as `lg.a_cos2phi` | **EDIT**: `return -(1.0 - y) / (y * y) * delta * DELTA_COS2PHI_NORM / phi_averaged_density(f1, f2, x, y);` on one line (97 columns; 156 lines of `src/core/*.cpp` already exceed 80) — or wrapped, then `--fix` (§9 step 5) |
| L3 | `include/lipolgen/constants.hpp:139-141` | `C_BAG = -0.012`, "Sather-Schmidt bag-model sum-rule coefficient" (line 139 carries an UNNAMED_ALLOW pin on exactly that phrase) | **NONE** to lines 139–141; **ADD** after line 141 the constant and its doc (§7 text, abridged): `inline constexpr double DELTA_COS2PHI_NORM = 0.5;` — shifts only `POLEMC_VALENCE_WINDOW_LO` (named, `--fix`) |
| L4 | `include/lipolgen/constants.hpp:31-40` | the `TENSOR_LL_SIGN` note: "the whole Delta (cos 2phi) sector [is] unaffected" by the sign | **NONE** (still true; the Δ sector depends on the new constant, not on the sign) |
| L5 | `include/lipolgen/xsec.hpp:13-32` | the master-formula header: `a_2 = -(1-y)/y^2 * c_eff(m) …` at line 21; "reproduces the Hoodbhoy-Jaffe-Manohar transcription digit for digit" at 28–29 (true of the c_m geometry) | **REWORD in place**: line 21 `a_2 = -(1/2)(1-y)/y^2 * c_eff(m) sin^2(theta_S) * Delta/D_phi [JM, D4]`; line 31–32 may add "(its normalisation is DELTA_COS2PHI_NORM)" only by rewording within the two lines; the range is fingerprinted (re-record, §9) |
| L6 | `include/lipolgen/xsec.hpp:209-216` | `tensor_gamma` doc: leakage "at most 0.109 % of the published amplitude" | **REWORD in place**: "0.109 % of the letter's `moment_A` amplitude as published (made with twice JM's coefficient); 0.22 % of the same Δ's amplitude since D4" — the leakage amplitude itself does not move |
| L7 | `include/lipolgen/xsec.hpp:296-304` | (ii) "Delta … keeps its massless kinematic factor in both paths" | **REWORD in place** (optional): "… and equals F[U T_TT] = −Δ (Pascalutsa–Vanderhaeghen Eq. (27b)) as γ → 0 since D4; it is not inserted there" |
| L8 | `include/lipolgen/asymmetries.hpp:9-27` | spin-1 master formula "(…; Hoodbhoy-Jaffe-Manohar NPB 312:571)" with the Δ line `- (1-y)/y^2 c_m sin^2(theta_m) Delta(x,Q2) cos(2 phi)`; `A2phi(…) = -(1-y)/y^2 Delta / D_phi` | **REWORD in place**: attribution (§8 text) and both coefficients to −½(1−y)/y²; range `include/lipolgen/asymmetries.hpp:9-21` is fingerprinted (re-record) |
| L9 | `include/lipolgen/asymmetries.hpp:123-124` | `a_cos2phi` doc "cos(2phi) amplitude for lambda_m = +-1, theta_m = 90 deg." | **REWORD in place**: "… = -DELTA_COS2PHI_NORM (1-y)/y^2 Delta/D_phi (Jaffe-Manohar, D4)" |
| L10 | `include/lipolgen/sf.hpp:428-430` | `toy_delta_gluon` doc: "scale is the peak Delta/F1. Use 1e-3, 3e-3, 1e-2 as discovery scenarios" | **REWORD in place** (line 429 only): "… (Delta in the kernel's Jaffe-Manohar normalisation since D4: the same scale gives half the pre-2026-10-27 amplitude)"; **the scale values do not move** |
| L11 | `include/lipolgen/sf.hpp:535-553`, `src/core/sf.cpp:472-503` | `make_moment_a` / `make_moment_b` with `c_moment = C_BAG` | **NONE** to code. Meaning: before D4 these fed SS90's −0.012 α_s into a kernel whose coefficient is twice JM's, so the amplitude they produced was that of Δ_JM = 2 C_BAG α_s (C2 §4.2, V6); after D4 the same functions produce SS90's bag model **at face value**. Optional one-line doc note in place |
| L12 | `src/core/pipeline.cpp:1262-1266` | the pipeline's Δ slot, `toy_delta_gluon(x, q2, f1, 1e-2)` | **NONE** (the scale keeps its value, D4) |
| L13 | `src/core/pipeline.cpp:1288-1292` | `rank2_input_report`: "delta_func = toy_delta_gluon at the 1e-2 discovery scale" → `meta["rank2_input"]` | **REWORD in place (recommended)**: "… at the 1e-2 discovery scale, Jaffe-Manohar normalisation of the cos 2phi term (D4)" — the one trace in an output file that says which convention its amplitude carries (no test pins the string) |
| L14 | `src/core/sampler.cpp:268-286` | `build_state`: `st.a2`, `a2_tensor`, `a2n`, `bound = 1 + |a1n| + |a2n|`, `density_min` refusal | **NONE** (reads the kernel). Consequence: the per-cell amplitudes halve; a user Δ that was refused as driving W negative may now be accepted (the refusal is on a₂) |
| L15 | `src/core/sampler.cpp:460-470` | φ by accept–reject on 1 + a1n cos φ′ + a2n cos 2φ′, last draw of the per-event stream | **NONE**. Consequence: on affected categories φ changes event by event; x, Q², m and cell (drawn earlier) do not; draws made after φ on the same stream change too — the struck-nucleon species draw `rng.uniform() < proton_fraction(…)` at `src/core/generator.cpp:144` |
| L16 | `src/core/sampler.cpp:539-612` | `mix_weights` → `weights_for` / `tensor_weights_for` (Mode-W, HepMC3 `spin_weight_k`) | **NONE**; the weights of tilted categories change |
| L17 | `src/core/rc.cpp:2141-2160`, `:2190-2205` | the RC band's τ = W_tensor/W uses `st.a2` (and `a2 − a2_tensor` under `TensorAll`) | **NONE**; `rc_tensor_lo/hi` on tilted fills move through W |
| L18 | `src/core/xsec.cpp:175-180`, `:182-219` | `tensor_moments` (Q_NN, 3Q_NN); `tensor_harmonics_gamma` (lam_tt = ε f_tt, tt harmonic → ½ sin²θ_S at γ → 0, scale −TENSOR_LL_SIGN/den) | **NONE** (§2.3) |
| L19 | `include/lipolgen/coherent.hpp:412-414`, `src/core/coherent.cpp:183-202` | coherent `amp` ("FLAT cos 2phi' modulation … at P_zz = 1", 0.01) and c₂ = 2 a₂(deformation) + amp·P_zz | **NONE**: an amplitude, not a Δ; nothing passes through the kernel's coefficient (open question 4) |
| L20 | `src/core/xsec.cpp:258-268` | `tensor_gamma_` branch adds `h.h2` | **NONE** |

### 3.2 Bindings and the Python package

| # | file:line | what | action |
|---|---|---|---|
| P1 | `python/bindings.cpp:1090-1091` | `m.attr("C_BAG") = C_BAG;` (the range's last line) | **ADD** after it `m.attr("DELTA_COS2PHI_NORM") = DELTA_COS2PHI_NORM;` (shifts 4 fingerprinted ranges and 3 named citations into the file — all `--fix`-able, §9) |
| P2 | `python/bindings.cpp:2102-2103` | `a_cos2phi` binding (no docstring) | **NONE** (passthrough); optional docstring |
| P3 | `python/bindings.cpp:4120-4140`, `python/lipolgen/__init__.py:224-236`, `python/lipolgen/cli.py:580-605`, `:200-212` | help text: the cos 2φ amplitude moves −8.33 / −6.73 / −3.14 % under `r1998`; −18.1195 % on δ(cos 2φ) | **NONE**: ratios, invariant under a common factor |
| P4 | `python/tests/test_module.py:63-80` | bound-name presence list | optional **ADD** `"DELTA_COS2PHI_NORM"` |

### 3.3 C++ tests

| # | file:line | what | action |
|---|---|---|---|
| T1 | `tests/test_xsec.cpp:441-461` | "the cos 2phi amplitude carries R and NOT F2": expected `-(eps_of(p.y) / 2.0) * d_over_f1 / (1.0 + …)` — the only C++ closed form of the coefficient | **RETAKE in place**: `/ 4.0`, comment "−(ε/4) c_eff sin²θ_S (Δ/F₁)/(1+εR), Jaffe–Manohar (D4, 2026-10-27; −(ε/2) before)" |
| T2 | `tests/test_xsec.cpp:805-827` | "production scenarios keep a healthy positivity margin" + "an oversized Delta drives the density negative" with Δ = 3·F₁ at m = +1 | **RETAKE in place**: `6.0 * f1` (reproduces the pre-D4 a₂ exactly; §2.4: 3·F₁ no longer goes negative at m = +1). The `> 0.9` production margin only grows |
| T3 | `tests/test_rc.cpp:205-210` | `delta_term = -(1.0 - y) / (y * y) * k.tensor_moments(m).second * …` written out | **RETAKE in place**: `-0.5 * (1.0 - y) / …` (the literal, not the constant: the test then pins the value) |
| T4 | `tests/test_tensor_gamma.cpp:338-378` | the leakage as an equivalent Δ: `d_fake = -h.h2 * p.y * p.y * dphi / (1.0 - p.y)`, want 0.1628248344088873 / 0.16006748610799001 / 0.15477697871375315, band `0.14 < coef < 0.17`, comment "the 1/6 of Eq. (17e)" | **RETAKE in place, same line count** (lines 152 and 382 carry UNNAMED_ALLOW pins; `docs/theory/SPIN32_FINITE_GAMMA.md` also names 406): `d_fake = -2.0 * h.h2 * …`, want **0.3256496688177746 / 0.32013497221598003 / 0.3095539574275063** (exactly ×2), band `0.28 < coef < 0.34`, comment "the 1/3 of Eq. (17e) — in JM units Δ_fake = −F[U T_TT] at leading order; the sibling's `tensor_gamma_leakage.py` prints 0.14–0.16, its pre-D4 units". `h2 < 0.0` unchanged |
| T5 | `tests/test_reference.cpp:346-402` | the port gate against `validation/reference/xsec.json` (polligen, LOI-Eq.-(9) coefficient): `a_cos2phi`, every state's `a2`, `dsigma_per_point` | **RETAKE** per §5 (helpers placed after `build_kernel`, never above line 96: `tests/test_reference.cpp:96-97`, the `kernel_scenario_rank2` `b1_32_func` / `delta_32_func` lines, is a cited range) |
| T6 | `tests/test_xsec.cpp:154-168`, `:280-293`, `:337-361`, `:463-512` | a₂ == `a_cos2phi`; the Δ sector vs the sign; the spin-3/2/spin-1 ratio; F₂-independence; the R1998 ratio | **NONE**: both sides move together or the ratio cancels — which is exactly why T1 and the new N1 are needed |
| T7 | `tests/test_sampler.cpp:238-260` | oversized Δ = 3·F₁ on `transverse_tensor_plan(0.6)` must throw | **NONE** (still throws via m = 0, §2.4); optional `6.0` for margin |
| T8 | `tests/test_sampler.cpp:503-519`, `:637-693` | closure truth via `a_cos2phi`; moment / binned-fit closures (fixed seeds) | **NONE** to code; outcomes are re-drawn (φ changes) — run them; errors are amplitude-independent, the "naive moment is biased" check is acceptance-driven |
| T9 | `tests/test_rc.cpp:3466-3501` | `TensorAll` buys back exactly `a2 − a2_tensor` | **NONE** (algebraic) |
| T10 | `tests/test_sf.cpp:367-425` (toy_delta_gluon and moment-model delta values), `tests/test_pipeline.cpp:1440-1445` and `:1559-1567` (t.delta equalities), `tests/test_b1_nuclear.cpp:1880-1890` and `:2375-2401` (delta tables, toy_delta_gluon), `tests/test_generator.cpp` (kinematics), `tests/test_rng.cpp:125-145` (batches identical at 1 and 8 threads), `tests/test_coherent.cpp` (c₂) | Δ *values*, determinism, kinematics, coherent c₂ | **NONE** |
| N1 | `tests/test_xsec.cpp` (append at the end) | — | **ADD** "D4: the cos 2phi coefficient is Jaffe-Manohar's": `CHECK(DELTA_COS2PHI_NORM == 0.5)`; the closed form −½(1−y)/y²·c_eff·sin²θ_S·Δ/D_φ written out with no reference to the constant, at m ∈ {+1, 0, −1}, θ_S ∈ {0.7, π/2}, for LI6 (`delta_func`) and LI7 (`delta_32_func`), on the analysis grid, kRtol; and `a_cos2phi(1, 1, 0, x, y) == -0.5 (1-y)/y^2` |

### 3.4 Python tests

| # | file:line | what | action |
|---|---|---|---|
| Y1 | `python/tests/test_bench_t5_delta_bound.py:283-300` | S2: `tree == approx(-(1.0 - y) / y ** 2)`; `ratio == approx(1.0 / mod.DS16_OVER_TREE, abs=2e-5)`; Q² = 5 check `approx(0.5, rel=0.02)` | **RETAKE**: `-0.5 * (1.0 - y) / y ** 2`; `DS16_OVER_TREE = 1.0` makes the second line hold at 3.3e−6 (§1 table); the Q² = 5 check becomes `approx(1.0, rel=0.02)` (0.99214); rename to `test_s2_the_kernel_coefficient_is_ds16_eq5` |
| Y2 | `python/tests/test_bench_t5_delta_bound.py:339-363` | headline pins | **RETAKE** `saturation_tight` only: **0.13596137455 / 0.16845554105 / 0.1968291457** (×½, V6 leaves the tight edge); `m2`, `bound`, `saturation`, `saturation_v6`, `f1ll_critical(_v6)`, statuses **unchanged** |
| Y3 | `python/tests/test_bench_t5_delta_bound.py:447-494` | letter rows | **RETAKE** `moment_A[letter]` `saturation_tight` → **1.5077958005 / 1.466479881 / 1.469045308**; docstring "~140 % under V6 (Δ_DS16 = 2Δ_tree)" → "~140 % for the letter's printed amplitudes (pre-D4 coefficient)"; everything else unchanged |
| Y4 | `python/tests/test_bench_t5_delta_bound.py:496-504` | positivity pins | **RETAKE** `w_min` → **0.9493997990** (rel 1e−9), `a_t_max` → **0.05060020096** (rel 1e−8); `a_t_min`, `w_min_at`, counts unchanged |
| Y5 | `python/tests/test_bench_t5_delta_bound.py:507-534` | `s6 == 2.0 * s`; synthetic rows; f1ll_critical_v6 "where twice it does" | **RETAKE** to `mod.PRE_D4_OVER_JM` (= 2.0, asserted) |
| Y6 | `python/tests/test_bench_t5_delta_bound.py:589-606` | the REPORT line text "(V6, Delta_DS16 = 2 Delta_tree: …)" and "FAIL under V6" | **RETAKE** to the new labels (§6.1) |
| Y7 | `python/tests/test_bench_t5_delta_bound.py:303-311` | kernel a₂ == `a_cos2phi` | **NONE** (docstring coefficient → ½) |
| Y8 | `python/tests/test_kernel_reference.py:98-174` | the port gate, Python side (`amplitudes` × 3 blocks, `dsigma`, `a_cos2phi`) | **RETAKE** per §5, mirroring T5 |
| Y9 | `python/tests/test_price_defaults.py:38-44` | `A2_TOY` | **RETAKE** to §2.1's halves; σ, A_zz, RATIOS unchanged (D7 retakes those) |
| Y10 | `python/tests/test_b1_model.py:1118-1129` (`a_cos2phi` as `c2p`) and `:1185-1209` (`toy_delta_gluon` carries no R; `R_COS2PHI_SHARED_PCT`, ratios); `python/tests/test_bench_t1_r_world.py:440-452`; `python/tests/test_rc.py:281-306`; `python/tests/test_bench_t4_t6.py:377-403` | `c2p` ratios; the kernel-factor ratio; `rc_scope` labels; T6 (plans `tensor-thirds`/`helicity-flip`, `n_spin_weights_seen == [0]`) | **NONE** |
| N2 | `python/tests/test_delta_convention.py` (new, SPDX) | — | **ADD**: the four §2.1 a₂ values of `default_inclusive_kernel(li6())` at rel 1e−12 (pre-D4 values in the comment), m = 0 = −2 × m = +1, `lg.a_cos2phi` closed form, `lg.DELTA_COS2PHI_NORM == 0.5` |

### 3.5 Validation harnesses, recipes, reference data

| # | file:line | what | action |
|---|---|---|---|
| V1 | `validation/benchmarks/t5_delta_moment_bound.py:134-199` | A7, S2 (`DS16_OVER_TREE = 2`), S3 (TIGHT includes V6) | **EDIT** per §6.1 (a D4 block in the header; constants; `saturation_row`, `measure`, `run`, `_print`, `plot`) |
| V2 | `validation/benchmarks/t5_delta_moment_bound.py:262-275` | "it does not decide which of A7 and V6 the tree's Delta follows" | **REWORD**: decided by the author (D4), not by this row |
| V3 | `validation/benchmarks/data/delta_moment_bound_sources.json` | LOI18 role/notes: "(9) is the tree's kernel term … the harness carries the other reading as V6" | **NONE** (a dated source record of 2026-10-26, sha256-pinned in the harness and the pytest); the header's D4 block says the note describes the pre-D4 kernel. Alternative: edit and re-record both digests together |
| V4 | `validation/benchmarks/t1_r_world.py:120-132` | docstring "In the kernel a_2 = -(1-y)/y^2 c_eff sin^2(theta_S) Delta/D_phi" | **REWORD** (½); the harness's number is a₂(R1998)/a₂(shipped), a **ratio, unchanged** (both halve exactly); re-run because the file changed (tiers rule) |
| V5 | `validation/benchmarks/t5_weighted_unweighted.py:81-85` | "the tensor a2 of the transverse plans is ~ 1e-3" | optional **REWORD** "~5e−4"; the row is θ_S = 0, unaffected |
| V6 | `validation/price_defaults.py:38-49` and `:294-318` (`effective_modulation`, `cos2phi_moment`) | A2(S) = a₂_eff/P_zz on `transverse-tensor` | **NONE** to code; absolute A2 halves, every ratio and the figure (ratios, `:1037-1147`) unchanged |
| V7 | `validation/dump_polligen_reference.py:280-297` (`asym.a_cos2phi`, `kern.amplitudes`), `:1368-1377` and `:1490-1491` (`_preflight_sd_phase`, run first in `main`) | how `xsec.json`'s `a_cos2phi`, `a2`, `dsigma` were made (polligen's `asym.a_cos2phi`, `kern.amplitudes`) | **ADD** a preflight `_preflight_delta_convention()` that measures polligen's coefficient (`asym.a_cos2phi(1, 1, 0, 0.1, 0.5)` against −2 / −1) and prints it; **REWORD** the README template's `xsec.json` paragraph (`:541-598`) and `validation/README.md` (generated; keep in sync by hand now) with the D4 note of §5 |
| V8 | `validation/reference/xsec.json` | 6Li: 24 of 30 states with a₂ ≠ 0 in each of 3 blocks; 7Li `kernel_scenario_rank2` and `kernel_tensor_gamma`: 32 of 40; `a_cos2phi` at 6 points per Δ-carrying block; the two θ_S = π/2 `dsigma` rows per block | **NONE** (never re-dumped from this tree; re-dumped only after the author updates the sibling, §5) |
| V9 | `validation/reference/b1_default_li6.json` | Δ *values* of the shipped toy (self-pin) | **NONE** |

### 3.6 Documents under the strict citation gate

| # | file:line | what | action |
|---|---|---|---|
| D1 | `docs/PHYSICS_CHANNELS.md:107-115` | §2 prose: "a₂ += −(1−y)/y²·c_eff·sin²θ_S·Δ/D_φ", "reducing to the HJM c_m … at J = 1" | **REWORD in place, 9 lines kept**: "a₂ += −½(1−y)/y²·c_eff·sin²θ_S·Δ/D_φ (Jaffe–Manohar's normalisation, D4)"; the HJM clause stays (geometry) |
| D2 | `docs/PHYSICS_CHANNELS.md:123` | "Rank-2 rate sector, massless HJM master formula": the formula with "− (1−y)/y²·c_m sin²θ_m·Δ cos 2φ" called "the HJM definition" | **REWORD in place** (one table line): ½, and "HJM's b-sector and c_m pattern with Jaffe–Manohar's Δ term — HJM has none"; references gain [JM89g] |
| D3 | `docs/PHYSICS_CHANNELS.md:124` | `tensor_gamma` row: "Δ_fake/(γ²b₁) ≈ 0.14–0.16" | **REWORD in place**: "0.28–0.32 in JM units (0.14–0.16 in the pre-D4 units)" |
| D4 | `docs/PHYSICS_CHANNELS.md:126` | the Δ row: coefficient, "Standalone identity: a_cos2φ … = −(1−y)/y²·Δ/D_φ", "**The repository carries no literature reference for Δ**" | **REWORD in place**: ½ twice; the references cell becomes "[JM89g] Eq. (16) (the normalisation adopted by D4), [SS90] (A7), [DS16] Eqs. (3)–(5), [LOI18] Eq. (6) — one coefficient; the LOI's Eq. (9), twice it, was the tree's until 2026-10-27. `toy_delta_gluon` stays an unsourced ansatz"; symbols cell gains `include/lipolgen/constants.hpp:<line>` `DELTA_COS2PHI_NORM` |
| D5 | `docs/PHYSICS_CHANNELS.md:209` | Δ models row | **REWORD in place**: "C_BAG enters at face value: the kernel's Δ is SS90's and DS16's (D4)" |
| D6 | `docs/PHYSICS_CHANNELS.md:568-575` | §13: "Δ_fake/(γ²b₁) ≈ 0.14–0.16"; "Δ (gluonometry) keeps its massless kinematic factor" | **REWORD in place** (numbers as D3) |
| D7 | `docs/PHYSICS_CHANNELS.md:988-1051` | references 1–58 (HJM89 = 1, JM89 = NPB 321 = 2, SS90 = 15) | **ADD at the end** (after 58; allow-table keys run to doc line 1041): 59 **[JM89g]** Jaffe, Manohar, *Nuclear gluonometry*, PLB 223 (1989) 218 (read as MIT preprint CTP#1712, Eq. (16)); 60 **[DS16]** Detmold, Shanahan, PRD 94 (2016) 014507, arXiv:1606.04505; 61 **[LOI18]** Maxwell *et al.*, arXiv:1803.11206, Eqs. (6) and (9) |
| D8 | `docs/theory/SPIN32_FINITE_GAMMA.md:919-925` | "a₂ ⊃ −(1−y)/y² · 3Q_NN · sin²θ_S · Δ / D_φ [ gluon transversity ]" | **REWORD in place** (½) |
| D9 | `docs/CONVENTIONS.md:75-96` | the `TENSOR_LL_SIGN` bullet, "What does NOT: |A_zz|, the whole Δ (cos 2φ) sector" | **NONE** (true of the sign) |
| D10 | `docs/CONVENTIONS.md:153-170` | `tensor_gamma` bullet: "Its size is Δ_fake/(γ² b₁) = 0.14–0.16" | **REWORD in place** (as D3) |
| D11 | `docs/CONVENTIONS.md` (end of file, after line 543) | — | **ADD** the §7 bullet at the END of the "Physics defaults that are a CHOICE" section — not after line 96, which would shift `docs/CONVENTIONS.md:207-211` (`EMC_VALENCE_DEPLETION_EPPS21`), `:325-327` ("X must be TIMELIKE") (RANGE_ALLOW keys) and `:415` (`KAPPA_NN_VIRTUAL`; an UNNAMED_ALLOW key) |
| D12 | `README.md:26` | inclusive row: "HJM spin-1 / spin-3/2 master formula: …, b₁, b₂, Δ (cos 2φ gluonometry)" | **REWORD in place** (§8) |

### 3.7 Other documents

| # | file:line | what | action |
|---|---|---|---|
| O1 | `docs/benchmarking/04_theory.md:110-113` | T-10 HJM "**CLOSED FORM** — b₁…b₄, the (1, −2, 1) helicity pattern, and Δ" | **REWORD**: drop "and Δ"; "(no Δ: A₊₋,₋₊ falls as 1/Q² there; Δ is T-16's)" |
| O2 | `docs/benchmarking/04_theory.md:127-152` | T-16 … T-19 and "The honest reading" | **REWORD**: T-16 relation "+ Eq. (16), the normalisation the kernel adopts (D4)"; T-17 "C_BAG in JM's normalisation — unconverted since D4"; T-18 "the kernel's Δ is DS16's (D4): the bound applies with no factor" |
| O3 | `docs/benchmarking/08_channel_priority.md:91` | rank 1: "18–28σ … (21–43σ on the toy …)"; the R price −15.6 / −15.1 / −7.1 / +4.6 % | **REWORD**: the σ are the letter's (twice-JM coefficient): 9–14σ / 10.5–21.5σ for the same Δ in JM units, still ≥ 5σ, **C2 = 3, no rank moves**; the R price is a RATIO and does **not** move |
| O4 | `docs/benchmarking/08_channel_priority.md:107-157` | §3.1: letter amplitudes (0.44–0.95) × 10⁻²; "Δ_fake/(γ²b₁) = 0.14–0.16"; the shipped A₂ = −1.06 / −1.30 / −1.43 / −1.33 × 10⁻³ "4–7× smaller" than `moment_A`; "the factor-2 normalisation … open" | **REWORD**: A₂ → −0.531 / −0.651 / −0.715 / −0.665 × 10⁻³; 8–14× smaller than the letter's *printed* amplitudes, 4–7× like for like; leakage 0.28–0.32 (JM units); "settled by D4" |
| O5 | `docs/benchmarking/08_channel_priority.md:327-328` | §4 row A: "xΔ = −a₂F₂[1 + y²/(2(1−y)(1+R))]" (the inversion with the pre-D4 coefficient); row C: "13–19 % … V6 … 6–9 % under … A7 … shown inconsistent"; development item "the factor-2 … a registry candidate" | **REWORD**: xΔ = **−2**a₂F₂[…]; row C per §6.1 (6.4–9.4 % is the generated amplitude; 13–19 % the pre-D4 one; letter 0.67–0.74 at face value, 1.35–1.48 for its printed amplitudes); "settled by D4" |
| O6 | `docs/benchmarking/08_channel_priority.md:392-399` | §7 revision log | **ADD** a D4 row ("no rank moved; C2 of rank 1 stays 3") |
| O7 | `docs/benchmarking/BENCHMARK_PLAN.md:89` | row 13 status text (headline at V6's 0.129 / 0.161 / 0.188; "registry candidate") | **REWORD** per §6.1 |
| O8 | `docs/benchmarking/BENCHMARK_PLAN.md:230` | U-5 "b₃, b₄, Δ … Detmold–Shanahan lattice A₂ for the φ meson — a different hadron" | **NONE** (optional: "same normalisation since D4; the distance is the hadron") |
| O9 | `validation/benchmarks/README.md:237` | row 13, `t5_delta_moment_bound.py` | **REWORD** per §6.1 |
| O10 | `docs/benchmarking/REPORT.md`, `REPORT.json`, `figures/t5_delta_moment_bound.png` | row 13 text, the figure (whiskers, legend) | **regenerate** (`make_report.py`), never by hand |
| O11 | `docs/references/REFERENCES.md:60-61` | HJM and JM-gluonometry entries "PLEASE DOWNLOAD … no preprint scan" | **REWORD**: preprint scans read in phase R1 (2026-10-27); journal texts not obtained |
| O12 | `docs/references/01_tensor-b1.md:86`, `:195-202`, `:621` | "NPB 312 defines b₁…b₄ and Δ"; "the b₁…b₄ and Δ definitions the tree transcribes" | **REWORD** (§8) |
| O13 | `docs/surveys/needs_survey.md:443` (the `test_recopseudo.py` row) and `:452` | survey of the sibling: "HJM NPB 312 Eq. (30)" for A = −[(1−y)/y²]…; the −(ε/2) closed form | dated note: "(HJM has no Δ term; the coefficient is the LOI's Eq. (9), twice JM's — plan_D4)" — or leave as the survey of the sibling it is |
| O14 | `docs/USAGE.md:150` ("rank-2 transfer + 1e-2 toy Delta"), `:336-352`, `:1176`, `:1233`, `:2754-2757` | "1e-2 toy Delta" comment; the R table (ratios); δ(cos 2φ) −18.1195 % (ratio); `tensor-all` label | **NONE** (optional in-place "(JM normalisation)" at 150; USAGE.md ranges are cited, so no added lines) |
| O15 | `docs/DEVELOPMENT_PLAN.md:250` | the leakage coefficient "Δ_fake/(γ² b₁) = 0.14–0.16" (a plan document) | optional dated note |
| O16 | `docs/HEPMC3_CONVENTION.md` | no example carries a₂; `spin_weight_k` of tilted categories change in files | **NONE** |
| O17 | `docs/open_items/run_2026-09-03/STATUS.md` registry rows 3, 12, 20 | the cos 2φ amplitude in ratio form (R; coherent ceiling; the P_zz divisor) | **NONE**: no registry row carries the coefficient; no label moves |
| O18 | `docs/open_items/run_2026-10-26/*` | dated records of the pre-D4 kernel (C2, A, D, E, F, SUMMARY, DECISIONS_PRICED item 4) | **NONE** (history); `--records` may ask for `(as of ca9b23a)` annotations where a cited line's text changes |

## 4. Pin decisions

**Re-take with a dated comment naming D4** (old value kept in the comment):
T1, T2, T3, T4 (C++); Y1–Y6, Y9 (Python); the port gate T5/Y8 in the form of
§5; the t5 harness's printed labels. Every new value is in §2 or §3 and is
exactly ½ (or ×2 for the leakage-as-Δ) of the old one where the old one is a
pure Δ-line quantity; w_min and A_T max are re-measured (§2.2).

**Keep comparing old numbers explicitly** (so the pre-D4 physics stays
visible, labelled):

- the **port gate** (T5, Y8): it keeps comparing to polligen's numbers, which
  carry the old coefficient, through the documented transformation of §5 — not
  re-dumped from this tree (that would make it a self-pin that passes by
  construction);
- the t5 row's **V6 column** (`saturation_v6`, `f1ll_critical_v6`, `status_v6`):
  kept, numerically unchanged, relabelled "the pre-D4 kernel / the letter's
  printed amplitudes (LOI Eq. (9), ×2)";
- T2: the oversized Δ doubles so that the *a₂* the test drives negative is the
  pre-D4 one, bit for bit;
- T4's comment keeps the sibling script's 0.14–0.16 as the pre-D4-unit value.

**New pins** (N1, N2): the constant's value, the closed form written out
independently, the four sweet-spot a₂ values, and the m-pattern. Without them
every surviving test would pass for either coefficient (T6, Y7 compare the
kernel with `a_cos2phi`, which move together).

## 5. The port gate after D4 — meaningful, and unable to pass for the wrong reason

`validation/reference/xsec.json` is dumped from the sibling's polligen
(`validation/dump_polligen_reference.py`), whose `a_cos2phi`
(`fastsim/polli_fastsim/asymmetries.py` lines 235–237) and kernel a₂
(`evgen/polligen/xsec.py` lines 732–733) carry −(1−y)/y². After D4 the tree
differs from it, in every Δ-carrying field, by a known amount. Design (both
`tests/test_reference.cpp` and `python/tests/test_kernel_reference.py`):

1. **Read the reference's own coefficient off its own numbers**, per block:
   c_ref = a_cos2phi_ref · D_φ,ref · y² / (−(1−y) Δ_ref) from the block's
   `tables` (f₁, f₂, delta) and `asymmetries` (y, a_cos2phi), at every point
   with Δ_ref ≠ 0. Require **all equal, and equal to 1 or to ½** (rtol 1e−12);
   anything else fails ("the reference carries an unknown Δ normalisation").
   Print which (doctest `MESSAGE`).
2. **Compare the tree to the literal Jaffe–Manohar target**, k = 0.5/c_ref
   with **0.5 written as a literal in the test** (not `DELTA_COS2PHI_NORM`, so
   a regression of the library constant fails here):
   - `a_cos2phi` == k · a_cos2phi_ref;
   - a₂ == a₂_ref − (1 − k) · c_eff(m) sin²θ_S · a_cos2phi_ref, c_eff from
     `tensor_moments(m).second` — this one formula covers the massless blocks
     (where a₂_ref is the Δ line) and the `kernel_tensor_gamma` blocks (where
     it is h₂ + the Δ line; the harmonics sub-block keeps its plain 1e−12
     comparison);
   - dsigma == dsigma_ref − (1 − k) · c_eff sin²θ_S · a_cos2phi_ref ·
     cos 2(φ − φ_S) · dsigma_unpol_ref/(2π) (W > 0 at every row: the
     modulation is O(10⁻³)).
   Everything Δ-free (`tables`, w_avg, a1, the harmonics, dsigma_unpol, the
   θ_S = 0 rows) keeps its unconditional 1e−12 comparison.
3. Effect: today c_ref = 1, k = ½, and the gate asserts that the tree is
   exactly JM's while polligen is exactly LOI-(9)'s. **It fails** if the tree
   regresses to the old coefficient (or to any other), if a re-dump carries a
   normalisation that is neither, or if a re-dump silently moves a Δ-free
   number. After the author updates the sibling and the JSON is re-dumped (the
   dump refuses moved numbers without `--accept-changed-numbers`, so that step
   is a human's), c_ref = ½, k = 1 and the gate is the plain 1e−12 comparison
   again **with no code change** — the MESSAGE then says so, and the D4
   comment can be retired at leisure.
4. Alternative, simpler, also safe: a literal factor ½ applied to the Δ part
   plus an assertion that c_ref == 1 (the gate then fails on the re-dump and
   the factor is removed in the same commit). Recommended: the data-driven
   form above.
5. The dump script gains `_preflight_delta_convention()` (prints polligen's
   coefficient before anything is written) and its README template — and the
   generated `validation/README.md`, kept in sync by hand until the next dump —
   gains: "`a_cos2phi`, `a2` and the θ_S ≠ 0 `dsigma` rows carry polligen's
   cos 2φ coefficient; the gates compare LiPolGen's (Jaffe–Manohar, D4)
   through the coefficient they read off the file".

The conftest `polligen` fixture (`python/tests/conftest.py:89-103`) imports
the sibling at test time for `test_hfs.py` only (HFS reading): no Δ path,
**NONE**.

## 6. Benchmark rows affected, and what their declared rules say

### 6.1 Row 13, `t5_delta_moment_bound` — PASS stays PASS; three things change

Declared rules (header `validation/benchmarks/t5_delta_moment_bound.py:212-233`):
PASS iff (a) the shipped Δ has 1 − |M₂|/B ≥ −1e−12 at Q² = 2, 5, 10 under A1–A8
and (b) every helicity-state rate ≥ −1e−12. **Neither the rule nor its
tolerance moves** (a bound, declared before the verdict). M₂ does not move
(Δ itself is unchanged), B does not move, so (a)'s numbers — 0.0644 / 0.0803 /
0.0941 — are bit-identical; (b) moves in the 8th decimal (§2.2). What D4
changes is the *meaning* of the columns, which the harness must say:

- **A7** ("the tree's Delta is the Delta of DS16",
  `validation/benchmarks/t5_delta_moment_bound.py:134-140`) stops being an
  assumption the second pass showed inconsistent and becomes the kernel's own
  convention. `DS16_OVER_TREE`
  (`validation/benchmarks/t5_delta_moment_bound.py:326-329`), "a derived ratio, re-derived
  by the pytest", becomes **1.0** — and the pytest re-derives 1.0 (Y1).
  "Never set from a result" stays true: it is set from the DS16 contraction.
- **V6** becomes a recorded reading of the *pre-D4 kernel and of the letter's
  printed amplitudes*: new constant `PRE_D4_OVER_JM = 2.0` drives
  `saturation_v6`, `status_v6`, `f1ll_critical_v6` (numbers unchanged);
  labels: "pre-D4 / letter as printed (LOI Eq. (9), ×2)".
- **TIGHT** (`validation/benchmarks/t5_delta_moment_bound.py:196-199`) loses V6 — back to the first pass's definition,
  V1 + V2(−) on the smaller gluon moment — because the normalisation is no
  longer an uncertainty of the kernel: toy 0.136 / 0.168 / 0.197, letter
  `moment_A` 1.508 / 1.466 / 1.469 (the first-pass edges recorded in C2 §4.3).
- the header gains a D4 block (text below); the REPORT generator string reads
  "|M2|/B = 0.0644 / 0.0803 / 0.0941 (Jaffe-Manohar normalisation, D4; the
  pre-D4 kernel's amplitude, x2: 0.1289 / 0.1605 / 0.1881); helicity-state
  rates min W = 0.949400 …"; the figure's open markers are relabelled and its
  whiskers end at the new tight edge; the letter rows' word "RECORDED: PASS
  centrally, FAIL under V6" becomes "PASS at face value, FAIL as printed".
- the vendored JSON is not touched (V3); the bound-side sub-row records both
  constants.

Proposed header block (to follow SECOND PASS):

    D4 (2026-10-27, author decision; docs/open_items/run_2026-10-27/
    DECISIONS.md D4, plan_D4_jm_convention.md).  The kernel's cos 2phi
    coefficient is Jaffe-Manohar's since this date, a2 = -(1/2) ((1-y)/y^2)
    c_eff sin^2(theta_S) Delta/D_phi (DELTA_COS2PHI_NORM = 1/2), which is
    DS16 Eqs. (3)-(5) contracted.  So A7 is now the kernel's own convention:
    DS16_OVER_TREE, re-derived by the pytest from Eqs. (3)-(5) against the
    kernel, is 1, and the central column is the shipped toy's GENERATED
    amplitude in DS16's units.  V6 is kept, numerically unchanged, as a
    RECORDED reading of the pre-D4 kernel and of the letter's printed
    amplitudes (LOI18 Eq. (9), twice Jaffe-Manohar: PRE_D4_OVER_JM = 2);
    it is in no edge and no verdict.  TIGHT loses V6 (V1 + V2(-) on the
    smaller gluon moment, the first pass's definition).  Nothing else moves:
    not M2 (the toy keeps 1e-2), not B, not the rule or its tolerance;
    sub-row (b)'s rates carry the halved Delta term (min W 0.9493997515 ->
    0.9493997990; A_T max 0.0506002485 -> 0.0506002010; A_T min, at
    theta_S = 0, unmoved).  The vendored file's LOI18 note is the second
    pass's reading of the pre-D4 kernel, dated 2026-10-26, and is not
    re-hashed.

**Order with D7.** `plan_D7_backend_defaults.md` re-points `shipped_kernel()`
to the D7 kernel (CT18NLO F₁, so the toy Δ ∝ F₁ moves M₂ and the saturation
|M₂|/B to 0.0760 / 0.0907 / 0.1023, and w_min to 0.9135) and keeps the letter rows on the toy kernel. The
two compose: D4's changes are structural (DS16_OVER_TREE = 1, V6 relabelled,
tight edge = |M₂|/B_tight, w_min re-measured with the halved Δ line) and must be
re-measured on whichever kernel the harness ships when D4 lands. §2.2's
numbers are for the ToyF2 kernel of `ca9b23a`.

### 6.2 The other rows

- **Row 12, `t1_r_world`**: its kernel factor is a₂(R1998)/a₂(shipped), a
  ratio of two Δ lines that both halve exactly — **bit-identical**; its rule
  (NMC R χ²) does not read a₂. Re-run only because its docstring changes (V4).
- **`t5_weighted_unweighted`**: the longitudinal `apar+` category, θ_S = 0,
  a₂ ≡ 0 — unaffected (docstring optional, V5).
- **T6 `t6_hepmc3_convention`**: plans `tensor-thirds` / `helicity-flip`, no
  tilt, `n_spin_weights_seen == [0]` — unaffected.
- **Every T1, T2, T3, T4 row**: reads F₂, g₁, b₁, form factors or tails, never
  a₂ — unaffected. `tiers.py --files src/core/xsec.cpp src/core/asymmetries.cpp
  include/lipolgen/xsec.hpp include/lipolgen/constants.hpp` prints all four
  **UNMAPPED** ("decide by hand"): by hand, D4 re-runs rows 12 and 13.
- **REPORT**: `make_report.py --all --plots` then `--compare` against the
  committed `REPORT.json` must show row 13's text and figure as the only
  differences attributable to D4 (row 12 identical).

## 7. The convention paragraph (for `docs/CONVENTIONS.md`, appended at the end)

> - **Δ normalisation — the cos 2φ coefficient.** `DELTA_COS2PHI_NORM = 0.5`
>   since 2026-10-27 (author decision D4, `open_items/run_2026-10-27/DECISIONS.md`):
>   Δ(x, Q²) is the double-helicity-flip amplitude A₊₋,₋₊ normalised like
>   F₁ = ⅓Σ_H A₊H,₊H — Jaffe and Manohar's Δ — and the kernel's term is
>
>       a₂ = −½ · (1−y)/y² · c_eff(m) sin²θ_S · Δ/D_φ,   c_eff = 3Q_NN (= 3m² − 2 at J = 1),
>
>   so a pure m = ±1 state with its axis transverse to the beam carries
>   a₂ = −½(1−y)/y² Δ/D_φ, and a fill carries Σ_m p_m c_eff(m) = P_zz. Every
>   source that derives the term prints this coefficient: Jaffe–Manohar,
>   *Nuclear gluonometry*, PLB 223 (1989) 218, Eq. (16) (read as the MIT
>   preprint CTP#1712); Sather–Schmidt, PRD 42 (1990) 1424, Eq. (A7);
>   Detmold–Shanahan, PRD 94 (2016) 014507, Eqs. (3)–(5) contracted with the
>   lepton tensor; the JLab LOI arXiv:1803.11206, Eq. (6); and, in
>   helicity-response language, F[U T_TT] = −Δ with T_TT cos 2φ_TT = ½ sin²θ
>   cos 2φ (Pascalutsa–Vanderhaeghen arXiv:2610.09862, Eqs. (26), (27b),
>   (29b), (A8a)). The coefficient without the ½ is the LOI's own Eq. (9),
>   which the LOI does not derive and which contradicts its Eq. (6), and a 2019
>   proceedings reprint of it; this tree carried it from its first port until
>   2026-10-27, transcribed from the sibling's `docs/Discussions.pptx` p. 5,
>   which it had called the Hoodbhoy–Jaffe–Manohar transcription.
>   Hoodbhoy–Jaffe–Manohar (NPB 312 (1989) 571) has no Δ term: there the
>   double-helicity-flip amplitude falls as 1/Q², and the leading-twist gluon
>   tower is Jaffe–Manohar's. The parton-model side — Δ = −(α_s/2π) Tr Q² x²
>   ∫ₓ¹ dy/y³ δG(y), its moments, h₁TT — is one normalisation in every source
>   up to a sign, so external Δ numbers enter UNCONVERTED: `C_BAG` = −0.012
>   (∫x Δ dx = C_BAG α_s) and the gluonic Soffer bound of row 13. What moved:
>   every Δ-induced cos 2φ amplitude halves per unit Δ — the shipped
>   `toy_delta_gluon(…, 1e−2)` keeps its scale, so its amplitude halves — and
>   the Δ-equivalent of the finite-γ b-sector leakage doubles, Δ_fake/(γ²b₁)
>   0.14–0.16 → 0.28–0.32. What did not: Δ itself, the b-sector and A_zz,
>   every R or F₂ ratio of the cos 2φ amplitude, the `tensor_gamma` harmonics,
>   the coherent channel's c₂, and any run at θ_S = 0 (the default
>   `tensor-thirds` plan). Setting the constant to 1.0 restores the pre-D4
>   kernel, which is the whole of the reversal; `tests/test_xsec.cpp` pins
>   both the closed form and the value, and the port gate reads polligen's
>   coefficient off `xsec.json` (it carries LOI Eq. (9)'s until the sibling
>   changes).

## 8. The corrected attribution (the texts)

- **b-sector**: "the spin-1 rate sector — F₁ + (2/3)a_m b₁ + (1−y)/(xy²)(F₂ +
  (2/3)a_m b₂), a_m = ¼c_m(3cos²θ_m − 1), c_m = 3m² − 2 → (1, −2, 1) — and the
  b₁…b₄ decomposition are Hoodbhoy–Jaffe–Manohar, NPB 312 (1989) 571 (b₃, b₄
  higher twist there); the overall sign of b₁ is the literature's
  (`TENSOR_LL_SIGN`)."
- **Δ term**: "the double-helicity-flip term −½(1−y)/y² c_m sin²θ_m Δ cos 2φ is
  Jaffe–Manohar, PLB 223 (1989) 218, Eq. (16) (normalisation adopted 2026-10-27,
  D4), the same in Sather–Schmidt PRD 42 (1990) 1424 Eq. (A7), Detmold–Shanahan
  PRD 94 (2016) 014507 Eqs. (3)–(5) and the JLab LOI arXiv:1803.11206 Eq. (6).
  HJM NPB 312 has no Δ term."
- Sites: the master-formula comments that name Hoodbhoy-Jaffe-Manohar,
  `include/lipolgen/asymmetries.hpp:9-27` and `include/lipolgen/xsec.hpp:13-32`;
  `docs/PHYSICS_CHANNELS.md:123` and `:126` (the double-helicity-flip row),
  `README.md:26` (the `inclusive` row's master formula),
  `docs/benchmarking/04_theory.md:112`, `docs/references/01_tensor-b1.md:86`,
  `docs/references/01_tensor-b1.md:200-201` ("NPB 312 defines b₁…b₄ and Δ"), `:621`.
  `include/lipolgen/spin.hpp:23` ("The J=1 conventions
  match Hoodbhoy-Jaffe-Manohar NPB 312:571 (c_m = 3m^2-2)") is right as it
  stands and carries an UNNAMED_ALLOW pin — leave it.

## 9. Ordered implementation steps, with the gate after each

0. **Baseline, before any edit.** `git rev-parse HEAD`; the strict gate
   tallies (`python3 validation/check_physics_channels_links.py`: 1240 refs /
   95 ranges / 7 external / 0 broken / 24 allow-listed; SPIN32 19 / 115 / 0 /
   0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0 / 0, measured 2026-10-27 13:3x) and
   `--records-only`'s tally; the pre-D4 dump
   `python3 d4_before_after.py dump pre_d4.json` (§10; one is already in session
   scratch, `d4/pre_d4_ca9b23a.json`, sha256 `0334237e…`). Decide the order
   with D5/D6/D7 (they touch `test_price_defaults.py`, the t5 harness,
   `PHYSICS_CHANNELS.md` and the REPORT): one decision per commit, the later
   one re-measures on top of the earlier.
1. **Library**: L1, L2, L3 (+ P1); comments L5–L10, L13 in place. Gate:
   `cmake --build build -j2 && ./build/lipolgen_tests` — **expected to fail in
   exactly** T1, T2, T3, T4 and the Δ-carrying checks of T5, and — the
   module is rebuilt with the library — on the Python side exactly Y1 (the
   S2 coefficient and ratio), Y4 (w_min, A_T max), Y8 and Y9 (and nothing
   else; any other failure means the inventory missed a site — stop). Then
   `python3 d4_before_after.py check pre_d4.json` must print "D4 contract
   holds" (5310 states exact half, 1593 `tensor_gamma` states half on the Δ
   line, every Δ-free quantity bit-identical).
2. **C++ tests**: T1–T5 re-taken, N1 added. Gate: `./build/lipolgen_tests`
   all pass; case count +1 (N1).
3. **Python tests and harness**: Y1–Y9, N2; V1–V5. Gate: the touched files
   (`python3 -m pytest python/tests/test_bench_t5_delta_bound.py
   python/tests/test_kernel_reference.py python/tests/test_price_defaults.py
   python/tests/test_delta_convention.py python/tests/test_module.py -q`),
   then the full `python3 -m pytest python/tests -q` and
   `LIPOLGEN_BENCH=1 python3 -m pytest python/tests -k bench -q` when the
   cores are free.
4. **Rows and report**: `python3 validation/benchmarks/t5_delta_moment_bound.py
   --plot` (PASS; numbers per §6.1), `python3 validation/benchmarks/t1_r_world.py`
   (identical kernel factor), then `make_report.py --all --plots` and
   `make_report.py --compare <committed REPORT.json>` (differences: row 13
   only), `make_report.py --check`.
5. **Documents**: D1–D12, O1–O18 as marked. Gate: the strict citation gate.
   Expected: EDITED fingerprints on the blocks whose text D4 changes
   (`src/core/xsec.cpp:260-284`, the tensor block of `tensor_amplitudes`;
   `include/lipolgen/xsec.hpp:13-21`, the header's `c_eff` master formula;
   `include/lipolgen/asymmetries.hpp:9-21`, the spin-1 master formula; and any
   other range the in-place edits touch) — re-read each citing row, then
   `--record-ranges` (it refuses a block the sentence no longer evidences);
   MOVED blocks and named citations from P1/L3's inserted lines — `--fix`.
   Never add a line to `docs/PHYSICS_CHANNELS.md` above line 1041, to
   `docs/CONVENTIONS.md` before its end, or to `tests/test_tensor_gamma.cpp`
   above line 407 (allow-table keys are doc lines and pinned target lines).
   Then `--audit-ranges` (0 refused), `--records-only` (no new R6 failures;
   a failure on a dated record is answered with `(as of ca9b23a)`, never by
   re-pointing history), `python3 -m pytest python/tests/test_doc_link_gate.py
   -q`, `python3 validation/check_spdx_headers.py`.
6. **Records**: `phase_I4_jm_convention.md` with every number of §2 re-taken
   on the new build, the gate tallies before/after, and the `--compare`
   output; the board; the D4 entry corrections (§14).

## 10. Verification plan

1. **The exact-half contract** — `d4_before_after.py` (session scratch
   `d4/`): on four kernels (shipped ⁶Li; the reference gate's toy ⁶Li; the same
   with `tensor_gamma` and b₃ = 0.05 f₁, b₄ = −0.02 f₁; a ⁷Li scenario with
   δ₃₂ = −10⁻² f₁) × 708 physical (x, Q², y) points × every m and θ_S ∈ {0,
   0.7, π/2}: a_cos2phi and every non-γ a₂ exactly ½, the γ-path a₂ half on
   its Δ line, w_avg / a₁ / h₀ / h₁ / h₂ / dsigma_unpol bit-identical. Its
   logic was checked here against an emulated D4 (passes), an emulated wrong
   factor ¼ (refused) and today's build (refused).
2. **The four cells** (§2.1) and the t5 positivity numbers (§2.2) re-measured
   on the new build, equal to the predictions (exactly, resp. to 1e−12).
3. **The finite-γ statement** (§2.3): the F[U T_TT]-equivalent of the Δ line
   is −1.000000001 at γ² = 3.5 × 10⁻¹⁰; the Δ line identical on both paths.
4. **Generated events**: a default run (`tensor-thirds`) bit-identical before
   and after (npz sha256); a `transverse-tensor` run of a few thousand events:
   x, Q², m, cell columns identical, φ different, `effective_modulation().a2`
   exactly ½.
5. **The port gate in both states**: green on the real `xsec.json`; red on a
   scratch copy whose `a_cos2phi` is scaled by 0.7; green on a scratch copy
   whose Δ-carrying fields are transformed to the JM coefficient (c_ref = ½,
   the post-sibling state); red on today's build (the pre-D4 library).
6. **Adversarial greps**: no spelled `(1.0 - y) / (y * y)` without the constant
   outside tests; no "Hoodbhoy" adjacent to Δ outside dated records; no bare
   "0.14–0.16" without "pre-D4"; no "Δ_DS16 = 2 Δ_tree" stated as current.
7. The gates of §9, with their tallies recorded.

## 11. Risks

1. **A half change**: one site edited, not the other. Caught by N1 (closed
   form) and the existing kernel-vs-`a_cos2phi` tests together; by §10.1.
2. **A gate passing for the wrong reason**: a re-dump from this tree, or a
   factor that tracks the library constant. §5 uses a literal JM target and
   reads the reference's coefficient from the file.
3. **Line shifts**: the strict gate's allow tables are keyed by document
   line and by pinned target line (`include/lipolgen/constants.hpp:139`
   "Sather-Schmidt bag-model sum-rule coefficient", include/lipolgen/xsec.hpp:19
   (written plain: no token of that line can pass the records gate's point
   rule), `src/core/asymmetries.cpp:85` `TENSOR_LL_SIGN`,
   `tests/test_tensor_gamma.cpp:152` and `:382` (two `tensor_gamma` TEST_CASEs),
   `docs/CONVENTIONS.md:415` `KAPPA_NN_VIRTUAL`);
   in-place edits and appending at the end avoid every one. The dated records
   (R6) cite `src/core/xsec.cpp:281` (`qc.second * st * st * t.delta`) and
   neighbours by line; in-place edits
   keep their shared tokens. (Incidental, not D4's: `docs/theory/SPIN32_FINITE_GAMMA.md`
   line 743 names ":406" of `tests/test_tensor_gamma.cpp` as the Table 1 row-3
   test, but line 406 is blank today and the TEST_CASE sits at 407; the strict
   gate does not parse a bare point continuation, so it never saw it.)
4. **Concurrent decisions on the same files**: D7 re-points the t5 shipped
   kernel and re-takes `test_price_defaults.py`; D5/D6/D7 all edit
   `PHYSICS_CHANNELS.md` and regenerate the REPORT. Sequential commits; the
   second re-measures; nobody adds lines to `PHYSICS_CHANNELS.md` above 1041.
5. **Statistical tests on new streams** (T8): expected to pass, must be run.
6. **Output provenance**: a pre-D4 and a post-D4 file with tilted categories
   are indistinguishable by `meta` unless L13 is done (open question 5).
7. **Tree vs letter**: until the sibling changes, the tree's cos 2φ amplitude
   for the letter's Δ is half the letter's printed one; any comparison needs
   the ×2, stated wherever both appear (08, BENCHMARK_PLAN, the t5 labels).
8. **The refusal threshold** for oversized Δ moves ×2 in Δ (L14).
9. **The b-basis question** (§12 item 1) is numerically large in the shipped
   configuration (§2.5); D4 neither creates nor resolves it, but halving the
   Δ line changes which reading a reader might infer.

## 12. Open questions for the author

1. **The b-basis and the separate Δ (possible double counting) — first.**
   Sather–Schmidt footnote 8 (p. 1428, read on the page image): "In the
   Bjorken limit they are related by b₁ = Δ₁ + ½Δ₄ and b₂ = Δ₂ + xΔ₄", with
   Δ₄ = Δ (`r1a/ss90.txt` line 345); their (A8) gives a longitudinally polarized
   target's rate with −⅓(xy²Δ₁ + λΔ₂), so A_zz measures SS90's Δ₁, while HJM's
   b₁ = Δ₁ + Δ/2. At finite Q², in the HJM basis, Pascalutsa–Vanderhaeghen
   (A9a): ⅓b₂ − b₃ = νΔ — so HJM's b₃ carries −νΔ, and the tree's own
   f_tt = −(γ²/x)(b₂/6 − b₃/2) = −(M/ν)(⅓b₂ − b₃) (inferred: the identity
   γ²/(2x) = M/ν) would produce F[U T_TT] = −Δ by itself on the `tensor_gamma`
   path. The tree adds Δ separately with b₃ = b₄ = 0 by default, which is
   consistent **only if its b-slots are read Δ-free** (SS90's Δ₁, Δ₂, …):
   true of every nucleonic b₁ model (Δ ≡ 0 for nucleons, JM's theorem), but
   the shipped toy Δ comes with no b-slot adjustment, and its Δ/2 is 66–114 % of
   the shipped b₁ at the sweet spots, 984 % at x = 0.3 (§2.5). Choices: state
   the convention "b-slots are Δ-free; Δ carries the double-flip alone"
   (documentation only), or treat the shipped b₁ as HJM's and subtract Δ/2
   (an O(1) change of the shipped A_zz at the sweet spots), or something else.
   **Recorded, not acted on.**
2. **J = 3/2**: JM state the analysis "applies without modification to a
   spin-3/2 target" (`r1a/jm89.txt` lines 93–94), so D4's ½ applies to the
   `delta_32_func` slot through the same line; the c_eff = 3Q_NN geometry at
   J = 3/2 remains the open J = 3/2 normalisation item (plans/04 #14). No
   shipped ⁷Li number moves (its slot is empty).
3. **The t5 TIGHT edge**: dropped V6 (§6.1). Keep V6 in it until the letter
   is updated instead?
4. **The coherent `amp`** (0.01 at P_zz = 1, band 3e−3–1e−2): an amplitude,
   untouched by D4. If it was chosen to mirror the inclusive toy's amplitude
   rather than as a coherent-channel scenario in its own right, say so.
5. **Provenance in files**: L13's banner text, a `meta` key, or a version
   bump, so post-D4 files say which coefficient their a₂ carries?
6. **A named constant** (`DELTA_COS2PHI_NORM`, recommended) **or** a literal
   ½ at the two sites?
7. **The vendored t5 sources file** (V3): leave as the dated record
   (recommended) or edit and re-record its digests?
8. The Δ dilution ⅓ vs the rank-2 transfer 0.3073 (C2 §9) — unchanged by D4.
9. **Δ inside the finite-γ path** (F[U T_TT] ⊃ −Δ, inserted as
   F_TT ⊃ TENSOR_LL_SIGN·Δ so that the Δ sector stays sign-independent) would
   give Δ the finite-γ kinematics SS90's (A3) carries; it moves the Δ term
   by O(γ²) (0.04 %–4 %, §2.3) and redefines the extracted function — a
   separate decision, not D4.

## 13. Sibling sites for the author (PolarizedLithiumSim, read-only here)

| file:line (sibling) | what | change for JM |
|---|---|---|
| `fastsim/polli_fastsim/asymmetries.py` lines 4–9 | master-formula docstring: "Hoodbhoy-Jaffe-Manohar NPB 312:571", Δ line without ½ | ½; attribution (§8) |
| `fastsim/polli_fastsim/asymmetries.py` lines 235–237 | `a_cos2phi` | ½ |
| `evgen/polligen/xsec.py` lines 17, 19–21 | docstring a_2; "Spin-1 (Hoodbhoy-Jaffe-Manohar NPB 312:571 …)" | ½; attribution |
| `evgen/polligen/xsec.py` lines 732–733 | the a_2 line | ½ |
| `evgen/scripts/tensor_gamma_leakage.py` lines 14, 34, 43, 77, 202 | Δ_fake = −a₂ y² D_φ/(1−y); 0.14–0.16; 0.23–0.26; 0.109 % → 0.17 % | −2a₂…; 0.28–0.32; 0.46–0.52; fractions ×2 |
| `evgen/tests/test_xsec_identity.py` lines 306, 328 | −(ε/2) closed form | −(ε/4) |
| `evgen/tests/test_recopseudo.py` lines 561, 576 | "HJM NPB 312:571 (1989) Eq. (30)"; −(1−y)/(y y) | attribution; ½ |
| `evgen/tests/test_tensor_gamma.py` | the leakage pins T4 ports | ×2 (JM units) |
| `fastsim/polli_fastsim/fom.py` line 286, `fastsim/scripts/money_delta*.py` (e.g. `money_delta.py` line 89), `evgen/scripts/money_cos2phi.py` | consumers of `a_cos2phi` / the kernel a₂ | amplitudes halve for the same Δ; significances halve |
| `paper/main.tex` lines 199–210 | Eq. (1) (`eq:master`), Δ term without ½, cited to `Jaffe:1989xy,Hoodbhoy:1988am` | ½; the Δ term is Jaffe–Manohar's |
| `paper/main.tex` lines 232–237 | Eq. (3): A = −(1−y)/y² Δ/D_φ, Δ̂ = −Â y² D_φ/(1−y) | −½…; Δ̂ = −2Â y² D_φ/(1−y) |
| `paper/main.tex` line 116; `paper/table1.tex`; `paper/figs/fig2_caption.tex` | 21–43σ; injected amplitudes; Â = 7.31 × 10⁻³, 8.88 × 10⁻³ | halve for the same Δ (or Δ doubles for the same amplitude — which the letter means is the author's) |
| `plans/07_plb_letter_gluonometry.md` lines 52, 105, 108 | sweet-spot amplitudes (0.44–0.95) × 10⁻²; 21–43σ / 32–64σ | as above |
| `docs/Discussions.pptx` p. 5 | LOI Eq. (9) transcribed and later called HJM's | annotate |
| `docs/code_review_2026-08-25.md` §3 | normalisation left open "for want of Jaffe–Manohar 1989" | settled |

After the update: LiPolGen re-dumps `validation/reference/xsec.json`
(`dump_polligen_reference.py --accept-changed-numbers`); the §5 gate turns
into the plain comparison by itself; `t5`'s V6 label "the letter as printed"
is then retired by a dated edit.

## 14. The D4 entry of `DECISIONS.md` — proposed corrections

1. "for a pure transverse m = ±1 state" next to "c_eff sin²θ_S" is redundant
   and reads as if the formula held only there: say "a₂ = −½(1−y)/y²·c_eff(m)
   sin²θ_S·Δ/D_φ for a pure state m on an axis at θ_S (c_eff = 3Q_NN; for
   m = ±1 transverse, a₂ = −½(1−y)/y² Δ/D_φ)".
2. "as in Jaffe–Manohar PLB 223 (1989) 218": add "read as its MIT preprint
   CTP#1712, Eq. (16)" — the journal text was not read.
3. "every cos 2φ amplitude the generator produces halves": "every
   **Δ-induced** cos 2φ amplitude — fills with sin θ_S ≠ 0; the default
   `tensor-thirds` output does not change; the O(γ²) b-sector leakage does
   not change, and its Δ-equivalent doubles (0.14–0.16 → 0.28–0.32 γ²b₁)".
4. "after which `moment_A` / `moment_B` represent the bag model at face
   value": in **this tree** they do so from D4 on (the kernel is JM's); what
   waits for the sibling is the **letter's printed amplitudes**, which encode
   2 × C_BAG until its `asymmetries.py`, `xsec.py`, Eqs. (1) and (3) and
   the inversion Δ̂ change.
5. Add: the attribution of the Δ term to HJM (`xsec.hpp`, `asymmetries.hpp`,
   `PHYSICS_CHANNELS.md`, `README.md`, `04_theory.md`, the reference survey)
   is corrected — HJM NPB 312 has no Δ term; and the b-basis question
   (§12 item 1) is filed for the author, not acted on.
6. Add: the port gate keeps comparing against polligen's (old-coefficient)
   numbers through the coefficient it reads off `xsec.json` (§5).

## 15. Commands, scratch, runtimes

All from `/home/cpeng/Projects/polli/PoliGen` after `source env.sh`; scripts in
the session scratch `d4/` (outside the repository), each under 1 s:
`measure_a2.py` (§2.1), `measure_pos.py` (§2.2, the harness's own
`helicity_positivity` and a halved re-run), `measure_gamma.py`,
`measure_gamma2.py` (§2.3), `measure_oversized.py` (§2.4), `measure_b1.py`
(§2.5), `ds16_ratio.py` (the t5 pytest's DS16 contraction, copied verbatim, §1),
`d4_before_after.py` (§10.1; `dump` wrote `pre_d4_ca9b23a.json`, 708 rows /
6903 states, sha256
`0334237e9416d6537b38240b7f8f5b34fd6c9673ceba955f319032305b2c8ba3`). Also
run: `validation/benchmarks/tiers.py --files …` (§6.2) and the strict
citation gate (tallies in §9 step 0). Not run: any build, the test suites,
`make_report.py`, any generator sample.

## Adversarial check (2026-10-27)

The adversarial check of this design, run on 2026-10-27 against `ca9b23a`
with nothing in the repository edited, reproduced in the checker's own words
(nothing trimmed, no number changed). Its file and line references are the
checker's, at `ca9b23a`, and are written as plain text, not as citations the
records gate checks.

### Verdict

SOUND IN PHYSICS AND NUMBERS, NOT READY TO IMPLEMENT AS WRITTEN. I checked the core claims against the papers and by re-running them on the ca9b23a build, and they hold. The ½ is printed in JM preprint Eq. (16), SS90 (A7) and LOI Eq. (6), and the DS16 contraction confirms it. The library evaluates the coefficient in only two places, src/core/xsec.cpp:281-282 and src/core/asymmetries.cpp:91, and applying the factor there gives exactly half, bit for bit. Every new pin value in the design reproduced exactly (sweet-spot a2, A2_TOY, leakage ×2, t5 w_min and A_T max, DS16 ratio 0.9999967–0.9999999), as did the finite-γ equivalence (−2Δ before D4, −Δ after) and the port gate's c_ref = 1. Fix these before implementing: (1) P1's one-line insertion into python/bindings.cpp. Contrary to the design, it is not fully --fix-able: SPIN32:986's range python/bindings.cpp:2187-2199 lands on a blank line and stays broken after --fix, and the records gate goes from 0 to 5 broken. Putting the new attr on the C_BAG line, or after the last cited line (5107), avoids both. (2) The step-1 expected-failure list is wrong: T4 does not fail. (3) One pin was missed (test_bench_t5_delta_bound.py:435), plus several text sites. (4) Six misstatements: the SS90 (A7) quote, the dump-script guard, spin_weight_k in files, the b-basis 'double counting', the claim of one parton-side normalisation, and N1's m values for J = 3/2.

### Sites the design missed (11)

| file | line | what |
|---|---|---|
| `python/tests/test_bench_t5_delta_bound.py` | 435 | `assert s["ds16_over_tree"] == 2.0` (the bound-side sub-row pin). It fails as soon as DS16_OVER_TREE becomes 1, and it is not among Y1–Y7. It must be retaken, and should also pin the new PRE_D4_OVER_JM field. |
| `validation/benchmarks/t5_delta_moment_bound.py` | 61-65 | The header's MODELS paragraph says, in the present tense: "The kernel builds the cos 2phi term as a2 += -((1-y)/y^2) c_eff sin^2(theta_S) Delta/D_phi". This is outside V1's listed ranges (134-199, 262-275, 326-329, …), so after D4 it keeps the old coefficient. The §10.6 grep for `(1.0 - y) / (y * y)` would not catch this spelling. |
| `tests/test_tensor_gamma.cpp` | 307 | "puts Delta_fake/(gamma^2 b1) at 0.14-0.16 -- the 1/6 of Eq. (17e)", in the TEST_CASE at 300-336. T4 covers only 338-378, so this comment keeps the pre-D4 unit (JM units: 0.28–0.32, the 1/3 of Eq. (17e)). |
| `docs/references/REFERENCES.md` | 424 | "#3 Hoodbhoy–Jaffe–Manohar 1989 (which *defines* b₁…b₄ and Δ …)". This is the same misattribution O11 fixes at :60-61, but this line is not listed. |
| `include/lipolgen/constants.hpp` | 26-28 | The TENSOR_LL_SIGN doc calls docs/Discussions.pptx p.5 "the transcription of Hoodbhoy-Jaffe-Manohar". The design's own §1 table names this line as the origin of the misattribution (R1: the slide is LOI Eq. (9)), but no rewording is listed; L4 covers only 31-40. Siblings of the same phrase: docs/CONVENTIONS.md:82-83 ("the repository's own transcription of Hoodbhoy–Jaffe–Manohar") and docs/PHYSICS_CHANNELS.md:123 ("flipped from the pre-2026-08-29 +1 HJM transcription"). |
| `docs/surveys/needs_survey.md` | 52, 59, 70, 198 | These carry the old coefficient ("a_2 = -(1-y)/y^2 * c_eff(m) …", "`a_cos2phi` … `-(1-y)/y² * Delta / D_phi`") and "HJM set {F₁, F₂, b₁…b₄, Δ, g₁, g₂}". O13 lists only :443 and :452; either annotate all of them or none. |
| `include/lipolgen/sf.hpp` | 529 | The `make_delta_toy` doc reads "`scale` is the peak Delta/F1". It is L10's twin and is not listed. It is also wrong before D4: scale·x^0.3(1−x)^4 peaks at 0.337·scale. |
| `docs/DEVELOPMENT_PLAN.md` | 28 | "xsec master formula dσ/dx dQ² dφ (HJM ⊕ CW)" credits the whole master formula to HJM. O15 lists only :250. |
| `docs/theory/SPIN32_FINITE_GAMMA.md` | 986 | Cites python/bindings.cpp:2187-2199. After P1's insertion, line 2187 is the blank old 2186 line. In a scratch emulation, --fix left this range broken ("the first line of the range is blank") and reported its fingerprint as orphaned. It needs a manual re-point to :2188-2200 plus --record-ranges, or no line shift at all. |
| docs/open_items/run_2026-09-03/phase_A_cdbonn.md, phase_D_sf_injection.md, phase_D_small_items.md | 482; 437 and 655; 279 and 450 | R6 breaks caused by P1's one-line shift, measured as `--records-only` going from 0 to 5 broken. They are: →bindings.cpp:1983 (no shared token); →:4079 (target blank); →:1535 (anchor-read pin 'arg("emc_ratio' moved off the line); →:5011-5015; →:5053 (anchor-read pin moved). The two anchor-read pins can only be re-dated, which collides with the design's 'never re-point history'. |
| `docs/PHYSICS_CHANNELS.md` | 126 | While D4 rewords this row, its symbols cell still cites src/core/xsec.cpp:266 `a2`. That line is the tensor_gamma h2 add, not the Δ line D4 edits (281-282). The gate passes only because the token `a2` is on the line. |

### Statements found wrong (11)

1. §3.2 P1 says the insertion "shifts 4 fingerprinted ranges and 3 named citations into the file — all --fix-able", and §9 step 5 expects "no new R6 failures". Both are false. I emulated L3 + P1 on a git-archive copy of ca9b23a. --fix re-pointed bindings.cpp :4158→:4159, :4404→:4405, :5107→:5108, constants.hpp :145→:149, and SPIN32 ranges 2341-2342, 2227-2229 and 2488-2501. It left python/bindings.cpp:2187-2199 broken, and the SPIN32 gate stayed at 1 broken. --records-only went from 0 to 5 broken.
2. §9 step 1 says the build "expected to fail in exactly T1, T2, T3, T4 and the Δ-carrying checks of T5". This is wrong for T4. tests/test_tensor_gamma.cpp:371 computes `d_fake = -h.h2 * p.y * p.y * dphi / (1.0 - p.y)` from the γ-path harmonic of a hand table with `t.delta = 0.0` (:60), so it never evaluates the Δ line and still passes after the library edit. Under the plan's 'any other deviation means a missed site — stop' rule, this is a false trip. The T4 retake is a relabel to JM units, not a catch.
3. §1 table, SS90 row, quotes (A7) as "[…+ ½(xy²Δ₁ + λΔ₂) − ½ xλΔ₄ cos 2φ …]". The p. 1427 page image prints "[xy²F₁ + λF₂ + ⅙(xy²Δ₁ + λΔ₂) − ½xλΔ₄cos2φ − ¼(xy²Δ₁ + λΔ₂ − xλΔ₄)sin²β(1+cos2φ) + 2x(1−y/2)Δ₃sin²βcosβ(1+cos2φ)]". The Δ₄ coefficient D4 rests on is right; the Δ₁ coefficient is misquoted (⅙, not ½).
4. §5 item 3 says "the dump refuses moved numbers without --accept-changed-numbers", and §13 gives the command `dump_polligen_reference.py --accept-changed-numbers`. That flag and guard belong to validation/dump_b1_default_li6.py (its docstring, lines 39 and 57-62). dump_polligen_reference.py parses no argv, runs only `_preflight_sd_phase()` (:1491), and rewrites all seven reference files atomically. The proposed `_preflight_delta_convention()` would only print; to be a guard it must refuse.
5. L16, O16 and the summary say that tilted-category "spin_weight_k" (HepMC3) change in files. No generator fills Event::spin_weights: only src/hepmc/hepmc_writer.cpp:134 reads it, and PHYSICS_CHANNELS.md:232 says it is "a declared, HepMC3-exported slot that no generator in this repository populates". `weights_for` changes only through the API. In files, what changes is φ, the e′ and X four-vectors, later draws of the per-event stream (the target hook at generator.cpp:138, the species draw at :144, the bridge's draws at pythia_bridge.cpp:435/465/656) and rc_tensor_lo/hi in slot 0.
6. §12 item 1 and risk 9 (b-basis "possible double counting … numerically O(1)") are misframed. SS90 (A2), on the page image, prints "Δ₁ = −½(A₊₊,₊₊ + A₊₋,₊₋ − 2A₊₀,₊₀)", which is the helicity-diagonal b₁ that A_zz measures; the design's own reading of (A8) says the same. The massless kernel, with its b-slot read as that b₁ and Δ ≡ Δ₄, reproduces SS90 (A7)/(A8) term by term at leading power: +⅙ for m_x = +1, −⅓ for m_z = +1 with TENSOR_LL_SIGN = −1, and −½ on Δ₄ after D4. Miller's and HERMES's b₁ are A_zz-type, so on the default path there is no double counting. The (Δ/2)/b₁ = 1.14/1.02/0.66/0.68 sizing (which I reproduced) has no bearing on consistency, and the offered option 'subtract Δ/2 (an O(1) change of the shipped A_zz)' would introduce an error. A genuine O(γ²Δ) basis question exists only on the opt-in tensor_gamma path. R1 §6 item 5 frames this correctly.
7. The §7 convention paragraph says "The parton-model side … is one normalisation in every source up to a sign, so external Δ numbers enter UNCONVERTED". R1 §2.3 and §6 item 2 contradict this. The JM preprint's Eq. (29) prints "∫₀¹dx xΔ(x,Q²) = (α_s(Q²)/8π) Tr Q² (1/(p⁺)²)⟨p,E_x|(E₂−B₁)² − (E₁−B₂)²|p,E_x⟩" (r1a/jm_eq29_crop.png), while its own Eqs. (20), (21) and (25), and SS90 (3.1), give α_s/16π; on the printed Eq. (29) the bag moment would be −0.023. Its Eq. (31) also has a = g_x̂ − g_ŷ against SS90 (2.11)'s −2[g_x̂ − g_ŷ]. The paragraph was written before R1 existed and must carry this as open. Its list of sources that derive the term also omits MWZ13 Eq. (23), R1's third independent derivation.
8. L10, and the design's site description "scale = peak Δ/F1", keep "`scale` is the peak Delta/F1" (sf.hpp:428-429). toy_delta_gluon = scale·F₁·x^0.3(1−x)⁴ peaks at Δ/F₁ = 0.337·scale (x = 0.0698), and the t5 harness measures scale as the prefactor. The line D4 rewrites should not repeat the error.
9. O3 gives the halved letter significances as "10.5–21.5σ" (toy) and "9–14σ" (grid). Per bin (R1 §4.2): 42.9/30.6/34.5/21.3σ → 21.5/15.3/17.3/10.7σ, and grid 27.7/18.4/27.1/19.4σ → 13.9/9.2/13.6/9.7σ, so the ranges are 10.7–21.5σ and 9.2–13.9σ. C2 = 3 (≥5σ) still holds.
10. N1 specifies the closed form "at m ∈ {+1, 0, −1} … for LI6 (delta_func) and LI7 (delta_32_func)". For J = 3/2 the projections are ±3/2 and ±1/2. amplitudes() checks only state.j against the ion spin, so the LI7 half would silently exercise the unphysical value c_eff = 3m² − 15/4.
11. The §1 PV row is marked 'read' for "F_UTTT = −Δ", and §7 and L7 cite Eq. (27b) for it. The equations are printed (r1b/arxiv_2610.09862.txt:551, 561, 578, 850), but −Δ needs (A8a) "∆ = W++,−− = (1/ν)(⅓b2 − b3)" in PV's Edelmann–Piller–Weise basis, and identifying that Δ with A₊₋,₋₊ is an inference (R1: 'inferred, medium confidence'). It also reproduces the post-D4 Δ line only because TENSOR_LL_SIGN = −1: the γ-path scale is −TENSOR_LL_SIGN/den.

### Underrated risks (9)

1. Line shifts from P1 (measured): one strict-gate break that --fix cannot repair (SPIN32:986), and five R6 breaks in run_2026-09-03 records, two of them anchor-read pins. The point citation plan_D7_backend_defaults.md:231 python/bindings.cpp:4177 would also drift. Binding on the C_BAG line, or just before the module's closing brace (cited lines end at 5107), costs nothing and removes all of it.
2. Gate mis-specification: T4 will not fail at step 1, and test_bench_t5_delta_bound.py:435 will fail at step 3 without being listed. A strict 'stop on deviation' implementer trips twice on non-physics causes. Meanwhile §10.6's adversarial greps are too narrow to find the surviving old-coefficient prose (t5 header:61 spelled -((1-y)/y^2); test_tensor_gamma.cpp:307 is a comment in tests).
3. Provenance of generated samples: post-D4 tilted-plan files differ from pre-D4 ones in φ, e′, X, the struck species, the target-hook draws and rc_tensor_*. With --hadronize, the bridge's per-event draws and PYTHIA's sequential engine (seeded once) diverge from the first tilted event on. There is no version bump or meta key, and L13's banner text is only 'recommended'. No test pins the generated φ modulation of a tilted plan independently: the closures' truth is a_cos2phi, which moves with the library. Make the meta/banner trace mandatory and pin it.
4. The 'only a human re-dumps xsec.json' safety does not exist: nothing in dump_polligen_reference.py refuses. A re-dump against a partially updated sibling rewrites all seven reference files. The c_ref logic would fail on inconsistent Δ coefficients (good) but accepts every other moved number as the new truth.
5. A factor applied twice across repositories. After D4, external Δ numbers must enter as printed. If the author also follows R1 §4.2's 'keep' recipe (C_BAG/2, 𝒩 = −0.146), or both halves the sibling's coefficient and doubles Δ to keep the letter's printed amplitudes, the bag scenario ends at ¼ or 1× instead of ½× the old amplitude. The convention paragraph should state the one-line rule explicitly.
6. The b-basis framing could lead the author to an O(1) wrong change of A_zz ('subtract Δ/2') for a double counting that does not exist on the default path.
7. Ordering with D7: both D4 and D7 retake test_price_defaults A2_TOY, t5's w_min (0.9493997990 here vs D7's 0.9135) and t5's saturations. The D4 numbers are ToyF2-specific, so whichever commit lands second must re-measure rather than reuse the plan's values.
8. Pre-registration: the design drops V6 from the TIGHT edge (open question 3) but already retakes the Y2/Y3 saturation_tight pins as if the author had answered. That changes a declared edge after the numbers were seen. Either keep TIGHT as declared and relabel it, or wait for the answer.
9. The step-5 'never add a line' rule omits docs/theory/SPIN32_FINITE_GAMMA.md, whose RANGE_ALLOW and UNNAMED_ALLOW keys are doc lines up to 1168 and 1167. Any added D4 sentence there breaks many keys.

### The checker's notes: what it ran and confirmed, and what it did not

No repository file was edited. Scratch work is in /tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/d4adv/ (check1.py/.out, check2.py/.out, gate_base.txt, gate_p1.txt, gate_p1_fix.txt, gate_p1_after.txt, rec_base.txt, rec_p1.txt, rec_inplace.txt, img/ page crops, and the git-archive scratch trees tree/, tree0/, tree2/).

CONFIRMED BY RUNNING on the ca9b23a build (source env.sh):
- Sweet-spot a2 = −1.0559482141520536e−3, −1.2999611414256306e−3, −1.4402307542212716e−3, −1.3264426806399593e−3. Their exact halves are N2's values. a2 equals a_cos2phi bit for bit, and m=0/m=+1 = −2.0.
- Finite-γ, by an independent route (analytic tt channel against the kernel's Δ line): the Δ line equals an F_TTT insertion of −2.000784754 / −2.000008946 / −2.000000001 / −2.007347840 times Δ before D4, and −1.000… after.
- Oversized Δ (1.5F₁ under the old coefficient, equivalent to 3F₁ after D4): the m=+1 margin is +0.277, so T2 fails as claimed. The sampler throws for m=0 at x = 0.6813, Q² = 242.1 with a2/den = 1.449, exactly as stated.
- t5 (b): w_min 0.9493997515410681 → 0.9493997990444439, A_T max 0.05060024845893199 → 0.05060020095555617, A_T min unchanged.
- S2: DS16 over the new tree coefficient = 0.999999603 / 0.999996731 / 0.999999873 at Q² = 1e5, and 0.992137 at Q² = 5.
- Port gate: c_ref = 1 (max |c−1| = 2.2e−16) at every Δ ≠ 0 point of the five Δ-carrying blocks. The 7Li default blocks are Δ-free, so the gate must skip them. The 24/30 and 32/40 state counts are as stated.
- Leakage ×2 and A2_TOY ×½ values are exactly as stated.
- Row 12 stores only ratios in REPORT.json, so it is unchanged. tiers.py prints the six code files UNMAPPED and re-runs t1_r_world and t5.
- The in-place comment edits (xsec.cpp:282, asymmetries.cpp:91, xsec.hpp:21, asymmetries.hpp:10/14/27) leave --records-only at 0 broken.

READ ON PAGE IMAGES / TEXT:
- JM preprint Eq. (16): \"lim … = e⁴ME/(4π²Q⁴){xy²F₁ + (1−y)F₂ − (x(1−y)/2)Δ cos 2φ}\", with E^x₊ = −(1/√2)(0,0,1,i).
- SS90 (A7): \"−½xλΔ₄cos2φ\", with Δ = Δ₄ in the Bjorken limit. SS90 fn 8: \"b₁ = Δ₁ + ½Δ₄ and b₂ = Δ₂ + xΔ₄\".
- LOI18 Eq. (6) (r1a/loi18.txt:138-143) carries \"− x(1−y)/2 Δ cos 2φ\"; Eq. (9) (:193-203) carries \"− (1−y)/y² b_m sin²θ_m Δ(x,Q²) cos(2φ)\".
- HJM (r1a/hjm89.txt:90-91, 209): \"b3 and b4 do not receive contributions at leading twist\"; \"A+−,−+ falls off like 1/Q²\".
- So the cross-section half of the convention paragraph is right. Its parton-side sentence is not (see 'wrong').

The sibling line numbers the design lists check out (asymmetries.py 4-9 and 235-237, xsec.py 17-21 and 732-733, main.tex 116, 199-210 and 232-237, plans/07 52/105/108, test_xsec_identity.py 306/328, test_recopseudo.py 561/576).

A separate D4 gate emulation already in the scratchpad (d4chk/gate_after_fix.txt) shows the same unfixable bindings.cpp:2187-2199 break.
