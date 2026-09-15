# Open items — solutions explored (2026-08-29)

Synthesis of four investigations (full reports in `docs/open_items/`:
`physics_literature.md`, `vmc_overlaps.md`, `code_designs.md`, `engineering.md`,
prototypes in `docs/open_items/prototypes/`). Each item below: what was found,
the recommended solution, effort, and status. Ordered by leverage.

## Run 2026-09-03 — where this file stands

Everything below was reopened, measured and revised by the 2026-09-03 open-items
run (phases A–F). **The run's own records — plan, per-phase measured-number
tables, the status board and the author-decision batch — are in
`docs/open_items/run_2026-09-03/`**: `PLAN.md`, `STATUS.md` (the ONE decision
registry; cite a decision by its ROW NUMBER), `AUTHOR_DECISIONS.md` (every
decision stated in full, with its options and their measured cost),
`SUMMARY.md` (one page: what the run measured, and what it did not establish),
and `phase_{A,B,C,D,E}_*.md` for the raw numbers. Where this file quotes a
number, that is where it was measured.

**Each item's state is one of five**, and the difference matters more than the
word "done" did:

| state | what it means here |
|---|---|
| **closed** | the question is answered and nothing about it is pending |
| **answered as a band** | the answer exists and is a BAND or an interval, never a point; quoting one edge as "the" answer is a misquotation |
| **opt-in shipped** | the code is in the tree behind a flag, the default is unchanged and bit for bit, and the flag's cost is measured |
| **deferred with design** | not implemented, deliberately; the design and the price are written down and the blocker is named |
| **author decision** | the tree will not choose; the question, its options and their measured cost sit in `run_2026-09-03/AUTHOR_DECISIONS.md` under the `STATUS.md` row number given |

An item can carry two: item 10's A = 2 gate is **closed** *and* it hands the
author six decisions (rows 1, 2, 3, 4, 18, 25). **No number in the table below stands without its window** —
the configuration, sample, energy or knob it was measured at.

| # | item | state | the deciding number, with its window | where |
|---|---|---|---|---|
| 1 | VMC α+d / α+t cluster wave functions | **opt-in shipped** (`--cluster-wave vmc`; Hulthén default bit for bit) **+ author decisions rows 9, 10, 11** | ⁶Li α-tag fraction **0.0264 → 0.0348** (×1.315; regenerated 2026-09-06 on the sign-fixed build — the pre-fix sample read 0.0249 → 0.0348, ×1.40) at 10 × 99.5 on the YR high-acceptance envelope, 0.2551 → 0.2486 (regenerated 2026-09-06; pre-fix sample 0.2530 → 0.2485) on the tagging optics; the shipped whole-nucleus ⁶Li polarization is **0.811228** against the ab-initio **0.848**, so **no inclusive ⁶Li polarization may be quoted without the band 0.81 … 0.91** | §1; `phase_C_numbers.md` §§C5.4–C5.5b |
| 2 | ePIC chain gate (HepMC3 → abconv → npsim) | **closed** | **10/10** events through `npsim`, directly and via `abconv -p ip6_hiacc_100x10`, with the writer fix that puts m_e = 0.51099895 MeV in `generated_mass` | §2; `T2_CHAIN.md` |
| 3 | Tensor sign `TENSOR_LL_SIGN` | **author decision — `STATUS.md` row 22** (carried forward, confirm only) | **−1**, the value four independent sources give for A_zz = −(2/3) b₁/F₁; the code has shipped it since 2026-08-29 while two documents still read "author to confirm" | §3–4; `AUTHOR_DECISIONS.md` §B1 |
| 4 | ⁶Li / ⁷Li effective polarizations | **author decision — row 23** | whole-nucleus VMC **0.85 ± 0.03** for ⁶Li (0.848, Wiringa 2014 Table I) against the shipped cluster-picture **0.811228**; ⁷Li 0.866 / −0.037 | §3–4; `AUTHOR_DECISIONS.md` §B6 |
| 5 | Coherent T2 final state | **opt-in shipped** (`--coherent-t2`, `--pom-set`) **+ author decision row 19** (`PomSet` 11 refused) | over all 15 sets × 20 000 events the T0 columns are **bit-identical across the fourteen sets that still run** (one md5), so the DPDF band on M_X, \|t\|, x_P and σ is **identically zero**; the systematic is the hadronic final state — ⟨n_charged⟩ **−2.6 % / +9.9 %** over the twelve DPDF fits, kaon fraction **×2.8** | §5.1–5.2; `phase_D_numbers.md` §D4 |
| 6 | Triton remnant (t* → N + …) | **opt-in shipped** (`--triton-sf ciofi-simula`) | the Ciofi–Simula n₀+n₁ model at **S₀ = 0.6525, untuned**; the sequential Hulthén default is unchanged bit for bit | §6 |
| 7 | Spectator FSI | **opt-in shipped** (`--fsi`) **+ a retraction** | the two Glauber variants are **not one**: **99.50 %** of events differ by more than **1 %** (\|w_nucleon/w_cluster − 1\| > 0.01), ratio to **68.5**, measured on ONE stream — `tagged-6Li-alpha`, 20 000 events, seed 1234, `tensor-thirds` at P_z 0.7 / P_zz 0.6 / P_e 0.7, **σ_XN = 40 mb**, one end of the mandatory 20–40 mb band. Σw/Σw_off **0.520954** vs **0.581605**; the model's own grid-integrated `survival()`, which `meta["fsi_survival"]` carries, **0.520239** vs **0.582899** | §7.1–7.2; `phase_D_numbers.md` §D3 |
| 8 | Spin-3/2 SF basis | **deferred with design** (theory note; no behaviour change) | the rank-≤2 truncation is **exact** for unpolarised-beam inclusive observables — a theorem, not an approximation; the finite-γ J = 3/2 decomposition and the full list of code changes for switching the rank-3 sector on are written down and unimplemented | §8–10; `theory/SPIN32_FINITE_GAMMA.md` |
| 9 | Tensor-sector RC | **answered as a band** (`--rc tensor-band` opt-in; `--rc off` byte-identical) **+ author decisions rows 6, 7, 17** | the tail ships as a **band of two models** and neither edge is "the" tail: POLRAD's "only the t-peak leads" holds **event-weighted (+0.61 %)** in the Q² ≥ 20 GeV², y ≤ 0.9 window and **fails per cell** — **331 of that window's 1356 accepted cells (24.4 %)** differ by > 1 %, worst **×6444** at x = 0.7943, y = 0.0088, 318 of them at y < 0.1. The ⁶Li C0 shape is a band too, in **both** sign and magnitude: σ^el_T/σ^el_U at x = 0.10, Q² = 5 is **+9.357e−04 on the `ho` edge** and **−2.442e−04 on `vmc-ft`** — below POLRAD's own deuteron elastic-tail value (+0.064) by a factor **68** resp. **≈ 260** — so the published SIGN is withdrawn as a band edge | §9; `phase_B_numbers.md` §§B1–B6 |
| 10 | b₁ for A > 2 | **A = 2 gate CLOSED — for a CONFIGURATION, not for a build** (`--b1-model li6-convolution`, opt-in) **+ author decisions rows 1, 2, 3, 4, 18, 25** | G3b peak ratio **0.843243** with MSTW2008 LO at CDKS Eq. (21)'s δ-function — inside the factor-2 window; **0.440 on the shipped `ToyF2` default, OUTSIDE it**; **1.000338** with CD-Bonn as well (a residual below the reference figure's own digitization error, specific to CD-Bonn *and* MSTW); **0.520** at Eq. (17)'s κ = 1. The gate is **A = 2**, so the mandatory **±100 %** band on every ⁶Li number stays | §10 |
| 11 | Coherent ⁶Li amplitude | **answered as a band — MARGINAL** **+ author decisions rows 8, 12, 16, 21** | the ⁶Li tensor a₂ in coherent J/ψ at one EIC year: **S = 2.63 σ at the band's low edge**, the top a **span 2.84 … 3.29 σ** (whether it crosses 3 σ is **not established**), **3 σ at 8.3 … 13.0 fb⁻¹/u** — inside {1, 10, 100} at both ends, quoted on the background-immune ⟨P_zz²⟩ = 0.81, from a **closed-form map and not a dipole-model amplitude**. `Optics::lumi_fraction` alone restores a NO (**0.73 … 0.92 σ**, 106–167 fb⁻¹/u) | §11.3–§11.3c; `phase_C_numbers.md` §§C2, C6 |
| 12 | Packaging | **closed / opt-in shipped** **+ author decisions rows 14, 15** | `pip install -e .` **66 s**; `auditwheel repair` **1 779 767 B → 7 232 971 B (×4.06)** under tag **`manylinux_2_35_x86_64`** (glibc 2.35). With the whole deps prefix replaced by an empty tmpfs, `import lipolgen` and the pure-C++ generator **work** and `LhapdfSF` and `--hadronize` **fail** — those two data trees are not in the wheel | §12–13; `PACKAGING.md`; `phase_E_numbers.md` §E4.6–E4.8 |
| 13 | License | **author decision — row 24** (carried forward; E2 built on it) | **GPL-3.0-or-later**, stamped as `SPDX-License-Identifier` on **101 of 101** files the gate checks; `AUTHORS` and `CITATION.cff` each carry one literal **placeholder** where the name goes — that is decision row 13, unanswered | §12–13; `AUTHOR_DECISIONS.md` §§B21, B22 |
| 14 | Structure-function backend injection (sweep **D2**) | **opt-in shipped** (`--unpol-sf`, `--pol-sf`; both default to the toys and are then bit for bit) | cost of having been on the toy, ⁶Li at config 1: accepted σ **×0.7985** (ct18nlo) / **×0.7934** (mstw), run-level A_zz **×1.2524 / ×1.2604**; and the shipped `ToyG1`'s g₁ⁿ has the **WRONG SIGN over roughly 0.25 < x < 0.6**, which bites hardest on the neutron-tagged `tagged-d-p` channel | §14; `phase_D_sf_injection.md` |
| 15 | ⁷Li rank-2 (tensor) input (sweep **D1**) | **deferred with design**, blocked on **decision row 20** (§15.5 D2) | a ⁷Li inclusive run's tensor term, cos 2φ amplitude and A_zz are **exactly 0** and now say so (banner block, `meta["rank2_input"]`, `b1_model = "none (spin 3/2: no rank-2 input)"`). The α–t convolution is **2.99 ± 0.02 ×** ⁶Li's orbital term **in the research note only, not in the code**, and its **sign flips with the unpolarised backend** the A = 2 gate tells you to use | §15–§15.6; `phase_D_li7_rank2.md` |

### The same fifteen rows in full, as the phases left them

The table below is the original inventory row per item, revised in place by each
phase as it measured. It is longer and it is the one to read for the reasoning;
the board above is the state and the deciding number.

| # | item | verdict | effort | status |
|---|---|---|---|---|
| 1 | VMC α+d / α+t cluster wave functions (plans/04 #15) | **closed & implemented** — ANL AV18 VMC tables (overlaps 2004; momentum distributions 2024) in `data/vmc/`, `VmcRadial` backend, `--cluster-wave vmc` | done | re-run on VMC at 5×41, 10×100, 18×275 (USAGE cluster-wave section; vmc_reconciliation.md "Impact on the tagged pipeline"); β band retired |
| 2 | ePIC chain gate (HepMC3 → abconv → npsim) | **passed** — 10/10 events through `npsim` directly and via `abconv -p ip6_hiacc_100x10` | done | writer fix in (m_e written as generated_mass; pinned by tests) |
| 3 | Tensor sign `TENSOR_LL_SIGN` (plans/08 D1) | **decided by literature: −1** (Cosyn Eq. 27, HERMES Eq. 6, HJM derivation, POLRAD Eqs. 9/10 all give A_zz = −(2/3) b₁/F₁) | 1 line + test | author to confirm; unblocks **plans/08 D2** (not the 2026-09-03 sweep's D2 of row 14, nor §15.5's D2). Status differs by site and that is the open question: `constants.hpp:17` records it as an author decision already TAKEN (−1 since 2026-08-29), this row and `docs/surveys/needs_survey.md` say "to confirm"; carried forward as `run_2026-09-03/AUTHOR_DECISIONS.md` §B1. The `TENSOR_LL_SIGN = +1.0` that survey quotes is the SIBLING `fastsim` tree's, not this one's |
| 4 | ⁶Li/⁷Li effective polarizations (plans/04 #6) | **closed** — 1/3 vs 0.81 was a convention mismatch; VMC (Wiringa 2014 Table I, Piarulli 2023) gives whole-nucleus P_p = P_n = 0.85 ± 0.03 (⁶Li), 0.87 / −0.03 (⁷Li) | 0.5 d | author to confirm — `STATUS.md` decision **row 23**, `AUTHOR_DECISIONS.md` §B6 |
| 5 | Coherent T2 final state | **implemented 2026-09-01** — γ*–Pomeron tier is the default coherent T2 (`CoherentT2::{Pomeron, Off}`, `docs/PYTHIA_BRIDGE.md` §12); ζ = β exact, 300-event chain conserves to 2.8e-14, M_had = M_X to 2.3e-11, veto 0 at M_X ≥ 1.4 GeV | done | `--coherent-t2`, `--pom-set`; **`PDF:PomSet` SCANNED 2026-09-04 (D4)** — all 15 sets × 20 000 events: the T0 columns are **bit-identical across the fourteen sets that still run** (1–10, 12–15; set 11 is refused by the constructor now, so fourteen is the count that reproduces — re-measured 2026-09-05, one md5 `ffd35a3a`), so the band on M_X, \|t\|, x_P and σ is *identically zero* and the systematic is on the **hadronic final state** — ⟨n_charged⟩ **−2.6 % / +9.9 %** and the **kaon fraction a factor 2.8** over the twelve DPDF fits. Set 11 refused (100 % e_q² fallback); the light-only fallback is **measured bit-identical**, not a systematic; `pom_set` now in the npz `meta` — §5.1–5.2 |
| 6 | Triton remnant (t* → N + …) | **implemented 2026-09-01** — `triton_sf.hpp` CS n₀+n₁ model, S₀ = 0.6525 untuned, third channel n → (pn); opt-in `--triton-sf ciofi-simula`, sequential Hulthén stays the default bit for bit | done | numbers in §6 below |
| 7 | Spectator FSI (plans/04 #16) | **implemented 2026-09-01** — `GlauberFsiWeight` per-event weight on `Event::weight` (`--fsi`, tagged channels; never a momentum shift); σ_Xα = 131.0 mb at σ_XN = 40, band 20–40 mb mandatory | done | numbers in §7 below; **D3 (2026-09-04): the inventory's "the two variants are one" is RETRACTED** — measured on ONE event stream (`--channel tagged-6Li-alpha --events 20000 --seed 1234`, `tensor-thirds` at P_z 0.7 / P_zz 0.6 / P_e 0.7, σ_XN = **40 mb**, one end of the mandatory 20–40 mb band), **99.50 %** of events differ by more than **1 %** (\|w_nucleon/w_cluster − 1\| > 0.01; ratio to 68.5, Σw/Σw_off **0.520954** vs **0.581605**; the model's own integrated `survival()`, which `meta["fsi_survival"]` carries, is 0.520239 vs 0.582899) — the whole claim is about that one stream and that one σ_XN (`CONVENTIONS.md` §"quote that window with the two numbers") — and `fsi` / `fsi_sigma_mb` are now in the npz `meta` — §7.1–7.2 |
| 8 | Spin-3/2 SF basis (plans/04 #14) | **implemented as a theory note** — Jaffe–Manohar NPB 321 (1989); explicit J=3/2 functions arXiv:2209.12161 Eqs. 19a–d; rank-≤2 truncation is *exact* for unpolarized-beam inclusive observables; the finite-γ J = 3/2 decomposition and the full list of code changes needed if the rank-3 sector is ever switched on are written up, with no default behaviour change | done | `docs/theory/SPIN32_FINITE_GAMMA.md` |
| 9 | Tensor-sector RC (plans/04 #10) | **implemented 2026-09-03** — `rc.hpp`/`rc.cpp`, opt-in `--rc tensor-band`: the two-sided band `rc_tensor_lo`/`rc_tensor_hi` on the tensor part of the rate (δ log-linear, 0.30 at x = 0.01 → 0.015 at x = 0.16) plus `rc_tail`, POLRAD's t-peak elastic tail with its tensor part and the unpolarised quasi-elastic tail. Weight-only, on `Event::rc_weights` and never on `Event::weight`; `--rc off` is byte-identical. **Phase B (2026-09-04):** the ⁶Li C0 shape ships as a BAND (`--rc-c0-shape ho\|vmc-ft`) and the published sign of (1/6)σ^el_T/σ^el_U at x = 0.10 is withdrawn as a band edge — the MAGNITUDE is a band edge too — σ^el_T/σ^el_U at x = 0.10, Q² = 5 is +9.357e−04 on `ho` against −2.442e−04 on `vmc-ft`, i.e. a factor 68 resp. ~260 below POLRAD's own deuteron value +0.064; `--rc-tail-model t-peak+ll` is the other edge of the tail band, and "< 1 % at Q² ≥ 20 GeV²" holds **event-weighted** (+0.61 %) and **fails per cell** (331 of 1356 = 24.4 %, worst ×6444 at y → 0); `--rc-qe-tensor-scale` (default 0) prices the uncomputed polarised quasi-elastic tail as a BORROWED magnitude; `--rc-a-transfer-frac` (default 0) prices the A = 2 → A = 6 transfer of δ(x); and the five ΔA_zz rows of §8.1c/§8.2 carried a ×3.253983 deuteron-b₁ A_zz — flagged, **not republished** | done | numbers in §9 below |
| 10 | b₁ for A > 2 (plans/04 #9) | **implemented 2026-09-03 as an OPT-IN backend; its A = 2 gate CLOSED the same day** — `b1_nuclear.hpp`/`b1_nuclear.cpp`, `--b1-model li6-convolution`: the **four**-term α–d convolution on the Cosyn–Dong–Kumano–Sargsian kernel (the struck-α orbital term is not optional, it is ≈ 0.5 × the struck-d one). The default stays `Li6B1(MillerB1)`, bit for bit. **The A = 2 validation gate PASSES**: read where design §5.4's checklist ends — MSTW2008 LO (CDKS's own PDF, now reachable as `MstwSF` from PYTHIA's own grid) at CDKS Eq. (21)'s δ-function — G3a passes on all three landmarks and G3b's peak ratio is **0.843243**, inside the factor-2 window; with the real CD-Bonn wave function as well, **1.000338**. Four conditions travel with that: **the lift is about a configuration, not about a build** — the shipped default `ToyF2` gives **0.440**, outside the window, so quote numbers made with `--b1-unpol mstw` and not the default backend's; it is marginal at Eq. (17)'s κ = 1 (0.520); the CD-Bonn agreement is below the reference's own digitization error and is not three-digit agreement with CDKS; and **the gate is A = 2** — so the **mandatory ±100 % band stays**. G3a's own pass is **qualified** by its counting ceiling (§10, "G3a's stated limitation"). The ⁶Li publication ban is lifted. Of the seven close conditions five are done and two are author decisions (Miller's normalisation; the ⁶Li R default) | 10–15 d, ~9 spent | **gate closed**, band mandatory — §10 below |
| 11 | Coherent ⁶Li amplitude (plans/04 #18) | **both halves in-tree** — eSTARlight ⁶Li unpolarized rates/slope (2026-09-02, `estarlight_li6.md`) settle `slope_b = 50 ± 10` GeV⁻² as citable and put a **lower bound** of 1.0e-3 on `f0` (the all-VM 3.0e-2 is not an upper bound — ρ/φ sit below the M_X floor); both recorded in `coherent.hpp`, with the cross sections themselves now in code as `estarlight_li6_coherent()`; the α+d configuration sampler (2026-09-03, `cluster_config.hpp`, `phase_G_numbers.md`) predicts ⁶Li's tensor a₂ ~10× smaller than the deuteron's and of **opposite sign**. **2026-09-04, §11.3: that a₂ is MARGINALLY measurable in coherent J/ψ — a BAND, S = 2.63 σ at the low edge and 2.84 … 3.29 σ at the top, 3 σ at 8.3 … 13.0 fb⁻¹/u, inside the {1, 10, 100} band at both ends; whether the top crosses 3 σ is NOT established (§11.3c)**, once the photoproduction region (Q² < 0.1 = 85 % of the rate, `estarlight_li6.md` §2f, `estarlight_li6_q2_floors()`) and the μ⁺μ⁻ channel are included. a₂ ∝ \|t\| on an e^{−B\|t\|} sample really is only a 0.25 % modulation — but on 4.2e5 events that is 2.6 σ. The photon-polarisation background separates for free and the separation survives Q² → 0 (the P_zz flip is a 1.50× *gain*). **Everything here is a closed-form map, not a dipole-model amplitude** (§11.3, §C2.0). *(§11.3a: the first pass said "NOT measurable — 0.75 σ, 160 fb⁻¹/u"; that used one lepton channel in one Q² window and is withdrawn.)* *(§11.3c: the fourth pass then found the band's TOP resting on an **unquantified factor presented as quantified** — a ⁷Li → ⁶Li efficiency substitution read off a list whose A-ordering is confounded with E/u and Z; read off the four entries that share a beam energy it **straddles 1**, so the top is a span, 2.84 … 3.29 σ, and the crossing of 3 σ is withdrawn.)* *(§11.3b: the second pass had shipped that verdict on an efficiency chain with **two defects of comparable size and opposite sign, neither in its own "complete" assumption list** — ε_det applied at ⁷Li's TOP energy 18 × 117.9 to a ⁶Li sample at 10 × 99.5, worth ×1.12–1.16 UP by the same paper's own ³He energy scan; and **no decay-lepton acceptance or reconstruction efficiency in the chain at all**, DOWN and unbounded below. They cancel to 0.7 %, which is why the point barely moved and why quoting either alone is worse than quoting neither.)* What is still open is not statistics: **no detection efficiency exists below Q² = 0.1 anywhere in this tree**; **no decay-lepton reconstruction efficiency exists in it at all** (only the geometry can be bounded, and it is 0.99 — not the problem); and the far-forward working point is unchosen — **that last is the single correction that on its own restores a NO**, taking the whole band to 0.73–0.92 σ at 106–167 fb⁻¹/u | done / done / **ask stands, on photoproduction** | §11 below — O1/O2/O4/O5 CLOSED/ANSWERED (§§11.2–11.4), O3 partially bounded offline at +18 %/+9 % (§11.5); the Mäntysaari-group ask is **reconciled into one canonical draft** (task C6, `docs/open_items/run_2026-09-03/mantysaari_collaboration_draft.md`), **not sent — sending is the author's decision**; **D5 (2026-09-04): the \|t\| ceiling STAYS 0.2 and its stated reason changed** — the anchor range (\|t\| ≤ 0.30, a property of the input table) is now primary and positivity is stated as secondary and contingent, because the positivity edge is `t_positivity_edge` = 0.245 GeV² only at the shipped `eps_b0` and would be **2.80** at the measured quadrupole; `--coherent-t-max` and the whole coherent block are now in the npz `meta` — §11.6 |
| 12 | Packaging | **implemented 2026-09-02** — `pyproject.toml` (scikit-build-core) in-tree, `pip install -e .` works (66 s); one copy of each `.so` in `lipolgen/`, `$ORIGIN`+deps-prefix RPATH, data/vmc vendored; portable wheel still needs `auditwheel` + GPL-3 terms | done | see §12–13 below |
| 13 | License | **GPL-3.0-or-later** (forced by HepMC3/LHAPDF; matches MCnet norms) | 0 | author to confirm — `STATUS.md` decision **row 24**, `AUTHOR_DECISIONS.md` §B21; E2 stamped the identifier onto 101 files on the strength of it |
| 14 | Structure-function backend injection (the 2026-09-03 sweep's **D2**, filed as `sf-backend-injection`) | **implemented 2026-09-04** — `--unpol-sf {toy,mstw,ct18nlo}` reaches **every kernel the pipeline builds**, the tagged struck-cluster kernel included, which had **no** structure-function slot of any kind before (and which `PipelineConfig::kernel` never reached either). **`--pol-sf {toy,nnpdfpol}` does not have that reach and never could** — it reaches the inclusive and tagged kernels only where the fill also carries `lam_e·P_e ≠ 0`, and not the coherent channel at all, whose rate is spin-independent; under `tensor-thirds`, this CLI's own default plan, it is read on **no** channel. It is *labelled, not credited*, wherever it did not run (`pol_sf_is_read(config, plan)`, §15.3). Both default to the toy backends and are then **bit for bit**, proved per channel by `np.array_equal`. Measured price of the toy on ⁶Li at config 1: accepted σ **×0.7985** (ct18nlo) / **×0.7934** (mstw), run-level A_zz ×1.2524 / ×1.2604 — and the shipped `ToyG1`'s g₁ⁿ has the **WRONG SIGN** over roughly 0.25 < x < 0.6. R and the EMC hook are deliberately **not** covered by *these two* selectors; the EMC hook stays hard-locked everywhere, and R did too until 2026-09-06, when `--r-source` opened it on ONE path — `--b1-model li6-convolution`, where it is the shared hook of registry row 3's option (iii) and is refused everywhere else (§10) | 1 session | done — §14 below |
| 15 | ⁷Li rank-2 (tensor) input (the 2026-09-03 sweep's **D1**) | **the ZERO is now LOUD; the b₁ is DEFERRED, deliberately** — a ⁷Li inclusive run's tensor term, cos 2φ amplitude and A_zz have always been *exactly* 0 (`default_inclusive_kernel` fills a rank-2 slot for spin 1 only) while `meta` recorded `b1_model = "miller"`, a backend that did not run, and no run-surface line said anything. Now: an unconditional banner block, `meta["rank2_input"]`, and `b1_model` / `b1_unpol` = `"none (spin 3/2: no rank-2 input)"`; plus four run-surface defects (F2 the spin-1 plans at J = 3/2, F3 the population-domain message, F4 a three-line segfault, F5 the silently ignored `--pzz` — no longer silent, and since 2026-09-06 no longer ignored either at the opt-in `--pzz-mode typed`, whose cost is measured in `docs/open_items/run_2026-09-06/phase_A_numbers.md` §A3). **No b₁(⁷Li) was implemented**: the α–t convolution is worked out and measured **in `phase_D_li7_rank2.md`, not here** — 2.99 ± 0.02 × ⁶Li's orbital term, and the ⟨r²⟩/P-wave character of the same wave function agrees with Q(⁷Li) to 13 % — but its **sign flips with the unpolarised backend** the A = 2 gate tells you to use, so shipping it would publish a tensor asymmetry whose direction is a flag. **The b₁ numbers are the research note's, quoted, not re-measured and not in the code. The QUADRUPOLE GATE, by contrast, is COMMITTED as of 2026-09-06** (task B3, `run_2026-09-06/phase_B_numbers.md` §B3, D11 paid): its reference is sourced — `LI7_QUADRUPOLE_FM2` = **−4.06 fm²** in `rc.hpp` beside `LI6_QUADRUPOLE_FM2`, from the **same TUNL A = 5, 6, 7 evaluation** (NPA 708 (2002) 3, `Q = −40.6 ± 0.8 mb`) that the ⁶Li constant comes from — and `li7_alpha_t_quadrupole` computes Q = **−3.485059 fm²** in shipped code, reproducing the research note bit for bit, for a **reported ratio of 0.858389 against −4.06 (0.871265 against the −4.00 the note quoted): 13–14 % low**, pinned by doctest T13 and pytest G8. It validates the α–t **wave function's** quadrupole only; **b₁(⁷Li) remains unimplemented** and ⁷Li's rank-2 sector is still exactly zero, which G8 re-asserts in the very test that runs the gate (§15.4, D11) | loud zero 1 session; b₁ ~1 week after the decision | **blocked on author decision D2** (one unpolarised backend, for both isotopes) — §15 below |

## 1. Cluster wave functions — the biggest physics correction — **OPT-IN SHIPPED** (rows 9, 10, 11 open)

The genuine two-cluster VMC overlaps live on ANL's *older* page (`overlap_old/`:
`li6.ad`, `li7.at`, AV18+UIX, 2004, r- and k-space with MC errors); the
`momenta/` page (2024, AV18+UX) has α–d and α–t relative-momentum distributions
with an S/D block split. Both were fetched via the Wayback mirror (the ANL site is
behind a Cloudflare challenge). Findings:

- ⁶Li α+d: P_D = 1.9–2.0 % (the scenario value 0.0867 is 4× too high — it was
  tuned to reproduce a vector dilution of 0.87 that is itself a convention error,
  see item 4). S_αd = 0.82.
- The k-space density has a Pauli node (α–d 2S relative state) at k ≈ 0.15–0.2 GeV
  that no nodeless Hulthén form carries; VMC has *more* strength than Hulthén
  β = 0.30 in the 0.2–0.35 GeV Roman-Pot window and much less above 0.45 GeV.
- Reconciled (`validation/vmc_reconcile.py`, `docs/open_items/vmc_reconciliation.md`):
  the "β band biased low" claim in `physics_literature.md` is withdrawn — it evaluated
  the S-wave Hulthén for the P-wave ⁷Li channel. Against the right form VMC α+t is
  softer than every β in the far tail; ⁶Li VMC is harder at 0.2–0.3 GeV (P(k>0.2)
  0.148 → 0.241) and softer above 0.45. Both ANL file families agree bin-for-bin;
  production input = `momenta/` magnitudes + `overlap_old/` S–D sign.
- **Implemented**: `VmcRadial` backend, `ClusterWaveSource::VmcAV18`
  (`--cluster-wave vmc`); default stays Hulthén. Measured at 10×99.5: ⁶Li α tag
  0.0264 → **0.0348** (YR high-acceptance, ×1.315), 0.255 → 0.249 (tagging optics)
  — regenerated 2026-09-06 on the fixed S–D sign; the Hulthén samples read
  0.0249 / 0.253 and the ratio ×1.40 before;
  ⁷Li 0.973 → 0.998. Tagged A_zz^tag(k) shrinks with P_D (0.087 → 0.019 plus
  shape); the S–D interference sign flips below the S node at 0.134 GeV, currently
  outside every Roman-Pot acceptance.
  **Corrected 2026-09-06** — this line read *"roughly halves"*, which was measured
  against the pre-fix curve. The tagged sector's own S–D interference sign was
  inverted until 2026-09-06 (`src/core/tagged.cpp` `build_amp2` summed ψ_L where
  the partial-wave amplitude needs φ_L = i^L ψ_L, so the relative phase
  `(-1)^floor(L/2)` — `+1` on L = 0, **`−1` on L = 2** — was missing;
  `docs/benchmarking/07_cw_sign_investigation.md`). Re-measured on the fixed
  library at the same optics, the VMC/Hulthén ratio of A_zz^tag runs
  **0.430 / 0.273 / 0.209 / 0.113 / −0.200** at k = 0.1979 / 0.2495 / 0.3012 /
  0.4001 / 0.4990 GeV (it read 0.535 / 0.398 / 0.348 / 0.234 / −0.554 before), so
  the shrink is a factor **2.3 at k = 0.20 falling to 8.8 at k = 0.40** — not a
  halving. **Both curves are now negative** across the accepted window; the tag
  ⟨k⟩ and P_D on this line are spin-blind and did not move; the tag fractions are NOT spin-blind for a tensor-polarised fill (they move with A_zz^tag, −24.5 % at the CLI's default fill) and their 40 000-event rows were regenerated on the fixed build — only the uniform-M mix is invariant.
- Consequence for the physics case: the ⁶Li α-tag acceptance is entirely a
  p_T-tail measurement, so the published tag fractions and the tagged A_zz
  curves were re-run on VMC (`validation/vmc_tag_fractions.py`,
  `docs/open_items/vmc_reconciliation.md` "Impact on the tagged pipeline",
  `docs/USAGE.md`): ⁶Li tag fraction 0.0264 → 0.0348 (YR high-acceptance),
  0.2551 → 0.2486 (tagging optics; the ⁶Li Hulthén samples were regenerated
  2026-09-06 on the fixed S–D sign and read 0.0249 / 0.2530 before, a re-drawn
  category average); ⁷Li 0.9730 → 0.9981; ⁶Li A_zz^tag at
  k = 0.20 GeV **−1.207 → −0.519** (**corrected 2026-09-06**; it read
  `+0.845 → +0.452`, the pre-fix pair, until the tagged sector's S–D
  interference sign was fixed — one `(-1)^floor(L/2)` in `src/core/tagged.cpp`
  `build_amp2`, the observable relative part of φ_L = i^L ψ_L, which was
  missing entirely; `docs/benchmarking/07_cw_sign_investigation.md`,
  `docs/open_items/run_2026-09-06/phase_CW_numbers.md`). Measured at the
  k = 0.1979 GeV cell — the grid point nearest 0.20 — on
  `azz_tensor_curve_weighted` at `default_configs("6Li")[1]`,
  `yr_optics(..., high_acceptance = true)`, `n_phi = 32`: **−1.2069** Hulthén
  β = 0.30 and **−0.5191** VMC AV18, against +0.8450 and +0.4518 before. **The
  Hulthén → VMC statement this bullet makes survives the fix and every other
  number on this line is unchanged except the tag-fraction rows, which were regenerated**: P_D is spin-blind,
  and the spin-blind accepted rate is **unmoved across the fix to ≤ 4.0e−16
  relative (1–2 ulp; the VMC channel is bit-identical)** — measured as the
  acceptance-weighted Σ_M n_M k² at `n_phi = 32`, giving an accepted fraction
  of **0.0246759321488** (Hulthén β = 0.30) and **0.0338102276258** (VMC
  AV18). Digits beyond those are summation-order dependent, not physics.
  P_D 0.0867 → 0.01935. Because no single β reproduces the VMC shape (Pauli node + window-dependent tail), the Hulthén
  β band is retired as the cluster-wave systematic rather than widened — the
  `--cluster-beta` knob itself stays (Hulthén is still the default radial
  form), but the quoted systematic is now the VMC-vs-Hulthén difference
  shown above, not a scan over β. Confirmed at the other two reference
  configurations (`validation/vmc_tag_fractions.py --configs 0,1,2`,
  `docs/USAGE.md`): the ⁶Li YR high-acceptance tag fraction moves
  0.0301 → 0.0365 at 5×41 and 0.0279 → 0.0348 at 18×275 (regenerated
  2026-09-06 on the fixed S–D sign; 0.0286 → 0.0365 and 0.0266 → 0.0349
  before), so the Hulthén→VMC shift is not a single-energy artifact.

## 2. Chain gate — **CLOSED** (passed)

`lipolgen-run … --hadronize` → `npsim --compactFile epic_craterlake_10x100.xml`
accepts the file directly (10-digit ion codes for beam ⁶Li status 4 and the α
spectator status 1); `abconv -p 1` cannot decode a ⁶Li ion (per-nucleon energy
0), `abconv -p ip6_hiacc_100x10` works and its output also runs through npsim.
EDM4hep `MCParticles` carry the full role chain. Cosmetic: write m_e = 0.511 MeV
as the electron `generated_mass` to silence DD4hep's ppm energy fix-ups.

**Done (2026-09-02)**: `src/hepmc/hepmc_writer.cpp` writes
`generated_mass = 0.51099895e-3` GeV for any massless electron (`|pdg| == 11`,
`p.mass == 0.0`); every other particle's generated mass is unchanged. Pinned
by a test in the doctest suite.

## 3–4. Conventions now decided by sources — **AUTHOR DECISIONS** (rows 22, 23)

- `TENSOR_LL_SIGN`: four independent conventions (Cosyn 2410.12764 Eq. 27, HERMES
  hep-ex/0506018 Eq. 6, a two-line HJM parton-model derivation, POLRAD 2.0 Eqs.
  9/10) all give A_zz = −(2/3) b₁/F₁; the program's +1 is opposite, exactly as its
  own test `test_the_program_sign_is_opposite_to_the_literature` says. Flipping
  the constant also flips κ and the O(γ²) subtraction (D2).
- ⁶Li effective polarization: whole-nucleus convention (Cloët–Bentz–Thomas Eq. 24,
  no ÷Z): VMC 0.848 (Wiringa 2014 Table I), 0.86 chiral (Piarulli 2023);
  adopt 0.85 ± 0.03 (per-nucleon 0.283 if the slot divides by Z/N — state which).
  ⁷Li 0.866 / −0.037 confirmed. Do **not** derive tensor inputs from these wave
  functions (Q(⁶Li) is off by 2.5×).

## 5. Coherent T2 — Pomeron beam through the existing bridge — **OPT-IN SHIPPED** (row 19 open)

`Diffraction:doHard` cannot be driven externally and never fires for virtual
photons, but `Beams:idA = 990` is a legal LHAup beam (`BeamSetup.cc:870`). With
P_IP = P_ion − P_recoil (massless), ζ = β = Q²/(M_X²+Q²) exactly; the Pomeron
remnant is a single antiquark, so charge and colour close automatically. Traps
found: colour-tag orientation for an incoming antiquark (else 50 % silent veto);
LO Pomeron grids have no quarks at Q² = 1, β < 0.1 (need the `q2_pdf_min` floor).
Raise `COHERENT_MX_MIN_DEFAULT` to 1.2 GeV; below it, an exclusive-VM channel
(ρ, φ, J/ψ) is a separate Phase-2 item. Fallback: hand-filled qq̄ string with
`ProcessLevel:all = off` (130k ev/s, 1e-15 conservation).

**IMPLEMENTED (2026-09-01), measured.** The design above is the shipped tier
(`docs/PYTHIA_BRIDGE.md` §12): `make_coherent` writes `Role::Pomeron`
(P_IP = P_ion − P_recoil, pdg 990, status 3), the bridge's third instance
(`Beams:idA = 990`, behind `PythiaBridgeOptions::coherent_t2`) hadronizes it
with W² → M_X², ζ = β verified exact to 10 digits, flavour drawn
e_q² × f_q(β, Q²) from `PDF:PomSet` (default 6) with the e_q² fallback. Both
predicted traps were hit and closed (antiquark colour-tag orientation;
quark-free LO grids → `q2_pdf_min` floor). 300-event chain: conservation
2.8e-14 relative, charge exact, recoil untouched, M_had = M_X to 2.3e-11,
veto 0 at M_X ≥ 1.4 GeV (the measured veto table sits in `coherent.hpp` and
set `COHERENT_MX_MIN_DEFAULT = 1.2`). One extra find: PYTHIA closes the
meson-like Pomeron-beam record on the electron beam's massive light-cone
minus — a constant Δ(m²) = −3.9e-6 GeV² the baryon-beam remnant path absorbs —
so the surrogate is built at a compensated w2_sur, keeping the mass-repair
rescale at ~1e-11.

### 5.1 `PDF:PomSet` — the tier's largest systematic, SCANNED (2026-09-04, D4)

`PDF:PomSet` was called the coherent T2 tier's largest systematic from the day
the tier was built and had never been scanned. It has been now: **every set
PYTHIA 8.317 ships (1–15), 20 000 coherent events each, ⁶Li config 1, seed
4242, `coherent_t2 = Pomeron`, single thread — 17 s wall for the whole scan**,
as run 2026-09-04; set 11 has been refused by the bridge constructor since.
Reproduction and the full per-set table are in
`docs/open_items/run_2026-09-03/phase_D_numbers.md` §D4.

**The observable had to move, and the reason is exact.** The T0 columns `t`,
`x_pom`, `q2`, `x` and the event `weight` come out **bit-identical across the
fourteen sets the constructor still admits** (1–10, 12–15) — one md5,
`ffd35a3a62b591c547e9ca2ac4301b5d`, re-measured 2026-09-05 in 14.3 s. Set 11
gave that md5 too in the 2026-09-04 scan, but it is refused now (below), so
**the reproducible claim is over fourteen, not fifteen**; this section said
"all fifteen" until 2026-09-05. The Pomeron PDF enters only the flavour
draw and PYTHIA's backward evolution, while |t|, x_P, M_X and the rate are
fixed upstream by `CoherentSampler` / `CoherentXpomModel`. So a `PomSet` band
on M_X (the observable this item originally proposed), on |t|, on x_P or on σ
is **identically zero by construction**, and quoting one would be meaningless.
The band is on the **hadronic final state**, which is the only thing that
moves.

**The band, measured over the twelve genuine diffractive-PDF fits (3–10,
12–15), about the default set 6:**

| observable | set 6 | min | max | band about the default |
|---|---|---|---|---|
| **⟨n_charged⟩** (primary) | 3.964 ± 0.016 | 3.860 (set 9) | 4.355 (set 5) | **−2.6 % / +9.9 %** |
| ⟨n_hadrons⟩ | 8.521 | 8.300 (10) | 9.394 (5) | −2.6 % / +10.2 % |
| ⟨p_T⟩ per particle [GeV] | 0.342 | 0.322 (5) | 0.379 (10) | −5.9 % / +10.9 % |
| **kaon fraction of the HFS** (secondary) | 0.0646 | 0.0359 (9) | 0.1002 (10) | **−44 % / +55 %, a factor 2.8** |

Statistical error on ⟨n_charged⟩ is 0.017 at 20 000 events and the
seed-to-seed scatter over 4242 / 777 / 31337 is ≤ 0.057, against a 0.50 spread
across sets — so 20 000 events per set already resolves the band by a factor
~30 and there is no case for more. The band itself is stable seed to seed:
−1.8…−2.6 % / +9.3…+9.9 % on ⟨n_charged⟩, factor 2.66–2.79 on the kaon
fraction, with set 9 always the low end and sets 5 / 10 always the high ones.

**How to quote it.** As an **envelope over re-runs, one npz per set**, never as
a per-event reweighting: the set changes the final state event by event and
there is no weight that maps one set onto another. That is why `meta["pom_set"]`
had to exist first (§5.2). And say which is which: **set 6 is the only LO H1
set**, so the band mixes LO and NLO DPDFs used in an LO Monte Carlo —
defensible for a *systematic envelope*, indefensible for a *central value*.

**What is excluded, and why — and what it costs.** Sets 1 (a Q²-independent
`N x^a (1−x)^b` toy) and 2 (π⁰ distributions) are not Pomeron fits and are
outside the band; set 1 is still worth running as a sanity floor and is in the
reproduction table. **Measured, so the exclusion is not a free choice**
(re-measured 2026-09-05, same recipe): including them moves the ⟨n_charged⟩
low edge from −2.6 % (set 9, 3.8600) to **−4.7 %** (set 2, 3.7773) and the
⟨n_hadrons⟩ low edge from −2.6 % to **−7.0 %** (set 2, 7.9261); the high edges
(set 5) and the ⟨p_T⟩ and kaon bands are unchanged, because sets 1 and 2 are
interior on those. So **−2.6 % / +9.9 % is a twelve-fit number and must never
be attached to "all 15 sets"**; **four** sites said "all fifteen" beside it
(`docs/DEVELOPMENT_PLAN.md`, `python/bindings.cpp`,
`docs/open_items/run_2026-09-03/PLAN.md` D4, all three corrected in the first
pass on 2026-09-05, and `docs/open_items/run_2026-09-03/STATUS.md` row D,
which read "all 15 sets × 20 000 events; … ⟨n_ch⟩ −2.6 %/+9.9 %" on one line
and was only found on the second pass the same day).
**Set 11 is now REFUSED by `PythiaBridge`'s constructor** when the Pomeron
instance is built: `PomHISASD` returns densities only after `setXPom(x_Pom)`,
which this bridge never calls (there is no Angantyr collision system here to
supply x_P), and measured, **100.00 % of 20 000 events took the charge-
democratic e_q² fallback** — the Pomeron PDF was not consulted on a single
event while `meta["pom_set"]` would have said 11. That is exactly the case
`PipelineConfig::validate`'s "a knob that did not run may not be recorded as
if it had" rule exists to prevent.

**Charm is a feature of the band and must be flagged.** Measured by toggling
`include_charm` and diffing the final states event by event (4 000 events,
seed 4242): set 5 (H1 2007 Jets) changes **29.98 %** of its events, the four
GKG18 sets 10.4–11.1 %, and sets 3, 4, 6, 7, 8, 9, 10 **exactly 0.00 %**. So
including set 5 in a scan produces an open-charm final state the default never
makes. Confirmed at the PYTHIA source, not merely observed: sets 3, 4 **and** 6
are all `PomH1FitAB`, whose `xfUpdate` sets `xc = xcbar = xb = xbbar = 0.`
unconditionally (`PartonDistributions.cc:2630`) — so the tree's old wording,
"the H1 **LO** grids carry no charm or bottom", was true but too narrow: it is
equally true of the two H1 **NLO** fits.

**The light-only e_q² fallback is NOT a systematic in its own right —
measured, it costs exactly zero.** Raising `q2_pdf_min` from 1.0 to 1.75 takes
the fallback share from 20.70 % → 0.00 % (set 6), 31.25 % → 0.00 % (3),
35.52 % → 29.32 % (4), 19.45 % → 2.65 % (12), 20.10 % → 3.23 % (13),
21.65 % → 5.12 % (15) — and leaves a **bit-identical** final state (`pid` and
`p4` arrays, ~34 000 hadrons over 4 000 events) in every one of those cases.
The reason is structural, not luck: every Pomeron DPDF PYTHIA ships carries a
single light-quark singlet (H1 sets it explicitly; the GKG18 LHAGrid1 grids
have columns −3…3 equal row by row), so e_q²·xf_q ∝ e_q² **exactly** over the
light flavours and the normalised "fallback" *is* the true draw. The fallback
fires only where every weight vanishes, which is below the charm threshold,
where charm is zero anyway.

*Corrected while measuring it:* `docs/PHYSICS_CHANNELS.md`'s remedy —
"raising `q2_pdf_min` to ≈1.75 removes it at the cost of clamping every
flavour weight to that Q²" — was true about the counter and **false about the
cost**; at 1.75 it changes nothing at all. There IS a cost at 3.0, and it is a
**charm** effect rather than a fallback one: sets 12, 13, 15 stop being
bit-identical there because the clamp lifts charm above threshold, and with
`include_charm = false` the same comparison is bit-identical again. Keep the
light-only restriction (it is the correct guard against the `2e404b8` bug) and
stop implying the fallback share is a modelling uncertainty; it is a
bookkeeping artefact.

Pinned at reduced statistics in `python/tests/test_pom_set_band.py` (4 sets ×
4 000 events, 4 s) so the band cannot rot silently, and carried at the run
surface by `--pom-set`'s help string and the coherent run banner.

### 5.2 The knob was absent from the npz `meta` (2026-09-04)

Two runs differing only in `--pom-set` had **identical `meta`** while their
entire hadronic final state differed. `python/bindings.cpp` states the rule on
the `b1_*` block — *without these keys three otherwise identical npz files are
indistinguishable, which is exactly the "never quote a single row" rule failing
silently* — and it had been applied to `b1_*` and to `rc_*` and to nothing
else. A conditional T2 block now records `t2_bridge`, `coherent_t2`,
`pom_set`, `pom_rescale`, `pom_q2_pdf_min`, `t2_include_charm`, `t2_n_ok`,
`t2_n_failed`, `t2_n_pomeron`, `n_pom_flavour_fallback` and
`pom_flavour_fallback_frac`. It is emitted only when the T2 tier was bound
through `set_pythia_hadronizer`, so a T0 npz keeps exactly today's key set.

The mechanics are worth recording: `PipelineConfig::hadronizer` is a type-
erased `std::function`, so the metadata writer could not see the bridge at all.
`set_pythia_hadronizer` now installs a **named** callable
(`PythiaHadronizerHook`) and `columns_to_dict` recovers the bridge with
`std::function::target`. `n_pom_flavour_fallback` was likewise counted and
surfaced nowhere; it is now in the `meta` and printed at the run banner beside
`n_ok` / `n_failed`.

## 6. Triton remnant — Ciofi–Simula A = 3 — **OPT-IN SHIPPED**

BeAGLE's `DT_KFERMI` carries n₀(k) = Σᵢ Aᵢ e^{−Bᵢk²}/(1+Cᵢk²)² with the CS
coefficients (A=3: 31.7/1.32/5.98 + 0.00266/0.365/0); its norm ∫n₀k²dk = 0.653 is
the ground-state-remnant (two-body) fraction — the split the breakup needs.
Design: `TritonSpectralFunction` backend with k-dependent branching
p₂(k) = n₀/(n₀+n₁), a third channel (struck n → (pn) continuum, absent today),
fixed uniform consumption per draw. Also recorded: BeAGLE renormalizes n₀ to 1 and
therefore drops the 35 % continuum for A = 3 (its tail is too soft by construction).

**IMPLEMENTED (2026-09-01), measured** (`triton_sf.hpp` / `triton_sf.cpp`,
opt-in via `BreakupOptions::triton_sf` / `PipelineConfig::triton_sf` /
`--triton-sf ciofi-simula`; the sequential Hulthén model stays the default
bit for bit). S₀ = 0.652548 untuned against the CS 0.6525; n₀ tail
P(k > 0.1/0.2/0.3/0.45 GeV) = 0.4267 / 0.06276 / 0.009621 / 0.002365 (all at
the transcription's own digits); sampled bound fraction 0.6524 at 10⁵ draws;
⟨k⟩ = 102 MeV on the bound-d channel, 126 MeV over all struck nucleons
(sequential model: 133/145). The (pn) continuum pair splits at the pn ¹S₀
pole `KAPPA_PN_SINGLET` beside the nn one; sample() consumes exactly 6
uniforms whatever the branch and the breakup exactly 9, verified empirically.
The Pipeline builds the model itself at the run's own `cluster_beta`, so no
physics number is defined twice. T1+T2 chain with the hadronizer on
conserves < 1e-9 on both options (`tests/test_t2.cpp`).

## 7. FSI — weight, not a shift — **OPT-IN SHIPPED**, and the "two variants are one" premise retracted

Eikonal rescattering transfers transverse momentum only, so the distortion is a
function of k_T at fixed k_z; pole dominance = small k_T, which is exactly the
Roman-Pot α-tag selection. Prototype (σ_XN = 40 mb): ≤ 20 % for k ≲ 0.02 GeV,
~27 % at 0.05, 40–60 % at 0.10; σ_XN(W) formation-length ramp is the weakest
input and must be banded (the 20 mb row is arguably realistic at EIC). Use σ_tot
in the absorptive term and σ_el in the gain term (74 % of α rescatterings destroy
the tag). Hook: `FsiWeight` interface on `TaggedSampler`, multiplied into
`Event::weight`; null = today's PWIA bit-for-bit.

**IMPLEMENTED (2026-09-01), measured** (`fsi.hpp` / `fsi.cpp`,
`TaggedSampler::set_fsi` → `TaggedEvent::weight` → `Event::weight`;
`PipelineConfig::fsi` / `--fsi {off,glauber-cluster,glauber-nucleon}` +
`--fsi-sigma-mb`, tagged channels only, `validate()` refuses it elsewhere).
Profile numbers at σ_XN = 40 mb: σ_Xα = 131.0 mb (Glauber-shadowed, not
4 × 40), σ_el = 35.2 mb, B_α = 27.2 GeV⁻²; at 20 mb: 72.5 mb. The pinned
FSI/IA table reproduces the `fsi_alpha.py` prototype to ≤ 3.3e-3 (0.787 at
k = 0.02; 0.607/0.381 at k = 0.10, θ = 0/90°; 0.398 at 0.20; the 20 mb row
0.877/0.844/0.764), production S+D channel 0.617 at k = 0.10, θ = 0;
survival 0.517 (cluster) / 0.582 (nucleon variant) on the ⁶Li α tag at
40 mb. Exact θ → π−θ symmetry, ratio → 1 as σ → 0, the wrong-spectator
guard, the ⁷Li P-wave channel and the σ_XN(W) formation ramp are all in
`tests/test_fsi.cpp`; the weight reaches HepMC3 `weights()[0]` and the npz
`weight` column, and every four-vector stays bit-identical to the PWIA run.
Quote it as an unpolarized-shape systematic banded over 20–40 mb, never as
a correction to A_zz.

### 7.1 The two variants are NOT one — D3, premise false (2026-09-04)

The inventory carried, as open item D3, *"the per-nucleon Glauber FSI variant
is algebraically identical to the cluster one — the two 'variants' are one."*
**Retracted: measured, 99.50 % of events differ by more than 1 %
(\|w_nucleon/w_cluster − 1\| > 0.01), and by construction rather than by
parameter choice.**

What is true is narrower and is what `fsi.hpp` already said: for an
**uncorrelated** cluster density the Ciofi degli Atti–Kaptari per-nucleon
product ⟨Π_i[1 − Γ_N(b − s_i)]⟩ factorises into [1 − (Γ_N ⊛ T_a)]^A, which is
the **cluster** form — so a literal per-nucleon code path would reproduce
variant (a), and that is precisely why the shipped variant (b) is deliberately
something else: the **single-scattering (optical, unshadowed) limit**
Γ_a = A·(Γ_N ⊛ T_a), i.e. σ_Xa = A σ_XN exactly.

**Measured on one event stream reweighted three ways** (`--channel
tagged-6Li-alpha --events 20000 --seed 1234 --fsi {off,glauber-cluster,
glauber-nucleon}`, plan `tensor-thirds` at P_z = 0.7 / P_zz = 0.6 / P_e = 0.7,
at the default σ_XN = 40 mb — one stream at one end of the mandatory 20–40 mb
band; the columns `k`, `cos_theta_k`, `phi_k`, `x`, `q2` are bit-identical
across the three files, so the weight is the whole difference):

| | `off` | `glauber-cluster` | `glauber-nucleon` |
|---|---|---|---|
| Σ weight (20 000 events) | 20000.0 | 10419.07 | 11632.10 |
| integrated survival Σw/Σw_off | 1 | **0.520954** | **0.581605** |
| per-event weight, min … max | 1 | 0.0208 … 1.4043 | 0.2139 … 6.9523 |
| σ_tot(X–α) [mb] | — | 131.045 | 160.006 |
| σ_el(X–α) [mb] | — | 35.224 | 68.186 |
| `survival()` | — | 0.520239 | 0.582899 |

Per-event ratio w_nucleon/w_cluster, percentiles [1, 5, 25, 50, 75, 95, 99] =
0.821, 0.827, 0.853, 0.897, 0.934, 3.337, 5.935; **maximum 68.52**;
`np.allclose(rtol = 1e-14)` is **False**; **99.50 %** of the 20 000 events
differ by more than 1 %. Pinned in `tests/test_fsi.cpp` ("the two variants
differ on essentially every event") and in
`python/tests/test_meta_provenance.py`.

**What a genuinely different per-nucleon variant would need**, and why it is
not a coding task. `fsi.hpp`'s standing TODO names both missing pieces and
they are the right two — the centre-of-mass constraint Σ_i **s**_i = 0, and
short-range NN correlations in the cluster density. Of the two:

* the **c.m. constraint** is a change of T_a alone, i.e. of
  `cluster_point_a2_fm2` (Gartenhaus–Schwartz on a Gaussian gives
  a²_int = a²(1 − 1/A) = 0.5250 fm² for the α). A better input to the **same**
  variant, not a different variant;
* **short-range correlations** need the **two-body** density ρ₂(**s**_i,
  **s**_j) — that is the whole content of ⟨Π_i(·)⟩ ≠ Π_i⟨·⟩ — and **the tree
  has none**. `data/vmc/density/{he4,li6}.density` are ONE-body point-proton
  densities, which is exactly the T_a the cluster form already convolves, and
  the ANL page the README fetches from publishes one-body densities only. So
  the SRC correction needs an input the repository does not have.

And the input is **not** σ_NN: the rescattering here is the DIS debris X off a
**spectator nucleon**, so the amplitude is X–N and the cross section is
`GlauberFsiOptions::sigma_xn_mb`. A σ_NN would enter only for the spectator
cluster's own *internal* absorption, which is not what this weight is.

### 7.2 The knob was absent from the npz `meta` (2026-09-04)

The three runs above produced **identical `meta` dicts**, key by key, while
their total rate differed by 48 % / 42 %. A conditional FSI block now records
`fsi`, `fsi_sigma_mb`, `fsi_sigma_cluster_mb`, `fsi_sigma_cluster_el_mb`,
`fsi_survival`, `fsi_clipped_grid_fraction` and `fsi_formation_ramp` (plus the
four ramp anchors when the ramp ran). It is emitted only when the run had an
FSI weight, so an `--fsi off` npz keeps exactly today's key set. This matters
most for `fsi_sigma_mb`: it carries a documented 20–40 mb band that the header
says must never be quoted as a single row, and until now the two ends of that
band produced indistinguishable files.

## 8–10. Theory notes — what each note became

- Spin-3/2: the complete basis is Jaffe–Manohar (1989); 2209.12161 writes the four
  leading-twist J = 3/2 functions (their "g₂" is the rank-3 partner of g₁ — rename
  `g1_rank3`); Cosyn–Weiss 2603.23699 App. D is the on-ramp. Time reversal makes
  odd-l multipoles unobservable with an unpolarized beam in inclusive DIS, so the
  generator's rank-≤2 truncation is a theorem, not an approximation. Gap: the
  finite-γ inclusive decomposition for J = 3/2 (5–10 d, ideal Cosyn/Weiss co-authorship).
  → **Written up 2026-09-03 as `docs/theory/SPIN32_FINITE_GAMMA.md`**: the
  rank-≤2 truncation is proved exact for unpolarised-beam observables (§4), the
  finite-γ J = 3/2 decomposition is derived there (§5 — two rank-3 structure
  functions, master formula (46)–(47); the repository's own derivation, still
  unpublished in the literature, with the Appendix listing every unsourced
  claim), and §6 lists the code changes if the rank-3 sector is ever switched
  on. **No code behaviour changed.**
- Tensor RC: adopt ISR shift (spin-blind) + POLRAD Eq. (A.4) tensor elastic tail
  with VMC ⁶Li form factors (Wiringa–Schiavilla 1998); quote 1.5 % (x ≳ 0.05) and a
  10–30 % band on the tensor part at x ≲ 0.01 (Gakh–Shekhovtsova, uncited).
  → **Implemented as §9, with three deliberate departures from this note**: the
  ISR *shift* is **not** applied (it would cost an extra uniform per event and
  break the "RC moves nothing" invariant — its residual is what the band
  prices); the tail is built on POLRAD's **t-peak closed forms Eqs. (37)–(39),
  (43)**, not the Eq. (A.4)/Appendix-B machinery, because POLRAD §2.1.3 B says
  the s- and p-peaks are *suppressed* for a tail; and the ⁶Li form-factor
  **normalisations are the measured moments, not VMC** (whose Q(⁶Li) is 3× the
  measured one). The 1.5 % anchor also moved to x = 0.16, E12-13-011's own
  lower kinematic edge.
- b₁(⁶Li): three-term α–d convolution (embedded deuteron b₁ ⊗ f_{d/Li}(z), α–d
  D-wave term with F₁ᵈ, CG depolarization); validate on A = 2 (Cosyn Figs. 4/5)
  first; 100 % band.
  → **Implemented as §10, with one structural correction to this note: it is
  FOUR terms, not three.** CDKS Eq. (10) sums the spectral function over
  *constituents*, and one level up that sum runs over {d, α}. The α is J = 0 so
  b₁^α ≡ 0, but its light-cone density carries the *same* (3cos²θ − 1) orbital
  alignment as the deuteron's — (3c² − 1) is even in k⃗ and the α carries −k⃗ —
  so the **struck-α orbital term (2α) exists and is ≈ 0.5 × the struck-d one
  (2d)**, fixed by counting (4/6 against 2/6) and by the 1/M² of the
  P₂-weighted density: 2(M_d/M_α)² = 0.5064 against a measured **0.494–0.504**
  (0.489–0.497 on the design's coarser y grid — a quadrature artefact, since
  fixed by the 2400/2400/3200 default; `phase_D_numbers.md`).
  A regression that silently drops it moves b₁ by ~30 %. The A = 2 validation
  was done and it **passes on magnitude in one configuration** — MSTW2008 LO
  as the unpolarised nucleon input (`--b1-unpol mstw`) at CDKS Eq. (21)'s
  δ-function, G3b = 0.843243. On the shipped default `ToyF2` the same clause
  is **0.440, outside** the [0.5, 2] window. See §10, which states the
  configuration and the limitation on G3a's counting window with it.

## 9. Tensor-sector RC — a band and a background, never a shift — **ANSWERED AS A BAND**, opt-in (rows 6, 7, 17 open)

**IMPLEMENTED (2026-09-03), measured** (`include/lipolgen/rc.hpp`,
`src/core/rc.cpp`; `PipelineConfig::rc` / `--rc {off,tensor-band}` plus
`--rc-delta-low-x`, `--rc-delta-high-x`, `--rc-fq-scale`,
`--rc-tail-tensor-scale`, `--rc-qe-suppression`; `docs/USAGE.md` §7b; the raw
run in `docs/open_items/run_2026-09-02/phase_C_numbers.md`).

> **Revised 2026-09-03 after review.** Every tail number below was
> regenerated: the elastic tail's per-nucleon reduction was **6× too large**
> (Eq. (38) is the whole-nucleus `d²σ/dx_A dy`, so the factor is `1/A²` and not
> `m_p/M_A`), the quasi-elastic **Pauli suppression `S(q)` is now on by
> default** (design Q9, closed), the nuclear map moved to this library's exact
> `x_A = x/A`, and the band is now **clamped** (`RcOptions::band_tau_max`)
> because `tau_tag` diverges at the nodes of the tagged spectator density and
> was publishing negative weights.

Three weights per event — `rc_tensor_lo = 1 − δ(x)τ`,
`rc_tensor_hi = 1 + δ(x)τ`, `rc_tail = 1 + σ_tail/σ_Born` — on
`Event::rc_weights` and **never** on `Event::weight`: the band is a systematic
variation and the tail a background, so neither may touch the Born sample.
No four-vector moves, no RNG is consumed (`RcModel::fill` takes no `Rng&`),
and an `--rc off` npz and HepMC3 file are **byte-identical** to what the same
seed wrote before `rc.hpp` existed. Applies on every channel except coherent
⁶Li (a φ-dependent tensor observable — nothing to cite); the tail is
identically 1 on the tagged channels, where the tag itself vetoes the elastic
recoil at x_L = 1. **That last clause is half a kinematic fact.** It is one for
the ELASTIC tail; for the QUASI-ELASTIC tail it is an **omission**, worked out
in the 2026-09-04 run: the A−1 remnant is unbound (⁵Li, ⁵He), so its α comes
out at x_L ≈ 2/3 — the tag window — and when the struck nucleon belongs to the
embedded deuteron the α is a **true spectator** carrying the same n_M(k, c) the
tagged Born does, so the tag neither vetoes it nor suppresses it in the
**ratio**. The omitted dilution is of the same order as the inclusive
quasi-elastic one (99.9 % of the inclusive tail at x = 0.30), and the design's
"probably small" is withdrawn. `exclusion_reason()` now says so in the banner
and in `meta` (`run_2026-09-03/phase_B_numbers.md` §B4).

**Configuration for every number below.** ⁶Li, `--config 1`
(s per nucleon = 3980 GeV²), `--channel inclusive`,
`--plan tensor-thirds --pzz 0.6`, θ_S = 0, defaults, `n_eta = 128`, tail node
grid 101 × 77 in (ln x, ln y).

**The band, at P_zz = +1, Q² = 5 GeV²** (`w_lo + w_hi == 2.0` bit-for-bit):

| x | δ(x) | A_zz | τ | `w_hi − 1` |
|---|---|---|---|---|
| 0.010 | 0.30000 | −1.4735e−03 | −7.3728e−04 | −2.2119e−04 |
| 0.063 | 0.11081 | −4.2848e−03 | −2.1470e−03 | **−2.3790e−04** |
| 0.100 | 0.06331 | −4.9751e−03 | −2.4937e−03 | −1.5789e−04 |
| 0.160 | 0.01500 | −4.1666e−03 | −2.0876e−03 | −3.1315e−05 |
| 0.300 | 0.01500 | −1.5852e−04 | −7.9268e−05 | −1.1890e−06 |

The band **peaks near x = 0.063, not at the lowest x**: δ(x) is still rising
there while |A_zz| has not yet fallen. The largest RC systematic on `A_zz`
sits in the *middle* of the low-x range.

> **⚠ The `A_zz` column of that table — and therefore its `τ` and `w_hi − 1`
> columns — carries the ×3.253983 deuteron-b₁ `A_zz`** of the WARNING further
> down (`run_2026-09-03/phase_B_numbers.md` §B3.2), exactly as
> `run_2026-09-02/phase_C_numbers.md` §8.2's does. The `δ(x)` column is
> **unaffected** — it is `rc_delta` alone.
>
> **REPUBLISHED BESIDE IT ON 2026-09-06 (registry row 6 / §B11).** The table
> above stays as published, labelled as what it was computed with; the same
> five points with the **shipped** kernel's ⁶Li b₁ (`Li6B1(MillerB1)`, i.e.
> `InclusiveKernel::tables`' own `b1`/`b2`) are:
>
> | x | δ(x) | `A_zz` | `τ` | `w_hi − 1` | half-width δ·\|A_zz\| |
> |---|---|---|---|---|---|
> | 0.010 | 0.30000 | −4.528242e−04 | −2.264634e−04 | −6.793901e−05 | 1.358472e−04 |
> | 0.063 | 0.11081 | −1.316774e−03 | −6.588209e−04 | **−7.300143e−05** | **1.459067e−04** |
> | 0.100 | 0.06331 | −1.528923e−03 | −7.650466e−04 | −4.843711e−05 | 9.680016e−05 |
> | 0.160 | 0.01500 | −1.280460e−03 | −6.406404e−04 | −9.609605e−06 | 1.920691e−05 |
> | 0.300 | 0.01500 | −4.871661e−05 | −2.435890e−05 | −3.653835e−07 | 7.307492e−07 |
>
> **The "peaks near x = 0.063" reading survives, measured**: |w_hi − 1| at
> 0.063 is ×1.0745 its x = 0.010 value (×1.0756 on the published column) and
> the half-width ×1.0740, and the peak stays at 0.063 on the shipped anchor
> and on both of registry row 17's alternatives
> (`run_2026-09-06/phase_A_numbers.md` §A2). The correction on `A_zz` is
> exactly ×(`LI6_B1_RANK2_TRANSFER` × `LI6_B1_PER_NUCLEON`) = ×0.30731567 =
> ÷3.253983147 at every x, because the shipped `b2` is `2·x·b₁` and `azz` is then linear in
> b₁; `τ` and `w_hi − 1` follow only approximately (τ's denominator is not).
> Full correction, and the eight-point §8.1c column beside it:
> `run_2026-09-02/phase_C_numbers.md` §8.1c-corr / §8.2-corr / §8.3-corr.

**The tails, `w_tail − 1`** (elastic + unpolarised quasi-elastic, as a
fraction of the Born, at m = ±1):

| x | Q² = 2 | Q² = 5 | Q² = 10 |
|---|---|---|---|
| 0.01 | 8.9654e−05 | 5.18398e−04 | 2.14577e−03 |
| 0.10 | 2.18453e−07 | 1.32298e−06 | 5.19642e−06 |
| 0.30 | — (y < 0.004) | 6.8624e−08 | 2.72302e−07 |

**and the headline that is not what a fixed-target intuition expects:** at EIC
collider kinematics these cells sit at small y (`y = Q²/(x s)`), where
`Y₊ = [1+(1−y)²]/(1−y) ≈ 2`. The radiative-tail dilution is therefore **five
to eight orders of magnitude smaller than at HERMES**. The same model run at
HERMES-like y gives 1.8 % at y = 0.5 and 33 % at y = 0.9. **`rc_tail` is a
small effect at the EIC and a large one at a fixed target, and the difference
is entirely `Y₊` and the Born's 1/Q⁴.**

**And the DEFAULT `rc_tail` is a LOWER BOUND — it is the t-peak ALONE — with
the other edge now shipped.** The leading-log s-/p-peak construction T8(c) used
to build privately is now `ll_peaks_spin1` / `ll_peaks_qe` in
`src/core/rc.cpp`, so T8(c) gates the **shipped** code path, and
`RcTailModel::TPeakPlusLL` (`--rc-tail-model t-peak+ll`) adds it to the tail.
**The two are a band and neither is "the" radiative tail.** The upper edge is a
**STATED MODEL of mixed approximation orders** — POLRAD's η_A quadrature plus a
single-z collinear leading log, good to the worse of the two (~5–10 %), with an
uncancelled soft 1/(1−z) that overshoots as y → 0 — and **not** a controlled
O(α) expansion. **Neither edge carries a tensor s/p peak**: POLRAD supplies
none, inventing one would be a second definition of a physics number, so the
s+p enter the **unpolarised** numerator only and `TPeakPlusLL` **lowers the
tensor fraction of the tail** — measured at Q² = 5 GeV², ⁶Li config 1,
production grid, 2026-09-06: `r_T/r_U` falls by **×0.66139 / ×0.0031829 /
×6.6076e−05** at x = 0.01 / 0.10 / 0.30. Since 2026-09-06 that fraction is
**bounded, not computed** — and the number has to be read with the word.
`--rc-sp-tensor-scale` (`RcOptions::sp_tensor_scale`, default 0, refused
unless `--rc-tail-model t-peak+ll`) lends the **elastic** s-/p-peaks the
elastic t-peak's own σ^el_T/σ^el_U, so scale 1 says *"the s/p tensor fraction
equals the elastic t-peak's"*. **It is a bound with no derivation**, and it is
**empty at the three standard points**: the s/p elastic vertex sits at
Q′²_s = 4.3728 / 4.9382 / 4.9801 GeV², where ⁶Li's coherent form factor is
F_c = −3.75e−45 / −6.44e−51 / −2.41e−51 against **+2.99260** at the t-peak's own
t_min = 8.7304e−05 GeV², so `u_sp/σ^el_U` = 5.25e−79 / 2.07e−86 / 4.06e−83 and
scale 1 is **bit-identical** to 0 there — it recovers **0.0 %** of the collapse
above. It bites on **77 of 3051 accepted cells** (x ≤ 7.94e−03, y ≥ 0.366),
where at scale 1 it reaches **582.9 % of the band half-width** (x = 4.169e−04,
Q² = 1.608, y = 0.9692) and recovers at most **21.72 %** of the collapse. The
**QUASI-ELASTIC s/p column** — `TailTriple::qe_sp`, which is essentially the
whole of the collapse at x ≥ 0.10 (`r_U` grows ×314.18 at x = 0.10 and ×15134
at x = 0.30 between the two tail models, none of it coherent) — is covered by
**neither** `sp_tensor_scale` nor `qe_tensor_scale` and stays exactly
tensor-blind. It is now the largest exactly-zero tensor term in `rc_tail`
(`open_items/run_2026-09-06/phase_B_numbers.md` §B1).

**AND SINCE 2026-09-06 THE ELASTIC HALF IS COMPUTED, NOT BOUNDED.**
`--rc-tail-model polrad-full` — POLRAD **Eq. (18) + Appendix B + Eq. (A.4)**,
the exact τ_A quadrature — carries the s- and p-peaks **with their own
Eq. (A.4) tensor content**, so *"POLRAD supplies no tensor s/p peak"* is true
of **Eq. (38)** and **false of the paper**, and `--rc-sp-tensor-scale` is
refused on that model because the term it stands in for **ran**. Scored against
the computed answer the bound was **not even one-sided**: `r_T/r_U` × the
t-peak is **×0.79395 / ×0.0020032 / ×0.00022591** computed against
×0.66139 / ×0.0031829 / ×6.6076e−05 bounded, at x = 0.01 / 0.10 / 0.30. **The
exact tail is NOT inside the t-peak pair** — between the two edges on 1725 of
3051 accepted cells (56.5 %), covering a median **0.6555** of the gap **over
the 3027 cells whose gap is nonzero** (0.6620 over all 3051; the other 24 have
`t-peak+ll` = `t-peak` exactly, so the fraction is undefined), with the **ratio
of the σ-weighted mean shifts** at **0.428** — a ratio of means, **not** an
event-weighted (nor a σ-weighted) mean of the per-cell fractions, which is
**6.483** — and **above both** in the Q² ≥ 20, y ≤ 0.9 window. It is still
**not** validated against Mo–Tsai (not in this tree). **It does not close the
quasi-elastic half**: a nucleon has no tensor structure function (Eq. (A.5) →
`Im₅…₈ ≡ 0`), so the quasi-elastic tail is tensor-blind at all three peaks on
all three models — and that is the column carrying ×202 and ×7123 of the
t-peak at x = 0.10 and 0.30
(`open_items/run_2026-09-06/phase_B_numbers.md` §B2).

POLRAD §2.1.3 B's "the s- and p-peaks are suppressed" is an **event-weighted**
statement about this generator's bulk and is **false cell by cell**. Both halves,
because either one alone misleads:

* **Event-weighted, it holds — and it is a statement about a P_z.** For ⁶Li
  at Q² ≥ 20 GeV² and y ≤ 0.9, at the CLI's **default fill P_z = 0.7**, the
  mean dilution ⟨w_tail − 1⟩ moves **8.719649e−03 → 8.773752e−03**, **+0.62 %**
  (5194 of 200 000 events, seed 1234), and **8.815031e−03** (+1.09 %) on the
  exact `polrad-full` tail — which is **above both**, so this window's
  event-weighted statement does not BRACKET the exact answer either. At
  **P_z = 0** (the plan both test suites use) the same build gives **5182**
  events and **8.489138e−03 → 8.540716e−03 → 8.581236e−03** (+0.61 %, +1.08 %).
  *(**CORRECTION, 2026-09-15 — re-measured at both P_z.** From 2026-09-06 this
  bullet said the "8.48914e−03 → 8.54072e−03, +0.61 %, 5182" absolutes "no
  longer reproduce" and blamed that run's Phase A for changing the ⁶Li kernel
  normalisation and with it the sampler's cell weights. **Both claims are
  withdrawn**: those digits reproduce exactly on this build at P_z = 0, and the
  2026-09-06 digits are the same run at P_z = 0.7 — the re-measurement had
  switched fill plans without saying so. The
  cross-section-weighted pair 8.31257e−03 → 8.35972e−03 was not re-measured
  and stays withdrawn.)*
* **Per cell, it fails.** **331 of that window's 1356 accepted cells — 24.4 %,
  28.2 % of its cross section — disagree by more than 1 %**, worst **×6444** at
  x = 0.7943, y = 0.0088, Q² = 27.8 (t-peak = 0.016 % of the total), and
  **318 of the 331 sit at y < 0.1**: the failure is at **y → 0**, not y → 1.
* **The surviving per-cell claim** is Q² ≥ 20 GeV² **and** 0.15 ≤ y ≤ 0.7,
  worst **0.55 %** (660 cells); ≤ 0.77 % to y ≤ 0.8, ≤ 2.0 % to y ≤ 0.9;
  **below y = 0.15, none**. The "0.16 % at y ≤ 0.7" published on 2026-09-03 was
  four table rows — three at y = 0.5, one at y = 0.7 — read as a window
  statement.

**Two corners, two mechanisms.** At y = 0.985 the two differ by **59 %**
because z_s = (1−y)/(1−x_A y) → 0 puts the elastic vertex at Q′² → 0 where the
form factor is 1 and the s-peak beats Y₊. At y → 0 it is the **uncancelled
soft radiator**: 1 − z_s = y(1−x_A)/(1−x_A y) and 1 − z_p = y(1−x_A), so
D(z_s) ∝ 1/y — 0.082 / 1.36 / 4.39 / 11.5 at x = 0.01 / 0.10 / 0.30 / 0.72
along Q² = 23.88 GeV², against ratios 1.0014 / 1.061 / 3.98 / 2244 — while the
t-peak has no matching growth (Y₊ → 2) and is crushed because Q² ≥ 20 forces
x·y ≥ 5.0e−3, so low y means **high x**, and its own vertex starts at
t_min = M_A²x_A²/(1−x_A) ∝ x² (0.63 GeV² at x = 0.79). The ratio leaves 1 %
exactly where D(z_s) passes 1 — **a breakdown of the upper edge, not evidence
that the t-peak is low by ×6444.** The excess there is entirely quasi-elastic,
both tails are ≤ 3.2e−04 of the Born, and the price is paid in the **tensor
fraction of the tail**: −1.41817e−08 → −2.20067e−12, **×1.55e−04**, inside the
Q² ≥ 20 window and below the smallest entry of the §B2.1 table (2.59e−04, at
Q² = 3). The absolute tensor term is unchanged (−6.975e−16 on both models).

It fails outright elsewhere too: at the HERMES deuteron point
(x = 0.012, y = 0.85, Q² = 0.53) the t-peak is only **23 %** of the total, low
by a factor **4.36**, and at the low-Q² corner of the generator window
(x = 0.01, y = 0.1, Q² ≈ 4 GeV²) the quasi-elastic **s+p** is **3.35×** the
t-peak unsuppressed and **7.09×** at the shipped Pauli k_F (the s-peak *alone*
is 2.33× / 4.94× — the earlier "s-peak … 3.3×" wording conflated the two).
Never quote either edge as "the" radiative tail.

**The derived quantities that settle three design estimates:**

| x | Q² | edge | (1/6)σ^el_T/σ^el_U | σ^q_U/σ^el_U (S(q) on / S = 1) | ΔA_zz (tail), as published | **ΔA_zz, shipped ⁶Li b₁ (2026-09-06)** |
|---|---|---|---|---|---|---|
| 0.01 | 5 | ho | −5.0943e−04 | 0.28237 / 0.5981 | +3.5167e−07 | **−1.769797e−07** |
| 0.01 | 5 | vmc-ft | −5.4705e−04 | 0.27942 / 0.5918 | +3.2291e−07 | **−2.100985e−07** |
| 0.03 | 5 | ho | −7.3583e−04 | 0.54628 / 0.8679 | — | — |
| 0.03 | 5 | vmc-ft | −8.0426e−04 | 0.53452 / 0.8492 | — | — |
| 0.10 | 5 | ho | **+1.5595e−04** | 2.7475 / 3.164 | +6.6751e−09 | **+2.127452e−09** |
| 0.10 | 5 | vmc-ft | **−4.0702e−05** | 2.4881 / 2.865 | +6.7163e−09 | **+2.042099e−09** |
| 0.30 | 5 | ho | +1.3287e−02 | **1000.2** / 1000.2 | +1.2665e−11 | **+5.130156e−12** |
| 0.30 | 5 | vmc-ft | +1.5556e−02 | 329.85 / 329.85 | +1.7334e−11 | **+9.784578e−12** |

**The last column is new on 2026-09-06 and the one beside it is left as
published** (registry row 6 / §B11): `ΔA_zz` carries `A_zz` through its
`−A_zz r_U` term, so it inherits the ×3.253983 of the WARNING below, while the
`σ^el_T` and `σ^q_U` columns carry no `A_zz` and are unchanged — re-measured,
identical to every printed digit. **`ΔA_zz` changes SIGN at x = 0.01 on both
C0 edges and at all three Q² of `phase_C_numbers.md` §8.1c** (the 2026-09-04
note recorded the flip at one point; it is six), and at x = 0.10 and 0.30 the
sign survives while the magnitude falls by ×0.32 and ×0.41 (`ho`). Recomputed
column, its recipe, and the eight-point version:
`run_2026-09-02/phase_C_numbers.md` §8.1c-corr.

The second row at every point is the other edge of the ⁶Li C0 **shape** band
(`--rc-c0-shape vmc-ft`, added 2026-09-04): the j₀ transform of the committed
ANL VMC point-proton density at the *same* measured ⟨r²⟩_point, so F_c(0),
F_q(0) and ⟨r²⟩ are identical on the two edges and only the unfitted shape
moves. σ^q_U itself is bit-identical between them.

1. **The tensor fraction of the elastic tail changes SIGN with x — but the
   crossing is not located, and at x = 0.10 the sign itself is a band edge.**
   Both edges are negative at x = 0.01 and 0.03 and positive at x = 0.30, so the
   sign *change* survives; the crossing sits between x = 0.03 and 0.10 on `ho`
   and between 0.10 and 0.30 on `vmc-ft`, and at x = 0.10 the table reads
   **+1.5595e−04 on one edge and −4.0702e−05 on the other**, at every Q².
   **The +1.5595e−04 published here on 2026-09-02 is withdrawn as a result and
   stands only as a band edge**; adding `fq_scale` 0…2 widens x = 0.10 to
   −9.0e−05 … +4.0e−04, still spanning zero on both edges. Unaffected: `σ^el_T`
   is a separate table and never a scale factor on `σ^el_U` — the ×(−0.29) the
   C0 band puts on σ^el_T against ×1.10 on σ^el_U is the sharpest demonstration
   of that yet — and design §2.1's `O(10⁻²)` per unit Q_N is met at x = 0.30 and
   missed by 20–60× at x ≤ 0.1 on **both** edges. Everything measured, with the
   two policy overrides that reading the VMC density required:
   `docs/open_items/run_2026-09-03/phase_B_numbers.md` §B1.
   > **Depends on the `A_zz(Born)` column? NO — re-derived 2026-09-06.** Every
   > quantity in this conclusion is `σ^el_T`, `σ^el_U` or `σ^q_U`; none of them
   > carries `A_zz`, and all three columns re-measure identical to every
   > printed digit. **It stands exactly as written.**
2. **`σ^q_U/σ^el_U` runs 0.28 → 2.75 → 1000** over x = 0.01 → 0.30 (0.60 →
   3.16 → 1000 with the Pauli suppression off; 0.28 → 2.49 → 330 on the
   `vmc-ft` edge, which does not change the reading). The design's "30–70 % of the
   ERT" estimate is right only at the bottom of the range: **`rc_tail` is
   quasi-elastic-dominated at every x ≳ 0.03** (22 % of the tail at x = 0.01,
   73 % at 0.1, **99.9 %** at 0.30), so `--rc-qe-suppression` and
   `RcOptions::qe_kf_gev` are the dominant tail knobs almost everywhere.
   > **Depends on the `A_zz(Born)` column? NO — re-derived 2026-09-06.** It is
   > a ratio of two tail cross sections. Unchanged, re-measured. (What the
   > correction *does* touch here is one step further on: the ΔA_zz **spread**
   > of `qe_suppression` and `qe_kf_gev`, because those two move `r_U` — see
   > the run-level rows below and `phase_C_numbers.md` §8.3-corr.)
3. **The tail is not the leading RC systematic at the EIC.** ΔA_zz from the
   whole tail is 3.5e−07 at x = 0.01, Q² = 5, against a **band half-width of
   4.4e−04** — a factor 1300. **CORRECTED 2026-09-04 (task B6): the factor is
   768, not 1300** — both numbers carried the ×3.253983 deuteron-b₁ `A_zz`
   (see the WARNING below). With this generator's ⁶Li b₁ the tail is
   −1.769797e−07 and the band half-width 1.358472e−04. **The ordering is
   unchanged and that is the point of the row**; the arithmetic is now right,
   and at `a_transfer_frac = 1` the band leads by ×1086.
   > **Depends on the `A_zz(Born)` column? YES, in both of its numbers — and
   > it was already corrected; re-derived and confirmed 2026-09-06.** The
   > ratio re-measures **×767.6** (1.358472e−04 ÷ 1.769797e−07), which is the
   > 768 above. What the 2026-09-06 recomputation adds is that the ordering
   > survives with **more** margin than the row claims at the other x too, and
   > that the tail entry it is compared against has changed **sign** since the
   > row was first written (+3.5167e−07 → −1.769797e−07): the row's conclusion
   > is about magnitudes and is unaffected, but the sign may not be quoted
   > from the published column. `phase_C_numbers.md` §8.1c-corr.

**Run-level:**

| quantity | value |
|---|---|
| `RcModel` construction (101 × 77 nodes, `n_eta = 128`) | **0.106 s** — the `Pipeline` setup difference, median of 7. ONE number, shared with `phase_C_numbers.md` §8.3 and `USAGE.md` §7b; the earlier 0.16 s / 0.134 s pair is withdrawn. |
| inclusive throughput, `--rc off` → `--rc tensor-band` (1 core, end to end) | 574 852 → 495 017 ev/s (**−13.9 %**), same build and machine as the row above |
| clipped tail nodes at `tail_max = 10`, global / y < 0.5 / 0.5–0.9 / y > 0.9 | **0 / 0 / 0 / 0** at the default `tail_model = TPeak` (was 1.71 % / 0 / 0 / 21.95 % before the per-nucleon fix). On `--rc-tail-model t-peak+ll`: **0.36 % / 0 / 0 / 4.62 %** on the 101 × 77 CLI grid (1.01 % / 0 / 0 / 5.85 % on T8(d)'s coarser 40 × 24 test grid) — the same y → 1 edge as everywhere else |
| clipped EVENTS — tail / band, default 2000-event inclusive run | **0 / 0** on BOTH tail models, and still 0 at 200 000 events (`meta["rc_clipped_tail_event_fraction"]`, `..._band_...`; a DIFFERENT quantity from the node fractions — nodes are not event-weighted, and no accepted cell centre lands in the nodes `t-peak+ll` clips) |
| whole-run mean `rc_tail`, 200 k inclusive ⁶Li config 1, `t-peak` → `t-peak+ll` | **1.021836 → 1.040369** at P_z = 0 (the suites' fill; at the CLI default P_z = 0.7: 1.021779 → 1.040274) (max 3.410153 → 5.981595; the tail dilution roughly doubles, and every event's kinematics and `weight` are bit-identical — `RcModel::fill` takes no `Rng&`) |
| clipped EVENTS on the BAND, 20 k tagged-alpha | **0.62 %** (⁶Li, seed 1 — the top of a ten-seed scatter 0.485–0.620 %, mean 0.528 %, sd 0.041 %), **2.6 %** (⁷Li, seed 11; seed 1 gives 2.51 %). `tau_tag` reaches 30.7 at the M = 0 density nodes; `band_tau_max = 1` holds every edge in [0.7, 1.3] instead of the −1.79 / −8.35 v0 published. |
| δ(A_zz) at x = 0.01, Q² = 5 from the band, δ_low = 0.19 vs 0.30 | 2.7996e−04 vs **4.4204e−04** — both on the ×3.253983 deuteron-b₁ `A_zz` of the WARNING below. Re-measured 2026-09-06 with this generator's ⁶Li b₁: **8.603659e−05 vs 1.358472e−04** |
| … from the LOW-x ANCHOR itself, δ_low = 0.30 / 0.266 / 0.113 (registry **row 17**, priced 2026-09-06) | band half-widths on A_zz at x = 0.01, Q² = 5 **1.358472e−04 / 1.204512e−04 / 5.116913e−05** (×1 / ×0.8867 / ×0.3767), at x = 0.063 **1.459067e−04 / 1.308566e−04 / 6.313127e−05** (×0.8969 / ×0.4327), at x = 0.10 **9.680016e−05 / 8.798804e−05 / 4.833349e−05** (×0.9090 / ×0.4993) and **1.920691e−05 unmoved** at x = 0.16. **0.113 is the value the source panel READS at x = 0.00966, the x nearest the 0.01 anchor**, so the shipped 0.30 errs **wide by ×2.65**. Two columns and one `meta` key move; the tagged clipped-event fraction does **not** (124/20 000 ⁶Li and 520/20 000 ⁷Li at all three). `run_2026-09-06/phase_A_numbers.md` §A2 |
| … from `fq_scale` 0 vs 2 | ΔA_zz +7.54154e−07 → −4.9925e−08, spread 8.041e−07 — **as published, on the ×3.253983 `A_zz`. Recomputed 2026-09-06 with the shipped ⁶Li b₁: +2.255019e−07 → −5.785810e−07, spread 8.040830e−07 — the SPREAD is unchanged to four digits** (this knob moves only `r_T`), **and the sign change, which is the claim, survives.** |
| … from `tail_tensor_scale` 0.5 vs 2 | +3.59193e−07 → +3.23118e−07, spread 3.61e−08 — **as published. Recomputed 2026-09-06: −1.694451e−07 → −2.055947e−07, spread 3.615e−08** (+0.21 % on the published 3.608e−08 — the third significant figure; this knob too moves only `r_T`). |
| … from `qe_suppression` 0 vs 1 | +1.83644e−07 → +3.51674e−07, spread 1.680e−07 — **as published. Recomputed 2026-09-06: −2.286503e−07 → −1.769797e−07, spread 5.167e−08, i.e. ×0.31** — this knob scales the quasi-elastic tail and so moves `r_U`, so its spread carries `A_zz` and does not survive the correction. |
| … from `qe_kf_gev` 0 (S = 1) vs 0.169 | +5.39491e−07 → +3.51674e−07, spread 1.878e−07 — **as published. Recomputed 2026-09-06: −1.192240e−07 → −1.769797e−07, spread 5.776e−08, ×0.31.** It stays a *larger* knob than `qe_suppression` (5.776e−08 against 5.167e−08), which is the reading this row carries. |
| … from `c0_shape` `ho` vs `vmc-ft` (the C0 SHAPE band, 2026-09-04) | +3.51683e−07 → +3.22911e−07, spread **2.877e−08** — **as published. Recomputed 2026-09-06: −1.769797e−07 → −2.100985e−07, spread 3.312e−08, i.e. it WIDENS by +15 %** (the two C0 edges have different `r_U`). It stays the smallest of the six. *(The published pair's absolute values do not re-measure exactly: this tree gives +3.516736e−07 → +3.229021e−07 for the same two settings — 2.7 × 10⁻⁵ and 2.9 × 10⁻⁵ relative away — while its **spread**, 2.8772e−08, reproduces. Deterministic here across processes; cause not established; `phase_C_numbers.md` §8.1c-corr records it.)* On ΔA_zz it is the *smallest* of these knobs; on `σ^el_T` itself it is the largest by far, because it **flips the sign** at x = 0.10 (§9's table). `RcModel` construction 0.111 s → 1.293 s on this edge. |

| … from `a_transfer_frac` 0 / 0.5 / 1 (the **A = 2 → A = 6 TRANSFER** price, 2026-09-04) | band half-widths on A_zz **1.358472e−04 / 1.518818e−04 / 1.921170e−04** — the band widens by √(1+f²) = 1 / 1.118 / **1.414**. Not a ΔA_zz row like the others: this knob moves the **band** and touches neither the tail nor τ (asserted bit-for-bit). Measured with this generator's ⁶Li b₁, so `f = 0` is the **corrected** 1.3585e−04, not the 4.4204e−04 row above. Default **0.0** and byte-identical output. `run_2026-09-03/phase_B_numbers.md` §B6 |

| … from `qe_tensor_scale` 0 vs 1 (the POLARISED QE stand-in, 2026-09-04) | ΔA_zz **−1.769797e−07 → −2.931806e−07**, spread **1.162e−07** at x = 0.01; +3.017677e−10 at x = 0.10 and +1.795956e−09 at x = 0.30 — **0.086 % / 0.0003 % / 0.25 %** of the band half-width there. Computed with **this generator's own ⁶Li b₁**; see the warning below. It is a BORROWED magnitude, exactly linear in the scale, and possibly ~10² too small at x ≤ 0.1. |

> **WARNING on the five ΔA_zz rows above, added 2026-09-04 (task B3), and its RESOLUTION dated 2026-09-06.** Their `A_zz(Born)` — and therefore the 4.4204e−04 band half-width and every ΔA_zz that carries the `−A_zz·r_U` term — was computed with the **deuteron** Miller b₁ table (`toy_b1`, `B1Mode::Digitized`), not with this generator's ⁶Li b₁ (`Li6B1(MillerB1)`), and is **×3.253983 too large**. The correct band half-width at x = 0.01, Q² = 5 is **1.3585e−04**, and the ΔA_zz rows move with it — the default point, measured both ways, reads +3.51674e−07 above and **−1.76980e−07** with the ⁶Li b₁, i.e. it changes **sign**. The `σ^el_T`, `σ^q_U`, `r_U` and `w_tail` rows are **unaffected** and reproduce exactly.
>
> **"the other rows in the block were not re-measured" is no longer true: all five were, on 2026-09-06 (registry row 6 / §B11), and each row now carries both numbers.** The published values stay in place; the recomputed ones are marked in the same cell. Two things came out of doing it. **(1) The correction is not a common factor.** A knob that leaves `r_U` alone has an `A_zz`-independent spread, so `fq_scale` and `tail_tensor_scale` keep theirs to four digits, while `qe_suppression` and `qe_kf_gev` shrink by ×0.31 and `c0_shape` **widens** by +15 %. **(2) The ORDER of the budget changes at second place, for a bookkeeping reason.** `qe_tensor_scale`'s row was already corrected when it was added in 2026-09-04, so the ladder above compared one corrected number against four uncorrected ones. All six on the same footing: **`fq_scale` 8.041e−07 ≫ `qe_tensor_scale` 1.162e−07 > `qe_kf_gev` 5.776e−08 > `qe_suppression` 5.167e−08 > `tail_tensor_scale` 3.615e−08 > `c0_shape` 3.312e−08** — `qe_tensor_scale` moves from fourth to **second**, because the three below it shrank and not because it grew. **The band still dominates every one of them** (168.9× the largest, on the corrected half-width). Full recomputation and its recipe: `run_2026-09-03/phase_B_numbers.md` §B3.2 and `run_2026-09-02/phase_C_numbers.md` §8.1c-corr / §8.2-corr / §8.3-corr.

**The `fq_scale` band changes the SIGN of ΔA_zz** and is not symmetric about
the nominal — σ^el_T is *quadratic* in F_q. That is why it must be **RUN**
(0, 1, 2) and never rescaled from one row, and why T12 fits a quadratic
through three runs instead of asserting linearity. **Nothing clips on the tail
any more** — dividing the elastic tail by the missing A = 6 removed the whole
21.95 % that used to sit in the `y > 0.9` band — but the four node fractions,
the two per-event fractions and `Event.rc_clipped` are all still reported,
because `Y₊ ~ 1/(1−y)` is still where a wider `y` window would bite. What does
clip is the **band on the tagged channels**, 0.6–2.6 % of events, where
`tau_tag = 1 − n̄/n_M` diverges at the nodes of `n_M`.

**Honest flags, mandatory wherever any of this is quoted.** `δ_low = 0.30` is
the **size of a correction this generator does not apply**, taken as a 1σ
band; its source (Gakh–Shekhovtsova, hep-ph/0403262) has **zero INSPIRE
citations** — **re-verified 2026-09-04** against the INSPIRE literature API
(record 647050, `citation_count: 0`); the flag stands. It is also, measured
from the paper's own figure data, a **Q² = 0.1 GeV² number quoted at Q² = 5**:
the paper states its "10 % to 30 %" for x ∼ 10⁻³–10⁻², and **exactly one** panel
of its Fig. 2 lies in that window — panel (a), Q² = 0.1, x = 0.00226–0.00966,
all 30 of its points inside it, where |δ| = 0.113–0.266. Panel (b) begins at
x = 0.01348, *above* 10⁻², and carries no point in the quoted range. So the
10 %/30 % pair is **that one panel's two x ends at a single Q²** — |δ| = 0.1133
at x = 0.00966 and 0.2662 at x = 0.00226, the paper scoping them to *x* and
never to Q² — and **not** "the spread of one calculation over Q², not two
independent edges", which `PHYSICS_CHANNELS.md`, `include/lipolgen/rc.hpp` and
`design_C_tensor_rc.md` all said until 2026-09-04. **The anchor's real basis,
written down here for the first time:** δ_low = 0.30 hangs at `RC_X_LOW` = 0.01,
*above* that panel's top x = 0.00966, where the paper reads |δ| = 0.113; the
0.266 is the panel's value at x = 0.00226, a factor 4.3 lower in x. "The
conservative end of a Q² spread" and "the panel's lowest-x value extrapolated
upward in x" are different bases, and the shipped documents asserted the first
while the anchor rests on the second. Conservative in magnitude either way
(0.30 > 0.266 > 0.113), which is all `w = 1 ∓ δτ` uses, so **nothing moves** —
recorded, not changed (`run_2026-09-03/phase_B_numbers.md` §B6.3).
**Since 2026-09-06 the two alternatives are PRICED and not merely named**
(registry row 17; `run_2026-09-06/phase_A_numbers.md` §A2), and the direction
of the error is the part to carry to any reader: **0.113 is the value the
panel actually reads at the x NEAREST the anchor**, so the shipped 0.30 is
**×2.65** it and **×1.13** the 0.266 at the panel's bottom — it errs **wide**,
which is the safe direction for a half-width and the wrong one for a quoted
precision. Running `--rc-delta-low-x` 0.30 / 0.266 / 0.113 at the
configuration of this section's own band table (⁶Li, `--config 1`,
`--plan tensor-thirds --pzz 0.6`, Q² = 5, P_zz = +1) but with **this
generator's ⁶Li b₁** — the corrected `A_zz` of the WARNING below, not the
table's own ×3.253983 column — the band half-width on `A_zz` is

| δ_low | x = 0.010 | x = 0.063 | x = 0.100 | x = 0.160 |
|---|---|---|---|---|
| **0.30 (shipped)** | **1.358472e−04** | **1.459067e−04** | **9.680016e−05** | 1.920691e−05 |
| 0.266 (panel bottom, x = 0.00226) | 1.204512e−04 | 1.308566e−04 | 8.798804e−05 | 1.920691e−05 |
| **0.113 (panel AT the nearest x, 0.00966)** | **5.116913e−05** | **6.313127e−05** | **4.833349e−05** | 1.920691e−05 |

— the anchor is carried in full only at x ≤ `RC_X_LOW` (×0.8867, ×0.3767),
**compressed** between the anchors by the log-linear interpolation (×0.8969 /
×0.4327 at 0.063, ×0.9090 / ×0.4993 at 0.100) and **identically zero** at
x ≥ `RC_X_HIGH` = 0.16, where the E12-13-011 anchor pins it. The band still
peaks at x = 0.063 on all three (and sharpens as the anchor drops), and still
leads the **whole** radiative tail at x = 0.01 by **×768 / ×681 / ×289**, so
the RC budget's ordering survives every choice. It moves **two npz columns
(`rc_tensor_lo`, `rc_tensor_hi`) and one `meta` key**, and — measured, not
assumed — **not** the tagged band's clipped-event fraction: the same 124 of
20 000 ⁶Li events (0.62 %, seed 1) and 520 of 20 000 ⁷Li (2.6 %, seed 11) clip
at all three anchors, because `clamp_tau` clips |τ| against `band_tau_max` and
δ never enters it. What the anchor sets **at** the clip is the width there:
|τ| = 1 exactly, so those events' edges are exactly 1 ∓ δ_low —
[0.700, 1.300] / [0.734, 1.266] / [0.887, 1.113]. **Nothing is decided and no
default moved**; `python/tests/test_rc_low_x_anchor.py` pins the table.
**The HIGH anchor has no alternative to run**: the record names none.
`RC_DELTA_HIGH_X = 0.015` is E12-13-011's 1.5 %, which lives in the
**unpublished proposal only** — the published companion arXiv:2506.04506 has
no radiative-correction discussion at all — and the only alternative the tree
states is an *action*, "cite it by page or **drop** the anchor"
(`docs/CONVENTIONS.md` (a), `run_2026-09-02/design_C_tensor_rc.md`), not a
second value. There is nothing to band it against, and that is the honest
statement of it. Throughout,
**δ ≡ (Δσ_RC − Δσ_Born)/Δσ_Born** — the RC as a fraction of the Born, which is
the base the paper's own sentence names ("as compared with the Born
contribution") and the base `w = 1 ∓ δτ` applies. It is **inferred** from the
paper's two plotted curves; the paper defines no δ. **Withdrawn 2026-09-04**:
the "Q² = 0.1–1 GeV²" above (panel (b) is out of range) and the "0.171–0.240 at
Q² = 1" this block used to print as a |δ| range — that pair is panel (b)'s two
*signed* extremes, −0.171 and +0.240, at opposite ends of its x window, while
its |δ| actually spans 0.006–0.240 (`run_2026-09-03/phase_B_numbers.md` §B6.3).
`RC_DELTA_LOW_X_OPTIMISTIC = 0.19` is HERMES's *measured*
fractional residual at its lowest-x bin (2×10⁻³ on A_zz = −1.06×10⁻²).
Everything cited is **deuteron** — whether the deuteron's fractional RC
transfers to ⁶Li is the largest unquantified assumption here. **Since
2026-09-04 that assumption has a price rather than only a sentence:**
`RcOptions::a_transfer_frac` / `--rc-a-transfer-frac` adds `f·δ(x)` in
quadrature, `δ_eff = δ√(1+f²)`, at **default 0.0** — which *is* the
unstated v0 assumption, "the deuteron fraction transfers exactly" — with
`f = 0.5` and `f = 1` as the priced edges. No measurement prefers any of
them; `meta["rc_a_transfer_frac"]` records which ran. The ⁶Li
form-factor shape parameters `(a, α, q₀) = (1.9069 fm, 0.13822,
3.0999 fm⁻¹)` and `(q_z, b) = (1.30 fm⁻¹, 1.85 fm)` are **unfitted starting
values** (the Suelzle–Yearian–Crannell and Li–Sick–Whitney–Yearian elastic
data are not in this repository in machine-readable form — **Q1 and Q10 stay
open**); their *normalisations* are the **measured** moments (μ = +0.822047
μ_N, Q = −0.0818 fm², TUNL A = 6) and never VMC, whose Q(⁶Li) = −0.23(9) fm²
is 3× the measured one. Since 2026-09-04 the **C0 shape** — `(a, α)`, i.e. the
monopole `F_c` and `F_q` share — is banded rather than left to look settled:
`--rc-c0-shape ho | vmc-ft`, both edges pinned to the same measured ⟨r²⟩_point,
F_c(0) and F_q(0), differing by ×2.5 in `F_point` at q = 2 fm⁻¹ and flipping
the sign of `(1/6)σ^el_T/σ^el_U` at x = 0.10. **Run both edges**; a band is a
price tag and not a fit, and Q1 is unchanged by it. The ⁶Li form-factor DIP
location is likewise a model number on either edge: q₀ = 3.0999 fm⁻¹ ⇒
\|t\| = 0.3742 GeV², gated over [2.9, 3.3] fm⁻¹ by T11, while the committed
VMC point-proton density has no C0 zero below q ≈ 4.3 fm⁻¹ at all (the uncited
"\|t\| ≈ 0.31 GeV²" that `PHYSICS_CHANNELS.md` carried until then is
withdrawn — `run_2026-09-03/phase_B_numbers.md` §B1.9). **No RC calculation exists for a tagged tensor
asymmetry** (the tagged band is an uncited extrapolation) and none for any
φ-dependent tensor observable at any axis. The **polarised** quasi-elastic
tail is **not computed anywhere** (POLRAD supplies no tensor partner to
Eq. (44); nothing exists for an A = 6 spin-1 nucleus) and is neglected by
default citing Zhou *et al.*, PRL **82** (1999) 687 — a *deuteron* statement,
and at x = 0.30 the tail is **99.9 %** quasi-elastic, so that zero is the
largest unpriced piece of `rc_tail`. `--rc-qe-tensor-scale` (default **0.0**)
prices the omission by lending the quasi-elastic tail the **elastic** tail's
own tensor fraction: a **borrowed magnitude, not a derived bound**, possibly
~10² too small at x ≤ 0.1, and with a meaningless sign
(`run_2026-09-03/phase_B_numbers.md` §B3). The **tensor fraction of the
leading-log s-/p-peaks** is the same shape of omission one level down, and
since 2026-09-06 it too is **bounded, not computed** — `--rc-sp-tensor-scale`
(default **0.0**, refused unless `--rc-tail-model t-peak+ll`) — but that bound
is **empty at x = 0.01 / 0.10 / 0.30, Q² = 5** (bit-identical to 0: the
coherent s/p vertex sits at Q′² ≈ Q², 45–51 decades below the form factor at
the t-peak's own t_min) and recovers at most **21.72 %** of the tensor-fraction
collapse anywhere on the grid. **The QUASI-ELASTIC s/p column is bounded by
neither knob and is now the largest exactly-zero tensor term in `rc_tail`**
(`run_2026-09-06/phase_B_numbers.md` §B1) — **and `--rc-tail-model
polrad-full`, the exact Eq. (18) tail shipped the same day, does NOT close it
either: a nucleon has no tensor structure function (Eq. (A.5) → `Im₅…₈ ≡ 0`),
so the quasi-elastic tail is tensor-blind at all three peaks on all three
models (§B2). What `polrad-full` DOES compute is the ELASTIC s/p tensor peak,
and scored against it the `sp_tensor_scale` bound was not even one-sided.**
**`rc_tail` at the DEFAULT is the t-peak only
and is a LOWER BOUND**, and "fine at Q² ≥ 20 GeV²" is **event-weighted and
never per cell** (T8(c), T8(d)(i)): the window mean moves only **+0.61 %**, but
over Q² ≥ 20 GeV² with **no y cut 346 of 1399 accepted cells — 24.7 %, 29.6 %
of that window's cross section — differ from the `t-peak+ll` edge by more than
1 %**, **318 of them at y < 0.1**, worst ×6444 — a breakdown of the *upper*
edge at y → 0, not a measured deficit of the t-peak. The **only** per-cell
agreement statement is **Q² ≥ 20 GeV² and 0.15 ≤ y ≤ 0.7** (660 cells, worst
0.55 %). What is unchanged is the lower bound elsewhere: the t-peak is **low by
4.4×** at the HERMES deuteron point **on the leading-log estimate, and by
2.36× on the exact `polrad-full` tail — the upper edge overshoots there by
1.85×** (the "POLRAD §2.1.3 B" block earlier in this section; §B2.4). The **tagged band is clamped** at
`|τ| ≤ band_tau_max = 1`, which is a *choice*: `n_M → 0` is exactly where the
fractional-rescale ansatz breaks down, and the clipped fraction is reported
rather than hidden.

**A RECORDED REJECTION, so it is not re-litigated: the Gakh–Shekhovtsova
*shape* is deliberately NOT adopted (2026-09-04, task B6).** The obvious
"better" move — replace the log-linear `rc_delta` with the shape digitised
from hep-ph/0403262's own Fig. 2 — was investigated and **refused on measured
grounds**, not overlooked. The arXiv source ships the four panels as EPS with
their `/x` and two `/y` arrays (RC-included and Born, `Q² = 0.1, 1, 4, 10`),
so this is read off the data, not off a picture:

1. **δ is SINGULAR, not a bounded fraction.** The Born Δσ **crosses zero** (the
   b₁ zero-crossing) inside the Q² = 4 and Q² = 10 panels and the RC **moves
   the crossing** — the paper's own stated result, reproduced here as
   x₀ = 0.20118 → 0.18467 (Q² = 4) and 0.20161 → 0.18331 (Q² = 10). So
   δ = (Δσ_RC − Δσ_Born)/Δσ_Born — whose pole is that Born zero, the
   denominator's — spans **−1.691 … +7.810** at Q² = 4 and
   **−0.887 … +4.091** at Q² = 10 (at Q² = 4: −0.170 at x = 0.111,
   **−0.537** at 0.169, **+1.026** at 0.214, +0.030 at 0.283), and changes
   sign already at Q² = 1 (+0.240 at the low-x end → −0.171 at the high).
   The band ansatz `w = 1 ± δτ` with
   `band_tau_max = 1` assumes a small, bounded, one-signed rescale and needs
   |δ| < 1 to keep both weights positive. A paper-shaped δ would emit
   negative weights.
2. **It is not monotone non-increasing in x**, which is the design's stated
   justification for the log-linear form and is *asserted* by
   `tests/test_rc.cpp:70-83` (T3: monotone on a 2000-point log grid, both
   anchors compared with `==`). The panels give **−0.128 at x = 0.386**
   (Q² = 4) and **−0.400 at x = 0.85** (Q² = 10) against δ_high = 0.015 for
   all x ≥ 0.16 — **8.5×** and **26.7×** past the E12-13-011 anchor in
   magnitude, and negative, which a half-width cannot be. The two
   anchors and the shape cannot both stand.
3. **Patchy, non-rectangular Q² support.** Four Q² values with x windows
   0.00226–0.00966 / 0.01348–0.09662 / 0.05394–0.38647 / 0.13097–0.85000 —
   the first two do not even overlap. This generator's own quoted point
   **x = 0.01, Q² = 5** falls *between* the Q² = 4 and 10 panels and outside
   the x window of both, with its x in the gap between panels a and b. There
   is nothing to interpolate from where the design publishes.
4. **Their Born is not our Born.** The ratio is against *their* tensor model —
   Eq. (60)'s HERMES A_zz parametrisation with ALLM97 F₂ and the
   Callan–Gross-type b₂ relation Eq. (59) — not this generator's b₁, so a
   fractional RC computed against a different b₁ shape is a fresh unpriced
   assumption on top of the A = 2 → A = 6 one.

**Kept:** the interpolation, the two visible anchors, and the *magnitude* with
its flags. Reopening the digitisation route requires a δ definition that
survives a zero crossing. `run_2026-09-02/design_C_tensor_rc.md` Q7;
`run_2026-09-03/phase_B_numbers.md` §B6.

## 10. b₁ of ⁶Li — a four-term α–d convolution, opt-in, and the A = 2 gate now PASSES — **CLOSED FOR ONE CONFIGURATION** (rows 1, 2, 3, 4, 18, 25 open)

**IMPLEMENTED (2026-09-03) AS AN OPT-IN BACKEND; ITS A = 2 GATE CLOSED THE SAME
DAY; THE ITEM ITSELF IS ALL BUT CLOSED — two author decisions and one
unquantified systematic remain, listed at the end**
(`include/lipolgen/b1_nuclear.hpp`, `src/core/b1_nuclear.cpp`;
`PipelineConfig::b1_model` / `--b1-model {miller,cdks,li6-convolution}` plus
`--b1-band-scale` and `--b1-alpha-d-dwave-weight`; `docs/USAGE.md` §2a; design
`docs/open_items/run_2026-09-02/design_D_b1_li6.md`; the gate as it now stands
in `docs/open_items/run_2026-09-03/phase_A_numbers.md` and `phase_A_cdbonn.md`;
the superseded failing verdict in `phase_D_gate.md`; the ⁶Li numbers in
`phase_D_numbers.md`).

> ### ⚠ READ THIS BEFORE QUOTING ANY NUMBER BELOW
>
> Design §5 makes the A = 2 validation gate **blocking**: the same kernel, fed
> the AV18 deuteron u(k), w(k) instead of the α–d waves, must reproduce the
> digitized CDKS Fig. 4 before any ⁶Li number is quoted. **Since 2026-09-03 it
> does.** Read at the end of §5.4's checklist — MSTW2008 LO, CDKS's own nucleon
> PDF (item 4), at CDKS Eq. (21)'s δ-function — G3a passes on all three
> landmarks and **G3b's peak ratio is 0.843243**, inside the factor-2 window
> with the lower edge cleared by a factor 1.69. With the real CD-Bonn wave
> function (item 5) as well it is **1.000338**. Design §5.4's *Escalation*
> clause is not triggered, so **the ⁶Li publication ban is LIFTED** — for that
> configuration. The tables below are **not** in it; see condition 1 and
> "The ⁶Li numbers" heading further down.
>
> **Four conditions travel with that, and none of them is optional.**
>
> 1. **The lift is about a CONFIGURATION, not about a build.** The gate passes
>    for the **MSTW2008 LO** nucleon input at CDKS Eq. (21)'s δ-function. On
>    the **shipped default** unpolarised backend, the library `ToyF2`, the same
>    gate gives **0.440** — *outside* the [0.5, 2] acceptance window. That row
>    is still pinned in the doctest as a regression guard and is *not* the
>    verdict: §5.4's Escalation clause reads the ratio *after* the checklist,
>    and item 4 is the nucleon PDF. Since 2026-09-04 the passing configuration
>    is emittable from the run surface — `PipelineConfig::b1_unpol` /
>    `--b1-unpol mstw` threads `MstwSF` into `Li6ConvolutionOptions::unpol`
>    through `default_inclusive_kernel` — so the rule is **quote numbers made
>    with `--b1-unpol mstw`; the default `toy` backend is outside the window
>    and its numbers are not covered by the lift**, and they cannot be rescaled
>    into covered ones (mstw/toy on `Li6ConvolutionB1::b1(x, 2.5)` is
>    **1.848 / 1.276 / 0.817** at x = 0.10 / 0.30 / 0.50 — a shape change).
>    Selecting `mstw` needs the optional PYTHIA tier and is **refused at
>    configuration time, never downgraded**; a build that cannot reproduce the
>    verdict row therefore cannot emit a number claiming it: `b1_nuclear T1v`
>    is skipped at REGISTRATION time, which doctest reports only in the
>    **skipped count** — never by name — so the name and the reason reach the
>    log through its companion case `T1r`, which prints in every build whether
>    the row was measured (`tests/test_b1_nuclear.cpp:881`).
> 2. **Comfortable at Eq. (21), marginal at Eq. (17).** MSTW at κ = 1 gives
>    **0.520** — inside by 4 % of its own value. "The gate passes" without "at
>    CDKS Eq. (21)'s δ-function" oversells it.
> 3. **The CD-Bonn agreement is better than the reference deserves.** The
>    target is a *digitization of a published figure* and this gate cannot
>    resolve the low-x zero better than a few per cent. Read 1.000338 as "the
>    residual is now below the digitization error", never as three-digit
>    agreement with CDKS — and it is specific to CD-Bonn **and** MSTW together.
> 4. **The gate is A = 2.** It validates the kernel on the **deuteron** and
>    tests nothing about the α–d step, for which no measurement exists at any
>    A > 2. **The mandatory ±100 % band on every ⁶Li number therefore stays** —
>    it comes from Q(⁶Li) against Q_d, not from the gate — and so does the rule
>    that a number is quoted with the configuration that made it.

### What was built

Per nucleon, ⁶Li (design §1.7):

> b₁^{⁶Li}(x,Q²) = (2/6) ∫ (dz/z) { [f_S(z) + w_CG f_D(z)] b₁ᵈ(x/z,Q²)
>   + w_αd δ_T f_αd(z) F₁ᵈ(x/z,Q²) } + (4/6) w_αd ∫ (dz/z) δ_T f_α(z) F₁^α(x/z,Q²)

with **w_CG = 1/10 exactly** — the Clebsch–Gordan tensor dilution of the
deuteron inside the L = 2 α–d component, which is the closed form quoted next
to `LI6_B1_RANK2_TRANSFER` — and w_αd = 1 nominal. The four terms are
`b1_embedded_s` (1), `b1_alpha_d_dwave_d` (2d), `b1_alpha_d_dwave_alpha` (2α)
and `b1_cg_dwave` (3); they sum to `b1` to 1e-12 (T7, P2). Inputs: α–d
magnitudes from `data/vmc/momenta/li6_ad1.momentum`, the S–D **sign structure**
from `data/vmc/li6_alpha_d/li6.ad`, the deuteron b₁ᵈ from whichever `TensorSF`
camp the caller injects (default: the **raw** digitized CDKS theory-1 column,
`cdks_b1_raw_per_nucleon()`), F₁ from **CDKS Eq. (22)**
`(1 + γ²) F₂/(2x(1+R))` — *not* `NuclearF2::f1a`, which is the massless form
and differs by 1 + γ² = 1.35 at x = 0.5, Q² = 2.5.

### The A = 2 gate, measured

| clause | result |
|---|---|
| G0a–G0g — CG algebra, Eq. (21)'s two coefficients, ∫δ_T f dy = 0, `k_range` against a brute-force scan | **PASS** (exact / 1e-12 / endpoints to the scan step) |
| G1a–G1f — ∫k²(u²+w²)dk = 0.99998, P_D = 0.05760, ∫f dy = 0.98766 against the moment identity's 0.98707, ⟨y⟩ = 0.99526, y_max = 1.4976, the F₁ᴰ/F₁ᴺ EMC shape | **PASS** |
| layer 2 — the α–d sign gate | **PASS**: Q(α–d) = **−0.3333 fm²** < 0, and the same code on the deuteron pair returns `fdeut.av18`'s own header qm = 0.269673 as **0.269362** (0.12 %) |
| **G3a** — exactly two sign changes in (0, 1.0], positions, peak position | **PASS, QUALIFIED.** It passes on every configuration the clause is *applied* to, and must be read with **"G3a's stated limitation — the counting ceiling"** below, which measures what the (0, 1.0] ceiling excludes. One measured configuration is **exempt by construction and is not a pass**: **ToyF2 + CD-Bonn** has **zero** sign changes in (0, 1.0] — `tests/test_b1_nuclear.cpp` checklist item 5 pins it as `CHECK(t.nz == 0)` and deliberately does *not* call `check_g3a` on it, because applying a shape clause there would gate the *wave function* on the *toy* PDF; both real PDFs put the two zeros back. Verdict row (MSTW2008 LO, Eq. (21)): zeros **0.0279** (falling) / **0.4952** (rising) and peak at **0.7716** against the digitized 0.0656 / 0.4572 / 0.7657 — misses of 0.038 / 0.038 / 0.0059 against tolerances 0.08 / 0.10 / 0.10. CDKS's own PDF reproduces CDKS's own figure better than either stand-in on all three landmarks, with nothing tuned. **Read the margin**: the low-x zero is resolved to a few per cent at best on any of these grids (0.0220 ToyF2, 0.0098 CT18NLO, 0.0279 MSTW) and must not be quoted to more than two significant figures; the second zero and the peak position are the discriminators. The counting window is **(0, 1.0]** and the scan floor **0.001**, both stated in the design since 2026-09-03 — the old [0.02, 1.0] was narrower than the ±0.08 tolerance it bracketed, and with CT18NLO the zero fell below it and aborted the clause |
| **G3b** — peak magnitude within a factor 2 | **PASS**: max\|x·b₁\| over [0.10, 0.80] = **9.15096e−4** against the digitized 1.08521e−3, **ratio 0.843243**, a factor **1.19 low**, inside [0.5, 2] with the lower edge cleared by a factor 1.69. With CD-Bonn as well: **1.000338**. Same clause on the other rows: ToyF2 at Eq. (21) 0.440 (outside), ToyF2 at κ = 1 0.272, CT18NLO at Eq. (21) 0.719, MSTW at κ = 1 0.520 |
| G3c — Close–Kumano, **reported, not enforced** | this kernel **+2.24896e−4** (MSTW, Eq. (21)) against the digitized CDKS **+4.59200e−4** — 0.49 of it, essentially where ToyF2 was (0.467). **The PDF that fixes the peak does not fix the integral.** With CD-Bonn as well it is +4.48580e−4 = 0.98, which is a *new* fact and not a restatement of the peak one. Miller's is +5.9147e−3. None is zero; none is enforced |

**G3a's stated limitation — the counting ceiling.** G3a counts sign changes in
**(0, 1.0]**. That upper edge is a **scope choice, not a derived tolerance**,
and it is load-bearing. Over the digitized reference's **own** domain
[0.010, 1.590] the reference has **two** sign changes and stays **positive**
from 0.4572 to its last point (+4.040e−6 at x = 1.590), while **every**
computed configuration crosses zero a **third** time and is **negative** at
the top. So over the reference's own domain the computed curve has **three**
sign changes and the reference has **two**; "exactly two" holds only because
of the ceiling. Measured on the gate's own grid `linspace(0.001, 1.59, 300)`
and pinned by `tests/test_b1_nuclear.cpp` **T17**:

| configuration | third zero (falling) | digitized x·b₁ there | % of the digitized peak |
|---|---|---|---|
| ToyF2, Eq. (21), AV18 — **the shipped default** | 1.220437 | +7.6061e−5 | 7.01 % |
| **MSTW2008 LO, Eq. (21), AV18 — the verdict row** | **1.217660** | +7.7860e−5 | 7.17 % |
| CT18NLO, Eq. (21), AV18 | 1.197722 | +9.1937e−5 | 8.47 % |
| ToyF2, Eq. (17) κ = 1, AV18 | 1.142281 | +1.4599e−4 | 13.45 % |
| MSTW2008 LO, Eq. (17) κ = 1, AV18 | 1.137076 | +1.5245e−4 | 14.05 % |
| CT18NLO, Eq. (17) κ = 1, AV18 | 1.133695 | +1.5680e−4 | 14.45 % |
| ToyF2, Eq. (21), CD-Bonn | 1.500803 | +8.3317e−6 | 0.77 % |
| MSTW2008 LO, Eq. (21), CD-Bonn | 1.493052 | +8.8438e−6 | 0.81 % |
| CT18NLO, Eq. (21), CD-Bonn | 1.473742 | +1.0256e−5 | 0.95 % |

**Widen the window to (0, 1.59] on exactly the argument used for the floor and
G3a fails on every configuration, MSTW + CD-Bonn included.** The floor was
widened in phase A on the argument that *"bookkeeping that brackets a physics
tolerance has to be at least as wide as the tolerance, or it is a second,
tighter, unstated cut."* Applied to the ceiling that argument does **not**
bite the same way — ±0.10 about 0.4572 reaches 0.5572 and the peak clause's
±0.10 about 0.766 reaches 0.866, so (0, 1.0] already brackets every tolerance
this clause writes down, with room. What the floor argument does not reach,
and what the ceiling actually is, is the **count**: a counting clause is a
statement about a domain, and the domain decides the answer. G3a's domain is
chosen, and it is narrower than the reference it compares against. **The phase
that re-opened the window argument moved the edge that cost it nothing and
stopped at the edge that would have cost it the clause.** That is recorded
here rather than left for the next reader to find.

**Two defences of the ceiling, and what survives of them.**

1. *"x > 1 is out of range."* **False, and it must not be used.** x per nucleon
   for a nucleus runs to A — to 2 for the deuteron this gate is about — so
   x ≈ 1.13–1.22 is kinematically allowed and physically meaningful (the Fermi
   motion / short-range-correlation region), and it is **inside** the digitized
   reference's own domain, which CDKS themselves plotted to 1.59.
2. *"the reference is unreadable that high."* **Partly true, and only for the
   CD-Bonn rows.** At the AV18 crossings the digitized x·b₁ is still **7.0–8.5 %**
   of its own peak, and at the κ = 1 crossings **13.5–14.5 %** — far above
   anything one can call a digitization floor, and the digitized column there
   is a smooth monotone decay, not jitter. So for those rows the third crossing
   is a **real sign disagreement in the several-per-cent-of-peak region**, not
   noise. The defence becomes honest only near x ≈ 1.47–1.50 (the CD-Bonn rows),
   where the reference is 0.8–1.0 % of its peak.

**Status: G3a is a QUALIFIED pass, and the clause is not changed here.**
Recorded, not resolved, because changing an acceptance criterion is the
design's and the author's call, not a correction pass's — and because moving
it in either direction after seeing the answer is exactly the failure mode
this note is about. What a reader must take from it:

* The ban lift is triggered by design §5.4's **Escalation** clause, which reads
  the **peak ratio** — that is G3b, and the ceiling does not touch it.
* But G3a is listed as *"shape, **hard**"* in an acceptance criterion whose own
  words are *"three parts, all must hold"*. **If the ceiling were widened to the
  reference's own domain, G3a would fail on every configuration and the
  criterion as a whole would not hold**, so item 10's close would have to be
  re-argued. Anyone quoting "the A = 2 gate passes" is quoting a pass that
  includes this scope.
* Nobody may quote "exactly two sign changes" without the window. The honest
  sentence is **"exactly two sign changes in (0, 1.0], and a third at
  x ≈ 1.22 that the window excludes."**

**The residual is attributed, and after this phase there is essentially none
left.** The §5.4 checklist was worked in order with no fudge factor tuned.
Starting point: the κ = 1 `ToyF2` column, ratio 0.271689. All effects measured
on the doctest grid `linspace(0.001, 1.59, 300)`.

| item | measured effect on the peak | direction |
|---|---|---|
| 0. finite-\|q⃗\| δ-function, κ = 1 → √(1+γ²) | **×1.6195** (ToyF2) / ×1.6206 (MSTW) | closes — **IN**: the gate's default since 2026-09-03 |
| 4. nucleon PDF, `ToyF2` → **MSTW2008 LO**, CDKS's own | **×1.9152** (at κ = 1) / ×1.9165 (at κ) | closes — **the verdict row**, opt-in through `MstwSF` |
| *(the retired CT18NLO stand-in, for reference)* | ×1.6680 at κ = 1; CT18NLO → MSTW is ×1.1721 | superseded |
| 1. target mass (CDKS Eq. 22, already in) | ×1.49 at x = 0.8 (1.52 at κ = 1) | in |
| 2. R (`r1998`, CDKS's own; in **on the gate** — the ⁶Li backend uses `r_sigma_lt`, see below) | ×0.98; swapping it back is +2.7 % on G3b | in |
| 5. wave function, AV18 → **the real CD-Bonn** | **×1.18631** (MSTW, Eq. (21)); ×1.198–1.202 on the other three rows | closes — opt-in through `Options::wave = kCdBonn` |
| *(the retired D-state rescaling proxy for item 5)* | ×0.88 — it had the **sign of the effect backwards** | superseded |
| **the gate's default wave function + MSTW + Eq. (21)** | **ratio 0.843243 (factor 1.19)** | **INSIDE G3b** |
| **CD-Bonn + MSTW + Eq. (21)** | **ratio 1.000338** | — |

Items 0 and 4 are independent to **0.07 %**: 0.271689 × 1.91519 × 1.61951 =
0.842682 predicted against 0.843243 measured, so there is no third unexplained
factor at the peak. The 1.186 that was left over after them is supplied, to
1.3 % across two PDFs and both δ-functions, by item 5 — **the wave function,
not a fudge**. `Li6ConvolutionB1`'s own `finite_q_delta` stays `false` (A10):
in ⁶Li the switch is worth −1 % to +7 %.

**Two things this budget does not say.** The convolution is not a pointwise
multiplication: MSTW/CT18NLO on the isoscalar F₁ is 1.40–1.44 at the peak x
but the peak x·b₁ moved only ×1.1721, and the peak itself moves
0.7344 → 0.7716. And **x = 0.50 got worse**: MSTW gives +1.65e−5 there against
the digitized +1.48e−4, a factor 9 low, because MSTW's second zero sits at
0.4952 and x = 0.50 is 0.005 past it while the digitized curve crossed at
0.4572. No clause of G3 sees that — G3a gates the zero's *position* (cleared by
0.062) and G3b the peak — but a reader comparing curves at x = 0.5 will.

### The ⁶Li numbers — REGENERATED AS A BAND on 2026-09-06

x·b₁ per nucleon, Q² = 2.5 GeV², default options — which include **ToyF2 F₁
with `r_sigma_lt`** (not the gate's `r1998`; see below) and **κ = 1** (not the
gate's Eq. (21)), and the 2400 / 2400 / 3200 y grid the review of 2026-09-03
set. Every row is the centre of a mandatory {0, 1, 2} × b₁ band:

> **THESE VALUES WERE REGENERATED AS A BAND ON 2026-09-06, and the band — not
> the row below — is what may be quoted.**
> `run_2026-09-06/phase_A_li6_tables.md` reruns **the same observables at the
> same (x, Q²) points with the same options**, over
> `--b1-unpol {toy, ct18nlo, mstw}` × `--unpol-sf {toy, ct18nlo, mstw}`. The
> table below is the **`toy`/`toy`** cell of that band and reproduces there to
> every printed digit; it is the band's **centre-of-record, not its answer**.
>
> **Which axis moves what — measured, not assumed.** `--b1-unpol` is the
> deuteron F₁ *inside* CDKS Eq. (22) and moves **terms (2d) and (2α) and
> nothing else**: terms (1) and (3) convolve b₁ᵈ and are **bit-identical**
> across all three backends, as are the densities, the T14 remainders and the
> `w_αd = 0` row. `--unpol-sf` is the **kernel's own** F₁/F₂ — the
> **denominator** of A_zz — and is **exactly flat on x·b₁** (bit-identical
> through the pipeline at all five x for both non-toy `b1_unpol` values),
> because `Li6ConvolutionB1::b1` ignores its `f1` argument. So the band on
> **these tables is three wide, not seven**; A_zz itself moves on both axes
> (a factor **2.23** across the seven legal cells at x = 0.10, Q² = 2.5,
> y = 0.5). Two of the nine cells — `--b1-unpol toy` under a non-toy
> `--unpol-sf` — are **refused by name** by `validate()`, because
> `meta["b1_unpol"]` would then record "toy" for something else.
>
> **The band on x·b₁ (Q² = 2.5), and then the ±100 % on top of it:**
>
> | x | `--b1-unpol` envelope | **× {0, 1, 2} — quotable** | spread |
> |---|---|---|---|
> | 0.05 | [+1.574759e−6, +3.650716e−6] | **[0, +7.301431e−6]** | ×2.318 |
> | 0.10 | [−6.212120e−6, −2.796612e−6] | **[−1.242424e−5, 0]** | ×2.221 |
> | 0.20 | [−2.803165e−5, −2.120522e−5] | **[−5.606329e−5, 0]** | ×1.322 |
> | 0.30 | [−5.324901e−5, −4.173248e−5] | **[−1.064980e−4, 0]** | ×1.276 |
> | 0.50 | [+4.226952e−5, +5.411263e−5] | **[0, +1.082253e−4]** | ×1.280 |
>
> **The ±100 % band does NOT bracket the PDF spread at the two lowest x**:
> ×2.318 and ×2.221 both exceed 2, so running `{0, 1, 2}` on one backend does
> not cover the others at x ≤ 0.1. It does at x = 0.20/0.30/0.50. **Both bands
> have to be quoted.** And the ±100 % band is not lifted by any of this: it
> comes from Q(⁶Li) = −0.0818(17) fm² against Q_d = +0.2859(3) fm², **not**
> from the A = 2 gate, and **the gate is A = 2 and tests nothing about the
> α–d step** — there is no b₁ measurement at any A > 2. `∫b₁ dx` over
> [0.01, 1.2] becomes **[+1.387049e−4, +1.458804e−4]** (a 5.17 % spread, the
> tightest quantity on the axis).
>
> **The 2026-09-03 lift is still a statement about a configuration**, and the
> `toy` column still sits outside it: `ToyF2`'s own G3b is **0.440**, outside
> the [0.5, 2] window the lift was read off, and the gap is a shape change and
> not a normalisation (mstw/toy on `Li6ConvolutionB1::b1(x, 2.5)` =
> **1.848 / 1.276 / 0.817** at x = 0.10 / 0.30 / 0.50; ×1.92 on the gate's own
> peak) — so no rescaling converts one column into another, which is why the
> rerun was a rerun and not a factor. **Two published readings do NOT survive
> it**, and they are named in place below: the orbital sector's *opposite
> sign*, and *"w = 2 nearly zeroes b₁ at x = 0.10"*. `docs/USAGE.md` §2a and
> `run_2026-09-02/phase_D_numbers.md` carry the same band and the same
> reasons; one row is pinned at rtol 1e-12 by
> `python/tests/test_li6_unpol_band.py`.

*(The `Li6B1(CdksB1)` column DOUBLED on 2026-09-03 — `b1_convolution()` stopped
halving a column that is already per nucleon; see Q1b below. Every other column
is unchanged, because `Li6ConvolutionB1` already reached the raw column through
`cdks_b1_raw_per_nucleon()`.)*

| x | (1) embedded d, S | (3) CG D-wave | (2d) struck d | (2α) struck α | **total** | `Li6B1(MillerB1)` | `Li6B1(CdksB1)` |
|---|---|---|---|---|---|---|---|
| 0.05 | +1.5686e−6 | +3.0979e−9 | +1.3838e−6 | +6.9528e−7 | **+3.6507e−6** | +3.6503e−4 | +1.7983e−6 |
| 0.10 | −4.6926e−6 | −9.2965e−9 | +1.2673e−6 | +6.3805e−7 | **−2.7966e−6** | +4.1222e−4 | −5.3132e−6 |
| 0.20 | −2.3664e−5 | −4.6736e−8 | +1.6669e−6 | +8.3886e−7 | **−2.1205e−5** | +2.6700e−4 | −2.6982e−5 |
| 0.30 | −4.4967e−5 | −8.8342e−8 | +2.2133e−6 | +1.1092e−6 | **−4.1732e−5** | +6.7167e−6 | −5.1928e−5 |
| 0.50 | +4.2448e−5 | +8.5642e−8 | +6.1632e−6 | +3.0422e−6 | **+5.1739e−5** | −1.8650e−4 | +4.5585e−5 |

Reading it:

* **Term (3) is 0.197 % of term (1)** at every x — exactly the CG algebra's
  0.1 × P_D^{αd}. It is kept because it is free and because it is the analytic
  bridge to `LI6_B1_RANK2_TRANSFER`, not because it matters numerically.
* **Terms (2d)+(2α) are 0.2–133 % of term (1)** over the band, and they are
  **the only part of these tables the unpolarised backend moves at all**
  (terms (1) and (3) convolve b₁ᵈ and are bit-identical across it). **This is
  the ⁶Li quadrupole puzzle showing up in b₁ and it is the quantitative
  reason the 100 % band is mandatory** — the PDF systematic and the A > 2 band
  live on the *same* term.
  > **⚠ The "OPPOSITE sign … over the whole window" this bullet published
  > until 2026-09-06 DOES NOT SURVIVE THE BAND, and is withdrawn as a
  > reading.** The signed ratio [(2d)+(2α)]/(1) at
  > x = 0.05/0.10/0.20/0.30/0.50/0.70 is
  > **+1.325 / −0.406 / −0.106 / −0.074 / +0.217 / +0.093** on `toy`
  > (opposite at three of the six, which the toy row's own signs already
  > showed), **+0.002 / +0.322 / +0.183 / +0.148 / +0.273 / +0.167** on
  > `ct18nlo` — *same* sign as term (1) at **every** x — and
  > **+0.272 / +0.099 / +0.147 / +0.182 / −0.006 / +0.166** on `mstw`, where
  > the one negative entry is a near-zero crossing, not an opposition. The
  > follow-on *"they nearly cancel term (1) at x ≈ 0.1 and dominate below"* is
  > toy-only: at x = 0.10 the total is **0.596** of |term (1)| on `toy` but
  > **×1.324** on `ct18nlo` and **×1.101** on `mstw` — the orbital sector
  > **adds** there. **The sign of the orbital sector relative to term (1) may
  > not be quoted without the unpolarised backend.**
  > `run_2026-09-06/phase_A_li6_tables.md` §A4.3.
* **The struck-α term is half the orbital sector**: (2α)/(2d) = 0.494–0.503
  on `toy` and **0.477–0.512 over the whole band** against the analytic
  2(M_d/M_α)² = 0.5064 (T15) — configuration-robust wherever it is readable.
  Dropping it is a 30 % error. (Two band entries read −1.185 and 13.08, at
  `ct18nlo`/x = 0.05 and `mstw`/x = 0.50; those are ratios taken across a zero
  of (2d) — ×0.0121 and ×0.0030 of its `toy` value — and are **not** a spread.)
  *(On the design's 600/2400/800 y grid this column read 0.489–0.497 and term
  (2d) itself was 2.3 % high at x = 0.10 — a quadrature artefact of the struck
  deuteron's wider z-distribution, fixed by the 2400/2400/3200 default and now
  guarded by a ⁶Li clause in T2.)*
* Against the production default `Li6B1(MillerB1)` the new model differs by two
  orders of magnitude **and by sign**, because Miller (HERMES-like) and CDKS
  (convolution) are different *camps* for b₁ᵈ — a pre-existing disagreement this
  work does not resolve. **The default does not change.**

Densities (struck d): raw `norm()` **0.8176766** against the identity
N_αd⟨E⟩/M_d = 0.817421 (3.1e−4 apart, T3's own tolerance); `mean_y()` =
⟨z²⟩/⟨z⟩ **0.9997871** (⟨z⟩ ≈ 1, *not* 1/3 — the fair-share normalisation of
design §1.6 is what makes this literally CDKS Eq. (16) one level up); `p_d()` 0.0193299, `p_d_momentum()` 0.0193516
against `VMC_P_D_LI6` = 0.0193549 (1.7e−4, the file's own quadrature spread).
∫b₁ dx over [0.01, 1.2]: **+1.4588e−4** on `toy` and
**[+1.387049e−4, +1.458804e−4]** over the band (5.17 %), against
`Li6B1(MillerB1)`'s +9.1243e−4
and `Li6B1(CdksB1)`'s +1.38596e−4 (that one doubled on 2026-09-03; it was
+6.9298e−5). *(The densities above are **backend-independent** — measured
identical on all three, since `LightConeDensities` never sees an `UnpolSF` —
so they were never a ToyF2 measurement and did not need the rerun.)*
`finite_q_delta` moves ⁶Li's x·b₁ by only
−1 % to +7 % on `toy` and **−1.1 % to +8.6 % over the band** (against
1.57–1.69 on the A = 2 gate), because in ⁶Li the term κ
multiplies is the small orbital one — which is why the two objects default
differently.

**The ⁶Li backend and the gate use different R — decided 2026-09-03, and both
defaults stay.** `Li6ConvolutionOptions::r_func` null ⇒ **`r_sigma_lt`** (the
kernel's own toy R); `DeuteronConvolutionB1::Options::r_func` null ⇒ **`r1998`**
(CDKS's SLAC world fit, because the gate reproduces *their* figure). So the R
the gate was validated with is not the R the shipped ⁶Li numbers above carry.

**What decides it is that the shipped observable is a RATIO.** The tensor
weight is K/D_φ with K ∋ b₁ and D_φ ∋ F₁, and that F₁ is `InclusiveKernel`'s,
built from the `ToyF2` the pipeline *shares* with `Li6ConvolutionOptions::unpol`
— whose R is `r_sigma_lt`, because nothing sets it. Give the b₁ backend a
different R and the (1 + R) in the numerator stops cancelling the one in the
denominator: (1 + r1998)/(1 + r_sigma_lt) at Q² = 2.5 is 1.088 at x = 0.1 and
1.039 at x = 0.3, the same size as the effect being chased. Fidelity to CDKS is
worth having on the *gate*, which has no denominator; inside the pipeline
consistency wins. **This is an author decision, not a derivation** — the third
option, and the better one, is to thread ONE R hook through
`default_inclusive_kernel` into both objects so they cannot disagree by
accident. **That wiring EXISTS since 2026-09-06, as an opt-in whose default
does not take it**: `--r-source {unset, sigma-lt, r1998}`
(`PipelineConfig::r_source`, `RSource`) puts ONE R object in
`Li6ConvolutionOptions::r_func` *and* `InclusiveKernel::Options::r_func`;
`unset` is the default and never enters that branch, so nothing here moved,
and the flag is refused off `--b1-model li6-convolution`. Its cost is priced
in `docs/open_items/run_2026-09-06/phase_A_numbers.md` §A1 and **it does not
simply undo the mismatch**: at y = 0.5, Q² = 2.5 the shared `r1998` moves
K/D_φ (and A_zz, which is −(2/3) of it) by −3.8357 % / +26.3988 % / +3.3518 %
at x = 0.05 / 0.10 / 0.30 against the numerator-only −5.5148 % / +24.6099 % /
+2.6629 % — **larger** at x = 0.10 — the cos 2φ amplitude moves by
−8.33 % / −6.73 % / −3.14 % where the numerator-only swap leaves it
bit-identical, the shared shift is **y-dependent** where the numerator-only
one is not, and the unpolarised rate moves too (σ −0.684681 %).

**Measured cost of the choice** (Q² = 2.5, default options), x·b₁:

| x | `r_sigma_lt` (the default) | `r1998` | change |
|---|---|---|---|
| 0.05 | +3.650716e−6 | +3.449387e−6 | −5.5 % |
| 0.10 | −2.796612e−6 | −3.484857e−6 | +24.6 % |
| 0.20 | −2.120522e−5 | −2.235704e−5 | +5.4 % |
| 0.30 | −4.173248e−5 | −4.284379e−5 | +2.7 % |
| 0.50 | +5.173934e−5 | +5.075777e−5 | −1.9 % |

**And it is not the (1 + R) prefactor that does it.** Term (1), the embedded
deuteron, is **bit-identical** under the swap — it carries b₁ᵈ from the injected
`TensorSF`, the raw digitized column, which has no R in it at all. The whole
effect sits in the two orbital terms, whose F₁ᵈ slot *is* `f1_cdks`, and there
it is **×0.54 to ×0.94**, far more than the ≤ 8 % a prefactor allows: those
terms convolve F₁ᵈ(x/z) against a density that **integrates to zero**, so they
respond to the *slope* of R — and `r_sigma_lt` is x-independent by construction
while `r1998` runs from 0.30 at x = 0.05 to 0.20 at x = 0.5. The +24.6 % at
x = 0.10 is that ×0.64 seen through the near-cancellation between term (1) and
terms (2d)+(2α). On the **gate**, where there is no such cancellation, the same
swap is worth **+2.7 %** on G3b (0.8432 → 0.8662, MSTW at Eq. (21), max over a
0.01 grid on [0.10, 0.80]) and, at κ = 1, moves the second zero from 0.392 to
0.365 — *away* from the digitized 0.457. **The gate's verdict does not rest on
this choice either way.**

**The truncation is MEASURED, not asserted.** The dropped P₂ and P₄ remainders
of the D-wave b₁ᵈ weight are **below 6e−5 of term (1) everywhere** (T14) — four
orders of magnitude inside the band. The reason is *not* "P_D × 1/10" (the P₄
coefficient is 15× the isotropic one term (3) keeps): they are small because the
P₂- and P₄-weighted densities integrate to zero, so they enter only through the
*curvature* of b₁ᵈ(x/z), and because they multiply b₁ᵈ, not F₁ᵈ. The other
dropped piece — the S–D interference in the b₁ᵈ sector — is term (2d)'s SD
structure with coefficient exactly 1/3 of it (G0f, exact) and b₁ᵈ in place of
F₁ᵈ: ~3e−4 of term (2d).

### Systematics that are stated, not hidden

* **The mandatory 100 % band.** Q(⁶Li) = −0.0818(17) fm² against
  Q_d = +0.2859(3) fm²: the α–d D wave enters the closest measured observable
  with the opposite sign to the deuteron's own D state and nearly cancels it.
  *(Q(⁶Li) is `LI6_QUADRUPOLE_FM2` in `include/lipolgen/rc.hpp`, TUNL's A = 6
  evaluation, 1998CE04; Pyykkö's compilation, Mol. Phys. **106** (2008) 1965,
  gives −0.0806(6) fm² from Cederberg et al., Phys. Rev. A **57** (1998) 2539
  — the repository constant is TUNL's. Q_d: Bishop & Cheung, Phys. Rev. A
  **20** (1979) 381.)* Every number is `--b1-band-scale 0/1/2`,
  and `--b1-alpha-d-dwave-weight 0/1/2` is the shape variant reported next to
  it.
* **±5 % on N_αd, unexplained — and now PROPAGATED (2026-09-04, §C5.3).**
  Three tabulations of the α–d spectroscopic factor span 5 %: **0.819481**
  (2014 `li6_ad1.momentum`, the default, and now the single home
  `VMC_N_ALPHA_D_LI6` that `VMC_P_D_LI6` is expressed through), 0.856 (2004
  `li6.ad`) and 0.863 (Wiringa et al., PRC **89** (2014) 024305 §III:
  0.846 + 0.017). Entries 1 and 3 are the same year and the same Hamiltonian
  family, so this is **not** a version difference anyone can name. b₁ is
  exactly linear in `norm_target`, so it is a flat ±5 % (T12).
  **Three numbers span three different amounts, and "N_αd 5 %, P_D 7 %" was
  quoting two of them about a third**: N_αd spans **5.311 %**, the D-wave
  *norm* spans **7.181 %**, and P_D — the ratio every tagged observable
  actually uses — spans only **2.729 %**. Outside b₁ the propagation is nearly
  nothing: `TaggedModel` renormalises each wave to its own P_L, so N_αd cancels
  exactly and the three readings move `tensor_dilution` by **0.048 %** and
  `vector_dilution` by **0.082 %** (T26, `phase_C_numbers.md` §C5.3).
* **Q4 — suppression or renormalisation?** The default gives the 18 % of ⁶Li
  that is not α+d **b₁ = 0**; `use_spectroscopic_factor = false` renormalises
  instead and is ×1/N_αd = **1.2203** on the whole answer. Both are defensible,
  nobody has published either, and the difference is well inside the band.
* **Q1b — RESOLVED on the CDKS half, an author decision on the Miller half
  (2026-09-03).** The two camps no longer share one applied constant.
  * **The arbiter is HERMES, not either theory paper.** Their published b₁ᵈ —
    the data both camps plot against — is **per nucleon**, because their Eq. (5)
    divides by an F₁ᵈ built from F₂ᵈ = (F₂ᵖ + F₂ⁿ)/2. Inverting their own
    Table II reproduces that F₁ᵈ in all six bins (mean ratio **0.946**, against
    **0.473** for the per-deuteron one) across five decades of b₁ᵈ and a factor
    9 in Q². *(Their own QPM definition table, read literally, is the
    per-deuteron object — and that inconsistency is exactly the trap that
    caught the two theory camps in opposite directions.)*
  * **CDKS: certain, and `CdksB1` stopped halving.**
    `B1_CDKS_TABLE_TO_PER_NUCLEON` = **1**. Eq. (10)'s spectral function carries
    an explicit 1/A, the text under Eq. (16) says in words *"the structure
    function b₁ is defined by the one per nucleon"*, f(y) is normalised to one
    nucleon and F₁ᴺ = (F₁ᵖ + F₁ⁿ)/2. **`b1_convolution` / `CdksB1` /
    `--b1-model cdks` therefore DOUBLED.** Nothing on the default path moved:
    the only b₁ reference JSON records `"b1_model": "miller"`, so the rtol-1e−12
    gate is untouched, and `Li6ConvolutionB1` already reached the raw column
    through `cdks_b1_raw_per_nucleon()`.
  * **Miller: likely, not certain — the 0.5 stays, as an author decision.**
    `B1_MILLER_TABLE_TO_PER_NUCLEON` = **0.5**, and it *is*
    `B1_PER_DEUTERON_TO_PER_NUCLEON` rather than a second copy of the number.
    His Eq. (1) densities are "in a target hadron", his Eq. (5) is a light-cone
    correlator in the normalised deuteron state, and Eq. (6)'s ½ is fully
    consumed by the quark-spin average of a *spinless* pion — every ½ in
    Eqs. (1)/(5)/(6)/(20) is spoken for and none of them is a 1/A. **Against
    that**, his Table I transcribes HERMES's per-nucleon numbers unrescaled,
    his Fig. 5 overlays them on the curve, and he tunes P₆q to one of them at
    face value. The paper is self-inconsistent by exactly this factor, and no
    third party adjudicates. Keeping the 0.5 is the status quo, moves no
    default and touches no gate; **the numerical consequence of the other
    reading is that every `MillerB1` number would double** — including the
    default inclusive tensor rate and its rtol-1e−12 reference. What would
    settle it: ask Miller, or evaluate his Eq. (20) at x = 0.012 and see
    whether it lands on 10.5 × 10⁻² (per nucleon, the 0.5 is wrong) or
    5.25 × 10⁻² (per deuteron, the 0.5 is right) — which needs the pion PDF set
    of his ref [29], model 1, that this tree does not have.
  * **A retraction while here.** This document previously said the answer to
    Miller's normalisation "is on the axis of Miller's Fig. 5". It is not: that
    ordinate is `100 b₁(x)` and nothing else, his caption is bare, and the
    strings "per nucleon" and "1/A" do not occur anywhere in the paper.
  * Full argument, with the six-bin inversion and its reproduction script:
    `docs/open_items/run_2026-09-03/phase_A_miller_normalisation.md`. One
    consequence to carry: `close_kumano_integral(false)` and
    `close_kumano_integral(true)` integrate the RAW columns and are therefore on
    **different scales** — never compare them to each other.
* **Q6.** `LI6_B1_RANK2_TRANSFER` = 0.921947 is tied to the Hulthén-scenario
  P_D = 0.0867; the VMC value 0.0193549 would give 0.982581. Changing it would
  move every published inclusive tensor number, so it stays — but this backend
  makes the inconsistency visible (the tagged channel already uses VMC at
  `--cluster-wave vmc`; the inclusive b₁ never does).
* **Q8 — there is nothing to validate against for A = 6.** No published b₁
  exists for any A > 2. The A = 2 gate plus the analytic limits (T5, T6, T7) are
  the *only* checks there are, which is why the gate is blocking and why the
  band is not a formality. **The gate closing on 2026-09-03 does not change
  this**: it says the *kernel* reproduces CDKS on the deuteron, and says nothing
  about the α–d step it is applied to one level up.

### The seven close conditions, and where each one landed

All seven were worked on 2026-09-03. Five are **done**, two are recorded as
**author decisions** — a decision, with its evidence and its measured cost, is
a close and not a deferral.

| # | condition | outcome |
|---|---|---|
| 1 | Decide item 0's default (the δ-function) | **DONE.** The A = 2 gate is quoted at CDKS Eq. (21), their exact definition; the ⁶Li backend keeps Eq. (17)'s κ = 1, where the switch is worth 1 %. It turned out to be load-bearing for the verdict, not only for the third digit: MSTW gives 0.843 at Eq. (21) against 0.520 at Eq. (17) |
| 2 | Get MSTW2008 LO (or Kumano's tabulated curves) and rerun checklist item 4 with CDKS's own PDF | **DONE.** MSTW2008 LO was on disk all along — the CENTRAL member ships with PYTHIA 8 as `pdfdata/mstw2008lo.00.dat`; only LHAPDF's set store lacked it. `MstwSF` (`include/lipolgen/mstw_sf.hpp`) reads it with the same charge weights as `LhapdfSF`, so an MSTW/CT18NLO ratio is a PDF comparison and not a convention one. **G3b 0.719 (CT18NLO stand-in) → 0.843243.** Kumano's curves were not needed |
| 3 | Get a real CD-Bonn u, w instead of the D-state rescaling proxy | **DONE.** Machleidt's Appendix-D parameterisation is `cdbonn_wave()` in `cluster.hpp`, with C₁₁ and D₉..D₁₁ **computed** from the r → 0 boundary conditions rather than typed; reached through `DeuteronConvolutionB1::Options::wave = kCdBonn`. It reproduces his Table XV (P_D 4.8562 % against 4.85, η 0.0255714 against 0.0256(4), Q_d +0.270178 fm² against 0.270) and is **×1.18631** on the G3b peak — it **closes**, which is the opposite sign to the retired proxy's ×0.88. The default stays AV18 |
| 4 | Decide whether `Li6ConvolutionOptions` should default to `r1998` like the gate | **DECIDED — it stays `r_sigma_lt`**, because the shipped observable is a ratio whose denominator carries `InclusiveKernel`'s R. Measured cost of the choice and the term-by-term reason are above ("The ⁶Li backend and the gate use different R"). The better third option — one R hook threaded through both — **is now BUILT and PRICED as the opt-in `--r-source`, and still not taken as a default** (2026-09-06): registry row 3 is open, and the question it now asks is *(iii) with which R?* (`docs/open_items/run_2026-09-06/phase_A_numbers.md` §A1) |
| 5 | Decide whether `CdksB1` should stop halving the per-nucleon column (Q5(a)) | **DECIDED — yes, and done.** `B1_CDKS_TABLE_TO_PER_NUCLEON` = 1; `--b1-model cdks` doubled. Certain, from CDKS's own words plus the HERMES inversion (Q1b above). Two pinned values in `tests/test_sf.cpp` doubled with it; no reference JSON moved |
| 6 | Check **Miller's** own normalisation independently (Q5(b)) | **DECIDED AS AN AUTHOR DECISION — the 0.5 stays, and it is *likely*, not certain.** His derivation is per deuteron and his presentation is per nucleon; the paper cannot be both. Evidence for each reading and the numerical consequence of the other (every `MillerB1` number, and the default inclusive tensor rate, would double) are in Q1b above. **The framing this document gave the question was wrong and is retracted**: the answer is not on the axis of his Fig. 5 |
| 7 | Only then re-open G3b | **RE-OPENED, AND IT PASSES.** Ratio **0.843243** with MSTW at Eq. (21); **1.000338** with CD-Bonn as well. Design §5.4's *Escalation* clause is not triggered, so **the ⁶Li publication ban is lifted** — read the four conditions in the warning block at the top of this section before quoting anything |

**What is still open, stated as what it is.**

* **The Miller normalisation (condition 6)** — `STATUS.md` decision **row 2**,
  drafted as `run_2026-09-03/AUTHOR_DECISIONS.md` §B12 — is a decision, not a
  measurement.
  Closing it properly needs the author or a numerical evaluation of his
  Eq. (20), which needs a pion PDF set this tree does not have.
* **The R wiring (condition 4)** — `STATUS.md` decision **row 3**,
  `AUTHOR_DECISIONS.md` §B7 — is a decision too; the clean fix is one R hook
  threaded through `default_inclusive_kernel` into both objects. **It was
  unpriced until 2026-09-06 and is not any more**: the hook is built as the
  opt-in `--r-source` (default `unset`, nothing moved) and measured against
  both alternatives in
  `docs/open_items/run_2026-09-06/phase_A_numbers.md` §A1. The decision is
  still open — sharing the hook does **not** shrink the shift (it is larger at
  x = 0.10), it makes the shift **y-dependent**, and it moves the unpolarised
  rate by −0.684681 % — so the row now asks *(iii) with which R?*
* **G3c is not enforced and did not improve with the PDF** (0.467 → 0.490 of
  the digitized integral; 0.98 only with CD-Bonn as well). The peak agreeing
  does not make the integral agree.
* **x = 0.50 is a factor 9 low** on the verdict row, for the zero-position
  reason recorded above. No clause of G3 sees it.
* **The ±5 % N_αd spread is still unexplained** (0.819481 / 0.856 / 0.863 from
  three tabulations, two of them the same year and Hamiltonian family) —
  unexplained, but since 2026-09-04 **propagated**: 5.311 % on N_αd, 7.181 % on
  the D-wave norm, 2.729 % on P_D, and only 0.048 % on the tagged tensor
  observable (§C5.3, T26). It is a real ±5 % in b₁ and a rounding error
  everywhere else.
* **Q6, Q4 and Q8 are unchanged** by this phase: `LI6_B1_RANK2_TRANSFER` is
  still tied to the Hulthén-scenario P_D, `use_spectroscopic_factor` is still a
  1.22 fork, and there is still nothing to validate an A = 6 b₁ against.
  *(2026-09-04, §C5.5: the size of the Q6 gap is now measured. Under
  `--cluster-wave vmc` the tagged model's `tensor_dilution` is 0.982576 against
  `LI6_B1_RANK2_TRANSFER` = 0.921947, i.e. **+6.58 %**, and the vector sector
  drifts **+11.61 %**. Deliberately NOT closed by substitution: the
  "consistent" reading of the inclusive polarization, 0.905427, sits **6.8 %
  above** the ab-initio six-body VMC 0.848 where the shipped 0.811228 is 4.3 %
  **below** it, so the formula carries the error, not the choice of P_D.)*
* **The default wave function is still AV18**, deliberately: making `kCdBonn`
  the default of `DeuteronConvolutionB1` is a strong argument — that object
  exists to reproduce CDKS's figure and CDKS used CD-Bonn — but it would move
  every pinned gate number, so it belongs to a task that says so.

## 11. Coherent ⁶Li amplitude — **ANSWERED AS A BAND (MARGINAL)** (rows 8, 12, 16, 21 open)

Sartre: nuclei hard-coded (`Nucleus.cpp` switch, no A = 6), spherical sampling
(no polarization axis), tables CPU-years (the 2026 speed-up code is unreleased) —
ruled out. Outdated project claims: 2605.00454 publishes coherent J/ψ down to
A = 3, 4 with α-clustering; 2511.05638 has e+⁷Li coherent J/ψ |t| distributions
with tagging efficiency. Path, both halves now **in-tree**: (11.1) eSTARlight
⁶Li for the unpolarized rate and slope baseline; (11.2) an α+d configuration
sampler (α core from the VMC density, p–n pair from AV18 u/w oriented by the
polarization axis, α–d separation from the VMC α–d overlap) whose output grafts
into `hejajama/subnucleondiffraction` (author asks to be contacted) for the
tensor cos 2φ. Warn in plans/06: photon-polarization cos 2φ (STAR 2204.01625)
is a distinct mechanism and a background at Q² > 0.

### 11.1 eSTARlight — the unpolarized rate and slope baseline

**In-tree since 2026-09-02.** eSTARlight (`github.com/eic/estarlight`, commit
`939b11a24499398392d959db81c7502aeec91046`, same code arXiv:2511.05638 cites)
runs e+⁶Li today with no code change — `nucleus::init()` has no Z = 3 case, so
it falls to the generic light-nucleus branch (`_Radius = 1.2·A^{1/3}`, Gaussian
form factor). Full run log, input files and an adversarial re-check are
`docs/open_items/run_2026-09-02/estarlight_li6.md`; nothing there was
committed. 2×10⁵ events per channel, `0.1 < Q² < 100 GeV²` (arXiv:2511.05638's
own window), e 10 GeV × ⁶Li 99.5 GeV/u:

| channel | σ_coh (default R = 2.1805 fm) | B fitted [GeV⁻²] | σ_coh (measured R = 2.589 fm) | B fitted |
|---|---|---|---|---|
| coherent ρ⁰ | **506.4 nb** | 38.6 | 379.5 nb | 54.1 |
| coherent φ | **30.16 nb** | 39.2 | 22.08 nb | 54.8 |
| coherent J/ψ | **1.773 nb** | 38.9 | 1.255 nb | 55.0 |

The fitted slope is VM-independent to 1.5 % (there is no VM-dependent slope in
a coherent Gaussian-form-factor calculation) and reproduces the analytic
B = R_G²/(3ħc²) to **3.7–5.2 %**, one-sided (φ −3.69 %, J/ψ −4.43 %, ρ −5.17 %
against 40.7026 — corrected 2026-09-04 from "4 %", which is the φ row alone;
pinned in `tests/test_coherent.cpp` T10a). The ⁷Li cross-check against arXiv:2511.05638 (identical
configuration, same commit-era code, same Q² window, same beam) agrees with
every *shape* statement the paper makes — no diffractive minimum, the |t| slope
and the Fourier-transform width read off its Figs. 6–7 — but the paper quotes
no absolute σ, so no number-to-number check is possible.

**Conclusion for `CoherentScenario::slope_b`.** `gaussian_slope(r_rms)` in
`coherent.cpp` (B = R_rms²/(3ħ²c²)) is algebraically the same object as
eSTARlight's light-nucleus form factor; only the radius differs. The scenario
band `slope_b ∈ {40, 60}` GeV⁻² corresponds to R_rms ∈ {2.162, 2.647} fm, which
brackets almost exactly the two defensible ⁶Li densities — eSTARlight's own
`1.2·A^{1/3}` (fitted **38.6–39.2**) at the bottom and the Angeli–Marinova
measured charge radius 2.589 fm (fitted **54.1–55.0**) at the top, i.e. a
citable **B = 39–55 GeV⁻²** band from simulation rather than a hand-tuned one.
**`slope_b = 50 ± 10` survives contact with eSTARlight unchanged, and the band
could be tightened to 45–58 GeV⁻² but there is no reason to** — the ±30 % rate
swing between the two radii shows the density assumption, not the dynamics, is
the leading systematic.

**Conclusion for `CoherentScenario::f0`.** eSTARlight generates *exclusive*
vector-meson production only — M_X is the meson mass exactly, with no
diffractive continuum — so **it cannot measure f0**, LiPolGen's coherent
fraction of the DIS rate at x → 0. What it *can* bound: running ⁶Li at
`MIN_GAMMA_Q2 = 0.7` against `generate_inclusive`'s σ_incl = 592 nb in the same
window gives σ(ρ+φ+J/ψ)/σ_incl = 3.0 × 10⁻² (all three channels; f0 = 0.04 is
not refuted at the order-of-magnitude level, but 90 % of this is ρ⁰, which sits
below LiPolGen's own M_X ≥ 1.2 GeV coherent floor) and
σ(J/ψ only)/σ_incl = **1.0 × 10⁻³** — the one channel actually inside
LiPolGen's window, and therefore a genuine **lower bound**, 40× below f0 = 0.04,
on the part of coherent diffraction eSTARlight's exclusive-VM model can see.
Determining f0 itself needs a coherent diffractive-DIS calculation (a
coherent-A analogue of the H1/ZEUS diffractive PDFs) that exists in neither
eSTARlight nor Sartre; **recommendation: leave `f0 = 0.04` and its {0.02, 0.08}
band as a scenario**, and do not treat eSTARlight's numbers as a
determination. **Done 2026-09-03**: the comment on `CoherentScenario::f0`
(`include/lipolgen/coherent.hpp`) now records a **lower bound of 1.0 × 10⁻³**
(exclusive J/ψ inside the M_X ≥ 1.2 GeV window). The all-VM figure 3.0 × 10⁻²
is *not* an upper bound on f0 — ~90 % of it is ρ⁰, which with φ sits below that
floor.

Caveats carried forward unchanged from `estarlight_li6.md` §6: one spherically
symmetric Gaussian density (no α+d clustering, no polarization axis, no
diffractive minimum — not an imaging baseline), and no saturation.

### 11.2 The α+d configuration sampler — the tensor cos 2Φ half

**In-tree since 2026-09-03**: `include/lipolgen/cluster_config.hpp` +
`src/core/cluster_config.cpp`, `lipolgen-configs`,
`tests/test_cluster_config.cpp`, `python/tests/test_cluster_config.py`, design
`docs/open_items/run_2026-09-02/design_G_cluster_config.md`, numbers
`.../phase_G_numbers.md`. **Opt-in and inert**: nothing in the generator calls
it, `CoherentScenario` is untouched, and its output is a file. It writes
`he3.dat`-compatible position tables (fm, ion rest frame, c.m. at the origin,
the polarization axis applied) plus a sidecar carrying the axis, the substate,
every option and every input table's md5. `docs/USAGE.md` §9 has the CLI,
the Python API and the output format.

**What it is.** ⁶Li(1⁺) as a rigid α(0⁺) core plus a deuteron with relative
L = 0, 2 coupled to S = 1: the α's four nucleons come from the ANL VMC ⁴He
one-body density (recentred so the configuration's own c.m. sits at the
origin, exactly), the α–d separation R and its orientation from the ANL VMC
α–d overlap (S+D wave, `li6.ad`/`li6.adr.fit`), and the p–n pair from AV18
u(r)/w(r) — all three sharing one drawn deuteron projection m_S per
configuration, so the R̂-to-r̂ correlation the design calls out is kept and
only the *relative azimuth* of the two is dropped (`tagged.hpp`'s own
truncation, in r space).

**The three m-state densities.** One (R, cos θ_R) table per (m, m_S) pair —
nine tables total, the same discipline `TaggedModel::build_amp2` uses in
k-space — because sampling R̂ from the m_S-*summed* density instead would
decorrelate it from the deuteron spin a moment test cannot see: for m = +1
the exact ⟨P₂(cos θ_R)⟩ is −2/7 at m_S = −1 and +1/7 at m_S = 0 —
reproduced by the sampler's grid to 1e-6 for m_S = −1 and to 2.3e-4 at the
default 96 cos θ cells for m_S = 0 (1e-6 at n_c = 3072: a midpoint-rule
residue at the simple zero of |Θ₂¹|², `phase_G_numbers.md` §3) — while an
m_S-marginal draw would give the m_S-**summed** ⟨P₂⟩ = **−0.0370** in
*every* branch (the design's own earlier estimate of it is −0.0355), which
the two conditioned branches sit 46σ and 17σ away from in a
2×10⁵-configuration test.
The three substates m = +1, 0, −1 (plus the interleaved unpolarized mix) all
read the *same* radial tables; only the Clebsch–Gordan recoupling changes.

**The quadrupole puzzle, and how the sampler handles it rather than hides
it.** This geometry's own point-matter quadrupole overshoots the measured
one by a factor ≈ 7.5: Q_charge(⁶Li) = **−0.615 fm²** (model range
−0.615…−0.730 across the three α–d sources) against the measured
**−0.0818 fm²** (`LI6_QUADRUPOLE_FM2`) and GFMC AV18+IL7's **−0.20(6) fm²**
(Pastore *et al.*, PRC 87, 035503 (2013)) — an independent-Hamiltonian
comparison point, not a check on the same tables. The asymptotic α–d D/S
ratio, correctly divided by the Coulomb Whittaker ratio (which is 3.3 at
R = 6 fm and 2.6 at R = 8 fm — the naive R₂/R₀ is *not* the asymptotic ratio),
gives **η = −0.048** against the measured **η = −0.025 ± 0.006 ± 0.010**
(George & Knutson, PRC 59, 598 (1999)): a real but *moderate* ≈ 2× D-wave
excess, not the 5–15× a naive ratio suggests — which demotes an ANL
normalization/phase-convention error as the leading explanation (design O1).

**The 7.5× decomposed, and its error bar (added 2026-09-03; closes O1).**
Numbers and reproduction in
`docs/open_items/run_2026-09-03/phase_C_numbers.md` §C1, pinned in
`tests/test_cluster_config.cpp` **T22b** and its pytest mirror. The (G9) dial
scales R₂ in place and `asymptotic_ds_ratio()` reads the mutated waves after,
so **η(s) = s·η(1) exactly** — which lets a *measured* observable, not the GFMC
number, supply the first leg:

| leg | ratio | value |
|---|---|---|
| model → the η-matched dial (Q = −0.18557 fm², η = −0.02507 ≈ GK) | −0.615448 / −0.18557 | **3.3165288362** |
| η-matched → measured (`LI6_QUADRUPOLE_FM2`) | −0.18557 / −0.0818 | **2.26858190709** |
| **product** | −0.615448 / −0.0818 | **7.52381731215** |

**Two factors, not three.** The missing ≈15 % enters as 1/S_αd = 1/0.8542 =
1.1706, and that is **already inside the −0.615**: `cluster_config.cpp:369-374`
divides *both* waves by √S_αd before the moments are recomputed. Writing it as a
third factor gives 8.808, which overshoots the measured 7.524 by **17 %**.
1.17× is a **ceiling** on what a coherent missing-component model could add, not
a multiplier.

**The band is the physics, not the 3.32.** σ_comb(GK) = 0.011662, and because η
is linear in the dial that maps straight onto model Q: **−0.4005 fm²** at one
edge (1.54× the model), **+0.0298 fm²** at the other — the interval passes
through **zero at η = −0.014971, 0.86 σ from the central value**. So the leg is
*3.32× (1 σ: 1.54× … sign change)*. Run backwards, η implies −0.019439 for the
measured Q (**+0.48 σ**) and −0.025854 for GFMC's (**−0.07 σ**): both inside
GK's 1 σ, so **η cannot discriminate between them**. What O1 gets is the
qualitative result only — the D-wave excess is moderate *and measured*, so the
residual 2.27× is cluster polarization plus the missing component.

The sampler does not paper over the gap: `quadrupole_band_fm2()` returns all
three numbers together and the writer stamps them on every output, and
**the rule stands — do not derive a published tensor input from this
geometry alone.** A `quadrupole_target_fm2` **deformation dial** (not a
wave-function fit) can rescale the α–d D-wave to land the geometry's own
Q_charge on any target in the reachable range [−0.615, +0.270] fm²; the root
for the measured Q(⁶Li) is s = 0.4032, which drags P_D(α–d) down to
3.3 × 10⁻³ from 0.02011 — the dial trades away the natural D-state
probability to match the tensor moment, and the sidecar always records which
choice was made.

**The measured numbers** (default `FitRescaled` source; full table, all
three sources and the grid-vs-analytic quadrature study in
`phase_G_numbers.md`): ⟨r²⟩(⁶Li) = **6.4447 fm²** (r_rms 2.5386 fm — **3 %**
above the measured point radius 2.4655 fm, `LI6_R2_POINT_FM2`, and **4 %**
above the VMC `li6.density` value 2.4433 fm, `LI6_R_POINT_VMC_FM`); the
point-matter quadrupole per
substate, Q_matter(m) = (3m² − 2)·[(4/3) Q[R₀,R₂] + 2 Q_d D_T] (eq. (G5)),
is **Q_matter(±1) = −1.2309 fm²** and **Q_matter(0) = +2.4618 fm²** — the two
non-zero substates carry opposite sign by construction, and Q_matter(+1)/2 is
exactly the Q_charge(+1) quoted above. A 2×10⁵-configuration Monte Carlo
closure reproduces both to within 5σ of the sampler's own per-configuration
variance (sampled ⟨r²⟩ = 6.4483 vs 6.4447, 5σ = 0.041; sampled
Q_matter(+1) = −1.2216 vs −1.2309, 5σ = 0.322), and Σᵢ r⃗ᵢ = 0 to 2.2 × 10⁻¹⁵ fm
worst component over 20 000 configurations. Feeding Q_matter(+1) through the
closed-form (G7)+(G8) quadrupole → a₂ map — which reproduces the only
published polarized coherent calculation, Mäntysaari *et al.*'s digitized
deuteron a₂, to 8 % at m = ±1 and 8–21 % at m = 0 with **zero free
parameters** — gives a₂(±1) = **+0.1976** at |t| = 0.3 GeV² for this
geometry's own (overshot) Q, and **+0.026** for the *measured* Q(⁶Li) — about
10× smaller than the deuteron's and of the **opposite sign**, because
Q(⁶Li) < 0 while Q_d > 0. **Timing**: 0.54 µs/config in-process, 1.31 µs/config
through the `lipolgen-configs` CLI, both well inside the design's < 5 µs/config
target; a 10⁵-configuration table with its sidecar writes in ≈ 1 s.

**How the graft would consume the table.** The writer emits one line per
configuration in `he3.dat`'s own layout, in fm, ion rest frame, c.m. at the
origin, with the polarization axis already applied — so the consumer needs no
knowledge of our conventions beyond "these are the six nucleon positions."
Their side then needs the `-configfile/-configid` generalization described
below. Their Good–Walker loop is unchanged; the amplitude uses only x and y.
One set per m ∈ {+1, 0, −1} plus one unpolarized set; a₂ comes out of the Φ
dependence of ⟨A⟩ exactly as in their Fig. 2.

**The ask itself is reconciled into ONE canonical draft** —
`docs/open_items/run_2026-09-03/mantysaari_collaboration_draft.md` (task C6,
2026-09-04) — which supersedes the paragraph that stood here (and the copy
in `design_G_cluster_config.md` §10, which disagreed with it by one number:
"3 %" here vs "4 %" there on the α+d model's point-radius agreement, both
correct against different reference radii, both now stated together in the
canonical file). **It is not sent; sending is the author's decision.** The
two editor's notes that stood here (2026-09-03 on the factor-7.5 decomposition,
2026-09-04 on O5) are folded into that file's own §§2–4 rather than repeated;
what they found, in one line each: the "factor ≈ 7.5" is two factors,
3.3165 × 2.2686, not three (§C1 below); and §11.3's O5 answer means the
letter's own request (c) — asking the group to run their polarised J/ψ setup
on our ⁶Li tables — asks for a channel that is **marginally measurable**
(the band of §11.3b–c: 2.63 σ at its low edge, 2.84–3.29 σ at its top, 3 σ at
8.3–13.0 fb⁻¹/u), so the canonical draft keeps that
request and moves it to **photoproduction**, Q² < 0.1, with the four things that are not established (no detection efficiency below Q² = 0.1, where 85 % of the coherent rate sits; no decay-lepton reconstruction efficiency anywhere in this tree, which is what leaves the band OPEN BELOW; the ⁷Li → ⁶Li efficiency substitution, which straddles 1 and makes the band's top a span; and the far-forward working point, unchosen — `cluster_config.hpp:525-534`) stated inside the ask. *(An intermediate version
said the channel was "blind" and withdrew the request; that came from one
lepton channel in one Q² window and is corrected in §11.3a.)*

**What remains open.** The tensor cos 2Φ prediction above (a₂(±1) ≈ +0.026 at
the measured Q, opposite sign and ~10× smaller than the deuteron's) is a
closed-form estimate from the geometry's own point-matter quadrupole — **it is
not yet a coherent-diffraction amplitude from a dipole-model run**, and it
cannot become one without the actual `subnucleondiffraction` graft (requests
(a)/(b) of the canonical draft above): the Good–Walker amplitude, its Φ-averaging and its own
statistical and saturation-model uncertainties are not reproduced by
`a2_from_quadrupole`, which only carries the target's quadrupole moment
through the deuteron's own published |t| dependence. **§11.3 now prices that
limitation**: a dipole-model run would have to move the signal **up by ×1.15**
to bring coherent J/ψ to 3 σ at 10 fb⁻¹/u (**down by ×1.31** to fall back to
2 σ, **up by ×1.91** for 5 σ) — the 4.0 / 6.7 / 16 quoted here until
2026-09-05 belonged to the **restricted** one-lepton-channel, one-Q²-window row
of §11.3a, which §11.3 withdrew. Also still open, in the
design's own numbering: ~~**O1**~~ **CLOSED 2026-09-03** — the 7.5× gap is
**two** factors, 3.3165 (model → η-matched dial, from the *measured* D/S ratio)
× 2.2686 (→ measurement), and the missing 15 % is *not* a third factor because
1/S_αd = 1.1706 is already inside the −0.615; the leg's own 1 σ runs from 1.54×
to a **sign change** in Q (see above and §C1 of the phase-C numbers);
~~**O2**~~ **CLOSED 2026-09-03** — the `li6.adr.fit` R₂ node **is real** (the
raw `li6.ad` R₂ over the fit's inner lobe is **3.3 σ below zero**; raw node
1.119 fm [1.089, 1.153] against the fit's 1.0648 fm) but it **does not decide
the default**: r < 1.5 fm is −0.09 % of q_int. `FitRescaled` **stays the
default by author decision** — the two sources differ by 8.7 % in Q inside a
band already a factor 7.5 wide, neither is closer to GK (−2.0 σ vs −2.4 σ), and
the raw block is rough at ≈1 σ per point (`phase_C_numbers.md` §C3); **O3
PARTIALLY BOUNDED 2026-09-04, §11.5 below** — how much an uncorrelated α core
(vs. GFMC ⁴He configurations with correlations) moves the incoherent/coherent
split is still unanswerable from this repository (no configuration table
exists here for any nucleus), but a genuinely correlated α → d+d overlap
already in the tree (`he4.dd`) bounds the SIZE of 4-body correlation on a
comparable position-space moment at +18 % (variance) / +9 % (rms) against the
uncorrelated-product null — a several-to-twenty-percent effect, not a factor
of several; ~~**O4**~~ **CLOSED 2026-09-04** — ΔB is now defined once, at the
`eps_b0` declaration (|F_m|² = exp(−|t|[B + ΔB_m cos2(Φ−Φ_S)]), ΔB_m = δ_m/2,
so a₂(m) = −(ΔB_m/2)|t| and **eps_b0 = δ_{±1}/B = 2ΔB_{±1}/B = −ΔB₀/B**: the
old label "ΔB₀/B" was off by a **sign**, not a factor 2), with
`CoherentScenario::delta_b_m()` as its single code home. And the band was
re-decided: `eps_b0` = −0.08 implies a ⁶Li charge quadrupole of **−0.9345 fm²**,
**11.42×** the measured −0.0818 and 1.52× even the α+d model's −0.615, so it is
a **deuteron** number; the honest ⁶Li band is **−(0.0070 … 0.0527)** at
`slope_b` = 50. **Author decision: the default stays** (it is a reference gate)
and the cost is recorded — every *generated* coherent tensor number is 11.4×
the measured-quadrupole expectation, and the positivity edge
`t_positivity_edge(−2)` = **0.245 GeV²** is a consequence of that oversized
`eps_b0` (it would sit at **2.80 GeV²** at the measured quadrupole) — the
ceiling `COHERENT_T_MAX_DEFAULT` = 0.2 itself rests on the **anchor range**
and does NOT move with this decision (D5, §11.6, `STATUS.md` decision row 12;
this sentence stated the dependency the other way round until 2026-09-05). `phase_C_numbers.md` §C4, T10b/T23;
~~**O5**~~ **ANSWERED 2026-09-04, and the answer is MARGINAL — as a BAND** — §11.3
below and `phase_C_numbers.md` §C2. Coherent J/ψ at one EIC year, over the
**whole** Q² range and with **both** lepton channels, gives **2.63 σ at the
band's low edge and 2.84 … 3.29 σ at its top** on the a₂ this section predicts,
with 3 σ at **8.3 … 13.0 fb⁻¹/u** — inside the {1, 10, 100} fb⁻¹/u band at both
ends. Whether the top crosses 3 σ is **not established** (§11.3c). The separation from the photon-polarisation cos 2φ
is **not** the obstacle: it is free, it survives Q² → 0, and the P_zz flip is
a 1.50× *gain*. a₂ ∝ |t| on an e^{−B|t|} sample really is only a **0.25 %**
modulation rather than the 2.6 % that a₂(0.3) suggests — but 0.25 % on
4.2 × 10⁵ events is 2.6 σ. *(The first pass of §11.3 read "NO for J/ψ,
0.75 σ, 160 fb⁻¹/u". That was computed from **one** lepton channel in
**one** Q² window and is withdrawn; §11.3a records what moved and by how
much.)* **This is the item that was said to decide whether the collaboration
is worth proposing; it no longer decides against it, and what it leaves open
is not statistics but an efficiency nobody has measured below Q² = 0.1.**

### 11.3 O5 — **MARGINAL: the a₂ is a 2.6 σ (low edge) to 2.8–3.3 σ (top) measurement in coherent J/ψ at one EIC year**

**Answered 2026-09-04; the answer was revised the same day, and §11.3a says
what moved.** Numbers, every assumption and the reproduction in
`docs/open_items/run_2026-09-03/phase_C_numbers.md` §C2; the arithmetic is
`validation/o5_a2_reach.py` (self-checking), pinned in
`python/tests/test_o5_reach.py` and `tests/test_coherent.cpp` **T10a** /
**T10c**. Every number it multiplies is read from code:
`a2_from_quadrupole` at the *measured* `LI6_QUADRUPOLE_FM2`; the
`estarlight_li6_coherent()` and `estarlight_li6_q2_floors()` tables in
`coherent.hpp` — the single code home of §11.1's cross sections and of the
2026-09-04 photoproduction scan; `COHERENT_JPSI_EFF_IR8_LI7`; `Scenario` /
`tensor_flip_plan`; and `tagging_optics` / `yr_optics`.

**THE DECIDING NUMBER.** Coherent J/ψ off tensor-polarised ⁶Li, at
`Scenario::lumi_fb_per_nucleon` = 10 fb⁻¹/u ("one EIC year"; band {1, 10, 100};
there is **no Li luminosity in any source this repository has seen**,
`needs_survey.md` §3.8), hence 10/6 = 1.667 fb⁻¹ of e+⁶Li; **σ_coh = 11.971 nb
over the whole Q² range** (`estarlight_li6.md` §2f — 6.75× the 1.773 nb of the
0.1 < Q² < 100 window, which was an acceptance study's kinematic range, not a
physics one); **J/ψ → e⁺e⁻ *and* μ⁺μ⁻ (0.11932)**;
`COHERENT_JPSI_EFF_IR8_LI7` = 0.1775 (a ⁷Li number, used as a stand-in);
`tensor_flip_plan(0.6)`, i.e. P_zz = {+0.6, −1.2} at equal shares, on the
**background-immune** ⟨P_zz²⟩ = 0.81 rather than the optimal 0.90 —

> **N = 4.22 × 10⁵ reconstructed events → δa₂(|t| = 0.3 GeV²) = 0.0100 against
> a predicted a₂ = +0.0263. S = 2.62 σ. Three sigma needs 13 fb⁻¹/u, INSIDE
> the {1, 10, 100} fb⁻¹/u band and close to one EIC year.** *(That row is the
> UNCORRECTED point and may not be quoted alone — see the band below.)* With a *perfect*
> detector it is 6.21 σ. On the measured-radius density, 1.58 σ and
> 36 fb⁻¹/u — still inside the band.

**AND THAT ROW IS THE MIDDLE OF A BAND, NOT THE ANSWER — §11.3b.** ε_det above
was applied at the wrong beam energy and carried no decay-lepton factor at all.
Correcting both:

> **S = 2.63 σ at the band's LOW EDGE and 2.84 … 3.29 σ at its TOP, 3 σ at
> 8.3 … 13.0 fb⁻¹/u — inside the {1, 10, 100} fb⁻¹/u band at BOTH ends, and
> OPEN below**, because the per-lepton reconstruction efficiency is unbounded
> in this tree. The 2.62 σ point lands 0.3 % under the band's low edge, because
> the two omissions cancel to 0.7 %. **The TOP is a span, not an edge, and
> whether it crosses 3 σ is NOT established** — §11.3c.

### 11.3a What the first pass got wrong, and by how much

The first pass of this section said **"NO … 0.75 σ … 3 σ needs 160 fb⁻¹/u …
the measurement does not exist at any luminosity the EIC is quoted at."** Two
restrictions, both artefacts rather than physics, produced that:

1. **One lepton channel.** The row used J/ψ → e⁺e⁻ only. The section's own
   prose said, three lines below its headline, that "adding μ⁺μ⁻ gains √2" —
   and 160/2 = **80 fb⁻¹/u is inside the band**, so the write-up refuted its
   own conclusion before any new run was made. `branching_all` = 0.11932
   (eSTARlight's own `JpsiBree` + `JpsiBrmumu`) is now in `coherent.hpp`.
2. **One Q² window.** σ = 1.773 nb is 0.1 < Q² < 100 GeV², which
   `estarlight_li6.md` §1 states is *"arXiv:2511.05638's range verbatim"* — an
   acceptance study's kinematic range. Q² < 0.1 was absent from the whole
   chain, and it is most of the rate. A fresh
   eSTARlight run of the same build, beams and seed (§2f; the Q² > 0.1 rows
   reproduce §2a/§2b to every printed digit) gives **11.971 nb with no floor
   at all — 6.75× — with the |t| slope unchanged to 0.4 %.** 85 % of the
   coherent J/ψ rate is below Q² = 0.1.

Together: **×13.5 in rate, 0.749 σ → 2.62 σ, 160 fb⁻¹/u → 13.1 fb⁻¹/u** (that
point is the band's uncorrected middle, not the answer; §11.3b–c). The
ladder, one correction at a time, is `phase_C_numbers.md` §C2.4a. Two smaller
things were wrong beside them: the self-check pinned the verdict on the
restricted row (`assert lumi_for_3sigma_jpsi > 100.0`, which could never
fail); and the headline ⟨P_zz²⟩ = 0.90 was the *optimal* two-fill combination
while the separation argument below requires the background-immune difference,
0.81. Both are fixed, and the pin now runs in both directions.

### 11.3b What the SECOND pass got wrong: the efficiency chain, in both directions

*(2026-09-04, third pass. `phase_C_numbers.md` §C2.3c and §C2.8 items 3a/3b/3c;
`validation/o5_a2_reach.py` §3b; `coherent.hpp`'s two new tables; pinned in
`tests/test_coherent.cpp` **T10d** and in `python/tests/test_o5_reach.py`.)*

The MARGINAL verdict survives. The chain it rested on had **two defects of
comparable size and opposite sign, and §11.3's own "every assumption" list —
`phase_C_numbers.md` §C2.8 — carried neither.**

**(1) ε_det was applied at the wrong beam energy, in the one direction the
write-up never hedged.** `COHERENT_JPSI_EFF_IR8_LI7` = 0.1775 is a ⁷Li number
at **18 × 117.9 GeV/u — that nucleus's own top energy** — and this tree's
configuration-identical eSTARlight row for those beams has ⟨W⟩ = 43.2 GeV. The
O5 rate sample is ⁶Li at **10 × 99.5**, ⟨W⟩ = 30.2. arXiv:2511.05638 measures
exactly that dependence in its §V.B, on ³He: **32.23 % at 18 × 183, 54.38 % at
10 × 100, 99.77 % at 5 × 41** (verified against the PDF; in code as
`chang26_he3_energy_scan()`). d ln ε/d ln E_ion = **−0.866** and **−0.681**.
**Direction UP; size ×1.122–1.158 for the 117.9 → 99.5 step**, i.e.
ε_det = 0.199–0.206 rather than 0.1775. ~~The ⁷Li → ⁶Li substitution beside it
is worth a further ×1.16–1.22 at fixed rigidity — the paper's own species list
is entirely at `A/Z × E` = 275 GeV/e (`chang26_species_efficiency()`), so its
A-ordering is a species lever on its own.~~ **RETRACTED 2026-09-04, fourth pass
— that is a non sequitur and the ×1.16–1.22 was an unquantified factor
presented as quantified. §11.3c.**

**It is not the paper's own 1.687.** 54.38/32.23 is their 183 → 100 ratio, a
step **3.56× larger in ln E** than this transfer needs; carried whole it gives
3.40 σ and 3 σ at 7.8 fb⁻¹/u, by using the wrong lever arm. *(The tree already
invoked the same Fig. 2 W-dependence for the far smaller 32.2 → 30.2 GeV
Q²-floor shift and called it "conservative on the W axis", while never invoking
it for the 43.2 → 30.2 GeV beam-energy shift; `coherent.hpp` recorded "none at
10 × 99.5" with neither the direction nor the size. Both now carry both.)*

**(2) There was no decay-lepton acceptance anywhere in the chain.** The chain
is σ × BR(ℓ⁺ℓ⁻) × ε_recoil, and ε_recoil is the fraction of scattered
**nuclei** in the far-forward acceptance — *"32.23 % of the scattered ³He
nuclei occur within a safe distance from the beam"* — under a simulation whose
own p. 4 says it *"only accounts for the acceptance effect and does not
incorporate the efficiencies of the detector. Additionally, we did not account
for the efficiency and acceptance of the reconstructed distribution"*. So there
was **no central-detector acceptance and no reconstruction efficiency for the
e⁺e⁻/μ⁺μ⁻ pair at all**. Direction **DOWN**.

**Bounded, in part.** At 10 × 99.5 a J/ψ at ⟨W⟩ = 30.2 comes from an
E_γ = 2.29 GeV lab photon and sits at y = −0.39 — nearly at rest — so its two
1.548 GeV leptons are central. Using this repository's **own**
`Scenario::eta_max` = 3.5 (documented there for the *scattered electron*,
because there is no decay-lepton acceptance in this tree), and SCHC transverse
decay as the pessimistic weighting: **A_geom = 0.994 at ⟨W⟩, never below 0.89
anywhere the beams can reach**. The geometry is not where the factor is. **The
per-lepton reconstruction efficiency is UNBOUNDED HERE** — nothing in this
repository supplies one. The 0.95/track used to draw the band's low end is the
sibling `../PolarizedLithiumSim`'s `HfsModel(eff_track=0.95)`, which that file
itself labels a stand-in.

**(3) They partly cancel, and that is the finding.** ×1.122 against ×0.897 is
**1.007**. A reader told only about the beam energy moves the verdict up; a
reader told only about the leptons moves it down; both together move it by
0.7 %. That is precisely why leaving both out of a list that claimed
completeness was a defect rather than a rounding.

**(4) A self-contradiction, printed at two of three sites.** *"Every single
conservatism on its own leaves 3 σ inside the band"* is false: the de-squeezed
optics row is a counterexample, and `validation/o5_a2_reach.py` printed that
sentence **three lines below** its own bullet saying the optics put 3 σ "back
outside the band". `mantysaari_collaboration_draft.md` repeated it uncorrected
**and omitted the optics row from its ladder entirely**. §C2.4a stated the
exception inside the same sentence as the claim. All three now carry the
exception as an exception, and the draft's ladder has the optics row.

**The verdict as the third pass left it, and the pin it wrote.** S =
2.63 … 3.15 σ, 3 σ at 9.1 … 13.0 fb⁻¹/u, inside {1, 10, 100} at both ends and
**open below**. `Optics::lumi_fraction` = 0.0781 remains **the single
correction that on its own restores a NO**, and it does so at both ends of the
band. The self-check no longer pins a point (`2.0 < S < 3.0`); it pins the
band, both directions, the cancellation, the fact that the beam factor is
*not* 1.687, and the optics exception — and nine mutations of that claim were
checked to make it fail. **The 3.15 σ top did not survive §11.3c: read the
corrected band there, 2.63 σ at the low edge and 2.84 … 3.29 σ at the top,
3 σ at 8.3 … 13.0 fb⁻¹/u, optics 0.73–0.92 σ at 106–167 fb⁻¹/u.**

**The band is asymmetric in what it protects.** *"3 σ inside the {1, 10, 100}
band"* survives a pair acceptance × efficiency all the way down to **0.117** —
a detector reconstructing about one J/ψ in eight. *"MARGINAL rather than NO"*
only survives down to **0.52**. Quote the two with that asymmetry attached.

**Why: the |t| slope, not the modelling.** a₂ is *linear* in |t| and the
coherent sample is e^{−B|t|} with B = 39–55, so the information-weighted
modulation is κ√⟨t²⟩ = **0.25 %** at B = 50, not the 2.6 % that a₂(0.3)
suggests — a factor **10.6**, and it follows from the coherent form factor
alone. The |t| window is *not* where anything is lost: |t| < 0.2 GeV²
(`COHERENT_T_MAX_DEFAULT`) keeps 99.995 % of the rate and 99.72 % of the
Fisher information.

**And it is §11.2's own factor 7.5 that sets the scale.** The same J/ψ sample
gives **19.7 σ** at this α+d geometry's own Q = −0.615 fm², **6.4 σ** at
GFMC's −0.20, and **2.6 σ** at the measured −0.0818. ⁶Li is attractive
precisely because Q(⁶Li) is a *near-null* against the deuteron's +0.286 — and
this is the price of that: **near-null is expensive.** What it is not, on
these numbers, is unmeasurable.

**The separation is free, and is not the problem.** Two independent handles.
(a) The tensor modulation is about the **spin axis**, the photon-polarisation
one about the **photon's linear-polarisation direction** (the lepton plane at
Q² > 0), which is uniform in the lab relative to a fixed vertical φ_S — so it
averages to zero in a spin-referenced histogram under a φ_γ-uniform
acceptance. (b) The tensor term is **odd** in P_zz (exactly: a₂(0) = −2a₂(±1)
makes Σ_m p_m a₂(m) = P_zz a₂(±1) identically) and the photon term is
**even**, so the two-fill difference is background-free at first order, and it
cancels a relative-luminosity offset because each fill's amplitude is a
self-normalised ratio. **Cost: negative.** Against the optimal use of the same
two fills the background-immune difference loses 5.4 % in δ; against putting
the whole luminosity into one +0.6 fill it **gains 1.50×**, because the
m = 0-rich fill carries |P_zz| = 1.2. The residual — the two fills' φ-averaged
|t| spectra differ only at second order in ε = |`eps_b0_equivalent()`| =
0.0067, i.e. by 1.13e−3 at |t| = 0.2 and 2.26e−5 averaged — matches √(2/N)
only for a background amplitude A_γ above **1.9 (J/ψ), 0.13 (φ), 0.017 (ρ⁰)**
on the whole-Q² samples (7.1 / 0.60 / 0.10 on the smaller 0.1 < Q² < 100
ones). **No magnitude for A_γ exists anywhere in this tree**
(`physics_literature.md` records the mechanism, STAR arXiv:2204.01625, not a
number), which is why the argument rests on parity and not on a size; for ρ⁰
it means binning in |t| is mandatory.

**BOTH HANDLES SURVIVE Q² → 0**, which is what lets the photoproduction region
be used at all. Below Q² = 0.1 the scattered electron is not detected and φ_γ
is unknown event by event; neither handle needs it — (a) is a statement about
a *spin-axis* histogram under a φ_γ-uniform acceptance, and (b) never
references the lepton plane. *(The premise that the photon effect "dies as Q²
rises" is not supported here and is not used: the virtual photon's
linear-polarisation degree is a function of **y**, not a decreasing function
of Q². What Q² > 0.7 buys is that the scattered electron is detected, so φ_γ
is known event by event — measurability, not suppression, and it is a control
this argument never needed.)*

**The light mesons have more statistics, and it does not change the ask.**
With the same assumptions and no Q² floor: **φ 39 σ** (8.3 σ over
0.1 < Q² < 100, 2.1 σ inside LiPolGen's own Q² > 0.7 window) and
**ρ⁰ 293 σ** (49 σ / 8.6 σ). But both sit **below**
`COHERENT_MX_MIN_DEFAULT` = 1.2 GeV, so the shipped coherent channel carries
no rate there at all; ρ⁰ is the channel in which the photon-polarisation
cos 2φ has actually been *observed* (STAR's ρ⁰ ultraperipheral measurement —
this tree records the reference, not its channel and not its magnitude), and
where the P_zz-flip residual first bites, now at A_γ ≈ 0.02; and it is where
the closed-form map is **least** defensible — a large, strongly absorbed
dipole at ⟨W⟩ = 14–18 GeV, where black-disc absorption saturates the very
anisotropy the map is linear in, validated only against a **J/ψ** calculation.
J/ψ being measurable is what makes that irrelevant to the ask.

**VERDICT: MARGINAL — as a BAND — and the ask stands.** Concretely:

1. **The ⁶Li tensor a₂ in coherent J/ψ is a 2.6 σ (band low edge) to
   2.8–3.3 σ (band top) measurement at one EIC year, with 3 σ at
   8.3–13.0 fb⁻¹/u** (§11.3b–c; the 2.62 σ / 13.1 fb⁻¹/u point is the
   uncorrected middle of that band, not the answer, and whether the top
   crosses 3 σ is **not established**). Propose the collaboration on that channel —
   which is also the only channel the published calculation covers, so the
   theory ask stays the easy one. *(This replaces "do not propose the
   collaboration on the ⁶Li tensor a₂ in coherent J/ψ; the measurement does
   not exist at any luminosity the EIC is quoted at", which was written from
   the restricted row of §11.3a and is withdrawn.)*
2. **Propose it as a photoproduction measurement.** Q² < 0.1 is 85 % of the
   rate, and it is also where [Mant24]'s own Fig. 4 lives — the digitised
   deuteron a₂ this whole map is validated against is a *photoproduction*
   calculation, so the request now matches the published kinematics better
   than the first draft did.
3. **Say what is missing, because it is not statistics.** No detection
   efficiency exists below Q² = 0.1 in any source this tree has seen; **no
   decay-lepton reconstruction efficiency exists in this tree at all** (§11.3b
   — the geometric half is bounded at 0.99 and is not the problem); and the
   far-forward working point is unchosen — at LiPolGen's own de-squeezed ⁶Li
   tagging optics (`lumi_fraction` = 0.0781 at 10 × 100) the whole band becomes
   0.73–0.92 σ with 3 σ at 106–167 fb⁻¹/u. **That last one is the single
   correction that on its own restores the NO**; the others widen the band.
   Say also that the ⁷Li → ⁶Li efficiency substitution the whole chain rests on
   is **undetermined in direction** (×0.99–1.33, §11.3c) and is what makes the
   band's top a span. All four are the collaboration's and the detector groups'
   to supply. Say also that the ε_det used is a **top-energy ⁷Li** number
   applied at 10 × 99.5, i.e. conservative by ×1.12–1.16 on the beam-energy
   axis alone — the recipient will find that themselves otherwise.
4. The cheap thing still worth doing is the **sign**: ⁶Li's a₂ is predicted
   **positive** where the deuteron's is negative, and a sign is one bit. A
   2.6 σ sample does deliver that bit.

**THE LIMITATION RIDES WITH THE NUMBER.** All of the above is built on
`a2_from_quadrupole`, a **closed form and not a dipole-model amplitude**: it
carries the target's quadrupole through [Mant24]'s published |t| dependence
(reproducing their a₂(m = ±1) to 8 % with zero free parameters) and contains
no Good–Walker average, no amplitude, no saturation model, none of their
uncertainties, and it uses the *matter* quadrupole as a proxy for the
transverse **gluon** anisotropy. A dipole-model run could move the signal by a
factor; it would now have to move it **up by 1.15** to bring coherent J/ψ to
3 σ at 10 fb⁻¹/u, **down by 1.31** to fall back to 2 σ, or up by 1.91 to reach
5 σ — equivalently, the luminosity would have to change by 1.31 in the first
case, which is the 13.1 fb⁻¹/u of the band's uncorrected middle. The first pass of this paragraph quoted
4.0 / 6.7 / 16; those belonged to the restricted row of §11.3a.

**The rule of this section is unchanged and is now quantified rather than
repealed: do not derive a tensor input for a published observable from these
wave functions.** And since 2026-09-04 there is a second rule beside it:
**do not quote a₂(|t| = 0.3) = +0.026 as "the ⁶Li tensor signal" without the
sample it would be measured on.** The coherent sample lives at |t| ≈ 1/B, so
the observable is κ√⟨t²⟩ = **0.25 %**, a factor 10.6 smaller. That factor is
real; what it is not is a verdict — the modulation on the sample that would
actually be taken is ā₂ = κ√⟨t²⟩ = **0.32 %** at that sample's own B = 38.8,
and 0.0032 × √(2 × 0.81 × 4.22 × 10⁵) = 2.6 σ, so
O5 comes out MARGINAL (§11.3, §11.3a). The α+d truncation reproduces the ⁶Li point radius to 3 %
(4 % against the VMC `li6.density` value) but
overshoots Q(⁶Li) by **7.5×** — **3.3165 × 2.2686, two factors and not three**,
the first of them anchored on the *measured* asymptotic D/S ratio and carrying
that measurement's error bar, which spans model Q from −0.4005 fm² to
**+0.0298 fm², through zero** (model −0.615…−0.730 fm² against the measured
−0.0818 and GFMC AV18+IL7's −0.20(6), Pastore *et al.*, PRC 87, 035503
(2013)). `quadrupole_band_fm2()` returns all three and the writer stamps them;
the `quadrupole_target_fm2` dial can put the geometry on the measured Q, but
it is a **deformation dial, not a wave function**, and the sidecar labels it
as one.

### 11.3c What the THIRD pass got wrong: the species leg, and the band's TOP

*(2026-09-04, **fourth pass**. `phase_C_numbers.md` §C2 fourth-pass changelog;
`validation/o5_a2_reach.py` §3b(b) and `species_scaling_same_energy`;
`coherent.hpp`'s `chang26_species_efficiency()` comment; pinned in
`tests/test_coherent.cpp` T10d and in `python/tests/test_o5_reach.py`.)*

**The MARGINAL verdict, its band's LOW edge and "inside {1, 10, 100} at both
ends" all survive. The band's TOP does not.**

**(1) A non sequitur, at every site that carried it.** *"Every entry is at the
same rigidity, so its A-ordering is a species lever on its own"* —
`o5_a2_reach.py` (twice), `coherent.hpp` (twice), `python/bindings.cpp`,
`phase_C_numbers.md` §C2.3c(b) and §C2.8 item 3a, `STATUS.md`, this section's
§11.3b, and the name of a `tests/test_coherent.cpp` TEST_CASE. **RETRACTED.** Eliminating rigidity
does not leave species alone. At fixed R = A E/Z the per-nucleon energy is
E/u = R Z/A and the total beam momentum is p_z = Z R, and **both still vary
down the list**: E/u from 118 GeV/u (⁷Li) to 183 (³He), Z from 1 to 8, p_z from
274 to 2192 GeV. The chain's *own other leg* puts d ln ε/d ln E at −0.68 … −0.87,
so the list's **×1.55 spread in E/u is worth ×1.35–1.46 in ε** — the size of
the entire claimed species gain. The list is a **joint (A, Z, E/u) lever**; the
chain asserted the decomposition and never tested it.

**(2) The one validation was mis-specified, and the direction it was read as
giving is withdrawn.** The form is ε = exp(−B(A)·p_T,cut²) with p_T,cut
inverted from ⁷Li (0.1958 GeV). The criterion arXiv:2511.05638 states is
*"within a safe distance **from the beam**"* — a cut on the **angle** — so
p_T,cut = θ·p_z = θ·Z·R and **it carries the charge**. ³He and ⁴He are Z = 2;
⁷Li is Z = 3. The shipped test applied ⁷Li's Z = 3 cut to the Z = 2 pair, and
its own comment said *"same Z"* without noticing it was the **wrong** Z. Three
readings of the one table, and **nothing in the table picks between them**:

| reading | on ³He/⁴He | on the list's absolute values |
|---|---|---|
| p_T,cut ∝ Z — the geometric criterion taken literally | pred **1.0967** vs measured **1.0955**: **0.1 %** | badly wrong — pred/meas 1.95 (²D), 2.00 (⁴He), 0.21 (⁹Be) |
| p_T,cut Z-independent | pred 1.2309 vs 1.0955 — this is the shipped one, "overstates by 12 %" | **fits the light end** — pred/meas **1.003** (²D), **1.034** (⁴He), **1.047** (⁹Be); ³He off by **16 %**, which is what its E/u = 183 against the family's 137 is worth on the chain's own scan (an **18–22 %** shortfall). It is a *local* form, not a description of the list: ¹²C off by 32 %, ¹⁶O by ×3.1 |

The shipped reading is the only one of the three that yields the directional
instruction **"Read the low end"**, printed at three sites. **That instruction
is deleted**, and `tests/test_coherent.cpp` now carries the Z-corrected
prediction beside the wrong-Z one.

**(3) The four entries the chain never used.** ⁶Li's own fixed-rigidity energy
is Z/A × 275 = **137.5 GeV/u**, and **four** of the seven entries — ²D, ⁴He,
¹²C, ¹⁶O, at A = 2, 4, 12, 16 — sit at **137 GeV/u**. They bracket A = 6 with
**no energy step at all**, and ⁴He → ¹²C is the adjacent pair. The same closed
forms on that same-energy set give

> ε(⁶Li, 137 GeV/u) = **0.1763 … 0.2366** against ⁷Li's 0.1775, i.e.
> **×0.993 … ×1.333 — it STRADDLES 1.** log R_G gives ×1.004, the p_T
> threshold anchored at ¹²C gives ×0.993, log A ×1.130, and linear-in-A
> ×1.333.

against the ×1.160 … ×1.219 the confounded reading gave. **The direction of the
⁷Li → ⁶Li substitution is not established, let alone its size.**

**The corrected verdict.** The band's low edge is untouched, and so is
everything that rests on it:

> **S = 2.63 σ at the band's LOW EDGE and 2.84 … 3.29 σ at its TOP, with
> 3 σ at 8.3 … 13.0 fb⁻¹/u — INSIDE the {1, 10, 100} fb⁻¹/u band at both ends
> on every form, and OPEN below.** MARGINAL stands.
>
> **"The band's TOP just crosses 3 σ" is NOT ESTABLISHED and is withdrawn.**
> The top straddles 3 σ: two of the five forms put it at 2.84 and 2.86, two at
> 3.03 and 3.04, one at 3.29.

`_self_check` and `test_o5_reach.py` no longer pin a scalar top. There is no
`significance_hi`, no `lumi_for_3sigma_lo` and no `band_factor_hi` key left in
the returned dict — a `KeyError` is the cheapest way to stop the point coming
back — and the pin is

```python
assert res["significance_hi_min"] < 3.0 < res["significance_hi_max"]
```

which **fails if "the top crosses 3 σ" is ever re-asserted as established**,
because that means the top's low reading has been pushed above 3.

The de-squeezed-optics rung moves with the top: the band there is
**0.73–0.92 σ with 3 σ at 106–167 fb⁻¹/u**, still outside {1, 10, 100} at both
ends, so `Optics::lumi_fraction` remains the single correction that on its own
restores a NO.

### 11.4 O4 and the VMC inputs — ΔB defined, and four opt-in knobs (2026-09-04)

Full record and every number: `docs/open_items/run_2026-09-03/phase_C_numbers.md`
§C4 and §C5, pinned in `tests/test_coherent.cpp` **T10b**,
`tests/test_cluster_config.cpp` **T23**, `tests/test_tagged.cpp`
**T24**/**T25**/**T26** and mirrored in `python/tests/test_module.py`. **No
shipped default moved and `validation/reference/*.json` is untouched.**

* **O4 — ΔB is defined, once**, at `CoherentScenario::eps_b0`:
  |F_m|² = exp(−|t|[B + ΔB_m cos 2(Φ−Φ_S)]) with ΔB_m = δ_m/2, hence
  a₂(m) = −(ΔB_m/2)|t| and **eps_b0 = δ_{±1}/B = 2ΔB_{±1}/B = −ΔB₀/B**. The old
  label was off by a **sign**. `delta_b_m()` / `slope_at_azimuth()` are its code
  home; `quadrupole_from_a2_slope()` inverts the map so a scenario can be
  *asked* what it assumes. **eps_b0 = −0.08 assumes Q_charge(⁶Li) = −0.9345 fm²,
  11.42× the measured value**, and a₂(±1, 0.3) = +0.300 — bigger than the
  *deuteron's* own −0.28. Honest ⁶Li band **−(0.0070 … 0.0527)** at B = 50.
  Author decision: default kept, cost recorded, and `eps_b0`/`slope_b` are
  flagged as **not independent** (band the product δ = eps_b0·B, never eps_b0
  alone).
* **C5.1 — no η dial.** η is exactly linear in `quadrupole_dial_s()`, so an
  `eta_target` option would be a second name for one knob. Added instead:
  `ClusterConfigSampler::quadrupole_for_eta()`, a converter that re-derives
  §11.2's whole GK band (−0.4005 / −0.1842 / **+0.0298**, through zero) from the
  tables instead of from typed numbers.
* **C5.2 — the MC errors are carried, and they are negligible.** The ANL 1σ
  columns were parsed and dropped; `VmcRadial` now carries `dpsi()`,
  `shifted_by_sigma()` and `norm2_error()` (correlated **and** quadrature, since
  one variational walk fixes neither), with `vmc_mc_sigma` /
  `PipelineConfig::cluster_vmc_mc_sigma` as the band knob (refused where no
  errors exist). Measured: **0.47 %** / 1.62 % on the S/D norms, **1.1 %** on
  P_D per σ, and **0.02 %** on the tagged tensor observable — against 6.6 % for
  the wave-function choice and a factor 7.5 for the quadrupole.
* **C5.3 — the 5 % is propagated.** N_αd spans 5.311 %, the D-wave *norm*
  7.181 %, P_D only **2.729 %**; outside b₁ (exactly linear, exactly ±5 %) it is
  worth **0.048 %** on `tensor_dilution`.
* **C5.4 — the deuteron control reaches AV18, opt-in.** `fdeut.av18`'s u(k),
  w(k) *are* the p–n relative waves, so `--cluster-wave vmc` was being silently
  ignored on that channel; it now selects P_D = **0.0576** (+28 % over the
  scenario 0.045), worth −2.0 % / −1.2 % on the channel's vector/tensor
  dilutions. The relative S–D sign does **not** flip (the **stored** ψ₂ = +W),
  unlike ⁶Li — but the amplitude sums φ₂ = i²ψ₂ = **−W**, which `build_amp2`
  did not apply until 2026-09-06 (`07_cw_sign_investigation.md`).
* **C5.5 — the split is real, the framing was not, and substitution is worse.**
  It is **not** two channels of one run (a `Pipeline` builds exactly one
  channel); it is a tagged row and an inclusive row of one programme, which is
  worse because nothing can notice it. Size: **+11.61 %** vector, **+6.58 %**
  rank-2. But the "consistent" inclusive value 0.905427 sits **6.8 % above** the
  ab-initio six-body VMC **0.848** where the shipped 0.811228 is 4.3 % below it,
  so the cluster product — not the choice of P_D — carries the error. Author
  decision: default kept, drift documented, and **no inclusive ⁶Li polarization
  is quoted without the band 0.81 … 0.91**.
* **C5.5b — the claim C5.5 shipped was false, and C5.4 is what made it false.**
  "Nothing *inside* one run is inconsistent" went to five sites and was wrong:
  `li6_alpha_channel` set `dis_target = DEUTERON()` **unconditionally**, and
  the T1 `ClusterBreakup` built its deuteron at `P_D_DEUTERON` = 0.045, so a
  `--cluster-wave vmc` α-tag run took the α–d **relative motion** from the ANL
  VMC AV18+UX overlap and the **embedded deuteron** from the 0.045 scenario —
  two deuteron wave-function families in one run, with every polarized
  tagged-α observable **+2.069 %** high (0.9325 against the AV18 deuteron's
  0.9136001620; exact, because g₁A is linear in the effective polarization).
  C5.4's own rationale forbids it: *switching only the control would make the
  control and the channel it controls two different wave-function families*.
  **Fixed:** `DEUTERON_AV18()` and `BreakupOptions::source` follow the flag, so
  it means one deuteron everywhere it is read; the opt-in path moved −2.027 %,
  the Hulthén default is bit for bit, and a `vmc` run's whole-nucleus reading
  is **0.887076** — one Hamiltonian end to end, and C5.5's own third table row.
  **And the secondary defect:** `validate()` accepted `cluster_wave` on
  `Inclusive` and `CoherentLi6`, where it is never read — the mechanism by
  which C5.5's 11.61 % reached a user who thought they had asked for a VMC
  ⁶Li. It is now refused, as `cluster_vmc_mc_sigma` (a 0.02 % effect) already
  was. `phase_C_numbers.md` §C5.5b, **T27**.

> **§11.4 checked against the 2026-09-06 S–D interference fix — no number in
> the five C5 bullets above inherits it.** The tagged sector's partial-wave
> amplitude was summing ψ_L where it needs φ_L = i^L ψ_L
> (`src/core/tagged.cpp` `build_amp2`, one missing `(-1)^floor(L/2)`;
> `docs/benchmarking/07_cw_sign_investigation.md`), which flips the S–D
> interference and therefore every *angle-differential* tagged quantity. Every
> quantity these bullets quote is **angle-integrated**, and `∫Θ₀Θ₂ dc = 0` by
> L-orthogonality, so the flipped cross term integrates away and only the
> 96-cell midpoint-quadrature residual of that zero survives. Re-measured:
> **C5.2**'s tagged tensor observable `0.982575817 → 0.982575872`
> (+5.6e−8 relative, against the 0.02 % the bullet reports); **C5.3**'s
> `tensor_dilution` sensitivities likewise, being differences of such values;
> **C5.4**'s `−2.0 % / −1.2 %` on the control's vector/tensor dilutions
> **unchanged to every printed digit** (−2.02683 % / −1.18214 %, from
> −2.02682 % / −1.18214 %); **C5.5**'s `0.905427`, `0.811228`, `0.848` and
> `0.887076` are closed forms in D-state probabilities (`beams.hpp`
> `li6_cluster_polarization`, `constexpr`) and never reach `build_amp2` at
> all; **C5.5b**'s **+2.069 %** re-measures as **1.020687624664** against
> **1.020687500612** — unmoved in the seventh figure. What *did* move is the
> spectator-differential observable these bullets do not quote: ⁶Li
> A_zz^tag(k = 0.20 GeV) went **+0.845 → −1.207** (Hulthén) and
> **+0.452 → −0.519** (VMC); see §1 and
> `docs/open_items/run_2026-09-06/phase_CW_numbers.md`.

### 11.5 O3 — what is bounded offline, and C6 — the collaboration ask reconciled (2026-09-04)

**O3, partially bounded.** The question — how much would the
incoherent/coherent split move if the α core were sampled from genuine
correlated GFMC ⁴He configurations instead of an uncorrelated product of
one-body densities — stays **unanswerable from this repository** in the form
it was asked: no configuration table for any nucleus is in this tree, and the
split is a property of a Good–Walker amplitude's configuration-to-
configuration fluctuation, which nothing here computes. What *can* be
bounded with data already committed: `data/vmc/bonus_other_clusters/he4.dd`
is a genuinely correlated VMC overlap of the same ⁴He wavefunction onto an
α → d+d channel (Forest *et al.*, PRC 54, 646 (1996), the same method as
`li6.ad`/`li7.at`). `validation/o3_alpha_correlation_bound.py` compares its
own relative-motion second moment to the analytic prediction of an
**uncorrelated** product of `he4.density` for the identical observable — the
separation between the centroids of the α's two 2-nucleon halves, ⟨D²⟩ =
(4/3)R₁² exactly, independent of recentring and of which 2-2 partition is
chosen (worked in the module's own docstring) — and measures the correlated
value **18.2 % above** that null in variance (**8.7 %** in rms; r_rms 1.809 fm
correlated vs. 1.664 fm uncorrelated). **Reading**: genuine 4-body correlation
is a several-to-twenty-percent effect on this kind of position-space moment
in this nucleus, not a factor of several — small next to the factor-7.5
quadrupole gap (§11.2) — but this is a different observable standing in for
the one O3 actually asks about, and it does not close O3. What genuinely
needs the configurations: the incoherent/coherent split itself.

**C6 — the collaboration ask, reconciled.** Two committed copies of the
Mäntysaari-group ask disagreed (this section's former paragraph and
`design_G_cluster_config.md` §10's, by one number: "3 %" vs "4 %" on the
α+d model's point-radius agreement — both correct, against different
reference radii). §11.2 above now points to the single reconciled draft,
`docs/open_items/run_2026-09-03/mantysaari_collaboration_draft.md`, which is
informed by §11.3's O5 answer: the letter's original request (c) asked the
group to run their polarised J/ψ machinery on this repository's ⁶Li
configuration tables, and §11.3 shows that channel is **marginally
measurable** — as a **BAND**, never as a point: **2.63 σ at the band's low
edge and 2.84–3.29 σ at its top, with 3 σ at 8.3–13.0 fb⁻¹/u**, inside the
{1, 10, 100} band at both ends and open below, and with the crossing of 3 σ at
the top **not established** (§11.3b–c; the 2.62 σ / 13.1 fb⁻¹/u point is the
band's uncorrected middle and may not be quoted alone). The
reconciled draft keeps requests (a) and (b) (the GFMC configurations and the
code generalization — both channel-independent, and (a) is the only route to
closing O3 above) and keeps (c) as a request to perform the calculation, moved
to **photoproduction** (Q² < 0.1, which is 85 % of the rate and where
[Mant24]'s own Fig. 4 lives), stating the deciding band and the four
unestablished factors plainly. *(An intermediate version of this paragraph
said the channel was "blind at any EIC luminosity this repository has ever
quoted" and rewrote (c) into a question about whether to bother; §11.3a is the
record of why that was withdrawn.)* **It is not sent. Sending it, in whatever form, is the author's
decision**, and the draft says so at its own top.

### 11.6 D5 — the |t| ceiling: 0.2 STAYS, and the reason it is written down changed (2026-09-04)

`COHERENT_T_MAX_DEFAULT` = 0.2 GeV² carried **two** stated reasons, presented
as agreeing: the Mäntysaari deformation input is digitised only to |t| ≤ 0.30,
*and* the linear c₂ crosses −1 at |t| = 0.245 for P_zz = −2. **They agree only
at the shipped `eps_b0`.** Written in closed form — with c₂ = A|t| + C,
A = −(P_zz/2)·eps_b0·B and C = amp·P_zz, the positivity edge is
|t|_pos = (1 − sign(A)·C)/|A|, i.e. 2(1/|P_zz| − amp)/(|eps_b0|·B) for the
shipped signs — the second reason moves with `eps_b0` and the first does not:

| eps_b0 | \|t\|_pos, P_zz = −2 | \|t\|_pos, P_zz = +1 |
|---|---|---|
| −0.08 (shipped, an EXACT input) | **0.2450** | **0.4950** |
| −0.0527 (α+d model Q, rounded) | 0.3719 | 0.7514 |
| −0.0171 (GFMC Q, rounded) | 1.1462 | 2.3158 |
| −0.0070 (**measured** Q, rounded) | **2.80** | 5.6571 |

**Read the digit count with the table.** The last three `eps_b0` are ROUNDED
and each row is arithmetic on the rounded value, so the bottom row is
**2.80 GeV²** — three figures, which is what a two-figure input supports —
and never "2.8000". On the DERIVED band (`ClusterConfigSampler::
quadrupole_band_fm2()` through `a2_from_quadrupole` at B = 50, amp = 0.01;
measured 2026-09-05) the three `eps_b0` are −0.0526846 / −0.0171207 /
−0.0070024 and the P_zz = −2 edges are **0.3720 / 1.1448 / 2.7990**;
`tests/test_coherent.cpp` T10b asserts the three derived `eps_b0` and the
0.37202 / 2.7990 edges.

So §11.2's conclusion is confirmed — the positivity edge at 0.245 is a
consequence of an `eps_b0` that is 11.4× the measured ⁶Li quadrupole, not a
property of ⁶Li — and **deriving the ceiling from positivity would be less
honest than the number it replaces**: at the measured quadrupole it would
license |t| up to 2.80 GeV², **9.3× outside the four digitised rows**
(|t| = 0.05, 0.10, 0.20, 0.30) the whole deformation term is scaled from.

**What was decided and implemented.** The ceiling **stays fixed at 0.2**, and
its stated primary reason is now the **anchor range**, which is a property of
the input table and does not move with any knob; the positivity edge is stated
as **secondary and contingent**, with the arithmetic on the page. Concretely:

1. `COHERENT_T_MAX_DEFAULT` unchanged — no reference gate moves.
2. Its comment rewritten in `coherent.hpp`; the same correction applied at
   `pipeline.hpp`'s `coherent_t_max`, `CoherentScenario::eps_b0`,
   `docs/CONVENTIONS.md`, `docs/PHYSICS_CHANNELS.md` and `docs/USAGE.md`.
3. `CoherentScenario::positivity_margin` and `CoherentSampler::check_positivity`
   **keep throwing, unchanged**. They are the GUARD on the truncated weight the
   sampler actually uses — an author who raises `t_max` or `eps_b0` past where
   it stops being a density is caught by nothing else — not the derivation of
   the ceiling.
4. **`CoherentScenario::t_positivity_edge(pzz)` is new**: the edge in closed
   form, so the number is *derived* rather than repeated in four docstrings,
   and so it **moves with `eps_b0`** for any caller who asks. It changes no
   default. Gated in `tests/test_coherent.cpp` (every sign combination against
   the zero of `positivity_margin`, plus the degenerate cases) and in
   `python/tests/test_meta_provenance.py`.
5. `coherent_t_max` is now **reachable as `--coherent-t-max`** and **recorded
   in the npz `meta`**, together with `coherent_slope_b`, `coherent_eps_b0`,
   `coherent_amp`, `coherent_f0`, `coherent_m_x_min`, `coherent_x_pom_max`,
   `coherent_weighted_azimuth` and the derived
   `coherent_t_positivity_edge_pzz_m2`. Before this a Python caller who moved
   `coherent_t_max` changed the entire |t| spectrum, the tag acceptance and
   every c₂ in the file, and the sidecar could not tell.

**The number that says the ceiling is not a rate question.** Measured on
**200 000 generated coherent events** at the shipped defaults (seed 99, plan
`tensor-thirds` at P_z = 0.7 / P_zz = 0.6, 4 threads — the plan splits the
events across categories, so a single-category run of the same seed gives
0.019974; max |t| and the zero count above 0.20 are the same either way,
`phase_D_small_items.md` §D5.3):
⟨|t|⟩ = 0.019963 GeV², max |t| = 0.197575, and the fractions above 0.05 /
0.10 / 0.15 / 0.20 are 0.082015 / 0.006550 / 0.000455 / 0. `sample_t`
**renormalises** on [0, t_max], so the ceiling does not lose rate — it
redistributes exp(−B·t_max) = **4.5e−5** of it. Moving the ceiling to 0.245 or
anywhere else inside the anchor range changes the generated sample at the
1e−5 level. The ceiling is a statement about **where the model is defined**,
and it is now justified as one.

**What this depends on, and what it does not.** Only the *secondary* reason
depends on the `eps_b0` author decision (STATUS.md decision table row 8): if
`eps_b0` is ever corrected to the ⁶Li band, positivity stops binding entirely
and the anchor range is the only ceiling left — the ceiling itself does not
move. The follow-on that *does* wait on that decision is un-truncating the
azimuthal modulation (sampling φ from Σ_m p_m e^{−ΔB_m|t|cos 2φ} directly
instead of from its O(ΔB|t|) truncation 1 + c₂ cos 2φ). That is the physically
right form — it is positive for every |t| by construction — but it moves every
pinned coherent azimuth, so it is filed as the follow-on to the `eps_b0`
decision rather than as an independent item.

## 12–13. Engineering — **CLOSED / OPT-IN SHIPPED** (rows 13, 14, 15, 24 open)

- **In-tree now** (2026-09-02; the prototype in `engineering.md` §B was
  redone for real, prototype files themselves are gone): `pyproject.toml`
  (`[build-system] requires = ["scikit-build-core>=0.9", "pybind11>=2.10"]` —
  CMakeLists.txt's own pybind11 discovery is `find_package(pybind11 CONFIG)`
  via `python -m pybind11 --cmakedir`, no vendored copy, so this is the only
  pybind11 declaration and it is what the isolated build env installs).
  `[project].version` is `dynamic`, read out of CMakeLists.txt's
  `project(LiPolGen VERSION 0.1.0 ...)` by
  `scikit_build_core.metadata.regex` — CMake's `project()` version is the
  one source of truth, not retyped in `pyproject.toml`.
- `wheel.packages = ["python/lipolgen"]`; `cmake.args =
  ["-DLIPOLGEN_BUILD_TESTS=OFF"]` (the option already existed, no new CMake
  option added); `build-dir = ".skbuild/{wheel_tag}"` (repo-local, gitignored,
  never the in-tree `build/`); `editable.mode = "redirect"` (`rebuild` left
  at its default/off — scikit-build-core's ninja/cmake are only guaranteed
  present in pip's *build*-isolation env, not in the venv doing the
  importing afterwards, so on-import rebuild is not reliable there).
  `LIPOLGEN_DEPS_PREFIX` honoured both via
  `--config-settings=cmake.define.LIPOLGEN_DEPS_PREFIX=…` (scikit-build-core,
  no extra wiring needed) and, as a convenience, from the environment
  (`CMakeLists.txt`: `if(NOT DEFINED CACHE{LIPOLGEN_DEPS_PREFIX} AND DEFINED
  ENV{LIPOLGEN_DEPS_PREFIX}) …`).
- `CMakeLists.txt`: `install()` is now branched on `SKBUILD`. The `SKBUILD`
  branch installs `_lipolgen` and all four `lipolgen_*` bridge libraries into
  **one** destination, `lipolgen/` (the package dir `wheel.packages` also
  populates) — the prototype's `lib/` + `lipolgen/` duplication is gone, one
  copy of each `.so`, confirmed via `python -m zipfile -l` on a built wheel.
  The non-`SKBUILD` branch keeps the previous `include/` +
  `${CMAKE_INSTALL_LIBDIR}` install for conventional `cmake --install
  --prefix …` use, untouched. `CMAKE_INSTALL_RPATH = "$ORIGIN:
  ${LIPOLGEN_DEPS_PREFIX}/lib"`, `CMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE` and
  `CMAKE_BUILD_WITH_INSTALL_RPATH=TRUE` all live inside `if(SKBUILD)` —
  **not** unconditional. That guard is load-bearing, not cosmetic: CMake
  pads a target's BUILD-tree RUNPATH with empty (`:`) entries at configure
  time for any target with an `install()` rule once `CMAKE_INSTALL_RPATH` is
  set at all, regardless of whether an install step ever runs, so setting it
  outside `if(SKBUILD)` — as an earlier revision briefly did — put
  CWD-lookup RUNPATH entries into `build/libLiPolGenCore.so` and the other
  in-tree `.so`s even though the in-tree flow never installs anything. With
  the fix, the non-`SKBUILD` configure sets no `CMAKE_INSTALL_RPATH` at all,
  so `build/lipolgen_tests` and the in-tree `build/python` module keep
  CMake's default build-tree RPATH exactly as before this packaging support
  existed (`readelf -d build/libLiPolGenCore.so` shows no RUNPATH entry, as
  at HEAD pre-packaging). `readelf -d` on the built wheel's `.so`s confirms
  `RUNPATH: $ORIGIN:<deps>/lib`; `ldd` (with `LD_LIBRARY_PATH` unset)
  resolves `libHepMC3.so.4`/`libpythia8.so`/`libLHAPDF.so` from the deps
  prefix alone.
- Data: `data/vmc` (1.1 MB) is installed to `lipolgen/data` inside the wheel
  under the `SKBUILD` branch; `python/lipolgen/__init__.py` points
  `$LIPOLGEN_DATA_DIR` at it when that env var is unset and the directory
  exists next to the installed module (no-op for the in-tree build, whose
  `build/python/lipolgen/` never has a `data/` sibling, so the compiled-in
  `$CMAKE_SOURCE_DIR/data` default keeps resolving exactly as before).
  PYTHIA8's `xmldoc` and LHAPDF's grids are **not** vendored (much larger,
  own licensing questions for the LHAPDF grids) — documented in
  `docs/USAGE.md` as coming from `$LIPOLGEN_DEPS_PREFIX` via
  `PYTHIA8DATA`/`LHAPDF_DATA_PATH`, same as `env.sh`.
- **Gate, measured 2026-09-02** (fresh venv, `numpy`+`pytest`+`pyhepmc`, all
  installed cleanly from PyPI):
  - `pip install -e .` with `LIPOLGEN_DEPS_PREFIX=…`: **66 s**, wheel
    `lipolgen-0.1.0-cp311-cp311-linux_x86_64.whl`.
  - From `/tmp` (no `env.sh`, no `PYTHONPATH`, only `PYTHIA8DATA`/
    `LHAPDF_DATA_PATH` exported): `import lipolgen; lipolgen.ion_spin('6Li')`
    → `1.0`; `python -m pytest python/tests -q` with the new
    `LIPOLGEN_TESTS_USE_INSTALLED=1` opt-out (`conftest.py`, default
    behaviour unchanged) → 144 passed, 2 skipped (missing optional `yaml`/
    `scipy` in the minimal venv — unrelated to packaging).
  - `pip wheel . --no-deps` + `python -m zipfile -l`: one `_lipolgen*.so`
    and one each of the four `libLiPolGen*.so`, all under `lipolgen/`, no
    duplicates.
  - In-tree flow re-verified bit for bit after all of the above:
    `build/lipolgen_tests` 276/276, `pytest python/tests` 146 passed.
- Portable wheel still needs `auditwheel repair` to vendor
  HepMC3/PYTHIA8/LHAPDF and rewrite RPATHs `$ORIGIN`-relative; documented in
  `docs/USAGE.md` together with the GPL-3 consequence of doing so (see §C in
  `engineering.md` and §13 below).
- License: GPL-3.0-or-later for LiPolGen's own code (the linked combination is
  GPL-3 regardless; permissive headers would mislead; MCnet-consistent).
  `pyproject.toml`'s `license = "GPL-3.0-or-later"` (SPDX string) matches.

## 14. Structure-function backends at the run surface — D2, `--unpol-sf` / `--pol-sf` — **OPT-IN SHIPPED**

**The finding (2026-09-03 sweep, D2).** Toy F₂/g₁/R were the shipped defaults
on every channel, and the tagged channels had no injection point for a real
backend. Confirmed in full on 2026-09-04 by reading every construction and
defaulting site of `InclusiveKernel`, `UnpolSF`, `PolSF`, `TensorSF` and
`RFunc` in the tree
(`docs/open_items/run_2026-09-03/phase_D_sf_injection.md`), with one
correction to the finding's own premise: **the documents had not missed it.**
`docs/PHYSICS_CHANNELS.md` §3 already said, verbatim, that on the tagged
channels `struck_cluster_kernel` exposes only `inclusive_b1` and `delta_func`
"so F₂, g₁, R and EMC there are hard-locked", and its EMC row already said
`Options::emc_ratio` "is empty in every default kernel". What was missing was
the *injection point*, not its documentation. D2 was a code task.

### What was built

Phase A's `B1UnpolSource` thread, copied piece for piece:

* **Two enums with a `Custom` value**, `UnpolSfSource {Toy, Mstw, Ct18Nlo,
  Custom}` and `PolSfSource {Toy, NnpdfPol, Custom}`, each with the whole
  "why it exists / what each value is / what the scope is and what it costs"
  block in the header (`include/lipolgen/pipeline.hpp`).
* **Two fields kept in step per selector**: the enum is the PROVENANCE (it is
  what `meta` records) and the object is the REALISATION. Assigning an object
  from Python sets `Custom`; clearing it sets `Toy`; `set_unpol_sf` /
  `set_pol_sf` go the other way.
* **The core library cannot build the named backends and does not try.**
  `MstwSF` is in the PYTHIA tier, `LhapdfSF`/`LhapdfG1` in the LHAPDF tier,
  and `sf.hpp`'s rule is that the core links neither, so the object is built
  one layer up in `python/bindings.cpp`, where the tiers are visible.
* **Six refusals in `validate()`**, all modelled on the b₁ ones: a named
  backend with an empty slot (naming the missing tier — never a silent
  fallback), `toy` with an object attached, either selector together with a
  caller-supplied `kernel`, a directly set `struck.f2_source`/`g1_model`
  under a `toy` selector, the one collision the two flag families
  create — `--b1-unpol toy` under a non-toy `--unpol-sf`, where "toy" would
  no longer name `ToyF2` while `meta["b1_unpol"]` still wrote `"toy"` — and,
  **added 2026-09-05**, a caller-supplied `kernel` on a **tagged** channel,
  which is never read there while `meta` wrote `"caller-supplied kernel"`
  into three keys (the `cluster_wave` rule applied to the escape hatch).
  A seventh refusal lives at the binding, where the two objects meet:
  `set_pythia_hadronizer` rejects a bridge whose `options.f2_source` is not
  the config's own `unpol_sf_obj`.
* **Four new `meta` keys, written unconditionally**: `unpol_sf`, `pol_sf`
  (both `"caller-supplied kernel"` when one is set), the pair
  `unpol_sf_grid_q2_min` / `unpol_sf_below_grid_frac`, and — **added
  2026-09-05** — `pol_sf_reach`. On the **coherent** channel `pol_sf` carries
  `not read on channel coherent-6Li` rather than a backend name, because
  nothing there evaluates g₁ (`pol_sf_is_read`); the below-grid fraction is
  taken against that channel's **own** rate (`cell_rate_weights_pb`), 0.447460
  and not the inclusive 0.361791. Both were knobs recorded as if they had run;
  see `docs/open_items/run_2026-09-03/phase_D_numbers.md` §D2.8.
* **Two kernels changed**: `default_inclusive_kernel` gained two trailing
  `nullptr` arguments (so the one-argument overload, which
  `validation/reference/b1_default_li6.json` pins at rtol 1e-12, is untouched)
  and `StruckClusterOptions` gained `f2_source` / `g1_model`. Nothing in
  `src/core/xsec.cpp` changed: `Toy` leaves both `Options` slots exactly as
  they were, so **the default path is unchanged by construction, not by
  inspection.**
* **The T2 bridge closed at the same time**, because it would otherwise have
  become a live inconsistency the moment the flag existed:
  `PythiaBridgeOptions::f2_source` was unbound and set by nothing, so the T2
  struck-nucleon species draw was on `ToyF2` while T0 and T1 would have moved.
  It is now bound, and both `lipolgen-run` and `lipolgen.run(hadronize=True)`
  hand it the same object they gave the kernel.

### What it costs to have been on the toy — measured 2026-09-04

⁶Li, `default_configs("6Li")[1]`, the shipped window, through the pipeline.
A_zz is the canonical thirds estimator on `tensor_thirds_plan(0, 0.6)`
divided by P_zz; A_∥ is `(σ₊−σ₋)/(σ₊+σ₋)/(P_e P_z)` on
`helicity_flip_plan(1, 0.7, 0.7)`.

| channel | backend | σ [pb] | A_zz | A_∥ |
|---|---|---|---|---|
| inclusive | `toy` | 591846.2 | −5.193231e−4 | −1.171716e−3 |
| inclusive | `ct18nlo` | 472571.9 (×0.79847) | −6.503970e−4 (×1.25239) | −9.668880e−4 (×0.82519) |
| inclusive | `mstw` | 469556.5 (×0.79338) | −6.545739e−4 (×1.26044) | −1.057030e−3 (×0.90212) |
| inclusive | `nnpdfpol` | 591846.2 (×1, exactly) | −5.193231e−4 (×1, exactly) | −1.270806e−4 (×0.10846) |
| tagged-6Li-α | `toy` | 591846.2 | **0 exactly** | −3.514730e−3 |
| tagged-6Li-α | `ct18nlo` | 472571.9 (×0.79847) | **0 exactly** | −2.900243e−3 (×0.82517) |
| tagged-6Li-α | `mstw` | 469556.5 (×0.79338) | **−1.3774e−16** (qual. 1) | −3.170626e−3 (×0.90210) |
| tagged-6Li-α | `nnpdfpol` | 591846.2 (×1, exactly) | **0 exactly** | −3.811967e−4 (×0.10846) |

**The A_zz column is at the shipped window, `Scenario::x_max` = 1.0** — it was
printed at `--x-max 0.95` until 2026-09-05, in the same rows whose σ and A_∥
were at 1.0, so one row carried two windows. At 0.95 the same three numbers are
−5.193192e−4 / −6.503956e−4 / −6.545724e−4; the ratios ×1.25239 / ×1.26044 are
the same to six digits in either window, which is why the error was invisible.

Three qualifications travel with that table and are not optional.

1. **The tagged A_zz is zero on every backend by construction — and on
   `mstw` that zero is a floating-point one, not a bit-level one.** At the
   shipped default the struck-cluster kernel carries no b₁, so the *total*
   per-category cross sections carry no tensor term. On `toy`, `ct18nlo` and
   `nnpdfpol` the three totals are then **equal to the last bit** and the
   estimator returns **0.0 exactly**. On `mstw` they are not: measured
   2026-09-05, `sigma_per_category_pb()` = 469556.45134814945 /
   469556.45134814945 / 469556.4513481495, i.e. σ₀ **one ULP** (5.82e−11 pb)
   above σ_±, so A_zz = **−1.3774e−16**. That is summation round-off —
   2.7e−13 of the inclusive channel's own toy A_zz — not a tensor signal, but
   "equal to the last bit" is false there and must not be written. The tagged
   tensor signal lives in the spectator-differential rate, not in σ_tot. With
   `--inclusive-b1 --x-max 0.95` A_zz is −1.557969e−3 (`toy`) → −1.951190e−3
   (×1.25239) → −1.963721e−3 (×1.26044), i.e. the same factors the inclusive
   channel shows.
2. **The A_∥ ratios are window integrals of a sign-changing integrand and
   are smaller than any of their own parts.** On the toy ⟨A_∥⟩ is −2.34618e−3
   for x < 0.01 (53.0 % of the rate) and positive above x = 0.05. Per x band
   the `nnpdfpol` swap is ×0.163 (x < 0.01), ×1.041, ×0.905, ×1.172, ×0.937 —
   the window-integrated **×0.108 is smaller than every band ratio** because
   the bands partly cancel. Quote the band, not the ×0.108, unless the window
   is the observable.
3. **The two selectors are not orthogonal.** `InclusiveKernel` builds its
   default `ToyG1` on its own base `UnpolSF`, so `--unpol-sf` alone moves g₁:
   at `--pol-sf toy`, `ct18nlo` moves ⁶Li F₁ by ×0.9942/×1.0711/×1.2038/×1.0014
   and g₁ by ×1.0544/×1.1346/×1.3050/×1.0551 at x = 0.05/0.10/0.30/0.50, so
   A₁ = g₁/F₁ moves 5–8 %. `--pol-sf toy` does **not** mean "g₁ unchanged".

**The neutron is a sign, not a factor.** `ToyG1`'s
a1n(x) = −0.07(1−x)² + 0.8x^2.2 crosses zero near x ≈ 0.25 and is positive
above it; NNPDFpol1.1's g₁ⁿ stays negative to x ≈ 0.6. At Q² = 10, g₁ⁿ is
−0.0800 / −0.0114 / **+0.00521** / +0.00779 (toy) against −0.1365 / −0.0608 /
**−0.02740** / −0.00037 (nnpdfpol) at x = 0.10/0.20/0.30/0.50. So **the
shipped toy g₁ⁿ has the wrong sign over roughly 0.25 < x < 0.6** — hidden on
isoscalar ⁶Li, not hidden at all on the neutron-tagged `TaggedDeuteronP`
channel. This is the single strongest argument for the polarised selector,
and it is a first measurement of this run.

**The grid clause.** CT18NLO's grid starts at Q² = 1.677 and the shipped
window's accepted cells start at 1.054; below its grid LHAPDF continues the
evolution downward rather than freezing, and F₂ᵖ at x = 3e−4 is a factor 2.34
below the toy at Q² = 0.7. **36.18 %** of a CT18NLO run's own accepted cell
cross section (42.33 % of the toy run's) is there. It is not refused —
refusing would make the flag unusable on the shipped scenario — but the
banner prints it and `meta` carries it, computed from the grid's own
`q2Min()` and written NaN (never 0) for a backend that reports no floor.

### What this deliberately does NOT cover, and why

* **R stays its own axis and no flag selects it.** Measured over the 3051
  accepted cells: `r_sigma_lt` spans 0.0046–0.1763 against `r1998`'s
  0.0124–0.4010, (1+r1998)/(1+r_sigma_lt) spans 0.9203–1.1910 with a
  σ-weighted mean of 1.1229, and **38.18 %** of the accepted cell cross
  section lies outside R1998's own stated support, where `r1998` returns a
  clipped boundary value. R1998 is therefore not a drop-in default, and
  shipping an `--r-model` in the same change as these two would make three
  independent 10–40 % movements unattributable. §10's close condition 4 (the
  ⁶Li R default) is untouched and still stands as decided.
* **The EMC hook is still empty in every kernel**, exactly as
  `PHYSICS_CHANNELS.md` already said. If one is ever wired, the transcribed
  EPPS21 depletion constant is quoted against **CT18ANLO** — which is not
  installed here — while `Epps21Ratio`'s own `proton_set` default is CT18NLO;
  referencing a ratio and a baseline to two different proton fits is the same
  4.2 %-class error that constant's own provenance note describes. That is
  also why `--unpol-sf` offers no `ct18anlo` row: the set is not on disk.
* **`--b1-unpol` stays a separate flag.** They are different quantities with
  different conventions (that one is the deuteron's per-nucleon F₁ inside
  CDKS Eq. (22) through `f1_cdks`, explicitly *not* `NuclearF2::f1a`), and it
  is a CDKS-comparability choice the §10 verdict rests on: a user who wants a
  realistic rate must not be forced to move the b₁ gate row off MSTW as a
  side effect. One refusal closes the ambiguity the two create.

### Adjacent defects found while mapping, and where they stand

* **`PipelineConfig::kernel` was silently ignored on tagged channels, and the
  `meta` then said `"caller-supplied kernel"` for a kernel that never ran** —
  **CLOSED 2026-09-05**. `Pipeline` reads `cfg_.kernel` on the non-tagged
  branch alone, and `validate()` now REFUSES a caller-supplied kernel on any
  tagged channel (`src/core/pipeline.cpp`, the sixth refusal listed under
  "What was built" above; `docs/USAGE.md` §7 carries the row). Nothing is
  filed for phase F on this. *(This bullet said the opposite — "This change
  does not close it … Filed for phase F" — until 2026-09-05, while the same
  section's "What was built" already recorded the refusal.)*
* **An `--isotope d --channel inclusive` run gets the ⁶Li rank-2 transfer**,
  because `default_inclusive_kernel`'s tensor branch keys on `ion.spin == 1`
  and the isotope guard lives only in `validate()`'s non-Miller branch.
  Documented behaviour (`PHYSICS_CHANNELS.md` says "for ANY spin-1 ion"), but
  a deuteron run is not a ⁶Li run. Unchanged here; a phase-F decision.
* **`PythiaBridgeOptions::f2_source` was unbound** — closed by this change,
  because it would have become a live inconsistency the moment the flag
  landed.


---

## 15. ⁷Li rank-2 (tensor) input — the zero made loud, the α–t b₁ deferred — **DEFERRED WITH DESIGN** (row 20 open)

**The finding (2026-09-03 sweep, D1; researched in
`docs/open_items/run_2026-09-03/phase_D_li7_rank2.md`, implemented and
measured in `phase_D_numbers.md` §D1).** A ⁷Li **inclusive** run's entire
rank-2 sector — the tensor term of the φ-averaged rate, the cos 2φ (gluon
transversity) amplitude, and therefore A_zz — is **exactly, bit-for-bit
zero**. `InclusiveKernel::tables` dispatches the rank-2 slots on the ion spin
(spin 1 → `b1_func / b2_func / delta_func`, spin 3/2 → `b1_32_func /
b2_32_func / delta_32_func`), an unset slot is 0.0 and an unset b₂ is 2x·b₁
and so 0 too; `default_inclusive_kernel` fills the **spin-1** slots only.
Measured: the two per-category cross sections of a T = +1 against T = −1 fill
are **the same double**, 590952.42641509 pb, and the asymmetry is exactly 0.0.

It was **silent**. The npz `meta` recorded `b1_model = "miller"` — a backend
that did not run — and the banner printed its b₁ block only when
`b1_model != Miller`, which a ⁷Li run can never reach. That is the
repository's own rule (*a knob that did not run may not be recorded as if it
had*, which `validate()` enforces for `b1_band_scale`) failing on the model
name itself.

### 15.1 What shipped (2026-09-04)

* **The banner**, unconditional on every ⁷Li inclusive run, tensor plan or
  not: the sentence, then what is unaffected (the unpolarised rate and the
  whole vector sector — measured, with a hand-supplied slot the sampler
  returns A_T = −0.043258 for b₁_32 = +0.05·F₁), then the two escape routes
  (`--isotope 6Li`; `--channel tagged-7Li-alpha`, whose α–t alignment **is**
  carried, in the event weight, gated at ⟨P₂⟩ = −T/5).
* **`meta["rank2_input"]`**, on every run and every channel, from the **one**
  C++ definition the banner prints (`rank2_input_report`) — so the claim
  cannot drift between the file and the log. Off the inclusive channel it says
  where that channel's tensor signal actually lives instead of borrowing the
  inclusive sentence.
* **`meta["b1_model"]` and `meta["b1_unpol"]` = `"none (spin 3/2: no rank-2
  input)"`** on such a run. The two **scales** stay numeric and stay at 1.0:
  `validate()`'s Miller branch already refuses any other value there, so
  neither can record a variation that did not run, and changing their type per
  isotope would break every consumer that reads them as floats.
* **Four run-surface defects.** F2: `make_plan` refuses the three spin-1-only
  tensor plans at J ≠ 1 (two of them used to *build* a spin-1 plan), and the
  `Pipeline` constructor now makes the run-plan-spin check on the inclusive
  branch that the tagged branch has always made. F3: `spin32_populations`
  prints the offending m and the **derived** domain. F4: a half-assigned
  `ClusterPartialWave` throws instead of segfaulting. F5: the run banner
  prints the fill's own moments, so `--plan helicity-flip --pzz 0.6` can no
  longer look as if it set T = 0.6 when the max-entropy ladder gave 0.4 —
  and since 2026-09-06 `--pzz-mode typed` will honour the typed value
  instead, refusing it (never clamping) where it is outside the plan's
  domain, which `--pzz 0.6` at `--pz 0.7` is at J = 3/2 by 0.02.
* **No b₁(⁷Li), in any form.** Nothing below is in the code.

### 15.2 Why the b₁ is deferred, and what would unblock it

The α–t construction is worked out, and it is **cleaner than ⁶Li's**: α is 0⁺
and the triton is ½⁺, so parity plus L ⊗ ½ ∋ 3/2 leaves **L = 1 alone** — one
partial wave, no interference, no node — and neither cluster carries rank-2, so
**b₁(⁷Li) is 100 % orbital**. Leading estimate at Q² = 2.5 with the shipped
`ToyF2`: b₁/nucleon = 1.246e−04 … 1.040e−04 over 0.05 ≤ x ≤ 0.70, i.e.
**2.99 ± 0.02 × ⁶Li's orbital term** — *larger*, not suppressed. **Every number
in this subsection was measured in `phase_D_li7_rank2.md` §§5–6 and is quoted
from it; none of it was re-measured when this item was written, and none of it
is in the code.**

**What blocks it is not the algebra. It is the sign.** Swapping the
unpolarised backend `ToyF2 → MSTW2008 LO` — the swap the A = 2 gate verdict
(§10) tells you to make before quoting a ⁶Li number — **flips the sign of
b₁(⁷Li) at four of six x points** and moves it by up to 348 %; a second
in-tree decomposition (⁷Li = ⁶Li + n) gives the opposite sign and up to 8× the
magnitude. Shipping that would publish a tensor asymmetry whose **direction is
a flag**, and no band expresses it: the band would have to contain zero and
both signs, at which point the number carries no information the loud zero
does not.

> **The unblocking decision is D2 below: one unpolarised backend, decided once
> for both isotopes.** It is the same question `--b1-unpol` asks for ⁶Li, where
> term (1) hides it; ⁷Li has nothing to hide it with.

### 15.3 The deferred design, in full

Write the α–t relative-momentum density in ion state M. One partial wave, so
the Clebsch–Gordan sum collapses to a single Legendre coefficient — **exact**,
for every k:

        n_M(k, k̂) = |φ₁(k)|² · A_M(k̂)
        A_M(k̂)   = Σ |⟨1 m_L ; ½ m_S | 3/2 M⟩|² |Y_{1 m_L}(k̂)|²
                  = (1/4π) [ 1 − Q_NN(M) · P₂(cos θ_k) ]                    (1)

with Q_NN(±3/2) = +1, Q_NN(±½) = −1 (`spin.hpp`; for J = 3/2 the plan's `pzz`
field carries **T**, not P_zz). Eq. (1) **is** the repository's own ⁷Li
polarimeter ⟨P₂⟩ = −T/5, which is a passing test today, and the identical CG
procedure reproduces the code's shipped ⁶Li coefficients 4.242641 and 1.500000
(`b1_nuclear.cpp`) — that is the load-bearing check on the normalisation.

Light-cone densities, `LightConeDensities`' own prefactor and δ-function:

        f_M(y)     = f₁(y) − Q_NN(M) · f₁^{P₂}(y)
        f₁(y)      = pre ∫ k dk |φ₁|²          (unpolarised α–t density)
        f₁^{P₂}(y) = pre ∫ k dk |φ₁|² P₂(c*)   (tensor density)             (2)

and, matching F₁^{(M)} = F₁ − Q_NN b₁ per nucleon,

> **b₁(⁷Li)(x, Q²) = (3/7) ∫ (dz/z) f₁^{P₂,t}(z) F₁^t(x/z, Q²)
>                  + (4/7) ∫ (dz/z) f₁^{P₂,α}(z) F₁^α(x/z, Q²)**           (3)

**Two terms, not ⁶Li's four**: there is no embedded spin-1 cluster, so ⁶Li's
dominant term (1) and its CG-depolarisation term (3) have no ⁷Li counterpart.
Kinematics {m_struck, m_recoil} = {M_t, M_α} and {M_α, M_t},
ε = 2.467 MeV (`LI7_ALPHA_TAG().separation_energy`); counting factors **3/7 and
4/7** replace ⁶Li's 2/6 and 4/6. Consistency check the implementation must
report: struck-t : struck-α = 1 : 0.756 against the analytic
[(4/7)/(3/7)]·(M_t/M_α)² = 0.757.

**What transfers from `b1_nuclear.hpp`**: `ConvolutionKinematics` unchanged
(only two masses and ε change; y_max = 1.6626 / 1.3761), the
`LightConeDensities` quadrature unchanged (converged: doubling every grid moves
b₁ by ≤ 0.05 %), `convolve` unchanged with `x_max_g = 1`, and Phase A's
`r_func` decision (`r_sigma_lt`) — which **matters more** here, worth
−9.66 % … −46.02 %, because *all* of b₁(⁷Li) is the orbital term that responds
to the slope of R. Since 2026-09-06 the shared-hook option (`--r-source`, §10)
exists for ⁶Li; a ⁷Li backend would inherit both the choice and its
consequences, including the **y dependence** the shared hook introduces and
the numerator-only swap does not.

**What must NOT transfer**: `LI6_B1_RANK2_TRANSFER` = 0.921947 — it is the
tensor dilution of an *embedded spin-1 deuteron*, ⁷Li has no spin-1
constituent, and the number that plays its role is the exact −1/5 of
⟨P₂⟩ = −T/5. Phase A's per-nucleon table constants
(`B1_CDKS_TABLE_TO_PER_NUCLEON`, `B1_MILLER_TABLE_TO_PER_NUCLEON`) **do not
arise**: they reconcile two published *deuteron* tables and there is no ⁷Li
table to reconcile. `LightConeDensities`' φ₀/φ₂ + SD/DD interface does not
transfer either — there is no S–D interference to split, and reusing `f_d` /
`f_d_p2` with the L = 1 wave in the φ₂ slot is numerically exact and
semantically a lie (decision D8).

**The nucleon input must be the run's own.** `--unpol-sf` (§14) already
reaches every ⁷Li kernel — measured: σ ×0.795936 (`ct18nlo`) and ×0.791318
(`mstw`) on the ⁷Li inclusive channel, ×0.792543 / ×0.788565 on tagged-⁷Li-α,
with `meta` recording which — so the unpolarised backend a ⁷Li b₁ folds
against **is** selectable today. `--pol-sf` reaches those same kernels, but
only where the fill carries `lam_e · P_e ≠ 0`: on a tensor plan it is
`not-read` and `meta` says so rather than naming a backend (§14, the
knob-provenance table). Since b₁ is a tensor observable, that is the plan a
⁷Li b₁ would be measured on — so a ⁷Li b₁ must take its polarised input from
the field, never from `meta["pol_sf"]`, which is a label there. A `Li7AtConvolutionB1` must therefore share
the kernel's own `f2_source` object exactly as the `Li6Convolution` branch
does, and `validate()` must refuse the same mislabelling collision `--b1-unpol
toy` under a non-toy `--unpol-sf` is refused for.

### 15.4 Step 2, before any b₁: the A = 7 quadrupole gate

**There is no A = 3 analogue of the A = 2 gate, and the reason is a selection
rule**: the A = 3 analogue of α + t would be ³H or ³He, which are **J = ½** and
have **no rank-2 structure function at all**. Nothing about the ⁷Li
convolution can be validated by running it on a lighter system. That must be
said in any header that ever ships this model.

What *does* exist is a measured rank-2 observable of the same nucleus that the
same wave function and the same CG coefficient predict with **no free
parameter**. Both clusters have zero intrinsic quadrupole, so the whole of
Q(⁷Li) is orbital:

        Q = Z_eff ⟨r²⟩ ⟨3cos²θ − 1⟩_{M=3/2} = − (2/5) Z_eff ⟨r²⟩
        Z_eff = Z_α (M_t/M₇)² + Z_t (M_α/M₇)² = 0.695075                    (4)

**COMMITTED 2026-09-06 (D11 paid; `run_2026-09-06/phase_B_numbers.md` §B3).**
It is now shipped code, not a research measurement:
`li7_alpha_t_quadrupole` (`b1_nuclear.hpp`) reads
`momenta/li7_at3.momentum` through `li7_alpha_channel` and returns
⟨r²⟩ = **12.534828961030 fm²**, r_rms = **3.540456038568 fm**,
Z_eff = **0.695075101362** and **Q = −3.485059004257 fm²** — reproducing the
research note's own numpy recipe (§9 of `phase_D_li7_rank2.md`) to the **last
bit** (relative difference 0.000e+00). Converged to **0.032 %** in r_max
(30 → 60 fm) and **±0.33 %** on the fully correlated VMC MC band.

The reference is now in the tree, with one home:
**`LI7_QUADRUPOLE_FM2` = −4.06 fm²** in `rc.hpp`, beside `LI6_QUADRUPOLE_FM2`
— sourced from **TUNL's A = 5, 6, 7 evaluation** (Tilley *et al.*, NPA 708
(2002) 3), whose A = 7 half prints `Q = −40.6 ± 0.8 mb (1988DI1B)` and whose
A = 6 half prints the `Q = −0.818(17) mb (1998CE04)` that IS
`LI6_QUADRUPOLE_FM2`. One document, one sign convention, one unit rule; the
−4.00(3) fm² the research note quoted from memory is one of the **eight**
⁷Li determinations N. J. Stone lists without recommending any, and it is
recorded in `rc.hpp` beside the adopted value with its price:

> **ratio 0.858389** against −4.06 (**14.16 % low**), **0.871265** against
> −4.00 (**12.87 % low**) — so "13–14 % low" is the verdict and the choice of
> compilation (1.5 %) does not decide it.

Against the ⁶Li contrast from the same code family — factor **4.0748**, i.e.
**307.5 %** — ⁷Li is better by a factor **21.7**, which is what "an order of
magnitude better" means here and is asserted with a threshold of 10.
Pinned by `tests/test_b1_nuclear.cpp` **T13** (26 assertions) and
`python/tests/test_li7_rank2.py` **G8** (4 tests).

**What it is NOT.** A **reported ratio**, never a pass/fail on b₁(⁷Li) — which
remains **unimplemented** and blocked on D2 / registry row 3, with ⁷Li's whole
rank-2 sector still exactly zero. G8's
`test_the_passing_wave_function_gate_ships_no_b1_at_all` runs the gate and
then re-asserts that zero in the same test, so a passing wave-function gate
cannot be read as licensing a b₁. It validates the α–t **wave function's**
quadrupole — its ⟨r²⟩ and P-wave character — and the wave function
`li7_vmc_waves` already ships to the tagged channel; nothing more.

Its honest limits, to travel with it: S_αt = 1.0084 > 1 proves the α–t overlap
is **not a probability**, and the missing 13 % is the expected size of cluster
distortion and non-α–t components; and it gates ⟨r²⟩ and the P-wave character,
**not** the light-cone convolution or the DIS input, which is where the real
uncertainty lives.

### 15.5 The author decisions this needs

*(Registry: these thirteen are carried in `run_2026-09-03/STATUS.md`'s decision
table as **row 20**, and drafted in `run_2026-09-03/AUTHOR_DECISIONS.md` §B3
(D2), §B18 (D13) and §B20 (the rest). Cite them as "§15.5 D<n>": "D1" and "D2"
each name three different decisions across this document set.)*

| # | decision | what is at stake |
|---|---|---|
| **D1** | **Which decomposition** — α + t (S = 1.008, purely orbital) or ⁶Li + n (S = 0.682, carries a core b₁) | They disagree in **sign** and by up to 8×, and they are non-orthogonal, so summing them double counts. α + t is the one already in the C++ (`li7_vmc_waves`, `li7_alpha_channel`), the one with a k-space table, and the one with a quadrupole gate |
| **D2** | **The unpolarised backend the convolution folds against** — `ToyF2` or MSTW2008 LO | **Decides the sign** at four of six x points. This is the blocker, and it should be settled **once, for both isotopes**, not twice |
| **D3** | **The normalisation target** — renormalise to S_αt = 1.0084, or to 1 | 0.83 % numerically, so a **provenance** decision, not a magnitude one; but "the file's own S" and "1" are different claims about ⁷Li, and S > 1 means the ⁶Li rule ("the non-α–d 18 % is given b₁ = 0") has no ⁷Li form |
| **D4** | **F₁ of the triton** — isoscalar (the ⁶Li shortcut) or the true Z = 1, N = 2 | 0.6 … 8.6 %. ⁶Li's `f1_alpha_ = f1_d_` is *correct* for a deuteron and **wrong** here; copying it would be a silent error, not a documented approximation. `triton_sf.hpp` already exists |
| **D5** | **κ: CDKS Eq. (17) or Eq. (21)** | +0.35 … **+69 %** here against −1 … +7 % on ⁶Li. The ⁶Li default's stated justification ("the term κ multiplies is the small orbital one") **inverts** for ⁷Li, where it is the only term |
| **D6** | **Which A_zz** — `A_T` = −b₁/F₁ or `A_zz^{(3/2)}` = −(2/3) b₁/F₁ | A factor 3/2, and the spin-1 `azz()`'s explicit `2.0/3.0` has **no J = 3/2 counterpart**. Must be chosen before any ⁷Li tensor number is published |
| **D7** | **`ClusterPartialWave` refuses odd L** | For a single wave the i^L phase is unobservable, so the refusal protects nothing — but the **type** must be taught that rather than worked around |
| **D8** | **`LightConeDensities`' interface** | Its φ₀/φ₂ + SD/DD shape is an A = 2 / ⁶Li shape. A first-class ⁷Li needs an L-generic alignment slot or an explicit `alignment_coefficient` (1 for L = 1 / S = ½ / J = 3/2; 1.5 and 6/√2 for ⁶Li's DD and SD), which would also let the ⁶Li coefficients be **derived** instead of hard-coded |
| **D9** | **b₂_32** | Silence means 2x·b₁ by default; state it rather than inherit it |
| **D10** | **Δ_32 (cos 2φ)** | Filling `b1_32_func` alone leaves cos 2φ at zero. There is **no ⁷Li Δ model**; reusing ⁶Li's toy means adopting an arbitrary 1e−2 scale for a second nucleus — and 3·Q_NN = ±3 at J = 3/2 against ±1/−2 at spin 1, so the same Δ gives a **larger** ⁷Li amplitude |
| **D11** | ~~**Q(⁷Li) is not in the tree**~~ — **PAID 2026-09-06** | §15.4. `LI7_QUADRUPOLE_FM2` = **−4.06 fm²** at `rc.hpp:669`, from TUNL's A = 5, 6, 7 evaluation (NPA 708 (2002) 3), the same document `LI6_QUADRUPOLE_FM2` comes from; the −4.00(3) vs −4.06 spread and all eight of Stone's ⁷Li entries are recorded in the constant's comment, and the gate is `li7_alpha_t_quadrupole` + doctest T13 + pytest G8 (ratio 0.858389 / 0.871265). It gates the wave function; **b₁(⁷Li) is still unimplemented and D1–D10, D12, D13 are untouched** |
| **D12** | **The run-plan surface** | There is no J = 3/2 tensor plan. The honest one is the **two-state T = +1 / T = −1 contrast**, not a thirds pattern — and at J = 3/2 the pure-alignment fill is rank-3 clean by construction |
| **D13** | **`--pzz` on `helicity-flip`** (F5; **priced 2026-09-06**) | At the default `--pzz-mode ladder` it is still not read — `use_explicit_pzz` stays false and the fill comes from the max-entropy ladder at `--pz` (T = 0.4, not the 0.6 typed) — but the run banner names both fills and the provenance table's `pzz_mode` row says which branch ran. **`--pzz-mode typed` is registry option (iii), built and priced**: it honours the typed value and refuses one outside the plan's domain with the edge named, never clamping. WHAT IT COSTS, measured at the standard configuration (inclusive, config 1, seed 20260713, `--pz 0.7 --pe 0.7`, `--pzz 0.5` against the ladder; `phase_A_numbers.md` §A3): on **⁷Li — the isotope this row is about — no observable moves.** The rank-2 sector is identically zero and both modes honour `--pz`, so the only fill moments the kernel sums are unchanged: σ agrees to 1 ulp (1.970 × 10⁻¹⁶ relative), the two per-category σ to 3.939 × 10⁻¹⁶ and 0, and A_∥ moves by 6.2 × 10⁻¹⁴ of its own statistical error. On **⁶Li** (J = 1, where the rank-2 sector is live) σ moves −0.0023527 % and A_∥ by −4.27 × 10⁻⁶ σ_stat at 100 000 events. What moves on both is the **recorded alignment** — 0.4 → 0.5 at J = 3/2, 0.409403 → 0.5 at J = 1 — which is the divisor of every tensor estimator (−20 % / −18.1195 % on δ(A_zz) and δ(cos 2φ) at fixed N), and the event sample, which is not bit-identical even on ⁷Li (0.1060 % of 100 000 events change (x, Q²) cell on last-bit arithmetic). The DEFAULT did not move |

### 15.6 One thing not to do

Do **not** reuse `LI6_B1_RANK2_TRANSFER` = 0.921947 for ⁷Li. It is
1 − (9/10)·P_D(α–d), the tensor dilution of an *embedded spin-1 deuteron*
inside an L = 2 component. ⁷Li has no spin-1 constituent; the number that
plays the analogous role is the −1/5 of ⟨P₂⟩ = −T/5, and it is **exact**
rather than model-dependent.
