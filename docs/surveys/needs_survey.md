I have comprehensive coverage. Writing the report now.

# Survey: PolarizedLithiumSim — event generation, tooling, and what a new dedicated generator must satisfy

Paths are relative to `/home/cpeng/Projects/polli/PolarizedLithiumSim`. Nothing was modified.

---

## 0. Orientation

Two Python packages, no C++ of the project's own:

- `fastsim/polli_fastsim/` — analytic fast simulation (structure functions, asymmetry formulas, beams, far-forward routing, FOM maps). 11 modules, ~3.3 kLOC.
- `evgen/polligen/` — the event generator, which **imports** `polli_fastsim` rather than duplicating it (`evgen/polligen/__init__.py:16-20` inserts `fastsim/` on `sys.path`). 12 modules, 6.6 kLOC, 305 tests.
- `tools/{pythia8,beagle,fullsim,analysis}` — third-party generator/full-sim glue.

Self-description (`README.md:29`): polligen is *"the doubly polarized e+⁶,⁷Li event generator (the first of its kind)"*. The build-vs-reuse rationale is `plans/05_doubly_polarized_generator.md:18-30`: *"**no public generator produces polarized e+A events for A > 1** — BeAGLE is unpolarized; DJANGOH's hadron polarization is nucleon-level and *longitudinal only*; PEPSI/CLASDIS are polarized nucleon-level LO … So we reinvent the wheel — but *only the wheel that does not exist*: the polarized-nucleus vertex … and the spin-correlated cluster-spectator sampler."*

---

## 1. Physics implemented by `evgen/polligen`

### 1.1 Module inventory

| file | lines | purpose (one line) |
|---|---|---|
| `polligen/__init__.py` | 28 | package doc + `sys.path` wiring to `fastsim/polli_fastsim`; declares scope "Step 5.A … tagged (5.B), reweighting (5.C) and HepMC3 output (5.D) build on these" |
| `polligen/spin.py` | 257 | spin-density matrices ρ(m,m′) for J=1, 3/2: Wigner-d, Clebsch–Gordan, multipole operators, populations ↔ (vector, tensor, octupole) moments, arbitrary quantization axis, spin-temperature fills |
| `polligen/xsec.py` | 395 | the **inclusive doubly polarized master cross section** (`InclusiveKernel`) — HJM spin-1 tensor sector, vector sector, Δ cos 2φ, spin-3/2 rank-0/1 exact + rank-2 scenario, optional exact finite-γ A∥ |
| `polligen/bookkeeping.py` | 207 | run plans: `SpinCategory` (fill populations × electron helicity × luminosity share), helicity-flip/tensor-thirds/transverse/tensor-flip plans, relative-luminosity bias formulas, polarimetry smearing, per-(run,bunch) RNG |
| `polligen/sample.py` | 366 | `InclusiveSampler` — grid inverse-CDF in (x,Q²) on the *polarized* rate, φ by accept–reject, per-spin-state Poisson counts, Mode-W weight matrix, exact φ-histogram pseudo-experiments |
| `polligen/estimators.py` | 124 | analysis-side estimators: helicity flip, tensor thirds, cos 2φ moment and binned LSQ fit (acceptance-hole robust), lumi-corrected yields, pulls |
| `polligen/tagged.py` | 574 | two-cluster light-front IA ⊗ spin: partial-wave amplitudes, m-dependent spectator densities n_M(k,k̂), struck-cluster spin populations, lab boost, far-forward routing, `TaggedSampler` |
| `polligen/coherent.py` | 425 | coherent intact-⁶Li channel: scenario coherent fraction/slope, exponential t, tag acceptance, deformation + flat-gluonic cos 2φ, breakup/veto tables |
| `polligen/reco.py` | 1596 | measured quantities: head-on↔lab frames (25 mrad crossing), covariant azimuths, electron/hadronic/mixed reconstruction, detector resolutions, spin-state harmonic ratio & likelihood estimators, Roman-Pot emulation |
| `polligen/recopseudo.py` | 1228 | reconstructed-level pseudo-experiments: `RecoResponse` (importance-sampled), `measure_inclusive`, forward-folded Δ response, `CoherentResponse` + `measure_coherent` |
| `polligen/hfs.py` | 886 | hadronic final state: `HFSSample` format, exact hadronic sums + JB/Σ/DA/mixed kinematics, `HadronResponse`, `ToyHFS` fragmenter, `HFSLibrary`/`HFSResponse` transfer |
| `polligen/radiative.py` | 732 | collinear ISR as a **migration bound**: exponentiated leading-log D(z,Q²), sampler, closed-form observed kinematics per method, E−p_z handle, `migration_bound` |
| `polligen/nearbeam.py` | 384 | thin-film energy deposit, Landau sampling, nanowire hot-spot threshold, Roman-Pot Z-ID fake rates (detector study, not event physics) |

### 1.2 The spin-density kernel

**Master formula** (`evgen/polligen/xsec.py:7-17`, docstring):

```
dsigma/(dx dQ2 dphi) = sigma_unpol(x,Q2)/(2 pi) * W(phi')
W = 1 + w_avg + a_1 cos(phi') + a_2 cos(2 phi'),   phi' = phi - phi_S

w_avg = t_geo(m, theta_S) * K / D_phi                      [tensor, J>=1]
      + lam_e P_e (m/J) cos(theta_S) * A_par(x,y)          [vector-L]
a_1   = lam_e P_e (m/J) sin(theta_S) * A_perp(x,y)         [vector-T, gT]
a_2   = -(1-y)/y^2 * c_eff(m) sin^2(theta_S) * Delta / D_phi  [gluonometry]
```
with `D_phi = F1 + (1-y)/(x y^2) F2` and `K = b1 + (1-y)/(x y^2) b2`.

Implementation: `InclusiveKernel.amplitudes` at `xsec.py:345-379`; `_dphi` at `:252-254`, `_tensor_kernel` at `:256-258`, and the actual lines

- tensor rate: `t_geo = TENSOR_LL_SIGN * q_nn * 0.5 * (3.0*ct*ct - 1.0)` (`xsec.py:365`), `w_avg += t_geo * kern / D_phi` (`:367`)
- double-helicity-flip: `a_2 = -(1.0-y)/(y*y) * c_eff * st*st * t["delta"] / D_phi` (`xsec.py:368-369`)
- vector: `w_avg += helicity * v * ct * a_parallel` (`:374`), `a_1 = helicity * v * st * a_perp` (`:378`)
- final density: `w = 1 + w_avg + a1 cos(phi') + a2 cos(2 phi')`, `dsigma_unpol/(2π) * max(w,0)` (`xsec.py:388-395`).

**Rank-2 geometry, one formula for both spins** (`xsec.py:311-341`), quoting the code comment: *"Q_NN(m) = [3 m² − J(J+1)] / 3 … T_LL = Q_NN P₂(cos Θ) and T_TT = (3/2) Q_NN sin²Θ. Returning (Q_NN, 3 Q_NN) therefore gives ONE geometry for both spins."* For J=1 this equals the HJM `c_m = 3m² − 2` transcription digit-for-digit; for J=3/2 it *changed the cos 2φ channel by a factor 3* relative to earlier code (`xsec.py:330-335`).

**Tensor sign is one constant.** In the SIBLING repository — `fastsim/` is a different tree, outside LiPolGen, and this whole survey is about it — `fastsim/polli_fastsim/asymmetries.py:41` reads `TENSOR_LL_SIGN = +1.0`, the OPPOSITE of the value this tree ships (`include/lipolgen/constants.hpp` `TENSOR_LL_SIGN = -1.0`, the Cosyn/HERMES convention, set 2026-08-29), with the comment at `:24-40` that `+1` gives `Azz = +(2/3) b1/F1` (the program's HJM transcription) and `-1` the Cosyn et al. EPJ A 61 (2025) 83 / HERMES convention `Azz = -(2/3) b1/F1`. **This is an unresolved author decision** (`plans/08_simulation_chain_completion.md:373`, item D1). The Δ sector does not depend on it. *(STALE against the sibling since 2026-08-29: `fastsim/polli_fastsim/asymmetries.py:66` now reads `TENSOR_LL_SIGN = -1.0`; :41 is a citation comment. Noted 2026-09-05, not rewritten — this is a dated survey.)*

**Analytic asymmetry formulas the kernel must reproduce** (`fastsim/polli_fastsim/asymmetries.py`):

- `azz` (`:101-108`): `num = TENSOR_LL_SIGN * (2/3) * (b1 + (1-y)/(x y²) b2) * ½(3cos²θ_m − 1)`, divided by `phi_averaged_density`.
- `a_cos2phi` (`:111-113`): `-(1-y)/y² * Delta / D_phi`.
- `a_parallel` (`:59-93`): `D(y) * g1/F1` with `depolarization_d = y(2−y)/(y² + 2(1−y)(1+R))` (`:44-56`).
- errors: `err_a_parallel = 1/(P_e P_z √N)` (`:119-122`), `err_azz = √(2/N)/P_zz` (`:125-129`), `err_cos2phi_amplitude = √(2/N)/P_zz` (`:132-136`).

**Finite-γ / target-mass**, default OFF (`xsec.py:118-146`, `:260-292`): `gamma² = 4M²x²/Q²`, E143 `eps`, `D_gamma`, `eta`; `A_par = D_gamma (A1 + eta A2)` with `A1 = (g1 − γ²g2)/F1`, `A2 = γ(g1+g2)/F1`. `g2` defaults to Wandzura–Wilczek, implemented in `g2_ww` (`xsec.py:85-102`). The W² ≥ 10 GeV² cut bounds `γ² ≤ M²/(W2_min − M²) = 0.0965` (`xsec.py:50-52`).

**Transverse vector amplitude** `a_perp` (`xsec.py:294-309`): `d(y)·γ·((y/2)g₁ + g₂)/F₁`, keeping both O(γ) pieces.

### 1.3 Structure-function inputs (all from `fastsim/polli_fastsim`)

- **b₁ (deuteron, per nucleon)**: two published camps, digitized, not toys — `polarized.toy_b1` (`polarized.py:451-470`) = Miller PRC 89:045203 Fig. 5 ("HERMES-like"); `polarized.b1_convolution` (`:473-492`) = Cosyn–Dong–Kumano–Sargsian PRD 95:074036 Fig. 4 (|b₁| < 10⁻³ at x ≳ 0.2). Close–Kumano sum rule `∫b₁dx = 0` is **reported, not enforced** (`:495-506`).
- **b₁ for ⁶Li**: `b1_li6_from_deuteron` (`polarized.py:525-533`): `b1(6Li)/nucleon = transfer × (2/6) × b1(d)/nucleon` with `LI6_B1_RANK2_TRANSFER = 0.921947` (`:515`) — which *is* `TaggedModel(li6_alpha_channel()).tensor_dilution()`, i.e. `1 − (9/10)P_D`. Comment at `:529-530`: *"No published 6Li b1 exists (plans/04 #9), so the embedded-deuteron scaling is our own inference."*
- **Δ(x,Q²)**: `polarized.toy_delta_gluon` (`:536-544`) = `scale·F1·x^0.3·(1−x)^4`; and the production registry `fastsim/polli_fastsim/delta_models.py` — `moment_A` (`Δ = A·α_s(Q²)·F1·shape`, A solved so `∫ x Δ dx = C_BAG·α_s` with `C_BAG = -0.012` from Sather–Schmidt PRD 42:1424, `delta_models.py:78-80`) and `moment_B` (same without F1). Shape variants `low_x/mid_x/high_x` = (0.3,4.0)/(0.7,3.0)/(1.5,2.0) (`:83-87`). ⁶Li per-nucleon dilution 1/3 (`:39-52`).
- **⁷Li polarized EMC**: `polarized.cbt_polarized_emc_ratio` (`:312-335`) = digitized Cloët–Bentz–Thomas PLB 642:210 Fig. 6 (⁷Li at Q²=5, Eq. 23 or 26); `tmt_polarized_emc_ratio` (`:338-...`) = Tronchin–Matevosyan–Thomas PLB 783:247 (nuclear matter at Q²=10) **transferred to ⁷Li** by one valence-window strength factor `TMT_VALENCE_SCALE = 0.397009` over `POLEMC_VALENCE_WINDOW = (0.35, 0.65)` (`:285-299`). Digitized range x = 0.028–0.871; outside it the curve is frozen (`:325-328`).
- **g₁**: `ToyG1` (A₁ shapes, `polarized.py:26-82`) or `PartonG1` on NNPDFpol11_100, 3 flavours (`:85-131`). Nucleus: `g1A = Z·P_p·g1p + N·P_n·g1n` (`:78-79`), with per-nucleon effective polarizations in `beams.Ion` — `LI6 = (1/3, 1/3)`, `LI7 = (0.866/3, −0.037/4)` (`beams.py:197-199`).
- **F₂**: `ToyF2` (`structure.py:41-60`) or `PartonF2` on CT18NLO, 5 flavours (`:224-300`); nucleus `F2A = Z F2p + N F2n` × optional EMC ratio (`NuclearF2.f2a`, `:322-326`); `F1 = F2/(2x(1+R))` (`:328-338`); cross section `d²σ/dxdQ² = 4πα²/(xQ⁴)[(1−y+y²/2)F₂ − (y²/2)F_L]` (`:341-359`).
- **Nuclear PDFs**: `EPPS21nlo_CT18Anlo_Li6` and `nNNPDF30_nlo_as_0118_A6_Z3` are installed and used **only in the dated `fastsim/scripts/money_delta_2026*.py` suite** (`money_delta_20260729.py:407` `EPPS21_SET = "EPPS21nlo_CT18Anlo_Li6"`), not in polligen's kernel path, which builds F2A from free-nucleon CT18NLO or the toy.

### 1.4 Nucleon / cluster momentum distributions and spectator tagging

**Kinematic layer — `fastsim/polli_fastsim/spectator.py`:**

- Two-cluster channels (`:187-209`): `LI6_ALPHA_TAG` (6Li → d struck + α spectator, S = 1.4743 MeV, L=0), `LI6_D_TAG`, `LI7_ALPHA_TAG` (7Li → t struck + α spectator, S = 2.4670 MeV, L=1), `LI7_T_TAG`, plus `DEUTERON_P_TAG`/`DEUTERON_N_TAG`/`HE3_P_TAG` controls.
- Momentum densities (`:212-219`): S-wave Hulthén `ψ ∝ 1/(k²+κ²) − 1/(k²+β²)`; P-wave `ψ ∝ k/((k²+κ²)(k²+β²))`. `β = 0.30` default, scan 0.20–0.40 as the model band (`:20-21`).
- `κ = √(2μS)` with μ the reduced mass of the two **free** clusters (`:168-184`): 60.66 MeV (⁶Li→α+d), 88.90 MeV (⁷Li→α+t).
- Masses are physical AME2020-derived nuclear masses (`NUCLEUS_MASS`, `:75-84`), not A·M_U; this moves the rigidity ratio at k=0 to `R = 0.99813` (⁶Li α) and `0.85571` (⁷Li α) (`:29-31`) — a 2×10⁻³ effect that nevertheless shifts the ⁶Li α Roman-Pot fraction by +15% relative because the R = 0.95 window edge is hard (`:294-299`).
- `sample_k` (`:222-235`): inverse-CDF of `k²·n(k)` on a 30 000-point grid, isotropic angles.
- `_boost_fragment` (`:238-277`) returns `{pT, theta, phi, p_lab, R, xL, k}`.
- Docstring is explicit that this is *"crude BY DESIGN — the high-k tail is the dominant model uncertainty and exactly what VMC densities/BeAGLE must eventually pin down"* (`:9-11`), with *"Upgrade path: replace n(k) with VMC two-cluster overlap densities (R.B. Wiringa et al., ANL) — same interface"* (`:35-36`).

**Spin layer — `evgen/polligen/tagged.py`:**

- The ion state is expanded in cluster relative partial waves; everything follows from the joint amplitude `A_{m_S}(M; k, k̂) = Σ_L a_L ψ̂_L(k) C_L(M,m_S) Y_L^{M−m_S}(k̂)` (`tagged.py:5-23`). From it: spectator density `n_M(k,k̂) = Σ|A|²`, struck-cluster spin populations `p(m_S|M,k,k̂)`, and pair decomposition via a second CG factor.
- `Wave` radial forms (`:126-149`): L=0 Hulthén, L=1 P-wave, L=2 D-wave `k²/((k²+κ²)(k²+β²)²)`.
- Channels (`:181-200`): `li6_alpha_channel` = S+D with `P_D_LI6 = 0.0867` (chosen so the vector dilution `1 − (3/2)P_D` reproduces 0.87, `:174-176`, flagged **SCENARIO**, VMC overlaps the scheduled replacement); `li7_alpha_channel` = pure P-wave onto `TRITON = Ion("t",3,1,0.5, eff_pol_p=0.86, eff_pol_n=-0.028)` (`:99`); `deuteron_channel` = S+D with `P_D_DEUTERON = 0.045`.
- `TaggedModel` (`:203-332`) builds `|A_{m_S}|²` on a (280 k-points × 96 cosθ) grid, exposes `n_of_kc`, `struck_populations`, `population_integrated`, `vector_dilution`, `tensor_dilution`, `p2_moment`, and `sample_kc`.
- **Documented physics consequences** (`:42-52`): ⁷Li α-tag `n_{3/2} ∝ sin²θ_k` with `⟨P₂(cosθ_k)⟩ = −T/5` (in-situ alignment polarimeter), triton polarization `P_t(M=3/2)=1`, `P_t(M=1/2)=1/3`; ⁶Li α-tag vector dilution `1 − (3/2)P_D`; deuteron S/D interference tensor structure.
- **Coherences are dropped** — *"diagonal truncation: coherences feed phi-dependent structures beyond the Step-5.A master formula and are dropped, documented"* (`:19-22`). **No FSI** — *"(IA, no FSI — quote at low spectator virtuality; plans/04 #16)"* (`:24-25`). **No light-cone α_s dependence of the struck-cluster SFs** — *"beyond this tier (Cosyn-Weiss upgrade path)"* (`:33-34`).
- `boost_spectator` (`:335-388`): spin-frame → lab rotation, longitudinal beam boost, returns `{pT, theta, p_lab, R, xL, kx, ky, kz, phi_spec}`. `phi_spec` is the spectator's **lab** azimuth, needed for the rectangular 10(σ_h,σ_v) Roman-Pot envelope; without it the cut degenerates to a circle and overstates the tag by 1.7× (`:349-351`).
- `TaggedSampler.sample_category` (`:507-569`): joint sampling M → m_S → (k,k̂) → (x,Q²,φ) from the inclusive sampler conditioned on m_S → lab boost → `route_charged`.
- The tagged cross section factorizes as `dσ^tag(λ_e, M) = |A_{m_S}(M;k)|² ⊗ dσ_struck(λ_e, m_S)` (`:27-28`), with the struck cluster carrying its **own** SFs (embedded deuteron: F₁d, F₂d, b₁d, g₁d — *"this is the embedded-b₁ observable"*; quasi-free triton: g₁t).
- Explicit reconciliation with the other sampler (`:54-70`): the D wave is the whole difference, and *"The D wave is not an optional tail: it IS the tensor observable — with P_D = 0 the alpha-d density is m-independent and A_zz^tag vanishes identically."*

### 1.5 How the hadronic final state is attached

There is **no** hadronization inside the generator. The attachment is a **response transfer**, not an event-level graft (`evgen/polligen/hfs.py:29-41`):

1. `HFSSample` (`hfs.py:68-164`) — a CSR-style container of an *unpolarized* generator's final states.
2. `hadronic_sums` (`:168-187`) → per-event `Σ = Σ(E − p_z)`, `p_Tx`, `p_Ty`; `hadronic_kinematics` (`:190-220`) → the HERA method set (JB, Σ, double-angle, mixed): `y_jb = Σ/2E_e`, `y_sigma = Σ/(Σ + E'(1−cosθ))`, `x_mixed = Q²_e/(s·y_Σ)`, etc.
3. `HadronResponse` (`:238-423`) — ePIC-like tracker/EMCal/HCal acceptance, thresholds, efficiencies, resolutions, calorimeter noise; resolutions are Yellow Report requirement values, everything else labelled *"this programme's own stand-ins"* (`:241-247`).
4. `HFSLibrary` (`:610-771`) bins the library in (log x, log Q²) and, for each pseudo-event, draws a library event from the same cell and transfers its **captured Σ fraction** `f = Σ_reco/Σ_true`, **p_T ratio** `r`, and azimuthal shift `dφ_pT`.
5. `HFSResponse.hadronic` (`:807-847`) applies `f`, `r` to the pseudo-event's exact truth Σ (including the target-mass term `m_N²/(E_N+p_N)`), adds noise, optionally divides by the **reconstructed-bin** calibration map `cell_means_reco`.
6. `ToyHFS` (`:462-585`) — a vectorized string-fragmentation stand-in (uniform-in-rapidity, Gaussian p_T, exact 4-momentum closure by longitudinal bisection, π⁰→γγ) *"so that the chain runs and is tested on a machine without PYTHIA"* (`:26-28`).

The production library is PYTHIA 8 via `tools/pythia8/gen_dis_hfs.py`. Because the tensor polarization enters only the inclusive azimuthal weight, an unpolarized HFS source is deemed correct (`docs/reproduction_manual.md:1280-1281`).

---

## 2. How PYTHIA 8 and BeAGLE are used today

### 2.1 PYTHIA 8 — `tools/pythia8/gen_dis_hfs.py` (307 lines)

- **Bindings**: native pybind11 built from PYTHIA 8.311 source against the analysis interpreter (`tools/pythia8/README.md:21-44`); `import pythia8` fails inside eic-shell, which has C++ library and headers only.
- **Settings**, verbatim (`gen_dis_hfs.py:180-211`): `Beams:frameType = 2`; `Beams:idA = 2212|2112`; `Beams:idB = 11`; `Beams:eA = <p_per_nucleon>`; `Beams:eB = <E_e>`; `WeakBosonExchange:ff2ff(t:gmZ) = on`; `PhaseSpace:Q2Min = 0.7`; `PhaseSpace:pTHatMinDiverge = 0.5`; `PhaseSpace:mHatMin = 0.5`; `SpaceShower:dipoleRecoil = on`; `PDF:lepton = off`; `TimeShower:QEDshowerByL = off`; `Random:setSeed/seed`.
- **The nucleus is faked**: ⁶Li is *"represented per nucleon by an equal mixture of e+p and e+n runs"* (`gen_dis_hfs.py:26-28`) — two separate free-nucleon runs merged downstream. No Fermi motion, no binding, no shadowing/EMC, no breakup, **no polarization at all** (`gen_dis_hfs.py:2`, `README.md:7,16`).
- **Output**: `np.savez_compressed` in polligen's own `HFSSample` format (`gen_dis_hfs.py:283-298`). Not HepMC, not ROOT.
- **Two silent-cut traps, discovered and fixed** (`docs/reproduction_manual.md:1311-1329`): `PhaseSpace:mHatMin` defaults to 4 GeV and, since `m̂² = x s`, silently removed everything below x = 16/s = 0.004 at 10×99.5 — **39% of the selected rate** (`docs/code_review_2026-08-28.md:49-54`, R2). `PhaseSpace:Q2Min` does nothing unless `Q2Min ≥ pTHatMinDiverge²`; the requested 0.7 GeV² was ignored until `pTHatMinDiverge = 0.5` was set.
- **Standing production**: 6 files, 8 M events, 3.7 GB, git-ignored, at the three γ-matched configurations × {p, n}, e.g. `evgen/samples/pythia8_e10_p99.5_dis.npz` (2 M events, `sigma_gen_mb = 9.47e-4`, seed 101, PYTHIA 8.311).

### 2.2 BeAGLE — `tools/beagle/`

- **Not installed, never run locally.** `tools/beagle/README.md:12-16`: *"**Nothing of the local build is present here**: no `~/Projects/eic/beagle_deps`, no BeAGLE checkout, no FLUKA."* The blocker is the personal FLUKA licence (`README.md:75-93`).
- `build_beagle.sh` (35 lines) checks `$FLUPRO/libflukahp.a`, LHAPDF, CERNLIB (`libmathlib`), RAPGAP-3.302, then invokes `make all` with `-std=legacy` Fortran flags. `env.sh` sets `BEAGLESYS`, `DEPS`, `LHAPATH`.
- **Cards** (`tools/beagle/cards/`): `eD_18x135_test.inp` (`TARPAR 2.0 1.0 0`, `FERMI 2 0.62 1 0`, `MOMENTUM 18.0 110.0`, `MODEL PYTHIA`, `START 1E6 0`); `eLi6_18x137_draft.inp` (`TARPAR 6.0 3.0 0`, `MOMENTUM 18.0 137.5`, `START 1E4 0`) and `eLi7_18x117_draft.inp` (`TARPAR 7.0 3.0 0`, `MOMENTUM 18.0 117.9`). Both Li cards carry `* FERMI  <-- intentionally left at BeAGLE defaults for A>4 (TODO)` at line 15 and a TODO block at lines 5-10 ("confirm the appropriate A>4 mode/defaults with the BeAGLE authors"; L-TAG ranges; USERSET max E*).
- **What BeAGLE is actually used for**: streaming *official* e+d samples over xrootd for a control study that calibrates the cluster model's p_T tail. `tools/analysis/dump_spectators.py` reads a HepMC3 `.tree.root` with `pyHepMC3.rootIO.ReaderRootTree` and writes CSV `ievt,kind,pdg,px,py,pz,e` (kind ∈ B/E/S). `tools/analysis/ed_control_analysis.py` parses that CSV, selects spectator protons, computes pT/xL/R/θ, overlays the Hulthén model and far-forward routing, and with `--beta-scan` fits β.
- **The control study's verdict** (`docs/reproduction_manual.md:1373-1382`):
  ```
  P(p_T >)                       0.10     0.20     0.30     0.45 GeV
  BeAGLE e+d (x_L ∈ [0.9,1.1))  0.1882   0.0494   0.0261   0.0144
  Hulthén β = 0.26 (default)    0.1522   0.0201   0.0037   0.0005
  Hulthén β = 0.40 (best fit)   0.1934   0.0360   0.0089   0.0015
  ```
  *"no β reproduces the shape — BeAGLE has a narrower core and a harder tail than a two-parameter Hulthén can have at once … its model uncertainty is one-sided upward."*
- **BeAGLE's own limitation for Li** (`plans/02_phase1_event_generation.md:257-264`): *"A>4 inherits the **C-12 Fermi-momentum parameterization** … geometry is Woods–Saxon with **no α+d / α+t cluster structure** (fragments come from FLUKA statistical de-excitation, not cluster knockout); code frozen since v1.03.02 (2023)."* There is no A = 6 or 7 BeAGLE sample anywhere (`docs/reproduction_manual.md:1356-1357`).

### 2.3 `tools/fullsim` — ePIC

- `ion_gun_hepmc.py` (112 lines) writes **HepMC3 Asciiv3** because *"npsim cannot shoot a nucleus"* — `--gun.particle Li6` returns "Bad particle type" (`:4-7`). The writer (`:58-83`) emits `HepMC::Version 3.02.05` / `E <i> 0 1` / `U GEV MM` / `P 1 0 <pdg> …` with the vertex deliberately omitted and the vertex count set to 0.
- `ff_gun_hits.py` reads `edm4hep.root` with uproot; `ff_gun_scan.sh` drives npsim directly for single-fragment scans. Hit level only; no EICrecon step, no beam divergence, no vertex spread (`tools/fullsim/README.md:518-525`).

### 2.4 What is missing / hacked in the tooling

- No polarization anywhere in `tools/`.
- No real nuclear target in the PYTHIA path (p/n mixture).
- BeAGLE unbuilt; A = 6,7 breakup entirely gated on FLUKA.
- BeAGLE Li cards are drafts with `FERMI` unset and TODOs unresolved.
- `pyHepMC3.rootIO.ReaderRootTree` segfaults in the current container; the legacy container is required for HepMC3 reading (`tools/fullsim/README.md:22-26`).
- **polligen itself has no HepMC3 writer** (`plans/03_phase2_full_simulation.md:231-239`: *"the generator has no HepMC3 writer (plans/05 step 5.D)"*; `eicrecon` is invoked in zero scripts).

---

## 3. Requirements a NEW dedicated generator must satisfy

Compiled from `plans/00`, `plans/02`–`plans/05`, `plans/08`, `plans/04`, `docs/code_review_2026-08-28.md`, `docs/reproduction_manual.md`.

### 3.1 Polarized beams and spin states

`plans/05:41-52` is the requirement table — five observable rows, all from one master cross section and one run-plan bookkeeper:

| observable | e beam | ion spin | axis | terms |
|---|---|---|---|---|
| A∥ → g₁ᴬ | λ_e = ±, P_e | vector P_z | longitudinal | F₁,F₂ + λ_e P_z (g₁,g₂) |
| A⊥ (g_T) | λ_e = ± | vector P_T | transverse | λ_e P_T cos(φ−φ_S), γ-suppressed |
| A_zz → b₁ | unpol | tensor P_zz | longitudinal (and ‖q) | (2/3)a_m (b₁,b₂) |
| A_cos2φ → Δ | unpol | tensor P_zz | **transverse** | c_m sin²θ_m Δ cos 2φ |
| tagged versions of all three | as above | as above | as above | tagged SFs vs (α_s, p_sT); spin ⊗ cluster wave function |

Plus (`plans/05:54-59`): *"bunch-by-bunch λ_e flips; ion fills with population patterns (p₊, p₀, p₋) [⁶Li] / (p₃/₂ … p₋₃/₂) [⁷Li]; relative-luminosity offsets between spin states at the 10⁻⁴ level; per-fill P_e, P_z, P_zz values with polarimetry uncertainty (δP/P ≈ 3%) … Every event carries its spin labels."*

- ρ(m,m′) for J = 1 **and** J = 3/2 from (P_z, P_zz, n̂(θ_S,φ_S)), rotated analytically to lab/photon frame, exposing vector/tensor/(rank-3) moments (`plans/05:129-134`).
- Spin-3/2 positivity: vector fills with zero tensor/octupole violate positivity above P_z ≈ 0.55 → spin-temperature populations as the default fill model (`plans/00:1028-1030`).
- φ-density positivity guard with the **exact** minimum over φ, not the `1+|A|+|B|` envelope, since a negative density silently samples `max(W,0)` (`plans/08:244-252`).
- Relative-luminosity tolerances differ by observable: 10⁻⁴ for the A_zz thirds estimator, 10⁻³ for the two-state cos 2φ ratio (`plans/05:345`).
- Transverse ion polarization at IP6 for Li is undefined machine-side (`plans/04:67-75`); clean b₁ prefers polarization **along q** (Cosyn arXiv:2410.12764) — a spin-direction systematic to design for.
- No HepMC3 convention exists for ion spin (`plans/04:208-210`, #17) — define one and propose it upstream.

### 3.2 Spin-dependent cross sections on spin-1 / spin-3/2

`plans/05:135-146`: spin-1 must carry the **HJM set {F₁, F₂, b₁…b₄, Δ, g₁, g₂} in the covariant Cosyn–Weiss classification**, *"whose inclusive limit fixes every sign and y-factor, including the transverse-vector λ_e P_T cos(φ−φ_S) g_T term the fastsim never needed."* Spin-3/2: rank-0/1 exact with P_p = 0.866, P_n = −0.037; rank-2 as scenario inputs; rank-3 dropped.

Open/blocking items:
- **b₁ sign** (`plans/08:373`, D1): code is opposite to Cosyn Eq. (27)/HERMES; what flips is κ and any O(γ²) subtraction built on it.
- **Exact finite-γ kernel** (`plans/08:381`, D2): Cosyn Eqs. 9/10/14/24, b₃/b₄ slots, Eqs. 17d/17e. Gated on D1. The currently quoted leakage γ²b₁/6 is *"≈ 7× low; the exact ratio a₂(full)/a₂(17e alone) ≈ 6.9"*.
- **Complete spin-3/2 basis does not exist in the literature** (`plans/04:192-200`, #14): *"rank-2 (b₁-analog) and rank-3 functions are not classified anywhere we can adopt; needed by the doubly polarized generator."* The rank-2 *kernel geometry* is already isotope-generic; what is blocked is the structure function.
- **No b₁ prediction for any A > 2; no EIC Δ projection for any target** (`plans/04:179-182`, #9).
- **⁶Li effective polarization value unresolved** (`plans/04:125-158`, #6): 1/3 per-nucleon (Cloët) vs cluster-model 0.81 whole-nucleus — a factor 1.23 in the g₁ FOM.
- **R = σ_L/σ_T**: switching from the toy to R1998 moves Δ/F₁ by +16.6/+18.0/+4.7/−4.4% (`plans/08:340-352`, C2); one `r_func` hook must reach **all four** consumers (F1 from F2, F_L, D(y), and g1/F1).
- Consistency gates: `ρ → P_z only` must reproduce `asymmetries.a_parallel`; tensor-only ρ must reproduce `azz` and `a_cos2phi` bin-by-bin (`plans/05:145-146`).

### 3.3 Cluster structure, Fermi motion, spectator tagging

- Two-cluster light-front IA ⊗ spin, transplanting Cosyn–Weiss (`plans/05:147-183`): ⁶Li(1⁺) = α⊗d at L = 0 (+small L = 2), with the struck "nucleon" replaced by the polarized deuteron carrying its own g₁/b₁ — *"this is precisely the embedded-b₁ observable"*; ⁷Li(3/2⁻) = α⊗t at L = 1, whose CG structure delivers **with no new parameters** both the per-m triton polarization and the m_L-dependent angular distribution |Y₁^{m_L}(k̂)|² of the tagged α.
- **Radial inputs must move to VMC two-cluster overlaps** (Wiringa, ANL) — `plans/04:201-204` (#15), `plans/05:166-172`. The two-parameter Hulthén tail is the dominant systematic and is demonstrably wrong against BeAGLE (§2.2 table).
- The β = 0.20/0.30/0.40 band rule *is not currently kept* — published α-tag numbers run β = 0.30 alone while the scan spans a factor 3.3, and the band needs restating as **one-sided upward** (`plans/05:173-179`).
- FSI not modelled; quote at low spectator virtuality / small |t′| (`plans/05:180-183`, `plans/04:205-207` #16).
- Physical nuclear masses, not A·M_U — the 2×10⁻³ mass error moved the ⁶Li α Roman-Pot fraction by 15% relative (`plans/08:324-338`, C1).
- t* remnant treatment (d vs nn) is *"too crude for double-tag studies"*; double-tag claims gated on a ³He control (`plans/05:342`).
- Off-shellness / light-cone α_s dependence of struck-cluster SFs is explicitly beyond the current tier (`tagged.py:33-34`).

### 3.4 Coherent / intact-recoil channel

- Process: `e + ⁶Li → e' + X + ⁶Li(g.s.)`, A/Z = 2 so the recoil is **exactly at beam rigidity** — the only handle is the Roman-Pot near-beam p_T tail (`coherent.py:1-11`).
- Everything is a scenario with explicit bands (`coherent.py:105-109`): `f0 = 0.04` (×2÷2), `x_coh = 0.01`, `slope_b = 50 GeV⁻²` (band 40–60), flat gluon-transversity `amp = 0.01` (band 3e-3…1e-2), deformation `eps_b0 = −0.08` (band −0.04…−0.13). *(Annotation 2026-09-04, and it applies to LiPolGen's inherited copy too — this line is a record of upstream and is not rewritten: `eps_b0 = −0.08` implies a ⁶Li charge quadrupole of −0.9345 fm², **11.4× the measured −0.0818**, and −(0.04…0.13) is a **deuteron** band. The ⁶Li band is −(0.0070…0.0527). `docs/open_items/run_2026-09-03/phase_C_numbers.md` §C4.)*
- `f_coh(x) = f0/(1+(x/x_coh)²)`, `dσ/dt ∝ exp(−B|t|)`, `acc = exp(−B p_T,cut²)`, `⟨|t|⟩_tagged = p_T,cut² + 1/B` (`coherent.py:111-130`).
- **No calculation exists for any polarized A > 2 nucleus** (`plans/04:211-230`, #18); ⁷Li has **no scenario at all** because ε_B0 cannot be linearly rescaled (it gives −2.2 to −4.5, |ΔB₀| > B).
- Blocks the coherent electron side and x_P/M_X binning (`plans/08:386`, D8); pairing the inclusive map with log-uniform x_P gives β = x/x_P > 1 for 58% of events.
- Published |t| binning: seven bins over 0.017–0.25 GeV² (`recopseudo.T_EDGES_PUBLISHED`).

### 3.5 Radiative corrections

- Mode R design (`plans/05:104-112`): not event-level at first — multiplicative RC bands per (x,y) attached as weights, from DJANGOH 4.6.22 on an effective polarized nucleon (the arXiv:2406.05591 ³He workaround).
- The **unpolarized collinear-ISR half is now done in-repo** (`polligen/radiative.py`), and two of four planned deliverables were void by algebra: collinear ISR fakes exactly zero cos φ′/cos 2φ′ (covariant azimuth invariant under k → (1−z)k for a massless target), and the spin-state ratio cancels anything common to the fills (`radiative.py:9-21`).
- Measured bound: +0.62/+0.50/+0.94/+1.22% of Δ̂ at the four sweet spots, ≤ 2.9% with low-Q² feed-in, against a 5% gate (`plans/08:382`, D3). Requires an **eight-seed average**: a single response seed is not publishable.
- **The tensor-sector RC has never been calculated by anybody** — default: *no band* (`plans/04:183-191`, #10; `docs/reproduction_manual.md:44`).
- Not modelled: polarized-lepton RC, non-collinear real emission, FSR, elastic/quasi-elastic radiative tails (removed by W² ≥ 10 GeV²) — `radiative.py:80-93`.

### 3.6 Event record / output format

- Target: **spin-labeled HepMC3** via `io_hepmc.py` (`plans/05:66-86, 282-294`): *"ASCII HepMC3; status-4 beams so `abconv` accepts it; 10-digit fragment PDG; spin labels as named attributes — no HepMC3 convention exists for ion spin states, so define one and propose it to the ePIC MC group."* Fragments as final-state particles with codes like α = 1000020040, t = 1000010030; excited-ion codes 10LZZZAAAI with I≠0 from BeAGLE (`plans/03:82-89`).
- **Three final-state fidelity tiers** (`plans/05:96-103`): T0 = (e′, spectator, X as one pseudo-particle) — enough for every Phase-1 FOM; T1 = + struck-cluster internal nucleon and partner spectators; T2 = + PYTHIA hadronization.
- Weighted mode (a vector of weights over all spin states per event) must be supported alongside unweighted, for FOM-efficiency studies; fixed RNG streams per (run, bunch) (`plans/05:184-190`).
- Practical HepMC3 gotchas already learned (`docs/reproduction_manual.md:1648-1649`): DD4hep routes `.hepmc` to `ReaderAscii`, so write Asciiv3 not HepMC2 IO_GenEvent; a vertex at the origin with nothing incoming must be **omitted entirely**.
- Chain gate: `HepMC3 → abconv → npsim [→ EICrecon]`, event-by-event 4-momentum and charge conservation, 100-event smoke (`plans/05:333`). **Never run** — Step 5.D is still ☐.

### 3.7 PDF / structure-function inputs demanded

- Unpolarized: CT18NNLO + **EPPS21 / nNNPDF3.0** nuclear ratios for A = 6, 7 (`plans/02:93-96`). Installed and used today: `CT18NLO`, `EPPS21nlo_CT18Anlo_Li6`, `nNNPDF30_nlo_as_0118_A6_Z3`, `NNPDFpol11_100` — via the pure-python `parton` package, **not LHAPDF** (`docs/reproduction_manual.md:84-100`).
- Polarized: JAM/DSSV14 or NNPDFpol1.1 for g₁p/g₁n (`plans/02:98-104`). Note the deliberate flavour asymmetry: F₂ over 5 flavours (d u s c b), g₁ over 3 (d u s), because NNPDFpol1.1 sets Δc = Δb = 0 (`structure.py:231-243`).
- **Polarized nuclear PDFs do not exist** (`plans/04:166-178`, #8) → effective-polarization convolution + CBT/TMT medium curves as scenarios.
- The unpolarized EMC baseline the whole polarized-EMC construction rests on is *"still the hand-written 12-point table awaiting EPPS21"* (`plans/04:174-176`).
- No tensor PDFs exist; b₁ is handled by two digitized theory curves, Δ by the sum-rule-constrained moment ansatz.

### 3.8 Kinematic coverage, beam energies, luminosity

- Scope (`plans/05:61-64`): *"inclusive DIS 10⁻⁴ < x < 1, Q² > 1 GeV², the three reference energies of `beams.default_configs` per isotope; spectator channels of `spectator.py` (⁶Li → α+d\*, d+α\*; ⁷Li → α+t\*, t+α\*), plus evaporation-n ZDC tags left to BeAGLE."*
- **Ions are γ-matched to the proton configuration, not rigidity-scaled** (`beams.py:1-31`), because HSR and ESR must share a revolution period. Verified against Yellow Report Table 10.2 (gold at 41 GeV/u, γ = 44.02). Confirmed by running `beams.default_configs`:

| | 5 × 41 | 10 × 100 | 18 × 275 |
|---|---|---|---|
| ⁶Li | 5.0 × **40.8** GeV/u (√s = 28.57) | 10.0 × **99.5** (63.09) | 18.0 × **137.5** (99.50, rigidity-limited) |
| ⁷Li | 5.0 × **40.8** (28.57) | 10.0 × **99.5** (63.09) | 18.0 × **117.9** (92.13) |

  `plans/00:21`: *"⁶Li sits at **40.8 / 99.5 / 137.5 GeV/u**, not 20.5 / 50 / 137.5."* The accessible menu is an isolated γ ≈ 43.7 point plus a band from γ ≈ 106.6 to the rigidity cap, *"with NOTHING in between"* (`beams.py:25-27`).
- Analysis cuts (`fom.Scenario`, `fastsim/polli_fastsim/fom.py:56-73`): `q2_min = 1.0`, `y_min = 0.01`, `y_max = 0.95`, `w2_min = 10.0`, `eta ∈ [−3.5, 3.5]`, `E'_min = 0.5 GeV`, `P_e = 0.70`, `P_z = 0.70`, `P_zz = 0.60`.
- Generator window must be **looser** than the analysis window so events can migrate in (`recopseudo.generator_scenario`, `:110-117`): `q2_min = 0.7`, `y ∈ [0.004, 0.985]`, `w2_min = 8.0`, `eta_pad = 0.3`, `E'_min = 0.3`.
- **Luminosity**: `lumi_fb_per_nucleon = 10.0` default, quoted ∈ {1, 10, 100}. *"**Li luminosity.** Confirmed gap — no Li number exists in any document (EPIOS included)"* (`plans/04:114-117`, #4). A `run_share` / `--lumi-fraction` knob prices run-plan splits; errors scale as 1/√f.
- The coherent tagging optics costs L/L_HA = 1/7.1 – 1/13.3 (`docs/reproduction_manual.md:1121-1123`).

### 3.9 Downstream consumers

The measurability audit (`docs/code_review_2026-08-28.md:131-143`) lists exactly what the analysis chain consumes: spin state + P_zz + L_f per state; the smeared scattered electron in the lab transformed to head-on → Q²_e, θ_e, φ′; Σ_h and p_T,h → y_Σ, x_mixed; the calibration factor from `cell_means_reco(x_mixed^uncal, Q²_e)`; reconstructed-quantity cuts and bin assignment; φ′ from the covariant formula on (k, k′_meas, P, S); purity/efficiency/D/K/fold (MC-derived, never fitted); the Roman-Pot angle pair → β, |t|; the β basis and t_ref; u₁, u₂ and luminosity shares (assumed, as fit arguments); truth references a_t(t_ref), a_e (plotted only).

---

## 4. Data structures and file formats a drop-in generator must produce

There is **exactly one persisted event format** in the repository. Everything else is in-memory dicts of NumPy arrays; the scripts write PNG figures, not data files (`grep` for `np.savez|to_csv|savetxt` over `evgen/` returns only `hfs.py:102`).

### 4.1 `HFSSample` — the on-disk `.npz` schema (`evgen/polligen/hfs.py:68-116`)

CSR/flat-array layout. This is what `tools/pythia8/gen_dis_hfs.py` writes and what `hfs_resolution.py`, `hfs_acceptance.py`, `money_cos2phi_reco.py --hfs-sample` read.

| key | shape | dtype | meaning |
|---|---|---|---|
| `offsets` | (n_ev+1,) | int64 | particle-index offsets; event i owns `[offsets[i], offsets[i+1])` |
| `pid` | (N,) | int64 | PDG id, **neutrinos kept** (dropped later by `HadronResponse`) |
| `charge` | (N,) | float64 | particle charge |
| `p4` | (N,4) | float64 | (E, px, py, pz) GeV, head-on frame, ion +z, electron −z |
| `x`, `q2`, `y` | (n_ev,) | float64 | per-event, **per-nucleon** truth kinematics |
| `kp` | (n_ev,4) | float64 | scattered-electron four-vector |
| `weight` | (n_ev,) | float64 | per-event weight (1.0 for unweighted) |
| `e_energy` | () | float64 | electron beam energy [GeV] |
| `p_per_nucleon` | () | float64 | ion momentum per nucleon [GeV] |
| `meta` | () | `<U…` | JSON string: generator, version, settings list, target, `sigma_gen_mb`, `sigma_err_mb`, n_events, n_tried, seconds, frame |

Notes for a replacement generator: `s = 4·e_energy·p_per_nucleon` (`hfs.py:92-94`); the scattered electron is **excluded** from the particle list; `meta["sigma_gen_mb"]` is required for cross-section-weighted p/n merging in `concatenate` (`hfs.py:119-164`); `HFSResponse.check_beams` refuses a library from different beam energies (`hfs.py:849-861`). Naming convention: `pythia8_e{E_e}_{p|n}{p_u}_dis.npz`.

Exact identities a new sample must satisfy (`hfs.truth_kinematics_check`, `:873-886`): `Σ(E − p_z) = 2 E_e y + m²/(E_N + p_N)` and `|Σp_T,h| = p_T,e`. The massless form holds exactly for the toy; PYTHIA satisfies the massive one to a few 10⁻³, and the p_T identity to ≤ 5×10⁻¹¹.

### 4.2 Inclusive event dict — `InclusiveSampler.sample_category` (`sample.py:187-224, 289-291`)

```python
{"x": f8[n], "q2": f8[n], "y": f8[n], "phi": f8[n],
 "m": f8[n],          # ion spin projection of the event
 "cell": i8[n],       # flat index into the accepted (x,Q2) cell grid
 "category": str,     # SpinCategory name
 "lam_e": int}        # electron helicity of the category
```
`phi` is the DIS azimuth in [0, 2π). Consumed by `estimators.binned_counts`, `closure_fom.py:38,65-67`, `reco_chain_figures.py:96-97`.

### 4.3 Tagged event dict — `TaggedSampler.sample_category` (`tagged.py:551-569`)

The inclusive dict minus `m`, plus:
```python
"m_ion", "m_struck",                       # spin labels
"k", "cos_theta_k", "phi_k",               # spectator in the SPIN frame
"pT", "theta", "p_lab", "R", "xL",         # spectator in the LAB
"kx", "ky", "kz", "phi_spec",              # lab components + lab azimuth
"route"                                    # far-forward routing code
```
`R` = rigidity ratio to the beam; `route` from `farforward.route_charged`; `rp_accepted(events)` = `(route == 1) | (route == 4)` (`tagged.py:572-574`). Consumed by `money_tagged_azz.py:160-170` (`ev["R"], ev["theta"], ev["pT"], ev["phi_spec"], ev["k"]`) and `tagged_polarimetry_7li.py:125-126,171,196` (`ev["cos_theta_k"], ev["x"]`).

### 4.4 `HadronResponse.reconstruct_sums` (`hfs.py:401-422`)

`{"sigma", "ptx", "pty", "sigma_true", "ptx_true", "pty_true"}` — all (n_ev,).

### 4.5 `hadronic_kinematics` (`hfs.py:190-220`)

`{"q2_e", "y_jb", "q2_jb", "x_jb", "y_sigma", "q2_sigma", "x_sigma", "y_da", "q2_da", "x_da", "x_mixed", "x_ejb", "x_eda"}`. The chain's production choice is `x_mixed = q2_e/(s·y_sigma)`.

### 4.6 `RecoResponse` attributes (`recopseudo.py:296-315`)

Flat arrays over pseudo-events, all length N: `cell`, `w` (pb weight), `x`, `q2`, `y` (hard y), `y_nominal`, `x_reco`, `q2_reco`, `y_reco`, `eta_reco`, `e_prime_reco`, `theta_reco`, `eff`, `dil` (cos 2Δφ′ dilution), `phi_true`, `isr_z`, `isr_dphi`, `sigma_capture`.

### 4.7 Bookkeeping

`SpinCategory(name, j, populations, lam_e, pe, theta_s, phi_s, lumi_fraction)` (`bookkeeping.py:36-49`) and `RunPlan(categories, pe_true, pz_true, pzz_true, delta_p_over_p, polarimetry_seed, measured)` (`:52-77`). **Populations are ordered m = +J … −J everywhere.**

---

## 5. Tests, and what "closure" would mean for a new generator

`evgen/tests/` — 17 files, 305 tests collected (`python -m pytest evgen/tests --collect-only -q`). `fastsim/tests/` adds 57. Named validation gates are in `plans/05:321-333`.

### 5.1 The gate matrix and its current status

| gate | reference | criterion | status |
|---|---|---|---|
| ρ moments, all axes | analytic Wigner rotations | exact, atol 1e-12 | ☑ `test_spin.py` (16 tests) |
| master formula ↔ `asymmetries.py` | code-to-code, both backends | bin-wise identity, rtol **1e-12** | ☑ `test_xsec_identity.py` (15) |
| pseudo-experiment estimators | `fom.py` maps | pulls unbiased; spreads within 15% of `err_a_parallel`/`err_azz`/`err_cos2phi_amplitude` | ☑ `test_pseudoexp.py` + `scripts/closure_fom.py` (~65 x-bins/isotope) |
| φ-modulation recovery | injected Δ | amplitude unbiased with uniform **and holey** φ acceptance | ☑ `test_cos2phi_fit_unbiased_with_holey_acceptance` |
| deuteron limit of tagged mode | Cosyn–Weiss arXiv:2603.23700 Eq. (6.12)–(6.14), TABLE II | analytic match | ☑ `test_cosyn_weiss_tensor_gate` |
| unpolarized spectator spectra | official BeAGLE e+d | bulk agreement, tail as model band | ☑ but **no β reproduces the tail** |
| forward limit of tagged ⁷Li | P_p = 0.866, P_n = −0.037 | within D-state band | **☐ only the proton half is tested**; model gives P_n = −0.028 vs the gate's −0.037, and no D-state band is defined for the neutron |
| ⁶Li embedded-d b₁ scaling | `b1_li6_from_deuteron` | 0.921947 × 2/6 | ☑ `test_li6_b1_rank2_transfer_constant_is_pinned_to_the_model` |
| conservation & chain | HepMC3 → abconv → npsim | 4-momentum/charge per event; 100-event smoke | **☐ never run** |

### 5.2 Genuine physics-closure tests (not code-to-code)

These are what a replacement generator must reproduce.

- **`test_tagged.py:730-785` `test_cosyn_weiss_tensor_gate`** — the strongest external gate. Against Cosyn–Weiss II Eq. (6.12): the P₂(cosθ_k) angular factorization is exact (`ratios.max() − ratios.min() < 1e-5`, mean `0.99940 ± 1e-4`); the radial quadratic form peaks at 1.000 at k = 0.3098 GeV against CW's 0.30 GeV for AV18 (`abs=0.02`); `A_T∥ = −2 A_zz^wf` gives **+0.9997** and **−2.000** against CW TABLE II's +1 and −2; the whole curve stays inside CW's stated [−2, 1].
  *(**CORRECTION 2026-09-06 — do NOT reproduce this gate.** It is wrong twice
  over, and LiPolGen's `docs/benchmarking/07_cw_sign_investigation.md` settled
  both. (i) **The mapping is `A_T∥ = +1 · A_zz^wf`, not −2** — CW's −2 is the
  value of their own angular factor at θ_k = 0, a node of the Λ = ±1 densities
  (CW Eq. 6.13), which `A_zz^wf` already carries; applying it again
  double-counts. (ii) The peak *"at 1.000 at k = 0.3098 GeV"* is on the
  **Hulthén** pair, whose f₂/f₀ never reaches √2 anywhere on the grid — it is
  matching CW Eq. (6.14)'s **minimum** at f₂/f₀ = 1/√2 while calling it
  Eq. (6.13)'s **maximum** at √2, and it agrees with CW's "0.30 GeV" only by a
  coincidence of β = 0.30 on a channel that is not AV18. Together those two
  errors hid an **inverted S–D interference sign** in the partial-wave sum
  (`polligen/tagged.py:243-248` `_amp2_table` omits the i^L phase, exactly as
  LiPolGen's `build_amp2` did until 2026-09-06). The three pass values quoted
  above — `0.99940`, `k = 0.3098 GeV`, `+0.9997 / −2.000` — are all pre-fix.
  **What a replacement generator must reproduce instead** is CW Eq. (6.12) as
  an *identity* on the **AV18** wave function CW quote TABLE II for:
  LiPolGen's rewritten gate measures `max|A_zz^wf − CW| = 8.881784e−16` over
  all 280 × 96 = 26 880 cells (it was **2.740499e+00**), CW's own k landmarks
  at 0.298121 and 1.034872 GeV, and TABLE II's three rows at
  −1.937124 / +0.999313 / +0.967340 against −2 / +1 / +1. LiPolGen's gate is
  `tests/test_tagged.cpp` "tagged: the Cosyn-Weiss deuteron tensor gate
  (CW TABLE II)"; the sibling was **not** modified.)*
- **`test_tensor_convention.py:49-69`** — the only literature anchor for the tensor sector: `A_zz(θ_S=0)·(1 + ε(y)R) == TENSOR_LL_SIGN·(2/3)·b₁/F₁` **exactly and at every y**, rel 1e-12, at six (x,Q²) points. Its docstring corrects a widely-quoted form that *"double-counts R … and misses by a factor 1.17"*. `:72-84` asserts `TENSOR_LL_SIGN == +1.0` with a message naming plans/08 D1 as the switch. `:87-107` checks the kernel's thirds combination `(w₊+w₋−2w₀)/(3+Σw)` reproduces `azz` including sign, and `w₊ = ½·azz`. `:110-131` proves the Δ sector is sign-independent. `:133-167` pins the unified rank-2 geometry for both spins.
- **`test_xsec_identity.py`** — sector-by-sector identity at rtol 1e-12: `test_vector_sector_matches_a_parallel` (`:41-53`, `(w₊−w₋)/(2+w₊+w₋) == P_e·a_parallel`), `test_vector_sector_spin32` (`:56`), `test_tensor_sector_matches_azz` (`:72`), `test_transverse_tensor_matches_a_cos2phi` (`:87`), `test_dsigma_reduces_to_unpolarized` (`:103`), `test_g2_ww_analytic_power_law` (`:119`), `test_identities_on_grid_backends` (`:181`, actually running CT18NLO + NNPDFpol11_100), and four R-sensitivity tests (`:273-379`).
- **`test_pseudoexp.py`** — estimator closure with tolerances: `test_apar_closure` and `test_azz_closure` assert the trial mean equals truth to `4·std/√ntrials` (`:88-89`, `:108-109`); `test_azz_relative_lumi_bias` / `test_apar_relative_lumi_bias` assert the **naive** estimator's bias equals the closed forms `−(2/3)δ/P_zz` and `δ/(2 P_e P_z)` to `5·se`, and that the lumi-corrected estimator is unbiased (`:133-134`, `:155-156`); `test_cos2phi_fit_unbiased_with_holey_acceptance` (`:190-218`) removes two asymmetric φ sectors and requires the binned fit unbiased at 5·se while the naive moment fails by >10·se.
- **`test_target_mass.py`** — E143 finite-γ closure: `test_finite_gamma_factors_match_the_lab_frame_definitions` (`:87`) pins ε, D, η to double precision against E143's (E, E′, θ) definitions; `test_target_mass_a_parallel_matches_a_lab_frame_construction` (`:115`); `test_gamma2_capped_by_the_w2_cut_at_every_configuration` (`:58`); `test_target_mass_off_is_bit_for_bit_the_published_kernel` (`:160`); `test_small_y_collapse_to_one_plus_gamma_squared` (`:211`).
- **`test_hfs.py`** — `test_toy_conserves_four_momentum_and_truth_sums` (`:32`); `test_kinematic_methods_are_exact_with_a_perfect_response` (`:48-56`): with a perfect response, `y_jb`, `y_sigma`, `y_da` and `x_mixed` recover truth to **< 5×10⁻³** and `q2_e` to **< 10⁻⁹**; `test_library_transfer_reproduces_the_response_statistics` (`:85`); `test_concatenate_refuses_to_merge_across_beam_energies` (`:163`); `test_reco_response_refuses_a_library_from_other_beams` (`:179`).
- **`test_reco.py`** — `test_covariant_azimuth_matches_explicit_collinear_frame` (`:140`); `test_phi_s_equals_lab_angle_for_massless_target` (`:73`) and `..._massive_target_deviation_is_order_gamma2` (`:87`); `test_covariant_azimuth_is_invariant_under_a_collinear_photon` (`:171`) — pinned at 3.6×10⁻¹⁵ rad over 2×10⁴ events; `test_ratio_estimator_unbiased_under_phi_dependent_efficiency` (`:262`) and `..._immune_to_relative_luminosity_offset` (`:291`); `test_lab_azimuth_about_detector_axis_only_odd_harmonics` (`:49`).
- **`test_radiative.py`** (30 tests) — `test_spectrum_normalisation_is_second_order_in_t` (`:31`); `test_sampler_reproduces_the_closed_form_moments` (`:63`); `test_observed_kinematics_match_a_four_vector_construction` (`:115`); `test_mixed_x_is_exact_and_only_the_q2_label_migrates` (`:128`); `test_isr_off_path_is_the_pre_hook_constructor` (`:210`).
- **`test_tagged.py`** (26 tests) — normalization of every channel (`:49`); `test_pure_s_wave_is_m_independent` (`:64`); `test_li7_p2_moments_and_polarimeter` (`:145`, the −T/5 relation); `test_li7_triton_polarization_forward_limit` (`:159`, P_p within 0.02 of 0.866); `test_li6_s_wave_reduces_to_inclusive_deuteron` (`:200`); `test_boost_matches_fastsim_spectator` (`:325`, agreement quantile-by-quantile at P_D = 0); `test_tagged_rate_asymmetry_matches_analytic` (`:343`).
- **`test_coherent.py`** — `test_sample_t_matches_slope_and_acceptance` (`:41-48`, ⟨t⟩ = 1/B to 2%, sampled acceptance = `exp(−B cut²)` to 2%); `test_veto_table_routing` (`:141-167`, α R = 0.99813, d R = 1.00452, ³He 0.75204, t 1.50437, all to 5×10⁻⁶); `test_mantysaari_table_m_state_relation` (`:109`, a₀ = −2a₁ to 25%).
- **`test_sampler.py`** — `test_mode_w_weights_reproduce_rate_ratio` (`:121`, the reweighting closure); `test_density_min_is_exact` (`:146`); `test_sampler_refuses_a_negative_phi_density` (`:159`); `test_sampled_events_respect_acceptance` (`:134`).
- **`test_recopseudo.py`** (largest, 993 lines) — reco-level closure: `test_expected_counts_match_generator_level_when_unsmeared` (`:94`); `test_ratio_fit_closure_on_reco_bin` (`:114`); `test_delta_bin_centering_exact_without_noise` (`:135`); `test_amplitude_is_exactly_linear_in_delta` (`:547`); `test_folded_fit_recovers_the_truth_from_a_wrong_prior` (`:711`); `test_pzz_scale_error_propagates_one_to_one_not_quadratically` (`:877`); `test_two_azimuth_fit_closure` / `_pulls` (`:194`, `:225`).
- **`test_likelihood_fit.py`** (16 tests) — `test_likelihood_is_exact_on_asimov_counts` (`:78`, Asimov closure 4×10⁻¹⁶); `test_likelihood_is_unbiased_where_the_ratio_fails` (`:135`).
- **`test_review_20260828.py`** (16 tests, the regression suite of the last code review) — `test_calibration_never_reads_the_true_cell` (`:82`, a monkeypatch that *forbids* the truth lookup); `test_hadron_acceptance_is_applied_in_the_lab_frame` (`:161`); `test_target_mass_term_enters_the_transferred_sigma` (`:129`); `test_beam_energy_spread_has_the_analysis_sign` (`:188`); `test_sigma_weighted_merge_of_p_and_n` (`:326`); `test_hadronic_methods_reproduce_the_exact_kinematics_by_hand` (`:248`, to 10⁻⁹); `test_selected_bins_of_the_r_script_use_measurable_criteria_only` (`:150`).
- **`test_bookkeeping.py`** — spin-temperature P_zz relations, `pzz_true` evenness in P_z, refusal of unsupported spins.

### 5.3 What a new generator would need for "closure"

A drop-in replacement would have to pass, in this order:

1. **Reproduce the analytic asymmetry sector-by-sector.** `A_par`, `A_zz`, `A_cos2φ` from the generated event rates must equal `polli_fastsim.asymmetries` bin-by-bin at rtol 1e-12, on both toy and PDF-grid backends, for J = 1 and J = 3/2, at arbitrary (θ_S, φ_S). This is the existing gate and is not negotiable — it is the only thing that pins signs, y-factors and the depolarization convention.
2. **Reproduce ρ moments exactly** for arbitrary axes and populations, J = 1 and 3/2, including the max-entropy fill model and its positivity limits.
3. **Estimator closure at pseudo-experiment scale**: pulls unbiased and spreads within 15% of the three analytic error formulas across ~65 x-bins per isotope, plus the two closed-form relative-luminosity biases and their exact removal.
4. **Azimuthal-amplitude recovery with holey acceptance**, since a naive moment estimator demonstrably fails there.
5. **The Cosyn–Weiss deuteron limit of the tagged mode**, quantitatively: the P₂(cosθ_k) factorization, the k = 0.30 GeV peak of the quadratic form, and A_T∥ ∈ [−2, +1] with the two TABLE II extrema. *(**Restated 2026-09-06** — as written this asks a replacement to reproduce a gate that is wrong; see the correction under §5.2. The requirement is CW Eq. (6.12) as an identity on **AV18** with the mapping `A_T∥ = +1 · A_zz^wf`, which LiPolGen now meets at 8.88e−16 over 26 880 cells.)*
6. **The unpolarized spectator spectrum against BeAGLE e+d** — bulk routing to better than 2 points, with the p_T tail carried as a one-sided upward band (a new generator with VMC overlaps would need to *close* this, not just bracket it).
7. **⁷Li forward limits**: P_p = 0.866 **and P_n = −0.037** — the neutron half is currently untested and the model disagrees (−0.028), so this is an open gate a new generator should actually pass.
8. **HFS truth identities** on any sample it produces: `Σ = 2E_e y + m²/(E_N+p_N)` and `|Σp_T,h| = p_T,e`; and the perfect-response kinematic methods recovering truth to 5×10⁻³ (Q²_e to 10⁻⁹).
9. **Frame/azimuth invariants**: covariant φ′ = lab azimuth for a massless target, O(γ²) with the physical ion mass; collinear-ISR invariance to ~10⁻¹⁵ rad; crossing-angle affecting only odd harmonics.
10. **The chain gate that has never run**: HepMC3 → `abconv` → `npsim` (→ EICrecon), event-by-event 4-momentum and charge conservation, 100-event smoke.
11. **No analysis-side truth leaks** — the `test_review_20260828.py` discipline: every correction keyed on reconstructed variables only, enforced by tests that forbid the truth lookup.
12. **Bit-for-bit reproduction of the published numbers** on the default paths, which the current suite enforces with several `..._is_bit_for_bit_the_published_kernel` tests; a new generator would need an equivalent regression anchor or an explicit, documented break.

---

## 6. Bottom line for the "new generator" decision

- The **physics kernel already exists and is validated to rtol 1e-12** against independent analytic formulas and one external literature relation each for the tensor sign and the tagged deuteron limit. A C++ rewrite inherits a working specification, not a blank page.
- What is genuinely missing on the **engineering** side is narrow and enumerated: Mode W (BeAGLE reweighting, step 5.C), the T1 final-state tier and `io_hepmc.py` (step 5.D), and the abconv→npsim→EICrecon smoke test. `plans/05:315-319`: *"what is left of the 7–9 weeks is 5.C + 5.D + 5.E."*
- What is missing on the **physics-input** side is not solvable by writing code: the spin-3/2 rank-2/3 structure-function basis is unclassified in the literature; no b₁ exists for A > 2; no coherent diffractive amplitude exists for any polarized A > 2 and none at all for ⁷Li; polarized nuclear PDFs do not exist; the tensor-sector radiative correction has never been calculated; VMC α+d/α+t overlaps are not yet in hand and the Hulthén stand-in is demonstrably wrong in the tail; and HepMC3 has no ion-spin convention to adopt.
- Two conventions remain **author decisions**, and both propagate: the b₁-sector sign (`TENSOR_LL_SIGN`, one constant, but it flips κ and any O(γ²) subtraction built on it) and the ⁶Li effective polarization value (1/3 vs 0.81, a factor 1.23 in the g₁ FOM). *(Both are carried forward for confirmation in `docs/open_items/run_2026-09-03/AUTHOR_DECISIONS.md` §B1 and §B6. In LiPolGen the sign has been **taken** — `constants.hpp` says "author decision, plans/08 D1", −1 since 2026-08-29 — while `OPEN_ITEMS_SOLUTIONS.md` row 3 and this survey still say "to confirm"; that difference is B1's whole question. The `+1.0` quoted above is the sibling `fastsim` tree's, not this one's.)*
- Forking BeAGLE would inherit the wrong nuclear model for Li (C-12 Fermi parameterization, Woods–Saxon, no α+d/α+t clustering, FLUKA-gated, code frozen since 2023). Forking PYTHIA 8 gives hadronization and a working build path — the project already builds its bindings natively in 2.5 minutes — but no nucleus, no spin, and no spectator. The repository's own conclusion (`plans/05:25-30`) is to reinvent only the polarized-nucleus vertex and the spin-correlated cluster-spectator sampler, and to delegate hadronization to PYTHIA, RC to DJANGOH/HERACLES, and evaporation/backgrounds to BeAGLE.

---

## Addendum to §5 — test-suite detail (verified)

A deeper pass over all 17 test files returned after my report. I spot-checked the load-bearing new numbers against the source; all confirmed. The additions below are new material only — nothing in §1–§4 changes.

### Corrections/additions to the gate picture

The suite defines expected behaviour at **three distinct tiers**, and this matters for a replacement generator because only the first is non-negotiable:

**Tier 1 — external physics anchors** (must reproduce regardless of implementation). Six independent published sources are pinned in code:

| anchor | file:line | pinned value |
|---|---|---|
| Cosyn et al. EPJ A 61 (2025) 83 Eq. (27) | `test_tensor_convention.py:61-69` | `A_zz·(1+ε(y)R) == ±(2/3)b₁/F₁`, rel 1e-12 |
| Cosyn–Weiss arXiv:2603.23700 Eq. 6.12/6.13, TABLE II | `test_tagged.py:730-785` | `A_T∥` extremes **+0.9997 / −1.9378** vs published +1 / −2. ***Superseded 2026-09-06** — both the −2 mapping and the sibling's missing i^L are wrong (see §5.2's correction). LiPolGen's replacement checks Eq. (6.12) as an identity on AV18: 8.881784e−16 over 26 880 cells, was 2.740499e+00.* |
| **Hoodbhoy–Jaffe–Manohar NPB 312:571 (1989) Eq. (30)** | `test_recopseudo.py:545-577` | `A = −[(1−y)/y²]⟨c_eff⟩sin²θ_S·Δ/D_φ`, written out *"from the paper rather than taken from the code"*, rtol 1e-12 |
| E143 PRD 58:112003 lab-frame ε/D/η | `test_target_mass.py:87-112` | rebuilt from random (E, E′, θ), rel 1e-12 over 64 draws |
| Kuraev–Fadin / Skrzypek LL spectrum | `test_radiative.py:40-46` | `t = (2α/π)(L−1)` with α = 1/137.035999; `S(t) = 1+3t/8+O(t²)` |
| Yellow Report Table 10.1 divergences | `test_reco.py:327-338` | proton σ_θ(h,v) = **[220, 380] / [180, 180] / [65, 65]** µrad, high-divergence [150, 150] |
| AME2020 ⁷Li breakup thresholds | `test_nearbeam.py` | 2.468 / 7.251 / 8.725 / 9.975 MeV, abs 0.002 |
| Textbook CG values | `test_spin.py:43-46` | `⟨1,1;1,−1\|0,0⟩ = 1/√3`, `⟨1,0;1,0\|2,0⟩ = √(2/3)` |

**Tier 2 — internal identities at rtol 1e-9…1e-14**, each rebuilt independently from raw kinematics rather than compared code-to-code. Beyond what §5.2 listed:

- `test_xsec_identity.py:297-321` — hand-derived closed form `a₂ = −(ε/2)·c_eff·sin²θ_S·(Δ/F₁)/(1+ε·R)` at rtol 1e-12, proving **F₂ cancels exactly** from the cos 2φ amplitude; `:324` then shows scaling F₂ by any k leaves a₂ bit-identical (rtol 1e-14, atol 0.0).
- `test_xsec_identity.py:241-272` — θ_S swept continuously; the **magic angle arccos(1/√3) kills the tensor rate shift to <1e-15**.
- `test_target_mass.py:58-81` — `cap = M²/(W2_min−M²) == 0.09654` (rel 1e-3), plus the cell-wise tighter bound `4M²x(1−x)/(W2_min−M²)`, and the six published per-configuration maxima pinned exactly: `[0.0258, 0.0332, 0.0577, 0.0577, 0.0854, 0.0854]`.
- `test_target_mass.py:184-236` — the γ² shift scales exactly as 1/Q² at fixed (x,y) (100× Q² divides it by 100 to 5%) and collapses to `(1+γ²)` to 8%/2% at y = 0.05/0.01.
- `test_hfs.py:198-231` — analytic energy-scale levers `dlnx/dlnE′ = 2−y` and `dlnx/dln(hadronic scale) = −(1−y)`, matched to abs 0.05; the companion test documents that the Gaussian stand-in collapses this lever to exactly 1.0.
- `test_review_20260828.py:51-81` — the ratio-inversion Jacobian `T = R(1+u)/(σ_P²−P̄R)` propagated variance vs exact analytic derivative, rtol 1e-9.
- `test_bookkeeping.py` — every `pzz_true` branch re-derived by a **hand-solved spin-temperature ladder** (bisection), deliberately *not* via `spin.populations_maxent`, i.e. not reusing the module under test. Pins the rational anchor `(P_z, P_zz) = (8/13, 4/13)` for J=1 and `(0.7, 0.4)` for J=3/2 at abs 1e-12, the geometric-ladder identity `p₀p₂ = p₁²`, and the small-P_z coefficients `(3/4)p_z²` (J=1) and `(18/25)p_z²` (J=3/2) at rel 1e-6.
- `test_recopseudo.py` — `test_pzz_scale_error_propagates_one_to_one_not_quadratically` pins `Â/A = P_zz,true/P_zz,assumed` at rel 1e-6, and `test_relative_luminosity_bias_is_one_third_per_unit_ratio_error` pins the analytic coefficient `−(P₁+P₂)/(P₁−P₂) = 1/3` exactly (rel 1e-9).
- `test_radiative.py` — `method_bias_table` pins the electron-method rows at `1.351` / `0.740` (abs 5e-3) and the four sweet-spot y and x distortions at `[9.156, 4.184, 8.450, 4.196]` and `[0.109, 0.239, 0.118, 0.238]`.
- `test_coherent.py:81-91` — the angular near-beam cut across the three configurations: acceptance `0.67 / 0.093 / <1e-7`; `test_angular_cut_kills_the_upper_energies` shows that at the machine's own divergence the coherent tag is dead everywhere (`max < 1e-5`).

**Tier 3 — statistical/estimator closure** at 5–35% tolerances (Gate 3): `test_pseudoexp.py`, `test_likelihood_fit.py`, most of `test_recopseudo.py` and `test_run_share.py`.

### Two tests that exist specifically to catch a *silent* regression

These are worth copying into any new generator's suite, because both guard failure modes that produce plausible-looking wrong numbers:

- `test_sampler.py:89-120` `test_vector_rate_sign_and_labels` — a **deterministic** guard on the populations↔m pairing. Its docstring notes the A∥ truths are smaller than the MC tolerances of every other test, so this is *the only* test that would catch a reversed vector-sector sign.
- `test_run_share.py:132-153` `test_the_sampler_cross_sections_are_share_invariant` — pins that a run-plan share moves event counts, never cross sections. The bug it guards reported **N_1yr = 1.17×10⁷ where 4.67×10⁷ was right**, with errors 4× rather than 2× the published ones, because the share was applied twice.

### Where the suite is deliberately *not* a physics assertion

`test_tensor_convention.py:167` `test_spin32_rate_and_cos2phi_channels_are_now_consistent` is explicitly labelled a **characterization** test: it pins the adopted spin-3/2 rank-2 normalization, not a derived result. Before 2026-08-25 the rate and cos 2φ channels disagreed by a factor 3. A new generator inherits this as a *convention choice to make*, not a target to hit — it sits under `plans/04` #14, the unclassified spin-3/2 basis.

Likewise `test_the_program_sign_is_opposite_to_the_literature` (`test_tensor_convention.py:72-84`) asserts `TENSOR_LL_SIGN == +1.0` with a failure message naming `plans/08` D1 — the test is written so that flipping the convention is a deliberate act with a visible consequence, and it doubles as the switch.