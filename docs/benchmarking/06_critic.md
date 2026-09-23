<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 06 — Completeness critic: what the survey missed, what it got wrong, and what can never be checked

**Critic date 2026-09-06.** Sixth and last file in the benchmarking series. It
does not add a domain; it audits the other five (`00`–`05`, 201 resources) against
four questions:

1. which LiPolGen observable has **no** external benchmark of any kind, even by
   decomposition, and how far the closest proxy actually is;
2. which surveyed resources were **refuted or downgraded** by the adversarial
   verification pass, and what replaces them;
3. which **search angle was never run** — and what running the cheap ones now
   adds (§4 is the answer, and it is the largest section, because five of them
   returned something real);
4. the **per-nucleon / per-nucleus normalisation trap**, resource by resource.

Row labels from the earlier files (`A-…`, `C-…`, `D-…`, `E-…`, `F-…`, `T-…`,
`CH-…`) are used without repeating their content. New resources found by this
file are labelled **`N-…`**.

Everything marked *verified* was fetched, downloaded, read or run on the critic
date. Nothing was written outside `docs/benchmarking/`; nothing was downloaded
into either repository (the two PDFs and two data files fetched here live in a
session scratchpad and are cited by URL, not shipped). Four read-only runs
against the existing `build/` are quoted in full where they carry a number.

---

## 1. The honest headline, restated after verification

`00` §10 said the tree's external validation is "real but narrow… strongest at
A = 2". The verification pass narrows it further, and the single most important
sentence in this file is:

> **After verification there is no surviving like-for-like external check of a
> tensor observable anywhere in the tree, at any A.** The one that came closest
> — `E-1`, Cosyn–Weiss II TABLE II, which `00` §5 calls "the single strongest
> external physics check in the tree" — came back **refuted at blocker severity
> on like-for-likeness**. The world's only tensor-polarized DIS measurement
> (`C-7`, HERMES b₁ᵈ) remains present in this tree **as a comment**, not as a
> gate.

What survives at full strength is what it always was: transcription gates
(`C-6` CD-Bonn, `E-5`/`E-6` Cosyn conventions, `E-7` HJM, `E-8` E143, `E-9`
R1998), one configuration-dependent A = 2 shape test (`E-2` CDKS), one ab-initio
wave-function read-back (`E-13`), unpolarized nuclear inputs, and the whole
chain layer. That is a respectable spine for a generator. It is not validation
of the physics the generator exists to produce, and §2 below is the list of
things that cannot be validated at all until somebody takes data.

---

## 2. Question 1 — the "cannot be benchmarked until data exists" list

Twelve observables. For each: the closest external proxy that exists anywhere,
and **how far** it is — stated as the specific transfer step that is unchecked,
not as a vague "different regime".

### 2.1 The eight with NO external anchor of any kind

| # | LiPolGen observable | symbol / channel | closest proxy in existence | the distance, precisely |
|---|---|---|---|---|
| **U-1** | **b₁(⁶Li), b₂(⁶Li)** — inclusive tensor structure function | `Li6B1`, `b1_li6_from_deuteron`, `LI6_B1_RANK2_TRANSFER`, `LI6_B1_PER_NUCLEON` | b₁ᵈ at A = 2: HERMES (`C-7`, data), CDKS (`E-2`), Miller (`E-3`), Kumano–Kuroki 2026 (`T-1`) | The proxy constrains the **input** of the convolution, never the convolution. Three unchecked steps sit between: (i) the α–d overlap replacing the free deuteron wave function; (ii) `LI6_B1_RANK2_TRANSFER` = 0.921947, a **single number with nothing behind it**; (iii) `LI6_B1_PER_NUCLEON` = 2/6, which assumes an **inert α**. INSPIRE full-text search run for this file (§4.2) confirms independently: **no calculation of b₁ for any A > 2 exists.** |
| **U-2** | **Tagged tensor asymmetry with an α or t spectator** | `TaggedModel`, `azz_tensor_curve` on `li6_alpha_channel` / `li7_alpha_channel` | Cosyn–Weiss II A_T∥ for a **nucleon** spectator on a deuteron (`E-1`, `T-6`) — and that gate is now **refuted as like-for-like** (§3), **repaired 2026-09-06** (`07_cw_sign_investigation.md` §8: it runs on the AV18 control CW quote TABLE II for, checks Eq. (6.12) as an identity to 8.88e−16 over 26 880 cells, and the code bug it exposed — an inverted S–D interference sign in `build_amp2` — is fixed) | Even taken at face value the proxy is A = 2, nucleon spectator, IA-only. The α-spectator case changes the spectator from a spin-½ point particle to a spin-0 composite with its own form factor and its own FSI, and no calculation of it exists. **Nothing external touches this channel.** |
| **U-3** | **Coherent ⁶Li tensor cos 2φ / a₂** | `coherent.hpp`, `a2_from_quadrupole`, `eps_b0` | Mäntysaari *et al.* polarized **deuteron** a₂(\|t\|) (`E-12`, `T-42`) | The only polarization-dependent coherent-diffraction calculation in existence, at A = 2. Re-verified in `04` §9: its 8 forward citations extend it to **no** other nucleus. The tree's own `T10b` prints that the shipped `eps_b0` = −0.08 is a **deuteron** number corresponding to a ⁶Li quadrupole **11.4× the measured one**. `T-46` = DOES NOT EXIST is correct. |
| **U-4** | **⁷Li rank-2 sector (J = 3/2)** | the `spin32_*` block; `tests/test_pipeline.cpp` D1 pins it **identically zero** | Jaffe–Manohar NPB 321 (`T-11`, library copy only); Fu–Sun–Dong GPDs (`T-12`); the spin-3/2 positivity bounds (`T-13`) | These give a **basis** and **inequalities**. There is no calculation, no measurement, and no second implementation of a J = 3/2 rank-2 structure function for anything, let alone ⁷Li. `T-13` would be satisfied trivially by zero, so it is not yet a gate. |
| **U-5** | **b₃, b₄, and Δ (double-helicity-flip gluon)** | `toy_delta_gluon`, `C_BAG`, the b₃/b₄ scenario shapes | Detmold–Shanahan lattice A₂ for the **φ meson** (`T-18`); Sather–Schmidt bag sum rule (`T-17`, 31 citations); new this pass: Xie & Lu T-odd gluon TMDs for a tensor-polarized deuteron (`N-9`) | `04` §4 already states the decisive point and it stands: the lattice number is a *valence-gluon* matrix element in a **meson**, `C_BAG` parametrises *exotic glue in a nucleus*. They are different objects. Only the **gluonic Soffer bound** transfers, and a bound is not a benchmark. |
| **U-6** | **Tensor-sector radiative corrections** | `RcScope::TensorAll`, `RC_DELTA_LOW_X`, `qe_tensor_scale` | Gakh–Shekhovtsova (`T-34`), **zero INSPIRE citations**, one panel at Q² = 0.1; RADGEN's `pnrun = ±2` branch (`01` §2.5), whose b₁ = −(3/2)F₁A₁ is a **placeholder computed from the vector asymmetry** | The RADGEN branch validates the **kernel** (how b₁, b₂ enter the POLRAD integrals) and nothing about the content. `T-38` confirms no modern successor exists. For the **φ-dependent** tensor observables (cos φ_TL, cos 2φ_TT) and for J = 3/2 there is not even a placeholder. |
| **U-7** | **Cluster-spectator Glauber FSI** | `GlauberFsiWeight` | Ciofi degli Atti–Kaptari ³He(e,e′d)X (`T-27`) — the **only** published cluster-spectator FSI, with a **d** spectator in A = 3; Cosyn's `physics-code` (`T-32`), deuteron/nucleon spectator | The published object is a deuteron spectator escaping A = 3. LiPolGen's is an α or t escaping A = 6/7, with a composite σ built from the cluster's own nucleons. `01` §5.4 verified that **no generator computes this**. `T-32` would give the first external FSI number, at A = 2 only. |
| **U-8** | **The α–d D-wave, from data** | `ClusterWaveSource`, `P_D^{αd}`, the quadrupole dial | George–Knutson η(⁶Li → α+d) (`C-1`) — **and this is exactly the entry the verification refuted as a benchmark** (§3) | Every one of the four measured α–d momentum distributions (`03` §3: Ent, Mitchell, (p,pd), (α,2α)) is **unpolarized**, so it constrains \|φ₀\|² + \|φ₂\|² summed and is blind to the S–D interference that b₁ is built from. The 1.93× η discrepancy has no better external anchor and this pass found none. |

### 2.2 The four with a *partial* anchor that is weaker than it looks

| # | observable | the anchor that exists | why it is only partial |
|---|---|---|---|
| **U-9** | **A_zz(⁶Li) inclusive** | `E-6`: A_zz = −(2/3) b₁/F₁ is agreed by four independent sources and gated at 1e-12 | This is a **convention gate on the relation**, not on either side of it. It would pass identically with a b₁ that is wrong by any factor. The verification found **four corrections** to the entry, including wrong HERMES equation numbers. |
| **U-10** | **⁶Li elastic form factors in the RC tail** | `C-2`/`C-3` anchor F_q(0), F_m(0) on the measured Q and μ; the C0 **shape** is an unfitted model band | **New this pass (§4.5, §4.6):** there is more external material than `03` concluded — a machine-readable ⁶Li Fourier–Bessel charge density (`N-5`) and a measured ⁶Li/⁷Li **magnetization**-density row (`N-6`). Both sit **outside** the tree's asserted [2.9, 3.3] fm⁻¹ C0-zero window. So this line moves from "no data" to "data exists and disagrees", which is worse for the tree and better for the survey. |
| **U-11** | **The spin bookkeeping (P_z, P_zz, the ladder, the flip plans)** | Nothing, in `00`–`05`. Table A-7 is an INTERNAL `polligen` pin. | **New this pass (§4.1):** this is the one gap the survey called unbenchmarkable that turns out to have a *clean external gate already sitting in `refs/`* — EPIOS Table II (`N-1`). See §4.1; it is now rank 1 of the ten. |
| **U-12** | **⁶Li quasi-elastic and elastic radiative tails, polarized** | `03` §4: 133 measured **unpolarized** ⁶Li QE points exist and are machine-readable | The tensor QE tail (`qe_tensor_scale`) is a "borrowed magnitude". No polarized QE or elastic measurement exists for any A = 6 spin-1 nucleus, and none was found this pass. |

### 2.3 One correction to the survey's own "does not exist" claims

`03` §8 item 4 states: *"A machine-readable ⁶Li elastic form factor. No HEPData,
no SOG in de Vries, no archive."* **The last clause is wrong** — see §4.5. The
tree's own `docs/open_items/physics_literature.md:260` already suspected it; this
file downloaded the file and confirmed it.

---

## 3. Question 2 — what the verification refuted or downgraded

Sixteen high-priority entries were adversarially checked. **Four came back
refuted**, two more carry major-severity corrections while surviving, and ten
survive with minor corrections. The finding texts handed to this file are
**truncated**; where a finding is cut off mid-sentence that is said, and the full
verification record must be read before acting.

### 3.1 Refuted — the four

| entry | severity | what was refuted | what is **not** refuted | what replaces it |
|---|---|---|---|---|
| **`C-1` George–Knutson η(⁶Li → α+d)** | major | **The claim that it benchmarks what it says it benchmarks.** The verdict is explicit: "REFUTED as a benchmark of the claimed thing." | The reference (INSPIRE 522531, PRC **59** (1999) 598–606), the value η = −0.025 ± 0.006 ± 0.010, its transcription into `LI6_ETA_DS_GK{,_STAT,_SYST}`, and the honesty of the "in-tree" obtainability label. The number is real and correctly typed. | **Nothing does.** `00` already knew the gate asserts `1.2 < \|η/η_GK\| < 3.0`, i.e. it pins a documented 1.93× discrepancy, and `T22b` shows η is **exactly linear** in the quadrupole dial — so it is a relabelling of that dial. The correct action is to **retitle the row** (it is a *consistency band on a dial*, not a D/S benchmark) and to demote it out of the DATA class. Candidate replacement, still absent: an ab-initio α–d D/S ratio (`T-23`, Hebborn–Brune–Phillips `arXiv:2510.19067`, "read the paper first"). |
| **`D-3` POLRAD 2.0 deuteron elastic** | major | Not the core: "CORE OF THE CLAIM STANDS; THE ENTRY AS WRITTEN **OVERSTATES ITS SCOPE IN FIVE PLACES** AND CARRIES ONE STALE NUMBER." | The three σ^el numbers are genuinely POLRAD's — the verifier re-downloaded CPC `ADGH_v1_0` (Mendeley `10.17632/37vgvzgr2w.1`, CPC licence, sha256 `83703668…`) and reproduced them. *(2026-09-23 re-drive with POLRAD's own `ffdeu`/`qunc8`: rows 1 and 3 reproduce to every printed digit; row 2 does not — σ_u = −1.0303e−03, σ_q/σ_u = +0.0623 — so the stale number is that row, `../open_items/run_2026-09-23/phase_B3_chain_rc.md` §3.)* | The entry, rewritten to its actual scope, plus the stale number corrected. The **substantive** replacement is `01` §4: DJANGOH/HERACLES as an independently written second one-loop engine — and its published Rad=1/Rad=0 table (`01` §2.4) needs no code at all. *(2026-09-23: that table was generated with the elastic radiative tail OFF, IEL2 = IEL31..33 = 0 in every log, so it shares no term with `rc_tail` and cannot be this comparison; `../open_items/run_2026-09-23/phase_B3_chain_rc.md` §1.2.)* |
| **`D-4` eSTARlight ⁶Li baseline** | major | **The "claimed to benchmark" line.** "It is not a like-for-like benchmark of `coherent.hpp`'s channel rate": `coherent.hpp` generates **coherent diffractive DIS with a continuum M_X ≥ 1.2 GeV**, while eSTARlight produces **exclusive vector mesons**. Different final state, different observable. | The run itself, the references, the commit, the tables and the O5 reach arithmetic all check out. | The table stays as what it always physically was — an **unpolarized rate and \|t\|-slope scale for exclusive VM production on ⁶Li**, and the input to the O5 J/ψ reach arithmetic. It must stop being described as a check on the coherent **channel rate**. Nothing else does ⁶Li (`01` §5.3: Sartre has no light nucleus); the nearest independent contact is Guzey *et al.* on ³He/⁴He (`T-44`). |
| **`E-1` Cosyn–Weiss II TABLE II** | **blocker** | **Like-for-likeness.** "The claimed benchmark is NOT like-for-like and its pass masks a physics …" — *the finding text handed to this file is truncated at exactly this point.* | The reference is real and cited accurately: `arXiv:2603.23700`, JLAB-THY-26-4661, 24 Mar 2026, 45 pp, CC BY 4.0; the local PDF p. 35 carries Eqs. (6.11)–(6.14) and TABLE II exactly as described; the "in-tree" obtainability label is honest. | **Nothing.** This was `00`'s strongest external physics check. Until the truncated finding is read in full and acted on, the tree should be described as having **no** like-for-like external check of a tensor observable. The nearest unexploited substitute is `T-8` (Cosyn–Weiss PRC 102 (2020) 065204), which gates the **vector** tagged sector — a different observable, and it would not repair this. |

### 3.2 Survived, but with major-severity corrections

| entry | what stands | what must change |
|---|---|---|
| **`C-2` TUNL Q(⁶Li) = −0.0818(17) fm²** | Not refuted. The reference is in-tree and says exactly this: TUNL A = 6 revised manuscript, ⁶Li ground-state block, "Q = −0.818(17) mb" = −0.0818(17) fm², 1998CE04 = Cederberg *et al.* PRA **57** (1998) 2539; also public at `nucldata.tunl.duke.edu`. | Severity is flagged major, so the entry's wording needs the same audit the others got. Note `03` §5 independently bounds the compilation spread at ±2 % (−0.0806 … −0.083), which is negligible beside the model's 7.52× discrepancy — the choice of compilation is **not** a live systematic. |
| **`E-12` Mäntysaari a₂** | "CORE CLAIM STANDS." Reference confirmed: PLB **858** (2024) 139053, DOI `10.1016/j.physletb.2024.139053`, `arXiv:2408.13213` (**only v1 exists**; the local copy is that). | "ONE OF ITS THREE **'VALIDATES' POINTS IS FALSE**, AND THE PRECISION / OBTAINABILITY / SCOPE WORDING OVERSTATES." One of the three things `00` E-12 claims this validates does not follow. Until the full finding is read, treat E-12 as validating **only** the closed-form quadrupole → a₂ map for the deuteron, and nothing about precision. |

### 3.3 Survived with minor corrections — the pattern worth naming

`C-6` (CD-Bonn: four sub-claim corrections), `C-7` (HERMES b₁: Table II has
exactly six rows, as claimed), `D-1` (PYTHIA: three descriptive claims to
tighten; manual is *SciPost Phys. Codebases* **8** (2022)), `D-2` (MSTW: "in-tree"
honest **with one precision — the grid is not inside the LiPolGen repository**,
it is under `deps/`), `E-2` (CDKS: `arXiv:1702.05337` = PRD **95** (2017)
074036, p. 9 Fig. 4 confirmed), `E-3` (Miller: scope "accurate and honestly
narrow"), `E-5` (Cosyn EPJ A: three statements to correct), `E-6` (**the HERMES
equation numbers in the entry are wrong**), `E-13` (VMC: like-for-like, the
ρ_L = A_L²/(4π) convention reconciles the 2004 amplitude family), `E-16` (Chang:
four wording corrections, none changing the O5 arithmetic).

**The pattern:** every single refutation and nearly every correction is about
**wording, scope and like-for-likeness**, never about a mistyped number. Not one
transcription failed. That is a real compliment to the tree and a real
indictment of the survey's *claim* lines — which is exactly the failure mode this
series was written to prevent, reproduced inside the series itself.

---

## 4. Question 3 — the search angles that were never run

Ten angles. Six were run for this file; four are named with a reason for not
running them. Five returned something the survey does not contain.

| angle | run? | result |
|---|---|---|
| polarized-target / ion-source community (P_z, P_zz) | **yes** | **§4.1 — the largest single find. `N-1`, `N-2`.** |
| INSPIRE **full-text** (`ft`) rather than title/metadata search | **yes** | §4.2 — corroborates the "does not exist" claims independently. `N-9`, `N-10`. |
| polarized-deuteron **beam** community (JINR/Dubna, JEDI, NICA) | **yes** | §4.3 — measured **tensor** observables on a spin-1 nucleus. `N-3`, `N-4`. |
| HEPData direct (not via INSPIRE's mirror) | **yes** | §4.4 — still 403. `02` §0's caveat re-confirmed. |
| the **1974** de Jager–de Vries–de Vries compilation (the survey read only the 1987 one) | **yes** | **§4.6 — a measured ⁶Li *and* ⁷Li magnetization density. `N-6`, `N-7`.** |
| the UVa nuclear-charge-density **archive files** (not the scan) | **yes** | **§4.5 — a machine-readable ⁶Li charge density. `N-5`.** Contradicts `03` §8 item 4. |
| EIC user group / ePIC software generator list | partly | §4.7 — the canonical URL 404s; the substantive answer is that the EICUG **MCEG** working group and the CFNS **EPIOS** community exist and are the right addressees, not a list. |
| CLAS Physics Database (not HEPData) | **yes, reachability only** | §4.7 — reachable, searchable by observable including F₁, F₂, g₁, g₂. Plausible route to EG1b / BONuS tables that INSPIRE does not index. **Queried 2026-09-23** (§4.7): 918 measurements / 19 experiment ids for a deuteron target; BONuS is one F₂ⁿ/F₂ᵖ-vs-W* measurement, no spectator spectrum; EG1b inclusive only; the only spectator-momentum distribution is Deeps (Klimenko 2006) at 0.30–0.53 GeV/c. |
| IAEA **EXFOR** | attempted | §4.7 — the two documented query endpoints 404 from here; EXFOR's electron-scattering coverage is thin and the QES archive (`03` §4) already supplies the ⁶Li (e,e′) points. **Low expected yield; not pursued.** |
| PDG / Durham structure-function review compilations | **no** | Deliberately not run. These are plots over the same HEPData records `02` already enumerates; they would add a citation, not a number. Stated so the next reader does not re-derive it. |

### 4.1 The polarized-target and ion-source community — the biggest miss

`00`–`05` contain **no** external resource for the spin-bookkeeping module at
all: `spin_temperature_ladder`, `spin_temperature_pzz`, `populations_maxent`,
`helicity_flip_plan`, `tensor_thirds_plan`, `tensor_flip_plan`,
`azz_rel_lumi_bias`. Table `A-7` pins them at rtol 1e-12 against `polligen`, and
that is an INTERNAL lock. Two external anchors exist and one of them is **already
in `refs/`**.

**`N-1` — EPIOS Table II: a real polarized-ion-source (P_z, P_zz) menu.**
`arXiv:2510.10794` (EPIOS Scientific Consortium, *Realizing the Scientific
Program with Polarized Ion Beams at EIC*; PRC **113**, 060501 (2026)) is
**already cited by this tree** for the γ-synchronisation windows (`C-12`,
`constants.hpp:80`). Its **Table II** — read out of the local copy
`PolarizedLithiumSim/refs/2510.10794.pdf` for this file — lists eight
configurations of the COSY/ANKE polarized **deuteron** source produced by three
RF transitions, with the ideal (P_z, P_zz) of each and a measured vector
polarization from the Low Energy Polarimeter at 539 MeV/c (November 2003).

Run read-only against this build (`spin1_populations`, the tree's own
(P_z, P_zz) → substate inversion):

```
 mode   Pz_ideal  Pzz_ideal    EST Pzz     lo=3|Pz|-2   where     Pz_meas/Pz_ideal
   0   +0.0000   +0.000   +0.00000    -2.0000   interior      n/a
   1   -0.6667   +0.000   +0.36701    +0.0000   LOW-EDGE     0.858
   2   +0.3333   +1.000   +0.08515    -1.0000   HI-EDGE      0.855
   3   -0.3333   -1.000   +0.08515    -1.0000   LOW-EDGE     0.906
   4   +0.5000   -0.500   +0.19722    -0.5000   LOW-EDGE     0.790
   5   -1.0000   +1.000     THROWS    +1.0000   LOW-EDGE     0.758
   6   +1.0000   +1.000     THROWS    +1.0000   LOW-EDGE     0.731
   7   -0.5000   -0.500   +0.19722    -0.5000   LOW-EDGE     0.834

 mode  (Pz,Pzz)          spin1_populations(Pz,Pzz)  ->  (p_+1, p_0, p_-1)
   0  (+0.0000,+0.000)   0.333333  0.333333  0.333333
   1  (-0.6667,+0.000)   0.000000  0.333333  0.666667
   2  (+0.3333,+1.000)   0.666667  0.000000  0.333333
   3  (-0.3333,-1.000)   0.000000  0.666667  0.333333
   4  (+0.5000,-0.500)   0.500000  0.500000  0.000000
   5  (-1.0000,+1.000)   0.000000  0.000000  1.000000
   6  (+1.0000,+1.000)   1.000000  0.000000  0.000000
   7  (-0.5000,-0.500)   0.000000  0.500000  0.500000
```

Three things follow, all new:

1. **Zero free parameters, exact rationals.** Every published source mode maps
   through the tree's own inversion onto exactly the pure two-substate (or
   single-substate) population pattern the three RF transitions physically
   produce. Every polarized mode (1–7) lands **on the boundary** of the physical domain
   (min p_m = 0); the unpolarized mode 0 is interior (1/3, 1/3, 1/3). That is a genuine external gate on the spin algebra, obtainable
   today from a PDF already in `refs/`, at an effort of hours.
2. **The tree's DEFAULT is the wrong model for an EIC beam.** `--pzz-mode
   ladder` (`use_explicit_pzz = false`) is the equal-spin-temperature fill.
   At mode 1's \|P_z\| = 2/3 the ladder gives P_zz = **+0.367**; the real source
   gives **0**. Only mode 0 is on the EST curve. EST describes a **solid target
   in thermal equilibrium**, not an atomic-beam ion source; `--pzz-mode typed`
   is the physically right default for an EIC ion beam and the docs do not say so.
   *(2026-09-23, measured — `../open_items/run_2026-09-23/phase_B1_spin_deuteron.md` §1: the EST
   fill reproduces the ideal P_zz of 3 of the 8 modes (0, and the two |P_z| = 1 pure states, where
   the tree throws) and cannot produce P_zz < 0 at any P_z. That an atomic-beam DEUTERON source is
   not an EST system is measured; which fill an EIC LITHIUM beam will have is not — no lithium
   source-mode table exists. No default was changed.)*
3. **An honesty correction to the paper, which must ride with any citation.**
   Table II's caption says "the **measured** vector and tensor polarizations";
   the table prints only a measured **P_z^LEP** column. **The P_zz column is the
   ideal RFT value, not a measurement.** The caption's "relative beam intensities" are not printed either. The measured/ideal vector ratio is
   0.73–0.91, i.e. the real source delivers 73–91 % of nominal (0.731–0.906 re-read 2026-09-23; "76–91 %" until then) — a realisation
   factor LiPolGen's run plans do not carry.

The same paper (§V D, read here) records that the **polarized ⁶Li/⁷Li source is
under development** at ANL + University of Kentucky, and, in §VIII B 4, lists a "Li-6/Li-7 beam polarimeter" (an atomic-beam target for
Li–Li CNI scattering with Breit–Rabi calibration) as **proposed** R&D — i.e.
**there is not yet a polarimeter of any kind for a lithium beam**, and the tensor
case is not addressed. (The sentence "tensor polarization presents additional
challenges, requiring specialized polarimeter development such as Lamb-shift-based
systems or BRP-type" is in §VIII A 2, the *deuteron* source, and ends
"…specifically designed for deuteron ions"; it was misattributed to lithium here
until 2026-09-23.) That is the correct citation for why U-11's measured
counterpart does not exist for lithium.

**`N-2` — the equal-spin-temperature relation, and COMPASS's ⁶LiD target.**
Run read-only:

```
j=1 ladder vs the polarized-target closed form P_zz = 2 - sqrt(4 - 3 Pz^2)
  Pz=0.10  tree=0.007514115483  EST=0.007514115483  diff= 6.1e-16
  Pz=0.30  tree=0.068679208417  EST=0.068679208417  diff=-7.4e-15
  Pz=0.50  tree=0.197224362268  EST=0.197224362268  diff= 1.6e-15
  Pz=0.70  tree=0.409402627941  EST=0.409402627941  diff= 7.7e-15
  Pz=0.90  tree=0.747003591386  EST=0.747003591386  diff=-5.6e-15
```

`spin_temperature_pzz(1, ·)` **is** the polarized-target community's
equal-spin-temperature relation, to the bisection's 1e-13-in-β tolerance: max 1.14e-14 over 19 999 points of (−1, 1) (the five points above are a sample, not the bound), and the tree throws at \|P_z\| = 1. That identity is
uncited anywhere in the tree (grep for `NMR`, `spin temperature`, `equal spin`,
`polarized target` returns nothing outside `bookkeeping.hpp`'s own prose).

The external anchor for it, fetched and read for this file: J. Koivuniemi *et
al.*, *Polarization Build Up in COMPASS ⁶LiD Target*, SPIN 2004 proceedings
(`wwwcompass.cern.ch/compass/publications/technical/koivuniemi_spin2004.pdf`).
Verbatim from that paper: the deuteron and ⁶Li both have spin 1; the nuclear
polarization follows the **Brillouin function** with J = 1 for D and ⁶Li and
**J = 3/2 for ⁷Li** — the same construction as `spin_temperature_ladder(j, pz)`;
and *"We have carefully checked that during the polarization all the nuclei share
the same spin temperature"*, so ⁶Li and ⁷Li polarizations are computed from the
measured deuteron polarization at a common spin temperature.

**And the limitation, which is the honest half:** *"In the absence of quadrupole
splitting, only one narrow nuclear magnetic resonance line, 3 kHz wide, is
seen"*, and the integrated signal is proportional to **p₊ − p₋** alone. So in the COMPASS ⁶LiD target the NMR line carries the **vector** polarization
only, and the tensor polarization was **estimated** from the EST relation, not
measured — Fig. 1's caption: "The tensor polarization can be estimated to be
T = 1/2(1 − 3p₀) = 11 %" (for the deuteron; note their T is **half** the
Cartesian P_zz = 1 − 3p₀ of this tree). No measurement of the tensor
polarization of ⁶Li was found by any search run for this survey: the INSPIRE
full-text search of §4.2 returned none, and EPIOS §VIII B 4 lists a ⁶Li/⁷Li beam
polarimeter only as proposed R&D (no lithium beam polarimeter exists; tensor
polarimetry is not addressed there). (Re-read and checked numerically
2026-09-23, `../open_items/run_2026-09-23/phase_B1_spin_deuteron.md` §2.) Solid
targets that *do* show a quadrupole-split (Pake) doublet — ND₃, deuterated
alcohols — are where both P_z and P_zz are extracted from the lineshape and
where the EST relation is deliberately broken by RF manipulation: Keller, Crabb
& Day, *Enhanced Tensor Polarization in Solid-State Targets*, NIM A **981**
(2020) 164504, `arXiv:2008.09515` (abstract verified; body not read, so the
(P_z, P_zz) pairs in it are **table-in-paper, unread**).

### 4.2 INSPIRE full-text search — an independent confirmation of three negatives

`02` and `04` used title/metadata queries. Full-text (`ft`) queries were run here:

* `ft "tensor polarized" and ft "lithium"` → **12 records, none of them a
  tensor-polarized-lithium calculation, measurement or target** (they are EIC
  and NICA programme documents, a long-range plan, conference proceedings).
  This independently corroborates `04` §0 item 1 by a different search route.
* `t b1 structure function and (t nucleus or t nuclear or t lithium or t helium)`
  → **zero hits**. No b₁ for A > 2.
* a sweep of everything tagged tensor-polarized since 2025 → the only new
  physics items are Kumano–Kuroki (`T-1`, already surveyed) and **`N-9`** below.

**`N-9`** — Xie & Lu, *T-odd transverse momentum dependent gluon distributions
for tensor polarized deuteron in a spectator model*, `arXiv:2607.23692`
(26 Jul 2026, 24 pp, 3 figures; abstract verified). Six T-odd gluon TMDs for a
tensor-polarized deuteron, with numerical x and k_T results. Relevant to `U-5`
as a second modern statement about tensor-polarized **gluons**; the abstract does
**not** say whether the gluon transversity Δ is among the six, so it is a
"read the body first" item, not a ready gate. A = 2, figures only.

**`N-10`** — the 2026 sweep also returned "Tensor analyzing power T₂₀ in
γd → pn at 350–680 MeV" (`arXiv:2608.09169`) — photodisintegration, outside the
generator's domain, listed so the negative is on the record.

### 4.3 The polarized-deuteron **beam** community — measured tensor observables

`02` searched the tensor-DIS literature and correctly found one datum. It did not
search the **relativistic polarized-deuteron beam** literature, where tensor
observables have been measured for thirty years on a spin-1 nucleus. Verified via
INSPIRE:

**`N-3`** — deuteron **tensor analyzing powers T₂₀ and A_yy in relativistic
deuteron breakup**, JINR/Dubna and SATURNE. Records confirmed: *Measurement of
the tensor analyzing power T₂₀ in inclusive deuteron breakup at 9 GeV/c*, PLB
**387** (1996) 37 (INSPIRE 328522); *Measurement of polarization transfer and the
tensor analyzing power in polarized deuteron breakup with deuteron momenta up to
9 GeV/c*, PLB **325** (1994) 327 (INSPIRE 363903); *Measurement of the
tensor-analyzing power A_yy in deuteron breakup at 4.5 GeV/c*, Few Body Syst.
**32** (2002) 127 (INSPIRE 602774); *New data on tensor analyzing power A_yy of
the relativistic deuteron breakup as additional test of deuteron structure at
small distances*, PLB **595** (2004) 151 (INSPIRE 1464537); Ladygin *et al.*,
`nucl-ex/0510050`, Phys. At. Nucl. **69** (2006) 852.

**Why this matters and exactly how far it is.** These are **measurements of a
rank-2 observable on a spin-1 nucleus, sensitive to the S–D interference, at
internal momenta 0.2–1.0 GeV/c** — the same physical object (`f₀`, `f₂` and their
interference) that drives LiPolGen's tagged tensor asymmetry, in the same
momentum region the tagged spectator spectrum populates, and far beyond the
≤ 500 MeV/c that `02` §1.3's BLAST entry reaches. It is the **closest measured
proxy in existence for `U-2`**, and `00`–`05` do not contain it.

The distance, stated so nobody over-reads it: (i) it is a **hadronic** probe —
the reaction mechanism (deuteron fragmentation on a nuclear target) carries its
own uncertainties and the extraction of n(k) from it is famously
model-dependent; (ii) it is A = 2, nucleon spectator; (iii) the light-cone
variable is the deuteron's, not a tagged α's; (iv) **the analysing-power
convention is not the target-polarization convention** — see §5, trap T-8.
Obtainability: **table-in-paper / digitizable-figure**; not in INSPIRE's HEPData
index.

**`N-4`** — the same community's 2026 review, which is the right citation for the
*flip* side of the bookkeeping module: B. Gou, V. P. Ladygin, N. N. Nikolaev,
F. Rathmann, A. J. Silenko, Y. N. Uzikov, *Fresh Look at Polarized Deuterons at
the Nuclotron/NICA & HIAF and Beyond*, Natural Science Review **3** (2026)
200802, `arXiv:2608.20945` (abstract verified). It treats **frequent tensor
polarization flips in storage rings**, the **decoherence of tensor polarization
under continuous RF spin flips**, and **spin-tensor dichroism as a tensor
polarimeter**. LiPolGen has `tensor_flip_plan` and `azz_rel_lumi_bias`; this is
the accelerator-physics literature behind them, and it is theory/technique, not
data. It also confirms EPIOS's point that tensor polarimetry is an open
instrument problem.

### 4.4 HEPData, retried directly

`02` §0 records `hepdata.net` behind a Cloudflare challenge. Retried for this
file with the JSON API and a browser user agent: `search/?q=…&format=json`,
`record/ins684394?format=json` and `record/ins394050?format=json` all return
**HTTP 403** (5.7 kB challenge page). `02` §0's method — INSPIRE's `data`
collection with the `literature` back-link resolved — remains the only route
from this environment — superseded 2026-09-23: the bare host `hepdata.net` (no `www.`) served ins394050 Table 1 with HTTP 200 to a plain user agent (`../open_items/run_2026-09-23/phase_B2_nuclear.md` §2.1) — and its qualifier ("not in INSPIRE's HEPData index", never
"does not exist on HEPData") stands.

### 4.5 `N-5` — a machine-readable ⁶Li charge density, and it moves the C0 zero

`03` §8 item 4 says no archive has one. **Downloaded and verified:**
`https://discovery.phys.virginia.edu/research/groups/ncd/dldata/FB_data.dat`
(HTTP 200, 18 306 bytes; the Nuclear Charge Density archive, maintained by
Donal Day, stated provenance "Atomic and Nuclear Data Tables, volumes 14, 36 and
60"). It contains a **⁶Li Fourier–Bessel row**: A = 6, Z = 3, seven coefficients

```
a1..a7 = 1.6353e-2  2.9603e-2  2.0807e-2  7.2731e-3  4.7580e-4  -8.7510e-4  1.1413e-3 ,  R = 6.0 fm
```

Evaluated here (ρ(r) = Σ aₙ j₀(nπr/R) for r ≤ R):

```
  4*pi*int rho r^2 dr = 2.991170          (should be Z = 3; a 7-term FB truncation)
  <r^2>^(1/2)         = 2.5206 fm         (e-scattering analyses 2.54-2.57; Angeli 2.589)
  first zero of F_C(q)= 2.6944 fm^-1   (closed form; 2.6950 in this file's first printing)
  measured |F_L|^2 diffraction minimum (Li et al. 1971, q^2 = 8 fm^-2) = 2.8284 fm^-1
  LiPolGen HoSpin1FF ships q0 = 3.0998 fm^-1, test window [2.9, 3.3]
```

So the row is self-consistent and physical, and it gives a **second, independent,
data-derived** number for the ⁶Li C0 zero. **Both external numbers (2.694 and
2.828 fm⁻¹) lie below the tree's asserted [2.9, 3.3] fm⁻¹ window.** `03` §2.4
found one of them and called it "a substantive finding"; there are two, and they
bracket each other rather than the model.

**Provenance is the one open item, and it must be closed before use.** I read
both compilations to trace it and **⁶Li is in neither Table IV**: de Jager, de
Vries & de Vries, ADNDT **14** (1974) 479, Table IV covers ³He, ³He, ⁴He, ⁴He,
¹²C, ³²S, ³⁹K, ⁴⁰Ca, ⁴⁰Ca, ⁴⁸Ca, ²⁰⁸Pb; de Vries, de Jager & de Vries, ADNDT
**36** (1987) 495, Table IV covers ³H, ³He, ¹²C, ¹²C, ¹⁵N and heavier — no
lithium in either (both PDFs fetched and their Table IV blocks read). The row
must therefore come from the archive's third stated source (ADNDT **60**) or from
the archive maintainers. The tree's own `docs/open_items/physics_literature.md:260`
already flags exactly this ("verify against the primary source"). Mark it
**UNVERIFIED provenance, verified content**. FB_data.dat re-fetched 2026-09-23, 18 306 bytes, sha256 `e9770460…3577c5d`; the ⁶Li row is vendored in `validation/benchmarks/data/uva_ncd_fb_data_li6.dat`. The archive states no licence.

Effort if the provenance holds: **hours.** `TabulatedSpin1FF::from_data_dir("ff/li6_elastic.csv", …)`
is implemented and its data directory does not exist; this fills its **C0**
column from seven numbers with no digitization at all — which is cheaper and
harder to argue with than `04`'s rank-1 item (digitizing WS98 Fig. 1), and it is
**data** rather than theory. It does **not** supply C2 (`03` §2's "what none of
this validates" stands: no experiment has separated F_C2 for ⁶Li).

### 4.6 `N-6`, `N-7` — the ⁶Li and ⁷Li **magnetization** densities nobody looked for

The survey read the **1987** compilation, which is charge-only. The **1974**
compilation has a fifth table the 1987 one does not: **Table V,
Magnetization-Density-Distribution Parameters**, p. 502. Fetched
(`https://www.cns.s.u-tokyo.ac.jp/~gunji/tmp/tmp2/AtomicData_NuclearData_14_485.pdf`,
HTTP 200, 2.0 MB) and read. Its lithium rows, verbatim:

| nucleus | μ [μ_N] | ⟨r²⟩^½_mag [fm] | R | a | (s.p. value) | q-range [fm⁻¹] | ref |
|---|---|---|---|---|---|---|---|
| ⁶Li | 0.822 | **3.42(6)** | 2.16(4) | 0.53(3) | 0.438 | 0.85–1.39 | Ra66 |
| ⁶Li | | 2.81 | 1.78 | 0.10(3) | 0.438 | 0.85–1.39 | Ra66, a₀ fixed from charge scattering |
| ⁷Li | 3.256 | **2.69(13)** | 1.70(8) | 0.09(4) | 0.117 | 0.70–1.97 | Ra66 |
| ⁷Li | | 2.72 | 1.72 | 0.08(4) | 0.117 | 0.70–1.97 | Ra66, a₀ fixed |
| ⁷Li | | **3.00(5)** | 1.90(3) | 0 | 0.117 | 0.25–0.90 | Ni71 |

`Ra66` and `Ni71` were resolved from the compilation's own reference list:

**`N-6` = R. E. Rand, R. F. Frosch, M. R. Yearian, Phys. Rev. 144 (1966)
859–873** — verified through INSPIRE: *"Elastic Electron Scattering from the
Magnetic Multipole Distributions of Li-6, Li-7, Be-9, B-10, B-11 and N-14"*,
180° backscatter geometry, electrons to 230 MeV. **The abstract's own reading of
its ⁶Li result: "The Li6 results may also be interpreted in terms of an 'alpha
plus deuteron' model."**

**`N-7` = G. J. C. van Niftrik, L. Lapikás, H. de Vries, G. Box, Nucl. Phys.
A174 (1971) 173** — the ⁷Li magnetization paper `03` §2 lists as "paywalled,
body not read". Its Table V row is now readable without the paywall.

**What these validate.** `00` row `C-3` says of the ⁶Li magnetic form factor:
"the q-dependence… is an unfitted **model band** ([2.9, 3.3] fm⁻¹) with **no
elastic ⁶Li data in this repository** to close it". These rows close the
**q → 0 slope**: F_m(q) ≈ 1 − q²⟨r²⟩_mag/6 with a measured ⟨r²⟩^½_mag. `HoSpin1FF`
anchors F_m(0) on μ (`C-3`) and fixes the slope with the unfitted (q_z, b) = (1.30 fm⁻¹, 1.85 fm): r_mag(⁶Li) = 2.9469 fm (measured 2026-09-23, `../open_items/run_2026-09-23/phase_B2_nuclear.md` §4), uncompared until RFY66 / ADNDT 14 is held; these rows would test it.

**What they do not.** (i) They are **magnetization-density model parameters**,
not tabulated F_M(q) — the underlying form-factor points are in the 1966 and 1971
journal articles (paywalled, **table-in-paper, unverified whether tabulated**).
(ii) The ⁶Li value is a **band, not a number**: free fit 3.42(6) fm vs 2.81 fm
when a₀ is constrained from charge scattering, a 22 % model spread. Any gate must
carry both rows. (iii) They are the **M1** form factor, and say nothing about
C0 or C2. (iv) The q-ranges are narrow (0.85–1.39 fm⁻¹ for ⁶Li).

**The physics point worth carrying into the docs regardless of any gate:** the
⁶Li **magnetization** rms (2.81–3.42 fm) is substantially **larger** than its
**charge** rms (2.54–2.57 fm). The magnetism of ⁶Li lives further out than its
charge — which is what an α core plus a spatially extended, spin-carrying
deuteron looks like, and it is what Rand *et al.* said in 1966. That is the only
**measured** electromagnetic signature of α+d clustering in ⁶Li this whole survey
series has turned up, and it is worth one paragraph in `PHYSICS_CHANNELS.md`
even if no test is ever written.

### 4.7 The remaining angles, briefly

**`N-8` — a second ab-initio ⁶Li charge form factor.** Dytrych, Hayes, Launey,
Draayer, Maris, Vary, Langr, Oberhuber, *Electron-scattering form factors for ⁶Li
in the ab initio symmetry-guided framework*, PRC **91** (2015) 024326,
`arXiv:1502.03066` — SA-NCSM, charge form factor to q ≈ 4 fm⁻¹, i.e. **past the
diffraction minimum**. It is **absent from `00`–`05`** but **not new to the
repository**: `docs/open_items/run_2026-09-02/design_C_tensor_rc.md:1105` already
cites it and correctly notes it is **C0 only**. Paired with WS98 (`T-20`) it
turns `04`'s rank-1 item from one digitized curve into a two-calculation
bracket on the C0 zero — which is exactly what §4.5's disagreement needs.
*A survey-hygiene point: the benchmarking series was not cross-checked against
`docs/open_items/physics_literature.md`, and two of this file's findings
(`N-5`, `N-8`) were already sitting there.*

**EIC generator list / ePIC wiki.** `https://eic.github.io/software/generators.html`
and `https://eic.github.io/software/` both **404** on the critic date, so there
is no canonical list to point at; `01` §2's per-repository verification remains
the better evidence. The organisational answer is that the addressees are the
EICUG **MCEG** working group (for a generator) and the **CFNS EPIOS** community
(`indico.cfnssbu.physics.sunysb.edu/event/343`, "Polarized Ion Sources and Beams
at EIC", March 2025) for anything about a polarized lithium beam. `05` §2.5's
conclusion — there is no ePIC generator-validation benchmark suite to plug into —
is unchanged.

**JLab tensor community.** It exists and has a name-able membership: Slifer &
Long (`arXiv:1311.4835`, *Novel Physics With Tensor Polarized Targets*, PSTP 2013
— verified: a workshop report, no EST relation and no P_zz values in the
abstract), Keller/Crabb/Day (`arXiv:2008.09515`), Farooq (JLab Indico event 949,
the E12-13-011 b₁ target experiment). Status as of the critic date: **E12-13-011
approved, not run**; the community is analysing CLAS12 **RG-C** data (which
carries a small incidental tensor polarization) as preparation, via a CAA
proposal. `02` §1.2's "does-not-exist (yet)" is correct and unchanged.

**CLAS Physics Database.** `https://clas.sinp.msu.ru/cgi-bin/jlab/db.cgi`
reachable (HTTP 200, 45 kB), searchable by observable including F₁, F₂, F_L, g₁,
g₂ and asymmetries. It is a **plausible machine-readable route to the EG1b and
BONuS numbers** that `02` F-9 / §5.5 and D-9 could only reach as digitizable
figures. **Queried 2026-09-23** (`../open_items/run_2026-09-23/phase_B1_spin_deuteron.md` §4): for a deuteron target, 918 measurements in 19 experiment ids. BONuS is one measurement (F₂ⁿ/F₂ᵖ vs W*), EG1b is inclusive (g₁, A₁, g₁/F₁), and the only spectator-momentum distribution is Deeps (Klimenko 2006, F₂ᴺ × P(p_s) at 0.30–0.53 GeV/c, 115 tables, vendored unwired). The BONuS p_s tables are in PRC 89 045206's Supplemental Material (HTTP 401 without an APS login).

**EXFOR.** Both documented query endpoints 404 from this environment. EXFOR's
electron-scattering coverage is thin and the QES archive (`03` §4, 133 ⁶Li
points, already machine-readable) covers the one ⁶Li electron-scattering dataset
this project needs. Not pursued; the negative is recorded.

---

## 5. Question 4 — the per-nucleon / per-nucleus normalisation trap

Every DATA resource in the series, with the convention it **publishes in**.
Marked **⚠** where getting it wrong changes a number by a factor and the error is
silent.

### 5.1 The table

| resource | published in | note |
|---|---|---|
| **HERMES b₁ᵈ, A_zz** (`C-7`, `02` §1.1) | **per nucleon** ⚠ | Their Eq. (5) divides by an F₁ᵈ built from F₂ᵈ = (F₂ᵖ + F₂ⁿ)/2. `constants.hpp:85–95` records this and says it was "confirmed by inverting their own Table II in all six bins". A_zz is a ratio and is convention-free; **b₁ is not**. |
| **Miller b₁ᵈ curve** (`E-3`) | **per deuteron** ⚠ | `constants.hpp:99–108`: his Eq. (1) densities are "in a target hadron". `B1_MILLER_TABLE_TO_PER_NUCLEON` = 0.5 is applied, and is labelled **LIKELY, NOT CERTAIN** — Miller's own paper is self-inconsistent by exactly this factor. **Open author decision.** |
| **CDKS b₁ᵈ curve** (`E-2`) | **per deuteron**, with x defined on the **nucleon** mass (x runs to 2) | The two halves of that sentence are independent and both matter. `E-4` records that the digitized column stops at x = 1.59 while the deuteron's x runs to 2, which is why `close_kumano_integral` reports rather than enforces. |
| **LiPolGen b₁(⁶Li)** | **per nucleon**, via `LI6_B1_PER_NUCLEON` = 2/6 ⚠⚠ | **The worst trap in the tree.** 2/6 is *not* a normalisation convention: it is the modelling assumption that the **α contributes nothing** to b₁. Anyone comparing LiPolGen's "per nucleon" b₁(⁶Li) with HERMES's "per nucleon" b₁ᵈ is off by **3** unless they know this. Say it in the same breath, every time. |
| **NMC F₂(⁶Li)/F₂(D)** (`03` §1) | a ratio of **per-nucleon** structure functions ⚠ | Both numerator and denominator are per nucleon, so the ratio is 1 for no nuclear effect. The trap is the **denominator's identity**: it is the *deuteron*, which carries its own few-percent EMC/binding effect at x ≳ 0.5, so the true ⁶Li-to-free-nucleon depletion is **larger** than the NMC ratio (`03` §1.4). Direction: **up** on 0.031. |
| **EPPS21 `…_Li6` LHAPDF grid** (`D-5`) | the **isoscalar-averaged nucleon** of ⁶Li, per nucleon ⚠⚠ | Measured in-tree (`03` §1.3): R_u = 0.578, R_d = 2.829 at x = 0.65 — exactly (1+d/u)/2 and (1+u/d)/2, the n/p isospin average, **not** a nuclear effect. Dividing by the **free proton** gives an apparent ⟨1−R⟩ of 0.268, **nine times** the real number, and it reproduces a plausible-looking EMC curve. The denominator must be ½(F₂ᵖ + F₂ⁿ). |
| **g₁ᵈ** — HERMES D-1, COMPASS D-2…D-5, E143 D-6, E155 D-7 | **per nucleon** ⚠ | And the second half: published g₁ᵈ is conventionally **not** corrected for the D state; the (1 − 1.5 ω_D) factor is applied when extracting g₁ⁿ. Applying it twice, or not at all, is a 7 % error on the isoscalar combination. Verify per paper before wiring D-2. |
| **F₂ᵈ** — NMC F-1, HERMES F-7, SLAC F-5/F-6, BCDMS F-4 | **per nucleon** (the deuteron's average nucleon) | This is the standard in all four collaborations, but it is the assumption behind `constants.hpp`'s "HERMES's b₁ is per nucleon" inversion, so wiring F-7 (`02` §10 rank 2) is precisely the way to make it checkable rather than argued. |
| **F₂ⁿ/F₂ᵖ** — NMC F-3, MARATHON F-8 | a ratio; convention-free | The trap is elsewhere: `ToyF2::f2n_over_f2p` is Q²-frozen at Q² = 10. |
| **H1 F_L / R** (R-5…R-7) | proton; per nucleon trivially | R = F_L/(F₂ − F_L) is a ratio. |
| **QES archive ⁶Li, ⁴He** (`03` §4) | **per nucleus**, σ in nb/sr/GeV ⚠ | Not per nucleon. A factor of 6 for ⁶Li, 4 for ⁴He. The archive supplies **cross sections only** — no L/T separation and explicitly no Coulomb corrections. |
| **eSTARlight coherent σ** (`D-4`) | **per nucleus** ⚠ | Coherent production scales ~A² in amplitude; a per-nucleon reading is meaningless. And after the verification (§3.1) this table is a rate scale for **exclusive VM**, not for `coherent.hpp`'s M_X ≥ 1.2 GeV continuum. |
| **Chang *et al.* far-forward efficiencies** (`E-16`) | dimensionless per-recoil-nucleus efficiency | Not per nucleon and not a cross section. Carries **no** central-detector acceptance and no decay-lepton efficiency; the missing factor is ≤ 1, direction **down**. |
| **AME2020 masses** (`C-4`) | **per nucleus**, atomic mass excess → nuclear mass | `mass_per_nucleon` is derived, not published. |
| **TUNL / Stone moments** μ, Q (`C-2`, `C-3`, `03` §5) | **per nucleus**, μ in μ_N and Q in b or fm² ⚠ | The b ↔ fm² factor is 100 and the compilations mix them: TUNL prints "−0.818(17) **mb**" for what the tree stores as −0.0818(17) **fm²** (verification, §3.2). Get this wrong and it is a factor 10 or 100, silently. |
| **de Vries / de Jager charge radii** (`03` §2) | **per nucleus**, ⟨r²⟩^½ in fm | Model-dependent (MI30 / HO / MI): 2.54(5), 2.56(5), 2.57(10) for ⁶Li. Angeli's 2.589(39) is ~1 % higher than all three. |
| **`N-5` UVa FB charge density** (§4.5) | ρ(r) normalised to **4π∫ρr²dr = Ze**, per nucleus | Verified numerically here: the row integrates to 2.9912 for Z = 3. |
| **`N-6`/`N-7` magnetization densities** (§4.6) | **per nucleus**, ⟨r²⟩^½_mag in fm, μ in μ_N | Model parameters, not F_M(q). Two ⁶Li rows differing by 22 % — carry both. |
| **CD-Bonn deuteron** (`C-6`) | wave function normalised **∫(u² + w²)dr = 1**, per deuteron | A_S in fm^{−1/2}; η dimensionless; Q_d in fm² per deuteron. |
| **ANL VMC overlaps** (`E-13`) | **per nucleus**, with the reconciling convention **ρ_L = A_L²/(4π)** ⚠ | The verification confirms this is the one normalisation that reconciles the 2004 amplitude family with the 2014/2024 momentum files. Two different ANL file generations, two conventions. |
| **`N-1` EPIOS Table II** (§4.1) | dimensionless (P_z, P_zz), **Cartesian** convention, per ion ⚠ | And: the P_zz column is **ideal**, only P_z is measured. |
| **`N-3` T₂₀ / A_yy** (§4.3) | **analysing powers**, spherical (T₂₀, t₂₀) or Cartesian (A_yy, A_zz) ⚠⚠ | **Trap T-8, and it is the reason this resource must not be wired casually.** The spherical tensor moments (t₂₀, T₂₀) and the Cartesian ones (P_zz, A_zz) differ by a factor of order √2 in the Madison convention, and the sign convention differs between beam-analysing-power and target-polarization usage. LiPolGen is Cartesian throughout (`P_zz ∈ [−2, 1]`). **Resolve the convention against the Madison-convention definition before any numerical comparison** — this file did not, and does not assert a factor. |
| **Blinov α–p σ_tot, σ_el** (`C-9`) | **per nucleus** (α), mb | Used as a ratio, so convention-safe as used. |
| **Moniz k_F(⁶Li)** (`C-10`) | a Fermi-gas fit parameter, per nucleus, GeV | An input, not a residual. |
| **HERMES ³He g₁ⁿ, E06-014** (`02` §6) | per **neutron** after effective-polarization correction ⚠ | The published g₁ⁿ already has P_n ≈ 0.86 divided out. Comparing it to a raw ³He number is a 16 % error. |
| **Lapikás–Wesseling–Wiringa ⁷Li(e,e′p)** (`03` §3) | summed **spectroscopic factor**, dimensionless (0.58 ± 0.05 vs 0.60 VMC) | Relative to the independent-particle expectation for the transitions summed — not a per-nucleon density. |

### 5.2 The three rules that fall out of the table

1. **Never compare two b₁ values without naming both normalisations *and* the
   A-scaling assumption.** The tree carries three distinct factors — `0.5`
   (per-deuteron → per-nucleon, an open decision), `2/6` (an inert-α model
   assumption), and CDKS's x-to-2 convention — and none of them is visible in a
   plotted curve.
2. **Nuclear PDF grids are the average nucleon, not the bound proton.** Divide by
   ½(F₂ᵖ + F₂ⁿ), always, and assert it in the test.
3. **Cross sections on a nucleus are per nucleus; structure functions on a
   nucleus are per nucleon.** The QES archive and eSTARlight are the former;
   every F₂ and g₁ row is the latter. Nothing in either convention is written on
   the file.

---

## 6. New resources found by this file

| id | resource | kind | in `00`–`05`? | obtainable | effort | priority |
|---|---|---|---|---|---|---|
| `N-1` | EPIOS Table II — eight measured polarized-deuteron source modes (`arXiv:2510.10794`, PRC 113 060501) | data | **no** (paper cited for beam energies only) | **in-tree** (`refs/2510.10794.pdf`) | hours | **high** |
| `N-2` | Equal-spin-temperature relation + COMPASS ⁶LiD common-spin-temperature check (Koivuniemi SPIN 2004); Keller–Crabb–Day NIM A 981 164504 | data + theory | **no** | free PDF / arXiv | hours (identity) – 1 d (Keller body) | **high** |
| `N-3` | Deuteron T₂₀ / A_yy in relativistic breakup, JINR/Dubna (PLB 387 37; PLB 325 327; FBS 32 127; PLB 595 151; `nucl-ex/0510050`) | data | **no** | table-in-paper / digitizable-figure | 3–5 d + a convention resolution | medium |
| `N-4` | Ladygin *et al.*, NICA/HIAF polarized-deuteron review, `arXiv:2608.20945` | theory | **no** | arXiv, free | 0.5 d (citation) | low |
| `N-5` | ⁶Li Fourier–Bessel charge density, UVa NCD archive `FB_data.dat` | data | **no** (`03` §8 says it does not exist) | **downloadable, verified** | hours | **high** |
| `N-6` | Rand–Frosch–Yearian ⁶Li/⁷Li magnetization densities, PR 144 (1966) 859, via de Jager ADNDT 14 Table V | data | **no** | compilation free / article paywalled | 0.5–1 d | **high** |
| `N-7` | van Niftrik *et al.* ⁷Li magnetization row, NPA 174 (1971) 173 | data | listed but body unread | compilation row free | hours | medium |
| `N-8` | Dytrych *et al.* SA-NCSM ⁶Li form factors, PRC 91 024326 | theory | **no** (but in `open_items/`) | arXiv, free | 1–2 d | medium |
| `N-9` | Xie & Lu T-odd gluon TMDs, tensor-polarized deuteron, `arXiv:2607.23692` | theory | **no** | arXiv, free | read first | low |
| `N-10` | CLAS Physics Database as a route to EG1b / BONuS tables | data | **no** | reachable, unqueried | hours to check | medium |

---

## 7. What is still true after all of it

1. **No polarized e + ⁶Li or ⁷Li DIS data.** Unchanged, and independently
   re-confirmed by INSPIRE full-text search (§4.2).
2. **No generator implements a tensor-polarized target of any species.**
   Unchanged.
3. **No calculation of b₁ for any A > 2.** Unchanged, re-confirmed.
4. **No measurement of P_zz for ⁶Li was found by any search run for this survey** — an absence, stated with its searches, not a proof (reworded 2026-09-23). New: the
   INSPIRE full-text search of §4.2 returned none; the one material in which ⁶Li has ever been polarized (COMPASS ⁶LiD) shows **no
   quadrupole splitting**, so its NMR line carries the vector polarization only
   (§4.1); and EPIOS §VIII B 4 lists a ⁶Li/⁷Li beam polarimeter only as proposed R&D, tensor polarimetry not addressed (the tensor-polarimeter
   sentence once attributed to it here is EPIOS §VIII A 2, on deuteron sources).
5. **No like-for-like external check of a tensor observable survives in the
   tree** (§1, §3.1).
6. Everything in `00` §9, `01` §5, `02` §11, `03` §8, `04` §11 and `05` §7 that
   this file has not explicitly corrected stands as written.

**One-line summary.** The survey's 201 resources are, after audit, honestly
labelled but over-claimed in four places and incomplete in five: the spin
bookkeeping had an external gate sitting unread in `refs/`, the ⁶Li
electromagnetic sector has more measured material than the survey concluded and
that material **disagrees** with the shipped model, a thirty-year body of
measured tensor observables on spin-1 nuclei was never searched for, and the
strongest external physics check in the tree did not survive verification.
