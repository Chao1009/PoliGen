<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-27 — phase R1: which normalisation of Δ the literature uses (the judge's record)

Judge R1, 2026-10-27, on `claude/run-2026-10-26` at `ca9b23a`. Inputs: the
structured findings of the two readers, R1A (1989–2016 and the LOI) and R1B
(2019–2026, the tree and the letter). I re-read each source where they
disagreed or were unsure. I also followed the LOI's and Sather–Schmidt's
citation trails to three sources neither reader had: Ma–Wang–Zhang 2013,
Edelmann–Piller–Weise 1998 and Piller–Weise 2000.

The question is the author's decision 1 of 2026-10-27, "we keep the factor
2, is this more common?". The factor is the kernel's cos 2φ coefficient per
unit Δ, which is twice Detmold–Shanahan's.

**Scope of what was done.** This record is the only repository file
written. No constant, tolerance, registry label or reference number moved,
and nothing was published or sent. Conventions used throughout:

- **read**: printed, and read from a page image unless "text" is said;
- **derived**: computed here or by a reader from printed equations;
- **inferred**: depends on an identification that the source does not print.

> **Which decision this record serves.** The task defines decision 1 as "the
> kernel's cos 2φ coefficient stays as coded"; §4.1–§4.3 and §5a answer
> that reading. This run's `DECISIONS.md` (13:04) records something else:
> once the literature had been read, the author revised D4 to the
> Jaffe–Manohar convention, which halves the coefficient. The D4 design
> (`plan_D4_jm_convention.md`) follows that revision. §4.4 and §5b give the
> same items for it. The literature verdict (§0–§3) is the same under both
> readings.

## 0. Verdict

Every source that derives the inclusive cross section carries the ½. That
covers Jaffe–Manohar 1989 (preprint Eqs. (1) and (16), and its spin-J
appendix independently), Sather–Schmidt 1990 (A7) and Ma–Wang–Zhang 2013
(Eq. (23)). Every source that prints the double-helicity-flip tensor carries
it too: Detmold–Shanahan Eq. (5), Edelmann–Piller–Weise Eq. (7),
Piller–Weise Eq. (3.4) and the LOI's own Eq. (5). For a pure transverse
m = ±1 state they all give a₂ = −½ (1−y)/y² · Δ/D_φ, so Δ_source = 2 Δ_LiPolGen
at equal cross section (k = 2).

The kernel's coefficient (k = 1) appears in exactly one printed equation: the
JLab LOI's Eq. (9), which has no derivation. Its only reprint is PoS
SPIN2018 174, by the same corresponding author. Eq. (9) contradicts the same
LOI's Eq. (6). It reproduces, symbol for symbol, Ma–Wang–Zhang's Eq. (23),
which the LOI cites, except that the ½ of its "(1−y)/(2y²)" is missing. That
makes it a transcription loss, not a second convention.

So the answer to "is this more common?" is **no**. The literature uses
Jaffe–Manohar's normalisation: Δ = A₊₋,₋₊ with F₁ = ⅓Σ_H A₊H,₊H.
LiPolGen's Δ as coded is half of that. Keeping the coefficient is a pure
relabelling, but it is consistent only if every external Δ number is halved
on entry: Sather–Schmidt's bag value, Detmold–Shanahan's relation and bound,
and lattice moments. The letter's moment_A / moment_B scenarios do not do
this.

## 1. The table

**The reference point** is `src/core/xsec.cpp:281`
(`out.a2 + (-(1.0 - y) / (y * y) * qc.second * st * st * t.delta`), in
`tensor_amplitudes`:

  a₂ = −(1−y)/y² · c_eff(m) sin²θ_S · Δ/D_φ

with c_eff = 3Q_NN (`src/core/xsec.cpp:179`, `std::make_pair(q_nn, 3.0 * q_nn)`) and D_φ = F₁ + (1−y)F₂/(xy²).
For a pure m = ±1 state with its axis ⊥ q this is a₂ = −(1−y)/y² · Δ/D_φ.
`a_cos2phi` (`src/core/asymmetries.cpp:90-92`) is the same.

**Definition of k.** Δ_source = k · Δ_LiPolGen for the same physical cross
section. So k = 2 means the source's Δ is twice LiPolGen's for the same
amplitude.

| # | source | polarization state | coefficient per unit Δ (quoted) | k | basis |
|---|---|---|---|---|---|
| 1 | **JM89** — Jaffe, Manohar, "Nuclear gluonometry", PLB 223 (1989) 218, read as MIT preprint CTP#1712 (March 1989, "Submitted to: Physics Letters B"), KEK scan 2000-36-731, sha256 `0873595123b2159a22ab3c908c014937863bc30c0c1d13343ab929b81236d9d1` | "We take E_μ = E′_μ ≡ E^x_{+μ} ≡ −(1/√2)(0,0,1,i)" (p. 5): pure m_x = +1, axis x̂ ⊥ q (rest frame, q along ẑ); "The same cross section is obtained if the target is polarized in the opposite sense, E = E^x_−" (p. 5) | p. 2, Eq. (1): "2π(dσ/dx dy dφ)/(dσ/dx dy) = 1 − ½ · x(1−y)/[xy²F₁(x,Q²) + (1−y)F₂(x,Q²)] · Δ(x,Q²) cos 2φ". p. 5, Eq. (16): "(e⁴ME/4π²Q⁴){xy²F₁ + (1−y)F₂ − (x(1−y)/2) Δ cos 2φ}". p. 9 appendix: "2(1−y)(A_{0H,0H′} − ½A_{+H,−H′}e^{−2iφ} − ½A_{−H,+H′}e^{2iφ}) + [(1−y) + y²/2](A_{+H,+H′} + A_{−H,−H′})", which with Eq. (15)'s ρ₊₋ = ¼ gives the same ½. Result: a₂ = −½(1−y)/y² · Δ/D_φ | **2** | read (pp. 2, 4, 5, 9). The published journal text was not seen (§2.3) |
| 2 | **SS90** — Sather, Schmidt, PRD 42 (1990) 1424, sha256 `63378c24…439dc37a` | "For a target with spin component m_J = 1 along the x axis, the density matrix can be written ρ = E_x^i E_x^{*j}, where E_x^i = (1/√2)(0,1,i)" (p. 1427). The beam defines ẑ; x ⊥ beam | (A7): "dσ/dx dy dφ = (4α²ME/Q⁴)[xy²F₁ + λF₂ + ⅙(xy²Δ₁ + λΔ₂) − ½xλΔ₄cos2φ − …]", with λ → 1−y. (A2): "Δ₄ = (1/κ)Δ = (1/κ)A₊₋,₋₊ … Also Δ = Δ₄ in the Bjorken limit". So a₂ = −½(1−y)/y² · Δ/D_φ. For any ρ, (A3)'s "2xλΔ₄Tr(I₃ρ)" with (A6)'s I₃ gives −½xλ·P_zz·Δ cos 2φ for populations along a transverse axis | **2** | read (p. 1427) |
| 3 | **MWZ13** — Ma, Wang, Zhang, "Azimuthal Dependence of DIS with Spin-1 Target", arXiv:1306.6693v2 (9 Dec 2013, unpublished; INSPIRE 1240510). Cited by the LOI as [2]. Fetched this phase | "H is polarized along the direction m⃗ in the sense that the projection of the spin vector along m⃗ is quantized, i.e. m⃗·S⃗ = λ_m = 0, ±1"; "the initial lepton moves in the −z-direction" (p. 2) | p. 8, Eq. (23): "dσ(λ_m)/dx_B dy dφ = (2yα²/Q²)[F₁ + a_mF̃₁ + ((1−y)/y²)(F₂ + a_mF̃₂) − ((1−y)/(2y²)) b_m sin²θ_m Ĝ_T cos 2φ − …]". Eq. (24): "b_m = 3\|λ_m\| − 2, a_m = ¼b_m(1 − 3cos²θ_m)". Their F₂ multiplies x p̂p̂/(p·q) (Eq. (11)), so F₂^MWZ = F₂/x. Ĝ_T multiplies T⊥^{μν} = ½(ε⊥^μ ε⊥^{*ν} + ε⊥^ν ε⊥^{*μ} − g⊥^{μν} ε⊥·ε⊥*) (Eqs. (5), (11)), which equals JM Eq. (14)'s tensor exactly (§7). Result: a₂ = −½(1−y)/y² · b_m sin²θ_m · Ĝ_T/D_φ | **2** (Ĝ_T ≡ Δ_JM) | read (pp. 3, 4, 8); the tensor identity is derived. The LOI itself writes "Ma et al. [2] call Δ(x,Q²), labeled as Ĝ_T in their paper" (LOI p. 2, text) |
| 4 | **LOI18 Eq. (6)** — Maxwell *et al.*, JLab PAC 44 LOI (June 2016), arXiv:1803.11206, sha256 `30d2bd48…236d` | "φ is the angle between the scattering plane and the target spin orientation. If the target is polarized in the opposite direction, the same cross section is obtained": pure m = ±1 along a transverse axis | p. 4, Eq. (6): "lim_{Q²→∞} dσ/dx dy dφ = (e⁴ME/4π²Q⁴)[xy²F₁ + (1−y)F₂ − (x(1−y)/2) Δ cos 2φ]" (JM Eq. (16) verbatim). Result: a₂ = −½(1−y)/y² · Δ/D_φ | **2** | read (p. 4) |
| 5 | **LOI18 Eq. (9)**, with Eq. (10) | "for the target spin in the m̂ direction λ_m = (1, 0, −1)", axis at θ_m from the beam: pure states | p. 5, Eq. (9): "(2yα²/Q²)(F₁ + ⅔a_mb₁ + ((1−y)/(xy²))(F₂ + ⅔a_mb₂) − ((1−y)/y²) b_m sin²θ_m Δ(x,Q²) cos(2φ))", with "b_m = 3\|λ_m\| − 2". Eq. (10), for N₊ + N₋ − 2N₀: "− 6 ((1−y)/y²) sin²θ_m Δ cos(2φ)". Result: −(1−y)/y² · Δ/D_φ at λ = ±1, θ_m = 90° | **1** | read (p. 5). No derivation is given. b₁/b₂ are attributed to "[6] R. Jaffe, Private communication (2016)" |
| 6 | **Maxwell 2019**, PoS(SPIN2018)174, sha256 `3b341df6…40bf` | as #5 | LOI Eq. (9) reprinted | **1** | read (R1B's page render, re-viewed). A reprint by the LOI's corresponding author |
| 7 | **DS16** — Detmold, Shanahan, PRD 94 (2016) 014507, arXiv:1606.04505v3, sha256 `6cea5915…b50a7` | the contraction uses pure m_x = ±1 ⊥ q | Eq. (5): "W^{Δ=2}_{μν} = ½{[(E′*_μ − …)(E_ν − …) + (μ↔ν)] − [g_μν − q_μq_ν/q² + …][E′*·E + …]}Δ(x,Q²)", which is JM Eq. (14). Δ ≡ A₊₋,₋₊ (Eq. (4)), with the projectors of Eqs. (2)–(3). No cross section is printed. Contracted with the lepton tensor it gives a₂ = −½(1−y)/y² · Δ/D_φ + O(M²/Q²) | **2** | derived, by four independent contractions (§7): phase C2, R1A, R1B and this judge. At x = 0.1, y = 0.5, Q² = 5 GeV² the ratio to the kernel is 0.4964 (0.4961 with F₂ = 0); it is 0.49982 at Q² = 100 GeV² and → ½ at large Q². m_x = 0 gives −2 × the m_x = ±1 value in all four |
| 8 | **EPW98** — Edelmann, Piller, Weise, PRC 57 (1998) 3392, hep-ph/9709455. Fetched | general 𝓔 | p. 4, Eq. (7): "+ (Δ/2){(−g_μν + (2x_T/κ_dν)p_μp_ν)((M_d²/κ_dν²) q·𝓔 q·𝓔* − 1) + [(𝓔_μ − (q·𝓔/κ_dν)p_μ)(𝓔*_ν − (q·𝓔*/κ_dν)p_ν) + (μ↔ν)]}". Eq. (8a): "F₁ = (1/6πe²)(Im A₊₊,₊₊ + Im A₊₋,₊₋ + Im A₊₀,₊₀)". Eq. (8g): "Δ = (1/2πe²) Im A₊₋,₋₊". These are SS90 (A1)–(A2)'s tensor and normalisation. No cross section is printed | **2** | derived (identity with SS90 (A1)). Eqs. (7) and (8a) read from the page image; (8g) from text |
| 9 | **PW00** — Piller, Weise, Phys. Rept. 330 (2000) 1, hep-ph/9908230. Fetched | general 𝓔 | p. 28, Eqs. (3.3)–(3.4): "M_A W^A_μν = … + t_μν Δ^A + …", with "t_μν = ½{(−g_μν + (2x_A/κ_A P·q)P_μP_ν)(… q·𝓔 q·𝓔* − 1) + [(𝓔_μ − …)(𝓔*_ν − …) + (μ↔ν)]}". No cross section is printed | **2** | derived (same tensor) |
| 10 | **PV26** — Pascalutsa, Vanderhaeghen, arXiv:2610.09862v1 (7 Oct 2026, unpublished) | ensemble P_zz along an axis at θ to q | p. 10: Eq. (26) "ε T_TT cos 2φ_TT F_UTTT"; Eq. (27b) "F_UTTT = −(M/3ν)(b₂ − 3b₃)"; Eq. (29b) "T_TT cos 2φ_TT = (P_zz/2) sin²θ cos 2φ"; fn. 5 "F_UTTT = −2M\|q\|τ_TT". App. A (text): (A6a) "W₊₊,₊₊ + W₊₋,₊₋ = 4\|q\|σ_TT", (A6b) "W₊₊,₋₋ = 2\|q\|τ_TT", (A7a) "2\|q\|σ_TT = F₁ − …", (A8a) "∆ = W₊₊,₋₋ = (1/ν)(⅓b₂ − b₃)" with M = 1. Hence F_UTTT = −Δ, and the cos 2φ term is −½(1−y)/y² · P_zz Δ/D_φ | **2** | inferred, medium confidence. Identifying their Δ with A₊₋,₋₊ uses HJM's "A₊₋,₋₊ = (M²/ν)a₃, a₃ = (b₂/3 − b₃)" (R1A, HJM preprint p. 5) and CRSZ25's F_UU,T = 2F₁ |
| 11 | **LiPolGen kernel** | pure m along n̂(θ_S, φ_S); c_eff → P_zz for a fill | `-(1.0 - y) / (y * y) * qc.second * st * st * t.delta / std::max(den, 1e-30)` | 1 | the reference |
| 12 | **The letter**, PolarizedLithiumSim `paper/main.tex` (DRAFT v0, sha256 `556c74f0…7169`): eq:master and eq:amp | "per spin projection m"; "per unit P_zz" at θ_m = 90° | "− ((1−y)/y²) c_m sin²θ_m Δ cos 2φ"; "A^{cos 2φ} = −((1−y)/y²)(Δ/D_φ)" | 1 | read. It cites [Jaffe:1989xy, Hoodbhoy:1988am] for this coefficient, which neither paper prints. Its eq:deltag is DS16 Eq. (13), a k = 2 relation (§4.2) |
| 13 | **The program's lineage**: the sibling's `docs/Discussions.pptx` slide 5 (R1A: "LOI Eq. (9) with b_m renamed c_m"); `refs/refs_dict.json` ("Master formula used throughout the repo (docs/Discussions.pptx p.5)"); `plans/01` §2.1; `fastsim/polli_fastsim/asymmetries.py` `a_cos2phi`; `evgen/polligen/xsec.py` | as #5 | LOI Eq. (9)'s coefficient | 1 | read (R1A, R1B). These are copies of #5, not independent sources |

**Sources with no Δ cross-section coefficient (k n/a).** Each was read in its
relevant sections:

- HJM NPB 312 (1989) 571, preprint CTP#1608 (R1A): "A₊₋,₋₊ = (M²/ν)a₃ … The
  double helicity flip amplitude A₊₋,₋₊ falls off like 1/Q²", so it has no Δ
  term. It does fix F₁ = ⅓Σ_H A₊H,₊H.
- Close–Kumano PRD 42 (1990) 2377: the b₁ sum rule only.
- Nzar–Hoodbhoy PRD 45 (1992) 2264: "dσ ∝ Δ(x,Q²)cos2φ", with no coefficient.
- Jaffe's 1995 lectures (hep-ph/9602236) and Jaffe–Saito PLB 382 (1996) 165
  (Drell–Yan): both fetched this phase.
- CvDM17, NPLQCD17, KS20, KS20b, Kumano24, Zhao25, FDKX26 (spin 3/2), EPIOS26
  (it quotes DS16's relation and SS90's moment only), EIC Yellow Report,
  Milner–Milner 2025, Xie–Lu 2026, Xie–Chen–Lu 2026, NICA-SPD 2021, Kaur 2024,
  Cosyn–Weiss 2026.
- CRSZ25: F_UTTT = O(γ²), with no Δ.
- Fetched this phase:
  - Herwig 7 polarized DIS, arXiv:2606.23845. It cites JM only for "For a
    spin-½ hadron there is no leading-twist gluon transversity".
  - USQCD "Hadrons and Nuclei", arXiv:1904.09512.
  - Song *et al.*, PRL 136 (2026) 131901 (jet energy correlators).
  - arXiv:2610.06483, 2610.32741 (TCS), 2607.23692 and 2410.13225 (Drell–Yan).
- arXiv:2410.00251 and 2410.20950 are not spin-1 papers.

**The factor-2 traps, checked against every k above.** None of them produces
the 2:

- *Amplitude vs moment.* JM's Eq. (1) is the cos 2φ amplitude of the
  normalised φ distribution, and so is the kernel's a₂. No source above
  quotes ⟨cos 2φ⟩.
- *Pure state vs P_zz weighting.* Every printed coefficient is for a pure
  state. SS90 (A3) gives −½xλ P_zz Δ cos 2φ, so the ½ survives the ensemble
  average, exactly as the kernel's c_eff → P_zz.
- *Δ vs xΔ.* The x in JM's "x(1−y)/2 Δ" is the xy²F₁ normalisation of the
  bracket. SS90's "lowest moment" is ∫dx x Δ, and so is the letter's.
- *Per nucleon vs per nucleus.* a₂ ∝ Δ/D_φ, and the n = 2 moment is invariant
  (C2's A8).
- *Lepton-tensor normalisation.* It cancels in a₂.
- *Transverse axis ⊥ q vs ⊥ beam.* SS90's axis is ⊥ beam and JM's ⊥ q; they
  agree at leading power (SS90: "sin β is O(1/Q)").

Where the 2 does come from is §3.3.

## 2. Which convention is more common

### 2.1 The count

The sources that fix the cross-section normalisation of Δ, by kind:

| kind | k = 2 (with the ½: Jaffe–Manohar) | k = 1 (LiPolGen as coded) |
|---|---|---|
| cross section printed, own derivation | JM89 Eqs. (1), (16) and the appendix; SS90 (A7); MWZ13 Eqs. (23)/(25) | — |
| cross section printed, quoted | LOI18 Eq. (6) | LOI18 Eqs. (9)/(10); Maxwell 2019 |
| tensor printed (coefficient by contraction) | DS16 Eq. (5); LOI18 Eq. (5); EPW98 Eq. (7); PW00 Eq. (3.4) | — |
| response-function framework | PV26 (inferred) | — |
| the program's own lineage | — | Discussions.pptx p. 5, refs_dict, plans/01, fastsim, evgen, the letter's Eqs. (1)/(3), LiPolGen |

- **Independent derivations:** 3 of 3 carry the ½.
- **Documents outside the program:** 8 with the ½, counting the LOI once for
  its Eqs. (5)–(6). 2 without: LOI Eq. (9) and Maxwell 2019. Those two are
  one underived equation by one corresponding author, plus its reprint.
- **The LOI is on both sides.**

### 2.2 Weighted by what the later literature cites

Citation counts on INSPIRE, queried 2026-10-27 through its API:

| paper | citations |
|---|---|
| JM89 | 111 |
| PW00 (a review) | 225 |
| SS90 | 31 |
| DS16 | 28 |
| EPW98 | 26 |
| MWZ13 | 17 |
| LOI18 | 12 |
| PV26 | 1 |
| Maxwell 2019 | 0 |

What later papers do with them:

- **DS16 defines Δ by JM:** "We follow that reference in defining and
  outlining the construction of Δ(x,Q²) below" (Sec. II, p. 1).
- **The 2017–2026 papers that define or use Δ cite JM89, DS16 and SS90 for it.**
  Examples: CvDM17, "h1TT(x) appears in the structure function Δ(x) defined
  in [9 = JM]"; KS20 Eq. (68)'s notation dictionary; Kumano24; EPIOS26, which
  takes DS16's relation and SS90's −0.012 α_s.
- **Of the 12 records citing the LOI, only Maxwell 2019 prints a Δ
  coefficient.** Kumano24 cites the LOI only as "a project to measure the
  gluon transversity at JLab … by observing the dependence on the angle"
  (R1B). I fetched and searched the rest (§1, last paragraph), except
  Cotogno's 2018 thesis (INSPIRE 1725056), which I did not read.
- **MWZ13's 17 citing records are mostly Kumano–Song papers.** I judged
  them from their titles, and read KS20, KS20b and Kumano24 (R1B); none
  prints a Δ coefficient.

No source outside the program has adopted Eq. (9)'s normalisation, apart
from its author's proceedings.

### 2.3 Jaffe–Manohar itself, and its secondary quotations

**The preprint was read.** R1A fetched the KEK scan. R1B's later plain request
to the same URL was refused (stub sha256 `f38815fe…6149`), but neither reader
changed client identity: R1A's `fetch.py` is a bare `requests.get`, which I
read. I re-read pp. 2, 4–7 and 9 from R1A's page renders.

**Its Eqs. (1) (p. 2), (14) (p. 4) and (16) (p. 5) and the spin-J appendix
(p. 9) agree on the ½.** The appendix gives it independently:

- Eq. (15) gives "E^x_{+μ} = ½(E_{+μ} + √2E_{0μ} + E_{−μ})", so both
  double-flip density-matrix elements are ¼.
- 2(1−y) · (−½) · (¼Δe^{−2iφ} + ¼Δe^{2iφ}) = −½(1−y)Δ cos 2φ, times the
  appendix's overall x.

**The published PLB text was not seen.** It is Elsevier, and INSPIRE record
277858 lists no document. The secondary quotations of the published paper
split as follows:

- **Agree (the ½):**
  - DS16 Eq. (5): "the double-helicity-flip term in W_μν(p,q) can be expressed
    as [7]".
  - LOI18 Eqs. (5)–(6): "Following reference [1]". Eq. (6) is JM Eq. (16)
    verbatim.
  - MWZ13: its Ĝ_T coefficient function "agrees with that in [6]" (JM), and
    its Eq. (23) carries the ½.
  - SS90 (A1)/(A7): an independent derivation ("We thank R. L. Jaffe for
    advising us throughout this work").
  - EPW98 Eq. (7) and PW00 Eq. (3.4), both citing HJM and SS90.
- **Do not agree:**
  - The letter's eq:master cites [Jaffe:1989xy, Hoodbhoy:1988am] for a
    coefficient neither paper prints.
  - LiPolGen attributes the master formula, Δ term included, to HJM, which
    has no Δ term:
    - `include/lipolgen/asymmetries.hpp` line 10 at `ca9b23a`: "Hoodbhoy-Jaffe-Manohar NPB
      312:571";
    - `include/lipolgen/xsec.hpp` lines 26-28 at `ca9b23a`: "reproduces the
      Hoodbhoy-Jaffe-Manohar transcription digit for digit";
    - `docs/PHYSICS_CHANNELS.md` §2: "[HJM89]" on the master-formula row.
  - The sibling's `refs_dict.json` cites JM and HJM for the master formula
    with `"file": null`.
  - LOI Eq. (9) and Maxwell 2019 attribute their Δ term to nobody.

**The preprint's Sec. III does not match SS90's use of it.** R1A found this
and I re-derived it (§7). This bears on the bag value's derivation, not on
k (open question 2):

- Its Eq. (29) prints α_s/(8π), where its own Eqs. (20), (21) and (25) give
  α_s/(16π), which is SS90 (3.1)'s value.
- Its Eq. (31) gives a = g_x̂ − g_ŷ for E^x₊, where SS90 (2.11) has
  −2[g_x̂ − g_ŷ].

### 2.4 The answer

**No, it is not more common.** The kernel's coefficient, twice DS16's per unit
Δ, appears in the literature in one underived equation, LOI18 Eq. (9), and in
that equation's reprint. Inside the program it appears in the program's own
copies of it.

Every derivation carries the ½: JM89, SS90 and MWZ13. So does every printed
tensor: DS16, EPW98, PW00 and LOI18 Eq. (5). So does the LOI's own
Bjorken-limit formula, Eq. (6). And so do the papers the later literature
cites for Δ.

The Jaffe–Manohar convention is the common one, and as far as this search
reaches it is the only one with a derivation behind it. Keeping the kernel's
coefficient means LiPolGen's Δ is half of everyone else's.

## 3. The LOI's Eq. (6) versus Eq. (9)

### 3.1 What is read

All of this is from LOI pp. 3–5, page images, read by both readers and by me.

- **Δ is defined once.** Eq. (4): "we can denote A₊₋,₋₊ = A₋₊,₊₋ as Δ(x,Q²),
  so that the double helicity contribution to W_μν,αβ is
  W^{Δ=2}_{μν,αβ} = [P(+−,−+) + P(−+,+−)]Δ(x,Q²)". Eq. (5) is JM Eq. (14),
  with a stray printed "v" after its ½.
- **Eq. (6)** (p. 4) is quoted in §1 #4, together with its sentence on φ and
  the opposite polarization.
- **Eq. (9)** (p. 5) is quoted in §1 #5, with "a_m = ¼b_m(3cos²θ_m − 1)" and
  "usual tensor structure functions b₁ and b₂ [6]". [6] is "R. Jaffe, Private
  communication (2016)".
- **The two equations share everything except the factor:**
  - the same prefactor, since e⁴ME/(4π²Q⁴) · xy² = 4α²MExy²/Q⁴ = 2yα²/Q² with
    Q² = 2MExy;
  - the same unpolarized bracket;
  - both are cross sections, i.e. amplitudes, not moments;
  - both are for pure states, and at λ_m = ±1, θ_m = 90° Eq. (9) describes
    the same vector-polarized transverse state as Eq. (6);
  - neither carries xΔ, a per-nucleus factor or a P_zz weighting.
- **The factor itself.** Eq. (9) there gives −(1−y)/y² · Δ; Eq. (6) gives
  −½(1−y)/y² · Δ.
- **Eq. (10) inherits the factor.** Its "− 6 (1−y)/y² sin²θ_m Δ cos(2φ)" for
  N₊ + N₋ − 2N₀ compares with −3 from Eq. (6).

### 3.2 Verdict: internally inconsistent

No sentence of the LOI reconciles the two equations or redefines Δ between
them. Eq. (6) agrees with JM (16), SS90 (A7), MWZ13 (23), and with the
contraction of DS16 (5) and of the LOI's own Eq. (5). Eq. (9) is the outlier.

The same pages are loose about overall factors elsewhere (read):

- After Eq. (9), "If we average over polarization: N₊ + N₋ + N₀ ⇒ σ̄" is
  given the value of the average of the three rates, not of their sum.
- Eq. (11)'s denominator is printed "N₊ + N₋ + 2N₀".
- Eq. (7) prints Tr𝒬 where Tr𝒬² belongs; Tr𝒬 is zero for the stated
  𝒬 = diag(2/3, −1/3, −1/3).

### 3.3 The explanation: Eq. (9) is MWZ13's Eq. (23) without its ½ (inferred)

The LOI does not say where Eq. (9) comes from. It does cite MWZ13 as [2]
(p. 2: "Ma et al. [2] call Δ(x,Q²), labeled as Ĝ_T in their paper, the most
interesting of the structure functions related to tensor polarization").
Eq. (9) is MWZ13's Eqs. (23)–(24) rewritten in the LOI's symbols, with three
departures:

| | MWZ13 Eqs. (23)–(24), p. 8 | LOI18 Eq. (9), p. 5 |
|---|---|---|
| frame, state | "the initial lepton moves in the −z-direction"; "polarized along the direction m⃗ … λ_m = 0, ±1"; θ_m | "electron incident from −z"; "the target spin in the m̂ direction λ_m = (1, 0, −1)"; θ_m |
| prefactor | 2yα²/Q² | 2yα²/Q² (same) |
| unpolarized | F₁ + ((1−y)/y²)F₂^MWZ, with F₂^MWZ = F₂/x | F₁ + ((1−y)/(xy²))F₂ (same) |
| b_m | 3\|λ_m\| − 2 | 3\|λ_m\| − 2 (same) |
| rank-2 rate term | a_mF̃₁, with a_m = ¼b_m(**1 − 3cos²θ_m**) | ⅔a_mb₁, with a_m = ¼b_m(**3cos²θ_m − 1**): angular factor reversed |
| cos 2φ term | −((1−y)/(**2**y²)) b_m sin²θ_m Ĝ_T cos 2φ | −((1−y)/y²) b_m sin²θ_m Δ cos(2φ): **the 2 is missing** |

MWZ13's Ĝ_T is JM's Δ in normalisation: its tensor T⊥ is JM Eq. (14)'s, term
by term (§7). Its Eq. (23) agrees with JM (16) and with the LOI's own
Eq. (6). The factor 2 is therefore most simply a ½ lost in transcription.

On the same reading, the reversed angular factor of a_m would explain why
Eq. (9)'s b₁ sign is opposite to HJM Sec. 6 (R1A). This part is inferred.
That MWZ13's F̃₁ = +⅔b₁, i.e. carries HJM's sign, follows from its
S̃_L = (−½, 1, −½) for λ_m = (±1, 0) along ẑ (Eqs. (5), (26)); MWZ13 does not
print it.

R1B listed four other mechanisms, each of which reproduces the factor:

- normalising per unit T_xx − T_yy = 2ρ₊₁,₋₁ instead of per unit ρ₊₁,₋₁
  (KS20's "E_x − E_y instead of (E_x − E_y)/2");
- reading Eq. (6)'s term as a ⟨cos 2φ⟩ moment;
- using the m_x = 0 magnitude for m_x = ±1;
- dropping the ½ of the tensor.

None can be excluded, but nothing in the LOI's text points to them, whereas
the MWZ13 match is notational. Either way, no source defines a Δ for which
Eq. (9) is correct.

The LOI makes no numerical Δ projection (R1A: its Table 1 lists HMS rates
only), so no LOI number depends on the error. Its measurement methods,
Eqs. (10)–(12), inherit it.

## 4. What "keep the factor 2" means operationally

"Keep the factor 2" is read here as the task defines it: the kernel's cos 2φ
coefficient stays as coded. §4.4 covers the other reading.

### 4.1 The conversion

- **The relation.** Δ_LiPolGen = ½Δ_JM = ½A₊₋,₋₊, in the normalisation
  F₁ = ⅓Σ_H A₊H,₊H. Δ is per nucleon in the per-nucleon x, like F₁.
- **Into the kernel**, Δ_LiPolGen = Δ_source/k:
  - **÷ 2** from JM89, SS90, MWZ13 (Ĝ_T), DS16, LOI Eq. (6), EPW98, PW00,
    PV26, and any lattice number for Δ's moments (DS16, NPLQCD17);
  - **× 1** from LOI Eq. (9), Maxwell 2019, the letter's Eqs. (1)/(3), and the
    sibling's `fastsim` and `evgen`.
- **The conversion is linear.** Every Mellin moment converts the same way and
  the x shape does not change. Per nucleon vs per nucleus adds no factor: a₂
  goes as Δ/D_φ, and M₂ is invariant (C2 A8).
- **Relations and bounds, read in the kernel's units:**
  - **DS16 Eq. (13)** becomes Δ_LiPolGen = −(α_s/4π)Tr𝒬² x² ∫_x^1 (dy/y³) δG.
  - **SS90 Eq. (3.5)** becomes ∫xΔ_LiPolGen dx = −0.006 α_s, i.e. C_BAG/2.
    `C_BAG` itself (`include/lipolgen/constants.hpp:139-141`, "int_0^1 x Delta
    dx = C_BAG * alpha_s(Q2)") stays −0.012. That is a correct quote of SS90
    in SS90's normalisation, so the Δ in that comment is Δ_SS90.
  - **DS16 Eq. (32)'s n = 2 moment form** reads |2M₂(Δ_LiPolGen)| ≤ B.
- **Extraction.** Δ̂ from the letter's inversion, −Â y²D_φ/(1−y), is in the
  kernel's units. In JM's units Δ̂_JM = 2Δ̂.

### 4.2 The letter's moment_A / moment_B in the kernel's convention

**What is encoded (read).** `make_moment_a`/`make_moment_b` fill the kernel's
Δ slot. They are in the sibling's `delta_models.py` (sha256 `c0729707…a93ee`)
and in the tree's port (`include/lipolgen/sf.hpp:546-552`). With
`c_moment = C_BAG = −0.012` and `dilution = 1/3`:

- ∫xΔ_LiPolGen dx = **−0.004 α_s**. For moment_B this holds at every Q²; for
  moment_A it holds at q²_ref = 4.474 GeV², the letter's ⟨Q²⟩ = 4.47.
- The letter's eq:model says the same thing: "𝒩 = −0.292 fixed by
  ∫₀¹ xΔ dx = −0.012 α_s for the undiluted ansatz".
- In JM/SS90 units that is ∫xΔ_JM dx = −0.008 α_s, which is
  **2 × ⅓ × SS90's bag value**.

**What would have to change.** For moment_A/B to represent SS90's bag model,
the letter's bag numbers must be entered **halved**: c_moment = C_BAG/2 =
−0.006 (𝒩 = −0.146), giving ∫xΔ_LiPolGen dx = −0.002 α_s at dilution ⅓.
Entered **as given**, which is the status quo, they represent twice SS90's
bag model.

**What halving does to the letter's numbers.** The amplitudes are linear in
Δ. The statistical errors are rate errors and do not move. The toy column of
plans/07 is at spots x = 0.0282 / 0.0112 / 0.0708 / 0.141.

| quantity | as printed | bag at face (halved) |
|---|---|---|
| injected amplitude, toy column (×10⁻³ per unit P_zz) | 7.416 / 4.346 / 9.491 / 9.517 (letter: "7.4, 4.4, 9.5 and 9.5") | 3.708 / 2.173 / 4.745 / 4.758 |
| significance, one year | 42.9 / 30.6 / 34.5 / 21.3σ (letter: "43, 31, 35 and 21σ") | 21.5 / 15.3 / 17.3 / 10.7σ |
| δA, one year (×10⁻⁴) | 1.728 / 1.422 / 2.750 / 4.464 | unchanged |
| injected amplitude, grid column (×10⁻³) | 5.217 / 3.040 / 7.320 / 8.197 | 2.608 / 1.520 / 3.660 / 4.098 |
| significance, grid column | 27.7 / 18.4 / 27.1 / 19.4σ | 13.9 / 9.2 / 13.6 / 9.7σ |
| Δ̂ at Q² = 1.14 GeV², x = 0.02 | −0.1352 ± 0.0034 | ≈ −0.0676 ± 0.0034 (δΔ/Δ doubles) |

The letter's model-free statement does not move: "a flat amplitude of
3×10⁻³ per unit P_zz would remain a 3–11σ effect per bin even at
P_zz = 0.3" is about an amplitude, not a Δ.

**The letter is then internally inconsistent by the same factor.**

- Its Eq. (2) (eq:deltag) is DS16 Eq. (13), in JM's normalisation, while its
  Eqs. (1) and (3) are in the kernel's.
- Its introduction's "∫xΔ dx = −0.012 α_s [Sather:1990bq]" is SS90's number
  in SS90's normalisation. In the normalisation of its Eq. (1) that is
  −0.006 α_s.

Under "keep", the letter needs one of two fixes, and both are the author's:

- either ×½ on Eq. (2) together with 𝒩 = −0.146;
- or a sentence saying its scenario is twice the bag value, which is beyond
  the central gluonic Soffer bound (§4.3).

### 4.3 The t5_delta_moment_bound row: the headline V6, "13–19 % of the bound"

- **Which reading is consistent.** The harness's derived ratio
  `DS16_OVER_TREE = 2` (`validation/benchmarks/t5_delta_moment_bound.py:329`)
  is exactly k for DS16. Under "keep", **V6 is the consistent reading and A7
  is not**: the bound is written for DS16's Δ, the kernel's Δ is half of it,
  so |2M₂| ≤ B.
- **What the headline means.** The shipped toy is
  `toy_delta_gluon(x, Q², F₁, 1e-2)` in kernel units, i.e. a DS16-normalised
  Δ of 2×10⁻² · F₁ x^0.3(1−x)⁴. It uses **12.9 / 16.1 / 18.8 %** of the
  gluonic Soffer room at Q² = 2 / 5 / 10 GeV² (|M₂|/B under V6 = 0.1289 /
  0.1605 / 0.1881). It is allowed: PASS.
- **The A7 column** (6.4 / 8.0 / 9.4 %) describes a Δ half the size of the
  generated one. It remains the row's pre-registered label, and the rules
  keep it recorded, but it is not the generator's Δ.
- **The letter's scenario, as coded:**
  - moment_A uses 1.397–1.429 of the room and moment_B 1.349–1.476. That is
    beyond the central bound, and allowed only if f₁LL^g ≥ +0.815 … +0.883
    f₁^g (moment_A) or +0.714 … +0.979 f₁^g (moment_B).
  - Entered halved (bag value at face), they use 0.699–0.715 and 0.674–0.738,
    which is allowed.

  Whether the letter's scenario is admissible is therefore the same question
  as §4.2.
- **Nothing else in the row moves.** The helicity-positivity sub-row reads
  the kernel's own rates.

### 4.4 If D4 as recorded in `DECISIONS.md` stands (the Jaffe–Manohar convention adopted)

- **The coefficient** becomes −½(1−y)/y² · c_eff sin²θ_S · Δ/D_φ, so
  Δ_LiPolGen = Δ_JM.
- **k** becomes 1 for every source that carries the ½. It becomes ½ for LOI
  Eq. (9), Maxwell 2019, the sibling's `fastsim`/`evgen`, and the letter's
  Eqs. (1)/(3) until the author changes them: a Δ quoted there is half of
  LiPolGen's for the same amplitude.
- **External numbers enter unconverted.** C_BAG = −0.012 is the bag value in
  the kernel's Δ. DS16 Eqs. (13) and (32) apply as printed.
- **The t5 row.** The derived ratio becomes 1 and A7 becomes the consistent
  reading. The toy at 1e-2 is now a JM Δ and uses 6.4 / 8.0 / 9.4 % of the
  room. Its generated amplitude halves, as `DECISIONS.md` says: "every cos 2φ
  amplitude the generator produces halves".
- **moment_A/B** at C_BAG and dilution ⅓ become SS90's bag model at face,
  using 0.67–0.74 of the room. A tree regenerated with them gives half the
  letter's printed amplitudes (3.7 / 2.2 / 4.7 / 4.8 ×10⁻³, 21.5 / 15.3 /
  17.3 / 10.7σ) until the author changes the sibling and the letter.
- **Cross-repository gates.** Any gate comparing the kernel's a₂ with the
  sibling's `a_cos2phi` will differ by exactly 2 on the Δ term.
- **The letter.** Its Eq. (2) and its −0.012 quote become consistent with a
  corrected Eq. (1)/(3).

## 5. The paragraph for the generator's documentation (`docs/CONVENTIONS.md` bullet)

### 5a. Under "keep the factor 2" (decision 1 as the task states it)

- **The normalisation of Δ (the double-helicity-flip cos 2φ term).**
  LiPolGen's Δ is normalised so that a target in a pure spin state m along
  n̂(θ_S, φ_S) has the cos 2φ′ amplitude a₂ = −(1−y)/y² · c_eff(m) sin²θ_S ·
  Δ/D_φ. Here D_φ = F₁ + (1−y)F₂/(xy²), c_eff = 3Q_NN (1 for a pure m = ±1
  state, P_zz for a spin-1 fill), and Δ is per nucleon in the per-nucleon x,
  like F₁. This is the coefficient of the JLab LOI's Eq. (9) (Maxwell et al.,
  arXiv:1803.11206).
  - **Relative to Jaffe–Manohar** (PLB 223 (1989) 218, Eqs. (1), (14), (16))
    **and Detmold–Shanahan** (PRD 94 (2016) 014507, Eqs. (3)–(5)) it is
    **half**: Δ_LiPolGen = ½Δ_JM = ½A₊₋,₋₊, in the normalisation
    F₁ = ⅓Σ_H A₊H,₊H. Any Δ, moment, relation or bound quoted there enters
    the kernel ×½. DS16's Δ–δG relation, Eq. (13), reads
    −(α_s/4π)Tr𝒬² x² ∫_x^1 dy/y³ δG here, and their gluonic Soffer bound reads
    |2∫xΔ dx| ≤ B (`t5_delta_moment_bound`'s V6).
  - **Relative to the LOI's Eq. (6)** it is half. The LOI's Eqs. (6) and (9)
    differ by this factor, and Eq. (6) is the one that agrees with
    Jaffe–Manohar.
  - **Relative to Sather–Schmidt** (PRD 42 (1990) 1424, Eq. (A7)) it is half.
    `C_BAG` = −0.012 is their bag value of ∫xΔ dx/α_s in their normalisation,
    −0.006 in LiPolGen's, so `make_moment_a/b` with `c_moment = C_BAG` encode
    twice the bag model.
  - **Ma–Wang–Zhang's Ĝ_T** (arXiv:1306.6693, Eq. (23)) is Jaffe–Manohar's Δ.

  The choice relabels Δ and moves no physics: a given amplitude is the same,
  only the Δ quoted for it differs by 2
  (`docs/open_items/run_2026-10-27/phase_R1_delta_convention.md`).

### 5b. If D4 as recorded stands (the Jaffe–Manohar convention)

- **The normalisation of Δ (the double-helicity-flip cos 2φ term).**
  LiPolGen's Δ is the Jaffe–Manohar structure function: Δ = A₊₋,₋₊ in the
  normalisation F₁ = ⅓Σ_H A₊H,₊H (Jaffe–Manohar, PLB 223 (1989) 218,
  Eqs. (10)–(11), (14), (16)), per nucleon in the per-nucleon x, like F₁. A
  target in a pure spin state m along n̂(θ_S, φ_S) therefore has the cos 2φ′
  amplitude a₂ = −½(1−y)/y² · c_eff(m) sin²θ_S · Δ/D_φ, with
  D_φ = F₁ + (1−y)F₂/(xy²) and c_eff = 3Q_NN (1 for a pure m = ±1 state, P_zz
  for a spin-1 fill).
  - **Relative to Jaffe–Manohar and Detmold–Shanahan** (PRD 94 (2016) 014507,
    Eqs. (3)–(5), (13), (32)) it is identical: their relation and bound apply
    as printed.
  - **Relative to the LOI** (Maxwell et al., arXiv:1803.11206) it is
    identical to Eq. (6). The LOI's Eq. (9) carries twice the coefficient, so
    a Δ quoted in that normalisation is half of LiPolGen's for the same
    amplitude. The same doubled coefficient appears in PoS(SPIN2018)174, the
    sibling program's slides, `polli_fastsim.asymmetries.a_cos2phi`,
    `polligen.xsec`, the letter's Eqs. (1) and (3), and LiPolGen before
    2026-10-27.
  - **Relative to Sather–Schmidt** (PRD 42 (1990) 1424, Eq. (A7)) it is
    identical: `C_BAG` = −0.012 is their bag value of ∫xΔ dx/α_s in
    LiPolGen's Δ.
  - **Ma–Wang–Zhang's Ĝ_T** (arXiv:1306.6693, Eq. (23)) is the same function.

  (`docs/open_items/run_2026-10-27/phase_R1_delta_convention.md`)

## 6. Open questions, each with what would settle it

1. **The published PLB 223 text.** Does the journal version keep the
   preprint's ½ in Eqs. (1), (14) and (16)? Yes is expected: DS16 (5) and LOI
   (5)–(6) reprint the preprint's forms while citing the journal, and MWZ13
   re-derives the cross section. *Settles it:* the journal pages (Elsevier,
   not reachable from here), e.g. a library copy.
2. **The preprint's Sec. III against SS90.** Eq. (29)'s α_s/(8π) disagrees
   with the α_s/(16π) of its own Eqs. (20), (21), (25) and of SS90 (3.1).
   Eq. (31)'s a = g_x̂ − g_ŷ (circular E^x₊) disagrees with SS90 (2.11)'s
   −2[g_x̂ − g_ŷ]. This does not touch k. It bears on:
   - the bag value's derivation: SS90's arithmetic with 1/(16π) gives
     −0.0115 α_s, i.e. −0.012; the preprint's Eq. (29) would give −0.023
     α_s;
   - the parton-model normalisation behind t5's B, which follows DS16/SS90 in
     C2 S1's linear-state reading.

   *Settles it:* the journal text, or a fresh calculation of the n = 2
   matrix element.
3. **The sign of Δ for a given gluon linear polarization.** The sources
   differ: JM (2) is +; DS16 (13) is −, with its fn. 1 "differs by a sign
   from that in Ref. [7]"; SS90 (2.11) has −2; KS20 has Eq. (68). This is not
   k. It decides whether the letter's moment_A sign, which is opposite to
   the shipped toy's (phase A §4.3), is SS90's. *Settles it:* a box-graph
   calculation in one fixed state convention.
4. **Where Eq. (9) came from**: an MWZ13 transcription, or the private
   communication [6]? It is not needed for the convention. *Settles it:*
   only the LOI's corresponding author.
5. **b₁ and Δ double counting.** SS90 fn. 8 (p. 1428): "In the Bjorken limit
   they are related by b₁ = Δ₁ + ½Δ₄ and b₂ = Δ₂ + xΔ₄". PV26 (A9a):
   "⅓b₂ − b₃ = νΔ". Once the double flip is leading twist, HJM's covariant
   b's contain it. A kernel that adds a separate Δ term must therefore be
   fed the helicity-diagonal b₁: CK90 Eq. (9b), b₁ = A₊₀,₊₀ − ½(A₊₊,₊₊ +
   A₊₋,₊₋), the A_zz-measured one. This matters because SS90 (A7)'s
   φ-independent tensor term at leading power carries Δ₁ and Δ₂ only. In the
   massless kernel b₁ enters only the φ-independent rate and Δ enters only
   a₂. The finite-γ path's cos 2φ harmonic reads b₁…b₄ (its O(γ²) leakage).
   *Settles it:* stating that the shipped b₁ models (Miller's convolution,
   `toy_b1`) are helicity-diagonal. This is a physics-lens item for the D4
   plan.
6. **Finite Q².** The kernel carries the Bjorken-limit coefficient. The exact
   contraction of JM (14) differs from it:
   - by −0.7 % at x = 0.1, y = 0.5, Q² = 5 GeV² (ratio 0.4964, not ½);
   - by 0.04–0.6 % at the four sweet spots (R1B: 0.49981 / 0.49879 /
     0.49718 / 0.49756).

   This is a modelling choice for D4, not a convention question. *Settles
   it:* a decision to implement JM (14) exactly.
7. **Spin 3/2.** JM p. 1 says the analysis "applies without modification to
   a spin-3/2 target as well". I did not check the normalisation of
   `delta_32_func` (c_eff = 3Q_NN at J = 3/2) against JM's spin-J appendix,
   and no spin-3/2 slot is filled by default. *Settles it:* contracting the
   appendix formula with J = 3/2 density matrices. The geometry is inferred,
   not computed: an unpolarized lepton sees only the symmetric W_μν, which
   couples to the even-rank (0, 2) part of ρ, so the J = 3/2 double-flip
   term is rank 2 and its m-dependence goes as Q_NN, as c_eff assumes. Its
   normalisation relative to J = 1 is what was not checked.
8. **PV26's identification** (medium confidence): their
   Δ = W₊₊,₋₋ = (1/ν)(b₂/3 − b₃), read as JM's A₊₋,₋₊ through HJM's Eq. (7).
   It is not needed for the count. *Settles it:* an explicit contraction of
   their responses, or their authors.

## 7. Provenance: what I re-read, ran and fetched

**Page images I re-read:**

- JM preprint pp. 2, 4, 5, 6, 7, 9: R1A's renders (`r1a/jmhi/`) of the PDF,
  whose sha256 I checked;
- SS90 p. 1427 and Eq. (3.1): `r1a/ss-4.png`, `r1a/ss_eq31_crop.png`; PDF
  sha256 checked;
- LOI pp. 4–5: `r1a/loi-04.png`, `loi-05.png`; PDF sha256 checked;
- Maxwell 2019: `r1b/maxwell_p3_eq.png`;
- PV26 p. 10: `r1b/pv26_p-10.png` (App. A from text);
- MWZ13 pp. 3, 4, 8, EPW98 p. 4 and PW00 p. 28, rendered here.

All of these are under the session scratchpad
`/tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/`.

**My own contraction.** `r1j/judge_contract.py` (output
`r1j/judge_contract.out`) uses numpy only and the printed equations only. The
setup:

- target rest frame, q along +ẑ (JM Sec. II), φ measured from the alignment
  axis x̂;
- JM Eq. (14) / DS16 Eq. (5) with E = E′ = E^x₊ = −(0,0,1,i)/√2, together with
  −F₁g̃ + (F₂/ν)p̃p̃;
- contracted with L^{μν} = 2(kk′ + k′k − g k·k′).

Results at x = 0.1, y = 0.5, F₂ = 2xF₁, as a₂(JM tensor)/a₂(kernel):

| Q² (GeV²) | ratio |
|---|---|
| 5 | 0.496416 |
| 100 | 0.499820 |
| 10⁴ | 0.499999 |
| 10⁶ | 0.500529 (float64 cancellation) |

m_x = −1 gives the same as m_x = +1. m_x = 0 (E = x̂) gives −2 × the m_x = ±1
amplitude. This agrees with C2, with R1A (`r1a/ds16_contract_mp.out`, 60
digits) and with R1B (`r1b_scripts/ds16_contract.py`).

**JM (14) = MWZ13's T⊥ Ĝ_T, term by term (derived).** Take p·E = 0, with
p̃ = p − (ν/q²)q, p̃² = κν²/Q² and g⊥ = g − qq/q² − p̃p̃/p̃². Then:

- JM's bracket [g − qq/q² + (q²/κν²)p̃p̃] = g⊥;
- E′*·E + (M²/κν²)(q·E′*)(q·E) = E⊥′*·E⊥;
- E − (q·E/κν)(p − (M²/ν)q) = E⊥.

So JM (14) = ½[E⊥′*_μE⊥_ν + E⊥′*_νE⊥_μ − g⊥_μν E⊥′*·E⊥]Δ = T⊥_μν Δ.
SS90 (A1)'s Δ₄ term (κΔ₄ = Δ), EPW98 (7) and PW00 (3.4) are the same tensor
with the q-proportional terms dropped.

**The JM (29) check (derived; agrees with R1A's `sumrule_n2_check`).**

1. Contract JM Eq. (21) at n = 2 with μ₁ = μ₂ = n (n·p = p⁺) and μ, ν
   transverse.
2. With Eq. (20)'s O = ½S{GG}, the left side is ⟨O_yynn − O_xxnn⟩ =
   ½⟨G_yn² − G_xn²⟩.
3. For E^x₊ the right side is ½(p⁺)²A₂. The transverse part of E^x₊ is
   −ŷ/√2; its ẑ part enters only at O(M/p⁺).
4. So A₂ = ⟨G_yn² − G_xn²⟩/(p⁺)² = ½⟨(E₂ ± B₁)² − (E₁ ∓ B₂)²⟩/(p⁺)², the
   signs depending on the light-cone direction.
5. Eq. (25) then gives ∫xΔ dx = (α_s/8π)Tr𝒬² A₂ = (α_s/16π)Tr𝒬²⟨…⟩/(p⁺)²,
   which is SS90 (3.1)'s coefficient and field combination. The preprint's
   Eq. (29) prints α_s/(8π).

**Fetched this phase.** Plain `python3 -I` with `requests`, no header
changes. Files are in `scratchpad/refs/r1j_fetch/`; the sha256 values are
full; none was put into a repository.

| file | sha256 |
|---|---|
| arXiv:1306.6693v2 (MWZ13) | `7404d85641c542f569959d25057eb37c1951810b2ee986f964593c90e7342bc7` |
| hep-ph/9709455 (EPW98) | `ebcdf53aa8e87c928794cb2bacbfb0a7773076392adb8b68e03d9ed4baa14f54` |
| hep-ph/9908230 (PW00) | `585c6746964ccedd3ddc30cfee4c114fc378a30a7083c6654d4ad81dd808eb84` |
| hep-ph/9602236 (Jaffe lectures) | `967fd52a9f9e64b4ec308037b73a9c44f68c168da23e988b61434b60c9df35ca` |
| hep-ph/9604220 (Jaffe–Saito) | `6200b5c8d60073032562e8388b365ca449fb25a58c13587ec7c1ce0f9cece4ee` |
| arXiv:2606.23845 (Herwig 7 polarized DIS) | `06a410f15097bfe1bcc2a339e158c812bc845cddf02e43f50eb20aab81275483` |
| arXiv:1904.09512 | `51ac7d0f4fa20947399b9ed0f2688ba0c292b34c9339cbcef1cddf8bc5b1fdf8` |
| arXiv:2509.14960 | `b496cf25f5b2f9db6e24a1852cc7e61bc177e25457bb3771ff339dc5d036068c` |
| arXiv:2610.06483 | `666ff387178ce953329d350fd52110d9618d5e70051debcf13a4ec91d04629b0` |
| arXiv:2610.32741 | `dad9703246e934fbe8aa76e06ea62038856618067b7e7a19e6267556fc6e3cac` |
| arXiv:2607.23692 | `89e2e5fe2c0105143c188654dff90a134249b8fcec8af58f35a4f9a3feaadf13` |
| arXiv:2410.13225 | `ced4afa846b62efd75c58891867b012f958c06b33f2779a99b96697e3b68182c` |

The readers' downloads (JM and HJM preprint scans, the DS16 erratum, PV26,
Maxwell 2019 and the rest) are listed with their sha256 in their findings.
They too are in the scratchpad, which is not permanent; filing any of them in
`PolarizedLithiumSim/refs/` is the author's call.

**Not done:**

- no repository file other than this record was edited;
- no test or harness was run against the tree;
- the sibling repository was only read, through the shell.
