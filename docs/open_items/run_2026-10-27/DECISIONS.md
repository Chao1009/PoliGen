<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Decisions of 2026-10-27

Taken by the author on the priced options of run 2026-10-26
(`../run_2026-10-26/DECISIONS_PRICED.md`, `../run_2026-10-26/SUMMARY.md`
"For the author"), in these words:

> 1. we keep the factor 2, is this more common? 2. PolradFull is wanted.
> 3. if that's the well accepted wave default, then do it, otherwise if there
> is a better model favored by the measured D/S ratio we can use it. 4. I
> don't quite understand the question.

and, on the fourth once restated (the shipped structure-function backends:
the toys, or the grid backends the letter's production already uses, with
the leading-order evaluation of the NLO sets filed as the development item):

> Your recommendation sounds the best approach, go with it.

and, on the first once the literature had been read (the kernel's cos 2φ
coefficient without the ½ appears only in the LOI's underived Eq. (9) and a
proceedings reprint of it; Jaffe–Manohar, Sather–Schmidt, Detmold–Shanahan
and the LOI's own Eq. (6) all carry the ½) and asked "Which convention is the
most common ones?":

> Then we go with the Jaffe-Manohar convention (your recommendation)

Each entry below records the decision, the evidence it rests on, what it
changes in the tree and what it does not. The records of this run
(`phase_R1_delta_convention.md`, `plan_D5_polradfull_default.md`,
`plan_D6_cluster_wave_default.md`, the implementation and verification
records) carry the detail; `STATUS.md` the board.

**D4 — the kernel adopts the Jaffe–Manohar convention for Δ: the cos 2φ
coefficient becomes −½ (1−y)/y² · c_eff sin²θ_S · Δ/D_φ for a pure
transverse m = ±1 state, so that Δ(x, Q²) in LiPolGen is the double-helicity-
flip amplitude A₊₋,₋₊ normalised like F₁, as in Jaffe–Manohar PLB 223 (1989)
218, Sather–Schmidt PRD 42 (1990) 1424, Detmold–Shanahan PRD 94 (2016)
014507 and the LOI's Eq. (6).** *(decision 1, revised by the author on the
literature reading of `phase_R1_delta_convention.md`: the first answer,
"we keep the factor 2", was conditional on the kernel's convention being the
common one, and it is not. Consequences: every cos 2φ amplitude the
generator produces halves; external Δ numbers — the bag value −0.012 α_s,
the lattice bound — enter unconverted; the t5_delta_moment_bound row's
first-pass identification A7 becomes the consistent one; the sibling's
`asymmetries.py`, `evgen/polligen/xsec.py` and the letter's Eqs. (1) and
(3) carry the old coefficient and are the author's to change, after which
`moment_A` / `moment_B` represent the bag model at face value. The
implementation follows `plan_D4_jm_convention.md`.)*

**D5 — `RcTailModel::PolradFull` becomes the shipped `rc_tail` default;
`TPeak` stays an option.** *(decision 2; the implementation follows
`plan_D5_polradfull_default.md`)*

**D6 — the VMC AV18+UX cluster waves (`--cluster-wave vmc`) become the
shipped default; the Hulthén forms stay an option. Supersedes D2 of
2026-09-27.** *(decision 3; the implementation follows
`plan_D6_cluster_wave_default.md`)*

**D7 — the shipped structure-function defaults become the grid backends
the letter's production uses: `--unpol-sf ct18nlo` (CT18NLO through LHAPDF)
and `--pol-sf nnpdfpol` (NNPDFpol1.1); `ToyF2` and `ToyG1` stay as opt-in
options; the leading-order evaluation of the NLO sets (no coefficient
functions, massless charm and bottom) is filed as the development item.**
*(decision 4, on the supervisor's recommendation; the implementation
follows `plan_D7_backend_defaults.md`)*
