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
common one, and it is not. Jaffe–Manohar PLB 223 (1989) 218 is read as its
MIT preprint CTP#1712, Eq. (16); in general a₂ = −½ (1−y)/y² c_eff(m)
sin²θ_S Δ/D_φ, of which the bold text is the m = ±1 case. Consequences:
every **Δ-induced** cos 2φ amplitude the generator produces halves —
everything Δ-free, the θ_S = 0 output among it, is bit-identical, and the
tensor-γ leakage's Δ-equivalent doubles; external Δ numbers — the bag value
−0.012 α_s, the lattice bound — enter unconverted; the
t5_delta_moment_bound row's first-pass identification A7 becomes the
consistent one; the attribution of the Δ term to Hoodbhoy–Jaffe–Manohar
(which has none) is corrected and the b-sector reading documented (no
double counting on the default path); the port gate reads polligen's
coefficient off `xsec.json` against a literal ½. `moment_A` / `moment_B`
represent the bag model at face value in this tree from D4 on; the
sibling's `asymmetries.py`, `evgen/polligen/xsec.py` and the letter's Eqs.
(1) and (3) carry the old coefficient and are the author's to change —
until they are, the letter's printed amplitudes encode 2 × C_BAG. The
implementation follows `plan_D4_jm_convention.md`; corrections (1)–(6) of
`phase_I4_jm_convention.md` §8 item 8 applied by the fix stage,
2026-11-03.)*

**D5 — `RcTailModel::PolradFull` becomes the shipped `rc_tail` default;
`TPeak` stays an option.** *(decision 2; the implementation follows
`plan_D5_polradfull_default.md`)*

**D6 — the VMC AV18+UX cluster waves (`--cluster-wave vmc`) become the
shipped default; the Hulthén forms stay an option. Supersedes D2 of
2026-09-27.** *(decision 3; implemented 2026-10-29 as the run-surface
default per `plan_D6_cluster_wave_default.md` and `phase_I6_cluster_wave.md`:
the channel factories and the polligen port gate keep the analytic forms.
Evidence beyond the deuteron's η_d: the Hulthén α–d pair's asymptotic D/S
ratio has the wrong sign, +0.070 against the VMC wave's −0.048 and the
measured −0.025 ± 0.012 — the VMC wave itself sits about 1.9 σ from that
measurement; the "well accepted" premise is the author's and its usage was
not surveyed; CD-Bonn, nearer the measurement, is the recorded alternative
for the deuteron only (the tree has no CD-Bonn α–d or α–t wave;
`plan_D6_cluster_wave_default.md`). Still the author's: whether `LI6_CLUSTER_POLARIZATION`
0.811228 and `LI6_B1_RANK2_TRANSFER` 0.921947, which the default inclusive
⁶Li run still reads, are re-derived on the VMC wave (0.887 and 0.983); they
were kept and the inconsistency recorded.)*

**D7 — the shipped structure-function defaults become the grid sets the
letter's production uses (free nucleon, without its EPPS21 ratio):
`--unpol-sf ct18nlo` (CT18NLO through LHAPDF)
and `--pol-sf nnpdfpol` (NNPDFpol1.1); `ToyF2` and `ToyG1` stay as opt-in
options; the leading-order evaluation of the NLO sets (no coefficient
functions, massless charm and bottom) is filed as the development item.**
*(decision 4, on the supervisor's recommendation; implemented 2026-10-29
per `plan_D7_backend_defaults.md` and `phase_I7_backends.md`. Two choices
made by the supervisor because the default must not depend on the machine
(in the design's check, an `lhapdf.conf` saying `Extrapolator: nearest`
instead of this machine's `continuation` moved the pre-freeze ⁶Li cross
section by +5.03 %; `plan_D7_backend_defaults.md`, underrated risk 1 — the
design had proposed continuation, its check proposed the freeze): the sets
are FROZEN at their grid edges in code (LHAPDF's nearest extrapolator set on
every PDF object the tier builds, plus a clamp of x and Q²) — the letter's
`--pdf grid` chain freezes Q² only, at Q₀, and on a nuclear set — and the
below-grid share of the rate is
recorded in meta (about 39 % of the default ⁶Li inclusive rate sits below
CT18NLO's Q² floor of 1.68 GeV²). Still the author's: the EPPS21 ⁶Li ratio by
default for lithium beams, which no pipeline path applies today, and R1998
in place of `r_sigma_lt` (registry row 3); the new default is free-nucleon
CT18NLO with the toy R, neither the letter's grid chain nor its published
toy footing. The development item carries a measured cost the run's records
first described too narrowly (verification of 2026-11-03): at LO the shipped
pair breaks |A₁| ≤ 1 from x ≳ 0.83 on the free neutron and from x ≳ 0.87 on
the deuteron and the proton (A₁ⁿ 1.22 at x = 0.85 and 1.40–1.45 at 0.871,
Q² = 10–1000; A₁ᵖ 1.01 at 0.871, Q² = 10; the shipped deuteron kernel's A₁ᵈ
1.15–1.22 at 0.90); `--x-max 0.95` prevents the positivity REFUSAL of the
top cell but not the violation — the x = 0.871 cell still generates with
|A₁| > 1 on `tagged-d-p` and the inclusive deuteron
(`phase_I7_backends.md` §1, §3.5).)*
