<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Decisions of 2026-11-03

Taken on the 17 open items of run 2026-10-27 (`../run_2026-10-27/STATUS.md`
"For the author"; `../run_2026-10-27/SUMMARY.md`). The author's words, in
full:

> Go with your recommendations for the 17 items, just record the choice and
> why.

The author delegated the choice on each item to the supervising session's
recommendation, asking for the choice and its reason to be recorded. Every
entry below is therefore the SUPERVISOR's choice under the author's
delegation — to be labelled so wherever it is quoted — not the author's own
words. The numbering follows run 2026-10-27's list (A1 is its item 1, …).
Each entry gives the choice (bold), why, what changes in the tree and what
does not. The phases that implement them, their records and gate tallies are
on the board, `STATUS.md`.

Ground rules carried over from the two previous runs: a default moves only
under an entry below; every number that depends on a moved default is
re-taken with a dated comment naming the entry (`A<n>, 2026-11-03`), or the
test passes the old value explicitly where its purpose is to pin the old
physics; tolerances are declared before verdicts; a benchmark that disagrees
is recorded, never tuned; `validation/reference/*.json` are not edited, with
one ADDITIVE exception (A13's dumped block); `docs/HEPMC3_CONVENTION.md` is
edited by A10 only, with the T6 row re-pinned; the hashed design blocks of
`t4_djangoh_rad_noRad.py` and `t2_boeglin_nk.py` are never edited; the
sibling repository `PolarizedLithiumSim` is read-only; nothing is published
or sent.

**A1 — the inclusive ⁶Li polarization factors follow the best-founded
calculation of each, within the ANL AV18+UX family the default wave belongs
to: `LI6_CLUSTER_POLARIZATION` 0.811228 → the six-body ab initio value 0.848
(R. B. Wiringa et al., PRC 89 (2014) 024305, Table I; the tree's
`LI6_POLARIZATION_VMC_SIX_BODY`); `LI6_B1_RANK2_TRANSFER` 0.921947 → the VMC
α–d value 0.982576 (the default wave's tensor dilution); the inclusive
deuteron species' `P_D_DEUTERON` 0.045 → AV18's 0.057600
(`deuteron_av18_p_d`; `DEUTERON_VECTOR_POLARIZATION` 0.9325 → 0.9136). The
Hulthén and scenario values stay named constants, reachable with
`--cluster-wave hulthen`, and the polligen port gates pass them explicitly,
so no reference JSON moves.** *Why: (a) since D6 the two default ⁶Li runs
describe the nucleus with different wave functions, +11.6 % in the vector
sector and +6.6 % in the rank-2 one (`../run_2026-10-27/phase_I6_cluster_wave.md`
§9 item 1). (b) The cluster product (1 − 1.5 P_D)(1 − 1.5 P_D) spans
0.811 … 0.905 depending on the wave it is fed, and the formula carries the
error (`include/lipolgen/beams.hpp`, C5.5); for the vector sector one ab
initio number for this very observable exists, 0.848, which the author's
decision of 2026-09-04 already named the reference point — 0.811 was kept
then as the nearer of the two cluster values to it, and the delegation now
allows the anchor itself. (c) For the rank-2 transfer no six-body number
exists, so the cluster value on the default wave is the best founded; the
Hulthén value is computed on an α–d wave whose asymptotic D/S ratio has the
wrong sign (D6's evidence); the non-cluster depletion that takes the vector
sector from the cluster 0.887 to 0.848 (4.4 %) is the unquantified
systematic of the rank-2 value and is stated beside it. (d) The deuteron:
C5.5b already made the embedded deuteron of a default tagged run follow AV18
("one Hamiltonian end to end"); the inclusive species' 0.045 is the
sibling's scenario value, chosen to reproduce the 0.87 of
`b1_li6_from_deuteron`, labelled "AV18-like" where AV18's is 0.0576; the
CD-Bonn–AV18 spread 0.0485–0.0576 and the experiments' ω_D = 0.05 ± 0.01
bound its band. What changes: every inclusive ⁶Li polarized observable in
the vector sector ×1.0453 (0.848 / 0.811228), b₁ and the inclusive A_zz
×1.0658; the inclusive deuteron's polarized observables ×0.9797; every pin
re-taken with a dated A1 comment; the kernels take the factors as parameters
so the port gates pass polligen's values explicitly (`beams.json`,
`xsec.json`, `b1_default_li6.json` untouched); the band around 0.848
(0.811 … 0.887, the two cluster readings) stays mandatory wherever an
inclusive ⁶Li polarization is quoted; the report regenerated; registry row
10 closed by this entry. What does not: the tagged channel (on the VMC wave
since D6), the Hulthén option, the sibling.*

*Supervisor's note, 2026-11-03, after phase P1 (its record
`phase_P1_polarization_factors.md` §5) — two consequences the entry above
does not name. (1) The default inclusive ⁶Li `--plan helicity-flip` and
`--plan apar` runs (P_z = P_e = 0.7, spin configuration 1) are now refused by
the sampler's positivity check at the top cell of the shipped window
(x = 0.955, Q² = 1897 GeV²): 1 + w_avg = −0.0185 (vector term −0.2578, rank-2
term −0.7607), where the pre-A1 factors gave +0.0396. The choice: the refusal
stands. It is the rank-2 (b₁) term of the toy tensor model at the top cell,
the same term behind the refused deuteron cases of
`../run_2026-10-27/phase_I7_backends.md` §3.5, and A4 keeps the rank-2
refusals; A3's clip acts on g₁ (A₁ = 0.616 at that cell) and will not lift
it; `x_max` stays 1.0 (A4). The tests whose subject is not that cell run that
one configuration at x_max = 0.95 through helpers carrying a dated reason;
USAGE and the CLI say so; P2 re-measures it beside the deuteron cases.
(2) The inclusive deuteron's b₁ / A_zz moved ×1.0658 as well, because the
class-level toy kernel fills the bare deuteron's b₁ with the ⁶Li transfer
(`Li6B1(MillerB1)` for every spin-1 ion). Pre-existing, recorded, not a
choice of this run; it goes to the author with the summary.*

*Correction to the note above, 2026-11-06, after phase P2 (its record
`phase_P2_positivity.md`): the clip DID lift that refusal. A₁ = 0.616 is
the nucleus's value at that cell; the nucleons sit at A₁ⁿ = 4.4 and
A₁ᵖ = 1.4 there, both clipped to 1 under A3, and the default inclusive ⁶Li
`--plan helicity-flip` / `--plan apar` run builds at x_max = 1.0
(1 + w_avg = +0.1210; vector −0.2578 → −0.1183, rank-2 −0.7607 unmoved).
Consequences: the x_max = 0.95 test helpers of P1 lose their reason and
are removed in the P2 commit (pins re-taken with `A3 (2026-11-03)`
comments, never loosened); the USAGE/CLI sentences say the run builds
since A3. What A4 keeps: the two inclusive deuteron refusals (spin
configurations 0 and 1 at the default), each a rank-2 one on the
class-level toy b₁ — the bare deuteron reading the ⁶Li transfer, item (2)
above, now the author's most consequential open item of this run.*

**A2 — EPPS21 and R1998 stay options; the shipped default stays the
free-nucleon CT18NLO with the toy R (`r_sigma_lt`).** *Why: EPPS21's Li6
ratio parametrizes the inclusive nuclear modification — shadowing, EMC,
Fermi smearing — fitted on heavier nuclei and interpolated in A, while the
generator builds its lithium cross section from an explicit cluster
convolution that already carries Fermi motion and the cluster structure, so
the ratio on top double-counts the Fermi region; no pipeline path applies it
today (`../run_2026-10-27/phase_I7_backends.md` §13 item 7), and a design
would first have to say which part of the ratio the convolution does not
already contain. R: `r_sigma_lt` against `r1998` is up to 19 % in (1 + R)
and 38 % of the shipped window lies outside R1998's own support
(`include/lipolgen/pipeline.hpp`'s note); both pass NMC's R (12.54/10
against 6.20/10, report row 12), and Whitlow and E140 are R1998's own
inputs, so the benchmark cannot separate them where the window needs it;
moving a 10–40 % lever without a measurement that favours it is what the
runs' rules forbid. The deciding measurement, recorded as the item: R below
Q² = 1 GeV² at x < 0.05, which no vendored table carries. What changes:
nothing in code; registry row 3 carries this entry; `--r-source r1998` and
the harnesses' EPPS21 configuration stay opt-in.*

**A3 — the positivity bound |g₁| ≤ F₁ is enforced where the nucleon
structure functions are evaluated, for every backend: g₁ is clipped to ±F₁
per nucleon and species, the clipped share of the rate is recorded in `meta`
and printed on the run banner; the top-cell positivity refusal stays for
what it still catches (the rank-2 term); the NLO coefficient functions stay
the filed development item.** *Why: a generator must not produce |A₁| > 1 —
one helicity state's cross section goes negative. At LO the shipped pair
breaks the bound from x ≳ 0.83 on the neutron and 0.87 on the deuteron and
the proton (`../run_2026-10-27/phase_I7_backends.md` §1), a mismatch between
CT18NLO's high-x quarks and the NNPDF2.3 reference NNPDFpol1.1's positivity
was imposed against, not a prediction. The bound is exact and
convention-free, the clip is visible (meta, banner) and must be bit-identical
wherever it does not bind — the toys, to be measured. The NLO coefficient
functions are a design of their own (F_L, charm, bottom) and would not by
themselves restore the positivity of a mixed pair. What changes: the
x ≳ 0.83 cells of the default run; the six top-cell vector refusals of I7
§3.5; the t5 positivity pins, re-taken.*

**A4 — follows from A3; `x_max` stays 1.0.** *Why: the two CLI plans
(`--isotope d --plan helicity-flip`, `--channel tagged-d-p --plan
helicity-flip`) are refused by the vector term at the top cell; A3 removes
the violation, where a smaller `x_max` would only hide it. Any refusal that
survives A3 (a rank-2 one) is the right answer and stays. Measured in P2's
record.*

**A5 — a report row's headline follows the library's shipped default
configuration; the configuration a harness's declaration names is a recorded
sub-row. Written once as a convention of `validation/benchmarks/README.md`;
the Deeps row re-pointed to the wave the default d–p run uses.** *Why: the
report's one question is what the generator SHIPS does against data; a
headline pinned to an opt-in configuration answers a question about an
option; both numbers stay visible either way; the two opposite readings of
one clause (`../run_2026-10-27/phase_I6_cluster_wave.md` §9 item 12,
`../run_2026-10-27/phase_I7_backends.md` §6) cannot recur once the rule is
written. What changes: the Deeps headline (FAIL on both waves, so the report
stays 14 / 9), the convention text, the README rows; the T1 rows already
follow it.*

*Supervisor's note, 2026-11-06, after phase P3's survey of every harness
against the rule (`phase_P3_report_rules.md` §1.2, §5) — the four rows the
entry did not foresee, each the supervisor's choice under the same
delegation. (1) `t2_bonus_spectator_shape`: the row's object is the opt-in
FSI weight; on the shipped `--fsi off` the weight is 1 and every block
passes trivially (12/12, 16/16), which tests nothing. The rule is refined,
not excepted: a row whose object is an opt-in component builds its headline
with that component ON in its shipped form and everything else at the
shipped default — here AV18 with `GlauberFsiWeight` at the shipped 40 mb
ceiling (3/12 and 0/16, FAIL → FAIL); the declared Hulthén-at-40 mb numbers
are the recorded sub-rows; the sub-rows built on the harness's default
configuration follow the headline, the declared ones stay recorded. The
refinement is written into convention 8. (2) `t2_deuteron_static`: its
declared verdict is a conjunction over the three shipped waves; the rule
applied as written puts the headline on AV18 (η −1.39σ, PASS) with CD-Bonn
(−0.07σ) and Hulthén (+17.2σ, the deficiency that decided D6) as recorded
sub-rows, and the declared all-waves verdict (FAIL) recorded beside it.
Status FAIL → PASS by the rule, not by physics; the row's text says so. The
choice is uniformity: the rule was written on 2026-11-03 before this
measurement and admits no exception by row; a headline FAIL driven by an
option misstates what ships. (3) `t2_hermes_b1_table2`: the headline moves
from `cdks_conv_mstw` to the shipped b₁ (`miller_shipped`), and, because the
rule names the whole shipped configuration, its A_zz numbers are built on
the shipped F₁ᵈ backend and R as `make_config()` builds them; the declared
`cdks_conv_mstw` on MSTW / R1990 stays a recorded sub-row with its A_zz
numbers unchanged; every moved number before → after in the record; status
FAIL → FAIL (χ²(b₁)/6 22.26 → 5.26). (4) `t2_boeglin_nk`: its headline wave
is named inside the hashed design block, which this run never edits; the
headline stays Hulthén, the AV18 value (11.54/21, p = 0.951, PASS by the
same rule) is recorded as a sub-row computed outside the block where the
harness allows it, otherwise in the README row, and convention 8 names
this one exception. (5) Deeps' Glauber and transverse-window sub-rows carry
no status and stay on the declared wave as records. Report tally after
(1)–(4): 15 pass / 8 fail, the one move labelled as the rule's.*

**A6 — V6 leaves the TIGHT edge of `t5_delta_moment_bound`, recorded as the
pre-D4 reading.** *Why: V6 is "every M₂ in DS16's normalisation is 2 M₂",
i.e. the pre-D4 kernel's DS16_OVER_TREE = 2; since D4 that factor is 1, so
V6 is the central value under another name and the TIGHT edge carries a
vacuous member; a verdict edge should state the convention the tree ships.
What changes: the edge's definition (its value unchanged, the factor being
1), the row's S3 text, a dated edit; no status moves (|M₂|/B 0.076–0.102).*

*Supervisor's amendment, 2026-11-06, from phase P3's measurement
(`phase_P3_report_rules.md` §1.3): the sentence above, "since D4 that factor
is 1 … the edge's value unchanged", does not describe the harness. The TIGHT
edge did not carry V6 through `DS16_OVER_TREE` (1 since D4) but through
`PRE_D4_OVER_JM = 2.0`, a ×2 on M₂ (the pre-D4 kernel / the letter as
printed), so removing V6 HALVES every TIGHT saturation: the verdict model
`toy_1e-2[shipped]` 0.3207 / 0.3808 / 0.4281 → 0.1603 / 0.1904 / 0.2141 at
Q² = 2 / 5 / 10 GeV², `moment_A[letter]` 3.016 / 2.933 / 2.938 → 1.508 /
1.466 / 1.469, `moment_B[letter]` 3.114 / 2.919 / 2.822 → 1.557 / 1.459 /
1.411, the dilution-1 rows 9.05 / 8.80 / 8.81 → 4.52 / 4.40 / 4.41 and
9.34 / 8.76 / 8.47 → 4.67 / 4.38 / 4.23. The conclusion stands as decided:
TIGHT = V1 + V2(−); the pre-A6 edge is kept as the recorded
`saturation_tight_v6`; V6's own sub-row fields are unchanged in value; the
TIGHT edge carries no status, so no verdict moved; the sources JSON is not
edited.*

**A7 — the grid row of `t5_delta_moment_bound` gets the harness's declared
quadrature convergence rule as its pass bound: the integrator's error
estimate ≤ `QUAD_REL_TOL` (1e−7) × |M₂|; the toy's 1e−15 absolute stays on
the analytic integrand.** *Why: a row with no bound cannot catch a degraded
interpolant; the relative rule pre-dates D7 as the harness's convergence
criterion for every integrand and is not a number fitted to the observed
error (2–4e−10 relative, 300× inside). The record says the choice was made
after the values were seen, as F7 asked (`../run_2026-10-27/phase_I7_backends.md`
§13 item 2).*

**A8 — a harness whose verdict configuration needs a missing set exits 2
(`validation/benchmarks/README.md` convention 6); the brief's "exit 0" was
wrong.** *Why: exit 0 is for a harness that ran and printed its row; a
harness whose verdict cannot be computed is broken in that environment;
`make_report` exits 0 with the row `blocked (environment)`, which is what
the report's reader sees. Records only.*

**A9 — a missing `lhapdf-config` at configure time is a CMake `FATAL_ERROR`
unless `-DLIPOLGEN_WITH_LHAPDF=OFF` is passed explicitly.** *Why: since D7
every default run needs the tier; a WARNING produces a build whose every
default run is refused at run time, a worse failure than a configure-time
stop that names the cure; the explicit OFF keeps the core-only build (the
regression lens built it after D7: 411 / 411). What changes: `CMakeLists.txt`,
`docs/PACKAGING.md`; CI always has the tier.*

**A10 — HepMC3 output records the run's physics provenance as `GenRunInfo`
attributes — the Δ convention, the rc mode and tail model, the cluster wave,
the structure-function backends and their extrapolation (the keys `meta`
already carries) — and `docs/HEPMC3_CONVENTION.md` documents them; the T6
row re-pins the document's sha256 with a dated comment naming A10.** *Why: a
HepMC3 file is the generator's deliverable; four defaults moved in run
2026-10-27, and files written before and after differ by ×0.43 on A_zz^tag
or ×0.5 on the cos 2φ amplitude with nothing in them to say why; the npz has
the keys, the HepMC3 file does not; the cost is one documented re-pin
(`../run_2026-10-27/phase_I4_jm_convention.md` §8 item 10,
`../run_2026-10-27/phase_I5_polradfull.md` §9 item 2,
`../run_2026-10-27/phase_I6_cluster_wave.md` §9 item 4,
`../run_2026-10-27/phase_I7_backends.md` §13 item 6). What changes: the
writer, the convention document (an appended section), the T6 harness's
pinned sha256 and any clause that reads the changed region, tests.*

**A11 — `cluster_vmc_mc_sigma` is bounded: |σ| ≤ 3, beyond which
`validate()` refuses.** *Why: the knob is "a number of the tables' printed
sigmas", a statistical band; 3σ is the widest reading of such a band; 1000 σ
builds a wave (α–d P_D 0.133) that prices nothing the tables say
(`../run_2026-10-27/phase_I6_cluster_wave.md` §9 item 11); the knob matrix
uses ±1. What changes: one refusal, its test, the USAGE paragraph.*

**A12 — the FSI weight-table ceiling stays at 40 mb; the clipped cell
fraction is printed on the banner of every FSI run beside its `meta` key,
and the ceiling is entered in the registry as an unpriced systematic with
the measured fractions (⁷Li 4.41 % / 9.52 % on the default wave); a band on
the ceiling is filed with the next FSI design.** *Why: the ceiling regulates
the low-momentum NN cross section inside the Glauber weight, so a band on it
would price a regulator, not a physics uncertainty; widening the grid moves
the clipped cells, not the regulator; printing it makes the bias visible on
every run while the design is pending (`../run_2026-10-27/phase_I6_cluster_wave.md`
§9 item 3).*

**A13 — a `li6_alpha_vmc` block joins the polligen port gate: the vector and
tensor dilutions of the VMC α–d channel, dumped from the sibling by the dump
script and compared at a declared relative 1e−8; the C++ and Python gates
read it.** *Why: D6 moved the default to a configuration no external gate
covers; the design measured the agreement at 5.4e−9 / 3.2e−9
(`../run_2026-10-27/phase_I6_cluster_wave.md` §9 item 2), so 1e−8 is a
declared tolerance a factor 2–3 above the measured residual — it catches a
sign or normalisation slip and nothing less; the block covers the α–d motion
only (the sibling keeps the scenario deuteron and refuses `vmc` on ⁷Li), and
says so. The dump is ADDITIVE: the diff of `tagged.json` must show the new
block and nothing else.*

**A14 — I5's two changed readings are recorded as facts of the default tail,
nothing re-tuned: the `rc_qe_tensor_scale` stand-in at scale 1 is 14.5× the
band half-width at x = 0.30 on `polrad-full`, and `qe_kf_gev` is the smaller
of the two quasi-elastic knobs there; the registry row and
`docs/PHYSICS_CHANNELS.md` state them; the quasi-elastic tensor term's own
calculation is the filed item.** *Why: re-deriving a band edge from the
model that now produces the central value would tune a systematic to itself;
the stand-in's size on the default tail is the honest statement of how
unpriced that systematic is (`../run_2026-10-27/phase_I5_polradfull.md` §9
item 5). Records only.*

**A15 — Δ is the physical input and the coefficient is what the sibling
changes: for the same Δ (the bag `C_BAG`, 𝒩 = −0.292) every amplitude of
`PolarizedLithiumSim` and of the letter halves, and so do the significances;
the Δ inputs stay; the change list is `../run_2026-10-27/phase_I4_jm_convention.md`
§7. The sibling is read-only here; the author applies it, and the dump
script's preflight then reads ½.** *Why: D4 settled that the tree's Δ is
Jaffe–Manohar's; the letter cites Jaffe–Manohar for its Δ term and so means
Δ in their units; letting the amplitude stand would redefine the letter's Δ
as 2 × C_BAG, a bag-model value no source carries. Records only.*

**A16 — USAGE §8's D7 throughput is re-measured on this machine while it is
idle, a fresh process each, median of 3, and the raw log is kept in this
run's phase record.** *Why: a number without a log is a claim; the
2026-10-30 values were taken under load (`../run_2026-10-27/phase_V_verification.md`
H6). Phase P6.*

**A17 — one rule: an explicit request for the DEFAULT value of a knob the run
never reads is accepted and labelled `not-read` (D7's behaviour); the D6
refusal of `--cluster-wave vmc` off the tags and the D5 "accepted, recorded
`refused`" of `--rc-tail-model polrad-full` on a tail-less run both become
that label. A non-default value keeps its refusal.** *Why: CONVENTIONS' own
criterion — REFUSE when a value would name a variation of a piece that did
not run — is about variations; an explicit default is not one and prices
nothing, so a refusal costs a scan a cell for no information, and an
`at_default` row that says `refused` misreads; three behaviours for one
situation were a trap for scan scripts (`../run_2026-10-27/phase_V_verification.md`
C4). What changes: `make_config`, the CLI and `validate()` for the wave; the
D5 knob row; the three-ways test becomes the one-way pin; the CONVENTIONS
paragraph; USAGE.*

The smaller points of run 2026-10-27's list (its "Also open, smaller"
paragraph: I4 §8 items 1–4, 6, 9; I5 §9 items 1, 3, 4, 7; I6 §9 items 6, 7,
9, 10; I7 §13 items 4, 5, 8 and §14) are outside this delegation and stand
as recorded.
