# polligen reference tables

Numeric reference tables dumped from the Python event generator `polligen`
(`PolarizedLithiumSim/evgen/polligen`) by `dump_polligen_reference.py`, for a
C++ reimplementation (LiPolGen) to be checked against at `rtol=1e-12`
unless a table's description says otherwise (a few quantities in polligen
are themselves iterative/quadrature results, e.g. `spin.populations_maxent`,
`TaggedModel.norm`, `p2_moment_mixture` -- those are called out below).

Every JSON file has an `"inputs"`-shaped section (grids, axes, kernel /
scenario configuration, constants) and an outputs section with the actual
computed numbers, organized per table below. Floats are written with Python's
`repr()`-precision JSON serialization, i.e. full IEEE-754 double precision
round-trips exactly.

Regenerate with:

```
python3 LiPolGen/validation/dump_polligen_reference.py
```

This script only *imports* `polligen`/`polli_fastsim` from `PolarizedLithiumSim`;
it does not modify anything there.

## spin.json

Spin-density matrices and multipole moments, from `polligen.spin`
(PolarizedLithiumSim/evgen/polligen/spin.py).

- `wigner_d`: `spin.wigner_d(j, beta)` (spin.py:68), package m-ordering
  (+J...-J, see `spin.m_values`, spin.py:40).
- `clebsch_gordan`: `spin.clebsch_gordan(j1,m1,j2,m2,j,m)` (spin.py:46),
  Racah formula, exact for these half-integer/integer spins.
- `spins.j1` / `spins.j3_2`: for J=1 and J=3/2, several population
  vectors p_m (m ordered +J...-J):
    * `pure_m=<m>`      -- one entry 1, rest 0
    * `uniform`         -- all populations equal
    * `maxent_pz=<pz>`  -- `spin.populations_maxent(j, pz)` (spin.py:210),
      spin-temperature p_m ~ exp(beta m); pz=8/13 (J=1) and pz=0.7 (J=3/2)
      are the t=3 rational anchors of evgen/tests/test_bookkeeping.py,
      giving tensor moments 4/13 and 0.4 respectively (bisection tol
      1e-13, so results are good to ~1e-13, not full double precision --
      compare at a looser tolerance than rtol=1e-12 for this table only).
    * `explicit_pz=..._pzz=...` / `..._t=...` -- `spin.spin1_populations`
      (spin.py:180) / `spin.spin32_populations` (spin.py:193), closed-form,
      exact to double precision.
  For each population vector: `moments_own_axis` is
  `spin.moments_along_axis(j, populations)` (spin.py:165) -- (vector,
  tensor_zz[, octupole_z]) along the population's OWN axis. Then for
  each of the 5 `axes` (theta, phi): `rho` is
  `spin.SpinDensity(j, populations, theta, phi).rho_lab()` (spin.py:247),
  a (2J+1)x(2J+1) complex matrix (m-ordering +J...-J), split into
  `re`/`im`; `vector_polarization` is `spin.vector_polarization(rho, j)`
  (spin.py:138); `tensor_polarization` is `spin.tensor_polarization(rho, j)`
  (spin.py:144, J=1: <3Jz^2-2>; J=3/2: <3Jz^2-J(J+1)>/3);
  `octupole_moment` (J=3/2 only) is `spin.octupole_moment(rho, j)`
  (spin.py:156); `lab_moments` is `SpinDensity.lab_moments()`
  (spin.py:251), which must equal (moment_own_axis * n_hat) for the
  vector and (tensor_own_axis * P2(cos theta)) for tensor_zz -- this is
  the analytic identity `test_rotated_moments_analytic` in
  evgen/tests/test_spin.py pins.

## xsec.json

`polligen.xsec.InclusiveKernel` (PolarizedLithiumSim/evgen/polligen/xsec.py)
and the underlying `polli_fastsim.asymmetries` analytic formulas
(PolarizedLithiumSim/fastsim/polli_fastsim/asymmetries.py), evaluated
on TOY structure functions only.

`constants`: `structure.ALPHA_EM` (structure.py:21, =1/137.036),
`structure.GEV2_TO_PB` (structure.py:38), `asymmetries.TENSOR_LL_SIGN`
(=-1.0 since 2026-08-29: the LITERATURE convention, Cosyn et al. Eq. 27
/ HERMES, A_zz = -(2/3) b1/F1), `beams.LI6_CLUSTER_POLARIZATION`
(=0.81123, the whole-nucleus 6Li vector polarization of the cluster
picture, from which beams.LI6's per-nucleon slots are a third each),
`asymmetries.M_NUCLEON` (=0.9383 GeV,
the FREE-nucleon mass used in gamma^2 = 4 M^2 x^2/Q^2), and per-ion
`beams.Ion` fields (A, Z, N=A-Z, spin, eff_pol_p, eff_pol_n,
mass_per_nucleon) for `beams.LI6`/`beams.LI7` (beams.py:197,199).

`toy_formulas`: the exact closed forms used everywhere below --
`ToyF2` (structure.py:41-59), `r_sigma_lt` (structure.py:62),
`NuclearF2.f2a`/`f1a` (structure.py:303-338), `ToyG1`
(polarized.py:26-82), `toy_b1(..., mode='toy')` (polarized.py:444,451 --
NOTE: the default mode='digitized' reads a CSV table and is NOT used
here; every b1 in this file is the pure formula), `toy_delta_gluon`
(polarized.py:536), the b2 default `2*x*b1` (xsec.py `tables()`,
xsec.py:234/238), and `g2_ww` (xsec.py:85, Wandzura-Wilczek by
quadrature -- reproduce it with the SAME npts=96 trapezoid rule to
match at 1e-12; it is not a closed form).

`results.<6Li|7Li>`: `beam_config` is `beams.default_configs(name)[1]`
(beams.py:224, the MID energy point) with its `sqrt_s_per_nucleon`
(beams.py:213) and `s = sqrt_s_per_nucleon^2`.

`kernel_default` (6Li: J=1, b1_func/delta_func = the toy formulas
above, target_mass=False; 7Li: J=3/2, all rank-2 kwargs None so
b1/b2/delta are identically zero per xsec.py's documented default) and
`kernel_default_target_mass` (same kernel with target_mass=True, the
exact finite-gamma E143 vector sector) each contain:
  * `tables`: `InclusiveKernel.tables(x, q2, with_g2=True)` (xsec.py:217)
    -- f1, f2, g1, g2, b1, b2, delta on `grid_points`.
  * `asymmetries`: `y = q2/(s x)`, `R = r_sigma_lt(x,q2)`,
    `asymmetries.depolarization_d` (asymmetries.py:44),
    `asymmetries.a_parallel` (asymmetries.py:59, the MASSLESS g1/F1
    form -- this is what `kernel_default`'s vector sector must equal
    exactly, per `test_vector_sector_matches_a_parallel`),
    `asymmetries.azz` (asymmetries.py:101, only when b1 is nonzero)
    and `asymmetries.a_cos2phi` (asymmetries.py:111, only when delta
    is nonzero).
  * `amplitudes`: `InclusiveKernel.amplitudes(tables, x, q2, s, state,
    with_perp=True)` (xsec.py:345) -- (w_avg, a1, a2) of
    W = 1 + w_avg + a1 cos(phi') + a2 cos(2 phi') -- for every
    (lam_e in {+1,-1}) x (5 `axes`) x (each spin projection m of the
    ion), at pe=0.7 (`PE` in this script).
  * `dsigma`: `InclusiveKernel.dsigma(x, q2, phi, s, state,
    with_perp=True)` (xsec.py:388) on `phi_grid` (8 uniform points in
    [0, 2pi)) for 3 representative states, at each of the 6
    `grid_points` (`dsigma_per_point[i]` is the phi-array at
    grid_points[i]).  `dsigma_unpol` is `InclusiveKernel.dsigma_unpol`
    (xsec.py:383), i.e. `structure.dsigma_dx_dq2` (structure.py:341) on
    the per-nucleon F2.
  * (tensor_gamma variant only) `tensor_gamma`: the pieces of the EXACT
    finite-gamma tensor sector (Cosyn et al. Eqs. 9/10/14/16/17/24,
    plans/08 D2) -- `xsec.theta_q_cos_sin` (Eq. 24),
    `xsec.cosyn_tensor_sfs` (Eqs. 17a-17e, on the block's own b1..b4),
    `xsec.cosyn_unpolarized_sfs` (Eq. 16), and the (h0, h1, h2)
    harmonics `InclusiveKernel._tensor_harmonics_gamma` returns for
    every (axis, m).  Those blocks are dumped from a kernel built with
    `tensor_gamma=True` and both higher-twist slots filled
    (b3_func=0.05*f1, b4_func=-0.02*f1 -- SCENARIO shapes, b3 and b4
    are unmeasured), so their `amplitudes` carry the exact b-sector
    while every other block carries the massless one.
  * (target_mass variant only) `target_mass`: `xsec.gamma_squared`
    (xsec.py:118), `xsec.epsilon_gamma` (xsec.py:129),
    `xsec.depolarization_gamma` (xsec.py:136), `xsec.eta_gamma`
    (xsec.py:142), `InclusiveKernel.a_parallel` with target_mass=True
    (xsec.py:260, the exact E143 D_gamma(A1+eta*A2) form) alongside the
    massless `asymmetries.a_parallel` for comparison.

`results.7Li.kernel_scenario_rank2`: the SAME structure with an
explicit rank-2 scenario (b1_32_func=0.05*f1, delta_32_func=-1e-2*f1,
matching `evgen/tests/test_tensor_convention.py`'s
`test_spin32_rate_and_cos2phi_channels_are_now_consistent`), so the
J=3/2 tensor geometry Q_NN=(3m^2-J(J+1))/3, t_geo=Q_NN*P2(cos theta_S),
c_eff=3*Q_NN (xsec.py:311-341) is exercised numerically too.

## beams.json

`polli_fastsim.beams` (PolarizedLithiumSim/fastsim/polli_fastsim/beams.py).

`constants`: `PROTON_TOP_MOMENTUM` (beams.py:35, =275.0 GeV),
`PROTON_MASS` (beams.py:36), `PROTON_CONFIG_ENERGIES` (beams.py:52,
the 41/100/275 GeV proton energies), `ELECTRON_ENERGIES` (beams.py:203,
5/10/18 GeV), `NUCLEUS_MASS` (beams.py:43, keyed `"<name>_<A>_<Z>"`).

`ions`: every `beams.IONS` entry (beams.py:194-201) -- A, Z, N=A-Z,
spin, eff_pol_p, eff_pol_n, `Ion.mass_per_nucleon` (beams.py:169, the
PHYSICAL nuclear mass / A, not amu), `Ion.momentum_per_nucleon_max`
(beams.py:176, = 275*Z/A GeV).

`configs.<6Li|7Li|d>`: `beams.default_configs(name)` (beams.py:224),
the low/mid/top (5x41-like / 10x100-like / 18x275-like) reference scan,
with `BeamConfig.sqrt_s_per_nucleon` (beams.py:213,
= sqrt(4*E_e*p_ion_per_nucleon)) and `gamma_ion = sqrt(p^2+m^2)/m` computed here from `ion_momentum_per_nucleon` and
`Ion.mass_per_nucleon` (not stored directly by beams.py, but a plain
function of the two fields above).

`gamma_of_proton_energy`: `beams.gamma_of(proton_energy)` (beams.py:86)
at the three reference proton energies -- the RING Lorentz factor that
every ion at that configuration is gamma-matched to.

## tagged.json

`polligen.tagged` (PolarizedLithiumSim/evgen/polligen/tagged.py) using
the DEFAULT channel constructors (all default beta=0.30, and for 6Li
p_d=tagged.P_D_LI6, for the deuteron control p_d=tagged.P_D_DEUTERON).

**THIS FILE NO LONGER TRACKS polligen FOR THE SPIN-1 MODEL BLOCKS (2026-09-06).**
`polligen.tagged._amp2_table` (tagged.py:243-248) sums the partial waves with
NO `i^L`: it feeds `psi_L` into an amplitude that needs `phi_L = i^L psi_L`.
For an S+D channel that is not a global phase -- it is +1 on L=0 and -1 on
L=2 -- so polligen's tagged sector carries the S-D interference sign
INVERTED, against Cosyn-Weiss II Eq. (6.12) (by up to 2.74 in an asymmetry
whose whole range is [-2, 1]), against LiPolGen's own deuteron quadrupole
sign gate, and against LiPolGen's own b1 sector, which applies the phase
explicitly.  LiPolGen fixed it on 2026-09-06 (`src/core/tagged.cpp`
`build_amp2`, one `(-1)^floor(L/2)`); polligen was NOT touched and still
carries the bug.  Regenerating `channels.li6_alpha.model` or
`channels.deuteron.model` from polligen would therefore re-bake the refuted
sign and the rtol-1e-12 gate would go on certifying it.

Those two blocks are instead RE-PINNED from the fixed C++ library by
`validation/repin_tagged_from_lipolgen.py` (provenance `"LiPolGen post-fix,
formerly polligen"`, recorded in the file's own `provenance` /
`provenance_note` keys), and `dump_polligen_reference.py` carries them
through unchanged rather than overwriting them.  What moved in the re-pin,
MEASURED 2026-09-06 as the re-pinned file against the pre-fix dump
(`git show HEAD:validation/reference/tagged.json`), not copied from the
investigation's own table: `n_of_kc` (up to +725% on 6Li, +19215% on the
AV18 deuteron control), `struck_populations` (up to 0.86 absolute),
`p2_moment` (sign flip, x1.28 to x2.03), `p2_moment_mixture_uniform`
(-5.3160743e-05 -> -5.2153460e-05 on 6Li, 1.9e-2 rel; -5.3694953e-05 ->
-5.3337762e-05 on the deuteron, 6.7e-3 rel -- it is gated at 1e-9, so it had
to be re-pinned too), `norm` (up to 5.8e-5 rel: 6Li M=0 1.000026689626 ->
0.999968571520; 4.0e-5 on the deuteron), `population_integrated` (up to
3.0e-6 abs: 6Li M=0 0.947966373524 -> 0.947963349333) and the dilutions
(<=4.3e-6 rel).  Until 2026-09-06 this list read `norm` "+2e-5" and
`population_integrated` "<=5e-7 abs" and did not mention
`p2_moment_mixture_uniform` at all -- three figures taken from
`07_cw_sign_investigation.md` section 6.1 rather than measured here, where
that section records the last one as not moving.  Reason, derivations and
the full before/after tables: `docs/benchmarking/07_cw_sign_investigation.md`
and `docs/open_items/run_2026-09-06/phase_CW_numbers.md`.

EVERYTHING ELSE IN THE FILE IS STILL polligen's, untouched: `waves`,
`base`, `beam_configs`, `boost_spectator`, `P_D_LI6`, `P_D_DEUTERON`, the
channel scalars, and the WHOLE of `channels.li7_alpha.model` -- 7Li alpha-tag
is a single L=1 wave, so the common `i` is a global phase, and the fix
moves its `n_of_kc` and `p2_moment` by exactly zero (measured, bit for bit).
The re-pin script re-checks that block against the live C++ at rtol 1e-12
instead of overwriting it.  No other reference JSON moved.

`channels.<li6_alpha|li7_alpha|deuteron>`: built by
`tagged.li6_alpha_channel()` / `li7_alpha_channel()` / `deuteron_channel()` (tagged.py:181,188,195). `base` is the underlying
`polli_fastsim.spectator.ClusterChannel` (spectator.py:111): `m_spec`
(spectator.py:124), `m_beam` (spectator.py:128, physical nuclear mass),
`m_partner` (spectator.py:143), `kappa` (spectator.py:168, sqrt(2*mu*S)).

`waves[i]`: one `tagged.Wave` (tagged.py:126) of `channel.waves` --
`l`, `prob`, `beta`, and `radial` = `Wave.radial(k_grid, kappa)`
(tagged.py:140; L=0 Hulthen 1/(k^2+kappa^2)-1/(k^2+beta^2), L=1
k/((k^2+kappa^2)(k^2+beta^2)), L=2 k^2/((k^2+kappa^2)(k^2+beta^2)^2)) on
`k_grid` = 40 points linspace(1e-4, 1.2, 40) GeV -- UNNORMALIZED (the
TaggedModel normalizes internally, see below).

`model`: `tagged.TaggedModel(channel)` (tagged.py:203) at its DEFAULT
grid (k_max=1.2, nk=280, nc=96 -- `grid` records the exact edges/sizes,
which the C++ port must reproduce bit-for-bit since `n_of_kc` looks up
the NEAREST cell at or below the query point, `np.clip(np.searchsorted(model.k, k)-1, 0, nk-2)`, tagged.py:256-258 -- this script replicates
that exact lookup in Python to pick `k_pts`/`c_pts` cell values, so the
reference values are themselves grid CELL values, not interpolated).
  * `n_of_kc[i].n_at_k_c`: `TaggedModel.n_of_kc(M)` (tagged.py:248) at
    `k_pts` x `c_pts` (6x5), indexed [k][c].
  * `struck_populations[i].p_at_k_c`: `TaggedModel.struck_populations(M)`
    (tagged.py:260), shape (n_mS, nk, nc), sliced at the same
    (k_pts, c_pts) cells; `m_s_order` gives the m_S ordering (+S_c...-S_c).
  * `population_integrated[i].p_m_s`: `TaggedModel.population_integrated(M)`
    (tagged.py:266), the k/khat-integrated channel-spin populations.
  * `norm[i].value`: `TaggedModel.norm(M)` (tagged.py:272), grid
    quadrature of the normalization integral (should be ~1, NOT exact
    -- it is itself a finite-grid quadrature, so match this one at a
    looser tolerance, e.g. 1e-6, not rtol=1e-12).
  * `vector_dilution`: `TaggedModel.vector_dilution()` (tagged.py:292).
  * `tensor_dilution` (only present when s_channel=1, i.e. li6_alpha and
    deuteron -- s_channel=0.5 for li7_alpha, where `tensor_dilution`
    raises `ValueError` by design, tagged.py:299-302; see the JSON's
    top-level `_could_not_call` list): `TaggedModel.tensor_dilution()`
    (tagged.py:299).
  * `p2_moment[i]`: `TaggedModel.p2_moment(M)` (tagged.py:307), per M.
  * `p2_moment_mixture_uniform`: `TaggedModel.p2_moment_mixture(pops)`
    (tagged.py:314) at a uniform fill (also a grid quadrature; same
    looser-tolerance caveat as `norm`).

`boost_spectator`: `tagged.boost_spectator(channel, k, c, phi_k,
p_per_nucleon, theta_s, phi_s)` (tagged.py:335), an EXACT closed-form
boost (no grid quadrature -- rtol=1e-12 applies), at a fixed set of
(k, c, phi_k) x (theta_s, phi_s) x each of the channel's 3
`beam_configs` (`beams.default_configs` for the channel's own beam
species: 6Li, 7Li, or d). Fields: pT, theta, p_lab, R (rigidity ratio),
xL, kx, ky, kz (spin-frame components after rotation), phi_spec (lab
azimuth).

## spectator.json

`polli_fastsim.spectator` (PolarizedLithiumSim/fastsim/polli_fastsim/spectator.py).

`constants`: `M_U` (spectator.py:43), `MASSES` (spectator.py:44, keyed
by particle name), `NUCLEUS_MASS` (spectator.py:75, keyed `"<Z>_<A>"`,
the AME2020-derived physical nuclear masses).

`channels.<name>`: one of `spectator.{DEUTERON_P_TAG, DEUTERON_N_TAG,
HE3_P_TAG, LI6_ALPHA_TAG, LI6_D_TAG, LI7_ALPHA_TAG, LI7_T_TAG}`
(spectator.py:187-207), each a `ClusterChannel` (spectator.py:111):
m_spec (spectator.py:124), m_beam (spectator.py:128), m_partner
(spectator.py:143), kappa (spectator.py:168, sqrt(2*mu*S) with mu the
reduced mass of the FREE spectator+partner). `R_at_k_zero` is the
closed form `(m_spec/spectator_Z)/(m_beam/beam_Z)` (spectator.py
module docstring, verified against `_boost_fragment` at k=0 --
independent of p_per_nucleon, see `boost_fragment` rows with
kx=ky=kz=0).

`momentum_density.n_of_k`: `spectator.momentum_density(k, kappa,
beta=0.30, l_wave)` (spectator.py:212, the UNNORMALIZED |psi(k)|^2 --
Hulthen^2 for l_wave=0, [k/((k^2+kappa^2)(k^2+beta^2))]^2 for l_wave=1)
on `k_grid` (40 points, linspace(1e-4, 1.2, 40) GeV).

`boost_fragment`: `spectator._boost_fragment(channel, p_per_nucleon,
kx, ky, kz, m=channel.m_spec, frag_Z=channel.spectator_Z,
frag_A=channel.spectator_A)` (spectator.py:238) -- i.e. exactly what
`spectator_lab_kinematics` (spectator.py:280) calls internally, but at
FIXED (kx, ky, kz) rather than a random `sample_k` draw: reproducing
`spectator_lab_kinematics`'s own random sampling bit-for-bit would
require replicating numpy's PCG64 bit generator in C++, which is out
of scope for a physics cross-check, so this table isolates the
boost/kinematics formula itself (the only thing `sample_k` feeds it).
Evaluated at p_per_nucleon in {40.8, 99.5, 137.5} GeV (representative
6Li/7Li values) x 4 fixed (kx, ky, kz) [GeV] triples including (0,0,0).
Fields: pT, theta, phi (lab azimuth), p_lab, R (rigidity ratio), xL, k.

## coherent.json

`polligen.coherent` (PolarizedLithiumSim/evgen/polligen/coherent.py).

`constants`: `M_LI6` (coherent.py:47), `GEV_PER_FM_INV` (coherent.py:48,
hbar*c), `RATE_WEIGHT_SYST` (coherent.py:213), `MANTYSAARI_A2_DEUTERON`
(coherent.py:221, keyed by |t| [GeV^2] as a string, value = (a2 m=0,
a2 m=+-1)).

`scenario_defaults`: the field defaults of `coherent.CoherentScenario`
(coherent.py:57): f0=0.04, x_coh=0.01, slope_b=50.0, amp=0.01,
eps_b0=-0.08.

`gaussian_slope`: `coherent.gaussian_slope(r_rms_fm)` (coherent.py:51,
= (r_rms_fm/GEV_PER_FM_INV)^2 / 3) at the matter/charge radii quoted in
the docstring.

`coherent_fraction`: `CoherentScenario.coherent_fraction(x)`
(coherent.py:111, = f0/(1+(x/x_coh)^2)) at the default scenario.

`tag_acceptance`: `CoherentScenario.tag_acceptance(pt_cut)`
(coherent.py:123, = exp(-B*pt_cut^2)) and `mean_t_tagged(pt_cut)`
(coherent.py:128, = pt_cut^2 + 1/B).

`tag_acceptance_angular`: `CoherentScenario.tag_acceptance_angular(
sigma_theta, p_per_nucleon, a_beam=6, n_sigma=10.0)` (coherent.py:132,
= exp(-B*(n_sigma*sigma_theta*a_beam*p_per_nucleon)^2)).

`a2_deformation`: `CoherentScenario.a2_deformation(t_abs, pzz)`
(coherent.py:155, = -(pzz/4)*eps_b0*B*|t|) and
`cos2phi_coefficient_deformation(t_abs, pzz)` (coherent.py:178,
= 2*a2_deformation).

`a2_tagged`: `CoherentScenario.a2_tagged(pt_cut, pzz)` (coherent.py:196,
= a2_deformation(mean_t_tagged(pt_cut), pzz)).

`recoil_lab`: `coherent.recoil_lab(t_abs, phi_t, p_per_nucleon,
x_pom=0.0)` (coherent.py:229, exact closed-form kinematics of the
intact-6Li recoil): pT, theta, R, xL.

`fragment_rigidity`: `coherent.fragment_rigidity(a, z, beam_a, beam_z)`
(coherent.py:352, = (nucleus_mass(z,a)/z) / (nucleus_mass(beam_z,
beam_a)/beam_z)) for every fragment of `coherent.LI6_BREAKUP` /
`LI7_BREAKUP` (coherent.py:306,331), whose raw contents are also given
verbatim in `breakup_tables`.

## bookkeeping.json

`polligen.bookkeeping` (PolarizedLithiumSim/evgen/polligen/bookkeeping.py).

`plans.<name>`: the four standard run-plan constructors, each a list of
`SpinCategory` (bookkeeping.py:36: name, j, populations, lam_e, pe,
theta_s, phi_s, lumi_fraction, and `moments()` = `spin.moments_along_axis(j, populations)`) plus the `RunPlan`
(bookkeeping.py:52) true polarizations and `lumi_shares(1e6)`
(bookkeeping.py:74, absolute lumi at an arbitrary reference total of
1e6 pb^-1 = 1 fb^-1 -- purely `lumi_fraction * total`, exact rational
arithmetic).
  * `helicity_flip_j1_maxent_anchor` /
    `helicity_flip_j32_maxent_anchor`: `bookkeeping.helicity_flip_plan(
    j, pz, pe)` (bookkeeping.py:80) with pzz=None (spin-temperature
    fill via `spin.populations_maxent`) at the t=3 rational anchors
    pz=8/13 (J=1, expect pzz_true=4/13) and pz=0.7 (J=3/2, expect
    pzz_true=0.4) -- see `spin_temperature_ladder_t3` below for the
    exact closed-form populations (bisection-derived in the module, so
    match this ladder at a looser tolerance than 1e-12, e.g. 1e-10).
  * `helicity_flip_j12`: J=1/2 branch (no rank-2 moment, populations
    (1+pz)/2, (1-pz)/2 exactly).
  * `helicity_flip_j1_explicit_pzz`: pzz=0.35 explicit (closed-form
    `spin.spin1_populations`, exact to double precision).
  * `helicity_flip_j1_tilted_offset`: theta_s=0.9, phi_s=0.4,
    rel_lumi_offset=1e-4 (boosts the lam_e=+1 share only, see
    bookkeeping.py:102).
  * `tensor_thirds_j1`: `bookkeeping.tensor_thirds_plan(pz=0.7,
    pzz=0.6, rel_lumi_offset=1e-4)` (bookkeeping.py:114): 3 categories
    (+pzz, -pzz, m0-enriched -2pzz).
  * `transverse_tensor_j1`: `bookkeeping.transverse_tensor_plan(
    pzz=0.6)` (bookkeeping.py:140): 1 category at theta_S=pi/2.
  * `tensor_flip_j1`: `bookkeeping.tensor_flip_plan(pzz=0.6,
    share_plus=0.5, rel_lumi_offset=1e-4)` (bookkeeping.py:150).
  * `tensor_thirds_j1_with_offset_azz0`: `bookkeeping.with_offset(
    tensor_thirds_plan(0.7, 0.6), "azz0", 1e-4)` (bookkeeping.py:178) --
    only the "azz0" category's lumi_fraction is scaled by (1+1e-4).

`spin_temperature_ladder_t3`: the closed-form geometric-ladder
populations p_m ~ t^m at t=3 for J=1 and J=3/2 (exact rationals, from
the docstring/tests of bookkeeping.py and test_bookkeeping.py), the
independent construction the maxent anchors above should match.

`rel_lumi_bias`: `bookkeeping.azz_rel_lumi_bias(offset, pzz)`
(bookkeeping.py:193, = -(2/3)*offset/pzz) and
`bookkeeping.apar_rel_lumi_bias(offset, pe, pz)` (bookkeeping.py:199,
= offset/(2*pe*pz)) at offset=1e-4 (the plans/05 reference delta).

## Functions this script could not call

- `TaggedModel.tensor_dilution channel=7Li alpha-tag (quasi-free t)`: s_channel=0.5 != 1: tensor_dilution raises ValueError by design (tagged.py:299-302)
- `SpinCategory.moments() name=apar+ (j=0.5)`: spin.moments_along_axis raises ValueError for j not in {1, 1.5} (spin.py:177); vector moment computed directly instead, tensor/octupole n/a
- `SpinCategory.moments() name=apar- (j=0.5)`: spin.moments_along_axis raises ValueError for j not in {1, 1.5} (spin.py:177); vector moment computed directly instead, tensor/octupole n/a

