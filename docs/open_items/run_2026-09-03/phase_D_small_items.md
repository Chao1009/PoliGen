# Phase D small items — D3 (the FSI "variants"), D4 (`PDF:PomSet`), D5 (the |t| ceiling)

Measured 2026-09-04 at commit `903fcc9`, tree clean, **read-only**: this file is
the only thing the task added, no source, doc, test or reference file was
touched. Baseline re-run on this machine before reporting and unchanged —
**390 doctest cases / 17 217 416 assertions / 1 skipped / 0 failed**,
**345 pytest passed** (50.3 s), gate `1094 references checked (strict), 97
ranges, 6 external, 0 broken, 6 allow-listed`.

Every number below was produced through `PYTHONPATH=build/python` at that
commit; §D3.6 lists the commands. Where a claim is *not* measured it says so.

Two of the three items turn out to rest on a premise that does not survive
measurement (D3 outright, D4 partly), and all three ran into **the same
third defect**: the knob is bound and reachable but is **not written to the
npz `meta`**, so two runs that differ in it are indistinguishable in the
record. §D6 collects that.

---

## D3 — the FSI variants

### D3.1 The premise is false as stated

`PLAN.md:332-333` says *"The per-nucleon Glauber FSI variant is algebraically
identical to the cluster one — the two 'variants' are one."* As implemented
they are not, and the header already says so correctly and in more detail than
the inventory line preserved.

`src/core/fsi.cpp:191-198` (anchor read 2026-09-15) — the whole of the difference:

```cpp
  if (opt_.variant == FsiVariant::GlauberNucleon) {
    return static_cast<double>(spec_a_) * g1;          // (b)  A * (Gamma_N conv T_a)
  }
  std::complex<double> keep(1.0, 0.0);
  const std::complex<double> one_minus = 1.0 - g1;
  for (int i = 0; i < spec_a_; ++i) keep *= one_minus;
  return 1.0 - keep;                                    // (a)  1 - (1 - Gamma_N conv T_a)^A
```

* **`GlauberCluster` (a)** — Γ_a(b) = 1 − [1 − (Γ_N ⊛ T_a)(b)]^A. Glauber
  shadowing over the cluster's own Gaussian point-nucleon profile.
* **`GlauberNucleon` (b)** — Γ_a(b) = A·(Γ_N ⊛ T_a)(b). The **single-scattering
  (optical, unshadowed) limit**, i.e. the first term of the same binomial, so
  σ_Xa = A·σ_XN exactly.

The *true* algebraic-identity statement, which `fsi.hpp:107-129` makes and
which the inventory line compressed into something else, is narrower: for an
**uncorrelated** cluster density the Ciofi degli Atti–Kaptari per-nucleon
product ⟨Π_i[1 − Γ_N(b − s_i)]⟩ factorises into Π_i⟨1 − Γ_N(b − s_i)⟩ =
[1 − (Γ_N ⊛ T_a)]^A — which is form (a). So *a literal per-nucleon code path
would reproduce the cluster variant*, and that is why the shipped (b) is
deliberately something else. The header is right; the inventory line is not.

**Which of the three the task offers:** the two shipped variants are
**genuinely different by construction** — not identical, and not
"identical only for the shipped parameters". The identity that *does* hold (CK
product ≡ cluster form on an uncorrelated density) is parameter-independent:
it is an algebraic factorisation, not a coincidence of σ_XN = 40 mb.

### D3.2 Measured — the two profiles

Both built on the production S+D ⁶Li α tag (`li6_alpha_channel()`), σ_XN = 40 mb,
ε = −0.5, B_XN = 6 GeV⁻², all other options default:

| quantity | `GlauberCluster` | `GlauberNucleon` |
|---|---|---|
| σ_tot(X–α) [mb] | 131.045 | 160.006 |
| σ_el(X–α) [mb] | 35.224 | 68.187 |
| σ_el/σ_tot | 0.269 | 0.426 |
| effective slope B_α [GeV⁻²] | 27.189 | 23.978 |
| `survival()` | 0.520239 | 0.582899 |

(σ_tot = 131.045 and σ_el = 35.224 reproduce `fsi.hpp:88-100` exactly; the
measured α–p ratio quoted there is 0.258.)

### D3.3 Measured — the two weights, on the same events

Direct evaluation on a 630-point (|k|, cos θ_k) grid, |k| ∈ [0.02, 0.6] GeV × 30,
cos θ_k ∈ [−1, 1] × 21, at W = 5 GeV, Q² = 10 GeV², x = 0.1:

* max |w_nucleon − w_cluster| = **6.001**
* w_nucleon/w_cluster ∈ **[0.817, 23.68]**, max |ratio − 1| = 22.68
* `np.allclose(rtol=1e-14)` → **False**

Pointwise rows (θ_k = 90°):

| \|k\| [GeV] | cluster | nucleon | ratio |
|---|---|---|---|
| 0.05 | 0.67707 | 0.61125 | 0.903 |
| 0.10 | 0.39135 | 0.32094 | 0.820 |
| 0.20 | 0.27747 | 0.92728 | 3.342 |
| 0.30 | 1.13927 | 4.07433 | 3.576 |
| 0.50 | 0.63510 | 6.54603 | 10.307 |

(The 0.10 GeV row reproduces the pair `docs/PHYSICS_CHANNELS.md:343` records
for the **production S+D tag**, 0.391 vs 0.321; `tests/test_fsi.cpp:258-281`
pins the *pure-S prototype* channel's 0.381 / 0.309 instead, which is why the
two pairs differ.)

**Through the pipeline, same events three times.** `--channel tagged-6Li-alpha
--events 20000 --seed 1234 --fsi {off,glauber-cluster,glauber-nucleon}`, plan
`tensor-thirds` at P_z = 0.7 / P_zz = 0.6 / P_e = 0.7, at the default
σ_XN = 40 mb (one stream at one end of the mandatory 20–40 mb band). The
kinematic columns `k`, `cos_theta_k`, `x`, `q2` are **bit-identical** across the
three files, so this is one event stream reweighted three ways:

| | `off` | `glauber-cluster` | `glauber-nucleon` |
|---|---|---|---|
| Σ weight (20 000 events) | 20000.0 | 10419.1 | 11632.1 |
| integrated survival Σw/Σw_off | 1 | **0.520954** | **0.581605** |
| per-event factor: mean | 1 | 0.520954 | 0.581605 |
| per-event factor: min … max | 1 | 0.0208 … 1.4043 | 0.2139 … 6.9523 |

Per-event ratio w_nucleon/w_cluster, percentiles [1, 5, 25, 50, 75, 95, 99] =
**0.821, 0.827, 0.853, 0.897, 0.934, 3.337, 5.935**; maximum **68.52**.
**99.50 %** of the 20 000 events differ by more than 1 %.

**Verdict: the two variants are different on essentially every event, and by
construction, not by parameter choice.** D3 should be closed as *premise
false*, and `PLAN.md:332-333` retracted where it stands (it is the only place
in the tree that carries the claim — `docs/PHYSICS_CHANNELS.md:343` and
`docs/CONVENTIONS.md:454` (anchor read 2026-09-15 "40); variant") already state it correctly).

### D3.4 What a genuinely NON-identical per-nucleon variant would have to include

`fsi.hpp:126-129` already names both missing pieces as an open TODO, and they
are the right two:

1. **the centre-of-mass constraint** Σ_i **s**_i = 0 — the one-body density is
   a *lab* density and the A nucleon positions in the product are not
   independent;
2. **short-range NN correlations** in the cluster density — i.e. the two-body
   density ρ₂(**s**_i, **s**_j), which is what makes ⟨Π_i(·)⟩ ≠ Π_i⟨·⟩.

**A correction to the task's phrasing, stated plainly.** The task asks for
"Glauber over the individual nucleons' positions … **with the nucleon–nucleon
σ** and the cluster density". σ_NN is **not** the input to this model. The
rescattering here is the DIS debris X off a **spectator nucleon**: the
amplitude is X–N, σ_XN (`GlauberFsiOptions::sigma_xn_mb`, band 20–40 mb,
`fsi.hpp:131-139`), and a σ_NN would enter only if one were computing the
spectator cluster's *internal* absorption, which is not what the weight is.
Substituting σ_NN would change the physics, not the variant.

**Do the in-tree densities allow it? No for (2); (1) is not a variant at all.**

* `data/vmc/density/{he4,li6}.density` are **one-body point-proton** densities
  (ANL VMC AV18+UX; `data/vmc/density/README.md`). A one-body density is
  *exactly* the T_a the cluster form already convolves. Feeding it in instead of
  the Gaussian would give a non-Gaussian T_a — a better input to the **same**
  variant, not a different variant.
* The **two-body** density needed for (2) is **not in the tree**, and the ANL
  density page the README fetches from publishes one-body densities only. So
  the SRC correction cannot be built from committed data; it needs an input the
  repository does not have.
* (1) is a change of T_a alone, and `cluster_point_a2_fm2` (`fsi.hpp:177`) is
  its single code home. Gartenhaus–Schwartz on a Gaussian gives
  a²_int = a²(1 − 1/A) = 0.7001 × 0.75 = **0.5250 fm²** for the α. Measured
  effect, recomputed analytically against the same formulas (and reproducing
  the shipped 131.04 / 35.22 / 0.269 exactly as a control):

  | a² [fm²] | σ_tot [mb] | σ_el [mb] | σ_el/σ_tot |
  |---|---|---|---|
  | 0.7001 (shipped, charge radius with the proton folded out) | 131.04 | 35.22 | **0.269** |
  | 0.5250 (c.m. removed) | 124.91 | 37.48 | **0.300** |

  σ_tot falls 4.7 %, and σ_el/σ_tot moves **away** from the measured α–p 0.258
  by about as much as it started away from it. Worth knowing before anyone
  implements the c.m. correction as an improvement: on the one number the model
  is checked against, it makes the agreement worse.

**Conclusion.** The item as filed is not open, and the thing that *is* open
(the `fsi.hpp:126-129` TODO) cannot be closed with in-tree data. Recommend:
retract the inventory claim, keep the TODO, and record in the header that the
two-body density is the missing input and where it would come from.

### D3.5 A real defect found on the way — FSI is absent from the npz `meta`

The three runs of §D3.3 produce **identical `meta` dicts**. Verified key by key:
`meta` for `--fsi off`, `--fsi glauber-cluster` and `--fsi glauber-nucleon`
compare **equal**, while their event weights differ by up to a factor 68.5 and
their total rate by 48 % / 42 %.

`python/bindings.cpp:440-556` writes 23 base keys, the four `b1_*` keys
**unconditionally**, and the whole `rc_*` block **conditionally**. There is no
`fsi` key of any kind. The comment that justifies the `b1_*` block —
`python/bindings.cpp:462-466` (as of adec442), *"without these three keys three otherwise
identical npz files are indistinguishable — which is exactly the 'never quote a
single row' rule failing silently"* — describes the FSI situation exactly, and
was never applied to it. `fsi_sigma_mb` (the documented 20–40 mb band, the one
knob the header says must never be quoted as a single row) is likewise absent.

The setting *is* printed at the run surface, `python/lipolgen/cli.py:1436-1443` (anchor read 2026-09-15),
but that is stdout, not the record.

Against the run's own ground rule — *exposed in bindings, reachable from
`cli.py`, recorded in the npz meta, gated by a pytest* — the first, second and
fourth hold (`python/tests/test_module.py:456,468`); **the third does not**.

Suggested minimal fix, following the **`rc` precedent** rather than the `b1`
one so that today's default (`--fsi off`) stays byte-identical: emit the block
only when `p.fsi_weight()` is non-null, with `fsi` (the
`pipeline_fsi_name`), `fsi_sigma_mb`, `fsi_survival`, `fsi_sigma_cluster_mb`,
`fsi_clipped_grid_fraction`, and the four ramp anchors when
`formation_ramp` is on. Plus one pytest asserting the three-way distinguishability.

### D3.6 Reproduction

```bash
source env.sh
for v in off glauber-cluster glauber-nucleon; do
  python -m lipolgen.cli --channel tagged-6Li-alpha --events 20000 --seed 1234 \
                         --fsi $v --npz /tmp/fsi_$v.npz --quiet
done
# then: load the three, check the kinematic columns are equal and diff `weight`
# and `meta`.  Profiles: lg._lipolgen.GlauberFsiWeight(li6_alpha_channel(), opt)
# with opt.variant = FsiVariant.{GlauberCluster,GlauberNucleon}.
```

---

## D4 — `PDF:PomSet`

### D4.1 What PYTHIA 8.317 ships

Deps live at `/home/cpeng/Projects/polli/deps/install/share/Pythia8` (**not**
`LiPolGen/deps` — the task's path does not exist). Version confirmed
`PYTHIA_VERSION 8.317`, `xmldoc/Version.xml:9` `versionNumber 8.317`,
`versionDate 20260120`.

`xmldoc/PDFSelection.xml:482-563` defines `PDF:PomSet`, **default 6**, with
integer options 1–15 plus the three external forms
(`LHAPDF5:`/`LHAPDF6:`/`LHAGrid1:`). Grids on disk in `pdfdata/`:
`pomH1FitA.data`, `pomH1FitB.data`, `pomH1FitBlo.data`, `pomH1Jets.data`,
`pomactw{b14,d14,d19,sg14}.pds`, `GKG18_DPDF_Fit{A,B}_{LO,NLO}_0000.dat`.

Flavour structure read out of the PYTHIA source
(`deps/src/pythia8317/src/PartonDistributions.cc`) and the dispatch out of
`BeamSetup.cc:1343-1392`:

| set | fit | PYTHIA class | order | light flavours | charm |
|---|---|---|---|---|---|
| 1 | `N x^a (1−x)^b` toy | `PomFix` | Q²-indep. | tunable mix | none |
| 2 | π⁰ distributions | (idIn → 111) | — | pion | pion |
| 3 | H1 2006 Fit A | `PomH1FitAB(iFit=1)` | NLO | **singlet, all six equal** | **≡ 0** |
| 4 | H1 2006 Fit B | `PomH1FitAB(iFit=2)` | NLO | **singlet** | **≡ 0** |
| 5 | H1 2007 Jets | `PomH1Jets` | NLO | singlet (`sn/6`) | **≠ 0** (`ch·9/8`) |
| **6** | **H1 2006 Fit B LO** | `PomH1FitAB(iFit=3)` | **LO** | **singlet** | **≡ 0** |
| 7 | ACTW B, ε = 0.14 | `CTEQ6pdf(iFit=11)` | NLO | grid | grid |
| 8 | ACTW D, ε = 0.14 | `CTEQ6pdf(iFit=12)` | NLO | grid | grid |
| 9 | ACTW SG, ε = 0.14 | `CTEQ6pdf(iFit=13)` | NLO | grid | grid |
| 10 | ACTW D, ε = 0.19 | `CTEQ6pdf(iFit=14)` | NLO | grid | grid |
| 11 | proton rescaled (Angantyr) | `PomHISASD` | — | **needs `xPomNow`** | — |
| 12 | GKG18 Fit A | `LHAGrid1` | **LO** | **singlet** (grid cols −3…3 equal) | ≠ 0 above threshold |
| 13 | GKG18 Fit B | `LHAGrid1` | **LO** | singlet | ≠ 0 above threshold |
| 14 | GKG18 Fit A | `LHAGrid1` | NLO | singlet | ≠ 0 above threshold |
| 15 | GKG18 Fit B | `LHAGrid1` | NLO | singlet | ≠ 0 above threshold |

Two things the tree does not record:

* `PartonDistributions.cc:2622-2631` — H1 2006 Fit A/B set
  `xu = xd = xubar = xdbar = xs = xsbar` and `xc = xcbar = xb = xbbar = 0.`.
  The tree says "the H1 **LO** grids carry no charm or bottom at ANY (β, Q²)"
  (`src/pythia/pythia_bridge.cpp:598-601`,
  `docs/PHYSICS_CHANNELS.md:497`) — true, and **also** true of the two H1 NLO
  fits (3, 4), which the wording does not cover.
* `BeamSetup.cc:1379-1387` — selecting an ACTW set (7–10) makes PYTHIA
  **overwrite** `SigmaDiffractive:PomFlux = 4` and set `PomFluxEpsilon`
  (0.14, or 0.19 for set 10) and `PomFluxAlphaPrime = 0.25`, with a warning.
  Harmless here — the bridge supplies the hard process through `LHAup`
  (strategy 3) and never runs PYTHIA's diffractive machinery — but it is a
  global settings mutation on a run that also hadronises T1 protons, and it is
  not written down anywhere in the tree.

### D4.2 What the code exposes

`include/lipolgen/pythia_bridge.hpp:212` `pom_set` (default 6) and `:271`
`pom_rescale` (default 1.0), bound at `python/bindings.cpp:5011-5015` and
reachable as `--pom-set` / `--pom-rescale` (`python/lipolgen/cli.py:328-337`).
The value is passed through as a raw integer, `src/pythia/pythia_bridge.cpp:329` (anchor read 2026-09-15 "to_string(opt").
No range check — which is fine, because an invalid value fails **hard and
early**: `--pom-set 99` → `RuntimeError: PythiaBridge: pythia.init() failed for
beam id 990` before any event is generated.

### D4.3 The scan, run. And the observable the task proposed is the wrong one

Pilot scan actually executed: **all 15 sets × 20 000 events**, `channel =
coherent`, beam config default, `coherent_t2 = Pomeron`, seed 4242, single
thread. **Whole scan 19.5 s wall, 238 MB RSS.** Cost is not a consideration.
Set 11 has been refused by the bridge constructor since §D4.6's guard, so
what re-runs today is the other **fourteen** (1–10, 12–15).

**The hadronised M_X spectrum cannot be a `PomSet` systematic.** Measured: the
T0 columns `t`, `x_pom`, `q2` and the event `weight` are **bit-identical across
those fourteen sets** (one md5 `ffd35a3a62b591c547e9ca2ac4301b5d`, re-measured
2026-09-05; set 11 gave it too before the guard, and this line said "all 15
sets" until 2026-09-05). The Pomeron PDF enters only the flavour draw
(`src/pythia/pythia_bridge.cpp:560-620`) and PYTHIA's backward evolution; the
diffractive mass, |t|, x_P and the rate are fixed upstream by
`CoherentSampler` / `CoherentXpomModel`. The hadronised ⟨M_X⟩ therefore
reproduces the T0 value by conservation — **4.3478 ± 0.0212 GeV for every
set** (set 5's 4.3204 is the mean over the 19 855 events that survived
hadronisation, not a moved spectrum). A `PomSet` band on M_X, on |t|, on x_P or
on σ is **identically zero by construction**; quoting one would be meaningless.

So the scan must be on the **hadronic final state**, which is the only thing
that moves.

### D4.4 What does move — the pilot scan, 20 000 events per set

`n_ch` = charged final-state hadrons per event; `n_had` = all final-state
particles; ⟨p_T⟩ = mean transverse momentum per final-state particle; fractions
are of all final-state particles. Errors are statistical (√N on the mean).

| set | n_ok | fail | fallback | ⟨n_ch⟩ | ⟨n_had⟩ | ⟨p_T⟩ [GeV] | π | K | p | γ |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 20000 | 0 | 0.00 % | 3.8515 ± 0.0149 | 8.341 | 0.3621 | 0.410 | 0.051 | 0.011 | 0.513 |
| 2 | 20000 | 0 | 0.00 % | **3.7899 ± 0.0147** | 7.944 | 0.3876 | 0.428 | 0.047 | 0.013 | 0.494 |
| 3 | 20000 | 0 | 31.04 % | 3.9186 ± 0.0156 | 8.365 | 0.3710 | 0.406 | 0.066 | 0.011 | 0.502 |
| 4 | 20000 | 0 | 35.38 % | 3.9663 ± 0.0162 | 8.477 | 0.3820 | 0.407 | 0.066 | 0.010 | 0.504 |
| 5 | 19855 | 145 | 0.00 % | **4.3457 ± 0.0192** | 9.366 | **0.3534** | 0.385 | 0.084 | 0.005 | 0.493 |
| **6** | 20000 | 0 | 20.09 % | **3.9407 ± 0.0160** | 8.457 | 0.3718 | 0.404 | 0.066 | 0.010 | 0.505 |
| 7 | 20000 | 0 | 0.00 % | 3.8930 ± 0.0154 | 8.491 | 0.3660 | 0.414 | 0.037 | 0.011 | 0.521 |
| 8 | 20000 | 0 | 0.00 % | 3.9556 ± 0.0165 | 8.311 | 0.4001 | 0.396 | 0.098 | 0.008 | 0.486 |
| 9 | 20000 | 0 | 0.00 % | **3.8635 ± 0.0150** | 8.460 | 0.3583 | 0.414 | **0.036** | 0.011 | 0.523 |
| 10 | 20000 | 0 | 0.00 % | 3.9457 ± 0.0165 | 8.251 | **0.4045** | 0.394 | **0.103** | 0.008 | 0.482 |
| 11 | 20000 | 0 | **100.00 %** | 3.8488 ± 0.0148 | 8.228 | 0.3620 | 0.406 | 0.065 | 0.011 | 0.502 |
| 12 | 20000 | 0 | 18.57 % | 4.0333 ± 0.0168 | 8.617 | 0.3677 | 0.403 | 0.068 | 0.010 | 0.500 |
| 13 | 20000 | 0 | 19.30 % | 4.0292 ± 0.0168 | 8.623 | 0.3677 | 0.402 | 0.068 | 0.010 | 0.501 |
| 14 | 20000 | 0 | 19.98 % | 4.0041 ± 0.0164 | 8.546 | 0.3664 | 0.404 | 0.068 | 0.010 | 0.500 |
| 15 | 20000 | 0 | 20.94 % | 4.0119 ± 0.0165 | 8.560 | 0.3662 | 0.404 | 0.068 | 0.010 | 0.500 |

Seed stability, 20 000 events at seeds 4242 / 777 / 31337:

| set | seed 4242 | seed 777 | seed 31337 |
|---|---|---|---|
| 2 | 3.7899 | 3.8186 | 3.8075 |
| 5 | 4.3457 | 4.3744 | 4.3574 |
| 6 | 3.9407 | 3.9941 | 3.9810 |
| 10 | 3.9457 | 3.9453 | 3.9597 |
| 12 | 4.0333 | 4.0649 | 4.0507 |

Seed-to-seed scatter ≤ 0.053 in ⟨n_ch⟩ against a 0.56 spread across sets, and
the ordering is stable, so **20 000 events per set is already enough** and
there is no case for more.

### D4.5 The scan design, and what the band should be recorded as

**Sets to scan.** The 13 that are real diffractive-PDF fits to H1/ZEUS data:
**3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15**, plus **1** kept as a
Q²-independent sanity floor and clearly labelled as not a fit.

* **Exclude 11 and say why.** `PomHISASD` needs `xPomNow`
  (`PartonDistributions.cc:2809-2823`); the bridge never calls `setXPom`.
  Measured: **100.00 % `n_pom_flavour_fallback` on 20 000 events** — the
  Pomeron PDF is not consulted on a single event, the whole run is the
  charge-democratic e_q² fallback, and **nothing at the run surface says so**.
  Set 11 should be refused by the same kind of guard `PipelineConfig::validate`
  uses elsewhere ("a knob that did not run may not be recorded as if it had"),
  or documented as unusable.
* **Exclude 2** from the systematic band: it is the π⁰ PDF, not a Pomeron fit.

**Observable.** Not M_X (§D4.3, identically zero). The band should be recorded
on the hadronic final state:

* **primary: ⟨n_charged⟩** — the cleanest single number, it is what a tracker
  sees, and its 0.4 % statistical error at 20 000 events resolves the spread by
  a factor 35;
* **secondary: the HFS strangeness fraction** — the *largest* mover, and the
  one most likely to matter for a PID-based analysis;
* **reported alongside:** ⟨p_T⟩ per particle and ⟨n_had⟩.

**The band, measured (13 fit sets, default = set 6):**

| observable | default (set 6) | min (set) | max (set) | band about the default |
|---|---|---|---|---|
| ⟨n_charged⟩ | 3.941 | 3.863 (9) | 4.346 (5) | **−2.0 % / +10.3 %** |
| ⟨n_had⟩ | 8.457 | 8.251 (10) | 9.366 (5) | −2.4 % / +10.8 % |
| ⟨p_T⟩ [GeV] | 0.3718 | 0.3534 (5) | 0.4045 (10) | −4.9 % / +8.8 % |
| K fraction of HFS | 0.066 | 0.036 (9) | 0.103 (10) | **−45 % / +56 %**, a **factor 2.9** end to end |

**How to record it.** As an **envelope over re-runs, one npz per set** — *not*
as a per-event reweighting. The set changes the final state event by event and
there is no weight that maps one set onto another. That makes §D4.7 a
precondition: without `pom_set` in the `meta`, an envelope built from 13 npz
files has nothing in the files themselves that says which is which.

Note for whoever writes the band down: **set 6 is the only LO H1 set**, so the
band mixes LO and NLO DPDFs used in an LO Monte Carlo. That is a defensible
choice for a *systematic* envelope and an indefensible one for a central value;
say which is which.

### D4.6 The charm finding — narrower than the tree states, and the fallback is not a systematic

1. **Confirmed at the PYTHIA source**, not merely observed: `xc = xcbar = 0.`
   on H1 2006 Fit A/B (`PartonDistributions.cc:2630`), i.e. sets 3, 4 **and** 6.
   The tree's wording covers the LO one only.
2. **But it is not a property of "the Pomeron grid" in general.** Measured by
   toggling `include_charm` and diffing the final states event by event
   (6 000 events, seed 4242, "changed" = the event's pid list differs):

   | set | events changed by `include_charm = false` |
   |---|---|
   | 5 (H1 2007 Jets) | **29.40 %** |
   | 12 / 13 / 14 / 15 (GKG18) | 10.97 / 11.23 / 11.20 / 11.53 % |
   | 3, 4, 6, 7, 8, 9, 10 | **0.00 %** |

   So the H1 2007 Jets set puts nearly **30 %** of a coherent sample on a
   charm initiator, and the ACTW grids give charm no support in this
   generator's (β, Q²) window even though `CTEQ6pdf` carries a charm slot.
   Any scan that includes set 5 will see an open-charm final state that the
   default never produces — that is a *feature* of the band, but it must be
   flagged, because it is the same signature the `2e404b8` bug faked.
3. **Is the light-only e_q² fallback a systematic in its own right?
   Measured: no — it costs exactly zero.** On set 6, `q2_pdf_min` = 1.0
   (20.09 % fallback) against 1.75 or 3.0 (0.00 % fallback) produces
   **bit-identical** final states — same `pid` array, same `p4` array, 33 689
   hadrons over 4 000 events. Same result for set 12 (18.57 % → 2.48 %
   fallback) and set 13 (19.30 % → 2.83 %): **bit-identical**.

   The reason is in the source, and it is structural rather than lucky: **every
   Pomeron DPDF in PYTHIA carries a single light-quark singlet** — H1 sets it
   explicitly (`PartonDistributions.cc:2624-2629`), and the GKG18 LHAGrid1 files
   have columns −3, −2, −1, 1, 2, 3 numerically equal row by row. So
   e_q²·xf_q ∝ e_q² **exactly** over the light flavours, and the normalised
   "fallback" *is* the true draw. And the fallback fires only where every weight
   vanishes, which is below the charm threshold, where charm is zero anyway —
   so the light-only restriction has never yet discarded anything.

   Consequence for the docs: `docs/PHYSICS_CHANNELS.md:450` (as of adec442)'s remedy —
   *"raising `q2_pdf_min` to ≈1.75 removes it at the cost of clamping every
   flavour weight to that Q²"* — is true about the **counter** and false about
   the **cost**. It changes nothing at all. The light-only restriction should
   stay (it is the correct guard against the `2e404b8` bug), and the tree should
   stop implying that the fallback is a modelling uncertainty. Measured, it is
   a bookkeeping artefact.

### D4.7 `pom_set` is absent from the npz `meta` too

Second instance of the D3.5 defect. The coherent npz `meta` carries **no**
`pom_set`, **no** `pom_rescale`, **no** `coherent_t2`, and **no**
`n_pom_flavour_fallback`. Two runs differing only in `--pom-set` have identical
`meta` while their entire hadronic final state differs — which is precisely the
band this item exists to establish.

And `docs/PHYSICS_CHANNELS.md:450` (as of adec442) says of the fallback *"Counted, so the
fallback share of a run is visible."* It is counted
(`include/lipolgen/pythia_bridge.hpp:341` (anchor read 2026-09-15 "n_pom_flavour_fallback = 0"), bound at
`python/bindings.cpp:5053` (anchor read 2026-09-15 `def_readonly("n_pom_flavour_fallback"`)) and **surfaced nowhere**: not in the `meta`, and
not in the run banner, which prints `n_ok`, `n_failed`, `n_retries` only
(`python/lipolgen/cli.py:712-715` (as of a7b3d18)). On set 11 that is the difference between a
run that used a Pomeron PDF and one that used none.

### D4.8 Reproduction

```bash
cd /home/cpeng/Projects/polli/LiPolGen && source env.sh
# per set: make_config(channel="coherent", events=20000, seed=4242);
# PythiaBridgeOptions with coherent_t2 = Pomeron and pom_set = N;
# set_pythia_hadronizer; Pipeline(...).generate(0, True, 1);
# lg.hfs_sample(cols["events"], ...) -> offsets/pid/p4/charge.
# Charm probe: the same, twice, with include_charm True/False, diffing pid per event.
# Fallback probe: the same, with q2_pdf_min 1.0 vs 1.75, diffing pid and p4 wholesale.
```

---

## D5 — the |t| ceiling

### D5.1 The arithmetic, reproduced

`src/core/coherent.cpp:179-202` (anchor read 2026-09-15):

    c₂(|t|, P_zz) = −(P_zz/2)·eps_b0·B·|t| + amp·P_zz
    positivity_margin(t_max, P_zz) = 1 − |c₂(t_max, P_zz)|      (:202)

Write c₂ = A·|t| + C with A = −(P_zz/2)·eps_b0·B and C = amp·P_zz. The
positivity edge |c₂| = 1 then has the closed form

    |t|_pos = (1 − sign(A)·C) / |A|

and for the shipped signs (eps_b0 < 0, amp > 0, so A and C always share a sign)
this is

    |t|_pos = 2·(1/|P_zz| − amp) / (|eps_b0|·B).

At B = 50, amp = 0.01 this reproduces `phase_C_numbers.md` §C4.5's table
**exactly**:

| eps_b0 | \|t\|_pos, P_zz = −2 | \|t\|_pos, P_zz = +1 |
|---|---|---|
| −0.08 (shipped, exact input) | **0.2450** | **0.4950** |
| −0.0527 (α+d model, rounded) | 0.3719 | 0.7514 |
| −0.0171 (GFMC, rounded) | 1.1462 | 2.3158 |
| −0.0070 (measured Q, rounded) | **2.80** | 5.6571 |

The last three `eps_b0` are ROUNDED and the rows are arithmetic on the rounded
values, so the bottom row is quotable as **2.80 GeV²** and not as "2.8000". On
the DERIVED band (measured 2026-09-05, B = 50, amp = 0.01) the `eps_b0` are
−0.0526846 / −0.0171207 / −0.0070024 and the P_zz = −2 edges 0.3720 / 1.1448 /
**2.7990**.

**Phase C's conclusion is confirmed on the arithmetic: the POSITIVITY EDGE is
a consequence of the oversized `eps_b0`, not a property of ⁶Li** — 0.245 GeV²
at the shipped value against 2.80 at the measured quadrupole. *(This sentence
read "`COHERENT_T_MAX_DEFAULT` = 0.2 is a consequence of the oversized eps_b0"
until 2026-09-05. **Superseded by D5 itself** — the section this line closes:
the ceiling stays 0.2 on the ANCHOR RANGE and does not move with `eps_b0`;
`STATUS.md` decision row 12, `OPEN_ITEMS_SOLUTIONS.md` §11.6.)*

### D5.2 "c₂ is linear and unbounded" — and whether a bounded form is motivated or a fudge

`coherent.hpp:496` is exact, and the header's **own definition** three hundred
lines earlier already contains the answer. `coherent.hpp:376`:

    |F_m(|t|, Φ)|² = exp(−|t|·[B + ΔB_m·cos 2(Φ − Φ_S)])

That is **positive for every |t| and every ΔB_m**. The linear c₂ is the
O(ΔB·|t|) truncation of this exponential's azimuthal Fourier expansion, and the
**truncation** is what goes negative — not the model. Exactly, with
z_m = ΔB_m·|t| and p_m the substate populations:

    W(φ) = Σ_m p_m e^{−z_m cos 2φ} / (2π Σ_m p_m I₀(z_m))          [a density, ∀|t|]
    c₂^exact(|t|) = −2 Σ_m p_m I₁(z_m) / Σ_m p_m I₀(z_m)  +  amp·P_zz

and as |t| → 0, I₁(z)/I₀(z) → z/2 and Σ_m p_m ΔB_m = (eps_b0·B/2)·P_zz, so
c₂^exact → −(P_zz/2)·eps_b0·B·|t| + amp·P_zz — the shipped linear form, recovered
analytically as the |t| → 0 limit.

**So a bounded form is physically motivated — but the motivated bounded form is
not a form factor bolted onto the modulation. It is the un-truncated
exponential the header already calls "THE definition".** A form factor on the
tensor modulation would be a fudge on two counts: nothing in the slope-splitting
mechanism supplies one, and it would break the ΔB ↔ quadrupole map that
`cluster_config.hpp`'s `quadrupole_from_a2_slope` inverts — the map §C4 used to
derive the whole eps_b0 band.

**Measured** (B = 50, amp = 0.01, p₀ = (1 − P_zz)/3, p_± = (2 + P_zz)/6;
`scipy.special.i0e`-stabilised, 60 001-point φ grid):

* The **pure deformation** density (amp = 0) never turns negative: min_φ W =
  2.85e−34 at |t| = 10 and underflows to 0 at |t| = 100 for eps_b0 = −0.08,
  P_zz = −2. It decays to 0⁺ and never crosses.
* With the flat `amp` term added, the **exact** edge sits at

  | eps_b0 | exact, P_zz = −2 | exact, P_zz = +1 | ratio to the linear edge |
  |---|---|---|---|
  | −0.08 | 0.6562 | 1.5055 | 2.678 / 3.041 |
  | −0.0527 | 0.9961 | 2.2854 | 2.678 / 3.041 |
  | −0.0171 | 3.0700 | 7.0432 | 2.678 / 3.041 |
  | −0.0070 | 7.4996 | 17.2056 | 2.678 / 3.041 |

  The ratio is eps_b0-independent because both edges scale as 1/|eps_b0|.
* At the point where the truncated form hits |c₂_lin| = 1, the exact
  coefficient is only |c₂| = 0.899 (P_zz = −2) / 0.896 (P_zz = +1), and the
  exact density is still min_φ W = 0.279 / 0.285 of its mean. **The shipped
  ceiling is set where the approximation fails, at a point where the physics is
  nowhere near failing.**

### D5.3 How much rate is at stake — the number that decides the question

Measured on **200 000 generated coherent events** at the shipped defaults
(`--channel coherent`, seed 99, **one P_zz = −2 category**, 4 threads):

    ⟨|t|⟩ = 0.019974 GeV²      max |t| = 0.197575 (the ceiling)
    frac |t| > 0.05  = 0.081960   (untruncated e^{−B|t|} = 8.208e−02)
    frac |t| > 0.10  = 0.006655   (6.738e−03)
    frac |t| > 0.15  = 0.000550   (5.531e−04)
    frac |t| > 0.20  = 0          (4.540e−05)

**Reconciled with `phase_D_numbers.md` §D5.3, which quotes 0.019963 for the
same quantity** (and which is the number `coherent.hpp`, `docs/USAGE.md` §4
and `OPEN_ITEMS_SOLUTIONS.md` §11.6 all carry): the two are DIFFERENT RUNS of
the same estimator, not a discrepancy. §D5.3 there uses the three-category
`tensor-thirds` plan at P_z = 0.7 / P_zz = 0.6; this block uses the
single-category P_zz = −2 plan. The run plan decides how the 200 000 events
are split across spin categories, so the two runs are not the same 200 000
|t| draws even at one seed, and the means differ in the fourth digit.
Re-measured 2026-09-05, both at seed 99 / 200 000 events / 4 threads:

| plan | ⟨\|t\|⟩ | max \|t\| | >0.05 | >0.10 | >0.15 | >0.20 |
|---|---|---|---|---|---|---|
| `tensor-thirds` (0.7, 0.6) — §D5.3 | 0.019963 | 0.197575 | 0.082015 | 0.006550 | 0.000455 | 0 |
| one P_zz = −2 category — here | 0.019974 | 0.197575 | 0.081960 | 0.006655 | 0.000550 | 0 |

**`max |t| = 0.197575` and the `> 0.20` count of 0 are identical in both**,
which is the part the ceiling argument rests on. Quote 0.019963 with the
`tensor-thirds` window, since that is the one every shipped docstring cites.

`sample_t` (`src/core/coherent.cpp:140-144`) **renormalises** on [0, t_max], so
the ceiling does not lose rate, it redistributes **4.5e−5** of it. Moving the
ceiling from 0.2 to 0.245, or to the exact-form 0.657, changes the generated
sample at the 1e−5 level, and the RP-tagged sample sits at
⟨|t|⟩ = pT_cut² + 1/B, far below it either way.

**That is the decisive fact.** The ceiling is not a rate question and never
was. It is a statement about **where the model is defined**, and it should be
justified as one.

### D5.4 The honest replacement — three options, priced

**Option A — derive the ceiling from positivity at the actual eps_b0.**
Add `CoherentScenario::t_positivity_edge(pzz)` returning the §D5.1 closed form
and make the ceiling `min(that, the input-range ceiling)`.
*For:* it moves with the author's eventual eps_b0 decision automatically, which
is what the task asks for. *Against, decisively:* at the measured
eps_b0 = −0.0070 it returns **2.80 GeV²**, which is **9.3×** outside the
Mantysaari digitisation the whole deformation term is scaled from
(`src/core/coherent.cpp:58-64`: four rows, |t| = 0.05, 0.10, 0.20, 0.30). A
positivity-only ceiling would license extrapolating a linear-in-|t| fit ten
times past its data — it is *less* honest than the number it replaces, not more.
*Cost:* `coherent_t_max` is a plain double read once (`pipeline.hpp:555`,
`pipeline.cpp:1016`); making it derived would move a shipped default the moment
eps_b0 moves, so it would have to be opt-in in any case.

**Option B — keep 0.2 fixed, and replace the stated reason.**
The honest "why" is **not** positivity: it is the **anchor range**.
`mantysaari_a2_deuteron()` carries four digitised rows spanning |t| ≤ 0.30, and
`coherent.hpp:45-50` already gives this as the **first** of its two reasons
("digitized over |t| <= 0.30 and LINEAR in |t| only as |t| -> 0, so 0.5 is
outside the input"). That reason is a property of the **input** and does not
move with eps_b0. The positivity reason is the artefact.

**Option C — un-truncate the modulation (the physics fix).**
Sample φ from W(φ) directly instead of from 1 + c₂ cos 2φ.
`src/core/coherent.cpp:361-369` already does rejection sampling, so only the
envelope changes (max_φ W is closed-form: it is the φ = 0 or φ = π/2 value).
Positivity then holds for every |t| and every P_zz in range, and the only
ceiling left is the anchor range. *Costs:* `CoherentEvent::c2`
(`coherent.hpp:620`) stops being the sampled weight's coefficient and becomes
its |t| → 0 limit, which is a documented output; and **every reference file
that pins a sampled azimuth moves**, which this run may not do.

### D5.5 The author decision, stated precisely

**Recommended: Option B's reason, Option A's machinery, Option C filed behind
the eps_b0 decision.** Concretely, five things:

1. **`COHERENT_T_MAX_DEFAULT` stays 0.2.** No reference gate moves; the rate at
   stake is 4.5e−5 (§D5.3); and 0.2 is inside the anchor range under every
   eps_b0 in the ⁶Li band.
2. **Rewrite `coherent.hpp:43-51` so the primary reason is the anchor's
   |t| ≤ 0.30 digitisation range**, and the positivity reason is stated as
   *secondary and contingent*, with the arithmetic on the page: the positivity
   edge is 0.245 GeV² **only because eps_b0 is 11.4× oversized**, and would be
   2.80 GeV² at the measured quadrupole — i.e. **positivity stops binding the
   moment eps_b0 is corrected, and the anchor range does not.** Today the
   header presents the two reasons as agreeing ("Two reasons, and they agree"),
   which is true at the shipped eps_b0 and false at every ⁶Li value of it.
3. **Keep `check_positivity()` exactly as it is** — it must still throw. It is
   not the derivation of the ceiling; it is the **guard** that catches an author
   who raises `t_max` or `eps_b0` past where the *truncated* weight the sampler
   actually uses stops being a density. Both failure modes are real and the
   guard is the only thing that catches them.
4. **Add `CoherentScenario::t_positivity_edge(double pzz)`** returning
   (1 − sign(A)·C)/|A|, so the number is *derived* rather than repeated in four
   docstrings and one open-items table, and so it **moves with eps_b0** for any
   caller who asks. Changes no default; one new pinned test value.
5. **Record `coherent_t_max` in the npz `meta`** (§D5.6).

**Why not Option C now — and the number that says when to take it.** It is the
physically right answer, and it is the one to take **when eps_b0 is corrected**.
Measured, max |c₂^exact − c₂^lin| over |t| ∈ (0, 0.30] (the anchor range):

| eps_b0 | P_zz = −2: abs / rel | P_zz = +1: abs / rel |
|---|---|---|
| −0.08 (shipped) | 0.1744 / **14.30 %** | 0.0255 / 4.18 % |
| −0.0527 (α+d) | 0.0559 / 6.90 % | 0.0075 / 1.86 % |
| −0.0171 (GFMC) | 0.0021 / 0.75 % | 0.0003 / 0.19 % |
| −0.0070 (measured Q) | 0.00014 / **0.12 %** | 0.00002 / 0.03 % |

So at the shipped eps_b0 the truncation is already a **14 %** error on c₂
*inside* the anchor range — Option C would change published numbers — while at
the measured quadrupole it is 0.12 % and Option C would change nothing
observable. **The truncation stops mattering exactly when the ceiling stops
mattering**, which is why the two are one decision: it moves every pinned
coherent azimuth, so file it as the follow-on to the eps_b0 decision, not as an
independent item.

### D5.6 `coherent_t_max` is absent from the npz `meta`, and from the CLI

Third instance of the D3.5 defect, and the worst of the three: `coherent_t_max`
is on `PipelineConfig` (`include/lipolgen/pipeline.hpp:873`) and bound
(`python/bindings.cpp:3864` (as of adec442)), but it is **not reachable from `cli.py` at all**
(there is no `--coherent-t-max`) and **not in the `meta`**. A Python caller who
moves it changes the entire |t| spectrum, the tag acceptance and every c₂ in
the file, and the sidecar cannot tell.

### D5.7 Reproduction

```bash
cd /home/cpeng/Projects/polli/LiPolGen && source env.sh
python -c "import lipolgen as lg; ..."     # 200k coherent events, seed 99, column 't'
# The exact/linear comparison is pure arithmetic on B, amp, eps_b0, p_m -- no
# library call -- using scipy.special.i0e / i1 for the Bessel ratio.
```

---

## D6 — one defect, three instances

| knob | bound | in `cli.py` | in npz `meta` | pytest |
|---|---|---|---|---|
| `fsi` (off / cluster / nucleon) | yes | yes | **no** | config only |
| `fsi_sigma_mb` (the 20–40 band) | yes | yes | **no** | config only |
| `pom_set` | yes | yes | **no** | — |
| `pom_rescale` | yes | yes | **no** | — |
| `coherent_t2` | yes | yes | **no** | — |
| `coherent_t_max` | yes | **no** | **no** | — |
| `n_pom_flavour_fallback` (a counter, not a knob) | yes | **not printed** | **no** | — |
| `b1_model`, `b1_band_scale`, `b1_alpha_d_dwave_weight`, `b1_unpol` | yes | yes | **yes, unconditional** | yes |
| the 19 `rc_*` keys | yes | yes | **yes, conditional** | yes |

The rule is already written down in the code, at `python/bindings.cpp:462-466` (as of adec442)
and `:474-479`: *without these keys, otherwise identical npz files are
indistinguishable*. It was applied to `b1_*` and to `rc_*` and to nothing else.
Measured consequences, all three from this task:

* three FSI runs of the same seed: `meta` identical, total rate differs by 48 %;
* fifteen `pom_set` runs: `meta` identical, the whole hadronic final state
  differs and the K fraction by a factor 2.9;
* set 11: `meta` identical to a real run, and 100 % of events did not use a
  Pomeron PDF at all — the exact case `PipelineConfig::validate`'s "a knob that
  did not run may not be recorded as if it had" rule exists to prevent.

**Recommend one Phase-E item**: extend the `meta` writer to the FSI, Pomeron and
coherent-|t| knobs, following the **`rc` precedent** (a conditional block, so
every current default stays byte-identical and no reference gate moves), print
`n_pom_flavour_fallback` in the run banner beside `n_ok`/`n_failed`, add
`--coherent-t-max`, and gate each with a pytest that asserts two runs differing
only in that knob are distinguishable in `meta`.
