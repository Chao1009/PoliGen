# The spin-3/2 inclusive sector at finite γ

**Status:** theory note, LiPolGen `docs/theory/`. Written to close
`docs/open_items/physics_literature.md` §1 ("Spin-3/2 inclusive
structure-function basis (rank-2 and rank-3)") and `docs/OPEN_ITEMS_SOLUTIONS.md`
row 8. Audience: the LiPolGen authors and a theory co-author.
No code is changed by this note; every statement about the generator carries a
`file:line` that was read to write it. Those references were **re-resolved by
hand on 2026-09-03**, after the Phase E code comments shifted several of them.
**As of 2026-09-05 this document IS machine-checked**:
`validation/check_physics_channels_links.py` covers it (alongside
`docs/PHYSICS_CHANNELS.md` and `docs/PYTHIA_BRIDGE.md`), strictly, with a
sha256 fingerprint on every cited range so a code edit that moves or changes a
cited block is caught rather than silently going stale the way the 15 ranges
above did between 2026-09-03 and 2026-09-05. `python/tests/test_doc_link_gate.py`
runs the gate over the real document in the suite, so a stale citation here
fails `pytest`, not just a hand run of the script.

Conventions of this note: `M` is the free nucleon mass (per-nucleon Bjorken `x`,
`include/lipolgen/asymmetries.hpp:57-61`); `J` is the ion spin; `m` a projection
on the fill axis `n̂(θ_S, φ_S)`; `φ' = φ − φ_S`; `λ_e = ±1` the electron helicity
sign and `P_e` its magnitude (`include/lipolgen/xsec.hpp:121-128`). Papers cited
as `[n]` are listed in §7 with equation numbers.

---

## 1. Scope and result

Today the generator treats a J = 3/2 ion (⁷Li) with **rank 0, 1 and 2 only**: the
unpolarized rate from `F1`/`F2`, a vector term `λ_e P_e (m/J) cos θ_S · A_∥`, and
a rank-2 (alignment) block that reuses the spin-1 machinery verbatim — the same
scalar `Q_NN(m) = [3m² − J(J+1)]/3` geometry (`include/lipolgen/xsec.hpp:23-29`,
`src/core/xsec.cpp:168-173`), the same `b1 + (1−y)/(xy²) b2` kernel and the same
optional exact finite-γ Cosyn kernel (`include/lipolgen/xsec.hpp:204-217`), fed
through the dedicated J = 3/2 slots `b1_32_func / b2_32_func / delta_32_func`
(`include/lipolgen/xsec.hpp:218-226`, dispatched at `src/core/xsec.cpp:108-118`); the
rank-3 (octupole) sector is computed in the spin bookkeeping
(`include/lipolgen/spin.hpp:94-102`, `src/core/spin.cpp:246-259`) but never
reaches the cross section. This note establishes two things. **(a)** The
rank-≤2 truncation is *exact*, not an approximation, for every inclusive
observable measured with an **unpolarized lepton beam** — A_zz, the cos 2φ
gluon-transversity amplitude and the rate itself — because parity and
time-reversal invariance of the one-photon-exchange hadronic tensor confine
odd-rank target multipoles to the *antisymmetric* part of W^{μν}, which is
contracted only with the λ_e-odd part of the lepton tensor (§4). The corollary
is sharper than the register's statement and cuts the other way as well: the
truncation is *not* exact for A_∥ — nor for any other λ_e-odd observable —
on a J = 3/2 target, because the rank-3 multipole sits in the same
beam-helicity channel as g1 and the code's default ⁷Li spin-temperature fill
carries R₃ ≠ 0 (§2, §6). **(b)** The finite-γ
decomposition of the J = 3/2 inclusive cross section: the rank-2 block is the
spin-1 block of [2] *unchanged in form*, so `cosyn_tensor_sfs` and
`tensor_harmonics_gamma` apply to ⁷Li verbatim with only `Q_NN` and the b-values
changed; the new rank-3 block contributes exactly **two** structure functions at
finite Q², both λ_e-odd, collapsing to a single leading-twist function at γ → 0
(§5).

---

## 2. The spin-3/2 density matrix and the code's parameters

### 2.1 Multipole content

A spin-J density matrix decomposes into multipoles l = 0, …, 2J of the rotation
group. Going from spin j − 1/2 to spin j adds exactly one new multipole, l = 2j,
carrying 2(2j)+1 parameters, collected in a symmetric traceless rank-2j Lorentz
tensor — *regular* for integer j, *pseudo*-tensor for half-integer j, orthogonal
to the target four-momentum in every index ([1] App. D, opening paragraph). For
J = 3/2:

        l = 0   1 parameter    (normalization)
        l = 1   3 parameters   vector       S^μ        (pseudovector)         (1)
        l = 2   5 parameters   alignment    T^{μν}     (regular tensor)
        l = 3   7 parameters   octupole     R^{μνρ}    (pseudotensor)

3 + 5 + 7 = 15 = (2J+1)² − 1 real parameters, as required.

The covariant form is [3] Eqs. (9)–(13): orthogonality

        P_μ S^μ = 0 ,   P_μ T^{μν} = 0 ,   P_μ R^{μνρ} = 0 ,                   (2)

([3] Eq. (9)) and the light-cone decomposition of S^μ, T^{μν}, R^{μνρ} into the
15 real parameters ([3] Eqs. (11), (12), (13)),

        l = 1 :  S_L , S_T^i                                (1 + 2)
        l = 2 :  S_LL , S_LT^i , S_TT^{ij}                  (1 + 2 + 2)        (3)
        l = 3 :  S_LLL , S_LLT^i , S_LTT^{ij} , S_TTT^{ijk} (1 + 2 + 2 + 2)

with the transverse metric g_T^{μν} = g^{μν} − n̄^μ n^ν − n^μ n̄^ν ([3] Eq. (14)).
The subscript pattern is the L/T index count; the number of T indices is |m_l|,
the multipole's azimuthal index. This is the labelling used throughout §4 and §5.

The three-dimensional (rest-frame) form is [4] Eq. (21):

        ρ = (1/4) [ 1 + (4/5) S^i Σ^i + (2/3) T^{ij} Σ^{ij}
                                      + (8/9) R^{ijk} Σ^{ijk} ] ,              (4)

with Σ^i, Σ^{ij}, Σ^{ijk} the rank-1/2/3 polarization operator basis for spin
3/2. Equation (4) is transcribed exactly as printed in [4]; the note does not
depend on the values 4/5, 2/3, 8/9 — every result below is anchored on the
*observable* per-state quark densities of §3, which are normalization-free.

### 2.2 What the generator actually stores

LiPolGen never stores ρ as 15 parameters. It stores an **axially symmetric fill**:
populations p_m along one axis n̂(θ_S, φ_S), m ordered +J … −J
(`include/lipolgen/spin.hpp:9-14`, `include/lipolgen/bookkeeping.hpp:59-68`).
An axially symmetric state has one non-zero parameter per multipole, so the 15
parameters collapse to three numbers plus two angles. The three numbers are
(`include/lipolgen/spin.hpp:18-21`, implemented at `src/core/spin.cpp:261-290`):

        P_z = ⟨J_z⟩ / J                                    vector, l = 1        (5)
        T   = ⟨3J_z² − J(J+1)⟩ / 3                         alignment, l = 2     (6)
        R₃  = ⟨J_z³ − (41/20) J_z⟩ / (3/10)                octupole, l = 3      (7)

all with respect to the fill axis, all dimensionless. Equation (6) is
`tensor_polarization` (`src/core/spin.cpp:230-244`, the `j = 1.5` branch divides
by 3); Eq. (7) is `octupole_moment` (`src/core/spin.cpp:246-259`, the operator
`J_z³ − (41/20)J_z` divided by 0.3). The operator in (7) is the standard rank-3
spherical tensor for J = 3/2: T_30 ∝ 5J_z³ − [3J(J+1) − 1] J_z = 5J_z³ −
(41/4)J_z, i.e. J_z³ − (41/20)J_z up to a factor 5.

Per-state values of the three weights, m = (3/2, 1/2, −1/2, −3/2):

        m/J           =  ( 1,   1/3,  −1/3,  −1 )
        Q_NN(m) ≡ T   =  ( 1,  −1,    −1,     1 )                              (8)
        o(m)          =  ( 0.3, −0.9,  0.9,  −0.3 )   [ = m³ − (41/20)m ]
        R₃(m) = o/0.3 =  ( 1,  −3,     3,    −1 )

Two normalization facts follow and matter later. **(i)** T is normalized so
that the stretched state |3/2, ±3/2⟩ has |T| = 1; its full physical range is
[−1, +1]. **(ii)** R₃ is normalized on the *stretched* state, R₃(+3/2) = +1, but
its physical range is [−3, +3], because the |m| = 1/2 states carry R₃ = ∓3. This
is a deliberate choice, pinned by
`tests/test_spin.cpp:138-141` ("pure m = +3/2 has all normalized moments = +1",
`spin32_populations(1,1,1) = (1,0,0,0)`); it is *not* a moment normalized to
unit range, and any positivity check written against |R₃| ≤ 1 would be wrong.

The J = 1 limit is a consistency check and is pinned in the code. For J = 1 the
code uses P_zz = ⟨3J_z² − 2⟩ with no division by 3
(`src/core/spin.cpp:269-275`), the Hoodbhoy–Jaffe–Manohar c_m = 3m² − 2 = (1, −2,
1) convention (`include/lipolgen/spin.hpp:23`). The cross-section kernel then
divides by 3 itself, `Q_NN = (3m² − J(J+1))/3 = c_m/3`
(`src/core/xsec.cpp:168-173`), giving **one** rank-2 geometry for both spins
(`include/lipolgen/xsec.hpp:23-29`, pinned at `tests/test_xsec.cpp:297-316` and
against the HJM transcription digit-for-digit at `tests/test_xsec.cpp:318-335`).
So Q_NN is the same object for J = 1 and J = 3/2 — but T is not: T ≡ P_zz for
J = 1 and T ≡ Q_NN for J = 3/2. That asymmetry is a live trap for anyone reading
`AxisMoments::tensor` across the two isotopes, and the two example CLIs are
already caught in it. For the hand-built `--plan pure` fill (the stretched
state, §2.4) `examples/generate_tagged.cpp:123` recorded `pzz_true = 1.0`, while
`examples/generate_full.cpp` recorded `pzz_true = (spin1 ? 1.0 : 0.0)`, i.e.
**0 for ⁷Li** — although `moments_along_axis(1.5, {1,0,0,0}).tensor = +1`. The
two CLIs therefore disagreed about the recorded rank-2 moment of the *same*
fill, and `generate_full` was the wrong one. **Fixed 2026-09-03 (commit
0961c60): both CLIs now write `j >= 1.0 ? 1.0 : 0.0`** — `generate_tagged.cpp:123`
and `examples/generate_full.cpp:135` — so a J = 3/2 pure fill records T = 1. No
cross section moved (the plan runs at `lam_e = 0` and the kernel reads
`state.m`, not the plan's bookkeeping), but any analysis that divides by
`measured_pzz` did.

The exact map from moments to populations is a 4×4 linear solve,
`spin32_populations(pz, t, o)` (`src/core/spin.cpp:306-324`), whose rows are
exactly the four weight vectors (1, m/J, Q_NN, R₃) of (8)
(`src/core/spin.cpp:308-313`), round-tripped at `tests/test_spin.cpp:129-137`.

### 2.3 The covariant tensors of an axially symmetric fill

For a fill with axis n̂ and moments (P_z, T, R₃), the rest-frame multipole tensors
are the unique symmetric traceless tensors built from n̂:

        S^i   = P_z n_i                                                        (9)
        t_ij  = (T/2) (3 n_i n_j − δ_ij)                                      (10)
        R_ijk = (R₃/2) (5 n_i n_j n_k − n_i δ_jk − n_j δ_ik − n_k δ_ij)       (11)

normalized so that with n̂ = ẑ one has S^z = P_z, t_zz = T, R_zzz = R₃.
Equation (10) is exactly the code's stated form
(`include/lipolgen/xsec.hpp:24-26`; line 23 carries the t_ij form itself, 24-25
the Q_NN definition and t_geo/c_eff) and is [2] Eq. (9) with the pure-state
scalar replaced by T: [2] Eq. (9) reads t^{αβ}(N, Λ_d) = (1/6)(g^{αβ} −
p^α p^β/p² + 3N^α N^β) × (+1 for Λ_d = ±1, −2 for Λ_d = 0), whose rest-frame
spatial part is (c_Λ/6)(3n_i n_j − δ_ij) with c_Λ/3 = Q_NN. Equation (11) is the
rank-3 continuation of the same construction; it is *not* taken from a source,
it is the unique normalization consistent with (10) and with R₃ of (7).

The projections that the cross section sees are the L/T components of these
tensors in the **virtual-photon** frame, following [2] Eqs. (14a)–(14c):

        T_LL = t_zz ,   T_LT cos φ_TL = t_xz ,   T_TT cos 2φ_TT = t_xx − t_yy  (12)

and, by the same rule applied to (11),

        R_LLL = R_zzz ,   R_LLT cos φ_LLT = R_xzz ,
        R_LTT cos 2φ_LTT = R_xxz − R_yyz ,   R_TTT cos 3φ_TTT = R_xxx − 3R_xyy (13)

with z along q. Section 4 shows that only the first two lines of (13) survive in
inclusive DIS.

### 2.4 ⁷Li fill states as worked examples

The **library-level** plan builders produce exactly one J = 3/2 plan,
`helicity_flip_plan(1.5, …)` (`include/lipolgen/bookkeeping.hpp:151-152`,
`src/core/bookkeeping.cpp:81-110`); `tensor_thirds_plan`,
`transverse_tensor_plan` and `tensor_flip_plan` all hard-code spin 1
(`src/core/bookkeeping.cpp:107-140`). The **example CLIs add one more by hand**:
`--plan pure` builds `pops = (1, 0, 0, 0)`, the stretched state |3/2, +3/2⟩,
with an *unpolarized* electron at θ_S = 0
(`examples/generate_full.cpp:128-134`, `examples/generate_tagged.cpp:119-124`).
That fill carries (P_z, T, R₃) = (1, +1, +1) — the **largest rank-3 moment any
J = 3/2 fill can carry** — so the programme does reach the extreme of the
octupole range today. It is protected by the theorem of §4 (`lam_e = 0`), so no
number it has produced moves; but it is the reason the rank-3 sector cannot be
waved away with "no run plan reaches it". Returning to `helicity_flip_plan`:
its default branch
(`use_explicit_pzz = false`, `include/lipolgen/bookkeeping.hpp:133-138`) is the
spin-temperature (maximum-entropy) ladder p_m ∝ u^m
(`src/core/bookkeeping.cpp:85`, `include/lipolgen/bookkeeping.hpp:198-208`). For
J = 3/2 with ladder ratio u the three moments have closed forms — derived here
and checked numerically against `populations_maxent`:

        P_z = (u−1)(3u² + 4u + 3) / [ 3(u+1)(u²+1) ]
        T   = (u−1)² / (u²+1)                                                 (14)
        R₃  = (u−1)³ / [ (u+1)(u²+1) ]

so that R₃ = T · (u−1)/(u+1) — the octupole is the alignment times tanh(β/2).
Worked values (u from bisection, moments from `moments_along_axis`):

| fill | p(+3/2), p(+1/2), p(−1/2), p(−3/2) | P_z | T | R₃ |
|---|---|---|---|---|
| pure \|3/2, +3/2⟩ — the CLIs' `--plan pure` | 1, 0, 0, 0 | **1** | **+1** | **+1** |
| pure \|3/2, +1/2⟩ | 0, 1, 0, 0 | **1/3** | **−1** | **−3** |
| unpolarized | ¼ each | 0 | 0 | 0 |
| spin-temperature, P_z = 0.30 | .40214 .27664 .19031 .13092 | 0.30 | 0.0661 | 0.0122 |
| spin-temperature, P_z = 0.50 | .52581 .26799 .13659 .06961 | 0.50 | 0.1909 | 0.0620 |
| spin-temperature, P_z = 0.70 (u = 3) | 27, 9, 3, 1 over 40 | **0.7** | **0.4** | **0.2** |
| spin-temperature, P_z = 0.80 | .76411 .18213 .04341 .01035 | 0.80 | 0.5489 | 0.3376 |
| spin-temperature, P_z = 0.90 | .86892 .11412 .01499 .00197 | 0.90 | 0.7418 | 0.5696 |
| "tensor-thirds" analogue, P_z = 0, alignment T | (1+T)/4, (1−T)/4, (1−T)/4, (1+T)/4 | 0 | T | **0** |

The u = 3 row is the rational anchor already pinned by
`tests/test_spin.cpp:220-229` (populations (27, 9, 3, 1)/40; P_z = 0.7,
T = 0.4, R₃ = 0.2 exactly) and quoted in
`include/lipolgen/bookkeeping.hpp:201-202`.

Three consequences for the run plans:

* **No ⁷Li fill the programme builds is rank-3 clean, except the aligned one.**
  Every spin-temperature row has R₃ ≠ 0, and R₃/P_z grows from 0.04 at
  P_z = 0.3 to 0.63 at P_z = 0.9; the CLIs' `--plan pure` sits at R₃ = +1, the
  stretched maximum. §6 quantifies what that costs A_∥.
* **A pure-alignment fill is rank-3 clean.** The last row (P_z = 0, symmetric
  populations) has R₃ = 0 identically, because o(m) is odd in m. So the ⁷Li
  analogue of the A_zz thirds plan and of the transverse cos 2φ plan are free of
  the rank-3 sector by construction, for any T — which is the run-plan
  statement of the theorem in §4.
* **The explicit-P_zz branch silently sets R₃ = 0.** `helicity_flip_plan` with
  `use_explicit_pzz = true` calls `spin32_populations(pz, opt.pzz)`
  (`src/core/bookkeeping.cpp:89`) and the third argument defaults to `o = 0.0`
  (`include/lipolgen/spin.hpp:118-119`). That is a *choice of fill*, not a
  property of the physics. No C++ example CLI parses a `--pzz` flag, so the
  branch is unreachable from *them*; the Python API does expose it
  (`HelicityFlipOptions.use_explicit_pzz` and `.pzz`,
  `python/bindings.cpp:1797-1798`), so a Python caller who sets `pzz` gets
  R₃ = 0 silently today. If the rank-3 sector is ever switched on,
  `HelicityFlipOptions` needs an explicit `o` alongside `pzz` — and `RunPlan`
  needs somewhere to record it (§6.1).

The in-situ ⁷Li polarimeter the pipeline uses is ⟨P₂(cos θ_k)⟩ = −T/5 for any
fill (`tests/test_tagged.cpp:337-353`, `tests/test_pipeline.cpp:486-512`),
i.e. it measures the rank-2 moment. **There is no rank-3 polarimeter in the
generator**; R₃ would have to come from the fill model (14), not from a
measurement.

---

## 3. The Jaffe–Manohar basis and the four leading-twist functions

### 3.1 The arbitrary-spin basis

Jaffe and Manohar [5] treat lepton scattering from a target of arbitrary spin J
and show that in the Bjorken limit the target is described by **2J+1 quark
distributions per flavour**, evolving logarithmically — one per helicity state of
the target. That is the basis; everything below is a change of basis on it, from
helicity states to multipoles. [5] is not on arXiv and was not read for this
note; it is cited through [6] (which attributes its Eqs. (19) to [5]) and [1].

### 3.2 The convention-free object

Let q_↑^λ(x) be the probability of finding a quark of momentum fraction x and
*positive* helicity in a target of helicity λ ([6], text under Eq. (19)), with
parity giving q_↑^λ = q_↓^{−λ}. Define, for each target projection m, the
per-state F1 and the per-state helicity difference:

        F1^{(m)}(x) ≡ (1/2) Σ_q e_q² [ q_↑^m + q_↓^m ] + {q → q̄}              (15)
        G1^{(m)}(x) ≡ (1/2) Σ_q e_q² [ q_↑^m − q_↓^m ] + {q → q̄}              (16)

F1^{(m)} = F1^{(−m)} and G1^{(m)} = −G1^{(−m)}. These are normalization-free:
F1^{(m)} is what an unpolarized beam measures on a pure |J, m⟩ target, G1^{(m)}
what the beam-helicity difference measures. For J = 3/2 there are exactly two
independent F1^{(m)} (m = 3/2, 1/2) and two independent G1^{(m)} — four functions,
which is [6]'s count: *"in the parton model for the spin-3/2 sector, there are
four independent structure functions in deep inelastic scattering at leading
twist and leading order in α_s"* ([6] §3.1).

### 3.3 The four functions of [6] Eq. (19)

[6] Eqs. (19a)–(19d), transcribed literally from the source (the ↑ superscripts
are quark helicity, the numerical superscripts target helicity):

        F1 = (1/2) Σ_q e_q² { [q_↑^{3/2} + q_↑^{−3/2}] + [q_↑^{1/2} + q_↑^{−1/2}] } / 2
                                                              + {q → q̄}      (17a)
        b1 = (1/2) Σ_q e_q² { [q_↑^{3/2} + q_↑^{−3/2}] − [q_↑^{1/2} + q_↑^{−1/2}] } / 2
                                                              + {q → q̄}      (17b)
        g1 = (1/2) Σ_q e_q² { 3[q_↑^{3/2} − q_↑^{−3/2}] + [q_↑^{1/2} − q_↑^{−1/2}] } / √20
                                                              + {q → q̄}      (17c)
        g2 = (1/2) Σ_q e_q² { [q_↑^{3/2} − q_↑^{−3/2}] − 3[q_↑^{1/2} − q_↑^{−1/2}] } / √20
                                                              + {q → q̄}      (17d)

with the sum rules ([6] Eq. (21))

        ∫₀¹ dx b1(x) = 0 ,      ∫₀¹ dx g2(x) = 0                              (18)

*"if the quark sea q − q̄ does not contribute to the integral"*.

**Naming.** [6]'s fourth function is called `g2` there and is described as *"the
new structure function as the spin goes from 1 up to 3/2"*. It is **not** the
twist-3 nucleon g2. In the language of (15)–(16) its weight vector over
(m = 3/2, 1/2) is (1, −3) ∝ (o(3/2), o(1/2)) = (0.3, −0.9): it is the **rank-3
octupole partner of g1**, the l = 3 member of the multipole basis, and it exists
for J ≥ 3/2 only. The nucleon g2 exists already for J = 1/2 and is the
twist-3 partner *of the same l = 1 multipole* as g1. Two different objects, and
they can appear side by side in the same J = 3/2 cross section (§5). **This note
names [6]'s fourth function `g1_rank3`** and recommends that name in the code
comments; §5 introduces its own twist-3 partner as `g2_rank3`.

### 3.4 The multipole weights, and the exact map to the code

Using q_↑^λ + q_↑^{−λ} = q_↑^λ + q_↓^λ, Eqs. (17a)–(17d) read, in terms of the
convention-free objects (15)–(16):

        F1 = [ F1^{(3/2)} + F1^{(1/2)} ] / 2 ,   b1^{[6]} = [ F1^{(3/2)} − F1^{(1/2)} ] / 2
        g1^{[6]} = [ 3 G1^{(3/2)} + G1^{(1/2)} ] / √20 ,
        g1_rank3^{[6]} = [ G1^{(3/2)} − 3 G1^{(1/2)} ] / √20                  (19)

Inverting the first pair:

        F1^{(m)} = F1 + Q_NN(m) · b1^{[6]} ,    Q_NN(m) = [3m² − J(J+1)]/3     (20)

**The alignment weight in [6] is exactly the code's Q_NN of Eq. (8), with the
same normalization.** This is the central mapping result of §3: the code's choice
to extend `Q_NN = (3m² − J(J+1))/3` from J = 1 to J = 3/2 with T = ±1 on the
stretched state (`src/core/xsec.cpp:168-173`) is precisely the normalization in
which [6]'s b1 is the coefficient — no rescaling is needed.

The same statement holds independently in [4]. [4] Eq. (23) gives the parton
densities in each target helicity state,

        M_{λ+,λ+}(λ = +3/2) = f1 + f1LL + (3/2) g1 + (3/10) g1LLL
        M_{λ+,λ+}(λ = +1/2) = f1 − f1LL + (1/2) g1 − (9/10) g1LLL             (21)

The λ < 0 rows are **not** these two repeated. [4] Eq. (23) rows 3-4 read

        M_{λ+,λ+}(λ = −3/2) = f1 + f1LL − (3/2) g1 − (3/10) g1LLL
        M_{λ+,λ+}(λ = −1/2) = f1 − f1LL − (1/2) g1 + (9/10) g1LLL

— only the *even-l* terms (f1, f1LL) are even in λ; the odd-l terms flip, which
is parity, and which is exactly what the m and o(m) weights of (22) encode. So
the weights of f1LL are (+1, −1) = Q_NN
and the weights of g1LLL are (3/10, −9/10) = o(m) of Eq. (8) — **the code's
unnormalized rank-3 operator J_z³ − (41/20)J_z**. Hence

        F1^{(m)} = Σ_q e_q² [ f1 + Q_NN(m) f1LL ] ,
        G1^{(m)} = Σ_q e_q² [ m g1 + o(m) g1LLL ]                             (22)

and, projecting (22) onto (19), the two papers are consistent:

        3 G1^{(3/2)} + G1^{(1/2)}
              = Σe² [ (9/2 + 1/2) g1 + (9/10 − 9/10) g1LLL ] = 5 Σe² g1
        G1^{(3/2)} − 3 G1^{(1/2)}
              = Σe² [ (3/2 − 3/2) g1 + (3/10 + 27/10) g1LLL ] = 3 Σe² g1LLL

(check: 3·(3/2) + 1/2 = 5 and 3·0.3 + 2.7 = 3). Dividing by the √20 of (19),

        g1^{[6]}       = (5/√20) Σe² g1^{[4]}     = (√5/2)   Σe² g1^{[4]}
        g1_rank3^{[6]} = (3/√20) Σe² g1LLL^{[4]}  = (3/2√5)  Σe² g1LLL^{[4]}

Only the *weights* are convention-free; the overall constants differ between [6]
and [4] and neither is the code's.

**The code's own normalization**, which §3.5 and §6.4 need, is fixed by the rate
terms themselves. The vector term the generator already computes is
`w = λ_e P_e (m/J) cos θ_S · A_∥` (`src/core/xsec.cpp:267-269`), and the rank-3
slot this note proposes is `w₃ = λ_e P_e R₃(m) P₃(cos θ_S) D g1_rank3/F1`
(Eq. (45), used again in (49)). So the per-state helicity difference *in the
code's slots* is (m/J) g1 + R₃(m) g1_rank3, and matching it term by term to
G1^{(m)} of (22) gives

        g1|code       = J · Σe² g1^{[4]}     = (3/2)  Σe² g1^{[4]}
        g1_rank3|code = 0.3 · Σe² g1LLL^{[4]} = (3/10) Σe² g1LLL^{[4]}        (22′)

i.e. the stretched state's helicity difference is G1^{(3/2)} = g1 + g1_rank3 and
the |m| = 1/2 state's is G1^{(1/2)} = g1/3 − 3 g1_rank3, both read straight off
the code's tables. Equation (22′) is the map that makes §3.5's positivity bound
and §6.4's test 7 evaluable on `SFTables` without rescaling.

**The one sign that must be fixed by hand.** The generator's rank-2 rate term is
`TENSOR_LL_SIGN · Q_NN · P₂(cos θ_S) · K / D_φ` (`src/core/xsec.cpp:257-260`)
with `TENSOR_LL_SIGN = −1` (`include/lipolgen/constants.hpp:40`), i.e.

        F1^{(m)}|_code = F1 − Q_NN(m) · b1|_code                              (23)

For J = 1 this is exactly the Hoodbhoy–Jaffe–Manohar convention: [7] Eq. (5)
defines δ_T f = f⁰ − (f^{+1} + f^{−1})/2 and b1 = (1/2) Σ e_i² (δ_T q_i + δ_T q̄_i),
so b1 = F1^{(0)} − F1^{(1)}, which with (23) and Q_NN(1) − Q_NN(0) = 1 is an
identity. It is also the sign of every published number: A_zz = −(2/3) b1/F1
([2] Eq. (27), HERMES [8]), which the code reproduces exactly at every y and is
pinned at `tests/test_xsec.cpp:223-240` and `tests/test_xsec.cpp:242-255`, and
which `docs/CONVENTIONS.md:18-25` states as "b₁ > 0 means the m = 0 state has the
LARGER cross section". Comparing (20) with (23),

        b1_32|_code  =  − b1^{[6]}  =  − Σ_q e_q² [ f1LL^q + f1LL^q̄ ]         (24)

**Anyone filling `b1_32_func` from [6] or [4] must flip the sign.** The clash is
real, not a transcription slip on our side: for spin 1 the community's b1 is
positive when the *planar* (m = 0) state has more quarks, while [6]'s spin-3/2 b1
and [4]'s f1LL are positive when the *stretched* state does. One reading note,
which is **not** a claim that any source contradicts itself: [2] Eq. (3),
"b1 ∝ [q^{+1}(x) − q⁰(x)]", states only a *proportionality* and fixes no sign.
[2] itself writes f1LL = −(2/3) δ_T q with δ_T q = q⁰ − (q^{+1} + q^{−1})/2, so
its f1LL is positive with a positive constant when the stretched state is
enhanced, while δ_T q — and hence b1 — is positive when the planar state is;
the two statements are consistent. The sign the generator uses is fixed by [2]
Eqs. (2) + (27) and independently by [7] Eq. (5), which both give
b1 = F1^{(0)} − F1^{(1)}; that pair is what the generator and its tests are
anchored on.

### 3.5 Positivity constraints for the toy inputs

[4] Eqs. (32a)–(32b) impose A_{a+,a+} ≥ 0 and A_{a+,a+}A_{b−,b−} ≥ A²_{b−,a+};
the consequence for the diagonal (helicity-conserving) PDFs is [4] Eqs. (33a),
(33b), transcribed literally:

        f1 + f1LL ≥ | (3/2) g1 + (3/10) g1LLL |                               (25a)
        f1 − f1LL ≥ | (1/2) g1 − (9/10) g1LLL |                               (25b)

In the convention-free objects of §3.2 this is simply

        F1^{(m)}(x) ≥ | G1^{(m)}(x) |    for every m and every x              (26)

— "you cannot have more polarized quarks than quarks" applied state by state.
For the code's toy inputs this is a **cheap, complete positivity test on any
(b1_32, g1, g1_rank3) triple**, and it is strictly stronger than the two checks
the generator has today (`InclusiveKernel::positivity_margin`,
`include/lipolgen/xsec.hpp:324-326`, only tests W(φ) ≥ 0 in the azimuth). In the
code's own slots — **all three** structure functions in the code normalization,
`b1_32` with the code's sign of (23) and (g1, g1_rank3) carrying the weights the
rate terms actually apply — it reads

        | (m/J)·g1 + R₃(m)·g1_rank3 |  ≤  F1 − Q_NN(m)·b1_32                   (27)

              (m/J, R₃, Q_NN)  =  ( 1,   +1, +1 )   at m = 3/2
                                  ( 1/3, −3, −1 )   at m = 1/2

Equation (27) is *identical* to [4] Eqs. (33a,b) under the normalization map
(22′) — g1|code = (3/2) Σe² g1^{[4]}, g1_rank3|code = (3/10) Σe² g1LLL^{[4]},
F1|code = Σe² f1, b1_32|code = −Σe² f1LL. Writing it instead with [4]'s
per-flavour weights (m, o(m)) = (3/2, 3/10), (1/2, −9/10) while feeding it the
code's tables would overweight g1 by 3/2 and g1_rank3 by 3/10 at m = 3/2, i.e.
test against the wrong bound — that is the whole reason (27) is stated twice
over. Note that (27) at m = 1/2 is the tight one whenever
b1_32 > 0: the |m| = 1/2 doublet is the depleted state. [4] also gives the full
transversity bounds ([4] Eqs. (30), (31), (33c)); those constrain h1, h1LLT and
h1TT, none of which the generator carries, and are not repeated here.

---

## 4. Theorem: the rank-≤2 truncation is exact for an unpolarized beam

> **Theorem.** In inclusive lepton–nucleus scattering on a target of any spin J,
> at Born level in the one-photon-exchange approximation, with parity *and*
> time-reversal invariance (the latter, with current hermiticity, being what
> forbids T-odd inclusive structure functions — §4.3), the cross section
> measured with an **unpolarized lepton beam** depends on the target density matrix only through its **even-rank**
> multipoles (l = 0, 2, 4, …). The odd-rank multipoles (vector l = 1, octupole
> l = 3, …) enter only multiplied by the lepton helicity λ_e.
>
> **Corollary (the one LiPolGen needs).** For J = 3/2 the generator's rank-≤2
> truncation is *exact*, not an approximation, for every unpolarized-beam
> observable: the rate, A_zz, and the cos 2φ gluon-transversity amplitude. The
> rank-3 sector contaminates only **λ_e-odd** observables — A_∥, A_⊥ and the
> λ_e-odd azimuthal harmonics — and no unpolarized-beam observable at all.

### 4.1 Step 1 — the lepton tensor splits by λ_e

The inclusive lepton tensor for a lepton of helicity λ_e is ([1] Eq. (4.5b))

        L^{μν} = 4 l^μ l^ν − Q² g^{μν} − 2 l^{{μ} q^{ν}} + (2λ_e) 2i ε^{μνρσ} q_ρ l_σ
               ≡ L_S^{μν} + (2λ_e) i L_A^{μν}                                 (28)

with L_S real symmetric and λ_e-independent, L_A real antisymmetric. This is
exact, kinematic, and carries no hadronic input.

**Helicity convention, because the two differ by a factor 2.** In (28)–(30), and
in every equation quoted from [1] in this section, λ_e is **[1]'s λ_e = ±1/2**
([1] §IV A, "a pure spin state described by the helicity λ_e = ±1/2"), *not* the
note's own λ_e = ±1 of the conventions paragraph (the code's `lam_e`). Written
in the note's convention the antisymmetric term of (28) is
λ_e · 2i ε^{μνρσ} q_ρ l_σ and (30) reads L_{μν}W^{μν} = L_S·W_S − λ_e L_A·W_A.
Nothing in §4 depends on the factor; §5's master formula and every code-facing
formula are in the note's convention throughout.

### 4.2 Step 2 — the hadronic tensor splits by symmetry

Hermiticity of the electromagnetic current gives W^{μν*} = W^{νμ}, so

        W^{μν} = W_S^{μν} + i W_A^{μν} ,                                       (29)

W_S real symmetric, W_A real antisymmetric — both real by **hermiticity of the
current alone** ([1] Eq. (4.12): ⟨W^{{μν}}⟩ real, ⟨W^{[μν]}⟩ imaginary, which
follows from the symmetric (antisymmetric) part of the lepton tensor being real
(imaginary) and the contraction having to produce a real number). Time-reversal
invariance is **not** used here; it enters only in Step 3, where it removes the
T-odd tensor structures. Contracting (28) with (29), the
symmetric×antisymmetric cross terms vanish index by index:

        L_{μν} W^{μν} = L_S · W_S − (2λ_e) L_A · W_A                          (30)

**An unpolarized beam sees W_S alone.** Everything else is the statement that
the surviving odd-rank structures live in W_A and the surviving even-rank ones
in W_S — which is §4.3, and which needs time reversal as well as parity.

### 4.3 Step 3 — parity and time reversal assign the multipoles

By [1] App. D, the l-th multipole is carried by a symmetric traceless rank-l
Lorentz tensor that is a *regular* tensor for even l and a *pseudo*-tensor for
odd l (the l = 1 case is the familiar axial spin vector S^μ; the l = 3 case is
R^{μνρ}). W^{μν} must be a regular tensor. Parity therefore fixes the
**ε-content** of the allowed structures — and nothing more:

* **Odd l.** The multipole tensor is a pseudo-tensor, so one ε is *required*.
  The natural structures ε^{μνρσ} q_ρ Y_σ, with Y the pseudovector obtained by
  contracting the multipole with q (Y = S for l = 1; Y^σ = R^{σαβ} q_α q_β for
  l = 3 — the P-contractions vanish by (2)), are antisymmetric and T-even: they
  live in W_A. But parity *also* allows the regular **symmetric** structure
  Z^{{μ} P̃^{ν}} with Z^μ = ε^{μρσλ} q_ρ P_σ Y_λ, which is current-conserving
  (q·Z = 0) — the T-odd single-spin structure, sitting in W_S.
* **Even l.** The multipole tensor is regular, so structures such as t^{μν},
  g^{μν} t^{αβ}q_α q_β, P̃^μ P̃^ν t^{αβ}q_α q_β, t^{{μ}_α q^α P̃^{ν}} can be built
  with no Levi-Civita symbol; all of these are symmetric and sit in W_S. **It is
  not true that an antisymmetric even-l structure would need an ε**:
  t^{[μ}{}_α q^α P̃^{ν]} — and, for l = 0, P̃^{[μ} q^{ν]} — is regular,
  antisymmetric and ε-free. Contracted with L_A it produces the
  t_{yz} ∝ sin φ_S structure, i.e. exactly the even-multipole/polarized-lepton
  sin terms of [1] rule (i)(b) and the T-odd λ_e-dependent tensor structures
  listed in [1] Eq. (4.40): F^{sin φ_h}_{L T_LL}, F^{sin φ_TL}_{L T_LT},
  F^{sin(φ_h−φ_TL)}_{L T_LT}, F^{sin(2φ_h−2φ_TT)}_{L T_TT}, and the rest. **So
  parity alone does not put even l in W_S.**

**Parity alone therefore closes neither direction, and the same statement closes
both.** Time-reversal invariance, combined with hermiticity of the
electromagnetic current, forbids T-odd structures in the inclusive cross section
at one-photon exchange: [1] Eq. (4.44) with its Ref. [62] (N. Christ and
T. D. Lee, Phys. Rev. **143**, 1310 (1966)), stated there as
∫dΓ_{P_h} F^{sin φ_S}_{U S_T} = 0 (odd multipole, unpolarized lepton) and
∫dΓ_{P_h} F^{sin φ_TL}_{L T_LT} = 0 (even multipole, polarized lepton) — one
example of each parity. That removes the odd-l Z^{{μ}P̃^{ν}} structure *and* the
even-l ε-free antisymmetric structure, and the theorem survives. This is the one
place the argument uses more than parity, it is used for **both** parities, and
it is the assumption that fails at two-photon exchange ([1] Ref. [63]:
A. Afanasev, M. Strikman, C. Weiss, Phys. Rev. D **77**, 014028 (2008)).

**The clean version, for m_l = 0, needs only parity.** For the azimuth-free
component the whole argument can be run in the photon-helicity basis of §4.4,
with no appeal to time reversal at all. Write the m_l = 0, rank-l part of a
diagonal density matrix as ρ_{MM}|_l = w_l(M), with
w_l(M) ∝ ⟨J, M; l, 0 | J, M⟩ a polynomial of degree l in M obeying
w_l(−M) = (−1)^l w_l(M) — l = 1: M; l = 2: 3M² − J(J+1); l = 3: M³ − (41/20)M
for J = 3/2, which is Eq. (8) read as a polynomial. Parity relates the
photon-helicity amplitudes as W_{−λ',−λ}(−M, −M) = W_{λ'λ}(M, M), so

        W_{−λ',−λ}[ρ_l]  =  (−1)^l  W_{λ'λ}[ρ_l]                             (30′)

Hence, channel by channel: **T** (W_{++} + W_{−−}) and **L** (W_{00}) are even
under λ → −λ and so vanish for odd l; **T′** (W_{++} − W_{−−}) is odd and so
vanishes for even l. T and L are the unpolarized-beam channels and T′ the
λ_e-odd one, so m_l = 0 odd-rank multipoles are λ_e-odd and m_l = 0 even-rank
multipoles are λ_e-even **by parity alone**. For |m_l| = 1 parity fixes only cos
versus sin (rule (i)), and it is time-reversal invariance that kills the sin
terms, as above.

### 4.4 Step 4 — the same result from photon helicities, with the counting

The rest-frame version is sharper and gives the structure-function count for
free. Project W^{μν} onto virtual-photon helicity states λ, λ' ∈ {+1, 0, −1}.
Angular momentum about the photon axis forces the target density-matrix element
ρ_{MM'} that a given (λ', λ) pair can reach to satisfy M − M' = λ' − λ, i.e. a
multipole component with azimuthal index

        m_l = λ' − λ ,   so   |m_l| ≤ 2 in inclusive DIS.                     (31)

Now decompose the lepton tensor in the same basis, ε*_{λ'μ} L^{μν} ε_{λν}:

* λ' = λ = ±1 : ε*_± ⊗ ε_± has both a symmetric and an antisymmetric part
  → **two** unpolarized channels are not available here, but one is: the
  T channel (λ'=λ=±1, summed) for L_S, and the T′ channel (the ± difference) for
  L_A.
* λ' = λ = 0 : ε_0 ⊗ ε_0 is symmetric → **L channel only** (unpolarized).
* |λ' − λ| = 1 : mixed → **LT** (unpolarized) *and* **LT′** (λ_e-odd).
* |λ' − λ| = 2 : ε*_{+} ⊗ ε_{−} = −(1/2)(ê_x − iê_y) ⊗ (ê_x − iê_y) is the outer
  product of a vector with **itself**, hence symmetric → **TT only**, no λ_e-odd
  partner.

Combining with (31), the inclusive structure-function count per multipole is:

        l even, unpolarized beam :  m_l = 0  → 2 functions (T and L)
                                    |m_l| = 1 → 1 function  (LT)
                                    |m_l| = 2 → 1 function  (TT)              (32)
        l odd,  λ_e-odd          :  m_l = 0  → 1 function  (T′)
                                    |m_l| = 1 → 1 function  (LT′)
                                    |m_l| ≥ 2 → 0

and every |m_l| ≥ 3 component decouples from inclusive DIS entirely. The
counting reproduces the known cases exactly:

        J = 1/2 :  l=0 → 2 (F1, F2)   l=1 → 2 (g1, g2)                  = 4
        J = 1   :  + l=2 → 4 (b1, b2, b3, b4)                           = 8   (33)
        J = 3/2 :  + l=3 → 2 (g1_rank3, g2_rank3)                       = 10
        J = 2   :  + l=4 → 4                                            = 14

The J = 1 line is exactly [2]'s four tensor structure functions F[U T_LL,T],
F[U T_LL,L], F[U T_LT], F[U T_TT] of its Eq. (10) mapped onto b1…b4 by its
Eqs. (17a)–(17e), plus F1, F2, g1, g2 — the Hoodbhoy–Jaffe–Manohar total. That
the same rule reproduces both the spin-1/2 and the spin-1 answers is the check on
(32).

### 4.5 Step 5 — the same result quoted from the source

[1] App. D rule (i), verbatim:

> (i) The azimuthal dependence alternates between sin and cos and reflects the
> parity of the terms. This follows the pattern
> (a) cos-dependence for even parity: l = even multipole combined with
> unpolarized lepton or l = odd multipole combined with polarized lepton.
> (b) sin-dependence for odd parity: l = odd multipole combined with unpolarized
> lepton or l = even multipole combined with polarized lepton.

Applied to the inclusive case, where the only azimuth is the target-spin azimuth
φ_S relative to the lepton plane: an odd-l multipole with an unpolarized lepton
can appear only as sin(|m_l| φ_S), which is T-odd and vanishes at Born level, and
for the m_l = 0 component sin(0) = 0 identically. Symmetrically, an even-l
multipole with a polarized lepton can appear only as a sin term and drops out for
the same reason. [2] states the spin-1 case of this outright, twice: *"the
inclusive cross section with an unpolarized electron beam is only sensitive to
tensor polarization Q, not vector polarization P"*, and footnote 2, *"In the case
of polarized electron beams, two additional structure functions appear (g1 and
g2), which couple to deuteron vector polarization only."* [1] Eq. (4.38) shows
the pattern in its full SIDIS form: S_T with an unpolarized lepton produces only
**sin**-type terms — five of them, sin(φ_h − φ_S) (carrying two structure
functions, T and L), sin(φ_h + φ_S), sin(3φ_h − φ_S), sin φ_S,
sin(2φ_h − φ_S) — and with (2λ_e) only **cos**-type terms, three of them,
cos(φ_h − φ_S), cos φ_S, cos(2φ_h − φ_S). After ∫dφ_h only sin φ_S and cos φ_S
survive, and the first of those is then set to zero by time-reversal invariance:
[1] Eqs. (4.41d) and (4.44a). That pair is the inclusive statement.

### 4.6 What is assumed, precisely

1. **One-photon exchange, Born level.** Two-photon exchange generates absorptive
   phases and hence T-odd inclusive structure functions; step 4.3 fails. The
   explicit statement, and the calculation of the resulting transverse
   single-spin asymmetry, is [1] Ref. [63]: **A. Afanasev, M. Strikman and
   C. Weiss**, *Transverse target spin asymmetry in inclusive DIS with
   two-photon exchange*, Phys. Rev. D **77**, 014028 (2008), arXiv:0709.0901.
   The effect is O(α) and is not in the generator at all.
2. **No T-odd inclusive structure functions at Born level.** This is not an
   extra assumption: it is a theorem, and it is the one §4.3 leans on. Time
   reversal combined with hermiticity of the electromagnetic current forbids
   T-odd structures in the inclusive cross section in the one-photon-exchange
   approximation — [1] Eq. (4.44) with its Ref. [62], **N. Christ and
   T. D. Lee**, Phys. Rev. **143**, 1310 (1966). [1] writes it as the vanishing
   of the φ_h-integrals ∫dΓ F^{sin φ_S}_{U S_T} = 0 and
   ∫dΓ F^{sin φ_TL}_{L T_LT} = 0. The *tagged* analogue, where FSI can supply
   the missing phases, is explicitly open — [10] §VII, item (i)(c): *"What is
   the size of T-odd structure functions, which are zero in the IA and provide a
   sensitive test of FSI?"*
3. **Pure electromagnetic current.** γ–Z interference gives a parity-violating
   term proportional to the target *vector* polarization with an unpolarized
   beam, breaking the theorem. Its scale is Q²/M_Z² ≈ 1.2 × 10⁻³ at
   Q² = 10 GeV², i.e. a PV asymmetry A_PV ≈ G_F Q²/(4√2 π α) ≈ 10⁻⁴ × Q²/GeV²
   ≈ 10⁻³ there — far below anything LiPolGen quotes, but a real exception, and
   it should be stated when the theorem is used in a paper.
4. **Nothing about the target's internal dynamics.** No convolution model, no
   nucleon-only assumption, no leading twist. The theorem is kinematic +
   discrete symmetries and therefore holds at every x, Q² and γ.

### 4.7 The two consequences for the generator

* **Exactness.** Every unpolarized-beam number the programme has published or
  will publish — the rate, A_zz, the cos 2φ amplitude — is *complete* at rank ≤ 2
  for J = 3/2. `docs/open_items/physics_literature.md` §1's "Default: rank ≤ 2
  truncation" can be restated as a theorem with a citation. The generator's own
  numerical statement of this is already in place at
  `tests/test_xsec.cpp:125-130`: an unpolarized electron on a J = 3/2 target with
  zero rank-2 slots gives `w_avg == a1 == a2 == 0.0` exactly.
* **Non-exactness where it matters.** *Every* λ_e-odd observable — A_∥, A_⊥ and
  the λ_e-odd azimuthal harmonics — receives the l = 1 **and** the l = 3
  multipole on a J = 3/2 target. At γ → 0 only the A_∥-like rate term does
  (Eq. (45)); at finite γ the LT′ channel of (42) feeds a₁ at O(γ), and the full
  cubic geometry of (43) feeds a₂ and a₃ as well (Eq. (47)). Since every ⁷Li
  vector fill the generator builds has R₃ ≠ 0 (§2.4), the rank-3 sector is a
  genuine, if small, systematic on any ⁷Li g1 extraction. §6 quantifies it.

---

## 5. Finite-γ inclusive decomposition for J = 3/2

### 5.1 Kinematics, in the code's variables

        γ² = 4 M² x² / Q²                    (`gamma_squared`, asymmetries.hpp:60)
        ε  = (1 − y − γ²y²/4) / (1 − y + y²/2 + γ²y²/4)     ([2] Eq. (11c))    (34)
        cos θ_q = (1 + γ²y/2)/√(1+γ²) ,
        sin θ_q = γ√(1 − y − γ²y²/4)/√(1+γ²)               ([2] Eqs. (24a,b))

θ_q is the target-rest-frame angle between the virtual photon and the incoming
electron; it vanishes as γ → 0. Equations (34) are `epsilon_gamma`
(`include/lipolgen/asymmetries.hpp:63`) and `theta_q_cos_sin`
(`src/core/xsec.cpp:12-18`). The spin axis at (θ_S, φ_S) to the **beam** has,
in the **photon** frame with x̂ along the incoming lepton's transverse projection
(`include/lipolgen/xsec.hpp:263-272`),

        N = ( c s_S cos φ' + s c_S ,  s_S sin φ' ,  c c_S − s s_S cos φ' )     (35)

with c, s = cos θ_q, sin θ_q and c_S, s_S = cos θ_S, sin θ_S, φ' = φ − φ_S.
Equation (35) is verbatim `include/lipolgen/xsec.hpp:266-267`; the frame choice
is the one that makes both finite-γ rows of [2] Table 1 come out, pinned at
`tests/test_tensor_gamma.cpp:382` and `:406`.

### 5.2 The rank-2 block: unchanged in form, and already in the code

By [1] App. D, passing from spin 1 to spin 3/2 *adds* the l = 3 multipole and
changes nothing about l ≤ 2. The alignment tensor T^{μν} of a spin-3/2 state is
the same object — symmetric, traceless, orthogonal to P — as the spin-1 t^{μν},
so it admits exactly the same independent basis tensors A_i^{μν}, the same four
inclusive structure functions and hence the same four b-functions. Concretely,
[2] Eq. (10) holds for J = 3/2 word for word:

        dσ/(dx dQ²) = π y² α² / [Q⁴(1 − ε)] × {
              F_[UU,T] + ε F_[UU,L]
            + T_LL ( F_[U T_LL,T] + ε F_[U T_LL,L] )
            + √(2ε(1+ε)) T_LT cos φ_TL F_[U T_LT]
            + ε T_TT cos 2φ_TT F_[U T_TT] }                                   (36)

with ([2] Eq. (16))

        F_[UU,T] = 2 F1 ,   F_[UU,L] = (1 + γ²) F2/x − 2 F1                    (37)

and ([2] Eqs. (17a)–(17e); the leading factor 2 in (17a) multiplies b1 alone).
**Variable note.** In (36)–(38) the note's `x` is **[2]'s x_d**, *not* [2]'s
x = 2x_d ([2] Eq. (12)), and F2, b2 are per-nucleon-normalized so that the
Callan–Gross relations read F2 = 2x F1 and b2 = 2x b1. [2] itself prints them as
F2 = 2x_d F1 = x F1 and b2 = 2x_d b1 = x b1 ([2] Eq. (18)) in its own x = 2x_d.
The generator uses the per-nucleon form throughout — x per nucleon
(`include/lipolgen/asymmetries.hpp:57-61`) and `t.b2 = 2·x·t.b1` as the default
(`src/core/xsec.cpp:116`) — so the note follows the code, and anyone comparing
against [2] must remember that the note's x is [2]'s x_d.

        F_[U T_LL,T] = −[ 2(1+γ²) b1 − (γ²/x)( b2/6 − b3/2 ) ]                (38a)
        F_[U T_LL,L] = (1/x)[ 2(1+γ²) x b1 − (1+γ²)²( b2/3 + b3 + b4 )
                              − (1+γ²)( b2/3 − b4 ) − ( b2/3 − b3 ) ]         (38b)
        F_[U T_LT]   = −(γ/2x)[ (1+γ²)( b2/3 − b4 ) + ( 2b2/3 − 2b3 ) ]       (38c)
        F_[U T_TT]   = −(γ²/x)( b2/6 − b3/2 )                                 (38d)

The γ → 0 collapse has to be stated carefully, because [2] Eqs. (25)–(26) take
the **full Bjorken limit** — γ → 0 *and* both Callan–Gross relations — whereas
the code takes γ → 0 at whatever b2 and R it was handed. At γ → 0 alone:

        F_[U T_LL,T] → −2 b1
        F_[U T_LL,L] → (2x b1 − b2)/x       (zero only at b2 = 2x b1)
        F_[U T_LT] = F_[U T_TT] → 0         (b3, b4 cancelling, for any values)
        F_[UU,L]    → F2/x − 2F1 = 2 F1 R   (zero only at R = 0)

which is exactly what the code's own header and its pinned test state
(`include/lipolgen/xsec.hpp:99-103`, `src/core/xsec.cpp:29-32`, and
`tests/test_tensor_gamma.cpp:152`, which pins the collapse explicitly "for any
b2"). The F_[UU,L] line is not a technicality: it is the same R that puts the
(1 + εR) into the pinned **spin-1** identity A_zz(1 + εR) = −(2/3) b1/F1
(`tests/test_xsec.cpp:223-240`; the J = 3/2 contrast has no 2/3, §5.4), and the
generator's default `r_sigma_lt` (`src/core/xsec.cpp:41`) is not zero. Imposing
the two Callan–Gross relations on top — b2 = 2x b1 ([2] Eq. (18)) and R = 0 —
recovers [2] Eqs. (25), (26) as printed.

**What differs for J = 3/2, and it is only two things.** (i) T_LL, T_LT, T_TT are
built from the *spin-3/2* alignment via Eqs. (10) and (12) — numerically this is
the single substitution Q_NN(m) = (3m² − 15/4)/3 = ±1 in place of the J = 1
values (1/3, −2/3, 1/3), which is exactly what
`InclusiveKernel::tensor_moments` already does off `ion().spin`
(`src/core/xsec.cpp:168-173`). (ii) b1…b4 are spin-3/2 matrix elements with
their own values and their own sign convention (§3.4, Eq. (24)).

**Therefore the generator's `b1_32_func / b2_32_func / delta_32_func` slots and
the whole `tensor_gamma` machinery apply to ⁷Li verbatim.** Verified in code:
`InclusiveKernel::tables` dispatches to the `_32` slots on `ion().spin == 1.5`
(`src/core/xsec.cpp:108-118`) and fills b3, b4 for *every* spin
(`src/core/xsec.cpp:119-122`); `tensor_harmonics_gamma` then calls
`tensor_moments(state.m)`, `cosyn_tensor_sfs` and `cosyn_unpolarized_sfs` with
no spin-dependent branch at all (`src/core/xsec.cpp:175-212`). Nothing in the
finite-γ path is spin-1-specific. The only thing the J = 3/2 case lacks is a
physics input for b1…b4 — which is `OPEN_ITEMS_SOLUTIONS.md` row 10 (b1 for
A > 2), not row 8.

Explicitly, in the code's variables, the rank-2 harmonics for a pure |J, m⟩ fill
are (this is `tensor_harmonics_gamma`, restated):

        w₂(φ') = −TENSOR_LL_SIGN ×
                 [ t_zz(φ') Λ_LL + t_xz(φ') Λ_LT + (t_xx − t_yy)(φ') Λ_TT ]
                 / ( F_[UU,T] + ε F_[UU,L] )                                  (39)
        Λ_LL = F_[U T_LL,T] + ε F_[U T_LL,L] ,
        Λ_LT = √(2ε(1+ε)) F_[U T_LT] ,   Λ_TT = ε F_[U T_TT]

with t_ij from (10) evaluated on N of (35) and T = Q_NN(m), each t_ij being a
quadratic in cos φ' and hence exactly a constant plus cos φ' plus cos 2φ'.

### 5.3 The rank-3 block

**Structures and count.** The l = 3 multipole enters W_A (§4.3). The only
pseudovector linear in R^{μνρ} is Y^σ = R^{σαβ} q_α q_β, since P-contractions
vanish by (2); with q it spans the two current-conserving antisymmetric
structures ε^{μνρσ} q_ρ Y_σ and ε^{μνρσ} q_ρ P̃_σ (q·Y), exactly parallel to the
pair that carries g1 and g2 for the spin vector. The photon-helicity count (32)
gives the same answer and is the one to quote: **two** structure functions,

        F_[L R_LLL, T′]   ( m_l = 0 ,  prefactor √(1 − ε²) )
        F_[L R_LLT, LT′]  ( |m_l| = 1, prefactor √(2ε(1 − ε)) )               (40)

and **zero** from R_LTT (|m_l| = 2: the λ_e-odd lepton tensor has no TT channel)
and from R_TTT (|m_l| = 3 > 2). Both prefactors are fixed, not free: [1]
App. D rule (iii) states that a projection's *"azimuthal and ε dependence follows
that from lower spin multipoles with the same m_l value"*, and [1]'s own spin-1
Eqs. (4.38), (4.39) exhibit √(1 − ε²) for the λ_e-odd m_l = 0 term and
√(2ε(1 − ε)) for the λ_e-odd |m_l| = 1 terms — **after integration over φ_h**.
The qualifier is needed: before that integration the S_L(2λ_e) block of (4.38)
also carries a √(2ε(1 − ε)) cos φ_h F^{cos φ_h}_{L S_L} term, which is m_l = 0
and nonetheless carries the |m_l| = 1 prefactor. What distinguishes it is its
φ_h dependence, not its multipole index, and it integrates away. That is the
whole content of "the kinematic prefactors" for the *inclusive* rank-3 sector:
**they are identical to those of g1 and g2.**

Name the two functions, in the code's convention,

        g1_rank3(x, Q²)  ≡ the leading-twist l = 3 function of (17d) / [4]'s g1LLL
        g2_rank3(x, Q²)  ≡ its twist-3 partner                                (41)

so that, by the identity of channels with the vector sector and with the
constants of [1]'s own **inclusive** vector reductions — [1] Eq. (4.41c),
∫dΓ_{P_h} F_{L S_L} = 2(g1 − γ² g2), and [1] Eq. (4.41d),
∫dΓ_{P_h} F^{cos φ_S}_{L S_T} = −2γ(g1 + g2), note the sign on the second —

        F_[L R_LLL, T′]  =  2 ( g1_rank3 − γ² g2_rank3 ) ,
        F_[L R_LLT, LT′] = −2γ ( g1_rank3 + g2_rank3 )                        (42)

Equation (42) is the **definition** of the pair (g1_rank3, g2_rank3): the
rank-3 functions are named so that they enter their channels with exactly the
constants the vector functions enter theirs, which is what makes the rank-3
block a literal copy of the vector block in the code. It is therefore also the
rank-3 statement of the generator's own A₁/A₂ decomposition, `a_parallel_exact`
(`include/lipolgen/asymmetries.hpp:69-81`): with A₁ = (g1 − γ²g2)/F1 and
A₂ = γ(g1 + g2)/F1 for the vector sector,

        A₁^{(3)} = ( g1_rank3 − γ² g2_rank3 ) / F1 ,
        A₂^{(3)} = γ ( g1_rank3 + g2_rank3 ) / F1

which §6.2 uses verbatim.

**Geometry.** With R_ijk of (11) evaluated on N of (35),

        R_LLL(φ')            = R₃ · P₃(N_z) = (R₃/2)( 5 N_z³ − 3 N_z )
        R_LLT cos φ_LLT (φ') = R_xzz = (R₃/2) N_x ( 5 N_z² − 1 )              (43)

Both are **cubic** in cos φ', so the rank-3 sector generates harmonics
cos 0φ′ … **cos 3φ′** — one more than the rank-2 sector. This is the single
structural surprise of the calculation and it is what forces a new amplitude slot
in the code (§6).

**γ → 0 limit.** At γ = 0, θ_q = 0, so N = (s_S cos φ', s_S sin φ', c_S) and

        R_LLL → R₃ P₃(cos θ_S)          (φ'-independent)
        R_LLT cos φ_LLT → (R₃/2) sin θ_S ( 5cos²θ_S − 1 ) cos φ'              (44)

while the LT′ prefactor of (40) carries an explicit γ through (42). So at γ = 0
**exactly one** rank-3 function survives, multiplying a φ-independent,
λ_e-odd, P₃(cos θ_S)-weighted rate shift:

        w₃(γ = 0) = λ_e P_e · R₃ · P₃(cos θ_S) · D(y) · g1_rank3 / F1          (45)

with D(y) the generator's own `depolarization_d`
(`include/lipolgen/asymmetries.hpp:39-40`). Compare the vector term the code
already computes, `w = λ_e P_e (m/J) cos θ_S · A_∥`
(`include/lipolgen/xsec.hpp:19`, `src/core/xsec.cpp:266-270`): the rank-3 term is
the same expression with (m/J) → R₃ and P₁(cos θ_S) → P₃(cos θ_S). The
leading-twist count checks out: F1, b1, g1, g1_rank3 — the four functions of [6]
§3.1 (Eq. (17) above).

### 5.4 Master formula

Collecting §5.2 and §5.3, the inclusive cross section for an electron of helicity
λ_e and polarization P_e on a J = 3/2 ion in an axially symmetric fill with
moments (P_z, T, R₃) along n̂(θ_S, φ_S), at finite γ, in the generator's own
variables and normalization
(`include/lipolgen/xsec.hpp:13-21`, `src/core/xsec.cpp:283-285`):

        dσ / (dx dQ² dφ)  =  [ σ_U(x, Q²) / 2π ] · W(φ')                       (46)
        W(φ') = 1 + w₀ + a₁ cos φ' + a₂ cos 2φ' + a₃ cos 3φ'

with, writing D_γ = F_[UU,T] + ε F_[UU,L] for the exact denominator ([2] Eq. (16)
and (37); `cosyn_unpolarized_sfs`, `src/core/xsec.cpp:39-42`),

        w₀, a₁, a₂  ⊃  the rank-2 harmonics of (39)      [ unpolarized beam ]
        a₂          ⊃  −(1−y)/y² · 3Q_NN · sin²θ_S · Δ / D_φ   [ gluon transversity ]
        w₀, a₁      ⊃  λ_e P_e × the rank-1 (vector) sector, i.e. A_∥ and A_⊥   (47)
        w₀ … a₃     ⊃  λ_e P_e × [ R_LLL(φ') √(1−ε²) F_[L R_LLL,T′]
                                 + R_LLT cos φ_LLT(φ') √(2ε(1−ε)) F_[L R_LLT,LT′] ] / D_γ

The rank-3 line of (47) is the *only* new term relative to what the generator
computes today. Its γ → 0 collapse is (45); its λ_e = 0 value is exactly zero;
its φ-average is R₃ times the constant harmonic of P₃(N_z), which at γ = 0 is
R₃ P₃(cos θ_S).

Two immediate uses of (46)–(47):

* **A_zz for J = 3/2 — the definition has to come before the number.** With an
  unpolarized beam, (47) reduces to the rank-2 lines alone, and for a fill at
  θ_S = 0 Eq. (23) gives σ^{(m)} = σ_U[1 − Q_NN(m) K/D_φ], with Q_NN(m) = ±1 and
  K = b1_32 + (1−y)/(xy²) b2_32. The natural ⁷Li observable is the **alignment
  contrast** between the stretched (T = +1) and the |m| = 1/2 (T = −1) fills,
  and it does *not* carry the spin-1 factor 2/3:

        A_T^{(3/2)}  ≡  [σ(T=+1) − σ(T=−1)] / [σ(T=+1) + σ(T=−1)]
        A_T^{(3/2)}(θ_S = 0) · (1 + εR)  =  − b1_32 / F1                       (48)

  The 2/3 of the spin-1 relation is not universal. It is the explicit `2.0/3.0`
  in `azz()` (`src/core/asymmetries.cpp:85`) and the 1/Q with Q = 3T_LL of [2]
  Eq. (19) — i.e. the spin-1 (1/3, −2/3) population weights of the thirds
  estimator, which have no J = 3/2 counterpart, because for J = 3/2
  Q_NN = ±1 exactly. If continuity with [2] Eq. (19) is wanted, *define*
  A_zz^{(3/2)} ≡ (2/(3T))(σ_T/σ_U − 1) and the −(2/3) b1_32/F1 comes back by
  construction; what must not happen is quoting −(2/3) b1/F1 for a contrast that
  actually returns −b1/F1. §6.4 test 10 pins whichever definition is adopted.
  Either way, the ⁷Li analogue of `tensor_thirds_plan` does not exist in the
  code yet (`src/core/bookkeeping.cpp:107-120` hard-codes spin 1).
* **A_∥ for J = 3/2.** From (45) and (47) the measured longitudinal asymmetry on
  a fill at θ_S = 0 is

        A_∥^meas / (P_e P_z D)  =  g1/F1  +  (R₃ / P_z) · g1_rank3/F1          (49)

  — an *additive* bias on every extracted g1/F1, with coefficient R₃/P_z read off
  §2.4: 0.041 at P_z = 0.3, 0.286 at P_z = 0.7, 0.633 at P_z = 0.9.

---

## 6. What changes in the code if the rank-3 sector is switched on

This section is a design sketch, not an instruction; nothing here is implemented.

### 6.1 New slots and fields

| where | add |
|---|---|
| `InclusiveKernel::Options` (`include/lipolgen/xsec.hpp:197-233`) | `SFFunc3 g1_rank3_func, g2_rank3_func;` and `bool rank3 = false;` — off by default, for the same reason `tensor_gamma` is (`include/lipolgen/xsec.hpp:204-217`): switching it on moves a published number with nothing on the analysis side to meet it |
| `SFTables` (`include/lipolgen/xsec.hpp:131-145`) | `double g1_rank3 = 0.0, g2_rank3 = 0.0;` filled in `tables()` only for `ion().spin == 1.5`, in the same branch as the `_32` slots (`src/core/xsec.cpp:108-118`) |
| `Amplitudes` (`include/lipolgen/xsec.hpp:155-160`) | `double a3 = 0.0;` — required by (43). `density()` and `positivity_margin()` (`src/core/xsec.cpp:311-318`) take it as one more term, but **`density_min` (`src/core/xsec.cpp:44-53`) changes algorithm, not just signature**: today it is an exact minimiser of a quadratic in c = cos φ′; with cos 3φ′ = 4c³ − 3c the stationary condition becomes a *cubic* in c, so it needs a cubic root solve (or a guarded grid scan over c ∈ [−1, 1] plus a Newton polish, with the endpoints kept) |
| `InclusiveSampler::StateTables` (`include/lipolgen/sampler.hpp:208-215`) | `a3` and `a3n` vectors; `bound` becomes `1 + \|a1n\| + \|a2n\| + \|a3n\|`; `margin` from the new `density_min`. All of it is filled at `src/core/sampler.cpp:234-240` |
| the per-event φ draw (`src/core/sampler.cpp:386-396`) | the accept–reject test is literally `u < 1 + a1n cos φ′ + a2n cos 2φ′` with `u` drawn against `bound`; both gain the cos 3φ′ term |
| `InclusiveSampler::effective_modulation` (`src/core/sampler.cpp:305-322`) | a third numerator accumulator and an `a3` field on `EffectiveModulation` |
| the per-category φ density (`src/core/sampler.cpp:480-482`) and `phi_histogram_pseudo` (`src/core/sampler.cpp:526-531`) | both evaluate `1 + w_avg + a1 cos φ′ + a2 cos 2φ′` by hand; the second also calls `density_min(a1, a2)` |
| `InclusiveKernel` | `double octupole_moments(double m) const` returning `(m³ − (41/20)m)/0.3` for spin 3/2 and `0.0` otherwise, the exact analogue of `tensor_moments` (`src/core/xsec.cpp:168-173`) and consistent with `src/core/spin.cpp:246-259`. **Plural, mirroring `tensor_moments`**: a member named `octupole_moment` would hide the free function `lipolgen::octupole_moment(const CplxMatrix&, double)` (`include/lipolgen/spin.hpp:102`) inside the class scope |
| `HelicityFlipOptions` (`include/lipolgen/bookkeeping.hpp:137-147`) | an explicit `o` beside `pzz`, since the explicit branch silently sets R₃ = 0 today (`src/core/bookkeeping.cpp:89`) |
| `RunPlan` (`include/lipolgen/bookkeeping.hpp:89-129`) | **there is nowhere to record R₃ today.** The class carries `pe/pz/pzz` true + measured only, and `helicity_flip_plan` records just `.tensor` (`src/core/bookkeeping.cpp:102-104`), so the `HelicityFlipOptions::o` above has no destination. Add `o_true_` / `measured_o_`, recorded from `moments_along_axis(j, pops).octupole`, and give the smear block (`src/core/bookkeeping.cpp:33-41`) a policy — noting §2.4 that **there is no rank-3 polarimeter**, so the honest default is *not* a fourth `rng.normal()` draw but R₃ taken from the fill model (14) with its own systematic |
| `python/bindings.cpp:1687-1689` | the two new `Options` members, beside the existing `b1_32_func` / `b2_32_func` / `delta_32_func` |
| `python/bindings.cpp:1647-1656` | the `Amplitudes` binding gains `a3` — and this is a **breaking** change, not an addition: its `__iter__` and `__repr__` are a fixed 3-tuple, so every Python caller doing `w, a1, a2 = amps` breaks the day `a3` appears |
| `python/bindings.cpp:1950-1951` | the `state_tables` dict export gains `a3`, `a3n` |

Naming: **use `g1_rank3`, never `g2`**, for [6]'s fourth function (§3.3); the
name `g2` is already taken in `SFTables` by the twist-3 nucleon g2
(`include/lipolgen/xsec.hpp:143`) and the collision would be silent and severe.

### 6.2 Where the λ_e term enters

`InclusiveKernel::amplitudes` (`src/core/xsec.cpp:214-276`). The rank-2 branch
sits at `src/core/xsec.cpp:248-264`; the λ_e-dependent vector block at
`src/core/xsec.cpp:266-274`. The rank-3 term
belongs **inside the λ_e block**, gated on `helicity != 0.0` exactly as the
vector term is, and structured like `tensor_harmonics_gamma`
(`src/core/xsec.cpp:175-212`): build N from (35), evaluate the two projections of
(43) as cubics in cos φ', and distribute the four harmonics into
`w_avg`, `a1`, `a2`, `a3`.

**Which finite-γ switch gates it: `target_mass`, not `tensor_gamma`.** This is
the one place where the obvious guess is wrong, and getting it wrong is silent.
`tensor_gamma` is the **b-sector (rank-2)** switch and defaults to *false*
(`include/lipolgen/xsec.hpp:204-217`). The **vector** sector's finite-γ switch is
`target_mass`, which defaults to **true** (`include/lipolgen/xsec.hpp:232`, with
the rationale at `include/lipolgen/xsec.hpp:183-186`) and routes `a_parallel` to `a_parallel_exact` =
D_γ(A₁ + η A₂) with a Wandzura–Wilczek g2 table (`src/core/xsec.cpp:145-152`,
`src/core/xsec.cpp:124-131`). By (42) the rank-3 block is the rank-3 *copy* of that A₁/A₂
decomposition, so it must follow the same switch:

* `target_mass == false` → Eq. (45) with `depolarization_d`, i.e.
  `w₃ = λ_e P_e R₃(m) P₃(cos θ_S) D(y) g1_rank3 / F1` — one addend to `w_avg`,
  nothing else, a₁ = a₃ = 0.
* `target_mass == true` (**the code's default**) → the geometry of (43) times
  D_γ(A₁^{(3)} + η A₂^{(3)}) with A₁^{(3)} = (g1_rank3 − γ² g2_rank3)/F1 and
  A₂^{(3)} = γ(g1_rank3 + g2_rank3)/F1 (§5.3). This needs a `g2_rank3` model: a
  `G2Mode`-like choice (a WW-type integral built from `g1_rank3`, or zero) with
  **zero as the stated default** until someone computes one, plus a
  `g2_rank3_scale` mirroring `g2_scale` so the twist-3 sensitivity can be
  measured the way the vector one already is
  (`include/lipolgen/xsec.hpp:188-191`, members at `:227-229`).

Gating the rank-3 block on `tensor_gamma` instead would, at the code's own
defaults, run the vector sector at finite γ while its rank-3 partner ran
massless — silently dropping the very (42) structure this note derives — and
would make a₃ appear only when an unrelated switch was set. **Leave
`tensor_gamma` untouched.** The LT′ piece (cos φ′, O(γ)) belongs beside `a_perp`
under the `with_perp` guard (`src/core/xsec.cpp:272-274`), which is where the
vector sector's own O(γ) transverse partner already lives.

### 6.3 Existing tests that must not move

Every one of these is either an unpolarized-beam identity (protected by the
theorem) or a J = 1 identity (protected by the multipole not existing):

* `tests/test_xsec.cpp:223-240` — `A_zz(θ_S=0)(1 + εR) = −(2/3) b1/F1`, exactly, at every y.
* `tests/test_xsec.cpp:242-255` — "the program sign IS the literature sign".
* `tests/test_xsec.cpp:257-278` — the thirds combination carries the sign of A_zz.
* `tests/test_xsec.cpp:170-193` — the population-averaged cross section is the unpolarized one.
* `tests/test_xsec.cpp:297-316`, `:318-335` — one rank-2 geometry for both spins; the J = 1 HJM transcription digit for digit.
* `tests/test_xsec.cpp:337-360` — the spin-3/2 rate / cos-2φ ratio is spin-independent (a characterization test, and one that (20) now upgrades to a physics statement).
* `tests/test_xsec.cpp:113-132` — the J = 3/2 vector sector, and `w_avg == a1 == a2 == 0.0` for an unpolarized electron on zero rank-2 slots. **This one is the theorem's own regression test** and must be extended, not replaced (see 6.4).
* `tests/test_xsec.cpp:649-683` — `amplitudes` refuses a spin state of the wrong J.
* the whole of `tests/test_tensor_gamma.cpp`, in particular `:151` (γ → 0 collapse of the Cosyn SFs), `:381` and `:406` ([2] Table 1 rows 2 and 3).
* `tests/test_spin.cpp:129-142`, `:204-229` — the J = 3/2 moment round trip and the (0.7, 0.4, 0.2) spin-temperature anchor.
* `tests/test_reference.cpp:79-80` — the reference `b1_32_func` / `delta_32_func` values.

### 6.4 New identity tests

1. **The theorem, numerically.** With `rank3 = true` and a non-zero
   `g1_rank3_func`, an unpolarized electron (`lam_e = 0` *and* `pe = 0`) on any
   J = 3/2 fill at any (x, Q², θ_S, φ') gives w_avg, a1, a2, a3 **bit-for-bit**
   what `rank3 = false` gives. This is the strongest possible statement of §4 and
   is a one-line test.
2. **γ → 0 collapse.** At **`target_mass = false`** (§6.2 — *not*
   `tensor_gamma`, which gates the b-sector and is a different switch) the
   rank-3 sector is exactly (45) at every Q²: a₁ = a₃ = 0 and
   w_avg − w_avg|_{rank3 off} = λ_e P_e R₃ P₃(cos θ_S) D(y) g1_rank3/F1. At
   `target_mass = true` — the code's default — the same difference must
   *approach* that value as γ → 0, with a₁ = O(γ) and a₃ = O(γ³). Mirror of
   `tests/test_tensor_gamma.cpp:152`.
3. **Pure-state weights.** `kern.octupole_moments(m)` = (1, −3, 3, −1) for
   m = (3/2, 1/2, −1/2, −3/2), and `moments_along_axis(1.5, p).octupole` agrees
   with the population-weighted sum. Mirror of `tests/test_xsec.cpp:297-316`.
4. **The population sum rule.** Averaging W(φ') over the 2J+1 pure states with
   equal weight returns the unpolarized cross section for **any** λ_e, because
   Σ_m m/J = Σ_m Q_NN(m) = Σ_m R₃(m) = 0. Extend
   `tests/test_xsec.cpp:170-193` to LI7 with all four sectors live.
5. **Magic angles.** The rank-3 rate term vanishes at P₃(cos θ_S) = 0, i.e.
   θ_S = 39.23°, 90° and 140.77° — the rank-3 analogue of the rank-2 magic
   angle (P₂ = 0 at θ_S = 54.74°) already tested at
   `tests/test_xsec.cpp:388-406`. The two sets are disjoint, so no single axis
   angle kills both rate terms: at 54.74° the rank-3 term survives at
   P₃ = −0.385, and at 90° the rank-2 term survives at P₂ = −1/2. The 90°
   case is nevertheless the useful one — it is the transverse-alignment
   geometry, whose fill has R₃ = 0 identically anyway (§2.4), so the rank-3
   sector is doubly absent from the cos 2φ run plan.
6. **Spin-1 rejection.** A J = 1 kernel given a `g1_rank3_func` must throw, in
   the spirit of `tests/test_xsec.cpp:649-683`: there is no l = 3 multipole for
   spin 1, and silently ignoring the slot would be worse than failing.
7. **Positivity.** Over the toy (x, Q²) grid, (27) holds for m = 3/2 and 1/2
   given the (b1_32, g1, g1_rank3) triple in use — evaluated in the **code's**
   normalization, |(m/J) g1 + R₃(m) g1_rank3| ≤ F1 − Q_NN(m) b1_32 with
   (m/J, R₃, Q_NN) = (1, +1, +1) and (1/3, −3, −1). Writing the bound with [4]'s
   per-flavour weights (m, o(m)) while feeding it the code's tables would
   overweight g1 by 3/2 and g1_rank3 by 3/10 at m = 3/2 — the test would fail on
   a physical input or pass on an unphysical one. (22′) is the map that makes
   the two forms identical. This is the first test in the programme that
   constrains a *physics input* rather than a kernel identity.
8. **[6]'s sum rules**, (18), for any toy `b1_32_func` / `g1_rank3_func` shape
   that is meant to be sea-free: ∫dx b1 = ∫dx g1_rank3 = 0.
9. **The A_∥ bias.** Reproduce (49) end-to-end from generated events: the
   extracted g1/F1 on a spin-temperature fill shifts by (R₃/P_z)·g1_rank3/F1.
10. **The J = 3/2 alignment contrast.** Whichever definition §5.4 adopts, pin it
    and pin the *absence* of the spin-1 2/3: for the ±1 contrast,
    A_T^{(3/2)}(θ_S = 0)(1 + εR) = −b1_32/F1 exactly at every y — the J = 3/2
    mirror of `tests/test_xsec.cpp:223-240`, which pins the spin-1
    −(2/3) b1/F1 and must keep doing so unchanged.

---

## 7. References

1. **W. Cosyn and C. Weiss**, *Semi-inclusive deep-inelastic scattering on a
   polarized spin-1 target. I. Cross section and spin observables*,
   arXiv:2603.23699 (2026). Local: `refs/2603.23699.pdf`. Used: **Appendix D**
   ("Extension to higher-spin targets", p. 25) rules (i)–(v) and Eq. (D2); the
   lepton tensor Eq. (4.5b) — where **λ_e = ±1/2**, half the note's own λ_e,
   §4.1; hermiticity Eq. (4.12); the spin-1 cross-section decomposition
   Eqs. (4.38), (4.39); the list of T-odd λ_e-dependent tensor structures
   Eq. (4.40); the inclusive reductions Eqs. (4.41a–d); the vanishing of the
   T-odd inclusive structures Eq. (4.44), with its **Ref. [62]** (N. Christ and
   T. D. Lee, *Possible tests of C_st and T_st invariances…*, Phys. Rev. **143**,
   1310 (1966)) and the two-photon-exchange exception **Ref. [63]** (A. Afanasev,
   M. Strikman and C. Weiss, *Transverse target spin asymmetry in inclusive DIS
   with two-photon exchange*, Phys. Rev. D **77**, 014028 (2008),
   arXiv:0709.0901); the multipole/pseudotensor assignment.
2. **W. Cosyn, B. Roldan Tomei, A. Sosa, A. Zec**, *Polarization options in
   inclusive DIS off tensor polarized deuteron*, Eur. Phys. J. A **61**, 83
   (2025), arXiv:2410.12764. Local: `refs/2410.12764v1.pdf`. Used: covariant
   spin-1 density matrix Eq. (7); pure-state alignment tensor Eq. (9);
   cross-section decomposition Eq. (10); kinematics Eqs. (11a–c), (12), (13);
   rest-frame projections Eqs. (14a–c); unpolarized SFs Eq. (16); the b-basis
   Eqs. (17a)–(17e); Callan–Gross Eq. (18); asymmetry Eq. (19); polarization
   directions Eqs. (20)–(22); θ_q Eqs. (24a,b); Bjorken limit Eqs. (25)–(27);
   **Table 1**. Reading notes: Eq. (3) states only a *proportionality* and fixes
   no sign — the sign comes from Eqs. (2)+(27) and independently from [7] Eq. (5)
   (§3.4); and **[2]'s x = 2x_d**, so this note's per-nucleon x is [2]'s x_d
   (§5.2).
3. **J. Zhao, Z. Zhang, Z.-t. Liang, T. Liu, Y.-j. Zhou**, *Inclusive and
   semi-inclusive production of spin-3/2 hadrons in e⁺e⁻ annihilation*,
   arXiv:2401.10031, Phys. Rev. D **109**, 074017 (2024). Used: covariant
   spin-3/2 polarization tensors and their orthogonality Eq. (9); light-cone
   decomposition Eqs. (10)–(13); transverse metric Eq. (14). The rank-1/2/3
   parameterization itself is from the companion Phys. Rev. D **106**, 094006
   (2022), arXiv:2206.11742 (their Ref. [84]).
4. **D. Fu, Y. Dong, S. Kumano, J.-J. Xie**, *Generalizing the Soffer Bound:
   Positivity Constraints on Parton Distributions of Spin-3/2 Particles*,
   arXiv:2602.11587. Used: spin-3/2 density matrix Eq. (21); the helicity-state
   PDF table Eq. (23); PDF definitions Eq. (29); positivity Eqs. (32a–c),
   (33a,b). Six quark PDFs (f1, g1, h1, f1LL, g1LLL, h1LLT) and five gluon PDFs.
5. **R. L. Jaffe and A. Manohar**, *Deep Inelastic Scattering from Arbitrary Spin
   Targets*, Nucl. Phys. B **321**, 343 (1989). The complete arbitrary-spin
   basis: 2J+1 quark distributions per flavour in the Bjorken limit. **Not on
   arXiv and not read for this note** — cited through [6] (which attributes its
   Eqs. (19) to it) and [1].
6. **D. Fu, B.-D. Sun, Y. Dong**, *Generalized parton distributions in spin-3/2
   particles*, Phys. Rev. D **106**, 116012 (2022), arXiv:2209.12161. Used: §3.1
   ("four independent structure functions … at leading twist"); the parton-model
   definitions Eqs. (19a)–(19d); the forward-limit GPD relations Eqs. (20a)–(20d);
   the sum rules Eq. (21).
7. **W. Cosyn, Y.-B. Dong, S. Kumano, M. Sargsian**, *Tensor-polarized structure
   function b1 in the standard convolution description of the deuteron*,
   Phys. Rev. D **95**, 074036 (2017), arXiv:1702.05337. Local:
   `refs/1702.05337.pdf`. Used: the HJM parton-model definition of the spin-1 b1,
   Eq. (5), δ_T f = f⁰ − (f^{+1} + f^{−1})/2, which fixes the sign in §3.4.
8. **HERMES Collaboration (A. Airapetian et al.)**, *First measurement of the
   tensor structure function b1 of the deuteron*, Phys. Rev. Lett. **95**, 242001
   (2005), arXiv:hep-ex/0506018. The A_zz = −(2/3) b1/F1 convention the
   generator's `TENSOR_LL_SIGN = −1` is anchored on.
9. **P. Hoodbhoy, R. L. Jaffe, A. Manohar**, *Novel effects in deep inelastic
   scattering from spin-one hadrons*, Nucl. Phys. B **312**, 571 (1989). The
   spin-1 b1…b4 basis and the c_m = 3m² − 2 convention the generator's J = 1 path
   transcribes (`include/lipolgen/spin.hpp:23`,
   `include/lipolgen/asymmetries.hpp:9-21`). Not read directly; used through [2]
   and [7].
10. **W. Cosyn and C. Weiss**, *Semi-inclusive deep-inelastic scattering on a
    polarized spin-1 target. II. Deuteron and spectator nucleon tagging*,
    arXiv:2603.23700 (2026), JLAB-THY-26-4661. Local: `refs/2603.23700.pdf`.
    Used: **§VII** ("Conclusions and extensions"), item (i)(c) — *"What is the
    size of T-odd structure functions, which are zero in the IA and provide a
    sensitive test of FSI?"* — the open question for the *tagged* case quoted in
    §4.6 item 2. (That quote was previously mis-attributed to [1] §VI C; it is
    not in [1].)

---

## Appendix: claims in this note that are *not* backed by a source

Listed so that the theory co-author can check them rather than inherit them.

* **Eq. (11)**, the rank-3 tensor R_ijk = (R₃/2)(5n_in_jn_k − n_iδ_jk − n_jδ_ik −
  n_kδ_ij), and the projections (13), (43). These are the unique symmetric
  traceless rank-3 tensor built from n̂, normalized to match Eq. (10) and the
  code's R₃ of Eq. (7). No source was found that writes them in this frame.
* **Eq. (32)–(33)**, the inclusive structure-function count per multipole. Derived
  here from photon-helicity conservation and the symmetry of ε*_{λ'} ⊗ ε_λ;
  checked against the known J = 1/2 (4) and J = 1 (8) answers. [5] presumably
  contains the general result; [5] was not read.
* **Eq. (40)**, "the rank-3 block has exactly two inclusive structure functions".
  This follows from (32) but is not quoted from any paper.
* **Eq. (30′)**, the parity relation W_{−λ',−λ}[ρ_l] = (−1)^l W_{λ'λ}[ρ_l] for
  the m_l = 0 component, and with it the conclusion that T and L select even l
  and T′ selects odd l by parity alone. Standard, but written out here rather
  than quoted; the only ingredient is w_l(−M) = (−1)^l w_l(M) for the m_l = 0
  multipole weight, which is Eq. (8) read as a polynomial in M.
* **Eq. (42)** is *no longer* on this list, and the earlier draft's disclaimer
  ("the absolute normalizations were not verified against any source") was too
  weak. The channel identity is [1] App. D rule (iii) with [1] Eqs. (4.38),
  (4.39), and the constants 2 and −2γ are those of [1] Eqs. (4.41c,d), adopted
  as the *definition* of (g1_rank3, g2_rank3). What remains unbacked is the
  assumption that rule (iii) fixes the constant and not merely the ε and
  azimuthal dependence — the rule as printed says "azimuthal and ε dependence".
  If the rank-3 pair is ever computed independently, the overall constant is the
  thing to check.
* **Eq. (24)**, b1_32|_code = −b1^{[6]}. Each half of the chain is verified
  ([7] Eq. (5) for the spin-1 sign; [6] Eq. (19b) for the spin-3/2 sign; the code
  and its tests for `TENSOR_LL_SIGN`), but the claim that [5] itself uses [6]'s
  sign for arbitrary spin could not be checked, [5] being unavailable.
* **Eq. (49)** and the bias coefficients R₃/P_z. The arithmetic is exact given
  (45); the *size* of the resulting bias depends on g1_rank3/g1, for which no
  calculation exists for any nucleus.
