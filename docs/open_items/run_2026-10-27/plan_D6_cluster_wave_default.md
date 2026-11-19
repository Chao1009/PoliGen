<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Plan D6 — the VMC AV18+UX cluster waves become the shipped default (design only)

Run 2026-10-27, task D6, written by an Opus worker on `claude/run-2026-10-26`
at `ca9b23a` (clean tree). **Design only: no repository file other than this
record was edited.** Everything measured here was measured on this machine
(4 shared cores, `source env.sh`, the in-tree build of `ca9b23a`), or on a
throw-away copy of `ca9b23a` exported with `git archive` into the session
scratch directory and built there (§4.6) — never in the working tree.

The author's decision (`DECISIONS.md` of this run, decision 3):

> if that's the well accepted wave default, then do it, otherwise if there
> is a better model favored by the measured D/S ratio we can use it.

The task's framing, taken as the authority for the default switch: AV18 with
the VMC cluster waves becomes the shipped default (`ClusterWaveSource::VmcAV18`
— the enum value is spelled `VmcAV18`, not `Vmc`; `--cluster-wave vmc`), the
Hulthén forms stay available as `--cluster-wave hulthen`, and CD-Bonn is the
recorded alternative for the deuteron.

Conventions of this file. **READ** marks what was read in a file or paper
(with the place); **MEASURED** what was run here (with the window);
**INFERRED** what follows from those and was not itself run or read. Numbers
are copied from the tool output; ratios are computed from the printed values.
**Every line number is at `ca9b23a`.** Line numbers are written as a file name
plus "line N", deliberately not in the backticked path-colon-line form that the
records citation gate (`validation/check_physics_channels_links.py
--records`) parses. The implementation this plan designs moves those lines,
and this record states the tree it read, not the tree I6 will leave.

---

## 0. In one paragraph

`--cluster-wave` today switches three things on the three tagged channels and
nothing else (§1): the cluster relative wave (⁶Li α–d, ⁷Li α–t, d–p), the
⁶Li α tag's embedded deuteron (`dis_target`) and the T1 breakup's deuteron
(`BreakupOptions::source`, forwarded by the `Pipeline`). The recommended design
(§3) flips the **run** default — `PipelineConfig::cluster_wave`, hence
`lipolgen-run`, `lg.run`, `lg.make_config` — through one named constant
`CLUSTER_WAVE_DEFAULT`, and **keeps** the analytic family as the argument
default of the three channel factories and of `BreakupOptions` (their β and
P_D are that family's parameters, and the polligen port gate `tagged.json` is
pinned on them; no `validation/reference` JSON depends on the run default). That is the same split the D7
plan chose for the structure-function objects ("Class-level null defaults stay
the toy, and say so", `plan_D7_backend_defaults.md` §2 refinement 1). A flip
experiment on a scratch copy measured the cost of both designs (§4.6): the
recommended one breaks **4 C++ cases and 15 pytest tests**, the
factories-too one **27 C++ cases** (the polligen port gates and the P_D = 0
limits among them — 22 267 assertions in one case, because `li6_alpha_channel(β,
0.0)` silently became VMC and ignored its `p_d`) and 24 + 5 pytest failures in
a subset where the recommended design fails 3. `validate()`'s off-channel
refusal must follow the default (otherwise every default inclusive run
throws); an explicit `--cluster-wave` off the tagged channels stays refused by
`make_config` and the CLI, which can tell an explicit request from the default
(§3.2). No benchmark verdict moves (§7): the four T2 deuteron rows keep their
statuses under either reading of their declared rules; Deeps' headline moves
627.85 → 579.12 on 187 dof because its rule names "the default wave". The
first risk, and a question for the author rather than a decision (§11.1):
after D6 the default tagged ⁶Li α run and the default inclusive ⁶Li run
describe ⁶Li with different wave functions — ×1.0935 on the whole-nucleus
polarization (0.887076 against `LI6_CLUSTER_POLARIZATION` = 0.811228) and
×1.0658 on the rank-2 transfer (0.982576 against 0.921947).

---

## 1. What the flag switches today — and what it does not (READ)

### 1.1 The three reads of `PipelineConfig::cluster_wave`

| channel | object | `hulthen` (today's default) | `vmc` (`VmcAV18`) | where |
|---|---|---|---|---|
| `tagged-6Li-alpha` | α–d relative wave | S+D analytic pair at β = 0.30 GeV, P_D = `P_D_LI6` = 0.0867 | ANL VMC tables: magnitudes `vmc/momenta/li6_ad1.momentum` (AV18+UX), S–D sign `vmc/li6_alpha_d/li6.ad` (AV18+UIX, 2004); P_D = `VMC_P_D_LI6` = 0.019355; **ignores β and `p_d`** | `src/core/tagged.cpp` lines 243-276 |
| `tagged-6Li-alpha` | embedded deuteron, struck-cluster g₁ (`dis_target`) | `DEUTERON()`, eff_pol 0.9325 (P_D_DEUTERON = 0.045) | `DEUTERON_AV18()`, eff_pol 0.913600 (P_D 0.057600 from `fdeut.av18`) | `src/core/tagged.cpp` lines 272-273 |
| `tagged-6Li-alpha` (T1) | the breakup's deuteron (struck-nucleon momentum **and** spin draw) | S+D Hulthén at 0.045 | AV18 `fdeut.av18` | `src/core/pipeline.cpp` line 1686 (`bo.source = cfg_.cluster_wave`); `include/lipolgen/breakup.hpp` lines 207-230 |
| `tagged-7Li-alpha` | α–t relative P wave | Hulthén P wave at β | `vmc/momenta/li7_at3.momentum` (3/2⁻ ground state); ignores β for the relative wave | `src/core/tagged.cpp` lines 278-305 |
| `tagged-d-p` | p–n relative S+D wave | analytic pair at β = 0.30, P_D = 0.045 | the exact AV18 deuteron `data/vmc/deuteron/fdeut.av18`, P_D 0.057600 | `src/core/tagged.cpp` lines 307-329; `src/core/pipeline.cpp` lines 1532-1533 |

READ, `include/lipolgen/pipeline.hpp` lines 766-768: "`Hulthen` is the default and
keeps every published number bit-for-bit; `VmcAV18` swaps in the ANL VMC
tables and then IGNORES `cluster_beta` and `p_d`". The knob table already
encodes the reach (`src/core/pipeline.cpp` lines 2800-2873): `cluster_wave` READ on
the three tagged channels and REFUSED elsewhere; `cluster_beta` NOT READ at
`vmc` except on ⁷Li at T1 (the triton breakup stays analytic at β); `p_d` NOT
READ at `vmc`.

### 1.2 What it does not reach (READ; each is out of D6)

- **The inclusive and coherent channels**: `validate()` refuses the field
  there (`src/core/pipeline.cpp` lines 482-506); the inclusive ⁶Li constants
  `LI6_CLUSTER_POLARIZATION` = 0.811228 (`include/lipolgen/beams.hpp` lines 137-138)
  and `LI6_B1_RANK2_TRANSFER` = 0.921947 (`include/lipolgen/constants.hpp` line 133)
  do not follow it — the subject of §11.1.
- **The ⁷Li triton's internal structure**: `TRITON()` on either setting
  (`src/core/tagged.cpp` lines 296-303, "this tree carries no AV18 triton wave
  function"); the T1 triton breakup (sequential Hulthén at β, or
  `--triton-sf ciofi-simula`).
- **`DeuteronConvolutionB1`'s wave** (`kFdeutFile`, AV18 already; registry
  row 5) and the `li6-convolution` α–d source (`FitRescaled`, VMC already;
  registry row 16).
- **The FSI weight** is built from the run's own `TaggedModel`
  (`src/core/pipeline.cpp` lines 1552-1562), so it follows the wave automatically —
  with a measured consequence under VMC (§4.4).

### 1.3 Sibling defaults the same reasoning would cover — reported, not folded in

| default | today | what the D6 reasoning would say | why not folded in |
|---|---|---|---|
| `--triton-sf` (`TritonSfChoice::Hulthen`, `include/lipolgen/pipeline.hpp` lines 335-337, `python/lipolgen/cli.py` line 805) | sequential two-body Hulthén triton decay, "CRUDE AND FLAGGED" (`docs/PHYSICS_CHANNELS.md` line 304) | a realistic spectral function (Ciofi degli Atti–Simula, S₀ = 0.6525 untuned) is the better-accepted model | the measured D/S ratio says nothing about the triton; the decision is about the cluster wave |
| the factories' and `BreakupOptions`' own `source` defaults | `Hulthen` | one default everywhere | §3.1: their β/P_D are the analytic family's parameters; measured cost §4.6 |
| `P_D_DEUTERON` = 0.045 → `DEUTERON()` / `DEUTERON_VECTOR_POLARIZATION` = 0.9325 (inclusive deuteron species, and a factor of `LI6_CLUSTER_POLARIZATION`) | scenario value; `beams.hpp` line 100 calls it "(AV18-like)" although AV18's is 0.057600 | one Hamiltonian would read 0.913600 | P_D is not an observable (READ, `docs/benchmarking/02_data_nucleon_deuteron.md` line 288: "P_D. The D-state probability is **not an observable**"); it is part of the §11.1 question, not of D6 |
| `DeuteronConvolutionB1::Options::wave` = AV18 (registry row 5) | AV18 | the measured η favours CD-Bonn numerically (−0.07σ against −1.39σ) | both pass the declared 2σ; already the D6 Hamiltonian |
| `--rc-c0-shape ho` (the ⁶Li C0 shape of the RC tail) | harmonic oscillator | `vmc-ft` is the VMC point-proton density | declared a **band edge**, not a default choice (`README.md` line 54) |

---

## 2. The decision's evidence, and the D/S question

READ:

- Rodning & Knutson, PRC 41 (1990) 898 (the copy fetched 2026-10-26,
  sha256 `4dbfb75f…ee0ff161`), abstract p. 898: "A weighted average of the
  individual values gives the result η = 0.0256 ± 0.0004, where the quoted
  uncertainty includes systematic and statistical contributions", and
  Eq. (10), p. 906: "η = 0.0256 ± 0.0004".
- R. B. Wiringa, V. G. J. Stoks, R. Schiavilla, PRC 51 (1995) 38
  (`PolarizedLithiumSim/refs/nucl-th_9408016.pdf`, sha256
  `44f551a2…10996f5`), text p. 16 of the arXiv version: "The asymptotic
  S-state normalization, A_S, the D/S-ratio, η, and the deuteron radius,
  r_d, all come out close to the experimental values"; Table X ("Static
  deuteron properties", p. 34): η experiment 0.0256(4), Argonne v18 0.0250.
- R. B. Wiringa, R. Schiavilla, S. C. Pieper, J. Carlson, PRC 89 (2014)
  024305 (`PolarizedLithiumSim/refs/1309.3794.pdf`, sha256
  `00fd24f1…a6f28e20b`), abstract p. 1: "The wave functions have been
  generated for a Hamiltonian containing the Argonne v18 two-nucleon and
  Urbana X three-nucleon potentials. … The nucleon-cluster momentum
  distributions include dp in ³He, tp and dd in ⁴He, αd in ⁶Li, αt in ⁷Li,
  and αα in ⁸Be."
- The tree's own row, `validation/benchmarks/t2_deuteron_static.py` lines 83-88
  (declared before the verdict): "eta: ASSERTED. |eta_model − 0.0256| <=
  2 sigma (sigma = 0.0004) for EVERY shipped wave function"; the committed
  REPORT row: "cdbonn 0.025571 (−0.07 sigma), av18 0.025045 (−1.39 sigma),
  hulthen 0.032478 (+17.20 sigma)". `fdeut.av18`'s own header prints
  `eta 0.025045`, `dstate 0.057599`.
- The sign reference is not the same three-body force:
  `include/lipolgen/cluster.hpp` line 407, "the ANL VMC tables (AV18+UX
  magnitudes, AV18+UIX signs)".

INFERRED (the reading the decision text in §12 carries): the measurement
rejects the Hulthén pair at +17.2σ and accepts both realistic NN potentials at
the declared 2σ. It favours CD-Bonn numerically (−0.07σ against −1.39σ) but
does not discriminate between the two at the declared tolerance. CD-Bonn
exists in the tree only as a deuteron (`cdbonn_wave()`, `kCdBonn`); there is no
CD-Bonn ⁶Li α–d or ⁷Li α–t table. The ANL tables share the AV18 two-nucleon
force across the α–d overlap, the α–t overlap and `fdeut.av18`. Choosing
CD-Bonn for the d–p control alone would put two wave-function families in one
programme, which the C5.4/C5.5b rule forbids (READ, registry row 9: "the
control's job is to be the Cosyn–Weiss tagged limit of the same analytic
family the ⁶Li channel is built from"; row 11: "§C5.4's own rationale forbids
it (*switching only the control would make the control and the channel it
controls two different wave-function families*)"). So AV18/VMC satisfies the
author's first branch, and the second branch ("a better model favored by the
measured D/S ratio") points at CD-Bonn only for a deuteron-only choice the
tree cannot make consistently. CD-Bonn stays the recorded alternative. The
"well accepted" premise is the author's. Beyond the two papers quoted, it was
not surveyed here.

---

## 3. Design

### 3.1 DD-1 — where the default lives: the run surface, not the building blocks (recommended)

- **A (recommended).** Flip `PipelineConfig::cluster_wave`
  (`include/lipolgen/pipeline.hpp` line 823) to a new
  `inline constexpr ClusterWaveSource CLUSTER_WAVE_DEFAULT =
  ClusterWaveSource::VmcAV18;` declared beside the enum in `cluster.hpp`
  (one home, the `cluster_wave_name` pattern of `cluster.hpp:409-415`). Keep
  `Hulthen` as the argument default of `li6_alpha_channel`,
  `li7_alpha_channel`, `deuteron_channel` (`include/lipolgen/tagged.hpp` lines
  216, 223 and 245-246; `python/bindings.cpp` lines 2826, 2832 and 2841) and of
  `BreakupOptions::source` (`include/lipolgen/breakup.hpp` line 230). Their docs
  say so: "the analytic family's constructor default; the RUN default is
  `CLUSTER_WAVE_DEFAULT`". `BreakupOptions::source` is overwritten by the
  `Pipeline` exactly as `beta` and `f2` are (`src/core/pipeline.cpp` lines 1680-1687),
  and the binding docstring (`python/bindings.cpp` lines 4474-4476) names only
  `beta` and `f2`. It gains `source`.
- **B (rejected).** Flip the factories and `BreakupOptions` too. MEASURED
  (§4.6): `li6_alpha_channel(β, p_d)` positional calls would silently return
  the VMC channel and ignore `p_d` and β (the P_D = 0 limits in
  `tests/test_tagged.cpp` lines 341, 514 and 834, `tests/test_fsi.cpp` line 48). The
  polligen port gates (`tagged.json`, polligen's own Hulthén) would compare
  VMC to Hulthén. The four T2 deuteron harnesses would label AV18 as
  "hulthen", and `t2_deuteron_static.hulthen()` raises "deuteron_channel() is
  no longer the analytic Hulthen S + D pair". The sibling prototype keeps
  the same split: READ, `PolarizedLithiumSim/evgen/polligen/tagged.py:365-367`
  "The tables are selectable (`li6_alpha_channel(wave='vmc')`) and are NOT
  the default: every published number is the Hulthen pair, bit for bit".
- Consistency: `plan_D7_backend_defaults.md` §2 chose the same split for
  the SF objects ("Class-level null defaults stay the toy, and say so … their
  comments change from 'the DEFAULT … every published number' to 'the
  class-level toy; the run's default is set by `PipelineConfig`'").

### 3.2 DD-2 — the off-channel refusal must follow the default

Today: `if (cluster_wave != ClusterWaveSource::Hulthen && !is_tagged(channel))
throw` (`src/core/pipeline.cpp` line 497). With `VmcAV18` the default this throws on
**every** default inclusive and coherent run (MEASURED, §4.6: the flip made
`CHECK_NOTHROW(cfg.validate())` throw on both).

- **Recommended:** refuse a value *moved from the default*,
  `cluster_wave != CLUSTER_WAVE_DEFAULT && !is_tagged(channel)`. That is the
  rule every RC sub-knob refusal already follows: each knob is compared with a
  default-constructed `const RcOptions od;` (`src/core/pipeline.cpp` lines
  718-748, `{"--rc-tail-model", rc_options.tail_model != od.tail_model}`).
  D7 adopts the same rule for its kernel check (§2 refinement 6(b): "compared
  with `PipelineConfig{}`"). The message is reworded: an explicit `hulthen` on
  inclusive is refused for being unread, and the inclusive constants are named
  as the Hulthén-era cluster product.
- **And, to keep registry row 11's protection for an explicit request**,
  which a value rule cannot see once `vmc` *is* the default
  (READ, `src/core/pipeline.cpp` lines 485-489: "they set the flag wanting a 'VMC 6Li'
  inclusive row, nothing refused it, and they got LI6_CLUSTER_POLARIZATION =
  0.811228 with no warning"), `make_config` refuses an explicitly passed
  `cluster_wave` (either value) on a non-tagged channel. It raises
  `RuntimeError` with a message that starts "PipelineConfig:". Then the CLI's
  type-based catch (`python/lipolgen/cli.py` lines 1270-1282, the `except` at
  line 1312) turns it into a one-line
  refusal, and `validation/price_defaults.py`'s `_attempt`/`_is_environment`
  classify it as a refusal unchanged. The CLI's `DEFAULTS["cluster_wave"]`
  becomes `None` (`cli.py:805`), so "not given" reaches `make_config` as `None`
  and leaves the C++ default. INFERRED from the failure mode measured in
  §4.6, not run with the guard: the inclusive probe
  `("cluster_wave", "vmc", "refused")` in `price_defaults.py` (lines 544, 746
  and 885) and the inclusive cells of the knob matrix then keep passing as
  written.
- Residual hole, stated: a C++ caller (or a Python caller assigning
  `cfg.cluster_wave` directly) who sets the default value on an inclusive
  config is indistinguishable from one who did not. It is recorded in the
  knob table as `cluster_wave = vmc (refused, at default)` with the reason
  text, the same as every other knob at its default.
- Alternative (not recommended, priced by reading): a tri-state
  (`std::optional<ClusterWaveSource>`) would make the explicit request visible
  in C++ too. It changes a public field's type, the binding (None ↔ nullopt)
  and about six tests that compare `cfg.cluster_wave` to an enum value.

### 3.3 DD-3 — one home for the default name

`m.attr("CLUSTER_WAVE_DEFAULT")` in the bindings. The CLI help and the
docstrings read `cluster_wave_name(CLUSTER_WAVE_DEFAULT)` instead of retyping
"vmc" (the rule `cli.py:806-808` already states for the RC anchors: "read them
rather than retyping"). Same arrangement D7 recommends for `unpol_sf`/`pol_sf`
(`plan_D7_backend_defaults.md` §2 item 4).

### 3.4 DD-4 — the C++ example generator

`examples/generate_tagged.cpp` builds `PipelineConfig cfg;` and sets
`cfg.cluster_beta = args.beta; cfg.p_d = args.p_d;` (lines 145-146) from `--beta`
/ `--pd` (lines 87-88). Under the VMC default both are silently not read, and the
banner prints `beta = %.2f` from the argument (lines 164-168). Recommended:
`--cluster-wave hulthen|vmc` (default from `CLUSTER_WAVE_DEFAULT`), refuse
`--beta`/`--pd` unless `hulthen`, and print the wave. Cost: lines move above
the citations `examples/generate_tagged.cpp` lines 119-124 (a range) and line 123
(an unnamed point), made by `docs/theory/SPIN32_FINITE_GAMMA.md` lines 208 and
151, and the named point at its line 79 (`--channel`), made by
`docs/PHYSICS_CHANNELS.md` line 520 (§9.3).

### 3.5 DD-5 — optional, recommended

1. A banner line on every tagged run naming the wave (today a default run
   names it nowhere; the KNOB PROVENANCE block lists only knobs set away from
   their default). Example: `cluster wave: vmc (ANL VMC AV18+UX; alpha-d P_D
   0.019355, embedded AV18 deuteron P_D 0.057600) -- --cluster-wave hulthen
   for the analytic pair`.
2. The missing-table error names the escape. MEASURED with an empty
   `$LIPOLGEN_DATA_DIR`: `RuntimeError: cannot open ANL table: <dir>/vmc/
   li6_alpha_d/li6.ad` (and `.../deuteron/fdeut.av18`,
   `.../momenta/li7_at3.momentum`). Before D6 a default tagged run needed no
   file.
3. A conventional `cmake --install` installs no `data/`
   (`CMakeLists.txt` lines 328-338; only the wheel branch copies it, lines 325-327). A
   C++ program run against an installed library outside the source tree, with
   no `$LIPOLGEN_DATA_DIR`, now throws at its first default tagged run.
   Options: install `data/` in that branch too (one line), or document it.

---

## 4. Measurements (this machine, 2026-10-27)

### 4.1 First-use load time and memory of the tables (MEASURED)

A fresh Python process per row (`lg._lipolgen` at `ca9b23a`), timing the first
and a second factory call and one `TaggedModel` build, with VmRSS deltas:

| channel | wave | 1st factory call | 2nd call | `TaggedModel` | ΔRSS factory | ΔRSS model | max RSS |
|---|---|---|---|---|---|---|---|
| ⁶Li α | hulthen | 0.12 ms | 0.05 ms | 4.5 ms | +0.1 MB | +4.6 MB | 39.3 MB |
| ⁶Li α | vmc | **6.7 ms** | 0.08 ms | 2.9 ms | **+2.4 MB** | +2.3 MB | 39.5 MB |
| ⁷Li α | hulthen | 0.12 ms | 0.05 ms | 3.6 ms | +0.1 MB | +3.8 MB | 38.5 MB |
| ⁷Li α | vmc | 0.17 ms | 0.04 ms | 2.7 ms | +0.1 MB | +3.8 MB | 38.6 MB |
| d–p | hulthen | 0.13 ms | 0.05 ms | 6.7 ms | +0.1 MB | +4.6 MB | 39.4 MB |
| d–p | vmc | **4.7 ms** | 0.13 ms | 2.6 ms | **+2.3 MB** | +2.5 MB | 39.6 MB |

CLI, 4 000 events (`lipolgen-run … --cluster-wave hulthen|vmc`): setup 0.07 s
on both waves; wall 0.22–0.24 s; max RSS 49.7 → 49.9 MB (⁶Li α) and 49.9 →
50.2 MB (d–p). **The tables cost ≲ 7 ms and ≈ 2.4 MB on first use; process
peaks differ by ≤ 0.3 MB.** Not a cost worth designing around.

### 4.2 The size of the change on the tagged channels (MEASURED)

The extraction is `validation/price_defaults.py`'s own (`_rp_fraction`,
`acceptance_weights(…, n_phi = 32)`, `azz_tensor_curve_weighted`), imported
read-only. Window: config 1 (e 10 GeV × ion 99.5 GeV/u), YR high acceptance,
`tensor-thirds` (⁶Li, d) / `helicity-flip` P_e 0.7 (⁷Li), seed 20261026.

**Confirmation of the price table (N = 200 000):** the hulthen and vmc rows of
`../run_2026-10-26/phase_A_baseline.md` §4.2 reproduce **exactly**. ⁶Li α:
A_zz^tag −1.2069113 → −0.5191289 at k = 0.19793 GeV (×0.43013) and −0.9532827
→ −0.1992630 at 0.30115 GeV (×0.20903); RP tag fraction 0.02466 ± 0.00035 →
0.03422 ± 0.00041 (×1.3877); vector dilution 0.8699432 → 0.9709661
(×1.11613), tensor 0.9219490 → 0.9825759 (×1.06576); ⟨k⟩ (all events)
0.12200 → 0.12243 GeV; σ 591846.1688784088 pb on both. ⁷Li α: ⟨k⟩ 0.28852 →
0.18608 GeV (×0.6449); RP 0.97386 → 0.99839; σ identical.

**The d–p control, not in the price table (N = 200 000, new):**

| d–p control | hulthen | vmc (AV18) | ratio |
|---|---|---|---|
| A_zz^tag at k = 0.19793 GeV | −0.8763449 | −0.9624645 | ×1.0983 |
| A_zz^tag at k = 0.30115 GeV | −0.5478468 | −0.6680525 | ×1.2194 |
| vector / tensor dilution | 0.9324961 / 0.9594889 | 0.9135960 / 0.9481463 | ×0.97973 / ×0.98818 |
| RP tag fraction | 0.02396 ± 0.00034 | 0.02107 ± 0.00032 | ×0.879 |
| ⟨k⟩ (all events) | 0.09849 ± 0.00019 GeV | 0.09530 ± 0.00021 GeV | ×0.968 |
| σ | 585044.2284977469 pb | 585044.228497747 pb | 1 ulp |

At N = 4 000 (the `lipolgen-run` CLI and `lg.run`, same seed): ⁶Li α tag
fraction 0.0225 → 0.0278 and ⟨k⟩ 0.1180 → 0.1200 GeV; d–p tag fraction 0.0232
→ 0.0200 and ⟨k⟩ 0.0962 → 0.0925 GeV (binomial σ on a fraction ≈ 0.0025,
on ⟨k⟩ ≈ 0.0015 GeV). A_zz^tag and the dilutions are model-level and do not
depend on N.

### 4.3 The Deeps row if its declared "default wave" is AV18 (MEASURED)

The committed `t2_deeps_spectator_tail.py` was imported read-only with its
module `DEFAULT_WAVE` set in memory, and `transverse(…, wave=…)` was called
explicitly. **Implementer's trap:** `def transverse(pts, tree,
wave=DEFAULT_WAVE)` (line 361) binds its default at import, so a monkeypatch does
not reach it. Editing the constant in the source does.

| | hulthen (committed headline) | av18 |
|---|---|---|
| PWIA-window LC shape χ²/ndf | 627.85/187 (p = 4.9e−49) | **579.12/187 (p = 1.1e−41)** — FAIL either way |
| data/(s·model) per p_s 0.30 … 0.53 GeV/c | 0.790 / 0.954 / 1.271 / 1.715 / 1.409 | 1.1395 / 1.1964 / 1.1737 / 1.0350 / 0.6408 |
| LC × GlauberFsiWeight (recorded) | 422.923/187 | 582.516/187 |
| transverse window data/PWIA | 0.760 / 1.198 / 1.760 / 2.642 / 2.957 | 1.1714 / 1.7026 / 1.9900 / 1.9391 / 1.5027 |
| transverse window tree Glauber | 1.681 / 3.115 / 6.048 / 12.404 / 18.936 | 1.9298 / 3.4634 / 5.4336 / 7.3071 / 7.1735 |
| same side of 1 | 4 of 5 | 5 of 5 |
| FSI weight table clipped fraction | 0.0145 | 0.0135 |

### 4.4 Published default-run numbers that ride on the wave (MEASURED)

**FSI, "the two variants are NOT one"** (`docs/USAGE.md` lines 1679-1697,
`README.md` line 49, `tests/test_fsi.cpp` lines 384-457, `python/tests/test_meta_provenance.py`):
⁶Li α, 20 000 events, seed 1234, `tensor-thirds` P_z 0.7 / P_zz 0.6 / P_e 0.7,
σ_XN = 40 mb.

| | published (2026-09-04) | hulthen, re-measured | **vmc** |
|---|---|---|---|
| Σw cluster / nucleon | 10419.07 / 11632.10 | 10403.75 / 11584.83 | **15989.35 / 16428.71** |
| `survival()` cluster / nucleon | 0.520239 / 0.582899 | 0.520239 / 0.582899 | **0.797094 / 0.819119** |
| per-event w, cluster | 0.0208 … 1.4043 | 0.0181 … 1.4138 | **0.0681 … 37.92** |
| per-event w, nucleon | 0.2139 … 6.9523 | 0.2124 … 6.9625 | **0.2026 … 50.0 (the `w_max` ceiling, 2 events)** |
| w_n/w_c median / max | 0.897 / 68.52 | 0.8975 / 69.89 | **0.9913 / 10.14** |
| events differing by > 1 % | 99.50 % | 99.425 % | **64.515 %** |

The published column is a pre-2026-09-06 sample (measured 2026-09-04, before
the S–D sign fix re-drew the spin-1 (M, cos θ_k) sample). Its model-level
survivals reproduce; its sample-level numbers do not, a pre-existing staleness
D6's re-measurement would fix. **Under VMC the FSI weight table clips:**
`clipped_grid_fraction()` 0.70 % (cluster, 40 mb), 5.25 % (nucleon, 40 mb),
0.34 % (nucleon, 20 mb), 0 (cluster, 20 mb); 0 everywhere on Hulthén. READ:
the weight is clipped at `w_max` = 50 because "the FSI/IA ratio is UNBOUNDED
wherever the impulse amplitude has a node" (`include/lipolgen/fsi.hpp` lines
324-330, worded so at `docs/PHYSICS_CHANNELS.md` line 347), and the VMC α–d S
wave has "one S node at 0.678 fm⁻¹ = 0.134 GeV" (`docs/CONVENTIONS.md` lines
363-364). INFERRED: that node is where the VMC weight clips; it was not
localised here.
The 20 mb survivals: hulthen 0.6712 / 0.6555, vmc 0.8636 / 0.8576.

**RC band clamp** (`docs/USAGE.md` lines 2226-2232 "1.2 % of ⁶Li tagged-alpha events
and 2.6 % of ⁷Li"; line 2242 "245 of 20 000 = 1.2250 % … seed 1";
`docs/PHYSICS_CHANNELS.md` line 855): `--rc tensor-band`, 20 000 events. ⁶Li α,
CLI-default fill, seed 1: **237 (1.185 %) → 14 (0.070 %)**; seed 1234: 1.095 %
→ 0.090 %. ⁷Li α (`helicity-flip`, `--pe 0`), seed 1: 698 (3.49 %) → 700
(3.50 %); seed 1234: 3.345 % → 3.48 %. The 245 at seed 1 was not reproduced
with the CLI-default fill (237); the published recipe is not stated.

**Tag fractions against the Python** (`docs/USAGE.md` lines 1994-1998;
`tests/test_pipeline.cpp` lines 398-465): under VMC, at 400 000 events, ⁶Li α 0.0334375
(YR) / 0.287445 (tagging, legacy levers) / 0.248635 (tagging) against the
Python's Hulthén 0.0247 / 0.3061 / 0.2545. ⁷Li α YR 0.99846 at all three
configurations, against the Python's 0.9655 / 0.9735 / 0.9788. These are port
comparisons with a Hulthén prototype, so they stay on Hulthén (§6).

### 4.5 The sibling prototype's VMC channel, read-only (MEASURED)

`polligen.tagged` at the sibling's HEAD (`4726812`) was imported with
`python3 -B` / `PYTHONDONTWRITEBYTECODE=1`. Nothing was written: `git status`
was clean and no new `*.pyc` appeared. Hulthén ⁶Li α: vector 0.8699431788603683,
tensor 0.9219489769684724 — identical to LiPolGen's to every digit. VMC ⁶Li α
(`wave='vmc'`): vector 0.970966091698 against LiPolGen's 0.970966086466
(5.4 × 10⁻⁹ relative), tensor 0.982575875485 against 0.982575872345
(3.2 × 10⁻⁹), P_D identical (0.019354933183319687). Its embedded deuteron
stays the scenario 0.9325 (LiPolGen's `VmcAV18`: 0.913600, the C5.5b
repair), and it refuses ⁷Li `vmc`. So an external port gate for the new
default's α–d relative motion is feasible at a declared tolerance of order
1e−8, not at the rtol 1e−12 the Hulthén blocks meet (§11.6).

### 4.6 The flip experiment — what breaks (MEASURED)

`git archive ca9b23a` into the scratch directory, built there (`cmake --build
-j2`, 2 min 35 s, peak RSS 4.4 GB), environment from `env.sh`. **The full
suites were not run** (task rule). The files run were chosen by grep: every
C++ file that builds a `PipelineConfig`/`BreakupOptions` or calls a factory,
and every pytest file that runs a tagged channel or names `cluster_wave`.

**Design A** — `pipeline.hpp` line 823 → `VmcAV18`; `pipeline.cpp` line 497
refuses `!= VmcAV18` off the tagged channels; `cli.py` line 805
`cluster_wave=None`; no
`make_config` guard:

- C++, 10 files (`test_pipeline`, `test_breakup`, `test_tagged`, `test_t2`,
  `test_fsi`, `test_b1_nuclear`, `test_rc_pipeline`, `test_rc`,
  `test_triton_sf`, `test_hepmc`): **4 failing cases, 11 assertions** (6 + 5) —
  `test_pipeline.cpp` "6Li alpha Roman-Pot tag fractions reproduce the Python"
  (0.0334375 vs 0.0247 at rtol 0.10), "A_zz^tag(k) closes on the
  acceptance-weighted truth" (`REQUIRE(ndf >= 8)`: 7 at YR; χ²/ndf 1.197,
  worst pull 1.75), "cluster_wave is refused on the channels that never read
  it" (4); `test_fsi.cpp` "the two variants differ on essentially every event"
  (5: survival 0.797094 vs 0.520239, …, 64.45 % vs > 98 %).
- pytest, 22 files: **15 failed** — `test_cli` 1, `test_knob_provenance` 7
  (the `cluster_wave: hulthen -> vmc` cell on all seven specs), `test_module` 3
  (`make_config` default, the C5.2 refusal via a default config, the C5.5b
  refusal), `test_meta_provenance` 2 (the FSI pair), `test_price_defaults` 1
  (the inclusive probe of `vmc` is no longer refused; the `make_config` guard
  of §3.2 would restore the refusal — INFERRED, not run), `test_pipeline.py` 1
  (tag fraction 0.03313 vs the documented 0.0248).
- Passing but with a changed meaning (found by reading, not by failure):
  `test_pipeline.cpp:438-465` (the ⁷Li tag fraction against the Python's
  Hulthén numbers stays inside its 5 %); the knob cells labelled "on
  hulthen" / "0.30 -> 0.33" / "0.0867 -> 0.10", now run on VMC.

**Design B** — A plus the three factory defaults, the bindings' defaults and
`BreakupOptions::source`:

- C++, 15 files: **27 failing cases** — `test_tagged` 14 (23 756 assertions:
  "channel construction against polligen", "the model grid and its tables
  against polligen", "a pure S wave is m-independent …" 22 267 alone, "the
  Hulthen path is untouched by the VMC backend", T25–T27, the S-wave limit,
  …), `test_fsi` 7 (incl. "the prototype FSI/IA table and profile are
  reproduced"), `test_pipeline` 3, `test_breakup` 2, `test_cluster_config` 1.
- pytest, 8 files: **24 failed + 5 errors** (`test_module` 8,
  `test_bench_static` 4, `test_bench_deeps` 2, `test_bench_boeglin` 4,
  `test_bench_bonus` 5 + 5 errors, `test_bench_t5_closures` 1), against 3 for
  design A in the same files.

`test_bench_t4_t6.py -k t6` on design A: 18 passed (`test_t6_measured_row_is_pinned`
included, with the default tagged files now VMC); the expensive T4 row was
skipped.

---

## 5. Inventory — every site (file, line, what, action)

Actions: **FLIP** the default changes · **REWORD** comment/doc text · **EXPLICIT**
pass `Hulthen`/`"hulthen"` explicitly, no number moves · **RETAKE** a pin re-taken with a
dated comment naming D6 · **LABEL** label/text only · **REMEASURE** a published
number re-measured · **AMEND** dated amendment · **KEEP** with the reason ·
**NEW** added. Line numbers are at `ca9b23a`. The historical run records
(`docs/open_items/run_2026-0*`, `vmc_reconciliation.md`, `vmc_overlaps.md`,
`07_cw_sign_investigation.md`) are history and are not inventoried, except
where an amendment is due (§5.7).

### 5.1 Library and run surface

| file | line | what | action |
|---|---|---|---|
| `include/lipolgen/cluster.hpp` | 405-408 | `enum class ClusterWaveSource`; `Hulthen = 0 ///< … THE DEFAULT` | REWORD (Hulthen: "an option, and the factories' argument default"; VmcAV18: "the RUN default since D6"); numbering kept (no value-initialised enum anywhere — grepped) |
| `include/lipolgen/cluster.hpp` | after 408 | — | NEW `inline constexpr ClusterWaveSource CLUSTER_WAVE_DEFAULT = ClusterWaveSource::VmcAV18;` |
| `include/lipolgen/pipeline.hpp` | 766-822 | doc of `cluster_wave`: "`Hulthen` is the default and keeps every published number bit-for-bit"; :799 "Hulthen default bit for bit (T27)" | REWORD (D6; reach unchanged; the inclusive drift now between defaults) |
| `include/lipolgen/pipeline.hpp` | 823 | `ClusterWaveSource cluster_wave = ClusterWaveSource::Hulthen;` | **FLIP** → `= CLUSTER_WAVE_DEFAULT;` |
| `include/lipolgen/pipeline.hpp` | 764-765, 824-833 | `cluster_beta`, `p_d`, `cluster_vmc_mc_sigma` | KEEP (not read at the default; the knob table says so) |
| `include/lipolgen/pipeline.hpp` | 335-337 | `TritonSfChoice::Hulthen` "the DEFAULT" | KEEP (sibling, §1.3) |
| `src/core/pipeline.cpp` | 482-506 | the off-channel refusal: `!= ClusterWaveSource::Hulthen`, "Leave it at Hulthen" | **FLIP** rule to `!= CLUSTER_WAVE_DEFAULT`; REWORD comment and message (§3.2) |
| `src/core/pipeline.cpp` | 1524-1534, 1679-1687 | factories and `bo.source = cfg_.cluster_wave` | KEEP (read the field); REWORD the line 1680 comment to name `source` beside `beta` |
| `src/core/pipeline.cpp` | 2800-2818 | knob row `cluster_wave`, the off-channel reason | REWORD ("refuses any value but the default here"; the inclusive constants are the Hulthén-era product) |
| `src/core/pipeline.cpp` | 2832-2873 | knob rows `cluster_beta`, `p_d` | KEEP (already conditional on `vmc`) |
| `include/lipolgen/tagged.hpp` | 52-75 | "on the HULTHEN DEFAULT … agree to 7.84e-6"; "UNDER VmcAV18 THEY DO DRIFT APART" | REWORD (true of the Hulthén option, false of the D6 default; registry row 10, §11.1) |
| `include/lipolgen/tagged.hpp` | 105-127 | `LI6_CLUSTER_POLARIZATION_VMC` doc | REWORD (a `VmcAV18` run — the default since D6 — reads 0.887076) |
| `include/lipolgen/tagged.hpp` | 147-177 | `DEUTERON_AV18` doc; :171 "the same approximation the Hulthen default makes" | REWORD ("the Hulthén option") |
| `include/lipolgen/tagged.hpp` | 193-217 | `li6_alpha_channel` doc :195 "`source` defaults to `Hulthen`, which keeps every existing number BIT-FOR-BIT"; arg :216 | KEEP the argument default (DD-1); REWORD the doc |
| `include/lipolgen/tagged.hpp` | 218-224, 225-246 | `li7_alpha_channel` (:223), `deuteron_channel` (:228 "the default, bit for bit what it always was", :245-246) | KEEP the arguments; REWORD the docs |
| `include/lipolgen/beams.hpp` | 46-91 | the C5.5 block: :50 "TRUE OF THE DEFAULT AND FALSE OF `ClusterWaveSource::VmcAV18`", :52 "On the Hulthen default the two ARE one wave function", :74 "AUTHOR DECISION 2026-09-04: the default stays 0.811228 …" | REWORD the first two (they invert under D6); the 2026-09-04 decision quoted unchanged + a D6 pointer to §11.1 |
| `include/lipolgen/beams.hpp` | 93-99 | `P_D_LI6` doc :96 "It is also the DEFAULT of the Hulthen tagged channel, which is the bit-compatibility path" | REWORD |
| `include/lipolgen/beams.hpp` | 100-101 | `P_D_DEUTERON` "(AV18-like)" — AV18's is 0.057600 | KEEP value; noted for the author (§1.3, §11.1) |
| `include/lipolgen/breakup.hpp` | 30-60, 207-230 | deuteron branch doc; `BreakupOptions::source = Hulthen` "(the default)" | KEEP the struct default (DD-1); REWORD (the Pipeline overwrites it with `cluster_wave`, whose default is VMC) |
| `include/lipolgen/constants.hpp` | 124-133 | `LI6_B1_RANK2_TRANSFER` = "`TaggedModel(li6_alpha_channel()).tensor_dilution()` evaluated at `P_D_LI6`" | KEEP value (stays literally true under DD-1); REWORD (+ "the Hulthén option's; the D6 default channel's is 0.982576 — registry row 10") |
| `src/core/tagged.cpp`, `src/core/breakup.cpp`, `src/core/b1_nuclear.cpp` | 258-273; 130-140; 289, 933 | comments; explicit `VmcAV18` | KEEP |
| `python/bindings.cpp` | 2730-2735 | enum docstring "`Hulthen` (the default everywhere) keeps every published number bit-for-bit" | REWORD |
| `python/bindings.cpp` | 2824-2851 | `py::arg("source") = ClusterWaveSource::Hulthen` ×3, docstrings | KEEP defaults; REWORD docstrings |
| `python/bindings.cpp` | 4473-4486 | `BreakupOptions` docstrings ("`beta` and `f2` are overwritten …"; "Hulthen (default)") | REWORD (add `source`) |
| `python/bindings.cpp` | 4672-4675 | `PipelineConfig.cluster_wave` "Hulthen by default (bit-compatible)" | REWORD; NEW `m.attr("CLUSTER_WAVE_DEFAULT")` |
| `python/lipolgen/__init__.py` | 98-108 | `CLUSTER_WAVES` comment "\"hulthen\" is the DEFAULT" | REWORD |
| `python/lipolgen/__init__.py` | 455-460, 555-557 | `make_config` docstring; the assignment | REWORD; NEW the explicit-request guard (§3.2) |
| `python/lipolgen/cli.py` | 223-257 | `--cluster-wave` help "'hulthen' (default, bit-compatible)" | REWORD (vmc default; the 2026-09-06 sign numbers kept, labelled) |
| `python/lipolgen/cli.py` | 805 | `DEFAULTS cluster_wave="hulthen"` | **FLIP** → `None` |
| `python/lipolgen/cli.py` | 1276-1277 | comment "--cluster-wave vmc off its channels" | REWORD |
| `python/lipolgen/cli.py` | banner | — | NEW (optional, DD-5.1) |
| `examples/generate_tagged.cpp` | 19-20, 42-56, 60-68, 71-96, 139-146, 164-174 | `--beta`/`--pd` into a default config; banner `beta = %.2f` | NEW `--cluster-wave` (DD-4); citations move (§9.3) |
| `examples/generate_full.cpp` | — | tagged runs on the default config | KEEP (now VMC; no β/P_D options) |
| `P_D_LI6` users outside the tagged factories: `tests/test_b1_nuclear.cpp` (12), `tests/test_beams.cpp` (5), `python/tests/test_b1_model.py` (2), `validation/b1_li6_table.py`, `include/lipolgen/constants.hpp`, `src/core/beams.cpp`, `validation/reference/beams.json`, `validation/reference/tagged.json` | — | the inclusive constants' and the Hulthén option's input, 0.0867 | KEEP — D6 moves no constant (`VMC_P_D_LI6` users likewise) |

### 5.2 C++ tests

| file | line | what | action |
|---|---|---|---|
| `tests/test_pipeline.cpp` | 73-94 | `pipe()` helper (cache key without the wave) | EXPLICIT: a `ClusterWaveSource` parameter in the key; call sites pass it |
| `tests/test_pipeline.cpp` | 398-436 | ⁶Li tag fractions vs the Python 0.0247 / 0.3061 / 0.2545 | EXPLICIT Hulthén (a port comparison) |
| `tests/test_pipeline.cpp` | 438-465 | ⁷Li optics-blind vs the Python 0.9655 … | EXPLICIT Hulthén (passes on VMC only because 5 % is wide) |
| `tests/test_pipeline.cpp` | 466-552 | A_zz^tag closure, `REQUIRE(ndf >= 8)` | EXPLICIT Hulthén; optional VMC twin, its ndf floor declared in the I6 record before it runs (VMC gives 7 at YR — known from §4.6) |
| `tests/test_pipeline.cpp` | 317-396, 554-606, 878-900 | the spectator-marginal closure, the rate identity, the HepMC3 round trip | KEEP — follow the default (identities; pass on VMC: χ²/ndf(k) 1.52 on 13, the identity at the tensor dilution 0.982576) |
| `tests/test_pipeline.cpp` | 1236-1293 | "runs on the VMC cluster waves" | KEEP (explicit VMC) |
| `tests/test_pipeline.cpp` | 1295-1336 | "cluster_wave is refused on the channels that never read it"; :1317 "Hulthen -- the default -- is legal everywhere" | **RETAKE** (the default legal everywhere; `Hulthen` legal on tagged, refused on inclusive/coherent) |
| `tests/test_fsi.cpp` | 384-457 | "the two variants differ on essentially every event" (0.520239 / 0.582899, > 98 %) | EXPLICIT Hulthén (the claim is the Hulthén α–d's; §4.4) |
| `tests/test_fsi.cpp` | 48, 206, 486 | `li6_alpha_channel(0.30, 0.0)`, `li6_alpha_channel()`, `li7_alpha_channel()` | KEEP (factory default); EXPLICIT recommended |
| `tests/test_tagged.cpp` | 74-85 | `li6_model()` / `li7_model()` / `deut_model()` — the rtol-1e-12 polligen port gates | KEEP (factory default); EXPLICIT recommended |
| `tests/test_tagged.cpp` | 341, 514, 834, 866, 888, 1118, 1150 | P_D = 0 limits and Hulthén model cases | KEEP; EXPLICIT recommended |
| `tests/test_tagged.cpp` | 1535-1576 (T25, :1574 "Hulthen stays bit for bit the default"), 1610-1650 (T26, :1617 "On the HULTHEN default"), 1667-1757 (T27, :1712 `bh.options().source == ClusterWaveSource::Hulthen`) | Hulthén sides of C5.4 / C5.5 / C5.5b | KEEP pins; REWORD comments ("the Hulthén option" / "the struct default"); EXPLICIT recommended |
| `tests/test_t2.cpp` | 341 | comment "already covers the Hulthen default" | REWORD |
| `tests/test_breakup.cpp`, `tests/test_cluster_config.cpp` line 551, `tests/test_rc.cpp` line 2200 | — | factory / `BreakupOptions` defaults | KEEP |

### 5.3 Python tests

| file | line | what | action |
|---|---|---|---|
| `python/tests/test_cli.py` | 131-134 | `resolve(["--channel","tagged-alpha"])["cluster_wave"] == "hulthen"` | **RETAKE** (resolves to `None`; the config's value is `VmcAV18`; NEW: an explicit `--cluster-wave` on `--channel inclusive` is refused in one line) |
| `python/tests/test_module.py` | 275-283 | `make_config(channel="tagged-alpha").cluster_wave == Hulthen` | **RETAKE** |
| `python/tests/test_module.py` | 594-600 | the band refused "on Hulthen" through a default config | EXPLICIT `cluster_wave="hulthen"` + RETAKE (the default accepts the band on the Li tags) |
| `python/tests/test_module.py` | 715-728 | C5.5b refusal matrix | **RETAKE** (default legal everywhere; `hulthen` refused off-tagged; an explicit request refused by `make_config`) |
| `python/tests/test_module.py` | 270-271, 610-711 | Hulthén sides through the factory defaults; :630 "on the Hulthen default", :694 `BreakupOptions().source == Hulthen` | KEEP pins; REWORD comments; EXPLICIT recommended |
| `python/tests/test_knob_provenance.py` | 332-342 | `V("cluster_wave","hulthen -> vmc", var=VMC)`; "0 -> 1 on hulthen"; β and p_d cells | **RETAKE**: `vmc -> hulthen` (var `HUL`); the "on hulthen" cells get `base=HUL`; the plain cells run on the default (a refused base on inclusive is skipped by the matrix, lines 546-549) |
| `python/tests/test_meta_provenance.py` | `fsi_runs` fixture; tests at :60-90 | > 98 % differ; survival 0.520239 | EXPLICIT Hulthén |
| `python/tests/test_pipeline.py` | 77-82 | tag fraction vs the documented 0.0248 (the Python's 0.0247) | EXPLICIT Hulthén (own fixture); optional VMC twin vs the documented expectation 0.033810 at the same rel 0.20 (VMC measured 0.03313) |
| `python/tests/test_rc.py` | 860-887 | the tagged band's clamp "bites here" (`rc_clipped_band_events > 0`); comment "-1.79" | KEEP (follows the default; passed on the flipped copy — VMC still clips, 14 of 20 000 at seed 1 with the CLI fill); REWORD the comment (the −1.79 is the Hulthén run's) |
| `python/tests/test_price_defaults.py` | — | li6-incl family | KEEP (fails without the guard, MEASURED; passes with it, INFERRED) |
| `python/tests/test_bench_deeps.py` | 153-226 | `in_verdict == ["hulthen/lc"]`; verdict head; `lc_times_glauber` 422.92334145; `by_ps`; transverse | **RETAKE** (§7.1: `["av18/lc"]`; 582.5163997; [1.1395, 1.1964, 1.1737, 1.0350, 0.6408]; data/PWIA [1.1714, 1.7026, 1.9900, 1.9391, 1.5027], tree [1.9298, 3.4634, 5.4336, 7.3071, 7.1735], same side all True) |
| `python/tests/test_bench_boeglin.py` | 245-254 | REPORT line "hulthen (default)" | **RETAKE** the text only |
| `python/tests/test_bench_bonus.py`, `test_bench_static.py` | :333-334, :354-399; :196 | Hulthén sub-row keys / FAIL `["hulthen"]` | KEEP (§7); REWORD the :196 comment |
| `python/tests/test_bench_spin_deuteron.py` | — | rows 1, 3, 9, 10-history | KEEP — **it pins no tagged-wave number**; the −1.2069 Hulthén reading is pinned only as −0.953283 (k = 0.30, weighted) in `tests/test_tagged.cpp` line 790 on the factory default; −0.5191 appears in no test (only in `validation/vmc_reconcile.py` line 845, `validation/vmc_tag_fractions.py` lines 92-94 output text and the docs) |

### 5.4 Harnesses and validation scripts

| file | line | what | action |
|---|---|---|---|
| `validation/benchmarks/t2_deeps_spectator_tail.py` | 60-70, 94, 149, 151, 237, 279, 361, 449, 462, 533 | rule "of the default wave"; `DEFAULT_WAVE = "hulthen"`; labels; FSI sub-row key hard-coded `"hulthen/lc x GlauberFsiWeight"`; generator text "Hulthen default"; `PLOT_WAVES` | **FLIP** `DEFAULT_WAVE = "av18"`; key and generator from `DEFAULT_WAVE`/`WAVE_LABELS`; dated note in the docstring (§7.1) |
| `validation/benchmarks/t2_bonus_spectator_shape.py` | 140-145, 181, 378, 380, 644, 1103 | rule (d) names "the DEFAULT configuration … the Hulthen deuteron … the tagged channel's shipped wave"; `DEFAULT_CONFIG = "hulthen-40mb"` | LABEL + dated note; `DEFAULT_CONFIG` unchanged (§7.2) |
| `validation/benchmarks/t2_boeglin_nk.py` | 42-45, 56-58 (hashed), 358-359, 367-369, 494, 662, 858 | hashed block "the row's status is the DEFAULT wave's (hulthen)"; `DEFAULT_WAVE`; generator "hulthen (default)" | LABEL + a post-hoc addendum OUTSIDE the block; the block and `DEFAULT_WAVE` unchanged (§7.3) |
| `validation/benchmarks/t2_deuteron_static.py` | 52-53, 144-146, 457-466, 504, 642 | "THE DEFAULT deuteron control"; "(tagged default)" | LABEL + dated note; rule unchanged (§7.4) |
| `validation/benchmarks/t5_sum_rules.py` | 514-518 | `li6_alpha_channel()` "defaults (Hulthen beta 0.3)" | KEEP (true under DD-1); EXPLICIT recommended |
| `validation/price_defaults.py` | 544, 746, 769, 774, 817-835, 885, 1090, 1140 | "shipped default (cluster_wave = hulthen)"; `cluster_wave = vmc` as the priced alternative; notes; figure series | REWORD: the default row named from the knob table; the alternative `cluster_wave = hulthen`; the refusal probes keep working (guard) |
| `validation/vmc_tag_fractions.py` | 181 | "`hulthen` the model default" | LABEL |
| `validation/repin_tagged_from_lipolgen.py` | 74-76 | factory defaults against `tagged.json` | KEEP (DD-1); EXPLICIT recommended |
| `validation/dump_polligen_reference.py`, `validation/vmc_reconcile.py` | 819-821, 1373-1374; — | polligen's own defaults; explicit sources | KEEP |
| `validation/reference/*.json` | — | model-level pins (`tagged.json` is polligen's Hulthén) | KEEP — **no reference JSON moves** under DD-1 |

### 5.5 Documents

| file | line | what | action |
|---|---|---|---|
| `README.md` | 27 | "radial forms: Hulthén S/P/D (default) or … (`--cluster-wave vmc`)" | REWORD |
| `README.md` | 38-41 | "Almost every physics ingredient beyond the toy backends is **opt-in**: the default run is the one every reference JSON is pinned against" | REWORD (shared with D5/D7; under DD-1 the tagged reference JSONs pin the Hulthén factory models, not the default run) |
| `README.md` | 46 | selector row: default `hulthen`, bit for bit | FLIP the row (default `vmc`; "the other value" = the analytic pair; the measured price kept with its window, direction stated) |
| `README.md` | 49 | FSI row "99.50 % … ratio to 68.5" | REMEASURE on the default (64.5 %, 10.14, §4.4) or label it Hulthén |
| `README.md` | 276, 321-326 | CLI example `--cluster-wave vmc` (now redundant); the ⁶Li polarization band | KEEP / REWORD (the default tagged run reads 0.887076; the band stays mandatory) |
| `python/README.md` | 377 | the CLI flag list has no `--cluster-wave` (pre-existing gap) | NEW line |
| `docs/USAGE.md` | 1314-1329 | "Cluster radial forms: Hulthen (default) …"; "stays the default, so every published number is unchanged bit-for-bit" | REWORD |
| `docs/USAGE.md` | 1343-1381 | the d-control table "Hulthén (default, bit for bit)"; :1378 "Hulthén stays the default because …" | REWORD (the rationale kept as dated history, superseded by D6) |
| `docs/USAGE.md` | 1429-1456, 1536-1538 | "The Hulthén default is untouched"; the inclusive drift; "the cluster-wave systematic is quoted as the difference …" | REWORD; KEEP the systematic rule (hulthen − vmc) |
| `docs/USAGE.md` | 1559 | T1 table "S + D Hulthén at `P_D_DEUTERON`" | REWORD (the run's deuteron: AV18 by default) |
| `docs/USAGE.md` | 1679-1697 | FSI two-variant table (pre-sign-fix sample) | REMEASURE (§4.4, both waves) |
| `docs/USAGE.md` | 1994-1998 | tag fractions vs the Python | LABEL Hulthén + add the VMC numbers |
| `docs/USAGE.md` | 2226-2232, 2242 | RC clamp 1.2 % / 2.6 %; "245 of 20 000 = 1.2250 %" | REMEASURE on the default (0.070 % / 3.50 % at seed 1, §4.4) |
| `docs/USAGE.md` | 2758-2761 | "a VMC and a Hulthén tagged file are no longer indistinguishable" | KEEP |
| `docs/PHYSICS_CHANNELS.md` | 49-51 | option lists `--cluster-wave vmc`; the d–p row lacks `--cluster-wave` (pre-existing) | REWORD (`--cluster-wave hulthen`); add it to the d–p row |
| `docs/PHYSICS_CHANNELS.md` | 87 | "`DEUTERON()` … on the Hulthén default, `DEUTERON_AV18()` under `--cluster-wave vmc`" | REWORD |
| `docs/PHYSICS_CHANNELS.md` | 264-272 | "Analytic Hulthén-family radial forms (**default**)"; "ANL VMC … (`--cluster-wave vmc`)"; the dilutions row | REWORD (single-line table rows: no line shift) |
| `docs/PHYSICS_CHANNELS.md` | 302, 347-349 | the T1 deuteron row; "⁷Li α tag under FSI … survival 0.365 … ⁶Li tag's 0.520" | REWORD (Hulthén-labelled; the default ⁶Li survival is 0.797) |
| `docs/PHYSICS_CHANNELS.md` | 660-700 | notes "The deuteron control channel **defaults to Hulthén**", "Hulthén stays the default", "nothing at all on the Hulthén default" | REWORD, keeping the line count (exemption keys at :996, :1002, :1017, :1036, :1041 sit below, §9.3) |
| `docs/PHYSICS_CHANNELS.md` | 855 | the clamp note | KEEP text; the number lives in USAGE |
| `docs/CONVENTIONS.md` | 95-145 | "**on the Hulthén default** the INCLUSIVE effective polarization and the TAGGED S/D interference … are the same wave function"; "the Hulthén default is bit for bit" | REWORD (+ §11.1) |
| `docs/CONVENTIONS.md` | 328-395 | "`ClusterWaveSource::Hulthen` is the DEFAULT everywhere"; :375 "Hulthén **by default and follows `--cluster-wave vmc`**" | REWORD (the RUN default vs the factories' argument default, DD-1) |
| `docs/OPEN_ITEMS_SOLUTIONS.md` | 55, 95, 116, 162, 2330-2345 | item 1 "opt-in shipped (… Hulthén default bit for bit)"; §1 heading "OPT-IN SHIPPED (rows 9, 10, 11 open)" | REWORD (default since D6; row 9 decided; row 10 re-posed) |
| `docs/benchmarking/08_channel_priority.md` | 40-42, 96-97, 202-203, 276-297, 311, 330, 362, 389, §7 log | ingredient G "Hulthén by default … D2 at M2"; dossiers 3.2/3.6/3.7; §5 item 3; §6 trigger | REWORD; G's development item "D2 at M2" → "taken, D6"; re-read C3 for ranks 6–7 under the file's own rule (expected unchanged: the shipped default still moves A_zz^tag ≥ 20 % against an opt-in, and α–d is unanchored, U-2) and log it |
| `docs/benchmarking/BENCHMARK_PLAN.md` | 86, 87, 92, 109 | rows 10/11/16 status texts; §6 T2 "FAIL recorded for the Hulthén pair the tagged channel ships" | LABEL (row 11's headline → AV18, §7.1) |
| `docs/benchmarking/BENCHMARK_PLAN.md` | §8 (258-267) | the three normalisation rules, written in each harness header | KEEP — no rule-3 statement of the four T2 deuteron rows depends on the wave (each compares per-deuteron shapes) |
| `docs/benchmarking/BENCHMARK_PLAN.md` | §7 (239-256) | the change → tier map: `cluster.*`, `tagged.*` → T2, T3; `breakup.*` → T3-cluster, T4-tagged; `cli.py` → T6; `pipeline.*` unmapped | KEEP; `tiers.py --since <pre-I6>` will print `pipeline.*` as UNMAPPED — `make_report.py --all` is the gate anyway |
| `docs/benchmarking/02_data_nucleon_deuteron.md` | 287, 312 | "the tagged channel's DEFAULT deuteron is the Hulthén pair" | dated note |
| `validation/benchmarks/README.md` | 206, 229, 230, 231 | the four rows' status/label texts | LABEL (row 11's numbers) |
| `docs/benchmarking/REPORT.md`, `REPORT.json`, `figures/` | — | generated | regenerate (`make_report.py --all --plots`), never hand-edited |

### 5.6 The citation gate's own files

| file | what | action |
|---|---|---|
| `validation/physics_channels_ranges.json` | fingerprints of cited blocks and pinned use sites | `--fix` (MOVED) / `--record-ranges` after a human re-read (EDITED) |
| `validation/check_physics_channels_links.py` | `ALLOW`, `RANGE_ALLOW` (46 keys on `PHYSICS_CHANNELS.md` lines), `UNNAMED_ALLOW` (18 on `PHYSICS_CHANNELS.md`; line 1103 `("SPIN32_FINITE_GAMMA.md", 151, "examples/generate_tagged.cpp:123")`) | update the keys whose citing line or target line moves (§9.3) |

### 5.7 Registry, decisions, records

| file | line | what | action |
|---|---|---|---|
| `docs/open_items/run_2026-09-03/STATUS.md` | row 9 | "The tagged deuteron CONTROL channel keeps Hulthén at P_D = 0.045" | AMEND: decided by D6 (text §12.3) |
| same | row 10 | "`LI6_CLUSTER_POLARIZATION` stays 0.811228 and the `--cluster-wave vmc` drift is DOCUMENTED, not closed" | AMEND: the drift is now between default runs; the question re-posed; no label moved (§12.3) |
| same | row 11 | "follows … REFUSED where it is never read" | AMEND: the refusal follows the default; the explicit request refused by `make_config`/CLI (§12.3) |
| same | row 23 | ⁶Li effective polarization vs 1/3 | KEEP (cross-reference row 10's note) |
| `docs/open_items/run_2026-09-03/AUTHOR_DECISIONS.md` | §B6 (305-356), §B9 (474-517), §B10 (519-545) | the drafted options | AMEND (dated boxes, §12.3) |
| `docs/open_items/run_2026-09-27/DECISIONS.md` | 22-33 | D2 | a dated pointer "superseded by D6" (§12.2) |
| `docs/open_items/run_2026-10-27/DECISIONS.md` | D6 entry | | the text of §12.1 |
| `docs/open_items/run_2026-10-27/STATUS.md` | I6 row | | the I6 phase record `phase_I6_cluster_wave.md` |

---

## 6. Pin decisions

The rule (this run's ground rules): a number that depends on the default is
**re-taken** with a dated comment naming D6. A test that pins Hulthén physics,
compares against the Hulthén-based sibling prototype, or reproduces a
Hulthén-era published number **passes `Hulthen` explicitly**.

| pin | decision | why |
|---|---|---|
| `test_pipeline.cpp:398-436` (Python 0.0247 / 0.3061 / 0.2545) | EXPLICIT | port comparison with `money_tagged_azz.py`, a Hulthén prototype |
| `test_pipeline.cpp:438-465` (Python 0.9655 …) | EXPLICIT | the same; passes on VMC only by tolerance |
| `test_pipeline.cpp:466-552` (ndf ≥ 8) | EXPLICIT (+ optional VMC twin) | the ndf floor was declared on the Hulthén spectrum; changing it now would be choosing it after the number |
| `test_pipeline.cpp:1295-1336` | RETAKE | the rule itself follows the default |
| `test_fsi.cpp:384-457` | EXPLICIT | "differ on essentially every event" is a Hulthén statement (64.5 % on VMC) |
| `test_cli.py:131-134` | RETAKE | pins the default |
| `test_module.py` lines 275-283 and 715-728 | RETAKE | pin the default / the refusal |
| `test_module.py:594-600` | EXPLICIT + RETAKE | the "refused on Hulthen" half needs Hulthén; the default half is new |
| `test_knob_provenance.py:332-342` | RETAKE | cells are variations from the default |
| `test_meta_provenance.py` FSI pair | EXPLICIT | the D3 claim (0.520239, > 98 %) is the Hulthén α–d's |
| `test_pipeline.py:77-82` | EXPLICIT (+ optional VMC twin) | reproduces a Hulthén-era documented number |
| `test_bench_deeps.py` pins | RETAKE | the rule follows "the default wave" (§7.1) |
| `test_bench_boeglin.py:245-254` | RETAKE (text) | label only |
| factory-level pins (`test_tagged.cpp` gates, T25–T27, P_D = 0 limits; `test_module.py` C5.x; harness `Tree`s; `t5_sum_rules`; `repin_tagged_from_lipolgen.py`) | KEEP, EXPLICIT recommended | bit for bit under DD-1; explicit makes them immune to a later factory flip and makes "hulthen" labels true by construction |

---

## 7. Benchmark rows and their declared rules

No verdict moves under D6 whichever reading is taken. That was established
before choosing a reading, from numbers already on record (rows 10/16/§5.2)
and from §4.3 (row 11).

### 7.1 `t2_deeps_spectator_tail` (plan row 11) — the rule follows the default

READ, lines 67-70: "TOLERANCE, declared BEFORE any chi2 of this row was computed
(2026-09-27): p(chi2, ndf) >= 0.01 on the COMBINED PWIA-window shape chi2 of
the default wave with the LC prescription." The identification is descriptive,
in another paragraph (lines 60-62): "Waves: the shipped default,
`deuteron_channel()` = the analytic Hulthen pair at P_D = 0.045 (the
tagged-channel default) -- THE HEADLINE VERDICT".
**Action:** the verdict row becomes `av18/lc` (579.12/187, p = 1.1 × 10⁻⁴¹,
FAIL), Hulthén a recorded sub-row (627.85/187), status unchanged; numbers §4.3.

### 7.2 `t2_bonus_spectator_shape` (plan row 10) — the rule names the Hulthén

READ, lines 140-145: "(d) TOLERANCE … THE ROW is PASS iff the DEFAULT configuration
passes at both beams: the Hulthen deuteron (deuteron_channel() default, P_D =
0.045, the tagged channel's shipped wave) at sigma_XN = 40 mb." The rule
sentence itself defines the default configuration *as* the Hulthén. No hashed
block. **Action:** keep `DEFAULT_CONFIG = "hulthen-40mb"`. Add a dated note
under (d): the shipped wave is AV18 since D6, and AV18 at 40 mb, already
recorded under (e), fails too ("AV18 40 mb: 3/12 and 0/16", lines 206-207; χ²
against the weight's shape 221.6/84 against 226.3/84). The row stays FAIL.
Re-pointing it to the shipped default belongs to the author's pending
re-declaration of this row (`../run_2026-10-26/SUMMARY.md`, "For the author"
item 6).

### 7.3 `t2_boeglin_nk` (plan row 16) — the hashed block names the Hulthén

READ, lines 56-58, inside the hashed block (`DESIGN_BLOCK_SHA256 = "437eab47…"`,
lines 358-359): "THE VERDICT, one per wave; the row's status is the DEFAULT wave's
(hulthen), av18 and cdbonn are sub-rows by the same rule." The block is never
edited. `design_block_sha256()` checks its text and the numeric rules, not
`DEFAULT_WAVE`, so a code edit would silently contradict the block's
parenthetical. **Action:** keep `DEFAULT_WAVE = "hulthen"`, relabel the
generator ("hulthen (the design block's verdict wave)"), and add a post-hoc
addendum OUTSIDE the block: since D6 the shipped wave is AV18, which passes by
the same rule (11.54/21, p = 0.951). PASS either way.

### 7.4 `t2_deuteron_static` (T2 §5.2) — the rule covers every shipped wave

READ, lines 83-85: "|eta_model − 0.0256| <= 2 sigma … for EVERY shipped wave
function; the row passes only if all three do." The Hulthén still ships (as
an option). **Action:** nothing moves but labels ("THE DEFAULT deuteron control"
→ "the Hulthén option"; the plot's "(tagged default)"). FAIL stays (+17.2σ).
Re-declaring it to "the default wave" would be a post-hoc re-declaration,
which is for the author, not D6.

### 7.5 The other rows

`t6_hepmc3_convention` writes default tagged files (now VMC):
`test_t6_measured_row_is_pinned` passed on the flipped copy, so no change is
expected. `t5_sum_rules`' FSI unitarity row uses the factory default:
unchanged under DD-1. Every other row reads no tagged wave. The REPORT summary
(statuses) is unchanged by D6. The plan-row status texts
(`BENCHMARK_PLAN.md` lines 86, 87, 92 and 109) and the README rows
(`validation/benchmarks/README.md` line 206, lines 229-231) take the labels of
§7.1–7.4. Only row 11 takes new numbers.

---

## 8. Published numbers that move, and what regenerates them

| number | where | regenerated by |
|---|---|---|
| the Deeps headline 627.85 → 579.12 and its sub-rows | REPORT, `BENCHMARK_PLAN.md` line 87, `validation/benchmarks/README.md` line 230, `08` §0/§3.6 | `make_report.py --all --plots` (REPORT); the prose by hand |
| FSI "two variants" table and README's "99.50 % / 68.5" | `docs/USAGE.md` lines 1679-1697, `README.md` line 49, `08` ingredient I "99.5 %" | by hand (recipe: `--channel tagged-6Li-alpha --events 20000 --seed 1234`, three `--fsi` values) |
| RC clamp 1.2 % / 2.6 %; 245 / 1.2250 % | `docs/USAGE.md` lines 2226-2232, line 2242 | by hand (`--rc tensor-band`, 20 000 events; state the seed and fill) |
| tag fractions vs the Python | `docs/USAGE.md` lines 1994-1998 | `tests/test_pipeline.cpp` MESSAGEs (400 000 events) |
| the `--cluster-wave` selector row (direction of the price) | `README.md` line 46 | `validation/price_defaults.py` (the li6-tag/li7-tag families) |
| the Hulthén-vs-VMC tables | `docs/USAGE.md` lines 1460-1530, `docs/open_items/vmc_reconciliation.md` | `validation/vmc_tag_fractions.py`, `validation/vmc_reconcile.py` — numbers do not move, only "(default)" labels |
| the ⁷Li FSI survival 0.365 vs ⁶Li 0.520 | `docs/PHYSICS_CHANNELS.md` line 349 | factory-level (Hulthén) — label |
| the price table (phase A) | `../run_2026-10-26/phase_A_baseline.md` | history — not regenerated; `price_defaults.py` would now print the inverse ratios (×1/0.430 …) with `hulthen` as the alternative |
| throughput tables | `docs/USAGE.md` §8 | machine-specific; not regenerated by D6 |

---

## 9. Ordered implementation steps (phase I6) and the gates after each

### 9.0 Order relative to the other implementations

The board orders I4, I5, I6, I7 (`STATUS.md`). D6 shares files with all three:
`include/lipolgen/pipeline.hpp`, `src/core/pipeline.cpp`, `python/lipolgen/cli.py`,
`python/lipolgen/__init__.py`, `python/bindings.cpp`, `README.md`,
`docs/USAGE.md`, `docs/PHYSICS_CHANNELS.md`, `docs/CONVENTIONS.md`,
`validation/price_defaults.py`, `python/tests/test_knob_provenance.py`,
`python/tests/test_cli.py`, `python/tests/test_module.py`, the citation
sidecar and exemption keys, and the REPORT. **Implement sequentially, one
decision per commit, each on top of the previous one's gates.** No physics pin
depends on both D6 and D7: the tagged σ moves with D7 only (MEASURED, σ
bit-identical across the waves, §4.2) and A_zz^tag with D6 only. The overlap
is in labels, the knob matrix and the price script. Before Step 1, save the
committed `docs/benchmarking/REPORT.json` (the previous phase's) to the
scratch directory as the `--compare` baseline.

### 9.1 Step 1 — make every Hulthén-meaning reliance explicit (no behaviour change)

`pipe()` gains the wave. These pass `Hulthen` explicitly:
`test_pipeline.cpp:398-552`, `test_fsi.cpp:384-457`, `test_pipeline.py:77-82`,
`test_meta_provenance.py`'s `fsi_runs`, the `test_module.py:594-600` refusal
half, and the knob cells "on hulthen" (`base=HUL`). Recommended in the same
step, all bit for bit: the factory-level Hulthén pins of §6's last row.

**Gate 1:** `cmake --build build -j2 && ./build/lipolgen_tests` (full) · `python3
-m pytest python/tests -q` · `LIPOLGEN_BENCH=1 python3 -m pytest python/tests -k
bench -q` · `python3 validation/benchmarks/make_report.py --all --compare <saved
REPORT.json>` → **0 changes** (bit for bit) · the citation gates (§9.3; test
files are cited, e.g. `tests/test_tagged.cpp` lines 396-413 and 723-726) ·
`python3 validation/check_spdx_headers.py`.

### 9.2 Step 2 — the flip

`cluster.hpp` (`CLUSTER_WAVE_DEFAULT`, enum docs), `pipeline.hpp` line 823 + doc,
`pipeline.cpp` lines 482-506 (rule + message) and 2800-2818 (knob reason), bindings
(attr + docstrings), `__init__.py` (guard + docs), `cli.py` (DEFAULTS `None`,
help, comment), `price_defaults.py` (labels/alternatives), the example
(DD-4), optional DD-5.

**Gate 2a — the predicted break list, and nothing else:** full C++ → exactly
`test_pipeline.cpp` "cluster_wave is refused …" (4 assertions). pytest →
exactly `test_cli::test_cluster_wave_flag_reaches_the_config`,
`test_module::test_make_config_takes_cluster_wave_by_name`,
`test_module::test_c5_5b_cluster_wave_is_refused_where_it_is_never_read`, and the
knob matrix's `cluster_wave` cell on the three tagged specs. The four
inclusive/coherent cells and `test_price_defaults` should pass through the
guard (INFERRED; §4.6 measured them failing without it).
Anything else is a reliance Step 1 missed: make it explicit or re-take it, and
record which. **Then RETAKE** those (dated "D6 (2026-10-27)").
**Gate 2b:** full C++, pytest, bench, the citation gates with `--fix` (§9.3),
SPDX, and the regression check of §10.1.

### 9.3 The citation gates and line shifts

`validation/check_physics_channels_links.py` (strict on
`docs/PHYSICS_CHANNELS.md`, `docs/theory/SPIN32_FINITE_GAMMA.md`,
`docs/PYTHIA_BRIDGE.md`; `--records` adds the relaxed R6 over
`docs/open_items/`). READ, its docstring: `--fix` "rewrites a broken citation:
a named point onto the nearest DECLARATION of the name … a range or a pinned
use site onto its new home when the fingerprint locates it". An EDITED block
("found NOWHERE") needs a human: "re-read the row, confirm the citation still
says something true, then re-record with `--record-ranges`", which refuses a
block its citing sentence does not evidence. Unnamed point citations (rule B)
are not fixed automatically. D6 edits files the strict documents cite
(`PHYSICS_CHANNELS.md` cites `pipeline.hpp` 98×, `cli.py` 71×,
`pipeline.cpp` 43×, `tagged.cpp` 40×, `cluster.hpp` 30×, `tagged.hpp` 26×,
`breakup.hpp` 20×, `beams.hpp` 18×, `bindings.cpp` 17×, `constants.hpp` 16×,
`cluster.cpp` 15×, `breakup.cpp` 15×, `beams.cpp` 14×, `__init__.py` 13×). Some
cited blocks are edited in content: `include/lipolgen/tagged.hpp` lines 54-62 ("the
'same wave function seen in two experiments' rule", cited at `PHYSICS_CHANNELS.md`
line 270), `include/lipolgen/beams.hpp` lines 99-100, `include/lipolgen/cluster.hpp`
lines 9-22 and the named points at its lines 406 and 407 (`Hulthen`, `VmcAV18`;
cited at lines 265-266), `include/lipolgen/pipeline.hpp` line 823 `cluster_wave`
(cited at line 270), `python/lipolgen/cli.py` line 221 (an unnamed point, cited
at line 302), `examples/generate_tagged.cpp` line 123 (an unnamed point, cited
at `SPIN32_FINITE_GAMMA.md` line 151; its `UNNAMED_ALLOW` key is at the gate's
line 1103) and its lines 119-124 (a range, cited at `SPIN32_FINITE_GAMMA.md`
line 208). Procedure after each code step:

1. `python3 validation/check_physics_channels_links.py` → read the report.
2. `python3 validation/check_physics_channels_links.py --fix` for MOVED named
   points and ranges; re-point unnamed citations by hand and update their
   `UNNAMED_ALLOW` key (the target line is part of the key).
3. For every EDITED block: re-read the row, then `--record-ranges`. A
   `RANGE_ALLOW` exemption is keyed (document, citing line, range), so a
   document edit that moves a citing line retires it loudly and the key must
   move with it. `PHYSICS_CHANNELS.md` table rows are single lines; keep the
   line count of the wrapped notes at lines 652-705 to leave the keys at lines
   996, 1002, 1017, 1036 and 1041 where they are.
4. `--records`: a dated record citing a line that moved is history. Annotate
   it `(as of <pre-I6 commit>)`, the honest repair per the gate's R6 note,
   rather than re-point it.
5. `python3 -m pytest python/tests/test_doc_link_gate.py -q` (strict, the
   extra documents, live fingerprints and exemptions, records ≥ 640 cited).

### 9.4 Step 3 — the harness rows

§7.1–7.4 actions; RETAKE `test_bench_deeps.py` and `test_bench_boeglin.py`.
**Gate 3:** `LIPOLGEN_BENCH=1 … -k bench` · `make_report.py --all --compare
<saved REPORT.json>`. Expected: `t2_deeps_spectator_tail` headline 627.85 →
579.12 and its sub-row set (`av18/lc` in the verdict, `av18/lc x
GlauberFsiWeight` 582.52 replacing `hulthen/…` 422.92, transverse values,
config label); label-only text on rows 10, 16 and T2 §5.2; **no status change
on any row**. Then `make_report.py --all --plots` and `make_report.py --check`.

### 9.5 Step 4 — documents and re-measured numbers

§5.5 and §8, with every re-measured number given its window. **Gate 4:** the
citation gates (§9.3), `test_doc_link_gate.py`, SPDX.

### 9.6 Step 5 — registry, decisions, records

§5.7 and §12. **Gate 5 (final):** full C++, full pytest, bench, all citation
modes, SPDX, `make_report.py --check`; `git diff --stat ca9b23a.. --
validation/reference` empty.

---

## 10. Verification plan (for phase V)

1. **The key regression check (code/regression lenses).** A pre-I6 default
   run and a post-I6 `--cluster-wave hulthen` run must be byte-identical in
   every column and every `meta` value except the `cluster_wave` knob row's
   `at_default`. Run it on each tagged channel (⁶Li α `tensor-thirds`, ⁷Li α
   `helicity-flip`, d–p), 20 000 events, two seeds, `_fingerprint` of
   `price_defaults.py`. And a pre-I6 `--cluster-wave vmc` run must be
   byte-identical to a post-I6 default run. Inclusive and coherent default
   runs must be byte-identical before and after, apart from the
   `cluster_wave` knob row's value (`hulthen` → `vmc`, refused, at default).
2. **Physics lens.** A default ⁶Li α run at N = 200 000, seed 20261026,
   reproduces §4.2's vmc column (A_zz^tag −0.5191289 / −0.1992630, RP 0.03422
   ± 0.00041, dilutions 0.9709661 / 0.9825759). The channel label reads "VMC
   AV18+UX", `dis_target.eff_pol_p` = 0.9136001620, and the T1 breakup's
   `deuteron_model()` P_D = 0.0575998919874. d–p: §4.2's vmc column.
   `LI6_CLUSTER_POLARIZATION`, `LI6_B1_RANK2_TRANSFER` and every
   `validation/reference/*.json` are sha256-unchanged.
3. **Silent-knob lens.** On the default, `--p-d 0.05` and `--cluster-beta 0.33`
   appear in the banner's loud "NOT READ, and set away from the default" block
   with the reason. `--cluster-wave vmc --channel inclusive` and
   `--cluster-wave hulthen --channel inclusive` are refused in one line;
   `lg.make_config(channel="inclusive", cluster_wave="vmc")` raises.
   `examples/generate_tagged --pd 0.05` refuses unless `--cluster-wave
   hulthen`.
4. **Honesty lens.** No document still says the Hulthén is the default or that
   a default run is "bit for bit" the pre-D6 one. Every re-measured number
   carries its window. The registry rows carry dated amendments, and row 10's
   question is posed, not answered.
5. **Benchmark lens.** `--compare` lists exactly §9.4's expected changes; rows
   10 and 16 carry their dated notes outside any hashed block;
   `DESIGN_BLOCK_SHA256` of `t2_boeglin_nk.py` and `t4_djangoh_rad_noRad.py`
   unchanged.

---

## 11. Risks

### 11.1 First: the ⁶Li effective polarization — a question for the author, not a decision

After D6 the two DEFAULT ⁶Li runs no longer share one wave function. MEASURED
here, on `ca9b23a`:

| quantity | default tagged ⁶Li α (D6) | default inclusive ⁶Li | ratio |
|---|---|---|---|
| α–d vector dilution | 0.970966 (`VMC_P_D_LI6` = 0.019355) | 0.869950 = 1 − 1.5 `P_D_LI6` | ×1.11613 |
| embedded-deuteron vector dilution | 0.913600 (AV18, P_D 0.057600) | 0.9325 (`P_D_DEUTERON` 0.045) | ×0.97973 |
| whole-nucleus polarization (the product) | **0.887076** | **0.811228** (`LI6_CLUSTER_POLARIZATION`) | **×1.09350** |
| α–d rank-2 transfer | **0.982576** | **0.921947** (`LI6_B1_RANK2_TRANSFER`) | **×1.06576** |
| against the six-body VMC 0.848 | +4.61 % | −4.34 % | |

What the tree says about re-deriving the constant (READ):

- Registry row 10 (`../run_2026-09-03/STATUS.md`): "**`LI6_CLUSTER_POLARIZATION`
  stays 0.811228 and the `--cluster-wave vmc` drift is DOCUMENTED, not
  closed.** … substitution is the wrong fix: the self-consistent reading
  0.905427 sits **6.8 % above** the ab-initio six-body VMC **0.848** where the
  shipped 0.811228 is 4.3 % below it, and with the AV18 deuteron too 0.887076
  is still 4.6 % above. The cluster **product** carries the error, not the
  choice of P_D — `P_D_LI6` = 0.0867 was never an α–d D-state probability but
  a fitted effective depolarisation. **Cost:** a tagged row and an inclusive
  row of one programme differ by 11.61 %, and **no inclusive ⁶Li polarization
  may be quoted without the band 0.81 … 0.91**."
- `include/lipolgen/beams.hpp` lines 74-79: "AUTHOR DECISION 2026-09-04: the default
  stays 0.811228 (it is the closer of the two to 0.848 and it is pinned bit for
  bit in validation/reference/), the drift under `--cluster-wave vmc` is
  DOCUMENTED at 11.61 % rather than "fixed", and nothing quotes an inclusive
  6Li polarization without the band 0.81 .. 0.91."
- `include/lipolgen/tagged.hpp` lines 67-71: "It is NOT closed by substitution --
  `LI6_CLUSTER_POLARIZATION_VMC` below overshoots the ab initio
  `LI6_POLARIZATION_VMC_SIX_BODY` = 0.848 by 6.8 % where the shipped constant
  undershoots it by 4.3 % -- so it is documented, banded and left alone", and
  lines 115-118: "It is a DIAGNOSTIC, not a default … must never be substituted
  for the shipped one on the strength of being more self-consistent."
- `AUTHOR_DECISIONS.md` §B6(b): options (i) 0.811228, (ii) 0.905427, (iii)
  0.887076, (iv) 0.848, (v) an opt-in override; §B6(d): "(ii)–(iv): every
  inclusive ⁶Li g₁ / A_∥ scales by new/0.811228 (×1.1161, ×1.0935, ×1.0453);
  `validation/reference/beams.json` and `xsec.json` (both carry 0.811228) move
  at rtol 1e−12".

**The consistent choice** (INFERRED): one wave function for the default tagged
and the default inclusive ⁶Li is option (iii) together with the rank-2 analogue.
That is `LI6_CLUSTER_POLARIZATION = li6_cluster_polarization(VMC_P_D_LI6,
deuteron_av18_p_d())` = 0.887076 — the same object a default tagged run reads
end to end — and `LI6_B1_RANK2_TRANSFER` = the VMC channel's `tensor_dilution()`
= 0.982576. Its price: every inclusive ⁶Li g₁ / A_∥ ×1.09350; b₁(⁶Li) and the
shipped inclusive A_zz ×1.06576 (A_zz is linear in b₁ through the transfer,
registry row 6); `beams.json`, `xsec.json` and `b1_default_li6.json` re-pinned
(a moved `validation/reference` file, which this run's rules reserve to the
author); and the inclusive number crosses the six-body anchor, from 4.34 %
below 0.848 to 4.61 % above it. That 0.27-point difference in distance is the
whole of the 2026-09-04 reason for keeping 0.811228.

**The question:** with VMC the tagged default, which reading should the
inclusive ⁶Li constants carry — (i) keep 0.811228 / 0.921947, the drift now
between the two default runs and the band mandatory; (iii) re-derive to
0.887076 / 0.982576, consistent; or (iv) 0.848 for the vector part (the rank-2
transfer then has no ab-initio anchor)? A sub-question: the inclusive deuteron
species' `DEUTERON_VECTOR_POLARIZATION` 0.9325 (`P_D_DEUTERON` 0.045, "AV18-like"
at `beams.hpp` line 100, AV18's being 0.057600). Under (i) the I6 implementation
changes no constant and states the drift in the D6 text; that is what §12 does.

### 11.2 Other risks

1. **The off-channel refusal** (§3.2): without the rule change every default
   inclusive and coherent run throws. Without the `make_config` guard,
   registry row 11's protection is lost for explicit requests.
2. **Silent knobs on the default**: `--p-d`, `--cluster-beta` (except ⁷Li T1)
   and the example's `--beta`/`--pd` become unread by default. The knob table
   flags the first two loudly; the example does not (DD-4).
3. **A runtime data dependency** where there was none: every default tagged
   run opens `data/vmc` (DD-5.2, DD-5.3). The error is loud but does not name
   the escape. The conventional `cmake --install` installs no `data/`. The
   wheel vendors the 20 tables (READ, `docs/PACKAGING.md` line 279).
4. **The FSI weight clips on the default** (`--fsi`, opt-in): 0.70 % / 5.25 % of
   the weight table at 40 mb (cluster / nucleon), w up to 37.9 and the 50
   ceiling (§4.4). Reported in `meta["fsi_clipped_grid_fraction"]`, not on the
   banner.
5. **The default loses its external port gate**: `tagged.json` stays
   polligen's Hulthén. The new default's checks are in-repo: the
   reconciliation tables (`validation/vmc_reconcile.py`) and the CW gate on
   AV18. The sibling's VMC α–d agrees to ~5 × 10⁻⁹ (§4.5). A `li6_alpha_vmc`
   block dumped from polligen is a candidate for the author (it would cover
   the α–d motion, not the embedded deuteron the sibling keeps at 0.9325).
6. **No wave-function band on the default row**: `fdeut.av18` prints no MC
   errors and `cluster_vmc_mc_sigma` is 0 by default (a 0.02 % effect). The
   quoted cluster-wave systematic stays hulthen − vmc (READ, `docs/USAGE.md` lines 1536-1538),
   which is the reason to keep the Hulthén option.
7. **Tables end at 5 fm⁻¹ = 0.9866 GeV**: no spectator is drawn past it on the
   default (READ, `docs/CONVENTIONS.md` lines 374-375); the Hulthén tail reached the
   1.2 GeV grid.
8. **Concurrent decisions**: I4/I5/I7 touch the same files (§9.0). Parallel
   implementation would collide in the citation sidecar and exemption keys.
9. **Two "defaults" in the API** (DD-1): `lg.li6_alpha_channel()` (Hulthén)
   and `lg.run(channel="tagged-6Li-alpha")` (VMC) differ. Mitigated by the
   docs and the named constant; the measured alternative is worse (§4.6).
10. **The sibling diverges**: its prototype's defaults are Hulthén, and its
    VMC path keeps the scenario embedded deuteron (2.069 % on vector
    observables, C5.5b). Any default-run tagged comparison with it now needs
    the wave named. The letter quotes no spectator-tagged number (READ,
    `PolarizedLithiumSim/paper/main.tex`: "tagged" there is the coherent
    recoil), so the letter is unaffected.

---

## 12. Texts

### 12.1 `docs/open_items/run_2026-10-27/DECISIONS.md` — D6 (supersedes D2)

> **D6 — the VMC AV18+UX cluster waves (`--cluster-wave vmc`) become the
> shipped default; the Hulthén forms stay an option (`--cluster-wave
> hulthen`). Supersedes D2 of 2026-09-27.** *(decision 3: "if that's the well
> accepted wave default, then do it, otherwise if there is a better model
> favored by the measured D/S ratio we can use it")*
>
> *Evidence.* The deuteron's asymptotic D/S ratio is measured, η_d = 0.0256 ±
> 0.0004 (Rodning & Knutson, PRC 41 (1990) 898, Eq. (10), p. 906). It rejects
> the Hulthén pair the tagged channels shipped (η = 0.032478, +17.2σ;
> `t2_deuteron_static`, |pull| ≤ 2 declared before the verdict). It accepts
> both realistic potentials in the tree: AV18 (0.025045, −1.39σ; the AV18
> paper prints 0.0250 against the same 0.0256(4), PRC 51 (1995) 38, Table X)
> and CD-Bonn (0.025571, −0.07σ). The measurement favours CD-Bonn numerically
> but does not reject AV18 at the declared tolerance. CD-Bonn exists in the
> tree only as a deuteron: it has no ⁶Li α–d or ⁷Li α–t table. The ANL VMC
> tables share the AV18 two-nucleon force — "a Hamiltonian containing the
> Argonne v18 two-nucleon and Urbana X three-nucleon potentials" (Wiringa
> et al., PRC 89 (2014) 024305, abstract), with the 2004 AV18+UIX overlap
> for the S–D sign — across the ⁶Li α–d and ⁷Li α–t overlaps, and through
> `fdeut.av18` across the d–p control and the embedded deuteron. CD-Bonn for
> the control alone would put two wave-function families in one programme,
> which the C5.4/C5.5b rule forbids. So the default is AV18/VMC, and CD-Bonn
> is the recorded alternative for the deuteron.
>
> *What changes.* `PipelineConfig::cluster_wave` — hence `lipolgen-run`,
> `lg.run`, `lg.make_config` — defaults to `VmcAV18` (`CLUSTER_WAVE_DEFAULT`,
> `cluster.hpp`) on the three tagged channels. That covers the ⁶Li α–d
> relative wave (P_D 0.0867 → 0.019355), its embedded deuteron (0.9325 →
> 0.913600) and its T1 breakup's deuteron; the ⁷Li α–t relative wave; and the
> d–p control's p–n wave (P_D 0.045 → 0.057600). Measured at config 1, YR
> high acceptance, N = 200 000, seed 20261026: ⁶Li A_zz^tag −1.2069 → −0.5191
> at k = 0.198 GeV (×0.430) and −0.9533 → −0.1993 at 0.301 GeV (×0.209); RP tag
> fraction 0.02466 → 0.03422; ⁷Li ⟨k⟩ 0.2885 → 0.1861 GeV; d–p A_zz^tag
> −0.8763 → −0.9625 and −0.5478 → −0.6681, tag fraction 0.02396 → 0.02107;
> every cross section unchanged. `validate()` refuses a value moved from the
> default off the tagged channels, and `make_config` and the CLI refuse an
> explicit `--cluster-wave` there.
>
> *What does not change.* The channel factories' and `BreakupOptions`' own
> argument default stays the analytic family: their β and P_D are its
> parameters, and the polligen port gate pins them. Every
> `validation/reference/*.json` stays as it is. The inclusive and coherent
> channels are unchanged. `LI6_CLUSTER_POLARIZATION` = 0.811228 and
> `LI6_B1_RANK2_TRANSFER` = 0.921947 are unchanged: the default tagged and
> inclusive ⁶Li now differ by ×1.0935 (whole-nucleus, 0.887076 against
> 0.811228) and ×1.0658 (rank-2), a question put to the author on registry
> row 10. The ⁷Li triton, `--triton-sf hulthen` and `DeuteronConvolutionB1`
> (AV18 already) are unchanged. No benchmark verdict moves: Deeps' headline
> follows its "default wave" (627.85 → 579.12 on 187 dof, FAIL either way);
> BONuS and Boeglin keep the Hulthén their declared rules name (FAIL / PASS
> on either wave); `t2_deuteron_static` stays FAIL, its rule covering every
> shipped wave. The Hulthén option is the partner of the quoted cluster-wave
> systematic (hulthen − vmc). D2's "physics numbers use AV18" is now the
> default, and D2's "revisit at M2" is done before M2.

### 12.2 `docs/open_items/run_2026-09-27/DECISIONS.md` — under D2

> *Superseded 2026-10-27 by D6 (`../run_2026-10-27/DECISIONS.md`): the VMC
> AV18+UX cluster waves are the shipped default and the Hulthén forms the
> `--cluster-wave hulthen` option. The `t2_deuteron_static` FAIL stays
> recorded — its rule covers every shipped wave function, and the Hulthén
> still ships.*

### 12.3 Registry and `AUTHOR_DECISIONS.md` amendments (dated boxes; no label moved)

- **Row 9 / §B9:** "*Decided 2026-10-27 (D6): option (ii) — the control
  follows the shipped default, now `VmcAV18`, the AV18 deuteron of
  `fdeut.av18` at P_D = 0.057600; the Hulthén pair at 0.045 is `--cluster-wave
  hulthen`. Leg (1): the CW gate already ran on AV18. Leg (2): `tagged.json`
  stays polligen's Hulthén port gate, untouched, since the factories keep the
  analytic default. Leg (3): `fdeut.av18` prints no MC errors; the default
  control carries no band of its own, and the quoted systematic is hulthen −
  vmc. Leg (4): the analytic family stays the option. Measured cost
  (`../run_2026-10-27/plan_D6_cluster_wave_default.md` §4.2): vector/tensor
  dilutions ×0.97973 / ×0.98818, A_zz^tag ×1.098 / ×1.219 at k = 0.198 / 0.301
  GeV, tag fraction ×0.879. Closed.*"
- **Row 10 / §B6:** "*2026-10-27 (D6): the drift this row documents is now
  between the two DEFAULT runs — a default tagged ⁶Li α run reads the VMC
  α–d overlap and the AV18 deuteron (0.887076 whole-nucleus, 0.982576
  rank-2), a default inclusive run 0.811228 and 0.921947. No label moved: the
  2026-09-04 decision stands until the author answers the question of
  `../run_2026-10-27/plan_D6_cluster_wave_default.md` §11.1 (keep / re-derive
  at 0.887076 and 0.982576, the consistent choice / the six-body 0.848). The
  band 0.81 … 0.91 stays mandatory.*"
- **Row 11 / §B10:** "*Amended 2026-10-27 (D6): with the default now
  `VmcAV18`, `validate()` refuses a `cluster_wave` moved from the default
  (`Hulthen`) on `inclusive`/`coherent`. An explicit `--cluster-wave` of
  either value off the tagged channels is refused by `make_config` and the
  CLI, which can tell an explicit request from the default. A C++ caller who
  sets the default value explicitly is indistinguishable from one who did not,
  as for every knob. 'Follow' (the embedded deuteron) is unchanged.*"
- **Row 23:** unchanged; cross-reference row 10's note.

---

## 13. What this design did not verify

- The whole C++ and pytest suites were not run (task rule): 10 + 15 C++ files
  and 22 + 8 pytest files under the two designs, chosen by grep (§4.6). The
  unrun files construct no `PipelineConfig`, call no factory and run no tagged
  channel (grepped). `make_report.py --all` was not run (task rule): the
  harness numbers of §4.3 come from importing the committed harnesses
  read-only. The expensive T4 row was skipped.
- The "well accepted microscopic default" premise is the author's. Two papers
  were quoted (§2); usage was not surveyed.
- The published RC-clamp number 245 / 1.2250 % (seed 1) was not reproduced
  (237 with the CLI-default fill); its recipe is not stated.
- The figures (`--plots`) were not drawn here.
- Whether `--compare` flags a label-only text change was not checked. Its
  documentation speaks of statuses and numbers only.
- This record itself adds nothing to the records citation gate: it carries
  no parsed `path:line` citation (`check_physics_channels_links.py
  --records-only` on 2026-10-27: 0 broken, 0 skipped from this file).

## 14. Files and commands (scratch, outside the repository)

Scratch directory: `/tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/`

- `d6scripts/load_cost.py` (§4.1), `d6scripts/tagged_prices.py` (§4.2:
  `python3 tagged_prices.py <repo> 200000 20261026`), `d6scripts/fsi_rc_prices.py`
  and `d6scripts/fsi_ceiling.py` (§4.4), `d6scripts/deeps_av18.py` (§4.3);
  outputs `runs/`.
- The flip copy `d6flip/` (left in the design-B state), logs
  `cxx_flip_*.log`, `py_flip_batch1.log`, `py_flip_batch2.log` (design A),
  `cxx_flip2_*.log`, `cxx_flip2_summary.txt`, `py_flip2_batch.log` (design B).
- Papers read (none fetched by this task): Rodning–Knutson
  (`refs/rodning_knutson1990_PhysRevC.41.898.pdf`, sha256
  `4dbfb75f5d67e437f9115a8427e6db0508b0f3984e59a2d1eebee79dee0ff161`);
  `PolarizedLithiumSim/refs/nucl-th_9408016.pdf` (AV18, sha256
  `44f551a2adb5899994ca73f51610bc803c289aa67227c5971118a130510996f5`);
  `PolarizedLithiumSim/refs/1309.3794.pdf` (Wiringa 2014, sha256
  `00fd24f199d2c0a899ad6cbbfd1501e73d14c7c477e3617c224ced4a6f28e20b`); text
  via `pdftotext -layout` into `txt/`.

## Adversarial check (2026-10-27)

The adversarial check of this design, run on 2026-10-27 against `ca9b23a`
with nothing in the repository edited, reproduced in the checker's own words
(nothing trimmed, no number changed). Its file and line references are the
checker's, at `ca9b23a`, and are written as plain text, not as citations the
records gate checks.

### Verdict

The direction is sound, but the plan is not ready to implement as written. What holds: I checked the core design against the code and by measurement. DD-1 (flip PipelineConfig::cluster_wave only; the factories and BreakupOptions keep Hulthen) and DD-2 (validate() refuses a value moved from the default; make_config refuses an explicit request off-channel) hold up. The design-A flip prediction matches the scratch logs: 4 C++ cases and 15 pytest failures. Every number I re-ran reproduces exactly: A_zz^tag (6Li: -1.2069113/-0.9532827 -> -0.5191289/-0.1992630; d-p: -0.8763449/-0.5478468 -> -0.9624645/-0.6680525), the dilutions, the FSI two-variant table on both waves (10403.75/11584.83, 99.425 % -> 15989.35/16428.71, 64.515 %, max ratio 10.14), and the Deeps AV18 row (579.1162/187; LC x Glauber 582.5164; by_ps and transverse as quoted). No validation/reference JSON depends on the run default. The T6 REPORT row is wave-independent (clause counts identical on both waves). The t4_djangoh and t2_boeglin hashed blocks are untouched by D6. The quotes from Rodning-Knutson (Eq. (10), p. 906), AV18 Table X and Wiringa 2014 are accurate. What does not hold: (1) DD-4 would refuse a knob the run reads (--beta on 7Li at the default T1 tier). (2) The knob-matrix retake instruction ("base=HUL") makes a cell pass for the wrong reason and drops coverage of the Hulthen read branch. (3) The RC-clamp re-measurement uses a different recipe from the published one, and the claimed staleness is false (245 and 520 reproduce exactly). (4) About 20 sites are missed: published Hulthen-default numbers stated with a default-run recipe (the FSI D3 numbers in fsi.hpp, CONVENTIONS, bindings.cpp, PHYSICS_CHANNELS 321/344, USAGE 1667/2559, OPEN_ITEMS_SOLUTIONS), plus citation-gate exemptions and fingerprints in the shift zones. (5) The citation-gate procedure omits that --fix never touches ALLOW entries and that RANGE_ALLOW/UNNAMED_ALLOW keys carry the target range or line. (6) The runtime cost misses a one-time first-generate cost of about 0.1-1.3 s on 7Li VMC. (7) Several risks are underrated: 7Li FSI clipping on the default (4.4-9.5 % of the grid, events at w_max = 50); HepMC3 files carry no record of the wave; I6's re-measured sample-level numbers move again at I7; the D/S evidence for the alpha-d channel itself is not examined (the Hulthen alpha-d pair has the opposite S-D sign, eta = +0.070 against George-Knutson -0.025); and re-pointing the Deeps verdict wave contradicts the record's own rule for T2 sec. 5.2. No benchmark verdict silently flips.

### Sites the design missed (21)

| file | line | what |
|---|---|---|
| `include/lipolgen/fsi.hpp` | 138-150 | D3 numbers stated for the default 6Li-tag recipe (--events 20000 --seed 1234): Sum w 10419.07 / 11632.10, ratio up to 68.52, '99.50 % of the 20 000 events differ by more than 1 %'. Under D6 that recipe gives 64.5 % and 10.14. Needs '--cluster-wave hulthen' in the recipe or a re-measurement. |
| `docs/CONVENTIONS.md` | 455-467 | The same D3 claim with its recipe: '99.50 % ... factor 68.52 ... 0.520239 vs 0.582899 -- the number meta["fsi_survival"] carries'. On the default run fsi_survival becomes 0.797094 / 0.819119. The record lists only CONVENTIONS 95-145 and 328-395. |
| `python/bindings.cpp` | 698-713 | Comments 'individual event weights by up to a factor 68.5' and '(131.0 vs 160.0 mb, 0.520 vs 0.583 on the 6Li alpha tag)'. These are default-run D3 numbers; the record lists only bindings 2730-2735, 2824-2851, 4473-4486 and 4672-4675. |
| `docs/PHYSICS_CHANNELS.md` | 321-322 | 'its event weights average to the survival probability (0.520 at 40 mb on the 6Li S+D tag)'. The default becomes 0.797. This is a wrapped paragraph above the exemption keys at 342-518, so any reword must keep its line count. |
| `docs/PHYSICS_CHANNELS.md` | 344 | The D3 table row: the recipe '--channel tagged-6Li-alpha --events 20000 --seed 1234' with Sum w 10419.07/11632.10, 0.520239/0.582899, 68.52 and 99.50 %. The row also cites tests/test_fsi.cpp:384 as an UNNAMED point, so a dated comment added above that TEST_CASE in Step 1 breaks the citation, and --fix does not repair unnamed citations. |
| `docs/USAGE.md` | 1665-1668 | '0.520 at 40 mb, 0.671 at 20 mb on the 6Li alpha tag (S+D)'. On VMC these are 0.797 and 0.864. This is outside every USAGE range the record lists. |
| `docs/USAGE.md` | 2559 | '1.2-2.6 % of tagged events clip on the band'. Measured on the published recipe, VMC gives 0.080 % (6Li) and 2.535 % (7Li). |
| `docs/OPEN_ITEMS_SOLUTIONS.md` | 61, 85, 414-444, 832-833, 858, 911-913 | The D3 FSI table (99.50 %, 68.52, 0.520239) stated with the default-run recipe, and the RC clamp '245 of 20 000 (1.2250 %, seed 1)' / '520 of 20 000 7Li (2.6 %, seed 11)' / '1.2-2.6 %'. Lines 832-833 also state the recipe the record calls unstated (P_z = 0, 'the suites' fill'). The record lists only lines 55, 95, 116, 162 and 2330-2345. |
| `docs/theory/SPIN32_FINITE_GAMMA.md` | 156 | '`generate_tagged.cpp:123`' is a citation without a directory. The gate's PATH regex does not parse it, so DD-4's insertion above line 123 leaves it silently stale. The record names only the parsed :151 and :208 citations. |
| `python/lipolgen/cli.py` | 219-222 | The --cluster-beta and --p-d help strings ('short-range scale of the cluster radial waves', 'D-state probability (6Li alpha tag)') carry no 'hulthen only' qualifier. Under D6 both are unread on the default 6Li-alpha run (beta is read only on 7Li at T1). The same two flags are listed as live in python/README.md:375. |
| `docs/benchmarking/figures/price_defaults.png` | - | A committed figure of 'the ratio to the shipped default' with the series 'cluster_wave = vmc' (price_defaults.py:1090). make_report --plots does not draw it ('price_defaults.png is not a harness's'), so after D6 it is stale and nothing regenerates it. |
| `docs/PHYSICS_CHANNELS.md` | 267 | 'VmcAV18 opens exactly four files relative to data_dir()' omits vmc/deuteron/fdeut.av18, which the d-p control, the embedded deuteron and the T1 breakup read under VmcAV18. This staleness predates D6, but the row now describes the default run's data dependency (DD-5). |
| `docs/DEVELOPMENT_PLAN.md` | 270-272 | 'the deuteron internal wave function is Hulthen S+D at P_D = 0.045 ... as good as the cluster model gets without VMC overlaps'. After D6 the default T1 deuteron is AV18 fdeut. |
| `docs/open_items/run_2026-09-03/AUTHOR_DECISIONS.md` | 503-517 (B9 d/e), 804-830 (B15) | B9(e) prices option (ii) as 'the deuteron-channel default in tagged.hpp / breakup.hpp (BreakupOptions::source), re-dump tagged.json, T25', with tagged.json moving at rtol 1e-12 (B9(d)). DD-1 does none of this, and the amendment box must say so. B15 lists cluster_vmc_mc_sigma 0 -> +-1 as opt-in physics left at its no-op default; under D6 the band applies to the default lithium tags. |
| `tests/test_pipeline.cpp` | 151-152 | Comments quote tagged asymmetries measured on the Hulthen default (-1.72e-03 6Li, +6.41e-04 7Li, -5.73e-03 d-p) as 'P_e IS read' evidence. |
| `validation/benchmarks/t2_boeglin_nk.py` | 817-819 | Plot docstring: 'the default `hulthen` ... around the default wave (grey)'. These are not in the record's label list (42-45, 56-58, 358-359, 367-369, 494, 662, 858). |
| `validation/dump_polligen_reference.py` | 819-821 | tagged.json is dumped from polligen's DEFAULT channels (li6_alpha_channel(), etc.). If the sibling mirrors D6 in its own defaults, the next dump silently re-pins VMC against LiPolGen's Hulthen factories. The record says KEEP; passing wave='hulthen' explicitly removes the coupling. |
| src/core/pipeline.cpp (cited in docs/PHYSICS_CHANNELS.md) | 847, 1877 (ALLOW); 1116 (UNNAMED_ALLOW) | Named-point exemptions ALLOW[(pipeline.cpp, 847, coherent)] and ALLOW[(pipeline.cpp, 1877, optics_lumi_factor)], and the unnamed pin 'kernel && is_tagged(channel)' at 1116, all sit below the validate() block D6 rewords (482-506). --fix 'never touch[es] an ALLOW entry' (gate docstring line 239; code at 2114), so a line-count change there needs a hand re-point plus a re-key. |
| include/lipolgen/constants.hpp; tests/test_pipeline.cpp; README.md; docs/CONVENTIONS.md (UNNAMED_ALLOW targets) | constants.hpp:139; test_pipeline.cpp:608; README.md:128; CONVENTIONS.md:415 | Unnamed exemptions whose key carries the target line, each in a D6 shift zone: the constants.hpp 124-133 reword, Step 1's pipe()/call-site edits at 73-552, the README 38-41 reword, and the CONVENTIONS 95-145/328-395 rewords. §9.3 names only generate_tagged.cpp:123 and cli.py:221. |
| validation/check_physics_channels_links.py RANGE_ALLOW (targets) | tests/test_pipeline.cpp:608-634; tests/test_tagged.cpp:396-412/396-413; include/lipolgen/pipeline.hpp:931-934; docs/USAGE.md:2795-2992; docs/CONVENTIONS.md:207-211, 325-327 | RANGE_ALLOW keys carry the TARGET range, so a D6 edit above any of these ranges (Step 1 test edits, the pipeline.hpp 766-822 doc reword, the USAGE and CONVENTIONS rewords) retires the exemption. --fix only prints 'RANGE_ALLOW[...] should become ...'. The record mentions only moves of the citing line. |
| docs/USAGE.md (fingerprinted block) | 1978-2169 | A fingerprinted block cited twice by PHYSICS_CHANNELS.md ('## 5. Far-forward routing'). It contains the 1994-1998 tag-fraction text D6 relabels, so it becomes EDITED and needs a human re-read plus --record-ranges, where R5 can refuse it. It is absent from §9.3's list of blocks edited in content. |

### Statements found wrong (12)

1. DD-4 (§3.4): 'Under the VMC default both are silently not read ... refuse `--beta`/`--pd` unless `hulthen`'. This is false for --beta on 7Li-alpha. generate_tagged.cpp never sets cfg.tier, so it runs at the default T1 (pipeline.hpp:885 `Tier tier = Tier::T1`). There cluster_beta is read through the triton breakup (pipeline.cpp knob row `beta_read = tag && (!vmc || (7Li && T1))`). Measured on VMC: the knob row says 'read' and the output hash moves between beta 0.30 and 0.33. Separately, --pd is never read on 7Li or d-p even on Hulthen. The refusal must follow the knob table, not the wave.
2. §5.1, pipeline.hpp 824-833: '`cluster_beta`, `p_d`, `cluster_vmc_mc_sigma` | KEEP (not read at the default; the knob table says so)'. After D6, cluster_vmc_mc_sigma IS read on the default lithium alpha tags (`mc_read = vmc && li_tag`), and validate() then accepts a non-zero band on a default config.
3. §5.3 knob provenance: 'the "on hulthen" cells get `base=HUL`; the plain cells run on the default'. V() does not copy `base` into `var` (test_knob_provenance.py:204, and the variant runs `_run(spec, variant.var, ...)` at :590); the existing VMC cells repeat VMC inside var. With base=HUL alone, the '0 -> 1 on hulthen' cell's variant runs on the VMC default. On the lithium tags it then passes because the WAVE moved the hash, which breaks the matrix's one-knob-per-cell rule. Running the plain beta / p_d cells on the default makes them duplicates of the 'on vmc' cells, so the Hulthen READ branch of cluster_beta and p_d is no longer tested. The fix is base=HUL together with var=dict(HUL, ...).
4. §4.4 / §13 / summary note (6): 'The 245 at seed 1 was not reproduced with the CLI-default fill (237); the published recipe is not stated.' The recipe is stated in OPEN_ITEMS_SOLUTIONS.md:832-833 and 911-913 (20 000 events, seed 1 / seed 11, P_z = 0 'the suites' fill'). tensor_thirds_plan(0.0, 0.6) at seed 1 gives exactly 245; the test_rc.py 7Li fill at seed 11 gives exactly 520 (2.600 %). The proposed replacements ('237 (1.185 %) -> 14 (0.070 %)'; '698 (3.49 %) -> 700 (3.50 %)') come from different recipes. On the published recipe VMC gives 16/20 000 = 0.080 % (6Li, seed 1; seeds 1-5 give 16/20/21/30/34) and 507/20 000 = 2.535 % (7Li, seed 11).
5. §4.4: 'Under VMC the FSI weight table clips ... 0 everywhere on Hulthen.' The d-p control clips 1.451 % on Hulthen at 40 mb (the record's own §4.3 table shows 0.0145). The statement also omits 7Li on VMC: 4.41 % / 9.52 % (cluster / nucleon, 40 mb) and 2.91 % / 4.98 % at 20 mb.
6. DD-5.3: 'A C++ program run against an installed library outside the source tree, with no $LIPOLGEN_DATA_DIR, now throws at its first default tagged run. Options: install data/ in that branch too (one line)'. CMakeLists.txt:74-80 compiles in LIPOLGEN_DATA_DIR_DEFAULT = ${CMAKE_SOURCE_DIR}/data, so an installed library keeps finding data/ while the source tree exists; it throws only once that tree is gone or moved. Installing data/ alone is also not enough: data_dir() would still point at the source tree unless the LIPOLGEN_DATA_DIR cache variable is set to the install location.
7. Summary note (3): 'The wave port gate is validation/reference/tagged.json ... in tests/test_tagged.cpp and tests/test_cluster.cpp, pinned on the factory defaults'. tests/test_cluster.cpp:55-70 builds `Wave{l, prob, beta}` from the JSON's own fields and pins no factory default. Only test_tagged.cpp's li6_model() / li7_model() / deut_model() (74-85) use the factory defaults.
8. §4.1: 'The tables cost ≲ 7 ms and ≈ 2.4 MB on first use ... Not a cost worth designing around.' Only the factory and TaggedModel were timed. The first generate() pass of a 7Li VMC pipeline carries an extra one-time cost. Measured in fresh processes, first call vs Hulthen: 50k events 0.217 vs 0.120 s; 200k 0.764 vs 0.449 s; 400k 2.035/2.051 vs 0.718/0.695 s; 800k 2.587 vs 1.905 s. Later VMC calls run at par (0.11 / 0.40 / 0.70-0.90 / 1.67 s). 6Li and d-p show no difference, and construction takes 60-73 ms on both waves.
9. §12.3 row 9: 'Decided 2026-10-27 (D6): option (ii) ... Leg (2): tagged.json stays polligen's Hulthen port gate, untouched, since the factories keep the analytic default.' AUTHOR_DECISIONS §B9(e) defines option (ii) as changing the factory and BreakupOptions defaults, re-dumping tagged.json and T25, and §B9(d) has tagged.json moving at rtol 1e-12. DD-1 is a different implementation, and the amendment calls it option (ii) without saying so.
10. §7.1: 'Deeps ... the rule follows the default' and 'the verdict row becomes av18/lc'. The declaration identifies the verdict wave as 'the shipped default, `deuteron_channel()` = the analytic Hulthen pair at P_D = 0.045 (the tagged-channel default)', and under DD-1 `deuteron_channel()` stays Hulthen. Re-pointing the headline to AV18 with both chi2 values already known (627.85 and 579.12) is a post-hoc re-reading. §7.4 reserves that same act to the author for t2_deuteron_static, and §7.2/§7.3 keep Hulthen for BONuS and Boeglin. The verdict is FAIL either way, but the treatment is inconsistent.
11. §9.3 steps 2-3: '`--fix` for MOVED named points and ranges; re-point unnamed citations by hand and update their UNNAMED_ALLOW key' and 'A RANGE_ALLOW exemption is keyed (document, citing line, range), so a document edit that moves a citing line retires it'. --fix never touches ALLOW-exempted named citations (gate docstring line 239; code at 2114 only reports 'allow-listed on this line but the line no longer contains its pin'), and two of them are pipeline.cpp:847 and :1877. RANGE_ALLOW keys also break when the TARGET range moves, not only the citing line.
12. §10.1: 'a pre-I6 --cluster-wave vmc run must be byte-identical to a post-I6 default run.' The cluster_wave knob row's at_default also flips in that comparison (False -> True). The exception is stated only for the first comparison.

### Underrated risks (10)

1. The evidence chapter answers the author's D/S criterion with the deuteron eta_d only. The observable that moves (6Li A_zz^tag x0.43 / x0.21) is driven by the alpha-d D state. Computed here with t2_deuteron_static.hulthen_forms' residue formula (which reproduces the deuteron's 0.032478), the Hulthen alpha-d pair (kappa = 0.06066 GeV, beta = 0.30, P_D = 0.0867) has eta(6Li -> alpha d) = +0.0701. George-Knutson's restricted phase-shift value is -0.025 +- 0.006 +- 0.010 (LI6_ETA_DS_GK). The VMC tables give about -0.048 (USAGE.md:2923), with the documented S-D sign opposite to the Hulthen forms (src/core/tagged.cpp:114-116). This supports D6 on the alpha-d channel too, but it also means the 'quoted cluster-wave systematic hulthen - vmc' the record keeps Hulthen for mixes a wrong S-D sign with P_D and shape. As stated, that is not a model band.
2. A default tagged HepMC3 file records no wave. GenRunInfo carries only the tool name and version (src/hepmc/hepmc_writer.cpp:68-93), and nothing bumps the version. A file made after D6 by the same command is indistinguishable by metadata from a pre-D6 Hulthen file, while A_zz^tag differs by x0.43. The optional DD-5.1 banner does not reach the file (the ePIC chain's input).
3. Sample-level numbers that I6 re-measures (USAGE FSI table, README FSI row, RC clamp counts, tag fractions at fixed seed, test_fsi's 4000-event shape) depend on the event stream. That stream changes again at I7: with unpol_sf=ct18nlo at seed 1234 the k, cos theta_k and x columns are not identical (Sum w 15989.349 -> 15989.374). plan_D7 treats the tag fraction as carrying no structure function. Either re-measure after I7, or put '--unpol-sf toy' in each recipe.
4. FSI on the default 7Li tag (opt-in --fsi) is much worse than the 6Li numbers the record quotes. The grid clips 4.41 % (cluster) and 9.52 % (nucleon) at 40 mb, and 2.91 % / 4.98 % at 20 mb. In 20 000 events (seed 1234, helicity-flip, 40 mb), glauber-nucleon has 4 events at the w_max = 50 ceiling, and 12 events with w > 10 carry 2.5 % of Sum w; on Hulthen the maximum weight is 1.28. The 7Li survival moves 0.365 -> 0.605 (PHYSICS_CHANNELS:349 quotes 0.365 as 'measured').
5. Step 1 trap. If the new pipe() wave parameter passes Hulthen for the inclusive and coherent call sites too (the four-momentum test at test_pipeline.cpp:232, the threaded test at :768, the throughput test at :957), Step 1 still passes. Step 2's new validate() rule then refuses Hulthen there, and Gate 2a reports failures that are not on the predicted list. The same holds for any Python fixture passing cluster_wave='hulthen' to an inclusive make_config once the guard lands. Only the tagged call sites may pass Hulthen.
6. Published numbers nobody regenerates: docs/benchmarking/figures/price_defaults.png (ratios to 'the shipped default', not drawn by make_report), and the D3 FSI recipe sites listed above. Each states a default-run command that will silently reproduce different numbers after D6.
7. Every default inclusive and coherent file's meta will carry knob_provenance cluster_wave value 'vmc' (status refused, at_default). Registry row 11 exists to prevent exactly the 'VMC 6Li inclusive' misreading, and any reader that takes the value without the status will make it.
8. Closures sit nearer their floors on VMC. The spectator-marginal closure passes with ndf_k = 13 against REQUIRE(ndf_k > 10), chi2/ndf 1.52 against < 2.0. The record handles only the A_zz closure (ndf 7 < 8). Any future seed, grid or N change on the default can trip the marginal test.
9. The sibling coupling runs both ways. dump_polligen_reference.py reads polligen's defaults, so a sibling-side flip to wave='vmc' would silently change the port gate's provenance on the next dump.
10. Proposed DD-4 refusal: '--pd unless hulthen' still accepts --pd on 7Li and d-p, where it is never read on either wave. The example writes cfg.p_d for every channel, so the refusal must also be channel-scoped.

### The checker's notes: what it ran and confirmed, and what it did not

For the supervisor. No repository file was edited; the git tree is clean apart from the existing untracked run_2026-10-27/. The sibling was only grep'd and pdftotext'd, never Read, and its git status is clean. Record sha256 3b30d171...76aa1c7a matches.

Verified by reading or running:
- Flip-experiment logs: design A gave 14 + 1 pytest failures, matching the record.
- Reproduced exactly: the Deeps AV18 numbers; the FSI table on both waves; A_zz^tag and the dilutions for 6Li and d-p.
- Line numbers and quotes of the listed sites are accurate: pipeline.cpp:497/2803-2818, cli.py:805/1312, bindings 2730-2735/2826/2832/2841/4473-4486/4672-4675, generate_tagged 87-88/145-146/164-168, gate line 1103, 46 RANGE_ALLOW and 18 UNNAMED_ALLOW keys on PHYSICS_CHANNELS.md.
- t2_boeglin's design_block_sha256 hashes only the marker-delimited text plus the _RULES constants, not DEFAULT_WAVE, as the record says. t4_djangoh runs an inclusive proton with no cluster_wave, so D6 cannot touch it.
- T6: running measure() on the three tagged specs with an explicit --cluster-wave gives identical clause counts on both waves (file bytes differ, but REPORT stores none).
- The Deeps test at test_bench_deeps.py:144-148 still passes on AV18 (13.23 > 5, 0.61 < 1).
- The citation-count table in the record's §9.3 is correct.

Scratch, all under /tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/:
- d6chk/cost.py, cost_tier.py, warm.py, warm2.py: runtime.
- d6chk/t6_wave.py: T6 clause counts per wave.
- d6chk/rk_layout.txt, d6chk/av18.txt, d6chk/w14p1.txt: paper text.
- d6chk_records.txt: --records-only output (exit 0; plan_D6 contributes nothing).
- Inline python runs (RC clamp, FSI tables, alpha-d eta) are reproducible from the numbers quoted above.

The repo-root validation/benchmarks/__pycache__ holds pre-existing .pyc files (make_report, _plots, _report); all my imports used -B.

Suggested minimal corrections before I6:
- DD-4: refuse only where the knob row says NotRead.
- Knob matrix: base=HUL together with var=dict(HUL, ...), and keep Hulthen beta/p_d cells.
- RC clamp: re-measure on the suites' fill (P_z = 0, seed 1; the 7Li J = 3/2 fill, seed 11).
- Add the missing sites to the inventory.
- §9.3: add the ALLOW / UNNAMED / RANGE_ALLOW re-keys and the USAGE.md:1978-2169 re-record.
- Decide the Deeps verdict wave consistently with T2 sec. 5.2, or put it to the author.
- Consider recording the wave in HepMC3 GenRunInfo.
- Sequence the sample-level re-measurements after I7.
