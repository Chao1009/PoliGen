# Phase D — measured numbers

Every number in this file was produced on **2026-09-04** against the build in
`build/` of this tree, with `source env.sh`, by the recipe printed next to it —
except §D2.8, a **second pass dated 2026-09-05**, §D6, a **third pass of the
same date**, the **claim-residue pass of the same date** (§D7) and the
**fourth pass of the same date** (§D8), whose numbers were all re-measured that
day and say so in place. Nothing here is
quoted from another document without saying so.

Suites, **after every section of this file** (D2 landed first; D1 added 4 C++
cases and 23 pytests on top of it; D3+D4+D5 added 2 C++ cases and 17 pytests on
top of that; §D2.8 added 26 pytests and no C++ case; §D6 added 415 pytests and
no C++ case; §D7 added neither and moved one citation count; §D8 added 58
pytests -- 40 passing and 18 refused-base skips -- and no C++ case):

| gate | after D2 | after D1 | after D3+D4+D5 | after §D2.8 | after §D6 | after §D7 | after §D8 (this tree) |
|---|---|---|---|---|---|---|---|
| `build/lipolgen_tests` | 394 cases / 17 217 500 assertions / 1 skipped / 0 failed | 398 cases / 17 219 756 assertions / 1 skipped / 0 failed | 400 cases / 17 240 262 assertions / 1 skipped / 0 failed | 400 cases / 17 240 262 assertions / 1 skipped / 0 failed | 400 cases / 17 240 262 assertions / 1 skipped / 0 failed | 400 cases / 17 240 262 assertions / 1 skipped / 0 failed | **400 cases / 17 240 262 assertions / 1 skipped / 0 failed** |
| `python -m pytest python/tests -q` | 367 passed | 390 passed | 407 passed | 433 passed | 848 passed, 94 skipped | 848 passed, 94 skipped | **888 passed, 112 skipped** (measured mid-§D8; the committed tree `de1a040` runs **909 passed, 112 skipped** — see the note under this table) |
| `python3 validation/check_physics_channels_links.py` | 1115 references checked (strict), 97 ranges, 6 external, 0 broken, 6 allow-listed | 1119 references checked (strict), 97 ranges, 6 external, 0 broken, 6 allow-listed | 1122 references checked (strict), 97 ranges, 7 external, 0 broken, 6 allow-listed | 1130 references checked (strict), 97 ranges, 7 external, 0 broken, 6 allow-listed | 1144 references checked (strict), 97 ranges, 7 external, 0 broken, 6 allow-listed | 1145 references checked (strict), 97 ranges, 7 external, 0 broken, 6 allow-listed | **1145 references checked (strict), 97 ranges, 7 external, 0 broken, 6 allow-listed** |

**The `888 passed, 112 skipped` cell is a measurement taken BEFORE the section
that follows it finished.** It was measured during §D8 and not re-run after the
§D7–F fourteen-site pass and §D9 added their tests. Re-measured 2026-09-05 on
the committed phase-D tree: `de1a040`'s `python/tests` **collects 1021 tests
(909 + 112)** and runs **909 passed, 112 skipped** — the figure `STATUS.md`
row D and `phase_E_numbers.md` §E0.2 carry, and the 21-test gap
`phase_E_numbers.md:59-66` records as "not independently traced" is exactly
this. `git diff --stat de1a040 HEAD -- python/tests` is phase E's 16 tests
only (`test_doc_link_gate.py`, `test_mstw_sf.py`, `test_release_metadata.py`,
`test_spdx_headers.py`). Both `888` cells in this file are left as the
historical measurement they are.

(The "after §D2.8" reference count above is **1130**, measured on this tree on
2026-09-05 before §D6 touched anything; §D2.8's own text below says 1128, which
is what it measured when it was written. The difference is two citations and
its cause was not tracked down — the number here is the one this run measured,
not a correction of that one.)

§D2.8's six new references are the entry points it added (`pol_sf_is_read`,
`pol_sf_reach_report`, `pol_sf_unread_label`, `cell_rate_weights_pb`); its C++
column is unchanged because the four defects it closes are all on the run
SURFACE — `validate()`, `meta`, the banner and one message — and the C++ suite
carries no test of any of those.

§D7's single new reference is `python/bindings.cpp:608` `unpol_sf`, added
where a row had said "nothing records it"; measured by ablation on 2026-09-05
(removing that one citation and re-running the gate gives 1144). §D7 also
re-anchored **71** citations that its own edits had displaced — counted
position by position against the pre-`--fix` document, 1243 citations in both,
**0 of them landing in a different file**: 31 in `python/lipolgen/cli.py`, 20
in `include/lipolgen/pythia_bridge.hpp`, 6 each in `python/bindings.cpp` and
`docs/USAGE.md`, 4 each in `docs/PYTHIA_BRIDGE.md` and
`src/pythia/pythia_bridge.cpp`. Most went through `--fix`; the six use sites
that `--fix` cannot re-point (a symbol its file declares nowhere) and the four
ranges were moved by hand, as was the `dis_parton_fraction` ALLOW pin, from
:694 to :700. Then 5 fingerprints re-recorded (5 new, 0 changed, 5 dropped):
the three moved blocks re-hashed to their recorded sha256 **exactly**
(41b6be95…, 1ec7fcc4…, 6d27ba15…, i.e. moved and not edited), and the one
block that CHANGED (`docs/PYTHIA_BRIDGE.md:500-561`, 215e460d… →
9adde79e…) changed because §D7 corrected the "all fifteen sets" T0 claim
inside it. It added no test and moved no generated array; its C++ and pytest
columns are unchanged.

§D6's 14 new references are the knob-provenance row's own entry points
(`KnobStatus`, `knob_status_name` ×2, `KnobProvenance`, `KnobRunContext`,
`plan_has_beam_helicity`, `knob_provenance` ×3, `rc_band_applies`,
`rc_tail_applies`, `cluster_wave_name`, `triton_sf_name`,
`knob_provenance_lines`). Its C++ column is unchanged for the same reason
§D2.8's was, and its pytest column moves by 415 because the new file is a
MATRIX — one case per (spec × knob) cell, where a "spec" is a
(isotope, channel, plan) triple: the 7 (isotope, channel) combinations carry
EVERY knob under ONE plan each, and 5 further specs carry the plan axis
(helicity-flip / transverse-tensor / tensor-flip / pe = 0 on inclusive-⁶Li and
tagged-⁷Li-alpha) for the 10 plan-sensitive knobs — 12 specs, 71 variants,
**547 cells** (582 collected tests, re-measured 2026-09-05). It is NOT the
full (channel × plan) product: coherent, tagged-⁶Li-alpha, tagged-d-p and
inclusive-d are never run under helicity-flip, transverse-tensor or
tensor-flip.  That is the point of it.
`python/tests/test_knob_provenance.py` alone is **409 passed, 94 skipped in
58 s**, every skip naming the refusal that made that cell's base configuration
unbuildable; the remaining 6 of the 415 are the two new
`python/tests/test_sf_backend.py` cases for the plan axis (one of them
parametrised over the five channels). The 94 skips are new — no pre-existing
test skipped.

The 7th external citation is new in D4: the Pomeron flavour row now cites
`PartonDistributions.cc:2630`, the line at which PYTHIA's one `PomH1FitAB`
class — sets 3, 4 **and** 6 — assigns `xc = xcbar = 0.` unconditionally.
`python/tests/test_doc_link_gate.py::test_real_document_external_citations`
carries the count and was updated with it.

The Phase D baseline this run started from (commit `903fcc9`) was 390 cases /
17 217 416 assertions / 345 pytest / 1094 references. **No file under
`validation/reference/` changed in either section** (`git status` on that
directory is empty), and every suite that reads one — `tests/test_reference.
cpp`, `test_tagged.cpp`, `test_bookkeeping.cpp`, `test_spectator.cpp`,
`test_coherent.cpp`, `test_pipeline.cpp`, `test_b1_nuclear.cpp` and
`python/tests/test_kernel_reference.py`, `test_b1_model.py` — passes unchanged
at its own rtol.

---

## D1 — the ⁷Li rank-2 (tensor) sector: the zero, out loud

The physics research is `phase_D_li7_rank2.md`; its recommendation is **path
(c)** — ship the loud zero and the run-surface defects now, the quadrupole
gate next, and the α–t b₁ itself only after the unpolarised-backend decision
(D2 of that note's §7) is made once for both isotopes. **This section is what
was implemented and measured: the loud zero and four defect fixes. No b₁(⁷Li)
was implemented, and none of the α–t numbers of that note entered the code.**
The deferred design is `docs/OPEN_ITEMS_SOLUTIONS.md` §15.

### D1.0 The claim, and its exact scope

> A ⁷Li **inclusive** run's rank-2 sector — the tensor term of the φ-averaged
> rate, the cos 2φ (gluon transversity) amplitude and therefore A_zz — is
> **exactly zero**, because `default_inclusive_kernel` fills a rank-2 slot for
> **spin 1 only** and `InclusiveKernel::tables` dispatches the slots on the ion
> spin. It is zero *before* this change and zero *after* it; what changed is
> that the run now **says so** — in the banner and in `meta` — instead of
> returning silent zeros under `meta["b1_model"] = "miller"`.

Where the claim does **not** apply, and neither does the banner: a ⁷Li
**tagged** run (`--channel tagged-7Li-alpha`) carries the α–t alignment in the
event weight and is gated today at ⟨P₂(cos θ_k)⟩ = −T/5; the **coherent**
channel's tensor signal is the recoil azimuth. Both say so in `meta`, in their
own words, rather than borrowing the inclusive sentence.

### D1.1 The zero, measured

Default ⁷Li beam config 1 (e 10 GeV × ⁷Li 99.5 GeV/u), shipped scenario, grid
and optics; the plan is the honest ⁷Li A_zz contrast, pure |m| = 3/2 (T = +1)
against pure |m| = ½ (T = −1), unpolarised beam, θ_S = 0, 60 000 events,
seed 11:

```
sigma per category (pb)  [590952.42641509, 590952.42641509]   <- the SAME double
A_T = (s+ - s-)/(s+ + s-) = 0.0                                <- exactly
kernel tables at x = 0.2, Q2 = 5:   7Li  b1 = b2 = delta = 0.0
                                    6Li  b1 = 1.33502308280647e-03
                                         b2 = 5.340092331225881e-04
                                         delta = 1.4130453571323692e-03
meta b1_model = b1_unpol = "none (spin 3/2: no rank-2 input)"   <- was "miller"
```

The ⁶Li row is the control: the same call on the spin-1 ion is not zero, so
this is a missing **input**, not a broken kernel. The two per-category cross
sections are the same double — there is no "small" and no rounding.

### D1.2 What is NOT broken, measured the same way

With a rank-2 slot supplied by hand on a caller-built kernel
(`b1_32_func = +0.05·F1`, `delta_32_func = −1e−2·F1`), same plan, 40 000
events, seed 11:

```
sigma+ = 565388.96105743 pb   sigma- = 616515.89177276 pb   A_T = -0.043258077
meta b1_model = "caller-supplied kernel"
```

−0.0433 against the −b₁/F₁ = −0.05 of `SPIN32_FINITE_GAMMA.md` Eq. (48) with
the sample's own 1/(1 + εR) ≈ 0.865. **The sampler, the tensor weight, the
per-category cross sections and the provenance line are already correct at
J = 3/2**; what is missing is a physics input. This reproduces §1.6 of the
research note exactly, and it is the number the run banner quotes.

### D1.3 The D2 selectors reach ⁷Li — which is why the b₁ is deferred, not why it is missing

Task requirement: the `--unpol-sf` / `--pol-sf` selectors must reach whatever
nucleon input a ⁷Li construction would use. Measured, 2000 events, seed 11:

| channel | backend | σ [pb] | ratio | tensor sector | `meta["unpol_sf"]` |
|---|---|---|---|---|---|
| inclusive | `toy` | 590952.426415 | 1 | **0 exactly** | `toy` |
| inclusive | `--unpol-sf ct18nlo` | 470360.061948 | **0.795936** | **0 exactly** | `ct18nlo` |
| inclusive | `--unpol-sf mstw` | 467631.418990 | **0.791318** | **0 exactly** | `mstw` |
| tagged-⁷Li-α | `toy` | 589760.769797 | 1 | (in the weight) | `toy` |
| tagged-⁷Li-α | `--unpol-sf ct18nlo` | 467410.920663 | **0.792543** | (in the weight) | `ct18nlo` |
| tagged-⁷Li-α | `--unpol-sf mstw` | 465064.709179 | **0.788565** | (in the weight) | `mstw` |

and on the polarised side, A_∥ on `helicity-flip` at J = 3/2 (P_z = P_e = 0.7):
**5.29188064e−04** (toy) → **2.81290908e−04** (`--pol-sf nnpdfpol`), a factor
**0.531552** (2.81290908e−04 / 5.29188064e−04 = 0.53155188; this line read
0.531557 until 2026-09-05, which is not the quotient of the two numbers beside
it — re-measured 2026-09-05, ⁷Li, 10 events, seed 11, `helicity-flip` at
j = 3/2, P_z = P_e = 0.7).

Read this table the right way round. The selectors reach the ⁷Li kernels —
inclusive and tagged — and move the rate by ~20 %, so the F₁ a ⁷Li b₁ would
fold against **is** selectable today. **The tensor sector stays exactly zero
under every one of them**, because nothing fills the slot that backend would
feed. And that same flag is what the research measured as deciding the **sign**
of the α–t estimate at four of six x points, which is why the b₁ is held
(§15 of `OPEN_ITEMS_SOLUTIONS.md`, decision D2).

### D1.4 Defect F2 — the spin-1 plans at J = 3/2, and one correction to the note

Three of the four standard plans hard-code j = 1 (`tensor_thirds_plan`,
`transverse_tensor_plan`, `tensor_flip_plan`, `src/core/bookkeeping.cpp`);
only `helicity_flip_plan` takes `j`. Before:

```
make_plan tensor-thirds      REFUSED, advising "helicity-flip or transverse-tensor"
make_plan transverse-tensor  BUILT a spin-1 plan
make_plan tensor-flip        BUILT a spin-1 plan
```

After: all three are refused at any J ≠ 1, with a message that names
`helicity-flip` (the only plan that takes `j`, and a **vector** plan), names
`--channel tagged-7Li-alpha` (which does carry the alignment) and quotes the
measured zero with its configuration. `helicity-flip` still builds at J = 3/2
(`j = [1.5, 1.5]`, T = 0.4 at P_z = 0.7), and all four still build at J = 1.

**Correction to `phase_D_li7_rank2.md` §1.5 F2, measured:**
`transverse_tensor_plan` builds **one** category, not three, and the run that
used it did **not** run silently — the `Pipeline` constructor threw two frames
down, out of `InclusiveKernel::amplitudes`, with a message that named neither
the plan nor the fix. `docs/PHYSICS_CHANNELS.md` carried the same wrong claim
("silently runs three spin-1 categories") and is corrected in the same change.
The inclusive branch of the `Pipeline` constructor now makes the run-plan-spin
check the tagged branch has always made, so the hand-built and C++ routes get
the named refusal too:

```
Pipeline: run-plan category "cos2phi" is J = 1 but the inclusive kernel's ion
7Li is spin 1.5 -- tensor_thirds_plan, transverse_tensor_plan and
tensor_flip_plan all hard-code j = 1 ...
```

### D1.5 Defect F3 — the J = 3/2 population domain, derived and checked

The four J = 3/2 populations are fixed **uniquely** by (1, P_z, T, R₃) — the
4×4 system in `spin32_populations` is square — so an unphysical request is a
**domain** statement, not a solver failure. Inverting the system by hand:

        p(+3/2) + p(−3/2) = (1 + T)/2 ,   p(+3/2) − p(−3/2) = 0.9 P_z + 0.1 R₃
        p(+½)  + p(−½)  = (1 − T)/2 ,   p(+½)  − p(−½)  = 0.3 (P_z − R₃)

so positivity is exactly

> **|0.9 P_z + 0.1 R₃| ≤ (1 + T)/2  AND  |0.3 (P_z − R₃)| ≤ (1 − T)/2**,
> which at R₃ = 0 collapses to **1.8|P_z| − 1 ≤ T ≤ 1 − 0.6|P_z|**
> (hence |P_z| ≤ 5/6 = 0.833333 at all).

**Checked against the function itself: 0 disagreements** on a 41 × 41 (P_z, T)
grid at R₃ = 0 and on a 21 × 21 × 11 (P_z, T, R₃) grid. The refusal now prints
the offending m, its negative population and both edges at this P_z; the C++
and Python suites both re-derive the boundary rather than assert the text.

**Correction to `phase_D_li7_rank2.md` §1.5 F3, measured: the repro line in
that note does not reproduce.** `lg.make_plan("helicity-flip", j=1.5, pz=0.7,
pzz=0.6)` does **not** raise: `make_plan` never forwards `pzz` to
`helicity_flip_plan`, whose `use_explicit_pzz` stays false, so the call
silently returns the **max-entropy** fill with **T = 0.4**, not 0.6. The
refusal is real but is reached through `spin32_populations(0.7, 0.6)` or
`helicity_flip_plan` with `use_explicit_pzz = true`. (0.6 is outside the
domain above by 0.02: the edge at P_z = 0.7 is T = 0.58, where p(−½) = 0
exactly.)

**And that non-reproduction is itself a defect — F5, new here.** `--pzz` /
`make_plan(pzz=…)` is read by the three spin-1 tensor plans **only**, and
`--plan helicity-flip --pzz 0.6` has always produced a different alignment
from the one typed with nothing saying so. **No fill was moved** (that would
change shipped numbers); instead the run banner now prints the plan's own
recorded moments on every run —

```
  fill helicity-flip: J = 1.5, P_z = 0.7, T = 0.4  (from the max-entropy ladder
                                     at --pz; --pzz is not read by this plan)
  fill tensor-thirds: J = 1, P_z = 0.7, P_zz = 0.6
```

— and the `--pzz` help says which plans read it.

### D1.6 Defect F4 — the three-line segfault

`ClusterPartialWave::k` and `::phi` are two public vectors that Python assigns
one at a time, so between the assignments they have different lengths;
`rebuild()` built the spline anyway and the constructor indexed `y_[i+1]` off
the end of an empty `phi`. Reproduced before the fix (exit 139, core dumped):

```python
w = _lipolgen.ClusterPartialWave(); w.k = [0.1, 0.2, 0.3]     # Segmentation fault
```

`rebuild()` now **defers** (drops the cached spline) so that either assignment
order works, and `operator()`, `eval_sorted` and `norm2` throw a message naming
both lengths. Both orders reach the same value, 1.6875 at k = 0.15 on the
(0.1, 0.2, 0.3) × (1, 2, 1) test wave.

### D1.7 What was NOT done, and is therefore not claimed

* **No b₁(⁷Li) exists in the code, in any form, opt-in or otherwise.** The
  α–t convolution of `phase_D_li7_rank2.md` §5 was not implemented; none of
  its numbers (b₁/nucleon ~ 1.2e−04 … 1.0e−04, the ratio 2.99 against ⁶Li's
  orbital term) entered any header, constant or test. They stay in that note,
  which states plainly that their **sign** is not predicted by the inputs in
  this tree.
* **The A = 7 quadrupole gate (§4.2 of the note) was NOT committed**, and this
  section reports no Q(⁷Li). It needs one constant that is **not in this
  tree** — the measured Q(⁷Li) = −4.00(3) fm², which that note quotes from
  memory of the standard compilations and itself says "must be sourced before
  it is committed" (its decision D11). A number nobody in this session
  verified against a source is not a gate; it is filed as step 2 of §15
  instead, with the −4.06 fm² alternative and the ratio it moves (0.871 →
  0.858) recorded there.
* **Nothing was measured about the ⁶Li + n decomposition** (§5.5 of the note)
  beyond reading it.
* **The rank-3 (octupole) sector was not touched.** `R₃` remains at 0 by
  default and reaches no unpolarised-beam observable
  (`SPIN32_FINITE_GAMMA.md` §4).
* **`--pzz` still does not reach `helicity-flip`** (F5 above). Making it do so
  would move a shipped fill; it is written up in §15 as an author decision.

### D1.8 Reproduction

From this tree, `source env.sh`, then:

```python
import numpy as np, lipolgen as lg
from lipolgen import _lipolgen as _l

def T_plan():                       # the honest 7Li A_zz contrast, T = +-1
    return _l.RunPlan([_l.SpinCategory("T+", 1.5, [.5,0,0,.5], 0, 0.,0.,0., .5),
                       _l.SpinCategory("T-", 1.5, [0,.5,.5,0], 0, 0.,0.,0., .5)],
                      0.0, 0.0, 1.0)

# D1.1
p = _l.Pipeline(lg.make_config(isotope="7Li", events=60000, seed=11), T_plan())
print(p.sigma_per_category_pb(), p.generate(0, False)["meta"]["b1_model"])
print(_l.default_inclusive_kernel(_l.ion_by_name("7Li")).tables(0.2, 5.0).b1)

# D1.2
o = _l.InclusiveKernel.Options()
o.b1_32_func = lambda x, q2, f1: 0.05 * f1
o.delta_32_func = lambda x, q2, f1: -1e-2 * f1
c = lg.make_config(isotope="7Li", events=40000, seed=11)
c.kernel = _l.InclusiveKernel(_l.ion_by_name("7Li"), o)
s = _l.Pipeline(c, T_plan()).sigma_per_category_pb()
print(s, (s[0] - s[1]) / (s[0] + s[1]))

# D1.3
for ch, pl in (("inclusive", T_plan()),
               ("tagged-7Li-alpha",
                lg.make_plan("helicity-flip", j=1.5, pz=0.7, pe=0.7))):
    for sel in ("toy", "ct18nlo", "mstw"):
        q = _l.Pipeline(lg.make_config(isotope="7Li", channel=ch, events=2000,
                                       seed=11, unpol_sf=sel), pl)
        print(ch, sel, q.sigma_pb(), q.generate(0, False)["meta"]["unpol_sf"])

# D1.5 -- the domain, against the function itself
for pz in np.linspace(-1, 1, 41):
    for t in np.linspace(-1, 1, 41):
        try: _l.spin32_populations(float(pz), float(t), 0.0); ok = True
        except RuntimeError: ok = False
        assert ok == (1.8*abs(pz) - 1 <= t + 1e-12 <= 1 - 0.6*abs(pz) + 2e-12)
```

and from the command line, for the banner:

```
python -m lipolgen.cli --isotope 7Li --plan helicity-flip --pzz 0.6 --events 100
```

---

## D2 — `--unpol-sf` / `--pol-sf`: what the toy backends have been costing

The design and the injection map are in `phase_D_sf_injection.md`; the shipped
surface is `docs/USAGE.md` §2b. This section is the measurement.

### D2.0 The default is bit for bit — the claim the whole change rests on

`unpol_sf = pol_sf = Toy` leaves `InclusiveKernel::Options::f2_source` at the
shared `ToyF2` the pipeline has always named and leaves `g1_model` **unset**,
so the kernel class takes exactly the null branch it always took. Asserted by
`np.array_equal` — not `approx` — on `cell_xsec_pb`, `sigma_per_category_pb()`
and the generated `x`/`q2`/`y`/`weight` columns, on the **inclusive**, the
**coherent** and the **tagged-6Li-alpha** channel
(`python/tests/test_sf_backend.py::test_the_default_is_bit_for_bit_on_every_channel`).
All eight rtol-1e-12 reference gates in `validation/reference/_manifest.json`
still pass unchanged; no reference file moved and none was re-recorded.

### D2.1 The run-level numbers

⁶Li, `default_configs("6Li")[1]` (e 10 GeV × ⁶Li 99.5 GeV/u), the shipped
scenario and grid. σ is `Pipeline::sigma_pb()`; A_zz is the canonical thirds
estimator `(σ₊ + σ₋ − 2σ₀)/(σ₊ + σ₋ + σ₀)` on `tensor_thirds_plan(0, 0.6)`'s
per-category cross sections divided by P_zz = 0.6; A_∥ is
`(σ₊ − σ₋)/(σ₊ + σ₋)/(P_e·P_z)` on `helicity_flip_plan(1, 0.7, 0.7)`.

**Inclusive**

| backend | σ [pb] | σ ratio | A_zz | A_zz ratio | A_∥ | A_∥ ratio |
|---|---|---|---|---|---|---|
| `toy` (shipped) | 591846.17 | 1 | −5.1932309e−4 | 1 | −1.1717156e−3 | 1 |
| `--unpol-sf ct18nlo` | 472571.92 | **0.798471** | −6.5039705e−4 | **1.252394** | −9.6688802e−4 | 0.825190 |
| `--unpol-sf mstw` | 469556.45 | **0.793376** | −6.5457386e−4 | **1.260437** | −1.0570300e−3 | 0.902122 |
| `--pol-sf nnpdfpol` | 591846.17 | 1 (exactly) | −5.1932309e−4 | 1 (exactly) | −1.2708063e−4 | **0.108457** |
| both | 472571.92 | 0.798471 | −6.5039705e−4 | 1.252394 | −1.5915927e−4 | 0.135834 |

**Corrected 2026-09-05 — the A_zz column had been measured in a DIFFERENT
WINDOW from the A_∥ column beside it.** Until then it read
−5.1932192e−4 / −6.5039558e−4 / −6.5457238e−4, which are this estimator at
`Scenario::x_max` = **0.95** — the window the tagged `--inclusive-b1`
sub-table below is measured in, and the one that section's numbers still
carry. The header says the shipped window and the A_∥ column is demonstrably
in it: at `x_max` = 1.0 A_∥ is −1.1717156e−3 (toy) /
−9.6688802e−4 (ct18nlo) / −1.2708063e−4 (nnpdfpol), exactly the printed
values, while at 0.95 the same three are −1.1717160e−3 / −9.6688807e−4 /
−1.2708071e−4. One row, two windows. The A_zz values above are the
shipped window, re-measured 2026-09-05, and they agree to six digits with
this section's own cell-level whole-window A_zz row (−5.193231e−4 /
−6.503970e−4 / −6.545740e−4), which was right all along.

**Second correction, same day, same class — the σ column was in the 0.95
window too, and the sentence above used to claim it was not.** Both tables'
`toy` and `--pol-sf nnpdfpol` rows read 591846.16; measured 2026-09-05,
`Pipeline::sigma_pb()` on config 1 is **591846.168878 pb** in the shipped
window against **591846.161405 pb** at `x_max` = 0.95, so the shipped window
rounds to **591846.17** and the column now says that. Those two rows are the
only place in either table where the two windows are distinguishable at the
printed precision: `ct18nlo` is 472571.917912 against 472571.917091 and `mstw`
469556.451348 against 469556.450465, and both pairs round to the printed
472571.92 / 469556.45 in either window — which is why the slip showed in one
digit of two rows and nowhere else. The `--inclusive-b1 --x-max 0.95`
sub-table below keeps 591846.16: that is right for ITS window (measured there
with `--inclusive-b1` on, 591846.162356 pb).

The **ratios are unchanged** by either correction: A_zz 1.252394 and 1.260437,
and σ 0.798471 and 0.793376, all to six digits in either window — which is
exactly why a window slip in an absolute column survived every consistency
check the ratios could give.

**Tagged ⁶Li α**

| backend | σ [pb] | σ ratio | A_zz | A_∥ | A_∥ ratio |
|---|---|---|---|---|---|
| `toy` | 591846.17 | 1 | **0 exactly** | −3.5147303e−3 | 1 |
| `--unpol-sf ct18nlo` | 472571.92 | **0.798471** | **0 exactly** | −2.9002425e−3 | 0.825168 |
| `--unpol-sf mstw` | 469556.45 | **0.793376** | **−1.3773675e−16** (qual. (a)) | −3.1706263e−3 | 0.902097 |
| `--pol-sf nnpdfpol` | 591846.17 | 1 (exactly) | **0 exactly** | −3.8119671e−4 | **0.108457** |
| both | 472571.92 | 0.798471 | **0 exactly** | −4.7740842e−4 | 0.135831 |

Three qualifications, all of them load-bearing.

**(a) The tagged A_zz is zero on every backend by construction, not by
accident — and on `mstw` it is a FLOATING-POINT zero, not a bit-level one.**
At the shipped default the struck-cluster kernel carries no b₁
(`--inclusive-b1` is off, deliberately: on a tagged channel the α–d density is
already in the event weight), so the per-category **total** cross sections
carry no tensor term at all. On `toy`, `ct18nlo` and `nnpdfpol` they come out
**equal to the last bit** and the thirds estimator returns **0.0 exactly**
(`st[0] == st[1] == st[2]` is `True`). On `mstw` it is `False`: measured
2026-09-05, `sigma_per_category_pb()` = 469556.45134814945 /
469556.45134814945 / 469556.4513481495 — σ₀ sits **exactly one ULP**
(5.820766e−11 pb) above σ_±, and the estimator returns **−1.3773675e−16**.
That is round-off in the summation order — 2.7e−13 of the inclusive channel's
own toy A_zz — and not a tensor signal; but the qualification as it stood
("equal to the last bit ... exactly 0") was false for one of the four
backends, so it is stated per backend from now on. The tagged tensor signal
lives in the spectator-momentum-differential rate, not in σ_tot. With
`--inclusive-b1 --x-max 0.95` the tagged channel does move, by the same
factors the inclusive one shows:

| backend | σ [pb] | A_zz | A_zz ratio |
|---|---|---|---|
| `toy` | 591846.16 | −1.5579685e−3 | 1 |
| `ct18nlo` | 472571.92 | −1.9511901e−3 | 1.252394 |
| `mstw` | 469556.45 | −1.9637206e−3 | 1.260437 |

(`--x-max 0.95` is needed there for a reason worth recording: with
`--inclusive-b1` **and** `--unpol-sf ct18nlo` the topmost default cell
(x = 0.955, Q² = 167.3) gives `1 + w_avg = −0.1302` and `InclusiveSampler`
refuses the run. On the toy backend the same configuration runs. So the
`--x-max 0.95` rule that `--b1-model cdks|li6-convolution` already carries now
also applies to `--inclusive-b1` off the toy unpolarised backend.)

**(b) The A_∥ ratios are window integrals of a sign-changing integrand and are
smaller than any of their own parts.** Cross-section-weighted ⟨A_∥⟩ per x band
over the accepted cells, inclusive ⁶Li:

| x band | σ share (toy) | ⟨A_∥⟩ toy | ⟨A_∥⟩ nnpdfpol | ratio | ⟨A_∥⟩ ct18nlo | ratio |
|---|---|---|---|---|---|---|
| [0, 0.01) | 0.5304 | −2.34618e−3 | −3.81601e−4 | **0.1626** | −2.31586e−3 | 0.9871 |
| [0.01, 0.05) | 0.3018 | −1.44253e−5 | −1.50146e−5 | 1.0409 | +1.25934e−6 | (sign flip, both ≈ 0) |
| [0.05, 0.1) | 0.1092 | +2.81976e−4 | +2.55328e−4 | 0.9055 | +2.92954e−4 | 1.0389 |
| [0.1, 0.3) | 0.0532 | +6.84532e−4 | +8.02597e−4 | 1.1725 | +6.72640e−4 | 0.9826 |
| [0.3, 1) | 0.0054 | +1.82940e−3 | +1.71330e−3 | 0.9365 | +1.57717e−3 | 0.8621 |
| **whole window** | 1 | **−1.171591e−3** | −1.270670e−4 | **0.1084** | −9.667590e−4 | 0.8251 |

The window-integrated ×0.108 is **smaller than every band ratio**, because the
negative low-x band and the positive high-x bands partly cancel. Quote the
band, not the window factor, unless the window is the observable. (The
window rows here are the cell-level σ-weighted means. **For A_∥** they differ
from the run-level table above in the fifth digit — the run-level estimator
carries the per-category spin weights, the cell-level one does not.
**For A_zz they do not**: at P_z = 0 the three thirds categories carry equal
luminosity and the two estimators coincide, so §D2.1's corrected run-level
A_zz column and the whole-window row of the A_zz table below agree to **six**
digits (the `mstw` row's seventh differs by one unit in the last place
printed). That agreement is the check that caught the window slip recorded
above.)

The same table for A_zz, whose integrand also changes sign but only in the
smallest band:

| x band | ⟨A_zz⟩ toy | ⟨A_zz⟩ ct18nlo | ratio | ⟨A_zz⟩ mstw | ratio |
|---|---|---|---|---|---|
| [0, 0.01) | −1.63842e−4 | −2.35149e−4 | 1.4352 | −2.21361e−4 | 1.3510 |
| [0.01, 0.05) | −7.26687e−4 | −8.71005e−4 | 1.1986 | −9.23380e−4 | 1.2707 |
| [0.05, 0.1) | −1.38538e−3 | −1.39169e−3 | 1.0046 | −1.52800e−3 | 1.1030 |
| [0.1, 0.3) | −1.34964e−3 | −1.15550e−3 | 0.8562 | −1.24461e−3 | 0.9222 |
| [0.3, 1) | +1.83740e−3 | +1.58274e−3 | 0.8614 | +1.50748e−3 | 0.8204 |
| **whole window** | −5.193231e−4 | −6.503970e−4 | 1.2524 | −6.545740e−4 | 1.2604 |

**(c) The two selectors are not orthogonal.** `InclusiveKernel` builds its
default `ToyG1` on its **own** base `UnpolSF`, so `--unpol-sf` alone moves g₁:
with `f2_source = CT18NLO` and `g1_model` left at the default, ⁶Li at y = 0.5:

| x | Q² | F₁ ratio | g₁ ratio | A₁ = g₁/F₁ ratio |
|---|---|---|---|---|
| 0.05 | 5 | 0.9942 | 1.0544 | 1.0605 |
| 0.10 | 10 | 1.0711 | 1.1346 | 1.0593 |
| 0.30 | 15 | 1.2038 | 1.3050 | 1.0841 |
| 0.50 | 25 | 1.0014 | 1.0551 | 1.0537 |

So `--pol-sf toy` does **not** mean "g₁ unchanged", and any plot must state
both settings.

### D2.2 Where the 20 % rate change comes from

Pointwise F₂ᵖ at Q² = 10:

| x | `toy` | `ct18nlo` (ratio) | `mstw` (ratio) |
|---|---|---|---|
| 0.01 | 0.666267 | 0.613967 (0.9215) | 0.543303 (0.8154) |
| 0.05 | 0.467304 | 0.472435 (1.0110) | 0.423859 (0.9070) |
| 0.10 | 0.385242 | 0.427421 (1.1095) | 0.390353 (1.0133) |
| 0.20 | 0.284398 | 0.366212 (1.2877) | 0.344625 (1.2118) |
| 0.30 | 0.207952 | 0.285077 (**1.3709**) | 0.283824 (1.3649) |
| 0.50 | 0.0915026 | 0.115664 (1.2641) | 0.132491 (**1.4479**) |
| 0.70 | 0.0233479 | 0.0225094 (0.9641) | 0.0301786 (1.2926) |

and F₂ⁿ/F₂ᵖ, which is what the T1 and T2 struck-nucleon **species** draws use:

| x | `toy` | `ct18nlo` | `mstw` |
|---|---|---|---|
| 0.05 | 0.9625 | 0.9218 | 0.9271 |
| 0.20 | 0.8500 | 0.7219 | 0.7427 |
| 0.50 | 0.6250 | 0.5035 | 0.4856 |

The toy's n/p ratio is the straight line `clip(1 − 0.75x, 0.25, 1)`; it is up
to **24 %** above the fits at x = 0.5. That is why the T2 PYTHIA bridge had to
be wired to the same backend in the same change (D2.5).

### D2.3 The polarised sector, and a sign

⁶Li per-nucleon g₁A, `g1_model` swapped from the default `ToyG1` to
`LhapdfG1("NNPDFpol11_100", 0)`:

| x | Q² | g₁A `toy` | g₁A `nnpdfpol` | ratio |
|---|---|---|---|---|
| 0.01 | 2.5 | −0.0847927 | −0.0567341 | 0.6691 |
| 0.05 | 5 | +0.0323251 | +0.0218549 | 0.6761 |
| 0.10 | 10 | +0.0343653 | +0.0382812 | 1.1139 |
| 0.20 | 10 | +0.0255593 | +0.0344831 | **1.3491** |
| 0.30 | 15 | +0.0184751 | +0.0239836 | 1.2982 |
| 0.50 | 25 | +0.00788548 | +0.00773991 | 0.9815 |
| 0.70 | 50 | +0.00196416 | +0.00124874 | **0.6358** |

×0.64 to ×1.35 over the generator window, and not monotone. The **neutron** is
worse, and it is a sign rather than a factor — at Q² = 10:

| x | g₁ⁿ `toy` | g₁ⁿ `nnpdfpol` | ratio |
|---|---|---|---|
| 0.05 | −0.242789 | −0.248012 | 1.0215 |
| 0.10 | −0.0800273 | −0.136452 | 1.7051 |
| 0.20 | −0.0113549 | −0.0607773 | **5.3525** |
| 0.30 | **+0.00520678** | **−0.0273959** | **−5.2616** |
| 0.50 | +0.00778817 | −0.000370697 | −0.0476 |
| 0.70 | +0.00247092 | +0.00221761 | 0.8975 |

`ToyG1`'s a1n(x) = −0.07(1−x)² + 0.8x^2.2 crosses zero near x ≈ 0.25 and is
positive above it; NNPDFpol1.1's g₁ⁿ stays negative to x ≈ 0.6. **The shipped
toy g₁ⁿ has the wrong sign over roughly 0.25 < x < 0.6.** On isoscalar ⁶Li the
proton term dominates the sum and this mostly hides; on the neutron-tagged
`TaggedDeuteronP` channel it does not. Neither the ⁷Li α tag (triton) nor the
d+p tag (free neutron) was measured end to end here — see "not done" below.

### D2.4 The grid clause

| set | XMin | Q range | Q²Min |
|---|---|---|---|
| CT18NLO | 1e−9 | Q ≥ 1.295 GeV | **1.677025** |
| NNPDFpol11_100 | 1e−5 | 1 ≤ Q ≤ 316.228 GeV | 1 |

The accepted cells of the shipped ⁶Li run span Q² = 1.0542 … 1897.2 (the
`generator_scenario` cut is `q2_min = 0.7`; the *grid* is what stops at
1.0542). Below its own grid LHAPDF does **not** freeze:

| x | Q² | `toy` F₂ᵖ | `ct18nlo` F₂ᵖ | toy/ct18 | `mstw` F₂ᵖ |
|---|---|---|---|---|---|
| 3e−4 | 0.7 | 0.681248 | 0.291574 | **2.336** | 0.492329 |
| 3e−4 | 1.0 | 0.681248 | 0.400770 | 1.700 | 0.549911 |

Fraction of the accepted cell cross section below Q² = 1.677025:

* weighted by the **toy** run's own cells: **42.328 %**;
* weighted by the **CT18NLO** run's own cells (which is what a CT18NLO run
  records): **36.179 %**.

That second number is what `meta["unpol_sf_below_grid_frac"]` carries and what
the banner prints; it is computed from `LhapdfSF::q2_min()`, i.e. the loaded
grid's own metadata, never from a typed-in constant.

**There is no `--q2-min` CLI flag** (`--x-max` is the only scenario knob on the
command line), so moving the window off the extrapolation means raising
`Scenario::q2_min` from Python, and it costs rate. Measured on the shipped ⁶Li
`ct18nlo` run:

| `Scenario::q2_min` | below-grid fraction | σ [pb] |
|---|---|---|
| 0.7 (shipped) | 0.361791 | 472571.92 |
| 1.7 | **0 exactly** | 301599.80 (×0.638) |
| 2.0 | 0 exactly | 247731.98 (×0.524) |

A backend that reports no floor writes **NaN**, never 0 — PYTHIA's `MSTWpdf`
keeps its `qsqmin` private, and a `Custom` object or a caller-supplied kernel
reports nothing. (The paragraph above was glued into the table's last row
until 2026-09-05, which made the `q2_min` = 2.0 row unreadable.)

### D2.5 R, and why it is not in this change

Over the 3051 accepted cells of the shipped ⁶Li window:

* `r_sigma_lt` spans **0.004622 – 0.17630**; `r1998` spans
  **0.012430 – 0.40100**;
* `(1 + r1998)/(1 + r_sigma_lt)` spans **0.9203 – 1.1910**, cross-section-
  weighted mean **1.1229**;
* **38.182 %** of the accepted cell cross section lies **outside** R1998's own
  stated support (`R1998_X_MIN/MAX` = 0.005/0.86, `R1998_Q2_MIN/MAX` =
  0.5/130), where `r1998(..., clip = true)` returns a clipped boundary value.

R1998 is therefore not a drop-in replacement default, and two independent
10–40 % movements landing in one change would be unattributable. R stays its
own axis with no CLI flag, exactly as before this change.

### D2.6 Reproduction

```bash
cd /path/to/LiPolGen && source env.sh
# the run-level table of D2.1
python - <<'PY'
import numpy as np, lipolgen as lg
from lipolgen import _lipolgen as _l
PZ, PZZ, PE = 0.7, 0.6, 0.7
def run(channel, plan, **kw):
    return _l.Pipeline(lg.make_config(channel=channel, config=1, events=10,
                                      seed=1, **kw), plan)
for ch in ("inclusive", "tagged-6Li-alpha"):
    for lbl, kw in (("toy", {}), ("ct18nlo", dict(unpol_sf="ct18nlo")),
                    ("mstw", dict(unpol_sf="mstw")),
                    ("nnpdfpol", dict(pol_sf="nnpdfpol"))):
        t = run(ch, lg.tensor_thirds_plan(0.0, PZZ), **kw)
        h = run(ch, lg.helicity_flip_plan(1.0, PZ, PE), **kw)
        st = np.asarray(t.sigma_per_category_pb())
        sh = np.asarray(h.sigma_per_category_pb())
        azz = ((st[0] + st[1] - 2*st[2]) / st.sum()) / PZZ
        apar = ((sh[0] - sh[1]) / sh.sum()) / (PE * PZ)
        print(ch, lbl, t.sigma_pb(), azz, apar)
PY
# the below-grid fraction of D2.4, as the run itself records it
python -c "
import lipolgen as lg
from lipolgen import _lipolgen as _l
p = _l.Pipeline(lg.make_config(events=10, unpol_sf='ct18nlo'),
                lg.tensor_thirds_plan(0.0, 0.6))
print(_l.unpol_sf_grid_report(p))"
```

### D2.7 What was NOT done, and is therefore not claimed

* **Only the ⁶Li inclusive and the ⁶Li α-tag channel were measured end to
  end.** The ⁷Li α tag (triton `dis_target`) and the d+p tag
  (`NEUTRON_TARGET`) were **not** — and those are exactly the two where the
  toy n/p ratio and the toy g₁ⁿ sign do the most damage. Nobody should quote
  a factor on those channels from this file. *(Partly answered by §D2.8: the
  d+p tag was RUN on 2026-09-05, and on the shipped window the pair
  `ct18nlo|mstw` + `nnpdfpol` under `--plan helicity-flip` at `--pe` ≠ 0 has
  no run to quote — it is refused. No factor for that channel is quoted here
  even now; what §D2.8 records is the refusal and its cure.)*
* **MSTW's grid boundary was not characterised** the way CT18NLO's was. Its
  below-grid F₂ᵖ was probed at the two points in D2.4 and looked milder; that
  is an observation, not a measured fraction, and `MstwSF` reports no floor
  so the pipeline writes NaN rather than guessing one.
* **No EMC hook was wired or measured**, because none exists to wire: no site
  in `src/`, `python/` or `examples/` sets `InclusiveKernel::Options::emc_ratio`.
* **No `--r-model` was added or measured** beyond the D2.5 statistics.
* **The A_∥ and A_zz estimators here are the run plans' own thirds/flip
  estimators**, not an unfolded extraction; they are comparable across
  backends because the plan is identical, and they are not a physics
  prediction of the asymmetries.

### D2.8 Second pass, 2026-09-05 — four reach defects, found by RUNNING the selectors

An adversarial review ran `--unpol-sf` / `--pol-sf` on every channel instead of
reading the wiring, and confirmed the unpolarised selector's reach by value on
all five kernels. What it refuted was **the claim as written at the run
surface**: in four places the run recorded, printed or documented a reach it
did not have. All four are now closed; every number below was re-measured on
this tree on 2026-09-05, at `default_configs("6Li")[1]` (config 1) and the
shipped window unless a row says otherwise.

#### (1) `--pol-sf` does not reach the coherent channel

The coherent yield is `f_coh(x)` times the **unpolarised** cell cross sections
and the channel's tensor signal is the recoil azimuth's
1 + c₂ cos 2(φ_t − φ_S), so no g₁ is evaluated anywhere in a coherent run —
not in the rate, not in a column, not in `RcModel` (`applies()` is false
there). Measured, `--channel coherent --events 400 --seed 11`, under
helicity-flip, tensor-thirds and transverse-tensor:

| quantity | `--pol-sf toy` vs `nnpdfpol` |
|---|---|
| `sigma_pb()` | equal to the last bit (12476.025113 pb at `toy` unpol) |
| `sigma_per_category_pb()` | equal, element by element |
| all 47 generated array columns | `array_equal(..., equal_nan=True)` — **0 differ** |
| `--unpol-sf toy` → `ct18nlo`, same channel | 12476.025113 → **8806.096207 pb, ×0.705841** |

The last row is why this is **labelled and not refused**: the *unpolarised*
half of the pair does reach this rate, so refusing the other half would break
`--unpol-sf X --pol-sf Y` on one channel of a three-channel scan and nowhere
else. Instead, on that channel:

* `meta["pol_sf"]` carries **`not read on channel coherent-6Li`**
  (`pol_sf_unread_label`) and not a backend name — and that is what a
  `--pol-sf toy` coherent run records too, because `ToyG1` did not run there
  either. This is `rank2_none_label`'s rule for the polarised selector.
* `meta["pol_sf_reach"]` carries the sentence (`pol_sf_reach_report`), the
  same one the banner prints.
* the banner's header line now reads `--unpol-sf reaches EVERY kernel this run
  builds …; --pol-sf reaches the inclusive and the tagged kernels, and NOT
  this run's`. It used to read `both reach EVERY kernel: inclusive, coherent
  and tagged`.

The kernel still *carries* the selected `g1_model` there — the coherent branch
builds the ordinary inclusive kernel — so `dis_sampler.kernel.tables(x, q2).g1`
does move. The predicate is about the **run**, not the object.

#### (2) The neutron-tagged pair is refused at the top of the window

The A_∥ measurement on a free neutron — `--isotope d --channel tagged-d-p
--unpol-sf ct18nlo --pol-sf nnpdfpol --plan helicity-flip --pz 0.7 --pe 0.7`,
at this CLI's own default `--pe` — is refused at configuration time, before an
event is drawn. Measured in the cell the sampler names, **x = 0.954993,
Q² = 1119.1**, on the `tagged-d-p` DIS target:

| backend pair | F₁ⁿ | g₁ⁿ | A₁ = g₁/F₁ | 1 + w_avg at P_z = P_e = 0.7 |
|---|---|---|---|---|
| `toy` + `toy` | 1.37602e−5 | 9.94559e−6 | 0.722778 | builds |
| `ct18nlo` + `toy` | 1.19311e−6 | 8.62353e−7 | 0.722778 | builds |
| `toy` + `nnpdfpol` | 1.37602e−5 | 5.25659e−6 | 0.382014 | builds |
| **`ct18nlo` + `nnpdfpol`** | 1.19311e−6 | 5.25659e−6 | **4.40579** | **−0.02414** |
| **`mstw` + `nnpdfpol`** | 1.10033e−6 | 5.25659e−6 | **4.77729** | **−0.1105** |

`ToyG1` is A₁(x)·F₁ on the kernel's *own* F₁, which is why its A₁ is the same
0.722778 under either F₂ backend: the ratio only runs away when a **grid** g₁
is divided by a **different** grid's F₁. That is a property of the PAIR, and
it is what makes it worth a message rather than a footnote.

**No clamp.** A positivity clamp on g₁ or on the density is a physics change
and would be a silent one; the sampler's own accept–reject would then draw
max(W, 0), diluting the modulation *and* skewing the (x, Q²) mixture. What
changed instead:

* `InclusiveSampler`'s existing configuration-time check — it already
  evaluates 1 + w_avg over **every** accepted cell before an event is drawn —
  now names the cell, that cell's A₁ with its F₁ and g₁, the spin state, and
  **the cure**: lower `Scenario::x_max` (`--x-max`, e.g. `--x-max 0.95`, the
  shipped grid's top cell being x = 0.955) or reduce P_e / P_z, which scale
  w_avg linearly;
* `lipolgen-run` catches it and exits 1 with that message as a **configuration
  error, not a Python traceback** — the `PythiaBridge` clause's rule;
* `--x-max`, `--unpol-sf` and `--pol-sf` all say so in `--help`, and
  `docs/USAGE.md` §2b carries the table above.

Verified cures, same command: `--x-max 0.95` builds; `--pe 0` builds; either
backend alone builds; `tensor-thirds` (P_z = 0) always built.

This is the **vector-sector** twin of the tensor-sector case §D2.1(a) records
for `--inclusive-b1` + `--unpol-sf ct18nlo` (1 + w_avg = −0.1302 at x = 0.955,
Q² = 167.3). Same window edge, same refusal, different sector and channel.

#### (3) The below-grid fraction is the channel's OWN rate

`meta["unpol_sf_below_grid_frac"]` and the banner used to take the **inclusive
sampler's** cell cross sections as the denominator on every channel. On the
coherent channel that is not the run's rate: `f_coh(x)` falls by a factor ≈26
across the window, so it is a large reweighting of the same cells.

| denominator (`--unpol-sf ct18nlo`, config 1) | Σ [pb] | below Q² = 1.677025 [pb] | fraction |
|---|---|---|---|
| inclusive cells (3051 accepted) | 472571.917912 | 170972.122640 | **0.361791** |
| coherent rate (1726 cells carry it) | **8806.096207** | 3940.376757 | **0.447460** |

The coherent sum is `Pipeline::sigma_pb()` **exactly**, which is what makes it
the rate and not an auxiliary weight. The fix is one accessor,
`Pipeline::cell_rate_weights_pb()`: the sampler's own `cell_xsec_pb()` on the
inclusive and tagged channels (so those numbers, the pinned 0.361791 included,
are unmoved) and σ_cell·f_coh(x) with the diffractive-mass gate on the
coherent one. Stated residual: the spin modulation (1 + w_avg) is per (cell,
spin state) and is in **none** of these weights on any channel — they are the
spin-blind rate, the only per-cell weight a run has that does not depend on
which category is being asked about.

#### (4) A caller-supplied kernel on a tagged channel

`Pipeline` reads `cfg.kernel` on the **ion-level** branch only. A tagged run
draws from the struck-cluster sampler and never looks at it — while `meta`
wrote `"caller-supplied kernel"` into `unpol_sf`, `pol_sf` **and** `b1_model`,
with NaN grid keys. Measured before the fix: a `tagged-6Li-alpha` config whose
`kernel` was an `InclusiveKernel` on CT18NLO passed `validate()` and ran at

| F₂A(0.3, 10), deuteron | value |
|---|---|
| what the run used (`ToyF2`) | **0.3691149345** |
| what `meta` implied (CT18NLO) | 0.4639703442 |

`validate()` now refuses it, naming the channel and the way in
(`struck.f2_source` / `struck.g1_model`, which the two selectors fill). **The
refusal is on the tagged channels, not "off the inclusive one":** the coherent
channel rides the ion-level sampler, so a hand-built kernel *does* reach its
rate there and stays accepted — the research note's one-clause phrasing would
have removed a working path.

#### (5) A hand-built T2 bridge on a different backend

The documented pybind path — `PythiaBridge(beams, PythiaBridgeOptions())` then
`set_pythia_hadronizer` — used to accept a bridge whose species draw was on
`ToyF2` while the config said `ct18nlo`, with no `meta` key recording the
bridge's backend at all. `set_pythia_hadronizer` now refuses any bridge whose
`options.f2_source` is not the config's own `unpol_sf_obj` (object
**identity**, the rule the b₁ convolution and `BreakupOptions::f2` already
follow; both unset is the same object, so the default run is untouched). The
price it protects is §D2.2's: the toy's F₂ⁿ/F₂ᵖ is up to 24 % from CT18NLO's,
on the species of every T2 event.

#### What was added

* `pol_sf_is_read` / `pol_sf_reach_report` / `pol_sf_unread_label`
  (`pipeline.hpp`, `pipeline.cpp`), bound to Python and read by
  `meta["pol_sf"]`, the new `meta["pol_sf_reach"]` and the CLI banner — one
  definition, three surfaces.
* `Pipeline::cell_rate_weights_pb()` and the `coh_cell_pb_` it hands out;
  `unpol_sf_grid_report` takes its denominator from it.
* Two refusals: `validate()` on a caller kernel off the ion-level channels,
  and `set_pythia_hadronizer` on a bridge that is not on the run's backend.
* `InclusiveSampler`'s negative-density message: the cell, its A₁/F₁/g₁, the
  spin state and the cure; and a CLI that turns every `Pipeline`-construction
  `RuntimeError` into a `SystemExit` carrying that text.
* 22 pytests (`python/tests/test_sf_backend.py` G9–G13), including both pinned
  fractions, the bit-identity of the coherent columns, the refusal messages by
  content, and the CLI's exit path.

#### Reproduction

The four probes are one file: `python/tests/test_sf_backend.py`, gates G9–G13,
which carry the numbers of this section as assertions rather than as prose. At
the command line the two user-visible ones are

```bash
cd /path/to/LiPolGen && source env.sh
# (2): refused, with the cell and the cure in the message; exit 1
python -m lipolgen.cli --isotope d --channel tagged-d-p --events 50 \
  --unpol-sf ct18nlo --pol-sf nnpdfpol --plan helicity-flip --pz 0.7 --pe 0.7
# ... and cured by the window the message names
python -m lipolgen.cli --isotope d --channel tagged-d-p --events 50 \
  --unpol-sf ct18nlo --pol-sf nnpdfpol --plan helicity-flip --pz 0.7 --pe 0.7 \
  --x-max 0.95 --quiet
# (1) and (3): the banner says --pol-sf did not run here, and prints 44.75 %
python -m lipolgen.cli --channel coherent --events 50 \
  --unpol-sf ct18nlo --pol-sf nnpdfpol
```

and the pinned values are `_l.unpol_sf_grid_report(p)` = (1.677025, 0.447460)
on that coherent pipeline against (1.677025, 0.361791) on the inclusive one.

#### What this pass did NOT do, and is therefore not claimed

* **No cross-section factor is quoted for the d+p or the ⁷Li α-tag channel.**
  The d+p tag was run; on the shipped window the interesting pair has no run
  to quote, and `--x-max 0.95` is a *different window*, so its σ is not
  comparable with §D2.1's table.
* **The x_max = 0.95 window was not characterised.** What it costs in rate,
  and how much of §D2.1's A_∥ lives in the cells it removes, is unmeasured.
* **The spin modulation is still not in the below-grid denominator** on any
  channel, and there is no per-category below-grid fraction. Stated above,
  not fixed.
* **The `pol_sf` label is per channel, not per plan.** On the inclusive and
  tagged channels the key records the backend name whether or not a given fill
  has polarised sensitivity — the kernel reads g₁ there either way, so nothing
  false is written, but the key is not a statement about the plan.
* **Nothing new was measured about `MstwSF`'s grid floor**, which still
  reports none and is still written NaN (§D2.7 stands).

---

## D3 — the FSI "variants": the premise is FALSE, and the real defect was the `meta`

The research note is `phase_D_small_items.md` §D3. Its headline reproduces
here exactly, so this section is a confirmation and not a re-derivation.

### D3.1 The claim that was retracted

`PLAN.md` carried, as open item D3: *"The per-nucleon Glauber FSI variant is
algebraically identical to the cluster one — the two 'variants' are one."*
**Retracted.** It confused two different statements:

* **True, and already in `fsi.hpp`:** for an *uncorrelated* cluster density the
  Ciofi degli Atti–Kaptari per-nucleon product ⟨Π_i[1 − Γ_N(b − s_i)]⟩
  factorises into [1 − (Γ_N ⊛ T_a)]^A — the **cluster** form. So a literal
  per-nucleon code path *would* reproduce variant (a).
* **False:** that the shipped `GlauberNucleon` *is* such a path. It is not
  (`src/core/fsi.cpp`, `gamma_profile`): it returns A·(Γ_N ⊛ T_a), the
  **single-scattering (optical, unshadowed) limit**, σ_Xa = A σ_XN exactly —
  the first term of the same binomial, chosen deliberately because the full
  product would print the cluster variant's numbers.

### D3.2 Measured — the two profiles

Both on the production S+D ⁶Li α tag (`li6_alpha_channel()`), σ_XN = 40 mb,
ε = −0.5, B_XN = 6 GeV⁻², everything else default. Read straight out of the
run's own `meta` (see §D3.4):

| quantity | `GlauberCluster` | `GlauberNucleon` |
|---|---|---|
| σ_tot(X–α) [mb] | **131.045** | **160.006** |
| σ_el(X–α) [mb] | 35.224 | 68.186 |
| σ_el/σ_tot | 0.269 | 0.426 |
| `survival()` | **0.520239** | **0.582899** |

131.045 and 35.224 reproduce `fsi.hpp`'s own header row exactly, and 160.006 is
A·σ_XN = 4 × 40 to 4e−5 relative — the unshadowed bound.

### D3.3 Measured — the two weights, on the same events

```bash
source env.sh
# make_config(channel="tagged-6Li-alpha", events=20000, seed=1234, fsi=V)
# Pipeline(cfg, make_plan("tensor-thirds", j=1, pz=0.7, pzz=0.6, pe=0.7))
#   .generate(0, False, 1)   for V in off / glauber-cluster / glauber-nucleon
```

The columns `k`, `cos_theta_k`, `phi_k`, `x`, `q2` are **bit-identical**
(`np.array_equal`) across the three runs, so this is one event stream
reweighted three ways and the weight column is the whole of the difference:

| | `off` | `glauber-cluster` | `glauber-nucleon` |
|---|---|---|---|
| Σ weight (20 000 events) | 20000.0000 | 10419.0712 | 11632.0979 |
| integrated survival Σw/Σw_off | 1 | **0.520954** | **0.581605** |
| per-event factor, min … max | 1 | 0.0207504 … 1.40428 | 0.213898 … 6.9523 |

Per-event ratio w_nucleon/w_cluster, percentiles [1, 5, 25, 50, 75, 95, 99] =
**0.8208, 0.8265, 0.8527, 0.8974, 0.9338, 3.3374, 5.9347**; **maximum 68.5222**,
minimum 0.8178; `np.allclose(rtol = 1e-14)` → **False**; **99.5000 %** of the
20 000 events differ by more than 1 %.

**Verdict: different on essentially every event, and by construction rather
than by parameter choice.** The identity that does hold is an algebraic
factorisation, not a coincidence of σ_XN = 40 mb, so no choice of parameters
makes the two shipped variants agree.

### D3.4 The real defect — FSI was absent from the npz `meta`

Before this task the three runs of §D3.3 produced **`meta` dicts that compared
equal, key by key**, while their total rate differed by 48 % / 42 % and
individual event weights by a factor up to 68.5. `python/bindings.cpp` states
the rule on the `b1_*` block — *without these keys three otherwise identical
npz files are indistinguishable* — and it had been applied to `b1_*` and to
`rc_*` and to nothing else.

Fixed, following the **`rc` precedent** (a conditional block, so an `--fsi off`
npz keeps exactly the old key set and no reference gate moves). Emitted when
`Pipeline::fsi_weight()` is non-null:

| key | `glauber-cluster` | `glauber-nucleon` |
|---|---|---|
| `fsi` | `"glauber-cluster"` | `"glauber-nucleon"` |
| `fsi_sigma_mb` | 40.0 | 40.0 |
| `fsi_sigma_cluster_mb` | 131.04542080791026 | 160.00623096559633 |
| `fsi_sigma_cluster_el_mb` | 35.223684907167275 | 68.18646124641093 |
| `fsi_survival` | 0.5202393836689237 | 0.5828989055760844 |
| `fsi_clipped_grid_fraction` | 0.0 | 0.0 |
| `fsi_formation_ramp` | False | False |

plus `fsi_ramp_w_lo` / `fsi_ramp_sigma_lo_mb` / `fsi_ramp_w_hi` /
`fsi_ramp_sigma_hi_mb` **only when the ramp ran** — a σ_XN(W) ramp is a
different run and `fsi_sigma_mb` alone does not describe it. This matters most
for `fsi_sigma_mb`: the header's mandatory 20–40 mb band was unquotable while
its two ends produced indistinguishable files.

### D3.5 What a genuinely different per-nucleon variant would need

`fsi.hpp`'s standing TODO names the two corrections and they are the right two.
Of them:

* **the centre-of-mass constraint Σ_i s_i = 0** is a change of T_a alone, i.e.
  of `cluster_point_a2_fm2` — Gartenhaus–Schwartz on a Gaussian gives
  a²_int = a²(1 − 1/A) = 0.7001 × 0.75 = **0.5250 fm²** for the α. A better
  *input to the same variant*, not a different variant. (The σ_tot / σ_el
  consequences of that substitution were **not measured here** — `a2_fm2_` is
  derived from the spectator (Z, A) and has no injection point — so the
  research note's analytic estimate is not repeated as a result of this task.)
* **short-range NN correlations** need the **two-body** density
  ρ₂(s_i, s_j) — that is the whole content of ⟨Π_i(·)⟩ ≠ Π_i⟨·⟩ — and **the
  tree has none**. `data/vmc/density/{he4,li6}.density` are ONE-body
  point-proton densities, which is exactly the T_a the cluster form already
  convolves, and the ANL page `data/vmc/density/README.md` fetches from
  publishes one-body densities only. **The SRC correction cannot be built from
  committed data.**

And the input is **not** σ_NN: the rescattering is the DIS debris X off a
*spectator nucleon*, so the amplitude is X–N and the cross section is
`GlauberFsiOptions::sigma_xn_mb`. σ_NN would enter only for the spectator
cluster's own internal absorption, which is not what this weight is.

### D3.6 What was added

* `include/lipolgen/fsi.hpp` — the retraction, the §D3.3 table, and the TODO
  extended with the two-body-density gap and the σ_NN correction.
* `tests/test_fsi.cpp`, *"the two variants differ on essentially every
  event"* — 4 000 generated events, both variants and `off` on one stream:
  bit-identical kinematics, the two profiles, Σw ratios, ratio range and the
  >98 % differing fraction (measured 99.375 % at 4 000, ratio [0.8179, 68.52]).
* `python/tests/test_meta_provenance.py` — the same statement through the npz,
  plus the three-way `meta` distinguishability and the 20 vs 40 mb band.

---

## D4 — `PDF:PomSet`: scanned, and the proposed observable was the wrong one

### D4.1 The scan, as run

```bash
source env.sh
# per set N in 1..15:
#   cfg = make_config(isotope="6Li", config=1, channel="coherent",
#                     events=20000, seed=4242)
#   po = PythiaBridgeOptions(); po.coherent_t2 = CoherentT2.Pomeron
#   po.pom_set = N;  bridge = PythiaBridge(default_configs("6Li")[1], po)
#   set_pythia_hadronizer(cfg, bridge)
#   cols = Pipeline(cfg, make_plan("tensor-thirds", j=1, pz=0.7, pzz=0.6,
#                                  pe=0.7)).generate(0, True, 1)
#   a = hfs_arrays(cols["events"], False, False)   # offsets/pid/p4/charge
```

**17.1 s wall for all fifteen sets at 20 000 events each**, single thread,
as the scan was run on 2026-09-04. Re-run 2026-09-05, after §D4.6's
constructor guard: **14.3 s for the fourteen sets that still run** (1–10,
12–15); set 11 now throws before the first event. Cost is not a
consideration. Sets are PYTHIA 8.317's, read off
`xmldoc/PDFSelection.xml` rather than remembered.

### D4.2 The T0 side does not move — so the proposed M_X band is identically zero

The md5 of `t ‖ x_pom ‖ q2 ‖ x ‖ weight`, over all five columns of all 20 000
events, is `ffd35a3a62b591c547e9ca2ac4301b5d` for **every one of the fourteen
sets that run today** — 1–10 and 12–15, re-measured 2026-09-05, one md5 over
the lot. It was that same md5 on set 11 when the scan was made on 2026-09-04,
but §D4.6's constructor guard has since put set 11 out of reach, so **the
reproducible claim is over fourteen, and "all 15 sets" is what the surface
should no longer say** (corrected 2026-09-05 at fourteen sites:
`python/lipolgen/cli.py` help and banner, `docs/USAGE.md` §4,
`docs/OPEN_ITEMS_SOLUTIONS.md` §5.1 and its item-5 summary row,
`docs/DEVELOPMENT_PLAN.md`, `docs/PYTHIA_BRIDGE.md` §11, `docs/T2_CHAIN.md`,
`docs/PHYSICS_CHANNELS.md` row 459, `include/lipolgen/pythia_bridge.hpp`,
`python/bindings.cpp`, `PLAN.md` D4, `STATUS.md` row D and
`phase_D_small_items.md` §D4.3 — see §D7 F).

The Pomeron PDF enters only the flavour draw
(`src/pythia/pythia_bridge.cpp`) and PYTHIA's backward evolution; |t|, x_P,
M_X and the rate are fixed upstream by `CoherentSampler` /
`CoherentXpomModel`. ⟨|t|⟩ = 0.020054 and ⟨x_P⟩ = 0.034111 on every set.

**A `PomSet` band on M_X, |t|, x_P or σ is identically zero by construction.**
The item as filed proposed M_X as the observable; quoting a band on it would
have been quoting the width of a constant. The band is on the **hadronic final
state**.

### D4.3 The full scan, 20 000 events per set, seed 4242

`n_ch` = charged final-state particles per event; `n_had` = all final-state
particles per event; ⟨p_T⟩ = mean transverse momentum per final-state particle;
fractions are of all final-state particles. `fallback` is
`n_pom_flavour_fallback / n_pomeron`. Errors on ⟨n_ch⟩ are √N on the mean.

| set | fit | n_ok | fail | fallback | ⟨n_ch⟩ | ⟨n_had⟩ | ⟨p_T⟩ [GeV] | π | K | p | γ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | toy N x^a(1-x)^b, Q²-indep. | 20000 | 0 | 0.00 % | 3.8559 ± 0.0150 | 8.367 | 0.3327 | 0.408 | 0.0511 | 0.0109 | 0.514 |
| 2 | π⁰ densities | 20000 | 0 | 0.00 % | 3.7773 ± 0.0148 | 7.926 | 0.3434 | 0.427 | 0.0462 | 0.0127 | 0.494 |
| 3 | H1 2006 Fit A NLO | 20000 | 0 | 30.80 % | 3.9223 ± 0.0156 | 8.391 | 0.3428 | 0.405 | 0.0653 | 0.0104 | 0.504 |
| 4 | H1 2006 Fit B NLO | 20000 | 0 | 35.00 % | 3.9846 ± 0.0165 | 8.509 | 0.3560 | 0.408 | 0.0646 | 0.0096 | 0.504 |
| 5 | H1 2007 Jets NLO | 19862 | 138 | 0.00 % | 4.3554 ± 0.0193 | 9.394 | 0.3220 | 0.384 | 0.0834 | 0.0056 | 0.494 |
| 6 | **H1 2006 Fit B LO (default)** | 20000 | 0 | 19.99 % | 3.9639 ± 0.0163 | 8.521 | 0.3420 | 0.404 | 0.0646 | 0.0103 | 0.507 |
| 7 | ACTW B NLO ε=0.14 | 20000 | 0 | 0.00 % | 3.8975 ± 0.0156 | 8.517 | 0.3401 | 0.413 | 0.0379 | 0.0110 | 0.522 |
| 8 | ACTW D NLO ε=0.14 | 20000 | 0 | 0.00 % | 3.9693 ± 0.0167 | 8.325 | 0.3763 | 0.398 | 0.0958 | 0.0081 | 0.486 |
| 9 | ACTW SG NLO ε=0.14 | 20000 | 0 | 0.00 % | 3.8600 ± 0.0151 | 8.452 | 0.3319 | 0.414 | 0.0359 | 0.0111 | 0.523 |
| 10 | ACTW D NLO ε=0.19 | 19999 | 1 | 0.00 % | 3.9601 ± 0.0166 | 8.300 | 0.3793 | 0.396 | 0.1002 | 0.0078 | 0.483 |
| 11 | PomHISASD (Angantyr) | 20000 | 0 | 100.00 % | 3.8377 ± 0.0148 | 8.252 | 0.3316 | 0.404 | 0.0637 | 0.0108 | 0.505 |
| 12 | GKG18 Fit A LO | 20000 | 0 | 18.52 % | 4.0356 ± 0.0168 | 8.648 | 0.3372 | 0.402 | 0.0671 | 0.0100 | 0.502 |
| 13 | GKG18 Fit B LO | 20000 | 0 | 19.28 % | 4.0272 ± 0.0168 | 8.661 | 0.3369 | 0.401 | 0.0668 | 0.0097 | 0.504 |
| 14 | GKG18 Fit A NLO | 20000 | 0 | 19.94 % | 4.0159 ± 0.0165 | 8.581 | 0.3366 | 0.403 | 0.0677 | 0.0098 | 0.501 |
| 15 | GKG18 Fit B NLO | 20000 | 0 | 20.83 % | 4.0105 ± 0.0166 | 8.580 | 0.3366 | 0.403 | 0.0679 | 0.0097 | 0.502 |

Set 11's row was measured **before** the guard of §D4.6 was added; that
configuration now throws, so this one row of the table cannot be reproduced
and is kept only as the record of what was seen. Everything quoted as a band
is over the twelve fits below, which never included it.

### D4.4 The band, and how to quote it

Over the **twelve genuine diffractive-PDF fits** — 3, 4, 5, 6, 7, 8, 9, 10, 12,
13, 14, 15 — about the default set 6:

| observable | set 6 | min | max | band about the default | max/min |
|---|---|---|---|---|---|
| **⟨n_charged⟩** (primary) | 3.9639 ± 0.0163 | 3.8600 (set 9) | 4.3554 (set 5) | **−2.62 % / +9.88 %** | 1.128 |
| ⟨n_hadrons⟩ | 8.5209 | 8.2997 (10) | 9.3938 (5) | −2.60 % / +10.24 % | 1.132 |
| ⟨p_T⟩ per particle [GeV] | 0.34199 | 0.32196 (5) | 0.37925 (10) | −5.86 % / +10.90 % | 1.178 |
| ⟨p_T⟩ per **charged** particle | 0.45405 | 0.42365 (5) | 0.49533 (10) | −6.70 % / +9.09 % | 1.169 |
| **kaon fraction of the HFS** (secondary) | 0.06456 | 0.03592 (9) | 0.10020 (10) | **−44.36 % / +55.20 %** | **2.790** |
| pion fraction | 0.40380 | 0.38402 (5) | 0.41367 (9) | −4.90 % / +2.45 % | 1.077 |
| proton fraction | 0.01030 | 0.00562 (5) | 0.01114 (9) | −45.50 % / +8.10 % | 1.984 |
| photon fraction | 0.50681 | 0.48334 (10) | 0.52270 (9) | −4.63 % / +3.14 % | 1.081 |

**Seed stability**, 20 000 events at seeds 4242 / 777 / 31337 — ⟨n_ch⟩:

| set | 4242 | 777 | 31337 | spread |
|---|---|---|---|---|
| 1 | 3.8559 | 3.8870 | 3.9011 | 0.0452 |
| 2 | 3.7773 | 3.8224 | 3.8342 | 0.0569 |
| 3 | 3.9223 | 3.9469 | 3.9249 | 0.0246 |
| 4 | 3.9846 | 3.9842 | 3.9796 | 0.0050 |
| 5 | 4.3554 | 4.3740 | 4.3712 | 0.0186 |
| 6 | 3.9639 | 4.0030 | 3.9865 | 0.0391 |
| 7 | 3.8975 | 3.9170 | 3.9351 | 0.0376 |
| 8 | 3.9693 | 3.9523 | 3.9928 | 0.0405 |
| 9 | 3.8600 | 3.9012 | 3.9147 | 0.0547 |
| 10 | 3.9601 | 3.9617 | 3.9716 | 0.0115 |
| 11 | 3.8377 | 3.8749 | 3.8776 | 0.0399 |
| 12 | 4.0356 | 4.0600 | 4.0694 | 0.0338 |
| 13 | 4.0272 | 4.0560 | 4.0762 | 0.0490 |
| 14 | 4.0159 | 4.0172 | 4.0304 | 0.0145 |
| 15 | 4.0105 | 4.0192 | 4.0363 | 0.0258 |

Seed-to-seed scatter ≤ **0.0569** against a **0.4954** spread across the fit
sets, and the band per seed is −2.62/−2.54/−1.80 % on the low side and
+9.88/+9.27/+9.65 % on the high side, with **set 9 always the minimum and set 5
always the maximum**; the kaon ratio is 2.790 / 2.747 / 2.659 with set 9 always
the minimum and set 10 always the maximum. **20 000 events per set already
resolves the band and there is no case for more.**

Four rules travel with the band, and they are in `pythia_bridge.hpp`,
`--pom-set`'s help, the run banner, `USAGE.md` §4 and
`OPEN_ITEMS_SOLUTIONS.md` §5.1:

1. **Envelope over re-runs, one npz per set** — never a per-event reweighting;
   there is no weight that maps one set onto another. `meta["pom_set"]` (§D4.7)
   is what tells the files apart.
2. **Sets 1 and 2 are outside the band** (a Q²-independent toy and π⁰
   densities; neither is a Pomeron fit). Set 1 is still worth running as a
   sanity floor.
3. **Set 6 is the only LO *H1* set**, so the band mixes LO and NLO DPDFs used
   in an LO Monte Carlo: fine for a systematic envelope, not for a central
   value. (The tree said "the only LO Q²-dependent set in the list" — false:
   12 and 13 are LO GKG18. Corrected at `pythia_bridge.hpp` and
   `PHYSICS_CHANNELS.md`.)
4. **Set 5 makes open charm the default never makes** — §D4.5.

### D4.5 Charm — a feature of the band, and narrower than the tree stated

Measured by toggling `PythiaBridgeOptions::include_charm` and diffing the
final states event by event (4 000 events, seed 4242; "changed" = the event's
`pid` list differs):

| set | events changed by `include_charm = false` |
|---|---|
| 5 (H1 2007 Jets) | **29.98 %** |
| 12 / 13 / 14 / 15 (GKG18) | 10.40 / 10.78 / 10.70 / 11.05 % |
| 3, 4, 6, 7, 8, 9, 10 | **0.00 %** |

So set 5 puts nearly 30 % of a coherent sample on a charm initiator, and the
ACTW grids give charm no support in this generator's (β, Q²) window even
though `CTEQ6pdf` carries a charm slot. Flag it in any scan that includes set
5: it is the same signature the `2e404b8` bug faked.

**Confirmed at the PYTHIA source, and the tree's wording was too narrow.**
`PartonDistributions.cc:2630` sets `xc = xcbar = 0.` and `xb = xbbar = 0.`
inside `PomH1FitAB::xfUpdate`, and `BeamSetup.cc:1371-1377` dispatches sets
**3, 4 and 6** to that one class. The tree said "the H1 **LO** grids carry no
charm or bottom at any (β, Q²)" — true, and equally true of the two H1 **NLO**
fits, which the wording did not cover. Corrected at
`include/lipolgen/pythia_bridge.hpp` and `docs/PHYSICS_CHANNELS.md`.

### D4.6 The light-only e_q² fallback is NOT a systematic — measured, it costs zero

4 000 coherent events per point at 6Li config 1, seed 4242, diffing the whole
final state (`pid` and `p4` arrays) wholesale by md5:

| set | q2_pdf_min 1.0 | → 1.75 | final state |
|---|---|---|---|
| 6 | 20.70 % fallback | 0.00 % | **bit-identical** |
| 3 | 31.25 % | 0.00 % | **bit-identical** |
| 4 | 35.52 % | 29.32 % | **bit-identical** |
| 12 | 19.45 % | 2.65 % | **bit-identical** |
| 13 | 20.10 % | 3.23 % | **bit-identical** |
| 15 | 21.65 % | 5.12 % | **bit-identical** |

**Set 4 is the row that must not be dropped when this table is summarised in a
sentence.** It is the only one where 1.75 does *not* collapse the share, and
**six** sites summarised the table as "0.00 % (sets 6, 3) or a few per cent
(12, 13, 15)" — leaving set 4 out. Four were corrected in the first pass on
2026-09-05 (`include/lipolgen/pythia_bridge.hpp` §`n_pom_flavour_fallback`,
`docs/USAGE.md` §4, `docs/PYTHIA_BRIDGE.md` §11, `docs/PHYSICS_CHANNELS.md`
row 459) and this note said "four sites" — but two more, both of them
CODE-doc, survived that pass and were corrected only in the second pass the
same day: the `n_pom_flavour_fallback` binding docstring in
`python/bindings.cpp` and the block comment above the light-only branch in
`src/pythia/pythia_bridge.cpp`, which had read "removes the fallback (or
nearly)" — the parenthesis was the whole of set 4. The exact
counters, re-measured 2026-09-05 on the same recipe, are 1421/4000 = 35.525 %
at 1.0 and 1173/4000 = 29.325 % at 1.75; the two-decimal rendering above
truncates, and it is the rendering all six sites now carry.

The reason is structural rather than lucky: every Pomeron DPDF PYTHIA ships
carries a **single light-quark singlet** — H1 sets it explicitly
(`PartonDistributions.cc:2623-2631`), and the GKG18 LHAGrid1 files have
columns −3 … 3 numerically equal row by row — so e_q²·xf_q ∝ e_q² **exactly**
over the light flavours and the normalised "fallback" *is* the true draw. And
it fires only where every weight vanishes, i.e. below the charm threshold,
where charm is zero anyway.

**The one real cost is the CLAMP, and it is a charm effect, not a fallback
one.** At `q2_pdf_min` = 3.0 sets 12, 13 and 15 stop being bit-identical (sets
3, 4, 6 stay identical) — and with `include_charm = false` the same 1.0 vs 3.0
comparison is **bit-identical again** on all three. So the clamp's only effect
is to lift charm above threshold. This is a **finding of this task that the
research note did not have**: it had checked 1.75 and 3.0 and reported both as
identical.

*Correction applied:* `docs/PHYSICS_CHANNELS.md`'s remedy — *"raising
`q2_pdf_min` to ≈1.75 removes it at the cost of clamping every flavour weight
to that Q²"* — was true about the counter and **false about the cost**. The
light-only restriction stays (it is the correct guard against the `2e404b8`
bug) and the tree no longer implies the fallback share is a modelling
uncertainty.

### D4.7 The `meta` and the run surface

Second instance of the D3.4 defect: two runs differing only in `--pom-set` had
**identical `meta`** while their whole hadronic final state differed.

The mechanics had to be solved first. `PipelineConfig::hadronizer` is a
type-erased `std::function`, so the metadata writer could not see the bridge at
all. `set_pythia_hadronizer` now installs a **named** callable
(`PythiaHadronizerHook`) and `columns_to_dict` recovers the bridge with
`std::function::target`. The block is emitted only when that hook is present,
so a T0 npz keeps exactly the old key set. Measured on three 800-event runs:

| key | set 6 | set 5 | set 11 (before the guard) |
|---|---|---|---|
| `t2_bridge` | `"pythia8"` | `"pythia8"` | `"pythia8"` |
| `coherent_t2` | `"pomeron"` | `"pomeron"` | `"pomeron"` |
| `pom_set` | 6 | 5 | 11 |
| `pom_rescale` | 1.0 | 1.0 | 1.0 |
| `pom_q2_pdf_min` | 1.0 | 1.0 | 1.0 |
| `t2_include_charm` | True | True | True |
| `t2_n_ok` / `t2_n_failed` | 800 / 0 | 794 / 6 | 800 / 0 |
| `t2_n_pomeron` | 800 | 794 | 800 |
| `n_pom_flavour_fallback` | 166 | 0 | 800 |
| `pom_flavour_fallback_frac` | 0.2075 | 0.0 | **1.0** |

`n_pom_flavour_fallback` was counted and **surfaced nowhere** — not in the
`meta`, not in the banner, which printed `n_ok` / `n_failed` / `n_retries`
only. It is now in both.

### D4.8 Set 11 is refused

`PomHISASD` (set 11) returns its densities only after `setXPom(x_Pom)`, which
this bridge never calls — there is no Angantyr collision system here to supply
x_P. Measured, 20 000 events: **`n_pom_flavour_fallback` = 20000 / 20000 =
100.00 %.** The Pomeron PDF is not consulted on a single event; the whole run
is the charge-democratic e_q² fallback, while `meta["pom_set"]` would say 11.
That is exactly the case `PipelineConfig::validate`'s *"a knob that did not run
may not be recorded as if it had"* rule exists to prevent, so
`PythiaBridge`'s constructor now **throws** on `pom_set == 11` when the Pomeron
instance is built (`coherent_t2 == Pomeron`), naming the measurement and the
alternatives. With `coherent_t2 = Off` it is accepted, because then no Pomeron
instance exists and the setting is inert.

### D4.9 What was added

* `include/lipolgen/pythia_bridge.hpp` — the corrected set list, the band, the
  seed stability, the charm finding, and the rewritten
  `n_pom_flavour_fallback` doc with the bit-identity measurement.
* `src/pythia/pythia_bridge.cpp` — the set-11 guard.
* `python/bindings.cpp` — `PythiaHadronizerHook` and the T2 `meta` block.
* `python/lipolgen/cli.py` — `--pom-set`'s help carries the band; the coherent
  run banner prints it; the fallback share is printed beside `n_ok`/`n_failed`;
  and a bridge refusal now comes out as a command-line error carrying the
  bridge's own message rather than as a traceback (the message is not
  restated at the CLI — one definition).
* `python/tests/test_pom_set_band.py` — the reduced re-run: 4 sets × 4 000
  events (3.9 s), pinning the T0 bit-identity, the HFS difference, the ⟨n_ch⟩
  ordering and spread, the kaon factor, the `meta` distinguishability, the
  fallback bit-identity and the set-11 refusal.

---

## D5 — the |t| ceiling: 0.2 stays, and the reason it is written down changed

### D5.1 The arithmetic, reproduced

With c₂(|t|, P_zz) = A|t| + C, A = −(P_zz/2)·eps_b0·B and C = amp·P_zz
(`src/core/coherent.cpp`, `cos2phi_coefficient`), the positivity edge |c₂| = 1
has the closed form

    |t|_pos = (1 − sign(A)·C) / |A|

and for the shipped signs (eps_b0 < 0, amp > 0, so A and C share a sign) this
is 2(1/|P_zz| − amp)/(|eps_b0|·B). Now implemented as
`CoherentScenario::t_positivity_edge(pzz)` and measured through it, at
B = 50, amp = 0.01:

| eps_b0 | \|t\|_pos, P_zz = −2 | \|t\|_pos, P_zz = +1 |
|---|---|---|
| −0.08 (shipped, an EXACT input) | **0.2450** | **0.4950** |
| −0.0527 (α+d model Q, rounded) | 0.3719 | 0.7514 |
| −0.0171 (GFMC Q, rounded) | 1.1462 | 2.3158 |
| −0.0070 (**measured** Q, rounded) | **2.80** | 5.6571 |

which reproduces `phase_C_numbers.md` §C4.5's table exactly, and confirms its
conclusion: **the positivity edge at 0.245 is a consequence of an `eps_b0`
that is 11.4× the measured ⁶Li quadrupole, not a property of ⁶Li.**

**Read the digit count with the table** (this section printed a four-figure
**2.800** with no qualification until 2026-09-05, which is neither of the two
numbers there are). The last three `eps_b0` are ROUNDED and every row is
arithmetic on the rounded value, so the bottom row is quotable as
**2.80 GeV²** — three figures, all a two-figure input supports — and never as
"2.8000". On the DERIVED band (`ClusterConfigSampler::quadrupole_band_fm2()`
through `a2_from_quadrupole`, at B = 50 and amp = 0.01; measured 2026-09-05)
the three `eps_b0` are −0.0526846 / −0.0171207 / **−0.0070024** and the
P_zz = −2 edges are 0.372025 / 1.144810 / **2.799047**, i.e. the same
**2.80**. `tests/test_coherent.cpp` asserts the derived 2.7990 and 0.37202,
which is where that digit comes from.

### D5.2 The decision, and why it is not Option A

The header carried **two** reasons for `COHERENT_T_MAX_DEFAULT` = 0.2 and said
*"Two reasons, and they agree"*. They agree only at the shipped `eps_b0`:

* the **anchor range** — `mantysaari_a2_deuteron()` carries four digitized
  rows, |t| = 0.05, 0.10, 0.20, 0.30, and the fit is linear in |t| and exact
  only as |t| → 0 — is a property of the **input table** and does not move
  with any knob;
* the **positivity edge** moves with `eps_b0`, and at the measured quadrupole
  it is **2.80 GeV²** (2.8000 on the rounded `eps_b0` = −0.0070, 2.7990 on the
  derived −0.0070024 — never "2.8000" as a published figure), **9.3× outside**
  that same |t| ≤ 0.30.

**So deriving the ceiling from positivity would be less honest than the fixed
number, not more**: it would license extrapolating a linear-in-|t| fit ten
times past its data. The ceiling therefore **stays fixed at 0.2**, the anchor
range becomes the stated primary reason, and positivity is stated as secondary
and contingent — with `t_positivity_edge` deriving the contingent number so it
is not retyped in four docstrings and cannot drift from `eps_b0`.

### D5.3 What the ceiling costs — measured

200 000 generated coherent events at the shipped defaults
(`make_config(channel="coherent", events=200000, seed=99)`, `tensor-thirds`
at pz = 0.7 / pzz = 0.6, 4 threads):

    ⟨|t|⟩    = 0.019963 GeV²        max |t| = 0.197575  (the ceiling)
    frac |t| > 0.05 = 0.082015      (untruncated e^{−B|t|} = 8.2085e−02)
    frac |t| > 0.10 = 0.006550      (6.7379e−03)
    frac |t| > 0.15 = 0.000455      (5.5308e−04)
    frac |t| > 0.20 = 0             (4.5400e−05)

`sample_t` **renormalises** on [0, t_max], so the ceiling does not lose rate —
it redistributes **4.54e−5** of it. Moving it to 0.245, or anywhere else
inside the anchor range, changes the generated sample at the 1e−5 level.
**The ceiling is not a rate question and never was**; it is a statement about
where the model is defined, and it is now justified as one.

### D5.4 What was implemented

1. **`COHERENT_T_MAX_DEFAULT` = 0.2, unchanged.** No reference JSON moved.
2. Its comment rewritten (`include/lipolgen/coherent.hpp`), with the primary /
   secondary split, the §D5.3 rate number and the dependency on the `eps_b0`
   decision stated as *contingent only*. The same correction applied at
   `PipelineConfig::coherent_t_max`, `CoherentScenario::eps_b0`,
   `CoherentScenario::positivity_margin`, `docs/CONVENTIONS.md`,
   `docs/PHYSICS_CHANNELS.md` (two sites) and `docs/USAGE.md`.
3. **`positivity_margin` and `CoherentSampler::check_positivity` unchanged and
   still throwing.** They are the GUARD on the truncated weight the sampler
   uses — an author who raises `t_max` or `eps_b0` past where it stops being a
   density is caught by nothing else — not the derivation of the constant. The
   header now says which job each does.
4. **`CoherentScenario::t_positivity_edge(pzz)` added** and bound to Python
   (as is `positivity_margin`, which had never been bound). Changes no default.
   `tests/test_coherent.cpp` gates it over 96 (eps_b0, B, amp, P_zz)
   combinations against the zero of `positivity_margin` — including the
   opposite-sign case the shipped scenario never reaches — plus the degenerate
   cases (A = 0 → +∞; |C| ≥ 1 → 0). The pre-existing hand-rolled closed form in
   that file's item (7) is replaced by a call, so the edge is defined once.
5. **`--coherent-t-max` added** (`make_config(coherent_t_max=…)` too) and the
   whole coherent block recorded in the npz `meta`: `coherent_t_max`,
   `coherent_slope_b`, `coherent_eps_b0`, `coherent_amp`, `coherent_f0`,
   `coherent_m_x_min`, `coherent_x_pom_max`, `coherent_weighted_azimuth` and
   the derived `coherent_t_positivity_edge_pzz_m2` = 0.245. Conditional on the
   channel, so no other npz moves. Before this, a Python caller who moved
   `coherent_t_max` changed the entire |t| spectrum, the tag acceptance and
   every c₂ in the file and the sidecar could not tell.
6. **The run banner** prints the ceiling with its reason and the derived edge
   at both P_zz signs, on every coherent run.
7. **STATUS.md decision row 12** records the decision, and row 8's cost clause
   — which said the ceiling *was* a consequence of `eps_b0` — is amended.

### D5.5 What was NOT done, and is therefore not claimed

* **The modulation was not un-truncated.** Sampling φ from
  Σ_m p_m e^{−ΔB_m|t| cos 2φ} directly, instead of from its O(ΔB|t|)
  truncation 1 + c₂ cos 2φ, is the physically right form — it is positive for
  every |t| by construction — but it moves **every pinned coherent azimuth**,
  so it is filed as the follow-on to the `eps_b0` decision (STATUS.md row 8),
  not taken here.
* **The exact-vs-truncated comparison of `phase_D_small_items.md` §D5.2 was
  not re-run**, so none of its numbers (the exact edges, the 14.30 % / 0.12 %
  truncation errors) entered a header, a constant or a test. They are the
  research note's and are cited as such.
* **`eps_b0` did not move.** That is STATUS.md decision row 8 and remains the
  author's.

---

## D6 — the knob-provenance table: one mechanism for "a knob that did not run"

**Third pass, 2026-09-05.** Every number in this section was measured on that
date against this tree, and every reach cell in the table below is printed by
the code itself (`Pipeline::knob_provenance`, `include/lipolgen/pipeline.hpp`)
rather than transcribed.

### D6.1 Why a mechanism, and not a sixth fix

The rule — *a knob that did not run may not be recorded in the `meta` or
printed in the banner as if it had* — was enforced knob by knob:
`PipelineConfig::validate()` refuses `b1_band_scale` on the Miller branch and
`rc_qe_tensor_scale` without the quasi-elastic tail; `rank2_input_report`
labels the 7Li rank-2 zero; `pol_sf_is_read` labelled the coherent channel
(§D2.8). Each of those closed **one cell of a two-dimensional table nobody had
written down**, and the rule then broke five times in this one run, each time
on an axis the previous fix had not looked at:

| round | axis | what was recorded as if it had run |
|---|---|---|
| 1 | the **channel** | `--pol-sf` on `coherent-6Li` (closed §D2.8) |
| 2 | the **run plan** | `--pol-sf` under every unpolarised-beam plan — including `tensor-thirds`, the CLI's own default — and `--pe` itself, which the three tensor plans never read |
| 3 | the **T2 tier** | `--pom-set` / `--pom-rescale` / `--coherent-t2` on the six non-coherent channels, where the Pomeron PYTHIA instance is built and hadronises zero events |
| 4 | the **rc sub-knobs** | `--rc-fq-scale` and five more on the tagged channels, where `rc_tail_applies` is false by construction; every rc knob on coherent, where `rc_applies` is |
| 5 | **silence** | `--cluster-wave`, `--triton-sf`, `--inclusive-b1`, `--cluster-beta`, `--p-d`, `--fsi-sigma-mb` at `--fsi off`, `--coherent-t-max` off the coherent channel — accepted, some of them **moving the output**, recorded nowhere at all |

The evidence is the adversarial matrix of 2026-09-05: 12 channel × plan specs ×
up to 38 knob cells = 336 cells, 600 events at seed 7, sha256 over all 47
ndarray columns plus `sigma_pb` and `sigma_per_category_pb` against a same-plan
baseline, plus 7 channels × {`pom_set`, `coherent_t2`} with `--hadronize`
hashing every particle (3278 of them, verified reproducible across processes).

### D6.2 The mechanism

`Pipeline::knob_provenance(KnobRunContext)` returns **one row per
user-settable knob** — 66 of them on the shipped default run — carrying
`{name, flag, value, status ∈ {read, not-read, refused}, reason, label,
at_default}`. Three consumers read that one table and nothing else:

* the npz/HFS `meta` — a new `knob_provenance` block with the whole table, and
  `pol_sf`, `rc_scope`, `coherent_t2`, `pom_set`, `pom_rescale` now written
  **through** it, so a not-read knob's key carries the label instead of a
  value (the `rank2_none_label` rule, generalised);
* the CLI run banner — one `KNOB PROVENANCE` block
  (`cli.knob_provenance_lines`), replacing the hand-written per-knob reach
  sentences, of which the `--unpol-sf` / `--pol-sf` one was **false on the
  default plan** and the PomSet one was printed on the coherent channel only;
* `python/tests/test_knob_provenance.py`, which rebuilds the matrix and
  asserts the table against the **output hash**.

`KnobRunContext` carries the two things the core cannot see: the T2 bridge
(the optional PYTHIA tier, recovered from `PythiaHadronizerHook` by the
metadata writer itself) and what the caller asked the plan factory for (a
`RunPlan` records its moments, not which flags produced them).

**The criterion, stated once** (`KnobProvenance`, pipeline.hpp). A knob this
run does not read is **REFUSED** when the value it would record names a
*variation of a piece that did not run* — a scale, a band edge or a shape on a
term the run computes as identically 1 — because such a value claims a
systematic was PRICED, and no label makes a priced systematic un-priced. It is
**LABELLED** when it names a *backend, an axis or a member of a family that a
channel-, plan- or set-scan sets uniformly across runs*, because refusing one
cell of such a scan costs more than it buys. Either way it is **written**.

### D6.3 The reach table, printed by the code

`R` = read, `-` = not read (labelled), `X` = refused. Plans: TT
tensor-thirds, HF helicity-flip, XT transverse-tensor. All at `--pz 0.7
--pzz 0.6 --pe 0.7` unless the column says otherwise, `--rc off`, no
`--hadronize`. Recipe:

```
# validation-style, run 2026-09-05:
#   cfg  = make_config(isotope=I, channel=C, events=50, seed=7)
#   plan = make_plan(P, j=J, pz=0.7, pzz=0.6, pe=PE)
#   ctx  = KnobRunContext(); ctx.plan_name = P; ctx.pz/pzz/pe = ...
#   Pipeline(cfg, plan).knob_provenance(ctx)
```

| knob | incl-6Li/TT | incl-6Li/HF | incl-7Li/HF | incl-d/TT | coh-6Li/TT | coh-6Li/HF | tag-6Li-a/TT | tag-6Li-a/HF | tag-7Li-a/HF | tag-d-p/TT | tag-d-p/HF | incl-6Li/XT | incl-6Li/HF pe0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `isotope` | R | R | R | R | R | R | R | R | R | R | R | R | R |
| `beam_config` | R | R | R | R | R | R | R | R | R | R | R | R | R |
| `channel` | R | R | R | R | R | R | R | R | R | R | R | R | R |
| `optics` | - | - | - | - | R | R | R | R | R | R | R | - | - |
| `n_sigma` | - | - | - | - | R | R | R | R | R | R | R | - | - |
| `pot_config` | - | - | - | - | R | R | R | R | R | R | R | - | - |
| `seed` | R | R | R | R | R | R | R | R | R | R | R | R | R |
| `run` | R | R | R | R | R | R | R | R | R | R | R | R | R |
| `events` | R | R | R | R | R | R | R | R | R | R | R | R | R |
| `lumi_pb` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `poisson` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `apply_optics_lumi_fraction` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `with_virtual_photon` | R | R | R | R | R | R | R | R | R | R | R | R | R |
| `pz` | R | R | R | R | R | R | R | R | R | R | R | - | R |
| `pzz` | R | - | - | R | R | - | R | - | - | R | - | R | - |
| `rel_lumi_offset` | R | R | R | R | R | R | R | R | R | R | R | - | R |
| `pe` | - | R | R | - | - | R | - | R | R | - | R | - | R |
| `x_max` | R | R | R | R | - | - | R | R | R | R | R | R | R |
| `kernel` | R | R | R | R | R | R | X | X | X | X | X | R | R |
| `unpol_sf` | R | R | R | R | R | R | R | R | R | R | R | R | R |
| `pol_sf` | - | R | R | - | - | - | - | R | R | - | R | - | - |
| `b1_model` | R | R | X | R | X | X | X | X | X | X | X | R | R |
| `b1_band_scale` | X | X | X | X | X | X | X | X | X | X | X | X | X |
| `b1_alpha_d_dwave_weight` | X | X | X | X | X | X | X | X | X | X | X | X | X |
| `b1_unpol` | X | X | X | X | X | X | X | X | X | X | X | X | X |
| `cluster_wave` | X | X | X | X | X | X | R | R | R | R | R | X | X |
| `cluster_vmc_mc_sigma` | X | X | X | X | X | X | X | X | X | X | X | X | X |
| `cluster_beta` | - | - | - | - | - | - | R | R | R | R | R | - | - |
| `p_d` | - | - | - | - | - | - | R | R | - | - | - | - | - |
| `triton_sf` | - | - | - | - | - | - | - | - | R | - | - | - | - |
| `tier` | - | - | - | - | - | - | R | R | R | R | R | - | - |
| `inclusive_b1` | - | - | - | - | - | - | R | R | - | - | - | - | - |
| `fsi` | X | X | X | X | X | X | R | R | R | R | R | X | X |
| `fsi_sigma_mb` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc` | R | R | R | R | R | R | R | R | R | R | R | R | R |
| `rc_delta_low_x` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_delta_high_x` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_x_low` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_x_high` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_a_transfer_frac` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_band_tau_max` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_scope` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_fq_scale` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_tail_tensor_scale` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_c0_shape` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_qe_suppression` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_qe_tensor_scale` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_qe_kf_gev` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_tail_model` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_with_qe_tail` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_with_tail` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_n_eta` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_tail_max` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `rc_m_lepton` | X | X | X | X | X | X | X | X | X | X | X | X | X |
| `coherent_t_max` | - | - | - | - | R | R | - | - | - | - | - | - | - |
| `coherent_f0` | - | - | - | - | R | R | - | - | - | - | - | - | - |
| `coherent_slope_b` | - | - | - | - | R | R | - | - | - | - | - | - | - |
| `coherent_amp` | - | - | - | - | R | R | - | - | - | - | - | - | - |
| `coherent_eps_b0` | - | - | - | - | R | R | - | - | - | - | - | - | - |
| `coherent_m_x_min` | - | - | - | - | R | R | - | - | - | - | - | - | - |
| `coherent_x_pom_max` | - | - | - | - | R | R | - | - | - | - | - | - | - |
| `coherent_weighted_azimuth` | - | - | - | - | R | R | - | - | - | - | - | - | - |
| `hadronize` | R | R | R | R | R | R | R | R | R | R | R | R | R |
| `coherent_t2` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `pom_set` | - | - | - | - | - | - | - | - | - | - | - | - | - |
| `pom_rescale` | - | - | - | - | - | - | - | - | - | - | - | - | - |

The rows that are the five rounds, in one line each:

* **`pol_sf`** — `-` under every TT and XT column and under `HF pe0`, `R`
  under HF at `--pe 0.7` off coherent, `-` on coherent under both. Two axes,
  one predicate (`pol_sf_is_read(config, plan)`).
* **`pe`** — `-` under every tensor plan, because
  `tensor_thirds_plan` / `transverse_tensor_plan` / `tensor_flip_plan` build
  every category at `lam_e = 0` (src/core/bookkeeping.cpp) and P_e multiplies
  nothing. `R` under `HF pe0`: the plan **did** consult `--pe`, and g₁ still
  did not run — two different questions, two rules.
* **`cluster_wave` / `triton_sf` / `inclusive_b1` / `cluster_beta` / `p_d` /
  `tier`** — every one of them now has a row, on every channel. Measured
  2026-09-05, 300 events seed 7: `--cluster-beta 0.33` moves
  tagged-7Li-alpha **even under `--cluster-wave vmc`** (the alpha-t relative
  motion is the ANL table, but the T1 triton breakup is still analytic at
  that β) and does **not** at `Tier::T0`; `--triton-sf ciofi-simula` moves
  tagged-7Li-alpha at T1 and nothing at T0; `--p-d` moves the 6Li alpha tag
  under Hulthén only; `--inclusive-b1` moves tagged-6Li-alpha only, because
  the other two tagged channels' struck clusters are spin ½ and
  `InclusiveKernel::tables` opens no rank-2 sector there.
* **`x_max`** — `-` on coherent, and the reason is **this run's own measured
  edge**, not a quoted ≈ 0.05: every cell carrying coherent rate lies below
  **x = 0.1**, the largest UPPER cell edge (`exp(logx_hi)`) with a non-zero
  `cell_rate_weights_pb` — the largest such cell's CENTRE is 0.0954992586 —
  so `--x-max 0.95` clips nothing there. A value below that edge would move
  it, and the label says so, with the number in it.
* **`optics` / `n_sigma` / `pot_config`** — `-` on the inclusive channel, and
  the claim is absolute rather than statistical: `route_of` returns
  `Route::Lost` **before** it looks at an envelope when the event carries
  neither a tagged spectator nor an intact recoil. Measured, 2000 events seed
  11: the whole `route` column is 0 on `--channel inclusive` (6Li and d), and
  `--optics yr-high-divergence`, `n_sigma` 3 and 1, and `--pot-config 18x275`
  each leave it bit-identical. On the read side the caveat is stated in the
  reason itself: all four move the route column on tagged-6Li-alpha and
  **none** of them does on tagged-d-p at config 1, whose proton spectator is
  outside the beam band by rigidity.

### D6.4 The rc families, under `--rc tensor-band`

Same symbols, at `--rc tensor-band` (the rc rows above are all `-` because
that table's runs are `--rc off`, where no `RcModel` is built at all).

| knob | incl-6Li/TT | incl-6Li/XT | coh-6Li/TT | tag-6Li-a/TT | tag-7Li-a/HF pe0 | tag-d-p/TT |
|---|---|---|---|---|---|---|
| `rc` | R | R | R | R | R | R |
| `rc_delta_low_x` | R | R | X | R | R | R |
| `rc_delta_high_x` | R | R | X | R | R | R |
| `rc_x_low` | R | R | X | R | R | R |
| `rc_x_high` | R | R | X | R | R | R |
| `rc_a_transfer_frac` | R | R | X | R | R | R |
| `rc_band_tau_max` | R | R | X | R | R | R |
| `rc_scope` | - | R | X | - | - | - |
| `rc_fq_scale` | R | R | X | X | X | X |
| `rc_tail_tensor_scale` | R | R | X | X | X | X |
| `rc_c0_shape` | R | R | X | X | X | X |
| `rc_qe_suppression` | R | R | X | X | X | X |
| `rc_qe_tensor_scale` | R | R | X | X | X | X |
| `rc_qe_kf_gev` | R | R | X | X | X | X |
| `rc_tail_model` | R | R | X | X | X | X |
| `rc_with_qe_tail` | R | R | X | X | X | X |
| `rc_with_tail` | R | R | - | - | - | - |
| `rc_n_eta` | R | R | X | X | X | X |
| `rc_tail_max` | R | R | X | X | X | X |
| `rc_m_lepton` | X | X | X | X | X | X |

Three things in that table are new, and two of them are **behaviour changes**:

1. **The tail sub-knobs are now REFUSED on the tagged channels** and every rc
   sub-knob on coherent. `validate()` gained the class rule, of which the
   pre-existing `qe_tensor_scale`-without-`with_qe_tail` clause is the first
   member; `rc_band_applies` / `rc_tail_applies` (rc.hpp) are `RcModel`'s own
   two predicates, read there rather than re-derived. `--rc tensor-band`
   itself stays accepted on every channel, and every sub-knob stays available
   wherever its piece runs — what is refused is the combination.
2. **`rc_scope` is not read under a θ_S = 0 fill**, and it is LABELLED rather
   than refused because that condition is a property of the RUN PLAN and
   `validate()` has no plan. `TensorAll` differs from `TensorRate` only in the
   cos 2φ amplitude, which `tensor_amplitudes` builds with a sin²(θ_S) factor.
   Measured 2026-09-05, 400 events seed 7, inclusive 6Li `--rc tensor-band`:
   bit-identical under `tensor-thirds`, and it **does** move under
   `transverse-tensor` and `tensor-flip`, whose categories sit at θ_S = π/2.
   `meta["rc_scope"]` therefore goes through the table.
   `python/tests/test_rc.py::test_meta_distinguishes_a_tensor_all_no_qe_run`
   used to assert "it really is a different file, not only a different label"
   while checking only the label; it now measures both sides.
3. **`rc_with_tail` is not a tail row**, and could not be: it is the knob that
   DECIDES whether the tail runs, so refusing it "because the tail did not
   run" would be circular. It is read wherever the tail would apply if it were
   true, and labelled on the channels that already force `rc_tail == 1`.

### D6.5 The T2 tier

| knob | coherent `--hadronize` | coherent, no `--hadronize` | inclusive `--hadronize` | inclusive, no `--hadronize` | tagged-6Li-alpha `--hadronize` | tagged-6Li-alpha, no `--hadronize` |
|---|---|---|---|---|---|---|
| `hadronize` | R | R | R | R | R | R |
| `coherent_t2` | R | - | - | - | - | - |
| `pom_set` | R | - | - | - | - | - |
| `pom_rescale` | R | - | - | - | - | - |

`meta["pom_set"]` stays the **int** it always was wherever it ran, and carries
the label where it did not — the `b1_model` / `pol_sf` rule. The context's T2
half is filled by the metadata writer from the bridge **actually bound** to
the config (`PythiaHadronizerHook`), never from what a caller says.

### D6.6 What was added

1. **`KnobStatus` / `KnobProvenance` / `KnobRunContext` and
   `Pipeline::knob_provenance`** (`include/lipolgen/pipeline.hpp`,
   `src/core/pipeline.cpp`), with the criterion stated once above
   `KnobStatus`. It is a method on the RUN and not on the config because three
   rules are properties of the run: the rc reaches are `RcModel::applies()` /
   `tail_applies()` on the model this run built, the tagged breakup rows key
   on `Pipeline::tier()` (which is T0 on channels with no struck cluster
   whatever the config says), and `x_max`'s coherent reach is decided against
   this run's own upper accepted-x edge.
2. **`pol_sf_is_read` / `pol_sf_reach_report` / `pol_sf_unread_label` take the
   RunPlan**, which is round 2's fix at its own site; `plan_has_beam_helicity`
   is the new one-line predicate they and the `pe` row share.
   `pol_sf_unread_label` is now SHORT on both axes, because it is the string a
   `meta` key carries and `KnobProvenance::label` derives the same string from
   the report's opening clause — two spellings of one label is how a file and
   a table start disagreeing.
3. **`cluster_wave_name` and `triton_sf_name`** (cluster.hpp, pipeline.hpp) —
   the two selectors that had no name function while every other one did.
4. **`rc_band_applies` / `rc_tail_applies`** (rc.hpp) — `RcModel`'s per-channel
   rule, extracted so `validate()` can read it.
5. **The `validate()` class rule** for rc sub-knobs whose piece did not run.
6. **`meta["knob_provenance"]`**, and `pol_sf` / `rc_scope` / `coherent_t2` /
   `pom_set` / `pom_rescale` routed through the table. Every other existing key
   keeps its name, its type and its meaning.
7. **The CLI `KNOB PROVENANCE` block**, and the deletion of the hand-written
   reach sentences it replaces: `--unpol-sf reaches EVERY kernel this run
   builds …; --pol-sf reaches the inclusive and the tagged kernels`, the
   `--pol-sf DID NOT RUN HERE` block, and the coherent-only PomSet reach line
   (its measured BAND stays — that is a physics number, not a reach claim).
8. **`python/tests/test_knob_provenance.py`** — the matrix, 409 cases in 58 s.

### D6.7 What was NOT done, and is therefore not claimed

* **`validate()` was not rewritten to read the table.** That would make one
  encoding instead of two, but it would also mean re-deriving every existing
  refusal message through a new path, and the brief's own rule is that no
  existing refusal may be weakened. The table's `refused` rows CITE
  `validate()`, and `test_knob_provenance.py` asserts the two agree in both
  directions: a `refused` row whose non-default value is accepted fails, and a
  `not-read` row whose non-default value is refused fails.
* **The HepMC3 writer gained no provenance block.** It carries per-event
  attributes and a `GenRunInfo` of weight names only; it has never carried the
  configuration, and adding one would move every HepMC file. The npz and HFS
  `meta` are the provenance surface.
* **`with_virtual_photon`, `lumi_pb`, `poisson` and
  `apply_optics_lumi_fraction` are in the table but not in the hash matrix**,
  each with the reason in `EXCUSED`: the first changes the event RECORD and no
  generated column, and the other three are reachable only by leaving
  fixed-count mode, which moves `events` at the same time — two knobs, not one.
* **No generated array moved.** `git status validation/reference/` is empty,
  the C++ reference gates pass at their own rtol, and the shipped defaults are
  unchanged: 6Li inclusive 591846.168878 pb, coherent 12476.025113 pb, 7Li
  inclusive 590952.426415, d inclusive 591260.012713, tagged-6Li-alpha
  591846.168878, tagged-7Li-alpha 589760.769797, tagged-d-p 585044.228498
  (2000 events, seed 20260713).

---

## D7 — the claim residues §D6 left behind, and the one §D6's own fix created

A fourth pass, 2026-09-05. §D6 closed the *reach* defect with a mechanism; the
review that followed it found **seven surviving claim residues** of the SAME
classes the earlier passes had been closing — one of them created by §D2.8's
own correction, and one of them (F) spanning fourteen sites on its own. All
seven are fixed in place, with the measurement beside each. Nothing here
moved a generated array: `build/lipolgen_tests` still 400 cases / 17 240 262
assertions / 1 skipped / 0 failed, `pytest` 848 passed / 94 skipped, and the
reference gates pass unchanged.

**A. A window slip in a σ column, in the paragraph that had just fixed one.**
§D2.1's correction note asserted "The header and the σ column say the shipped
window" while the σ column read 591846.16 — the `x_max` = 0.95 value.
Re-measured 2026-09-05: `sigma_pb()` at config 1 is 591846.168878 pb at the
shipped `x_max` = 1.0 against 591846.161405 at 0.95, so the column now reads
**591846.17** on the `toy` and `--pol-sf nnpdfpol` rows of both tables. Those
are the only two rows either window can be told apart in: `ct18nlo` is
472571.917912 / 472571.917091 and `mstw` 469556.451348 / 469556.450465, both
rounding to the printed 472571.92 / 469556.45 either way. The
`--inclusive-b1 --x-max 0.95` sub-table keeps 591846.16 (measured there,
591846.162356). The A_zz column had already been moved into the shipped window
by §D2.8 and the A_∥ column was in it all along; neither moved again here — the
σ column was the one clause of that correction note that was still false when
it was written. See §D2.1.

**B. A citation the gate could not see.** `docs/PHYSICS_CHANNELS.md`'s
coherent-scenario bullet cited `python/lipolgen/cli.py:421` for
`--coherent-t-max`; that line is inside the `--x-max` help string. The flag is
declared at cli.py:445 (444 before §D7's own `--b1-model` edit), and the
citation now carries the flag NAME, which moves it from the gate's shape B
(unnamed: only "a line somebody could have meant") to shape A, where a drift
is caught.

**C. Set 4's surviving 29.32 % at two more sites.** §D4.6's "four sites
summarised the table without set 4" was itself two short: the
`n_pom_flavour_fallback` binding docstring in `python/bindings.cpp` said
"removes the fallback and leaves the final state bit-identical", and the block
comment above the light-only branch in `src/pythia/pythia_bridge.cpp` said
"removes the fallback (or nearly)" — that parenthesis was the whole of set 4.
Both now carry 35.52 % → 29.32 % and the "set 4 INCLUDED" bit-identity, and
§D4.6 says **six sites**.

**D. Two bare claims on one line of `STATUS.md` row D.** D3's "they differ on
99.50 % of events" now carries its window (one event stream, `--channel
tagged-6Li-alpha --events 20000 --seed 1234`, `tensor-thirds` at
P_z = 0.7 / P_zz = 0.6 / P_e = 0.7, at the default σ_XN = 40 mb, one end of the
mandatory 20–40 mb band) and the "by more than 1 %" it was always about; D4's
−2.6 %/+9.9 % is attached to the **twelve DPDF fits**, not to "all 15 sets"
(over all 15 the low edge is −4.7 %). `OPEN_ITEMS_SOLUTIONS.md` §5.1's count of
sites that had made the second mistake goes from three to **four**.

**E. Four figures on a two-figure input, at four sites.** §D5.1's table cell and
its prose, and `docs/USAGE.md`'s code block and quadrupole-band table, printed
**2.800** — neither of the two numbers there are. Measured 2026-09-05:
`t_positivity_edge(-2.0)` is 2.8000000000 on the ROUNDED `eps_b0` = −0.0070 and
**2.799047** on the derived −0.0070024 (`quadrupole_band_fm2()` through
`a2_from_quadrupole` at B = 50, amp = 0.01; the other two band members are
−0.0171207 → 1.144810 and −0.0526846 → 0.372025). All four sites now say
**2.80**, with the derived value beside them and "never 2.8000" — the wording
every other site in the tree already carried.

**F. "Bit-identical across all 15 sets" reproduces over 14.** Re-measured
2026-09-05 on §D4.1's own recipe: the md5 of `t ‖ x_pom ‖ q2 ‖ x ‖ weight` over
20 000 events is `ffd35a3a62b591c547e9ca2ac4301b5d` on **every one of sets
1–10 and 12–15** — one md5, 14.3 s wall — and `pom_set = 11` throws in
`PythiaBridge`'s constructor before the first event, so its row cannot be
re-run at all. Set 11 gave the same md5 in the 2026-09-04 scan; the scan as run
was over fifteen and the claim that REPRODUCES is over fourteen. **Fourteen
sites outside this file now say fourteen**: `python/lipolgen/cli.py` help and
banner (2), `docs/USAGE.md` §4, `docs/OPEN_ITEMS_SOLUTIONS.md` §5.1 and its
item-5 summary row (2), `docs/DEVELOPMENT_PLAN.md`, `docs/PYTHIA_BRIDGE.md`
§11, `docs/T2_CHAIN.md`, `docs/PHYSICS_CHANNELS.md` row 459,
`include/lipolgen/pythia_bridge.hpp`, `python/bindings.cpp`, `PLAN.md` D4,
`STATUS.md` row D and `phase_D_small_items.md` §D4.3 — plus §D4.2 here.

**Nit.** `PHYSICS_CHANNELS.md`'s `LhapdfSF` row said that handing an object
straight to a kernel's `Options::f2_source` means "nothing records it".
Measured 2026-09-05 on a config carrying a `default_inclusive_kernel`:
`meta["unpol_sf"]`, `["pol_sf"]`, `["b1_model"]` and `["b1_unpol"]` all read
`caller-supplied kernel`. What is missing is the BACKEND, not the record, and
the row says that now — with `python/bindings.cpp:608` `unpol_sf` cited, which
is §D7's one new reference.

### D7.1 The reach sentences §D6 made obsolete

§D6 replaced the hand-written reach claims in the BANNER; the same claims
survived in prose. Re-checked against `Pipeline::knob_provenance` on four
channels × two plans (2026-09-05):

* `unpol_sf` is `read` on **every** channel and every plan — so `docs/USAGE.md`
  §2b's and `PHYSICS_CHANNELS.md` row 163's "every kernel the pipeline builds"
  survive as written.
* `pol_sf` is `read` only on `inclusive` and the tagged channels **and** only
  under `helicity-flip` at `--pe ≠ 0`; it is `not-read` on `coherent` under
  either plan and `not-read` on every channel under `tensor-thirds`, the CLI's
  own default. `docs/USAGE.md` §2b said only the channel half; it now says both,
  as does the "what is not refused" list further down.
* `--b1-model cdks|li6-convolution` is **refused** off `inclusive` and off ⁶Li,
  not merely unread: `python/lipolgen/cli.py`'s help said it "would reach the
  metadata but not the rate", which had stopped being true.
* `OPEN_ITEMS_SOLUTIONS.md` §15.3's "`--unpol-sf` / `--pol-sf` already reach
  every ⁷Li kernel" is true of the unpolarised half only; a ⁷Li b₁ is a TENSOR
  observable and is measured on exactly the plans where `pol_sf` is not read,
  so that paragraph now says where the polarised input has to come from.

## D8 — the fourth pass (2026-09-05): six cells the mechanism still slipped, and eleven text residues

Every number in this section was measured on **2026-09-05** on this tree, after
§D7, by the recipe printed beside it. §D6 built the mechanism, §D7 reconciled
the prose §D6 made obsolete; an independent re-run of the 336-cell matrix found
that **the mechanism itself still had six edge cells**, and a second reviewer
found **eleven surviving hand-written `--pol-sf` reach claims** plus two number
cells and one citation. All of them are closed here.

Suites after this section, against the "after §D7" column at the top of this
file:

| gate | after §D7 | after §D8 (this tree) |
|---|---|---|
| `build/lipolgen_tests` | 400 cases / 17 240 262 assertions / 1 skipped / 0 failed | **400 cases / 17 240 262 assertions / 1 skipped / 0 failed** |
| `python -m pytest python/tests -q` | 848 passed, 94 skipped | **888 passed, 112 skipped** (measured before the §D7–F pass and §D9 added their tests; `de1a040` as committed collects 1021 and runs **909 passed, 112 skipped**, re-measured 2026-09-05 — the figure `STATUS.md` row D and `phase_E_numbers.md` §E0 carry) |
| `python3 validation/check_physics_channels_links.py` | 1145 references checked (strict), 97 ranges, 7 external, 0 broken, 6 allow-listed | **1145 references checked (strict), 97 ranges, 7 external, 0 broken, 6 allow-listed** |

The C++ column does not move: every one of the six defects is on the run
SURFACE. The pytest column moves by **+58 tests** — 56 new matrix cells
(3 route knobs + 2 coherent `x_max` edges + 3 quasi-elastic sub-knobs, each on
the 7 full specs), one new test for the coherent `pot_config` cell and one more
channel on the route test — of which **+18 land as skips, every one of them a
cell whose BASE is refused** (`rc_options.with_qe_tail = false` is itself a
tail knob, so on the coherent and the three tagged channels the base of a
quasi-elastic cell does not run and the cell says so). `python/tests/test_knob_
provenance.py` alone is **449 passed, 112 skipped** (526 matrix cells over 12
specs and 68 variants, plus 35 other tests).

**Bit for bit at every default.** σ and the sha256 over every generated array
column, 2000 events at seed 20260713 on all seven (isotope, channel)
combinations, are unchanged: 591846.1688784091 (inclusive ⁶Li),
590952.426415097 (⁷Li), 591260.0127125566 (d), 12476.025113185518 (coherent),
591846.1688784088 (tagged-⁶Li-α), 589760.7697973427 (tagged-⁷Li-α),
585044.2284977469 (tagged-d-p). Nothing under `validation/reference/` moved.

### D8.1 `meta["optics"]` and `meta["pot_config"]` bypassed the table

`python/bindings.cpp` wrote both keys unconditionally —
`meta["optics"] = p.optics().name; meta["pot_config"] = p.pot_config();` —
so on the INCLUSIVE channel, where the table says not-read and the banner's
provenance block prints the label, `--optics yr-high-divergence` recorded
`meta["optics"] = "10x100 high-divergence"` and `meta["pot_config"] = "10x100"`
as bare values while all 47 generated columns were bit-identical (16 cells over
the five inclusive specs of the §D6 matrix). That is the `pol_sf` defect on two
top-level keys, and it survived the fix that was meant to close the class
because these two never went through the mechanism.

Both keys are `knob("optics").meta_value()` / `knob("pot_config").meta_value()`
now, and so is the banner's **header line**, which printed `p.optics.name` on
every channel — announcing the envelope as a run property three lines above a
block saying it is not read. Measured after the fix:

```text
$ python -m lipolgen.cli --events 200 --channel inclusive --optics yr-high-divergence
inclusive  e(10) x 6Li(99.5/u)  not read on channel inclusive
$ python -m lipolgen.cli --events 200 --channel coherent
coherent-6Li  e(10) x 6Li(99.5/u)  10x100 high-acceptance
```

Nothing is lost: the bare name stays in
`meta["knob_provenance"]["optics"]["value"]`.

> **AMENDED, phase F, 2026-09-05.** "Nothing is lost" was true of the inclusive
> channel, where the classifier is never reached, and FALSE of every channel
> that HAS a route. There the label was the bare scope clause "not read on this
> run", so `meta["optics"]` on a `tagged-d-p` run — 95.75 % of whose events
> carry an optics-classified route — went from `10x100 high-acceptance` at
> `a94fd6e` to `not read on this run`, and the banner printed that three lines
> above a tag fraction quoted AT that envelope. The status is right (the labels
> are insensitive to the one alternative envelope tabulated for d at 10x100);
> the scalar said more than the row it came from. Fixed in the ONE place the
> label is built: on a channel with a route the reason now opens
> `not read on this run at 10x100 high-acceptance (the classifier IS consulted
> here; this run's route labels are insensitive to yr-high-divergence)`, and
> `n_sigma` and `pot_config` carry their own value the same way, so the name is
> in the scalar, in the banner header and in the row — with no per-key branch
> in `bindings.cpp`. Measured after that fix:
>
> ```text
> $ python -m lipolgen --isotope d --channel tagged-d-p --plan tensor-thirds --events 400 --seed 4242
> tagged-d-p  e(10) x d(100/u)  not read on this run at 10x100 high-acceptance (the classifier IS
>     consulted here; this run's route labels are insensitive to yr-high-divergence)
> $ python -m lipolgen --events 200 --channel inclusive --optics yr-high-divergence
> inclusive  e(10) x 6Li(99.5/u)  not read on channel inclusive
> ```

### D8.2 The coherent `x_max` rule was self-referential, and wrong at both edges

The rule read this run's own upper rate-carrying x edge off
`cell_rate_weights_pb()` — i.e. off the ALREADY-CLIPPED cell set — and compared
it exactly. Both halves were wrong, in opposite directions. Measured (⁶Li
coherent, config 1, `tensor-thirds` at P_z = 0.7 / P_zz = 0.6):

| `--x-max` | cells | carrying rate | σ [pb] | row, before | row, now |
|---|---|---|---|---|---|
| 1.0 (default) | 3051 | 1726 | 12476.025113185518 | not-read | **not-read** |
| 0.95 | 3027 | 1726 | 12476.025113185518 | not-read | **not-read** |
| 0.10 | 1822 | 1726 | 12476.025113185518 | **read** | **not-read** |
| 0.09 | 1770 | 1702 | 12475.91392063426 | read | **read** |
| 0.05 | 1407 | 1406 | 12446.10902999575 | **not-read** | **read** |

At 0.05 the clipped set lies below 0.05 by construction, so the edge came back
**0.047863009232263845**, `x_max >= edge` held, and the row said NOT READ of a
value that had just removed **320 of the 1726** rate-carrying cells and moved σ
by 29.916 pb. At 0.10 the edge is `exp(logx_hi)` of a grid boundary and comes
back **0.10000000000000002**, one ulp above the value that produced it, so a
run bit-identical to the default said READ.

`Pipeline::coherent_rate_x_edge()` now measures the edge on the UNCLIPPED cell
set — a scratch `InclusiveSampler` at the shipped `Scenario::x_max`, with this
run's own kernel, beams and grid, built only when `x_max` is away from its
default (at the default the run's own cells ARE the unclipped set) — and the
comparison carries a relative tolerance of 1e-12. The matrix gained `x_max`
cells at **0.05** and **0.10**; its only one was 0.95, which is on the same
side of the edge as the default and could see neither error.

The same cell exposed one thing about the matrix test itself: it judged the
BASELINE's row against a measurement of the VARIANT. At `x_max = 1.0` the row
honestly says not-read ("a value BELOW the edge would move this channel") while
0.05 moves the file, so the two are both true and only the variant's row is a
statement about the value that moved. The cell reads the variant's row now.

### D8.3 The three route knobs are measured on the run, not tabulated by channel

`--optics`, `n_sigma` and `pot_config` reach exactly one quantity, the per-event
`route` column, and through exactly two branches of `route_charged`:
`Optics::clears` only for a NEAR-BEAM fragment (|R − 1| < 0.05, θ < 5 mrad),
`over_rigid_route` only for an OVER-RIGID one (R > 1.05). The table said `read`
on all three wherever a route exists and hung a hand-written caveat about
tagged-d-p on the read side; the matrix flagged **22 "did not move but read"**
cells, and `python/tests/test_knob_provenance.py` excused all three knobs from
the hash test as "statistics-dependent at 60 events".

Measured (config 1; `route` column compared cell by cell against the same run's
baseline; the movers are named with the number of labels they move):

| channel | 60 events, seed 7 | 2000 events, seed 11 |
|---|---|---|
| inclusive ⁶Li / ⁷Li / d | none of the eight alternatives moves anything | none |
| coherent-⁶Li | optics `tagging` 13, `tagging-legacy` 17; n_sigma 1 → 31; **pot_config nothing** | optics 529 / 681; n_sigma 1 → 1007, 3 → 12; **pot_config nothing** |
| tagged-⁶Li-α | optics 18 / 20; n_sigma 29 / 3; pot_config 18x275 → 2 | optics 1 / 462 / 560; n_sigma 790 / 88 / 1; pot_config 7 / 25 |
| tagged-⁷Li-α | **optics nothing**; n_sigma 30 → 1; **pot_config nothing** | optics 3 / 32 / 32; n_sigma 37 / 26 / 6; pot_config 5x41 → 2 |
| tagged-d-p | optics `yr-high-divergence` → 1; n_sigma 30 → 1; **pot_config nothing** | **none of the eight moves anything** |

So the answer is not a property of the channel, and it is not a property of the
statistics either: it flips BOTH ways between the two columns (tagged-⁷Li-α's
`--optics` moves nothing at 60 events and does at 2000; tagged-d-p's moves at 60
and nothing at 2000). A tabulated rule would be wrong in one column or the
other whichever way it was written.

`Pipeline::RouteReach` measures it instead, by the idiom the header above
`route_of` already states: **the sample is drawn ONCE and re-routed per
optics**. The probe walks this run's own events (without the rc weights and
without the T2 hadronizer — neither moves a four-vector, and PYTHIA is not
re-entrant and its stats are in `meta`), re-routes each at the other three
tabulated envelopes, at n_sigma ∈ {1, 3, 30} and at the other machine
configurations, and stops as soon as all three knobs have been seen to move.
It costs one extra generation pass on the channels that write a far-forward
fragment (0.20–0.26 s per 100 000 events, memoized per `Pipeline`, and on
tagged-⁶Li-α at 60 events it exits after **5**), and nothing at all on the
inclusive channel. Each row states what was measured and against what.

`--optics` has one non-route reader and it is now tested first: in LUMINOSITY
mode `apply_optics_lumi_fraction` multiplies every category's COUNT by
`Optics::lumi_fraction`, which reaches every column on every channel.

The three `EXCUSED` entries are gone. The matrix cells for these knobs are
ANY-OF cells — the row is `read` iff at least one value the run HAS moves the
hash — over exactly the value lists the core re-routes at, so the cell and the
table are the same statement measured twice and independently; a value the run
surface refuses (the tagging working point is tabulated for ⁶Li and ⁷Li only,
so `--optics tagging` throws on a deuteron beam) is skipped, not counted as
"did not move". `test_the_route_knobs_are_read_where_a_route_moves` runs the
same rule at 2000 events on six channels, and
`test_pot_config_is_never_read_on_the_coherent_channel` pins the one cell that
is structural rather than statistical: the coherent fragment is the INTACT beam
nucleus, which can only have lost momentum (measured max R ≤ 1), so
`over_rigid_route` is never reached at any statistics.

### D8.4 The quasi-elastic sub-knobs under `with_qe_tail = false`

`rc_qe_suppression` is a flat multiplier ON POLRAD Eq. (44)'s σ^q_U and
`rc_qe_kf_gev` the Fermi momentum of the Pauli factor S(q) INSIDE it. With
`rc_options.with_qe_tail = false` that term is not computed at all, so neither
is a factor of anything — and both were **accepted, bit-identical, and recorded
`read` / non-default** (inclusive ⁶Li, `--rc tensor-band`, 300 events seed 7:
`meta["rc_qe_suppression"] = 0.5`, `meta["rc_qe_kf_gev"] = 0.25`). The stated
criterion — *a scale on a term the run computes as identically 1 is REFUSED,
not labelled* — had exactly one member there, `qe_tensor_scale`, which has been
refused since the knob existed.

All three are refused now, in `PipelineConfig::validate()` and in `RcModel`'s
constructor (the C++ API path), and the rows are a `qe_row` rather than a
`tail_row`: `Read` only where the radiative tail applies AND the quasi-elastic
piece is on, `Refused` where `rc` is on and either fails, `NotRead` at `rc off`.
What is NOT refused, and the reason is in the code: either knob away from its
default with the tail ON — `qe_suppression = 0` degenerates the tail the same
way, but it is then a numeric edge of a term that RAN. The matrix gained a
`with_qe_tail = false` base for all three.

### D8.5 The `hadronize` row looked at the context, not at the config

The row keyed on `KnobRunContext::t2_bound`, which only a caller who went
through `set_pythia_hadronizer` (or the CLI) sets. A hadronizer bound as a plain
callable — `cfg.hadronizer = f`, the route `python/README.md:246` documents, and
the C++ lambda of `docs/USAGE.md:1741` — was called on **every event** (measured:
20 of 20) while the row said `hadronize = off`, "no T2 tier is bound: every
record stops at T0", and the three Pomeron rows said "not read without
`--hadronize`: no PYTHIA instance of any kind is built".

The core can see `PipelineConfig::hadronizer` is bound and the row looks at it
now. What the context adds is not WHETHER a hook is bound but WHICH bridge it
is — the core does not link the PYTHIA tier and cannot ask a type-erased hook
for `PDF:PomSet` — so the Pomeron rows have three cases instead of two, and the
middle one says what it is: *not read through this run's hadronizer: it is a
callable bound to `PipelineConfig::hadronizer` and not a `PythiaBridge`
attached by `set_pythia_hadronizer`, so the value shown is
`PythiaBridgeOptions`' own default, not something this run configured*.

### D8.6 `lg.run` held the plan name and never passed it

`python/lipolgen/__init__.py`'s `run()` called `p.generate(0, keep_events,
nthreads)` with no `KnobRunContext`, so its `meta["knob_provenance"]` came back
with **63 rows and no `pz` / `pzz` / `rel_lumi_offset` row at all**:
`lg.run(plan="helicity-flip", pzz=0.6)` recorded the typed `pzz` nowhere. That
is the silence class the table exists to close, on the library's own entry
point. It passes the context now (**66 rows**; `pzz` = 0.6 `not-read` under
`helicity-flip`, with the max-entropy-ladder reason). The T2 half is left
alone: the metadata writer overwrites it from the bridge actually bound.

### D8.7 Eight more `--pol-sf` reach claims, and where the correct wording is

`include/lipolgen/pipeline.hpp:590-597` has said since §D6 that *the flag pair
reaching "EVERY kernel" is true of `--unpol-sf` alone*, and the table agrees:
`pol_sf` is `not-read` on the coherent channel under both plans and on EVERY
channel under `tensor-thirds`, the CLI's own default. §D7.1 fixed the claim in
`docs/USAGE.md` §2b; it survived at eight more sites, all corrected here to the
two-axis wording of `USAGE.md:628-638` / `PHYSICS_CHANNELS.md:164` /
`OPEN_ITEMS_SOLUTIONS.md` §15.3 / `cli.py:757-765`:

1. `docs/open_items/run_2026-09-03/STATUS.md` row D — "`--unpol-sf` / `--pol-sf`
   reach every kernel (§14)".
2. `docs/OPEN_ITEMS_SOLUTIONS.md` item-14 row — "reach **every kernel the
   pipeline builds**" said of both flags.
3. `docs/PHYSICS_CHANNELS.md` §3's opening — "the unpolarised and polarised
   backends of **every kernel the pipeline builds** … every channel".
4. `python/bindings.cpp` `PolSfSource` docstring — "to every kernel a Pipeline
   builds", with no scope caveat anywhere in it; its C++ twin carries the SCOPE
   paragraph, and the binding now carries it too.
5. `python/bindings.cpp` `PipelineConfig.pol_sf` — "which POLARISED backend
   every kernel takes its g1 … from".
6. `include/lipolgen/pipeline.hpp` the `pol_sf` FIELD comment — the same bare
   sentence, seven lines below the enum that contradicts it.
7. `python/lipolgen/__init__.py` `POL_SF` ("to every kernel a `Pipeline`
   builds") and `make_config`'s docstring ("the structure-function backends of
   EVERY kernel the run builds — inclusive, coherent and tagged").
8. `docs/open_items/run_2026-09-03/PLAN.md` D2 — the channel half only
   ("labelled, not credited, on the coherent one"), with no plan axis: the
   exact one-axis defect §D7.1 fixed in USAGE §2b.

### D8.8 Two number cells, and one citation that had drifted onto a doc comment

`t_positivity_edge(-2.0)` on the ROUNDED `eps_b0` = −0.0070 is 2.8000000000,
and the derived band gives 2.799047 — so the quotable figure is **2.80**, three
figures, which is all a two-figure input supports, and the "never 2.8000"
caveat sits one paragraph below each table. Three table CELLS still printed
`2.8000`: `phase_D_numbers.md` §D5's, `docs/USAGE.md:1578`'s, and — the same
table, a third copy the review did not name —
`docs/OPEN_ITEMS_SOLUTIONS.md`'s. All three now read 2.80. The prose that
quotes "2.8000" as the figure NOT to publish is left alone, as is the code
comment at `USAGE.md:1509` that shows what the function returns.

`docs/PHYSICS_CHANNELS.md` cited `include/lipolgen/pythia_bridge.hpp:338` for
`n_pom_flavour_fallback`; :338 is the doc comment (`/// Recorded per run since
2026-09-04 as meta["n_pom_flavour_fallback"]`) and :340 the declaration. At
`903fcc9` the citation was on the declaration, so `--fix` had re-anchored it
onto a comment — **and the gate could not see it**, because it agreed with
itself: `cxx_declares` calls `literal_defines` on the RAW line, `mask_cxx`
having erased the string literals it needs, and the raw line carries the
comments too. So any `///` line spelling the name inside a bracketed string
literal counted as a declaration.

`mask_cxx` grew a `keep_strings` mode — comments blanked, string and character
literals kept — and both `cxx_declares` and `Source.in_code` read THAT line.
A doc comment now offers nothing and the declaration is the only target left.
Re-running the gate on this tree found **100 broken** references, of which the
comment-line rule accounts for the ones `--fix` re-anchored; the rest are line
drift from this section's own code. All are re-anchored, 13 use-site
fingerprints re-recorded with **identical line text** (`--record-ranges`: 13
new, 0 changed, 11 dropped) and three allow-list entries moved with their pinned
lines, ending at **1145 references checked (strict), 97 ranges, 7 external,
0 broken, 6 allow-listed** — the §D7 counts exactly.

## D9 — fifth pass, 2026-09-05: the coherent `x_max` edge is a cell CENTRE, and two latent residues recorded

**Fixed.** `Pipeline::coherent_rate_x_edge` returned the top of the last
rate-carrying cell, `exp(logx_hi)` = 0.1 + 1 ulp; the sampler admits a cell on
its CENTRE and the coherent branch never re-tests x per event, so every
`x_max` in (0.0954992586, 0.1] — a 4.5 %-wide window — was bit-identical to
the default (1822 cells, 1726 carrying, σ = 12476.025113185518 pb) while the
row said READ with a reason that named cells it did not remove. The edge is
now the largest rate-carrying centre; 0.097 and 0.0955 are `not-read`, 0.09549
is `read` (1770 / 1702, σ = 12475.913921), and three matrix cells at those
values pin the interval (`test_knob_provenance.py`).

**Recorded, not fixed (latent — zero live instances, no shipped number
affected):**

- (B) In LUMINOSITY mode the `optics` row is asserted `read` through
  `Optics::lumi_fraction` rather than measured against the other envelopes: on
  a deuteron run (`--lumi 1e-3 --optics yr-high-divergence`) `yr-high-acceptance`
  is bit-identical (counts [14, 25, 21] both, lumi factor 1.0, route column all
  Lost) yet the row is `read`. On ⁶Li/⁷Li the tagging working point does move
  the counts, so `read` is right there. The matrix excuses luminosity mode; a
  `--lumi` cell on the deuteron would expose it.
- (C) `check_physics_channels_links.py`'s doc-comment exclusion is C++-only:
  `py_declares` and `Source.in_code` read the RAW line for `.py` files, so a
  `#` comment carrying a flag name is accepted as a declaration and `--fix`
  can re-anchor onto it (constructed in a scratch copy at `cli.py:472`). A scan
  of all named point citations finds no Python citation on a `#` line today.
