<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-27 — phase I5: `PolradFull` becomes the shipped `rc_tail` default (decision D5)

Worker I5, 2026-10-29, on `claude/run-2026-10-26` at `9779edd` (the D4
implementation committed; clean tree apart from the supervisor's `STATUS.md`
edit that set this phase RUNNING). Uncommitted: the supervising session
commits. Decision: `DECISIONS.md` D5. Design: `plan_D5_polradfull_default.md`,
AMENDED by its adversarial check (the plan's last section; session scratch
`checks_1027/D5.json`) as the task's nine amendments name. Every measurement
below was taken on this machine (4 cores, nearly idle, load 0.6–1.5 unless a
line says otherwise) with `source env.sh`; the scripts are in session scratch
`i5/meas/`, the edit scripts (every edit an exact old/new pair, applied by
`i5/edits/apply.py`, which refuses an ambiguous or missing old string and
prints each file's line-count change) in `i5/edits/`.

Conventions: **pre** = `9779edd` before any I5 edit (measured on its build);
**default** = `polrad-full` after this phase; a number is labelled with its
tail model wherever both exist. Line numbers below are the tree as this record
leaves it; they are written as plain text except where a citation is gated.

## 0. What D5 changed, in one paragraph

The default is ONE line, `include/lipolgen/rc.hpp:1503` (`RcTailModel tail_model = RcTailModel::PolradFull;`),
and the CLI's second definition of it is now READ from that line instead of
retyped: `python/lipolgen/cli.py:817` (`rc_tail_model=_l.rc_tail_model_name(_l.RcOptions().tail_model)`).
Two code consequences, both the check's amendment 5, both line-neutral:
(i) the `PolradFull` preconditions (`long double` wider than `double`;
`n_eta >= 64`) now fire only where a tail is BUILT — `src/core/rc.cpp:1391` and
`src/core/pipeline.cpp:616` gate them on `rc_tail_applies(channel, with_tail)`
— so a tagged or coherent `--rc tensor-band` run, which builds no tail, is not
refused by the default on a narrow-`long double` platform (where the refusal's
own advice, `--rc-tail-model t-peak`, is a refused tail knob); (ii)
`rc_options.m_lepton` joins the tail sub-knobs `validate()` refuses where the
tail does not apply (`src/core/pipeline.cpp:745`, the `m_lepton` entry): under the TPeak default the
t-peak `m_lepton` rule refused it there, under PolradFull nothing did, and the
knob row (`refused`) would have contradicted an accepted run — found by reading
`validate()` against the knob rows before the flip, and gated since by the C++
D5 pipeline case and by the knob matrix's `rc_m_lepton … at polrad-full` cells
on the tagged and coherent specs, which run since D5 (§5). A third, of the same kind, came out
of the verification (F5, §12): `validate()`'s `sp_tensor_scale` clause binds only where a tail
is built, so off it the tail sub-knob rule refuses `--rc-sp-tensor-scale` with that cause, and
`RcModel`'s refusal there says "no tail is BUILT" instead of "double-count a term that ran"
(no accept/refuse outcome moved). **Nothing else moves in a number**: `--rc off` (the
generator's default) is bit for bit; every explicit-model run is bit for bit;
no reference JSON carries an rc field; the t4 row's rules name the models.
**What moves** is every `--rc tensor-band` run on the inclusive channel that
does not name a model: it now builds POLRAD Eq. (18)'s exact tail — 12.0 s of
table build instead of 0.16 s (⁶Li config 1, unloaded), mean `rc_tail`
1.020341 → 1.028010 on the 5000-event CLI file, 0.386 % of the table's nodes
clipped at `tail_max` (0 events) — and the label `meta["rc_tail_model"]`, which
the bindings write on every rc-on run, tagged and coherent included.

## 1. The decisive check (design V2): A/B bit identity, and the knob rows

The commands of the design's §2.4 (`python3 -m lipolgen --isotope 6Li
--channel inclusive --events 5000 --seed 20261027 --rc tensor-band
[--rc-tail-model …] --npz …`; and the same with no rc flag), column sha256 over
every array column (the design's `npz_hash.py`, `meta` excluded):

| run | pre (`9779edd`) | after I5 |
|---|---|---|
| `--rc tensor-band --rc-tail-model t-peak` | `2a312c26af9f18d1…` (then the default) | `2a312c26af9f18d1…` (explicit) |
| `--rc tensor-band` default | `2a312c26af9f18d1…` (t-peak) | **`00c012b459e25cfe…`** = pre's explicit `polrad-full` |
| `--rc tensor-band --rc-tail-model polrad-full` | `00c012b459e25cfe…` | — (it is the default run) |
| no rc flag (`--rc off`) | `2c48b4e516fea455…` | `2c48b4e516fea455…` |

All three pre hashes equal the design's `ca9b23a` ones: D4 does not move a
θ_S = 0 run. In `meta` the after/pre pairs differ only in
`knob_provenance["rc_tail_model"]`: `at_default` (true ↔ false) and the row's
reason (it now names the default first), and at `--rc off` its `value`
(`t-peak` → `polrad-full`, status `not-read`). The default run's wall time is
12.72 s against 0.37 s for the explicit t-peak run.

**The knob rows on the three channel classes** (measured on the I5 build,
`lg.make_config(…, rc="tensor-band")`, 200 events): on `tagged-alpha` and
`coherent` `meta["rc_tail_model"]` is now `polrad-full` (it was `t-peak`;
`rc_tail_applies` false), and the rows `rc_tail_model`, `rc_m_lepton`,
`rc_n_eta`, `rc_sp_tensor_scale` are all `refused` at their defaults, with the
"the radiative tail does not apply on this run" reason — `rc_n_eta`'s
description is the tanh-sinh one, `rc_m_lepton` is now refused by the
tail sub-knob rule as well. On `inclusive`, `rc_tail_model` and `rc_m_lepton`
are `read` (the latter was `refused` under the t-peak default), `rc_n_eta`
reads "tanh-sinh nodes PER PANEL", and `rc_sp_tensor_scale` is `refused` for
the OPPOSITE reason (the term ran).

## 2. Every re-taken number, before → after

Tests (each with a dated comment naming D5; t-peak pins name the model):

| file:line | what | before | after |
|---|---|---|---|
| `tests/test_rc.cpp` 3519–3601 (new; F5 appended a SUBCASE, §12) | D5 case: the default is `PolradFull` (`"polrad-full"`, enum values 0/1/2 unchanged), a default-built model equals an explicit one node for node and gives T19(i)'s pinned σ^el_U = 2.2301252856028388e−05 at (0.01, 5); the preconditions bind only where a tail is built (n_eta = 32 refused on the inclusive channel, accepted on tagged-6Li-alpha and coherent, and on a band-only run) | — | new |
| `tests/test_rc_pipeline.cpp` 277–341 (new; F5 appended its `sp_tensor_scale` checks, §12) | D5 case: `validate()` on tagged — a non-default n_eta is refused by the tail sub-knob rule ("RADIATIVE TAIL"), `m_lepton` refused, an explicit t-peak refused; coherent `m_lepton` refused; inclusive `m_lepton` READ, refused at `with_tail = false` | — | new |
| `tests/test_rc_pipeline.cpp` 247 | the reset at the end of the m_lepton SUBCASE | `RcTailModel::TPeak` | `RcOptions().tail_model` |
| `python/tests/test_rc.py` 263 | `meta["rc_tail_model"]` of the default run | `"t-peak"` | `"polrad-full"` |
| `python/tests/test_rc.py` 611–612 | the seven-file chain's default run | `"t-peak"` | `"polrad-full"`; the t+ll comparisons now run against an explicit t-peak base |
| `python/tests/test_rc.py` 650–651 | `cli.DEFAULTS` / `cli.resolve(["--rc", "tensor-band"])` | `"t-peak"` | `== rc_tail_model_name(RcOptions().tail_model) == "polrad-full"` |
| `python/tests/test_rc.py` 481–484 | the s/p refusal (amendment 6) | one clause on the default | clause 1 on t-peak BY NAME (`match="t-peak-only tail does not compute"`), clause 2 the default's double-count branch (`match="double-count"`) |
| `python/tests/test_knob_provenance.py` 395, 400 | the two `rc_tail_model` cells | `t-peak -> t-peak+ll`, `t-peak -> polrad-full` | `polrad-full -> t-peak+ll`, `polrad-full -> t-peak` (as written the second would vary nothing) |

The C++ t-peak pins (§4) and the Python runtime pins (§5) change no number: they
name the model the number was always computed on.

## 3. The published numbers, re-measured on BOTH models (design V9, check amendment 3)

`i5/meas/m_points.py` and `i5/meas/m_extra.py` (⁶Li inclusive config 1, the
CLI's 101 × 77 grid, `tensor_thirds_plan(0.0, 0.6)`, Q² = 5 GeV², every model
named explicitly, run on the pre build; the default rows re-checked on the I5
build by the A/B identity of §1). **Every t-peak number the documents publish
reproduces to every printed digit but two** (corrected by F5, §12): `w_tail − 1`,
whose published values sit 0.01–0.26 % above the table route (§9 item 7), and
the C0 ΔA_zz ratios −8 / +0.6 / +37 %, which reproduce only on the withdrawn
deuteron-b₁ `A_zz` (−8.2 / +0.6 / +36.9 % with ×3.253983, against +18.7 / −4.0 /
+90.7 % on the shipped ⁶Li b₁, the table's last row). The rest is what licenses
printing the default's beside them with the same recipe. ΔA_zz is the recipe of
`run_2026-09-02/phase_C_numbers.md` §8.1c-corr, (A_zz + 2 r_T)/(1 + r_U) − A_zz
with r_U the tail ratio at q_n = 0 and r_T its Q_N coefficient, A_zz(Born) from
the shipped kernel's tables (`test_rc_low_x_anchor.py`'s `AZZ` recipe).

| quantity (x = 0.01 / 0.03 / 0.10 / 0.30 unless stated) | `t-peak` (published, reproduced) | `polrad-full` (the default) |
|---|---|---|
| (1/6)σ^el_T/σ^el_U, `ho` | −5.0943e−04 / −7.3583e−04 / +1.5595e−04 / +1.3287e−02 | −5.1540e−04 / −7.4709e−04 / **+4.4562e−05** / +1.0947e−02 |
| (1/6)σ^el_T/σ^el_U, `vmc-ft` | −5.4705e−04 / −8.0426e−04 / −4.0702e−05 / +1.5556e−02 | −5.5461e−04 / −8.1840e−04 / **−1.6369e−04** / +1.1293e−02 |
| w_tail − 1 (table route) | 5.18222e−04 / 3.54397e−05 / 1.31963e−06 / 6.86176e−08 (the published 5.18398e−04 / 1.32298e−06 / 6.8624e−08 sit +0.034 / +0.254 / +0.009 % above this route: §9 item 7) | 6.62870e−04 / 1.75050e−04 / 1.95963e−04 / 4.88356e−04 |
| σ^q_U/σ^el_U | 0.28237 / 0.54628 / 2.7475 / 1000.2 | 0.63146 / 6.5693 / 533 / 3.6846e+06 |
| quasi-elastic share of the tail | 22.0 / 35.3 / 73.3 / 99.90 % | 38.7 / 86.8 / 99.81 / 100.00 % |
| ΔA_zz of the whole tail | −1.769797e−07 / −6.848432e−09 / +2.127452e−09 / +5.130156e−12 | −1.178978e−07 / +9.821438e−08 / +2.995857e−07 / +2.378230e−08 |
| band lead at x = 0.01, anchors 0.30 / 0.266 / 0.113 | ×767.6 / ×680.6 / ×289.1 (the ×768 / ×681 / ×289) | **×1152.2 / ×1021.7 / ×434.0** |
| band lead at x = 0.01, `a_transfer_frac` 0 / 0.5 / 1 | ×768 / ×858.2 / ×1085.5 | ×1152 / ×1288.3 / ×1629.5 |
| band lead at x = 0.10 | ×4.55e+04 | ×323.1 |
| QRT(S)/QRT(1) at x = 0.01 / 0.10 / 0.30 | 0.4721 / 0.8684 / 1.0000 | 0.6665 / 0.9992 / 1.0000 |
| spread at x = 0.01: `qe_suppression` 0 vs 1 | −2.286503e−07 → −1.769797e−07 (5.167e−08) | −2.343600e−07 → −1.178978e−07 (1.165e−07) |
| … `qe_kf_gev` 0 vs 0.169 | −1.192240e−07 → −1.769797e−07 (5.776e−08) | −5.970809e−08 → −1.178978e−07 (5.819e−08) |
| … `c0_shape` `ho` vs `vmc-ft` | −1.769797e−07 → −2.100985e−07 (3.312e−08) | −1.178978e−07 → −1.525695e−07 (3.467e−08) |
| … `qe_tensor_scale` 0 vs 1 | −1.769797e−07 → −2.931806e−07 (1.162e−07) | −1.178978e−07 → −3.828754e−07 (2.650e−07) |
| … `fq_scale` 0 vs 2 | +2.255019e−07 → −5.785810e−07 (8.041e−07) | +2.908174e−07 → −5.257069e−07 (8.165e−07) |
| … `tail_tensor_scale` 0.5 vs 2 | −1.694451e−07 → −2.055947e−07 (3.615e−08) | −1.102803e−07 → −1.468396e−07 (3.656e−08) |
| `qe_tensor_scale` 1: ΔA_zz moves at x = 0.01 / 0.10 / 0.30 (% of the band) | −1.162e−07 / +3.018e−10 / +1.796e−09 (0.086 / 0.0003 / 0.25 %) | −2.650e−07 / +1.741e−08 / +1.0625e−05 (**0.195 / 0.018 / 1454 %**) |
| … and the tail's tensor part r_T | ×1.282 / ×3.747 / ×1006 | ×1.634 / ×533.8 / ×3.699e+06 |
| σ^el_U(`vmc-ft`)/σ^el_U(`ho`) | ×1.011 / ×1.022 / ×1.104 / ×3.032 | ×1.011 / ×1.023 / ×1.105 / ×2.761 |
| ΔA_zz(`vmc-ft`)/ΔA_zz(`ho`) − 1 at x = 0.01 / 0.10 / 0.30 | +18.7 / −4.0 / +90.7 % (the published −8 / +0.6 / +37 % are the withdrawn deuteron-b₁ ones) | +29.4 / −0.0 / +0.0 % |

Build cost (`Pipeline` construction, `--rc tensor-band` minus `--rc off`,
`i5/meas/m_c0cost.py` and the checker's `build_time.py`): `t-peak` 0.160–0.166 s
(1.67–1.70 s on `vmc-ft`), `polrad-full` 11.89–12.08 s (16.08–16.35 s on
`vmc-ft`). Node clip at `tail_max` = 10 on the default
(`clip_configs.py`, 50 000 events each): 0.3986 / 0.3858 / 0.3600 % of the nodes
on ⁶Li configs 0 / 1 / 2 (by y band 0.471 / 0.124 / 0, 0.456 / 0.124 / 0,
0.440 / 0 / 0 %), **0 events** on each; 0 % on t-peak. Event level, the
5000-event CLI files at P_z 0.7: mean `rc_tail` 1.020341 (t-peak) → 1.028010
(default), max 3.2380 → 4.0053, 0 clipped events.

Two readings change on the default, and both are recorded where the documents
state them, not tuned: **(a) the `qe_tensor_scale` stand-in at x = 0.30** moves
ΔA_zz by 14.5× the band half-width at scale 1 (the quasi-elastic tail is
3.7 × 10⁶ σ^el_U there and is lent the elastic tensor fraction +0.066), where on
t-peak it was 0.25 % — "the band dominates every one of them" holds at x = 0.01
(166.4× the largest, `fq_scale`) and not for that stand-in at x = 0.30; **(b)
`qe_kf_gev` is no longer the larger of the two quasi-elastic knobs at x = 0.01**
(5.819e−08 against `qe_suppression`'s 1.165e−07 on the default).

## 4. C++: the flip-only failures, each a t-peak pin

Step 1a was the switch ALONE (`rc.hpp` line 1503), rebuilt (1 min 58 s, no
warning) and the whole suite run: **417 cases, 410 passed, 7 failed, 30 of
17 240 342 assertions** (528.15 s; the baseline on `9779edd` was 417 / 417,
17 244 191 assertions, 456.88 s). Every failure is a pin of t-peak physics:

| case (test_rc.cpp) | failing assertion | why it is t-peak's | pin |
|---|---|---|---|
| T8(0) | line 817 ×11, `sq[...] < 1e-30` at the x = 1 edge (3.12e+08 on the default) | the QRT vanishes at the nucleon elastic point only in the t-peak's eta_A range; Eq. (18)'s s-/p-peaks reach it (design §2.8) | `rc_model()` → `tpk()` |
| T8(a') | lines 981–982 ×4 each, 1007–1008 | it divides the model's `tail_sigma_at` by the raw Eq. (38) t-peak quadratures (`polrad_sigma_el_u/_t`, `polrad_sigma_qe_u`) | same |
| B3 | line 1524 ×3, `f < r.hi` (1.6562 < 1.6, 532.5 < 5, …) | "the SIZE of what was unpriced": ×1.28 / ×3.75 / ×1006 are the t-peak's (`phase_B_numbers.md` §B3); on the default ×1.63 / ×534 / ×3.7e6 | every local `RcOptions` → `tpk()` (it compares with `rc_model()`) |
| Q9 | line 1837, QRT(S)/QRT(1) = 0.4721 | the t-peak's Pauli cut (0.6665 on the default) | `on`, `off` → `tpk()` |
| T14 | lines 2609, 2612 (perp/par = −0.5 at rel 1e−12 / 1e−8) | the Q_N = 600 trick is calibrated on the t-peak's tensor fraction (~1e−3 per unit Q_N; ~1e−7 on the default at x = 0.1, so the cancellation eats the tolerance) | via `rc_model()` |
| T8(d) | line 2742, `REQUIRE(mp.sigma_tail_u()[k] == mt.sigma_tail_u()[k])` (0.758 vs 51434.6) | "t-peak+ll leaves the THREE t-peak tables untouched" — `mt` must be the t-peak | via `rc_model()` |
| T19 | line 3301 (`tt.u > mt...u`: equal), line 3324 (0 of 407 cells outside [t-peak, t-peak+ll]) | `mt` is the t-peak the exact tail is compared with | via `rc_model()` |

`rc_model()` (line 347, the base of the sampler's tensor-weights test, T8(0),
T8(a'), B3, B1, B6, T1, T4, T5, T14's θ_S subcase, T18, T8(d), T19 and the
`tail_ratio(cell)` case) now takes `tpk()`, defined on line 346 in the line the old doc
comment held. T8' and T12/T12' (whose comments read Eq. (38)'s σ^el_T, the
t-peak's) take `tpk()` as the design names them, and B1's two TPeak-refusal
lines too. **One more mixed pair was found by the step-1 gate itself**: B6 compares
its local model with `rc_model()`, so once that was the t-peak B6 failed 111
assertions at line 2017 (`m.tail_ratio_at == base.tail_ratio_at` across
models); its two locals take `tpk()` (the band is B6's subject). The
model-independent cases stay on the default and so cover it: T14's channel
table, T17, the TensorAll case, the `y_max` refusal, and `test_rc_pipeline.cpp`
T6 / T18b; plus the two new D5 cases. All 24 + 2 test edits above line 1877
are line-neutral (README, `PHYSICS_CHANNELS.md` and several records cite
`test_rc.cpp` by line there); the new cases are appended.

**Gate after step 1** (the whole suite): **419 cases, 419 passed, 0 skipped,
17 244 225 assertions, 461.28 s** (+4.4 s on the baseline: the default-path
cases cost T14 +7.2 s, T17 +3.6, TensorAll +3.6, T18b +2.1, T6 +2.0 and the
new D5 case 3.6; B3, Q9, T12/T12', T8', B6 are back at their t-peak times).

## 5. Python: the flip-state failures, the split, and the pass/skip comparison

**The flip state** — the library at its D5 state (the switch, the derived CLI
default, the two code consequences of §0) with the Python tests not yet edited
— run on the design's four files (`test_rc.py`, `test_rc_low_x_anchor.py`,
`test_knob_provenance.py`, `test_cli.py`): **8 failed, 579 passed, 145
skipped, 948 s** (they took 210 s on the baseline). Every failure is on the
design's re-take / pin list: `test_rc.py` — the default run's
`meta["rc_tail_model"]` (twice: line 263 and the seven-file chain), the
quasi-elastic s/p column test (its `tpk` base was the default), and the
polrad-full test (its `base` was the default run, so "the exact tail is above
the t-peak" compared the default with itself); `test_knob_provenance.py` — the
`t-peak -> polrad-full` cell on inclusive-6Li (a variation of nothing: the row
says `read`, the hash does not move) and on coherent, tagged-6Li-alpha and
tagged-d-p (polrad-full is the default there now, so the run is accepted while
the row says `refused`). The runtime is the other half of the reason for the
pins: in the flip state `test_knob_provenance.py` took 545 s, `test_rc.py`
253 s, `test_rc_low_x_anchor.py` 120 s and `test_cli.py` 30 s.

**The knob-matrix split (amendment 1).** An explicit `--rc-tail-model t-peak` is a
non-default TAIL knob since D5, refused wherever the tail does not apply, so the
pin cannot be a plain base dict: the design's `RC_TPK = dict(RC,
rc_tail_model="t-peak")` would have refused the base of 57 cells on the tagged
and coherent specs and skipped them in silence. It is a SENTINEL instead
(`python/tests/test_knob_provenance.py:231` `TPK`), resolved per spec in `_run`
(`python/tests/test_knob_provenance.py:507` `inclusive`): t-peak on the inclusive specs — where the tail is
built and the ~13 s polrad-full table is what the pin saves — and the default
everywhere else, where the cell runs exactly as it did before D5. Pinned
(`RC_TPK`): the band knobs, the tail sub-knobs (`fq_scale`,
`tail_tensor_scale`, `c0_shape`, the three quasi-elastic ones and their
`with_qe_tail = false` variants, `with_qe_tail`, `tail_max`) and the two t-peak
refusal cells (`rc_sp_tensor_scale` "at t-peak", the plain `rc_m_lepton`).
On the default: `rc` itself, the two `rc_tail_model` cells (re-based:
`polrad-full -> t-peak+ll`, `polrad-full -> t-peak`), `rc_with_tail` (with an
explicit t-peak its `false` variant would be refused for the wrong reason, the
check's point), `rc_n_eta` (now a polrad-full variant, legal at 64) and the
explicit `polrad-full` cells. Every edit is line-neutral (records cite this file
at 204, 298–301, 312–315, 332–342 and 586–597).

The other pins: `test_rc.py` — t-peak BY NAME where the purpose is t-peak
physics or a comparison with a t-peak number (`test_c0_shape…`'s 1.1043, the
quasi-elastic s/p column test, the bases of the seven-file chain and of the
polrad-full test), for RUNTIME where the subject is model-free
(`meta_distinguishes`, the two `qe_tensor_scale` tests, `a_transfer_frac`,
`optimistic_delta`, the Python-subclass form factor), and the polrad-full test
now reuses the default `rc_run` instead of building its own exact tail;
`test_rc_low_x_anchor.py` — t-peak on the INCLUSIVE builds only (amendment 2:
A5's tagged-alpha builds would refuse it); `test_cli.py`'s untyped-flags run —
`--rc-tail-model t-peak` in both runs (typed, identical); and
`test_bench_djangoh_eltail.py` — wording before line 267 (a record cites it)
and, appended, the label assertion (`"polrad-full (shipped)"`,
`"t-peak (option)"` in the row's generator string).

**The four files after the edits**: 587 passed, 145 skipped, 0 failed, 262.7 s;
per file 32 / 0 (`test_rc.py`, 59.4 → 58.3 s), 5 / 0 (anchor, 2.3 → 2.3 s),
517 / 145 (knob matrix, 511 / 151 before; 144.0 → 197.4 s), 33 / 0
(`test_cli.py`, 4.3 → 4.0 s). With the three relabelled cells mapped to their
new ids (`i5/cell_map.py`), **no test that passed before is skipped, failed or
missing after**, and six that were skipped now pass: `rc_m_lepton … at
polrad-full` and `rc_sp_tensor_scale 0 -> 1 at polrad-full` on coherent,
tagged-6Li-alpha and tagged-d-p — their explicit `polrad-full` base is the
default now, so it runs there, and the `m_lepton` refusal of §0 (ii) is what
makes the first of them pass (by construction: without it the muon variant
would be accepted while its row says `refused`; not run without it).

## 6. Gates, in the design's order (amended), with tallies and wall times

**Baseline on `9779edd`, before any edit:** `./build/lipolgen_tests` 417 cases,
417 passed, 17 244 191 assertions, 456.88 s; full pytest 1463 passed, 155
skipped, 0 failed, 489.65 s (wall 490.67 s); `LIPOLGEN_BENCH=1 … -k bench` 351
passed, 1 skipped (the nnpdf-data wheel is not on this machine), 305.73 s;
strict gate 1241 references / 95 ranges / 7 external / 0 broken / 24
allow-listed (SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0 / 0);
`--records-only` 1423 citations in 75 records, 0 broken (924 checked live, 61
into another record, 110 annotated historical, 129 read by hand, 199 skipped);
`--audit-ranges` 209 range citations on 179 blocks, 0 REFUSED; SPDX 133 files;
`make_report.py --check` OK (64 s); T6 in process 4.66 / 4.48 / 4.46 s; the
A/B hashes of §1.

1. **cmake build** after the switch: 1 min 58 s, no warning; after the code and
   test edits: 1 min 41 s, no warning.
2. **`./build/lipolgen_tests`**: flip-only 410 / 417 (§4); final **419 / 419,
   17 244 225 assertions, 461.28 s** (+4.4 s).
3. **The targeted pytest files** (§5): flip state 8 failed / 579 passed / 145
   skipped (948 s); after the edits **587 passed, 145 skipped, 0 failed**
   (262.7 s), no lost test.
4. **The full pytest**: **1469 passed, 149 skipped, 0 failed, 534.93 s (wall
   535.84 s)** against 1463 / 155 / 0, 489.65 s: +45 s (+9 %), of which the knob
   matrix +52 s and T6's row +11 s; `test_rc.py` −1 s. Per file — the only file
   whose counts moved is the knob matrix, 511 / 151 → 517 / 145, and the six
   gained passes are the six cells of §5; with the renames mapped no test that
   passed before fails, skips or disappears (`i5/cell_map.py`, 42 files):

| file | before passed / skipped / failed | after | time before → after |
|---|---|---|---|
| `test_b1_model.py` | 51 / 0 / 0 | 51 / 0 / 0 | 47.7 → 43.2 s |
| `test_bench_boeglin.py` | 13 / 0 / 0 | 13 / 0 / 0 | 2.1 → 1.8 s |
| `test_bench_bonus.py` | 26 / 0 / 0 | 26 / 0 / 0 | 11.9 → 11.6 s |
| `test_bench_chain_rc.py` | 4 / 0 / 0 | 4 / 0 / 0 | 21.8 → 20.2 s |
| `test_bench_deeps.py` | 16 / 0 / 0 | 16 / 0 / 0 | 2.5 → 2.5 s |
| `test_bench_djangoh_eltail.py` | 18 / 1 / 0 | 18 / 1 / 0 | 44.1 → 41.3 s |
| `test_bench_li6_ff_theory.py` | 18 / 0 / 0 | 18 / 0 / 0 | 2.1 → 2.0 s |
| `test_bench_nuclear.py` | 10 / 0 / 0 | 10 / 0 / 0 | 1.9 → 1.7 s |
| `test_bench_report.py` | 25 / 1 / 0 | 25 / 1 / 0 | 3.3 → 3.2 s |
| `test_bench_rfy.py` | 16 / 0 / 0 | 16 / 0 / 0 | 1.0 → 0.9 s |
| `test_bench_spin_deuteron.py` | 24 / 0 / 0 | 24 / 0 / 0 | 13.9 → 12.5 s |
| `test_bench_static.py` | 18 / 0 / 0 | 18 / 0 / 0 | 5.4 → 5.1 s |
| `test_bench_t1_nucleon.py` | 49 / 1 / 0 | 49 / 1 / 0 | 4.0 → 3.5 s |
| `test_bench_t1_r_world.py` | 51 / 0 / 0 | 51 / 0 / 0 | 6.7 → 6.3 s |
| `test_bench_t4_t6.py` | 18 / 1 / 0 | 18 / 1 / 0 | 8.8 → 19.9 s |
| `test_bench_t5_closures.py` | 19 / 0 / 0 | 19 / 0 / 0 | 6.1 → 5.8 s |
| `test_bench_t5_delta_bound.py` | 23 / 0 / 0 | 23 / 0 / 0 | 10.2 → 9.5 s |
| `test_cli.py` | 33 / 0 / 0 | 33 / 0 / 0 | 4.3 → 4.1 s |
| `test_cluster_config.py` | 27 / 0 / 0 | 27 / 0 / 0 | 3.4 → 3.2 s |
| `test_delta_convention.py` | 6 / 0 / 0 | 6 / 0 / 0 | 0.1 → 0.1 s |
| `test_doc_link_gate.py` | 112 / 0 / 0 | 112 / 0 / 0 | 24.8 → 24.1 s |
| `test_hepmc.py` | 6 / 0 / 0 | 6 / 0 / 0 | 0.1 → 0.1 s |
| `test_hfs.py` | 7 / 0 / 0 | 7 / 0 / 0 | 0.6 → 0.6 s |
| `test_kernel_reference.py` | 11 / 0 / 0 | 11 / 0 / 0 | 0.0 → 0.0 s |
| `test_knob_provenance.py` | 511 / 151 / 0 | 517 / 145 / 0 **(changed)** | 144.0 → 196.0 s |
| `test_li6_unpol_band.py` | 12 / 0 / 0 | 12 / 0 / 0 | 18.6 → 17.7 s |
| `test_li7_rank2.py` | 27 / 0 / 0 | 27 / 0 / 0 | 7.8 → 7.4 s |
| `test_meta_provenance.py` | 8 / 0 / 0 | 8 / 0 / 0 | 2.8 → 2.7 s |
| `test_module.py` | 89 / 0 / 0 | 89 / 0 / 0 | 3.4 → 3.4 s |
| `test_mstw_sf.py` | 25 / 0 / 0 | 25 / 0 / 0 | 0.3 → 0.3 s |
| `test_o5_reach.py` | 18 / 0 / 0 | 18 / 0 / 0 | 0.7 → 0.7 s |
| `test_pipeline.py` | 18 / 0 / 0 | 18 / 0 / 0 | 2.2 → 2.1 s |
| `test_pom_set_band.py` | 9 / 0 / 0 | 9 / 0 / 0 | 4.8 → 4.7 s |
| `test_price_defaults.py` | 1 / 0 / 0 | 1 / 0 / 0 | 3.1 → 3.1 s |
| `test_pzz_mode.py` | 25 / 0 / 0 | 25 / 0 / 0 | 1.2 → 1.2 s |
| `test_rc.py` | 32 / 0 / 0 | 32 / 0 / 0 | 59.4 → 58.2 s |
| `test_rc_low_x_anchor.py` | 5 / 0 / 0 | 5 / 0 / 0 | 2.3 → 2.4 s |
| `test_release_metadata.py` | 8 / 0 / 0 | 8 / 0 / 0 | 0.0 → 0.0 s |
| `test_sampler.py` | 9 / 0 / 0 | 9 / 0 / 0 | 0.1 → 0.1 s |
| `test_sf_backend.py` | 54 / 0 / 0 | 54 / 0 / 0 | 9.1 → 8.9 s |
| `test_spdx_headers.py` | 6 / 0 / 0 | 6 / 0 / 0 | 0.0 → 0.0 s |
| `test_throughput.py` | 5 / 0 / 0 | 5 / 0 / 0 | 2.0 → 1.9 s |

5. **The bench subset** (`LIPOLGEN_BENCH=1 … -k bench`): **350 passed, 1
   skipped, 1 failed, 328.63 s** (wall 329.49 s) — the one failure is
   `test_bench_report.py::test_the_committed_report_is_what_regeneration_gives`,
   whose diff is exactly row 7's two strings: the committed `REPORT.json` still
   read "t-peak (shipped)" until the report was regenerated (item 8), after which
   the same file gave **26 passed** (77.8 s). Row 7's harness itself: 19 / 19, its R pins at
   rel 2e−7 unchanged, the design block's sha256 unchanged, the new label
   assertion passing; T6 (`test_bench_t4_t6.py`) 19 / 19, `runtime_s < 30.0`:

| file | before passed / skipped / failed | after | time before → after |
|---|---|---|---|
| `test_bench_boeglin.py` | 13 / 0 / 0 | 13 / 0 / 0 | 1.9 → 2.0 s |
| `test_bench_bonus.py` | 26 / 0 / 0 | 26 / 0 / 0 | 11.5 → 11.8 s |
| `test_bench_chain_rc.py` | 4 / 0 / 0 | 4 / 0 / 0 | 20.1 → 20.0 s |
| `test_bench_deeps.py` | 16 / 0 / 0 | 16 / 0 / 0 | 2.4 → 2.5 s |
| `test_bench_djangoh_eltail.py` | 19 / 0 / 0 | 19 / 0 / 0 | 131.7 → 131.8 s |
| `test_bench_li6_ff_theory.py` | 18 / 0 / 0 | 18 / 0 / 0 | 2.1 → 1.9 s |
| `test_bench_nuclear.py` | 10 / 0 / 0 | 10 / 0 / 0 | 1.9 → 1.8 s |
| `test_bench_report.py` | 26 / 0 / 0 | 25 / 0 / 1 **(changed)** | 67.1 → 79.1 s |
| `test_bench_rfy.py` | 16 / 0 / 0 | 16 / 0 / 0 | 1.0 → 1.0 s |
| `test_bench_spin_deuteron.py` | 24 / 0 / 0 | 24 / 0 / 0 | 13.0 → 12.5 s |
| `test_bench_static.py` | 18 / 0 / 0 | 18 / 0 / 0 | 5.0 → 4.9 s |
| `test_bench_t1_nucleon.py` | 49 / 1 / 0 | 49 / 1 / 0 | 3.7 → 3.8 s |
| `test_bench_t1_r_world.py` | 51 / 0 / 0 | 51 / 0 / 0 | 6.3 → 6.1 s |
| `test_bench_t4_t6.py` | 19 / 0 / 0 | 19 / 0 / 0 | 22.4 → 33.7 s |
| `test_bench_t5_closures.py` | 19 / 0 / 0 | 19 / 0 / 0 | 5.8 → 5.8 s |
| `test_bench_t5_delta_bound.py` | 23 / 0 / 0 | 23 / 0 / 0 | 9.4 → 9.4 s |

6. **The citation gates** after every edit set and at the end: strict **1241 /
   95 / 7 / 0 / 24** (SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0 /
   0) — unchanged, every edit line-neutral in the cited files, no `--fix`, no
   `--record-ranges`; `--audit-ranges` **209 on 179 blocks, 0 REFUSED** (the
   rewritten PHYSICS rows keep their range evidence); `--records-only` 0 broken
   with this record in place (1430 citations in 76 records: 931 checked against a live file — this record's seven gated citations among them — 61 into another record, 110 annotated historical, 129 read by hand, 199 skipped) — no record citation moved, so no
   `(as of 9779edd)` annotation was needed and the census counts in
   `test_doc_link_gate.py` are untouched (110 annotated, ≤ 200 skipped);
   `test_doc_link_gate.py` 112 passed (23.6 s) with this record in place. Line-neutrality, checked with
   `git diff --numstat`: every file a document or a record cites by line has
   equal insertions and deletions, except the appends — the two C++ test files'
   new cases, `08_channel_priority.md`'s revision-log row,
   `test_bench_djangoh_eltail.py`'s label assertion (below its cited line 267)
   — and the t4 harness's addendum above `run()`, below its last cited line;
   `test_cli.py` gains two lines below its cited ones and
   `test_rc_low_x_anchor.py` (not cited) six.
7. **SPDX**: 133 files, all carry the identifier (no new file but this record,
   which carries it).
8. **The report** (§7) and **`make_report.py --check`**: OK (74.6 s).
9. **The grep audit (design V8)** for "t-peak … default / shipped" and
   "polrad-full … opt-in / not the default" over `include src python README.md
   docs/USAGE.md docs/PHYSICS_CHANNELS.md docs/OPEN_ITEMS_SOLUTIONS.md
   docs/benchmarking validation/benchmarks`: every hit is dated (D5, "until
   2026-10-27"), names `shipped ⁶Li b₁`, or is inside t4's hashed block.

**T6 (amendment 7)**: 16.41 / 16.22 / 16.20 s in process on this machine
(4.66 / 4.48 / 4.46 before), 23 / 23 clauses, PASS; the report's row reads
16.289 s. It stays on the default tail — 13.6 s of margin on its 30 s pin —
so T6 now exercises the shipped `polrad-full` table on its `--rc tensor-band`
file. The fallback (t-peak in T6's rc file) was not needed.

## 7. The report (`make_report.py --all --plots --compare <pre REPORT.json>`)

`LIPOLGEN_BENCH=1 python3 validation/benchmarks/make_report.py --all --plots
--compare <session scratch i5/pre_REPORT.json, sha256 618e3170…>
--out-dir docs/benchmarking` (86 s): **13 pass, 10 fail, 0 blocked — as
before — and `compare: 0 changes`** (no status, no headline number, no sub-row
moved). A structured diff of the two `REPORT.json` finds exactly: row 7's
`generator` ("polrad-full 1.0199 / … ; t-peak (shipped) 0.6118 / …" →
"polrad-full (shipped) 1.0199 / … ; t-peak (option) 0.6118 / …") and `reason`
("polrad-full (two-sided) … the shipped t-peak passes its one-sided CEILING
only …" → "polrad-full (the shipped default since D5; two-sided) … the t-peak
option (shipped until D5) passes its one-sided CEILING only …"); T6's
`runtime_s` 4.608 → 16.289; 21 other `runtime_s` (wall-time noise; row 7
20.207 → 19.987); and the generation stamp (date, revision `9779edd`); `REPORT.md` changes in the
same places (its table carries the runtimes). Nothing else — as the design
expected. **Figures**: `t4_djangoh_rad_noRad.png` changed
in content (polrad-full now in the palette's "default" slot as "LiPolGen
polrad-full (shipped default)", t-peak in `alt2` as "t-peak (option; a lower
bound)"); the other 18 regenerated figures differed from the committed ones
only in their 8-row footer stamp (`i5/fig_diff.py`, pixel by pixel) and were
restored from `HEAD`, as I4 did.

## 8. The check's items: applied, and not

**Applied — every missing site, every "wrong" statement, every underrated risk:**

*Missing sites (15).* `rc.hpp` 30–35 and 52–57 (the header now says Eq. (18) by
default; σ^q_U/σ^el_U on both models) and 241–244 (the band's lead, ×1152 /
×1022 / ×434 on the default beside ×768 / ×681 / ×289); `PHYSICS_CHANNELS.md`
rows 384, 385 (both orderings, both models), 388 (both share triples; the
stand-in's price on both models) and 390 (100 % / 99.9 %); `USAGE.md`
2273–2289 (the paragraph and its correction box name the t-peak tail and print
the default's −1.178978e−07, ×1152), 2250 (both refusal reasons) and 2251 (the
Pauli cut on both models); `OPEN_ITEMS_SOLUTIONS.md` 806–822 (item 3), the
836–844 spread rows (each now carries the default's numbers), 853–855 ("Nothing
clips on the `t-peak` tail", the default's 0.399 / 0.386 / 0.360 % and 0
events) and 907; `__init__.py` 337–341 (both models); `cli.py` 353–355 (both
models), 291–302 (the whole t-peak paragraph through line 302), 458–460 (both
refusal reasons) and 1651–1653 (the banner's clip line names both corners);
`bindings.cpp`'s `qe_kf_gev` docstring (both models; the task's
"docs/USAGE.md:3637-3640" does not exist — the file has 2992 lines — and is the
check's "bindings 3637-3640; USAGE 2251", both done); `test_rc.cpp` 1330, 1355,
1550, 2716, 2727, 2878 (comments and SUBCASE names now say t-peak);
`test_rc.py` 233–234, 342–343 (labelled); the t4 harness's run note (12.0 s,
dated). The design's own inventory (§3.1–§3.5) was applied as well, with the
numbers re-measured instead of copied.

*Wrong (8).* (1) the knob split is per spec (§5); (2) the anchor helper pins
t-peak on the inclusive builds only; (3) no document says "t-peak reproduces
every earlier number" — `USAGE.md` §7b's tail-model bullet says all of what is
true, README's `--rc-tail-model` row, `PHYSICS_CHANNELS.md` row 386 and
`rc.hpp` point (6) its first part: every RC-tail number published before 2026-10-27 without a model label is
a t-peak number, made on the structure-function backends shipped then; event-level
numbers on a tilted plan cannot be regenerated after D4; a HepMC3 file records
`rc_tail` without the rc mode, the tail model or the backend (filed, §9); (4)
the tagged / coherent `meta` and knob-row change is measured and recorded (§1);
(5) the deadlock is removed in code (§0 (i)) and the platform restriction
documented — `rc.hpp` point (5), `rc.cpp`'s refusal comment, `USAGE.md` §7b
(the tail-model row and the `polrad-full` paragraph), the bindings' enum
docstring, `__init__.py`; (6) line-neutrality is not taken for green: the
strict gate, `--audit-ranges` (R5) and `--records-only` (R6) were run after
every edit set (§6); (7) the s/p refusal test is two clauses (§2); (8) the
`rc.hpp` lines are 1225–1228 and were rewritten there; `test_rc_pipeline.cpp`
line 247's reset (no "back to the default" text) now reads
`RcOptions().tail_model`; the stale `test_rc.cpp:374` citations are the range
374–397 and are not D5's (§9).

*Underrated risks (6).* (1) the per-file pass/skip comparison is a gate here
(§5, §6); (2) narrow `long double`: the tagged / coherent band no longer
depends on it, MSVC and Apple arm64 are named, aarch64 Linux is recorded as
unmeasured (§9); (3) the reproduction recipe: as (3) above; (4) the whole-tail
orderings are re-measured on the default and printed beside the t-peak ones
everywhere the check found them (§3); (5) T6: 16.2–16.4 s here, a 13.6 s margin
on its 30 s pin, the runner speed unmeasured (§6, §9); (6) the (x, y) table's
+1.9 … +2.5 % and the 1 + `tail_max` = 11 cap are stated on the default path in
`USAGE.md` §7b (the clip paragraph and the `polrad-full` paragraph),
`PHYSICS_CHANNELS.md` rows 386 and 391 and `08_channel_priority.md` row F;
`rc.hpp`'s `tail_max` comment names the default's clip corner.

**Not applied, with the reason:**
- the design's README sentence "`--rc-tail-model t-peak` reproduces every
  earlier number" and its rc.hpp point (6) "reproduce with tail_model = TPeak" —
  replaced by the truth (amendment 4);
- the design's `RC_TPK` base on every spec — replaced by the per-spec sentinel
  (amendment 1);
- the design's §3.5 optional registry note — the design's own decision text says
  "Registry row 17 untouched; no registry row is added", so the registry is out
  of scope; the text, re-measured, is in §10 for the supervisor;
- T6's fallback (t-peak in its rc file) — not triggered: 16.4 s < 25 s, so T6
  keeps exercising the default tail (amendment 7);
- the optional note in `docs/references/01_radiative.md` (a dated survey, not
  rewritten) and any edit of the run records (history; none of their citations
  broke, §6);
- the design's §6 "label only, not re-measured" for the C0-shape ratios, the
  fq-scale sign change, the `qe_tensor_scale` price and the spreads — they WERE
  re-measured (§3), which the amendment allows and which turned up the two
  reading changes of §3.

## 9. Open points

1. **No tail-table cache — filed as a development item.** The tables are
   member vectors of `RcModel`, built serially in its constructor, one model
   per `Pipeline`; nothing is cached on disk or across pipelines. Every
   inclusive `--rc tensor-band` run, `lg.run(rc=…)` and `Pipeline` now pays
   ~12 s (16 s on `--rc-c0-shape vmc-ft`; the t4 harness's proton with its
   Python form factor ~23 s), and a batch of jobs pays it per job. A cache keyed
   on (ion, beam, scenario, grid, every tail knob) or a threaded build is a
   separate design, not D5.
2. **HepMC3 provenance — filed, out of scope.** A HepMC3 file carries `rc_tail`
   without the rc mode, the tail model or the structure-function backend (the
   npz `meta` carries all three), so the same command line writes a different
   `rc_tail` before and after 2026-10-27 with nothing in the file to say why.
   A run-level attribute would change `docs/HEPMC3_CONVENTION.md`, which T6
   pins by sha256 — a follow-up with its own re-pin, like D4's §8 item 10.
3. **Platforms.** On MSVC and Apple arm64 (`long double` = `double`) the
   default is refused on the inclusive channel with a message that names
   `--rc-tail-model t-peak`; the tagged and coherent band runs (§0 (i)). On
   aarch64 Linux `long double` is a software binary128 and the Eq. (18) build
   (the check counted 67 `long double` uses in `rc.cpp`) is UNMEASURED — it could take minutes.
   CI is x86_64 only.
4. **T6's margin on CI.** 16.2–16.4 s here against `runtime_s < 30.0`; the
   4-core GitHub runner's speed is an extrapolation (`ci.yml`'s header).
5. **Two readings that change on the default (§3), for the author:** the
   `qe_tensor_scale` stand-in at scale 1 is 14.5× the band half-width at
   x = 0.30 (0.25 % on t-peak), and `qe_kf_gev` is the smaller of the two
   quasi-elastic knobs at x = 0.01 (`PHYSICS_CHANNELS.md` row 387's "which is
   why it is the default rather than a band edge" holds on t-peak only). Both
   are recorded where the documents state them; nothing was tuned.
6. **The order with I7.** Every number in §3 is on the toy structure-function
   backends. When D7's grid backends become the default, the Born-dependent ones
   (w_tail − 1, the ΔA_zz, the band's leads, ⟨rc_tail⟩, the clip fractions)
   move again; the ratios of tail pieces (σ^q_U/σ^el_U, the shares,
   (1/6)σ^el_T/σ^el_U, the Pauli cut) do not.
7. **The published w_tail − 1 = 5.18398e−04** (`PHYSICS_CHANNELS.md` row 386,
   t-peak) is 0.03 % above the table-route value 5.18222e−04; its recipe is not
   reproduced (the design's §2.5 note). Recorded, not resolved. The whole
   published table of `OPEN_ITEMS_SOLUTIONS.md` §9 (eight cells, Q² = 2 / 5 / 10)
   sits 0.01–0.26 % above the `t-peak` table route, at q_n = 0 and at q_n = 1
   alike (F5, §12): a `t-peak` table, and the residual is not the tensor term.
8. **Two stale citations not D5's**: the range 374–397 of `tests/test_rc.cpp`
   cited from `docs/benchmarking/02_data_nucleon_deuteron.md` line 299 and
   `docs/references/01_nucleon-deuteron-data.md` line 429, stale since before
   this run (the design's §7 and its check note them; no D5 edit moved them).
9. **For the supervisor, DECISIONS.md D5.** The entry is the one-line decision;
   if it is expanded with the design's §12 text, its clause "and reproduces
   with `--rc-tail-model t-peak`" must not go in (amendment 4); this record's
   §8 (3) has the sentence to use instead.

## 10. The registry (row 17) — text only, not applied

The design's decision text keeps the registry out of scope ("Registry row 17
untouched; no registry row is added"), so nothing was written there. Its
optional note, re-measured, for `docs/open_items/run_2026-09-03/STATUS.md`
(append after the last line, so nothing shifts):

> 2026-10-29 — D5 taken 2026-10-27 (`../run_2026-10-27/DECISIONS.md`) and
> implemented (`../run_2026-10-27/phase_I5_polradfull.md`): `PolradFull` is
> the shipped `rc_tail`. Not a registry row. Row 17's ×768 / ×681 / ×289 was
> measured on the t-peak tail; on the default exact tail the band leads the
> whole tail at x = 0.01, Q² = 5 by ×1152 / ×1022 / ×434 (its ΔA_zz
> −1.178978e−07 against −1.769797e−07), so the row's ordering survives every
> anchor with more margin.

## 11. The files, and the provenance

**Library** (line-neutral): `include/lipolgen/rc.hpp` (the switch at line 1503;
23 comment edits — the header, HONEST FLAGS, "FIVE THINGS THIS IS NOT", the
`RcTailModel` enumerators and the `PolradFull` points (3), (5), (6), the
`TPeakPlusLL` point (3), and the `qe_tensor_scale`, `qe_kf_gev`,
`sp_tensor_scale`, `tail_max` and clip-fraction field comments);
`src/core/rc.cpp` (the gated preconditions and their comment; the tagged
exclusion reason, written to `meta["rc_exclusion_reason"]`, names both models'
shares; one comment); `src/core/pipeline.cpp` (the gated preconditions, the
`m_lepton` tail sub-knob, the `m_lepton` refusal's text, the `rc_tail_model`
knob row's reason, one comment). **Bindings and Python package**:
`python/bindings.cpp` (two comments, the `RcTailModel` enum, `tail_model`,
`sp_tensor_scale` and `qe_kf_gev` docstrings); `python/lipolgen/cli.py` (the
derived default; the `--rc`, `--rc-c0-shape`, `--rc-tail-model`,
`--rc-qe-tensor-scale`, `--rc-sp-tensor-scale` help; the banner's clip line, QE
shares and both model paragraphs); `python/lipolgen/__init__.py`
(`RC_C0_SHAPES`, `RC_TAIL_MODELS` comments). **Tests**: `tests/test_rc.cpp`,
`tests/test_rc_pipeline.cpp` (§4, §2); `python/tests/test_rc.py`,
`test_knob_provenance.py`, `test_rc_low_x_anchor.py`, `test_cli.py`,
`test_bench_djangoh_eltail.py` (§5, §2). **Harnesses**:
`validation/benchmarks/t4_djangoh_rad_noRad.py` (labels, plot slots, the run
note, the D5 addendum above `run()` and its pointer in the ADDENDA header — the
hashed block untouched, `design_block_sha256()` unchanged);
`validation/benchmarks/_plots.py` (the `alt2` slot's comment).
**Documents** (line-neutral but one appended log row): `README.md` (the
`--rc-tail-model` and `--rc-c0-shape` rows; the "still open" sentence);
`docs/USAGE.md` §7b (the tagged-omission box, the mandatory-band table's
fq-scale, c0-shape, tail-model, qe-suppression, qe-tensor-scale,
sp-tensor-scale and `qe_kf_gev` rows, the C0 box, the whole-tail paragraph and
its correction box, the honest-flags bullets, the banner / clip paragraph, the
cost paragraph); `docs/PHYSICS_CHANNELS.md` (rows 384–388, 390, 391, 393 and
the §13 bullets at 747–751, 815 and 840–841); `docs/OPEN_ITEMS_SOLUTIONS.md`
(§9 and the summary rows, dated notes inside existing lines);
`docs/benchmarking/BENCHMARK_PLAN.md` (row 7's status cell — read by the report
only for rows without a harness, so the report does not see it);
`docs/benchmarking/08_channel_priority.md` (row F; a revision-log row
appended); `validation/benchmarks/README.md` (rows 7 and T6);
`docs/benchmarking/REPORT.md`, `REPORT.json`, `figures/t4_djangoh_rad_noRad.png`
(regenerated, §7); this record; `STATUS.md`'s I5 row. **Not touched**: the
hashed design block, `docs/HEPMC3_CONVENTION.md`, `validation/reference/*.json`,
`docs/theory/SPIN32_FINITE_GAMMA.md`, the run records, the registry, the
sibling repository, every tolerance, census and `ASOF_READ_BY_HAND` entry.

**Session scratch `i5/`**: `edits/` (`s1a_flip.py`, `s1b_code.py`,
`s1c_cpp_tests.py`, `s1d_cpp_new.py`, `s1e_b6.py`, `s2_*.py`, `s3_*.py`,
`s4_harness.py`, `s5_*.py`, `s5b_extra.py`, `apply.py`); `meas/` (`run_pre.sh`
and its output `pre.out`, `m_points.py`, `m_extra.py` / `extra_pre.out`,
`m_qrt_cut.py`, `m_c0cost.py` / `c0cost.out`, `t6_time.py` / `t6_post.out`,
`npz_hash.py`); `ab/` (the six 5000-event npz files and their logs);
`base_cpp.txt`, `flip_cpp.txt`, `rc_subset_cpp*.txt`, `step1_cpp.txt`;
`base_full.xml` / `.txt`, `base_bench.xml`, `flip_py_targeted.xml`,
`post_py_targeted.xml`, `post_full.xml`, `post_bench.xml` and their `.txt`;
`junit_counts.py`, `cell_map.py`, `fig_diff.py`; `gate_*.txt`;
`report_regen.txt`, `report_check_*.txt`; `pre_REPORT.json` (sha256
`618e3170b9ed2246…`, the committed one). Nothing was committed, published or
sent; the sibling repository was not opened.

## 12. Verification

The adversarial verification of the tree I5 left (session scratch `v5/`).
Its verdict: "SOUND, with label fixes needed". The one-place switch, the
derived CLI default, the two amendment-5 code changes, the test split, the
gates, the report and row 7 check out when run, and every number of this
record the verifier re-ran reproduces. It found one medium and seven low
findings, all in labels, messages or wording; no number, tolerance or verdict
in code was wrong. Worker F5 fixed them on 2026-10-29 (session scratch `f5/`).
Each finding was re-checked against the tree before anything was touched, and
none was refuted; three of the verifier's own figures were corrected on the way
(V-1, V-3). F5 is line-neutral in every file a record or a gated document cites
by line; the C++ checks it adds sit at the end of the two D5 cases, at the end
of their files. The verifier labelled no finding with a lens; the lens column
is F5's assignment.

| # | lens | finding | F5's re-check | what was done |
|---|---|---|---|---|
| V-1 (medium) | honesty | `OPEN_ITEMS_SOLUTIONS.md` §9's `w_tail − 1` table (x = 0.01 / 0.10 / 0.30 × Q² = 2 / 5 / 10) and its headline ("five to eight orders of magnitude smaller than at HERMES … entirely `Y₊` and the Born's 1/Q⁴") were unlabelled `t-peak` data, and the design, its check and I5 had all missed them. | Confirmed (`f5/meas/f5_wtail.py`, both models, table route). On `t-peak` the route reproduces every printed cell to 0.01–0.26 %; on the default, at q_n = 0, it gives 1.71661e−03 / 6.62870e−04 / 2.19742e−03 at x = 0.01, 2.52430e−03 / 1.95963e−04 / 2.51627e−05 at 0.10 and — / 4.88356e−04 / 5.04628e−05 at 0.30 (x = 0.30, Q² = 2 lies below the table's y support, as the printed "—" says), the verifier's figures to every digit. Two of its readings were corrected: the ratio to the printed values is ×4.9 to ×1.2e4 at x ≥ 0.10 and at Q² = 2, not ×19 to ×1.2e4 (x = 0.10, Q² = 10 is ×4.85), and the default's tail at those cells is 94.4–100 % quasi-elastic (`f5_qeshare.py`), not 99.8–100 % (that is Q² = 5 alone). | Line 615 labelled "on `t-peak` (the default until D5)". A dated D5 note appended to line 630 gives the default's eight values, the ratios (×1.02–×1.28 at x = 0.01, Q² ≥ 5; ×4.9 to ×1.2e4 elsewhere), the factor 7.1 between 2.5e−03 (x = 0.10, Q² = 2) and the HERMES-like 1.8 %, and says that the "five to eight orders" comparison and the "entirely `Y₊`" reading are the `t-peak`'s and were not re-made on `polrad-full`. Line-neutral. |
| V-2 (low) | code | Under the new default, `--rc-sp-tensor-scale 1` on a tagged or coherent `--rc tensor-band` run, or at `with_tail = false`, was refused with `polrad-full`'s "would double-count a term that ran", where no tail runs; `docs/CONVENTIONS.md` requires a refusal to name its cause. The C++ path (`RcModel`'s constructor) gave the same text. | Confirmed: the CLI on `tagged-6Li-alpha` and `coherent` printed the double-count text (I5 build), and `RcModel`'s constructor reads the same way. | `validate()`'s `sp_tensor_scale` clause now binds only where a tail is built (`src/core/pipeline.cpp:640-641`, `rc_tail_applies`), so off the tail the tail sub-knob rule refuses `--rc-sp-tensor-scale` with its own reason; `RcModel` keeps the double-count text where a tail is built and says "no tail is BUILT on this run" where none is (`src/core/rc.cpp:1483-1494`; line 1481, a pinned use site, untouched; the first string reflowed from three lines to two to keep the count). No accept/refuse outcome moves. Gated by checks appended to the two D5 cases: in `tests/test_rc_pipeline.cpp`, `validate()` on tagged, coherent and inclusive `with_tail = false` names the RADIATIVE TAIL rule, and the inclusive tail still gives "double-count"; in `tests/test_rc.cpp`, a new SUBCASE does the same on the C++ path ("no tail is BUILT"). By construction both fail on the I5 build, whose CLI runs above printed the double-count text; they were not run against it. After the fix the same CLI runs print "--rc-sp-tensor-scale is a knob on the RADIATIVE TAIL, and rc_tail == 1 exactly on this run", and the inclusive runs keep their two reasons (the default's double count; t-peak's "never paid"). |
| V-3 (low) | honesty | More unlabelled `t-peak` numbers in living text: `rc.hpp`'s `qe_tensor_scale` point (3) (σ^el_T/σ^el_U on both C0 edges, the factors 35 / 67 / 1.5 and 67 / ~255), `PHYSICS_CHANNELS.md` rows 388 (the same triple) and 389 (the C0 (1/6)σ^el_T/σ^el_U pairs, σ^el_U ×1.011 / ×1.022 / ×1.104 / ×3.03, the 1.33 s / 0.146 s cost, the `fq_scale` and `tail_tensor_scale` spreads), `OPEN_ITEMS_SOLUTIONS.md` line 87 and the budget ladder at line 848. | Confirmed (`f5/meas/f5_points.py`; every `t-peak` number reproduces). The default gives σ^el_T/σ^el_U = −3.0924e−03 / +2.6737e−04 / +6.5680e−02 (`ho`) and −3.3277e−03 / −9.8212e−04 / +6.7761e−02 (`vmc-ft`); the verifier's +2.673e−04, −9.822e−04 and +6.777e−02 are off in the last digit. Factors below the deuteron: 34 / 233 / 1.8, and 63 at x = 0.10 on `vmc-ft`, taken against the re-driven +0.0623 like the published 67 / 255 (66.6 / 255.1). σ^el_U(`vmc-ft`)/σ^el_U(`ho`) = ×1.011 / ×1.023 / ×1.105 / ×2.761. The default ladder: `fq_scale` 8.165e−07 ≫ `qe_tensor_scale` 2.650e−07 > `qe_suppression` 1.165e−07 > `qe_kf_gev` 5.819e−08 > `tail_tensor_scale` 3.656e−08 > `c0_shape` 3.467e−08, the band 166.4× the largest. Construction (`--rc tensor-band` minus `--rc off`): 11.9 / 16.2 s on `ho` / `vmc-ft` (default), 0.165 / 1.68 s (`t-peak`). | Each site labelled `t-peak`, with the default's value printed beside it, appended inside existing lines: `rc.hpp` point (3) (13 lines → 13), `PHYSICS_CHANNELS.md` rows 388 and 389 (four places; the spreads' row also says the pair is as published on the withdrawn ×3.253983 `A_zz`, which neither spread depends on), `OPEN_ITEMS_SOLUTIONS.md` lines 87 and 848, and line 63, which carries the same triple and was found by grep. README's I5 text gave the default's factor as ×232, against the +0.062 printed beside it, while the 67 / 255 next to it are against +0.0623; it now reads ×233, on the same basis as every other site. |
| V-4 (low) | docs | `constants.hpp`'s `M_ELECTRON` comment said `validate()` refuses a non-default `m_lepton` "ON THOSE TWO ONLY" (the `t-peak` tails), false since I5 made `m_lepton` a tail sub-knob; `PHYSICS_CHANNELS.md` row 393's three-way row left out the same case. | Confirmed by running `validate()` with the muon mass: refused by the tail sub-knob rule on `tagged-alpha`, on `coherent` and on inclusive at `with_tail = false`; accepted on inclusive with the tail. | Reworded, line-neutral: "on those two, and (D5) on ANY model where no tail is built"; row 393 gains the same clause ("since D5, on every model wherever no tail is built"). |
| V-5 (low) | test | The knob cell `rc_sp_tensor_scale: 0 -> 1 at t-peak (an option since D5)` takes `RC_TPK`, which resolves to the default off the inclusive channel, so on `coherent`, `tagged-6Li-alpha` and `tagged-d-p` it ran on `polrad-full` under a `t-peak` label. It passed on status. | Confirmed by reading `_run`'s resolution of the `TPK` sentinel and the cell. | Relabelled, line-neutral: "0 -> 1 at t-peak (inclusive specs; the default elsewhere)". Its test ids change on every spec; the per-file counts do not (the gates below). |
| V-6 (low) | honesty | §3 said "Every t-peak number the documents publish reproduces to every printed digit", which its own table contradicts twice. | Confirmed: the published `w_tail − 1` sits +0.034 / +0.254 / +0.009 % above the table route at x = 0.01 / 0.10 / 0.30 (the §3 cell said 0.03 % for all three), and the published C0 ΔA_zz ratios reproduce only on the withdrawn `A_zz`: `f5/meas/f5_c0azz.py` (the verifier's script) gives −8.2 / +0.6 / +36.9 % with ×3.253983 and +18.7 / −4.0 / +90.7 % on the shipped ⁶Li b₁. | §3 names both exceptions; the §3 cell gives the three percentages; §9 item 7 is extended to the whole published table. |
| V-7 (low) | docs | `PHYSICS_CHANNELS.md` row 386 and `USAGE.md` §7b's `polrad-full` paragraph read "1.0004 / 1.0010 / 0.9995 / 0.9866 from the quadrature, +1.9 … +2.5 % through the (x, y) table", which reads as R − 1 through the table; R through the table is 1.0199 / 1.0261 / 1.0219 / 1.0098. | Confirmed by arithmetic: table over quadrature is 1.0195 / 1.0251 / 1.0224 / 1.0235, so +1.9 … +2.5 % is the table's excess, while R − 1 through the table is +1.0 … +2.6 %. A grep finds no third site: `08_channel_priority.md` row F, `BENCHMARK_PLAN.md` row 7 and the t4 harness already say "adds". | Both sites now read "the (x, y) table adds +1.9 … +2.5 % — R = 1.0199 / 1.0261 / 1.0219 / 1.0098 through it". Line-neutral. |
| V-8 (low) | docs | `PHYSICS_CHANNELS.md` row 393 called t-peak's 0.106 s "the one t-peak construction-time number, and the earlier 0.134 s / 0.16 s pair is withdrawn", while I5's new texts (`rc.hpp` `PolradFull` point (5), the `--rc-tail-model` help, README's row) quote 0.16 s. | Confirmed; `OPEN_ITEMS_SOLUTIONS.md`'s cost table (line 828) says the same, and `USAGE.md` §7b's cost paragraph quotes 0.106 s. Re-measured the same way as the 12.0 s: 0.165 s median on `t-peak` (0.146–0.180 s, five runs) against 11.9 s on the default. | The 0.106 s dated 2026-09-02 ("the one number of that date") and the 2026-10-29 0.16 s added beside the 12.0 s at all three sites (`PHYSICS_CHANNELS.md` row 393, `OPEN_ITEMS_SOLUTIONS.md` line 828, `USAGE.md` line 2582); I5's texts keep 0.16 s, the like-for-like partner of the 12.0 s. |

**One more, found by F5** while running the D5 case: its comment said the
default clips "0 on this grid", while its own MESSAGE prints 0.00925 of the
40 × 24 test grid's nodes on the default (0 on `t-peak`; the I5 and V5 logs
print the same). The comment now says 0.925 %, line-neutral, at the end of the
file. Nothing asserts on it.

**Gates after F5** (the tree F5 leaves; uncommitted):

- **cmake build**: after the code edits 2 min 13 s (near-full: `constants.hpp`
  changed), after the `rc.hpp` comment 1 min 37 s, after the D5-case comment
  10 s; no warning, no error.
- **`./build/lipolgen_tests`**: 419 cases, 419 passed, 0 skipped, **17 244 240
  assertions** (+15 on I5's 17 244 225: the 4 checks appended to the pipeline
  D5 case, the 4 of the new SUBCASE, and the D5 case's 7 top-level checks,
  which doctest runs once more for the new SUBCASE), 446 s wall.
- **The full pytest**: **1469 passed, 149 skipped, 0 failed, 532.8 s** (I5:
  1469 / 149 / 0, 534.9 s). The per-file counts are I5's in all 42 files
  (`f5/junit_counts.py` against I5's `post_full.xml`); per test, the only
  differences are the 7 ids of the relabelled cell (V-5), which map one to one
  onto the old ids with the same outcome (4 passed, 3 skipped;
  `f5/cell_map_f5.py`: 0 lost, 0 gained, 0 unmapped).
- **The bench subset** (`LIPOLGEN_BENCH=1 … -k bench`): **351 passed, 1 skipped
  (the nnpdf-data wheel, as before), 0 failed, 327.4 s**, per test identical to
  the pre-I5 baseline (`i5/base_bench.xml`: 0 outcome changes). Against I5's
  run, the report test that failed until I5 regenerated the report now passes.
- **`make_report.py --check`**: OK (75 s). F5 changed no harness, no `REPORT.*`
  and no figure.
- **The citation gates**, after the document edits and again with this section
  in place: strict **1241 / 95 / 7 / 0 / 24** (SPIN32 19 / 115 / 0 / 0 / 10;
  PYTHIA_BRIDGE 0 / 0 / 8 / 0 / 0), unchanged, with no `--fix` and no
  `--record-ranges`; `--audit-ranges` **209 on 179 blocks, 0 REFUSED** (the R5
  evidence of the edited rows, row 393's `qe_suppression` for its USAGE §7b range
  and row 389's `tests/test_rc.cpp` for its `rc.hpp` range, sits after every F5
  insertion, so their 200-character windows are unchanged); `--records-only`
  **1432 citations in 76 records, 0 broken** (933 checked against a live file,
  this section's two among them; 61 into another record, 110 annotated
  historical, 129 read by hand, 199 skipped). No `(as of …)` annotation was needed and no census was
  re-pinned; `test_doc_link_gate.py` 112 passed (23.3 s) with this section in
  place.
- **SPDX**: 133 files. **Line neutrality** (`git diff --numstat` against the I5
  tree's copies in `f5/pre_f5/`): F5 adds lines only at the ends of
  `tests/test_rc.cpp` (+22) and `tests/test_rc_pipeline.cpp` (+21); every other
  F5 edit keeps its file's line count, and line 1481 of `src/core/rc.cpp` (a
  pinned use site) is byte-identical.

Scratch `f5/`: `edits/` (`apply.py`; `step_a_code.py`, `step_a2_comment.py`,
`step_a3_d5comment.py`, `step_b_docs.py`, `step_c_record.py`,
`step_d_record.py`, `step_e_status.py`, every edit an exact old/new pair) and `pre_f5/` (each edited
file as I5 left it); `meas/` (`f5_wtail.py`, `f5_qeshare.py`, `f5_points.py`,
`f5_c0azz.py` and their `.out`); `build_a.txt`, `build_b.txt`, `cpp.txt`,
`full.xml` / `.txt`, `bench.xml` / `.txt`, `gate_*.txt`, `report_check.txt`;
`cell_map_f5.py`, `junit_counts.py`. Nothing was committed, published or sent,
and the sibling repository was not opened.
