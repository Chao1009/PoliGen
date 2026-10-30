<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-27 — phase I6: the VMC AV18+UX cluster waves become the shipped run default (decision D6)

Worker I6, 2026-10-29, on `claude/run-2026-10-26` at `cb75ec7` (D4 and D5
committed; clean tree apart from the supervisor's `STATUS.md` edit that set
this phase RUNNING). Uncommitted: the supervising session commits. Decision:
`DECISIONS.md` D6. Design: `plan_D6_cluster_wave_default.md`, AMENDED by its
adversarial check (the plan's last section; session scratch
`checks_1027/D6.json`) as the task's eleven amendments name. Every
measurement below was taken on this machine (4 cores; the heavy gates one at a
time) with `source env.sh`; the scripts are in session scratch
`i6/meas/`, the edit scripts — every edit an exact old/new pair, applied by
`i6/edits/apply.py`, which refuses a missing or ambiguous old string and prints
each file's line-count change — in `i6/edits/`.

Conventions: **pre** = `cb75ec7` before any I6 edit; **default** = the run
default after this phase (`vmc`); **hulthen** = `--cluster-wave hulthen`,
the option. Every sample-level number this record re-measures (event
weights, clamp counts, tag fractions at a fixed seed) is **on the
structure-function backends shipped at the time — the toys**: D7 moves the
event stream again. File and line references are written as plain text
("src/core/pipeline.cpp line 497", or file:line without backticks) and
describe the tree this record leaves; none of them is a gated citation.

## 0. What D6 changed, in one paragraph

The default lives in ONE place, include/lipolgen/cluster.hpp line 419,
`inline constexpr ClusterWaveSource CLUSTER_WAVE_DEFAULT =
ClusterWaveSource::VmcAV18;`, and `PipelineConfig::cluster_wave` reads it
(include/lipolgen/pipeline.hpp line 823) — hence `lipolgen-run`, `lg.run` and
`lg.make_config`. Every second definition reads it too: the bindings export it
with `cluster_wave_name` (python/bindings.cpp lines 5163-5169, at the end of
the module), the CLI's `DEFAULTS["cluster_wave"]` is `None` and its help is
built from it (python/lipolgen/cli.py lines 225-227 and 805), the example
generator's `--cluster-wave` defaults to it, and `validation/price_defaults.py`
names its tagged default rows from it. The design's DD-1 is implemented as
written: the channel FACTORIES (`li6_alpha_channel`, `li7_alpha_channel`,
`deuteron_channel`) and `BreakupOptions::source` keep `Hulthen` as their own
argument default, so `validation/reference/tagged.json`, its port gates and every
factory-level pin are untouched (`git diff --stat -- validation/reference`
empty). Three code consequences: `validate()` refuses a `cluster_wave` MOVED
from the default off the tagged channels (src/core/pipeline.cpp line 497;
otherwise every default inclusive and coherent run would throw), `make_config`
refuses an EXPLICIT `cluster_wave` of either value there (python/lipolgen/__init__.py
line 556 and its helper at line 877; registry row 11's protection, which a value
rule can no longer give), and the off-channel knob row withholds the value — it
reads "not read on channel inclusive", never a bare `vmc` (src/core/pipeline.cpp
line 2806). Provenance: `meta["cluster_wave"]` on every generated npz (through the
knob table, python/bindings.cpp line 641) and a run-banner line on every tagged
run naming the wave, its label, its P_L and the struck cluster's polarizations;
a HepMC3 file carries no such trace (its `GenRunInfo` has the tool name and
version only, src/hepmc/hepmc_writer.cpp) — filed, §9. A missing table now
throws with the escape named (src/core/cluster.cpp lines 79-89). **What moves**
is every default tagged run: the ⁶Li α–d relative wave (P_D 0.0867 → 0.019355),
its embedded deuteron and T1 breakup deuteron (eff_pol 0.9325 → 0.913600), the
⁷Li α–t wave and the d–p control (P_D 0.045 → 0.057600). **Nothing else
moves**: the explicit-wave runs are bit for bit, the inclusive and coherent
columns are bit for bit, no constant moved — `LI6_CLUSTER_POLARIZATION` =
0.811228 and `LI6_B1_RANK2_TRANSFER` = 0.921947 are KEPT, pending the author's
answer (§9 item 1) — and no benchmark verdict moved (§6).

## 1. The decisive checks: the A/B regression and the physics lens

**A/B (the design's §10 item 1, amended by the check's wrong-item 12).** The CLI, 20 000
events, seeds 1234 and 20261029, on the three tagged channels (⁶Li α
`tensor-thirds`, ⁷Li α `helicity-flip`, d–p `tensor-thirds`) at each of
default / `hulthen` / `vmc`, and the inclusive and coherent defaults — 22 files
made on the pre build, 22 on the I6 build (`i6/meas/ab_runs.sh`), compared by a
sha256 over every array column (`i6/meas/ab_compare.py`, log
`i6/logs/ab_compare.txt`): **every one of the 28 pairs is column-identical**.

| pre file | post file | columns | `meta` |
|---|---|---|---|
| `--cluster-wave vmc` | default | identical (6 / 6) | the new `cluster_wave` key; the `cluster_wave` knob row's `at_default` false → **true** and its reason; on the lithium tags the `cluster_vmc_mc_sigma` row's reason (it says "READ on the default lithium tags (D6)") |
| `--cluster-wave hulthen` | `--cluster-wave hulthen` | identical (6 / 6) | `cluster_wave` key; the row's `at_default` true → false |
| default (then Hulthén) | `--cluster-wave hulthen` | identical (6 / 6) | the same |
| `--cluster-wave vmc` | `--cluster-wave vmc` | identical (6 / 6) | as row 1 |
| inclusive / coherent default | the same | identical (4 / 4) | `cluster_wave` key; the row's value `hulthen` → "not read on channel inclusive" / "… coherent-6Li" (status `refused`, at default) |

So the pre-I6 `vmc` run equals the post-I6 default run in every column, and
the only flip in the knob row there is `at_default` (the check's wrong-item 12).

**Physics lens (the design's §10 item 2).** The design's own extraction
(`d6scripts/tagged_prices.py`, a copy with a third pass for the default:
`i6/meas/m_tagged_prices.py`, N = 200 000, seed 20261026, config 1, YR high
acceptance) on the I6 build: the `hulthen` and `vmc` rows are **identical in
every printed field** to the design's table (`runs/tagged_prices_200k.jsonl`),
and the default rows equal the `vmc` rows. ⁶Li α: A_zz^tag −1.2069113 /
−0.9532827 (hulthen) → **−0.5191289 / −0.1992630** (default) at k = 0.19793 /
0.30115 GeV; RP tag fraction 0.02466 ± 0.00035 → 0.03422 ± 0.00041; vector /
tensor dilution 0.8699432 / 0.9219490 → 0.9709661 / 0.9825759; the label reads
"6Li alpha-tag (embedded d, VMC AV18+UX)", `dis_target.eff_pol_p` =
0.9136001620. d–p: −0.8763449 / −0.5478468 → −0.9624645 / −0.6680525; RP
0.02396 → 0.02107; dilutions 0.9324961 / 0.9594889 → 0.9135960 / 0.9481463.
⁷Li α: ⟨k⟩ 0.28852 → 0.18608 GeV, RP 0.97386 → 0.99839. σ unchanged on all
three. The T1 breakup deuteron of the default ⁶Li run is the AV18 one, P_D =
`deuteron_av18_p_d()` (0.0575998919874, the design's §10 item 2; the new C++ D6
case pins the equality).
`LI6_CLUSTER_POLARIZATION`, `LI6_B1_RANK2_TRANSFER` and every
`validation/reference/*.json` are byte-unchanged.

## 2. Every number that moved or was re-taken, before → after

**Tests** (each re-take with a dated comment naming D6):

| file:line | what | before | after |
|---|---|---|---|
| tests/test_pipeline.cpp:1291-1337 | "cluster_wave is refused on the channels that never read it" (RETAKE) | `Hulthen` — the default — legal everywhere; `VmcAV18` refused off the tags | the default (`CLUSTER_WAVE_DEFAULT` = VmcAV18) legal everywhere; `Hulthen` (moved from the default) refused off the tags, legal on them |
| tests/test_pipeline.cpp:1652-1777 (new) | three D6 cases: the one home and the enum values; the factories' and `BreakupOptions`' defaults still analytic; a default tagged run reads VMC (label, AV18 embedded deuteron, T1 breakup source and P_D); the knob rows (`cluster_beta` READ on ⁷Li at T1, not at T0, not on ⁶Li or d–p; `p_d` not read on the default; `cluster_vmc_mc_sigma` READ on the lithium tags, refused on d–p; the off-channel row's value "not read on channel …", `meta_value` = value, no "vmc" in it); the MC band accepted on the default lithium tags and refused on Hulthén, d–p and inclusive | — | new |
| python/tests/test_cli.py:131-134 | the `--cluster-wave` default (RETAKE, line-neutral) | `resolve(...)["cluster_wave"] == "hulthen"` | `is None`: the config keeps the C++ `CLUSTER_WAVE_DEFAULT` (an explicit `hulthen` is exercised by the banner test below) |
| python/tests/test_cli.py:349-391 (new, at the END of the file) | an explicit `--cluster-wave` (either value) on inclusive / coherent is a one-line refusal; the tagged banner names the wave; `meta["cluster_wave"]` is `vmc` / `hulthen` on the tags and "not read on channel inclusive" off them | — | new |
| python/tests/test_module.py:274-284 | `make_config` default (RETAKE, line-neutral) | `Hulthen` | `CLUSTER_WAVE_DEFAULT` == `VmcAV18`; `cluster_wave="hulthen"` gives `Hulthen` |
| python/tests/test_module.py:596-616 | the MC band | refused through a DEFAULT config | the refusal half on `cluster_wave="hulthen"` by name; NEW the positive half: a default config accepts the band on ⁶Li and ⁷Li α and refuses it on d–p |
| python/tests/test_module.py:729-756 | C5.5b refusal matrix (RETAKE) | Hulthén legal everywhere, VMC refused off the tags | the default legal everywhere, Hulthén refused off the tags; an explicit `make_config(cluster_wave=…)` of `vmc`, `hulthen` or the enum refused off the tags (`match="^PipelineConfig: cluster_wave = '(vmc\|hulthen)' was asked"`), accepted on them |
| python/tests/test_knob_provenance.py:236, 331-342 | the seven tagged-cluster cells (RETAKE, §4) | `hulthen -> vmc`; `0 -> 1 on hulthen` (no base); `0.30 -> 0.33`, `0.0867 -> 0.10` (no base); their `on vmc` twins | `vmc -> hulthen`; every Hulthén cell carries `HUL` in base AND var; the plain cells run `on the default (vmc)` |
| python/tests/test_pipeline.py:99-108 (new) | the tag fraction on the default wave | — | 0.03313 at 10⁵ events, seed 1, against the model's uniform-M 0.033810 (−2.0 %) at the Hulthén test's own rel 0.20 |
| python/tests/test_bench_boeglin.py:249-252 | the REPORT line (text only, line-neutral) | "hulthen (default) chi2 = 21.98/21 …" | "hulthen (the design block's verdict wave) chi2 = 21.98/21 …" |

**Made explicit, no number moved** (the step-1 gate, §5, proves it): the
`pipe()` helper of tests/test_pipeline.cpp takes the wave in its cache key
(lines 42-93; the inclusive, coherent, threaded and throughput call sites keep
the default — the check's step-1 trap) and the three port comparisons with the
Hulthén prototype pass `kHulthen` (lines 408, 450 and 498): the tag fractions
still print 0.024885 / 0.305123 / 0.25432 and 0.96623 … 0.99441, the A_zz^tag
closure 17.84/9 (YR) and 15.15/11 (tagging); tests/test_fsi.cpp lines 48 (the
prototype's pure-S wave), 206 (the S+D row, relabelled) and 394 (the D3 pair,
whose TEST_CASE line 384 — an unnamed citation of `PHYSICS_CHANNELS.md` — stays
put); tests/test_tagged.cpp lines 73-85 (the polligen port-gate models) and 341,
514, 834 (the P_D = 0 limits); `python/tests/test_pipeline.py`'s new
`tagged_run_hulthen` fixture for the documented 0.0248; `test_meta_provenance.py`'s
`fsi_runs` fixture. The comments of T25 – T27 (tests/test_tagged.cpp lines
1574, 1617, 1667) now say "the factory default" and "the Hulthén path". T25 – T27
themselves keep the factory-default calls on purpose: with DD-1 they pin that
the factories ARE still the analytic family.

**Published numbers, re-measured** (toy SFs; scripts `i6/meas/m_fsi.py`,
`m_fsi_pct.py`, `m_rc.py`, `m_tagfrac.py`, run on the pre build and again on
the I6 build — identical JSON with the wave named, `meta["cluster_wave"]` aside):

| quantity (recipe) | published | `hulthen`, 2026-10-29 | **default `vmc`** | where |
|---|---|---|---|---|
| ⁶Li α D3 FSI, 20 000 events, seed 1234, `tensor-thirds` 0.7 / 0.6 / 0.7, 40 mb: Σw cluster / nucleon | 10419.07 / 11632.10 (2026-09-04, pre-S–D-fix sample) | 10403.75 / 11584.83 | **15989.35 / 16428.71** | fsi.hpp 128-150; USAGE §3; CONVENTIONS; PHYSICS_CHANNELS row 344; OPEN_ITEMS §7.1, rows 61, 85; README row; bindings comments |
| … Σw/Σw_off | 0.520954 / 0.581605 | 0.520188 / 0.579241 | **0.799467 / 0.821436** | same |
| … `survival()` = `meta["fsi_survival"]` | 0.520239 / 0.582899 | 0.520239 / 0.582899 (unmoved: model-level) | **0.797094 / 0.819119** | same |
| … per-event w, min … max | 0.0208 … 1.4043 / 0.2139 … 6.9523 | 0.0181 … 1.4138 / 0.2124 … 6.9625 | **0.0681 … 37.92 / 0.2026 … 50.0** (2 events at `w_max`) | same |
| … w_n/w_c percentiles 1 / 5 / 25 / 50 / 75 / 95 / 99, max | 0.821 / 0.827 / 0.853 / 0.897 / 0.934 / 3.337 / 5.935, 68.52 | 0.821 / 0.826 / 0.853 / 0.897 / 0.934 / 3.322 / 5.935, **69.89** | **0.777 / 0.805 / 0.980 / 0.991 / 1.006 / 1.220 / 6.600, 10.14** | fsi.hpp; OPEN_ITEMS §7.1 |
| … events differing by > 1 % | 99.50 % | 99.425 % | **64.515 %** | same, README, test comments |
| … weight-table cells clipped at `w_max` = 50, 40 / 20 mb | — | 0 / 0 | **0.70 % / 5.25 %; 0 / 0.34 %** | PHYSICS_CHANNELS row 347; USAGE; README |
| ⁶Li α survivals at 20 mb | 0.671 (cluster) | 0.6712 / 0.6555 | **0.8636 / 0.8576** | USAGE §3; PHYSICS_CHANNELS 321-322 |
| ⁷Li α FSI (`helicity-flip`, seed 1234), 40 mb: survival, clipped cells, w max | 0.365 (cluster, "measured, not test-pinned") | 0.3648 / 0.4677; 0 / 0; 0.888 / 1.280 | **0.6055 / 0.6290; 4.41 % / 9.52 %; 49.84 / 50.0 (4 nucleon events at `w_max`; 12 with w > 10 carry 2.46 % of Σw)** | PHYSICS_CHANNELS row 349 |
| … at 20 mb | — | 0.5541 / 0.5326; 0 / 0 | **0.7377 / 0.7244; 2.90 % / 4.98 %; w max 24.06 / 47.39** | same |
| d–p control FSI, 40 / 20 mb: survival, clipped cells, w max | — | 0.8339; **1.451 %** / 0; **50.0** / 13.0 | 0.8785; 1.346 % / 0.031 %; 18.06 / 6.01 | PHYSICS_CHANNELS row 347 |
| RC band clamp, the PUBLISHED recipe (OPEN_ITEMS lines 832-833, 911-913: 20 000 events, `--rc tensor-band`, ⁶Li seed 1 at `tensor_thirds_plan(0.0, 0.6)`, ⁷Li seed 11 at the J = 3/2 fill of `test_rc.py`) | 245 (1.2250 %) / 520 (2.600 %) | **245 / 520 — reproduced exactly**; ⁶Li seeds 1-5: 245 / 223 / 269 / 219 / 235; ⁷Li seed 1: 502 | **16 (0.080 %) / 507 (2.535 %)**; ⁶Li seeds 1-5: 16 / 20 / 21 / 30 / 34; ⁷Li seed 1: 499 | USAGE §7b (2229-2231, 2242, 2559); OPEN_ITEMS 832, 858, 911-913 |
| … on the CLI-default fill (the design's recipe, labelled) | — | 237 (6Li, P_z 0.7, seed 1); 698 (7Li `helicity-flip --pe 0`, seed 1) | 14; 700 | record only |
| ⁶Li α tag fraction, 400 k events, the C++ test's setup (YR / tagging, legacy levers / tagging) | 0.0248 / 0.3046 / 0.2547 (stale: pre-2026-09-06 sample) | **0.024885 / 0.305123 / 0.25432** (the C++ MESSAGE, reproduced in Python) | **0.0334375 / 0.287445 / 0.248635** | USAGE §5 |
| ⁷Li α tag fraction, configs 0 / 1 / 2, YR (tagging) | 0.966 / 0.974 / 0.979 (0.981 / 0.993 / 0.994) | 0.96623 / 0.97387 / 0.97925 (0.98065 / 0.99305 / 0.99441) | **0.99846 at all three (0.99955 / 0.99938 / 0.99951)** | USAGE §5 |
| spectator-marginal closure (tests/test_pipeline.cpp line 317; a default run: it follows the wave) | — | χ²/ndf(k) 20.43/24 = 0.85, (cos θ_k) 22.84/16 = 1.43 | **19.79/13 = 1.52 (floor ndf_k > 10, χ²/ndf < 2.0: margins 3 bins and 0.48), 17.81/16 = 1.11** | record |
| T1 mean struck-nucleon virtuality ⟨p_N² − M_N²⟩, 50 k events (tests/test_breakup.cpp's setup: config 1, seed 20260830, tagging optics on the lithium tags), ⁶Li α / ⁷Li α / d–p — **found by this phase**: neither the design nor the check lists it | −0.057 / −0.124 / −0.035 GeV² | −0.05675 / −0.12375 / −0.03494 | **−0.05817 / −0.08468 / −0.03609** | docs/T2_CHAIN.md 109-113 |
| … the draws T1's own timelike-X rejection redraws (the T0 and T1 `x` columns compared index by index; `i6/meas/m_t1_virt.py`) | 0.024 / 0.022 / 0.004 % | 0.038 / 0.022 / 0.002 % (19 / 11 / 1 of 50 k) | **0.054 / 0.014 / 0.006 % (27 / 7 / 3)** | docs/T2_CHAIN.md 112-113; USAGE line 1578 (was "0.02 %"); tests/test_breakup.cpp 416-417 (its 2000-index case: 2 redraws on Hulthén, 0 on the default, §5) |
| `price_defaults.py`, li6-tag / li7-tag families (N = 200 000, seed 20261026) | default rows `cluster_wave = hulthen` | — | default `cluster_wave = vmc`; the priced `hulthen` row: A_zz^tag ×2.32488 / ×4.78404, RP ×0.72074; ⁷Li ⟨k⟩ ×1.55053, RP ×0.97543 | `docs/benchmarking/figures/price_defaults.png` |

Two pre-existing stalenesses were found on the way and are now labelled, not
silently replaced: the D3 table's 68.52 / 99.50 % are printed TODAY by the
reduced-statistics C++ case (4000 events on its 40 × 28 test grid: max 68.5222,
99.5 %), while the 20 000-event CLI recipe the text names gives 69.89 / 99.425 %;
and USAGE §5's 0.3046 / 0.2547 predate the 2026-09-06 sign fix (0.305123 /
0.25432 since, as the C++ test prints).

The two T1 rows were labelled line-neutrally in their three files, outside
every fingerprinted range (T2_CHAIN.md's 41-45, 79-96 and 102-106; USAGE's
ranges) and on lines nothing cites; tests/test_breakup.cpp's twelve cases
re-run after the comment edit: 12 / 12, 6 104 706 assertions. T2_CHAIN.md
§1a's other T1 sample numbers were NOT re-measured — §9 item 10.

**The D/S criterion on the α–d channel itself** (the check's risk 1; the design
answered the author's "measured D/S ratio" with the deuteron's η_d only).
Computed with `t2_deuteron_static`'s own residue formula η = −N₂A/N₀
(`i6/meas/m_eta_ad.py`, which reproduces the deuteron's +0.032478): the
Hulthén α–d pair (κ = 0.06066 GeV, β = 0.30 GeV, P_D = 0.0867) has
**η(⁶Li → α + d) = +0.070098**, against George–Knutson's −0.025 ± 0.006 ± 0.010
(`LI6_ETA_DS_GK`) and the VMC tables' −0.0482 (USAGE.md line 2923, `OverlapRaw`
−0.0538 ± 0.0021). INFERRED: the Hulthén α–d is the wrong SIGN (+8.2 σ on GK's
errors in quadrature) where VMC is −2.0 σ; the S–D interference this sign enters
is what drives A_zz^tag ×0.43 / ×0.21 between the waves. So the measured D/S
supports D6 on the α–d channel too — and it makes the quoted cluster-wave
systematic, hulthen − vmc, a difference of two models one of which has the
wrong asymptotic D/S sign, not a band around the default. USAGE.md line 1536
now says so.

**Runtime** (the check's wrong-item 8; its own `d6chk/warm2.py`, fresh
processes, output `i6/meas/warm_7li.txt`): the first `generate()` of a ⁷Li α VMC pipeline carries a one-time cost —
50 000 events 0.163 s against 0.092 s on Hulthén, 800 000 events 1.941 s against
1.444 s; later calls at par on both waves (0.082 – 0.089 s; 1.30 – 1.39 s). ⁶Li α shows
no difference (0.357 – 0.395 s for 200 000 on both). The check's own run gave 0.217 / 0.120 s and
2.587 / 1.905 s. The table load itself is the design's
≲ 7 ms and ≈ 2.4 MB.

## 3. The flip-only failures (Gate 2a): each a Hulthén pin

The library edits of the flip (`s2a` – `s2g`: the constant, the field, the
`validate()` rule and the knob row, the bindings, `make_config`, the CLI, the
example, `price_defaults.py`, the header documentation), built and run before
any re-take, on top of the step-1 tree:

- **C++, full suite**: 419 cases, **418 passed, 1 failed — 4 of 17 244 258
  assertions** (452.18 s): "pipeline: cluster_wave is refused on the channels
  that never read it", whose comment read "Hulthen -- the default -- is legal
  everywhere": `CHECK_NOTHROW(validate())` with `Hulthen` now throws on
  inclusive and coherent (a value moved from the default), and the `VmcAV18`
  refusal no longer fires there (it is the default). A pin of the Hulthén
  DEFAULT → re-taken. The +18 assertions against 17 244 240 are counted, not
  failed (attributed in §5).
- **pytest, full suite**: **14 failed**, 1455 passed, 149 skipped (546.07 s).
  Six are the design's predicted list, each a pin of the Hulthén default:
  `test_cli::test_cluster_wave_flag_reaches_the_config` (the default resolved to
  `"hulthen"`), `test_module::test_make_config_takes_cluster_wave_by_name` (the
  default was `Hulthen`), `test_module::test_c5_5b_cluster_wave_is_refused_where_it_is_never_read`
  (Hulthén legal everywhere, VMC refused), and the knob cell
  `cluster_wave:hulthen -> vmc` on the three tagged specs (a variation from the
  Hulthén default; the variant IS the default now, so the hash does not move and
  the row says `read`). The other eight are `test_doc_link_gate.py`: the strict
  gate saw the 15 named citations the flip's two growing files moved
  (include/lipolgen/cluster.hpp +12 lines below line 408, src/core/cluster.cpp +10
  from line 79) and the one cited block it edited (include/lipolgen/tagged.hpp
  54-62) — repaired in §7. Nothing else failed: `test_price_defaults` passed
  through the `make_config` guard (the design's INFERRED point, now measured),
  `test_rc.py`'s tagged band still clips on the VMC default at its seed 4, and
  the step-1 explicit pins held.

## 4. The knob matrix, and the per-file pass/skip tables

**The re-take (amendment 3).** `V()` does not copy `base` into `var`, so
`base=HUL` alone would have run the "on hulthen" variants on the VMC default.
Every Hulthén cell now carries the wave in BOTH (`base=HUL, var=dict(HUL, …)`),
and the plain cells are named for what they run on:

| cell (python/tests/test_knob_provenance.py lines 332-341) | tagged specs | inclusive / coherent specs |
|---|---|---|
| `cluster_wave: vmc -> hulthen` (`var=HUL`) | read, the hash moves | the explicit `hulthen` is refused by `make_config`; the row says `refused` |
| `cluster_vmc_mc_sigma: 0 -> 1 on hulthen` (`base=HUL, var=HUL`) | refused (needs VmcAV18); the row says `refused` | base refused → skipped |
| `cluster_vmc_mc_sigma: 0 -> 1 on the default (vmc)` | read on ⁶Li / ⁷Li (the band runs), refused on d–p | refused (needs a lithium tag) |
| `cluster_beta: 0.30 -> 0.33 on hulthen` | **read on all three — the Hulthén READ branch stays covered** | base refused → skipped |
| `cluster_beta: 0.30 -> 0.33 on the default (vmc)` | not read on ⁶Li and d–p, **read on ⁷Li at T1** | not read (off channel) |
| `p_d: 0.0867 -> 0.10 on hulthen` | read on ⁶Li α, not read on ⁷Li and d–p | base refused → skipped |
| `p_d: 0.0867 -> 0.10 on the default (vmc)` | not read | not read (off channel) |

Line-neutral (the records cite the file at 204, 231, 298-301, 312-315,
332-342, 507 and 586-597; lines 332-341 hold these cells, as they held the old ones). **Proof by the per-test comparison**
(`i6/cell_map.py`, mapping the old ids to the new ones per spec class: on the
tagged specs the old plain cell — default Hulthén — is the new "on hulthen" one
and the old "on vmc" the new "on the default (vmc)"; on the untagged specs the
other way round): **517 passed / 145 skipped before and after; under that map 0 tests
that ran stopped running, 0 gained, 0 unmapped**. Counted literally (F6, §12): 42 ids
renamed, and 4 kept their id and went passed → skipped — the coherent and three
inclusive `cluster_vmc_mc_sigma: 0 -> 1 on hulthen` cells, whose coverage moved to the
renamed "… on the default (vmc)" cells of the same specs, all four passing.

**Full pytest, per file** (`i6/logs/base_full.xml` against `i6/logs/g2_full.xml`):

| file | before passed / skipped / failed | after | time before → after |
|---|---|---|---|
| `test_b1_model.py` | 51 / 0 / 0 | 51 / 0 / 0 | 42.9 → 43.3 s |
| `test_bench_boeglin.py` | 13 / 0 / 0 | 13 / 0 / 0 | 1.8 → 1.9 s |
| `test_bench_bonus.py` | 26 / 0 / 0 | 26 / 0 / 0 | 11.5 → 11.7 s |
| `test_bench_chain_rc.py` | 4 / 0 / 0 | 4 / 0 / 0 | 20.3 → 20.1 s |
| `test_bench_deeps.py` | 16 / 0 / 0 | 16 / 0 / 0 | 2.4 → 2.4 s |
| `test_bench_djangoh_eltail.py` | 18 / 1 / 0 | 18 / 1 / 0 | 41.7 → 41.2 s |
| `test_bench_li6_ff_theory.py` | 18 / 0 / 0 | 18 / 0 / 0 | 2.2 → 2.0 s |
| `test_bench_nuclear.py` | 10 / 0 / 0 | 10 / 0 / 0 | 1.8 → 1.7 s |
| `test_bench_report.py` | 25 / 1 / 0 | 25 / 1 / 0 | 3.0 → 3.1 s |
| `test_bench_rfy.py` | 16 / 0 / 0 | 16 / 0 / 0 | 1.0 → 0.9 s |
| `test_bench_spin_deuteron.py` | 24 / 0 / 0 | 24 / 0 / 0 | 12.5 → 12.4 s |
| `test_bench_static.py` | 18 / 0 / 0 | 18 / 0 / 0 | 5.1 → 5.1 s |
| `test_bench_t1_nucleon.py` | 49 / 1 / 0 | 49 / 1 / 0 | 3.7 → 3.8 s |
| `test_bench_t1_r_world.py` | 51 / 0 / 0 | 51 / 0 / 0 | 6.0 → 6.3 s |
| `test_bench_t4_t6.py` | 18 / 1 / 0 | 18 / 1 / 0 | 20.1 → 20.0 s |
| `test_bench_t5_closures.py` | 19 / 0 / 0 | 19 / 0 / 0 | 5.9 → 5.8 s |
| `test_bench_t5_delta_bound.py` | 23 / 0 / 0 | 23 / 0 / 0 | 9.6 → 9.5 s |
| `test_cli.py` | 33 / 0 / 0 | **35 / 0 / 0** (two new D6 tests) | 4.0 → 4.2 s |
| `test_cluster_config.py` | 27 / 0 / 0 | 27 / 0 / 0 | 3.1 → 3.1 s |
| `test_delta_convention.py` | 6 / 0 / 0 | 6 / 0 / 0 | 0.1 → 0.1 s |
| `test_doc_link_gate.py` | 112 / 0 / 0 | 112 / 0 / 0 | 24.4 → 23.8 s |
| `test_hepmc.py` | 6 / 0 / 0 | 6 / 0 / 0 | 0.1 → 0.1 s |
| `test_hfs.py` | 7 / 0 / 0 | 7 / 0 / 0 | 0.6 → 0.6 s |
| `test_kernel_reference.py` | 11 / 0 / 0 | 11 / 0 / 0 | 0.0 → 0.0 s |
| `test_knob_provenance.py` | 517 / 145 / 0 | 517 / 145 / 0 (42 ids renamed and 4 same-id swaps, mapped one to one) | 198.3 → 196.3 s |
| `test_li6_unpol_band.py` | 12 / 0 / 0 | 12 / 0 / 0 | 17.9 → 17.3 s |
| `test_li7_rank2.py` | 27 / 0 / 0 | 27 / 0 / 0 | 7.5 → 7.5 s |
| `test_meta_provenance.py` | 8 / 0 / 0 | 8 / 0 / 0 | 2.7 → 2.8 s |
| `test_module.py` | 89 / 0 / 0 | 89 / 0 / 0 | 3.4 → 3.4 s |
| `test_mstw_sf.py` | 25 / 0 / 0 | 25 / 0 / 0 | 0.3 → 0.3 s |
| `test_o5_reach.py` | 18 / 0 / 0 | 18 / 0 / 0 | 0.7 → 0.7 s |
| `test_pipeline.py` | 18 / 0 / 0 | **19 / 0 / 0** (the VMC twin) | 2.1 → 2.5 s |
| `test_pom_set_band.py` | 9 / 0 / 0 | 9 / 0 / 0 | 4.9 → 4.8 s |
| `test_price_defaults.py` | 1 / 0 / 0 | 1 / 0 / 0 | 3.0 → 3.1 s |
| `test_pzz_mode.py` | 25 / 0 / 0 | 25 / 0 / 0 | 1.2 → 1.2 s |
| `test_rc.py` | 32 / 0 / 0 | 32 / 0 / 0 | 58.7 → 59.3 s |
| `test_rc_low_x_anchor.py` | 5 / 0 / 0 | 5 / 0 / 0 | 2.4 → 2.4 s |
| `test_release_metadata.py` | 8 / 0 / 0 | 8 / 0 / 0 | 0.0 → 0.0 s |
| `test_sampler.py` | 9 / 0 / 0 | 9 / 0 / 0 | 0.1 → 0.1 s |
| `test_sf_backend.py` | 54 / 0 / 0 | 54 / 0 / 0 | 9.0 → 8.9 s |
| `test_spdx_headers.py` | 6 / 0 / 0 | 6 / 0 / 0 | 0.0 → 0.0 s |
| `test_throughput.py` | 5 / 0 / 0 | 5 / 0 / 0 | 2.0 → 1.9 s |

Per test: 91 id changes = the 3 new tests + the matrix's 88 (42 ids gone, 42 new,
and the 4 same-id swaps above); once mapped, no test that passed before fails,
skips or disappears.

**Bench subset, per file** (`LIPOLGEN_BENCH=1 … -k bench`; `base_bench.xml`
against `g2_bench.xml`): identical in every file but `test_bench_report.py`,
26 / 0 / 0 → 25 / 0 / 1 — the committed-report test, whose diff is exactly the
four label strings of §6 — and 26 passed again once the report was regenerated:

| file | before | after | time |
|---|---|---|---|
| `test_bench_boeglin.py` | 13 / 0 / 0 | 13 / 0 / 0 | 1.9 → 1.8 s |
| `test_bench_bonus.py` | 26 / 0 / 0 | 26 / 0 / 0 | 11.6 → 11.7 s |
| `test_bench_chain_rc.py` | 4 / 0 / 0 | 4 / 0 / 0 | 20.6 → 20.5 s |
| `test_bench_deeps.py` | 16 / 0 / 0 | 16 / 0 / 0 | 2.4 → 2.4 s |
| `test_bench_djangoh_eltail.py` | 19 / 0 / 0 | 19 / 0 / 0 | 133.4 → 131.6 s |
| `test_bench_li6_ff_theory.py` | 18 / 0 / 0 | 18 / 0 / 0 | 1.9 → 1.9 s |
| `test_bench_nuclear.py` | 10 / 0 / 0 | 10 / 0 / 0 | 1.8 → 1.7 s |
| `test_bench_report.py` | 26 / 0 / 0 | 25 / 0 / 1 → **26 / 0 / 0 after the regeneration** | 79.2 → 78.2 s |
| `test_bench_rfy.py` | 16 / 0 / 0 | 16 / 0 / 0 | 1.0 → 0.9 s |
| `test_bench_spin_deuteron.py` | 24 / 0 / 0 | 24 / 0 / 0 | 12.3 → 12.3 s |
| `test_bench_static.py` | 18 / 0 / 0 | 18 / 0 / 0 | 4.9 → 5.1 s |
| `test_bench_t1_nucleon.py` | 49 / 1 / 0 | 49 / 1 / 0 | 3.8 → 3.7 s |
| `test_bench_t1_r_world.py` | 51 / 0 / 0 | 51 / 0 / 0 | 6.2 → 6.0 s |
| `test_bench_t4_t6.py` | 19 / 0 / 0 | 19 / 0 / 0 | 34.0 → 34.2 s |
| `test_bench_t5_closures.py` | 19 / 0 / 0 | 19 / 0 / 0 | 5.8 → 5.9 s |
| `test_bench_t5_delta_bound.py` | 23 / 0 / 0 | 23 / 0 / 0 | 9.3 → 9.5 s |

## 5. Gates, in the design's order, with tallies and wall times

**Baseline on `cb75ec7`, before any edit:** C++ 419 / 419, 17 244 240
assertions, 449.31 s; full pytest 1469 passed / 149 skipped / 0 failed, 539.06 s
(wall 539.95 s); bench 351 / 1 / 0, 330.65 s (the skip is the nnpdf-data wheel,
an environment fact); strict 1241 references / 95 ranges / 7 external / 0
broken / 24 allow-listed (SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8
/ 0 / 0); `--records-only` 1432 citations in 76 records: 933 checked live, 61
into another record, 110 annotated historical, 129 read by hand, 199 skipped,
0 broken; `--audit-ranges` 209 range citations on 179 blocks (1 by name, 66 by
quoted text, 39 by quoted symbol, 103 by exemption), 0 REFUSED; SPDX 133 files.

1. **Gate 1 — step 1, every Hulthén reliance made explicit, the dump preflight
   added; no behaviour change:** build 12.3 s; C++ 419 / 419, **17 244 240**
   assertions (unchanged), 448.81 s; full pytest **1469 / 149 / 0, 0 per-test
   outcome changes**, 539.12 s; the citation gates and SPDX identical to the
   baseline. The dump preflight was exercised alone (`python3 -B`, no file
   written, no `.pyc`; no file in the sibling is newer than the design record): it passes on the
   sibling's `4726812` and REFUSES a monkeypatched VMC `li6_alpha_channel`
   default and a tabulated deuteron default.
2. **cmake build** of the flip: 2 min 14 s, no warning; of the re-takes 7.7 s;
   of the post-gate comment edits 7.9 s; no warning anywhere.
3. **`./build/lipolgen_tests`**: the flip alone, §3; final **422 cases, 422
   passed** (419 + the three D6 cases), **17 244 309 assertions**, 448.12 s.
   The +69 against the baseline: +18 from the flip itself (Gate 2a, every one
   passing; the per-case attribution is in the paragraph below) and +51 from the
   new D6 cases.
4. **The targeted pytest files** (the twelve the change touches, with
   `LIPOLGEN_BENCH=1`): **886 passed, 145 skipped, 0 failed**, 319.17 s.
5. **The full pytest**: **1472 passed, 149 skipped, 0 failed, 537.06 s** (wall
   538.06 s) — §4's table.
6. **The bench subset**: 350 passed, 1 skipped, 1 failed (the committed-report
   test, §4), 328.45 s; after the regeneration that test and
   `test_doc_link_gate.py` together: 138 passed, 101.8 s.
7. **The three citation-gate modes**, final: strict **1241 / 95 / 7 / 0 / 24**
   (SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0 / 0);
   `--records-only` **1432 citations in 76 records, 933 / 61 / 110 / 129 / 199,
   0 broken** before this record, which adds none (it carries no gated
   citation); `--audit-ranges` **209 on 179 blocks, 0 REFUSED** (1 / 66 / 39 /
   103, the baseline's).
8. **SPDX**: 133 files.
9. **The report** (§6) and **`make_report.py --check`**: OK (75.3 s).
10. **`price_defaults.png`** regenerated (§6).
11. **After this record and the last edits** (the example's banner, the T1
    labels in T2_CHAIN.md, USAGE and tests/test_breakup.cpp, the board row):
    build 3.2 s, no warning; tests/test_breakup.cpp's cases 12 / 12 (6 104 706
    assertions, 1.7 s); the example re-run on all three channels (§8, §9 item
    9); strict 1241 / 95 / 7 / 0 / 24 (2.4 s); `--records-only` **1432
    citations in 77 records** — this record adds none — 933 / 61 / 110 / 129 /
    199, 0 broken; `--audit-ranges` 209 on 179 blocks, 0 REFUSED; SPDX 133;
    `test_doc_link_gate.py`, `test_meta_provenance.py` and
    `test_spdx_headers.py`: 126 passed, 26.7 s. Then the relocations of §7
    (the two new CLI tests to the end of test_cli.py; the early re-take
    comments of test_module.py, test_rc.py and test_bench_boeglin.py folded
    to their original line counts; no assertion changed but the CLI
    re-take's two extra asserts, dropped — the banner test covers an explicit
    `hulthen`), each file re-run: test_cli.py 35 passed (4.3 s);
    test_module.py with test_bench_boeglin.py (`LIPOLGEN_BENCH=1`) 102 passed
    (5.5 s); the edited test_rc.py case passed; the gates again identical
    (strict 1241 / 95 / 7 / 0 / 24; R6 1432 in 77 records, 0 broken). Test ids
    are unchanged, so §4's per-test comparison stands.
12. **The final tree, both full suites again, one at a time**: pytest
    **1472 passed, 149 skipped, 0 failed, 536.92 s** (wall 537.97 s; 0
    per-test outcome changes against item 5's run, per file as §4's table);
    `./build/lipolgen_tests` **422 / 422, 17 244 309 assertions**, wall
    449.17 s (the tier-switch case prints "T1 timelike-X redraws: 0 of 2000").
    Then `PHYSICS_CHANNELS.md` row 347's 2.91 → 2.90 % (§8) and this record's
    last words: strict 1241 / 95 / 7 / 0 / 24; `--records-only` 1432 citations
    in 77 records, 933 / 61 / 110 / 129 / 199, 0 broken; `--audit-ranges` 209
    on 179 blocks, 0 REFUSED; SPDX 133; `test_doc_link_gate.py`,
    `test_meta_provenance.py` and `test_spdx_headers.py` 126 passed, 26.4 s.

**The +18 assertions of the flip, attributed.** The final test tree was run
twice with doctest's XML reporter — as built, and as a scratch copy of the
tracked files (`i6/attr_copy`, never the working tree) with
`CLUSTER_WAVE_DEFAULT` set back to `Hulthen` (`i6/attr_run.sh`;
`i6/logs/attr_vmc.xml`, `attr_hul.xml`; compared case by case by
`i6/attr_cmp.py`): 17 244 309 against 17 244 291 assertions, **+18, Gate 2a's
number**. Three cases carry all of it, each passing on both waves and each
counting assertions per event of a DEFAULT tagged run, so the count follows
the event stream: "T1 pipeline: T0 is unchanged by the tier switch, and T1
adds only the partners and the struck nucleon" (tests/test_breakup.cpp line
399) 19 990 → 20 006, +16 — eight checks are skipped per index that T1's
timelike-X rejection redraws, 2 of its 2000 indices on Hulthén and 0 on VMC;
"T1 pipeline: whole-record 4-momentum and charge conservation on every tagged
channel" (line 301) 766 817 → 766 818, +1; and "triton sf: T1 pipeline with
--triton-sf ciofi-simula …" (tests/test_triton_sf.cpp line 175) 65 155 →
65 156, +1 (its per-event branches differ in length by breakup channel). On the
Hulthén copy the four D6-aware cases fail, as they must — the re-taken refusal
case (2 of its 12 assertions) and the three new D6 cases (2 of 12, 14 of 34,
2 of 5): the new pins catch a default moved back.

## 6. The report, the figures and `price_defaults.png`

`LIPOLGEN_BENCH=1 python3 validation/benchmarks/make_report.py --all --plots
--compare <session scratch i6/pre_REPORT.json, sha256 6ae3d5f6…ff5e64> --out-dir
docs/benchmarking` (86.5 s): **13 pass, 10 fail, 0 blocked — as before.**
`--compare` printed **2 changes, both `MOVED` generator strings with every
number unchanged** (`t2_bonus_spectator_shape`: "Hulthen 40 mb GlauberFsiWeight …"
→ "Hulthen (rule (d)'s configuration; AV18 ships since D6 and fails too) 40 mb
GlauberFsiWeight …"; `t2_deeps_spectator_tail`: "…, Hulthen default, LC …" →
"…, Hulthen (the declared verdict wave, deuteron_channel(); an option since
D6), LC …"). A structured diff of the two `REPORT.json` (`i6/logs/report_json_diff.txt`)
finds exactly the label text the task expected on rows **10, 11, 16 and
`t2_deuteron_static`**: the BONuS and Deeps generators, Deeps' tolerance string
("default wave" → "the declared wave (deuteron_channel(): Hulthen)"), Boeglin's
generator ("hulthen (default)" → "hulthen (the design block's verdict wave)" —
`--compare` does not flag it: its headline numbers are unchanged), and the
`hulthen` sub-row label of `t2_deuteron_static`; plus 21 wall-time `runtime_s`
values and the generation stamp. **No status, verdict or number moved. T6**
(`t6_hepmc3_convention`): 23 / 23 clauses, PASS, unchanged — runtime 16.289 →
16.446 s, `docs/HEPMC3_CONVENTION.md` untouched (its sha256 is T6's pin).

**Figures**: `t2_deeps_spectator_tail.png` (the legend: "Hulthen (the verdict
wave; an option since D6)") and `t2_deuteron_static.png` (the row label
"Hulthen (tagged option)") changed in content; the other 17 regenerated figures
differed from the committed ones only in their 8-row footer stamp
(`i6/fig_diff.py`, pixel by pixel) and were restored from `HEAD`, as I4 and I5
did. **`price_defaults.png`**: `validation/price_defaults.py --plot` run once
(N = 200 000, six families, 166.6 s): its default rows are now named from
`CLUSTER_WAVE_DEFAULT` ("shipped default (cluster_wave = vmc)",
"shipped default (vmc, hulthen, toy, toy)"), the priced alternative is
`cluster_wave = hulthen`, and the figure's series follows — the ratios are the
inverses of the 2026-10-26 ones (§2) — but the plot's fixed 0.1 … 2.6 axis left
the ×4.78 A_zz^tag(0.30) point undrawn, which this phase did not notice (F6
widened the axis to the drawn points, §12). Its refusal probes of `cluster_wave = vmc`
on the inclusive and coherent families now quote the `make_config` refusal.

**The four rows, read against their declared rules** (no verdict re-pointed
after its numbers were seen — the ground rule):

- **Row 11, Deeps**: the declaration names "the shipped default,
  `deuteron_channel()` = the analytic Hulthen pair … (the tagged-channel
  default)". Under DD-1 `deuteron_channel()` stays Hulthén, and every χ² was
  known when D6 was taken, so the verdict wave stays Hulthén (the design's §7.1
  re-pointing is NOT applied — the check's wrong-item 10); AV18, now the
  tagged runs' shipped wave, stays the recorded sub-row 579.12/187; FAIL either
  way. A dated D6 note in the harness header says this.
- **Row 10, BONuS**: rule (d) defines the verdict configuration as the Hulthén
  deuteron at 40 mb; `DEFAULT_CONFIG` kept; AV18 at 40 mb, recorded under (e),
  fails too (3/12, 0/16). A dated note under the MEASURED block.
- **Row 16, Boeglin**: the hashed block's "(hulthen)" names the verdict wave;
  the block and `DEFAULT_WAVE` are untouched (`design_block_sha256()`
  unchanged — the bench pins it); a D6 addendum OUTSIDE the block: AV18 passes
  by the same rule (11.54/21, p = 0.951). PASS either way.
- **`t2_deuteron_static`**: the rule covers every shipped wave function and
  the Hulthén still ships, as the option: FAIL stays (+17.2σ); labels and a
  dated note only.

## 7. Citations: `--fix`, the re-read EDITED blocks, the records

- **`--fix`** re-pointed the 15 named citations of `docs/PHYSICS_CHANNELS.md`
  that the two growing files moved (include/lipolgen/cluster.hpp
  `cluster_wave_name`, `Wave`, `theta_lm`; src/core/cluster.cpp `data_dir`,
  `data_path`, `read_anl_overlap`, `read_anl_momentum`, `read_anl_momentum_norms`,
  `read_fdeut_k`, `vmc_sign_steps`, `vmc_from_overlap_k`, `vmc_from_overlap_r`,
  `vmc_from_momentum`, `radial`, `theta_lm`) — each onto its declaration, +10 or
  +12 lines. No ALLOW-exempted citation moved (`--fix` never touches those); the
  ALLOW / UNNAMED_ALLOW / RANGE_ALLOW targets the check named — src/core/pipeline.cpp
  lines 847, 1116, 1877; include/lipolgen/constants.hpp line 139;
  tests/test_pipeline.cpp line 608 and 608-634; README.md line 128;
  docs/CONVENTIONS.md line 415, 207-211, 325-327; tests/test_tagged.cpp
  396-412/413; include/lipolgen/pipeline.hpp 931-934; docs/USAGE.md 2795-2992 —
  all kept line-neutral (every edit above them in the same file is line-neutral;
  `git diff --numstat` shows equal insertions and deletions in each such file).
- **Two EDITED blocks, re-read, then `--record-ranges`** (0 new, 2 changed, 0
  dropped, R5 0 refused): include/lipolgen/tagged.hpp 54-62, cited by
  `PHYSICS_CHANNELS.md` row 270 as the "same wave function seen in two
  experiments" rule — the block now says that rule holds on the Hulthén wave
  (the run default until D6), which is what the row says of "the Hulthén ⁶Li
  model", and the quoted phrase is kept on ONE line so R5 still finds it;
  docs/USAGE.md 1978-2169, cited by rows 390 and 393 for §7b's RC rules — the
  edit inside it is §5's tag-fraction sentence (relabelled, re-measured), which
  neither row cites, so both still say something true.
- **Records**: one record citation broke on the way — `plan_D7_backend_defaults.md`
  line 813 → README.md 39-41 ("pinned against at rtol 1e-12") — because the first
  draft of the README paragraph dropped the phrase. The paragraph was rewritten to
  keep it (and to say that the tagged JSON pins the factories' Hulthén family, not
  the default run), and the gate is back at 0 broken. **No citation R6 parses
  moved** (`i6/r6_moved.py` replays R6's own path resolution and maps every
  parsed record citation through `git diff -U0 cb75ec7`): the parsed ones into
  the changed files land on unshifted lines; ten land in hunks this phase
  relabelled — six checked live (USAGE.md 2761-2768; price_defaults.py 769 and
  817; bindings.cpp 640-660; README.md 39-41, twice), all still passing, and
  four annotated historical or read by hand — so no `(as of cb75ec7)` annotation
  was needed and no census in `test_doc_link_gate.py` was re-pinned (110
  annotated, 199 skipped ≤ 200, 62 read-pins). Bare-filename citations, which R6
  does not parse, were mapped the same way (`i6/moved_cites.py`), and those my
  shifts had moved while they were still right at `cb75ec7` were put back by
  making the edits above them line-neutral: the two new CLI tests went to the
  END of python/tests/test_cli.py (plan_D7's test_cli.py 163-178, which I7 will
  follow, lands on its test again), and the re-take comments of test_module.py
  (phase_D_small_items.md's 456 / 468), test_rc.py (plan_D5's
  999-1004) and test_bench_boeglin.py (plan_D6's 245-254) were folded back to
  their original line counts. Four remain moved, and are named here: plan_D6's
  own cluster.hpp 409-415 (the `cluster_wave_name` pattern, now 421-427 — the
  constant sits above it, beside the enum, where the design put it); plan_D6's
  test_module.py 594-600 (its end now 602 — the block it names is the one
  re-taken); docs/code_review_2026-08-29.md line 86's
  examples/generate_tagged.cpp 243 (the optics menu's lumi_fraction column, now
  line 275, below the example's D6 code; an ungated review); and plan_D6's
  test_pipeline.py 77-82 at its lines 554, 653 and 762 (the Hulthén tag-fraction
  test, now 89-96, below the new `tagged_run_hulthen` fixture — missed by this
  phase's map, named by the verification, §12). The four 2026-09-02
  design-record citations of src/core/cluster.cpp the scan also lists were
  already off their targets at `cb75ec7`. One ungated reference my shift had made stale was
  corrected: `docs/benchmarking/02_data_nucleon_deuteron.md` line 243's
  src/core/cluster.cpp 305–441 → 315–451 (CD-Bonn). The directory-less
  generate_tagged.cpp:123 (cited without its directory) in `docs/theory/SPIN32_FINITE_GAMMA.md` line 156 —
  which the gate does not parse — still lands on the line it meant: every cited
  line of the example keeps its number (SPIN32 cites 119-124 and 123,
  `PHYSICS_CHANNELS.md` 79, plan_D7 139), because the header's and `parse()`'s
  edits balance before each (only the header's old lines 15-18, now 14-17, and
  `parse()`'s old 87-97, now 88-98, sit one line off, and nothing cites them). `docs/theory/SPIN32_FINITE_GAMMA.md` gained no line (it was
  not edited).

## 8. The check's items: applied, and not

**Applied — the eleven amendments, every missing site, every "wrong" statement,
every underrated risk, with these resolutions:**

*The amendments.* (1) The example's `--beta` / `--pd` refusals read the
library's knob table for THIS run (examples/generate_tagged.cpp, after the
`Pipeline` is built: refused iff moved from the default and the row is not
`read`) — so `--beta` is accepted on ⁷Li α at T1 even on VMC (the triton
breakup reads it) and `--pd` is refused on ⁷Li and d–p on either wave; measured
on the I6 build: `--channel 7Li-alpha --plan apar --beta 0.33` runs and prints
"wave vmc (the run default, D6); beta = 0.33 (read)"; `--pd 0.1` on 7Li-alpha (on
either wave) and on d-p `--cluster-wave hulthen` are one-line usage errors quoting
the row's reason; `--channel 6Li-alpha --beta 0.33` is refused on the default and
runs on `--cluster-wave hulthen` ("hulthen (an option; the run default is vmc,
D6); beta = 0.33 (read)"); an unknown `--cluster-wave` is a usage error. The banner prints the wave, whether β was
read and the run's true tier — "tagged T1 generator" (it said T0 while running
at the library's default T1, which is why `--beta` IS read on ⁷Li). (2)
`cluster_vmc_mc_sigma`: `validate()` already handles it (a non-zero band is
accepted on a default lithium α tag, refused on Hulthén, d–p and untagged);
pinned in C++ and Python; its knob row is `read` there and its reason says so;
pipeline.hpp, the bindings and USAGE say so. (3) The matrix re-take, §4. (4) The
RC clamp on the published recipe, §2 (245 and 520 reproduce exactly; the
design's 237 / 698 are a different fill, kept only as a labelled record row). (5)
The FSI numbers, every site the check listed (fsi.hpp 129-150 rewritten with
both waves; CONVENTIONS 455-467; bindings 697-714; `PHYSICS_CHANNELS.md` 321-322,
344 (line-neutral; its unnamed citation of tests/test_fsi.cpp line 384 still on
the TEST_CASE line, whose name it quotes), 347 and 349; USAGE 1665-1668,
1679-1697 and 2559; OPEN_ITEMS 61, 85, 414-444, 832-833, 858, 911-913; README's
row), and the other missed sites: `--cluster-beta` / `--p-d` help;
`PHYSICS_CHANNELS.md` 267 (five files, the install caveat, the escape);
DEVELOPMENT_PLAN 270-272; `AUTHOR_DECISIONS.md` §B9 and §B15 (and §B6, §B10);
tests/test_pipeline.cpp 151-152; the Boeglin plot docstring; the dump
preflight; the fingerprinted USAGE 1978-2169 re-read and re-recorded. (6) The
statements: the error message and the docs say an installed library finds the
compiled-in source-tree `data/` while that tree exists and name
`$LIPOLGEN_DATA_DIR` or `--cluster-wave hulthen`; only `test_tagged.cpp`'s three
port-gate models use the factory defaults (now by name) — `test_cluster.cpp`
builds `Wave` from the JSON; the ⁷Li first-generate cost is recorded (§2);
"`--fix` never touches ALLOW-exempted named citations" is how §7 was done; the
pre-I6 `vmc` run equals the post-I6 default run with the `at_default` flip (§1).
(7) `meta["cluster_wave"]` and the banner line; HepMC3's missing trace stated
(USAGE line 2761, the CLI comment, §9). (8) The α–d D/S statement (§2) and the
row-11 misreading closed in the knob row itself (§0; pinned in C++ and
Python). (9) Only the tagged call sites receive Hulthén (no inclusive /
coherent site of `pipe()` or of a Python fixture does); the closure margins are
recorded (§2: ndf_k = 13 against > 10, χ²/ndf 1.52 against < 2.0); every
sample-level number is labelled "toy SFs". (10) The gates in order, §5. (11) This
record, the board row, the D2 pointer, the registry.

*Wrong statements and risks, not already covered.* DD-5.3: no CMake change
(installing `data/` alone would not repoint `data_dir()`, the check's point);
documented instead. Wrong-item 9: the registry and §B9 say truthfully that D6
is implemented as the run-surface default (DD-1), not as §B9(e)'s option (ii),
and that the factories and `tagged.json` are unchanged. Risk 2 (HepMC3) and
risk 3 (I7 re-moves the event stream): §9.

*Found by this phase, beyond the design's inventory and the check's lists.*
The T1 virtuality and timelike-X redraw numbers of docs/T2_CHAIN.md §1a,
USAGE line 1578 and tests/test_breakup.cpp's comment — default-run numbers
measured on the Hulthén — re-measured on both waves and labelled (§2); the
example's banner ("vmc (default; the run default is vmc, D6)" became "vmc
(the run default, D6)"); the ungated src/core/cluster.cpp range of
02_data_nucleon_deuteron.md (§7); `PHYSICS_CHANNELS.md` row 302's two
pre-existing inaccuracies (§9 item 8); the example's d–p run, which stops
at its optics menu on `cb75ec7` already (§9 item 9); and one rounding of my
own first pass — `PHYSICS_CHANNELS.md` row 347 carried the check's "2.91 %" for
the ⁷Li cluster cells clipped at 20 mb, where this phase's measurement is
0.0290498, i.e. 2.90 % (corrected, line-neutral).

**Not applied, with the reason:**

- the design's §7.1 Deeps re-pointing (`DEFAULT_WAVE = "av18"`, `test_bench_deeps.py`
  re-takes): superseded by the check's wrong-item 10 and the task's ground rule;
  labels only (§6);
- the optional A_zz^tag closure twin on VMC (tests/test_pipeline.cpp 473-552):
  its ndf floor would be chosen after the number is known (7 at YR, the design's
  §4.6) — left to the author, §9;
- making T25 – T27 and the wave-agnostic factory calls of the C++ tests
  explicit: they keep the factory-default calls on purpose (they pin DD-1);
  the polligen port gates and the P_D = 0 limits ARE explicit;
- the explicit `source=Hulthen` in the harness `Tree`s, `t5_sum_rules.py` and
  `repin_tagged_from_lipolgen.py` (the design's "recommended"): the factory
  default stays the Hulthén under DD-1 and `t2_deuteron_static.hulthen()`
  already refuses a non-analytic default; no code moved there;
- a CMake install of `data/` (DD-5.3), for the reason above;
- the D6 entry of `DECISIONS.md`: the author's decision record, which this task
  did not name — the corrected text is in §10 for the supervisor;
- `tests/test_t2.cpp` line 341 ("covers the Hulthen default"): read and left —
  it is the `--triton-sf` default, not the cluster wave;
- `beams.hpp`'s "(AV18-like)" for `P_D_DEUTERON` = 0.045: inside a fingerprinted
  range and a pre-existing statement; §9 item 1's sub-question;
- the knob rows' `cluster_beta` / `p_d` reasons: unchanged — they were already
  conditional on the wave and are exactly the rule the example now reads.

## 9. Open points

1. **For the author — the ⁶Li cluster constants (the design's §11.1).** With
   VMC the tagged default, the two DEFAULT ⁶Li runs describe ⁶Li with different
   wave functions: a default tagged α run reads 0.887076 whole-nucleus and
   0.982576 rank-2, a default inclusive run `LI6_CLUSTER_POLARIZATION` = 0.811228
   and `LI6_B1_RANK2_TRANSFER` = 0.921947 (×1.0935 / ×1.0658). Both constants
   are KEPT, as instructed. Which reading should the inclusive constants carry —
   (i) keep (the drift now between the two default runs, the band 0.81 … 0.91
   mandatory); (iii) re-derive to 0.887076 / 0.982576 — one wave function for
   both default runs, at the price of every inclusive ⁶Li g₁ / A_∥ ×1.0935, b₁ and
   the shipped inclusive A_zz ×1.0658, `beams.json` / `xsec.json` /
   `b1_default_li6.json` re-pinned, and the six-body 0.848 crossed from −4.34 %
   to +4.61 %; or (iv) 0.848 for the vector part (no ab-initio anchor for the
   rank-2 transfer)? Sub-question: the inclusive deuteron species'
   `DEUTERON_VECTOR_POLARIZATION` 0.9325 (`P_D_DEUTERON` 0.045, called
   "AV18-like" in beams.hpp, where AV18's is 0.057600). The inconsistency is
   stated where it lives: beams.hpp, tagged.hpp, constants.hpp, pipeline.hpp, the
   knob row's reason, CONVENTIONS, USAGE, `PHYSICS_CHANNELS.md` rows 270 and the
   §C notes, README, registry row 10, §B6.
2. **A polligen VMC port-gate block (a candidate).** `tagged.json` stays
   polligen's Hulthén (now NAMED in the gates, and the dump refuses a sibling-side
   flip of its defaults); the new default has no external port gate. The design
   measured the sibling's VMC α–d at 5.4 × 10⁻⁹ / 3.2 × 10⁻⁹ of LiPolGen's vector
   / tensor dilutions, so a `li6_alpha_vmc` block at a declared ~1e−8 is feasible
   — it would cover the α–d motion, not the embedded deuteron (the sibling keeps
   the scenario 0.9325) and not ⁷Li (the sibling refuses `vmc` there).
3. **The ⁷Li FSI weights on the default.** On VMC the ⁷Li weight table clips
   4.41 % / 9.52 % of its cells at 40 mb (2.90 % / 4.98 % at 20 mb), 4 of 20 000
   nucleon-variant events sit at `w_max` = 50 and 12 with w > 10 carry 2.46 % of
   Σw; the ⁶Li default clips 0.70 % / 5.25 % with 2 events at the ceiling; the
   d–p control clipped 1.45 % on Hulthén already (and reached 50). The clipped
   fraction is in `meta["fsi_clipped_grid_fraction"]`, not on the banner. Whether
   to band the ceiling, widen the grid or report it on the banner is open.
4. **HepMC3 provenance** — a HepMC3 file carries no record of the cluster wave
   (nor of the D4 convention or the D5 tail model): the same command writes
   files that differ by ×0.43 on A_zz^tag before and after D6 with nothing in
   `GenRunInfo` to say why. A run-level attribute changes
   `docs/HEPMC3_CONVENTION.md`, which T6 pins by sha256 — the same follow-up as
   I4 §8 item 10 and I5 §9 item 2.
5. **Order with I7.** Every sample-level number above is on the toy backends;
   D7's grids move the event stream (the check measured Σw 15989.349 →
   15989.374 at seed 1234 on `ct18nlo`), so the FSI table, the clamp counts and
   the fixed-seed tag fractions move again; the model-level ones (dilutions,
   survivals, A_zz^tag, η) do not.
6. **The optional VMC twin of the A_zz^tag closure** (§8): VMC gives ndf 7 at
   YR against the declared floor 8 (χ²/ndf 1.197, worst pull 1.75 — the design's
   flip copy); a twin would need a floor declared by the author.
7. **The spectator-marginal closure now runs on VMC** with smaller margins
   (ndf_k = 13 against > 10; χ²/ndf 1.52 against < 2.0): a seed, grid or N change
   on the default can trip it.
8. **Two pre-existing inaccuracies corrected in passing**: `PHYSICS_CHANNELS.md`
   row 302 said `--p-d` was the T1 deuteron's lever (it never reaches
   `BreakupOptions::p_d`; the `Pipeline` forwards `beta` and `source` only) and
   omitted `source` from what the `Pipeline` overwrites; `examples/generate_tagged`
   called itself a T0 generator while running T1.
9. **The example's d–p run stops at its optics menu** — pre-existing, not D6's.
   `./build/generate_tagged --channel d-p` generates, prints its category table
   and exits 1 with "tagging_optics: no tabulated working point for d at 10x100
   (the scan is only published for 6Li and 7Li)": the same on
   `--cluster-wave hulthen`, and on a build of `cb75ec7` (session scratch
   `i6/pre_copy`, `i6/logs/pre_copy_dp.txt`). The example prices the whole
   optics menu, tagging included, on every channel; the menu wants a d-control
   skip. Not fixed here — outside D6.
10. **T2_CHAIN.md §1a's other T1 sample numbers** — the 300-event residual
   table (inside the fingerprinted 79-96) and the "1.4 % (⁶Li α; 0.7 % for ⁷Li
   and the d control)" level of the collinear `hfs_sigma_empz_truth` (inside
   102-106) — were taken on the then-default Hulthén waves and are not
   re-measured. The residuals are numerical (10⁻¹³); the 1.4 % / 0.7 % follow
   the struck nucleon's Fermi p_T and so the wave (⁷Li's ⟨p_N² − M_N²⟩ moves
   −0.124 → −0.085 GeV² on the default). Re-measuring them edits two
   fingerprinted blocks: a human re-read and `--record-ranges`.
11. **A finite bound on `cluster_vmc_mc_sigma`** (F6, §12 V-8). `validate()` now
   refuses a non-finite band, but 1000 σ still builds (the α–d P_D becomes
   0.133). No bound is declared anywhere, and choosing one is a new rule, so it
   is left to the author.

## 10. Texts for the supervisor

**The D2 pointer** — APPLIED, one dated line under D2 in
`../run_2026-09-27/DECISIONS.md`: *"Superseded 2026-10-27 by D6 … the
`t2_deuteron_static` FAIL stays recorded — its rule covers every shipped wave
function, and the Hulthén still ships."*

**Registry rows 9, 10, 11** (`../run_2026-09-03/STATUS.md`) — APPLIED in place,
the design's §5.7 putting them in scope, as dated italic amendments at the end of
each row's text cell (no row grew a line, so nothing below moved): row 9
*decided and closed by D6 — implemented as the run-surface default (DD-1), NOT
as §B9(e)'s option (ii); the factories, `BreakupOptions::source`, `tagged.json`
and T25 unchanged; legs (1)–(4) re-read; the measured cost*; row 10 *the drift
is now between the two DEFAULT runs; no label and no constant moved; the
question re-posed (§9 item 1); the band mandatory*; row 11 *the refusal follows
the default; an explicit request refused by `make_config` / the CLI; the
off-channel knob row withholds the value*. Row 23 gained a pointer to row 10's
note. `AUTHOR_DECISIONS.md` got dated boxes at the ends of §B6, §B9 (saying
truthfully that this is not option (ii) and that (d)'s "`tagged.json` moves at
rtol 1e−12" does not happen) and §B10, and §B15's `cluster_vmc_mc_sigma` entry
says the band is read on the default lithium tags.

**`DECISIONS.md` D6, proposed** (not applied — the author's record; the
design's §12.1 with the check's corrections): keep the decision line, and in
the evidence and change paragraphs replace "Deeps' headline follows its
'default wave' (627.85 → 579.12 on 187 dof, FAIL either way)" by **"Deeps keeps
the Hulthén verdict its declaration names (`deuteron_channel()`, which D6 leaves
analytic; AV18 579.12/187 recorded; FAIL either way)"**; add after the η_d
sentence **"On the α–d channel itself the Hulthén pair's asymptotic D/S ratio is
+0.070, the opposite sign to George–Knutson's −0.025 ± 0.012 and the VMC −0.048
(`phase_I6_cluster_wave.md` §2)"**; add to "What changes" **"`validate()`
refuses a value moved from the default off the tagged channels, `make_config`
and the CLI an explicit `--cluster-wave` there, and the off-channel knob row
reads 'not read on channel …'; `meta["cluster_wave"]` and the run banner name
the wave (HepMC3 does not)"**; and keep "`LI6_CLUSTER_POLARIZATION` = 0.811228
and `LI6_B1_RANK2_TRANSFER` = 0.921947 are unchanged … a question put to the
author on registry row 10" verbatim.

## 11. The files, and the provenance

**Library** — include/lipolgen/cluster.hpp (the constant, the enum comment;
+12 lines below line 408); include/lipolgen/pipeline.hpp (the field and its
documentation, line-neutral); src/core/pipeline.cpp (the `validate()` rule and
message, the breakup comment, the `cluster_wave` and `cluster_vmc_mc_sigma` knob
rows; line-neutral); src/core/cluster.cpp (the missing-table message, +10);
include/lipolgen/tagged.hpp, beams.hpp, breakup.hpp, constants.hpp, fsi.hpp
(documentation, line-neutral). **Bindings and Python package** —
python/bindings.cpp (`meta["cluster_wave"]`, docstrings, the D3 comments,
`CLUSTER_WAVE_DEFAULT` and `cluster_wave_name` at the end); python/lipolgen/__init__.py
(the `CLUSTER_WAVES` comment, the `make_config` docstring and guard, the helper
at the end); python/lipolgen/cli.py (`DEFAULTS`, three help strings, one comment,
the banner line after line 1777). **Example** — examples/generate_tagged.cpp
(+32 lines net, all below line 144: the header's and `parse()`'s edits balance
by line 99, so the cited lines 79, 119-124 and 139 keep their numbers).
**Tests** — tests/test_pipeline.cpp, tests/test_fsi.cpp, tests/test_tagged.cpp,
tests/test_breakup.cpp (one comment, line-neutral);
python/tests/test_cli.py, test_module.py, test_knob_provenance.py,
test_pipeline.py, test_meta_provenance.py, test_rc.py, test_bench_boeglin.py,
test_bench_static.py. **Validation** — validation/price_defaults.py,
validation/dump_polligen_reference.py (the preflight, below every cited line),
validation/vmc_tag_fractions.py (one docstring word), the four T2 harnesses
(labels and dated notes; no rule, verdict wave or hashed block moved),
validation/physics_channels_ranges.json (two fingerprints). **Documents** —
README.md, python/README.md, docs/USAGE.md, docs/PHYSICS_CHANNELS.md, docs/T2_CHAIN.md,
docs/CONVENTIONS.md, docs/OPEN_ITEMS_SOLUTIONS.md, docs/DEVELOPMENT_PLAN.md,
docs/benchmarking/BENCHMARK_PLAN.md, 02_data_nucleon_deuteron.md,
08_channel_priority.md (labels; a revision-log row: C3 re-read, rank 6 keeps 3,
rank 7 keeps 2, **no rank moved**), validation/benchmarks/README.md;
docs/benchmarking/REPORT.md, REPORT.json, figures/t2_deeps_spectator_tail.png,
figures/t2_deuteron_static.png, figures/price_defaults.png (regenerated).
**Records** — `../run_2026-09-03/STATUS.md` (rows 9, 10, 11, 23),
`../run_2026-09-03/AUTHOR_DECISIONS.md`, `../run_2026-09-27/DECISIONS.md`, this
record, `STATUS.md`'s I6 row. **Not touched**: `validation/reference/*.json`
(empty diff), `docs/HEPMC3_CONVENTION.md`, `docs/theory/SPIN32_FINITE_GAMMA.md`,
every hashed design block, `LI6_CLUSTER_POLARIZATION`, `LI6_B1_RANK2_TRANSFER`
and every other constant, tolerance, registry label, census and
`ASOF_READ_BY_HAND` entry, the sibling repository (only grep'd, sed'd and
imported with `python3 -B`; no file in it is newer than the design record).

**Session scratch `i6/`**: `edits/` (`apply.py`; `s1a_cpp_tests.py`,
`s1b_py_tests.py`, `s1c_dump.py`; `s2a_core.py` … `s2j_knob.py`;
`s3_harness.py`; `s4a_readme.py` … `s4g_test_comments.py`; `s5_records.py`;
`s6_post.py`), `meas/` (`ab_runs.sh`, `ab_compare.py`, `m_fsi.py`,
`m_fsi_pct.py`, `m_rc.py`, `m_tagfrac.py`, `m_eta_ad.py`,
`m_tagged_prices.py`, `m_t1_virt.py`, `dump_preflight.py`, and their `*_pre.json` /
`*_post.json` / `.jsonl` outputs, `price_defaults_d6.json`, `warm_7li.txt`),
`ab/` (the 44 A/B `.npz` files and their logs), `logs/` (every gate's `.txt` / `.xml`, `ab_compare.txt`,
`report_regen.txt`, `report_check.txt`, `report_json_diff.txt`,
`price_defaults.txt`), `cell_map.py`, `junit_counts.py`, `fig_diff.py`,
`cites.py` (the citation inventory used to keep every cited line in place),
`gate.sh`, `baseline.sh`, `attr_run.sh`, `attr_cmp.py`, `moved_cites.py` and
`r6_moved.py` (the citation shift maps of §7), `pre_copy/` (a build of
`cb75ec7` for the example's d–p check), `smoke/`, `pre_REPORT.json` (sha256
`6ae3d5f6213be062…`, the committed one). Nothing was committed, published or
sent.

## 12. Verification

The adversarial verification of the tree I6 left (session scratch `v6/`). Its
verdict: "SOUND. D6 is implemented as designed." VmcAV18 is the run default,
defined in one place; the factories, `BreakupOptions`, `tagged.json` and the
⁶Li constants are bit-identical to `cb75ec7`; explicit-wave, inclusive and
coherent runs are column-identical before and after, and the pre-I6 `vmc` run
equals the post-I6 default. Every physics number, gate count and suite tally of
this record that the verifier re-ran reproduces. It found one medium and eight
low findings, all in documents, one regenerated figure, one error message and
this record; no defect in code correctness, physics numbers, tolerances or
verdicts. Worker F6 fixed them on 2026-10-29 (session scratch `f6/`). Each
finding was re-checked against the tree before anything was touched, and none
was refuted. F6 is line-neutral in every file a record or a gated document
cites by line. Its two additions sit where nothing is cited: the end of
tests/test_pipeline.cpp (the D6 MC-band case is the file's last), and inside
`draw()` of validation/price_defaults.py, from line 1098, below every line a
record cites there (plan_D7's 1068-1074, plan_D6's 1090). The verifier
labelled no finding with a lens; the lens column is F6's assignment.

| # | lens | finding | F6's re-check | what was done |
|---|---|---|---|---|
| V-1 (medium) | figures | The regenerated `docs/benchmarking/figures/price_defaults.png` silently dropped the largest D6 price. The priced alternative is now `hulthen`, whose A_zz^tag(0.30) ratio to the default is ×4.78404, above the plot's fixed log limit of 2.6, so that category had no marker. §6 said only that the series "follows". | Confirmed by viewing I6's PNG: the bottom panel's A_zz^tag(0.30) column is empty, and the ×2.32488 point at 0.20 sits just under the limit (−0.9532827 / −0.1992630 and −1.2069113 / −0.5191289). Before D6 the alternative was `vmc`, at ×0.209 and ×0.430, inside the range. | `draw()` now takes its y range from the points it draws and the f0 band, never narrower than 0.1 … 2.6, with ticks at 1-2-5 and a dated comment. `price_defaults.py --plot` was run once (N = 200 000, seed 20261026, 168.4 s). Its table equals I6's in every field but the wall times (f6/meas/price_defaults_f6.json against i6/meas/price_defaults_d6.json). The bottom panel now spans 0.1 … 5.98 with a tick at 5 and draws the ×4.78 point. Compared pixel by pixel with I6's figure, only rows 429–662, the bottom panel's plot area, differ; the top panel and the footer are identical. §6's sentence now says what was missed. |
| V-2 (low) | honesty | Three D6 statements in `docs/benchmarking/08_channel_priority.md`. (a) Ingredient I's "the two variants differ on 99.5 % of events" was unlabelled, although the design's §8 listed it. (b) Ingredient G still gave the lever as "×2.3–4.8 … under `--cluster-wave vmc`", and its development item still read "D2 at M2". (c) The D6 log row said "the shipped default moves A_zz^tag ×2.3 / ×4.8 … against the opt-in", where the default's is ×0.430 / ×0.209 of the opt-in's. | Confirmed by reading and on the record's numbers, which the verifier re-ran: 99.425 % on `hulthen` and 64.515 % on the default for the D3 recipe; hulthen over default ×2.32488 / ×4.78404, default over hulthen ×0.43013 / ×0.20903; the d–p dilutions hulthen over default ×1.02069 / ×1.01196. So the board's "every default-run number the design, the check and this phase found … labelled" did not hold for (a). | Line-neutral, each a single table row. (a) "… on the Hulthén α–d (2026-09-04; 99.4 % re-measured 2026-10-29) and on 64.5 % on the VMC default since D6 (20 000 events, seed 1234, 40 mb, toy SFs)". (b) The lever is "`--cluster-wave hulthen` against the VMC default (×0.43 / ×0.21 the other way, the `--cluster-wave vmc` lever until D6)". P_D is named per wave, and the d–p dilutions are given in the same direction (+2.07 % / +1.20 %; "D6 moved them −2.03 % / −1.18 %"). The development item reads "D2 at M2 — taken by D6 (2026-10-27), before M2". (c) "the shipped default moves A_zz^tag ×0.43 / ×0.21 … against the opt-in — the opt-in ×2.3 / ×4.8 against it"; the row's stage cell names F6. C3's ≥ 20 % criterion holds in either direction, so no rank moves. The board's completeness sentence now names (a) as F6's. |
| V-3 (low) | code | The missing-table message sits in the shared low-level reader `read_whole`, so readers that are not tagged channels also offered the tagged escape: an inclusive `--b1-model li6-convolution` run with an empty `LIPOLGEN_DATA_DIR` suggested `--cluster-wave hulthen`, and following that advice was refused. | Reproduced: the run printed "… or select the analytic forms, which read no file: --cluster-wave hulthen", and the same run with `--cluster-wave hulthen` was refused ("cluster_wave = 'hulthen' was asked for on channel inclusive …"). The readers built on `read_whole` are the tagged VMC waves and the AV18 deuteron, the `li6-convolution` α–d source and the §10 sampler (src/core/cluster_config.cpp), `DeuteronConvolutionB1`'s fdeut (src/core/b1_nuclear.cpp) and the `vmc-ft` C0 shape (src/core/rc.cpp). | Reworded, line-neutral (src/core/cluster.cpp lines 79-88). The `data_dir()` explanation comes first, then "On a TAGGED channel, whose run default reads the cluster-wave tables since D6 (CLUSTER_WAVE_DEFAULT), --cluster-wave hulthen (ClusterWaveSource::Hulthen) opens none of them; no other table (an inclusive --b1-model's, --rc-c0-shape vmc-ft's) has a stand-in". The condition names the channel the user runs, not the file. F6's first wording keyed it on the file, and was replaced before any gate: the inclusive run's missing file is `li6.ad`, itself a cluster-wave table. Re-run with an empty data dir: the inclusive `li6-convolution` run and the three default tagged channels print it and exit 1, and `--cluster-wave hulthen` runs on all three tagged channels (f6/meas/v3_emptydata.txt). No test matches the message. |
| V-4 (low) | honesty | The T1 timelike-X redraw fractions. docs/T2_CHAIN.md lines 112-113 gave "0.024 % / 0.022 % / 0.004 %" beside its own "19 / 11 / 1 of 50 k draws", which is 0.038 / 0.022 / 0.002 %. tests/test_breakup.cpp lines 416-417 gave the old values as the current Hulthén ones. USAGE line 1578's "0.01–0.05 %" excluded values it lists itself (0.054, 0.006 and 0.002 %). | Confirmed by arithmetic on the counts, which the verifier reproduced (50 k events, seed 20260830, the T0 and T1 `x` columns compared index by index: Hulthén 19 / 11 / 1, VMC 27 / 7 / 3). The 0.024 / 0.022 / 0.004 % were first published on 2026-08-29 (commit 7d08e04). | Line-neutral. T2_CHAIN.md gives 0.038 / 0.022 / 0.002 % on Hulthén with the counts, the 2026-08-29 values as "first published … on an earlier sample", and 0.054 / 0.014 / 0.006 % (27 / 7 / 3) on VMC. The test comment now carries 0.038 / 0.022 / 0.002 %, and USAGE reads "0.002–0.054 %". tests/test_breakup.cpp's twelve cases pass, 6 104 706 assertions, as under I6. |
| V-5 (low) | records | §7's "Three remain moved" missed plan_D6's citations of test_pipeline.py 77-82 at its lines 554, 653 and 762. §4's "88 ids renamed, mapped one to one" and "0 tests that ran stopped running" were loose. | Confirmed. At `cb75ec7` the Hulthén tag-fraction test sits at 77-82; it is now at 89-96, below the 12-line `tagged_run_hulthen` fixture, with two more docstring lines. plan_D6 declares its line numbers to be at `ca9b23a`, and test_pipeline.py did not change between `ca9b23a` and `cb75ec7`. Counted literally (f6/meas/knob_raw_ids.py on I6's base_full.xml and g2_full.xml): 42 ids only before (30 passed, 12 skipped) and 42 only after (34 passed, 8 skipped). Four ids kept their name and went passed → skipped: the coherent, inclusive-6Li, inclusive-7Li and inclusive-d `cluster_vmc_mc_sigma: 0 -> 1 on hulthen` cells. Their coverage moved to the same specs' renamed "… on the default (vmc)" cells, all four passing. 42 + 42 + 4 = 88 counts id changes, not renames. | §7 names the fourth moved citation. This is a record change only: the fixture was not moved, because plan_D6 is history and declares its tree. §4's sentence and table cell now say "42 ids renamed and 4 same-id swaps", and that "0 stopped running" holds under the per-spec map. |
| V-6 (low) | docs | `docs/PHYSICS_CHANNELS.md` row 520 listed `generate_tagged`'s flags without `--cluster-wave`, and did not say that `--beta` / `--pd` are now refused where the knob table says the run does not read them. | Confirmed against examples/generate_tagged.cpp (usage lines 18-20, `parse()` line 87). The verifier ran the example (`--beta` accepted on ⁷Li at T1 on either wave). | Line-neutral, inside the row: "`--cluster-wave hulthen\|vmc` (default `CLUSTER_WAVE_DEFAULT`, vmc since D6), `--beta`, `--pd` (each refused where the run's knob table says it is not read)". The row's citations are untouched, and the strict gate and R5 are unchanged (below). |
| V-7 (low) | code | The default VALUES derive from `CLUSTER_WAVE_DEFAULT`, but the TEXTS retype vmc, contrary to DD-3: the `--cluster-beta` / `--p-d` help, the `--cluster-wave` help's description, the bindings' docstrings, `validate()`'s message, the knob reasons and the example's usage. Nothing tied those texts to the constant; I4 added a `static_assert` for the same problem. | Confirmed by reading every listed site. cluster.hpp's own comment said "nothing retypes it", and the bindings' comment "READ by cli.py and the docstrings, never retyped". | Added `static_assert(CLUSTER_WAVE_DEFAULT == ClusterWaveSource::VmcAV18, "texts say vmc")` beside the constant (include/lipolgen/cluster.hpp line 420), line-neutral. The constant now fits on its line 419, and the comment above it (409-418, reflowed) lists the texts the assert ties, the Python help among them. The bindings comment (python/bindings.cpp 5163-5165) now says that cli.py READS the value and that the assert ties the texts. Checked in scratch (f6/sa/): a copy of the headers with the constant set to `Hulthen` fails to compile with "static assertion failed: texts say vmc" at line 420, and the tree's header compiles. The help strings were not rebuilt from `cluster_wave_name()`: their sentences are true only of vmc, so a derived token would not make them true after a revert. The assert forces the review instead, which is I4's pattern. |
| V-8 (low) | code | `validate()` accepted any `cluster_vmc_mc_sigma` on the default lithium α tags, NaN and 1000 included. NaN then failed at `Pipeline` construction with the misleading "the plan has no rate", and 1000 built silently. Pre-existing on explicit `vmc`, but reached by a default config since D6. | Reproduced through the bindings (the field assigned on a `make_config` result; neither `make_config` nor the CLI exposes it). NaN and ±inf passed `validate()` and failed at construction with "Pipeline: the plan has no rate". 1000 built, with the factory's α–d P_D at 0.13299. | `validate()` now refuses a non-finite value first, "PipelineConfig: cluster_vmc_mc_sigma must be finite (a number of the tables' printed sigmas)", on any channel and either wave. The change is line-neutral in src/core/pipeline.cpp: the C5.2 comment and the two existing throws were reflowed with their message text byte-identical, so the pinned use site at line 557 and the ALLOW / UNNAMED_ALLOW lines 847, 1116 and 1877 do not move. Pinned at the end of the D6 MC-band case, the file's last lines: NaN, +inf and −inf, each on `VmcAV18` and on `Hulthen`, six checks. Re-run: NaN and ±inf are refused by name on ⁶Li and ⁷Li α, and 1 and 1000 are accepted (f6/meas/v8_nan.txt). A finite bound was NOT added: none is declared anywhere, and choosing one is a new rule (§9 item 11). |
| V-9 (low) | docs | Labels left over after the relabelling: the Deeps docstring's "Where the default misses: data / fitted model runs 0.79, 0.95, 1.27, 1.72, 1.41" (the Hulthén's numbers) and its printout's "default wave, per (Q2, W*)". Registry row 11's "MOVED from the default (`Hulthen`)" read as "the default, Hulthén". | Confirmed. The printout's `head` is `m["waves"][DEFAULT_WAVE]["lc"]`, with `DEFAULT_WAVE = "hulthen"`; 0.79 … 1.41 is plan_D6 §4.3's hulthen column (AV18's is 1.1395 … 0.6408). Neither text is in a hashed block: the harness hashes only its vendored file. | Line-neutral. The docstring (lines 99-101) reads "Where the verdict wave (Hulthen) misses: … -- its tail falls too fast …"; the printout reads "verdict wave (%s), per (Q2, W*)" % DEFAULT_WAVE; row 11 reads "(i.e. `Hulthen`)". The REPORT is unchanged (`--check` OK), and test_bench_deeps.py passes 16 / 16. |

**Gates after F6** (the tree F6 leaves; uncommitted):

- **cmake build**: after the code edits 1 min 48.5 s (near-full, since
  `cluster.hpp` changed), after the message rewording 18.1 s; no warning and no
  error.
- **`./build/lipolgen_tests`**: 422 cases, 422 passed, 0 skipped, **17 244 315
  assertions**, 473.6 s wall. That is +6 on I6's 17 244 309, the six NaN / ±inf
  checks of V-8: the four D6-aware cases now count 69 assertions against I6's 63.
- **The targeted pytest files** (with `LIPOLGEN_BENCH=1`): test_price_defaults,
  test_bench_deeps, test_doc_link_gate, test_cli, test_module,
  test_spdx_headers and test_meta_provenance, 267 passed in 41.9 s, each file
  at I6's count.
- **The full pytest**: **1472 passed, 149 skipped, 0 failed, 543.7 s** (wall
  544.6 s). Per file it is identical to I6's final run in all 42 files, and per
  test there are 0 outcome changes (f6/logs/full.xml against i6/logs/g5_full.xml).
- **The bench subset** (`LIPOLGEN_BENCH=1 … -k bench`): **351 passed, 1 skipped
  (the environment skip, as before), 0 failed, 340.9 s**, per test identical to
  the `cb75ec7` baseline (i6/logs/base_bench.xml). The report test passes:
  F6 moved nothing in the report.
- **`make_report.py --check`**: OK (75.2 s). REPORT.md, REPORT.json and the two
  label-changed figures are byte-identical to I6's. Of the generated files, only
  `price_defaults.png` changed (V-1).
- **The citation gates**, after the code and document edits and again with this
  section in place. Strict: **1241 / 95 / 7 / 0 / 24** (SPIN32 19 / 115 / 0 / 0 / 10;
  PYTHIA_BRIDGE 0 / 0 / 8 / 0 / 0), unchanged, with no `--fix` and no
  `--record-ranges`. `--audit-ranges`: **209 on 179 blocks**, 1 / 66 / 39 / 103,
  **0 REFUSED**. `--records-only`: **1432 citations in 77 records**, 933 / 61 /
  110 / 129 / 199, **0 broken**; this section adds no gated citation. No
  `(as of …)` annotation was needed and no census was re-pinned.
  test_doc_link_gate.py, test_meta_provenance.py and test_spdx_headers.py:
  126 passed (26.5 s) with this section in place.
- **SPDX**: 133 files. **Line neutrality**, each edited file against
  f6/pre_f6/: F6 adds lines only at the end of tests/test_pipeline.cpp (+14),
  inside `draw()` of validation/price_defaults.py (+10, from line 1098) and in
  this record. Every other edit keeps its file's line count. The pinned and
  allow-listed lines (src/core/pipeline.cpp 557, 847, 1116 and 1877) and the
  cited names of include/lipolgen/cluster.hpp (406, 427, 431) hold their text.

Scratch `f6/`: `edits/` (`apply.py`, which refuses a missing or ambiguous old
string; `step_a_code.py`, `step_a2_message.py`, `step_b_docs.py`,
`step_c1_record.py`, `step_c2_record.py`, `step_d_status.py`, every edit an
exact old/new pair) and `pre_f6/` (every file as I6 left it, and I6's diff,
sha256 `8102c8de…`, the verifier's); `meas/` (`knob_raw_ids.py`,
`price_defaults_f6.json`, `v3_emptydata.txt`, `v8_nan.txt`); `sa/` (the
`static_assert` check); `logs/` (`build_a.txt`, `build_a2.txt`, `cpp_full.txt`,
`full.xml` / `.txt`, `bench.xml` / `.txt`, `targeted.*`, `gate_*.txt`,
`report_check.txt`, `price_defaults.txt`); `emptydata/`. Nothing was
committed, published or sent, and the sibling repository was not opened.
