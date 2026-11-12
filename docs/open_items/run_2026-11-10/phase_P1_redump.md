<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-10 — phase P1: the re-dump (decision J2)

Written by the P1 implementer (an Opus worker under the supervising session)
on 2026-11-10, on base commit `dd7f4d4` with a clean tree. It implements
**J2, the author's decision** ("re-dump"; `DECISIONS.md` J2). Scratch files
live in the phase's scratch directory (`run_1110/p1/`) and are named below by
file.

## 0. What changed, in one paragraph

`validation/dump_polligen_reference.py` was RUN in place with the tree's
`env.sh`, and nothing was edited by hand. All three preflights passed. The
reference tables in the tree are now exactly what that run wrote. Six JSONs
moved: 3170 numbers at the rounding level. `_manifest.json` lost its one
stale string. `beams.json`, `b1_default_li6.json` and `validation/README.md`
are byte-identical to `dd7f4d4`. A second run reproduced the first byte for
byte, and both are byte-identical to the 2026-11-06 dump of
`../run_2026-11-03/phase_P5_port_gate.md` (first pass). Every gate passes at
its present tolerance, and no tolerance, test, source or document moved. The
largest move, 9.4e−13 relative in `xsec.json`, has a measured mechanism (§1).
One INPUT of the dump changed by one ULP: the magic-angle spin axis
θ = arccos(1/√3). The kernel rows at that axis follow it. LiPolGen reads the
angle off the file, so the C++ and Python gates track it, and the re-dumped
tables use less of the 1e−12 tolerance than the committed ones did (§2).

## 1. What moved, before -> after

The recipe was `source env.sh && PYTHONDONTWRITEBYTECODE=1 python3 -B
validation/dump_polligen_reference.py` (exit 0, 0.51 s, `dump1.log`). Every
reference JSON was compared, pre (`pre/`, copied before the run) against
post, by `diff_stats.py` (`diff_stats.log`). "Numbers" are int and float
leaves. ULP is |Δ| / ulp(pre).

| file | numbers compared | changed | max \|Δ\| | max rel, \|v\| > 1e−6 | max rel overall (on value) | > 1 ULP | > 2 ULP |
|---|---|---|---|---|---|---|---|
| `xsec.json` | 8655 | 2878 | 1.4e−9 | 9.36e−13 (`results/7Li/kernel_scenario_rank2/amplitudes[14]/w_avg[1]`, −1.0086e−5) | 3.8 (on 1.28e−19, a `kernel_tensor_gamma` `w_avg`) | 1761 | 999 (877 among \|v\| > 1e−6, max 5570 ULP) |
| `spin.json` | 3350 | 202 | 1.85e−16 | 1.41e−15 (12 ULP, `j1` `lab_moments/tensor_zz` = −0.0589) | ∞ (4 values that were exactly 0) | 117 | 96 (8 among \|v\| > 1e−6) |
| `tagged.json` | 2906 | 51 | 4.4e−16 | 2.19e−16 (`boost_spectator` of `li6_alpha`, `li7_alpha`, `deuteron`) | 2.19e−16 | 0 | 0 |
| `spectator.json` | 1583 | 28 | 4.3e−19 | 2.2e−16 (`boost_fragment`, seven channels) | 2.2e−16 | 0 | 0 |
| `coherent.json` | 444 | 9 | 2.8e−17 | 2.05e−16 (`tag_acceptance*`, `fragment_rigidity`) | 2.1e−16 (on 5.0e−7) | 0 | 0 |
| `bookkeeping.json` | 282 | 2 | 6.9e−18 | 2.78e−16 (`helicity_flip_j32_maxent_anchor` populations) | 2.78e−16 | 2 | 0 |
| `beams.json` | 101 | 0 | – | – | – | – | – |
| `b1_default_li6.json` | 2607 | 0 (not this script's file) | – | – | – | – | – |
| `_manifest.json` | 0 | 1 string: `could_not_call[0].reason` `tagged.py:299-302` → `tagged.py:849-852` | – | – | – | – | – |
| total | 19928 | 3170 | | | | | |

No key was added or removed. No other non-numeric value changed.
`git diff --numstat` gives 2878 / 202 / 51 / 15 / 6 / 2 / 1 lines. The
`li6_alpha_vmc` block of `tagged.json` (A13's splice) is equal to the
committed one both parsed and as text (15 lines). `validation/README.md` is
byte-identical to `dd7f4d4`.

**Against the picture in `phase_P5_port_gate.md` §1.** It matches: 3170
numbers in the same six files with the same per-file counts and maxima, and
the same manifest string. That section gave only an absolute bound for
`spin.json`. Measured here: 96 of the 202 moves exceed 2 ULP, and 8 of
those sit on values above 1e−6, up to 12 ULP (1.41e−15 relative).

**Where the 9.4e−13 comes from.** The axis list `AXES` of the dump script
holds `np.arccos(1.0 / np.sqrt(3.0))`. It was 0.9553166181245091 in the
committed `spin.json` and `xsec.json` and is 0.9553166181245092 today, one
ULP higher. Against a 40-digit `mpmath.acos`, the old value is 0.965 ULP
from the true one and today's is 0.035 ULP from it: today's is the correctly
rounded result. In `xsec.json` the rows at θ_s on that axis carry 776 of the
2878 changes and the worst relative move (9.36e−13). Every other row of the
file changes by at most 5.6e−14 among |v| > 1e−6 (2101 changes). Which
library changed the arccos (the libc6 2.35-0ubuntu3.13 → .15 upgrade of
2026-10-26 in `dpkg.log`, numpy 2.2.6 being unchanged since 2026-06-12) is
the likely reading. It is not proven, because the old library is not on
the machine to A/B.

**The environment of the dump** (`environment.txt`). Python 3.10.12
(`/usr/bin/python3`). numpy 2.2.6 (user site, OpenBLAS 0.3.29
scipy-openblas DYNAMIC_ARCH; SIMD baseline SSE/SSE2/SSE3, found up to AVX2,
no AVX512). scipy 1.15.3, which polligen's import pulls in. glibc
2.35-0ubuntu3.15, libpython3.10 3.10.12-1~22.04.18, kernel
6.8.0-138-generic, AMD Ryzen 9 5900HX. The sibling is PolarizedLithiumSim
`472681248b705a61c3132d784bea75fe2ecde09d` with 0 lines of
`status --short`, before and after. The manifest carries `polligen_version`
and no tool versions, so the environment is recorded here alone. The
script is unchanged.

## 2. What was verified and how

- **Determinism.** The tree's output after run 1 was copied to `dump1/`,
  and the script was run a second time (`dump2.log`, exit 0). `cmp` found
  all nine JSONs, `.keep` and `validation/README.md` identical to `dump1/`.
  The two logs differ only by the wall-time line that `/usr/bin/time`
  added to run 1. All nine JSONs and the README are also byte-identical to
  the 2026-11-06 dump (`run_1103/p5/dumped/`), so the output is stable
  across processes and days on this machine. File modes are unchanged
  (664), and no staged temporary file was left in `validation/reference/`.
- **Headroom: how much of the 1e−12 the port gates use.** This was measured
  pre against post with two read-only scratch harnesses.
  - C++ (`headroom_cpp/`): the six reference-reading test files,
    `test_reference`, `test_spectator`, `test_coherent`, `test_bookkeeping`,
    `test_tagged` and `test_cluster`, were compiled unmodified against the
    tree's libraries. A copy of `tests/check_close.hpp` records
    |got − want| / (atol + rtol·|want|) for every comparison at
    rtol ≤ 1e−12, and the reference directory comes from the environment.
    The harness runs the same 93 cases / 527 783 assertions as the tree
    binary's subset, all passing on both file sets.
  - Python (`headroom_py/`): copies of `test_kernel_reference.py` and
    `test_port_gate_li6_vmc.py` were run with `pytest.approx` at
    rel ≤ 1e−12 and the module's `_close` instrumented. 15 passed on each
    set.

  The largest fraction of the tolerance used, before → after:

  | gate (reference file) | C++ | Python |
  |---|---|---|
  | `xsec.json` | **0.056 → 0.0022** (5.6e−14 → 2.2e−15 relative) | 0.056 → 0.0022 (`test_tensor_gamma_block` 0.056 → 0; amplitudes ≤ 0.0089 → ≤ 0.0022) |
  | `spin.json` (1e−12 and 1e−13 maxent, absolute floors 1e−14 / 1e−15) | 0.042 → 0.042 | – |
  | `tagged.json` (model grid / `boost_spectator`) | 0.0061 → 0.0061 / 0.00022 → 0 | – |
  | `bookkeeping.json` (plans / ladder) | 0.0028 → 0 / 0.00019 → 0.00019 | – |
  | `spectator.json` (`_boost_fragment`) | 0.00022 → 0 | – |
  | `coherent.json` (scenario formulas) | 0.00021 → 0 | – |
  | `beams.json` | 0 → 0 | – |

  `validation/repin_tagged_from_lipolgen.py` prints the same tagged worst
  ratios on the re-dumped file: 6.057e−3 / 9.591e−4 / 2.977e−3, and the VMC
  block 5.389e−9 / 3.196e−9 against its 1e−8. The cases are excluded from
  the Python column. `test_the_port_gate_reads_the_coefficient_off_the_file`
  probes `approx(known)` for both 1 and 0.5 by design, so one probe always
  misses. The VMC dilutions are compared at 1e−8, outside this filter.
  **Why the axis move costs nothing:** the gates rebuild every row from the
  file's own θ_s. At the worst point LiPolGen returns the pre value exactly
  against the pre file and the post value exactly against the post file
  (`pre_all.log`, `post_all.log`).
- **No harness reads `validation/reference/`.** `grep -rn "reference/"
  validation/benchmarks/` matches one comment in `_report.py`, a JSON key
  path and not this directory. The report was not regenerated.
- **Documents.** In `docs/*.md`, `docs/benchmarking/*.md`, `README.md` and
  `validation/README.md` no sentence states, as a present fact, that the
  dump does not reproduce or that the tables landed by splice. The mentions
  are in dated run records and in the 2026-09-23 history paragraph of
  `docs/benchmarking/BENCHMARK_PLAN.md`. So no `J2 (2026-11-10)` clause was
  needed, and none was added.

## 3. Files

- `validation/reference/xsec.json`, `spin.json`, `tagged.json`,
  `spectator.json`, `coherent.json`, `bookkeeping.json`, `_manifest.json`:
  as the dump wrote them.
- This record, and `STATUS.md` (the P1 row).

Nothing else changed, including `validation/README.md`, the dump script, the
tests, any tolerance and the sibling.

## 4. Gates (against the re-dumped files)

| gate | result | wall |
|---|---|---|
| the dump's three preflights (twice) | passed (cos 2φ coefficient 1, the LOI's; defaults Hulthén) | 0.51 s per dump |
| `cmake --build build -j2` | nothing rebuilt (the JSONs are read at test time), 0 warnings | 0.4 s |
| `lipolgen_tests` full | **435/435**, 17 244 595 assertions | 8 m 28 s |
| pytest A / B / C | 358 + 3 skip / 136 + 1 skip / 186 | 191 s / 112 s / 83 s |
| pytest D / E (`test_knob_provenance` halves 339 + 338 of 677) | 249 + 90 skip / 277 + 61 skip | 425 s / 293 s |
| pytest F (18 files) | 317 | 179 s |
| pytest sum | **1523 passed, 155 skipped, 0 failed** | ≈ 21 m 23 s |
| bench subset (`LIPOLGEN_BENCH=1 … -k bench`) | **368 passed, 1 skipped** | 7 m 16 s |
| `repin_tagged_from_lipolgen.py` (check-only) | 3 Hulthén blocks + the VMC block agree; nothing written | 0.17 s |
| citation strict | 1243 / 95 / 7 / 0 broken / 24 | 3.2 s |
| citation `--records-only` | 1432 citations in 92 records, 0 broken | 0.3 s |
| citation `--audit-ranges` | 209 on 179 blocks, 0 REFUSED | 2.9 s |
| `check_spdx_headers.py` | 136 files, all carry the line | 0.04 s |
| `test_doc_link_gate.py` | 112 passed | 51 s |
| report | not regenerated: no harness reads `validation/reference/` | – |
| re-run after this record and the board row were written: the three citation modes, SPDX, `test_doc_link_gate.py` | the same, except `--records-only` now counts 1432 citations in 93 records (this record added), 0 broken; 112 passed | 33 s (doc-link) |

## 5. Open points

1. **The cause is located, but its library is not proven.** The one-ULP
   arccos change on the magic-angle axis is measured, and today's value is
   the correctly rounded one. The library that changed is inferred from
   `dpkg.log`.
2. **Upstream robustness (the sibling's, not this tree's).** A dump input
   that depends on libm's last bit of `arccos` will move again on the next
   libm change. The gates are immune, because they read the angle off the
   file. Only the byte-reproducibility of the tables is exposed. No change
   is proposed here.
3. `phase_P5_port_gate.md` open point 5 (the splice not reproducible from
   one run) is now moot. `tagged.json` is one dump run's output, and its
   `li6_alpha_vmc` block is unchanged.
