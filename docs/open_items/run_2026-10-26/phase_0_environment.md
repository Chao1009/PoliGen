<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — phase 0: the environment, built and measured

Written by the supervising session before any worker started. Everything
here is a measurement on this machine (Ubuntu 22.04, gcc 11.4.0, cmake
3.22.1, Python 3.10.12, numpy 2.2.6, scipy 1.15.3, 4 cores, 15 GB).

## The dependency stack

`../deps/install` did not exist (a fresh checkout at `0cea6c0`). Built with
the project's own `.github/scripts/build_deps.sh` (unmodified, `jobs = 4`)
into `/home/cpeng/Projects/polli/deps/install`: HepMC3 3.3.0 (sha256
verified), LHAPDF 6.5.5 (sha256 verified) with **CT18NLO, NNPDFpol11_100 and
EPPS21nlo_CT18Anlo_Li6 installed**, PYTHIA 8.317 (its tarball has no pin in
the script, as the script says) with `--with-python-config`; **439 MB**;
started 17:16 UTC, done 18:02 UTC (46 min wall on 4 cores, against 337 s on
the author's 8 cores on 2026-09-05 — this machine was also running the
documentation reading in parallel). `pythia8.so` (the Python module) is in
the prefix's `lib/`.

## LiPolGen on it

`cmake -S . -B build -DCMAKE_BUILD_TYPE=Release` + `cmake --build build -j4`:
2 min 1 s, all three optional tiers ON (`lipolgen 0.1.0 LHAPDF True HepMC3
True PYTHIA True`).

| gate | measured 18:04–18:12 UTC | the 2026-09-27 machine (no LHAPDF set, PYTHIA 8.312) |
|---|---|---|
| `build/lipolgen_tests` | **416 cases / 416 passed / 0 failed / 0 skipped; 17 242 712 assertions, 0 failed** (7 min 28 s) | 416 / 403 / 13 failed — all 13 "Info file not found for PDF set" |
| `python3 -m pytest python/tests -q` | **1269 passed / 158 skipped / 0 failed** (7 min 47 s; the first run of the supervising chain called `python`, which this machine does not have, and was re-run with `python3`; the `import_pythia8` edit below landed while it ran, so the T4/T6 tests may have seen either state — both pass) | 39 failed / 1209 passed / 163 skipped — all 39 the same cause |
| `make_report.py --all --compare docs/benchmarking/REPORT.json` | **5 pass / 8 fail / 5 blocked / 0 blocked (environment) / 0 error**, 38 changes — the LHAPDF sub-rows measured for the first time (below); T4 and T6-C14 BLOCKED by the PYTHIA python module (below, fixed) | 6 / 8 / 3 / 1 / 0 |

So the 13 + 39 LHAPDF failures the previous report carried are gone, as
PLAN.md phase A1 expected, and no other C++ failure appeared.

## What the LHAPDF sets showed, first numbers ever on a LiPolGen build

From the `--compare` output (`phase_A_baseline.md` carries the full list):

- `t1_nmc_f2d[ct18nlo]`: BLOCKED → **FAIL**, χ²/ndf = 3401.8/146 (p = 0),
  mean model/data **1.086** on the DIS cut, worst pull 12.3 at x = 0.05,
  Q² = 2.5 GeV². CT18NLO as this tree evaluates F₂ sits **8.6 % above** the
  NMC deuteron data — worse than MSTW2008 LO (416.5/155, mean ratio 1.025).
  A reading to be checked by the workers, not a verdict: `LhapdfSF` builds
  F₂ from an NLO PDF with a leading-order formula, so the excess is of the
  size of the missing NLO coefficient functions and target-mass / higher-twist
  terms at Q² ≈ 2–5 GeV², not necessarily of the PDF. For
  `08_channel_priority.md` ingredient A this means **no shipped unpolarised
  backend passes the absolute F₂ᵈ data**; the ratio does.
- `t1_nmc_f2d_over_f2p[ct18nlo]`: BLOCKED → **PASS**, 215.5/189 (p = 0.090).
- `t1_g1d_world[ct18nlo]` (ToyG1 on CT18NLO's F₁): FAIL — COMPASS 20.2/11
  (p = 0.043), HERMES 64.3/11, E143 59.0/24; `t1_g1d_world[nnpdfpol]`
  (NNPDFpol1.1 g₁): **FAIL** — COMPASS 45.0/15 (p = 7.7e−5), HERMES 52.6/15,
  E143 51.0/28 (p = 0.005). A fitted polarised PDF failing the data it was
  fitted to at p ≥ 0.01 points at the harness's deuteron construction (the
  D-state factor the 2026-09-26 run left UNVERIFIED — PLAN.md phase D1) or at
  the F₁ it divides by, before it points at the PDF.
- `t5_sum_rules[bjorken[nnpdfpol11]]`: BLOCKED → **FAIL**: within 9.6 % /
  12.4 % / 19.3 % of g_A/6·C_Bj at Q² = 10 / 5 / 2 against the declared
  "converged and within 10 %" — converged (last decade ≤ 0.4 %), outside 10 %
  at the two lower Q².
- `t3_nmc_li6_over_d`: blocked (environment) → **PASS** again (the 2026-09-23
  numbers).
- `t2_deuteron_static`: the two `q_d_quad_error` sub-values moved in their
  15th digit (7.38e−15 → 7.44e−15; 1.22347e−13 → 1.22291e−13) — numpy
  2.2.6 against 2.4.6; no verdict or headline moved.

## PYTHIA 8.317's Python module cannot be imported from the repository root

`import pythia8` from `/home/cpeng/Projects/polli/PoliGen` raised
`ImportError: initialization failed`, with the chained cause `SyntaxError:
invalid decimal literal` at `File "<string>", line 183: /// 10LZZZAAAI carry A
in digits 4..6.` — a line of **`python/bindings.cpp`**. The module's
initialiser, `Pythia8::loadExtraPython` in
`plugins/python/src/PythiaPython.cpp` of the 8.317 tarball, opens the
directory literally named `python` **relative to the current working
directory** and `pybind11::exec`s every file in it as Python. From `/tmp`
or from `validation/` the same module imports; from any directory with a
`python/` entry it does not. `make_report.py` runs every harness with
`cwd = root` and pytest runs from the root, so the T4 row and T6's clause C14
(the PDG charge table) read BLOCKED (environment) on this machine's first
report. A PYTHIA defect (the 2026-09-26 run's conda-forge 8.312 did not show
it); worked around in ONE place, `validation/benchmarks/_report.py`
`import_pythia8()` — a chdir to a fresh temporary directory around the one
import, restored on every exit path — and used at the four import sites
(`t4_pythia_ep_closure.py` ×2, `t6_hepmc3_convention.py` ×2) and in
`python/tests/test_bench_t4_t6.py`. No number, tolerance or verdict text
changed; the measurement after the fix is in `phase_A_baseline.md`. Worth
reporting upstream to the PYTHIA authors.

## Reachability that the previous runs did not have

- HEPData's JSON API answers the standard Python `requests` client
  (`https://www.hepdata.net/record/<ins>?format=json`); `curl` without a
  client identity receives 403 and **nothing here spoofs one** (the auto-mode
  classifier refused the attempt; the rule is in PLAN.md).
- APS journals serve PDFs and supplemental files to `curl`'s own identity
  from this network (institutional licence): the BONuS Supplemental Material
  (`BONuSprcSupMat2.txt`, 88 845 bytes), Phys. Rev. 144 (1966) 859, Phys.
  Rev. 162 (1967) 992, PRC 41 (1990) 898, PRA 20 (1979) 381, PRD 42 (1990)
  2377 and 1424, PRL 26 (1971) 445 were fetched into the session scratch
  directory, **outside the repository**; only numbers are vendored.
  Elsevier (ADNDT 14) returns a 2.7 kB redirect page to a non-browser client
  and was not pursued.
- `gfortran` 11 is installed; `github.com/spiesber/DJANGOH` at `e356b84`
  (the commit the 2026-09-26 review read) was cloned into the scratch
  directory for the row-7 worker.
- `docker` (snap) and `singularity` exist; no eic-shell image is mounted;
  the ePIC chain legs are not attempted (PLAN.md, "Not in this run").
