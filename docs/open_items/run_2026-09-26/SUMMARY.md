<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-09-26 — one page

**What this run is.** Two pieces of work on `5e5140a` (run 2026-09-23,
committed): a code review of that commit, applied; and the remainder of M1
(`../../benchmarking/BENCHMARK_PLAN.md` §3) implemented — the T1 harnesses,
four more tiers' first harnesses, and the `REPORT.md` generator. Two build
commits came first and are committed (`302b57a`, `f0a8f1e`); everything else
is **one commit** on base `f0a8f1e`, made by seven agents in one shared tree
and integrated by a closing stage, then checked by four adversarial verifiers
and a fix stage (dates in UTC; integration, verification and fixes ran on
2026-09-27; see "Verification" below). **No physics default, registry status label, constant,
tolerance or reference-JSON number moved, and no decision was taken.** No C++
source or binding changed except one comment-only header edit
(`include/lipolgen/coherent.hpp`, line count kept); nothing was rebuilt,
pushed, published or sent.

The environment it ran in is **not** the development machine: PYTHIA 8.312
from conda-forge (the project pins 8.317), HepMC3 3.3.1, LHAPDF 6.5.6 with
**no PDF set installed** (none downloadable: hepdata.net, arXiv, INSPIRE,
CERN and hepforge are blocked; PyPI and conda-forge are reachable), a
shallow clone of 51 commits, no `/usr/bin/time` (every runtime here is
`time.perf_counter` around a subprocess, interpreter start-up included).

## The review of `5e5140a`

Six lenses (t2, t3, t4/t5 harnesses; the bench tests; the dump scripts; core
and build), each finding verified independently twice: **30 findings**, 25
distinct (0/13/27, 1/10, 5/14 and 8/17 are one issue found two or three
times); one major confirmed by both verifiers (26), one refuted by both (16),
five with a split verdict (4, 7, 11, 24, 29). **27 applied, 3 skipped.**

| # | area | what | disposition |
|---|---|---|---|
| 0, 13, 27 | bench tests | the HERMES χ² pin failed whenever the MSTW grid was absent (χ²(A_zz) depends on the F₂ backend) | fixed: χ²(b₁) pinned unconditionally, χ²(A_zz) only under `MstwSF` (else a skip naming it); a new test forces ToyF2 and finds every χ²(b₁) identical |
| 1, 10 | harnesses | "an altered vendored file exits 2" was implemented only for row 4; an altered EPIOS table still passed | fixed: every harness now records a sha256 of each vendored file (taken at `f0a8f1e`) and re-checks it on every read; a pytest per file flips one byte. The digest covers every byte **except in row 4** (`t3_nmc_li6_over_d`: the body after its provenance marker; not a file this run edited). Row 5 hashed only its ⁶Li row until the verification stage, which added the whole-file digest (the first wording here said "every byte" of every harness; a verifier refuted it) |
| 2 | t2 | BONuS crashed (exit 2) on a non-UTF-8 stdout | fixed; exit 0 under latin-1, ascii and the C locale, a subprocess test |
| 3 | t2 | the Deeps census in the docstring was wrong | fixed in the docstring and printout (60 blocks at Q² 1.8 / 2.8; 55 with unverified labels, 25 triples duplicated; mid 78 empty; 60 value = stat = 0 rows are missing, not zeros); the file's numbers not touched |
| 4 | t3 | a C0 with no zero crashed the FB harness | fixed: a recorded FAIL with its reason |
| 5, 14 | bench tests | the NMC table guard never ran without the LHAPDF sets, and any `Tree()` error became a skip | fixed: a module-only fixture; the LHAPDF fixture skips only on "Info file not found", naming each set |
| 7 | t4 | the DJANGOH unblock recipe ignored the samples' W_h ≥ 3 GeV cut | fixed (wording): the cut removes the elastic tail whatever IEL is; DJANGOH's `HSPRLG` zeroes IEL2, IEL31..33 when WMIN > M_p (public source `github.com/spiesber/DJANGOH` at `e356b84`, read, not vendored), so the row unblocks only with IEL31..33 ≠ 0 **and** WMIN ≤ M_p |
| 8, 17 | t4 | the ROW line printed `np.float64(...)` | fixed |
| 9 | t4 | a missing numpy exited 1 | fixed: 2, with a subprocess test |
| 12 | t5 | EPIOS counted modes 5, 6 (where the tree throws) as EST matches | fixed: `est_modes_matched = [0]`, 1 of the 6 modes it computes; "3 of 8" recorded as what it was |
| 15 | bench tests | a tautological assertion | replaced by `pzz_cartesian == est.closed_form(0.530)` |
| 18 | bench tests | the as-of tests failed in a shallow clone | fixed: they skip, naming the missing commits, only when the clone is shallow **and** a named commit is absent |
| 19, 20 | dump scripts | the check-only re-pin used 1e-12 where `tests/test_tagged.cpp` uses 1e-9, and missed NaN / missing / extra entries | fixed: the C++ gate's tolerances and pass rule; any NaN, infinity, key or length mismatch fails |
| 21 | dump scripts | the `tagged.json` README gave movement figures from channels not in the file | re-measured on the committed files (`n_of_kc` up to +729 %, `struck_populations` 0.83 absolute, `p2_moment` sign flip ×1.38) |
| 22, 23 | dump scripts | the b1 self-pin's label was hard-coded; a failing build truncated the pin | fixed: the label is kept only while the numbers are byte-identical, moved numbers need `--accept-changed-numbers`; `--check`; atomic write |
| 24, 25 | dump scripts | the polligen dump wrote file by file; two `tagged.py` lines cited for one error | fixed: the whole set or nothing; one citation |
| 29 | core | the ESTARLIGHT SCOPE comment named the wrong symbol for the 0.7 GeV² floor | fixed, comment only, 792 lines kept |
| 26 | build | pybind11 < 2.12 with NumPy 2 returns stride-0 event columns | fixed in `f0a8f1e` (floor 2.12 in `pyproject.toml`, `README.md` and `find_package`) |
| 28 | build | scikit-build-core 0.9 / 0.10 reject `pyproject.toml` | fixed in `f0a8f1e` (floor 0.11) |
| 6 | t3 | the UVa FB file carries an SPDX GPL tag over third-party content that states no licence | **skipped** — a licence-tag nit; its header already says "Licence: NONE STATED"; a `LicenseRef-` tag is the author's call. *(Corrected at verification: the skip first also said "the file's bytes are now sha256-pinned by its harness", which was false — the harness hashed only the ⁶Li row, so the tag was pinned only by `REPORT.json`'s `data_sha256` through `--check`. Since the fix stage the harness hashes every byte, the tag included, and a test flips it.)* |
| 11 | t5 | the EST row would FAIL if the tree ever returned the correct limit at \|P_z\| = 1 | **skipped** — refuted by one verifier: the row pins the tree's documented throw (`populations_maxent` needs \|p_z\| < 1), and no such tree exists |
| 16 | bench tests | row 5's "verdict follows from the numbers" asserts "cannot fail" | **skipped** — refuted by both verifiers |

Build and portability fixes found by running here, committed before the
review was applied: `302b57a` — LHAPDF 6.5.6 rejects `lhapdf-config --incdir`
(CMake captured the help text as an include directory); `f0a8f1e` — the test
target builds with any optional tier off (375/375 cases with all four off),
the pybind11 and scikit-build-core floors above, and numpy 2.4's removal of
`np.trapz` (one test and three scripts use the `vmc_reconcile.py` shim).

## The M1 remainder: what was measured

Every harness below re-ran three times on 2026-09-27 with identical rows;
none was tuned; each tolerance is in its header, with a note wherever the
number had been previewed before it was chosen. Details and sub-rows:
`validation/benchmarks/README.md` and `../../benchmarking/REPORT.md`.

| row | harness | result | the number, with its tolerance | runtime |
|---|---|---|---|---|
| T1 F-1 | `t1_nmc_f2d` | **FAIL, recorded** | NMC F₂ᵈ, 155 points on Q² ≥ 1, W² ≥ 4 GeV²: shipped ToyF2 χ² = 8307.9/155 against p ≥ 0.01 (diagonal); MSTW 416.5/155; CT18NLO blocked (environment) | 0.78–0.85 s |
| T1 F-2 | `t1_nmc_f2d_over_f2p` | **FAIL, recorded** | NMC F₂ᵈ/F₂ᵖ, 211 points: ToyF2 2161.4/211; MSTW 247.2/211 (p = 0.044; 260.7 via the Q²-frozen hook); CT18NLO blocked | 0.78–0.88 s |
| T1 D-2/D-1/D-6 | `t1_g1d_world` | **FAIL, recorded** | g₁ᵈ per nucleon, ToyG1 as g₁_nucleus(d)/2: COMPASS 107.8/15, HERMES 69.4/15, E143 54.9/28; MSTW-F₁ 53.2 / 47.1 / 50.2; NNPDFpol11 and CT18NLO blocked; the data's D-state convention UNVERIFIED | 0.70–0.81 s |
| T2 §5.2 | `t2_deuteron_static` | **FAIL, recorded** | η_d = 0.0256(4) within 2σ for each shipped deuteron: CD-Bonn −0.07σ, AV18 −1.39σ, **Hulthén (the tagged default) +17.2σ**; Q_d −5.4 / −5.7 / +1.6 % recorded, not asserted | 1.3–1.5 s |
| T3 | `t3_li6_rms_devries` | **PASS** | r_ch = 2.57100 fm against de Vries 2.54(5) / 2.56(5) / 2.57(10), within 2σ of each — the point radius is a tree input derived from Angeli–Marinova, so this is a consistency check of its folding, not a prediction | 0.10–0.13 s |
| T4 D-1 | `t4_pythia_ep_closure` | **PASS** (opt-in) | bridge vs stock PYTHIA 8.312: worst \|ratio − 1\| = 0.023 against 0.20 over 15 sub-rows | 9.2–10.2 s |
| T5 | `t5_sum_rules` | **FAIL, recorded** | BC closure 45/45 cells, worst 0.667 of its bound (PASS); Bjorken on ToyG1 0.6430 against 0.1854 ± 0.0023 at Q² = 5, not converged (FAIL); Bjorken NNPDFpol11 blocked (environment); FSI unitarity blocked (not definable) | 1.4–1.6 s |
| T5 | `t5_weighted_unweighted` | **PASS** | same point both routes: pull +0.359 against \|pull\| ≤ 3 (89.1934 ± 0.0512 vs 89.1674 ± 0.0516 pb) | 2.2–2.5 s |
| T6 | `t6_hepmc3_convention` | **FAIL, recorded** | 21 of 23 clauses on 800 events in 8 files; C13, C14 fail on the inclusive records (the four-momentum 1e-9 read relative to the summed beam energy, looser on p_x, p_y than `tests/test_hepmc.cpp`'s rule; the same verdict under either) | 3.7–4.2 s |

The ten rows of 2026-09-23 re-ran with the same verdicts and numbers (the
EPIOS reading is now 1 of 6, above); row 4 is `blocked (environment)` here
and keeps its 2026-09-23 PASS. `make_report.py --all` over all 17 harnesses:
**5 PASS, 8 FAIL recorded, 3 BLOCKED, 1 blocked (environment), 0 error**,
43.97 s wall. `tiers.py` turns §7 into the harness list a diff re-runs.

**The integration stage** (2026-09-27) changed no number: `t2_deuteron_static`
now exits 2, not 1, when numpy is missing (the DJANGOH pattern);
`test_doc_link_gate.py`'s external-citation count is all-or-none over both
documents' externals (tallies below); the plan's U-3 row escapes its `|t|`
(it split the table row); and it wrote the README, the plan's status notes,
dated notes where a sentence had gone stale (`PHYSICS_CHANNELS.md`'s EPIOS
count; `01_generators.md` §5.6 and §6, and two 2026-09-23 records, on
DJANGOH; `02_data_nucleon_deuteron.md` §5.2 on η; `05_chain.md` §1.1 and
CH-25 on T5), this page, and the report.

**Findings beyond the plan.**
* The shipped `ToyG1` neutron behaves as x^−(1+λ) at small x (λ = 0.176 /
  0.217 / 0.248 at Q² = 2 / 5 / 10): ∫g₁ⁿ, the Bjorken integral and the
  neutron and deuteron BC integrals do not exist for the default polarized
  model. A property of the model, not of the WW code.
* The tagged channel's default deuteron (Hulthén, β = 0.30 GeV, P_D = 0.045)
  has η = 0.0325, +17.2σ from the measured 0.0256(4); `02_data_nucleon_deuteron.md`
  §5.2's "the tagged channel rides on" holds for CD-Bonn and AV18 only (a
  dated note there now).
* The inclusive HepMC3 records do not balance four-momentum or charge: they
  leave out the (A−1) remnant by design (`docs/USAGE.md:162-166`, the
  remnant deliberately not written), while the conservation section of the
  convention (`docs/HEPMC3_CONVENTION.md:123-138`) states the identity with
  no exception. The missing four-momentum is 0.8196 (⁶Li) / 0.8450 (⁷Li) of
  the beam energy to 7.5e-14; the tagged and coherent records balance to
  1.0e-13 with exact charge.
* The polarized rate shift on the full default window is −6.8e-4: a
  weighted-vs-unweighted closure there cannot tell a working Mode W from
  w = 1. The closure runs in x ≥ 0.1, y ≥ 0.5, chosen from the exact tables
  before any pull there was computed.

## What this run did NOT establish

1. **Nothing that needs an LHAPDF set was run.** Row 4, the CT18NLO and
   NNPDFpol11 sub-rows of T1, and Bjorken on NNPDFpol11 read blocked
   (environment); those code paths never executed here. `REPORT.md` was made
   on this machine and records that in its header; **it is to be regenerated
   on the reference machine** (`make_report.py --all`, then `--check`), where
   `--compare` against this file will show row 4 and those sub-rows changing.
2. **PYTHIA 8.312, not the pinned 8.317.** T4 and T6 ran on 8.312 and say
   so; no number was taken on 8.317.
3. **The T1 references are a second-hand vendoring.** The tables come from
   the nnpdf-data 4.1.5 wheel (sha256 `6ed209ad…3490269a`), not from
   HEPData directly; the papers could not be read. Unverified: whether the
   published g₁ᵈ is D-state corrected, which R NMC used, whether NMC's
   normalisation is inside "sys", E143's beam energy, E143's arXiv field
   (hep-ex/9705012 against the survey's hep-ph/9802357), SMC's paper
   (recorded, never in the verdict). The verdict uses the diagonal χ²; the
   full covariance, built on NNPDF's CORR labels, is recorded beside it.
4. **Two reference inputs are the least verified.** η_d and Q_d are a
   transcription of the tree's own text (02 §5.2, `b1_nuclear.hpp`); the de
   Vries radii a survey transcription; α_s(M_Z) = 0.1180(9) in the Bjorken
   reference was written from memory (not in the PDG 2026 database the other
   constants came from).
5. **Not done:** a Deeps harness (a new plan row — the author's call; the
   file stays vendored, unwired); rows 6, 7, 10 still BLOCKED; T2's
   wave-function data (A(Q²), t₂₀, (e,e′p)) and every R / F_L dataset are
   not reachable offline or not in the wheel; T6's chain legs (abconv, npsim,
   EICrecon, Rivet) not run; the `spin_weight_<k>` names checked with zero
   entries only (the CLI never fills them).
6. **Shallow clone.** The two as-of tests skip here, naming the four missing
   commits; they pass in a full clone. CI's `actions/checkout` default is
   also shallow (a review verdict suggests `fetch-depth: 0`; not changed).
7. **Left for their owners** (outside what this run could edit):
   `docs/HEPMC3_CONVENTION.md`'s conservation clause (scope it to tagged and
   coherent records, or write the remnant — T6 pins the document's sha256,
   so either edit re-reads the clauses); the `tests/test_cluster.cpp`
   comment that should attribute 0.0256(4) to Rodning–Knutson as a
   measurement (C++); `t3_nmc_li6_over_d.py` builds `Tree()` before reading
   its table, so without the sets a tampered CSV is not hashed (it exits 2
   either way), and it hashes only the body after its provenance marker, so
   its header (source, licence) is pinned only by `REPORT.json`'s
   `data_sha256` through `--check`; `validation/reference/_manifest.json`
   still cites the old `tagged.py` lines (299–302) where the script and
   `validation/README.md` now say 849–852 — two outputs of one script that
   disagree, so a `polligen` re-dump **with unchanged inputs will rewrite
   that one string** (`could_not_call[0].reason`; no number). It was not
   re-dumped (no `polligen` here), not hand-edited (not a file this run
   touched; its sha256 `a0ff5177…` is quoted by
   `../run_2026-09-23/phase_A_port_gate.md`), and the script's string was not
   reverted (that would reopen review finding 25); `LI6_R2_POINT_FM2` = 6.0788 fm² sits 0.00048 fm² above
   its stated recipe (recorded, no action proposed).
8. **27 decisions pending; this run took none.**

## For the user

**Run on the reference machine:** `source env.sh && python3
validation/benchmarks/make_report.py --all`, then `--check`, with the
LHAPDF sets installed — the only way row 4 and the LHAPDF sub-rows get a
2026-09-26 number.

**Decide:** whether the inclusive HepMC3 record should carry the (A−1)
remnant or the convention should scope its conservation identity; whether
the Hulthén deuteron's η (+17.2σ) matters for the tagged channel's default
(the static observable says nothing about w(k) at high k); whether Deeps
becomes a plan row.

**Run or request:** a DJANGOH production with IEL31..IEL33 ≠ 0 **and**
WMIN ≤ M_p — the only thing that turns row 7 into a benchmark; the source is
public.

## Verification (2026-09-26)

Four adversarial verifiers (physics, honesty, code and regression lenses)
read the uncommitted tree on 2026-09-27 (UTC): **14 findings** (0 major,
6 minor, 8 nits), 12 distinct (the row-5 hash and the T4 / T6 exit code were
each found twice). A fix stage re-checked each one against the tree: **all
14 confirmed, none rejected; 13 fixed in code or text, 1 recorded** (the
manifest). No number, status, default, tolerance or reference value moved:
`make_report.py --compare` of the regenerated `REPORT.json` against the
pre-fix one reports 0 changes, and every harness prints the rows it printed
before (T6's tolerance text excepted, below).

| lens | finding | what was done |
|---|---|---|
| physics | `t2_deuteron_static`'s Hulthén closed forms take the S wave's β for the D wave, unchecked (true today: both 0.30 GeV) | `hulthen()` raises unless the two β agree; the transform test reads β from the channel and checks both `Wave.radial`s against the transformed pair; a new test feeds a D wave with β = 0.35 GeV and sees the refusal |
| physics | T6's four-momentum "1e-9" is relative to the summed beam energy (607 GeV for ⁶Li, 706.5 for ⁷Li): 200–700× looser on p_x, p_y than the `tests/test_hepmc.cpp` check it calls the document's | the rule was **not** changed after the verdict; the docstring, the row's tolerance text and the README row now say what is applied and how it differs. Measured on the same 800 events, not asserted: under doctest's per-component rule the tagged and coherent records balance to ≤ 4.3e-12 and every inclusive event misses by ≥ 0.846 — the same verdict |
| honesty, code | README convention 2 and this page said every harness hashes every byte; row 5 hashed only its ⁶Li row and row 4 the body after its marker, so review finding 6's skip reason was false | row 5 now also checks a whole-file `FILE_SHA256` (a test refuses one digit of the A/Z note and the SPDX tag); row 4 (not a file this run edited) is named as the exception in convention 2 and in item 7 above; review rows 1/10 and 6 corrected |
| honesty | the opt-in skip reasons said "costs more than 30 s"; T4, the only opt-in row, takes ≈ 10 s | `make_report.py`'s reason, comment and REPORT header line, and `conftest.py`'s skip reason and marker text, now say "> 30 s, or a generator-vs-generator sample" |
| honesty, code | T4 and T6 turned a missing or broken `lipolgen` into a BLOCKED row with exit 0, against README convention 6 | the `lipolgen` import now propagates: exit 2, `error` in the report. BLOCKED stays for a build tier that is off and a missing `pythia8` / `pyhepmc`; a subprocess test per harness; convention 6 records the change |
| honesty | the pytest baseline was labelled "at `f0a8f1e`" but was taken before `f0a8f1e` was committed | relabelled in the tallies below |
| honesty | the T3 PASS rows dropped the harness's own caveat that the radius is a tree input | the caveat is now in the README row, the plan's status note and the table above |
| code | `--compare` scaled by max(1, \|x\|, \|y\|), so a move below ≈ 1e-9 (the EST residual 1.138e-14 doubling; a 4e-16 p-value becoming 9e-10) read "0 changes" | the rule is purely relative (0 → 0 is no move; to or from nan / inf is one); a test replays both moves |
| code | `np.trapezoid` (NumPy ≥ 2.0 only) in `t2_deuteron_static` and `test_bench_t5_closures.py`, while `pyproject.toml` allows numpy ≥ 1.20 | the repo's `vmc_reconcile.py` shim in both |
| code | `tiers.py` gave a `_report.py`-only diff no "regenerate the report" note | the note now comes with any `_report.py` change; asserted |
| code | the no-scipy `chi2_sf` fallback in the three T1 harnesses took 1 − P: p = 0 for 4.17e-16, a domain error at χ² = 0 | the continued fraction for Q itself above χ²/2 = ndf/2 + 1 (Numerical Recipes' `gammq`), p = 1 at χ² = 0; with scipy hidden the three harnesses now print exactly what they print with it; a test against scipy (rel 1e-10) |
| regression | `validation/reference/_manifest.json` cites the old `tagged.py` lines while the script and `validation/README.md` cite the new ones | **recorded, not fixed** (item 7 above): no `polligen` here to re-dump, the file is outside this run's edits, and reverting the script's string would reopen review finding 25 |

The regression verifier found the rest of its lens clean: the ten harnesses
of 2026-09-23 give the same numbers and exit codes as at `f0a8f1e`, the
vendored data and `validation/reference` are unchanged, `coherent.hpp` is
comment-only at 792 lines, no test was removed or went from pass to fail or
skip, and the doc gates and both dump scripts behave as before.

## Tallies

Measured 2026-09-27 (UTC) on the working tree (base `f0a8f1e` + this run):
at the end of the integration stage, and every row but the C++ one again
after the verification fixes (the C++ sources and `build/` did not change).
Baseline = the same machine before this run's uncommitted work: for
`build/lipolgen_tests`, the build of `f0a8f1e`'s sources (not rebuilt
since); for pytest, **the Python tests of `302b57a`**, run at 21:29 UTC on
2026-09-26 against the same `build/python` module (built 21:16 UTC, not
rebuilt since) — i.e. *before* `f0a8f1e` was committed at 22:43 UTC
(corrected at verification: this said "at `f0a8f1e`"). Its one `np.trapz`
failure is fixed by `f0a8f1e` itself (`python/tests/test_b1_model.py`), so
at `f0a8f1e` it would read 42 failed (inferred, not measured; a verifier's
run on a non-git copy of `f0a8f1e`, where the as-of tests skip, gave
40 failed / 1086 passed / 161 skipped).

| gate | measured | baseline |
|---|---|---|
| `build/lipolgen_tests` (not rebuilt) | **416 cases / 403 passed / 13 failed**; 17 241 417 assertions, 0 failed; 228.5 s. All 13 throw "Info file not found for PDF set" (CT18NLO 9, NNPDFpol11_100 2, EPPS21nlo_CT18Anlo_Li6 2) | 416 / 403 / 13, the same 13 |
| `python -m pytest python/tests -q` | **39 failed / 1209 passed / 163 skipped**, 288.3 s, after the fix stage (39 / 1204 / 163, 286.9 s, at integration; the +5 passes are the fix stage's new tests: the Hulthén β refusal, T4 / T6 exit 2 without `lipolgen` (×2), the tiny-number `--compare`, the no-scipy χ² fallback). **All 39** throw "Info file not found for PDF set" (CT18NLO 27, NNPDFpol11_100 12): `test_sf_backend` 14, `test_knob_provenance` 9, `test_li6_unpol_band` 7, `test_mstw_sf` 6, `test_li7_rank2` 2, `test_b1_model` 1 — the same 39 as the baseline. Skips +5: the two as-of tests (shallow clone, naming the four commits), the two `bench_expensive` tests (naming `LIPOLGEN_BENCH=1`), the wheel-reproducibility test, Bjorken on `NNPDFpol11_100`, less one (row 4's table guard now runs without the sets). `test_real_document_external_citations` failed on the first run here (15 skip lines where it allowed 0 or 7: the 8 externals of `PYTHIA_BRIDGE.md` also print one, since 2026-09-05, whenever the PYTHIA source tree is absent) and was corrected to all-or-none over 7 + 8 | 43 failed / 1086 passed / 158 skipped, `302b57a`'s Python tests (above): 39 missing LHAPDF sets, 1 `np.trapz` (fixed by `f0a8f1e`), the two as-of tests (shallow clone), 1 external-citation test (no PYTHIA source tree) |
| `LIPOLGEN_BENCH=1 python -m pytest python/tests -q -k bench` | **146 passed / 7 skipped** (1258 deselected), 98.6 s, after the fix stage (141 / 7, 96.4 s, at integration; +5 new tests, as above). The two `bench_expensive` tests ran and passed (T4's full row; `make_report.py --check` of the committed report). Skips, each naming what is missing: row 4 × 3 (`EPPS21nlo_CT18Anlo_Li6`, `CT18NLO`), row 6 (the RFY / ADNDT document), row 10 (the BONuS supplement), the nnpdf-data wheel reproducibility test (`$LIPOLGEN_NNPDF_WHEEL` unset), Bjorken on `NNPDFpol11_100` | — |
| strict link gate | **1240 / 95 / 7 external / 0 broken / 24 allow-listed** (the 7 externals, and `PYTHIA_BRIDGE.md`'s 8, skipped by name: no PYTHIA source tree here); the same after the fix stage | 1240 / 95 / 7 / 0 / 24 |
| `--audit-ranges` | **209 range citations, 179 blocks, 0 REFUSED**; the same after the fix stage | 209 / 179 / 0 REFUSED |
| `--records-only` | **801 citations in 50 dated run records, 0 broken** (+1 record, this page, with its 2 citations); the same after the fix stage | 799 / 49 / 0 broken |
| `check_spdx_headers.py` | **125/125** (+9 harnesses, +4 tools under `validation/benchmarks/`); the same after the fix stage | 112/112 at `f0a8f1e` (a `git archive` of it) |
| the 17 harnesses | each run three times at integration: exit 0 and identical rows for 16; `t3_nmc_li6_over_d` exits 2 (its only configuration needs `EPPS21nlo_CT18Anlo_Li6`). Run once more after the fix stage: the same exit codes and rows (T6's tolerance text reworded, no number). With numpy hidden, or `lipolgen` shadowed by a package that raises ImportError: **16 exit 2**, BONuS (needs neither) exits 0 — T4 / T6 printed a BLOCKED row and exited 0 here until the verification fix | — |
| `make_report.py --all` (`LIPOLGEN_BENCH=1`) | 17 rows, 43.97 s wall; then `--check`: OK, 42.9 s. Regenerated after the fix stage (T6's tolerance text and the header's opt-in line changed): 17 rows, 5 PASS / 8 FAIL / 3 BLOCKED / 1 blocked (environment) / 0 error, 44.46 s; `--compare` against the pre-fix `REPORT.json`: 0 changes; `--check`: OK, 43.88 s | — |
