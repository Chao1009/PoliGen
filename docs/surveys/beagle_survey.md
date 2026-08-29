## 1. Repository layout, languages, build system

### Layout
`/home/cpeng/Projects/polli/BeAGLE/` contains only `docs/` (one PDF) and `source/` (a clone of `https://github.com/eic/BeAGLE.git`, branch `master`, last commit **2024-01-17**, `eba1b85`). Repo 84 MB (18 MB `.git`).

```
/home/cpeng/Projects/polli/BeAGLE/source/
  Makefile                       (5.2 kB)
  README.md
  user-eA3.0-5.f                 409 lines   main program + DT_USRHIS stub
  src/                           79,799 lines of .f
     dpmjet3.0-5F-new.f          29,983   DPMJET-3.0-5F + all BeAGLE nuclear physics
     phojet1.12-35c4.f           42,458   PHOJET (essentially dead code in BeAGLE)
     dpm_pythia.f                 2,411   PYTHIA6 <-> DPMJET interface  (the real BeAGLE core)
     dpm_gcf.f                    2,410   GCF external-event-file reader mode
     dpm_rapgap.f                   631
     + 18 small helper files (anear, azmass, dcalc, deutfix, fixestar, kickit,
       nbary, nrbinom, pfshift, reinit, rmkick, sclsubevt, shdmap, ...)
     Old/                        17 MB of dated backup copies (dead weight)
  PYTHIA-6.4.28/                 82,178 lines (vendored, locally patched)
  radgen-6.4.28/                  3,632 lines (POLRAD 2.0 derived — see §3/§4)
  PyQM/                           2,411 lines (nucdens.f, qweight.f, density.f)
  RAPGAP-3.302/                  include/ only — lib/ contains ONLY a README
  include/                       521 lines, 32 .inc files + LHAPDF C headers
  nuclear.bin  (43 MB)  dpmjet.dat (581 kB)  fitpar.dat  shadowmapAu/Ca.dat
  Examples/    AAAREADME, eD_...inp, muXe_pyqm_qhat_g1.inp, eAt1dfJn, S1ALL003
```

Total Fortran: **168,483 lines**. Everything is fixed-form FORTRAN 77 (`IMPLICIT DOUBLE PRECISION (A-H,O-Z)`, `SAVE`, computed `GOTO`, 467 `COMMON` blocks and 427 `GOTO`s in `dpmjet3.0-5F-new.f` alone). **Zero C, C++, or Python.** The only `.h` files are LHAPDF's own headers under `include/LHAPDF/`.

BeAGLE-authored (as opposed to vendored DPMJET/PHOJET/PYTHIA) code is roughly **7,800 lines** across `src/dpm_pythia.f`, `src/dpm_gcf.f`, `src/dpm_rapgap.f`, the 18 helpers and `user-eA3.0-5.f`, **plus** several thousand lines of heavily-modified DPMJET routines inside `dpmjet3.0-5F-new.f`.

### Build system
Hand-written GNU Makefile, `gfortran` only, no autotools/CMake (`/home/cpeng/Projects/polli/BeAGLE/source/Makefile`):

- L12–16: `CXXFLAGS = -c -g -Wall -m64 -fno-inline -fno-automatic -fno-align-commons` (the `-fno-automatic` and `-fno-align-commons` are load-bearing: the code relies on static locals and misaligned commons).
- Builds four archives (`libpythia6.a`, `libradgen.a`, `libpyqm.a`, `libdpmjet.a`) then links `BeAGLE`.

### External dependencies and hardcoded paths (all four link paths are hardcoded to BNL AFS/CERN paths that do not exist here)

| Dep | Makefile line | Status locally |
|---|---|---|
| **FLUKA** (`libflukahp.a` + `$FLUPRO/flukapro` includes) | L59 `FLUKA = /afs/rhic/eic/PACKAGES/fluka-64`, L60 `LIB2`, L82 `FLUINC` | **ABSENT.** `$FLUPRO` unset, `/afs` does not exist. Note L59 *overrides* `FLUKA=$(FLUPRO)` on L57 — a "temporary kludge" that must be edited. |
| **CERNLIB** (`mathlib kernlib packlib_noshift`) | L64 `-L/cern64/pro/lib` | **ABSENT.** |
| **LHAPDF 5** (`libLHAPDF.a`, and `$LHAPATH` at runtime) | L65 `-L/afs/rhic/eic/lib` | **ABSENT** (only headers vendored). LHAPDF **5**, not 6 — the code calls `STRUCTM`/`STRUCTA`/`SETLHAPARM`. |
| **RAPGAP 3.302** (`librapgap33.a libar4.a libbases.a`) | L52 `LIBRAP` | **ABSENT** — `RAPGAP-3.302/lib/` holds only `AAAREADME`. Must be built externally and `ar d librapgap33.a sfecfe.o` applied. |
| PYTHIA 6.4.28, DPMJET, PHOJET, PyQM, radgen | vendored | present |
| ROOT | — | not used at all |

FLUKA include files referenced with DEC-style parenthesised names in `src/dpmjet3.0-5F-new.f`: `(DIMPAR) (EVAFLG) (FHEAVY) (FRBKCM) (GENSTK) (NUCDAT) (PAREVT) (RESNUC)`. These ship only with a licensed FLUKA distribution.

### Buildable on Ubuntu 22.04 / gfortran 11?
**Not as-is; feasible with effort.** Assessment:

- Fortran itself: fixed-form F77 + a few GNU extensions (e.g. the `$` non-advancing format edit descriptor, 16 occurrences in `src/dpm_pythia.f`, format label 34 at L2283). gfortran 11 still accepts these (with warnings). The `-fno-automatic -fno-align-commons` flags are already present, which is what usually breaks legacy DPMJET on modern gfortran. So compilation of the vendored code is plausible.
- Blockers, in order of difficulty:
  1. **FLUKA.** You must register at fluka.cern / infn.it, accept the license (redistribution prohibited, source not open), get the 64-bit `libflukahp.a` + `flukapro` includes. The 2022 paper explicitly states *"Since FLUKA is not an open-source program, the BeAGLE event generator has no handle on changing the evaporation process."* Modern FLUKA releases (FLUKA4-x, and CERN's separate FLUKA.CERN fork) have a different library name/ABI than the 2011-era `fluka-64` this expects; matching an old-enough release is the single biggest practical risk. There is an unmerged branch `origin/copilot/beagle-fluka-integration` that suggests upstream is aware.
  2. **CERNLIB.** Ubuntu 22.04 has no packaged CERNLIB; you need a source build (`cernlib` 2005/2006 or the Debian `cernlib-*` backport). Only `mathlib`, `kernlib`, `packlib_noshift` symbols are used.
  3. **LHAPDF 5** (deprecated since 2014) with the **EPS09LO/EPS09NLO** sets installed under `$LHAPATH`.
  4. **RAPGAP 3.302** — optional in practice; you could stub out `LIBRAP` if you never use `MODEL RAPGAP`.
- Runtime also needs `$BEAGLESYS` pointing at the tree (`src/dpmjet3.0-5F-new.f:253-266`) and `nuclear.bin` (FLUKA evaporation data) in CWD.

**Verdict:** a determined week of work for someone with FLUKA/CERNLIB experience; a multi-week rabbit hole otherwise. Nothing about the Fortran itself is the hard part — the closed-source FLUKA link is.

---

## 2. How BeAGLE models the nucleus

### 2a. Nucleon spatial configuration — `DT_CONUCL` / `DT_COORDI`
`/home/cpeng/Projects/polli/BeAGLE/source/src/dpmjet3.0-5F-new.f:6817` (`DT_CONUCL`), `:6891` (`DT_COORDI`). Three regimes:

- **A = 1** (L6928): nucleon at origin.
- **A = 2** (L6931–6942): a hardcoded **sum of three Gaussians**, back-to-back:
  `DATA WD / 0.0, 0.178, 0.465, 1.0/` (L6920) and `DATA RD /2.09, 0.935, 0.697/` (L6921) fm. This is the original Shmakov deuteron parameterisation, not Hulthén/Reid.
- **A = 3, 4** (L6943–6968): isotropic Gaussian with `SIGMA = R/sqrt(2)`, then CM-shifted.
- **A ≥ 5** (L6970 onwards): rejection sampling from `DT_DENSIT`, with a hard-core exclusion `R2MIN = 0.16 fm²` (i.e. ≥ 0.4 fm nucleon separation), optional **deformed 3D Woods–Saxon** with β₂/β₄ (`Y20`, `Y40`, L7010–7025) and random nuclear orientation. Then CM-shifted.

**Consequence for ⁶Li/⁷Li: they take the A ≥ 5 branch — a plain spherical (or β-deformed) Woods–Saxon ball. No α+d, no α+t, no cluster correlations whatsoever.**

### 2b. Density — `DT_DENSIT` and PyQM `nucdens.f`
`src/dpmjet3.0-5F-new.f:7098` — `DT_DENSIT` is now just a table lookup into `density_table` filled by PyQM. The original DPMJET shell-model (A ≤ 18) / Woods–Saxon (A > 18) code is **commented out** (L7137–7148).

The table is built in `/home/cpeng/Projects/polli/BeAGLE/source/PyQM/density.f:42` `GenNucDens`, with **`idist = 1` hardwired at L57** — i.e. always Woods–Saxon, never the Reid soft-core deuteron density. `PyQM/nucdens.f:395` gives the 3-parameter form ρ ∝ (1 + c₀ r²/R²)/(1 + exp((r−R)/a₀)).

**Critical for Li:** `PyQM/nucdens.f:146-147` sets `Rws(3) = aws(3) = cws(3) = 0` for **Z = 3 (lithium)**. The code therefore falls through to the Bialas generic parameterisation at L376:
```
RR = (0.978 + 0.0206*A**(1/3)) * A**(1/3)   ;  a0 = 0.523 ;  c0 = 0
```
→ for A = 6: R ≈ 2.44 fm, a = 0.523 fm. Purely generic. (Real measured Li has Rws ≈ 1.77–2.0 with a large surface diffuseness and a two-body halo-ish tail; and of course no cluster substructure.)

`GLAUB-3D` card (`src/dpmjet3.0-5F-new.f:2501`) lets you override R₀, a, w, β₂, β₄ from the input file — the **only** existing user handle on Li geometry, and it cannot express α+d.

`DT_RNCLUS` (`:7158`): `RADNUC(1..8) = 0`, so A = 6, 7 get `R = 1.12 × A^(1/3)` = 2.04 / 2.14 fm for Glauber purposes.

### 2c. Fermi motion — `DT_FER4M`, `DT_KFERMI`, `DT_DFERMI`
`DT_FER4M` at `src/dpmjet3.0-5F-new.f:5533` dispatches (L5583–5588):

- **A = 2, 3, 4** → `DT_KFERMI` (`:17219`): Ciofi degli Atti & Simula PRC 53 (1996) 1689 parameterisation, three-Gaussian/pole form. Separate hardcoded parameter sets for **d** (A0=157.4, B0=1.24, C0=18.3, …), **³He** (L17311), **⁴He** (L17321). Also implements the SRC/high-k-tail sampling variants via `KRANGE` = `IFMDIST`.
- **A > 4** → `DT_DFERMI` (`:17101`): coarse A-bins (>208, 56–208, 40–56, 16–40, 12–16, 4–12).
- `IFMDIST = -1` → `DT_DFERMIO` (`:17069`), the ancient "largest of three random numbers" DPMJET hack.

> **Bug worth knowing about.** `DT_DFERMI`'s branch conditions use `.OR.` where `.AND.` is clearly intended:
> `src/dpmjet3.0-5F-new.f:17125`: `ELSE IF( (ANUCLEUS .LE. 208) .OR. (ANUCLEUS .GT. 56) ) THEN`
> (same pattern at 17132, 17139, 17146, 17153). For **any** A ≤ 208 the second branch is unconditionally true. So **⁶Li and ⁷Li are given the Fe–Pb (56 < A ≤ 208) momentum distribution** — A0=1.80, B0=4.77, D0=25.5, F0=40.3 — not the "4 < A ≤ 12" set the author intended. This is a real, silent physics error for any A between 5 and 56.

Momentum is then scaled by `FERMOD` (the `FERMI` card's 2nd number, D = 0.62) and given an isotropic direction — i.e. **no spin–momentum correlation and no D-state/tensor structure**.

### 2d. Off-shellness / energy conservation
This is BeAGLE's acknowledged weak point. PYTHIA6 generates the γ*N sub-event on an **on-shell** nucleon; Fermi motion is grafted on afterwards. Options via `FERMI` what(1) = `IFERPY` (`src/dpmjet3.0-5F-new.f:719-745`):

- `-1` no Fermi motion; `1` DPMJET has p_F but the PYTHIA sub-event does not; `2` p_F applied to the PYTHIA skeleton after the fact (`PFSHIFT`, `src/pfshift.f`, 277 lines, called at `src/dpm_pythia.f:1218,1232,1235,1239`); `3` **A = 2, 3 or SRC kinematics** (A > 3 non-SRC reverts to 2).
- `FERMI` what(4) = `IFMPOST`: `1` = ad-hoc IRF energy correction via **`DEUTFIX`** (`src/deutfix.f`, 301 lines; called at `src/dpm_pythia.f:1284`, guarded by `STOP "FATAL: CAN ONLY POST-FIX DEUTERON KINEMATICS"`); `2` = **light-front kinematics** via `DT_SPECTRALFUNC` (`src/dpmjet3.0-5F-new.f:5204`), which implements the Strikman–Weiss LF α_spec/α_active prescription with `Md` **hardcoded** at 1.87561 GeV and a `STOP "FATAL: CAN ONLY DO DEUTERON KINEMATICS"` guard at `src/dpm_pythia.f:1176`.
- Excess energy otherwise goes into the residual-nucleus excitation E*; when that goes negative the event used to be re-rolled — the 2023 **`ESTARFIX`** card (`:2615`, plus `src/fixestar.f`, 234 lines, and `src/getestar.f`) now resamples E* from a truncated Gaussian/flat/histogram instead.

Note the 2022 paper says `IFERPY=3` and the post-processing flag are "not yet implemented"; **the code is newer than the paper** — both exist now, as do `ESTARFIX` and `CHARMPOT` (mid-2023).

### 2e. Short-range correlations
`DT_PICKSRC` at `src/dpmjet3.0-5F-new.f:5288`, called from `src/dpm_pythia.f:1112`. Comment at L5311: *"for now only A > 12 assign SRC pairs and bring them half way"*; L5334 hardcodes a **20 % SRC probability** (`IF( D00 .LE. 0.2D0 )`). Selected via `FERMI` what(3) = 4 or 5. Returns `SRC_PARTNER_INDEX`, which is then treated like the deuteron spectator. **⁶Li/⁷Li (A < 12) get no SRC pairs at all.**

### 2f. The eD mode with spectator tagging — where it lives
This is the single most reusable piece for your purposes:

- `src/dpm_pythia.f:1180-1196` — finds the spectator (`ISTHKK == 14`) in A = 2, Z = 1, with `STOP 'FATAL: Could not find a spectator nucleon in deuteron.'`
- `src/dpm_pythia.f:1300-1305` — *"If it is a spectator nucleon from e+D, leave it alone. If it is a spectator nucleon from e+A, put in the event record, but don't ptkick it."* The deuteron spectator is deliberately **not** kicked or recoil-corrected, so its momentum is the raw wave-function momentum. This is exactly the design choice you would want to preserve for α/d/t spectators from Li.
- `src/deutfix.f` (301 lines) and `src/pfshift.f` (277 lines) are the 4-momentum repair machinery.
- `src/dpmjet3.0-5F-new.f:5175-5182` — for A = 2 the two nucleons are forced exactly back-to-back in the rest frame.
- **No FSI at all in eD** (paper, §II).

### 2g. INC and evaporation/fragmentation
- **INC**: `DT_FOZOCA` at `src/dpmjet3.0-5F-new.f:9391` ("FOrmation ZOne suppressed intranuclear CAscade"). Controlled by `TAUFOR` (`:783`): τ₀ in fm/c, generations followed (D = 25), p_T-dependent vs constant formation zone. **Explicitly disabled for A = 3** (`:800-810`, "Intranuclear cascade disallowed for A=3", `KTAUGE` forced to 0). It is *not* disabled for A = 6, 7 — but a formation-zone cascade in a 6-nucleon α+d system is physically dubious.
- **Residual nucleus / E***: `DT_RESNCL` (`:11077`) accumulates recoil from particles leaving the nuclear potential; `DT_NCLPOT` (`:10897`) sets the potential.
- **Evaporation**: `DT_FICONF` (`:12049`) computes A, Z of the remnant (`IBRES`, `ICRES`, L12659-12667), fills FLUKA's `FKFINU` common, and calls **`CALL EVEVAP(WE)`** at `:12706` — FLUKA's evaporation/fission/Fermi-break-up/γ-deexcitation. `DT_EVA2HE` (`:12737`) copies results back.
- **Can BeAGLE produce intact d / t / ³He / α spectators?** **Yes** — but only as FLUKA *evaporation / Fermi-break-up products*, never as pre-formed clusters. `src/dpmjet3.0-5F-new.f:12849-12862` is the "heavy fragments" loop reading FLUKA's `FHEAVY` common (`NPHEAV`, `KHEAVY(i)`, `IBHEAV`, `ICHEAV`, `PHEAVY`, `TKHEAV`) and emitting them as `IDHKK = 80000`, `ISTHKK = -1`, with A in `IDRES` and Z in `IDXRES`. FLUKA's `FHEAVY` species set is d, t, ³He, α (+ break-up fragments via `(FRBKCM)`). Fermi break-up is enabled by the `EVAP` card's what(2) `if` sub-flag (`:960-1001`, sets `LFRMBK`).
- **The final residual nucleus** is written as `ISTHKK = 1001` (`:12866-12870`).
- **Coherent / intact-recoil: NOT POSSIBLE.** A nucleon is always struck and removed, so A_remnant ≤ A−1 by construction. The paper says so directly: *"Due to the structure of the BeAGLE generator coherent diffraction is currently not included… for any nuclear beam, the target nucleus will break up or at least be excited."* For diffraction it simply assumes σ_eA = A·σ_eN.

---

## 3. Input card format, output format, polarization parameters

### Input card format
Fixed format `A10, 6E10.0, A8` (codeword, `WHAT(1..6)`, `SDUM`); `*` = comment; the `PHOINPUT…ENDINPUT` block is format-free. Parser: `DT_INIT` at `src/dpmjet3.0-5F-new.f:36`; the full codeword table is at **`src/dpmjet3.0-5F-new.f:220-233`** — **69 cards**:

```
TITLE PROJPAR TARPAR ENERGY MOMENTUM CMENERGY EMULSION FERMI TAUFOR PAULI
COULOMB HADRIN EVAP EMCCHECK MODEL PHOINPUT GLAUBERI FLUCTUAT CENTRAL RECOMBIN
COMBIJET XCUTS INTPT CRONINPT SEADISTR SEASU3 DIQUARKS RESONANC DIFFRACT SINGLECH
NOFRAGME HADRONIZE POPCORN PARDECAY BEAM LUND-MSTU LUND-MSTJ LUND-MDCY LUND-PARJ
LUND-PARU OUTLEVEL FRAME L-TAG L-ETAG ECMS-CUT VDM-PAR1 HISTOGRAM XS-TABLE
GLAUB-PAR GLAUB-INI PY-INPUT INPFILE OUTFILE VDM-PAR2 XS-QELPRO RNDMINIT
LEPTO-CUT LEPTO-LST LEPTO-PARL FSEED OUTPUT PYQM-CTRL USERSET PYVECTORS
GLAUB-3D CHARMPOT ESTARFIX START STOP
```

The physics-relevant ones (documented in-line at the cited lines):

| Card | Line | Parameters |
|---|---|---|
| `PROJPAR` | 446 | `SDUM` = `ELECTRON` / `MUON+` |
| `TARPAR` | 499 | A, Z, n/p handling (0=seq n then p, 1=all en, 2=all ep, 3=random mix), target-mass method (`AZMASS` IMETH), hypernucleus handling |
| `MOMENTUM` | 610 | p_lepton, p_ion, **crossing angle (mrad)**, c.a. orientation |
| `FERMI` | 719 | `IFERPY` (−1/1/2/3), `FERMOD` scale (D=0.62), `IFMDIST` (−1,0..5,11..14 — n(k) variants and SRC), `IFMPOST` (0/1/2 = none / ad-hoc / light-front), extra unbinding energy (MeV), light-cone potential flag |
| `TAUFOR` | 783 | τ₀ (fm/c), # generations (D=25), p_T-dep vs const formation zone, cascade-nucleus selection mode |
| `EVAP` | 885 | level-density option i1+10·i2+100·i3+10000·i4; what(2) = `ig + 10·if` → deexcitation γ on/off, **Fermi break-up on/off** |
| `MODEL` | 1032 | `PYTHIA` (D) / `RAPGAP` / **`GCF-FT`** (read external event file) |
| `OUTLEVEL` | 1860 | 6 verbosity flags (last = PyQM) |
| `L-TAG` | 1900 | y_min, y_max, Q²_min, Q²_max, θ_min, θ_max (lab, rad) |
| `L-ETAG` | ~1930 | E'_min, E_γ min/max |
| `PY-INPUT` | — | name of the PYTHIA6 control file (**≤ 8 characters**) |
| `FSEED` | 2304 | `INSEED, IJKLIN, ISEED1, ISEED2` |
| `OUTPUT` | 2318 | output unit (D=89), format: −1 old / +1 DPMJET nuclei / +2 GEANT PDG ion codes (D) |
| `PYQM-CTRL` | 2338 | recoil fraction to E*, p_T model `iPtF` 0–3, gluon-emission mode 0–3, SupFactor, heavy-quark SW |
| `USERSET` | 2375 | userset # 0–17 (defines `USER1/2/3`), **E*_max for evaporation (D = 9.0 GeV)** |
| `PYVECTORS` | 2485 | 0 = all VMs, 1–4 = ρ/ω/φ/J/ψ only |
| `GLAUB-3D` | 2501 | R₀, a, w, **β₂, β₄**, orientation (0 = random) |
| `CHARMPOT` | 2574 | charm-meson nuclear potential on/off, ratio |
| `ESTARFIX` | 2615 | E* fix distribution type (−1/0,1/2/3), correct 4-mom?, centre/min, sigma/max |
| `START` | ~2660 | # events, Glauber init mode 0–3 |

The **PYTHIA control file** (referenced by `PY-INPUT`; examples `Examples/eAt1dfJn`, `Examples/S1ALL003`) is a separate free-format file: line 1 output filename, line 2 `xmin,xmax`, line 3 F2 model + R parameterisation, line 4 radiative-correction switch, line 5 Pythia model, line 6 `A-Tar, Z-Tar`, line 7 `genShd` (1/2/3), line 8 nPDF order+error set (`x*100+y`), line 9 quenching on/off, line 10 `QHAT`, then `SHDFAC`, `FSEED`, and arbitrary `MSTP/PARP/PARJ/MSTJ/MSTU/CKIN/MSUB/MDCY` assignments.

### Output format
ASCII, written from `src/dpm_pythia.f` (`DT_PYOUTEP`). Header `" BEAGLE EVENT FILE "` at `src/dpm_pythia.f:1878`; column names in format statement **31** at `src/dpm_pythia.f:1903-1911` (current, `OLDOUT=.FALSE.`):

**Per-event line (58 fields):**
```
I, ievent, genevent, lepton, Atarg, Ztarg, pzlep, pztarg, pznucl, crang, crori,
subprocess, nucleon, targetparton, xtargparton, beamparton, xbeamparton,
thetabeamprtn, truey, trueQ2, truex, trueW2, trueNu, leptonphi,
s_hat, t_hat, u_hat, pt2_hat, Q2_hat, F2, F1, R, sigma_rad, SigRadCor, EBrems,
photonflux, b, Phib, Thickness, ThickScl, Ncollt, Ncolli, Nwound, Nwdch,
Nnevap, Npevap, Aremn, NINC, NINCch, d1st, davg, pxf, pyf, pzf, Eexc, RAevt,
User1, User2, User3, nrTracks
```
(format 30 at `:1892-1901` is the older variant without `crang, crori, Nwound, Nwdch`.)

**Per-track line (18 fields)**, header at `src/dpm_pythia.f:1915-1918`, written by format 34 at `src/dpm_pythia.f:2283`:
```
I  ISTHKK(I)  IDHKK(I)  JMOHKK(2,I)  JMOHKK(1,I)  JDAHKK(1,I)  JDAHKK(2,I)
PHKK(1..3,I)  PHKK(4,I)  PHKK(5,I)  VHKK(1..3,I)  IDRES(I)  IDXRES(I)  NOBAM(I)
```
`IDRES` = fragment A, `IDXRES` = fragment Z. With `OUTPUT` what(2) = 2 (default), nuclei are re-coded to PDG ion IDs `1000000000 + 10000*Z + 10*A` (`src/dpm_pythia.f:2246`). Status remapping for output (`:2262-2270`): evaporation products `ISTHKK=-1` → `KSOUT=1, NOBAM=3`; residual nucleus `ISTHKK=1001` → `KSOUT=1, NOBAM=4`. Terminated by `=============== Event finished ===============`. This is the format `eic-smear` parses.

### Polarization / spin — grep results

**In BeAGLE's own code (`src/*.f`, `user-eA3.0-5.f`, `include/*.inc`): nothing.** Grepping `polariz|helicit|spin|tensor pol|g1|A1|b1|asymmetr` over `src/dpm_pythia.f`, `src/dpm_gcf.f`, `user-eA3.0-5.f` and all `include/*.inc` returns **zero hits**. The only hits in `src/dpmjet3.0-5F-new.f` are two irrelevant comments (L2165 longitudinal photon flag; L26590 polarized *lepton decay* i.e. μ→eνν). PHOJET's hits are all γ helicity/VMD, unrelated to beam/target spin.

**BUT — `radgen-6.4.28/` is a POLRAD 2.0 derivative and contains a complete polarized-DIS radiative-correction skeleton.** This is the single most interesting discovery for your project:

- `/home/cpeng/Projects/polli/BeAGLE/source/include/polcom.inc` — `common/pol/ as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks`
- `/home/cpeng/Projects/polli/BeAGLE/source/include/tailcom.inc` — `common/tail/ un,pl,pn,qn,ita,isf1,isf2,isf3,ire` (`pl` = lepton polarization, `pn` = nucleon polarization, `qn` = tensor)
- `/home/cpeng/Projects/polli/BeAGLE/source/include/radgenkeys.inc:4` — `common/radgenkeys/ plrun,pnrun,ixytest,kill_elas_res`
- **`radgen-6.4.28/radgen.f:1410-1435`** — the live branch:
  ```
  if ((plrun.ne.0).and.(pnrun.ne.0)) then
     call mkasym (t,aks,itara,itarz,asym1,asym2)
     ga=sqrt(t)/anu
     g1=(asym1+ga*asym2)/(1+ga**2)*f1*fdilut(t,aks,itara)     ! line 1417
     g2=(asym2-ga*asym1)/ga/(1+ga**2)*f1*fdilut(t,aks,itara)
     b1=0.d0 ... b4=0.d0
  elseif((plrun.eq.0).and.(abs(pnrun).eq.2)) then
     call mkasym(t,aks,itara,itarz,asym1,asym2)
     b1=-3./2.*f1*asym1                                        ! line 1425  ← TENSOR
     b2=2.d0*aks*b1
  ```
  i.e. **g₁, g₂ and the tensor structure functions b₁–b₄ are already plumbed through the cross-section sum** (`radgen.f:790-794`: `if(isf.eq.3.or.isf.eq.4) ppol=pl*pn` for g₁/g₂; `if(isf.ge.5) ppol=qn/6` for b₁–b₄). There is also elastic-tail tensor code at `radgen.f:1514` (`b1=2.*tau**2*fm**2`).
- **It is deliberately switched off**: `radgen-6.4.28/radgen_init.f:80-83` — *"Using radgen with Pythia the target and beam polarisation can be hardcoded to zero"* → `plrun = 0.` / `pnrun = 0.`
- **`mkasym` is an empty stub**: `radgen-6.4.28/pythia_radgen_extras.f:254-268` — *"Routine is empty because Pythia is unpolarised, but radgen expects it"*, returns `dA1 = dA2 = 0`.
- `fdilut` (dilution factor) only special-cases ³He.
- radgen's nuclear list (`radgen_init.f:55-78`) is n, p, D, ³He, ⁴He, ¹⁴N, ²⁰Ne, ⁴⁰Ca, ⁸⁴Kr, ¹³¹Xe, ¹⁹⁷Au, ²³⁸U — **no A = 6 or 7**.

So: BeAGLE ships a dormant, ~50-line-away hook for A₁/A₂ → g₁/g₂ and a tensor b₁ branch, but only inside the *radiative-correction* engine, not in the event-generation matrix element. It would give you weighted radiative corrections, not spin-dependent event generation.

---

## 4. Where the DIS cross section is computed; nPDF/shadowing; spin-injection points

**Not DPMJET, not PYSIGH.** For `MODEL PYTHIA` (the only mode anyone uses), the γ*N interaction is generated by **PYTHIA 6.4.28**, driven from `src/dpm_pythia.f`. The DIS/photoproduction subprocess machinery is PYTHIA's `PYINIT`/`PYEVNT` → `PYSIGH` path with the low-Q² VMD/GVMD extensions (`MSTP(14)`, `MSTP(17)`, `PARP(161-166)`). `radgen` is bolted on for QED radiative corrections. `PHOJET` is present but unreachable (`user-eA3.0-5.f:52-53` refuses hadronic/photoproduction collisions).

**Nuclear PDFs (EPS09):** applied as a *cross-section-level reweight*, not inside the matrix element. `src/dpm_pythia.f:847-899`:
```
CALL STRUCTM(XX,Q2,...,USEAP,...)                 ! free-proton u-sea
Nprm='EPS09LO,'//chpset ;  CALL SETLHAPARM(Nprm)   ! L875 (or EPS09NLO L881)
ATEMP = ANEAR(INUMOD)                              ! L878
CALL STRUCTA(XX,Q2,ATEMP,...,USEAA,...)            ! L885 nuclear
RAVAL = USEAA/USEAP                                ! L893 — u-sea R_A for ALL processes
IF (RAVAL.GE.1.D0) RAVAL=1.D0                      ! only shadowing kept, no antishadowing/EMC
```
Applied only for x ≤ 0.1 (L851). Flavour-blind (comment L890: *"Proper flavor dependence will take a bit of care…"*).

`ANEAR` (`src/anear.f`) snaps A to the nearest EPS09 grid value; the valid list is `4, 6, 9, 12, 16, 27, 40, 56, 64, 108, 115, 117, 184, 195, 197, 208, 238` — **A = 6 is exact for ⁶Li, and ⁷Li (A = 7) snaps to 6** (edges 4.9–7.3 → 6). So EPS09 is usable for Li out of the box.

**Multi-nucleon shadowing:** `genShd` 1/2/3. `src/dpmjet3.0-5F-new.f:8020`: `SIGEFF = SHDFAC*SHDMAP(RAVAL)` — the R_A value is mapped to an effective *dipole* cross section via a quintic spline (`src/shdmap.f`, knots in `shadowmapAu.dat`/`shadowmapCa.dat`), giving a transverse radius D = sqrt(σ_eff/π) inside which additional nucleons participate. `SHDFAC` is a user tuning factor (e.g. 1.32 in the example card).

**Where spin-dependence could be injected (in rough order of increasing effort):**
1. **Event weight** — compute A_∥/A₁ externally per event from (x, Q², struck-nucleon flavour, nucleon polarization) and store in `USER1/2/3` (`USERSET` card, `src/dpmjet3.0-5F-new.f:2375`) or in the free `RAevt` slot. Zero physics inside the generator; you just reweight. Cheapest, and probably right for a first study.
2. **`RAVAL` analogue** — mirror the `src/dpm_pythia.f:847-899` block with polarized nPDFs (there is no EPS09-equivalent for Δf^A; you would need e.g. a model of Δq_Li from α+d). Same flavour-blindness caveat.
3. **`mkasym` in `radgen-6.4.28/pythia_radgen_extras.f:257`** plus un-hardcoding `plrun/pnrun` in `radgen_init.f:82-83` — gets you g₁, g₂ *and* the b₁ tensor branch (`radgen.f:1425`) for the radiative-correction cross sections. Structural, ~100 lines, but it only affects radiative weights.
4. **PYTHIA6 `PYSIGH`/`PYSTFU`** — would require polarized parton densities and helicity-dependent hard matrix elements throughout PYTHIA 6. Effectively rewriting PYTHIA6. **Not viable.**
5. **cos 2φ double-helicity-flip (Δ(x) / the b-type gluon transversity term)** — there is no azimuthal-modulation infrastructure anywhere in BeAGLE; `leptonphi` is an output column only. Would be entirely new.

---

## 5. Assessment

### License
**There is no LICENSE, COPYING or COPYRIGHT file anywhere in the tree.** No SPDX headers, no license statement in `README.md`. The only copyright notices are PYTHIA's (`PYTHIA-6.4.28/pydata.f:51`, `pylogo.f:66` — "Copyright Torbjörn Sjöstrand"). The GitHub repo `eic/BeAGLE` likewise carries no license as of this snapshot. Practically:
- Redistributing a fork is legally murky (all-rights-reserved by default under Berne).
- **FLUKA is definitively closed-source and non-redistributable**; you can link against it but not ship it, and you cannot modify the evaporation model.
- DPMJET-3 and PHOJET have their own restrictive terms; PYTHIA 6 is GPL-2-ish (Sjöstrand's terms); LHAPDF 5 is GPL.
- If you plan to release a generator publicly, inheriting this dependency stack is a real problem, and the FLUKA dependency alone forecloses a clean open-source release.

### Code quality / maintainability
- **Fixed-form FORTRAN 77 throughout.** 467 `COMMON` blocks and 427 `GOTO`s in one 30 k-line file. Global mutable state everywhere; `SAVE` + `-fno-automatic` reliance.
- No unit tests, no CI, no `configure`, hardcoded AFS paths, 17 MB of dated backup copies in `src/Old/`.
- Commented-out code is pervasive (e.g. the entire event-writing block in `user-eA3.0-5.f:198-397` is dead).
- Physics comments are actually **good** — the card documentation blocks in `DT_INIT` are thorough, and Mark Baker's `MDB` annotations date and explain most changes. This is much better than typical legacy Fortran.
- The `.OR.`/`.AND.` bug at `src/dpmjet3.0-5F-new.f:17125-17153` is evidence that light/medium-A configurations are not exercised.
- Development is essentially dormant: last commit Jan 2024; ~7 branches, one from a Copilot FLUKA-integration experiment.

### (a) Fork BeAGLE and add ⁶Li/⁷Li cluster + spin
**Hard, and the payoff is bounded.** What you would have to build:

| Need | Status in BeAGLE | Work |
|---|---|---|
| α+d / α+t cluster configuration | absent; A ≥ 5 → spherical WS | new `DT_COORDI` branch + cluster wave function. Moderate (200–400 lines), but touches Glauber, `dcalc`, INC positions. |
| Cluster-resolved Fermi motion (2-body α–d relative momentum, not single-nucleon n(k)) | absent; and the A>4 dispatcher is buggy | new `DT_KFERMI` branch. Moderate. |
| Intact α/d/t spectator emitted as a *pre-formed* cluster | absent — clusters only come out of FLUKA evaporation | **Structural.** Requires a new branch parallel to the eD spectator path (`src/dpm_pythia.f:1180-1305`), plus telling FLUKA the remnant is (A−1−A_cluster) rather than (A−1). Feasible but invasive. |
| Coherent / intact-nucleus recoil | impossible by construction | would require bypassing the whole struck-nucleon architecture. Effectively a different generator. |
| g₁ / A₁ | absent (except the radgen stub) | reweighting only, realistically |
| Tensor b₁ | `radgen.f:1425` skeleton only | reweighting only |
| cos 2φ double-helicity-flip | absent entirely | new from scratch |
| C++/Python interface | none; ASCII text out | wrap with `f2py`/ISO_C_BINDING, or parse the ASCII (what eic-smear does) |

Plus you inherit the FLUKA + CERNLIB + LHAPDF5 build problem permanently, an unlicensed codebase, and a "Fermi motion doesn't enter the DIS kinematics" limitation the authors themselves call unreasonable for few-nucleon systems.

### (b) New C++ generator on PYTHIA 8, mining BeAGLE for ideas/numbers
**This is the better path, and it is not close.** Reasons:
- Everything you actually need from BeAGLE is *physics parameters and algorithms*, not code: the deuteron 3-Gaussian configuration (`:6920-6921`), the Ciofi–Simula n(k) coefficient tables (`DT_KFERMI` `:17219-17330`), the light-front spectator prescription (`DT_SPECTRALFUNC` `:5204`), the "don't recoil-correct the light-nucleus spectator" rule (`src/dpm_pythia.f:1300-1305`), the formation-time law τ = τ₀(E/m)m²/(m²+p_T²), the EPS09 reweight recipe (`src/dpm_pythia.f:847-899`), the `d`/thickness geometry observables (`src/dcalc.f`), and the output column set (which you should keep, for eic-smear compatibility).
- PYTHIA 8 gives you a real C++/Python (`pythia8mc`) API, `Pythia::forceHadronLevel`, user hooks, HepMC3 output, and a maintained DIS treatment.
- For ⁶Li/⁷Li you are dealing with *two-body* breakup dominated by α+d / α+t, not a 200-nucleon cascade. INC and FLUKA evaporation are largely overkill; a Fermi-break-up/two-body-decay model for the α+d and α+t channels plus a small excited-⁶Li/⁷Li level scheme is tractable and *more* accurate than FLUKA's statistical evaporation at these A.
- The spin content (g₁, b₁, cos 2φ) has to be written from scratch either way — there is nothing to reuse. Writing it in C++ against a clean spin-dependent structure-function interface is far easier than grafting it onto 467 COMMON blocks.
- You avoid the FLUKA license entirely, which means you can actually release the thing.
- **Hybrid option worth noting:** BeAGLE's `MODEL GCF-FT` mode (`src/dpm_gcf.f`, `DT_GCFINITQE` L8, `DT_GCFEVNTQE` reader at L370-436) already **reads an external ASCII event file in BeAGLE's own format** and then runs INC + FLUKA evaporation + output on it. That is a ready-made "use BeAGLE as a nuclear-breakup back end" hook. Currently hardcoded to `nrTracks.NE.9` (QE-specific, `src/dpm_gcf.f:384`) and a fixed track layout, but relaxing that is a ~100-line change. If you ever want FLUKA-quality de-excitation of a ⁶Li* remnant produced by your own C++ front end, this is how you would do it without forking the physics.

**Recommendation:** build new in C++ on PYTHIA 8; treat BeAGLE as a *reference implementation and parameter source*, and keep its ASCII output format so existing EIC far-forward analysis tooling (eic-smear, the ZDC/B0/Roman-Pot acceptance studies in the paper) works unchanged. Optionally keep the `GCF-FT` path as an escape hatch for FLUKA de-excitation.

---

## 6. Contents of `/home/cpeng/Projects/polli/BeAGLE/docs`

Exactly one file: **`2204.11998v1.pdf`** (1.35 MB, 20 pages). No manual, no user guide.

**Paper:** *BeAGLE: Benchmark eA Generator for LEptoproduction in high energy lepton-nucleus collisions*, W. Chang, E.-C. Aschenauer, M. D. Baker, A. Jentsch, J.-H. Lee, Z. Tu, Z. Yin, L. Zheng. arXiv:2204.11998v1 [physics.comp-ph], dated April 2022. **No journal-ref, no DOI, no license or repository statement** — the only pointer is the user's guide at `https://eic.github.io/software/beagle.html` (and the BNL wiki `https://wiki.bnl.gov/eic/index.php/BeAGLE`). Version described: 1.01.03. The code in this tree is ~1.5 years newer than the paper.

### Key statements — light nuclei
- *"DPMJet is not designed for light nuclei, so substantial changes had to be made for the case when the nucleus is a deuteron."* The light-nucleus special-casing is **deuteron-only**.
- Fermi momentum: Ciofi degli Atti–Simula, with **A = 2, 3, 4** using one functional form and **A > 4** another. For A > 4: *"the parametrizations are based on a few typical nuclei A_typical, e.g., Carbon-12, Oxygen-16, Calcium-40, Iron-56, Lead-208… Any nucleus between them, A_select, will use one of the nearest typical nuclei."* → ⁶Li/⁷Li are meant to get the ¹²C distribution (and, per the bug above, actually get the Fe–Pb one).
- **Lithium is never mentioned. "Cluster", "α+d", "α+t" never appear.** "Alpha" appears once, describing evaporation products.
- Only n₀(k) is implemented; **SRC tail n₁(k) is absent for A > 2**.
- GCF is listed as *future work*, not current capability.
- Deuteron-only extras: Strikman–Weiss light-front extension; ad-hoc energy-conservation post-fix.

### Key statements — spectator tagging
- Motivation: *"tagging of the spectator nucleon in eD scattering to allow for the extraction of free nucleon structure, as well as to study Short-Range Correlations in the deuteron."*
- The design rule you should copy: *"in the case of deuteron (or light nuclei in general) this correction will not be reasonable because there is only one spectator nucleon… The correction would artificially distort the spectator momentum distribution. Therefore, for the deuteron, we leave the spectator unmodified… By using this approach, the spectator tagging and related physics topics can be studied with the genuine information from the wave function. **No Final-State Interactions (FSI) are present in the BeAGLE generator for lepton-deuteron collisions.**"*
- Far-forward fragment kinematics (ePb 18×110, Fig. 10): *"The nuclei distributed within 7 < η < 10 are light nuclei, e.g., deuterons and alpha particles. The large remnant nuclei are distributed within 10 < η < 15."* Protons η ∈ [−4, 10] via B0 / off-momentum / Roman Pots; evaporation neutrons at few-mrad angles (≲6 mrad at 50 GeV/n, ~3× smaller at 110 GeV/n); few low-angle evaporation protons because of the Coulomb barrier.
- **Coherent/intact recoil: explicitly not available** (quoted in §2g above), and listed as future work.
- FLUKA closed-source caveat: *"BeAGLE… has no handle on changing the evaporation process and can only adjust the INC in the previous step."*

### Key statements — polarization
**None.** Every occurrence of "spin"/"polarized" is EIC boilerplate (*"the first collider to scatter polarized electrons off polarized light ions"*, *"highly polarized (70%) electron, proton, and light-ion beams"*) except one technical use of R_VMD = ratio of longitudinally to transversely polarized *vector meson* cross sections (Eq. A1, PARP(165)=0.47679, PARP(166)=0.67597) — photon/VMD polarization, unrelated to beam or target spin. **No g₁, A₁, b₁, tensor polarization, or helicity anywhere in the model description.** This matches the code exactly.

### Other useful numbers from the paper
τ₀ default = 10 fm (tuned to E665 μPb neutron multiplicity with a 24 % coherent-fraction correction, ⟨N_n⟩ = 4.7 ± 0.5), but μXe charged multiplicities prefer τ₀ ≈ 2–3 fm — an acknowledged tension. PARJ(170) = 0.32 is a **non-standard parameter the authors added to PYTHIA** to decouple beam-remnant hadron p_T from PARJ(21). MSTP(17)=6 and PARP(166) are likewise non-standard. ⟨d⟩ = 4.402 fm for min-bias ePb. BeAGLE underestimates charged multiplicity everywhere, especially for negatives and in the current-fragmentation region.