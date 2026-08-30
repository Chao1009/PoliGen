# LiPolGen engineering investigation — findings

All prototypes are under `/tmp/claude-1000/-home-cpeng-Projects-polli/ad581d75-f90b-4f90-9f7b-196cb9e9bda5/scratchpad/` (`eng/` for the smoke-test artifacts, `pkg/LiPolGen/` for the packaging copy). Nothing was committed or written back into `/home/cpeng/Projects/polli/LiPolGen`.

## A. ePIC chain smoke test — result: it passes, with one concrete writer-side fix to make and one known-broken external step

**Environment**: `apptainer` is not installed; `singularity-ce 3.11.2` is (`/usr/bin/singularity`). Both images exist: `~/Projects/eic-2026/local/lib/eic_xl-nightly.sif` (4.2 GB, epic-main git `9aaa2969`, ships `npsim`, `abconv`, `eicrecon`) and `~/Projects/eic/local/lib/jug_xl-nightly.sif` (3.0 GB, legacy). `./build/lipolgen-run` already existed; no rebuild needed.

**Generation** (as specified, from `/home/cpeng/Projects/polli/LiPolGen`):
```
source env.sh && ./build/lipolgen-run --isotope 6Li --config 1 --channel tagged-alpha \
  --plan tensor-thirds --events 100 --seed 1 --hepmc <scratch>/smoke.hepmc --hadronize
```
→ `PYTHIA: 100 ok, 0 failed, 0 retries`; wrote 100-event HepMC3 Asciiv3 (`<scratch>/eng/smoke.hepmc`, 396 KB).

**npsim, fed the raw LiPolGen HepMC3 file directly** (no `abconv`), `epic_craterlake_10x100.xml`, `-N 10`: **all 10 events completed**, 0 errors. It does accept the 10-digit ion PDG codes (beam 6Li `1000030060` status 4, alpha spectator `1000020040` status 1) straight off the file, confirming `T2_CHAIN.md` §2's claim.

**`abconv`**: `-p 1` ("auto read energy from source") **fails**: `0 is not a valid ion Beam Energy!! ... terminate called ... Ion energy setting is incorrect` — it cannot decode a 6Li (`Z=3,A=6`) ion PDG code at all, computing per‑nucleon energy as 0 rather than ~99.5 GeV/u (this is the documented `plans/03_phase2_full_simulation.md` limitation — abconv only snaps to `{41,100,110,130,166,250,275}` GeV/nucleon brackets keyed to known species — now empirically confirmed). Using a **manual preset** instead, `abconv smoke.hepmc -o smoke_manual -p ip6_hiacc_100x10 -f hepmc3`, works: applies the 25 mrad/100 µrad crossing angle + vertex smearing to all 100 events, correctly reading `pdg 1000030060 e:597.0` (=6×99.5 GeV/u) for identification purposes even though it can't auto-derive the energy bracket. Feeding **that** file through `npsim` also completes cleanly, 10/10 events, 0 errors/warnings. EDM4hep `MCParticles` were verified (via `uproot` inside the container) to carry the full expected role chain: beam e⁻(11,status4)/6Li(1000030060,status4) → scattered e⁻(11,1) + virtual γ(22,3) + α spectator(1000020040,1) + StruckCluster deuteron(1000010020,3) + PartnerSpectator p(2212,1)/StruckNucleon n(2112,3) → HadronicX(92,3) → real PYTHIA hadrons(status1). Status/PDG mapping survives npsim→EDM4hep intact, exactly per `HEPMC3_CONVENTION.md`.

**One real writer-side issue found**: `npsim`/DD4hep logs, several times per 10 events, `createG4Primary INFO Change particle e- energy from <E> MeV by <~15-40> ppm to avoid negative Energy^2`. Root cause traced to source: `src/core/generator.cpp:114` (`beam_e.mass = 0.0;`) and `:132` (`escat.mass = 0.0;`) — LiPolGen's electron 4-vectors are built exactly massless (the standard DIS approximation) and that `0.0` is written through verbatim as HepMC3 `generated_mass`. Geant4's particle table assigns the electron its real PDG mass (0.511 MeV) regardless, so at floating-point precision `E²-p²` for the primary comes out (very slightly) below the real `m_e²`, and DD4hep silently bumps `E` by O(10) ppm to compensate. **Non-fatal** (all events still process and conserve to Geant4's working precision) but worth fixing on the writer side if a clean chain-gate log matters: either write `m_e = 0.000511` GeV as `generated_mass` for `BeamElectron`/`ScatteredElectron` (cosmetic, since the momentum components are still built massless) or just document/accept the ppm-level correction as harmless.

**Recommendation for the real chain**: run `abconv` with an explicit `ip6_hiacc_100x10`-style manual preset (not `-p 1` auto-detect) until either abconv gains a 6Li species entry or LiPolGen's own writer applies the crossing-angle/vertex afterburner upstream (the ⁶Li equivalent of `ion_gun_hepmc.py`'s hardcoded `XING_IP6`). This mirrors `plans/03`'s already-documented options (a)/(b) for the missing `ip6_eLi6_*` preset — now confirmed as the actual failure mode rather than a predicted one.

Artifacts: `eng/smoke.hepmc`, `eng/run/smoke.edm4hep.root` (direct), `eng/abconv2/smoke_manual.hepmc`, `eng/run3/smoke_manual.edm4hep.root` (abconv→npsim), `eng/run3/npsim.log`.

## B. Packaging prototype — `pip install -e .` works in ~96 s; real gaps are RPATH portability and vendored data

Copied the repo (minus `build/`, `.git`) to `pkg/LiPolGen/`. Added to `pyproject.toml`:
- `[build-system] requires = ["scikit-build-core>=0.9", "pybind11>=2.10"]`
- `[tool.scikit-build] wheel.packages = ["python/lipolgen"]`, `cmake.args = ["-DLIPOLGEN_BUILD_TESTS=OFF"]`

And three small `CMakeLists.txt` additions (scratchpad copy only): `CMAKE_INSTALL_RPATH = "$ORIGIN:${LIPOLGEN_DEPS_PREFIX}/lib"` (+ `CMAKE_BUILD_WITH_INSTALL_RPATH`/`_USE_LINK_PATH`), and `install(TARGETS _lipolgen lipolgen_core lipolgen_lhapdf lipolgen_hepmc lipolgen_pythia LIBRARY DESTINATION lipolgen)`.

```
pip3 install scikit-build-core[pyproject]   # pybind11 3.1.0 / numpy 1.25.1 already present
pip3 install --no-build-isolation -e . \
  --config-settings=cmake.define.LIPOLGEN_DEPS_PREFIX=/home/cpeng/Projects/polli/deps/install
```
→ **succeeded in 1m36s**, built the full C++ core + `_lipolgen` extension, produced `lipolgen-0.1.0-cp311-cp311-linux_x86_64.whl` (1.23 MB). Verified from `/tmp` (outside the repo, no `env.sh`, no `PYTHONPATH`/`LD_LIBRARY_PATH`/`PYTHIA8DATA`/`LHAPDF_DATA_PATH` set):
```python
import lipolgen as lg
lg._lipolgen.HAVE_LHAPDF, lg._lipolgen.HAVE_HEPMC3, lg._lipolgen.HAVE_PYTHIA8  # -> True True True
lg.run(channel='inclusive', isotope='6Li', events=1000, seed=1)['x'].mean()     # works
lg.PythiaBridge(...); lg.LhapdfSF('CT18NLO')                                    # both work
```
PYTHIA8's xmldoc and LHAPDF's grid path resolved with **no env vars at all** — both libraries bake in their own build-time default data paths, so that part is already robust.

**Remaining issues, confirmed by inspection**:
1. **RPATH is an absolute, non-relocatable path.** `readelf -d _lipolgen*.so` shows `RUNPATH: $ORIGIN:/home/cpeng/Projects/polli/deps/install/lib` — `$ORIGIN` correctly finds LiPolGen's own bridge libs (installed alongside the extension), but `ldd` shows `libHepMC3.so.4`, `libpythia8.so`, `libLHAPDF.so` resolve only via that hardcoded machine-local absolute path. The wheel built here (`pip wheel .`, inspected with `python3 -m zipfile -l`) contains **only** LiPolGen's own 4 `.so`s (dupliacted under `lib/` and `lipolgen/` — a redundancy to clean up) plus headers and the `.py`/extension — **it does not vendor HepMC3/PYTHIA8/LHAPDF at all**. This wheel is not installable/runnable on any other machine as-is.
2. **Fix path**: run `auditwheel repair` (Linux) after `pip wheel .` to copy the third-party `.so`s into the wheel and rewrite RPATHs to `$ORIGIN`-relative — standard for compiled-extension PyPI wheels, but raises the GPL redistribution question directly (see C).
3. **Data files not addressed by packaging at all**: PYTHIA8's `xmldoc` and LHAPDF's grid files (CT18NLO etc., can be 100s of MB) are referenced via `LIPOLGEN_DEPS_PREFIX`-derived compiled-in absolute paths, not shipped in the wheel; a redistributable wheel needs to either vendor them (large, and its own licensing question for LHAPDF grids) or document that users need a system/conda install of HepMC3+PYTHIA8+LHAPDF (the more conventional route for compiled HEP packages — c.f. why most of this ecosystem ships via `conda-forge`, not PyPI wheels).
4. `data/vmc/*` (VMC tables) are **not** referenced anywhere in `src`/`include`/`python` yet (confirmed by grep) — not a packaging concern today, matches DEVELOPMENT_PLAN §6 ("VMC overlaps ... external").
5. `sf_tables.cpp`'s digitized curves are compiled in verbatim — no runtime data-file dependency there, good.

Files: `pkg/LiPolGen/pyproject.toml`, `pkg/LiPolGen/CMakeLists.txt` (diffed from the repo's), `eng/pip_install.log`, `eng/pip_wheel.log`, `eng/wheelout/lipolgen-0.1.0-cp311-cp311-linux_x86_64.whl`.

## C. Licensing

PYTHIA 8 = GPL-2.0-or-later, HepMC3 = GPL-3.0 (confirmed directly from the Sartre bundle's convention as a data point, and standard knowledge), LHAPDF = GPL-3.0. Because GPL is viral through linking (not just distribution of the GPL'd code itself), **any option weaker than GPL-3.0-or-later for LiPolGen's own code does not avoid the constraint** — the moment LiPolGen is *distributed* as a combined/derivative work with PYTHIA8+HepMC3+LHAPDF statically or dynamically linked, the combination as a whole is bound by the union of all three licenses, whose effective floor is GPL-3.0 (GPL-3 is not compatible with distributing under GPL-2-only, so the union's practical license is GPL-3-or-later regardless of what MIT/BSD/LGPL header LiPolGen's own files carry).

- **MIT/BSD/LGPL on LiPolGen's own code**: legally fine to put on the LiPolGen-authored files themselves (permissive licenses are just weaker requirements layered under a stronger copyleft from the deps), but gives **no practical freedom** to a downstream binary/wheel consumer — anyone who receives the *linked* binary still receives it under GPL-3 obligations (source availability, same-license redistribution of the combination) because of PYTHIA8/HepMC3/LHAPDF. It mainly matters if someone later strips those deps out and links LiPolGen's core against something else — the permissive header travels with just the LiPolGen-authored source in that scenario.
- **GPL-2-or-later**: incompatible as the *sole* license for the combined work, since HepMC3/LHAPDF are GPL-3-only-compatible-or-later and GPL-2-only code cannot be relicensed into a GPL-3 combination; "GPL-2-or-later" (not "-only") sidesteps this by permitting upgrade to GPL-3, effectively converging on the same result as choosing GPL-3 outright.
- **GPL-3-or-later**: the license that actually matches what's already true of any distributed binary/wheel — simplest, no compatibility gymnastics, and what's expected for MCnet-community glue code around PYTHIA8.
- **(i) binaries/wheels**: any distributed compiled artifact (including the wheel prototyped in B, once it vendors the third-party `.so`s) must ship under GPL-3-or-later terms — source availability on request, no additional restrictions, same-license propagation for anyone who redistributes it further; this is unavoidable regardless of LiPolGen's own header, so choosing GPL-3-or-later for LiPolGen's own files is honest rather than restrictive.
- **(ii) paper-code release**: fully compatible with a public GitHub release under GPL-3-or-later (the norm for HEP generator codebases — PYTHIA8 itself, Herwig, Sherpa are all GPL); no tension with journal/arXiv code-availability requirements, which don't mandate permissive licensing.
- **(iii) MCnet guidelines**: the MCnet Guidelines (governing PYTHIA8 and sibling generators) explicitly expect derivative/linked tools to credit and cite the underlying generators and to respect their licenses; GPL-3-or-later on LiPolGen's own code is the path of least friction with that expectation and with the generator community's norms, while a permissive license would create a mismatch between what LiPolGen's own LICENSE file says and what a recipient's actual rights are once PYTHIA8/HepMC3/LHAPDF are linked in.

**Recommendation** (three sentences): License LiPolGen's own code GPL-3.0-or-later — it's the only choice that isn't misleading about what a recipient of the linked binary actually gets, since PYTHIA8/HepMC3/LHAPDF already force the combined work to GPL-3 terms regardless of what's stamped on LiPolGen's own files. This costs nothing for the paper-code release (GPL-3 is standard and expected in this exact niche — PYTHIA8/Herwig/Sherpa are all GPL) and is the friction-free choice under the MCnet Guidelines' citation/license-respect norms. The only reason to reach for something else is a hypothetical future core (physics kernels, spin formalism, sampler) meant to be reused independently of the T2/PYTHIA bridge — that part alone could be dual-licensed or split into a separately-licensed permissive module, but the moment PYTHIA8/HepMC3/LHAPDF are linked in, GPL-3.0-or-later is the honest label for the result.

## D. Sartre

**Not installed system-wide, not on PyPI** (`pip3 index versions sartre` / `pip3 download sartre-eic` → no matching distribution; there is no `sartre-eic` GitHub org — `github.com/sartre-eic` is a 404). No pure-`pip install` path exists or has ever existed for it. `github.com/eic/SARTREdataset` exists but is pre-generated event/lookup-table **data** releases, not the generator source. Historically the source was hosted at Google Code (`code.google.com/p/sartre-mc`, defunct) and distributed as a tarball from the Yale RHIG group's site (`rhig.physics.yale.edu/~ullrich/sartre-doc/`); no current git mirror was found via GitHub/GitLab search.

A copy of the source (v1.39) was available in this session's scratchpad (`scratchpad/s139/sartre/`) and was inspected directly:
- **License**: GPL, `gpl.txt` bundled is **GPLv3** text; every source file header reads "GNU General Public License as published by the Free Software Foundation" (copyright Tobias Toll & Thomas Ullrich, 2010–2019).
- **Dependencies** (from `cmake/modules/`): **ROOT** and **GSL** are required (own `FindROOT.cmake`/`FindGSL.cmake`), **Boost** is optional (only for multithreaded table generation), **Cuba** and **Gemini++** are bundled in-tree (`cuba/`, `gemini/` subdirs) rather than external deps.
- **Nuclei supported — hardcoded, not user-defined**: `Nucleus::init(unsigned int A)` in `src/Nucleus.cpp` is a `switch(mA)` with Woods-Saxon (3-parameter Fermi) parameters hardcoded for exactly `A = 208 (Pb), 197 (Au), 110 (Cd), 90 (Zr), 63 (Cu), 40 (Ca), 27 (Al), 16 (O), 2 (D, Hulthén), 1 (p)`. **No case for A=6 (or any other light nucleus besides the deuteron)** — the `default:` branch prints `"Nucleus::init(): Error, cannot handle A=<N>, no parameters defined for this mass number."` and calls `exit(1)`. So Sartre does **not** support 6Li out of the box, and there is no generic/runtime nuclear-parameter input — adding it means editing and recompiling `Nucleus.cpp`/`TableGeneratorNucleus.cpp` with a physically appropriate 6Li density (a Woods-Saxon fit is a poor model at A=6; an α+d cluster form, closer to what LiPolGen's own T1 tier already does, would be more defensible), **then** running the separate table-generator step to build new Q²/W²/t dipole-amplitude lookup tables — a step the current (2606.14633, June-2026) Sartre paper describes as costing "a few CPU-years... per combination" per species/meson, i.e. adding 6Li is a substantial standalone project, not a config change.