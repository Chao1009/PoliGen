# Packaging LiPolGen — what was actually run, and what it produced

Everything on this page was executed on **2026-09-05** on one machine —
Ubuntu 22.04.5 LTS, glibc 2.35, gcc 11.4.0, cmake 3.22.1, Python 3.11.4,
8 cores — and every number is the output of the command printed above it.
Nothing here is an estimate unless it says so in the same sentence, and
nothing was uploaded anywhere.

Three separate things live here because they are three separate promises:

1. **§1 the dependency stack** — one script that builds HepMC3 + LHAPDF
   (+ its PDF sets) + PYTHIA 8 into one prefix, and how long that takes.
2. **§2 the wheel** — `pip wheel`, `auditwheel repair`, and exactly how far
   the repaired wheel travels (further than the tree used to claim for the
   libraries, not as far as it claims for the *data*).
3. **§3 the licence consequence** — vendoring GPL'd `.so`s into a wheel.
   This is a statement about what is already true, not a change.

Continuous integration is `.github/workflows/ci.yml`. That file has **never
been run on GitHub**; its own header says so and lists what that leaves
unverified. §1 and §2 below are the parts of it that *were* run, by hand,
here.

---

## 1. The dependency stack — `.github/scripts/build_deps.sh`

`env.sh` expects one prefix at `<repo>/../deps/install` containing all three
dependencies. The script that builds it is
`.github/scripts/build_deps.sh <prefix> [jobs]`. It is the two ad-hoc scripts
the author's own prefix was built by — `../deps/src/build_deps.sh` and
`../deps/src/build_pythia.sh` — merged, parameterised by prefix, given the
download and PDF-set steps they assumed away, and given one correction (see
"one real defect", below).

```bash
git clone <this repo> LiPolGen
LiPolGen/.github/scripts/build_deps.sh "$PWD/deps/install" 8
```

### 1.1 What it fetches, and from where

Every URL below was checked to return HTTP 200 on 2026-09-05 (`curl -sSI -L`).
The one the tree's older notes would have used for PYTHIA —
`https://pythia.org/download/pythia83/pythia8317.tgz` — **404s**, as does
that whole `/download/` path for every 8.3 release tried (8310 … 8318);
`pythia.org/releases.html` now links `/releases/pythia83/…`, which is what
the script uses.

| what | URL | bytes (`content-length`) | sha256 pinned in the script? |
|---|---|---|---|
| HepMC3 3.3.0 | `https://hepmc.web.cern.ch/hepmc/releases/HepMC3-3.3.0.tar.gz` | 9 341 637 | yes — `6f876091…5769d6` |
| LHAPDF 6.5.5 | `https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.5.tar.gz` | (302 → downloader) | yes — `641d5ea0…f5e3a5` |
| PYTHIA 8.317 | `https://pythia.org/releases/pythia83/pythia8317.tgz` | 31 291 506 | **no — see below** |
| set `CT18NLO` | `https://lhapdfsets.web.cern.ch/current/CT18NLO.tar.gz` | 18 756 939 | no |
| set `NNPDFpol11_100` | `https://lhapdfsets.web.cern.ch/current/NNPDFpol11_100.tar.gz` | 20 128 894 | no |
| set `EPPS21nlo_CT18Anlo_Li6` | `https://lhapdfsets.web.cern.ch/current/EPPS21nlo_CT18Anlo_Li6.tar.gz` | 26 352 280 | no |

The two pinned checksums are the sha256 of the tarballs still sitting in
`../deps/src/` on the author's machine, i.e. hashed from the artefact this
tree was actually built against. The HepMC3 `content-length` above is
**equal** to that local file's size, so the URL serves the same artefact.
PYTHIA's tarball had already been deleted after unpacking, so there is
nothing here to hash it against and it is left **unpinned on purpose** rather
than given a checksum copied off a web page. The PDF sets are unpinned for the
same reason (`/current/` is a moving target by construction).

Total download for a cold run: **≈ 105 MiB** of sources and **62.2 MiB** of
PDF sets. Measured throughput from this machine to `lhapdfsets.web.cern.ch`:
5.6 and 8.9 MB/s on two consecutive full fetches of `CT18NLO.tar.gz`, i.e.
the downloads are ~15-25 s in total and not the bottleneck.

### 1.2 One real defect this found, and fixed

HepMC3 3.3.0's `HEPMC3_ENABLE_PYTHON=ON` defaults its install path for the
Python bindings to **the interpreter's own `site-packages`**, not to
`CMAKE_INSTALL_PREFIX`, and warns about it at configure time:

```
-- HepMC3 python: WARNING: The installation path of the python modules is
   outside of the global instalation path CMAKE_INSTALL_PREFIX=...
```

`../deps/src/build_deps.sh` did not set it, so `make install` wrote
`pyHepMC3` outside the prefix — verified on the author's machine
2026-09-05: `pyHepMC3.so` sits in the interpreter's `site-packages` dated
`2026-08-29 10:43`, the day the deps were built, while
`<prefix>/lib/python` is an **empty directory**. Harmless there; wrong for a
cache whose whole content is supposed to be the prefix. The script therefore
passes `-DHEPMC3_Python_SITEARCH311="$PREFIX/lib"` (the `311` suffix is
computed from the running interpreter), which is exactly where `env.sh`
already points `PYTHONPATH`.

Worth knowing while reading that: **nothing in LiPolGen's suite imports
`pyHepMC3`.** `python/tests/test_hepmc.py` uses the unrelated PyPI package
`pyhepmc`. The bindings are built for parity with `env.sh`, not for a test.

### 1.3 Measured: a cold build of the whole stack

```bash
$ time .github/scripts/build_deps.sh <scratch>/deps/install 8
```

**337.28 s wall (5 min 37 s)**, exit 0, on 8 cores. That figure is short by
exactly one download: HepMC3's 9 341 637-byte tarball was already in the
source directory from an earlier aborted run of the same script and was
re-used rather than re-fetched (≈1–2 s at the rate measured in §1.1).
Everything else is inside it — all three `configure`s, all three builds, all
three installs, PYTHIA's 31.3 MiB tarball and the whole 62.2 MiB of PDF sets.

Stage boundaries reconstructed from the mtimes of the artefacts each stage
wrote (`hepmc3-build/install_manifest.txt`, `LHAPDF-6.5.5/src/libLHAPDF.la`,
the set tarballs, `pythia8317/Makefile.inc`, `pythia8317/lib/*.so`). The
script's own log markers are **not** used for this: the pipeline that
timestamped them buffers, so they cluster and are wrong by up to ~15 s.

| stage | wall | what dominates it |
|---|---:|---|
| HepMC3 3.3.0 — cmake configure, build, install | **57.7 s** | its Python bindings. `libHepMC3.so.4` was linked **14 s** in; the other 44 s is `pyHepMC3` |
| LHAPDF 6.5.5 — download, extract, configure, build, install | **49.6 s** | |
| the three PDF sets — download 62.2 MiB, untar | **≈12 s** | network, not CPU |
| PYTHIA 8.317 — download 31.3 MiB, extract, configure, build, install | **≈218 s** | its Python bindings again. `libpythia8.so` was linked **70 s** in; the other ~146 s is the single enormous `pythia8.so` binding translation unit, which does not parallelise at all |
| **total** | **337.28 s** | |

Read that table twice: **more than half the wall clock is Python bindings that
nothing in this tree imports** (§1.2). `-DHEPMC3_ENABLE_PYTHON=OFF` plus
dropping `--with-python-config` would cut a cold build to roughly two minutes.
That is deliberately *not* done, because `env.sh` puts `$LIPOLGEN_DEPS/lib` on
`PYTHONPATH` and clearly intends them to exist — but it is the obvious lever
if CI time ever matters more than that parity.

**Extrapolating to a GitHub runner** (`ubuntu-22.04`, 4 cores instead of 8):
the parallel stages roughly double, the two serial binding TUs barely move, so
a cache MISS should land near **8–12 minutes** plus the cache save. That is an
**extrapolation**, not a measurement — this workflow has never run on GitHub.
A cache HIT skips the whole job and costs one ~114 MiB restore (§1.5).

### 1.4 Measured: LiPolGen built and tested against that fresh prefix

§1.3 proves the prefix builds. What matters is whether LiPolGen is *the same
program* against it. So the repository was copied (minus `.git`, `build/`,
`.skbuild/`) to `<scratch>/ci-clone/LiPolGen`, with §1.3's prefix as its
sibling `<scratch>/ci-clone/deps` — the exact layout
`.github/workflows/ci.yml` produces with `actions/checkout` and
`path: LiPolGen`. **`source env.sh` then needed no CI-specific patching at
all**: it resolved `LIPOLGEN_DEPS`, `PYTHIA8DATA`, `LHAPDF_DATA_PATH`,
`LD_LIBRARY_PATH` and `PYTHONPATH` onto the new prefix by itself.

`ldd build/lipolgen_tests` confirms the binary is linked against the *new*
prefix, not the author's:

```
libLHAPDF.so   => <scratch>/ci-clone/deps/install/lib/libLHAPDF.so
libHepMC3.so.4 => <scratch>/ci-clone/deps/install/lib/libHepMC3.so.4
libpythia8.so  => <scratch>/ci-clone/deps/install/lib/libpythia8.so
```

| step | command | result |
|---|---|---|
| configure | `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release` | **3.08 s**, exit 0 |
| build | `cmake --build build -j8` | **95.93 s**, **0** lines matching `warning\|error`, exit 0 |
| C++ suite | `./build/lipolgen_tests` | **400 cases / 17 240 262 assertions / 1 skipped / 0 failed**, 154.81 s |
| Python suite | `python -m pytest python/tests -q` | **923 passed, 114 skipped, 0 failed**, 128.63 s |
| docs gate | `python3 validation/check_physics_channels_links.py` | `1145 … 97 ranges, 7 external, 0 broken, 6 allow-listed` + `0 broken` in both extra documents |
| SPDX gate | `python3 validation/check_spdx_headers.py` | 101 files, 0 missing |

The doctest line is **bit-for-bit the tallies the author's own prefix
produced at the commit this was run against** (`7f68339`) —
400 / 17 240 262 / 1 skipped / 0 failed; phase F then added one case, so both
trees read 401 / 17 240 286 at `5af0427`.  **The tally has moved twice since
and this page's table is not re-measured with it** — the 2026-09-06 run's
close-out measures **408 cases / 17 241 975 assertions / 0 skipped / 0
failed**, **1085 pytest passed / 151 skipped**, and the two gate rows above
read **1240 strict / 95 ranges / 7 external / 0 broken / 24 allow-listed**
(+ 19/115 and 8 external in the extra documents) and **SPDX 103/103**.  The
CI clone has not been re-run at those commits, so the bit-for-bit claim below
stands **as of `7f68339` only**. That is the strongest single statement
on this page: a dependency stack rebuilt from four URLs by
the script in §1 reproduces, to the digit, the environment this generator's
rtol-1e-12 pinned numbers were measured in.

The Python suite read **923 / 114** where the repository itself read
**925 / 112** at that same commit (phase F then added one pytest: the tree
reads **926 / 112** now). Both extra skips were identified rather than
assumed:

* `test_release_metadata.py:121` — *no 'origin' remote configured in this
  checkout*. The copy has no `.git`. **On GitHub this one would RUN**:
  `actions/checkout` creates a `.git` with an `origin`.
* `test_hfs.py:125` — *polligen is not importable*. `conftest.py`'s `polligen`
  fixture looks for a sibling `PolarizedLithiumSim/` checkout. **On GitHub
  this one would still skip**, for the same reason.

No test failed, and no test that ran gave a different answer.

### 1.5 What that costs a CI cache

| | bytes | note |
|---|---:|---|
| `deps/install` on disk | 448 336 030 (**431 MiB**) | 273 MiB of it is the three PDF sets |
| the same as `tar \| zstd -3` | 119 039 665 (**113.5 MiB**) | what `actions/cache` stores, near enough |
| — of which, everything but the PDF sets | 47 652 839 | |
| — of which, the PDF sets alone | 71 390 593 | |

One cache entry of ~114 MiB against a repository budget of 10 GiB. The
workflow keys it on the three versions **and** the set names, so bumping any
one of them is what forces a rebuild.

---

## 2. The wheel

### 2.1 Building it

```bash
export LIPOLGEN_DEPS_PREFIX=/home/cpeng/Projects/polli/deps/install
pip wheel . --no-deps -w wheelhouse
```

**87.77 s**, producing `lipolgen-0.1.0-cp311-cp311-linux_x86_64.whl`,
**1 779 767 bytes** (36 entries, 5 489 876 bytes uncompressed). It contains
LiPolGen's own five `.so`s (`_lipolgen`, `libLiPolGenCore`, `…HepMC`,
`…LHAPDF`, `…Pythia`), the five `.py` files and the 20 `data/vmc` tables —
and **not one** of HepMC3, LHAPDF or PYTHIA8.

```bash
auditwheel show wheelhouse/lipolgen-0.1.0-cp311-cp311-linux_x86_64.whl
```

> `…is consistent with the following platform tag: "linux_x86_64".`
> … `This constrains the platform tag to "manylinux_2_35_x86_64".`

`manylinux_2_35` is glibc 2.35 = this machine's Ubuntu 22.04. To get a wider
tag the wheel has to be built inside an older `manylinux` image; `auditwheel`
says so itself and nothing here changes that.

The reason the raw wheel is not portable, straight from the ELF headers:

```
$ readelf -d lipolgen/_lipolgen.cpython-311-x86_64-linux-gnu.so
 (RUNPATH)  $ORIGIN:/home/cpeng/Projects/polli/deps/install/lib
 (NEEDED)   libLiPolGenLHAPDF.so  libLiPolGenPythia.so  libLiPolGenHepMC.so …
$ readelf -d lipolgen/libLiPolGenPythia.so
 (NEEDED)   libpythia8.so         # resolved only by that absolute RUNPATH
```

### 2.2 Repairing it

`auditwheel` 6.8.2 was installed into a throwaway venv; it needs the
`patchelf` **binary** on `PATH` and does not vendor one, so the PyPI
`patchelf` wheel (0.19.1) supplies it:

```bash
python3 -m venv venv && ./venv/bin/pip install auditwheel patchelf
./venv/bin/auditwheel repair -w repaired wheelhouse/lipolgen-0.1.0-cp311-cp311-linux_x86_64.whl
```

**5.35 s.** Output: `lipolgen-0.1.0-cp311-cp311-manylinux_2_35_x86_64.whl`,
**7 232 971 bytes** (53 entries, 21 955 667 uncompressed) — **4.06×** the raw
wheel. What it vendored, into a new `lipolgen.libs/` directory:

| vendored library | uncompressed | in the wheel | upstream licence |
|---|---:|---:|---|
| `libpythia8-6b6416ce.so` | 13 999 009 | 4 629 974 | GPL-2.0-or-later |
| `libLHAPDF-b6160c34.so` | 1 238 777 | 421 568 | GPL-3.0 |
| `libHepMC3-f2bce5c3.so.4` | 1 175 953 | 398 671 | GPL-3.0 |

and it rewrote the RPATHs of the four `.so`s that need them:

```
lipolgen/_lipolgen.cpython-311-…so   RPATH $ORIGIN:$ORIGIN/../lipolgen.libs
lipolgen/libLiPolGenHepMC.so         RPATH $ORIGIN:$ORIGIN/../lipolgen.libs
lipolgen/libLiPolGenLHAPDF.so        RPATH $ORIGIN:$ORIGIN/../lipolgen.libs
lipolgen/libLiPolGenPythia.so        RPATH $ORIGIN:$ORIGIN/../lipolgen.libs
```

`auditwheel show` on the repaired wheel reports it *consistent with*
`manylinux_2_35_x86_64`, and the filename and `WHEEL` tag both carry that tag.
All 20 `data/vmc` tables survive the repair.

**One thing it does not rewrite**, and this is worth knowing before publishing
anything: `lipolgen/libLiPolGenCore.so` keeps its original
`RUNPATH $ORIGIN:/home/cpeng/Projects/polli/deps/install/lib`, because that
library needs nothing from the prefix and so `auditwheel` never touches it.
The path is inert — but the **builder's home directory is a literal string in
a published artefact**. Four other strings of the same kind survive too:

```
lipolgen/libLiPolGenPythia.so        …/deps/install/share/Pythia8/xmldoc
lipolgen.libs/libpythia8-….so        …/deps/install/share/Pythia8/xmldoc
lipolgen.libs/libLHAPDF-….so         …/deps/install/share
lipolgen.libs/libHepMC3-….so.4       …/deps/install/lib
```

The first three are not cosmetic. They are §2.4.

### 2.3 Measured: does the repaired wheel import somewhere else?

Installed into a clean venv, run from `/tmp`, with `env -i` — no
`LD_LIBRARY_PATH`, no `PYTHONPATH`, no `PYTHIA8DATA`, no `LHAPDF_DATA_PATH`:

```
LD_LIBRARY_PATH = None      PYTHONPATH       = None
PYTHIA8DATA     = None      LHAPDF_DATA_PATH = None
version 0.1.0
tiers HAVE_LHAPDF/HEPMC3/PYTHIA8 = True True True
run ok: <x> = 0.027023337580082865
```

and the full T2 chain through the console script, same empty environment:

```
$ lipolgen-run --isotope 6Li --config 1 --channel tagged-alpha \
      --plan tensor-thirds --events 25 --seed 1 --hepmc smoke.hepmc --hadronize
  PYTHIA: 25 ok, 0 failed, 0 retries
  wrote smoke.hepmc (0.1 MB)
```

That is not yet proof of portability, because
`/home/cpeng/Projects/polli/deps/install` still **exists on this machine** and
the strings in §2.2 still resolve. So the same commands were re-run with that
directory replaced by an empty tmpfs inside a `bwrap` user namespace — the
closest thing to "a machine that never had the deps" available without a
second machine.

**With the prefix gone and nothing exported** (`os.path.isdir(prefix)` →
`False`, verified in-process):

| | result |
|---|---|
| `import lipolgen`, three tiers True | **works** — all three vendored `.so`s load via `$ORIGIN/../lipolgen.libs` |
| `lg.run(channel='inclusive', …)` | **works**, `<x> = 0.027023337580082865` — *bit-identical* to the unmasked run |
| `lg.run(channel='tagged-alpha', cluster_wave='vmc')` | **works** — the wheel's own `data/vmc` tables are found with nothing set |
| `LhapdfSF('CT18NLO')` | **FAILS**: `RuntimeError: Couldn't find required lhapdf.conf system config file`, then LHAPDF's own `terminate` |
| `lipolgen-run --hadronize` | **FAILS**: `PYTHIA Error in Settings::init: settings file /home/cpeng/Projects/polli/deps/install/share/Pythia8/xmldoc/Index.xml not found` → `PYTHIA Abort from Pythia::Pythia: settings unavailable`, exit 1 |

**With the prefix gone but the two data trees present at a completely
different path** and the three variables exported at it:

```bash
export PYTHIA8DATA=$DATA/Pythia8/xmldoc
export LHAPDF_DATA_PATH=$DATA/LHAPDF
export LIPOLGEN_PYTHIA8_PDFDATA=$DATA/Pythia8/pdfdata
```

| | result |
|---|---|
| `LhapdfSF('CT18NLO').f2p(0.1, 5.0)` | **0.4238786262338316** — identical to the unmasked value |
| `lipolgen-run … --hadronize`, 25 events | **works**, and the HepMC3 file is **byte-identical** to the unmasked run (sha256 `68ed6f06…e405ac`, 93 662 bytes) |

### 2.4 So: what does the repaired wheel actually promise?

- **Libraries: solved.** HepMC3, LHAPDF and PYTHIA8 travel inside the wheel.
  Nothing needs `LD_LIBRARY_PATH`; the deps prefix can be deleted.
- **Data: not solved, and not solvable by `auditwheel`.** PYTHIA's `xmldoc`
  (3.7 MiB) and `pdfdata` (52 MiB), and LHAPDF's `lhapdf.conf` + set store
  (273 MiB here), are read at *run* time from paths compiled into
  `libpythia8`/`libLHAPDF`/`libLiPolGenPythia` at *build* time, and those
  paths are the build machine's. A wheel handed to a second machine imports,
  and does everything that needs neither library's data — the whole pure-C++
  generator, including the VMC tables, which **are** in the wheel — and stops
  at the first `LhapdfSF` or `PythiaBridge` with the two messages quoted
  above.
- The fix is three environment variables and a copy of the two data trees, at
  any path. That is measured, not asserted, in the last table of §2.3.
- Shipping the data trees inside the wheel instead would mean a ~300 MiB wheel
  and raises a separate redistribution question for the LHAPDF grids. Not
  done, not recommended here without a decision.

---

## 3. The licence consequence of a repaired wheel

`auditwheel repair` copies **libpythia8 (GPL-2.0-or-later), libLHAPDF
(GPL-3.0) and libHepMC3 (GPL-3.0)** into the distributable artefact. That makes
the wheel a combined/linked work of all three, and the union of those terms is
**GPL-3.0-or-later** — the licence this project already carries
(`LICENSE`, `CITATION.cff`, and the `SPDX-License-Identifier: GPL-3.0-or-later`
header on all 101 covered source files).

So this is a **documentation statement, not a change**: nothing about
LiPolGen's licensing moves because of packaging. What it does mean concretely,
for whoever publishes such a wheel:

- the recipient gets it under GPL-3.0-or-later **for the combination**,
  including the vendored libraries;
- corresponding source must be available for all of it, not only LiPolGen's
  own half;
- no additional restrictions may be layered on top;
- the same obligations propagate to anyone who redistributes it further.

The reasoning behind choosing GPL-3.0-or-later in the first place, and why a
permissive licence on LiPolGen's own files would not have avoided any of this,
is `docs/open_items/engineering.md` §C and `docs/OPEN_ITEMS_SOLUTIONS.md`
§12–13.

**Nothing has been published.** No wheel from this tree has been uploaded to
PyPI or anywhere else, and `.github/workflows/ci.yml` deliberately contains no
upload step at all — not even `actions/upload-artifact` — so that pushing the
workflow is not, by itself, an act of distribution. Publishing is the author's
decision; it is recorded as a row in
`docs/open_items/run_2026-09-03/STATUS.md`'s author-decision table.
