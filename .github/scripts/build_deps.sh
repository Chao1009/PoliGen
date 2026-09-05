#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Build LiPolGen's dependency stack into ONE prefix, the layout `env.sh`
# expects (`<repo>/../deps/install`):
#
#   HepMC3 3.3.0      -- the HepMC3 writer tier, + its Python bindings
#   LHAPDF 6.5.5      -- the LHAPDF structure-function backend
#     + the three PDF sets the suite actually reads:
#       CT18NLO, NNPDFpol11_100, EPPS21nlo_CT18Anlo_Li6
#   PYTHIA 8.317      -- the T2 hadronization tier, linked against both of the
#                        above, with its Python bindings
#
# This is the transcription of the two ad-hoc scripts the deps on the author's
# machine were actually built by -- `../deps/src/build_deps.sh` and
# `../deps/src/build_pythia.sh` -- into one re-runnable, prefix-parameterised
# file, so CI and a fresh developer machine take the same path.  Every
# `configure`/`cmake` flag below is copied from those two scripts verbatim;
# the only additions are the download step (they assumed the tarballs were
# already unpacked), the PDF-set download (they did not do it -- `lhapdf
# install` failed there, see `../deps/lhapdf_install.log`, and the sets were
# fetched by hand), the checksum pinning, and ONE deliberate correction:
#
#   `-DHEPMC3_Python_SITEARCH<NN>="$PREFIX/lib"`.  Without it HepMC3 3.3.0
#   defaults that path to the INTERPRETER'S OWN site-packages and says so at
#   configure time -- "WARNING: The installation path of the python modules
#   is outside of the global instalation path CMAKE_INSTALL_PREFIX=..." --
#   so `make install` writes `pyHepMC3` outside the prefix entirely.  That is
#   what happened on the author's machine (verified 2026-09-05: `pyHepMC3.so`
#   sits in the interpreter's site-packages, dated the day the deps were
#   built, while `<prefix>/lib/python` is an empty directory).  Harmless
#   there; wrong for CI, where the cached artifact must be the prefix and
#   nothing else.  `$PREFIX/lib` is exactly where `env.sh` already points
#   PYTHONPATH, so the bindings stay importable.  NOTE: nothing in LiPolGen's
#   own suite imports `pyHepMC3` -- `python/tests/test_hepmc.py` uses the
#   unrelated PyPI package `pyhepmc` -- so these bindings are built for
#   parity with `env.sh`, not because a test needs them.
#
# Usage:  build_deps.sh <prefix> [jobs]
#   <prefix>  where to install (created if absent), e.g. .../deps/install
#   [jobs]    parallelism for make (default: nproc)
#
# Sources are downloaded into <prefix>/../src and left there; a second run
# with the tarballs already present re-uses them and does not re-download.
#
# MEASURED 2026-09-05 on Ubuntu 22.04.5 (glibc 2.35, gcc 11.4.0, cmake 3.22.1,
# Python 3.11.4), 8 cores, `jobs` = 8 -- see docs/PACKAGING.md for the table.

set -euo pipefail

PREFIX="${1:-}"
JOBS="${2:-$(nproc)}"

if [ -z "$PREFIX" ]; then
  echo "usage: $0 <prefix> [jobs]" >&2
  exit 2
fi

HEPMC3_VERSION="${HEPMC3_VERSION:-3.3.0}"
LHAPDF_VERSION="${LHAPDF_VERSION:-6.5.5}"
PYTHIA_VERSION="${PYTHIA_VERSION:-8317}"          # 8.317, as PYTHIA names it
LHAPDF_SETS="${LHAPDF_SETS:-CT18NLO NNPDFpol11_100 EPPS21nlo_CT18Anlo_Li6}"

# SHA-256 of the two tarballs that were still on the author's machine on
# 2026-09-05 and could therefore be hashed from the real thing.  PYTHIA's
# tarball had already been deleted after unpacking, so there is NOTHING to
# pin it against here and it is deliberately left unverified rather than
# given a checksum copied off a web page -- said plainly instead of quietly.
HEPMC3_SHA256="6f876091edcf7ee6d0c0db04e080056e89efc1a61abe62355d97ce8e735769d6"
LHAPDF_SHA256="641d5ea0942b79e4447e15e5a33491ff3c7032d71d618119935e14ad27f5e3a5"

mkdir -p "$PREFIX"
PREFIX="$(cd "$PREFIX" && pwd)"
SRC="$(dirname "$PREFIX")/src"
mkdir -p "$SRC"

PYVER="$(python3 -c 'import sys; print("%d.%d" % sys.version_info[:2])')"
PYVER_TAG="${PYVER//./}"        # 3.11 -> 311, the suffix HepMC3 names its
                                # per-version Python install variables with

say() { printf '\n=== %s\n' "$*"; }

fetch() {  # fetch <url> <dest> [sha256]
  local url="$1" dest="$2" want="${3:-}"
  if [ ! -f "$dest" ]; then
    say "download $(basename "$dest")"
    curl -fsSL --retry 3 --retry-delay 5 -o "$dest.part" "$url"
    mv "$dest.part" "$dest"
  fi
  if [ -n "$want" ]; then
    local got
    got="$(sha256sum "$dest" | cut -d' ' -f1)"
    if [ "$got" != "$want" ]; then
      echo "checksum mismatch for $dest: got $got, want $want" >&2
      exit 1
    fi
    echo "sha256 ok: $(basename "$dest")"
  else
    echo "sha256 NOT PINNED for $(basename "$dest") -- see the comment above"
  fi
}

# ------------------------------------------------------------------ HepMC3
say "HepMC3 $HEPMC3_VERSION"
fetch "https://hepmc.web.cern.ch/hepmc/releases/HepMC3-${HEPMC3_VERSION}.tar.gz" \
      "$SRC/HepMC3-${HEPMC3_VERSION}.tar.gz" "$HEPMC3_SHA256"
[ -d "$SRC/HepMC3-${HEPMC3_VERSION}" ] || tar -xzf "$SRC/HepMC3-${HEPMC3_VERSION}.tar.gz" -C "$SRC"
mkdir -p "$SRC/hepmc3-build"
cmake -S "$SRC/HepMC3-${HEPMC3_VERSION}" -B "$SRC/hepmc3-build" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DHEPMC3_ENABLE_ROOTIO=OFF \
  -DHEPMC3_ENABLE_PYTHON=ON \
  -DHEPMC3_PYTHON_VERSIONS="$PYVER" \
  -DHEPMC3_Python_SITEARCH${PYVER_TAG}="$PREFIX/lib" \
  -DHEPMC3_ENABLE_TEST=OFF \
  -DHEPMC3_INSTALL_INTERFACES=ON \
  -DHEPMC3_BUILD_STATIC_LIBS=OFF \
  -DCMAKE_BUILD_TYPE=Release
cmake --build "$SRC/hepmc3-build" -j "$JOBS"
cmake --install "$SRC/hepmc3-build"

# ------------------------------------------------------------------ LHAPDF
say "LHAPDF $LHAPDF_VERSION"
fetch "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-${LHAPDF_VERSION}.tar.gz" \
      "$SRC/LHAPDF-${LHAPDF_VERSION}.tar.gz" "$LHAPDF_SHA256"
[ -d "$SRC/LHAPDF-${LHAPDF_VERSION}" ] || tar -xzf "$SRC/LHAPDF-${LHAPDF_VERSION}.tar.gz" -C "$SRC"
(
  cd "$SRC/LHAPDF-${LHAPDF_VERSION}"
  # --disable-python: nothing in this tree imports LHAPDF from Python (the
  # backend is reached through `lipolgen_lhapdf`), and LHAPDF's own bindings
  # need Cython.  Same flag the author's build used.
  ./configure --prefix="$PREFIX" --disable-python
  make -j "$JOBS"
  make install
)

# ---------------------------------------------------- LHAPDF sets (grids)
# `lhapdf install <SET>` is NOT used: it needs the index/config the install
# does not have yet at this point (it failed exactly that way on the author's
# machine -- `No PDFs known matching patterns: CT18NLO, ...`), so the sets are
# taken straight from the same CERN store `lhapdf` would have used.
say "LHAPDF sets: $LHAPDF_SETS"
mkdir -p "$PREFIX/share/LHAPDF"
for set_name in $LHAPDF_SETS; do
  if [ -d "$PREFIX/share/LHAPDF/$set_name" ]; then
    echo "already present: $set_name"
    continue
  fi
  fetch "https://lhapdfsets.web.cern.ch/current/${set_name}.tar.gz" \
        "$SRC/${set_name}.tar.gz"
  tar -xzf "$SRC/${set_name}.tar.gz" -C "$PREFIX/share/LHAPDF"
  echo "installed set: $set_name"
done

# ------------------------------------------------------------------ PYTHIA
say "PYTHIA $PYTHIA_VERSION"
fetch "https://pythia.org/releases/pythia83/pythia${PYTHIA_VERSION}.tgz" \
      "$SRC/pythia${PYTHIA_VERSION}.tgz"
[ -d "$SRC/pythia${PYTHIA_VERSION}" ] || tar -xzf "$SRC/pythia${PYTHIA_VERSION}.tgz" -C "$SRC"
(
  cd "$SRC/pythia${PYTHIA_VERSION}"
  # --with-hepmc3/--with-lhapdf6 point at the prefix just built, so
  # libpythia8hepmc3 / libpythia8lhapdf6 exist and the T2 tier can hand
  # events straight to the HepMC3 writer.  --cxx-common is copied verbatim
  # from the author's build: -std=c++17 is what LiPolGen compiles at, -fPIC
  # is required for the shared build.
  ./configure --prefix="$PREFIX" \
    --with-hepmc3="$PREFIX" \
    --with-lhapdf6="$PREFIX" \
    --with-python-config="$(command -v python3-config)" \
    --with-gzip \
    --cxx-common='-O2 -std=c++17 -pedantic -W -Wall -Wshadow -fPIC -pthread'
  make -j "$JOBS"
  make install
)

say "done -- prefix $PREFIX"
du -sh "$PREFIX"
