cmake_minimum_required(VERSION 3.16)
project(LiPolGen VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(LIPOLGEN_WITH_LHAPDF  "Build the LHAPDF6 structure-function backend" ON)
option(LIPOLGEN_WITH_HEPMC3  "Build the HepMC3 writer"                       ON)
option(LIPOLGEN_WITH_PYTHIA  "Build the PYTHIA 8 hadronization tier (T2)"     ON)
option(LIPOLGEN_WITH_PYTHON  "Build the pybind11 Python module"               ON)
option(LIPOLGEN_BUILD_TESTS  "Build doctest unit tests"                       ON)

# Default dependency prefix: ../deps/install (see env.sh), or $LIPOLGEN_DEPS_PREFIX
# from the environment when the cache variable has not already been set by
# -DLIPOLGEN_DEPS_PREFIX=... / --config-settings=cmake.define.LIPOLGEN_DEPS_PREFIX=...
# (either of those pre-seeds the cache before this file runs, so this branch
# is skipped and the explicit value always wins).
if(NOT DEFINED CACHE{LIPOLGEN_DEPS_PREFIX} AND DEFINED ENV{LIPOLGEN_DEPS_PREFIX})
  set(LIPOLGEN_DEPS_PREFIX "$ENV{LIPOLGEN_DEPS_PREFIX}" CACHE PATH "Prefix with HepMC3/LHAPDF/PYTHIA8")
else()
  set(LIPOLGEN_DEPS_PREFIX "${CMAKE_SOURCE_DIR}/../deps/install" CACHE PATH "Prefix with HepMC3/LHAPDF/PYTHIA8")
endif()
list(APPEND CMAKE_PREFIX_PATH "${LIPOLGEN_DEPS_PREFIX}")

# ---------------------------------------------------------------- RPATH
# Wheel packaging only (`pip install`/`pip wheel`, SKBUILD).  Both
# CMAKE_INSTALL_RPATH and CMAKE_INSTALL_RPATH_USE_LINK_PATH must stay inside
# this `if(SKBUILD)` guard: CMake pads a target's BUILD-tree RUNPATH with
# empty (":") entries at configure time for every target that has an
# install() rule, once CMAKE_INSTALL_RPATH is set at all -- regardless of
# whether an install step ever actually runs -- so setting them
# unconditionally would put CWD-lookup entries into `build/libLiPolGenCore.so`
# and friends even for the in-tree, non-SKBUILD flow.  Keeping them under
# SKBUILD means the non-SKBUILD configure has no INSTALL_RPATH at all, so
# `build/lipolgen_tests` and `build/python` keep CMake's untouched default
# build-tree RPATH, exactly as before this packaging support existed.
# $ORIGIN lets an installed `_lipolgen*.so` find LiPolGen's own bridge
# libraries installed alongside it; the deps prefix covers libHepMC3/
# libpythia8/libLHAPDF, which are not shipped in the wheel (see docs/USAGE.md).
if(SKBUILD)
  set(CMAKE_INSTALL_RPATH "$ORIGIN:${LIPOLGEN_DEPS_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  # Belt-and-suspenders for the editable/"redirect" build: make the
  # build-dir copies (under `.skbuild/`) carry the same RPATH as the
  # installed ones too, so they stay independently loadable (e.g. for
  # debugging straight out of the build dir) rather than relying solely on
  # install-time rewriting.  The copy actually imported by `import lipolgen`
  # after `pip install -e .` is the one under site-packages/lipolgen/,
  # written by the install() rules below -- this flag does not affect that.
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()

# ---------------------------------------------------------------- core library
file(GLOB LIPOLGEN_CORE_SOURCES CONFIGURE_DEPENDS src/core/*.cpp)
add_library(lipolgen_core SHARED ${LIPOLGEN_CORE_SOURCES})
target_include_directories(lipolgen_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_options(lipolgen_core PRIVATE -Wall -Wextra -Wshadow -O2)
set_target_properties(lipolgen_core PROPERTIES OUTPUT_NAME LiPolGenCore)

# Single source of truth for the version string seen at runtime (Python
# `__version__`, HepMC3 GenRunInfo tool version, ...): CMake's own
# project(LiPolGen VERSION ...) above, not retyped anywhere downstream.
# PUBLIC so it propagates to every target that links lipolgen_core
# (lipolgen_lhapdf/hepmc/pythia, _lipolgen).
target_compile_definitions(lipolgen_core PUBLIC
  LIPOLGEN_VERSION="${PROJECT_VERSION}")

# Where the packaged data tree (VMC tables, ...) lives.  `cluster.cpp`'s
# `data_dir()` prefers $LIPOLGEN_DATA_DIR at RUNTIME and falls back to this,
# so an in-tree build finds data/ with no environment set at all.
set(LIPOLGEN_DATA_DIR "${CMAKE_SOURCE_DIR}/data" CACHE PATH
    "Compiled-in default for the runtime data tree (override at run time with \$LIPOLGEN_DATA_DIR)")
target_compile_definitions(lipolgen_core PRIVATE
  LIPOLGEN_DATA_DIR_DEFAULT="${LIPOLGEN_DATA_DIR}")

# ---------------------------------------------------------------- LHAPDF
if(LIPOLGEN_WITH_LHAPDF)
  find_program(LHAPDF_CONFIG lhapdf-config HINTS ${LIPOLGEN_DEPS_PREFIX}/bin)
  if(LHAPDF_CONFIG)
    execute_process(COMMAND ${LHAPDF_CONFIG} --incdir OUTPUT_VARIABLE LHAPDF_INCDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${LHAPDF_CONFIG} --libdir OUTPUT_VARIABLE LHAPDF_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    file(GLOB LIPOLGEN_LHAPDF_SOURCES CONFIGURE_DEPENDS src/lhapdf/*.cpp)
    if(LIPOLGEN_LHAPDF_SOURCES)
      add_library(lipolgen_lhapdf SHARED ${LIPOLGEN_LHAPDF_SOURCES})
      target_include_directories(lipolgen_lhapdf SYSTEM PUBLIC ${LHAPDF_INCDIR})
      target_compile_options(lipolgen_lhapdf PRIVATE -Wall -Wextra -Wshadow -O2)
      target_link_directories(lipolgen_lhapdf PUBLIC ${LHAPDF_LIBDIR})
      target_link_libraries(lipolgen_lhapdf PUBLIC lipolgen_core LHAPDF)
      target_compile_definitions(lipolgen_lhapdf PUBLIC LIPOLGEN_HAVE_LHAPDF=1)
      set_target_properties(lipolgen_lhapdf PROPERTIES OUTPUT_NAME LiPolGenLHAPDF)
    endif()
  else()
    message(WARNING "lhapdf-config not found; LHAPDF backend disabled")
    set(LIPOLGEN_WITH_LHAPDF OFF)
  endif()
endif()

# ---------------------------------------------------------------- HepMC3
if(LIPOLGEN_WITH_HEPMC3)
  find_package(HepMC3 QUIET HINTS ${LIPOLGEN_DEPS_PREFIX}/share/HepMC3/cmake)
  if(HepMC3_FOUND)
    file(GLOB LIPOLGEN_HEPMC_SOURCES CONFIGURE_DEPENDS src/hepmc/*.cpp)
    if(LIPOLGEN_HEPMC_SOURCES)
      add_library(lipolgen_hepmc SHARED ${LIPOLGEN_HEPMC_SOURCES})
      target_link_libraries(lipolgen_hepmc PUBLIC lipolgen_core HepMC3::HepMC3)
      target_compile_definitions(lipolgen_hepmc PUBLIC LIPOLGEN_HAVE_HEPMC3=1)
      set_target_properties(lipolgen_hepmc PROPERTIES OUTPUT_NAME LiPolGenHepMC)
    endif()
  else()
    message(WARNING "HepMC3 not found; HepMC3 writer disabled")
    set(LIPOLGEN_WITH_HEPMC3 OFF)
  endif()
endif()

# ---------------------------------------------------------------- PYTHIA 8
if(LIPOLGEN_WITH_PYTHIA)
  find_program(PYTHIA8_CONFIG pythia8-config HINTS ${LIPOLGEN_DEPS_PREFIX}/bin)
  if(PYTHIA8_CONFIG)
    execute_process(COMMAND ${PYTHIA8_CONFIG} --includedir OUTPUT_VARIABLE PYTHIA8_INCDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${PYTHIA8_CONFIG} --libdir     OUTPUT_VARIABLE PYTHIA8_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${PYTHIA8_CONFIG} --datadir    OUTPUT_VARIABLE PYTHIA8_DATADIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    # Where PYTHIA's own PDF grids live.  `xmldoc` was already compiled in;
    # `pdfdata` was not, and `mstw_sf.cpp` needs it -- the MSTW2008 LO central
    # grid CDKS used is `${PYTHIA8_DATADIR}/pdfdata/mstw2008lo.00.dat`, which
    # is on disk even though LHAPDF's set store has no MSTW2008 set.  Derived
    # from pythia8-config like the rest of the deps prefix, a CACHE PATH so a
    # build can point elsewhere, and overridden at RUN time by
    # $LIPOLGEN_PYTHIA8_PDFDATA -- exactly the LIPOLGEN_DATA_DIR arrangement.
    set(LIPOLGEN_PYTHIA8_PDFDATA "${PYTHIA8_DATADIR}/pdfdata" CACHE PATH
        "Compiled-in default for PYTHIA 8's pdfdata grid directory (override at run time with \$LIPOLGEN_PYTHIA8_PDFDATA)")
    file(GLOB LIPOLGEN_PYTHIA_SOURCES CONFIGURE_DEPENDS src/pythia/*.cpp)
    if(LIPOLGEN_PYTHIA_SOURCES)
      add_library(lipolgen_pythia SHARED ${LIPOLGEN_PYTHIA_SOURCES})
      target_include_directories(lipolgen_pythia SYSTEM PUBLIC ${PYTHIA8_INCDIR})
      target_compile_options(lipolgen_pythia PRIVATE -Wall -Wextra -Wshadow -O2)
      target_link_directories(lipolgen_pythia PUBLIC ${PYTHIA8_LIBDIR})
      target_link_libraries(lipolgen_pythia PUBLIC lipolgen_core pythia8 dl)
      target_compile_definitions(lipolgen_pythia PUBLIC LIPOLGEN_HAVE_PYTHIA8=1
        LIPOLGEN_PYTHIA8_XMLDOC="${PYTHIA8_DATADIR}/xmldoc"
        LIPOLGEN_PYTHIA8_PDFDATA="${LIPOLGEN_PYTHIA8_PDFDATA}")
      if(TARGET lipolgen_hepmc)
        target_link_libraries(lipolgen_pythia PUBLIC lipolgen_hepmc)
      endif()
      set_target_properties(lipolgen_pythia PROPERTIES OUTPUT_NAME LiPolGenPythia)
    endif()
  else()
    message(WARNING "pythia8-config not found; PYTHIA tier disabled")
    set(LIPOLGEN_WITH_PYTHIA OFF)
  endif()
endif()

# ---------------------------------------------------------------- examples
file(GLOB LIPOLGEN_EXAMPLES CONFIGURE_DEPENDS examples/*.cpp)

# An example that USES a tier cannot link when that tier is off -- the
# conditional `target_link_libraries` below only stops the library being
# added, it cannot stop the example's own calls into it.  Measured
# 2026-09-05: with -DLIPOLGEN_WITH_PYTHIA=OFF, `hadronize_example` and
# `generate_full` failed at link with `undefined reference to
# lipolgen::PythiaBridge::...`, so that whole configuration did not build at
# all.  Pair each tier with the public header only that tier's target
# provides and skip an example that includes one of an absent tier; the
# pairing is by header, not a hardcoded example list, so a new example needs
# no edit here.  With all three tiers on -- the shipped default -- nothing
# below is skipped and the generated build is unchanged.
#
# `file(READ)` is a configure-time read with no dependency edge: adding or
# removing one of these includes in an existing example needs a reconfigure
# (touching CMakeLists.txt or `cmake -S . -B build` again) to take effect.
# Adding or deleting an example file is picked up by CONFIGURE_DEPENDS above.
foreach(ex ${LIPOLGEN_EXAMPLES})
  get_filename_component(name ${ex} NAME_WE)
  file(READ ${ex} LIPOLGEN_EX_SOURCE)
  set(LIPOLGEN_EX_SKIP "")
  if(NOT TARGET lipolgen_lhapdf AND LIPOLGEN_EX_SOURCE MATCHES "lipolgen/lhapdf_sf\\.hpp")
    set(LIPOLGEN_EX_SKIP "LHAPDF")
  elseif(NOT TARGET lipolgen_hepmc AND LIPOLGEN_EX_SOURCE MATCHES "lipolgen/hepmc_writer\\.hpp")
    set(LIPOLGEN_EX_SKIP "HepMC3")
  elseif(NOT TARGET lipolgen_pythia AND LIPOLGEN_EX_SOURCE MATCHES "lipolgen/pythia_bridge\\.hpp")
    set(LIPOLGEN_EX_SKIP "PYTHIA 8")
  endif()
  if(LIPOLGEN_EX_SKIP)
    message(STATUS
      "example ${name}: NOT BUILT -- it uses the ${LIPOLGEN_EX_SKIP} tier, which is off")
    continue()
  endif()
  add_executable(${name} ${ex})
  target_link_libraries(${name} PRIVATE lipolgen_core)
  foreach(t lipolgen_lhapdf lipolgen_hepmc lipolgen_pythia)
    if(TARGET ${t})
      target_link_libraries(${name} PRIVATE ${t})
    endif()
  endforeach()
endforeach()

# ---------------------------------------------------------------- tests
if(LIPOLGEN_BUILD_TESTS)
  enable_testing()
  file(GLOB LIPOLGEN_TEST_SOURCES CONFIGURE_DEPENDS tests/test_*.cpp)
  if(LIPOLGEN_TEST_SOURCES)
    add_executable(lipolgen_tests tests/main.cpp ${LIPOLGEN_TEST_SOURCES})
    target_include_directories(lipolgen_tests PRIVATE third_party/doctest)
    target_compile_definitions(lipolgen_tests PRIVATE
      LIPOLGEN_REFERENCE_DIR="${CMAKE_SOURCE_DIR}/validation/reference")
    target_link_libraries(lipolgen_tests PRIVATE lipolgen_core)
    foreach(t lipolgen_lhapdf lipolgen_hepmc lipolgen_pythia)
      if(TARGET ${t})
        target_link_libraries(lipolgen_tests PRIVATE ${t})
      endif()
    endforeach()
    add_test(NAME lipolgen_tests COMMAND lipolgen_tests)
  endif()
endif()

# ---------------------------------------------------------------- python
if(LIPOLGEN_WITH_PYTHON)
  find_package(Python3 COMPONENTS Interpreter Development.Module QUIET)
  execute_process(COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
                  OUTPUT_VARIABLE PYBIND11_CMAKEDIR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  find_package(pybind11 CONFIG QUIET HINTS ${PYBIND11_CMAKEDIR})
  if(pybind11_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/python/bindings.cpp)
    pybind11_add_module(_lipolgen python/bindings.cpp)
    target_link_libraries(_lipolgen PRIVATE lipolgen_core)
    foreach(t lipolgen_lhapdf lipolgen_hepmc lipolgen_pythia)
      if(TARGET ${t})
        target_link_libraries(_lipolgen PRIVATE ${t})
      endif()
    endforeach()
    target_compile_options(_lipolgen PRIVATE -O2)
    set_target_properties(_lipolgen PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python/lipolgen)

    # The importable package lives in the build tree next to the extension --
    # env.sh already puts ${CMAKE_BINARY_DIR}/python on PYTHONPATH, so
    # `import lipolgen` works straight after a build with no install step.
    # The pure-Python half is copied on every build (copy_if_different, so a
    # no-op when nothing changed).
    file(GLOB LIPOLGEN_PY_SOURCES CONFIGURE_DEPENDS python/lipolgen/*.py)
    add_custom_target(lipolgen_python ALL
      COMMAND ${CMAKE_COMMAND} -E make_directory
              ${CMAKE_BINARY_DIR}/python/lipolgen
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LIPOLGEN_PY_SOURCES}
              ${CMAKE_BINARY_DIR}/python/lipolgen
      COMMENT "Staging the lipolgen Python package")
    add_dependencies(lipolgen_python _lipolgen)

    # `./build/lipolgen-run` -- the console entry point, spelled the way the
    # C++ example generators already are (./build/generate_tagged ...).  There
    # is no `pip install` step: the shim puts the build tree's package on
    # PYTHONPATH and runs `python -m lipolgen.cli`.
    file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/lipolgen-run CONTENT
"#!/bin/sh
# generated by CMake -- LiPolGen command-line generator
PYTHONPATH=\"${CMAKE_BINARY_DIR}/python:\${PYTHONPATH}\" \\
  exec \"${Python3_EXECUTABLE}\" -m lipolgen.cli \"$@\"
" FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

    # `./build/lipolgen-configs` -- the sibling shim for the Phase G
    # configuration-table writer.  In-tree runs get their console entry points
    # from `file(GENERATE ...)`, not from [project.scripts], so without this
    # line `lipolgen-configs` simply does not exist in a build tree.
    file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/lipolgen-configs CONTENT
"#!/bin/sh
# generated by CMake -- LiPolGen configuration-table writer
PYTHONPATH=\"${CMAKE_BINARY_DIR}/python:\${PYTHONPATH}\" \\
  exec \"${Python3_EXECUTABLE}\" -m lipolgen.configs \"$@\"
" FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif()
endif()

# ---------------------------------------------------------------- install
include(GNUInstallDirs)
if(SKBUILD)
  # Wheel packaging (`pip install`/`pip wheel`, scikit-build-core): the
  # extension module and LiPolGen's own shared libraries go ONCE into the
  # "lipolgen" package directory, next to each other and next to the
  # pure-Python files `wheel.packages` (pyproject.toml) already copies there
  # -- that is what makes the plain `$ORIGIN` RPATH above sufficient.  (An
  # earlier prototype installed the core libs a second time under a
  # top-level `lib/`, which just duplicated every .so in the wheel for no
  # benefit; don't repeat that here.)
  foreach(t _lipolgen lipolgen_core lipolgen_lhapdf lipolgen_hepmc lipolgen_pythia)
    if(TARGET ${t})
      install(TARGETS ${t} LIBRARY DESTINATION lipolgen)
    endif()
  endforeach()

  # VMC data tables (~1.1 MB): small enough to vendor, so a wheel installed
  # on a machine with no LiPolGen checkout still finds them.  `python/
  # lipolgen/__init__.py` points $LIPOLGEN_DATA_DIR at this directory
  # whenever it exists and the caller has not already set the variable.
  # PYTHIA8's xmldoc and LHAPDF's grids are NOT vendored here -- both are
  # much larger and already resolve from the deps prefix at run time with no
  # help from LiPolGen (see docs/USAGE.md); LiPolGen has no data of its own
  # in that category.
  #
  # In an editable ("redirect") install this copy is dead weight: __init__.py
  # is redirected to python/lipolgen in the source tree, whose data_dir()
  # then resolves against the SOURCE tree (which already has data/ right
  # there), never against this site-packages copy. Harmless, just ~1.1 MB
  # unused per editable install -- not worth branching on SKBUILD_STATE for.
  if(EXISTS "${CMAKE_SOURCE_DIR}/data")
    install(DIRECTORY "${CMAKE_SOURCE_DIR}/data/" DESTINATION lipolgen/data)
  endif()
else()
  # Conventional `cmake --install [--prefix ...]`: headers and shared
  # libraries at their usual locations.  Unused by, and unaffected by, both
  # the in-tree build/test flow and the wheel packaging above.
  install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  foreach(t lipolgen_core lipolgen_lhapdf lipolgen_hepmc lipolgen_pythia)
    if(TARGET ${t})
      install(TARGETS ${t} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
  endforeach()
endif()
