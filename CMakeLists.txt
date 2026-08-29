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

# Default dependency prefix: ../deps/install (see env.sh)
set(LIPOLGEN_DEPS_PREFIX "${CMAKE_SOURCE_DIR}/../deps/install" CACHE PATH "Prefix with HepMC3/LHAPDF/PYTHIA8")
list(APPEND CMAKE_PREFIX_PATH "${LIPOLGEN_DEPS_PREFIX}")

# ---------------------------------------------------------------- core library
file(GLOB LIPOLGEN_CORE_SOURCES CONFIGURE_DEPENDS src/core/*.cpp)
add_library(lipolgen_core SHARED ${LIPOLGEN_CORE_SOURCES})
target_include_directories(lipolgen_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_options(lipolgen_core PRIVATE -Wall -Wextra -Wshadow -O2)
set_target_properties(lipolgen_core PROPERTIES OUTPUT_NAME LiPolGenCore)

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
    file(GLOB LIPOLGEN_PYTHIA_SOURCES CONFIGURE_DEPENDS src/pythia/*.cpp)
    if(LIPOLGEN_PYTHIA_SOURCES)
      add_library(lipolgen_pythia SHARED ${LIPOLGEN_PYTHIA_SOURCES})
      target_include_directories(lipolgen_pythia SYSTEM PUBLIC ${PYTHIA8_INCDIR})
      target_compile_options(lipolgen_pythia PRIVATE -Wall -Wextra -Wshadow -O2)
      target_link_directories(lipolgen_pythia PUBLIC ${PYTHIA8_LIBDIR})
      target_link_libraries(lipolgen_pythia PUBLIC lipolgen_core pythia8 dl)
      target_compile_definitions(lipolgen_pythia PUBLIC LIPOLGEN_HAVE_PYTHIA8=1
        LIPOLGEN_PYTHIA8_XMLDOC="${PYTHIA8_DATADIR}/xmldoc")
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
foreach(ex ${LIPOLGEN_EXAMPLES})
  get_filename_component(name ${ex} NAME_WE)
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
  endif()
endif()

# ---------------------------------------------------------------- install
include(GNUInstallDirs)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
foreach(t lipolgen_core lipolgen_lhapdf lipolgen_hepmc lipolgen_pythia)
  if(TARGET ${t})
    install(TARGETS ${t} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()
endforeach()
