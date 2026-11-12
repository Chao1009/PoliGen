// SPDX-License-Identifier: GPL-3.0-or-later
// HepMC3 Asciiv3 writer for lipolgen::Event records.
//
// Layout and the ion-spin attribute convention are documented in
// docs/HEPMC3_CONVENTION.md (proposed convention, PolarizedLithiumSim
// plans/04 #17). The public interface here is HepMC3-free (pimpl) so that
// including this header never requires HepMC3 include paths; only
// src/hepmc/hepmc_writer.cpp needs them.
#pragma once
#include "lipolgen/event.hpp"
#include <memory>
#include <string>
#include <map>
// Defined PUBLIC on the lipolgen_core CMake target from project()'s VERSION
// (CMakeLists.txt) -- the fallback below only matters for a translation
// unit that includes this header without linking lipolgen_core at all.
#ifndef LIPOLGEN_VERSION
#define LIPOLGEN_VERSION "0.0.0"
#endif

namespace lipolgen {

enum class HepMC3Format {
  Asciiv3,     // HepMC3::WriterAscii — the only format implemented.
  HepMC2Ascii  // not implemented; constructing with this throws.
};

class HepMC3Writer {
 public:
  explicit HepMC3Writer(const std::string& filename,
                         HepMC3Format format = HepMC3Format::Asciiv3,
                         std::string generator_name = "LiPolGen",
                         std::string generator_version = LIPOLGEN_VERSION);
  ~HepMC3Writer();

  HepMC3Writer(const HepMC3Writer&) = delete;
  HepMC3Writer& operator=(const HepMC3Writer&) = delete;

  // Writes one event. The GenRunInfo (weight names, tool info) is finalized
  // from the first event written: "nominal" plus one "spin_weight_<k>" name
  // per entry of Event::spin_weights (1-based). Later events with a
  // different spin_weights size are still written correctly (the per-event
  // weight vector is set explicitly every time), but their weights will not
  // all be reachable by name through GenRunInfo::weight_index.
  void write(const Event& ev);

  // Flushes and closes the underlying file. Safe to call more than once;
  // the destructor also closes if this was not called.
  void close();

  // A10 (2026-11-03; docs/open_items/run_2026-11-03/DECISIONS.md A10, the
  // supervisor's choice under the author's delegation): the run's physics
  // provenance, written into the file's one GenRunInfo as string attributes
  // named "lipolgen.<key>" with the values exactly as given -- pass
  // `run_provenance(pipeline, ctx)` (pipeline.hpp), the map the npz `meta`
  // reads its provenance keys from, so a file and its npz say the same thing
  // (docs/HEPMC3_CONVENTION.md, "Run provenance attributes").  Must be called
  // before the first `write`, which is when the GenRunInfo is serialized;
  // throws otherwise.  Not calling it writes no such attribute at all.
  void set_run_attributes(const std::map<std::string, std::string>& attributes);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace lipolgen
