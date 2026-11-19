// SPDX-License-Identifier: GPL-3.0-or-later
// D7 (2026-10-27): the run's ONE unpolarised backend, enforced where a
// PythiaBridge meets a PipelineConfig -- see pythia_bridge.hpp, "D7".

#include <stdexcept>
#include <string>
#include <utility>

#include "lipolgen/pipeline.hpp"
#include "lipolgen/pythia_bridge.hpp"

namespace lipolgen {

namespace {

struct BridgeHook {
  std::shared_ptr<PythiaBridge> bridge;
  void operator()(Event& ev, Rng& rng) const { bridge->hadronize(ev, rng); }
};

}  // namespace

void check_bridge_species_backend(const PipelineConfig& cfg,
                                  const PythiaBridge& bridge) {
  if (bridge.options().f2_source == cfg.unpol_sf_obj) return;
  throw std::runtime_error(
      std::string("set_pythia_hadronizer: the bridge's species draw is on a "
                  "different unpolarised backend from this run's kernel -- "
                  "bridge f2_source is ") +
      (bridge.options().f2_source ? "an attached object"
                                  : "unset (the bridge's own ToyF2)") +
      " while the config's unpol_sf is " + unpol_sf_name(cfg.unpol_sf) +
      (cfg.unpol_sf_obj ? " with an object attached" : " with no object") +
      ".  The T2 struck-nucleon species draw is P(p) = Z F2p/(Z F2p + N F2n) "
      "on PythiaBridgeOptions::f2_source, so the run would draw its T2 "
      "species from one backend while meta[\"unpol_sf\"] recorded the other, "
      "and the toy's F2n/F2p is up to 24 % away from CT18NLO's.  Hand the "
      "bridge THE SAME OBJECT before constructing it: since D7 the shipped "
      "grid attaches lazily, so call attach_shipped_sf_backends(config) "
      "(Python _lipolgen.attach_shipped_sf_backends) and then set "
      "opts.f2_source = config.unpol_sf_obj -- what lipolgen-run and "
      "lipolgen.run do -- or name the toy (unpol_sf = toy) and leave both "
      "unset");
}

void set_pythia_hadronizer(PipelineConfig& cfg,
                           std::shared_ptr<PythiaBridge> bridge) {
  if (!bridge) {
    throw std::invalid_argument("set_pythia_hadronizer: bridge must not be null");
  }
  fill_grid_sf_backends(cfg);
  check_bridge_species_backend(cfg, *bridge);
  cfg.hadronizer = BridgeHook{std::move(bridge)};
}

}  // namespace lipolgen
