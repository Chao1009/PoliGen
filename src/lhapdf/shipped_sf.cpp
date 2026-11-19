// SPDX-License-Identifier: GPL-3.0-or-later
// D7 (2026-10-27): the SHIPPED structure-function defaults, built by the one
// tier that can build them and attached to the configs that name them.  The
// rules are in include/lipolgen/lhapdf_sf.hpp ("D7") and pipeline.hpp ("D7").

#include <mutex>
#include <stdexcept>
#include <string>

#include "lipolgen/lhapdf_sf.hpp"

namespace lipolgen {

namespace {

std::mutex& shipped_mutex() {
  static std::mutex m;
  return m;
}

}  // namespace

const char* lhapdf_set_name(UnpolSfSource source) {
  return source == UnpolSfSource::Ct18Nlo ? "CT18NLO" : nullptr;
}

const char* lhapdf_set_name(PolSfSource source) {
  return source == PolSfSource::NnpdfPol ? "NNPDFpol11_100" : nullptr;
}

std::shared_ptr<const UnpolSF> make_lhapdf_unpol_sf(UnpolSfSource source) {
  const char* set = lhapdf_set_name(source);
  if (set == nullptr) {
    throw std::invalid_argument(
        std::string("make_lhapdf_unpol_sf: unpol_sf = ") + unpol_sf_name(source) +
        " names no LHAPDF grid (ct18nlo is the only one)");
  }
  std::lock_guard<std::mutex> lock(shipped_mutex());
  static std::shared_ptr<const UnpolSF> ct18;
  if (!ct18) {
    ct18 = std::make_shared<const LhapdfSF>(set, 0);
    register_shipped_grid_sf(ct18.get());  // validate() accepts THIS object only
  }
  return ct18;
}

std::shared_ptr<const PolSF> make_lhapdf_pol_sf(PolSfSource source) {
  const char* set = lhapdf_set_name(source);
  if (set == nullptr) {
    throw std::invalid_argument(
        std::string("make_lhapdf_pol_sf: pol_sf = ") + pol_sf_name(source) +
        " names no LHAPDF grid (nnpdfpol is the only one)");
  }
  std::lock_guard<std::mutex> lock(shipped_mutex());
  static std::shared_ptr<const PolSF> nnpdfpol;
  if (!nnpdfpol) {
    nnpdfpol = std::make_shared<const LhapdfG1>(set, 0);
    register_shipped_grid_sf(nnpdfpol.get());  // validate() accepts THIS object only
  }
  return nnpdfpol;
}

void attach_shipped_sf_backends(PipelineConfig& cfg) {
  if (lhapdf_set_name(cfg.unpol_sf) != nullptr && !cfg.unpol_sf_obj) {
    cfg.unpol_sf_obj = make_lhapdf_unpol_sf(cfg.unpol_sf);
  }
  if (lhapdf_set_name(cfg.pol_sf) != nullptr && !cfg.pol_sf_obj) {
    cfg.pol_sf_obj = make_lhapdf_pol_sf(cfg.pol_sf);
  }
}

void install_lhapdf_sf_filler() { set_grid_sf_filler(&attach_shipped_sf_backends); }

}  // namespace lipolgen
