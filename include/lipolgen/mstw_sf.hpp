// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LIPOLGEN_MSTW_SF_HPP
#define LIPOLGEN_MSTW_SF_HPP

/// \file mstw_sf.hpp
/// MSTW2008 LO unpolarized structure functions, read from the grid PYTHIA 8
/// already ships.
///
/// WHY THIS EXISTS.  Every document in this tree used to say MSTW2008 LO "is
/// not installed" -- `docs/open_items/run_2026-09-02/phase_D_gate.md` (three
/// places), this repository's `b1_nuclear.hpp` header block, `docs/USAGE.md`
/// section 2a and `python/lipolgen/cli.py`'s run banner.  All of those were
/// corrected on 2026-09-03 when the gate was reread; the claim was true only
/// of LHAPDF's set store.  The MSTW2008 LO CENTRAL grid is
/// on disk as `<pythia8 datadir>/pdfdata/mstw2008lo.00.dat`, and PYTHIA ships
/// `Pythia8::MSTWpdf` to read it (`iFit = 3` selects exactly that member).
/// CDKS's b1 predictions were computed with MSTW2008 LO, so the CT18NLO
/// stand-in was the dominant remaining systematic of the 6Li b1 gate; this
/// backend removes the stand-in.
///
/// WHY IT IS NOT IN `sf.hpp`.  `sf.hpp` is on phase D's do-not-touch list,
/// and `lhapdf_sf.hpp` is the standing precedent for a backend-specific
/// extra that the core headers do not declare.  Unlike `LhapdfSF` (declared
/// in `sf.hpp`, implemented in `src/lhapdf/`), `MstwSF` is declared here and
/// implemented in `src/pythia/mstw_sf.cpp`; the core library still links
/// neither LHAPDF nor PYTHIA.
///
/// COMPARABILITY.  The charge weights and the F2 construction are copied
/// VERBATIM from `src/lhapdf/lhapdf_sf.cpp` (`kF2ProtonE2`, `kF2NeutronE2`,
/// `f2_from_weights`, including the `0 < x < 1` guard and the
/// `std::max(tot, 0.0)` clamp), so `MstwSF` and `LhapdfSF` differ ONLY in
/// which grid they read.  Any MSTW/CT18NLO ratio quoted from these two
/// classes is therefore a PDF comparison and not a convention comparison --
/// which is the entire point of the exercise.
///
/// LANDMINE, measured: `sizeof(Pythia8::MSTWpdf)` is 7 988 336 bytes,
/// because its `c[13][64][48][5][5]` interpolation grid is a by-value
/// member.  A stack-local `MSTWpdf` overflows an 8 MB stack DURING
/// CONSTRUCTION and segfaults.  The pimpl below keeps it on the heap
/// (`std::make_unique`), following `LhapdfSF::Impl`'s
/// `std::unique_ptr<LHAPDF::PDF>`.

#include <memory>
#include <string>

#include "lipolgen/sf.hpp"

namespace lipolgen {

/// The compiled-in PYTHIA 8 `pdfdata` directory (`LIPOLGEN_PYTHIA8_PDFDATA`,
/// set by CMake beside `LIPOLGEN_PYTHIA8_XMLDOC`), overridden at RUN time by
/// `$LIPOLGEN_PYTHIA8_PDFDATA` when that is set and non-empty -- the same
/// env-beats-compiled-in rule `cluster.cpp`'s `data_dir()` uses for the VMC
/// tables.  Empty only in a build where PYTHIA was not configured.
const std::string& pythia8_pdfdata_dir();

/// MSTW2008 LO F2p / F2n over `Pythia8::MSTWpdf` -- the same 5-flavour
/// (d u s c b) e_q^2 combination `LhapdfSF` uses, so the two are directly
/// comparable.  `i_fit` is PYTHIA's own fit selector: 1 = MRST LO*,
/// 2 = MRST LO**, 3 = MSTW 2008 LO (the default and the only member on
/// disk that CDKS used), 4 = MSTW 2008 NLO.  `pdfdata_path` empty means
/// `pythia8_pdfdata_dir()`.
///
/// Throws `std::runtime_error` if the grid file cannot be read (PYTHIA
/// reports that through `PDF::isSetup()`, not through an exception).
class MstwSF : public UnpolSF {
 public:
  explicit MstwSF(int i_fit = 3, std::string pdfdata_path = "");
  ~MstwSF() override;

  double f2p(double x, double q2) const override;
  double f2n(double x, double q2) const override;
  double f2n_over_f2p(double x) const override;

  /// The fit member actually loaded, and the directory it was read from.
  int i_fit() const;
  const std::string& pdfdata_path() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_MSTW_SF_HPP
