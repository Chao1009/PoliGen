// SPDX-License-Identifier: GPL-3.0-or-later
// MSTW2008 LO unpolarized structure functions over `Pythia8::MSTWpdf`
// (declared in include/lipolgen/mstw_sf.hpp).
//
// This is the SAME F2 construction as `src/lhapdf/lhapdf_sf.cpp`'s
// `LhapdfSF`, over a different grid.  The charge weights `kF2ProtonE2` /
// `kF2NeutronE2` and the body of `f2_from_weights` below are copied verbatim
// from `src/lhapdf/lhapdf_sf.cpp:40-70` -- same 5 flavours (d u s c b), same
// proton weights, same isospin u<->d swap for the neutron, same `0 < x < 1`
// guard (`PartonF2._f2p_scalar` returns 0.0 outside the open unit interval
// before touching the grid) and same `std::max(tot, 0.0)` clamp.  The ONLY
// difference is the grid accessor: LHAPDF's `PDF::xfxQ2(pid, x, q2)` (const)
// becomes PYTHIA's `PDF::xf(pid, x, Q2)` (non-const -- it memoizes the last
// (x, Q2) point in the object, which is why the reference below is not
// const).  Both return x*f(x, Q2) with pid on the PDG convention
// (1 = d, 2 = u, 3 = s, 4 = c, 5 = b, negative = antiquark), and for a
// proton beam (idBeam = 2212, beamType = 1) `PDF::xf` performs no isospin
// rearrangement at all -- it just hands back the stored xd/xu/xs/xc/xb and
// their antiquark partners.  So MSTW and CT18NLO are compared like with
// like, which is the whole point of this class.
//
// The one behavioural difference inherited from PYTHIA, stated rather than
// hidden: `PDF::xf` clamps EACH flavour at max(0, .) before we weight it,
// whereas LHAPDF's `xfxQ2` may return a small negative value in a
// low-x/low-Q2 corner and let the sum absorb it.  Both are then clamped
// again as a sum by `std::max(tot, 0.0)`.  Nothing in the DIS region this
// library samples is anywhere near that corner.
//
// LANDMINE (measured, gdb): sizeof(Pythia8::MSTWpdf) == 7 988 336 bytes --
// its c[13][64][48][5][5] grid is a by-value member -- so a stack-local
// instance segfaults during construction on a default 8 MB stack.  It is
// heap-allocated here, exactly as `LhapdfSF::Impl` heap-allocates its
// `LHAPDF::PDF`.

#include "lipolgen/mstw_sf.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>

#include "Pythia8/PartonDistributions.h"

#ifndef LIPOLGEN_PYTHIA8_PDFDATA
#define LIPOLGEN_PYTHIA8_PDFDATA ""
#endif

namespace lipolgen {

namespace {

// F2 (PartonF2._E2): five flavours d u s c b, proton e_q^2 weights.
constexpr std::array<std::pair<int, double>, 5> kF2ProtonE2 = {{
    {1, 1.0 / 9.0}, {2, 4.0 / 9.0}, {3, 1.0 / 9.0}, {4, 4.0 / 9.0}, {5, 1.0 / 9.0}}};
// F2 neutron (PartonF2.f2n's e2n dict): isospin u<->d swap of the charges above.
constexpr std::array<std::pair<int, double>, 5> kF2NeutronE2 = {{
    {1, 4.0 / 9.0}, {2, 1.0 / 9.0}, {3, 1.0 / 9.0}, {4, 4.0 / 9.0}, {5, 1.0 / 9.0}}};

double f2_from_weights(Pythia8::PDF& pdf, double x, double q2,
                       const std::array<std::pair<int, double>, 5>& weights) {
  // 0 < x < 1 guard: PartonF2._f2p_scalar / f2n's scalar both return 0.0
  // outside the open unit interval before touching the grid.
  if (!(x > 0.0 && x < 1.0)) return 0.0;
  double tot = 0.0;
  for (const auto& pe : weights) {
    tot += pe.second * (pdf.xf(pe.first, x, q2) + pdf.xf(-pe.first, x, q2));
  }
  return std::max(tot, 0.0);
}

// PYTHIA's own name for the file `i_fit` selects, for the error message --
// MSTWpdf::init reports a missing grid only through isSetup() and a line on
// stdout, so the throw below has to say which file was looked for.
const char* grid_file_name(int i_fit) {
  switch (i_fit) {
    case 1: return "mrstlostar.00.dat";
    case 2: return "mrstlostarstar.00.dat";
    case 3: return "mstw2008lo.00.dat";
    case 4: return "mstw2008nlo.00.dat";
    default: return "<no grid for this iFit>";
  }
}

}  // namespace

const std::string& pythia8_pdfdata_dir() {
  // Same env-beats-compiled-in rule as cluster.cpp's data_dir().
  static const std::string dir = [] {
    const char* env = std::getenv("LIPOLGEN_PYTHIA8_PDFDATA");
    if (env != nullptr && *env != '\0') return std::string(env);
    return std::string(LIPOLGEN_PYTHIA8_PDFDATA);
  }();
  return dir;
}

// -------------------------------------------------------------------- MstwSF

struct MstwSF::Impl {
  std::unique_ptr<Pythia8::MSTWpdf> pdf;  // 8 MB: NEVER a by-value member.
  int i_fit = 3;
  std::string path;
};

MstwSF::MstwSF(int i_fit_in, std::string pdfdata_path) : impl_(new Impl) {
  impl_->i_fit = i_fit_in;
  impl_->path = pdfdata_path.empty() ? pythia8_pdfdata_dir()
                                     : std::move(pdfdata_path);
  if (impl_->path.empty()) {
    throw std::runtime_error(
        "MstwSF: no PYTHIA 8 pdfdata directory (LIPOLGEN_PYTHIA8_PDFDATA is "
        "unset at build and at run time)");
  }
  // Heap, not stack -- see the file header.  MSTWpdf appends its own '/' if
  // the path lacks one, so no separator handling is needed here.
  impl_->pdf = std::make_unique<Pythia8::MSTWpdf>(2212, impl_->i_fit,
                                                  impl_->path);
  if (!impl_->pdf->isSetup()) {
    throw std::runtime_error("MstwSF: MSTWpdf failed to initialise from " +
                             impl_->path + " (iFit=" +
                             std::to_string(impl_->i_fit) + ", expected " +
                             grid_file_name(impl_->i_fit) + ")");
  }
}

MstwSF::~MstwSF() = default;

double MstwSF::f2p(double x, double q2) const {
  return f2_from_weights(*impl_->pdf, x, q2, kF2ProtonE2);
}

double MstwSF::f2n(double x, double q2) const {
  return f2_from_weights(*impl_->pdf, x, q2, kF2NeutronE2);
}

double MstwSF::f2n_over_f2p(double x) const {
  // Verbatim LhapdfSF::f2n_over_f2p: PartonF2.f2n_over_f2p(x, q2=10.0), the
  // Q2 default the Python signature carries (the abstract interface's
  // f2n_over_f2p is single-argument, exactly as ToyF2's is), so it is fixed
  // here, not threaded through.
  constexpr double kQ2 = 10.0;
  const double fp = f2p(x, kQ2);
  if (!(fp > 0.0)) return 1.0;  // PartonF2.f2n_over_f2p: np.where(f2p > 0, ..., 1.0)
  return f2n(x, kQ2) / std::max(fp, 1e-30);
}

int MstwSF::i_fit() const { return impl_->i_fit; }

const std::string& MstwSF::pdfdata_path() const { return impl_->path; }

}  // namespace lipolgen
