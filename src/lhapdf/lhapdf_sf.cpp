// SPDX-License-Identifier: GPL-3.0-or-later
// LHAPDF6 backends: `LhapdfSF` / `LhapdfG1` (declared in include/lipolgen/sf.hpp)
// and `Epps21Ratio` (declared in include/lipolgen/lhapdf_sf.hpp).
//
// Ports, formula for formula, of the LHAPDF-grid classes in
// PolarizedLithiumSim/fastsim/polli_fastsim/{structure,polarized}.py:
//   structure.PartonF2   -> LhapdfSF   (CT18NLO default; F2 = x Sum e_q^2 (q+qbar),
//                                       5 flavours d u s c b, neutron by u<->d swap)
//   polarized.PartonG1   -> LhapdfG1   (NNPDFpol11_100 default; g1 = (1/2) Sum e_q^2
//                                       (Dq+Dqbar), 3 flavours d u s -- Delta c =
//                                       Delta b = 0 is NNPDFpol1.1's own convention,
//                                       not an omission here)
// Nucleus combinations (F2A = Z F2p + N F2n, g1A = Z P_p g1p + N P_n g1n) are NOT
// duplicated here: wrap them in sf.hpp's `NuclearF2` / `PolSF::g1_nucleus`.
// SINCE D7 (2026-10-27) LhapdfSF("CT18NLO", 0) and LhapdfG1("NNPDFpol11_100", 0)
// are the SHIPPED run defaults (`UNPOL_SF_DEFAULT` / `POL_SF_DEFAULT`, pipeline.hpp),
// evaluated at LEADING ORDER -- the two formulas above, no coefficient functions,
// massless c and b -- which is the development item of DECISIONS.md D7.  Every
// grid this file opens is FROZEN at its own knot edges (`Grid::open`, at the end:
// LHAPDF's `nearest` extrapolator, the `logcubic` interpolator and ForcePositive 0
// set ON THE OBJECT, never read from an lhapdf.conf), and a missing set or
// lhapdf.conf is refused by `lhapdf_set_problem` before any LHAPDF Config exists.
// `Epps21Ratio` ports `money_delta_20260729.py`'s `NuclearF2FromGrid` (lhapdf_sf.hpp).

#include "lipolgen/lhapdf_sf.hpp"
#include "lipolgen/sf.hpp"
#include <algorithm>
#include <array>
#include <cstdlib>
#include <mutex>
#include <stdexcept>
#include <utility>
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
namespace lipolgen {

namespace {

bool ensure_quiet_impl() {  // D7: LHAPDF::setVerbosity(0) only BESIDE an lhapdf.conf
  return !LHAPDF::findFile("lhapdf.conf").empty() && (LHAPDF::setVerbosity(0), true); }
// F2 (PartonF2._E2): five flavours d u s c b, proton e_q^2 weights.
constexpr std::array<std::pair<int, double>, 5> kF2ProtonE2 = {{
    {1, 1.0 / 9.0}, {2, 4.0 / 9.0}, {3, 1.0 / 9.0}, {4, 4.0 / 9.0}, {5, 1.0 / 9.0}}};
// F2 neutron (PartonF2.f2n's e2n dict): isospin u<->d swap of the charges above.
constexpr std::array<std::pair<int, double>, 5> kF2NeutronE2 = {{
    {1, 4.0 / 9.0}, {2, 1.0 / 9.0}, {3, 1.0 / 9.0}, {4, 4.0 / 9.0}, {5, 1.0 / 9.0}}};
// g1 (PartonG1._E2): three flavours d u s only -- Delta c = Delta b = 0.
constexpr std::array<int, 3> kG1Pids = {1, 2, 3};

// proton e_q^2 for pid in {1 (d), 2 (u), 3 (s)}, matching PartonG1._E2.
double e2_dus(int pid) { return pid == 2 ? 4.0 / 9.0 : 1.0 / 9.0; }
// D7 (2026-10-27): one LHAPDF grid, FROZEN at its own knot edges; `open` is at the end.
struct Grid {
  std::unique_ptr<LHAPDF::PDF> pdf;
  double x_lo = 0.0, x_hi = 1.0, q2_lo = 0.0, q2_hi = 0.0;  // the knot array's own edges
  static Grid open(const std::string& set, int member);  // pre-flight, mkPDF, freeze
  double xfx(int pid, double x, double q2) const { return pdf->xfxQ2(pid, std::min(std::max(x, x_lo), x_hi), std::min(std::max(q2, q2_lo), q2_hi)); }
};
double f2_from_weights(const Grid& pdf, double x, double q2,
                       const std::array<std::pair<int, double>, 5>& weights) {
  // 0 < x < 1 guard: PartonF2._f2p_scalar / f2n's scalar both return 0.0
  // outside the open unit interval before touching the grid.
  if (!(x > 0.0 && x < 1.0)) return 0.0;
  double tot = 0.0;
  for (const auto& pe : weights) {
    tot += pe.second * (pdf.xfx(pe.first, x, q2) + pdf.xfx(-pe.first, x, q2));
  }
  return std::max(tot, 0.0);
}

double g1_scalar(const Grid& pdf, double x, double q2, bool swap_ud) {
  // PartonG1._g1_scalar: 0 < x < 1 guard, xfxQ2 already returns x*Delta-q,
  // so g1 = 0.5 * sum / x.
  if (!(x > 0.0 && x < 1.0)) return 0.0;
  double tot = 0.0;
  for (int pid : kG1Pids) {
    double e2 = e2_dus(pid);
    if (swap_ud && (pid == 1 || pid == 2)) e2 = e2_dus(3 - pid);  // isospin swap
    tot += e2 * (pdf.xfx(pid, x, q2) + pdf.xfx(-pid, x, q2));
  }
  return 0.5 * tot / x;
}

}  // namespace

void lhapdf_quiet() {  // D7: retried until an lhapdf.conf exists, never Config without one
  static std::mutex mu; static bool done = false; std::lock_guard<std::mutex> lock(mu);
  if (!done) done = ensure_quiet_impl();
}

// ------------------------------------------------------------------ LhapdfSF

struct LhapdfSF::Impl {
  Grid grid;  // D7: the set, frozen at its own knot edges (Grid::open)
};

LhapdfSF::LhapdfSF(std::string setname, int member) : impl_(new Impl) {
  // D7: lhapdf_set_problem() BEFORE lhapdf_quiet() and mkPDF, then frozen
  impl_->grid = Grid::open(setname, member);
}

LhapdfSF::~LhapdfSF() = default;

double LhapdfSF::q2_min() const { return impl_->grid.pdf->q2Min(); }

double LhapdfSF::f2p(double x, double q2) const {
  return f2_from_weights(impl_->grid, x, q2, kF2ProtonE2);
}

double LhapdfSF::f2n(double x, double q2) const {
  return f2_from_weights(impl_->grid, x, q2, kF2NeutronE2);
}

double LhapdfSF::f2n_over_f2p(double x) const {
  // PartonF2.f2n_over_f2p(x, q2=10.0): the Q2 default the Python signature
  // carries: q2 is optional there and unused by any call site in the
  // Python (there isn't one -- the interface's abstract f2n_over_f2p is
  // single-argument, exactly as ToyF2's is), so it is fixed here, not
  // threaded through.
  constexpr double kQ2 = 10.0;
  const double fp = f2p(x, kQ2);
  if (!(fp > 0.0)) return 1.0;  // PartonF2.f2n_over_f2p: np.where(f2p > 0, ..., 1.0)
  return f2n(x, kQ2) / std::max(fp, 1e-30);
}

// ------------------------------------------------------------------ LhapdfG1

struct LhapdfG1::Impl {
  Grid grid;  // D7: the set, frozen at its own knot edges (Grid::open)
};

LhapdfG1::LhapdfG1(std::string setname, int member) : impl_(new Impl) {
  // D7: lhapdf_set_problem() BEFORE lhapdf_quiet() and mkPDF, then frozen
  impl_->grid = Grid::open(setname, member);
}

LhapdfG1::~LhapdfG1() = default;

double LhapdfG1::g1p(double x, double q2) const {
  return g1_scalar(impl_->grid, x, q2, false);
}

double LhapdfG1::g1n(double x, double q2) const {
  return g1_scalar(impl_->grid, x, q2, true);
}

// -------------------------------------------------------------- Epps21Ratio

struct Epps21Ratio::Impl {
  Grid nuclear;  // D7: both grids frozen at their OWN knot edges (Grid::open)
  Grid proton;
};

Epps21Ratio::Epps21Ratio(std::string nuclear_set, int nuclear_member,
                         std::string proton_set, int proton_member)
    : impl_(new Impl) {
  // D7: each set through Grid::open -- lhapdf_set_problem() BEFORE
  // lhapdf_quiet() and mkPDF -- and then frozen at its own knot edges, so the
  // ratio below a grid's Q0 is f_A at the nuclear grid's edge over f_p at the
  impl_->nuclear = Grid::open(nuclear_set, nuclear_member);
  impl_->proton = Grid::open(proton_set, proton_member);  // proton grid's edge
}

Epps21Ratio::~Epps21Ratio() = default;

double Epps21Ratio::ratio(int pid, double x, double q2) const {
  const double fp = impl_->proton.xfx(pid, x, q2);
  if (fp == 0.0) return 0.0;
  const double fa = impl_->nuclear.xfx(pid, x, q2);
  return fa / fp;
}

double Epps21Ratio::f2_per_nucleon(double x, double q2) const {
  if (!(x > 0.0 && x < 1.0)) return 0.0;
  double tot = 0.0;
  for (const auto& pe : kF2ProtonE2) {
    const int pid = pe.first;
    const double e2 = pe.second;
    const double xfp_q = impl_->proton.xfx(pid, x, q2);
    const double xfp_qbar = impl_->proton.xfx(-pid, x, q2);
    const double r_q = ratio(pid, x, q2);
    const double r_qbar = ratio(-pid, x, q2);
    tot += e2 * (r_q * xfp_q + r_qbar * xfp_qbar);
  }
  return std::max(tot, 0.0);
}

double Epps21Ratio::f2a(double x, double q2, int A) const {
  return f2_per_nucleon(x, q2) * static_cast<double>(A);
}

// ------------------------------------------------- D7 (2026-10-27) additions
//
// Everything below this line was added by D7 and sits after every line the
// documents cite, so nothing above it moved.

namespace {

std::string escape_for(const std::string& set) {
  if (set == "CT18NLO")
    return "it is the shipped --unpol-sf default (and --b1-unpol ct18nlo): run "
           "the toy explicitly with --unpol-sf toy (Python unpol_sf='toy', C++ "
           "UnpolSfSource::Toy or lipolgen::toy_config()) and say so";
  if (set == "NNPDFpol11_100")
    return "it is the shipped --pol-sf default: run the toy explicitly with "
           "--pol-sf toy (Python pol_sf='toy', C++ PolSfSource::Toy or "
           "lipolgen::toy_config()) and say so";
  return "no shipped default reads it; it was asked for by name";
}

}  // namespace

std::string lhapdf_set_problem(const std::string& set, int member) {
  // ONLY LHAPDF::findFile / paths(), which read the environment and the
  // install prefix and never construct LHAPDF's Config -- so a missing
  // lhapdf.conf cannot reach Config::get(), whose function-static throws
  // again from its own destructor at exit (std::terminate, exit 134,
  // measured 2026-10-27).  The set's .info and member file are checked
  // first, so a set that is simply not installed is reported in LHAPDF's
  // own words, "Info file not found for PDF set '<name>'", which
  // validation/benchmarks/_report.py and the harnesses read as "blocked
  // (environment)" (`LHAPDF_MISSING`).
  const std::string conf = LHAPDF::findFile("lhapdf.conf");
  const std::string info = LHAPDF::findpdfsetinfopath(set);
  const std::string mem =
      info.empty() ? std::string() : LHAPDF::findpdfmempath(set, member);
  if (!conf.empty() && !info.empty() && !mem.empty()) return std::string();
  std::string searched;
  for (const std::string& p : LHAPDF::paths())
    searched += (searched.empty() ? "" : ", ") + p;
  const char* env = std::getenv("LHAPDF_DATA_PATH");
  const std::string what =
      info.empty()
          ? "Info file not found for PDF set '" + set + "'"
          : mem.empty() ? "Couldn't find a PDF data file for " + set + " #" +
                              std::to_string(member) + " (the set's info file is " +
                              info + ")"
                        : "LHAPDF's own lhapdf.conf was not found (the set itself is "
                          "at " + info + ")";
  // F7 (2026-10-30): the ADVICE follows what is missing.  With the set found
  // and only lhapdf.conf absent, "install the set" was the wrong cure.
  const std::string cure =
      (!info.empty() && !mem.empty())
          ? "Put LHAPDF's own lhapdf.conf (share/LHAPDF/lhapdf.conf of the "
            "LHAPDF install this build links) in a directory on "
            "LHAPDF_DATA_PATH -- the set itself is there"
          : "Install the set (.github/scripts/build_deps.sh fetches CT18NLO, "
            "NNPDFpol11_100 and EPPS21nlo_CT18Anlo_Li6, sha256-pinned) and point "
            "LHAPDF_DATA_PATH at it";
  return "LHAPDF cannot open the set '" + set + "': " + what +
         ".  LiPolGen's pre-flight (D7, 2026-10-27), run before any LHAPDF "
         "call: searched " + (searched.empty() ? "<nothing>" : searched) +
         "; LHAPDF_DATA_PATH=" + (env ? env : "<unset>") + "; lhapdf.conf " +
         (conf.empty() ? "NOT found" : "found at " + conf) + ".  " + cure +
         " -- or " + escape_for(set) +
         ".  It is never substituted by a toy in silence.";
}

namespace {

Grid Grid::open(const std::string& set, int member) {
  const std::string problem = lhapdf_set_problem(set, member);
  if (!problem.empty()) throw std::runtime_error(problem);
  lhapdf_quiet();
  Grid g;
  g.pdf.reset(LHAPDF::mkPDF(set, static_cast<std::size_t>(member)));
  auto* grid = dynamic_cast<LHAPDF::GridPDF*>(g.pdf.get());
  if (grid == nullptr) {
    throw std::runtime_error("LHAPDF set '" + set + "' is not a grid PDF "
                             "(lhagrid1), so it cannot be frozen at its grid "
                             "edges the way every set LiPolGen opens is (D7)");
  }
  // THE THREE SETTINGS NO lhapdf.conf MAY CHANGE (D7).  None of the three
  // installed sets' .info files carries them, so without this block they came
  // from whichever lhapdf.conf LHAPDF found first: measured 2026-10-27 by
  // changing nothing but that file, `Extrapolator: nearest` moved the 6Li
  // default cross section +5.03 % and the coherent one +7.02 % against
  // `continuation`, and `ForcePositive: 1` moved the 7Li A_par x1.488.
  if (grid->interpolator().type() != LHAPDF_INTERPOLATOR) {
    grid->setInterpolator(std::string(LHAPDF_INTERPOLATOR));
  }
  grid->setExtrapolator(std::string(LHAPDF_EXTRAPOLATOR));
  grid->setForcePositive(LHAPDF_FORCE_POSITIVE);
  // The edges `nearest` freezes at are the KNOT ARRAY's (what
  // NearestPointExtrapolator clamps to), read from the grid, not typed.
  // `xfx` clamps to them itself, which is bit-identical to LHAPDF's own
  // `nearest` -- that extrapolator interpolates at the same clamped knot --
  // without its per-call copy of the whole knot array (measured: a 5000-event
  // inclusive 6Li run took 3.7 s through it, 0.016 s through the clamp).
  const std::vector<double>& xs = grid->xKnots();
  const std::vector<double>& q2s = grid->q2Knots();
  g.x_lo = xs.front();
  g.x_hi = xs.back();
  g.q2_lo = q2s.front();
  g.q2_hi = q2s.back();
  // Fill LHAPDF's lazily cached members (flavour list, positivity mode) now,
  // single-threaded, so the object is read-only by the time it is shared.
  (void)grid->flavors();
  (void)grid->forcePositive();
  (void)g.xfx(2, 0.1, std::max(10.0, g.q2_lo));
  return g;
}

}  // namespace

double LhapdfG1::q2_min() const { return impl_->grid.pdf->q2Min(); }

}  // namespace lipolgen
