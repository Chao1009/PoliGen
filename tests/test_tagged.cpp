// SPDX-License-Identifier: GPL-3.0-or-later
// Tagged mode: the two-cluster spin (x) spectator model.
//
// Structural identities that must hold for any such model (normalization,
// isotropy sums, pure-wave CG limits), the analytic 7Li polarimetry and
// forward-limit predictions, the 6Li embedded-deuteron dilutions and S-wave
// reduction to the inclusive master formula, the Cosyn-Weiss deuteron-limit
// tensor gate, the boost against the fast simulation's own spectator sampler,
// and the reference tables of validation/reference/tagged.json.

#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "json_min.hpp"
#include "lipolgen/asymmetries.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/breakup.hpp"
#include "lipolgen/cluster_config.hpp"   // VMC_DEUTERON_WAVE
#include "lipolgen/constants.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/numerics.hpp"
#include "lipolgen/spin.hpp"
#include "lipolgen/b1_nuclear.hpp"
#include "lipolgen/tagged.hpp"
#include "lipolgen/xsec.hpp"

using namespace lipolgen;

#ifndef LIPOLGEN_REFERENCE_DIR
#define LIPOLGEN_REFERENCE_DIR "validation/reference"
#endif

namespace {

constexpr double kRtol = 1e-12;
/// `TaggedModel::norm` and `p2_moment_mixture` are themselves finite-grid
/// quadratures; validation/README.md asks for ~1e-6 on those two.
constexpr double kQuadRtol = 1e-9;

bool load_tagged(jsonmin::Value& out) {
  return jsonmin::load_file(std::string(LIPOLGEN_REFERENCE_DIR) + "/tagged.json",
                            out);
}

/// Is the reference blob on disk?  Called at doctest REGISTRATION time by the
/// `doctest::skip()` decorators below, so that a checkout without
/// `validation/reference/` reports these cases as SKIPPED in the tally
/// instead of PASSED WITH ZERO ASSERTIONS.  The pattern is
/// `tests/test_coherent.cpp`'s (phase F, 2026-09-05), which introduced it for
/// exactly this reason and was applied to test_coherent / test_cluster /
/// test_b1_nuclear and to no other file -- these twelve rtol-1e-12 port gates
/// kept the in-case `MESSAGE(...); return;` and went on passing, silently and
/// vacuously, until 2026-09-16.  A blob that IS there and does not parse is a
/// FAILURE, not a skip: that is what the `REQUIRE` on the loader inside each
/// case is for.
bool tagged_blob_present() {
  return std::ifstream(std::string(LIPOLGEN_REFERENCE_DIR) +
                       "/tagged.json")
      .good();
}


const TaggedModel& li6_model() {
  static const TaggedModel m(li6_alpha_channel());
  return m;
}
const TaggedModel& li7_model() {
  static const TaggedModel m(li7_alpha_channel());
  return m;
}
const TaggedModel& deut_model() {
  static const TaggedModel m(deuteron_channel());
  return m;
}
const TaggedModel& model_by_key(const std::string& key) {
  if (key == "li6_alpha") return li6_model();
  if (key == "li7_alpha") return li7_model();
  return deut_model();
}

/// `np.argmin(np.abs(v - target))`, first index on a tie -- NumPy's rule, and
/// the Cosyn-Weiss gate depends on it (the c grid is symmetric).
std::size_t argmin_abs(const std::vector<double>& v, double target) {
  std::size_t best = 0;
  double bd = std::fabs(v[0] - target);
  for (std::size_t i = 1; i < v.size(); ++i) {
    const double d = std::fabs(v[i] - target);
    if (d < bd) { bd = d; best = i; }
  }
  return best;
}
/// `np.argmax(np.abs(v))`, first index on a tie.
std::size_t argmax_abs(const std::vector<double>& v) {
  std::size_t best = 0;
  double bd = std::fabs(v[0]);
  for (std::size_t i = 1; i < v.size(); ++i) {
    if (std::fabs(v[i]) > bd) { bd = std::fabs(v[i]); best = i; }
  }
  return best;
}

/// Is the AV18 deuteron k-space table on disk?  The Cosyn-Weiss gate runs on
/// `data/vmc/deuteron/fdeut.av18` -- CW's OWN wave function -- which is a
/// different file from the `vmc/momenta/` overlaps `vmc_data_present` guards.
bool fdeut_av18_present() {
  std::ifstream f(data_path(VMC_DEUTERON_WAVE));
  return static_cast<bool>(f);
}

/// The k at which `radial_table(2)/radial_table(0)` first crosses `target`,
/// by linear interpolation between grid cells.  Restricted to
/// |ratio - target| < 1 so the blow-up at the S-wave node (k = 0.4130 GeV on
/// AV18, where f0 changes sign) is never mistaken for a crossing.  NaN if the
/// ratio never crosses.
double k_where_ratio(const TaggedModel& m, double target) {
  const std::vector<double>& f0 = m.radial_table(0);
  const std::vector<double>& f2 = m.radial_table(2);
  double prev_r = 0.0, prev_k = 0.0;
  bool have = false;
  for (std::size_t i = 0; i < m.nk(); ++i) {
    const double r = f2[i] / f0[i];
    if (!std::isfinite(r) || std::fabs(r - target) > 1.0) { have = false; continue; }
    if (have && (prev_r - target) * (r - target) <= 0.0 && prev_r != r) {
      const double t = (target - prev_r) / (r - prev_r);
      return prev_k + t * (m.k()[i] - prev_k);
    }
    prev_r = r;
    prev_k = m.k()[i];
    have = true;
  }
  return std::nan("");
}

double quantile(std::vector<double> v, double q) {
  std::sort(v.begin(), v.end());
  const double pos = q * static_cast<double>(v.size() - 1);
  const std::size_t i = static_cast<std::size_t>(pos);
  const double f = pos - static_cast<double>(i);
  return i + 1 < v.size() ? v[i] * (1.0 - f) + v[i + 1] * f : v.back();
}

}  // namespace

// --------------------------------------------------------- reference tables

TEST_CASE("tagged: channel construction against polligen" *
          doctest::skip(!tagged_blob_present())) {
  jsonmin::Value ref;
  REQUIRE(load_tagged(ref));
  CHECK_CLOSE(P_D_LI6, ref["P_D_LI6"].num(), kRtol);
  CHECK_CLOSE(P_D_DEUTERON, ref["P_D_DEUTERON"].num(), kRtol);
  for (const auto& kv : ref["channels"].obj()) {
    const TaggedChannel& ch = model_by_key(kv.first).channel();
    const jsonmin::Value& r = kv.second;
    CAPTURE(kv.first);
    CHECK_CLOSE(ch.j_ion, r["j_ion"].num(), kRtol);
    CHECK_CLOSE(ch.s_struck, r["s_struck"].num(), kRtol);
    CHECK_CLOSE(ch.s_spec, r["s_spec"].num(), kRtol);
    CHECK_CLOSE(ch.s_channel, r["s_channel"].num(), kRtol);
    CHECK(ch.label == r["label"].str());
    CHECK_CLOSE(ch.base.kappa(), r["base"]["kappa"].num(), kRtol);
    CHECK(ch.waves.size() == r["waves"].size());
    for (std::size_t i = 0; i < ch.waves.size(); ++i) {
      CHECK(ch.waves[i].l == static_cast<int>(r["waves"][i]["l"].num()));
      CHECK_CLOSE(ch.waves[i].prob, r["waves"][i]["prob"].num(), kRtol);
      CHECK_CLOSE(ch.waves[i].beta, r["waves"][i]["beta"].num(), kRtol);
    }
  }
}

TEST_CASE("tagged: the model grid and its tables against polligen" *
          doctest::skip(!tagged_blob_present())) {
  // 2026-09-06 RE-PIN.  `li6_alpha` and `deuteron` no longer read polligen's
  // numbers: polligen's `tagged._amp2_table` sums psi_L with no i^L, i.e. the
  // INVERTED S-D interference sign this file's `build_amp2` fix removed, so
  // those two `model` blocks are dumped from the FIXED C++ by
  // validation/repin_tagged_from_lipolgen.py and the file's own
  // `provenance` key reads "LiPolGen post-fix, formerly polligen".
  // `li7_alpha`'s block is still polligen's -- one L = 1 wave, the common i
  // is a global phase, and the fix moves its n_of_kc and p2_moment by
  // exactly zero (measured, bit for bit) -- so for THAT channel this is
  // still the port gate it always was.  Everything outside `model` (waves,
  // base, beam_configs, boost_spectator, P_D_*) is untouched polligen.
  // docs/benchmarking/07_cw_sign_investigation.md section 7.2;
  // docs/open_items/run_2026-09-06/phase_CW_numbers.md.
  jsonmin::Value ref;
  REQUIRE(load_tagged(ref));
  for (const auto& kv : ref["channels"].obj()) {
    const TaggedModel& model = model_by_key(kv.first);
    const jsonmin::Value& m = kv.second["model"];
    CAPTURE(kv.first);
    // the grid IS part of the model: n_of_kc(k, c) takes the nearest cell
    CHECK(model.nk() == static_cast<std::size_t>(m["grid"]["nk"].num()));
    CHECK(model.nc() == static_cast<std::size_t>(m["grid"]["nc"].num()));
    CHECK_CLOSE(model.k().front(), m["grid"]["k_first"].num(), kRtol);
    CHECK_CLOSE(model.k().back(), m["grid"]["k_last"].num(), kRtol);
    CHECK_CLOSE(model.c().front(), m["grid"]["c_first"].num(), kRtol);
    CHECK_CLOSE(model.c().back(), m["grid"]["c_last"].num(), kRtol);

    const std::vector<double> k_pts = m["k_pts"].flat();
    const std::vector<double> c_pts = m["c_pts"].flat();

    for (const jsonmin::Value& row : m["n_of_kc"].arr()) {
      const double mm = row["M"].num();
      CAPTURE(mm);
      for (std::size_t a = 0; a < k_pts.size(); ++a) {
        for (std::size_t b = 0; b < c_pts.size(); ++b) {
          CHECK_CLOSE(model.n_of_kc(mm, k_pts[a], c_pts[b]),
                      row["n_at_k_c"][a][b].num(), kRtol);
        }
      }
    }
    for (const jsonmin::Value& row : m["struck_populations"].arr()) {
      const double mm = row["M"].num();
      const std::vector<double> order = row["m_s_order"].flat();
      CAPTURE(mm);
      CHECK(order.size() == model.m_struck_values().size());
      for (std::size_t i = 0; i < order.size(); ++i) {
        CHECK_CLOSE(model.m_struck_values()[i], order[i], kRtol);
      }
      const std::vector<std::vector<double>> pops = model.struck_populations(mm);
      for (std::size_t i = 0; i < order.size(); ++i) {
        for (std::size_t a = 0; a < k_pts.size(); ++a) {
          for (std::size_t b = 0; b < c_pts.size(); ++b) {
            // the reference is sliced at the same nearest-cell lookup points
            const std::size_t ik = static_cast<std::size_t>(
                std::lower_bound(model.k().begin(), model.k().end(), k_pts[a])
                - model.k().begin()) - 1;
            const std::size_t ic = static_cast<std::size_t>(
                std::lower_bound(model.c().begin(), model.c().end(), c_pts[b])
                - model.c().begin()) - 1;
            CHECK_CLOSE(pops[i][ik * model.nc() + ic],
                        row["p_at_k_c"][i][a][b].num(), kRtol);
          }
        }
      }
    }
    for (const jsonmin::Value& row : m["population_integrated"].arr()) {
      const std::vector<double> want = row["p_m_s"].flat();
      const std::vector<double> got = model.population_integrated(row["M"].num());
      CAPTURE(row["M"].num());
      for (std::size_t i = 0; i < want.size(); ++i) {
        CHECK_CLOSE_AT(got[i], want[i], kRtol, 1e-300);
      }
    }
    for (const jsonmin::Value& row : m["norm"].arr()) {
      // a finite-grid quadrature: the README asks for a looser tolerance
      CHECK_CLOSE(model.norm(row["M"].num()), row["value"].num(), kQuadRtol);
    }
    for (const jsonmin::Value& row : m["p2_moment"].arr()) {
      CHECK_CLOSE(model.p2_moment(row["M"].num()), row["p2_moment"].num(), kRtol);
    }
    CHECK_CLOSE(model.vector_dilution(), m["vector_dilution"].num(), kRtol);
    if (m.has("tensor_dilution")) {
      CHECK_CLOSE(model.tensor_dilution(), m["tensor_dilution"].num(), kRtol);
    } else {
      // s_channel = 1/2: the tensor dilution is undefined BY DESIGN
      CHECK_THROWS(model.tensor_dilution());
    }
    const std::vector<double> ms_ion = m_values(model.channel().j_ion);
    const std::vector<double> uniform(ms_ion.size(), 1.0 / ms_ion.size());
    CHECK_CLOSE(model.p2_moment_mixture(uniform),
                m["p2_moment_mixture_uniform"].num(), kQuadRtol);
  }
}

TEST_CASE("tagged: boost_spectator against polligen (closed form, rtol 1e-12)" *
          doctest::skip(!tagged_blob_present())) {
  jsonmin::Value ref;
  REQUIRE(load_tagged(ref));
  for (const auto& kv : ref["channels"].obj()) {
    const TaggedChannel& ch = model_by_key(kv.first).channel();
    CAPTURE(kv.first);
    for (const jsonmin::Value& row : kv.second["boost_spectator"].arr()) {
      const SpectatorLab lab = boost_spectator(
          ch, row["k"].num(), row["c"].num(), row["phi_k"].num(),
          row["ion_momentum_per_nucleon"].num(), row["theta_s"].num(),
          row["phi_s"].num());
      CHECK_CLOSE_AT(lab.pT, row["pT"].num(), kRtol, 1e-300);
      CHECK_CLOSE_AT(lab.theta, row["theta"].num(), kRtol, 1e-300);
      CHECK_CLOSE(lab.p_lab, row["p_lab"].num(), kRtol);
      CHECK_CLOSE(lab.R, row["R"].num(), kRtol);
      CHECK_CLOSE(lab.xL, row["xL"].num(), kRtol);
      CHECK_CLOSE_AT(lab.kx, row["kx"].num(), kRtol, 1e-300);
      CHECK_CLOSE_AT(lab.ky, row["ky"].num(), kRtol, 1e-300);
      CHECK_CLOSE_AT(lab.kz, row["kz"].num(), kRtol, 1e-300);
      CHECK_CLOSE_AT(lab.phi_spec, row["phi_spec"].num(), kRtol, 1e-300);
    }
  }
}

// ---------------------------------------------------- structural identities

TEST_CASE("tagged: every channel is normalized") {
  for (const TaggedModel* model : {&li6_model(), &li7_model(), &deut_model()}) {
    for (double m : m_values(model->channel().j_ion)) {
      CAPTURE(model->channel().label);
      CAPTURE(m);
      CHECK_CLOSE_AT(model->norm(m), 1.0, 0.0, 2e-3);
    }
  }
}

TEST_CASE("tagged: sum_M n_M is isotropic about the axis") {
  // unpolarized-beam isotropy: the L interference cancels in the M sum
  for (const TaggedModel* model : {&li6_model(), &li7_model(), &deut_model()}) {
    const std::size_t nk = model->nk(), nc = model->nc();
    std::vector<double> tot(nk * nc, 0.0);
    for (double m : m_values(model->channel().j_ion)) {
      const std::vector<double>& n = model->n_of_kc(m);
      for (std::size_t j = 0; j < tot.size(); ++j) tot[j] += n[j];
    }
    double worst = 0.0;
    for (std::size_t ik = 0; ik < nk; ++ik) {
      double lo = tot[ik * nc], hi = tot[ik * nc], sum = 0.0;
      for (std::size_t ic = 0; ic < nc; ++ic) {
        const double v = tot[ik * nc + ic];
        lo = std::fmin(lo, v);
        hi = std::fmax(hi, v);
        sum += v;
      }
      worst = std::fmax(worst, (hi - lo) / std::fmax(sum / nc, 1e-300));
    }
    CAPTURE(model->channel().label);
    CHECK(worst < 1e-10);
  }
}

TEST_CASE("tagged: a pure S wave is m-independent and diagonal in the spin") {
  const TaggedModel model(li6_alpha_channel(BETA_DEFAULT, 0.0));
  const std::vector<double>& n1 = model.n_of_kc(1.0);
  for (double m : {0.0, -1.0}) {
    const std::vector<double>& n = model.n_of_kc(m);
    for (std::size_t j = 0; j < n.size(); ++j) CHECK_CLOSE(n[j], n1[j], kRtol);
  }
  // S wave: the struck deuteron's spin IS the ion's, exactly
  const std::vector<double> p = model.population_integrated(1.0);
  CHECK_CLOSE_AT(p[0], 1.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(p[1], 0.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(p[2], 0.0, 0.0, 1e-12);
  // and with P_D = 0 the alpha-d density is m-independent, so A_zz^tag
  // vanishes identically -- the D wave IS the tensor observable
  const std::size_t ic = argmin_abs(model.c(), 0.0);
  const std::vector<double> a = azz_tensor_curve(model, ic);
  for (double v : a) CHECK_CLOSE_AT(v, 0.0, 0.0, 1e-12);
}

TEST_CASE("tagged: the pure D-wave forward limit is CG-exact") {
  // theta_k = 0: only m_L = 0 survives, so
  // n_{+-1}/n_0 = CG(2 0 1 1|1 1)^2 / CG(2 0 1 0|1 0)^2 = (1/10)/(4/10) = 1/4,
  // and the thirds combination gives A_zz^wf(0) = -1 exactly.
  TaggedChannel ch = li6_alpha_channel();
  ch.waves = {Wave{2, 1.0, BETA_DEFAULT}};
  const TaggedModel model(ch, 1.2, 280, 4001);
  const std::size_t ic = argmax_abs(model.c());
  const double expected = std::pow(clebsch_gordan(2, 0, 1, 1, 1, 1), 2)
                          / std::pow(clebsch_gordan(2, 0, 1, 0, 1, 0), 2);
  CHECK_CLOSE(expected, 0.25, 1e-12);
  const std::vector<double>& n1 = model.n_of_kc(1.0);
  const std::vector<double>& n0 = model.n_of_kc(0.0);
  const std::vector<double> a = azz_tensor_curve(model, ic);
  for (std::size_t ik = 10; ik < model.nk(); ik += 37) {
    CHECK_CLOSE(n1[ik * model.nc() + ic] / n0[ik * model.nc() + ic], expected, 5e-3);
    CHECK_CLOSE_AT(a[ik], -1.0, 0.0, 3e-3);
  }
}

// ------------------------------------------ 7Li polarimetry + forward limit

TEST_CASE("tagged: 7Li stretched state is |Y_1^1|^2 ~ sin^2 theta_k") {
  const TaggedModel& m = li7_model();
  const std::vector<double>& n = m.n_of_kc(1.5);
  const std::size_t ik = 50;
  double peak = 0.0, epeak = 0.0;
  for (std::size_t ic = 0; ic < m.nc(); ++ic) {
    peak = std::fmax(peak, n[ik * m.nc() + ic]);
    epeak = std::fmax(epeak, 1.0 - m.c()[ic] * m.c()[ic]);
  }
  for (std::size_t ic = 0; ic < m.nc(); ++ic) {
    CHECK_CLOSE_AT(n[ik * m.nc() + ic] / peak,
                   (1.0 - m.c()[ic] * m.c()[ic]) / epeak, 0.0, 1e-12);
  }
}

TEST_CASE("tagged: 7Li P2 moments and the <P2> = -T/5 polarimeter") {
  const TaggedModel& m = li7_model();
  const double tol = 3e-4;   // midpoint-rule c grid, (dc)^2 ~ 1e-4
  CHECK_CLOSE_AT(m.p2_moment(1.5), -0.2, 0.0, tol);
  CHECK_CLOSE_AT(m.p2_moment(0.5), +0.2, 0.0, tol);
  CHECK_CLOSE_AT(m.p2_moment(-0.5), +0.2, 0.0, tol);
  CHECK_CLOSE_AT(m.p2_moment(-1.5), -0.2, 0.0, tol);
  // the in-situ alignment polarimeter: <P2(cos theta_k)> = -T/5 for ANY fill
  const std::array<double, 4> pops = spin32_populations(0.5, 0.4, 0.1);
  const std::vector<double> pv(pops.begin(), pops.end());
  const AxisMoments mom = moments_along_axis(1.5, pv);
  CHECK_CLOSE_AT(m.p2_moment_mixture(pv), -mom.tensor / 5.0, 0.0, tol);
  // and on a second, unrelated fill
  const std::array<double, 4> p2 = spin32_populations(-0.3, -0.2, 0.0);
  const std::vector<double> pv2(p2.begin(), p2.end());
  CHECK_CLOSE_AT(m.p2_moment_mixture(pv2),
                 -moments_along_axis(1.5, pv2).tensor / 5.0, 0.0, tol);
}

TEST_CASE("tagged: 7Li triton polarization and the forward-limit gate") {
  const TaggedModel& m = li7_model();
  const std::vector<double> p32 = m.population_integrated(1.5);
  CHECK_CLOSE_AT(p32[0], 1.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(p32[1], 0.0, 0.0, 1e-12);
  const std::vector<double> p12 = m.population_integrated(0.5);
  CHECK_CLOSE_AT(p12[0], 2.0 / 3.0, 0.0, 3e-4);
  CHECK_CLOSE_AT(p12[1], 1.0 / 3.0, 0.0, 3e-4);
  CHECK_CLOSE_AT(p12[0] - p12[1], 1.0 / 3.0, 0.0, 6e-4);
  // khat-resolved: P_t(theta; M = 1/2) = (5 c^2 - 1)/(3 c^2 + 1)
  const std::vector<std::vector<double>> pop = m.struck_populations(0.5);
  for (std::size_t ic = 0; ic < m.nc(); ++ic) {
    const double c = m.c()[ic];
    const double got = pop[0][50 * m.nc() + ic] - pop[1][50 * m.nc() + ic];
    CHECK_CLOSE_AT(got, (5 * c * c - 1.0) / (3 * c * c + 1.0), 0.0, 1e-10);
  }
  // the plans/05 forward-limit gate: a stretched fill leaves the triton fully
  // polarized, so P_p(7Li) is the triton's own effective proton polarization
  const double p_p = 1.0 * TRITON().eff_pol_p;
  CHECK(std::fabs(p_p - 0.866) < 0.02);
  // The NEUTRON half is an OPEN gate.  The model's own value is 2 N P_n of the
  // per-nucleon triton slot; the plans/05 target is -0.037 and this model does
  // not reach it.  It is reported, not fudged.
  const double p_n = TRITON().eff_pol_n;
  CHECK_CLOSE(p_n, -0.028, 1e-12);
  MESSAGE("7Li forward limit: P_p = " << p_p << " (gate 0.866 +- 0.02, met); "
          "P_n = " << p_n << " against the plans/05 gate -0.037 -- OPEN, the "
          "model's own value is asserted and the gate is not");
}

TEST_CASE("tagged: the A = 3 slots are BISSEY's convention, with the explicit "
          "factor 2 on the proton") {
  // [Bissey02] (Bissey, Guzey, Strikman, Thomas, hep-ph/0109069, PRC 65
  // (2002) 064317) Eq. (2), quoted verbatim off the source PDF (re-fetched
  // 2026-09-16):
  //
  //     g1He(x, Q2) = Pn g1n(x, Q2) + 2 Pp g1p(x, Q2)
  //
  //   with "Pn = 0.86 +- 0.02 and Pp = -0.028 +- 0.004 [10]".
  //
  // THE FACTOR 2 IS THE WHOLE POINT.  It says P_p is PER PROTON, so the
  // whole-nucleus proton sum is 2 P_p = -0.056 -- which is what `HE3()` and
  // `TRITON()` store and what `PolSF::g1_nucleus`'s Z * eff_pol_p then
  // reproduces.  `docs/PHYSICS_CHANNELS.md` recorded this as "an unresolved
  // convention conflict" with `docs/open_items/physics_literature.md:95`
  // until 2026-09-16; the paper resolves it, the literature note was the site
  // in error, and only `TRITON().eff_pol_n == -0.028` was pinned, so nothing
  // guarded the COMBINATION a reader following that note would have halved.
  //
  // Halving them flips the sign of g1(3He): on these backends the shipped
  // code gives -2.74450912e-03 and the whole-nucleus misreading gives
  // +8.31711655e-04.
  const ToyF2 f2;
  const ToyG1 g(std::make_shared<ToyF2>(f2));
  const double x = 0.3, q2 = 5.0;
  const double g1p = g.g1p(x, q2), g1n = g.g1n(x, q2);

  CHECK(HE3().Z == 2);
  CHECK(HE3().N() == 1);
  CHECK(HE3().eff_pol_p == -0.028);
  CHECK(HE3().eff_pol_n == 0.86);
  CHECK_CLOSE(g.g1_nucleus(HE3(), x, q2),
              0.86 * g1n + 2.0 * (-0.028) * g1p, 1e-12);

  // The mirror nucleus, the same equation with p <-> n.
  CHECK(TRITON().Z == 1);
  CHECK(TRITON().N() == 2);
  CHECK_CLOSE(g.g1_nucleus(TRITON(), x, q2),
              0.86 * g1p + 2.0 * (-0.028) * g1n, 1e-12);

  // And the misreading is FAR away, not a rounding: opposite sign on 3He.
  const double whole_nucleus = 0.86 * g1n + (-0.028) * g1p;
  CHECK(g.g1_nucleus(HE3(), x, q2) < 0.0);
  CHECK(whole_nucleus > 0.0);
}

// -------------------------------------- 6Li dilutions + inclusive reduction

TEST_CASE("tagged: 6Li embedded-deuteron dilutions") {
  const TaggedModel& m = li6_model();
  const double v = m.vector_dilution();
  CHECK_CLOSE_AT(v, 1.0 - 1.5 * P_D_LI6, 0.0, 1e-3);
  CHECK_CLOSE_AT(v, 0.87, 0.0, 2e-3);
  const double t = m.tensor_dilution();
  CHECK(t > 0.7);
  CHECK(t < 1.0);
  // plans/08 D9: the RANK-2 transfer constant is the model's own tensor
  // dilution, and the closed form 1 - (9/10) P_D
  CHECK_CLOSE_AT(LI6_B1_RANK2_TRANSFER, t, 0.0, 1e-4);
  CHECK_CLOSE_AT(t, 1.0 - 0.9 * P_D_LI6, 0.0, 1e-4);
  CHECK_CLOSE_AT(LI6_B1_LEGACY_TRANSFER, v, 0.0, 2e-3);
  // deuteron control: the standard 1 - (3/2) w_D
  CHECK_CLOSE_AT(deut_model().vector_dilution(), 1.0 - 1.5 * P_D_DEUTERON,
                 0.0, 1e-3);
}

TEST_CASE("tagged: the 6Li S-wave limit reduces to the inclusive deuteron") {
  // P_D = 0: p(m_d|M) = delta_{m_d,M}, so the k-integrated tagged tensor
  // asymmetry is the inclusive polarized-deuteron azz, bin by bin.
  const TaggedModel model(li6_alpha_channel(BETA_DEFAULT, 0.0));
  InclusiveKernel::Options opt;
  opt.b1_func = [](double x, double q2, double f1) { return toy_b1(x, q2, f1); };
  const InclusiveKernel kern(DEUTERON(), opt);
  const BeamConfig cfg = default_configs("6Li")[1];
  const double s = cfg.s_per_nucleon();
  const double xs[3] = {0.05, 0.15, 0.35};
  const double q2s[3] = {5.0, 12.0, 30.0};
  for (int i = 0; i < 3; ++i) {
    const double x = xs[i], q2 = q2s[i], y = q2 / (s * x);
    const SFTables t = kern.tables(x, q2);
    double sig[3];
    const double ms_ion[3] = {1.0, 0.0, -1.0};
    for (int a = 0; a < 3; ++a) {
      const std::vector<double> p = model.population_integrated(ms_ion[a]);
      double w_mix = 0.0;
      for (int b = 0; b < 3; ++b) {
        if (p[b] <= 0.0) continue;
        EventSpinState st;
        st.lam_e = 0;
        st.pe = 0.0;
        st.j = 1.0;
        st.m = ms_ion[b];
        w_mix += p[b] * (1.0 + kern.amplitudes(t, x, q2, s, st).w_avg);
      }
      sig[a] = w_mix;
    }
    const double measured = (sig[0] + sig[2] - 2 * sig[1])
                            / (sig[0] + sig[2] + sig[1]);
    CHECK_CLOSE(measured, azz(t.b1, t.f1, t.f2, x, y), 1e-10);
  }
}

// ------------------------------------------------ deuteron control (CW gate)

TEST_CASE("tagged: the deuteron S/D interference shape") {
  // RE-PINNED 2026-09-06 with the `build_amp2` i^L fix
  // (docs/benchmarking/07_cw_sign_investigation.md).  The OLD gate asserted
  // an interior peak, `0.2 < k[argmax|A_zz|] < 0.45`, and it was an artefact
  // of the inverted S-D sign: with psi_2 summed in place of phi_2 the curve
  // was CW Eq. (6.12) at MINUS f2/f0, which turns over at f2/f0 = 1/sqrt(2)
  // (k = 0.3098 GeV on this channel).  With the phase applied the Hulthen
  // pair's |A_zz| is MONOTONE in k across the whole grid -- measured: 279 of
  // 279 steps non-decreasing, zero decreasing, worst drop 0.000e+00 -- and
  // saturates at CW's +1 as f2/f0 -> 1.2866, so the maximum sits at the top
  // cell k = 1.2 and there is no interior peak to pin.  The threshold
  // suppression and the O(1) CW window are unchanged statements and survive.
  const TaggedModel& m = deut_model();
  const std::size_t ic = argmin_abs(m.c(), 0.0);
  const std::vector<double> a = azz_tensor_curve(m, ic);
  double thresh = 0.0, best = 0.0;
  std::size_t jbest = 0, ndec = 0;
  double worst_drop = 0.0;
  for (std::size_t ik = 0; ik < m.nk(); ++ik) {
    if (m.k()[ik] < 0.02) thresh = std::fmax(thresh, std::fabs(a[ik]));
    if (std::fabs(a[ik]) > best) { best = std::fabs(a[ik]); jbest = ik; }
    if (ik > 0 && std::fabs(a[ik]) < std::fabs(a[ik - 1])) {
      ++ndec;
      worst_drop = std::fmax(worst_drop, std::fabs(a[ik - 1]) - std::fabs(a[ik]));
    }
  }
  CHECK(thresh < 0.02);                      // D-wave threshold suppression
                                             // MEASURED 0.006418
  double window = 0.0, window_lo = 1e300;
  for (std::size_t ik = 0; ik < m.nk(); ++ik) {
    if (m.k()[ik] > 0.25 && m.k()[ik] < 0.5) {
      window = std::fmax(window, std::fabs(a[ik]));
      window_lo = std::fmin(window_lo, std::fabs(a[ik]));
    }
  }
  CHECK(window > 0.45);                      // O(1) in the CW window
  CHECK_CLOSE_AT(window, 0.957798, 0.0, 1e-5);      // MEASURED
  CHECK_CLOSE_AT(window_lo, 0.731102, 0.0, 1e-5);   // MEASURED
  // monotone: the whole point of the corrected sign on this channel
  CHECK(ndec == 0);
  CHECK(worst_drop == 0.0);
  CHECK(jbest == m.nk() - 1);
  CHECK_CLOSE_AT(m.k()[jbest], 1.2, 0.0, 1e-12);
  CHECK_CLOSE_AT(a[jbest], +0.996607, 0.0, 1e-5);   // MEASURED, -> CW's +1
  CHECK(a[jbest] < 1.0 + 1e-12);             // CW's stated ceiling
  // integrated over all angles the wf tensor asymmetry vanishes: the
  // observable lives in the angular structure (midpoint residual ~ dc^2)
  for (double mm : {1.0, 0.0}) {
    for (std::size_t ik = 0; ik < m.nk(); ik += 29) {
      double s0 = 0.0, s1 = 0.0;
      for (std::size_t ii = 0; ii < m.nc(); ++ii) {
        s0 += m.n_of_kc(mm)[ik * m.nc() + ii];
        s1 += m.n_of_kc(1.0)[ik * m.nc() + ii];
      }
      CHECK_CLOSE(s0, s1, 3e-4);
    }
  }
}

TEST_CASE("tagged: the struck-neutron polarization of the triplet") {
  const TaggedModel& m = deut_model();
  const double cases[3][2] = {{1.0, +1.0}, {0.0, 0.0}, {-1.0, -1.0}};
  for (const auto& c : cases) {
    const std::vector<std::vector<double>> pair = m.pair_populations(c[0]);
    const std::vector<double> mn = m_values(0.5);
    double pol = 0.0;
    for (std::size_t i = 0; i < pair.size(); ++i) {
      double marg = 0.0;
      for (double v : pair[i]) marg += v;
      pol += 2.0 * marg * mn[i];
    }
    CHECK_CLOSE_AT(pol, c[1], 0.0, 1e-12);
  }
}

TEST_CASE("tagged: the Cosyn-Weiss deuteron tensor gate (CW TABLE II)"
          * doctest::skip(!fdeut_av18_present())) {
  REQUIRE(fdeut_av18_present());   // the decorator already guaranteed it
  // Cosyn-Weiss II (arXiv:2603.23700) p. 35, Eqs. (6.11)-(6.14) and TABLE II.
  // REWRITTEN 2026-09-06 -- docs/benchmarking/07_cw_sign_investigation.md
  // section 8, alongside the `build_amp2` i^L fix.
  //
  // THE MAPPING IS +1, NOT -2.  A_T|| = sqrt(2/3) P_[T_LL,U]/P_[U,U], and with
  // CW Eq. (2.30b) {T_LL,T_LT,T_TT} = W(Lambda) x {1/3,0,0}, W = (1,-2,1), the
  // (+1,+1,-2) combination of P_[U](T_D) gives sqrt(6)/3 = sqrt(2/3) times the
  // same ratio.  A_T|| IS our A_zz^wf.  The (1 - 3cos^2 theta_k) = -2 P2 factor
  // is INSIDE A_zz^wf already; CW's "-2" is the value at theta_k = 0, where the
  // Lambda = +-1 densities have a node (Eq. 6.13), not a conversion factor.
  // The old gate applied it a second time and so double-counted.
  //
  // TABLE II is quoted for the AV18 radial wave functions, so the CW rows run
  // on the AV18 control.  The Hulthen pair carries only the two statements that
  // do not depend on which wave function it is -- the old gate ran CW's k
  // landmarks on Hulthen, whose f2/f0 never reaches sqrt(2) anywhere on the
  // grid (it tops out at 1.2866), and its "peak at 0.3098 GeV vs CW's 0.30"
  // was a coincidence at f2/f0 = 0.7053, i.e. CW's Eq. (6.14) MINIMUM read
  // through a flipped f2 and called Eq. (6.13)'s maximum.
  const TaggedModel v(deuteron_channel(BETA_DEFAULT, P_D_DEUTERON,
                                       ClusterWaveSource::VmcAV18));
  const TaggedModel& h = deut_model();

  // ---- (a) the mapping, as an IDENTITY on every cell of the grid.
  // A_zz^wf(k,c) == [(2 f0 + f2/sqrt2)(f2/sqrt2)/(f0^2+f2^2)] (1 - 3 cos^2),
  // with f_L the model's own normalized radial tables.  Machine precision:
  // this is the whole of Eq. (6.12), not a shape comparison.
  for (const TaggedModel* m : {&v, &h}) {
    double worst = 0.0;
    const std::vector<double>& n1 = m->n_of_kc(1.0);
    const std::vector<double>& n0 = m->n_of_kc(0.0);
    const std::vector<double>& nm = m->n_of_kc(-1.0);
    for (std::size_t ik = 0; ik < m->nk(); ++ik) {
      const double f0 = m->radial_table(0)[ik], f2 = m->radial_table(2)[ik];
      const double r = f2 / f0;
      if (!std::isfinite(r)) continue;
      const double q = (2.0 + r / std::sqrt(2.0)) * (r / std::sqrt(2.0))
                       / (1.0 + r * r);
      for (std::size_t ic = 0; ic < m->nc(); ++ic) {
        const std::size_t j = ik * m->nc() + ic;
        const double a = (n1[j] + nm[j] - 2.0 * n0[j]) / (n1[j] + nm[j] + n0[j]);
        const double c = m->c()[ic];
        worst = std::fmax(worst, std::fabs(a - q * (1.0 - 3.0 * c * c)));
      }
    }
    CHECK(worst < 1e-12);          // MEASURED 8.882e-16 on BOTH channels
                                   // (2.740e+00 before the fix)
  }

  // ---- (b) the P2 factorization is exact: A_zz^wf/P2 is a function of k alone.
  // (The old gate's (a); it is a real property and it survives untouched, only
  // moved onto AV18 and re-pinned to the corrected value.)
  std::vector<double> p2(v.nc());
  for (std::size_t i = 0; i < v.nc(); ++i) p2[i] = 0.5 * (3.0 * v.c()[i] * v.c()[i] - 1.0);
  const std::size_t ik30 = argmin_abs(v.k(), 0.30);
  CHECK_CLOSE_AT(v.k()[ik30], 0.3012, 0.0, 5e-4);
  double rmin = 1e300, rmax = -1e300;
  for (std::size_t ic = 0; ic < v.nc(); ++ic) {
    if (std::fabs(p2[ic]) <= 1e-3) continue;
    const double r = azz_tensor_curve(v, ic)[ik30] / p2[ic];
    rmin = std::fmin(rmin, r); rmax = std::fmax(rmax, r);
  }
  CHECK(rmax - rmin < 1e-5);                                  // MEASURED 1.2e-14
  CHECK_CLOSE_AT(0.5 * (rmin + rmax), -1.99928, 0.0, 1e-4);   // MEASURED -1.999276
                                                              // (+0.636821 before)

  // ---- (c) CW's OWN k landmarks, on CW's own wave function.
  // "Eq. (6.13) [f2/f0 = +sqrt2] is satisfied for k = 0.30 GeV"
  // "Eq. (6.14) [f2/f0 = -1/sqrt2] is satisfied only at k ~ 1 GeV"
  //   MEASURED on the model grid: 0.298121 GeV and 1.034872 GeV
  //   (raw fdeut.av18 k-block: 0.2984 GeV and 1.0257 GeV)
  CHECK_CLOSE_AT(k_where_ratio(v, +std::sqrt(2.0)), 0.30, 0.0, 0.01);
  CHECK_CLOSE_AT(k_where_ratio(v, -1.0 / std::sqrt(2.0)), 1.00, 0.0, 0.05);

  // ---- (d) TABLE II's two extremes, through the P2 factorization, with the
  // +1 mapping and NO extra factor anywhere.
  const std::size_t ic0 = argmax_abs(v.c());          // |c| = 0.989583
  const std::vector<double> env = [&] {
    std::vector<double> e(v.nk());
    const std::vector<double> curve = azz_tensor_curve(v, ic0);
    for (std::size_t i = 0; i < v.nk(); ++i) e[i] = curve[i] / p2[ic0];
    return e;
  }();
  std::size_t jlo = 0, jhi = 0;
  for (std::size_t i = 1; i < v.nk(); ++i) {
    if (env[i] < env[jlo]) jlo = i;
    if (env[i] > env[jhi]) jhi = i;
  }
  CHECK_CLOSE_AT(env[jlo], -2.0, 0.0, 1e-3);          // MEASURED -1.999864, CW: -2
  CHECK_CLOSE_AT(v.k()[jlo], 0.2968, 0.0, 5e-4);      //  at f2/f0 = +1.3942 (-> sqrt2)
  CHECK_CLOSE_AT(env[jhi], +1.0, 0.0, 1e-3);          // MEASURED +0.999997, CW: +1
  CHECK_CLOSE_AT(v.k()[jhi], 1.0366, 0.0, 5e-4);      //  at f2/f0 = -0.7056 (-> -1/sqrt2)
  CHECK_CLOSE_AT(env[jlo] * (-0.5), +1.0, 0.0, 1e-3); // the theta = pi/2 row of TABLE II

  // ---- (e) TABLE II row by row, at the actual cell centres.
  const std::size_t ic90 = argmin_abs(v.c(), 0.0);    // |c| = 0.010417
  CHECK_CLOSE_AT(azz_tensor_curve(v, ic0)[ik30],  -1.93712, 0.0, 1e-4);  // CW -2 at the cell centre; measured -1.9371243623, residual 4.4e-6 (x23 slack)
  CHECK_CLOSE_AT(azz_tensor_curve(v, ic90)[ik30], +0.99931, 0.0, 1e-4);  // CW +1 at the cell centre; measured +0.9993127695, residual 2.8e-6 (x36 slack; 2e-3 would have been x722)
  const std::size_t ik100 = argmin_abs(v.k(), 1.00);
  CHECK_CLOSE_AT(azz_tensor_curve(v, ic0)[ik100], +0.96734, 0.0, 1e-4);  // CW +1; measured +0.9673403636, residual 3.6e-7 (x275 slack; 2e-3 here was x5501, tightened 2026-09-15)

  // ---- (f) the whole curve stays inside CW's stated range [-2, 1].
  for (std::size_t ic = 0; ic < v.nc(); ++ic)
    for (double x : azz_tensor_curve(v, ic)) {
      if (std::isnan(x)) continue;
      CHECK(x > -2.001);
      CHECK(x <  1.001);
    }

  // ---- (g) REGRESSION GUARD.  The pre-2026-09-06 amplitude summed psi_2 = +W
  // with no i^L, which is CW Eq. (6.12) evaluated at MINUS f2/f0.  Pin the sign
  // that distinguishes them, so the old behaviour cannot come back silently:
  // at k = 0.30 GeV, theta_k ~ 0, on AV18, CW is NEGATIVE (-1.937); the old
  // code gave +0.617.
  CHECK(azz_tensor_curve(v, ic0)[ik30] < -1.5);
}

TEST_CASE("tagged: the acceptance-weighted curve reduces to the cell curve") {
  const TaggedModel& m = li6_model();
  const std::size_t ic = argmin_abs(m.c(), 0.0);
  const std::vector<double> ref = azz_tensor_curve(m, ic);
  std::vector<double> w(m.nk() * m.nc(), 0.0);
  for (std::size_t ik = 0; ik < m.nk(); ++ik) w[ik * m.nc() + ic] = 1.0;
  const std::vector<double> got = azz_tensor_curve_weighted(m, w);
  for (std::size_t i = 0; i < ref.size(); ++i) CHECK_CLOSE(got[i], ref[i], kRtol);
  // a uniform weight is the 4pi average: the L cross terms integrate away and
  // CG completeness makes the c-integral of n_M the same for every M
  const std::vector<double> flat =
      azz_tensor_curve_weighted(m, std::vector<double>(w.size(), 1.0));
  for (double v : flat) CHECK(std::fabs(v) < 1e-3);
  // a real acceptance table is a probability and, at the Yellow Report optics,
  // sculpts the sample into the LONGITUDINAL half of the sphere
  const BeamConfig cfg = default_configs("6Li")[1];
  const Optics o = yr_optics("6Li", cfg.ion_momentum_per_nucleon, true);
  const std::vector<double> eps = acceptance_weights(
      m, cfg.ion_momentum_per_nucleon, o,
      yr_config_key("6Li", cfg.ion_momentum_per_nucleon), 32);
  double lo = 1e300, hi = -1e300, wsum = 0.0, wabs = 0.0;
  for (std::size_t ik = 0; ik < m.nk(); ++ik) {
    for (std::size_t ii = 0; ii < m.nc(); ++ii) {
      const double e = eps[ik * m.nc() + ii];
      lo = std::fmin(lo, e);
      hi = std::fmax(hi, e);
      const double ww = e * m.n_of_kc(1.0)[ik * m.nc() + ii] * m.k()[ik] * m.k()[ik];
      wsum += ww;
      wabs += ww * std::fabs(m.c()[ii]);
    }
  }
  CHECK(lo >= 0.0);
  CHECK(hi <= 1.0);
  CHECK(wsum > 0.0);
  CHECK(wabs / wsum > 0.5);        // longitudinal at the YR optics
  // RE-PINNED 2026-09-06: BOTH signs flip with the `build_amp2` i^L fix
  // (docs/benchmarking/07_cw_sign_investigation.md section 7.1).  The
  // STATEMENT is unchanged and is the one worth keeping -- the YR optics
  // sculpt the sample into the longitudinal half of the sphere, where P2 has
  // the opposite sign to the theta_k = 90 deg cell, so the acceptance-weighted
  // curve is the cell curve's mirror.  Only the polarity of the pair moved:
  // was ref -0.4915 / weighted +0.5226, now ref +0.8966 / weighted -0.9533.
  const std::size_t jk = argmin_abs(m.k(), 0.30);
  CHECK(ref[jk] > 0.4);                                    // the 90 deg curve
  CHECK(azz_tensor_curve_weighted(m, eps)[jk] < -0.4);     // opposite sign
  CHECK_CLOSE_AT(ref[jk], +0.896640, 0.0, 1e-5);                     // MEASURED
  CHECK_CLOSE_AT(azz_tensor_curve_weighted(m, eps)[jk], -0.953283,
                 0.0, 1e-5);                                         // MEASURED
}

// ----------------------------------------------------- sampling and records

TEST_CASE("tagged: sample_kc reproduces its own density") {
  const TaggedModel& m = li6_model();
  Rng rng(2026, 0, 0, 21);
  const std::size_t n = 200000;
  std::vector<double> k, c, phi;
  m.sample_kc(1.0, 1.0, n, rng, k, c, phi);
  // cos-theta marginal, 16 bins aligned with 6 cells each of the 96-cell grid
  const std::vector<double>& a2 = m.amp2_table(1.0)[0];
  std::vector<double> prob_c(m.nc(), 0.0);
  double tot = 0.0;
  for (std::size_t ik = 0; ik < m.nk(); ++ik) {
    const double k2 = m.k()[ik] * m.k()[ik];
    for (std::size_t ic = 0; ic < m.nc(); ++ic) {
      prob_c[ic] += a2[ik * m.nc() + ic] * k2;
      tot += a2[ik * m.nc() + ic] * k2;
    }
  }
  std::vector<double> expect(16, 0.0);
  for (std::size_t ic = 0; ic < m.nc(); ++ic) expect[ic / 6] += prob_c[ic] / tot * n;
  std::vector<double> counts(16, 0.0);
  for (double v : c) {
    std::size_t b = static_cast<std::size_t>((v + 1.0) * 8.0);
    if (b > 15) b = 15;
    counts[b] += 1.0;
  }
  for (std::size_t b = 0; b < 16; ++b) {
    CHECK(std::fabs(counts[b] - expect[b]) / std::sqrt(std::fmax(expect[b], 1.0)) < 5.0);
  }
  // and the azimuth is uniform
  double cs = 0.0, sn = 0.0;
  for (double v : phi) { cs += std::cos(v); sn += std::sin(v); }
  CHECK(std::fabs(cs / n) < 0.01);
  CHECK(std::fabs(sn / n) < 0.01);
}

TEST_CASE("tagged: the P_D = 0 boost matches spectator.py quantile by quantile") {
  // Same density, same boost algebra, different (and unmatched) random
  // streams: the comparison is quantile by quantile, as in the Python.
  const TaggedModel model(li6_alpha_channel(BETA_DEFAULT, 0.0));
  const double p_u = default_configs("6Li")[1].ion_momentum_per_nucleon;
  const std::size_t n = 100000;
  Rng rng_a(31337, 0, 0, 23);
  std::vector<double> k, c, phi;
  model.sample_kc(1.0, 1.0, n, rng_a, k, c, phi);
  std::vector<double> pt_new, r_new, xl_new;
  for (std::size_t i = 0; i < n; ++i) {
    const SpectatorLab lab = boost_spectator(model.channel(), k[i], c[i], phi[i], p_u);
    pt_new.push_back(lab.pT);
    r_new.push_back(lab.R);
    xl_new.push_back(lab.xL);
  }
  const MomentumSampler ms(LI6_ALPHA_TAG());
  Rng rng_b(4242, 0, 0, 24);
  std::vector<double> pt_ref, r_ref, xl_ref;
  for (std::size_t i = 0; i < n; ++i) {
    double kx, ky, kz;
    ms.sample(rng_b, kx, ky, kz);
    const FragmentLab lab = boost_spectator_fragment(LI6_ALPHA_TAG(), p_u, kx, ky, kz);
    pt_ref.push_back(lab.pT);
    r_ref.push_back(lab.R);
    xl_ref.push_back(lab.xL);
  }
  for (double q : {0.25, 0.50, 0.75}) {
    CHECK_CLOSE(quantile(pt_new, q), quantile(pt_ref, q), 0.02);
    CHECK_CLOSE(quantile(r_new, q), quantile(r_ref, q), 0.02);
    CHECK_CLOSE(quantile(xl_new, q), quantile(xl_ref, q), 0.02);
  }
}

TEST_CASE("tagged: phi_spec is the LAB azimuth, not the DIS one") {
  const TaggedChannel ch = li6_alpha_channel();
  Rng rng(9, 0, 0, 31);
  double drift = 0.0;
  const int n = 5000;
  for (int i = 0; i < n; ++i) {
    const double k = 0.02 + 0.58 * rng.uniform();
    const double c = 2.0 * rng.uniform() - 1.0;
    const double phi_k = 2.0 * kPi * rng.uniform();
    const SpectatorLab lab = boost_spectator(ch, k, c, phi_k, 99.5);
    CHECK_CLOSE_AT(lab.phi_spec, std::atan2(lab.ky, lab.kx), 0.0, 1e-12);
    // the longitudinal boost does not touch the transverse plane, so with an
    // untilted axis the lab azimuth IS the spin-frame one
    CHECK_CLOSE_AT(std::cos(lab.phi_spec), std::cos(phi_k), 0.0, 1e-12);
    CHECK_CLOSE_AT(lab.pT * std::cos(lab.phi_spec), lab.kx, 0.0, 1e-12);
    const SpectatorLab tilt = boost_spectator(ch, k, c, phi_k, 99.5, 0.5);
    CHECK_CLOSE_AT(tilt.phi_spec, std::atan2(tilt.ky, tilt.kx), 0.0, 1e-12);
    drift += std::fabs(std::cos(tilt.phi_spec) - std::cos(phi_k));
  }
  CHECK(drift / n > 0.05);   // a tilted axis is rotated BEFORE the boost
}

TEST_CASE("tagged: the off-shell struck cluster and the light-front variables") {
  const TaggedChannel ch = li6_alpha_channel();
  const double p_u = 99.5;
  // at k = 0 the struck cluster carries exactly m_beam - m_spec = m_free - S
  const SpectatorLab rest = boost_spectator(ch, 0.0, 0.0, 0.0, p_u);
  const StruckCluster s0 = struck_cluster(ch, rest, p_u);
  CHECK_CLOSE(std::sqrt(s0.m2), ch.base.m_beam() - ch.base.m_spec(), 1e-12);
  // m_beam - m_spec and m_partner - S are the two routes to the same number;
  // they agree to 3.6 keV, which is the rounding of `separation_energy`, not
  // an inconsistency (spectator.py's m_partner docstring pins the same gap).
  CHECK_CLOSE(std::sqrt(s0.m2), s0.m_free - ch.base.separation_energy, 1e-7);
  CHECK(s0.virtuality < 0.0);
  CHECK_CLOSE(s0.alpha_s, ch.base.beam_A * ch.base.m_spec() / ch.base.m_beam(), 1e-12);
  CHECK_CLOSE_AT(s0.alpha_s, 4.0, 0.0, 0.01);
  CHECK_CLOSE(s0.alpha_s + s0.alpha_x, static_cast<double>(ch.base.beam_A), 1e-12);
  CHECK_CLOSE(s0.pt_s, 0.0, 0.0);
  // four-momentum is conserved exactly, the spectator is exactly on shell,
  // and the light-front fraction is invariant under the beam boost
  Rng rng(5, 0, 0, 1);
  for (int i = 0; i < 200; ++i) {
    const double k = 0.6 * rng.uniform();
    const double c = 2.0 * rng.uniform() - 1.0;
    const double phi = 2.0 * kPi * rng.uniform();
    const SpectatorLab lab = boost_spectator(ch, k, c, phi, p_u);
    const StruckCluster sc = struck_cluster(ch, lab, p_u);
    const Vec4 sum = sc.p + sc.p_spectator;
    CHECK_CLOSE(sum.e, sc.p_ion.e, 1e-13);
    CHECK_CLOSE_AT(sum.px, sc.p_ion.px, 0.0, 1e-12);
    CHECK_CLOSE_AT(sum.py, sc.p_ion.py, 0.0, 1e-12);
    CHECK_CLOSE(sum.pz, sc.p_ion.pz, 1e-13);
    CHECK_CLOSE(sc.p_spectator.m2(), ch.base.m_spec() * ch.base.m_spec(), 1e-9);
    CHECK(sc.virtuality < 0.0);          // off shell, and always below it
    CHECK(sc.alpha_s > 0.0);
    CHECK(sc.alpha_s < ch.base.beam_A);
    // alpha_s is the light-front ratio of the minus components, so it must be
    // the same computed from the LAB four-vectors
    const double alpha_lab = ch.base.beam_A * (lab.e_lab - lab.pz_lab)
                             / (sc.p_ion.e - sc.p_ion.pz);
    CHECK_CLOSE(sc.alpha_s, alpha_lab, 1e-9);
    CHECK_CLOSE(sc.pt_s, lab.pT, 1e-12);
  }
}

namespace {

/// A toy DIS side, so the sampler can be exercised without the P3 sampler.
struct ToyKinematics : KinematicsSource {
  void sample(double m_struck, int lam_e, double pe, std::size_t n, Rng& rng,
              std::vector<double>& x, std::vector<double>& q2,
              std::vector<double>& y, std::vector<double>& phi) override {
    (void)m_struck; (void)lam_e; (void)pe;
    x.resize(n); q2.resize(n); y.resize(n); phi.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      x[i] = 1e-3 * std::pow(500.0, rng.uniform());
      q2[i] = 1.5 + 98.5 * rng.uniform();
      y[i] = 0.5;
      phi[i] = 2.0 * kPi * rng.uniform();
    }
  }
  double sigma_tot_pb(double m_struck, int lam_e, double pe) const override {
    (void)lam_e; (void)pe;
    return 1.0 + 0.05 * m_struck;   // a mild m_S dependence, as a real one has
  }
};

}  // namespace

TEST_CASE("tagged: the sampler joins spin, spectator and DIS") {
  const TaggedModel& m = li6_model();
  const BeamConfig cfg = default_configs("6Li")[1];
  ToyKinematics dis;
  const TaggedSampler s(m, cfg.ion_momentum_per_nucleon, &dis);
  CHECK(s.pot_config() == "10x100");
  CHECK(s.optics().sigma_theta_v > 0.0);   // anisotropic-capable by default

  IonFill fill;
  fill.name = "t+";
  fill.j = 1.0;
  fill.populations = {1.0, 0.0, 0.0};
  Rng rng(11, 0, 0, 22);
  const std::vector<TaggedEvent> ev = s.sample_category(fill, 20000, rng);
  CHECK(ev.size() == 20000);
  double r_sum = 0.0;
  std::vector<double> r_all;
  std::size_t n_tail = 0, n_rp = 0, n_lost = 0;
  for (const TaggedEvent& e : ev) {
    CHECK(e.m_ion == 1.0);
    CHECK((e.m_struck == 1.0 || e.m_struck == 0.0 || e.m_struck == -1.0));
    CHECK(e.lab.pT >= 0.0);
    CHECK(e.lab.p_lab > 0.0);
    CHECK(e.x > 0.0);
    CHECK(e.q2 > 0.0);
    CHECK(e.route != kRouteZDC);
    CHECK(e.route >= 0);
    CHECK(e.route <= 6);
    r_sum += e.lab.R;
    r_all.push_back(e.lab.R);
    if (e.route == kRouteRPNearBeam) ++n_tail;
    if (e.route == kRouteRomanPots) ++n_rp;
    if (e.route == kRouteLost) ++n_lost;
  }
  // the 6Li alpha spectator is BEAM-BLIND: R ~ 1, so it reaches the pots
  // through the off-rigidity R < 0.95 slice and, at the per-configuration
  // Yellow Report envelope, hardly at all through the near-beam tail
  const double r_med = quantile(r_all, 0.5);
  CHECK(r_med > 0.97);
  CHECK(r_med < 1.03);
  CHECK(static_cast<double>(n_tail) / ev.size() < 0.005);
  CHECK(static_cast<double>(n_rp) / ev.size() > 0.01);
  CHECK(static_cast<double>(n_rp) / ev.size() < 0.05);
  CHECK(static_cast<double>(n_lost) / ev.size() > 0.5);
  // the rates are the fill populations x P(m_S|M) x the DIS rate
  const std::vector<double> rt = s.rates(fill);
  const std::vector<double> pms = m.population_integrated(1.0);
  for (std::size_t b = 0; b < 3; ++b) {
    CHECK_CLOSE(rt[b], pms[b] * dis.sigma_tot_pb(m.m_struck_values()[b], 0, 0.0),
                1e-12);
    CHECK(rt[3 + b] == 0.0);
  }
  CHECK_CLOSE(s.sigma_tot_pb(fill), rt[0] + rt[1] + rt[2], 1e-12);
}

TEST_CASE("tagged: fill_event adds the spectator and the struck cluster") {
  const TaggedModel& m = li7_model();
  const double p_u = default_configs("7Li")[1].ion_momentum_per_nucleon;
  const TaggedSampler s(m, p_u, nullptr);
  IonFill fill;
  fill.j = 1.5;
  fill.populations = {1.0, 0.0, 0.0, 0.0};
  Rng rng(3, 0, 0, 4);
  const TaggedEvent te = s.sample_one(fill, s.rate_cdf(fill), rng);

  Event ev;
  Particle beam_e;
  beam_e.pdg = 11;
  beam_e.status = Status::Beam;
  beam_e.role = Role::BeamElectron;
  beam_e.p = Vec4{10.0, 0.0, 0.0, -10.0};
  ev.particles.push_back(beam_e);
  Particle beam_i;
  beam_i.pdg = 1000030070;
  beam_i.status = Status::Beam;
  beam_i.role = Role::BeamIon;
  const double m_beam = m.channel().base.m_beam();
  const double p_beam = 7 * p_u;
  beam_i.p = Vec4{std::sqrt(p_beam * p_beam + m_beam * m_beam), 0.0, 0.0, p_beam};
  beam_i.mass = m_beam;
  beam_i.charge = 3;
  ev.particles.push_back(beam_i);

  s.fill_event(ev, te);
  CHECK(ev.particles.size() == 4);
  const Particle* spec = ev.find(Role::Spectator);
  const Particle* x = ev.find(Role::StruckCluster);
  REQUIRE(spec != nullptr);
  REQUIRE(x != nullptr);
  CHECK(spec->pdg == 1000020040);          // alpha
  CHECK(spec->charge == 2);
  CHECK(spec->status == Status::Final);
  CHECK(spec->mother1 == 1);
  CHECK(x->pdg == 1000010030);             // triton
  CHECK(x->charge == 1);
  CHECK(x->status == Status::Intermediate);
  // the two ADD UP to the beam ion, exactly
  const Vec4 sum = spec->p + x->p;
  CHECK_CLOSE(sum.e, beam_i.p.e, 1e-13);
  CHECK_CLOSE(sum.pz, beam_i.p.pz, 1e-13);
  CHECK_CLOSE_AT(sum.px, 0.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(sum.py, 0.0, 0.0, 1e-12);
  CHECK(spec->charge + x->charge == 3);
  // and the tagging block of the kinematics is filled
  CHECK_CLOSE(ev.kin.k, te.k, 1e-12);
  CHECK_CLOSE(ev.kin.cos_theta_k, te.cos_theta_k, 1e-12);
  CHECK_CLOSE(ev.kin.phi_k, te.phi_k, 1e-12);
  CHECK(ev.kin.alpha_s > 0.0);
  CHECK(ev.kin.alpha_s < 7.0);
  CHECK_CLOSE(ev.kin.pt_s, te.lab.pT, 1e-12);
  CHECK(ev.spin.m_ion == te.m_ion);
  CHECK(ev.spin.m_struck == te.m_struck);
  CHECK(ev.channel == Channel::TaggedLi7Alpha);
  CHECK(nuclide_pdg(1, 1) == 2212);
  CHECK(nuclide_pdg(0, 1) == 2112);
  CHECK(nuclide_pdg(3, 6) == 1000030060);
}

TEST_CASE("tagged: single-core throughput of sample_kc + boost") {
  // The T0 tagged inner loop without the DIS side: one (k, khat) draw from the
  // (280 x 96) cell CDF plus one lab boost and one far-forward route.
  const TaggedModel& m = li6_model();
  const double p_u = default_configs("6Li")[1].ion_momentum_per_nucleon;
  const TaggedSampler s(m, p_u, nullptr);
  IonFill fill;
  fill.j = 1.0;
  fill.populations = {1.0 / 3, 1.0 / 3, 1.0 / 3};
  const std::vector<double> cdf = s.rate_cdf(fill);
  Rng warm(1, 0, 0, 0);
  for (int i = 0; i < 1000; ++i) (void)s.sample_one(fill, cdf, warm);

  Rng rng(2, 0, 0, 0);
  const std::size_t n = 400000;
  const clock_t t0 = std::clock();
  double acc = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const TaggedEvent e = s.sample_one(fill, cdf, rng);
    acc += e.lab.pT + e.route;
  }
  const double secs = static_cast<double>(std::clock() - t0) / CLOCKS_PER_SEC;
  CHECK(acc > 0.0);
  MESSAGE("tagged sample_kc + boost + route throughput: "
          << (n / secs) / 1e6 << " Mevents/s single core ("
          << 1e9 * secs / n << " ns/event)");
  CHECK(n / secs > 1e5);   // the plans/08 P8 target for the T0 tier
}


// C3.  `TaggedModel`'s amplitude / density / cell-CDF grids used to be
// `mutable` maps filled lazily on first use, with no lock: two threads
// sampling different (M, m_S) states raced on the same std::map.  Thread
// safety rested on the caller having warmed every state first, which only
// `Pipeline`'s constructor did.  They are built by the constructor now, so a
// FRESH model driven directly from two threads has to be safe.
//
// This test drives TaggedSampler directly -- no Pipeline, no warm-up -- with
// two threads on a model constructed inside the test, and demands that the
// per-thread event streams are bit-identical to the single-threaded ones.
// The pre-fix code fails it under ThreadSanitizer (data races in
// std::_Rb_tree insert/find on all three maps) and segfaults outright in a
// -O1 build; at -O2 without a sanitizer the window is narrow, so the
// SANITIZER, not the plain run, is the gate this test is written for.
TEST_CASE("tagged: a fresh TaggedModel is safe on two threads with no warm-up") {
  const clock_t tc0 = std::clock();
  const TaggedModel model(li6_alpha_channel());
  const double build_ms =
      1e3 * static_cast<double>(std::clock() - tc0) / CLOCKS_PER_SEC;
  MESSAGE("TaggedModel eager grid build: " << build_ms << " ms");

  const double p_u = default_configs("6Li")[1].ion_momentum_per_nucleon;
  const TaggedSampler s(model, p_u, nullptr);
  IonFill fill;
  fill.j = 1.0;
  fill.populations = {0.5, 0.3, 0.2};   // all three M states carry rate
  const std::vector<double> cdf = s.rate_cdf(fill);

  // the reference streams, single threaded, on a model of their own
  const std::size_t n = 20000;
  auto stream = [&](const TaggedSampler& smp, std::uint64_t seed) {
    std::vector<double> out;
    out.reserve(4 * n);
    Rng rng(seed, 0, 0, 0);
    for (std::size_t i = 0; i < n; ++i) {
      const TaggedEvent e = smp.sample_one(fill, cdf, rng);
      out.push_back(e.k);
      out.push_back(e.cos_theta_k);
      out.push_back(e.m_ion);
      out.push_back(e.m_struck);
    }
    return out;
  };
  const std::vector<double> ref_a = stream(s, 11);
  const std::vector<double> ref_b = stream(s, 22);

  // ... and the same two streams from a model built fresh in this test and
  // driven concurrently.  Under the old lazy caches this is the race.
  const TaggedModel fresh(li6_alpha_channel());
  const TaggedSampler sf(fresh, p_u, nullptr);
  std::vector<double> got_a, got_b;
  std::thread ta([&] { got_a = stream(sf, 11); });
  std::thread tb([&] { got_b = stream(sf, 22); });
  ta.join();
  tb.join();

  REQUIRE(got_a.size() == ref_a.size());
  REQUIRE(got_b.size() == ref_b.size());
  bool same = true;
  for (std::size_t i = 0; i < ref_a.size(); ++i) {
    if (got_a[i] != ref_a[i] || got_b[i] != ref_b[i]) same = false;
  }
  CHECK(same);

  // the grids themselves are identical objects, built once
  for (double m : {1.0, 0.0, -1.0}) {
    CHECK(&model.amp2_table(m) == &model.amp2_table(m));   // pure lookup
    CHECK(model.n_of_kc(m).size() == fresh.n_of_kc(m).size());
  }
  // an M that is not a projection of this ion is refused rather than built:
  // no accessor may write to the object any more
  CHECK_THROWS(model.amp2_table(2.0));
  CHECK_THROWS(model.n_of_kc(2.0));
}

// ------------------------------------------- the VMC cluster-wave backend
//
// `ClusterWaveSource::VmcAV18` swaps the analytic radial forms for the ANL
// VMC tables.  Reference numbers: validation/vmc_reconcile.py and
// validation/vmc_tag_fractions.py, tabulated in
// docs/open_items/vmc_reconciliation.md.

namespace {

/// int k^2 sum_L psihat_L^2 dk and its moments, on the model's OWN grid --
/// the grid is part of the model (`TaggedModel`'s docstring), so this is the
/// distribution the sampler actually draws from.
struct GridMoments {
  double mean_k = 0.0;
  double p_gt_02 = 0.0, p_gt_03 = 0.0, p_gt_045 = 0.0;
  double p_d = std::nan("");
};

GridMoments grid_moments(const TaggedModel& m) {
  const std::vector<double>& k = m.k();
  std::vector<double> dens(k.size(), 0.0), d2(k.size(), 0.0);
  bool has_d = false;
  for (const Wave& w : m.channel().waves) {
    const std::vector<double>& r = m.radial_table(w.l);
    for (std::size_t i = 0; i < k.size(); ++i) {
      const double v = k[i] * k[i] * r[i] * r[i];
      dens[i] += v;
      if (w.l == 2) {
        d2[i] += v;
        has_d = true;
      }
    }
  }
  const double tot = trapezoid(dens, k);
  std::vector<double> kd(k.size());
  for (std::size_t i = 0; i < k.size(); ++i) kd[i] = k[i] * dens[i];
  GridMoments g;
  g.mean_k = trapezoid(kd, k) / tot;
  if (has_d) g.p_d = trapezoid(d2, k) / tot;
  const double cuts[3] = {0.20, 0.30, 0.45};
  double* into[3] = {&g.p_gt_02, &g.p_gt_03, &g.p_gt_045};
  for (int c = 0; c < 3; ++c) {
    std::vector<double> kk, yy;
    for (std::size_t i = 0; i < k.size(); ++i) {
      if (k[i] >= cuts[c]) {
        kk.push_back(k[i]);
        yy.push_back(dens[i]);
      }
    }
    *into[c] = trapezoid(yy, kk) / tot;
  }
  return g;
}

bool vmc_data_present() {
  std::ifstream f(data_path("vmc/momenta/li6_ad1.momentum"));
  return static_cast<bool>(f);
}

}  // namespace

TEST_CASE("tagged: the S-D interference SIGN agrees with b1_nuclear, and only "
          "the sign" * doctest::skip(!vmc_data_present())) {
  // `TaggedModel::build_amp2` applies the CDKS relative phase
  // `(-1)^floor(L/2)` and its comment claimed, until 2026-09-16, that this is
  // "identical to ClusterPartialWave::from_vmc / from_uw ... WHICH IS WHY THE
  // b1 AND TAGGED SECTORS THEN AGREE".  No test compared the two modules, and
  // they do NOT agree beyond the sign: this class renormalises each wave to
  // sqrt(P_L) on its own 280-point grid through `VmcRadial`'s LINEAR
  // interpolation, while `b1_nuclear` splines the file's own nodes with one
  // common unit factor, so the relative S/D MAGNITUDE differs by +0.78 % /
  // +0.36 % / -4.4 % at k = 0.197 / 1.003 / 2.529 fm^-1 (the last just above
  // the D node).  What DOES hold, and what this case pins, is the phase rule.
  //
  // MEASURED 2026-09-16, both modules on the VMC (AV18) path:
  //
  //   k [fm^-1]  k [GeV]    phi0        phi2        sign  n(c=+1) - n(c=0)
  //   0.197      0.03887   +6.605e+01  +8.888e-01    +     +1.159e+01
  //   1.003      0.19792   -6.984e+00  +1.374e+00    -     -2.055e+00
  //   2.529      0.49905   -2.174e-01  -2.907e-02    +     +1.134e-03
  const std::pair<ClusterPartialWave, ClusterPartialWave> pw =
      li6_alpha_d_partial_waves();
  const TaggedModel m(li6_alpha_channel(BETA_DEFAULT, P_D_LI6,
                                        ClusterWaveSource::VmcAV18));
  const std::vector<double>& k = m.k();
  const std::vector<double>& c = m.c();
  const std::size_t nc = c.size();
  const std::vector<double>& n1 = m.n_of_kc(1.0);

  // The forward-most and the equatorial cell of the M = +1 density: their
  // difference carries the S-D interference, whose sign is what the common
  // phase rule fixes.
  const std::size_t ic_fwd = nc - 1;
  const std::size_t ic_eq = nc / 2;
  int compared = 0;
  for (const double k_fm : {0.197, 1.003, 2.529}) {
    const double k_gev = k_fm * HBARC_GEV_FM;
    CAPTURE(k_fm);
    // nearest grid k of the tagged model
    std::size_t ik = 0;
    for (std::size_t i = 1; i < k.size(); ++i) {
      if (std::fabs(k[i] - k_gev) < std::fabs(k[ik] - k_gev)) ik = i;
    }
    const double phi0 = pw.first(k_gev);
    const double phi2 = pw.second(k_gev);
    const double interference = n1[ik * nc + ic_fwd] - n1[ik * nc + ic_eq];
    REQUIRE(phi0 != 0.0);
    REQUIRE(phi2 != 0.0);
    REQUIRE(interference != 0.0);
    CHECK((phi0 * phi2 > 0.0) == (interference > 0.0));
    ++compared;
  }
  CHECK(compared == 3);

  // ... and the RADIAL tables are where the two part company, in two separate
  // ways that the old comment's "the b1 and tagged sectors then agree"
  // covered up.
  //
  // (i) SIGN.  `TaggedModel` keeps its radial tables UNSIGNED-by-L and applies
  //     `(-1)^floor(L/2)` in `build_amp2`, on the ANGULAR factor; `b1_nuclear`
  //     folds the same phase into `phi` itself.  So psi2/psi0 = -(phi2/phi0)
  //     by construction -- which is why the sign comparison above has to be
  //     made on `n_of_kc`, where the phase has been applied, and not here.
  // (ii) MAGNITUDE.  Even then the two differ: this class renormalises each
  //      wave to sqrt(P_L) on its own 280-point grid with LINEAR
  //      interpolation, `b1_nuclear` splines the file's nodes with one common
  //      unit factor.  MEASURED 2026-09-16 at k = 1.003 fm^-1:
  //      phi2/phi0 = -0.195594 against psi2/psi0 = +0.196269, i.e. +0.35 % on
  //      the magnitude.  Bounded, not pinned: the bound is what the comment
  //      now claims (<= 0.8 % away from the nodes), and a bit-for-bit pin
  //      here would break on any re-spline.
  const double k_gev = 1.003 * HBARC_GEV_FM;
  const double ratio_b1 = pw.second(k_gev) / pw.first(k_gev);
  std::size_t ik = 0;
  for (std::size_t i = 1; i < k.size(); ++i) {
    if (std::fabs(k[i] - k_gev) < std::fabs(k[ik] - k_gev)) ik = i;
  }
  const double r2 = m.radial_table(2)[ik] / m.radial_table(0)[ik];
  MESSAGE("S-D at k = 1.003 fm^-1: b1_nuclear phi2/phi0 = " << ratio_b1
          << ", tagged psi2/psi0 = " << r2 << " (opposite sign by the "
          "build_amp2 phase convention; |ratio| differs by "
          << 100.0 * (std::fabs(r2 / ratio_b1) - 1.0) << " %)");
  CHECK((ratio_b1 > 0.0) != (r2 > 0.0));          // (i)
  const double mag = std::fabs(r2 / ratio_b1);
  CHECK(mag != 1.0);                              // (ii): NOT bit for bit ...
  CHECK(std::fabs(mag - 1.0) < 0.008);            // ... but <= 0.8 %
}

TEST_CASE("tagged: the VMC channels are built from the tables, and normalize" *
          doctest::skip(!vmc_data_present())) {
  REQUIRE(vmc_data_present());   // the decorator already guaranteed it
  const TaggedChannel c6 = li6_alpha_channel(BETA_DEFAULT, P_D_LI6,
                                             ClusterWaveSource::VmcAV18);
  REQUIRE(c6.waves.size() == 2);
  for (const Wave& w : c6.waves) {
    REQUIRE(w.vmc);
    CHECK(w.vmc->l() == w.l);
    CHECK(!w.vmc->provenance().empty());
  }
  // `p_d` is IGNORED on the VMC path: the D-state probability is a property
  // of the wave function, and it is the file's own printed value.
  CHECK_CLOSE_AT(c6.waves[1].prob, VMC_P_D_LI6, 0.0, 1e-15);
  CHECK_CLOSE_AT(c6.waves[0].prob, 1.0 - VMC_P_D_LI6, 0.0, 1e-15);
  c6.validate();

  const TaggedChannel c7 = li7_alpha_channel(BETA_DEFAULT,
                                             ClusterWaveSource::VmcAV18);
  REQUIRE(c7.waves.size() == 1);            // Aat11 is a selection-rule zero
  REQUIRE(c7.waves[0].vmc);
  CHECK(c7.waves[0].l == 1);

  for (const TaggedChannel& ch : {c6, c7}) {
    const TaggedModel m(ch);
    for (double mi : m_values(ch.j_ion)) {
      CHECK_CLOSE_AT(m.norm(mi), 1.0, 0.0, 1e-3);
    }
  }
}

TEST_CASE("tagged: VMC P_D(6Li) = 1.94% within the file's own MC error" *
          doctest::skip(!vmc_data_present())) {
  REQUIRE(vmc_data_present());   // the decorator already guaranteed it
  // The constant is the file's printed 0.015861 / (0.80362 + 0.015861).  The
  // 1-sigma MC errors on those two integrals are ~1e-4 relative on the S
  // block and ~1% on the D block (DRHOKA2/RHOKA2 near the D peak), so 1.94%
  // is good to about +-0.02% absolute -- and the INDEPENDENT 2004 AV18+UIX
  // overlap file gives 2.01%, a Hamiltonian difference, not an error.
  CHECK_CLOSE_AT(VMC_P_D_LI6, 0.0194, 0.0, 2e-4);
  const TaggedModel m(li6_alpha_channel(BETA_DEFAULT, P_D_LI6,
                                        ClusterWaveSource::VmcAV18));
  const GridMoments g = grid_moments(m);
  CHECK_CLOSE_AT(g.p_d, VMC_P_D_LI6, 0.0, 1e-6);
  // and it is 4.5x SMALLER than the scenario placeholder it replaces
  CHECK_CLOSE_AT(P_D_LI6 / VMC_P_D_LI6, 4.48, 0.0, 0.01);
}

TEST_CASE("tagged: the VMC moments reproduce the reconciled table" *
          doctest::skip(!vmc_data_present())) {
  REQUIRE(vmc_data_present());   // the decorator already guaranteed it
  // docs/open_items/vmc_reconciliation.md, "Moments on the TaggedModel grid".
  const GridMoments g6 = grid_moments(TaggedModel(
      li6_alpha_channel(BETA_DEFAULT, P_D_LI6, ClusterWaveSource::VmcAV18)));
  CHECK_CLOSE_AT(g6.mean_k, 0.1225, 0.0, 5e-4);
  CHECK_CLOSE_AT(g6.p_gt_02, 0.2410, 0.0, 5e-4);
  CHECK_CLOSE_AT(g6.p_gt_03, 0.0581, 0.0, 5e-4);
  CHECK_CLOSE_AT(g6.p_gt_045, 0.0019, 0.0, 5e-4);

  const GridMoments g7 = grid_moments(TaggedModel(
      li7_alpha_channel(BETA_DEFAULT, ClusterWaveSource::VmcAV18)));
  CHECK_CLOSE_AT(g7.mean_k, 0.1864, 0.0, 5e-4);
  CHECK_CLOSE_AT(g7.p_gt_03, 0.2111, 0.0, 5e-4);
  CHECK_CLOSE_AT(g7.p_gt_045, 0.0114, 0.0, 5e-4);

  // The independent Python reconciliation integrates the SAME tables on a
  // 4000-point grid capped at 5 fm^-1 = 0.98663 GeV instead of the model's
  // 280-point grid capped at 1.2 GeV.  Agreement at the percent level is the
  // statement that the grid is not doing the physics.
  CHECK_CLOSE_AT(g6.mean_k, 0.1224, 0.02, 0.0);
  CHECK_CLOSE_AT(g7.mean_k, 0.1860, 0.02, 0.0);
  CHECK_CLOSE_AT(g7.p_gt_03, 0.2120, 0.02, 0.0);

  // 7Li alpha+t is SOFTER than every beta in the model band -- against the
  // P-WAVE form the channel really uses.  (Thread B's "the band is biased
  // low" compared it against the S-wave form; see the reconciliation doc.)
  for (double beta : {BETA_BAND_LO, BETA_DEFAULT, BETA_BAND_HI}) {
    const GridMoments h = grid_moments(TaggedModel(li7_alpha_channel(beta)));
    CHECK(g7.mean_k < h.mean_k);
    CHECK(g7.p_gt_045 < h.p_gt_045);
  }
}

TEST_CASE("tagged: the Hulthen path is untouched by the VMC backend") {
  // The bit-compatibility guarantee.  The default channels carry NO table,
  // the deuteron control channel can never carry one (the Cosyn-Weiss gate
  // above runs on it), and `Wave::radial` on a null `vmc` is the analytic
  // switch to the last bit.
  for (const TaggedChannel& ch : {li6_alpha_channel(), li7_alpha_channel(),
                                  deuteron_channel()}) {
    for (const Wave& w : ch.waves) CHECK(!w.vmc);
  }
  // ... and the DEFAULT control channel carries no table either.  (This
  // comment used to read "there is no VmcAV18 overload of `deuteron_channel`
  // at all -- the deuteron IS the cluster, there is no d -> p + n two-cluster
  // table".  C5.4 disproved that on 2026-09-04: `fdeut.av18`'s u(k) and w(k)
  // ARE the p-n relative S and D waves, and `deuteron_channel(.., VmcAV18)`
  // now builds from them.  What this loop checks is unchanged and is the
  // point -- the DEFAULT is the analytic pair to the last bit.)
  const TaggedChannel d = deuteron_channel();
  const double kappa = d.base.kappa();
  for (const Wave& w : d.waves) {
    for (double k = 0.01; k < 1.2; k += 0.05) {
      Wave bare;
      bare.l = w.l;
      bare.beta = w.beta;
      CHECK(w.radial(k, kappa) == bare.radial(k, kappa));
    }
  }
}

TEST_CASE("tagged: the VMC S-D interference flips sign below the S node" *
          doctest::skip(!vmc_data_present())) {
  REQUIRE(vmc_data_present());   // the decorator already guaranteed it
  // With the global phase fixed by psi_0(k -> 0) > 0 the model's own radial
  // tables carry the sign, and it is NOT the sign the positive-definite
  // Hulthen forms assume.  The S node is at 0.678 fm^-1 = 0.1338 GeV.
  const TaggedModel m(li6_alpha_channel(BETA_DEFAULT, P_D_LI6,
                                        ClusterWaveSource::VmcAV18));
  const std::vector<double>& s = m.radial_table(0);
  const std::vector<double>& d = m.radial_table(2);
  const std::size_t below = argmin_abs(m.k(), 0.08);
  const std::size_t above = argmin_abs(m.k(), 0.20);
  CHECK(s[below] > 0.0);
  CHECK(d[below] < 0.0);
  CHECK(s[below] * d[below] < 0.0);     // opposite to Hulthen
  CHECK(s[above] < 0.0);
  CHECK(d[above] < 0.0);
  CHECK(s[above] * d[above] > 0.0);     // same as Hulthen, above the S node
  // the Hulthen forms have no node at all
  const TaggedModel h(li6_alpha_channel());
  for (std::size_t i = 0; i < h.nk(); ++i) {
    CHECK(h.radial_table(0)[i] > 0.0);
    CHECK(h.radial_table(2)[i] > 0.0);
  }
}

// ================================================== open items C5.2 - C5.5
//
// T24.  The ANL Monte Carlo errors are CARRIED, not parsed and dropped, and
// the band they buy is measured rather than asserted.
// docs/open_items/run_2026-09-03/phase_C_numbers.md sec. C5.2.
TEST_CASE("T24 the VMC tables carry their Monte Carlo band" *
          doctest::skip(!vmc_data_present())) {
  REQUIRE(vmc_data_present());   // the decorator already guaranteed it
  const TaggedChannel c0 =
      li6_alpha_channel(BETA_DEFAULT, P_D_LI6, ClusterWaveSource::VmcAV18);
  REQUIRE(c0.waves.size() == 2);
  for (const Wave& w : c0.waves) {
    REQUIRE(w.vmc);
    CHECK(w.vmc->has_errors());
    CHECK(w.vmc->dpsi().size() == w.vmc->psi().size());
    // dpsi is SIGNED and carries psi's sign, so psi + n dpsi is covariant
    // under the global-phase flip the pair is fixed with.
    for (std::size_t i = 0; i < w.vmc->psi().size(); ++i) {
      if (w.vmc->psi()[i] != 0.0) {
        CHECK(w.vmc->psi()[i] * w.vmc->dpsi()[i] > 0.0);
      }
    }
    double cor = 0.0, quad = 0.0;
    w.vmc->norm2_error(&cor, &quad);
    CHECK(cor > quad);                     // correlated is the envelope
    CHECK(cor / w.vmc->norm2() < 0.02);    // and it is a 2 % effect at most
  }
  // The two limits, pinned: S wave 0.47 % / 0.11 %, D wave 1.62 % / 0.40 %.
  {
    double cor = 0.0, quad = 0.0;
    c0.waves[0].vmc->norm2_error(&cor, &quad);
    CHECK_CLOSE(cor / c0.waves[0].vmc->norm2(), 0.004705, 1e-2);
    CHECK_CLOSE(quad / c0.waves[0].vmc->norm2(), 0.001085, 1e-2);
    c0.waves[1].vmc->norm2_error(&cor, &quad);
    CHECK_CLOSE(cor / c0.waves[1].vmc->norm2(), 0.016176, 1e-2);
    CHECK_CLOSE(quad / c0.waves[1].vmc->norm2(), 0.003974, 1e-2);
  }
  // n_sigma = 0 is bit for bit the cached table.
  CHECK(li6_alpha_channel(BETA_DEFAULT, P_D_LI6, ClusterWaveSource::VmcAV18,
                          0.0).waves[1].prob == c0.waves[1].prob);
  // THE BAND, on the observables.  P_D moves ~1.1 % per sigma; the tagged
  // tensor dilution moves 0.02 %.  The ANL statistics are NOT the systematic
  // that matters -- the wave-function choice, 6.6 %, is (T26).
  const TaggedModel m0(c0);
  const TaggedModel mp(li6_alpha_channel(BETA_DEFAULT, P_D_LI6,
                                         ClusterWaveSource::VmcAV18, +1.0));
  const TaggedModel mm(li6_alpha_channel(BETA_DEFAULT, P_D_LI6,
                                         ClusterWaveSource::VmcAV18, -1.0));
  CHECK_CLOSE(mp.channel().waves[1].prob / VMC_P_D_LI6, 1.011346, 1e-4);
  CHECK_CLOSE(mm.channel().waves[1].prob / VMC_P_D_LI6, 0.988851, 1e-4);
  CHECK(std::fabs(mp.tensor_dilution() / m0.tensor_dilution() - 1.0) < 3e-4);
  CHECK_CLOSE(mp.tensor_dilution() / m0.tensor_dilution() - 1.0, -2.012e-4,
              2e-2);
  CHECK_CLOSE(mm.tensor_dilution() / m0.tensor_dilution() - 1.0, +1.977e-4,
              2e-2);
  // A band that would be silently ignored is REFUSED.
  CHECK_THROWS(li6_alpha_channel(BETA_DEFAULT, P_D_LI6,
                                 ClusterWaveSource::Hulthen, 1.0));
  CHECK_THROWS(li7_alpha_channel(BETA_DEFAULT, ClusterWaveSource::Hulthen,
                                 1.0));
  // 7Li is one wave: the band moves the SHAPE and never a probability.
  CHECK(li7_alpha_channel(BETA_DEFAULT, ClusterWaveSource::VmcAV18, 1.0)
            .waves[0].prob == 1.0);
}

// T25.  The deuteron control on the exact AV18 wave function (C5.4).
TEST_CASE("T25 the deuteron control channel, Hulthen against AV18" *
          doctest::skip(!vmc_data_present())) {
  REQUIRE(vmc_data_present());   // the decorator already guaranteed it
  // The k-space block's own D fraction reproduces the file's r-space header
  // `dstate` = 0.057599 to 1.55e-5 relative -- the reader and the convention
  // validated by the file against itself.  The residual is the file's own
  // r-space/k-space quadrature spread and NOT print rounding: 0.05759989
  // would print as 0.057600 at the header's six figures, not as 0.057599.
  CHECK_CLOSE(deuteron_av18_p_d(), 0.057599, 2e-5);
  CHECK(std::fabs(deuteron_av18_p_d() / 0.057599 - 1.0) > 1e-6);
  CHECK_CLOSE(deuteron_av18_p_d(), 0.0575998919874, 1e-9);
  // ... and it is 28 % ABOVE the scenario P_D_DEUTERON the control is pinned
  // to.
  CHECK_CLOSE(deuteron_av18_p_d() / P_D_DEUTERON, 1.28000, 1e-4);

  const TaggedModel h(deuteron_channel());
  const TaggedModel v(deuteron_channel(BETA_DEFAULT, P_D_DEUTERON,
                                       ClusterWaveSource::VmcAV18));
  CHECK(v.channel().waves[1].prob == deuteron_av18_p_d());
  // RE-PINNED 2026-09-06 with the `build_amp2` i^L fix
  // (docs/benchmarking/07_cw_sign_investigation.md).  The dilutions are
  // ANGLE-INTEGRATED, so L-orthogonality kills the S-D cross term and the
  // physics does not move; what moves is the 96-cell midpoint quadrature
  // residual of the integral int Theta_0 Theta_2 dc = 0, whose SIGN the
  // phase flips.  That is 1.4e-6 relative on the vector dilutions -- just
  // past this gate's 1e-6 -- and 8e-7 on the tensor ones.  Pre-fix, for the
  // record: vector 0.932494769 / 0.913594777, tensor 0.959488074 /
  // 0.948145618.  The v/h RATIOS are unchanged to every digit they pin.
  CHECK_CLOSE(h.vector_dilution(), 0.932496109, 1e-6);   // MEASURED 0.932496109312
  CHECK_CLOSE(v.vector_dilution(), 0.913595979, 1e-6);   // MEASURED 0.913595978560
  CHECK_CLOSE(v.vector_dilution() / h.vector_dilution() - 1.0, -0.020268, 1e-3);
  CHECK_CLOSE(h.tensor_dilution(), 0.959488878, 1e-6);   // MEASURED 0.959488878164
  CHECK_CLOSE(v.tensor_dilution(), 0.948146339, 1e-6);   // MEASURED 0.948146339359
  CHECK_CLOSE(v.tensor_dilution() / h.tensor_dilution() - 1.0, -0.011822, 1e-3);
  // THE RELATIVE S-D SIGN DOES NOT FLIP.  CDKS fix phi_2 = -W and
  // phi_L = i^L psi_L, so the physical deuteron has psi_2 = +W > 0 at low k,
  // which is what the positive-definite Hulthen forms already assume OF THE
  // STORED TABLE -- the two CHECKs below are on `vmc->psi()`, which is what
  // this paragraph is about.  It does NOT license summing psi_2 into the
  // amplitude: `build_amp2` consumes phi_L = i^L psi_L, i.e. phi_2 = -W, and
  // applied no phase at all until 2026-09-06 (the re-pin above).  This is the
  // OPPOSITE of the 6Li alpha-d case (see the S-node test above).
  CHECK(v.channel().waves[0].vmc->psi()[1] > 0.0);
  CHECK(v.channel().waves[1].vmc->psi()[1] > 0.0);
  // fdeut.av18 prints no MC error column, so this table carries no band.
  CHECK(!v.channel().waves[1].vmc->has_errors());
  // Hulthen stays bit for bit the default.
  CHECK(h.channel().waves[1].prob == P_D_DEUTERON);
}

// T26.  C5.3 (the N_ad / P_D spreads) and C5.5 (the inclusive-tagged drift),
// in one place because they are the same question asked of two constants.
TEST_CASE("T26 the alpha-d normalisation spread, and the inclusive drift" *
          doctest::skip(!vmc_data_present())) {
  REQUIRE(vmc_data_present());   // the decorator already guaranteed it
  // --- C5.3.  THREE NUMBERS SPAN THREE DIFFERENT AMOUNTS, and the document
  // that says "N_ad 5 %, P_D 7 %" is quoting two of them about a third: the
  // 7 % belongs to the D-wave NORM, and P_D -- the RATIO the tagged sector
  // actually uses -- spans only 2.7 %.
  struct Row { double s, d; };
  const Row r14{0.80362, 0.015861};   // 2014 li6_ad1.momentum, the default
  const Row r04{0.838, 0.017};        // 2004 li6.ad, as printed
  const Row rw{0.846, 0.017};         // Wiringa PRC 89 (2014) 024305
  CHECK_CLOSE(r14.s + r14.d, VMC_N_ALPHA_D_LI6, 1e-15);
  CHECK_CLOSE((rw.s + rw.d) / (r14.s + r14.d) - 1.0, 0.05311, 1e-3);
  CHECK_CLOSE(r04.d / r14.d - 1.0, 0.07181, 1e-3);
  const double p14 = r14.d / (r14.s + r14.d);
  const double p04 = r04.d / (r04.s + r04.d);
  const double pw = rw.d / (rw.s + rw.d);
  CHECK_CLOSE(p14, VMC_P_D_LI6, 1e-15);
  CHECK_CLOSE(p04 / p14 - 1.0, 0.02729, 1e-3);
  // ... and IT PROPAGATES TO ALMOST NOTHING here: `TaggedModel` renormalises
  // each wave to its own P_L, so N_ad drops out of every tagged observable
  // and only the RATIO survives.  0.05 % on the tensor dilution.
  TaggedChannel c = li6_alpha_channel(BETA_DEFAULT, P_D_LI6,
                                      ClusterWaveSource::VmcAV18);
  const TaggedModel m14(c);
  c.waves[0].prob = 1.0 - p04;
  c.waves[1].prob = p04;
  const TaggedModel m04(c);
  c.waves[0].prob = 1.0 - pw;
  c.waves[1].prob = pw;
  const TaggedModel mw(c);
  CHECK_CLOSE(m04.tensor_dilution() / m14.tensor_dilution() - 1.0, -4.84e-4,
              5e-2);
  CHECK_CLOSE(m04.vector_dilution() / m14.vector_dilution() - 1.0, -8.16e-4,
              5e-2);
  CHECK(std::fabs(mw.tensor_dilution() / m14.tensor_dilution() - 1.0) < 1e-3);

  // --- C5.5.  On the HULTHEN default the inclusive constant and the tagged
  // model ARE one wave function, to 7.84e-6 (0.869950 closed form against the
  // measured 0.8699431789; 1.22e-5 before the 2026-09-06 S-D sign fix moved
  // the quadrature residual).  The CHECK below stays at 2e-5.
  const TaggedModel hul(li6_alpha_channel());
  CHECK_CLOSE(hul.vector_dilution(), ALPHA_D_VECTOR_POLARIZATION, 2e-5);
  CHECK_CLOSE(LI6_CLUSTER_POLARIZATION,
              li6_cluster_polarization(P_D_LI6, P_D_DEUTERON), 1e-15);
  // Under VmcAV18 they DRIFT: +11.61 % vector, +6.58 % rank-2.
  CHECK_CLOSE(m14.vector_dilution() / hul.vector_dilution() - 1.0, 0.116131,
              1e-4);
  CHECK_CLOSE(m14.tensor_dilution() / LI6_B1_RANK2_TRANSFER - 1.0, 0.065762,
              1e-4);
  CHECK_CLOSE(LI6_CLUSTER_POLARIZATION_VMC / LI6_CLUSTER_POLARIZATION - 1.0,
              0.116119, 1e-4);
  // AND SUBSTITUTION IS NOT THE FIX.  The "consistent" reading overshoots the
  // ab initio six-body VMC number by 6.8 % where the shipped one undershoots
  // it by 4.3 %; adding the AV18 deuteron's own P_D leaves 4.6 % over.
  CHECK_CLOSE(LI6_CLUSTER_POLARIZATION / LI6_POLARIZATION_VMC_SIX_BODY - 1.0,
              -0.04336, 1e-3);
  CHECK_CLOSE(LI6_CLUSTER_POLARIZATION_VMC / LI6_POLARIZATION_VMC_SIX_BODY
                  - 1.0,
              +0.06772, 1e-3);
  const double both = li6_cluster_polarization(VMC_P_D_LI6,
                                               deuteron_av18_p_d());
  CHECK_CLOSE(both, 0.8870761569, 1e-8);
  CHECK_CLOSE(both / LI6_POLARIZATION_VMC_SIX_BODY - 1.0, +0.04608, 1e-3);
  CHECK(std::fabs(LI6_CLUSTER_POLARIZATION / LI6_POLARIZATION_VMC_SIX_BODY
                  - 1.0)
        < std::fabs(LI6_CLUSTER_POLARIZATION_VMC
                        / LI6_POLARIZATION_VMC_SIX_BODY - 1.0));
}


// T27.  C5.5b -- ONE RUN, ONE DEUTERON.
//
// `--cluster-wave vmc` selects the ANL VMC AV18+UX alpha-d overlap for the
// alpha-d RELATIVE motion.  Until 2026-09-04 the two places a tagged-alpha run
// reads the EMBEDDED deuteron's own wave function did not follow it:
//
//   * `TaggedChannel::dis_target` was `DEUTERON()` unconditionally, so the
//     struck cluster's g1 (`InclusiveKernel::g1a` -> `PolSF::g1_nucleus`) used
//     the SCENARIO `P_D_DEUTERON` = 0.045, i.e. 1 - 1.5 P_D = 0.9325;
//   * the T1 `ClusterBreakup` drew the struck-nucleon spin from the Hulthen
//     deuteron at the same 0.045.
//
// The AV18 deuteron belonging to that overlap has P_D = 0.0575998919874
// (`deuteron_av18_p_d`, T25), i.e. 0.9136 -- so every polarized tagged-alpha
// observable under the flag was 2.069 % HIGH against the wave function the flag
// claims to select, and the run held two deuteron wave-function families.  This
// pins that it now holds one, and that the Hulthen default is untouched.
TEST_CASE("T27 the embedded deuteron follows --cluster-wave" *
          doctest::skip(!vmc_data_present())) {
  REQUIRE(vmc_data_present());   // the decorator already guaranteed it
  // --- what the two DIS targets are.
  const TaggedChannel h = li6_alpha_channel();
  const TaggedChannel v = li6_alpha_channel(BETA_DEFAULT, P_D_LI6,
                                            ClusterWaveSource::VmcAV18);
  // the default is `DEUTERON()` itself, bit for bit
  CHECK(h.dis_target.eff_pol_p == DEUTERON().eff_pol_p);
  CHECK(h.dis_target.eff_pol_n == DEUTERON().eff_pol_n);
  CHECK(h.dis_target.eff_pol_p == DEUTERON_VECTOR_POLARIZATION);
  // and the VMC one differs from it in the effective polarizations ALONE
  CHECK(v.dis_target.name == DEUTERON().name);
  CHECK(v.dis_target.A == DEUTERON().A);
  CHECK(v.dis_target.Z == DEUTERON().Z);
  CHECK(v.dis_target.spin == DEUTERON().spin);
  CHECK(v.dis_target.eff_pol_p == vector_dilution_of(deuteron_av18_p_d()));
  CHECK(v.dis_target.eff_pol_n == v.dis_target.eff_pol_p);
  CHECK_CLOSE(v.dis_target.eff_pol_p, 0.9136001620189, 1e-12);

  // --- THE SIZE THAT WAS WRONG.  2.0687 %, and it is EXACT rather than
  // approximate: g1A = Z P_p g1p + N P_n g1n is linear in the effective
  // polarization, so the whole polarized sector carried one factor.
  const double was_high = DEUTERON().eff_pol_p / v.dis_target.eff_pol_p;
  CHECK_CLOSE(was_high - 1.0, 0.0206872095, 1e-8);
  const InclusiveKernel k_now(v.dis_target);
  const InclusiveKernel k_before(DEUTERON());
  const double xs[4] = {0.05, 0.1, 0.3, 0.5};
  const double q2s[4] = {2.0, 5.0, 10.0, 20.0};
  for (int i = 0; i < 4; ++i) {
    CHECK_CLOSE(k_before.tables(xs[i], q2s[i]).g1
                    / k_now.tables(xs[i], q2s[i]).g1,
                was_high, 1e-12);
    CHECK_CLOSE(k_before.a_parallel(k_before.tables(xs[i], q2s[i]), xs[i],
                                    q2s[i], 0.5)
                    / k_now.a_parallel(k_now.tables(xs[i], q2s[i]), xs[i],
                                       q2s[i], 0.5),
                was_high, 1e-12);
  }

  // --- THE T1 BREAKUP FOLLOWS TOO, and the gate is the identity breakup.hpp
  // states: the dilution the spin DRAW implies is the one the RATE uses.
  BreakupOptions bo;
  const ClusterBreakup bh(bo);
  CHECK(bh.options().source == ClusterWaveSource::Hulthen);
  CHECK_CLOSE(bh.deuteron_model().vector_dilution(), h.dis_target.eff_pol_p,
              1e-5);
  bo.source = ClusterWaveSource::VmcAV18;
  const ClusterBreakup bv(bo);
  CHECK_CLOSE(bv.deuteron_model().vector_dilution(), v.dis_target.eff_pol_p,
              1e-5);
  CHECK(bv.deuteron_model().channel().waves[1].prob == deuteron_av18_p_d());
  // ... and the two are 2.0688 % apart, which is what a run that mixed them
  // would have been wrong by.
  CHECK_CLOSE(bh.deuteron_model().vector_dilution()
                  / bv.deuteron_model().vector_dilution() - 1.0,
              0.0206879, 1e-4);

  // --- WHAT ONE `--cluster-wave vmc` RUN NOW SAYS 6Li's POLARIZATION IS.
  // alpha-d from the VMC overlap x the AV18 deuteron = C5.5's THIRD table row,
  // 0.887076 -- not the mongrel 0.905427 (VMC alpha-d x scenario deuteron)
  // it used to be.  Both are still readings of the cluster PRODUCT and neither
  // reproduces the ab initio 0.848; that is C5.5 and it is unchanged.
  const TaggedModel mv(v);
  CHECK_CLOSE(mv.vector_dilution() * v.dis_target.eff_pol_p,
              li6_cluster_polarization(VMC_P_D_LI6, deuteron_av18_p_d()),
              2e-5);
  CHECK_CLOSE(mv.vector_dilution() * v.dis_target.eff_pol_p, 0.887075, 1e-5);
  // ... i.e. 2.027 % BELOW the mongrel, which is the whole size of the repair.
  CHECK_CLOSE(mv.vector_dilution() * v.dis_target.eff_pol_p
                  / LI6_CLUSTER_POLARIZATION_VMC - 1.0,
              -0.020270, 2e-3);
  // the Hulthen run's own reading is `LI6_CLUSTER_POLARIZATION`, as it was
  const TaggedModel mh(h);
  CHECK_CLOSE(mh.vector_dilution() * h.dis_target.eff_pol_p,
              LI6_CLUSTER_POLARIZATION, 2e-5);

  // --- AND THE 7Li CHANNEL DELIBERATELY DOES NOT MOVE: there is no AV18
  // A = 3 wave function in this tree to switch the triton's internal spin
  // structure to, so the flag selects the alpha-t relative motion alone.
  CHECK(li7_alpha_channel().dis_target.eff_pol_p == TRITON().eff_pol_p);
  CHECK(li7_alpha_channel(BETA_DEFAULT, ClusterWaveSource::VmcAV18)
            .dis_target.eff_pol_p == TRITON().eff_pol_p);
  CHECK(li7_alpha_channel(BETA_DEFAULT, ClusterWaveSource::VmcAV18)
            .dis_target.eff_pol_n == TRITON().eff_pol_n);
  // the control channel's struck object is a free neutron on both settings
  CHECK(deuteron_channel().dis_target.eff_pol_n == NEUTRON_TARGET().eff_pol_n);
  CHECK(deuteron_channel(BETA_DEFAULT, P_D_DEUTERON,
                         ClusterWaveSource::VmcAV18)
            .dis_target.eff_pol_n == NEUTRON_TARGET().eff_pol_n);
}
