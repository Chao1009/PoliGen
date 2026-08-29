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
#include <string>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "json_min.hpp"
#include "lipolgen/asymmetries.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/spin.hpp"
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

double quantile(std::vector<double> v, double q) {
  std::sort(v.begin(), v.end());
  const double pos = q * static_cast<double>(v.size() - 1);
  const std::size_t i = static_cast<std::size_t>(pos);
  const double f = pos - static_cast<double>(i);
  return i + 1 < v.size() ? v[i] * (1.0 - f) + v[i + 1] * f : v.back();
}

}  // namespace

// --------------------------------------------------------- reference tables

TEST_CASE("tagged: channel construction against polligen") {
  jsonmin::Value ref;
  if (!load_tagged(ref)) {
    MESSAGE("tagged.json not found -- skipping");
    return;
  }
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

TEST_CASE("tagged: the model grid and its tables against polligen") {
  jsonmin::Value ref;
  if (!load_tagged(ref)) {
    MESSAGE("tagged.json not found -- skipping");
    return;
  }
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

TEST_CASE("tagged: boost_spectator against polligen (closed form, rtol 1e-12)") {
  jsonmin::Value ref;
  if (!load_tagged(ref)) {
    MESSAGE("tagged.json not found -- skipping");
    return;
  }
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
  const TaggedModel& m = deut_model();
  const std::size_t ic = argmin_abs(m.c(), 0.0);
  const std::vector<double> a = azz_tensor_curve(m, ic);
  double thresh = 0.0, best = 0.0;
  std::size_t jbest = 0;
  for (std::size_t ik = 0; ik < m.nk(); ++ik) {
    if (m.k()[ik] < 0.02) thresh = std::fmax(thresh, std::fabs(a[ik]));
    if (std::fabs(a[ik]) > best) { best = std::fabs(a[ik]); jbest = ik; }
  }
  CHECK(thresh < 0.02);                      // D-wave threshold suppression
  double window = 0.0;
  for (std::size_t ik = 0; ik < m.nk(); ++ik) {
    if (m.k()[ik] > 0.25 && m.k()[ik] < 0.5) window = std::fmax(window, std::fabs(a[ik]));
  }
  CHECK(window > 0.45);                      // O(1) in the CW window
  CHECK(m.k()[jbest] > 0.2);
  CHECK(m.k()[jbest] < 0.45);
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

TEST_CASE("tagged: the Cosyn-Weiss deuteron tensor gate (CW TABLE II)") {
  // Cosyn-Weiss II (arXiv:2603.23700) p. 35, Eq. (6.12)-(6.13): a ratio of
  // quadratic forms in the S and D radials times (1 - 3 cos^2 theta_k) =
  // -2 P2(cos theta_k), taking values in [-2, 1], the quadratic form peaking
  // at +1 where f2/f0 = sqrt(2) (k = 0.30 GeV for AV18).  Our A_zz^wf maps
  // onto theirs as A_T|| = -2 A_zz^wf.
  const TaggedModel& m = deut_model();
  std::vector<double> p2(m.nc());
  for (std::size_t i = 0; i < m.nc(); ++i) p2[i] = 0.5 * (3.0 * m.c()[i] * m.c()[i] - 1.0);

  // (a) the P2 factorization is EXACT: A_zz^wf / P2 does not depend on the
  //     angle bin at fixed k.  Cells within 1e-3 of the P2 zero are excluded,
  //     where the ratio is unbounded and says nothing.
  const std::size_t ik = argmin_abs(m.k(), 0.30);
  CHECK_CLOSE_AT(m.k()[ik], 0.3012, 0.0, 5e-4);
  double rmin = 1e300, rmax = -1e300, rsum = 0.0;
  std::size_t nr = 0;
  for (std::size_t ic = 0; ic < m.nc(); ++ic) {
    if (std::fabs(p2[ic]) <= 1e-3) continue;
    const double r = azz_tensor_curve(m, ic)[ik] / p2[ic];
    rmin = std::fmin(rmin, r);
    rmax = std::fmax(rmax, r);
    rsum += r;
    ++nr;
  }
  CHECK(rmax - rmin < 1e-5);
  CHECK_CLOSE_AT(rsum / nr, 0.99940, 0.0, 1e-4);

  // (b) the k envelope reaches its maximum 1 at f2/f0 = sqrt(2)
  const std::size_t ic0 = argmax_abs(m.c());   // nearest cell to theta_k = 0
  const std::vector<double> curve0 = azz_tensor_curve(m, ic0);
  std::vector<double> envelope(m.nk());
  for (std::size_t i = 0; i < m.nk(); ++i) envelope[i] = curve0[i] / p2[ic0];
  std::size_t j = 0;
  for (std::size_t i = 1; i < m.nk(); ++i) if (envelope[i] > envelope[j]) j = i;
  CHECK_CLOSE_AT(envelope[j], 1.0, 0.0, 1e-3);
  CHECK_CLOSE_AT(m.k()[j], 0.3098, 0.0, 5e-4);
  CHECK_CLOSE_AT(m.k()[j], 0.30, 0.0, 0.02);   // against CW's 0.30 GeV

  // (c) A_T|| = -2 A_zz^wf against CW TABLE II's +1 and -2.  The outermost
  //     cos theta_k cell is 0.9896, not 1, so the exact extremes come through
  //     the P2 factorization pinned in (a).
  const std::size_t ic90 = argmin_abs(m.c(), 0.0);
  const double a_par_90 = -2.0 * azz_tensor_curve(m, ic90)[j];
  const double a_par_0 = -2.0 * azz_tensor_curve(m, ic0)[j];
  CHECK_CLOSE_AT(a_par_90, 0.9997, 0.0, 1e-3);       // CW: +1
  CHECK_CLOSE_AT(a_par_0, -1.9378, 0.0, 2e-3);       // cell centre
  CHECK_CLOSE_AT(-2.0 * envelope[j] * 1.0, -2.0, 0.0, 3e-3);   // CW: -2
  CHECK_CLOSE_AT(-2.0 * envelope[j] * (-0.5), 1.0, 0.0, 3e-3); // CW: +1
  // and the whole curve stays inside CW's stated range [-2, 1]
  for (std::size_t ic = 0; ic < m.nc(); ++ic) {
    for (double v : azz_tensor_curve(m, ic)) {
      if (std::isnan(v)) continue;
      CHECK(-2.0 * v > -2.001);
      CHECK(-2.0 * v < 1.001);
    }
  }
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
  const std::size_t jk = argmin_abs(m.k(), 0.30);
  CHECK(ref[jk] < -0.4);                                   // the 90 deg curve
  CHECK(azz_tensor_curve_weighted(m, eps)[jk] > 0.4);      // opposite sign
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
