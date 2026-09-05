// SPDX-License-Identifier: GPL-3.0-or-later
// The cross-language gate: every number in validation/reference/*.json was
// produced by the Python (PolarizedLithiumSim) through
// validation/dump_polligen_reference.py, and the C++ port must reproduce it
// at rtol 1e-12.
//
// The tables are produced by another task in this repository.  When the
// directory is empty the test SKIPS with a message rather than failing, so
// the suite is green before and after they land.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "json_min.hpp"
#include "lipolgen/asymmetries.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/spin.hpp"
#include "lipolgen/xsec.hpp"

using namespace lipolgen;

#ifndef LIPOLGEN_REFERENCE_DIR
#define LIPOLGEN_REFERENCE_DIR "validation/reference"
#endif

namespace {

constexpr double kRtol = 1e-12;
/// The maxent populations come out of a bisection with tol 1e-13, so the
/// reference README asks for a looser comparison on that table only.  The
/// measured agreement is in fact 1-2 ulps (this bisection is transcribed
/// step for step, including the exp(beta (m - m_max)) stabilisation); the
/// margin here is for other libm implementations, not for this one.
constexpr double kMaxentRtol = 1e-13;

std::string ref_path(const char* name) {
  return std::string(LIPOLGEN_REFERENCE_DIR) + "/" + name;
}

bool load(const char* name, jsonmin::Value& out) {
  return jsonmin::load_file(ref_path(name), out);
}

/// Rebuild the kernel a named reference block was dumped with.
InclusiveKernel build_kernel(const std::string& ion_name,
                             const std::string& block) {
  const Ion& ion = ion_by_name(ion_name);
  InclusiveKernel::Options opt;
  // C2: `target_mass` DEFAULTS to true (xsec.py:149), so the finite-gamma
  // block is built from the DEFAULT options -- this is the assertion that the
  // C++ default-constructed kernel IS the Python default-constructed one --
  // and the massless block has to ask for target_mass = false by name.
  if (block == "kernel_default_target_mass") {
    CHECK(opt.target_mass);
  } else {
    opt.target_mass = false;
  }
  if (block == "kernel_tensor_gamma") {
    // the EXACT finite-gamma tensor sector, with both higher-twist slots
    // filled by the dump's scenario shapes (`B3_SCENARIO`, `B4_SCENARIO`)
    CHECK(!opt.tensor_gamma);   // ... and it is OFF in a default kernel
    opt.tensor_gamma = true;
    opt.b3_func = [](double, double, double f1) { return 0.05 * f1; };
    opt.b4_func = [](double, double, double f1) { return -0.02 * f1; };
  }
  if (ion_name == "6Li") {
    // b1_func = toy_b1(..., mode="toy"), delta_func = toy_delta_gluon(scale=1e-3)
    opt.b1_func = [](double x, double q2, double f1) {
      return toy_b1(x, q2, f1, B1Mode::kToy);
    };
    opt.delta_func = [](double x, double q2, double f1) {
      return toy_delta_gluon(x, q2, f1, 1e-3);
    };
  } else if (block == "kernel_scenario_rank2" || block == "kernel_tensor_gamma") {
    opt.b1_32_func = [](double, double, double f1) { return 0.05 * f1; };
    opt.delta_32_func = [](double, double, double f1) { return -1e-2 * f1; };
  }
  return InclusiveKernel(ion, opt);
}

}  // namespace

TEST_CASE("reference tables: beams.json") {
  jsonmin::Value doc;
  if (!load("beams.json", doc)) {
    MESSAGE("validation/reference/beams.json absent -- skipping");
    return;
  }
  const jsonmin::Value& c = doc["constants"];
  CHECK_CLOSE(PROTON_TOP_MOMENTUM, c["PROTON_TOP_MOMENTUM"].num(), kRtol);
  CHECK_CLOSE(PROTON_MASS, c["PROTON_MASS"].num(), kRtol);
  for (std::size_t i = 0; i < 3; ++i) {
    CHECK(PROTON_CONFIG_ENERGIES[i] == c["PROTON_CONFIG_ENERGIES"][i].num());
    CHECK(ELECTRON_ENERGIES[i] == c["ELECTRON_ENERGIES"][i].num());
  }
  // the 6Li cluster wave function: the two D-state probabilities live in
  // beams.hpp (as they do in beams.py since 2026-08-29), and the inclusive
  // eff_pol slots below are LI6_CLUSTER_POLARIZATION/3 built from them
  CHECK(P_D_LI6 == c["P_D_LI6"].num());
  CHECK(P_D_DEUTERON == c["P_D_DEUTERON"].num());
  CHECK_CLOSE(ALPHA_D_VECTOR_POLARIZATION,
              c["ALPHA_D_VECTOR_POLARIZATION"].num(), kRtol);
  CHECK_CLOSE(DEUTERON_VECTOR_POLARIZATION,
              c["DEUTERON_VECTOR_POLARIZATION"].num(), kRtol);
  CHECK_CLOSE(LI6_CLUSTER_POLARIZATION,
              c["LI6_CLUSTER_POLARIZATION"].num(), kRtol);
  CHECK(LI6_NAIVE_ONE_THIRD == c["LI6_NAIVE_ONE_THIRD"].num());

  for (const auto& kv : c["NUCLEUS_MASS"].obj()) {
    // keys are "<name>_<A>_<Z>"
    const std::size_t p2 = kv.first.rfind('_');
    const std::size_t p1 = kv.first.rfind('_', p2 - 1);
    const std::string name = kv.first.substr(0, p1);
    const int a = std::atoi(kv.first.substr(p1 + 1, p2 - p1 - 1).c_str());
    const int z = std::atoi(kv.first.substr(p2 + 1).c_str());
    CHECK_CLOSE(nucleus_mass(name, a, z), kv.second.num(), kRtol);
  }

  for (const auto& kv : doc["ions"].obj()) {
    const Ion& ion = ion_by_name(kv.first);
    const jsonmin::Value& v = kv.second;
    CHECK(ion.A == static_cast<int>(v["A"].num()));
    CHECK(ion.Z == static_cast<int>(v["Z"].num()));
    CHECK(ion.N() == static_cast<int>(v["N"].num()));
    CHECK(ion.spin == v["spin"].num());
    CHECK_CLOSE(ion.eff_pol_p, v["eff_pol_p"].num(), kRtol);
    CHECK_CLOSE(ion.eff_pol_n, v["eff_pol_n"].num(), kRtol);
    CHECK_CLOSE(ion.mass_per_nucleon(), v["mass_per_nucleon"].num(), kRtol);
    CHECK_CLOSE(ion.momentum_per_nucleon_max(),
                v["momentum_per_nucleon_max"].num(), kRtol);
  }

  for (const auto& kv : doc["configs"].obj()) {
    const std::vector<BeamConfig> cfgs = default_configs(kv.first);
    REQUIRE(cfgs.size() == kv.second.size());
    for (std::size_t i = 0; i < cfgs.size(); ++i) {
      const jsonmin::Value& row = kv.second[i];
      CHECK(cfgs[i].electron_energy == row["electron_energy"].num());
      CHECK(cfgs[i].ion_momentum_per_nucleon
            == row["ion_momentum_per_nucleon"].num());
      CHECK_CLOSE(cfgs[i].sqrt_s_per_nucleon(),
                  row["sqrt_s_per_nucleon"].num(), kRtol);
      const double m = cfgs[i].ion.mass_per_nucleon();
      const double p = cfgs[i].ion_momentum_per_nucleon;
      CHECK_CLOSE(std::sqrt(p * p + m * m) / m, row["gamma_ion"].num(), kRtol);
      CHECK(cfgs[i].label() == row["label"].str());
    }
  }

  for (const jsonmin::Value& row : doc["gamma_of_proton_energy"].arr()) {
    CHECK_CLOSE(gamma_of(row["proton_energy"].num()), row["gamma"].num(), kRtol);
  }
}

TEST_CASE("reference tables: spin.json") {
  jsonmin::Value doc;
  if (!load("spin.json", doc)) {
    MESSAGE("validation/reference/spin.json absent -- skipping");
    return;
  }

  for (const jsonmin::Value& e : doc["wigner_d"].arr()) {
    const double j = e["j"].num();
    const RealMatrix d = wigner_d(j, e["beta"].num());
    const std::vector<double> ms = m_values(j);
    for (std::size_t i = 0; i < ms.size(); ++i) {
      CHECK(ms[i] == e["m_order"][i].num());
      for (std::size_t k = 0; k < ms.size(); ++k) {
        CHECK_CLOSE_AT(d(i, k), e["d"][i][k].num(), kRtol, 1e-15);
      }
    }
  }

  for (const jsonmin::Value& e : doc["clebsch_gordan"].arr()) {
    const jsonmin::Value& a = e["args"];
    const double got = clebsch_gordan(a[0].num(), a[1].num(), a[2].num(),
                                      a[3].num(), a[4].num(), a[5].num());
    CHECK_CLOSE_AT(got, e["value"].num(), kRtol, 1e-15);
  }

  for (const auto& kv : doc["spins"].obj()) {
    const jsonmin::Value& blk = kv.second;
    const double j = blk["j"].num();
    for (const jsonmin::Value& set : blk["population_sets"].arr()) {
      const std::string name = set["name"].str();
      std::vector<double> pops;
      for (const jsonmin::Value& v : set["populations"].arr()) {
        pops.push_back(v.num());
      }

      // the generators themselves must reproduce the population vector
      if (name.rfind("maxent_pz=", 0) == 0) {
        // The name carries pz at "%.10g", which is NOT the double that was
        // bisected on: 8/13 prints as 0.6153846154.  Recover the exact
        // argument by matching the two anchors against that formatting.
        const std::string tag = name.substr(10);
        double pz = std::atof(tag.c_str());
        for (double cand : {8.0 / 13.0, 0.7}) {
          char buf[32];
          std::snprintf(buf, sizeof(buf), "%.10g", cand);
          if (tag == buf) pz = cand;
        }
        const std::vector<double> got = populations_maxent(j, pz);
        REQUIRE(got.size() == pops.size());
        for (std::size_t i = 0; i < pops.size(); ++i) {
          CHECK_CLOSE_AT(got[i], pops[i], kMaxentRtol, 0.0);
        }
      } else if (name == "explicit_pz=0.7_pzz=0.4") {
        const std::array<double, 3> got = spin1_populations(0.7, 0.4);
        for (std::size_t i = 0; i < 3; ++i) CHECK_CLOSE(got[i], pops[i], kRtol);
      } else if (name == "explicit_pz=8/13_pzz=4/13") {
        const std::array<double, 3> got =
            spin1_populations(8.0 / 13.0, 4.0 / 13.0);
        for (std::size_t i = 0; i < 3; ++i) CHECK_CLOSE(got[i], pops[i], kRtol);
      } else if (name == "explicit_pz=0.7_t=0.4") {
        const std::array<double, 4> got = spin32_populations(0.7, 0.4, 0.0);
        for (std::size_t i = 0; i < 4; ++i) {
          CHECK_CLOSE_AT(got[i], pops[i], kRtol, 1e-15);
        }
      } else if (name == "explicit_pz=8/13_t=4/13") {
        const std::array<double, 4> got =
            spin32_populations(8.0 / 13.0, 4.0 / 13.0, 0.0);
        for (std::size_t i = 0; i < 4; ++i) {
          CHECK_CLOSE_AT(got[i], pops[i], kRtol, 1e-15);
        }
      }

      const AxisMoments own = moments_along_axis(j, pops);
      const jsonmin::Value& mo = set["moments_own_axis"];
      CHECK_CLOSE_AT(own.vector, mo["vector"].num(), kRtol, 1e-14);
      CHECK_CLOSE_AT(own.tensor, mo["tensor_zz"].num(), kRtol, 1e-14);
      if (mo.has("octupole_z")) {
        CHECK_CLOSE_AT(own.octupole, mo["octupole_z"].num(), kRtol, 1e-14);
      }

      for (const jsonmin::Value& ax : set["axes"].arr()) {
        const double theta = ax["theta"].num(), phi = ax["phi"].num();
        const SpinDensity dens{j, pops, theta, phi};
        const CplxMatrix rho = dens.rho_lab();
        for (std::size_t r = 0; r < rho.size(); ++r) {
          for (std::size_t c = 0; c < rho.size(); ++c) {
            CHECK_CLOSE_AT(rho(r, c).real(), ax["rho"]["re"][r][c].num(), kRtol,
                           1e-14);
            CHECK_CLOSE_AT(rho(r, c).imag(), ax["rho"]["im"][r][c].num(), kRtol,
                           1e-14);
          }
        }
        const std::array<double, 3> vp = vector_polarization(rho, j);
        for (std::size_t i = 0; i < 3; ++i) {
          CHECK_CLOSE_AT(vp[i], ax["vector_polarization"][i].num(), kRtol, 1e-14);
        }
        CHECK_CLOSE_AT(tensor_polarization(rho, j),
                       ax["tensor_polarization"].num(), kRtol, 1e-14);
        const SpinDensity::LabMoments lm = dens.lab_moments();
        for (std::size_t i = 0; i < 3; ++i) {
          CHECK_CLOSE_AT(lm.vector[i], ax["lab_moments"]["vector"][i].num(),
                         kRtol, 1e-14);
        }
        CHECK_CLOSE_AT(lm.tensor_zz, ax["lab_moments"]["tensor_zz"].num(), kRtol,
                       1e-14);
        if (ax.has("octupole_moment")) {
          CHECK_CLOSE_AT(octupole_moment(rho, j), ax["octupole_moment"].num(),
                         kRtol, 1e-14);
          CHECK_CLOSE_AT(lm.octupole_z, ax["lab_moments"]["octupole_z"].num(),
                         kRtol, 1e-14);
        }
      }
    }
  }
}

TEST_CASE("reference tables: xsec.json") {
  jsonmin::Value doc;
  if (!load("xsec.json", doc)) {
    MESSAGE("validation/reference/xsec.json absent -- skipping");
    return;
  }

  const jsonmin::Value& consts = doc["constants"];
  CHECK_CLOSE(ALPHA_EM, consts["ALPHA_EM"].num(), kRtol);
  CHECK_CLOSE(GEV2_TO_PB, consts["GEV2_TO_PB"].num(), kRtol);
  CHECK(TENSOR_LL_SIGN == consts["TENSOR_LL_SIGN"].num());
  CHECK(M_NUCLEON == consts["M_NUCLEON"].num());
  CHECK_CLOSE(LI6_CLUSTER_POLARIZATION,
              consts["LI6_CLUSTER_POLARIZATION"].num(), kRtol);

  std::vector<double> xs, q2s;
  for (const jsonmin::Value& p : doc["grid_points"].arr()) {
    xs.push_back(p["x"].num());
    q2s.push_back(p["q2"].num());
  }
  const std::size_t np = xs.size();

  for (const auto& kv : doc["results"].obj()) {
    const std::string ion_name = kv.first;
    const jsonmin::Value& entry = kv.second;
    const double s = entry["beam_config"]["s"].num();
    CHECK_CLOSE(default_configs(ion_name)[1].s_per_nucleon(), s, kRtol);

    for (const auto& bkv : entry.obj()) {
      if (bkv.first.rfind("kernel", 0) != 0) continue;
      const std::string block = bkv.first;
      const jsonmin::Value& b = bkv.second;
      const InclusiveKernel kern = build_kernel(ion_name, block);
      if (block == "kernel_default_target_mass") {
        CHECK(kern.target_mass());
        CHECK(kern.g2_scale() == 1.0);
        // ... and so is the kernel built with no Options at all.
        CHECK(InclusiveKernel(ion_by_name(ion_name)).target_mass());
      }

      std::vector<SFTables> tabs(np);
      for (std::size_t i = 0; i < np; ++i) {
        tabs[i] = kern.tables(xs[i], q2s[i], true);
      }

      // --- structure-function tables
      const jsonmin::Value& jt = b["tables"];
      for (std::size_t i = 0; i < np; ++i) {
        CHECK_CLOSE_AT(tabs[i].f1, jt["f1"][i].num(), kRtol, 1e-300);
        CHECK_CLOSE_AT(tabs[i].f2, jt["f2"][i].num(), kRtol, 1e-300);
        CHECK_CLOSE_AT(tabs[i].g1, jt["g1"][i].num(), kRtol, 1e-300);
        CHECK_CLOSE_AT(tabs[i].g2, jt["g2"][i].num(), kRtol, 1e-300);
        CHECK_CLOSE_AT(tabs[i].b1, jt["b1"][i].num(), kRtol, 1e-300);
        CHECK_CLOSE_AT(tabs[i].b2, jt["b2"][i].num(), kRtol, 1e-300);
        CHECK_CLOSE_AT(tabs[i].b3, jt["b3"][i].num(), kRtol, 1e-300);
        CHECK_CLOSE_AT(tabs[i].b4, jt["b4"][i].num(), kRtol, 1e-300);
        CHECK_CLOSE_AT(tabs[i].delta, jt["delta"][i].num(), kRtol, 1e-300);
      }

      // --- the analytic asymmetries the master formula must reproduce
      const jsonmin::Value& ja = b["asymmetries"];
      for (std::size_t i = 0; i < np; ++i) {
        const double y = q2s[i] / (s * xs[i]);
        CHECK_CLOSE(y, ja["y"][i].num(), kRtol);
        CHECK_CLOSE(r_sigma_lt(xs[i], q2s[i]), ja["R"][i].num(), kRtol);
        CHECK_CLOSE(depolarization_d(y, xs[i], q2s[i]),
                    ja["depolarization_d"][i].num(), kRtol);
        CHECK_CLOSE_AT(a_parallel(tabs[i].g1, tabs[i].f1, y, xs[i], q2s[i]),
                       ja["a_parallel"][i].num(), kRtol, 1e-300);
        CHECK_CLOSE_AT(azz(tabs[i].b1, tabs[i].f1, tabs[i].f2, xs[i], y,
                           &tabs[i].b2),
                       ja["azz"][i].num(), kRtol, 1e-300);
        CHECK_CLOSE_AT(a_cos2phi(tabs[i].delta, tabs[i].f1, tabs[i].f2, xs[i], y),
                       ja["a_cos2phi"][i].num(), kRtol, 1e-300);
      }

      // --- the modulation amplitudes, every (lam_e, axis, m)
      for (const jsonmin::Value& amp : b["amplitudes"].arr()) {
        EventSpinState st;
        st.lam_e = static_cast<int>(amp["lam_e"].num());
        st.pe = amp["pe"].num();
        st.j = kern.ion().spin;
        st.m = amp["m"].num();
        st.theta_s = amp["theta_s"].num();
        st.phi_s = amp["phi_s"].num();
        for (std::size_t i = 0; i < np; ++i) {
          const Amplitudes a =
              kern.amplitudes(tabs[i], xs[i], q2s[i], s, st, true);
          CHECK_CLOSE_AT(a.w_avg, amp["w_avg"][i].num(), kRtol, 1e-300);
          CHECK_CLOSE_AT(a.a1, amp["a1"][i].num(), kRtol, 1e-300);
          CHECK_CLOSE_AT(a.a2, amp["a2"][i].num(), kRtol, 1e-300);
        }
      }

      // --- the full differential cross section at eight azimuths
      for (const jsonmin::Value& ds : b["dsigma"].arr()) {
        EventSpinState st;
        st.lam_e = static_cast<int>(ds["lam_e"].num());
        st.pe = ds["pe"].num();
        st.j = kern.ion().spin;
        st.m = ds["m"].num();
        st.theta_s = ds["theta_s"].num();
        st.phi_s = ds["phi_s"].num();
        for (std::size_t i = 0; i < np; ++i) {
          for (std::size_t k = 0; k < ds["phi_grid"].size(); ++k) {
            const double phi = ds["phi_grid"][k].num();
            CHECK_CLOSE_AT(kern.dsigma(xs[i], q2s[i], phi, s, st, true),
                           ds["dsigma_per_point"][i][k].num(), kRtol, 1e-300);
          }
        }
      }

      for (std::size_t i = 0; i < np; ++i) {
        CHECK_CLOSE(kern.dsigma_unpol(xs[i], q2s[i], s),
                    b["dsigma_unpol"][i].num(), kRtol);
      }

      // --- the finite-gamma block
      if (b.has("target_mass")) {
        const jsonmin::Value& tm = b["target_mass"];
        for (std::size_t i = 0; i < np; ++i) {
          const double y = q2s[i] / (s * xs[i]);
          const double g2v = gamma_squared(xs[i], q2s[i]);
          const double r = r_sigma_lt(xs[i], q2s[i]);
          CHECK_CLOSE(g2v, tm["gamma_squared"][i].num(), kRtol);
          CHECK_CLOSE(epsilon_gamma(y, g2v), tm["epsilon_gamma"][i].num(), kRtol);
          CHECK_CLOSE(depolarization_gamma(y, g2v, r),
                      tm["depolarization_gamma"][i].num(), kRtol);
          CHECK_CLOSE(eta_gamma(y, g2v), tm["eta_gamma"][i].num(), kRtol);
          CHECK_CLOSE_AT(kern.a_parallel(tabs[i], xs[i], q2s[i], y),
                         tm["a_parallel_finite_gamma"][i].num(), kRtol, 1e-300);
          CHECK_CLOSE_AT(a_parallel(tabs[i].g1, tabs[i].f1, y, xs[i], q2s[i]),
                         tm["a_parallel_massless"][i].num(), kRtol, 1e-300);
        }
      }

      // --- the exact finite-gamma TENSOR sector (Cosyn, plans/08 D2)
      if (b.has("tensor_gamma")) {
        CHECK(kern.tensor_gamma());
        const jsonmin::Value& tg = b["tensor_gamma"];
        for (std::size_t i = 0; i < np; ++i) {
          const double y = q2s[i] / (s * xs[i]);
          const double g2v = gamma_squared(xs[i], q2s[i]);
          const std::pair<double, double> cs = theta_q_cos_sin(y, g2v);
          CHECK_CLOSE(cs.first, tg["cos_theta_q"][i].num(), kRtol);
          CHECK_CLOSE_AT(cs.second, tg["sin_theta_q"][i].num(), kRtol, 1e-300);
          const CosynTensorSFs f = cosyn_tensor_sfs(
              tabs[i].b1, tabs[i].b2, tabs[i].b3, tabs[i].b4, xs[i], g2v);
          CHECK_CLOSE_AT(f.f_t, tg["F_TLL_T"][i].num(), kRtol, 1e-300);
          CHECK_CLOSE_AT(f.f_l, tg["F_TLL_L"][i].num(), kRtol, 1e-300);
          CHECK_CLOSE_AT(f.f_lt, tg["F_TLT"][i].num(), kRtol, 1e-300);
          CHECK_CLOSE_AT(f.f_tt, tg["F_TTT"][i].num(), kRtol, 1e-300);
          const std::pair<double, double> fu =
              cosyn_unpolarized_sfs(tabs[i].f1, tabs[i].f2, xs[i], g2v);
          CHECK_CLOSE_AT(fu.first, tg["F_UU_T"][i].num(), kRtol, 1e-300);
          CHECK_CLOSE_AT(fu.second, tg["F_UU_L"][i].num(), kRtol, 1e-300);
        }
        for (const jsonmin::Value& hh : tg["harmonics"].arr()) {
          EventSpinState st;
          st.lam_e = 0;
          st.pe = 0.0;
          st.j = kern.ion().spin;
          st.m = hh["m"].num();
          st.theta_s = hh["theta_s"].num();
          st.phi_s = hh["phi_s"].num();
          for (std::size_t i = 0; i < np; ++i) {
            const double y = q2s[i] / (s * xs[i]);
            const TensorHarmonics h =
                kern.tensor_harmonics_gamma(tabs[i], xs[i], q2s[i], y, st);
            CHECK_CLOSE_AT(h.h0, hh["h0"][i].num(), kRtol, 1e-300);
            CHECK_CLOSE_AT(h.h1, hh["h1"][i].num(), kRtol, 1e-300);
            CHECK_CLOSE_AT(h.h2, hh["h2"][i].num(), kRtol, 1e-300);
          }
        }
      }
    }
  }
}
