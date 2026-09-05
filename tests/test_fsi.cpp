// SPDX-License-Identifier: GPL-3.0-or-later
// FSI as a per-event weight (fsi.hpp): the Glauber X-cluster distortion of
// the tagged spectator spectrum.
//
// The physics gates: the prototype FSI/IA table of
// docs/open_items/code_designs.md section 3 (fsi_alpha.py) at sigma_XN = 40
// and 20 mb on the prototype's own pure-S Hulthen wave, the production S+D
// row next to it so the difference is on the record, the X-alpha profile
// numbers (sigma_tot / sigma_el / B_alpha and the Glauber shadowing against
// the unshadowed A sigma_XN), the exact theta_k -> pi - theta_k symmetry, the
// sigma -> 0 limit, weight_normalised = weight / survival with unit mean over
// the model's own density, construction determinism, the wrong-spectator /
// wrong-channel guard, the sigma_XN(W) formation ramp, the 7Li P-wave
// channel, and -- through the Pipeline -- the one rule everything rests on:
// THE WEIGHT NEVER MOVES A FOUR-VECTOR.

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <memory>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "lipolgen/fsi.hpp"
#include "lipolgen/numerics.hpp"
#include "lipolgen/pipeline.hpp"
#include "lipolgen/tagged.hpp"

using namespace lipolgen;

namespace {

/// The weight at (k [GeV], theta_k [deg]) with the default (alpha) spectator.
double wt(const GlauberFsiWeight& f, double k, double theta_deg) {
  FsiKinematics kin;
  kin.k = k;
  kin.cos_theta_k = std::cos(theta_deg * kPi / 180.0);
  kin.spectator_z = 2;
  kin.spectator_a = 4;
  return f.weight(kin);
}

/// The PROTOTYPE's wave: pure S Hulthen at beta = 0.30 (p_d = 0), which is
/// exactly what fsi_alpha.py tabulated.  The production S+D channel is pinned
/// separately below.
const TaggedChannel& pure_s_channel() {
  static const TaggedChannel ch = li6_alpha_channel(0.30, 0.0);
  return ch;
}

const GlauberFsiWeight& f40() {
  static const GlauberFsiWeight f(pure_s_channel(), GlauberFsiOptions());
  return f;
}

const GlauberFsiWeight& f20() {
  static const GlauberFsiWeight f = [] {
    GlauberFsiOptions o;
    o.sigma_xn_mb = 20.0;
    return GlauberFsiWeight(pure_s_channel(), o);
  }();
  return f;
}

}  // namespace

// ------------------------------------------------- the prototype table pins

TEST_CASE("fsi: the prototype FSI/IA table and profile are reproduced") {
  const GlauberFsiWeight& f = f40();

  // The X-alpha profile at sigma_XN = 40 mb, eps = -0.5, B_XN = 6 GeV^-2,
  // a^2 = 0.700 fm^2: Glauber shadowing gives 131.0 mb, NOT 4 x 40 = 160.
  CHECK_CLOSE_AT(f.sigma_cluster_mb(40.0), 131.0, 0.0, 0.1);
  CHECK_CLOSE_AT(f.sigma_cluster_el_mb(40.0), 35.2, 0.0, 0.1);
  CHECK_CLOSE_AT(f.slope_cluster_gev2(40.0), 27.2, 0.0, 0.1);
  CHECK_CLOSE_AT(f.sigma_cluster_mb(20.0), 72.5, 0.0, 0.1);
  // sigma_el / sigma_tot = 0.269 against the measured alpha-p 0.258
  // (Blinov et al.) -- the split is automatic, not a knob.
  CHECK_CLOSE_AT(f.sigma_cluster_el_mb(40.0) / f.sigma_cluster_mb(40.0),
                 0.269, 0.0, 2e-3);

  // The pinned FSI/IA rows (fsi_alpha.py; measured C++ deviations <= 3.3e-3).
  struct Row { double k, theta_deg, want; };
  const Row rows40[] = {
      {0.02, 0.0, 0.787},
      {0.10, 0.0, 0.607},
      {0.10, 90.0, 0.381},
      {0.20, 0.0, 0.398},
  };
  for (const Row& r : rows40) {
    CAPTURE(r.k);
    CAPTURE(r.theta_deg);
    CHECK_CLOSE_AT(wt(f, r.k, r.theta_deg), r.want, 0.0, 5e-3);
  }
  const Row rows20[] = {
      {0.02, 0.0, 0.877}, {0.05, 0.0, 0.844}, {0.10, 0.0, 0.764},
  };
  for (const Row& r : rows20) {
    CAPTURE(r.k);
    CHECK_CLOSE_AT(wt(f20(), r.k, r.theta_deg), r.want, 0.0, 5e-3);
  }

  MESSAGE("40 mb: survival = " << f.survival()
          << ", 20 mb: survival = " << f20().survival()
          << ", clipped fraction = " << f.clipped_grid_fraction());
}

TEST_CASE("fsi: the production S+D channel row is on the record") {
  // The default 6Li alpha channel (P_D = 0.0867) SHIFTS the pinned numbers;
  // this row is what a production run at sigma_XN = 40 mb actually applies,
  // pinned so the S-vs-S+D difference is loud rather than folklore.
  const TaggedChannel ch = li6_alpha_channel();
  const GlauberFsiWeight f(ch, GlauberFsiOptions());
  CHECK_CLOSE_AT(wt(f, 0.10, 0.0), 0.617, 0.0, 5e-3);
  CHECK_CLOSE_AT(wt(f, 0.10, 90.0), 0.391, 0.0, 5e-3);
  CHECK_CLOSE_AT(wt(f, 0.20, 0.0), 0.446, 0.0, 5e-3);
}

// ----------------------------------------------------- structural identities

TEST_CASE("fsi: the weight is EXACTLY symmetric under theta_k -> pi - theta_k") {
  // The residue approximation drops the principal-value piece, so the kernel
  // depends on |k_z| alone and the symmetry is exact by construction (the
  // grid is read at |k c|) -- bitwise, not approximately.
  const GlauberFsiWeight& f = f40();
  for (double k : {0.02, 0.11, 0.23, 0.47, 0.9}) {
    for (double c : {0.15, 0.5, 0.77, 0.99}) {
      FsiKinematics kin;
      kin.k = k;
      kin.cos_theta_k = c;
      double up = f.weight(kin);
      kin.cos_theta_k = -c;
      CHECK(up == f.weight(kin));
    }
  }
}

TEST_CASE("fsi: the FSI/IA ratio -> 1 as sigma_XN -> 0") {
  GlauberFsiOptions o;
  o.sigma_xn_mb = 1e-6;
  const GlauberFsiWeight f(pure_s_channel(), o);
  for (double k : {0.02, 0.1, 0.2, 0.5, 1.0}) {
    for (double theta : {0.0, 30.0, 60.0, 90.0}) {
      CAPTURE(k);
      CAPTURE(theta);
      CHECK_CLOSE_AT(wt(f, k, theta), 1.0, 0.0, 1e-4);
    }
  }
  CHECK_CLOSE_AT(f.survival(), 1.0, 0.0, 1e-4);
}

TEST_CASE("fsi: grid weights are non-negative, bounded and unclipped") {
  const GlauberFsiWeight& f = f40();
  const double w_max = f.options().w_max;
  for (double v : f.grid(0)) {
    REQUIRE(v >= 0.0);
    REQUIRE(v <= w_max);
  }
  // At the default grid nothing reaches the ceiling: the > 1 region past the
  // amplitude sign flip at k_T ~ 0.25 GeV stays far below w_max = 50.
  CHECK(f.clipped_grid_fraction() == 0.0);
}

TEST_CASE("fsi: weight_normalised = weight / survival and averages to 1") {
  const GlauberFsiWeight& f = f40();

  // The identity, pointwise (one ladder layer: survival is W-independent).
  for (double k : {0.05, 0.1, 0.3}) {
    for (double c : {0.0, 0.6, -0.9}) {
      FsiKinematics kin;
      kin.k = k;
      kin.cos_theta_k = c;
      CHECK_CLOSE(f.weight_normalised(kin), f.weight(kin) / f.survival(),
                  1e-12);
    }
  }

  // Unit mean over the model's OWN density n(k) k^2 dk dc, on an INDEPENDENT
  // midpoint quadrature (different resolution from the survival's internal
  // 240 x 96 one, so this is a statement about the integral, not an echo).
  const TaggedChannel& ch = pure_s_channel();
  const double kappa = ch.base.kappa();
  const std::vector<double> kn = linspace(1e-4, 1.2, 280);
  // The unpolarized density n(k) = sum_L P_L psihat_L(k)^2 from the channel's
  // own waves, normalized on the SAME grid the weight's build uses.
  std::vector<const Wave*> waves;
  std::vector<double> inv_norm2;
  for (const Wave& wv : ch.waves) {
    if (!(wv.prob > 0.0)) continue;
    const std::vector<double> psi_n = wv.radial(kn, kappa);
    std::vector<double> integ(kn.size());
    for (std::size_t i = 0; i < kn.size(); ++i) {
      integ[i] = psi_n[i] * psi_n[i] * kn[i] * kn[i];
    }
    const double norm2 = trapezoid(integ, kn);
    REQUIRE(norm2 > 0.0);
    waves.push_back(&wv);
    inv_norm2.push_back(1.0 / norm2);
  }

  const std::size_t nk = 193, nc = 77;
  double sw = 0.0, sn = 0.0;
  for (std::size_t ik = 0; ik < nk; ++ik) {
    const double k = (ik + 0.5) * 1.2 / static_cast<double>(nk);
    double nk_dens = 0.0;
    for (std::size_t wi = 0; wi < waves.size(); ++wi) {
      const double psi = waves[wi]->radial(k, kappa);
      nk_dens += waves[wi]->prob * psi * psi * inv_norm2[wi];
    }
    const double meas = nk_dens * k * k;
    for (std::size_t ic = 0; ic < nc; ++ic) {
      FsiKinematics kin;
      kin.k = k;
      kin.cos_theta_k = -1.0 + (ic + 0.5) * 2.0 / static_cast<double>(nc);
      sw += meas * f.weight_normalised(kin);
      sn += meas;
    }
  }
  MESSAGE("independent quadrature of <weight_normalised> = " << sw / sn);
  CHECK_CLOSE_AT(sw / sn, 1.0, 0.0, 1e-3);
}

TEST_CASE("fsi: construction is deterministic") {
  // Two independently built weights are BIT-identical: nothing in the build
  // touches an RNG or global state.
  const GlauberFsiWeight a(pure_s_channel(), GlauberFsiOptions());
  const GlauberFsiWeight b(pure_s_channel(), GlauberFsiOptions());
  const std::vector<double>& ga = a.grid(0);
  const std::vector<double>& gb = b.grid(0);
  REQUIRE(ga.size() == gb.size());
  for (std::size_t i = 0; i < ga.size(); ++i) REQUIRE(ga[i] == gb[i]);
  CHECK(a.survival() == b.survival());
}

TEST_CASE("fsi: a wrong spectator or wrong channel is left undistorted") {
  const GlauberFsiWeight& f = f40();
  FsiKinematics kin;
  kin.k = 0.10;
  kin.cos_theta_k = 0.0;
  // A deuteron spectator through an alpha-built weight: EXACTLY 1.
  kin.spectator_z = 1;
  kin.spectator_a = 2;
  CHECK(f.weight(kin) == 1.0);
  // Right spectator, wrong channel label: EXACTLY 1.
  kin.spectator_z = 2;
  kin.spectator_a = 4;
  kin.channel = Channel::TaggedDeuteronP;
  CHECK(f.weight(kin) == 1.0);
  // The matching channel (and the "unspecified" default) distort.
  kin.channel = Channel::TaggedLi6Alpha;
  CHECK(f.weight(kin) < 1.0);
  kin.channel = Channel::Inclusive;
  CHECK(f.weight(kin) < 1.0);
}

// ------------------------------------------------------------- the variants

TEST_CASE("fsi: GlauberNucleon is the POINTWISE bracket, not the integrated one") {
  GlauberFsiOptions o;
  o.variant = FsiVariant::GlauberNucleon;
  const GlauberFsiWeight fn(pure_s_channel(), o);
  const GlauberFsiWeight& fc = f40();

  // Unshadowed single-scattering limit: sigma_Xa = A sigma_XN = 160 mb, the
  // absorption bound ABOVE the shadowed 131.0.
  CHECK_CLOSE_AT(fn.sigma_cluster_mb(40.0), 160.0, 0.0, 0.1);

  // Pointwise at low k it IS more absorptive than the cluster variant ...
  CHECK_CLOSE_AT(wt(fn, 0.10, 90.0), 0.309, 0.0, 5e-3);
  CHECK(wt(fn, 0.10, 90.0) < wt(fc, 0.10, 90.0));
  CHECK(wt(fn, 0.10, 0.0) < wt(fc, 0.10, 0.0));

  // ... but the INTEGRATED survival comes out ABOVE the shadowed default:
  // the unshadowed quadratic (gain) term scales by A^2 -- second order in
  // Gamma_N, formally beyond single scattering -- and feeds strength back
  // into the tag (fsi.hpp header).  Pinned so a change of behaviour is loud.
  CHECK_CLOSE_AT(fn.survival(), 0.582, 0.0, 1e-2);
  CHECK_CLOSE_AT(fc.survival(), 0.517, 0.0, 1e-2);
  CHECK(fn.survival() > fc.survival());
  MESSAGE("survival: nucleon " << fn.survival() << " vs cluster "
                               << fc.survival());
}

// D3.  The inventory carried "the per-nucleon Glauber FSI variant is
// algebraically identical to the cluster one -- the two 'variants' are one".
// It is not: (a) shadows and (b) does not, and the identity that DOES hold
// (the Ciofi-Kaptari per-nucleon product on an uncorrelated density collapses
// onto (a)) is about a code path this file does not have.  Pinned on GENERATED
// EVENTS, not on the profile alone, so the claim cannot come back.
TEST_CASE("fsi: the two variants differ on essentially every event") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::TaggedLi6Alpha;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = 4000;
  cfg.seed = 1234;
  cfg.grid.nx = 40;
  cfg.grid.nq2 = 28;
  const RunPlan plan = tensor_thirds_plan(0.7, 0.6);

  PipelineConfig cfg_c = cfg, cfg_n = cfg;
  cfg_c.fsi = PipelineFsi::GlauberCluster;
  cfg_n.fsi = PipelineFsi::GlauberNucleon;
  const Pipeline off(cfg, plan), pc(cfg_c, plan), pn(cfg_n, plan);

  // The profiles: shadowed 131.0 mb against the unshadowed A sigma_XN = 160.0,
  // and the elastic (gain) term differs by nearly a factor two.
  const GlauberFsiWeight& wc = *pc.fsi_weight();
  const GlauberFsiWeight& wn = *pn.fsi_weight();
  CHECK_CLOSE_AT(wc.sigma_cluster_mb(40.0), 131.045, 0.0, 1e-2);
  CHECK_CLOSE_AT(wn.sigma_cluster_mb(40.0), 160.006, 0.0, 1e-2);
  CHECK_CLOSE_AT(wc.sigma_cluster_el_mb(40.0), 35.224, 0.0, 1e-2);
  CHECK_CLOSE_AT(wn.sigma_cluster_el_mb(40.0), 68.186, 0.0, 1e-2);
  CHECK_CLOSE_AT(wc.survival(), 0.520239, 0.0, 1e-5);
  CHECK_CLOSE_AT(wn.survival(), 0.582899, 0.0, 1e-5);

  // ONE event stream reweighted three ways: the kinematics are bit-identical
  // (the weight moves nothing), so the weight difference is the whole of it.
  Event ea, eb, ec;
  std::uint64_t n = 0, n_differ = 0;
  double sum_off = 0.0, sum_c = 0.0, sum_n = 0.0;
  double r_lo = 1e300, r_hi = 0.0;
  for (std::uint64_t i = 0; i < off.size(); ++i) {
    off.event(i, ea);
    pc.event(i, eb);
    pn.event(i, ec);
    REQUIRE(ea.kin.k == eb.kin.k);
    REQUIRE(ea.kin.k == ec.kin.k);
    REQUIRE(ea.kin.cos_theta_k == ec.kin.cos_theta_k);
    REQUIRE(ea.kin.x == ec.kin.x);
    REQUIRE(ea.kin.q2 == ec.kin.q2);
    sum_off += ea.weight;
    sum_c += eb.weight;
    sum_n += ec.weight;
    const double r = ec.weight / eb.weight;
    r_lo = std::min(r_lo, r);
    r_hi = std::max(r_hi, r);
    if (std::fabs(r - 1.0) > 0.01) ++n_differ;
    ++n;
  }
  REQUIRE(n == 4000);
  MESSAGE("sum w: off " << sum_off << ", cluster " << sum_c << ", nucleon "
                        << sum_n << "; ratio in [" << r_lo << ", " << r_hi
                        << "], differing by >1 % on "
                        << 100.0 * static_cast<double>(n_differ) / n << " %");
  // The integrated survivals bracket in the direction the header says they do
  // (the unshadowed quadratic term feeds strength BACK into the tag), and the
  // total rate differs by tens of percent, not by rounding.
  CHECK(sum_c / sum_off < 0.56);
  CHECK(sum_n / sum_off > 0.56);
  CHECK(sum_n > 1.1 * sum_c);
  // ... and the per-event ratio spans more than an order of magnitude, with
  // essentially every event moved.  Measured through the pipeline at
  // tagged-6Li-alpha, 20 000 events, seed 1234, tensor-thirds at
  // pz = 0.7 / pzz = 0.6 / pe = 0.7, sigma_XN = 40 mb (the default, one end
  // of the mandatory 20-40 mb band): 99.50 % differ by more than 1 % and the
  // ratio reaches 68.52.  This case reproduces the SHAPE of that at reduced
  // statistics; the quoted percentages belong to the 20 000-event run.
  CHECK(r_lo < 0.9);
  CHECK(r_hi > 5.0);
  CHECK(static_cast<double>(n_differ) / n > 0.98);
}

TEST_CASE("fsi: the formation ramp interpolates sigma_XN in W") {
  GlauberFsiOptions o;
  o.formation_ramp = true;   // defaults: 40 mb at W <= 2, 20 mb at W >= 10
  const GlauberFsiWeight f(pure_s_channel(), o);
  CHECK(f.sigma_eff_mb(2.0) == 40.0);
  CHECK_CLOSE(f.sigma_eff_mb(6.0), 30.0, 1e-12);
  CHECK(f.sigma_eff_mb(10.0) == 20.0);
  CHECK(f.sigma_eff_mb(0.5) == 40.0);    // clamped below
  CHECK(f.sigma_eff_mb(50.0) == 20.0);   // clamped above
  CHECK(f.sigma_ladder().size() == 5);
  // At the ladder ends the layers are built at exactly 20 / 40 mb, so the
  // ramped survival closes on the single-sigma builds.
  CHECK_CLOSE(f.survival(10.0), f20().survival(), 1e-9);
  CHECK_CLOSE(f.survival(2.0), f40().survival(), 1e-9);
  CHECK_CLOSE_AT(f.survival(10.0), 0.668, 0.0, 1e-2);
  // In between, the weight interpolates monotonically between the rows.
  FsiKinematics kin;
  kin.k = 0.10;
  kin.cos_theta_k = 1.0;
  kin.w = 6.0;
  const double w_mid = f.weight(kin);
  CHECK(w_mid > wt(f40(), 0.10, 0.0));
  CHECK(w_mid < wt(f20(), 0.10, 0.0));
}

TEST_CASE("fsi: the 7Li P-wave channel builds and stays physical") {
  // Pure P wave: the m = +-1 azimuthal component exercises the order-1
  // Hankel pair, which no 6Li (S+D, m = 0 dominated) build reaches at m = 1.
  const TaggedChannel ch = li7_alpha_channel();
  const GlauberFsiWeight f(ch, GlauberFsiOptions());
  for (double v : f.grid(0)) {
    REQUIRE(std::isfinite(v));
    REQUIRE(v >= 0.0);
  }
  CHECK(f.survival() > 0.0);
  CHECK(f.survival() < 1.0);
  FsiKinematics kin;
  kin.k = 0.10;
  kin.cos_theta_k = 0.0;
  kin.channel = Channel::TaggedLi7Alpha;
  CHECK(f.weight(kin) > 0.0);
  CHECK(f.weight(kin) < 1.0);
  MESSAGE("7Li alpha (P wave): w(0.1, 90 deg) = " << f.weight(kin)
          << ", survival = " << f.survival());
}

// ------------------------------------------------------------- the pipeline

TEST_CASE("fsi: pipeline weights never move a four-vector") {
  // THE rule of fsi.hpp, stated on generated events: an FSI-on run is
  // BIT-IDENTICAL to the FSI-off run in every particle four-vector, every
  // residual and every kinematic label -- only Event::weight differs.  The
  // weight consumes no RNG, so this is exact, not statistical.
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::TaggedLi6Alpha;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = 4000;
  cfg.grid.nx = 40;
  cfg.grid.nq2 = 28;

  // One UNPOLARIZED category: every ion projection equally filled and pe = 0,
  // so the drawn spectator density is the isotropic unpolarized one and
  // <Event::weight> estimates exactly the survival the model logs.
  SpinCategory cat("unpol", 1.0, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}, 0, 0.0,
                   0.0, 0.0, 1.0);
  const RunPlan plan(std::vector<SpinCategory>{cat}, 0.0, 1.0, 1.0);

  const Pipeline off(cfg, plan);
  PipelineConfig cfg_on = cfg;
  cfg_on.fsi = PipelineFsi::GlauberCluster;
  const Pipeline on(cfg_on, plan);

  REQUIRE(off.fsi_weight() == nullptr);
  REQUIRE(on.fsi_weight() != nullptr);
  const GlauberFsiWeight& fw = *on.fsi_weight();
  CHECK(fw.options().variant == FsiVariant::GlauberCluster);
  CHECK(fw.options().sigma_xn_mb == 40.0);

  const double survival = fw.survival();
  double sum_w = 0.0;
  std::uint64_t n = 0, n_tag = 0;
  const int nbin = 6;
  const double k_edges[nbin + 1] = {0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 1.2};
  double bin_w[nbin] = {0.0};
  std::uint64_t bin_n[nbin] = {0};

  Event ea, eb;
  for (std::uint64_t i = 0; i < off.size(); ++i) {
    off.event(i, ea);
    on.event(i, eb);
    // Same particles, bit for bit -- the weight moved NOTHING.
    REQUIRE(ea.particles.size() == eb.particles.size());
    for (std::size_t j = 0; j < ea.particles.size(); ++j) {
      REQUIRE(ea.particles[j].p.e == eb.particles[j].p.e);
      REQUIRE(ea.particles[j].p.px == eb.particles[j].p.px);
      REQUIRE(ea.particles[j].p.py == eb.particles[j].p.py);
      REQUIRE(ea.particles[j].p.pz == eb.particles[j].p.pz);
      REQUIRE(ea.particles[j].pdg == eb.particles[j].pdg);
    }
    REQUIRE(ea.kin.k == eb.kin.k);
    REQUIRE(ea.kin.cos_theta_k == eb.kin.cos_theta_k);
    // ... including the conservation residual, event by event.
    const Vec4 ra = momentum_residual(ea), rb = momentum_residual(eb);
    REQUIRE(ra.e == rb.e);
    REQUIRE(ra.pz == rb.pz);

    // Off carries weight 1; on carries the model's own weight at the drawn
    // spectator, reproducible from the finished record.
    REQUIRE(ea.weight == 1.0);
    REQUIRE(eb.weight >= 0.0);
    FsiKinematics kin;
    kin.k = eb.kin.k;
    kin.cos_theta_k = eb.kin.cos_theta_k;
    kin.channel = eb.channel;
    REQUIRE(eb.weight == fw.weight(kin));

    sum_w += eb.weight;
    ++n;
    if (rp_tagged(eb, on.optics(), on.pot_config())) ++n_tag;
    for (int b = 0; b < nbin; ++b) {
      if (eb.kin.k >= k_edges[b] && eb.kin.k < k_edges[b + 1]) {
        bin_w[b] += eb.weight;
        ++bin_n[b];
        break;
      }
    }
  }
  REQUIRE(n == 4000);

  // <w> estimates the logged survival probability (se ~ 0.006 at n = 4000).
  const double mean_w = sum_w / static_cast<double>(n);
  MESSAGE("tag fraction at " << on.optics().name << " = "
          << static_cast<double>(n_tag) / static_cast<double>(n)
          << "; <w> = " << mean_w << " vs survival() = " << survival);
  CHECK_CLOSE_AT(mean_w, survival, 0.0, 0.02);
  for (int b = 0; b < nbin; ++b) {
    char line[128];
    std::snprintf(line, sizeof(line),
                  "  k in [%.2f, %.2f): <w> = %.3f  (n = %llu)", k_edges[b],
                  k_edges[b + 1], bin_n[b] ? bin_w[b] / bin_n[b] : 1.0,
                  static_cast<unsigned long long>(bin_n[b]));
    MESSAGE(line);
  }
}

TEST_CASE("fsi: the configuration refuses FSI off the tagged channels") {
  PipelineConfig cfg;
  cfg.n_events = 10;
  cfg.fsi = PipelineFsi::GlauberCluster;
  cfg.channel = PipelineChannel::Inclusive;
  CHECK_THROWS_AS(cfg.validate(), std::runtime_error);
  cfg.channel = PipelineChannel::CoherentLi6;
  CHECK_THROWS_AS(cfg.validate(), std::runtime_error);
  cfg.channel = PipelineChannel::TaggedLi6Alpha;
  cfg.isotope = "6Li";
  CHECK_NOTHROW(cfg.validate());
  cfg.fsi_sigma_mb = -1.0;
  CHECK_THROWS_AS(cfg.validate(), std::runtime_error);
}
