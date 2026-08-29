// Integration layer: the `Pipeline` end-to-end gates (plans/05 step 5.D).
//
// Per-event four-momentum and charge conservation for every channel, the
// spectator spectrum against the TaggedModel marginals, the Roman-Pot tag
// fractions against the published Python numbers, the acceptance-weighted
// A_zz^tag closure (the C++ shape of
// evgen/tests/test_tagged.py::test_tagged_rate_asymmetry_matches_analytic plus
// evgen/scripts/money_tagged_azz.py's folded closure), the 7Li
// <P2(cos theta_k)> = -T/5 polarimeter, the coherent <|t|> and tag acceptance,
// bit-identity of the inclusive path with InclusiveGenerator::run_n, thread
// determinism and a HepMC3 round trip.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "lipolgen/asymmetries.hpp"
#include "lipolgen/pipeline.hpp"

#if __has_include(<HepMC3/ReaderAscii.h>)
#define LIPOLGEN_TEST_HEPMC 1
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/ReaderAscii.h>
#include "lipolgen/hepmc_writer.hpp"
#endif

using namespace lipolgen;

namespace {

constexpr double kPZ = 0.7;
constexpr double kPZZ = 0.6;

/// One built pipeline per (channel, config, optics, n_events) key: the tagged
/// model's amplitude tables and the sampler grids cost more than the events.
struct Key {
  PipelineChannel ch;
  int config;
  OpticsChoice optics;
  std::uint64_t n;
  int plan_id;
  bool operator<(const Key& o) const {
    if (ch != o.ch) return ch < o.ch;
    if (config != o.config) return config < o.config;
    if (optics != o.optics) return optics < o.optics;
    if (n != o.n) return n < o.n;
    return plan_id < o.plan_id;
  }
};

RunPlan plan_of(int plan_id) {
  switch (plan_id) {
    case 0: return tensor_thirds_plan(kPZ, kPZZ);
    case 1: return helicity_flip_plan(1.5, kPZ, kPZZ);
    case 2: {  // one PURE ion state M = +1, unpolarized electron
      SpinCategory c("pure+1", 1.0, {1.0, 0.0, 0.0}, 0, 0.0, 0.0, 0.0, 1.0);
      return RunPlan({c}, 0.0, 1.0, 1.0);
    }
    default: break;
  }
  throw std::runtime_error("unknown plan id");
}

const Pipeline& pipe(PipelineChannel ch, int config, OpticsChoice optics,
                     std::uint64_t n, int plan_id) {
  static std::map<Key, std::shared_ptr<Pipeline>> cache;
  const Key key{ch, config, optics, n, plan_id};
  const auto it = cache.find(key);
  if (it != cache.end()) return *it->second;
  PipelineConfig cfg;
  cfg.channel = ch;
  const std::string iso = channel_isotope(ch);
  cfg.isotope = iso.empty() ? "6Li" : iso;
  cfg.beam_config = config;
  cfg.optics_choice = optics;
  cfg.n_events = n;
  // A coarser grid than the 100x72 default: these tests are about the
  // integration layer, and the grid is P3's gate.
  cfg.grid.nx = 40;
  cfg.grid.nq2 = 28;
  auto p = std::make_shared<Pipeline>(cfg, plan_of(plan_id));
  cache.emplace(key, p);
  return *p;
}

/// Worst |residual| / |p_in| over the four components.
double relative_residual(const Event& ev) {
  const Vec4 r = momentum_residual(ev);
  const Vec4 s = momentum_scale(ev);
  const double scale = std::fabs(s.e);
  return std::max(std::max(std::fabs(r.e), std::fabs(r.px)),
                  std::max(std::fabs(r.py), std::fabs(r.pz))) / scale;
}

double p2(double c) { return 0.5 * (3.0 * c * c - 1.0); }

}  // namespace

// ------------------------------------------------- conservation, all channels

TEST_CASE("pipeline: every channel conserves four-momentum and charge") {
  struct Row { PipelineChannel ch; int cfg; const char* what; int plan; };
  const Row rows[] = {
      {PipelineChannel::Inclusive, 1, "inclusive 6Li (per-nucleon balance)", 0},
      {PipelineChannel::TaggedLi6Alpha, 1, "6Li alpha tag", 0},
      {PipelineChannel::TaggedLi7Alpha, 1, "7Li alpha tag", 1},
      {PipelineChannel::TaggedDeuteronP, 1, "d(e,e'p) control", 0},
      {PipelineChannel::CoherentLi6, 1, "coherent 6Li", 0},
  };
  for (const Row& row : rows) {
    const std::string what = row.what;
    CAPTURE(what);
    const Pipeline& p = pipe(row.ch, row.cfg,
                             OpticsChoice::YellowReportHighAcceptance, 20000,
                             row.plan);
    double worst_p = 0.0, worst_q = 0.0;
    std::uint64_t n = 0;
    p.for_each([&](const Event& ev) {
      ++n;
      worst_p = std::max(worst_p, relative_residual(ev));
      worst_q = std::max(worst_q, std::fabs(charge_residual(ev)));
      // structural: the beams are first and the record is self-consistent
      REQUIRE(ev.particles.size() >= 4);
      CHECK(ev.particles[0].role == Role::BeamElectron);
      CHECK(ev.particles[1].role == Role::BeamIon);
      REQUIRE(ev.find(Role::ScatteredElectron) != nullptr);
      REQUIRE(ev.find(Role::HadronicX) != nullptr);
      // C1: X is the whole hadronic final state and must be TIMELIKE on
      // EVERY channel.  `Pipeline::add_hadronic_x` used to clip a spacelike
      // residual to a massless X; it throws now, so this is a second, direct
      // statement of the same invariant.
      const Particle* px = ev.find(Role::HadronicX);
      CHECK(px->p.m2() >= 0.0);
      CHECK(px->mass >= 0.0);
      CHECK_CLOSE_AT(px->mass * px->mass, px->p.m2(), 1e-9, 1e-12);
      CHECK(ev.kin.q2 > 0.0);
      CHECK(ev.kin.y > 0.0);
      CHECK(ev.kin.y < 1.0);
    });
    CHECK(n == 20000);
    MESSAGE(what << ": worst relative |dP| = " << worst_p
                     << ", worst |dQ| = " << worst_q);
    CHECK(worst_p < 1e-9);
    CHECK(worst_q < 1e-12);
  }
}

TEST_CASE("pipeline: a tagged event carries the off-shell struck cluster T2 needs") {
  const Pipeline& p = pipe(PipelineChannel::TaggedLi6Alpha, 1,
                           OpticsChoice::YellowReportHighAcceptance, 5000, 0);
  REQUIRE(p.tagged_channel() != nullptr);
  Event ev = p.event(0);
  const Particle* spec = ev.find(Role::Spectator);
  const Particle* clus = ev.find(Role::StruckCluster);
  REQUIRE(spec != nullptr);
  REQUIRE(clus != nullptr);
  CHECK(spec->pdg == 1000020040);        // alpha
  CHECK(clus->pdg == 1000010020);        // embedded deuteron
  CHECK(clus->status == Status::Intermediate);
  CHECK(spec->status == Status::Final);
  // The spectator is on shell and the struck cluster is not.  The tolerance
  // is 1e-9, not 1e-12: E^2 - p^2 at E = 400 GeV against m^2 = 14 GeV^2 is a
  // cancellation of five orders of magnitude, so double precision leaves
  // ~1e-12 of relative accuracy in m^2 no matter how exactly the boost was
  // done (the boost itself is checked below, on the four-vectors).
  CHECK_CLOSE(spec->p.m2(), spec->mass * spec->mass, 1e-9);

  StruckCluster sc;
  REQUIRE(struck_cluster_of(ev, *p.tagged_channel(), sc));
  CHECK(sc.virtuality < 0.0);            // impulse approximation
  CHECK(sc.alpha_s > 0.0);
  CHECK(sc.alpha_s < 6.0);
  CHECK_CLOSE(sc.alpha_s + sc.alpha_x, 6.0, 1e-12);
  CHECK(sc.p_per_nucleon_eff > 0.0);
  // P_ion - p_spec is the record's own struck cluster, to the last bit
  const Vec4 d = (ev.particles[1].p - spec->p) - clus->p;
  CHECK(std::fabs(d.e) < 1e-9);
  CHECK(std::fabs(d.pz) < 1e-9);
  MESSAGE("6Li alpha tag, event 0: M_X^2 = " << sc.m2 << " GeV^2, virtuality "
          << sc.virtuality << " GeV^2, alpha_s = " << sc.alpha_s
          << ", p_per_nucleon_eff = " << sc.p_per_nucleon_eff << " GeV");
}

// ------------------------------------------------------ spectator spectrum

TEST_CASE("pipeline: the tagged spectator spectrum is the TaggedModel marginal") {
  // ONE PURE ion state and an unpolarized electron, so the per-(M, m_S) rates
  // are pop_M * P(m_S|M) exactly (the struck-cluster DIS cross section carries
  // no m_S dependence without an inclusive b1 and without a helicity term) and
  // the k, cos(theta_k) marginals of the generated sample are the model's own
  // n_M(k, c) k^2 to the last cell.
  const std::size_t n_ev = 200000;
  const Pipeline& p = pipe(PipelineChannel::TaggedLi6Alpha, 1,
                           OpticsChoice::YellowReportHighAcceptance, n_ev, 2);
  const TaggedModel& m = *p.tagged_model();
  const std::vector<double>& n1 = m.n_of_kc(1.0);

  // The sampler draws a (k, c) CELL by inverse CDF and then places the event
  // UNIFORMLY inside it (`TaggedModel::sample_kc_one`: k_[ik] +- dk/2,
  // c_[ic] +- dc/2), so the marginal is exact cell by cell and the histogram
  // bins have to be unions of whole cells or the comparison measures the
  // binning, not the sampler.  280 k cells -> 28 bins of 10, 96 c cells ->
  // 16 bins of 6.
  const int kper = 10, cper = 6;
  const int nkb = static_cast<int>(m.nk()) / kper;
  const int ncb = static_cast<int>(m.nc()) / cper;
  REQUIRE(static_cast<std::size_t>(nkb * kper) == m.nk());
  REQUIRE(static_cast<std::size_t>(ncb * cper) == m.nc());
  const double k_edge0 = m.k()[0] - 0.5 * m.dk();
  const double kbw = kper * m.dk();

  std::vector<double> hk(nkb, 0.0), hc(ncb, 0.0);
  std::uint64_t n = 0;
  p.for_each([&](const Event& ev) {
    ++n;
    // k is folded through fabs at the bottom cell, which is why the lowest
    // bin has to catch everything below its own upper edge.
    int bk = static_cast<int>((ev.kin.k - k_edge0) / kbw);
    if (bk < 0) bk = 0;
    if (bk < nkb) hk[static_cast<std::size_t>(bk)] += 1.0;
    const int bc = static_cast<int>((ev.kin.cos_theta_k + 1.0) * 0.5 * ncb);
    if (bc >= 0 && bc < ncb) hc[static_cast<std::size_t>(bc)] += 1.0;
  });
  REQUIRE(n == n_ev);

  std::vector<double> ek(nkb, 0.0), ec(ncb, 0.0);
  double tot = 0.0;
  for (std::size_t ik = 0; ik < m.nk(); ++ik) {
    const double k2 = m.k()[ik] * m.k()[ik];
    const std::size_t bk = ik / kper;
    for (std::size_t ic = 0; ic < m.nc(); ++ic) {
      const double w = n1[ik * m.nc() + ic] * k2;
      tot += w;
      ek[bk] += w;
      ec[ic / cper] += w;
    }
  }

  auto chi2 = [&](const std::vector<double>& h, const std::vector<double>& e,
                  int& ndf) {
    double c2 = 0.0;
    ndf = 0;
    for (std::size_t i = 0; i < h.size(); ++i) {
      const double mu = static_cast<double>(n_ev) * e[i] / tot;
      if (mu < 20.0) continue;
      const double d = h[i] - mu;
      c2 += d * d / mu;
      ++ndf;
    }
    return c2;
  };
  int ndf_k = 0, ndf_c = 0;
  const double c2k = chi2(hk, ek, ndf_k);
  const double c2c = chi2(hc, ec, ndf_c);
  MESSAGE("spectator marginals vs TaggedModel::n_of_kc over " << n_ev
          << " events: chi2/ndf(k) = " << c2k << "/" << ndf_k
          << " = " << c2k / ndf_k << ", chi2/ndf(cos theta_k) = " << c2c << "/"
          << ndf_c << " = " << c2c / ndf_c);
  REQUIRE(ndf_k > 10);
  REQUIRE(ndf_c > 10);
  CHECK(c2k / ndf_k < 2.0);
  CHECK(c2c / ndf_c < 2.0);
}

// ------------------------------------------------------------ tag fractions

TEST_CASE("pipeline: 6Li alpha Roman-Pot tag fractions reproduce the Python") {
  // The published numbers (evgen/scripts/money_tagged_azz.py, plans/09 B2) at
  // 10 x 99.5 GeV/u: 0.0247 at the Yellow Report high-acceptance optics of the
  // configuration, and 0.3061 at the lithium tagging optics as it was priced
  // before 2026-08-29 (the 18x275 pot levers everywhere).  The per-
  // configuration levers, the default since then, give 0.2545 for the same
  // sample -- both are checked, because the pair is the whole content of the
  // lever change.
  const std::uint64_t n_ev = 400000;
  const Pipeline& p = pipe(PipelineChannel::TaggedLi6Alpha, 1,
                           OpticsChoice::YellowReportHighAcceptance, n_ev, 0);
  const double p_u = p.beam_config().ion_momentum_per_nucleon;
  const Optics ha = optics_for(OpticsChoice::YellowReportHighAcceptance, "6Li", p_u);
  const Optics tg = optics_for(OpticsChoice::Tagging, "6Li", p_u);
  const Optics tl = optics_for(OpticsChoice::TaggingLegacyLevers, "6Li", p_u);
  const std::string pc = p.pot_config();

  std::uint64_t a = 0, b = 0, c = 0, n = 0;
  p.for_each([&](const Event& ev) {
    ++n;
    a += rp_tagged(ev, ha, pc);
    b += rp_tagged(ev, tl, pc);
    c += rp_tagged(ev, tg, pc);
  });
  const double f_ha = static_cast<double>(a) / n;
  const double f_tl = static_cast<double>(b) / n;
  const double f_tg = static_cast<double>(c) / n;
  MESSAGE("6Li alpha tag at 10 x 99.5 over " << n << " events: "
          << "YR high-acceptance " << f_ha << " (Python 0.0247), "
          << "tagging/18x275 levers " << f_tl << " (Python 0.3061), "
          << "tagging/per-config levers " << f_tg << " (Python 0.2545)");
  CHECK_CLOSE(f_ha, 0.0247, 0.10);
  CHECK_CLOSE(f_tl, 0.3061, 0.10);
  CHECK_CLOSE(f_tg, 0.2545, 0.10);
  // the tagging optics buys acceptance with luminosity: a factor 13 at this
  // configuration, which is the programme statement of Report 4
  CHECK(tg.lumi_fraction < 0.09);
  CHECK(ha.lumi_fraction == 1.0);
}

TEST_CASE("pipeline: the 7Li alpha tag is optics-blind") {
  // R = 0.856 puts the 7Li alpha in the MIDDLE of the Roman-Pot momentum
  // window, so it never has to clear the near-beam envelope and its
  // acceptance barely moves with the optics (plans/09 B3).  Python, same
  // fills and the same 2e5 events:
  //     YR   0.9655 / 0.9735 / 0.9788
  //     tag  0.9803 / 0.9929 / 0.9944
  const double py_ha[3] = {0.9655, 0.9735, 0.9788};
  const double py_tg[3] = {0.9803, 0.9929, 0.9944};
  for (int ci = 0; ci < 3; ++ci) {
    CAPTURE(ci);
    const Pipeline& p = pipe(PipelineChannel::TaggedLi7Alpha, ci,
                             OpticsChoice::YellowReportHighAcceptance, 100000, 1);
    const double p_u = p.beam_config().ion_momentum_per_nucleon;
    const Optics ha = optics_for(OpticsChoice::YellowReportHighAcceptance, "7Li", p_u);
    const Optics tg = optics_for(OpticsChoice::Tagging, "7Li", p_u);
    std::uint64_t a = 0, b = 0, n = 0;
    p.for_each([&](const Event& ev) {
      ++n;
      a += rp_tagged(ev, ha, p.pot_config());
      b += rp_tagged(ev, tg, p.pot_config());
    });
    const double f_ha = static_cast<double>(a) / n;
    const double f_tg = static_cast<double>(b) / n;
    MESSAGE("7Li alpha tag, config " << ci << " (" << p.beam_config().label()
            << "): YR " << f_ha << " (Python " << py_ha[ci] << "), tagging "
            << f_tg << " (Python " << py_tg[ci] << ")");
    CHECK_CLOSE(f_ha, py_ha[ci], 0.05);
    CHECK_CLOSE(f_tg, py_tg[ci], 0.05);
    CHECK(f_tg / f_ha < 1.05);   // optics-blind: at most +5 %
  }
}

// -------------------------------------------------------------- A_zz^tag

TEST_CASE("pipeline: A_zz^tag(k) closes on the acceptance-weighted truth") {
  // The C++ shape of money_tagged_azz.py's folded closure.  Tensor-thirds
  // fills, events routed through the far-forward windows, the Roman-Pot
  // (main window + near-beam tail) thirds estimator vs the spectator momentum
  // k, against
  //     A_zz^acc(bin) = sum_{k in bin, c} eps(k,c) [n_+1 + n_-1 - 2 n_0] k^2
  //                   / sum_{k in bin, c} eps(k,c) [n_+1 + n_-1 +   n_0] k^2,
  // i.e. `tagged.azz_tensor_curve_weighted` averaged over the bin the way the
  // marker is, which is the comparison the Python's own 8e6-event run makes.
  //
  // That the thirds estimator returns the WAVE-FUNCTION asymmetry exactly is
  // an identity of the fill pattern: for the three tensor-thirds fills,
  // p^+ + p^- - 2 p^0 = P_zz (1, -2, 1) and p^+ + p^- + p^0 = (1, 1, 1) over
  // m = +1, 0, -1, so the P_zz divides out and only the m-dependence of the
  // spectator density is left.  This is the tagged-mode counterpart of
  // evgen/tests/test_tagged.py::test_tagged_rate_asymmetry_matches_analytic,
  // which pins the same statement on the k-INTEGRATED rates.
  const std::uint64_t n_ev = 1200000;
  const int nb = 12;
  const double klo = 0.0, khi = 0.6, kw = (khi - klo) / nb;

  for (int which = 0; which < 2; ++which) {
    const OpticsChoice oc = which == 0
                                ? OpticsChoice::YellowReportHighAcceptance
                                : OpticsChoice::Tagging;
    const Pipeline& p = pipe(PipelineChannel::TaggedLi6Alpha, 1, oc, n_ev, 0);
    REQUIRE(p.plan().categories().size() == 3);

    std::vector<std::vector<double>> cnt(3, std::vector<double>(nb, 0.0));
    p.for_each([&](const Event& ev) {
      if (!rp_tagged(ev, p.optics(), p.pot_config())) return;
      const int b = static_cast<int>((ev.kin.k - klo) / kw);
      if (b < 0 || b >= nb) return;
      cnt[static_cast<std::size_t>(ev.spin.bunch)][static_cast<std::size_t>(b)] += 1.0;
    });

    const TaggedModel& m = *p.tagged_model();
    const std::vector<double> eps = acceptance_weights(
        m, p.beam_config().ion_momentum_per_nucleon, p.optics(),
        p.pot_config(), 64);
    const std::vector<double>& np1 = m.n_of_kc(1.0);
    const std::vector<double>& n00 = m.n_of_kc(0.0);
    const std::vector<double>& nm1 = m.n_of_kc(-1.0);
    std::vector<double> num(nb, 0.0), den(nb, 0.0);
    for (std::size_t ik = 0; ik < m.nk(); ++ik) {
      const int b = static_cast<int>((m.k()[ik] - klo) / kw);
      if (b < 0 || b >= nb) continue;
      const double k2 = m.k()[ik] * m.k()[ik];
      for (std::size_t ic = 0; ic < m.nc(); ++ic) {
        const std::size_t j = ik * m.nc() + ic;
        const double w = eps[j] * k2;
        num[static_cast<std::size_t>(b)] += w * (np1[j] + nm1[j] - 2.0 * n00[j]);
        den[static_cast<std::size_t>(b)] += w * (np1[j] + nm1[j] + n00[j]);
      }
    }

    double chi2 = 0.0, worst = 0.0;
    int ndf = 0;
    for (int b = 0; b < nb; ++b) {
      const std::size_t i = static_cast<std::size_t>(b);
      const double n = cnt[0][i] + cnt[1][i] + cnt[2][i];
      if (n < 400.0 || !(den[i] > 0.0)) continue;
      const double a = estimators::azz_thirds(cnt[0][i], cnt[1][i], cnt[2][i], kPZZ);
      const double truth = num[i] / den[i];
      const double err = err_azz(n, kPZZ);
      const double pull = (a - truth) / err;
      chi2 += pull * pull;
      worst = std::max(worst, std::fabs(pull));
      ++ndf;
      MESSAGE("  " << p.optics().name << "  k = " << (klo + (b + 0.5) * kw)
              << "  n = " << n << "  A_zz = " << a << " +- " << err
              << "  truth " << truth << "  pull " << pull);
    }
    MESSAGE(p.optics().name << ": chi2/ndf = " << chi2 << "/" << ndf << " = "
            << chi2 / ndf << ", worst |pull| = " << worst);
    REQUIRE(ndf >= 8);
    CHECK(chi2 / ndf < 3.0);
    CHECK(worst < 4.0);
  }
}

TEST_CASE("pipeline: the tagged RATE asymmetry is the diluted inclusive Azz") {
  // The direct port of test_tagged.py::test_tagged_rate_asymmetry_matches_
  // analytic: with an INCLUSIVE b1 in the struck-cluster kernel (the
  // `li6_sampler` fixture's `b1_func=toy_b1`), the per-pure-ion-state accepted
  // cross sections satisfy
  //     (s_+1 + s_-1 - 2 s_0)/(s_+1 + s_-1 + s_0)
  //         = tensor_dilution * <Azz>_{sigma-weighted over the DIS cells}.
  // The wave-function part integrates to zero over 4 pi, so only the DIS-side
  // b1 survives here; the k-RESOLVED asymmetry of the test above is the other
  // half of the same physics and is why the generator's default kernel has no
  // inclusive b1 (money_tagged_azz.py, 2026-08-29).
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::TaggedLi6Alpha;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = 1;
  cfg.struck.inclusive_b1 = true;
  cfg.grid.nx = 24;
  cfg.grid.nq2 = 18;
  cfg.grid.x_min = 1e-3;
  cfg.grid.x_max = 0.5;
  cfg.grid.q2_min = 1.5;
  cfg.grid.q2_max = 100.0;

  std::vector<SpinCategory> cats;
  const std::vector<double> pops[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  for (int i = 0; i < 3; ++i) {
    cats.push_back(SpinCategory("s" + std::to_string(i), 1.0, pops[i], 0, 0.0,
                                0.0, 0.0, 1.0 / 3.0));
  }
  const Pipeline p(cfg, RunPlan(cats, 0.0, 0.0, 1.0));
  const std::vector<double>& s = p.sigma_per_category_pb();
  const double measured = (s[0] + s[2] - 2.0 * s[1]) / (s[0] + s[2] + s[1]);

  // sigma-weighted <Azz> over the struck cluster's accepted cells
  const InclusiveSampler& inner = p.dis_sampler();
  double wsum = 0.0, asum = 0.0;
  for (std::size_t c = 0; c < inner.n_cells(); ++c) {
    const double x = inner.x_cells()[c], q2 = inner.q2_cells()[c];
    const double y = y_from_xq2(x, q2, inner.s());
    const SFTables& t = inner.tables()[c];
    const double w = inner.cell_xsec_pb()[c];
    wsum += w;
    asum += w * azz(t.b1, t.f1, t.f2, x, y);
  }
  const double expected = p.tagged_model()->tensor_dilution() * (asum / wsum);
  MESSAGE("tagged rate asymmetry " << measured << " vs diluted inclusive "
          << expected << " (tensor dilution "
          << p.tagged_model()->tensor_dilution() << ")");
  CHECK_CLOSE(measured, expected, 2e-2);
}

// ------------------------------------------------------------ 7Li polarimeter

TEST_CASE("pipeline: 7Li generated events give <P2(cos theta_k)> = -T/5") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::TaggedLi7Alpha;
  cfg.isotope = "7Li";
  cfg.beam_config = 1;
  cfg.n_events = 300000;
  cfg.grid.nx = 40;
  cfg.grid.nq2 = 28;
  for (double t : {-0.6, -0.3, 0.0, 0.3, 0.6}) {
    CAPTURE(t);
    const std::array<double, 4> pops = spin32_populations(0.0, t, 0.0);
    const SpinCategory cat("t", 1.5,
                           {pops[0], pops[1], pops[2], pops[3]}, 0, 0.0, 0.0,
                           0.0, 1.0);
    const Pipeline p(cfg, RunPlan({cat}, 0.0, 0.0, t));
    double sum = 0.0;
    std::uint64_t n = 0;
    p.for_each([&](const Event& ev) {
      ++n;
      sum += p2(ev.kin.cos_theta_k);
    });
    const double got = sum / n;
    // sigma of the mean: P2 on the P-wave density has variance ~ 0.2
    const double err = 0.45 / std::sqrt(static_cast<double>(n));
    MESSAGE("7Li T = " << t << ": <P2> = " << got << " +- " << err
            << ", -T/5 = " << (-t / 5.0));
    CHECK(std::fabs(got + t / 5.0) < 4.0 * err);
  }
}

// --------------------------------------------------------------- coherent

TEST_CASE("pipeline: the coherent channel reproduces <|t|> = 1/B and the tag") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::CoherentLi6;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = 300000;
  cfg.grid.nx = 40;
  cfg.grid.nq2 = 28;
  // The lithium TAGGING optics: at the Yellow Report envelope the near-beam
  // cut on the intact 6Li is 1.07 GeV and exp(-B pT_cut^2) is 1e-25, so the
  // coherent recoil is only measurable at the de-squeezed working point.
  cfg.optics_choice = OpticsChoice::Tagging;
  const Pipeline p(cfg, tensor_thirds_plan(kPZ, kPZZ));

  const double p_u = p.beam_config().ion_momentum_per_nucleon;
  const CoherentScenario& sc = cfg.coherent;
  double sum_t = 0.0, t_max_seen = 0.0;
  std::uint64_t n = 0, acc = 0;
  p.for_each([&](const Event& ev) {
    ++n;
    sum_t += ev.kin.t;
    t_max_seen = std::max(t_max_seen, ev.kin.t);
    acc += rp_tagged(ev, p.optics(), p.pot_config());
    CHECK(ev.channel == Channel::CoherentLi6);
    CHECK(ev.find(Role::IntactRecoil) != nullptr);
  });

  // <|t|> of the exponential truncated at t_max
  const double b = sc.slope_b, tm = cfg.coherent_t_max;
  const double e_tm = std::exp(-b * tm);
  const double mean_t = 1.0 / b - tm * e_tm / (1.0 - e_tm);
  const double err_t = (1.0 / b) / std::sqrt(static_cast<double>(n));
  const double got_t = sum_t / n;

  // The tag acceptance.  `CoherentScenario::tag_acceptance_angular` is the
  // ISOTROPIC closed form exp(-B pT_cut^2), i.e. the recoil clearing a CIRCLE
  // of radius n_sigma sigma_h.  The pots are planar and `Optics::clears` is a
  // RECTANGLE -- |theta_x| > env_x OR |theta_y| > env_y -- whose inside
  // strictly contains the inscribed circle, so the true acceptance is lower.
  // The exact rectangle acceptance is a one-dimensional quadrature: at fixed
  // recoil azimuth the event is inside the rectangle iff
  //     |t| <= [pz tan(min(env_x/|cos phi|, env_y/|sin phi|))]^2,
  // and the |t| integral of B exp(-B|t|) truncated at t_max is closed form.
  const Optics& o = p.optics();
  const double pz = 6.0 * p_u;
  const int nphi = 2048;
  double inside = 0.0;
  for (int i = 0; i < nphi; ++i) {
    const double phi = (i + 0.5) * 2.0 * kPi / nphi;
    const double ca = std::fabs(std::cos(phi)), sa = std::fabs(std::sin(phi));
    const double tx = ca > 1e-12 ? std::tan(o.envelope_x() / ca) : 1e300;
    const double ty = sa > 1e-12 ? std::tan(o.envelope_y() / sa) : 1e300;
    const double pt_thresh = pz * std::min(tx, ty);
    const double t_thresh = std::min(pt_thresh * pt_thresh, tm);
    inside += (1.0 - std::exp(-b * t_thresh)) / (1.0 - e_tm);
  }
  const double want_acc = 1.0 - inside / nphi;
  const double got_acc = static_cast<double>(acc) / n;
  const double err_acc =
      std::sqrt(want_acc * (1.0 - want_acc) / static_cast<double>(n));
  const double circle =
      sc.tag_acceptance_angular(o.sigma_theta, p_u, 6, o.n_sigma);
  MESSAGE("coherent: sigma = " << p.sigma_pb() << " pb, <|t|> = " << got_t
          << " +- " << err_t << " (truncated exponential " << mean_t
          << ", 1/B = " << 1.0 / b << "); tag at " << o.name << " = "
          << got_acc << " +- " << err_acc << " (rectangle " << want_acc
          << ", inscribed circle exp(-B pT_cut^2) = " << circle
          << ", pT_cut = " << o.pt_cut_for(pz) << " GeV)");
  CHECK(t_max_seen <= tm);
  CHECK(std::fabs(got_t - mean_t) < 4.0 * err_t);
  CHECK(std::fabs(got_acc - want_acc) < 4.0 * err_acc);
  CHECK(got_acc < circle);           // the rectangle is the tighter cut
  CHECK(p.sigma_pb() > 0.0);
  // ... and at the Yellow Report envelope nothing at all is tagged
  const Optics yr = optics_for(OpticsChoice::YellowReportHighAcceptance, "6Li", p_u);
  std::uint64_t yr_acc = 0;
  p.for_each([&](const Event& ev) { yr_acc += rp_tagged(ev, yr, p.pot_config()); });
  CHECK(yr_acc == 0);
  MESSAGE("coherent tag at " << yr.name << ": " << yr_acc << " / " << n
          << " (analytic exp(-B pT_cut^2) = "
          << sc.tag_acceptance_angular(yr.sigma_theta, p_u, 6, yr.n_sigma) << ")");
}

// ------------------------------------------------- determinism and identity

TEST_CASE("pipeline: the inclusive path is bit-identical to InclusiveGenerator::run_n") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::Inclusive;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = 20000;
  cfg.grid.nx = 40;
  cfg.grid.nq2 = 28;
  const RunPlan plan = tensor_thirds_plan(kPZ, kPZZ);
  const Pipeline p(cfg, plan);

  // The generator on the pipeline's own sampler; the aliasing deleter keeps
  // ownership with the pipeline.
  const std::shared_ptr<const InclusiveSampler> shared(
      &p.dis_sampler(), [](const InclusiveSampler*) {});
  const InclusiveGenerator gen(shared);

  std::vector<Event> a, b;
  a.reserve(20000);
  b.reserve(20000);
  p.for_each([&](const Event& ev) { a.push_back(ev); });
  gen.run_n(plan, cfg.n_events, cfg.seed, [&](const Event& ev) { b.push_back(ev); },
            cfg.run);
  REQUIRE(a.size() == b.size());
  REQUIRE(a.size() == 20000);
  for (std::size_t i = 0; i < a.size(); ++i) {
    REQUIRE(a[i].particles.size() == b[i].particles.size());
    CHECK(a[i].number == b[i].number);
    CHECK(a[i].spin.category == b[i].spin.category);
    CHECK(a[i].kin.x == b[i].kin.x);      // bit-identical, not "close"
    CHECK(a[i].kin.q2 == b[i].kin.q2);
    CHECK(a[i].kin.phi == b[i].kin.phi);
    CHECK(a[i].spin.m_ion == b[i].spin.m_ion);
    for (std::size_t j = 0; j < a[i].particles.size(); ++j) {
      CHECK(a[i].particles[j].pdg == b[i].particles[j].pdg);
      CHECK(a[i].particles[j].p.e == b[i].particles[j].p.e);
      CHECK(a[i].particles[j].p.px == b[i].particles[j].p.px);
      CHECK(a[i].particles[j].p.py == b[i].particles[j].p.py);
      CHECK(a[i].particles[j].p.pz == b[i].particles[j].p.pz);
    }
  }
}

TEST_CASE("pipeline: threaded generation is bit-identical to single-threaded") {
  const PipelineChannel chans[3] = {PipelineChannel::Inclusive,
                                    PipelineChannel::TaggedLi6Alpha,
                                    PipelineChannel::CoherentLi6};
  for (PipelineChannel ch : chans) {
    const std::string chname = pipeline_channel_name(ch);
    CAPTURE(chname);
    const Pipeline& p = pipe(ch, 1, OpticsChoice::YellowReportHighAcceptance,
                             30000, 0);
    std::vector<Event> ref;
    ref.reserve(30000);
    p.for_each([&](const Event& ev) { ref.push_back(ev); });
    for (unsigned nt : {2u, 3u, 8u}) {
      std::vector<Event> got;
      got.reserve(30000);
      p.for_each([&](const Event& ev) { got.push_back(ev); }, nt, 777);
      REQUIRE(got.size() == ref.size());
      bool same = true;
      for (std::size_t i = 0; i < ref.size() && same; ++i) {
        if (got[i].number != ref[i].number ||
            got[i].particles.size() != ref[i].particles.size() ||
            got[i].kin.x != ref[i].kin.x || got[i].kin.phi != ref[i].kin.phi ||
            got[i].kin.k != ref[i].kin.k || got[i].kin.t != ref[i].kin.t ||
            got[i].weight != ref[i].weight) {
          same = false;
          break;
        }
        for (std::size_t j = 0; j < ref[i].particles.size(); ++j) {
          const Vec4& u = got[i].particles[j].p;
          const Vec4& v = ref[i].particles[j].p;
          if (u.e != v.e || u.px != v.px || u.py != v.py || u.pz != v.pz) {
            same = false;
            break;
          }
        }
      }
      CHECK(same);
    }
    // ... and so is the pull API, and a direct index lookup
    Pipeline& mut = const_cast<Pipeline&>(p);
    mut.rewind();
    Event ev;
    std::size_t i = 0;
    while (mut.next(ev)) {
      REQUIRE(i < ref.size());
      CHECK(ev.number == ref[i].number);
      CHECK(ev.kin.x == ref[i].kin.x);
      ++i;
    }
    CHECK(i == ref.size());
    const Event mid = p.event(12345);
    CHECK(mid.kin.x == ref[12345].kin.x);
    CHECK(mid.number == ref[12345].number);
  }
}

TEST_CASE("pipeline: cross sections are share-invariant and counts follow shares") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::TaggedLi6Alpha;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.lumi_pb = 20.0;
  cfg.poisson = false;
  cfg.grid.nx = 40;
  cfg.grid.nq2 = 28;
  const RunPlan plan = tensor_thirds_plan(kPZ, kPZZ);
  const Pipeline p1(cfg, plan);
  // the same plan with the m0-enriched share boosted by 20 %
  const Pipeline p2(cfg, with_offset(plan, "azz0", 0.2));
  REQUIRE(p1.sigma_per_category_pb().size() == 3);
  for (std::size_t k = 0; k < 3; ++k) {
    CHECK_CLOSE(p2.sigma_per_category_pb()[k], p1.sigma_per_category_pb()[k],
                1e-15);
  }
  // ... and a Scenario run_share moves counts, never pb
  PipelineConfig cfg3 = cfg;
  cfg3.scenario.run_share = 0.25;
  const Pipeline p3(cfg3, plan);
  for (std::size_t k = 0; k < 3; ++k) {
    CHECK_CLOSE(p3.sigma_per_category_pb()[k], p1.sigma_per_category_pb()[k],
                1e-15);
  }
  // counts are lumi_fraction * lumi * sigma
  for (std::size_t k = 0; k < 3; ++k) {
    const double mu = p1.lumi_per_category_pb()[k] * p1.sigma_per_category_pb()[k];
    CHECK(p1.counts()[k] == static_cast<std::uint64_t>(std::llround(mu)));
  }
  CHECK(p2.counts()[2] > p1.counts()[2]);
  MESSAGE("tagged 6Li sigma per category [pb]: " << p1.sigma_per_category_pb()[0]
          << ", " << p1.sigma_per_category_pb()[1] << ", "
          << p1.sigma_per_category_pb()[2] << "; mixture " << p1.sigma_pb()
          << " pb, " << p1.size() << " events at 20 pb^-1");
}

// -------------------------------------------------------------- HepMC3 I/O

#ifdef LIPOLGEN_TEST_HEPMC
namespace {

double pdg_charge(int pdg) {
  if (pdg == 11) return -1.0;
  if (pdg == -11) return 1.0;
  if (pdg == 22 || pdg == 2112) return 0.0;
  if (pdg == 2212) return 1.0;
  if (pdg > 1000000000) return static_cast<double>((pdg / 10000) % 1000);
  return 0.0;
}

}  // namespace

TEST_CASE("pipeline: 100 tagged events survive a HepMC3 round trip") {
  const Pipeline& p = pipe(PipelineChannel::TaggedLi6Alpha, 1,
                           OpticsChoice::YellowReportHighAcceptance, 100, 0);
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "lipolgen_pipeline_tagged.hepmc";
  {
    HepMC3Writer w(path.string());
    p.for_each([&](const Event& ev) { w.write(ev); });
    w.close();
  }

  std::vector<Event> ref;
  p.for_each([&](const Event& ev) { ref.push_back(ev); });
  REQUIRE(ref.size() == 100);

  HepMC3::ReaderAscii reader(path.string());
  std::size_t n = 0;
  double worst_p = 0.0, worst_q = 0.0;
  for (;; ++n) {
    HepMC3::GenEvent ge(HepMC3::Units::GEV, HepMC3::Units::MM);
    if (!reader.read_event(ge) || reader.failed()) break;
    REQUIRE(n < ref.size());
    CHECK(ge.particles().size() == ref[n].particles.size());

    // conservation checked on what came BACK off disk, from status codes and
    // PDG ids alone -- nothing of the in-memory record is consulted
    // The T0 hadronic system is written as DOCUMENTATION status 3 by the
    // writer's convention (`hepmc_status`, docs/HEPMC3_CONVENTION.md) so it
    // can never be double-counted against real T2 hadrons, so the outgoing
    // side is "status 1, plus the status-3 pdg-92 pseudo-particle".  With the
    // PYTHIA tier attached the pdg-92 X carries no status-1 twin and the same
    // sum still closes.
    HepMC3::FourVector in(0, 0, 0, 0), out(0, 0, 0, 0);
    double q_in = 0.0, q_out = 0.0;
    for (const auto& gp : ge.particles()) {
      const HepMC3::FourVector& v = gp->momentum();
      const bool final_state =
          gp->status() == 1 || (gp->status() == 3 && gp->pid() == 92);
      if (gp->status() == 4) {
        in.set(in.px() + v.px(), in.py() + v.py(), in.pz() + v.pz(),
               in.e() + v.e());
        q_in += pdg_charge(gp->pid());
      } else if (final_state) {
        out.set(out.px() + v.px(), out.py() + v.py(), out.pz() + v.pz(),
                out.e() + v.e());
        // pdg 92 is a container code and carries no charge of its own; the
        // charge it stands for is the struck cluster's, which is what closes
        // the nuclear balance.  Everything else is read from its PDG id.
        q_out += gp->pid() == 92 ? ref[n].find(Role::HadronicX)->charge
                                 : pdg_charge(gp->pid());
      }
    }
    worst_p = std::max(worst_p,
                       std::max(std::fabs(in.e() - out.e()),
                                std::max(std::fabs(in.px() - out.px()),
                                         std::max(std::fabs(in.py() - out.py()),
                                                  std::fabs(in.pz() - out.pz())))) /
                           std::fabs(in.e()));
    worst_q = std::max(worst_q, std::fabs(q_in - q_out));

    // the momenta themselves survived
    for (std::size_t i = 0; i < ge.particles().size(); ++i) {
      CHECK_CLOSE(ge.particles()[i]->momentum().e(), ref[n].particles[i].p.e,
                  1e-9);
      CHECK(ge.particles()[i]->pid() == ref[n].particles[i].pdg);
    }
  }
  reader.close();
  std::filesystem::remove(path);
  CHECK(n == 100);
  MESSAGE("HepMC3 round trip of 100 tagged events: worst relative |dP| = "
          << worst_p << ", worst |dQ| = " << worst_q);
  CHECK(worst_p < 1e-9);
  CHECK(worst_q < 1e-12);
}
#endif  // LIPOLGEN_TEST_HEPMC

// ------------------------------------------------------------- throughput

TEST_CASE("pipeline: single-core throughput per channel") {
  struct Row { PipelineChannel ch; int plan; };
  const Row rows[] = {{PipelineChannel::Inclusive, 0},
                      {PipelineChannel::TaggedLi6Alpha, 0},
                      {PipelineChannel::TaggedLi7Alpha, 1},
                      {PipelineChannel::CoherentLi6, 0}};
  for (const Row& row : rows) {
    const Pipeline& p = pipe(row.ch, 1, OpticsChoice::YellowReportHighAcceptance,
                             200000, row.plan);
    double acc = 0.0;
    const clock_t t0 = std::clock();
    p.for_each([&](const Event& ev) { acc += ev.kin.x; });
    const double secs = static_cast<double>(std::clock() - t0) / CLOCKS_PER_SEC;
    const double rate = static_cast<double>(p.size()) / secs;
    CHECK(acc > 0.0);
    MESSAGE(std::string(pipeline_channel_name(row.ch)) << ": " << rate / 1e3
            << " kevents/s single core (" << 1e9 * secs / p.size()
            << " ns/event)");
    CHECK(rate > 1e5);   // the plans/08 P8 target for the T0 tier
  }
}


// ------------------------------------------------- C1: the pomeron and M_X

// C1.  `PipelineConfig` had no x_P knob and `CoherentSampler::set_x_pom` was
// never called, so x_P was 0 in every event: the recoil was the beam ion
// itself and X = k + P_ion - k' - P_recoil was just the virtual photon,
// M_X^2 = -Q^2 < 0 in 100 % of events, silently clipped to a massless X.
TEST_CASE("pipeline: the coherent channel draws x_P and gives a timelike X") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::CoherentLi6;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = 200000;
  cfg.grid.nx = 40;
  cfg.grid.nq2 = 28;
  cfg.optics_choice = OpticsChoice::Tagging;
  const Pipeline p(cfg, tensor_thirds_plan(kPZ, kPZZ));
  const CoherentXpomModel& xm = cfg.coherent_xpom;
  const double m_x_min2 = xm.m_x_min * xm.m_x_min;

  double xp_lo = 1e300, xp_hi = 0.0, mx2_lo = 1e300, beta_hi = 0.0;
  double r_lo = 1e300, r_hi = 0.0, worst_mx = 0.0;
  std::uint64_t n = 0, above_max = 0;
  p.for_each([&](const Event& ev) {
    ++n;
    const Particle* x = ev.find(Role::HadronicX);
    const Particle* rec = ev.find(Role::IntactRecoil);
    REQUIRE(x != nullptr);
    REQUIRE(rec != nullptr);

    // --- X is timelike and carries EXACTLY the drawn diffractive mass
    CHECK(x->p.m2() > 0.0);
    CHECK(ev.kin.m_x2 >= m_x_min2 * (1.0 - 1e-12));
    worst_mx = std::max(worst_mx,
                        std::fabs(x->p.m2() / ev.kin.m_x2 - 1.0));
    mx2_lo = std::min(mx2_lo, ev.kin.m_x2);
    // ... and X is NEUTRAL and carries no charge of the intact nucleus
    CHECK(x->charge == 0.0);

    // --- x_P is inside the model's own window at this event's kinematics
    const double lo = xm.x_pom_min(ev.kin.q2, ev.kin.w2);
    CHECK(ev.kin.x_pom >= lo * (1.0 - 1e-12));
    if (ev.kin.x_pom > xm.x_pom_max * (1.0 + 1e-9)) ++above_max;
    xp_lo = std::min(xp_lo, ev.kin.x_pom);
    xp_hi = std::max(xp_hi, ev.kin.x_pom);

    // --- beta = x/x_P is a momentum fraction.  With the per-nucleon
    // W^2 = Q^2 (1-x)/x + M_N^2 the identity is exact as written below and
    // collapses to the textbook Q^2/(M_X^2 + Q^2) when the nucleon mass is
    // dropped (they differ by the factor 1 + x M_N^2/Q^2, ~1e-3 here).
    CHECK(ev.kin.beta_pom > 0.0);
    CHECK(ev.kin.beta_pom <= 1.0);
    CHECK_CLOSE(ev.kin.beta_pom,
                ev.kin.x * (ev.kin.w2 + ev.kin.q2)
                    / (ev.kin.m_x2 + ev.kin.q2), 1e-12);
    CHECK_CLOSE(ev.kin.beta_pom, ev.kin.q2 / (ev.kin.m_x2 + ev.kin.q2), 0.05);
    beta_hi = std::max(beta_hi, ev.kin.beta_pom);

    // --- the recoil is ON SHELL at the 6Li mass and still near-beam
    CHECK_CLOSE(rec->p.m2(), rec->mass * rec->mass, 1e-9);
    const double pt = rec->p.pt();
    const double rr = std::sqrt(pt * pt + rec->p.pz * rec->p.pz)
                      / (6.0 * p.beam_config().ion_momentum_per_nucleon);
    r_lo = std::min(r_lo, rr);
    r_hi = std::max(r_hi, rr);
    CHECK(std::fabs(rr - 1.0) < NEAR_BEAM_BAND);
    // pT^2 IS the sampled |t|; the exact Mandelstam t adds |t_min|
    CHECK_CLOSE(pt * pt, ev.kin.t, 1e-9);
  });
  MESSAGE("coherent x_P in [" << xp_lo << ", " << xp_hi << "] (x_P,max = "
          << xm.x_pom_max << ", " << above_max << "/" << n
          << " above it from the intra-cell redraw), min M_X^2 = " << mx2_lo
          << " (M_X,min^2 = " << m_x_min2 << "), max beta = " << beta_hi
          << ", recoil R in [" << r_lo << ", " << r_hi
          << "], worst |M_X^2(record)/M_X^2(model) - 1| = " << worst_mx);
  CHECK(n == 200000);
  // the recoil construction is EXACT, not a small-x_P expansion
  CHECK(worst_mx < 1e-9);
  // the window is really scanned, not pinned
  CHECK(xp_lo < 1e-2);
  CHECK(xp_hi > 5e-2);
  // only the intra-cell redraw can leave the window, and only just
  CHECK(above_max * 200 < n);
  CHECK(xp_hi < 1.5 * xm.x_pom_max);
}

// P4.  c2 = -(P_zz/2) eps_B0 B |t| + amp P_zz is linear and UNBOUNDED in |t|,
// so the azimuthal weight 1 + c2 cos 2(phi_t - phi_S) goes negative at large
// |t|.  The default t_max was 0.5, past the |t| = 0.245 zero at P_zz = -2 and
// past the |t| = 0.495 zero at P_zz = +1, and past the |t| <= 0.30 range of
// the Mantysaari input the deformation term is scaled from.
TEST_CASE("coherent: the azimuthal weight is positive over the whole t range") {
  const CoherentScenario sc;
  CHECK(COHERENT_T_MAX_DEFAULT == 0.2);

  // the analytic zeros the new default avoids
  CHECK_CLOSE_AT(sc.positivity_margin(0.495, 1.0), 0.0, 0.0, 1e-3);
  CHECK_CLOSE_AT(sc.positivity_margin(0.245, -2.0), 0.0, 0.0, 1e-3);
  CHECK(sc.positivity_margin(0.5, 1.0) < 0.0);
  CHECK(sc.positivity_margin(0.5, -2.0) < 0.0);
  for (double pzz : {1.0, 0.6, 0.0, -0.6, -2.0}) {
    CHECK(sc.positivity_margin(COHERENT_T_MAX_DEFAULT, pzz) > 0.0);
  }

  // the SAMPLER refuses the range at setup, the way the inclusive kernel does
  const double p_u = default_configs("6Li")[1].ion_momentum_per_nucleon;
  CoherentSampler s(sc, p_u, -2.0, 0.0);
  CHECK(s.t_max() == COHERENT_T_MAX_DEFAULT);
  CHECK_THROWS(s.set_t_max(0.5));
  CHECK(s.t_max() == COHERENT_T_MAX_DEFAULT);   // and leaves it alone
  CHECK_NOTHROW(s.set_t_max(0.24));
  // ... and so does the Pipeline, at every category's own P_zz
  PipelineConfig bad;
  bad.channel = PipelineChannel::CoherentLi6;
  bad.isotope = "6Li";
  bad.beam_config = 1;
  bad.n_events = 100;
  bad.coherent_t_max = 0.5;
  CHECK_THROWS(Pipeline(bad, tensor_thirds_plan(kPZ, kPZZ)));

  // no generated event carries a negative weight
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::CoherentLi6;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = 50000;
  cfg.grid.nx = 30;
  cfg.grid.nq2 = 20;
  // the m0-enriched fill sits at P_zz = -2 kPZZ, the worst case
  const Pipeline p(cfg, tensor_thirds_plan(kPZ, kPZZ));
  double w_lo = 1e300, w_hi = 0.0, t_hi = 0.0;
  p.for_each([&](const Event& ev) {
    w_lo = std::min(w_lo, ev.weight);
    w_hi = std::max(w_hi, ev.weight);
    t_hi = std::max(t_hi, ev.kin.t);
  });
  MESSAGE("coherent event weights in [" << w_lo << ", " << w_hi
          << "] over |t| <= " << t_hi);
  CHECK(w_lo > 0.0);
  CHECK(t_hi <= COHERENT_T_MAX_DEFAULT);
}

// C4.  `Pipeline::event` called `cfg_.hadronizer` on EVERY channel, coherent
// included.  A coherent event carries neither a struck nucleon nor a struck
// cluster, so PythiaBridge falls through to its inclusive branch and invents
// a nucleon that is not in the record's balance.
TEST_CASE("pipeline: the hadronizer hook is refused on the coherent channel") {
  struct Row { PipelineChannel ch; const char* iso; int plan; bool coherent; };
  const Row rows[] = {
      {PipelineChannel::Inclusive, "6Li", 0, false},
      {PipelineChannel::TaggedLi6Alpha, "6Li", 0, false},
      {PipelineChannel::TaggedDeuteronP, "d", 0, false},
      {PipelineChannel::CoherentLi6, "6Li", 0, true},
  };
  for (const Row& row : rows) {
    const std::string what = pipeline_channel_name(row.ch);
    CAPTURE(what);
    std::uint64_t calls = 0;
    PipelineConfig cfg;
    cfg.channel = row.ch;
    cfg.isotope = row.iso;
    cfg.beam_config = 1;
    cfg.n_events = 500;
    cfg.grid.nx = 24;
    cfg.grid.nq2 = 16;
    cfg.hadronizer = [&calls](Event&, Rng&) { ++calls; };

    if (row.coherent) {
      // refused at configuration time, before a single event is built
      CHECK_THROWS_AS(Pipeline(cfg, tensor_thirds_plan(kPZ, kPZZ)),
                      std::runtime_error);
      CHECK(calls == 0);
      // ... and the escape hatch is explicit and documented
      cfg.hadronize_coherent = true;
      const Pipeline opt_in(cfg, tensor_thirds_plan(kPZ, kPZZ));
      opt_in.for_each([](const Event&) {});
      CHECK(calls == opt_in.size());
      // no hook at all is always fine
      PipelineConfig plain = cfg;
      plain.hadronizer = nullptr;
      plain.hadronize_coherent = false;
      CHECK_NOTHROW(Pipeline(plain, tensor_thirds_plan(kPZ, kPZZ)));
    } else {
      const Pipeline p(cfg, tensor_thirds_plan(kPZ, kPZZ));
      p.for_each([](const Event&) {});
      CHECK(calls == p.size());
      CHECK(calls > 0);
    }
  }
}

// C6.  `Optics::lumi_fraction` -- the share of the machine luminosity a
// far-forward working point actually delivers -- was computed, carried on
// every `Optics`, and never applied to anything.  The lithium TAGGING optics
// buy their acceptance by de-squeezing beta*_x, and that costs luminosity.
TEST_CASE("pipeline: the optics luminosity fraction reaches the event count") {
  auto counts = [](OpticsChoice oc, bool apply) {
    PipelineConfig cfg;
    cfg.channel = PipelineChannel::TaggedLi6Alpha;
    cfg.isotope = "6Li";
    cfg.beam_config = 0;              // 5x41, where the de-squeeze costs most
    cfg.lumi_pb = 20.0;
    cfg.poisson = false;              // the expectation, not a draw
    cfg.optics_choice = oc;
    cfg.apply_optics_lumi_fraction = apply;
    cfg.grid.nx = 24;
    cfg.grid.nq2 = 16;
    return Pipeline(cfg, tensor_thirds_plan(kPZ, kPZZ));
  };
  const Pipeline yr = counts(OpticsChoice::YellowReportHighAcceptance, true);
  const Pipeline tag = counts(OpticsChoice::Tagging, true);
  const Pipeline tag_off = counts(OpticsChoice::Tagging, false);

  const double f = tag.optics().lumi_fraction;
  MESSAGE("6Li 5x41: YR lumi_fraction = " << yr.optics().lumi_fraction
          << ", tagging lumi_fraction = " << f << "; events "
          << yr.size() << " (YR) vs " << tag.size() << " (tagging) vs "
          << tag_off.size() << " (tagging, knob off)");

  CHECK(yr.optics().lumi_fraction == 1.0);
  CHECK(f < 0.2);
  CHECK(yr.optics_lumi_factor() == 1.0);
  CHECK(tag.optics_lumi_factor() == f);
  CHECK(tag_off.optics_lumi_factor() == 1.0);

  // the tagging point delivers lumi_fraction x the events at the same lumi_pb
  REQUIRE(yr.size() > 0);
  CHECK_CLOSE(static_cast<double>(tag.size())
                  / static_cast<double>(yr.size()), f, 2e-3);
  // ... and the knob is what does it
  CHECK_CLOSE(static_cast<double>(tag_off.size())
                  / static_cast<double>(yr.size()), 1.0, 2e-3);
  // THE SHARE RULE: cross sections never see it (bookkeeping.hpp)
  REQUIRE(tag.sigma_per_category_pb().size()
          == tag_off.sigma_per_category_pb().size());
  for (std::size_t k = 0; k < tag.sigma_per_category_pb().size(); ++k) {
    CHECK(tag.sigma_per_category_pb()[k] == tag_off.sigma_per_category_pb()[k]);
  }
  CHECK(tag.sigma_pb() == tag_off.sigma_pb());
  // it is the LUMINOSITY that carries it
  for (std::size_t k = 0; k < tag.lumi_per_category_pb().size(); ++k) {
    CHECK_CLOSE(tag.lumi_per_category_pb()[k],
                f * tag_off.lumi_per_category_pb()[k], 1e-12);
  }
}
