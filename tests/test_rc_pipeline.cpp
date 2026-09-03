// The PIPELINE half of the tensor-sector radiative corrections (rc.hpp):
// the wiring, and the one gate everything else rests on.
//
// T6 -- `--rc tensor-band` IS `--rc off`, bit for bit, in every four-vector,
// in `Event::weight`, in every kinematic and spin label AND in the random
// stream.  RC is a weight family that hangs off the FINISHED record; it must
// move nothing.  The structural guarantee is `RcModel::fill(Event&) const`
// taking no `Rng&` (rc.hpp); this file is the runtime one.
//
// The stream check uses a `HadronizerHook` and not a "counting wrapper":
// `Rng` (rng.hpp) exposes only uniform()/normal()/next_u64() and
// `Pipeline::event` constructs its own `Rng` internally, so the stream
// position is not observable from outside.  Installing
// `[&](Event&, Rng& r){ log.push_back(r.uniform()); }` on BOTH pipelines and
// comparing the two logs with `==` is.  Comparing T0 four-vectors alone does
// NOT prove it, because the RC fill runs AFTER every T0 draw and before the
// hook.
//
// Also here: the config-level guards (`validate()` accepts the coherent
// channel, refuses a bad knob), and the `Pipeline::rc_model()` accessor.

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "doctest.h"
#include "lipolgen/pipeline.hpp"
#include "lipolgen/rc.hpp"

using namespace lipolgen;

namespace {

/// The tensor-thirds plan `--plan tensor-thirds --pzz 0.6` uses, with an
/// UNPOLARISED beam (lam_e = 0), which is what the whole A_zz programme --
/// and every tau identity in rc.hpp -- assumes.
RunPlan thirds_plan() { return tensor_thirds_plan(0, 0.6); }

PipelineConfig base_cfg(PipelineChannel ch, std::uint64_t n) {
  PipelineConfig cfg;
  cfg.channel = ch;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = n;
  cfg.grid.nx = 40;
  cfg.grid.nq2 = 28;
  cfg.seed = 20260903u;
  cfg.run = 7;
  return cfg;
}

/// Bit-for-bit equality that also holds for the NaN sentinels the record
/// carries (`SpinLabels::m_struck` is NaN off the tagged channels), where
/// `==` is false against itself.
bool same_bits(double a, double b) {
  return (a == b) || (std::isnan(a) && std::isnan(b));
}

/// Every scalar of a finished record that is NOT `rc_weights`, flattened so a
/// mismatch names itself.  `==`, never `Approx`: this is a bit-for-bit gate.
void require_same_record(const Event& a, const Event& b) {
  REQUIRE(a.number == b.number);
  REQUIRE(a.spin.category == b.spin.category);
  REQUIRE(a.channel == b.channel);
  REQUIRE(a.weight == b.weight);
  REQUIRE(a.xsec_pb == b.xsec_pb);
  REQUIRE(a.xsec_err_pb == b.xsec_err_pb);
  REQUIRE(a.spin_weights == b.spin_weights);
  REQUIRE(a.spin.j == b.spin.j);
  REQUIRE(a.spin.m_ion == b.spin.m_ion);
  REQUIRE(same_bits(a.spin.m_struck, b.spin.m_struck));
  REQUIRE(a.spin.lam_e == b.spin.lam_e);
  REQUIRE(a.spin.pe == b.spin.pe);
  REQUIRE(a.spin.theta_s == b.spin.theta_s);
  REQUIRE(a.spin.phi_s == b.spin.phi_s);
  REQUIRE(a.kin.x == b.kin.x);
  REQUIRE(a.kin.q2 == b.kin.q2);
  REQUIRE(a.kin.y == b.kin.y);
  REQUIRE(a.kin.w2 == b.kin.w2);
  REQUIRE(a.kin.nu == b.kin.nu);
  REQUIRE(a.kin.phi == b.kin.phi);
  REQUIRE(a.kin.cell == b.kin.cell);
  REQUIRE(same_bits(a.kin.k, b.kin.k));
  REQUIRE(same_bits(a.kin.cos_theta_k, b.kin.cos_theta_k));
  REQUIRE(a.particles.size() == b.particles.size());
  for (std::size_t j = 0; j < a.particles.size(); ++j) {
    const Particle& pa = a.particles[j];
    const Particle& pb = b.particles[j];
    REQUIRE(pa.pdg == pb.pdg);
    REQUIRE(pa.status == pb.status);
    REQUIRE(pa.role == pb.role);
    REQUIRE(pa.mass == pb.mass);
    REQUIRE(pa.charge == pb.charge);
    REQUIRE(pa.mother1 == pb.mother1);
    REQUIRE(pa.mother2 == pb.mother2);
    REQUIRE(pa.p.e == pb.p.e);
    REQUIRE(pa.p.px == pb.p.px);
    REQUIRE(pa.p.py == pb.p.py);
    REQUIRE(pa.p.pz == pb.p.pz);
  }
}

}  // namespace

TEST_CASE("T6: --rc tensor-band is --rc off bit for bit, RNG included") {
  const std::uint64_t n = 5000;
  const RunPlan plan = thirds_plan();

  for (PipelineChannel ch : {PipelineChannel::Inclusive,
                             PipelineChannel::TaggedLi6Alpha}) {
    PipelineConfig cfg_off = base_cfg(ch, n);
    PipelineConfig cfg_on = cfg_off;
    cfg_on.rc = PipelineRc::TensorBand;

    // The stream probe: one uniform per finished event, drawn from the
    // event's OWN counter-based stream at the point the hadronizer would see
    // it -- i.e. AFTER the RC fill.  If RC consumed a single random number,
    // these two logs would diverge from the first event.
    std::vector<double> log_off, log_on;
    cfg_off.hadronizer = [&log_off](Event&, Rng& r) {
      log_off.push_back(r.uniform());
    };
    cfg_on.hadronizer = [&log_on](Event&, Rng& r) {
      log_on.push_back(r.uniform());
    };

    const Pipeline off(cfg_off, plan);
    const Pipeline on(cfg_on, plan);
    REQUIRE(off.rc_model() == nullptr);
    REQUIRE(on.rc_model() != nullptr);
    REQUIRE(on.rc_model()->mode() == RcMode::TensorBand);
    REQUIRE(off.size() == on.size());

    Event ea, eb;
    std::size_t n_nonunit = 0;
    for (std::uint64_t i = 0; i < off.size(); ++i) {
      off.event(i, ea);
      on.event(i, eb);
      require_same_record(ea, eb);
      // The one thing that DOES differ.
      REQUIRE(ea.rc_weights.empty());
      REQUIRE(eb.rc_weights.size() == kRcWeightCount);
      // w_lo + w_hi == 2 exactly (T4a), and the tail never subtracts events.
      REQUIRE(eb.rc_weights[0] + eb.rc_weights[1] == 2.0);
      REQUIRE(eb.rc_weights[2] >= 1.0);
      if (eb.rc_weights[1] != 1.0) ++n_nonunit;
    }
    REQUIRE(log_off.size() == n);
    REQUIRE(log_on == log_off);   // zero extra uniforms consumed
    // The gate would pass trivially if the band were identically 1.
    CHECK(n_nonunit > n / 2);
    MESSAGE("T6 " << std::string(pipeline_channel_name(ch)) << ": " << n_nonunit << " / "
            << n << " events carry a non-unit band");
  }
}

TEST_CASE("T6b: the coherent channel is ALLOWED, and prices nothing") {
  // Unlike FSI, `validate()` must NOT refuse a channel.  `RcModel::applies()`
  // decides, the weights are exactly 1.0 and the run prints the reason, so a
  // scan over channels does not have to special-case `--rc`
  // (design_C_tensor_rc.md sec. 1.5.3).
  PipelineConfig cfg = base_cfg(PipelineChannel::CoherentLi6, 400);
  cfg.rc = PipelineRc::TensorBand;
  REQUIRE_NOTHROW(cfg.validate());
  const RunPlan plan = thirds_plan();
  const Pipeline p(cfg, plan);
  REQUIRE(p.rc_model() != nullptr);
  CHECK_FALSE(p.rc_model()->applies());
  CHECK_FALSE(p.rc_model()->tail_applies());
  CHECK_FALSE(p.rc_model()->exclusion_reason().empty());
  Event ev;
  for (std::uint64_t i = 0; i < p.size(); ++i) {
    p.event(i, ev);
    REQUIRE(ev.rc_weights.size() == kRcWeightCount);
    REQUIRE(ev.rc_weights[0] == 1.0);
    REQUIRE(ev.rc_weights[1] == 1.0);
    REQUIRE(ev.rc_weights[2] == 1.0);
  }
}

TEST_CASE("T6c: PipelineConfig::validate refuses only a bad RC knob") {
  PipelineConfig cfg = base_cfg(PipelineChannel::Inclusive, 100);
  REQUIRE_NOTHROW(cfg.validate());   // rc = Off: nothing is checked
  cfg.rc_options.n_eta = 2;
  REQUIRE_NOTHROW(cfg.validate());   // ... still nothing

  cfg.rc = PipelineRc::TensorBand;
  CHECK_THROWS_AS(cfg.validate(), std::runtime_error);
  cfg.rc_options.n_eta = 128;
  REQUIRE_NOTHROW(cfg.validate());

  SUBCASE("delta") {
    cfg.rc_options.delta_low_x = -1.0;
    CHECK_THROWS_AS(cfg.validate(), std::runtime_error);
  }
  SUBCASE("anchors") {
    cfg.rc_options.x_high = cfg.rc_options.x_low;
    CHECK_THROWS_AS(cfg.validate(), std::runtime_error);
  }
  SUBCASE("scales") {
    cfg.rc_options.fq_scale = -0.5;
    CHECK_THROWS_AS(cfg.validate(), std::runtime_error);
  }
  SUBCASE("qe") {
    cfg.rc_options.qe_suppression = -1.0;
    CHECK_THROWS_AS(cfg.validate(), std::runtime_error);
  }
  SUBCASE("tail model") {
    cfg.rc_options.tail_model = RcTailModel::PolradFull;
    CHECK_THROWS_AS(cfg.validate(), std::runtime_error);
  }
  SUBCASE("a POLARISED beam is refused -- by RcModel, one frame later") {
    // design sec. 3.2 puts this loop in `validate()`; `PipelineConfig` has no
    // run plan (RunPlan is the Pipeline's second constructor argument), so
    // RcModel's own constructor is where it lives.  The run still throws.
    std::vector<SpinCategory> cats;
    cats.push_back(SpinCategory("plus", 1.0, {1.0, 0.0, 0.0}, +1, 0.7, 0.0,
                                0.0, 1.0));
    const RunPlan pol(cats, 0.0, 1.0, 1.0);
    REQUIRE_NOTHROW(cfg.validate());
    CHECK_THROWS_AS(Pipeline(cfg, pol), std::runtime_error);
  }
}

TEST_CASE("T18b: threaded generation reproduces the RC weights bit for bit") {
  PipelineConfig cfg = base_cfg(PipelineChannel::Inclusive, 2000);
  cfg.rc = PipelineRc::TensorBand;
  const RunPlan plan = thirds_plan();
  const Pipeline p(cfg, plan);
  std::vector<std::vector<double>> one(p.size()), many(p.size());
  p.for_each([&one](const Event& ev) { one[ev.number] = ev.rc_weights; });
  p.for_each([&many](const Event& ev) { many[ev.number] = ev.rc_weights; }, 4);
  REQUIRE(one.size() == many.size());
  for (std::size_t i = 0; i < one.size(); ++i) {
    REQUIRE(one[i].size() == kRcWeightCount);
    REQUIRE(one[i] == many[i]);
  }
}
