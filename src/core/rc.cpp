#include "lipolgen/rc.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace lipolgen {

double rc_delta(double x, double delta_high, double delta_low, double x_high,
                double x_low) {
  if (!(x_low > 0.0) || !(x_high > x_low)) {
    throw std::runtime_error(
        "rc_delta: the band anchors must satisfy 0 < x_low < x_high (got "
        "x_low = " + std::to_string(x_low) + ", x_high = " +
        std::to_string(x_high) + "); the interpolation divides by "
        "ln(x_high/x_low)");
  }
  // Both clamp branches return the anchor UNMODIFIED, which is what makes
  // `rc_delta(x_high) == delta_high` and `rc_delta(x_low) == delta_low` exact
  // (T3 compares them with ==, legitimately).
  if (x >= x_high) return delta_high;
  if (x <= x_low) return delta_low;
  const double f = std::log(x_high / x) / std::log(x_high / x_low);
  return delta_high + (delta_low - delta_high) * f;
}

const char* rc_mode_name(RcMode m) {
  switch (m) {
    case RcMode::Off:        return "off";
    case RcMode::TensorBand: return "tensor-band";
  }
  throw std::runtime_error("rc_mode_name: unhandled RcMode");
}

const char* pipeline_rc_name(PipelineRc r) { return rc_mode_name(r); }

const char* rc_weight_name(std::size_t i) {
  switch (i) {
    case 0: return "rc_tensor_lo";
    case 1: return "rc_tensor_hi";
    case 2: return "rc_tail";
    default: break;
  }
  throw std::runtime_error("rc_weight_name: index " + std::to_string(i) +
                           " is outside the kRcWeightCount = 3 block");
}

std::string rc_weight_name(std::size_t i, std::size_t slot) {
  // Slot 0 is the event's OWN spin state and keeps the bare name, so an
  // unweighted run's HepMC3 weight names are exactly
  // {"nominal", "rc_tensor_lo", "rc_tensor_hi", "rc_tail"}.  Slot 1 + k is
  // spin category k and carries the "_k" suffix of "spin_weight_k".
  const std::string base = rc_weight_name(i);
  return slot == 0 ? base : base + "_" + std::to_string(slot);
}

}  // namespace lipolgen
