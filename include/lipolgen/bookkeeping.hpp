#ifndef LIPOLGEN_BOOKKEEPING_HPP
#define LIPOLGEN_BOOKKEEPING_HPP

/// \file bookkeeping.hpp
/// Run-plan bookkeeping: spin categories, luminosity shares, polarimetry.
/// C++17 port of `evgen/polligen/bookkeeping.py`.
///
/// A run plan is a list of SpinCategory: each combines an ion fill state
/// (populations along an axis n(theta_S, phi_S)) with an electron helicity
/// sign and carries its absolute share of the integrated luminosity,
/// including relative-luminosity offsets.  Every generated event is labeled
/// by its category, and the analysis-side estimators (`estimators` in
/// sampler.hpp) see exactly the per-category luminosities the "experiment"
/// would measure.
///
/// THREE LUMINOSITY FACTORS MULTIPLY AND NONE IS THE OTHER
/// (`polli_fastsim.fom` module docstring, evgen/tests/test_run_share.py):
///
///   programme share            how a year divides between observables,
///                              isotopes and configurations -- carried as
///                              `Scenario::run_share` (sampler.hpp);
///   optics fraction            what a de-squeezed beta*_x costs at fixed
///                              wall time (1/7 - 1/13 for lithium tagging);
///   SPIN-STATE share           how ONE measurement's own luminosity divides
///                              between its fills -- `lumi_fraction` here.
///                              It sums to one within a run plan.
///
/// A share moves COUNTS, never CROSS SECTIONS: `lumi_shares()` multiplies a
/// total luminosity, and nothing in this header or in InclusiveSampler lets
/// a share reach a pb number (test_run_share.py
/// `test_the_sampler_cross_sections_are_share_invariant`).
///
/// This module also owns the two Step-5.A systematics knobs (plans/05 5.0):
///
///  * relative-luminosity offsets between spin states.  First-order biases of
///    the naive (equal-share) estimators:
///        tensor thirds:   bias(Azz)   = -(2/3) delta / P_zz
///                         (offset delta on the m0-enriched share)
///        helicity flips:  bias(A_par) = +delta / (2 P_e P_z)
///                         (offset delta on the lam_e=+1 share)
///    The lumi-corrected estimators remove them.
///  * polarimetry scale: measured = true * (1 + delta_p_over_p * gauss),
///    drawn once per plan from a fixed-seed counter-based stream
///    (delta P/P ~ 3 %, plans/04 #5).

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "lipolgen/constants.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/spin.hpp"

namespace lipolgen {

/// One (ion fill state, electron helicity) luminosity category.
struct SpinCategory {
  std::string name;
  double j = 1.0;
  /// Populations p_m ordered m = +J ... -J along the axis (library-wide).
  std::vector<double> populations;
  int lam_e = 0;             ///< +1/-1 electron helicity; 0 = unpolarized
  double pe = 0.0;           ///< electron polarization magnitude
  double theta_s = 0.0;      ///< quantization-axis polar angle [rad]
  double phi_s = 0.0;        ///< quantization-axis azimuth [rad]
  double lumi_fraction = 1.0;  ///< absolute share of the total luminosity

  SpinCategory() = default;
  SpinCategory(std::string name_, double j_, std::vector<double> populations_,
               int lam_e_ = 0, double pe_ = 0.0, double theta_s_ = 0.0,
               double phi_s_ = 0.0, double lumi_fraction_ = 1.0)
      : name(std::move(name_)), j(j_), populations(std::move(populations_)),
        lam_e(lam_e_), pe(pe_), theta_s(theta_s_), phi_s(phi_s_),
        lumi_fraction(lumi_fraction_) {}

  /// `spin.moments_along_axis(j, populations)` -- (vector, tensor[, octupole])
  /// along the fill's OWN axis.  Throws for j outside {1, 3/2}, exactly as
  /// the Python does (a spin-1/2 fill has no rank-2 moment at all; use
  /// `vector_moment()` there).
  AxisMoments moments() const { return moments_along_axis(j, populations); }

  /// <J_z>/J along the fill axis, defined for every j (including 1/2).
  double vector_moment() const;
};

/// Spin categories + true / measured polarization bookkeeping.
class RunPlan {
 public:
  RunPlan() = default;
  /// Throws std::runtime_error if the categories carry no luminosity.
  RunPlan(std::vector<SpinCategory> categories, double pe_true = 0.0,
          double pz_true = 0.0, double pzz_true = 0.0,
          double delta_p_over_p = 0.0,
          std::uint64_t polarimetry_seed = 20260713);

  const std::vector<SpinCategory>& categories() const { return categories_; }
  double pe_true() const { return pe_true_; }
  double pz_true() const { return pz_true_; }
  double pzz_true() const { return pzz_true_; }
  double delta_p_over_p() const { return delta_p_over_p_; }
  std::uint64_t polarimetry_seed() const { return polarimetry_seed_; }

  /// The polarimeter's numbers -- what the analysis is allowed to divide by.
  /// Equal to the true values when delta_p_over_p == 0.  A true value of
  /// exactly zero stays exactly zero under the smear (a smeared 1e-17
  /// divisor would be worse than a zero one, test_bookkeeping.py).
  double measured_pe() const { return measured_pe_; }
  double measured_pz() const { return measured_pz_; }
  double measured_pzz() const { return measured_pzz_; }

  /// {category name -> integrated luminosity [pb^-1]} at `total_lumi_pb`.
  /// Pure `lumi_fraction * total`: the ONE place a share becomes counts.
  std::map<std::string, double> lumi_shares(double total_lumi_pb) const;

  /// Positional variant (same order as `categories()`), for hot loops.
  std::vector<double> lumi_share_vector(double total_lumi_pb) const;

  /// Index of a named category, or -1.
  int index_of(const std::string& name) const;

 private:
  std::vector<SpinCategory> categories_;
  double pe_true_ = 0.0, pz_true_ = 0.0, pzz_true_ = 0.0;
  double delta_p_over_p_ = 0.0;
  std::uint64_t polarimetry_seed_ = 20260713;
  double measured_pe_ = 0.0, measured_pz_ = 0.0, measured_pzz_ = 0.0;
};

// --------------------------------------------------------- standard plans

/// Options of `helicity_flip_plan`.  `use_explicit_pzz = false` (the default)
/// is the Python's `pzz=None`: SPIN-TEMPERATURE (max-entropy) populations for
/// the requested pz -- the physical vector-fill model, which for j >= 1 drags
/// a non-zero rank-2 moment along with it (recorded in `pzz_true`).
struct HelicityFlipOptions {
  bool use_explicit_pzz = false;
  /// P_zz in [-2, 1] for j = 1, the normalized T in [-1, 1] for j = 3/2.
  double pzz = 0.0;
  double theta_s = 0.0;
  double phi_s = 0.0;
  /// Boosts the lam_e = +1 share only: L+ = L/2 (1+offset), L- = L/2 --
  /// the convention `apar_rel_lumi_bias` is written for.
  double rel_lumi_offset = 0.0;
  std::string name = "apar";
};

/// A_par run plan: one vector fill, electron helicity flipped +-.
/// Throws for j outside {1/2, 1, 3/2} or an unphysical (pz, pzz).
RunPlan helicity_flip_plan(double j, double pz, double pe,
                           const HelicityFlipOptions& opt = {});

/// Spin-1 A_zz run plan: equal-thirds fills (pz, +pzz), (-pz, +pzz) and the
/// m0-enriched (0, -2 pzz) -- the HERMES-style pattern that makes the
/// canonical thirds estimator exact.  Electron unpolarized.
/// `rel_lumi_offset` boosts the 0-enriched share only: L0 = L/3 (1+offset),
/// L+- = L/3 -- the convention of `azz_rel_lumi_bias`.
RunPlan tensor_thirds_plan(double pz, double pzz, double rel_lumi_offset = 0.0,
                           double theta_s = 0.0, double phi_s = 0.0,
                           const std::string& name = "azz");

/// Gluonometry run plan: tensor-polarized fill, axis in the transverse plane
/// (theta_S = 90 deg), electron unpolarized.  Single category; the cos(2 phi')
/// amplitude is extracted from the phi distribution.
RunPlan transverse_tensor_plan(double pzz, double phi_s = 0.0,
                               const std::string& name = "cos2phi");

/// Gluonometry run plan with the acceptance-cancelling spin-state pattern:
/// m = +-1-rich bunches at (P_z, P_zz) = (0, +pzz) and m = 0-rich bunches at
/// (0, -2 pzz) -- the same source purity as `tensor_thirds_plan` -- both with
/// the alignment axis transverse at azimuth phi_s (default pi/2, the vertical
/// stable-spin direction), unpolarized electrons.  The two states are meant to
/// alternate BUNCH BY BUNCH: the ratio estimator cancels only an acceptance
/// common to both samples.  `share_plus` is the luminosity share of the +pzz
/// sample; `rel_lumi_offset` boosts that share by (1+offset).
RunPlan tensor_flip_plan(double pzz, double phi_s = kPi / 2.0,
                         double share_plus = 0.5, double rel_lumi_offset = 0.0,
                         const std::string& name = "flip");

/// Copy of `plan` with ONE category's lumi share scaled by (1+offset), the
/// others untouched -- the convention the bias formulas below are written
/// for.  Throws if the name is not in the plan.
RunPlan with_offset(const RunPlan& plan, const std::string& category_name,
                    double offset);

// ------------------------------------------------- relative-luminosity bias

/// First-order naive-thirds bias from an offset on the m0-enriched share:
/// bias(Azz) = -(2/3) offset / pzz.
double azz_rel_lumi_bias(double offset, double pzz);

/// First-order helicity-flip bias: bias(A_par) = offset / (2 pe pz).
double apar_rel_lumi_bias(double offset, double pe, double pz);

// --------------------------------------------------- spin-temperature ladder

/// The geometric (spin-temperature) population ladder p_m ~ t^m with
/// <J_z>/J = pz, solved for t by bisection in log t.  This is the closed-form
/// construction behind the `pzz=None` branch of `helicity_flip_plan`; at
/// t = 3 it is rational -- (9,3,1)/13 for J = 1 and (27,9,3,1)/40 for J = 3/2,
/// giving (P_z, rank-2) = (8/13, 4/13) and (7/10, 2/5).
struct SpinTemperatureLadder {
  std::vector<double> populations;  ///< m = +J ... -J
  double t = 1.0;                   ///< the ladder ratio exp(beta)
};
SpinTemperatureLadder spin_temperature_ladder(double j, double pz,
                                              int iterations = 400);

/// Rank-2 moment of the spin-temperature fill at vector polarization pz --
/// the number `helicity_flip_plan` records as `pzz_true` when no explicit
/// pzz is given.  Zero (exactly) below spin 1.
double spin_temperature_pzz(double j, double pz);

// ------------------------------------------------------------ RNG streams

/// Fixed stream per (run, bunch): the reproducibility discipline of
/// plans/05 5.2.4.  `event` selects the counter within the bunch, so the
/// same (seed, run, bunch, event) is the same event on any thread count.
inline Rng bunch_rng(std::uint64_t seed, std::uint64_t run, std::uint64_t bunch,
                     std::uint64_t event = 0) {
  return Rng(seed, run, bunch, event);
}

/// The counter reserved for per-bunch decisions that are not per-event
/// (the Poisson event count of a category).  Kept out of the event range so
/// that adding events never changes the count draw.
inline constexpr std::uint64_t kCountStreamEvent = ~static_cast<std::uint64_t>(0);

}  // namespace lipolgen

#endif  // LIPOLGEN_BOOKKEEPING_HPP
