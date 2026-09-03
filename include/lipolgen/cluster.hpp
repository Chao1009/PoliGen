#ifndef LIPOLGEN_CLUSTER_HPP
#define LIPOLGEN_CLUSTER_HPP

/// \file cluster.hpp
/// Cluster relative partial waves: the radial shapes and the azimuth-stripped
/// spherical harmonics the two-cluster amplitude of `tagged.hpp` is built from.
/// C++17 port of `tagged.Wave` / `tagged.theta_lm`
/// (evgen/polligen/tagged.py:126,101).
///
/// The radial forms are UNNORMALIZED; `TaggedModel` normalizes them on its own
/// grid so that integral psihat^2 k^2 dk = 1.  kappa comes from the channel's
/// separation energy (spectator.hpp), beta is the short-range scale scanned as
/// the model band (0.20 - 0.40, default 0.30):
///
///   L = 0  Hulthen   1/(k^2+kappa^2) - 1/(k^2+beta^2)
///   L = 1  P-wave    k / ((k^2+kappa^2)(k^2+beta^2))
///   L = 2  D-wave    k^2 / ((k^2+kappa^2)(k^2+beta^2)^2)
///
/// The L = 0 form is EXACTLY `spectator.momentum_density`'s psi, which is what
/// makes the P_D = 0 limit of the tagged sampler agree with the fast
/// simulation's spectator sampler quantile by quantile.
///
/// The VMC (R. B. Wiringa et al., ANL) two-cluster tables land as `VmcRadial`
/// behind the same `radial()` interface -- see `ClusterWaveSource` below and
/// `docs/open_items/vmc_reconciliation.md` for which file feeds which wave.

#include <memory>
#include <string>
#include <vector>

#include "lipolgen/spectator.hpp"

namespace lipolgen {

// ------------------------------------------------------------- data lookup

/// Root of the packaged data tree.  `$LIPOLGEN_DATA_DIR` when it is set and
/// non-empty, else the compiled-in default (CMake passes
/// `-DLIPOLGEN_DATA_DIR_DEFAULT="${CMAKE_SOURCE_DIR}/data"`, so an in-tree
/// build finds `data/` with no environment at all).
const std::string& data_dir();
/// `data_dir() + "/" + relative`.
std::string data_path(const std::string& relative);

// ---------------------------------------------------- the ANL table readers

/// One block of an ANL Fortran-output table: an abscissa plus N value columns
/// and (where the file carries them) their 1-sigma Monte Carlo errors.
struct AnlTable {
  std::vector<double> x;                      ///< k [fm^-1] or r [fm]
  std::vector<std::vector<double>> col;       ///< value columns
  std::vector<std::vector<double>> err;       ///< 1-sigma, empty if absent
};

/// The two blocks of an `overlap_old/` raw VMC output (`li6.ad`, `li7.at`):
/// `[0]` is the k-space amplitude block (header `k(fm-1)`), `[1]` the r-space
/// amplitude block (header `rij`).  Both carry two signed columns with MC
/// errors.  Throws if either block is missing.
std::vector<AnlTable> read_anl_overlap(const std::string& path);

/// The blocks of a `momenta/` momentum-distribution file
/// (`li6_ad1.momentum`, `li7_at3.momentum`): `[0]` is the total rho(K),
/// `[1]` (only where the file splits partial waves) is `rho_0, rho_2`.
/// Columns are DENSITIES rho_L = A_L^2 / (4 pi) [fm^3], never negative.
std::vector<AnlTable> read_anl_momentum(const std::string& path);

/// The file's own printed `4*PI*TOTINT(RHO*K**2:K)/(2*PI)**3 = S` lines, in
/// order (total first, then the per-wave ones).
std::vector<double> read_anl_momentum_norms(const std::string& path);

/// The k-space block of an ANL deuteron wave-function file
/// (`deuteron/fdeut.av18`) plus the binding energy its own header carries.
///
/// UNITS.  `k_gev` is the file's k [fm^-1] times `HBARC_GEV_FM`, and `u`, `w`
/// are the file's columns divided by `HBARC_GEV_FM^(3/2)`, so that the pair
/// is in ONE consistent system: integral k^2 (u^2 + w^2) dk = 1 with k in GeV
/// (0.999976 on the file's own 0.1 fm^-1 grid).  Converting the abscissa
/// without the fm^(3/2) of the ordinate would leave the norm off by
/// hbar c^3 = 7.7e-3, which is why the reader does both or neither.
///
/// `ebind_gev` is the header's `ebind` (2.224574 MeV -> 2.224574e-3), so the
/// convolution's epsilon_d has ONE source and no literal is re-typed;
/// `DEUTERON_P_TAG().separation_energy` = 2.2246e-3 is the same number
/// rounded (2.6e-5 relative).
struct FdeutTable {
  std::vector<double> k_gev;   ///< strictly increasing, 0 .. 20 fm^-1
  std::vector<double> u;       ///< S-wave phi_0(k) = u(k) >= 0 at low k
  std::vector<double> w;       ///< D-wave W(k); the CDKS phi_2 is -w
  double ebind_gev = 0.0;      ///< the header's `ebind` [GeV]
};

/// Read the `k  u(k)  w(k)` block and the `ebind` header of a `fdeut.*` file.
///
/// hbar c.  This reader converts with `HBARC_GEV_FM` (constants.hpp, 0.19733),
/// which `docs/CONVENTIONS.md` names as THE fm <-> GeV conversion, and so does
/// everything in `b1_nuclear.cpp`.  Its neighbours in this file
/// (`vmc_from_overlap_k`, `vmc_from_momentum`) use a file-local
/// 0.1973269804 instead; the two disagree at 1.5e-5 relative (20 keV at
/// k = 5 fm^-1), which is a PRE-EXISTING inconsistency recorded here rather
/// than silently propagated -- it is far inside every tolerance in
/// docs/open_items/run_2026-09-02/design_D_b1_li6.md section 5.
FdeutTable read_fdeut_k(const std::string& path);

// ----------------------------------------------------------------- VmcRadial

/// A tabulated, L-specific radial amplitude psi_L(k), linearly interpolated
/// and ZERO outside the tabulated range.
///
/// SIGN.  psi_L is signed and the sign is load-bearing: the S-D interference
/// term of `tagged.hpp`'s n_M(k, khat) goes as psi_0 psi_2, so flipping the
/// relative sign of one wave flips the tensor structure.  The ANL momentum
/// files carry |psi_L|^2 only and cannot supply it; the `overlap_old`
/// amplitude files can and do.  See `vmc_sign_steps` below.
///
/// The table is UNNORMALIZED, like every other radial form in this header --
/// `TaggedModel` normalizes each wave on its own grid.  That normalization is
/// by a POSITIVE factor sqrt(P_L)/norm, so it preserves signs.
class VmcRadial {
 public:
  VmcRadial() = default;
  VmcRadial(std::vector<double> k_gev, std::vector<double> psi, int l,
            std::string provenance = std::string());

  /// psi_L(k), linearly interpolated; 0 outside [k.front(), k.back()].
  double operator()(double k) const;

  int l() const { return l_; }
  bool empty() const { return k_.empty(); }
  const std::vector<double>& k() const { return k_; }
  const std::vector<double>& psi() const { return psi_; }
  const std::string& provenance() const { return provenance_; }

  /// The same table times `factor` -- used to fix the (unobservable) global
  /// phase of a multi-wave channel to psi_{Lmin}(k -> 0) > 0.
  VmcRadial scaled(double factor) const;
  /// integral k^2 psi^2 dk over the table [GeV^3 x psi^2].
  double norm2() const;

 private:
  std::vector<double> k_, psi_;
  int l_ = 0;
  std::string provenance_;
};

/// Sign changes of a tabulated amplitude, used to give a momentum-file
/// magnitude its phase.  Returns the interpolated zero crossings of
/// `amp` below `x_max` (in the abscissa's own units).  Crossings above
/// `x_max` are deliberately NOT returned: past ~3 fm^-1 both VMC columns are
/// at the Monte Carlo noise floor and their sign wanders, while the strength
/// there is ~1e-4 of the norm, so a noise-driven flip would be all cost and
/// no signal.
std::vector<double> vmc_sign_steps(const std::vector<double>& x,
                                   const std::vector<double>& amp,
                                   double x_max);

/// psi_L(k) straight from an `overlap_old` k-space column (already a signed
/// amplitude).  `column` is 0 or 1.
VmcRadial vmc_from_overlap_k(const std::string& path, int column, int l);

/// psi_L(k) = 4 pi integral A_L(r) j_L(k r) r^2 dr from an `overlap_old`
/// r-space column.  This is the ANL cluster-overlap convention, which is
/// (2 pi)^(3/2) times the reduced-radial sqrt(2/pi) one; it is validated to
/// 1.000 against `deuteron/fdeut.av18`'s own tabulated u(k)/w(k) and to 1e-3
/// against the file's independent k-space block (validation/vmc_reconcile.py).
VmcRadial vmc_from_overlap_r(const std::string& path, int column, int l,
                             double k_max_fm = 5.0, std::size_t nk = 251);

/// psi_L(k) = s(k) sqrt(rho_L(K)) from a `momenta` column, with the sign
/// structure s(k) taken from `sign_from` (an overlap-derived table of the
/// same L) through `vmc_sign_steps`.  A null `sign_from` leaves psi positive,
/// which is correct only for a single-wave channel, where the overall phase
/// is unobservable.
VmcRadial vmc_from_momentum(const std::string& path, int block, int column,
                            int l, const VmcRadial* sign_from,
                            double node_search_max_fm = 3.0);

/// Which family of radial forms a cluster channel is built from.
enum class ClusterWaveSource {
  Hulthen = 0,   ///< the analytic two-parameter forms above -- THE DEFAULT
  VmcAV18 = 1,   ///< the ANL VMC tables (AV18+UX magnitudes, AV18+UIX signs)
};

/// One cluster relative partial wave: orbital L, probability P_L = a_L^2, and
/// the short-range scale beta of its radial shape.
struct Wave {
  int l = 0;
  double prob = 1.0;
  double beta = BETA_DEFAULT;
  /// Null => the analytic Hulthen switch; non-null => the table.  Shared
  /// because a `TaggedChannel` is copied freely and the tables are ~kB each.
  std::shared_ptr<const VmcRadial> vmc;

  /// Unnormalized radial shape psi_L(k).  Dispatches to `vmc` when it is set
  /// (and `kappa` is then unused -- the table already encodes the physical
  /// falloff), else to the analytic switch.  Throws for L > 2 on the
  /// analytic branch.
  double radial(double k, double kappa) const;
  /// The same on a grid.
  std::vector<double> radial(const std::vector<double>& k, double kappa) const;
};

/// |m|-azimuth-stripped spherical harmonic Theta_l^m(theta) with
/// Condon-Shortley signs: Y_l^m = Theta_l^m(theta) exp(i m phi).  `c` is
/// cos(theta).  Implemented for l <= 2 (throws otherwise), which is every
/// wave the two-cluster expansion carries.
double theta_lm(int l, int m, double c);

}  // namespace lipolgen

#endif  // LIPOLGEN_CLUSTER_HPP
