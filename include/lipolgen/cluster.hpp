// SPDX-License-Identifier: GPL-3.0-or-later
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

#include <array>
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

/// A plain `x v1 v2 ... vN` table with NO MC-error columns, introduced by the
/// same `****  *********  *********` column rule `read_anl_momentum` keys on
/// (`li6_alpha_d/li6.adr.fit`, the `*.rho1` fits).
///
/// It exists because `read_anl_momentum` CANNOT read those files: their rows
/// also carry three numbers, so that reader would silently land the second
/// VALUE column (`R2LI6FIT`) in `err[0]` and report one wave where there are
/// two.  `AnlTable::err` is left empty here -- these files print no errors.
/// Throws if the block is missing or if a row is short of `1 + ncol` numbers.
AnlTable read_anl_plain(const std::string& path, std::size_t ncol);

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

// -------------------------------------------------------- the CD-Bonn deuteron

/// CD-Bonn's own deuteron binding energy: R. Machleidt, Phys. Rev. C 63,
/// 024001 (2001) (arXiv:nucl-th/0006014), Table XV (LaTeX label `tab_deu`).
///
/// This is NOT a second copy of any binding energy already in the library.
/// `fdeut.av18`'s header `ebind` = 2.224574 MeV is AV18's and stays where it
/// is; `DEUTERON_P_TAG().separation_energy` is the tagged channel's rounded
/// 2.2246e-3 GeV.  A wave function's epsilon_d must belong to the wave
/// function it is convolved with, so `cdbonn_fdeut_table` puts THIS one in
/// `FdeutTable::ebind_gev`.
///
/// It is TYPED, not derived from `CD_BONN_GAMMA_FM`.  gamma does imply it, to
/// 5e-7, but only through Machleidt's relativistic Eq. (D11); the
/// non-relativistic gamma = sqrt(M_N B_d)/hbar c misses the published gamma
/// by 3.0e-4 (docs/open_items/run_2026-09-03/phase_A_cdbonn.md section 5.2),
/// and reproducing (D11) would need Machleidt's own hbar c, M_p and M_n --
/// a third hbar c, which docs/CONVENTIONS.md forbids.
inline constexpr double CD_BONN_BINDING_MEV = 2.224575;

/// THE PUBLISHED CD-BONN COEFFICIENTS, verbatim.
///
/// Source: R. Machleidt, "The high-precision, charge-dependent Bonn
/// nucleon-nucleon potential (CD-Bonn)", Phys. Rev. C 63, 024001 (2001);
/// e-print arXiv:nucl-th/0006014, Table XX (LaTeX label `tab_dwpar`,
/// "Coefficients for the parametrized deuteron wave functions (n = 11)"),
/// with the parameterisation in Appendix D, Eqs. (D19)-(D25) and gamma in
/// Eq. (D6).  Read out of the e-print's LaTeX source and cross-read out of
/// `pdftotext` of the PDF, digit for digit -- see
/// docs/open_items/run_2026-09-03/phase_A_cdbonn.md sections 1 and 3, which
/// also reproduces both listings.
///
/// The table publishes C_1..C_10 and D_1..D_8 ONLY.  C_11 and D_9..D_11 are
/// determined by the r -> 0 boundary conditions and are computed in
/// `cdbonn_wave()`; see `CdBonnWave` below.
inline constexpr double CD_BONN_GAMMA_FM = 0.2315380;   ///< gamma [fm^-1], Eq. (D6)
inline constexpr double CD_BONN_M0_FM = 0.9;            ///< m_0 [fm^-1], Eq. (D25)
inline constexpr int CD_BONN_N = 11;                    ///< n, Table XX's caption
inline constexpr double CD_BONN_C[10] = {                            // fm^-1/2
    0.88472985e+00, -0.26408759e+00, -0.44114404e-01, -0.14397512e+02,
    0.85591256e+02, -0.31876761e+03,  0.70336701e+03, -0.90049586e+03,
    0.66145441e+03, -0.25958894e+03};
inline constexpr double CD_BONN_D[8] = {                             // fm^-1/2
    0.22623762e-01, -0.50471056e+00,  0.56278897e+00, -0.16079764e+02,
    0.11126803e+03, -0.44667490e+03,  0.10985907e+04, -0.16114995e+04};

/// The CD-Bonn deuteron wave function, r- and momentum-space parameterisation
/// of Machleidt, PRC 63, 024001 (2001), **Appendix D**, Eqs. (D19)-(D25),
/// with the coefficients of **Table XX** (LaTeX label `tab_dwpar`).
///
/// WHY IT IS HERE.  Gate condition 3 of open item 10
/// (docs/OPEN_ITEMS_SOLUTIONS.md) asks the A = 2 gate for a REAL CD-Bonn u, w
/// rather than the D-state RESCALING PROXY of item 5 (AV18's w scaled to
/// P_D = 4.85 %), because CDKS built the Fig. 4 curve the gate compares
/// against on CD-Bonn.  The proxy has the wrong SIGN of the effect: measured
/// through this repository's own gate, CD-Bonn moves the G3b peak ratio the
/// opposite way from the proxy.  The reason is that CD-Bonn's D wave is not a
/// rescaled AV18 D wave: w(CD-Bonn)/w(AV18) runs 1.02 -> 0.99 -> 0.36 over
/// p = 0.1, 1.0, 5.0 fm^-1, so ONE factor (which is all the proxy is) is
/// ~11 % too small where the D wave is largest and 2.6x too large in the
/// tail.  On top of that CD-Bonn's first S node sits 13 % higher in k
/// (between 2.3 and 2.4 fm^-1 against AV18's 2.0-2.1) and its u is half as
/// big beyond it, which nearly removes the NEGATIVE lobe of the S-D
/// interference that partly cancels the positive one in AV18 -- and a
/// rescaling can never see that, because it keeps AV18's node.  The
/// measurements are in
/// docs/open_items/run_2026-09-03/phase_A_cdbonn.md (sections 6-8) and the
/// gate rows in .../phase_A_numbers.md section 8.
///
/// A TABLE-NUMBERING WARNING.  Several planning documents in this repository
/// say "Tables XVII/XVIII".  In the e-print those are the scalar-isoscalar
/// boson parameters; the deuteron material is Table XIX (numerical u, w) and
/// **Table XX** (these coefficients).  The LaTeX labels `tab_dwpar` /
/// `tab_dwaves` are the unambiguous handles and are what is cited here; the
/// published journal's numbering has not been checked against the e-print's.
///
/// THE PARAMETERISATION, n = 11, m_0 = 0.9 fm^-1, gamma from Eq. (D6):
///
///   u_a(r)     = sum_j C_j exp(-m_j r)                              (D19)
///   w_a(r)     = sum_j D_j exp(-m_j r) [1 + 3/(m_j r) + 3/(m_j r)^2](D20)
///   psi_0^a(q) = sqrt(2/pi) sum_j C_j/(q^2 + m_j^2)                 (D21)
///   psi_2^a(q) = sqrt(2/pi) sum_j D_j/(q^2 + m_j^2)                 (D22)
///   m_j        = gamma + (j - 1) m_0                                (D25)
///
/// ONLY C_1..C_10 and D_1..D_8 ARE PUBLISHED.  C_11 and D_9..D_11 are fixed
/// by the boundary conditions u_a(r) -> r and w_a(r) -> r^3 as r -> 0
/// (Eqs. (D23)/(D24)) and are COMPUTED in `cdbonn_wave()`, never typed --
/// typing them would be a second definition of a physics quantity, which
/// docs/CONVENTIONS.md forbids, and the four residuals are gated in
/// tests/test_cluster.cpp.  This file uses the equivalent form of those
/// constraints, obtained by expanding
/// e^{-x}(1 + 3/x + 3/x^2) = 3/x^2 - 1/2 + x^2/8 - ... (the 1/x term cancels
/// identically), which is three SUM RULES rather than Eq. (D24)'s three
/// circular permutations and so cannot be mis-permuted:
///
///   sum_j C_j = 0 ,  sum_j D_j/m_j^2 = 0 ,  sum_j D_j = 0 ,
///   sum_j D_j m_j^2 = 0 .
///
/// Machleidt's own warning, quoted because it is a real trap: "The
/// constraints Eqs. (D23) and (D24) must be enforced by double precision
/// (i.e., to about 15 decimal digits), otherwise the wave function is not
/// reproduced correctly for r <= 0.5 fm.  This applies, particularly, to the
/// D wave."
///
/// THE SIGN -- the (-i)^L trap that b1_nuclear.hpp's `ClusterPartialWave`
/// already warns about, and CD-Bonn walks straight into it.  Machleidt's
/// Eq. (D13) is printed with a BARE j_L kernel for both L; taken together
/// with (D20)/(D22) that is inconsistent at L = 2 -- the bare-j_2 transform
/// of psi_2^a comes back as MINUS w_a(r), while the L = 0 pair round-trips
/// (measured, phase_A_cdbonn.md section 4).  `fdeut.av18`'s own r and k
/// blocks ARE related by the bare kernel for BOTH L.  So the drop-in
/// convention, and what `psi_s`/`psi_d` return, is
///
///   u(p) = + psi_0^a(p) ,     w(p) = - psi_2^a(p) ,
///
/// i.e. exactly `FdeutTable`'s `u`, `w`, hence the CDKS phi_0 = u,
/// phi_2 = -w.  Three independent confirmations are recorded in
/// phase_A_cdbonn.md sections 4 and 8.3; the cheapest is
/// `alpha_d_quadrupole_fm2`, which returns +0.2702 fm^2 with this sign
/// (CD-Bonn's published Q_d = 0.270 fm^2) and -0.3037 with the other.
///
/// UNITS are the paper's throughout this struct: p in fm^-1, C_j, D_j in
/// fm^-1/2, psi in fm^3/2, and the normalisation is
/// (2/pi) int_0^inf dp p^2 (u^2 + w^2) = 1 (Eq. D14).  The single conversion
/// to the library's GeV happens in `cdbonn_fdeut_table`, with
/// `HBARC_GEV_FM`, exactly as `read_fdeut_k` does it.
struct CdBonnWave {
  std::array<double, 11> m{};   ///< m_j [fm^-1], Eq. (D25)
  std::array<double, 11> c{};   ///< C_j [fm^-1/2]; c[10] is COMPUTED
  std::array<double, 11> d{};   ///< D_j [fm^-1/2]; d[8..10] are COMPUTED

  /// u(p) = +psi_0^a(p) [fm^3/2], p in fm^-1.
  double psi_s(double p_fm) const;
  /// w(p) = -psi_2^a(p) [fm^3/2] -- SIGNED as above, p in fm^-1.
  double psi_d(double p_fm) const;

  /// The momentum-space moments in CLOSED FORM.  With
  /// (2/pi) int_0^inf dp p^2/[(p^2+a^2)(p^2+b^2)] = 1/(a+b) the two norms are
  /// the double sums sum_ij C_i C_j/(m_i+m_j) and sum_ij D_i D_j/(m_i+m_j),
  /// so they carry NO quadrature error at all -- which is what makes them a
  /// sharp test of the coefficients rather than of an integrator.
  double norm_s() const;   ///< (2/pi) int dp p^2 u^2
  double norm_d() const;   ///< (2/pi) int dp p^2 w^2 -- the D-state probability
  double norm() const;     ///< the sum; = 1 by Eq. (D14)

  /// The asymptotics of Eq. (D15).  m_1 = gamma is the slowest-decaying mass,
  /// so u_a -> C_1 e^{-gamma r} and w_a -> D_1 e^{-gamma r}[1 + 3/(gamma r) +
  /// 3/(gamma r)^2], i.e. A_S = C_1 and A_D = D_1 with no integral needed.
  double a_s() const;      ///< A_S = C_1
  double a_d() const;      ///< A_D = D_1
  double eta() const;      ///< A_D/A_S, the asymptotic D/S ratio

  /// The four boundary-condition sums, in the order
  /// {sum C_j, sum D_j/m_j^2, sum D_j, sum D_j m_j^2}.  All four are zero by
  /// construction; the test measures HOW zero, which is Machleidt's
  /// "about 15 decimal digits" demand made mechanical.
  std::array<double, 4> constraint_residuals() const;
};

/// The one CD-Bonn wave function: Table XX with C_11 and D_9..D_11 solved
/// for.  Built once, on first use.
const CdBonnWave& cdbonn_wave();

/// CD-Bonn sampled on `k = 0, dk, .. k_max` [fm^-1] and returned as an
/// `FdeutTable`, i.e. converted to the library's GeV units with
/// `HBARC_GEV_FM` in BOTH the abscissa and the fm^3/2 ordinate (the same
/// "both or neither" rule `read_fdeut_k` states), and carrying
/// `CD_BONN_BINDING_MEV` as `ebind_gev`.  So a caller can swap
/// `read_fdeut_k(...)` for this and change the WAVE FUNCTION and nothing
/// else.
///
/// The defaults are `fdeut.av18`'s OWN grid -- 0 to 20 fm^-1 in steps of 0.1,
/// 201 rows -- deliberately, so that an AV18 row and a CD-Bonn row of the
/// gate differ in the wave function alone: same spline, same node count, same
/// 20 fm^-1 truncation.  Refining to dk = 0.02 moves the gate's peak by 5e-4
/// relative (phase_A_cdbonn.md section 8.4), so the coarse grid is not what
/// drives any of it.  Throws if `dk_fm <= 0` or `k_max_fm < dk_fm`.
FdeutTable cdbonn_fdeut_table(double k_max_fm = 20.0, double dk_fm = 0.1);

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
  /// With the table's own 1-sigma Monte Carlo error on psi.  `dpsi` must be
  /// the same length as `psi` or empty; it is SIGNED and transforms exactly
  /// like psi (see `scaled`), so psi + n dpsi is covariant under the global
  /// phase convention.
  VmcRadial(std::vector<double> k_gev, std::vector<double> psi,
            std::vector<double> dpsi, int l, std::string provenance);

  /// psi_L(k), linearly interpolated; 0 outside [k.front(), k.back()].
  double operator()(double k) const;

  int l() const { return l_; }
  bool empty() const { return k_.empty(); }
  const std::vector<double>& k() const { return k_; }
  const std::vector<double>& psi() const { return psi_; }
  /// The 1-sigma Monte Carlo error of `psi`, point by point -- EMPTY when the
  /// source file printed none (`li6.adr.fit`, `fdeut.av18`, CD-Bonn).
  ///
  /// THE ANL FILES CARRY THIS AND THE LIBRARY THREW IT AWAY until 2026-09-04
  /// (open item C5.2): `AnlTable::err` was parsed by `read_anl_momentum` /
  /// `read_anl_overlap` and then dropped on the floor by the two factories
  /// below, so `ClusterWaveSource::VmcAV18` had no error band of any kind.
  /// For a momentum file psi = s sqrt(rho), hence dpsi = s drho / (2 sqrt rho)
  /// and 0 wherever rho <= 0.
  ///
  /// IT IS NOT AN INDEPENDENT-POINT ERROR.  Every point of an ANL table comes
  /// from the SAME variational Monte Carlo walk, so the point-to-point
  /// correlation is unknown and is neither 0 nor 1.  `shifted_by_sigma` takes
  /// the FULLY CORRELATED reading (every point moved by n sigma in the same
  /// direction), which is the conservative envelope on any smooth functional;
  /// adding the same errors in quadrature gives the other extreme.  Quote
  /// both, or quote the correlated one and say so.
  const std::vector<double>& dpsi() const { return dpsi_; }
  bool has_errors() const { return !dpsi_.empty(); }
  const std::string& provenance() const { return provenance_; }

  /// The same table times `factor` -- used to fix the (unobservable) global
  /// phase of a multi-wave channel to psi_{Lmin}(k -> 0) > 0.  `dpsi` is
  /// scaled by the SAME factor, so a phase flip carries the band with it.
  VmcRadial scaled(double factor) const;
  /// psi -> psi + n_sigma * dpsi, fully correlated across the table.  Returns
  /// *this unchanged at n_sigma == 0 or with no errors, so the band's zero row
  /// is bit for bit today's table.
  VmcRadial shifted_by_sigma(double n_sigma) const;
  /// integral k^2 psi^2 dk over the table [GeV^3 x psi^2].
  double norm2() const;
  /// The 1-sigma Monte Carlo error of `norm2()`, in the two limits the errors
  /// admit: `correlated` = |d/dn integral k^2 (psi + n dpsi)^2 dk| at n = 0 =
  /// |2 integral k^2 psi dpsi dk|, `quadrature` = sqrt(sum of the same
  /// integrand's per-cell contributions squared).  Both are 0 with no errors.
  void norm2_error(double* correlated, double* quadrature) const;

 private:
  std::vector<double> k_, psi_, dpsi_;
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
/// The run-surface name of a wave source -- the same two strings the CLI's
/// `--cluster-wave` takes (`python/lipolgen/__init__.py`'s `CLUSTER_WAVES`),
/// defined ONCE here so a `meta` key, a `KnobProvenance` row and a banner line
/// cannot spell them differently (docs/CONVENTIONS.md).  Every other selector
/// in the library already has one (`b1_model_name`, `unpol_sf_name`, ...);
/// this enum and `TritonSfChoice` were the two that did not.
const char* cluster_wave_name(ClusterWaveSource s);

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
