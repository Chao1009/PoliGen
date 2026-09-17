// SPDX-License-Identifier: GPL-3.0-or-later
// pybind11 bindings for LiPolGen (task P7).
//
// Module `_lipolgen`, imported by the `lipolgen` package
// (python/lipolgen/__init__.py).  The layout mirrors the C++ headers one
// section per header; nothing physical is defined here, every number comes
// from the core library.
//
// NUMPY CONVENTIONS
//   * Every array returned by a *generation* entry point (`Pipeline.generate`,
//     `InclusiveSampler.sample_n/sample_lumi`, `hfs_arrays`) is a fresh
//     `numpy.ndarray` that OWNS its buffer: the C++ `std::vector` is moved into
//     a capsule and handed to numpy, so there is no copy and no lifetime tie to
//     the generator object.  This is the only place zero-copy is used.
//   * Every array returned by an *accessor* (`InclusiveSampler.x_cells`, ...)
//     is a copy, so it can outlive the object it came from.
//
// GIL
//   `Pipeline.generate`, `Pipeline.for_each_count`, `InclusiveSampler.sample_*`
//   and `PythiaBridge.hadronize` release the GIL around the C++ work.  A
//   hadronizer hook set from Python re-acquires it per event (pybind11 does
//   that itself); use `PipelineConfig.set_pythia_hadronizer` to bind the C++
//   bridge directly and stay GIL-free.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/operators.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "lipolgen/asymmetries.hpp"
#include "lipolgen/b1_nuclear.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/cluster.hpp"
#include "lipolgen/cluster_config.hpp"
#include "lipolgen/coherent.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/fsi.hpp"
#include "lipolgen/generator.hpp"
#include "lipolgen/numerics.hpp"
#include "lipolgen/pipeline.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/spectator.hpp"
#include "lipolgen/spin.hpp"
#include "lipolgen/tagged.hpp"
#include "lipolgen/triton_sf.hpp"
#include "lipolgen/xsec.hpp"

#ifndef LIPOLGEN_HAVE_LHAPDF
#define LIPOLGEN_HAVE_LHAPDF 0
#endif
#ifndef LIPOLGEN_HAVE_HEPMC3
#define LIPOLGEN_HAVE_HEPMC3 0
#endif
#ifndef LIPOLGEN_HAVE_PYTHIA8
#define LIPOLGEN_HAVE_PYTHIA8 0
#endif

#if LIPOLGEN_HAVE_HEPMC3
#include "lipolgen/hepmc_writer.hpp"
#endif
#if LIPOLGEN_HAVE_LHAPDF
#include "lipolgen/lhapdf_sf.hpp"
#endif
#if LIPOLGEN_HAVE_PYTHIA8
#include "lipolgen/mstw_sf.hpp"
#include "lipolgen/pythia_bridge.hpp"
#endif

namespace py = pybind11;
using namespace lipolgen;

namespace {

// ------------------------------------------------------------ numpy helpers

/// Move a std::vector into a 1-D numpy array (no copy; numpy owns the buffer).
template <typename T>
py::array_t<T> move_array(std::vector<T>&& v) {
  auto* held = new std::vector<T>(std::move(v));
  py::capsule owner(held, [](void* p) {
    delete reinterpret_cast<std::vector<T>*>(p);
  });
  return py::array_t<T>(static_cast<py::ssize_t>(held->size()), held->data(),
                        owner);
}

/// Move a flat std::vector into an (n, ncol) numpy array.
template <typename T>
py::array_t<T> move_array2(std::vector<T>&& v, py::ssize_t ncol) {
  auto* held = new std::vector<T>(std::move(v));
  py::capsule owner(held, [](void* p) {
    delete reinterpret_cast<std::vector<T>*>(p);
  });
  const py::ssize_t n = ncol > 0 ? static_cast<py::ssize_t>(held->size()) / ncol : 0;
  return py::array_t<T>({n, ncol},
                        {static_cast<py::ssize_t>(sizeof(T)) * ncol,
                         static_cast<py::ssize_t>(sizeof(T))},
                        held->data(), owner);
}

/// Move a flat std::vector into an (n, d1, d2) numpy array.  `move_array2`
/// is 2-D only, and `ClusterConfigSet::positions()` is (N, 6, 3).
template <typename T>
py::array_t<T> move_array3(std::vector<T>&& v, py::ssize_t d1, py::ssize_t d2) {
  auto* held = new std::vector<T>(std::move(v));
  py::capsule owner(held, [](void* p) {
    delete reinterpret_cast<std::vector<T>*>(p);
  });
  const py::ssize_t stride = static_cast<py::ssize_t>(sizeof(T));
  const py::ssize_t n =
      d1 * d2 > 0 ? static_cast<py::ssize_t>(held->size()) / (d1 * d2) : 0;
  return py::array_t<T>({n, d1, d2},
                        {stride * d1 * d2, stride * d2, stride},
                        held->data(), owner);
}

/// Copy a std::vector into a 1-D numpy array.
template <typename T>
py::array_t<T> copy_array(const std::vector<T>& v) {
  return py::array_t<T>(static_cast<py::ssize_t>(v.size()), v.data());
}

py::array_t<double> matrix_to_array(const RealMatrix& m) {
  const py::ssize_t n = static_cast<py::ssize_t>(m.size());
  py::array_t<double> out({n, n});
  auto r = out.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < n; ++i)
    for (py::ssize_t j = 0; j < n; ++j)
      r(i, j) = m(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
  return out;
}

py::array_t<std::complex<double>> matrix_to_array(const CplxMatrix& m) {
  const py::ssize_t n = static_cast<py::ssize_t>(m.size());
  py::array_t<std::complex<double>> out({n, n});
  auto r = out.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < n; ++i)
    for (py::ssize_t j = 0; j < n; ++j)
      r(i, j) = m(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
  return out;
}

using CplxArray =
    py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast>;

CplxMatrix array_to_cplx_matrix(const CplxArray& a) {
  if (a.ndim() != 2 || a.shape(0) != a.shape(1))
    throw std::runtime_error("expected a square 2-D complex array");
  CplxMatrix m(static_cast<std::size_t>(a.shape(0)));
  auto r = a.unchecked<2>();
  for (py::ssize_t i = 0; i < a.shape(0); ++i)
    for (py::ssize_t j = 0; j < a.shape(1); ++j)
      m(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) = r(i, j);
  return m;
}

double pseudorapidity(const Vec4& p) {
  const double pmag = p.p();
  if (pmag <= std::abs(p.pz)) return p.pz >= 0.0 ? 1e30 : -1e30;
  return 0.5 * std::log((pmag + p.pz) / (pmag - p.pz));
}

/// Mass number of a PDG code: nucleon codes are A = 1, 10-digit ion codes
/// 10LZZZAAAI carry A in digits 4..6.
int pdg_mass_number(int pdg) {
  const int a = std::abs(pdg);
  if (a == 2212 || a == 2112) return 1;
  if (a > 1000000000) return (a / 10) % 1000;
  return 0;
}

// ------------------------------------------------- columnar event generation

/// Every column `Pipeline.generate` fills.  Preallocated to the event count so
/// that worker threads can write disjoint row ranges without any locking.
struct Columns {
  std::vector<double> x, q2, y, phi, w2, nu, s, weight;
  std::vector<double> m_ion, m_struck, pe, theta_s, phi_s, pz, pzz, j;
  std::vector<double> k, cos_theta_k, phi_k, alpha_s, pt_s, t, x_pom;
  std::vector<double> e_prime, theta_e, eta_e, kp;   // kp flat (n, 4)
  std::vector<double> spec_pT, spec_theta, spec_p_lab, spec_R, spec_xL;
  std::vector<double> spec_kx, spec_ky, spec_kz, spec_phi;
  std::vector<double> struck_pol, struck_virtuality;
  std::vector<std::int64_t> number, lam_e, category_index, route;
  std::vector<std::int64_t> cell, n_partner, struck_pdg;
  /// rc.hpp, slot 0 (the event's own pure spin state).  Sized only when the
  /// run has RC on -- `columns_to_dict` then emits three keys that an
  /// `--rc off` npz does not have at all.
  std::vector<double> rc_lo, rc_hi, rc_tail;
  /// `Event::rc_clipped`, per event: bit kRcClipTail | bit kRcClipBand.
  /// Preallocated so worker threads write disjoint rows with no locking; the
  /// two run-level fractions in `meta` are reduced from it afterwards.
  std::vector<std::uint8_t> rc_clip;
  bool with_rc = false;

  void resize(std::size_t n) {
    for (auto* v : {&x, &q2, &y, &phi, &w2, &nu, &s, &weight, &m_ion,
                    &m_struck, &pe, &theta_s, &phi_s, &pz, &pzz, &j, &k,
                    &cos_theta_k, &phi_k, &alpha_s, &pt_s, &t, &x_pom,
                    &e_prime, &theta_e, &eta_e, &spec_pT, &spec_theta,
                    &spec_p_lab, &spec_R, &spec_xL, &spec_kx, &spec_ky,
                    &spec_kz, &spec_phi, &struck_pol, &struck_virtuality})
      v->assign(n, 0.0);
    kp.assign(n * 4, 0.0);
    for (auto* v : {&number, &lam_e, &category_index, &route})
      v->assign(n, 0);
    cell.assign(n, -1);
    n_partner.assign(n, 0);
    struck_pdg.assign(n, 0);
    if (with_rc) {
      // 1.0 is the right fill: it is what every RC weight means when the
      // model does not apply (the coherent channel, a tagged tail).
      rc_lo.assign(n, 1.0);
      rc_hi.assign(n, 1.0);
      rc_tail.assign(n, 1.0);
      rc_clip.assign(n, 0);
    }
  }
};

void fill_row(Columns& c, std::size_t i, const Event& ev, const Pipeline& p,
              bool with_route) {
  c.x[i] = ev.kin.x;
  c.q2[i] = ev.kin.q2;
  c.y[i] = ev.kin.y;
  c.phi[i] = ev.kin.phi;
  c.w2[i] = ev.kin.w2;
  c.nu[i] = ev.kin.nu;
  c.s[i] = ev.kin.s;
  c.weight[i] = ev.weight;
  c.m_ion[i] = ev.spin.m_ion;
  c.m_struck[i] = ev.spin.m_struck;
  c.pe[i] = ev.spin.pe;
  c.theta_s[i] = ev.spin.theta_s;
  c.phi_s[i] = ev.spin.phi_s;
  c.pz[i] = ev.spin.pz;
  c.pzz[i] = ev.spin.pzz;
  c.j[i] = ev.spin.j;
  c.k[i] = ev.kin.k;
  c.cos_theta_k[i] = ev.kin.cos_theta_k;
  c.phi_k[i] = ev.kin.phi_k;
  c.alpha_s[i] = ev.kin.alpha_s;
  c.pt_s[i] = ev.kin.pt_s;
  c.t[i] = ev.kin.t;
  c.x_pom[i] = ev.kin.x_pom;
  c.number[i] = static_cast<std::int64_t>(ev.number);
  c.lam_e[i] = ev.spin.lam_e;
  c.cell[i] = ev.kin.cell;
  // rc.hpp, slot 0.  Untouched (and left at 1.0) when the run had --rc off,
  // in which case these three vectors are empty and never reach the dict.
  if (c.with_rc && ev.rc_weights.size() >= kRcWeightCount) {
    c.rc_lo[i] = ev.rc_weights[0];
    c.rc_hi[i] = ev.rc_weights[1];
    c.rc_tail[i] = ev.rc_weights[2];
    c.rc_clip[i] = static_cast<std::uint8_t>(ev.rc_clipped);
  }

  // The T1 block: the struck nucleon and how many partner spectators came
  // out of the cluster with it (0 for a channel whose struck object already
  // was a nucleon, 1 for d* -> N + N, 1 or 2 for t* -> n + d / p + nn).
  if (const Particle* pn = ev.find(Role::StruckNucleon)) {
    c.struck_pdg[i] = pn->pdg;
    c.struck_pol[i] = pn->pol;
    c.struck_virtuality[i] = pn->p.m2() - M_NUCLEON * M_NUCLEON;
  } else {
    c.struck_virtuality[i] = std::nan("");
    c.struck_pol[i] = 9.0;
  }
  std::int64_t n_part = 0;
  for (const Particle& q : ev.particles)
    if (q.role == Role::PartnerSpectator) ++n_part;
  c.n_partner[i] = n_part;

  if (const Particle* e1 = ev.find(Role::ScatteredElectron)) {
    c.e_prime[i] = e1->p.e;
    c.theta_e[i] = std::atan2(e1->p.pt(), e1->p.pz);
    c.eta_e[i] = pseudorapidity(e1->p);
    c.kp[4 * i + 0] = e1->p.e;
    c.kp[4 * i + 1] = e1->p.px;
    c.kp[4 * i + 2] = e1->p.py;
    c.kp[4 * i + 3] = e1->p.pz;
  }

  const Particle* frag = ev.find(Role::Spectator);
  if (!frag) frag = ev.find(Role::IntactRecoil);
  if (frag) {
    // READ OFF THE RECORD, not re-derived: `Kinematics` now stores
    // `boost_spectator`'s own lab block (and `Pipeline::make_coherent`
    // stores the recoil's), so the four lines of inverse-boost and rigidity
    // algebra this used to repeat have exactly one home.
    c.spec_pT[i] = ev.kin.spec_pt;
    c.spec_theta[i] = ev.kin.spec_theta;
    c.spec_p_lab[i] = ev.kin.spec_p_lab;
    c.spec_phi[i] = ev.kin.phi_spec;
    c.spec_kx[i] = ev.kin.spec_kx;
    c.spec_ky[i] = ev.kin.spec_ky;
    c.spec_kz[i] = ev.kin.spec_kz;
    c.spec_R[i] = ev.kin.spec_r;
    c.spec_xL[i] = ev.kin.spec_xl;
    if (with_route) c.route[i] = route_of(ev, p.optics(), p.pot_config());
  } else {
    c.spec_R[i] = std::nan("");
    c.spec_xL[i] = std::nan("");
    c.route[i] = kRouteLost;
  }
}

/// Generate `n` events (from index 0) into preallocated columns, optionally
/// keeping the `Event` records too.  Called with the GIL released.
void generate_columns(const Pipeline& p, std::uint64_t n, Columns& cols,
                      std::vector<Event>* events, unsigned nthreads,
                      std::size_t chunk) {
  const bool with_route = is_tagged(p.config().channel)
                          || p.config().channel == PipelineChannel::CoherentLi6;
  auto worker = [&](std::uint64_t lo, std::uint64_t hi) {
    if (!events) {
      // Streaming: ONE Event reused for the whole range, i.e.
      // `Pipeline::for_each_range` written out so that `fill_row` can keep
      // its row index without a captured counter.  Materializing a
      // `generate_range` chunk here would buy nothing and cost cache:
      // re-measured 2026-08-30, after `Pipeline::event` learned to
      // reconstruct in place, 2.19 M ev/s through a 4096-record buffer
      // against 2.87 M streaming (inclusive 6Li, mid configuration).
      Event ev;
      for (std::uint64_t i = lo; i < hi; ++i) {
        p.event(i, ev);
        fill_row(cols, static_cast<std::size_t>(i), ev, p, with_route);
      }
      return;
    }
    std::vector<Event> buf;
    for (std::uint64_t i = lo; i < hi; i += chunk) {
      const std::uint64_t j = std::min<std::uint64_t>(i + chunk, hi);
      p.generate_range(i, j, buf);
      for (std::size_t r = 0; r < buf.size(); ++r) {
        const std::size_t row = static_cast<std::size_t>(i) + r;
        fill_row(cols, row, buf[r], p, with_route);
        (*events)[row] = std::move(buf[r]);
      }
    }
  };
  if (nthreads <= 1 || n < 4 * chunk) {
    worker(0, n);
    return;
  }
  std::vector<std::thread> pool;
  const std::uint64_t block = (n + nthreads - 1) / nthreads;
  for (unsigned t = 0; t < nthreads; ++t) {
    const std::uint64_t lo = std::min<std::uint64_t>(t * block, n);
    const std::uint64_t hi = std::min<std::uint64_t>(lo + block, n);
    if (lo < hi) pool.emplace_back(worker, lo, hi);
  }
  for (auto& th : pool) th.join();
}

/// The selected unpolarised backend's own grid floor, and how much of THIS
/// run's own accepted RATE sits below it.
///
/// ONE DEFINITION, read by both `meta["unpol_sf_grid_q2_min"]` /
/// `meta["unpol_sf_below_grid_frac"]` and the CLI run banner, so the fraction
/// is never a typed-in number (docs/CONVENTIONS.md).  The floor comes from
/// the loaded grid's own metadata (`LhapdfSF::q2_min`, LHAPDF's `q2Min()`).
///
/// THE DENOMINATOR IS THE CHANNEL'S OWN (`Pipeline::cell_rate_weights_pb`),
/// not the inclusive sampler's cells.  On the inclusive and the tagged
/// channels the two are the same vector; on the COHERENT one they are not,
/// because the coherent yield is `sigma_cell * f_coh(x)` over the cells that
/// admit a diffractive mass, and `f_coh` falls by a factor ~26 across the
/// window.  Measured on the shipped 6Li ct18nlo run (2026-09-05): the
/// coherent rate has 0.447460 below CT18NLO's floor where the inclusive
/// cells it rides on have 0.361791 -- and the coherent weights sum to
/// 8806.096207 pb, which is `Pipeline::sigma_pb()` exactly.  Reporting the
/// inclusive number on a coherent run was a quantified `meta` key that did
/// not describe the run it was attached to.
///
/// Three cases, and NaN is used where the honest answer is "not known":
///   * a caller-supplied `cfg.kernel` -- the backend is the caller's and this
///     function cannot see its grid: (NaN, NaN);
///   * `unpol_sf = Toy` -- `ToyF2` is closed form and has no grid at all, so
///     nothing can be below one: (0, 0);
///   * an `LhapdfSF`: its own `q2Min()` and the measured fraction;
///   * anything else (`MstwSF`, a `Custom` object) -- PYTHIA's `MSTWpdf`
///     keeps its `qsqmin` private and a custom backend reports nothing, so
///     the floor is NOT known and must not be guessed: (NaN, NaN).
std::pair<double, double> unpol_sf_grid_report(const Pipeline& p) {
  const double kNan = std::numeric_limits<double>::quiet_NaN();
  const PipelineConfig& cfg = p.config();
  if (cfg.kernel) return {kNan, kNan};
  if (cfg.unpol_sf == UnpolSfSource::Toy) return {0.0, 0.0};
#if LIPOLGEN_HAVE_LHAPDF
  if (const auto lh =
          std::dynamic_pointer_cast<const LhapdfSF>(cfg.unpol_sf_obj)) {
    const double q2_min = lh->q2_min();
    const std::vector<double>& sc = p.cell_rate_weights_pb();
    const std::vector<double>& q2c = p.dis_sampler().q2_cells();
    double tot = 0.0, below = 0.0;
    for (std::size_t i = 0; i < sc.size(); ++i) {
      tot += sc[i];
      if (q2c[i] < q2_min) below += sc[i];
    }
    return {q2_min, tot > 0.0 ? below / tot : kNan};
  }
#endif
  return {kNan, kNan};
}

#if LIPOLGEN_HAVE_PYTHIA8
/// `PipelineConfig::hadronizer` when the T2 tier was bound by
/// `set_pythia_hadronizer` -- a NAMED callable rather than a lambda, so
/// `columns_to_dict` can recover the bridge with `std::function::target` and
/// write the Pomeron block of the `meta`.  Without it the bridge is invisible
/// to the metadata writer (the hook is type-erased), which is why two runs
/// differing only in `--pom-set` used to produce identical `meta` while their
/// whole hadronic final state differed (D4.7).
struct PythiaHadronizerHook {
  std::shared_ptr<PythiaBridge> bridge;
  void operator()(Event& ev, Rng& rng) const { bridge->hadronize(ev, rng); }
};
#endif

py::dict columns_to_dict(Columns& c, const Pipeline& p, std::uint64_t n,
                         KnobRunContext ctx = KnobRunContext()) {
  py::dict d;
  d["x"] = move_array(std::move(c.x));
  d["q2"] = move_array(std::move(c.q2));
  d["y"] = move_array(std::move(c.y));
  d["phi"] = move_array(std::move(c.phi));
  d["w2"] = move_array(std::move(c.w2));
  d["nu"] = move_array(std::move(c.nu));
  d["s"] = move_array(std::move(c.s));
  d["weight"] = move_array(std::move(c.weight));
  d["m_ion"] = move_array(std::vector<double>(c.m_ion));
  d["m"] = move_array(std::move(c.m_ion));   // polligen's inclusive key
  d["m_struck"] = move_array(std::move(c.m_struck));
  d["pe"] = move_array(std::move(c.pe));
  d["theta_s"] = move_array(std::move(c.theta_s));
  d["phi_s"] = move_array(std::move(c.phi_s));
  d["pz"] = move_array(std::move(c.pz));
  d["pzz"] = move_array(std::move(c.pzz));
  d["j"] = move_array(std::move(c.j));
  d["k"] = move_array(std::move(c.k));
  d["cos_theta_k"] = move_array(std::move(c.cos_theta_k));
  d["phi_k"] = move_array(std::move(c.phi_k));
  d["alpha_s"] = move_array(std::move(c.alpha_s));
  d["pt_s"] = move_array(std::move(c.pt_s));
  d["t"] = move_array(std::move(c.t));
  d["x_pom"] = move_array(std::move(c.x_pom));
  d["e_prime"] = move_array(std::move(c.e_prime));
  d["theta_e"] = move_array(std::move(c.theta_e));
  d["eta_e"] = move_array(std::move(c.eta_e));
  d["kp"] = move_array2(std::move(c.kp), 4);
  d["pT"] = move_array(std::move(c.spec_pT));
  d["theta"] = move_array(std::move(c.spec_theta));
  d["p_lab"] = move_array(std::move(c.spec_p_lab));
  d["R"] = move_array(std::move(c.spec_R));
  d["xL"] = move_array(std::move(c.spec_xL));
  d["kx"] = move_array(std::move(c.spec_kx));
  d["ky"] = move_array(std::move(c.spec_ky));
  d["kz"] = move_array(std::move(c.spec_kz));
  d["phi_spec"] = move_array(std::move(c.spec_phi));
  d["number"] = move_array(std::move(c.number));
  d["lam_e"] = move_array(std::move(c.lam_e));
  d["route"] = move_array(std::move(c.route));
  // polligen's Mode-W reweighting column, and the T1 block.
  d["cell"] = move_array(std::move(c.cell));
  d["n_partner"] = move_array(std::move(c.n_partner));
  d["struck_pdg"] = move_array(std::move(c.struck_pdg));
  d["struck_pol"] = move_array(std::move(c.struck_pol));
  d["struck_virtuality"] = move_array(std::move(c.struck_virtuality));
  // rc.hpp: the three RC columns exist ONLY when the run had RC on, so an
  // `--rc off` npz carries exactly today's key set and the reference gates do
  // not have to move.  `export.columns_from_events` applies the same rule on
  // the records path -- both paths or neither, or `export.rc_columns` would
  // silently return ones for a run that had RC on.
  if (c.with_rc) {
    d["rc_tensor_lo"] = move_array(std::move(c.rc_lo));
    d["rc_tensor_hi"] = move_array(std::move(c.rc_hi));
    d["rc_tail"] = move_array(std::move(c.rc_tail));
  }

  py::list names;
  for (const SpinCategory& cat : p.plan().categories()) names.append(cat.name);
  d["category_names"] = names;
  py::array_t<std::int64_t> idx = move_array(std::move(c.category_index));
  d["category_index"] = idx;
  // Per-event category name, so a consumer sees polligen's "category" key.
  py::object np = py::module_::import("numpy");
  d["category"] = np.attr("asarray")(names).attr("__getitem__")(idx);

  // THE KNOB-PROVENANCE TABLE, built once and read by every key below that
  // has a reach rule (`Pipeline::knob_provenance`, pipeline.hpp).  The
  // context's T2 half is filled from the bridge ACTUALLY bound to this
  // config, not from what a caller says: `set_pythia_hadronizer` wraps the
  // bridge in a NAMED callable exactly so the metadata writer can recover it
  // (`PythiaHadronizerHook`), and taking `--pom-set` from the caller instead
  // would be one more place for the file and the run to disagree.  The
  // run-plan half is the caller's, because a `RunPlan` records its moments
  // and not which flags produced them.
#if LIPOLGEN_HAVE_PYTHIA8
  if (const auto* hook0 = p.config().hadronizer.target<PythiaHadronizerHook>();
      hook0 != nullptr && hook0->bridge) {
    const PythiaBridgeOptions& po0 = hook0->bridge->options();
    ctx.t2_bound = true;
    ctx.t2_pomeron = (po0.coherent_t2 == CoherentT2::Pomeron);
    ctx.pom_set = po0.pom_set;
    ctx.pom_rescale = po0.pom_rescale;
  } else {
    ctx.t2_bound = false;
  }
#else
  ctx.t2_bound = false;
#endif
  const std::vector<KnobProvenance> knobs = p.knob_provenance(ctx);
  auto knob = [&knobs](const char* name) -> const KnobProvenance& {
    for (const KnobProvenance& r : knobs)
      if (r.name == name) return r;
    throw std::runtime_error(std::string("columns_to_dict: no knob_provenance "
                                         "row named ") + name);
  };

  py::dict meta;
  meta["generator"] = "LiPolGen";
  meta["version"] = LIPOLGEN_VERSION;
  meta["channel"] = pipeline_channel_name(p.config().channel);
  meta["isotope"] = p.config().isotope;
  meta["beam_config"] = p.config().beam_config;
  meta["electron_energy"] = p.beam_config().electron_energy;
  meta["p_per_nucleon"] = p.beam_config().ion_momentum_per_nucleon;
  meta["s_per_nucleon"] = p.beam_config().s_per_nucleon();
  meta["seed"] = p.config().seed;
  meta["run"] = p.config().run;
  // THE PLAN'S NAME.  `--plan` was the one physics-relevant CLI flag with no
  // provenance of any kind until 2026-09-16: `RunPlan` records its MOMENTS
  // and not which pattern produced them, so the plan was recoverable from an
  // npz only by reading `category_names` and the pz/pzz/pe/lam_e columns back
  // out.  Written only when the caller supplied one (`cli.main` and `lg.run`
  // both do), so a `Pipeline.generate()` with no context is unchanged.
  if (!ctx.plan_name.empty()) meta["plan"] = ctx.plan_name;
  meta["n_events"] = n;
  meta["n_events_total"] = p.size();
  meta["sigma_pb"] = p.sigma_pb();
  meta["sigma_gen_mb"] = p.sigma_pb() * 1e-9;   // 1 mb = 1e9 pb
  meta["tier"] = (p.tier() == Tier::T1) ? "T1" : "T0";
  meta["sigma_per_category_pb"] = p.sigma_per_category_pb();
  meta["lumi_per_category_pb"] = p.lumi_per_category_pb();
  // THE TWO ROUTE KEYS GO THROUGH THE TABLE, like `pol_sf` and the three
  // Pomeron keys below and for the identical reason.  They were written
  // unconditionally until 2026-09-05, so on the INCLUSIVE channel -- where
  // the table says not-read and the banner's provenance block labels them --
  // `--optics yr-high-divergence` recorded `meta["optics"] =
  // "10x100 high-divergence"` and a `pot_config` of "18x275" as bare values
  // while all 47 columns were bit-identical (16 cells over the five inclusive
  // specs of the 2026-09-05 matrix).  That is the `pol_sf` defect on two
  // top-level keys, and it survived the fix that was supposed to close the
  // class because these two bypassed the mechanism entirely.
  //
  // AND THE NAME IS NOT LOST BY LABELLING THEM -- but it WAS, between
  // 2026-09-05 and the phase-F fix of the same day.  The first form of this
  // block wrote the bare scope clause "not read on this run" into both keys
  // on every channel that HAS a route (tagged-d-p, tagged-6Li-alpha, the
  // coherent channel), where the classifier IS consulted and only the
  // SENSITIVITY to the envelope is missing; the a94fd6e files carried
  // "10x100 high-acceptance" / "10x100" there, and the CLI banner printed
  // "not read on this run" three lines above a tag fraction quoted AT that
  // envelope.  The scalar said more than the provenance row it came from.
  //
  // The fix is in the ONE place the label is built (`Pipeline::knob_provenance`,
  // src/core/pipeline.cpp): on a channel with a route the label now opens
  // "not read on this run at 10x100 high-acceptance (...)", so the name is in
  // the scalar, in the banner and in the provenance row, and the scope clause
  // is still there.  Nothing is special-cased HERE: these two keys go through
  // `meta_value()` exactly like `pol_sf` and the three Pomeron keys, which is
  // the whole point of the mechanism (CONVENTIONS.md "one mechanism").
  // The bare name also stays available in
  // `meta["knob_provenance"]["optics"]["value"]`, which is where a consumer
  // that wants to re-route the sample reads it.
  meta["optics"] = knob("optics").meta_value();
  meta["pot_config"] = knob("pot_config").meta_value();
  meta["frame"] = "head-on: ion +z, electron -z; per-nucleon x, y, Q2";
  // Which b1 backend and which band row made this file (design_D_b1_li6.md
  // sec. 4.2).  UNCONDITIONAL, unlike the `rc` block below: the band is
  // produced by re-running with `--b1-band-scale 0/1/2`, and without these
  // three keys three otherwise identical npz files are indistinguishable --
  // which is exactly the "never quote a single row" rule failing silently.
  // On a caller-supplied `cfg.kernel` the flag is not what produced the
  // numbers, so say so rather than printing a model name that did not run.
  //
  // THE 7Li ZERO IS RECORDED, NOT LABELLED WITH A BACKEND THAT DID NOT RUN
  // (`inclusive_rank2_is_empty`, pipeline.hpp).  `default_inclusive_kernel`
  // fills a rank-2 slot for spin 1 only, so on `--isotope 7Li --channel
  // inclusive` nothing filled one and `"miller"` here named a backend the run
  // never used -- the same defect `validate()` refuses for `b1_band_scale`,
  // applied to the model name.  The two SCALES stay numeric and stay at 1.0:
  // `validate()`'s Miller branch already refuses any other value on that path,
  // so neither can record a variation that did not run, and changing their
  // TYPE per isotope would break every consumer that reads them as floats.
  // THROUGH THE TABLE, like every other key in this block (2026-09-05).  The
  // same three-way ternary -- caller-supplied kernel / rank2_none_label() /
  // the backend name -- used to be typed here AND in the `b1_model` row of
  // `Pipeline::knob_provenance`, which is the one thing CONVENTIONS.md's
  // "the three surfaces read that table and nothing else" forbids: two
  // spellings of one label drift.  The row's status is Read or Refused, never
  // NotRead, so `meta_value()` is its `value` and this is bit-identical.
  meta["b1_model"] = knob("b1_model").meta_value();
  meta["b1_band_scale"] = p.config().b1_band_scale;
  meta["b1_alpha_d_dwave_weight"] = p.config().b1_alpha_d_dwave_weight;
  // Which unpolarised PDF the b1 convolution folded against.  Same
  // provenance rule and the same unconditional write as the three keys above:
  // "toy" and "mstw" differ by up to a factor 1.85 on Li6ConvolutionB1::b1
  // and by nothing at all in any other npz key, so without this one they are
  // indistinguishable files.
  meta["b1_unpol"] = knob("b1_unpol").meta_value();
  // Which R = sigma_L/sigma_T the b1 convolution's numerator AND the kernel
  // that carries its denominator were BOTH built from (`RSource`; CLI
  // --r-source).  Same provenance rule and the same unconditional write:
  // "unset" and "r1998" differ by -0.6847 % in sigma_pb and by up to
  // +26.4 % in the tensor weight at the standard points (y = 0.5, Q2 = 2.5,
  // x = 0.10; docs/open_items/run_2026-09-06/phase_A_numbers.md sec. A1), so
  // without this key the two are indistinguishable files.  The row is Read
  // or Refused at every value but `sigma-lt`, where it is NotRead and this
  // key carries the label -- which is honest: that run IS bit-identical to
  // the unset one.
  meta["r_source"] = knob("r_source").meta_value();
  // ... and the SENTENCE, on every run and every channel, from the one
  // definition the CLI run banner also prints (`rank2_input_report`).  This
  // is the key that makes a 7Li tensor run self-describing: before it, an npz
  // whose whole tensor sector was exactly 0.0 was indistinguishable in `meta`
  // from one whose b1 merely happened to be small.
  meta["rank2_input"] = rank2_input_report(p.config(), p.plan());
  // Which structure-function backend EVERY kernel of this run actually used
  // (`UnpolSfSource`, `PolSfSource`; CLI --unpol-sf / --pol-sf).  Same
  // provenance rule and the same UNCONDITIONAL write as the b1 keys above,
  // and for a bigger reason: `toy` and `ct18nlo` differ by a FACTOR 0.80 on
  // the 6Li rate and move every column of the file, so without these two
  // keys two npz files that disagree by 20 % in cross section would be
  // indistinguishable in `meta`.  On a caller-supplied `cfg.kernel` the
  // selectors are not what produced the numbers (and on a tagged channel the
  // kernel itself is not even read), so say that rather than print a backend
  // name that did not run.
  meta["unpol_sf"] = p.config().kernel
                         ? std::string("caller-supplied kernel")
                         : std::string(unpol_sf_name(p.config().unpol_sf));
  // ... and the POLARISED one, which has one channel the unpolarised one
  // does not: `CoherentLi6` never evaluates g1 at all (`pol_sf_is_read`), so
  // a coherent run records the "not read on channel ..." LABEL and not a
  // backend name -- the `rank2_none_label` rule.  "toy" would be as wrong as
  // "nnpdfpol" there: `ToyG1` did not run on that channel either.  Measured
  // 2026-09-05: every generated column of a coherent run is bit-identical
  // between the two settings, while `meta["pol_sf"]` used to print whichever
  // one was typed.
  //
  // AND IT IS TWO AXES NOW, NOT ONE.  The label went in on the CHANNEL axis
  // on 2026-09-04 and the identical defect on the RUN-PLAN axis survived it:
  // under `tensor-thirds`, `transverse-tensor`, `tensor-flip` or
  // `helicity-flip --pe 0` the fill carries lam_e * P_e = 0 in every category
  // and g1 is multiplied by zero, so `--pol-sf nnpdfpol` is bit-identical to
  // `toy` -- on the CLI's own DEFAULT plan -- while this key wrote
  // "nnpdfpol".  Both axes now come from the ONE `knob_provenance` row, which
  // is `pol_sf_is_read(config, plan)`.
  meta["pol_sf"] = knob("pol_sf").meta_value();
  // ... with the reason beside it, on every run and every channel, from the
  // one definition the CLI run banner also prints (`pol_sf_reach_report`) --
  // the `rank2_input` arrangement.  This is the key that makes a coherent
  // run self-describing: the label above says the selector did not run, this
  // one says why, and says that `--unpol-sf` DID reach the same run.
  meta["pol_sf_reach"] = knob("pol_sf").reason;
  // ... and how much of THIS run's accepted rate sits below the selected
  // grid's own floor.  `--unpol-sf ct18nlo` on the shipped 6Li window puts
  // 36 % of its own cell cross section below CT18NLO's Q2 = 1.677, where
  // LHAPDF keeps evolving downward instead of freezing; that is a property
  // of the run, not of the flag, so it is measured per run rather than
  // quoted.  NaN means "the backend reports no grid floor" (MstwSF, a Custom
  // object, or a caller-supplied kernel) -- never a silent 0.
  const std::pair<double, double> grid = unpol_sf_grid_report(p);
  meta["unpol_sf_grid_q2_min"] = grid.first;
  meta["unpol_sf_below_grid_frac"] = grid.second;
  // fsi.hpp -- CONDITIONAL, following the `rc` precedent below rather than
  // the unconditional `b1` one, so an `--fsi off` npz carries exactly
  // today's key set and no reference gate moves.  Measured on the three runs
  // that made this block necessary (D3.3, --channel tagged-6Li-alpha
  // --events 20000 --seed 1234, tensor-thirds at pz 0.7 / pzz 0.6 / pe 0.7,
  // at the default sigma_XN = 40 mb): the kinematic columns are BIT-IDENTICAL
  // across `off` / `glauber-cluster` / `glauber-nucleon`, the total rate
  // differs by 48 % / 42 %, individual event weights by up to a factor 68.5
  // -- and the `meta` dicts compared EQUAL, key by key.  That is exactly the
  // "never quote a single row" rule failing silently, one level worse than
  // the b1 case it was written for, because `fsi_sigma_mb` carries a
  // documented 20-40 mb band that must never be quoted as a single row.
  if (const GlauberFsiWeight* fw = p.fsi_weight()) {
    meta["fsi"] = std::string(pipeline_fsi_name(p.config().fsi));
    meta["fsi_sigma_mb"] = p.config().fsi_sigma_mb;
    // What the weight actually applied, not only what was asked for: the
    // cluster profile and the integrated survival ARE the difference between
    // the two variants (131.0 vs 160.0 mb, 0.520 vs 0.583 on the 6Li alpha
    // tag), so a file that records the name records the numbers too.
    const double sig = fw->sigma_eff_mb(0.0);
    meta["fsi_sigma_cluster_mb"] = fw->sigma_cluster_mb(sig);
    meta["fsi_sigma_cluster_el_mb"] = fw->sigma_cluster_el_mb(sig);
    meta["fsi_survival"] = fw->survival();
    meta["fsi_clipped_grid_fraction"] = fw->clipped_grid_fraction();
    // The sigma_XN(W) ramp is a DIFFERENT run: `fsi_sigma_mb` alone does not
    // describe it, so the anchors go in only when the ramp ran.
    meta["fsi_formation_ramp"] = fw->options().formation_ramp;
    if (fw->options().formation_ramp) {
      meta["fsi_ramp_w_lo"] = fw->options().ramp_w_lo;
      meta["fsi_ramp_sigma_lo_mb"] = fw->options().ramp_sigma_lo_mb;
      meta["fsi_ramp_w_hi"] = fw->options().ramp_w_hi;
      meta["fsi_ramp_sigma_hi_mb"] = fw->options().ramp_sigma_hi_mb;
    }
  }
  // coherent.hpp -- CONDITIONAL on the channel, same precedent.  `t_max` is
  // the one that had to be here: it moves the whole |t| spectrum, the tag
  // acceptance and every c_2 in the file, and until 2026-09-04 it was
  // neither reachable from the CLI nor recorded anywhere (D5.6).  The
  // positivity edge is DERIVED (`t_positivity_edge`), so the file says how
  // far the run was from the edge rather than leaving a reader to recompute
  // it from eps_b0 -- and it is quoted at P_zz = -2, the worst case of the
  // shipped scenario.
  if (p.config().channel == PipelineChannel::CoherentLi6) {
    const CoherentScenario& sc = p.config().coherent;
    meta["coherent_t_max"] = p.config().coherent_t_max;
    meta["coherent_slope_b"] = sc.slope_b;
    meta["coherent_eps_b0"] = sc.eps_b0;
    meta["coherent_amp"] = sc.amp;
    meta["coherent_f0"] = sc.f0;
    meta["coherent_t_positivity_edge_pzz_m2"] = sc.t_positivity_edge(-2.0);
    meta["coherent_m_x_min"] = p.config().coherent_xpom.m_x_min;
    meta["coherent_x_pom_max"] = p.config().coherent_xpom.x_pom_max;
    meta["coherent_weighted_azimuth"] = p.config().coherent_weighted_azimuth;
  }
#if LIPOLGEN_HAVE_PYTHIA8
  // pythia_bridge.hpp -- CONDITIONAL on the T2 tier having been bound
  // through `set_pythia_hadronizer`, which is the only route by which the
  // metadata writer can see the bridge at all.  Measured (D4; 20 000 coherent
  // events per set at 6Li config 1, seed 4242): the `--pom-set` runs of the
  // same seed produce BIT-IDENTICAL T0 columns (t, x_pom, q2, x, weight --
  // one md5, ffd35a3a62b591c547e9ca2ac4301b5d, over the FOURTEEN sets the
  // constructor still admits, 1-10 and 12-15; set 11 gave that md5 too in
  // the 2026-09-04 scan but is refused now, so fourteen is the count that
  // reproduces, re-measured 2026-09-05) and a hadronic final state that
  // moves by
  // -2.6 % / +9.9 % in <n_charged> and by a factor 2.8 in the kaon fraction,
  // while their `meta` compared equal.  THAT BAND IS OVER THE TWELVE DPDF
  // FITS (3-10, 12-15) ABOUT SET 6, not over all fifteen: sets 1 (a toy) and
  // 2 (pi^0 densities) are not Pomeron fits and 11 is refused, and set 2
  // would take the <n_charged> low edge to -4.7 % (3.7773 against 3.9639).
  // `pom_set`'s own declaration in `pythia_bridge.hpp` is where the band and
  // its window are stated once.  The
  // fallback COUNTER is here for the reason PipelineConfig::validate refuses
  // a knob that did not run: on `pom_set = 11` it reads 1.000, i.e. the
  // Pomeron PDF was not consulted on a single event.
  if (const auto* hook = p.config().hadronizer
                             .target<PythiaHadronizerHook>();
      hook != nullptr && hook->bridge) {
    const PythiaBridge& br = *hook->bridge;
    const PythiaBridgeOptions& po = br.options();
    meta["t2_bridge"] = std::string("pythia8");
    // THE THREE POMERON KNOBS GO THROUGH THE TABLE, and the reason is the
    // measurement above turned on its head.  The Pomeron PYTHIA instance is
    // built whenever `coherent_t2 == Pomeron` REGARDLESS of channel
    // (pythia_bridge.cpp) and then hadronizes ZERO events off the coherent
    // one, so on the six non-coherent channels these three moved NOTHING --
    // measured 2026-09-05 with a deterministic hadron hash, 150 events seed
    // 4242: not a T0 column and not one of 3278 particles -- while this block
    // wrote `pom_set = 5` and `coherent_t2 = "pomeron"` as if they had run.
    // `pom_set` stays an INT where it ran and carries the label where it did
    // not, the `b1_model` / `pol_sf` rule; `meta["knob_provenance"]` carries
    // the whole sentence beside it.
    const KnobProvenance& kt2 = knob("coherent_t2");
    const KnobProvenance& kps = knob("pom_set");
    const KnobProvenance& kpr = knob("pom_rescale");
    meta["coherent_t2"] = kt2.meta_value();
    if (kps.status == KnobStatus::Read) meta["pom_set"] = po.pom_set;
    else meta["pom_set"] = kps.meta_value();
    if (kpr.status == KnobStatus::Read) meta["pom_rescale"] = po.pom_rescale;
    else meta["pom_rescale"] = kpr.meta_value();
    meta["pom_q2_pdf_min"] = po.q2_pdf_min;
    meta["t2_include_charm"] = po.include_charm;
    const PythiaBridgeStats& st = br.stats();
    meta["t2_n_ok"] = st.n_ok;
    meta["t2_n_failed"] = st.n_failed;
    meta["t2_n_pomeron"] = st.n_pomeron;
    meta["n_pom_flavour_fallback"] = st.n_pom_flavour_fallback;
    meta["pom_flavour_fallback_frac"] =
        st.n_pomeron > 0
            ? static_cast<double>(st.n_pom_flavour_fallback)
                  / static_cast<double>(st.n_pomeron)
            : 0.0;
  }
#endif
  // rc.hpp -- the WHOLE block, `meta["rc"]` included, is emitted only when
  // the run has RC on, so an `--rc off` npz is byte-identical to today's and
  // not merely key-compatible with it (T6).
  if (const RcModel* rc = p.rc_model()) {
    meta["rc"] = rc_mode_name(p.config().rc);
    py::list rc_names;
    for (std::size_t i = 0; i < kRcWeightCount; ++i)
      rc_names.append(rc_weight_name(i, 0));
    meta["rc_weight_names"] = rc_names;
    meta["rc_applies"] = rc->applies();
    meta["rc_tail_applies"] = rc->tail_applies();
    meta["rc_exclusion_reason"] = rc->exclusion_reason();
    meta["rc_ff_provenance"] = rc->ff_provenance();
    meta["rc_delta_low_x"] = rc->options().delta_low_x;
    meta["rc_delta_high_x"] = rc->options().delta_high_x;
    // The A = 2 -> A = 6 transfer price (design Q8).  Recorded for the same
    // reason as `rc_qe_tensor_scale` below: it widens EVERY rc_tensor_* column
    // by sqrt(1 + f^2) without changing any other key, so without it a
    // 1.41x-wider band is indistinguishable in `meta` from the shipped one.
    meta["rc_a_transfer_frac"] = rc->options().a_transfer_frac;
    meta["rc_fq_scale"] = rc->options().fq_scale;
    // The C0 SHAPE edge.  Same provenance rule as `b1_unpol` above: "ho" and
    // "vmc-ft" move F_c and F_q TOGETHER and flip the sign of
    // sigma^el_T/sigma^el_U at x = 0.1 (phase_B_numbers.md sec. B1), so
    // without this key the two npz files are indistinguishable in `meta`.
    meta["rc_c0_shape"] =
        std::string(c0_shape_name(rc->options().c0_shape));
    meta["rc_tail_tensor_scale"] = rc->options().tail_tensor_scale;
    meta["rc_qe_suppression"] = rc->options().qe_suppression;
    // The POLARISED quasi-elastic stand-in.  Recorded for the same reason as
    // every other knob here: it is 0 by default (the tail is tensor-blind on
    // its dominant piece), and a run that turned it on must be
    // distinguishable in `meta` from one that did not -- it is a BORROWED
    // magnitude, not a computed tail, and an analysis has to be able to see
    // that it was applied (rc.hpp, RcOptions::qe_tensor_scale).
    meta["rc_qe_tensor_scale"] = rc->options().qe_tensor_scale;
    // The s-/p-PEAK stand-in, recorded for the same reason and under the same
    // warning: 0 by default (the leading-log peaks are tensor-blind, and
    // under the shipped `t-peak` tail they do not exist at all), and a run
    // that turned it on must be distinguishable in `meta` from one that did
    // not -- it is a BOUND WITH NO DERIVATION, not a computed tensor s/p peak
    // (rc.hpp, RcOptions::sp_tensor_scale).
    meta["rc_sp_tensor_scale"] = rc->options().sp_tensor_scale;
    meta["rc_qe_kf_gev"] = rc->options().qe_kf_gev;
    meta["rc_band_tau_max"] = rc->options().band_tau_max;
    // The NON-numeric knobs, recorded for the same provenance reason as
    // `b1_band_scale` above: `scope = TensorAll` and `with_qe_tail = false`
    // are API-reachable and change EVERY rc_tensor_*/rc_tail column (the QRT
    // is 73-99.9 % of the tail at x >= 0.1, phase_C_numbers.md), so without
    // these keys such an npz is indistinguishable in `meta` from a default
    // run.  `rc_tail_applies` above covers `with_tail` only.
    // ... and `rc_scope` goes THROUGH THE TABLE, because it is the one rc
    // knob whose reach depends on the RUN PLAN and so cannot be refused by
    // `validate()`, which has no plan.  `tensor-all` differs from
    // `tensor-rate` only in the cos 2phi amplitude, which carries a
    // sin^2(theta_S): under `tensor_thirds_plan`, whose categories sit at
    // theta_S = 0, the two runs are BIT-IDENTICAL (measured 2026-09-05) and
    // "tensor-all" here would name a variation that did not happen.  Under
    // `transverse-tensor` / `tensor-flip` (theta_S = pi/2) it does run, and
    // the name goes in.
    meta["rc_scope"] = knob("rc_scope").meta_value();
    meta["rc_with_qe_tail"] = rc->options().with_qe_tail;
    meta["rc_tail_model"] =
        std::string(rc_tail_model_name(rc->options().tail_model));
    meta["rc_n_eta"] = rc->options().n_eta;
    meta["rc_tail_max"] = rc->options().tail_max;
    meta["rc_clipped_cell_fraction"] = rc->clipped_cell_fraction();
    const std::array<double, 3> by_y = rc->clipped_fraction_by_y();
    meta["rc_clipped_fraction_by_y"] =
        std::vector<double>(by_y.begin(), by_y.end());
    // ... and the EVENT-level clipping, which is a DIFFERENT quantity: the
    // node fractions above are not event-weighted, and the node statistic is
    // computed at q_n = 0 while the per-event clip carries the (q_n/6) tensor
    // term.  In the default 2000-event inclusive run the tail ceiling IS hit
    // by real events while the global node fraction reads 1.7 %.
    std::size_t n_band = 0, n_tail = 0;
    for (std::size_t i = 0; i < c.rc_clip.size() && i < n; ++i) {
      if (c.rc_clip[i] & kRcClipBand) ++n_band;
      if (c.rc_clip[i] & kRcClipTail) ++n_tail;
    }
    const double denom = (n > 0) ? static_cast<double>(n) : 1.0;
    meta["rc_clipped_band_events"] = n_band;
    meta["rc_clipped_tail_events"] = n_tail;
    meta["rc_clipped_band_event_fraction"] =
        static_cast<double>(n_band) / denom;
    meta["rc_clipped_tail_event_fraction"] =
        static_cast<double>(n_tail) / denom;
  }
  // ... AND THE WHOLE TABLE, on every run and every channel, unconditionally.
  // This is the block that closes the fourth defect class: `--cluster-wave`,
  // `--triton-sf`, `--inclusive-b1`, `--cluster-beta`, `--p-d`,
  // `--fsi-sigma-mb` at `--fsi off`, `--coherent-t-max` off the coherent
  // channel and `--x-max` on it were each accepted, some of them MOVED the
  // output, and NONE of them was recorded anywhere -- a VMC and a Hulthen
  // tagged file were indistinguishable in `meta`.  Every knob is here now,
  // with its status and the reason, so silence is no longer reachable: a knob
  // added without a `knob_provenance` row fails
  // python/tests/test_knob_provenance.py.
  //
  // It is a NESTED dict and `export.write_columns_npz` writes `meta` as JSON,
  // so it survives the round trip unchanged.  The existing keys keep their
  // names and their meaning; this one is added beside them.
  py::dict prov;
  for (const KnobProvenance& r : knobs) {
    py::dict row;
    row["value"] = r.value;
    row["status"] = std::string(knob_status_name(r.status));
    row["reason"] = r.reason;
    row["at_default"] = r.at_default;
    if (!r.flag.empty()) row["flag"] = r.flag;
    if (!r.label.empty()) row["label"] = r.label;
    prov[py::str(r.name)] = row;
  }
  meta["knob_provenance"] = prov;

  d["meta"] = meta;
  return d;
}

py::object pipeline_generate(const Pipeline& p, std::uint64_t n, bool events,
                             unsigned nthreads, std::size_t chunk,
                             const KnobRunContext& ctx = KnobRunContext()) {
  if (n == 0 || n > p.size()) n = p.size();
  Columns cols;
  cols.with_rc = p.rc_model() != nullptr;
  cols.resize(static_cast<std::size_t>(n));
  // category index is cheap and needs no event: recover it from the index
  for (std::uint64_t i = 0; i < n; ++i)
    cols.category_index[static_cast<std::size_t>(i)] =
        static_cast<std::int64_t>(p.category_of(i));
  std::vector<Event> evs;
  if (events) evs.resize(static_cast<std::size_t>(n));
  {
    py::gil_scoped_release unlock;
    generate_columns(p, n, cols, events ? &evs : nullptr, nthreads, chunk);
  }
  py::dict d = columns_to_dict(cols, p, n, ctx);
  if (events) {
    py::list lst;
    for (auto& e : evs) lst.append(py::cast(std::move(e)));
    d["events"] = lst;
  }
  return d;
}

// ------------------------------------------------------------- HFS exporter

/// Flat particle arrays (polligen `HFSSample` layout) from a sequence of
/// events.  `roles` selects which particles enter the list; the scattered
/// electron is never in it (it is the `kp` column), exactly as
/// tools/pythia8/gen_dis_hfs.py writes its samples.
py::dict hfs_arrays(const std::vector<Event>& evs, bool include_spectators,
                    bool include_hadronic_x) {
  const std::size_t n = evs.size();
  std::vector<std::int64_t> offsets(n + 1, 0), pid;
  std::vector<double> charge, p4;
  std::vector<double> x(n, 0.0), q2(n, 0.0), y(n, 0.0), weight(n, 1.0);
  std::vector<double> kp(4 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const Event& ev = evs[i];
    x[i] = ev.kin.x;
    q2[i] = ev.kin.q2;
    y[i] = ev.kin.y;
    weight[i] = ev.weight;
    if (const Particle* e1 = ev.find(Role::ScatteredElectron)) {
      kp[4 * i + 0] = e1->p.e;
      kp[4 * i + 1] = e1->p.px;
      kp[4 * i + 2] = e1->p.py;
      kp[4 * i + 3] = e1->p.pz;
    }
    for (const Particle& pa : ev.particles) {
      if (pa.status != Status::Final) continue;
      const bool keep =
          pa.role == Role::Hadron
          || (include_spectators
              && (pa.role == Role::Spectator || pa.role == Role::PartnerSpectator
                  || pa.role == Role::IntactRecoil))
          || (include_hadronic_x && pa.role == Role::HadronicX);
      if (!keep) continue;
      pid.push_back(pa.pdg);
      charge.push_back(pa.charge);
      p4.push_back(pa.p.e);
      p4.push_back(pa.p.px);
      p4.push_back(pa.p.py);
      p4.push_back(pa.p.pz);
    }
    offsets[i + 1] = static_cast<std::int64_t>(pid.size());
  }
  py::dict d;
  d["offsets"] = move_array(std::move(offsets));
  d["pid"] = move_array(std::move(pid));
  d["charge"] = move_array(std::move(charge));
  d["p4"] = move_array2(std::move(p4), 4);
  d["x"] = move_array(std::move(x));
  d["q2"] = move_array(std::move(q2));
  d["y"] = move_array(std::move(y));
  d["kp"] = move_array2(std::move(kp), 4);
  d["weight"] = move_array(std::move(weight));
  return d;
}

// -------------------------------------------------------- sampler -> numpy

py::dict batch_to_dict(EventBatch&& b) {
  py::dict d;
  d["x"] = move_array(std::move(b.x));
  d["q2"] = move_array(std::move(b.q2));
  d["y"] = move_array(std::move(b.y));
  d["phi"] = move_array(std::move(b.phi));
  d["m"] = move_array(std::move(b.m));
  d["weight"] = move_array(std::move(b.weight));
  d["cell"] = move_array(std::move(b.cell));
  d["category"] = b.category;
  d["lam_e"] = b.lam_e;
  return d;
}

EventBatch dict_to_batch(const py::dict& d) {
  EventBatch b;
  auto get = [&](const char* key) {
    return d.contains(key)
               ? py::cast<std::vector<double>>(d[key])
               : std::vector<double>();
  };
  b.x = get("x");
  b.q2 = get("q2");
  b.y = get("y");
  b.phi = get("phi");
  b.m = get("m");
  b.weight = get("weight");
  if (d.contains("cell")) b.cell = py::cast<std::vector<int>>(d["cell"]);
  if (d.contains("category")) b.category = py::cast<std::string>(d["category"]);
  if (d.contains("lam_e")) b.lam_e = py::cast<int>(d["lam_e"]);
  return b;
}

}  // namespace

// =========================================================== module sections

static void bind_event(py::module_& m);
static void bind_beams(py::module_& m);
static void bind_spin(py::module_& m);
static void bind_sf(py::module_& m);
static void bind_xsec(py::module_& m);
static void bind_bookkeeping(py::module_& m);
static void bind_sampler(py::module_& m);
static void bind_spectator(py::module_& m);
static void bind_tagged(py::module_& m);
static void bind_cluster_config(py::module_& m);
static void bind_fsi(py::module_& m);
static void bind_rc(py::module_& m);
static void bind_coherent(py::module_& m);
static void bind_pipeline(py::module_& m);
static void bind_io(py::module_& m);

PYBIND11_MODULE(_lipolgen, m) {
  m.doc() = "LiPolGen: doubly polarized e + 6Li / 7Li DIS event generator";
  m.attr("__version__") = LIPOLGEN_VERSION;

  m.attr("HAVE_LHAPDF") = static_cast<bool>(LIPOLGEN_HAVE_LHAPDF);
  m.attr("HAVE_HEPMC3") = static_cast<bool>(LIPOLGEN_HAVE_HEPMC3);
  m.attr("HAVE_PYTHIA8") = static_cast<bool>(LIPOLGEN_HAVE_PYTHIA8);

  // constants.hpp -- one definition, mirrored here read-only
  m.attr("TENSOR_LL_SIGN") = TENSOR_LL_SIGN;
  m.attr("ALPHA_EM") = ALPHA_EM;
  m.attr("GEV2_TO_PB") = GEV2_TO_PB;
  m.attr("M_NUCLEON") = M_NUCLEON;
  m.attr("PROTON_MASS") = PROTON_MASS;
  // "THE one conversion between fm^-1 and GeV" (constants.hpp).  Mirrored
  // here so a Python test can build t = (q hbar c)^2 from the LIBRARY's
  // own value instead of retyping 0.19733 -- the C0 shape band is stated
  // in q [fm^-1] and every form factor takes t [GeV^2].
  m.attr("HBARC_GEV_FM") = HBARC_GEV_FM;
  m.attr("PROTON_TOP_MOMENTUM") = PROTON_TOP_MOMENTUM;
  m.attr("B1_PER_DEUTERON_TO_PER_NUCLEON") = B1_PER_DEUTERON_TO_PER_NUCLEON;
  m.attr("B1_MILLER_TABLE_TO_PER_NUCLEON") = B1_MILLER_TABLE_TO_PER_NUCLEON;
  m.attr("B1_CDKS_TABLE_TO_PER_NUCLEON") = B1_CDKS_TABLE_TO_PER_NUCLEON;
  m.attr("LI6_B1_RANK2_TRANSFER") = LI6_B1_RANK2_TRANSFER;
  m.attr("LI6_B1_LEGACY_TRANSFER") = LI6_B1_LEGACY_TRANSFER;
  m.attr("LI6_B1_PER_NUCLEON") = LI6_B1_PER_NUCLEON;
  m.attr("C_BAG") = C_BAG;
  m.attr("M_U") = M_U;
  m.attr("BETA_DEFAULT") = BETA_DEFAULT;
  // beams.hpp -- the 6Li cluster wave function, one source of truth for the
  // inclusive effective polarization and the tagged S/D interference
  m.attr("P_D_LI6") = P_D_LI6;
  m.attr("P_D_DEUTERON") = P_D_DEUTERON;
  m.attr("ALPHA_D_VECTOR_POLARIZATION") = ALPHA_D_VECTOR_POLARIZATION;
  m.attr("DEUTERON_VECTOR_POLARIZATION") = DEUTERON_VECTOR_POLARIZATION;
  m.attr("LI6_CLUSTER_POLARIZATION") = LI6_CLUSTER_POLARIZATION;
  m.attr("LI6_POLARIZATION_VMC_SIX_BODY") = LI6_POLARIZATION_VMC_SIX_BODY;
  m.def("vector_dilution_of", &vector_dilution_of, py::arg("p_d"),
        "1 - (3/2) P_D -- the ONE home of the vector-dilution expression "
        "(C5.5b).  ALPHA_D_VECTOR_POLARIZATION, DEUTERON_VECTOR_POLARIZATION, "
        "li6_cluster_polarization and the AV18 embedded deuteron "
        "(DEUTERON_AV18, tagged.hpp) are all this function.");
  m.def("li6_cluster_polarization", &li6_cluster_polarization,
        py::arg("p_d_alpha_d"), py::arg("p_d_deuteron"),
        "(1 - 1.5 P_D^{alpha-d})(1 - 1.5 P_D^d) -- the ONE home of the "
        "cluster-picture product.  LI6_CLUSTER_POLARIZATION is this at "
        "(P_D_LI6, P_D_DEUTERON).");
  m.attr("LI6_NAIVE_ONE_THIRD") = LI6_NAIVE_ONE_THIRD;

  // numerics.hpp -- the NumPy primitives, exposed for the bit-level tests
  m.def("pairwise_sum", [](const std::vector<double>& v) {
    return pairwise_sum(v);
  }, py::arg("a"));
  m.def("trapezoid", &trapezoid, py::arg("y"), py::arg("x"));
  m.def("np_interp", [](double x, const std::vector<double>& xp,
                        const std::vector<double>& fp) {
    return np_interp(x, xp, fp);
  }, py::arg("x"), py::arg("xp"), py::arg("fp"));

  // rng.hpp
  py::class_<Rng>(m, "Rng",
                  "Counter-based stream: the same (seed, run, bunch, event)\n"
                  "gives the same numbers on any thread count.")
      .def(py::init<std::uint64_t, std::uint64_t, std::uint64_t, std::uint64_t>(),
           py::arg("seed"), py::arg("run") = 0, py::arg("bunch") = 0,
           py::arg("event") = 0)
      .def("uniform", &Rng::uniform)
      .def("normal", &Rng::normal)
      .def("next_u64", &Rng::next_u64)
      .def("uniforms", [](Rng& r, std::size_t n) {
        std::vector<double> v(n);
        for (std::size_t i = 0; i < n; ++i) v[i] = r.uniform();
        return move_array(std::move(v));
      }, py::arg("n"), "n uniform deviates as a numpy array");

  bind_event(m);
  bind_beams(m);
  bind_spin(m);
  bind_sf(m);
  bind_xsec(m);
  bind_bookkeeping(m);
  bind_sampler(m);
  bind_spectator(m);
  bind_tagged(m);
  bind_cluster_config(m);
  bind_fsi(m);
  bind_rc(m);
  bind_coherent(m);
  bind_pipeline(m);
  bind_io(m);

  m.def("hfs_arrays", &hfs_arrays, py::arg("events"),
        py::arg("include_spectators") = false,
        py::arg("include_hadronic_x") = false,
        "Flat polligen `HFSSample` particle arrays (offsets, pid, charge, p4,\n"
        "x, q2, y, kp, weight) from a sequence of Event records.  The\n"
        "scattered electron is excluded from the particle list and carried in\n"
        "`kp`, exactly as tools/pythia8/gen_dis_hfs.py writes its samples.");
}

// ------------------------------------------------------------------- event

static void bind_event(py::module_& m) {
  py::class_<Vec4>(m, "Vec4", "(E, px, py, pz) in GeV, head-on frame")
      .def(py::init<>())
      .def(py::init([](double e, double px, double py, double pz) {
        return Vec4{e, px, py, pz};
      }), py::arg("e"), py::arg("px"), py::arg("py"), py::arg("pz"))
      .def_readwrite("e", &Vec4::e)
      .def_readwrite("px", &Vec4::px)
      .def_readwrite("py", &Vec4::py)
      .def_readwrite("pz", &Vec4::pz)
      .def("m2", &Vec4::m2)
      .def("pt", &Vec4::pt)
      .def("p", &Vec4::p)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def("as_array", [](const Vec4& v) {
        return move_array(std::vector<double>{v.e, v.px, v.py, v.pz});
      }, "numpy array (E, px, py, pz)")
      .def("__repr__", [](const Vec4& v) {
        char b[128];
        std::snprintf(b, sizeof b, "Vec4(e=%.6g, px=%.6g, py=%.6g, pz=%.6g)",
                      v.e, v.px, v.py, v.pz);
        return std::string(b);
      });

  py::enum_<Status>(m, "Status")
      .value("Beam", Status::Beam)
      .value("Final", Status::Final)
      .value("Decayed", Status::Decayed)
      .value("Intermediate", Status::Intermediate);

  py::enum_<Role>(m, "Role")
      .value("BeamElectron", Role::BeamElectron)
      .value("BeamIon", Role::BeamIon)
      .value("ScatteredElectron", Role::ScatteredElectron)
      .value("VirtualPhoton", Role::VirtualPhoton)
      .value("StruckCluster", Role::StruckCluster)
      .value("StruckNucleon", Role::StruckNucleon)
      .value("Spectator", Role::Spectator)
      .value("PartnerSpectator", Role::PartnerSpectator)
      .value("HadronicX", Role::HadronicX)
      .value("Hadron", Role::Hadron)
      .value("IntactRecoil", Role::IntactRecoil)
      .value("Pomeron", Role::Pomeron)
      .value("Other", Role::Other);

  py::enum_<Channel>(m, "Channel")
      .value("Inclusive", Channel::Inclusive)
      .value("TaggedLi6Alpha", Channel::TaggedLi6Alpha)
      .value("TaggedLi6D", Channel::TaggedLi6D)
      .value("TaggedLi7Alpha", Channel::TaggedLi7Alpha)
      .value("TaggedLi7T", Channel::TaggedLi7T)
      .value("TaggedDeuteronP", Channel::TaggedDeuteronP)
      .value("TaggedDeuteronN", Channel::TaggedDeuteronN)
      .value("TaggedHe3P", Channel::TaggedHe3P)
      .value("CoherentLi6", Channel::CoherentLi6);

  py::class_<Particle>(m, "Particle")
      .def(py::init<>())
      .def_readonly("pdg", &Particle::pdg)
      .def_readonly("status", &Particle::status)
      .def_readonly("role", &Particle::role)
      .def_readonly("p", &Particle::p)
      .def_readonly("mass", &Particle::mass)
      .def_readonly("charge", &Particle::charge)
      .def_readonly("mother1", &Particle::mother1)
      .def_readonly("mother2", &Particle::mother2)
      .def_readonly("pol", &Particle::pol)
      .def("__repr__", [](const Particle& p) {
        char b[192];
        std::snprintf(b, sizeof b,
                      "Particle(pdg=%d, status=%d, E=%.6g, pz=%.6g, q=%.3g)",
                      p.pdg, static_cast<int>(p.status), p.p.e, p.p.pz,
                      p.charge);
        return std::string(b);
      });

  py::class_<SpinLabels>(m, "SpinLabels")
      .def(py::init<>())
      .def_readonly("j", &SpinLabels::j)
      .def_readonly("m_ion", &SpinLabels::m_ion)
      .def_readonly("m_struck", &SpinLabels::m_struck)
      .def_readonly("lam_e", &SpinLabels::lam_e)
      .def_readonly("pe", &SpinLabels::pe)
      .def_readonly("theta_s", &SpinLabels::theta_s)
      .def_readonly("phi_s", &SpinLabels::phi_s)
      .def_readonly("pz", &SpinLabels::pz)
      .def_readonly("pzz", &SpinLabels::pzz)
      .def_readonly("category", &SpinLabels::category)
      .def_readonly("run", &SpinLabels::run)
      .def_readonly("bunch", &SpinLabels::bunch);

  py::class_<Kinematics>(m, "Kinematics")
      .def(py::init<>())
      .def_readonly("x", &Kinematics::x)
      .def_readonly("q2", &Kinematics::q2)
      .def_readonly("y", &Kinematics::y)
      .def_readonly("phi", &Kinematics::phi)
      .def_readonly("w2", &Kinematics::w2)
      .def_readonly("nu", &Kinematics::nu)
      .def_readonly("s", &Kinematics::s)
      .def_readonly("cell", &Kinematics::cell,
                    "Accepted-cell index of the sampler the (x, Q2) was "
                    "drawn from -- polligen's Mode-W `cell` column.")
      .def_readonly("k", &Kinematics::k)
      .def_readonly("cos_theta_k", &Kinematics::cos_theta_k)
      .def_readonly("phi_k", &Kinematics::phi_k)
      .def_readonly("alpha_s", &Kinematics::alpha_s)
      .def_readonly("pt_s", &Kinematics::pt_s)
      // the spectator's lab block, as `tagged.SpectatorLab` spells it
      .def_readonly("spec_pt", &Kinematics::spec_pt)
      .def_readonly("spec_theta", &Kinematics::spec_theta)
      .def_readonly("spec_p_lab", &Kinematics::spec_p_lab)
      .def_readonly("spec_r", &Kinematics::spec_r)
      .def_readonly("spec_xl", &Kinematics::spec_xl)
      .def_readonly("spec_kx", &Kinematics::spec_kx)
      .def_readonly("spec_ky", &Kinematics::spec_ky)
      .def_readonly("spec_kz", &Kinematics::spec_kz)
      .def_readonly("phi_spec", &Kinematics::phi_spec)
      .def_readonly("t", &Kinematics::t)
      .def_readonly("x_pom", &Kinematics::x_pom)
      .def_readonly("beta_pom", &Kinematics::beta_pom)
      .def_readonly("m_x2", &Kinematics::m_x2);

  py::class_<Event>(m, "Event")
      .def(py::init<>())
      .def_readonly("number", &Event::number)
      .def_readonly("channel", &Event::channel)
      .def_readonly("weight", &Event::weight)
      .def_readonly("spin_weights", &Event::spin_weights)
      .def_readonly("rc_weights", &Event::rc_weights,
                    "rc.hpp's opt-in RC weight block, row-major "
                    "(n_slot x 3) with n_slot = 1 + len(spin_weights): "
                    "[0:3] is (rc_tensor_lo, rc_tensor_hi, rc_tail) of the "
                    "event's OWN pure spin state, slot 1+k is spin category "
                    "k's population mixture.  EMPTY when the run had "
                    "--rc off, and then nothing downstream changes.")
      .def_readonly("rc_clipped", &Event::rc_clipped,
                    "Which RC ceilings SLOT 0 hit: bit 1 = the tail hit "
                    "RcOptions.tail_max, bit 2 = |tau| hit "
                    "RcOptions.band_tau_max.  0 with --rc off.  The npz "
                    "carries the run-level fractions as "
                    "meta['rc_clipped_tail_event_fraction'] and "
                    "meta['rc_clipped_band_event_fraction'].")
      .def_readonly("xsec_pb", &Event::xsec_pb)
      .def_readonly("xsec_err_pb", &Event::xsec_err_pb)
      .def_readonly("spin", &Event::spin)
      .def_readonly("kin", &Event::kin)
      .def_readonly("particles", &Event::particles)
      .def("find", [](const Event& ev, Role r) -> py::object {
        const Particle* p = ev.find(r);
        return p ? py::cast(*p) : py::none();
      }, py::arg("role"), "The first particle with this role, or None")
      .def("total_final", &Event::total_final)
      .def("total_charge_final", &Event::total_charge_final)
      .def("__len__", [](const Event& ev) { return ev.particles.size(); })
      .def("__repr__", [](const Event& ev) {
        char b[192];
        std::snprintf(b, sizeof b,
                      "Event(number=%llu, x=%.5g, Q2=%.5g, y=%.5g, n=%zu)",
                      static_cast<unsigned long long>(ev.number), ev.kin.x,
                      ev.kin.q2, ev.kin.y, ev.particles.size());
        return std::string(b);
      });
}

// ------------------------------------------------------------------- beams

static void bind_beams(py::module_& m) {
  m.def("nucleus_mass", &nucleus_mass, py::arg("name"), py::arg("a"),
        py::arg("z"));
  m.def("gamma_of", &gamma_of, py::arg("proton_energy"));
  m.def("epios_window_of", &epios_window_of, py::arg("proton_energy"),
        py::arg("shift_gamma") = 2.7, py::arg("tol") = 0.005);

  py::class_<Ion>(m, "Ion")
      .def(py::init<>())
      .def(py::init([](std::string name, int a, int z, double spin,
                       double pp, double pn) {
        Ion i;
        i.name = std::move(name);
        i.A = a;
        i.Z = z;
        i.spin = spin;
        i.eff_pol_p = pp;
        i.eff_pol_n = pn;
        return i;
      }), py::arg("name"), py::arg("A"), py::arg("Z"), py::arg("spin") = 0.5,
          py::arg("eff_pol_p") = 0.0, py::arg("eff_pol_n") = 0.0)
      .def_readwrite("name", &Ion::name)
      .def_readwrite("A", &Ion::A)
      .def_readwrite("Z", &Ion::Z)
      .def_readwrite("spin", &Ion::spin)
      .def_readwrite("eff_pol_p", &Ion::eff_pol_p)
      .def_readwrite("eff_pol_n", &Ion::eff_pol_n)
      .def_property_readonly("N", &Ion::N)
      .def("mass", &Ion::mass)
      .def("mass_per_nucleon", &Ion::mass_per_nucleon)
      .def("momentum_per_nucleon_max", &Ion::momentum_per_nucleon_max)
      .def("momentum_per_nucleon_at", &Ion::momentum_per_nucleon_at,
           py::arg("proton_energy"))
      .def("__repr__", [](const Ion& i) {
        return "Ion(" + i.name + ", A=" + std::to_string(i.A) + ", Z="
               + std::to_string(i.Z) + ")";
      });

  m.def("proton", []() { return PROTON(); });
  m.def("deuteron", []() { return DEUTERON(); });
  m.def("he3", []() { return HE3(); });
  m.def("li6", []() { return LI6(); });
  m.def("li7", []() { return LI7(); });
  m.def("ion_by_name", [](const std::string& n) { return ion_by_name(n); },
        py::arg("name"));
  m.attr("PROTON_CONFIG_ENERGIES") =
      py::make_tuple(PROTON_CONFIG_ENERGIES[0], PROTON_CONFIG_ENERGIES[1],
                     PROTON_CONFIG_ENERGIES[2]);
  m.attr("ELECTRON_ENERGIES") =
      py::make_tuple(ELECTRON_ENERGIES[0], ELECTRON_ENERGIES[1],
                     ELECTRON_ENERGIES[2]);

  py::class_<BeamConfig>(m, "BeamConfig")
      .def(py::init<>())
      .def(py::init([](double ee, Ion ion, double pu) {
        BeamConfig c;
        c.electron_energy = ee;
        c.ion = std::move(ion);
        c.ion_momentum_per_nucleon = pu;
        return c;
      }), py::arg("electron_energy"), py::arg("ion"),
          py::arg("ion_momentum_per_nucleon"))
      .def_readwrite("electron_energy", &BeamConfig::electron_energy)
      .def_readwrite("ion", &BeamConfig::ion)
      .def_readwrite("ion_momentum_per_nucleon",
                     &BeamConfig::ion_momentum_per_nucleon)
      .def("sqrt_s_per_nucleon", &BeamConfig::sqrt_s_per_nucleon)
      .def("s_per_nucleon", &BeamConfig::s_per_nucleon)
      .def("label", &BeamConfig::label)
      .def("__repr__", [](const BeamConfig& c) {
        return "BeamConfig(" + c.label() + ")";
      });

  m.def("default_configs", &default_configs, py::arg("ion_name") = "7Li");
}

// -------------------------------------------------------------------- spin

static void bind_spin(py::module_& m) {
  m.def("m_values", [](double j) { return move_array(m_values(j)); },
        py::arg("j"), "Spin projections ordered +J ... -J");
  m.def("clebsch_gordan", &clebsch_gordan, py::arg("j1"), py::arg("m1"),
        py::arg("j2"), py::arg("m2"), py::arg("j"), py::arg("m"));
  m.def("wigner_d", [](double j, double beta) {
    return matrix_to_array(wigner_d(j, beta));
  }, py::arg("j"), py::arg("beta"));
  m.def("rotation_matrix", [](double j, double theta, double phi) {
    return matrix_to_array(rotation_matrix(j, theta, phi));
  }, py::arg("j"), py::arg("theta"), py::arg("phi"));
  m.def("angular_momentum_ops", [](double j) {
    AngularMomentumOps o = angular_momentum_ops(j);
    return py::make_tuple(matrix_to_array(o.jx), matrix_to_array(o.jy),
                          matrix_to_array(o.jz));
  }, py::arg("j"));
  m.def("multipole_operator", [](double j, int k, int q) {
    return matrix_to_array(multipole_operator(j, k, q));
  }, py::arg("j"), py::arg("k"), py::arg("q"));
  m.def("rho_from_populations", [](double j, const std::vector<double>& p,
                                   double theta, double phi) {
    return matrix_to_array(rho_from_populations(j, p, theta, phi));
  }, py::arg("j"), py::arg("populations"), py::arg("theta") = 0.0,
     py::arg("phi") = 0.0);
  m.def("vector_polarization",
        [](const CplxArray& rho, double j) {
          return vector_polarization(array_to_cplx_matrix(rho), j);
        }, py::arg("rho"), py::arg("j"));
  m.def("tensor_polarization",
        [](const CplxArray& rho, double j) {
          return tensor_polarization(array_to_cplx_matrix(rho), j);
        }, py::arg("rho"), py::arg("j"));
  m.def("octupole_moment",
        [](const CplxArray& rho, double j) {
          return octupole_moment(array_to_cplx_matrix(rho), j);
        }, py::arg("rho"), py::arg("j") = 1.5);

  py::class_<AxisMoments>(m, "AxisMoments")
      .def_readonly("vector", &AxisMoments::vector)
      .def_readonly("tensor", &AxisMoments::tensor)
      .def_readonly("octupole", &AxisMoments::octupole)
      .def_readonly("has_octupole", &AxisMoments::has_octupole)
      .def("__repr__", [](const AxisMoments& a) {
        char b[128];
        std::snprintf(b, sizeof b, "AxisMoments(vector=%.10g, tensor=%.10g)",
                      a.vector, a.tensor);
        return std::string(b);
      });
  m.def("moments_along_axis", &moments_along_axis, py::arg("j"),
        py::arg("populations"));
  m.def("spin1_populations", [](double pz, double pzz) {
    auto p = spin1_populations(pz, pzz);
    return std::vector<double>(p.begin(), p.end());
  }, py::arg("pz"), py::arg("pzz"));
  m.def("spin32_populations", [](double pz, double t, double o) {
    auto p = spin32_populations(pz, t, o);
    return std::vector<double>(p.begin(), p.end());
  }, py::arg("pz"), py::arg("t") = 0.0, py::arg("o") = 0.0);
  m.def("populations_maxent", &populations_maxent, py::arg("j"), py::arg("pz"),
        py::arg("beta_max") = 60.0, py::arg("tol") = 1e-13,
        "Spin-temperature (maximum-entropy) populations p_m ~ exp(beta m)");

  py::class_<SpinDensity::LabMoments>(m, "LabMoments")
      .def_readonly("vector", &SpinDensity::LabMoments::vector)
      .def_readonly("tensor_zz", &SpinDensity::LabMoments::tensor_zz)
      .def_readonly("octupole_z", &SpinDensity::LabMoments::octupole_z)
      .def_readonly("has_octupole", &SpinDensity::LabMoments::has_octupole);

  py::class_<SpinDensity>(m, "SpinDensity")
      .def(py::init([](double j, std::vector<double> pops, double theta,
                       double phi) {
        SpinDensity s;
        s.j = j;
        s.populations = std::move(pops);
        s.theta = theta;
        s.phi = phi;
        return s;
      }), py::arg("j"), py::arg("populations"), py::arg("theta") = 0.0,
          py::arg("phi") = 0.0)
      .def_readwrite("j", &SpinDensity::j)
      .def_readwrite("populations", &SpinDensity::populations)
      .def_readwrite("theta", &SpinDensity::theta)
      .def_readwrite("phi", &SpinDensity::phi)
      .def("rho_lab", [](const SpinDensity& s) {
        return matrix_to_array(s.rho_lab());
      })
      .def("lab_moments", &SpinDensity::lab_moments);
}

// ---------------------------------------------------------------------- sf

static void bind_sf(py::module_& m) {
  m.def("r_sigma_lt", &r_sigma_lt, py::arg("x"), py::arg("q2"));
  py::enum_<R1998Form>(m, "R1998Form")
      .value("Average", R1998Form::kAverage)
      .value("A", R1998Form::kA)
      .value("B", R1998Form::kB)
      .value("C", R1998Form::kC);
  m.def("r1998", &r1998, py::arg("x"), py::arg("q2"),
        py::arg("form") = R1998Form::kAverage, py::arg("clip") = true);
  m.def("r1998_spread", &r1998_spread, py::arg("x"), py::arg("q2"),
        py::arg("clip") = true);
  m.def("r1998_fit_error", &r1998_fit_error, py::arg("x"), py::arg("q2"));

  py::class_<UnpolSF, std::shared_ptr<UnpolSF>>(m, "UnpolSF")
      .def("f2p", &UnpolSF::f2p, py::arg("x"), py::arg("q2"))
      .def("f2n", &UnpolSF::f2n, py::arg("x"), py::arg("q2"))
      .def("f2n_over_f2p", &UnpolSF::f2n_over_f2p, py::arg("x"))
      .def("f1p", &UnpolSF::f1p, py::arg("x"), py::arg("q2"))
      .def("f1n", &UnpolSF::f1n, py::arg("x"), py::arg("q2"))
      .def("flp", &UnpolSF::flp, py::arg("x"), py::arg("q2"))
      .def("fln", &UnpolSF::fln, py::arg("x"), py::arg("q2"))
      .def("r", &UnpolSF::r, py::arg("x"), py::arg("q2"))
      .def("set_r_func", &UnpolSF::set_r_func, py::arg("f"));
  py::class_<ToyF2, UnpolSF, std::shared_ptr<ToyF2>>(m, "ToyF2")
      .def(py::init<>());

  py::class_<NuclearF2>(m, "NuclearF2")
      .def(py::init([](Ion ion, std::shared_ptr<UnpolSF> base,
                       std::function<double(double)> emc, RFunc r) {
        return NuclearF2(std::move(ion), std::move(base), std::move(emc),
                         std::move(r));
      }), py::arg("ion"), py::arg("base"), py::arg("emc_ratio") = nullptr,
          py::arg("r_func") = nullptr)
      .def("f2a", &NuclearF2::f2a, py::arg("x"), py::arg("q2"))
      .def("f1a", &NuclearF2::f1a, py::arg("x"), py::arg("q2"))
      .def_property_readonly("ion", &NuclearF2::ion);

  m.def("dsigma_dx_dq2", [](double x, double q2, double s, double f2,
                            std::optional<double> fl, RFunc r) {
    const double v = fl ? *fl : 0.0;
    return dsigma_dx_dq2(x, q2, s, f2, fl ? &v : nullptr, r);
  }, py::arg("x"), py::arg("q2"), py::arg("s"), py::arg("f2"),
     py::arg("fl") = py::none(), py::arg("r_func") = nullptr);

  m.def("g2_ww", [](const std::function<double(double, double)>& g1, double x,
                    double q2, int npts) {
    return g2_ww(g1, x, q2, npts);
  }, py::arg("g1"), py::arg("x"), py::arg("q2"), py::arg("npts") = 96);

  py::class_<PolSF, std::shared_ptr<PolSF>>(m, "PolSF")
      .def("g1p", &PolSF::g1p, py::arg("x"), py::arg("q2"))
      .def("g1n", &PolSF::g1n, py::arg("x"), py::arg("q2"))
      .def("g1_nucleus", &PolSF::g1_nucleus, py::arg("ion"), py::arg("x"),
           py::arg("q2"), py::arg("medium_ratio") = nullptr)
      .def("g2p", &PolSF::g2p, py::arg("x"), py::arg("q2"), py::arg("npts") = 96)
      .def("g2n", &PolSF::g2n, py::arg("x"), py::arg("q2"), py::arg("npts") = 96);
  py::class_<ToyG1, PolSF, std::shared_ptr<ToyG1>>(m, "ToyG1")
      .def(py::init([](std::shared_ptr<UnpolSF> base, RFunc r) {
        return std::make_shared<ToyG1>(std::move(base), std::move(r));
      }), py::arg("base") = nullptr, py::arg("r_func") = nullptr)
      .def("a1p", &ToyG1::a1p, py::arg("x"))
      .def("a1n", &ToyG1::a1n, py::arg("x"));

#if LIPOLGEN_HAVE_LHAPDF
  m.def("lhapdf_quiet", &lhapdf_quiet);
  py::class_<LhapdfSF, UnpolSF, std::shared_ptr<LhapdfSF>>(m, "LhapdfSF")
      .def(py::init<std::string, int>(), py::arg("setname") = "CT18NLO",
           py::arg("member") = 0)
      .def_property_readonly("q2_min", &LhapdfSF::q2_min,
          "The loaded grid's OWN lower Q2 edge (LHAPDF's q2Min()), read from "
          "the set metadata and never transcribed: 1.677 for CT18NLO "
          "(QMin 1.295 GeV), 1 for NNPDFpol11_100.  Below it LHAPDF does NOT "
          "freeze -- it keeps evolving downward and F2p falls fast -- which "
          "is why unpol_sf_grid_report() exists and why the run banner and "
          "meta['unpol_sf_below_grid_frac'] are computed from this.");
  py::class_<LhapdfG1, PolSF, std::shared_ptr<LhapdfG1>>(m, "LhapdfG1")
      .def(py::init<std::string, int>(), py::arg("setname") = "NNPDFpol11_100",
           py::arg("member") = 0);
  py::class_<Epps21Ratio>(m, "Epps21Ratio")
      .def(py::init<std::string, int, std::string, int>(),
           py::arg("nuclear_set") = "EPPS21nlo_CT18Anlo_Li6",
           py::arg("nuclear_member") = 0, py::arg("proton_set") = "CT18NLO",
           py::arg("proton_member") = 0)
      .def("ratio", &Epps21Ratio::ratio, py::arg("pid"), py::arg("x"),
           py::arg("q2"))
      .def("f2_per_nucleon", &Epps21Ratio::f2_per_nucleon, py::arg("x"),
           py::arg("q2"))
      .def("f2a", &Epps21Ratio::f2a, py::arg("x"), py::arg("q2"), py::arg("A"));
#endif

  // MSTW2008 LO, read from PYTHIA's own pdfdata grid store -- bound HERE,
  // inside bind_sf, and not in the `#if LIPOLGEN_HAVE_PYTHIA8` block further
  // down, because bind_sf runs before that block and `UnpolSF` (the base it
  // is registered against) is declared in this function.  Same placement
  // rule LhapdfSF follows for LHAPDF just above.
#if LIPOLGEN_HAVE_PYTHIA8
  m.def("pythia8_pdfdata_dir", &pythia8_pdfdata_dir,
        "PYTHIA 8's pdfdata grid directory (compiled-in default, overridden "
        "by $LIPOLGEN_PYTHIA8_PDFDATA).");
  py::class_<MstwSF, UnpolSF, std::shared_ptr<MstwSF>>(m, "MstwSF")
      .def(py::init<int, std::string>(), py::arg("i_fit") = 3,
           py::arg("pdfdata_path") = "")
      .def_property_readonly("i_fit", &MstwSF::i_fit)
      .def_property_readonly("pdfdata_path", &MstwSF::pdfdata_path);
#endif

  py::enum_<B1Mode>(m, "B1Mode")
      .value("Digitized", B1Mode::kDigitized)
      .value("Toy", B1Mode::kToy);
  py::enum_<EmcMode>(m, "EmcMode")
      .value("Digitized", EmcMode::kDigitized)
      .value("Constant", EmcMode::kConstant);
  py::enum_<EmcBaseline>(m, "EmcBaseline")
      .value("LegacyTable", EmcBaseline::LegacyTable)
      .value("Epps21", EmcBaseline::Epps21);
  m.attr("EMC_BASELINE_DEFAULT") = EMC_BASELINE_DEFAULT;
  m.attr("EMC_VALENCE_DEPLETION_EPPS21") = EMC_VALENCE_DEPLETION_EPPS21;

  m.def("toy_b1_shape", &toy_b1_shape, py::arg("x"), py::arg("q2"), py::arg("f1"));
  m.def("toy_b1", &toy_b1, py::arg("x"), py::arg("q2"), py::arg("f1"),
        py::arg("mode") = B1Mode::kDigitized);
  m.def("b1_convolution", &b1_convolution, py::arg("x"), py::arg("q2"),
        py::arg("f1"), py::arg("mode") = B1Mode::kDigitized);
  m.def("close_kumano_integral", &close_kumano_integral, py::arg("cdks") = true);
  m.def("b1_li6_from_deuteron", &b1_li6_from_deuteron, py::arg("b1_d"),
        py::arg("transfer") = LI6_B1_RANK2_TRANSFER,
        py::arg("per_nucleon") = LI6_B1_PER_NUCLEON);
  m.def("toy_delta_gluon", &toy_delta_gluon, py::arg("x"), py::arg("q2"),
        py::arg("f1"), py::arg("scale") = 1e-3);
  m.def("unpolarized_emc_ratio", &unpolarized_emc_ratio, py::arg("x"));
  m.def("cbt_unpolarized_emc_ratio", &cbt_unpolarized_emc_ratio, py::arg("x"));
  m.def("cbt_polarized_emc_ratio", &cbt_polarized_emc_ratio, py::arg("x"),
        py::arg("mode") = EmcMode::kDigitized, py::arg("eq") = 23,
        py::arg("baseline") = EMC_BASELINE_DEFAULT);
  m.def("tmt_polarized_emc_ratio", &tmt_polarized_emc_ratio, py::arg("x"),
        py::arg("mode") = EmcMode::kDigitized,
        py::arg("baseline") = EMC_BASELINE_DEFAULT);
  m.def("cbt_published_emc_ratio", &cbt_published_emc_ratio, py::arg("x"),
        py::arg("eq") = 23);
  m.def("tmt_published_emc_ratio", &tmt_published_emc_ratio, py::arg("x"));
  m.def("cbt_ratio_of_effects", &cbt_ratio_of_effects, py::arg("x"),
        py::arg("eq") = 23);
  m.def("tmt_ratio_of_effects", &tmt_ratio_of_effects, py::arg("x"));
  m.def("emc_valence_depletion", &emc_valence_depletion, py::arg("baseline"));
  m.def("cbt_valence_scale", &cbt_valence_scale,
        py::arg("baseline") = EMC_BASELINE_DEFAULT);
  m.def("tmt_valence_scale", &tmt_valence_scale,
        py::arg("baseline") = EMC_BASELINE_DEFAULT);

  py::class_<TensorSF, std::shared_ptr<TensorSF>>(m, "TensorSF")
      .def("b1", &TensorSF::b1, py::arg("x"), py::arg("q2"), py::arg("f1"))
      .def("b2", &TensorSF::b2, py::arg("x"), py::arg("q2"), py::arg("f1"))
      .def("delta", &TensorSF::delta, py::arg("x"), py::arg("q2"), py::arg("f1"))
      .def("b1_func", &TensorSF::b1_func)
      .def("b2_func", &TensorSF::b2_func)
      .def("delta_func", &TensorSF::delta_func);
  py::class_<MillerB1, TensorSF, std::shared_ptr<MillerB1>>(m, "MillerB1")
      .def(py::init<B1Mode>(), py::arg("mode") = B1Mode::kDigitized);
  py::class_<CdksB1, TensorSF, std::shared_ptr<CdksB1>>(m, "CdksB1")
      .def(py::init<B1Mode>(), py::arg("mode") = B1Mode::kDigitized);
  py::class_<Li6B1, TensorSF, std::shared_ptr<Li6B1>>(m, "Li6B1")
      .def(py::init([](std::shared_ptr<TensorSF> d, double transfer,
                       double per_nucleon) {
        return std::make_shared<Li6B1>(std::move(d), transfer, per_nucleon);
      }), py::arg("deuteron"), py::arg("transfer") = LI6_B1_RANK2_TRANSFER,
          py::arg("per_nucleon") = LI6_B1_PER_NUCLEON);

  // ================================================= b1_nuclear.hpp (design D)
  // The four-term alpha-d convolution b1 of 6Li and its A = 2 validation gate.
  // The gate PASSES since 2026-09-03 -- for ONE unpolarised nucleon input,
  // not for the shipped default: G3b ratio 0.843 with CDKS's own MSTW2008 LO
  // at their Eq. (21), against 0.440 on the DEFAULT ToyF2, which is OUTSIDE
  // the gate's [0.5, 2] window.  Quote numbers made with b1_unpol = Mstw
  // (CLI --b1-unpol mstw); the default toy backend's are not covered by the
  // lift.  With CD-Bonn as well the ratio is 1.000338 -- a residual below the
  // error of digitizing a published figure, NOT three-digit agreement with
  // CDKS, and specific to CD-Bonn AND MSTW together.  See the header block of
  // b1_nuclear.hpp for the four conditions that come with that, and
  // docs/open_items/run_2026-09-03/phase_A_numbers.md for the measurements.
  // Still opt-in and still BAND-ONLY: the band comes from Q(6Li) vs Q_d, not
  // from the gate, so never quote a single row.
  m.def("f1_cdks", &f1_cdks, py::arg("f2"), py::arg("x"), py::arg("q2"),
        py::arg("r_func") = nullptr,
        "CDKS Eq. (22) F1 = (1 + gamma^2) F2 / (2 x (1 + R)).  NOT "
        "`UnpolSF.f1_from_f2` and NOT `NuclearF2.f1a`, which are the massless "
        "form; the extra factor is 1 + gamma_squared(x, q2) and is 1.90 at "
        "x = 0.8, Q2 = 2.5.");
  m.def("cdks_b1_raw_per_nucleon", []() {
    return std::const_pointer_cast<TensorSF>(cdks_b1_raw_per_nucleon());
  }, "The RAW digitized CDKS theory-1 column as a TensorSF.  Since 2026-09-03 "
     "`CdksB1` carries the same normalisation -- the column is per nucleon by "
     "CDKS Eq. (10) and the text under their Eq. (16), and "
     "`B1_CDKS_TABLE_TO_PER_NUCLEON` is 1 -- so this accessor is now a second "
     "name for it rather than a route around a double halving.  This is "
     "`Li6ConvolutionOptions.deuteron_b1 = None`.");

  py::class_<CdBonnWave>(m, "CdBonnWave",
      "The CD-Bonn deuteron wave function (Machleidt, PRC 63 (2001) 024001, "
      "Appendix D, Table XX).  Units are the PAPER's: p in fm^-1, psi in "
      "fm^3/2, (2/pi) int p^2 (u^2 + w^2) dp = 1.  `psi_d` carries the "
      "`fdeut.av18` / CDKS sign w = -psi_2^a.")
      .def_property_readonly("m", [](const CdBonnWave& w) {
        return copy_array(std::vector<double>(w.m.begin(), w.m.end()));
      }, "m_j [fm^-1], Eq. (D25).")
      .def_property_readonly("c", [](const CdBonnWave& w) {
        return copy_array(std::vector<double>(w.c.begin(), w.c.end()));
      }, "C_j [fm^-1/2]; C_11 is COMPUTED from Eq. (D23), not published.")
      .def_property_readonly("d", [](const CdBonnWave& w) {
        return copy_array(std::vector<double>(w.d.begin(), w.d.end()));
      }, "D_j [fm^-1/2]; D_9..D_11 are COMPUTED from Eq. (D24).")
      .def("psi_s", &CdBonnWave::psi_s, py::arg("p_fm"))
      .def("psi_d", &CdBonnWave::psi_d, py::arg("p_fm"))
      .def("norm_s", &CdBonnWave::norm_s)
      .def("norm_d", &CdBonnWave::norm_d,
           "The D-state probability, CLOSED FORM.  CD-Bonn publishes 4.85 %.")
      .def("norm", &CdBonnWave::norm, "= 1 by Eq. (D14).")
      .def("a_s", &CdBonnWave::a_s)
      .def("a_d", &CdBonnWave::a_d)
      .def("eta", &CdBonnWave::eta,
           "A_D/A_S, the asymptotic D/S ratio.  CD-Bonn publishes 0.0256.")
      .def("constraint_residuals", [](const CdBonnWave& w) {
        const std::array<double, 4> r = w.constraint_residuals();
        return copy_array(std::vector<double>(r.begin(), r.end()));
      }, "{sum C_j, sum D_j/m_j^2, sum D_j, sum D_j m_j^2}, all zero by "
         "construction -- Machleidt's \"about 15 decimal digits\".");
  m.def("cdbonn_wave", &cdbonn_wave,
        py::return_value_policy::reference,
        "The one CD-Bonn wave function: Table XX with C_11 and D_9..D_11 "
        "solved for.");
  m.def("cdbonn_fdeut_table", [](double k_max_fm, double dk_fm) {
    const FdeutTable t = cdbonn_fdeut_table(k_max_fm, dk_fm);
    py::dict out;
    out["k_gev"] = copy_array(t.k_gev);
    out["u"] = copy_array(t.u);
    out["w"] = copy_array(t.w);
    out["ebind_gev"] = t.ebind_gev;
    return out;
  }, py::arg("k_max_fm") = 20.0, py::arg("dk_fm") = 0.1,
     "CD-Bonn sampled on 0 .. k_max_fm in dk_fm [fm^-1] and converted to the "
     "library's GeV units exactly as `read_fdeut_k` converts the AV18 file.  "
     "The defaults ARE `fdeut.av18`'s grid, so an AV18 row and a CD-Bonn row "
     "differ in the wave function alone.");

  py::class_<ClusterPartialWave>(m, "ClusterPartialWave")
      .def(py::init<>())
      .def_readwrite("l", &ClusterPartialWave::l)
      .def_property("k",
          [](const ClusterPartialWave& w) { return copy_array(w.k); },
          [](ClusterPartialWave& w, std::vector<double> v) {
            w.k = std::move(v); w.rebuild();
          })
      .def_property("phi",
          [](const ClusterPartialWave& w) { return copy_array(w.phi); },
          [](ClusterPartialWave& w, std::vector<double> v) {
            w.phi = std::move(v); w.rebuild();
          })
      .def_readwrite("provenance", &ClusterPartialWave::provenance)
      .def("__call__", &ClusterPartialWave::operator(), py::arg("k"))
      .def("norm2", &ClusterPartialWave::norm2)
      .def("scaled", &ClusterPartialWave::scaled, py::arg("factor"))
      .def("rebuild", &ClusterPartialWave::rebuild)
      .def_static("from_uw", &ClusterPartialWave::from_uw, py::arg("k_gev"),
                  py::arg("uw"), py::arg("l"));
  m.def("li6_alpha_d_partial_waves", &li6_alpha_d_partial_waves,
        "(phi_0, phi_2) of the 6Li alpha-d relative motion in the CDKS i^L "
        "convention: magnitudes from the ANL momentum file, sign structure "
        "from the overlap file.");
  m.def("alpha_d_quadrupole_fm2", &alpha_d_quadrupole_fm2, py::arg("phi0"),
        py::arg("phi2"), py::arg("r_max_fm") = 20.0, py::arg("n_r") = 2000,
        py::arg("n_k") = 4001,
        "THE SIGN GATE (design 2.2): the quadrupole moment of the unit-"
        "normalised relative wave function.  Must be NEGATIVE for alpha-d and "
        "positive for the deuteron, from the same code path.");

  py::class_<AlphaTQuadrupole>(m, "AlphaTQuadrupole",
      "The A = 7 gate's four numbers: <r^2> [fm^2], r_rms [fm], Z_eff and "
      "Q = -(2/5) Z_eff <r^2> [fm^2] of the alpha-t relative motion.  Compare "
      "with LI7_QUADRUPOLE_FM2; it validates the alpha-t WAVE FUNCTION's "
      "quadrupole only, and b1(7Li) remains unimplemented (open item 15).")
      .def_readonly("r2_fm2", &AlphaTQuadrupole::r2_fm2)
      .def_readonly("r_rms_fm", &AlphaTQuadrupole::r_rms_fm)
      .def_readonly("z_eff", &AlphaTQuadrupole::z_eff)
      .def_readonly("q_fm2", &AlphaTQuadrupole::q_fm2);
  m.def("alpha_t_quadrupole", &alpha_t_quadrupole, py::arg("phi1"),
        py::arg("r_max_fm") = 30.0, py::arg("n_r") = 4000,
        py::arg("n_k") = 8001,
        "Q = -(2/5) Z_eff <r^2> for an arbitrary L = 1 relative wave, read "
        "for its (k, psi) table only (j_1 is applied here).");
  m.def("li7_alpha_t_quadrupole", &li7_alpha_t_quadrupole,
        py::arg("vmc_mc_sigma") = 0.0, py::arg("r_max_fm") = 30.0,
        py::arg("n_r") = 4000, py::arg("n_k") = 8001,
        "THE A = 7 GATE on 7Li's own shipped alpha-t wave "
        "(momenta/li7_at3.momentum, the 3/2- ground state).  It validates the "
        "WAVE FUNCTION's quadrupole against the measured LI7_QUADRUPOLE_FM2 "
        "and nothing else: b1(7Li) is NOT implemented (open item 15) and "
        "7Li's rank-2 sector is exactly zero by construction.");

  py::class_<ConvolutionKinematics>(m, "ConvolutionKinematics")
      .def(py::init<>())
      .def_readwrite("m_struck", &ConvolutionKinematics::m_struck)
      .def_readwrite("m_recoil", &ConvolutionKinematics::m_recoil)
      .def_readwrite("separation", &ConvolutionKinematics::separation)
      .def_readwrite("kappa", &ConvolutionKinematics::kappa)
      .def("k_range", [](const ConvolutionKinematics& k, double y)
                          -> std::optional<std::pair<double, double>> {
        double lo = 0.0, hi = 0.0;
        if (!k.k_range(y, &lo, &hi)) return std::nullopt;
        return std::make_pair(lo, hi);
      }, py::arg("y"), "(lo, hi) in GeV, or None when the interval is empty.")
      .def("cos_star", &ConvolutionKinematics::cos_star, py::arg("k"),
           py::arg("y"))
      .def("y_max", &ConvolutionKinematics::y_max);

  py::class_<LightConeDensities> lcd(m, "LightConeDensities");
  py::class_<LightConeDensities::Options>(lcd, "Options")
      .def(py::init<>())
      .def_readwrite("n_y_low", &LightConeDensities::Options::n_y_low)
      .def_readwrite("n_y_mid", &LightConeDensities::Options::n_y_mid)
      .def_readwrite("n_y_high", &LightConeDensities::Options::n_y_high)
      .def_readwrite("y_mid_lo", &LightConeDensities::Options::y_mid_lo)
      .def_readwrite("y_mid_hi", &LightConeDensities::Options::y_mid_hi)
      .def_readwrite("n_k", &LightConeDensities::Options::n_k,
                     "k points per y.  MUST BE ODD -- the inner k integral is "
                     "composite Simpson and an even count would fall back to "
                     "the trapezoid, which is 22 % low at x = 0.05.  An even "
                     "value is bumped to the next odd one at construction and "
                     "`options()` then reports what actually ran.")
      .def_readwrite("y_min", &LightConeDensities::Options::y_min)
      .def_readwrite("renormalize", &LightConeDensities::Options::renormalize)
      .def_readwrite("norm_target", &LightConeDensities::Options::norm_target);
  lcd.def(py::init<ClusterPartialWave, ClusterPartialWave,
                   ConvolutionKinematics, LightConeDensities::Options>(),
          py::arg("phi0"), py::arg("phi2"), py::arg("kin"), py::arg("options"))
      .def(py::init<ClusterPartialWave, ClusterPartialWave,
                    ConvolutionKinematics>(),
           py::arg("phi0"), py::arg("phi2"), py::arg("kin"))
      .def_static("delta_limit", &LightConeDensities::delta_limit,
                  py::arg("p_d"), py::arg("norm"))
      .def("f_s", &LightConeDensities::f_s, py::arg("y"))
      .def("f_d", &LightConeDensities::f_d, py::arg("y"))
      .def("f_unpol", &LightConeDensities::f_unpol, py::arg("y"))
      .def("delta_t_f", &LightConeDensities::delta_t_f, py::arg("y"))
      .def("delta_t_f_sd", &LightConeDensities::delta_t_f_sd, py::arg("y"))
      .def("delta_t_f_dd", &LightConeDensities::delta_t_f_dd, py::arg("y"))
      .def("f_d_p2", &LightConeDensities::f_d_p2, py::arg("y"))
      .def("f_d_p4", &LightConeDensities::f_d_p4, py::arg("y"))
      .def("norm", &LightConeDensities::norm)
      .def("mean_y", &LightConeDensities::mean_y)
      .def("p_d", &LightConeDensities::p_d)
      .def("p_d_momentum", &LightConeDensities::p_d_momentum)
      .def("renormalization", &LightConeDensities::renormalization)
      .def("y_grid", [](const LightConeDensities& d) {
        return copy_array(d.y_grid());
      })
      .def_property_readonly("kinematics", &LightConeDensities::kinematics)
      .def_property_readonly("options", &LightConeDensities::options)
      .def("convolve", &LightConeDensities::convolve, py::arg("dens_at"),
           py::arg("g"), py::arg("x"), py::arg("x_max_g") = 1.0,
           py::arg("n") = 4001,
           "int (dy/y) dens_at(y) g(x/y), starting at x / x_max_g.  "
           "`x_max_g` is the x above which `g` is ZERO -- 1.0 for a "
           "free-nucleon F1 (the default), tables.kB1CdksQ2p5().x_max() = "
           "1.59 for the raw digitized CDKS b1_d column.  Cutting at "
           "x/y = 1 instead drops 5 % of term (1) at x = 0.95.");

  py::class_<Li6ConvolutionOptions>(m, "Li6ConvolutionOptions")
      .def(py::init<>())
      // The two `shared_ptr<const T>` fields go through def_property with a
      // const_pointer_cast: `def_readwrite` on a shared_ptr<const T> does not
      // round-trip with the registered shared_ptr<TensorSF> / <UnpolSF>
      // holders.  Same pattern as `PipelineConfig::kernel`.
      .def_property("deuteron_b1",
          [](const Li6ConvolutionOptions& o) {
            return std::const_pointer_cast<TensorSF>(o.deuteron_b1);
          },
          [](Li6ConvolutionOptions& o, std::shared_ptr<TensorSF> v) {
            o.deuteron_b1 = std::move(v);
          })
      .def_property("unpol",
          [](const Li6ConvolutionOptions& o) {
            return std::const_pointer_cast<UnpolSF>(o.unpol);
          },
          [](Li6ConvolutionOptions& o, std::shared_ptr<UnpolSF> v) {
            o.unpol = std::move(v);
          })
      .def_readwrite("deuteron_b1_x_max",
                     &Li6ConvolutionOptions::deuteron_b1_x_max,
                     "The x above which `deuteron_b1` is ZERO (its support, "
                     "which convolve needs and a TensorSF cannot be asked "
                     "for).  0 => automatic: 1.59 (the digitized table's own "
                     "x_max) when deuteron_b1 is None, else 1.0.")
      .def_readwrite("alpha_f1", &Li6ConvolutionOptions::alpha_f1)
      .def_readwrite("r_func", &Li6ConvolutionOptions::r_func)
      .def_readwrite("w_embedded_s", &Li6ConvolutionOptions::w_embedded_s)
      .def_readwrite("w_alpha_d_dwave",
                     &Li6ConvolutionOptions::w_alpha_d_dwave)
      .def_readwrite("w_cg_dwave", &Li6ConvolutionOptions::w_cg_dwave)
      .def_readwrite("use_spectroscopic_factor",
                     &Li6ConvolutionOptions::use_spectroscopic_factor)
      .def_readwrite("norm_target", &Li6ConvolutionOptions::norm_target)
      .def_readwrite("finite_q_delta", &Li6ConvolutionOptions::finite_q_delta)
      .def_readwrite("quad", &Li6ConvolutionOptions::quad);

  py::class_<Li6ConvolutionB1, TensorSF, std::shared_ptr<Li6ConvolutionB1>>(
      m, "Li6ConvolutionB1",
      "b1 of 6Li per nucleon, four-term alpha-d convolution (CDKS Eqs. 16/17/"
      "21 one level up).  OPT-IN.  Its A = 2 magnitude gate PASSES since "
      "2026-09-03 for ONE unpolarised nucleon input: ratio 0.843 with "
      "MSTW2008 LO at CDKS Eq. (21) (`Li6ConvolutionOptions.unpol = "
      "MstwSF()`, CLI --b1-unpol mstw).  On the DEFAULT `unpol` -- the "
      "kernel's own ToyF2 -- it is 0.440, OUTSIDE the [0.5, 2] window, and "
      "those numbers are NOT covered by the lift.  With CD-Bonn as well the "
      "ratio is 1.000338, a residual below the error of digitizing a "
      "published figure rather than three-digit agreement with CDKS.  And "
      "that gate is A = 2 and tests nothing about the alpha-d step: band "
      "every number, {0, 1, 2} x b1.")
      .def(py::init([](Li6ConvolutionOptions o) {
        return std::make_shared<Li6ConvolutionB1>(std::move(o));
      }), py::arg("options") = Li6ConvolutionOptions())
      .def("b1_embedded_s", &Li6ConvolutionB1::b1_embedded_s, py::arg("x"),
           py::arg("q2"), py::arg("f1") = 0.0)
      .def("b1_alpha_d_dwave_d", &Li6ConvolutionB1::b1_alpha_d_dwave_d,
           py::arg("x"), py::arg("q2"), py::arg("f1") = 0.0)
      .def("b1_alpha_d_dwave_alpha",
           &Li6ConvolutionB1::b1_alpha_d_dwave_alpha, py::arg("x"),
           py::arg("q2"), py::arg("f1") = 0.0)
      .def("b1_cg_dwave", &Li6ConvolutionB1::b1_cg_dwave, py::arg("x"),
           py::arg("q2"), py::arg("f1") = 0.0)
      .def("b1_alpha_d_dwave", &Li6ConvolutionB1::b1_alpha_d_dwave,
           py::arg("x"), py::arg("q2"), py::arg("f1") = 0.0)
      .def("b1_alpha_d_sd", &Li6ConvolutionB1::b1_alpha_d_sd, py::arg("x"),
           py::arg("q2"), py::arg("f1") = 0.0)
      .def("b1_alpha_d_dd", &Li6ConvolutionB1::b1_alpha_d_dd, py::arg("x"),
           py::arg("q2"), py::arg("f1") = 0.0)
      .def("b1_cg_dwave_p2_remainder",
           &Li6ConvolutionB1::b1_cg_dwave_p2_remainder, py::arg("x"),
           py::arg("q2"), py::arg("f1") = 0.0)
      .def("b1_cg_dwave_p4_remainder",
           &Li6ConvolutionB1::b1_cg_dwave_p4_remainder, py::arg("x"),
           py::arg("q2"), py::arg("f1") = 0.0)
      .def("banded", [](const Li6ConvolutionB1& b, double scale) {
        return std::const_pointer_cast<TensorSF>(b.banded(scale));
      }, py::arg("scale"),
         "A copy with the WHOLE b1 scaled: the mandatory 100 % band.")
      .def("band_scale", &Li6ConvolutionB1::band_scale)
      .def("f1_deuteron", &Li6ConvolutionB1::f1_deuteron, py::arg("x"),
           py::arg("q2"))
      .def("f1_alpha", &Li6ConvolutionB1::f1_alpha, py::arg("x"), py::arg("q2"))
      .def("deuteron_b1", [](const Li6ConvolutionB1& b) {
        return std::const_pointer_cast<TensorSF>(b.deuteron_b1());
      })
      .def("close_kumano_integral", &Li6ConvolutionB1::close_kumano_integral,
           py::arg("lo") = 0.01, py::arg("hi") = 1.2, py::arg("q2") = 2.5,
           py::arg("n") = 241)
      // reference_internal: these return references INTO an object the
      // shared_ptr owns, so the returned view must keep the parent alive.
      .def("densities", &Li6ConvolutionB1::densities,
           py::return_value_policy::reference_internal)
      .def("densities_alpha", &Li6ConvolutionB1::densities_alpha,
           py::return_value_policy::reference_internal)
      .def("options", &Li6ConvolutionB1::options,
           py::return_value_policy::reference_internal);

  py::enum_<DeuteronWaveSource>(m, "DeuteronWaveSource",
      "Which deuteron wave function `DeuteronConvolutionB1` convolves.  The "
      "DEFAULT IS `kFdeutFile` and it is a CHOICE: every published gate "
      "number was measured on AV18, so moving it would silently move all of "
      "them.  `kCdBonn` is what CDKS actually used for Fig. 4.")
      .value("kFdeutFile", DeuteronWaveSource::kFdeutFile)
      .value("kCdBonn", DeuteronWaveSource::kCdBonn);

  py::class_<DeuteronConvolutionB1, TensorSF,
             std::shared_ptr<DeuteronConvolutionB1>> dcb1(
      m, "DeuteronConvolutionB1",
      "THE A = 2 VALIDATION GATE: the same kernel on the AV18 deuteron u(k), "
      "w(k) and an isoscalar nucleon F1, to be compared with the digitized "
      "CDKS Fig. 4 (`tables::kB1CdksQ2p5`).");
  py::class_<DeuteronConvolutionB1::Options>(dcb1, "Options")
      .def(py::init<>())
      .def_readwrite("fdeut_path", &DeuteronConvolutionB1::Options::fdeut_path)
      .def_readwrite("wave", &DeuteronConvolutionB1::Options::wave,
                     "Which deuteron wave function: `DeuteronWaveSource."
                     "kFdeutFile` (the default, `fdeut_path`) or `kCdBonn`.")
      .def_readwrite("cdbonn_k_max_fm",
                     &DeuteronConvolutionB1::Options::cdbonn_k_max_fm)
      .def_readwrite("cdbonn_dk_fm",
                     &DeuteronConvolutionB1::Options::cdbonn_dk_fm)
      .def_property("unpol",
          [](const DeuteronConvolutionB1::Options& o) {
            return std::const_pointer_cast<UnpolSF>(o.unpol);
          },
          [](DeuteronConvolutionB1::Options& o, std::shared_ptr<UnpolSF> v) {
            o.unpol = std::move(v);
          })
      .def_readwrite("r_func", &DeuteronConvolutionB1::Options::r_func)
      .def_readwrite("target_mass",
                     &DeuteronConvolutionB1::Options::target_mass)
      .def_readwrite("finite_q_delta",
                     &DeuteronConvolutionB1::Options::finite_q_delta)
      .def_readwrite("quad", &DeuteronConvolutionB1::Options::quad);
  dcb1.def(py::init([](DeuteronConvolutionB1::Options o) {
             return std::make_shared<DeuteronConvolutionB1>(std::move(o));
           }), py::arg("options") = DeuteronConvolutionB1::Options())
      .def("b1_sd", &DeuteronConvolutionB1::b1_sd, py::arg("x"), py::arg("q2"),
           py::arg("f1") = 0.0)
      .def("b1_dd", &DeuteronConvolutionB1::b1_dd, py::arg("x"), py::arg("q2"),
           py::arg("f1") = 0.0)
      .def("f1_nucleon", &DeuteronConvolutionB1::f1_nucleon, py::arg("x"),
           py::arg("q2"))
      .def("densities", &DeuteronConvolutionB1::densities,
           py::return_value_policy::reference_internal)
      .def("options", &DeuteronConvolutionB1::options,
           py::return_value_policy::reference_internal);

  py::class_<B1Landmarks>(m, "B1Landmarks")
      .def(py::init<>())
      .def_property_readonly("zeros", [](const B1Landmarks& l) {
        return copy_array(l.zeros);
      })
      .def_readonly("zero_slope", &B1Landmarks::zero_slope)
      .def_readonly("x_min", &B1Landmarks::x_min)
      .def_readonly("xb1_min", &B1Landmarks::xb1_min)
      .def_readonly("x_max", &B1Landmarks::x_max)
      .def_readonly("xb1_max", &B1Landmarks::xb1_max)
      .def_readonly("integral_b1", &B1Landmarks::integral_b1);
  m.def("b1_landmarks", &b1_landmarks, py::arg("xb1"), py::arg("x_grid"),
        "Sign changes (with slope), extrema of x*b1 and int b1 dx -- COMPUTED, "
        "never typed.");
  m.def("b1_landmarks_of_table", &b1_landmarks_of_table,
        "The same, applied to the digitized `tables::kB1CdksQ2p5` raw column, "
        "so the gate compares like with like.");

  py::class_<DeltaVariant>(m, "DeltaVariant")
      .def_readonly("alpha", &DeltaVariant::alpha)
      .def_readonly("beta", &DeltaVariant::beta);
  m.def("delta_variant", &delta_variant, py::arg("name"));
  m.def("xab_peak_value", &xab_peak_value, py::arg("alpha"), py::arg("beta"));
  m.def("shape_normalized", &shape_normalized, py::arg("x"), py::arg("variant"));
  m.def("alpha_s_lo", &alpha_s_lo, py::arg("q2"), py::arg("lambda_qcd") = 0.22,
        py::arg("n_f") = 4.0);
  m.def("solve_A_interp_a", &solve_A_interp_a, py::arg("f1_func"),
        py::arg("q2_ref"), py::arg("variant") = "mid_x",
        py::arg("c_moment") = C_BAG);
  m.def("solve_A_interp_b", &solve_A_interp_b, py::arg("variant") = "mid_x",
        py::arg("c_moment") = C_BAG);

  py::class_<DeltaModel, TensorSF, std::shared_ptr<DeltaModel>>(m, "DeltaModel")
      .def(py::init<std::string, SFFunc3, std::string>(), py::arg("name"),
           py::arg("func"), py::arg("info") = "")
      .def("__call__", &DeltaModel::operator(), py::arg("x"), py::arg("q2"),
           py::arg("f1"))
      .def_property_readonly("name", &DeltaModel::name)
      .def_property_readonly("info", &DeltaModel::info);
  m.def("make_delta_toy", &make_delta_toy, py::arg("scale") = 1e-3);
  m.def("make_moment_a", &make_moment_a, py::arg("f1_func"), py::arg("q2_ref"),
        py::arg("variant") = "mid_x", py::arg("c_moment") = C_BAG,
        py::arg("dilution") = 1.0, py::arg("alphas") = nullptr);
  m.def("make_moment_b", &make_moment_b, py::arg("variant") = "mid_x",
        py::arg("c_moment") = C_BAG, py::arg("dilution") = 1.0,
        py::arg("alphas") = nullptr);
}

// -------------------------------------------------------- xsec, asymmetries

static void bind_xsec(py::module_& m) {
  // asymmetries.hpp
  m.def("depolarization_d", &depolarization_d, py::arg("y"), py::arg("x"),
        py::arg("q2"), py::arg("r_func") = nullptr);
  // `a_parallel(..., g2=None)`: the massless limit with no g2, exactly
  // `a_parallel_exact` with one -- the Python's own split.
  m.def("a_parallel", [](double g1, double f1, double y, double x, double q2,
                         const RFunc& r_func, std::optional<double> g2) {
    return g2 ? a_parallel(g1, f1, y, x, q2, r_func, *g2)
              : a_parallel(g1, f1, y, x, q2, r_func);
  }, py::arg("g1"), py::arg("f1"), py::arg("y"), py::arg("x"), py::arg("q2"),
     py::arg("r_func") = nullptr, py::arg("g2") = py::none());
  m.def("a_parallel_exact", &a_parallel_exact, py::arg("g1"), py::arg("g2"),
        py::arg("f1"), py::arg("y"), py::arg("x"), py::arg("q2"),
        py::arg("r_func") = nullptr);
  m.def("depolarization_effective",
        [](double y, double x, double q2, std::optional<double> g2_over_g1,
           const RFunc& r_func) {
          return depolarization_effective(y, x, q2,
                                          g2_over_g1 ? *g2_over_g1 : 0.0,
                                          static_cast<bool>(g2_over_g1),
                                          r_func);
        }, py::arg("y"), py::arg("x"), py::arg("q2"),
           py::arg("g2_over_g1") = py::none(), py::arg("r_func") = nullptr);
  m.def("a_perp", &a_perp, py::arg("g1"), py::arg("g2"), py::arg("f1"),
        py::arg("y"), py::arg("x"), py::arg("q2"), py::arg("r_func") = nullptr);
  m.def("phi_averaged_density", &phi_averaged_density, py::arg("f1"),
        py::arg("f2"), py::arg("x"), py::arg("y"));
  m.def("azz", [](double b1, double f1, double f2, double x, double y,
                  std::optional<double> b2, double theta_m) {
    const double v = b2 ? *b2 : 0.0;
    return azz(b1, f1, f2, x, y, b2 ? &v : nullptr, theta_m);
  }, py::arg("b1"), py::arg("f1"), py::arg("f2"), py::arg("x"), py::arg("y"),
     py::arg("b2") = py::none(), py::arg("theta_m") = 0.0);
  m.def("a_cos2phi", &a_cos2phi, py::arg("delta"), py::arg("f1"), py::arg("f2"),
        py::arg("x"), py::arg("y"));
  m.def("err_a_parallel", &err_a_parallel, py::arg("n"), py::arg("pe"),
        py::arg("pz"));
  m.def("err_azz", &err_azz, py::arg("n"), py::arg("pzz"));
  m.def("err_cos2phi_amplitude", &err_cos2phi_amplitude, py::arg("n"),
        py::arg("pzz"));

  // the finite-gamma kinematics -- asymmetries.hpp since 2026-08-29, exposed
  // under both names exactly as the Python re-exports them from xsec
  m.def("gamma_squared", &gamma_squared, py::arg("x"), py::arg("q2"),
        py::arg("m") = M_NUCLEON);
  m.def("epsilon_gamma", &epsilon_gamma, py::arg("y"), py::arg("gamma2"));
  m.def("depolarization_d_gamma", &depolarization_d_gamma, py::arg("y"),
        py::arg("gamma2"), py::arg("r"));
  m.def("depolarization_gamma", &depolarization_gamma, py::arg("y"),
        py::arg("gamma2"), py::arg("r"));
  m.def("eta_gamma", &eta_gamma, py::arg("y"), py::arg("gamma2"));

  // xsec.hpp -- the exact finite-gamma tensor sector (Cosyn, plans/08 D2)
  m.def("theta_q_cos_sin", &theta_q_cos_sin, py::arg("y"), py::arg("gamma2"),
        "(cos theta_q, sin theta_q): Cosyn Eq. (24)");
  m.def("cosyn_tensor_sfs", [](double b1, double b2, double b3, double b4,
                               double x, double gamma2) {
    const CosynTensorSFs f = cosyn_tensor_sfs(b1, b2, b3, b4, x, gamma2);
    return py::make_tuple(f.f_t, f.f_l, f.f_lt, f.f_tt);
  }, py::arg("b1"), py::arg("b2"), py::arg("b3"), py::arg("b4"), py::arg("x"),
     py::arg("gamma2"),
     "(F_TLL_T, F_TLL_L, F_TLT, F_TTT): Cosyn Eqs. (17a)-(17e)");
  m.def("cosyn_unpolarized_sfs", [](double f1, double f2, double x,
                                    double gamma2) {
    const std::pair<double, double> u = cosyn_unpolarized_sfs(f1, f2, x, gamma2);
    return py::make_tuple(u.first, u.second);
  }, py::arg("f1"), py::arg("f2"), py::arg("x"), py::arg("gamma2"),
     "(F_UU_T, F_UU_L): Cosyn Eq. (16)");
  m.def("density_min", &density_min, py::arg("a1n"), py::arg("a2n"));

  py::class_<EventSpinState>(m, "EventSpinState")
      .def(py::init([](int lam_e, double pe, double j, double mm,
                       double theta_s, double phi_s) {
        EventSpinState s;
        s.lam_e = lam_e;
        s.pe = pe;
        s.j = j;
        s.m = mm;
        s.theta_s = theta_s;
        s.phi_s = phi_s;
        return s;
      }), py::arg("lam_e") = 0, py::arg("pe") = 0.0, py::arg("j") = 1.0,
          py::arg("m") = 0.0, py::arg("theta_s") = 0.0, py::arg("phi_s") = 0.0)
      .def_readwrite("lam_e", &EventSpinState::lam_e)
      .def_readwrite("pe", &EventSpinState::pe)
      .def_readwrite("j", &EventSpinState::j)
      .def_readwrite("m", &EventSpinState::m)
      .def_readwrite("theta_s", &EventSpinState::theta_s)
      .def_readwrite("phi_s", &EventSpinState::phi_s);

  py::class_<SFTables>(m, "SFTables")
      .def(py::init<>())
      .def_readwrite("f1", &SFTables::f1)
      .def_readwrite("f2", &SFTables::f2)
      .def_readwrite("g1", &SFTables::g1)
      .def_readwrite("b1", &SFTables::b1)
      .def_readwrite("b2", &SFTables::b2)
      .def_readwrite("b3", &SFTables::b3)
      .def_readwrite("b4", &SFTables::b4)
      .def_readwrite("delta", &SFTables::delta)
      .def_readwrite("g2", &SFTables::g2)
      .def_readwrite("has_g2", &SFTables::has_g2)
      .def("as_dict", [](const SFTables& t) {
        py::dict d;
        d["f1"] = t.f1; d["f2"] = t.f2; d["g1"] = t.g1; d["g2"] = t.g2;
        d["b1"] = t.b1; d["b2"] = t.b2; d["b3"] = t.b3; d["b4"] = t.b4;
        d["delta"] = t.delta;
        return d;
      });

  py::class_<TensorHarmonics>(m, "TensorHarmonics")
      .def_readonly("h0", &TensorHarmonics::h0)
      .def_readonly("h1", &TensorHarmonics::h1)
      .def_readonly("h2", &TensorHarmonics::h2)
      .def("__iter__", [](const TensorHarmonics& h) {
        return py::iter(py::make_tuple(h.h0, h.h1, h.h2));
      });

  py::class_<Amplitudes>(m, "Amplitudes")
      .def_readonly("w_avg", &Amplitudes::w_avg)
      .def_readonly("a1", &Amplitudes::a1)
      .def_readonly("a2", &Amplitudes::a2)
      .def("__iter__", [](const Amplitudes& a) {
        return py::iter(py::make_tuple(a.w_avg, a.a1, a.a2));
      })
      .def("__repr__", [](const Amplitudes& a) {
        char b[160];
        std::snprintf(b, sizeof b,
                      "Amplitudes(w_avg=%.12g, a1=%.12g, a2=%.12g)",
                      a.w_avg, a.a1, a.a2);
        return std::string(b);
      });

  py::enum_<G2Mode>(m, "G2Mode")
      .value("WandzuraWilczek", G2Mode::kWandzuraWilczek)
      .value("Zero", G2Mode::kZero);

  py::class_<InclusiveKernel, std::shared_ptr<InclusiveKernel>> kern(
      m, "InclusiveKernel");
  py::class_<InclusiveKernel::Options>(kern, "Options")
      .def(py::init<>())
      .def_property("f2_source",
          [](const InclusiveKernel::Options& o) {
            return std::const_pointer_cast<UnpolSF>(o.f2_source);
          },
          [](InclusiveKernel::Options& o, std::shared_ptr<UnpolSF> v) {
            o.f2_source = std::move(v);
          })
      .def_property("g1_model",
          [](const InclusiveKernel::Options& o) {
            return std::const_pointer_cast<PolSF>(o.g1_model);
          },
          [](InclusiveKernel::Options& o, std::shared_ptr<PolSF> v) {
            o.g1_model = std::move(v);
          })
      .def_readwrite("b1_func", &InclusiveKernel::Options::b1_func)
      .def_readwrite("b2_func", &InclusiveKernel::Options::b2_func)
      .def_readwrite("delta_func", &InclusiveKernel::Options::delta_func)
      .def_readwrite("b1_32_func", &InclusiveKernel::Options::b1_32_func)
      .def_readwrite("b2_32_func", &InclusiveKernel::Options::b2_32_func)
      .def_readwrite("delta_32_func", &InclusiveKernel::Options::delta_32_func)
      .def_readwrite("b3_func", &InclusiveKernel::Options::b3_func)
      .def_readwrite("b4_func", &InclusiveKernel::Options::b4_func)
      .def_readwrite("tensor_gamma", &InclusiveKernel::Options::tensor_gamma)
      .def_readwrite("g2_mode", &InclusiveKernel::Options::g2_mode)
      .def_readwrite("g2_scale", &InclusiveKernel::Options::g2_scale)
      .def_readwrite("emc_ratio", &InclusiveKernel::Options::emc_ratio)
      .def_readwrite("r_func", &InclusiveKernel::Options::r_func)
      .def_readwrite("target_mass", &InclusiveKernel::Options::target_mass)
      .def_readwrite("g2_npts", &InclusiveKernel::Options::g2_npts);

  kern.def(py::init([](Ion ion) {
             return std::make_shared<InclusiveKernel>(std::move(ion));
           }), py::arg("ion"))
      .def(py::init([](Ion ion, InclusiveKernel::Options opt) {
             return std::make_shared<InclusiveKernel>(std::move(ion),
                                                      std::move(opt));
           }), py::arg("ion"), py::arg("options"))
      .def_property_readonly("ion", &InclusiveKernel::ion)
      .def_property_readonly("target_mass", &InclusiveKernel::target_mass)
      .def_property_readonly("tensor_gamma", &InclusiveKernel::tensor_gamma)
      .def_property_readonly("g2_scale", &InclusiveKernel::g2_scale)
      .def("tables", &InclusiveKernel::tables, py::arg("x"), py::arg("q2"),
           py::arg("with_g2") = false)
      .def_static("dphi", &InclusiveKernel::dphi, py::arg("t"), py::arg("x"),
                  py::arg("y"))
      .def_static("tensor_kernel", &InclusiveKernel::tensor_kernel,
                  py::arg("t"), py::arg("x"), py::arg("y"))
      .def("a_parallel", &InclusiveKernel::a_parallel, py::arg("t"),
           py::arg("x"), py::arg("q2"), py::arg("y"))
      .def("a_perp", &InclusiveKernel::a_perp, py::arg("t"), py::arg("x"),
           py::arg("q2"), py::arg("y"))
      .def("tensor_moments", &InclusiveKernel::tensor_moments, py::arg("m"))
      .def("tensor_harmonics_gamma", &InclusiveKernel::tensor_harmonics_gamma,
           py::arg("t"), py::arg("x"), py::arg("q2"), py::arg("y"),
           py::arg("state"),
           "(h0, h1, h2) of the exact finite-gamma b-sector")
      .def("amplitudes", &InclusiveKernel::amplitudes, py::arg("tables"),
           py::arg("x"), py::arg("q2"), py::arg("s"), py::arg("state"),
           py::arg("with_perp") = false)
      .def("amplitudes_at", [](const InclusiveKernel& k, double x, double q2,
                               double s, const EventSpinState& st,
                               bool with_perp) {
             return k.amplitudes(k.tables(x, q2, with_perp), x, q2, s, st,
                                 with_perp);
           }, py::arg("x"), py::arg("q2"), py::arg("s"), py::arg("state"),
              py::arg("with_perp") = false,
              "amplitudes() with the SF tables computed on the spot")
      .def("dsigma_unpol", &InclusiveKernel::dsigma_unpol, py::arg("x"),
           py::arg("q2"), py::arg("s"))
      .def("dsigma", &InclusiveKernel::dsigma, py::arg("x"), py::arg("q2"),
           py::arg("phi"), py::arg("s"), py::arg("state"),
           py::arg("with_perp") = false)
      .def_static("density", &InclusiveKernel::density, py::arg("a"),
                  py::arg("phip"))
      .def_static("positivity_margin", &InclusiveKernel::positivity_margin,
                  py::arg("a"))
      .def_property_readonly("nuclear_f2", &InclusiveKernel::nuclear_f2,
                             py::return_value_policy::reference_internal);
}

// ------------------------------------------------------------- bookkeeping

static void bind_bookkeeping(py::module_& m) {
  py::class_<SpinCategory>(m, "SpinCategory")
      .def(py::init<std::string, double, std::vector<double>, int, double,
                    double, double, double>(),
           py::arg("name"), py::arg("j"), py::arg("populations"),
           py::arg("lam_e") = 0, py::arg("pe") = 0.0, py::arg("theta_s") = 0.0,
           py::arg("phi_s") = 0.0, py::arg("lumi_fraction") = 1.0)
      .def_readwrite("name", &SpinCategory::name)
      .def_readwrite("j", &SpinCategory::j)
      .def_readwrite("populations", &SpinCategory::populations)
      .def_readwrite("lam_e", &SpinCategory::lam_e)
      .def_readwrite("pe", &SpinCategory::pe)
      .def_readwrite("theta_s", &SpinCategory::theta_s)
      .def_readwrite("phi_s", &SpinCategory::phi_s)
      .def_readwrite("lumi_fraction", &SpinCategory::lumi_fraction)
      .def("moments", &SpinCategory::moments)
      .def("vector_moment", &SpinCategory::vector_moment)
      .def("__repr__", [](const SpinCategory& c) {
        return "SpinCategory('" + c.name + "', j=" + std::to_string(c.j) + ")";
      });

  py::class_<RunPlan>(m, "RunPlan")
      .def(py::init<>())
      .def(py::init<std::vector<SpinCategory>, double, double, double, double,
                    std::uint64_t>(),
           py::arg("categories"), py::arg("pe_true") = 0.0,
           py::arg("pz_true") = 0.0, py::arg("pzz_true") = 0.0,
           py::arg("delta_p_over_p") = 0.0,
           py::arg("polarimetry_seed") = 20260713)
      .def_property_readonly("categories", &RunPlan::categories)
      .def_property_readonly("pe_true", &RunPlan::pe_true)
      .def_property_readonly("pz_true", &RunPlan::pz_true)
      .def_property_readonly("pzz_true", &RunPlan::pzz_true)
      .def_property_readonly("delta_p_over_p", &RunPlan::delta_p_over_p)
      .def_property_readonly("measured_pe", &RunPlan::measured_pe)
      .def_property_readonly("measured_pz", &RunPlan::measured_pz)
      .def_property_readonly("measured_pzz", &RunPlan::measured_pzz)
      .def("lumi_shares", &RunPlan::lumi_shares, py::arg("total_lumi_pb"))
      .def("lumi_share_vector", &RunPlan::lumi_share_vector,
           py::arg("total_lumi_pb"))
      .def("index_of", &RunPlan::index_of, py::arg("name"))
      .def("__len__", [](const RunPlan& p) { return p.categories().size(); });

  py::class_<HelicityFlipOptions>(m, "HelicityFlipOptions")
      .def(py::init<>())
      // `--pzz-mode` IS this bool: false = "ladder" (the default max-entropy
      // fill at pz), true = "typed" (the caller's pzz through
      // spin1_populations / spin32_populations).  `lipolgen.make_plan` is the
      // one place the mode string becomes this flag.
      .def_readwrite("use_explicit_pzz", &HelicityFlipOptions::use_explicit_pzz)
      .def_readwrite("pzz", &HelicityFlipOptions::pzz)
      .def_readwrite("theta_s", &HelicityFlipOptions::theta_s)
      .def_readwrite("phi_s", &HelicityFlipOptions::phi_s)
      .def_readwrite("rel_lumi_offset", &HelicityFlipOptions::rel_lumi_offset)
      .def_readwrite("name", &HelicityFlipOptions::name);

  m.def("helicity_flip_plan", [](double j, double pz, double pe,
                                 const HelicityFlipOptions& o) {
    return helicity_flip_plan(j, pz, pe, o);
  }, py::arg("j"), py::arg("pz"), py::arg("pe"),
     py::arg("options") = HelicityFlipOptions());
  m.def("tensor_thirds_plan", &tensor_thirds_plan, py::arg("pz"),
        py::arg("pzz"), py::arg("rel_lumi_offset") = 0.0,
        py::arg("theta_s") = 0.0, py::arg("phi_s") = 0.0,
        py::arg("name") = "azz");
  m.def("transverse_tensor_plan", &transverse_tensor_plan, py::arg("pzz"),
        py::arg("phi_s") = 0.0, py::arg("name") = "cos2phi");
  m.def("tensor_flip_plan", &tensor_flip_plan, py::arg("pzz"),
        py::arg("phi_s") = kPi / 2.0, py::arg("share_plus") = 0.5,
        py::arg("rel_lumi_offset") = 0.0, py::arg("name") = "flip");
  m.def("with_offset", &with_offset, py::arg("plan"), py::arg("category_name"),
        py::arg("offset"));
  m.def("azz_rel_lumi_bias", &azz_rel_lumi_bias, py::arg("offset"),
        py::arg("pzz"));
  m.def("apar_rel_lumi_bias", &apar_rel_lumi_bias, py::arg("offset"),
        py::arg("pe"), py::arg("pz"));

  py::class_<SpinTemperatureLadder>(m, "SpinTemperatureLadder")
      .def_readonly("populations", &SpinTemperatureLadder::populations)
      .def_readonly("t", &SpinTemperatureLadder::t);
  m.def("spin_temperature_ladder", &spin_temperature_ladder, py::arg("j"),
        py::arg("pz"), py::arg("iterations") = 400);
  m.def("spin_temperature_pzz", &spin_temperature_pzz, py::arg("j"),
        py::arg("pz"));
  m.def("bunch_rng", &bunch_rng, py::arg("seed"), py::arg("run"),
        py::arg("bunch"), py::arg("event") = 0);
}

// ----------------------------------------------------------------- sampler

static void bind_sampler(py::module_& m) {
  m.def("y_from_xq2", &y_from_xq2, py::arg("x"), py::arg("q2"), py::arg("s"));
  m.def("w2_from_xq2", &w2_from_xq2, py::arg("x"), py::arg("q2"),
        py::arg("m_n") = M_NUCLEON);

  py::class_<ScatteredElectron>(m, "ScatteredElectronKin")
      .def_readonly("e_prime", &ScatteredElectron::e_prime)
      .def_readonly("theta", &ScatteredElectron::theta)
      .def_readonly("eta", &ScatteredElectron::eta);
  m.def("scattered_electron", &scattered_electron, py::arg("x"), py::arg("y"),
        py::arg("s"), py::arg("electron_energy"));

  py::class_<LogGrid>(m, "LogGrid")
      .def_readonly("x_edges", &LogGrid::x_edges)
      .def_readonly("q2_edges", &LogGrid::q2_edges)
      .def_readonly("x_c", &LogGrid::x_c)
      .def_readonly("q2_c", &LogGrid::q2_c);
  m.def("log_grid", &log_grid, py::arg("x_min"), py::arg("x_max"),
        py::arg("q2_min"), py::arg("q2_max"), py::arg("nx"), py::arg("nq2"));

  py::class_<Scenario>(m, "Scenario")
      .def(py::init<>())
      .def_readwrite("lumi_fb_per_nucleon", &Scenario::lumi_fb_per_nucleon)
      .def_readwrite("run_share", &Scenario::run_share)
      .def_readwrite("pol_electron", &Scenario::pol_electron)
      .def_readwrite("pol_ion_vector", &Scenario::pol_ion_vector)
      .def_readwrite("pol_ion_tensor", &Scenario::pol_ion_tensor)
      .def_readwrite("q2_min", &Scenario::q2_min)
      .def_readwrite("y_min", &Scenario::y_min)
      .def_readwrite("y_max", &Scenario::y_max)
      .def_readwrite("w2_min", &Scenario::w2_min)
      .def_readwrite("x_max", &Scenario::x_max)
      .def_readwrite("eta_min", &Scenario::eta_min)
      .def_readwrite("eta_max", &Scenario::eta_max)
      .def_readwrite("e_prime_min", &Scenario::e_prime_min)
      .def("lumi_effective_fb_per_nucleon",
           &Scenario::lumi_effective_fb_per_nucleon)
      .def("lumi_effective_pb_per_nucleon",
           &Scenario::lumi_effective_pb_per_nucleon)
      .def("validate", &Scenario::validate);
  m.def("generator_scenario", &generator_scenario, py::arg("analysis"),
        py::arg("q2_min") = 0.7, py::arg("y_min") = 0.004,
        py::arg("y_max") = 0.985, py::arg("w2_min") = 8.0,
        py::arg("eta_pad") = 0.3, py::arg("e_prime_min") = 0.3);
  m.def("rng_poisson", &rng_poisson, py::arg("rng"), py::arg("mean"));

  py::class_<SamplerGrid>(m, "SamplerGrid")
      .def(py::init<>())
      .def_readwrite("nx", &SamplerGrid::nx)
      .def_readwrite("nq2", &SamplerGrid::nq2)
      .def_readwrite("x_min", &SamplerGrid::x_min)
      .def_readwrite("x_max", &SamplerGrid::x_max)
      .def_readwrite("q2_min", &SamplerGrid::q2_min)
      .def_readwrite("q2_max", &SamplerGrid::q2_max);

  py::class_<EventDraw>(m, "EventDraw")
      .def_readonly("x", &EventDraw::x)
      .def_readonly("q2", &EventDraw::q2)
      .def_readonly("y", &EventDraw::y)
      .def_readonly("phi", &EventDraw::phi)
      .def_readonly("m", &EventDraw::m)
      .def_readonly("cell", &EventDraw::cell);

  py::class_<InclusiveSampler::EffectiveModulation>(m, "EffectiveModulation")
      .def_readonly("sigma_pb", &InclusiveSampler::EffectiveModulation::sigma_pb)
      .def_readonly("a1", &InclusiveSampler::EffectiveModulation::a1)
      .def_readonly("a2", &InclusiveSampler::EffectiveModulation::a2);

  py::class_<InclusiveSampler>(m, "InclusiveSampler")
      .def(py::init([](std::shared_ptr<InclusiveKernel> kern, BeamConfig cfg,
                       Scenario sc, SamplerGrid g, bool wp) {
        return std::make_unique<InclusiveSampler>(std::move(kern),
                                                  std::move(cfg), std::move(sc),
                                                  g, wp);
      }), py::arg("kernel"), py::arg("config"), py::arg("scenario") = Scenario(),
          py::arg("grid") = SamplerGrid(), py::arg("with_perp") = false)
      .def_property_readonly("kernel", [](const InclusiveSampler& s) {
        return std::const_pointer_cast<InclusiveKernel>(s.kernel_ptr());
      })
      .def_property_readonly("config", &InclusiveSampler::config)
      .def_property_readonly("scenario", &InclusiveSampler::scenario)
      .def_property_readonly("s", &InclusiveSampler::s)
      .def_property_readonly("n_cells", &InclusiveSampler::n_cells)
      .def_property_readonly("x_cells", [](const InclusiveSampler& s) {
        return copy_array(s.x_cells());
      })
      .def_property_readonly("q2_cells", [](const InclusiveSampler& s) {
        return copy_array(s.q2_cells());
      })
      .def_property_readonly("cell_xsec_pb", [](const InclusiveSampler& s) {
        return copy_array(s.cell_xsec_pb());
      })
      .def_property_readonly("logx_lo", [](const InclusiveSampler& s) {
        return copy_array(s.logx_lo());
      })
      .def_property_readonly("logx_hi", [](const InclusiveSampler& s) {
        return copy_array(s.logx_hi());
      })
      .def_property_readonly("logq2_lo", [](const InclusiveSampler& s) {
        return copy_array(s.logq2_lo());
      })
      .def_property_readonly("logq2_hi", [](const InclusiveSampler& s) {
        return copy_array(s.logq2_hi());
      })
      .def("in_acceptance", &InclusiveSampler::in_acceptance, py::arg("x"),
           py::arg("q2"))
      .def("state_tables", [](const InclusiveSampler& s,
                              const SpinCategory& cat, double mm) {
        const InclusiveSampler::StateTables& t = s.state_tables(cat, mm);
        py::dict d;
        d["w_avg"] = copy_array(t.w_avg);
        d["a1"] = copy_array(t.a1);
        d["a2"] = copy_array(t.a2);
        d["a1n"] = copy_array(t.a1n);
        d["a2n"] = copy_array(t.a2n);
        d["bound"] = copy_array(t.bound);
        d["cdf"] = copy_array(t.cdf);
        d["sigma_pb"] = t.sigma_pb;
        d["margin"] = t.margin;
        return d;
      }, py::arg("category"), py::arg("m"))
      .def("sigma_state_pb", &InclusiveSampler::sigma_state_pb,
           py::arg("category"), py::arg("m"))
      .def("sigma_tot_pb", &InclusiveSampler::sigma_tot_pb, py::arg("category"))
      .def("expected_events", &InclusiveSampler::expected_events,
           py::arg("category"), py::arg("lumi_pb"))
      .def("effective_modulation", [](const InclusiveSampler& s,
                                      const SpinCategory& cat) {
        return s.effective_modulation(cat, nullptr);
      }, py::arg("category"))
      .def("sample_n", [](const InclusiveSampler& s, const SpinCategory& cat,
                          std::size_t n, std::uint64_t seed, std::uint64_t run,
                          std::uint64_t bunch, std::uint64_t event0,
                          unsigned nthreads) {
        EventBatch b;
        {
          py::gil_scoped_release unlock;
          b = s.sample_n(cat, n, seed, run, bunch, event0, nthreads);
        }
        return batch_to_dict(std::move(b));
      }, py::arg("category"), py::arg("n"), py::arg("seed"), py::arg("run") = 1,
         py::arg("bunch") = 0, py::arg("event0") = 0, py::arg("nthreads") = 1,
         "dict of numpy arrays: x, q2, y, phi, m, weight, cell + category,"
         " lam_e (polligen `sample_category` schema)")
      .def("sample_lumi", [](const InclusiveSampler& s, const SpinCategory& cat,
                             double lumi_pb, std::uint64_t seed,
                             std::uint64_t run, std::uint64_t bunch,
                             bool poisson, unsigned nthreads) {
        EventBatch b;
        {
          py::gil_scoped_release unlock;
          b = s.sample_lumi(cat, lumi_pb, seed, run, bunch, poisson, nthreads);
        }
        return batch_to_dict(std::move(b));
      }, py::arg("category"), py::arg("lumi_pb"), py::arg("seed"),
         py::arg("run") = 1, py::arg("bunch") = 0, py::arg("poisson") = true,
         py::arg("nthreads") = 1)
      .def("weights_for", [](const InclusiveSampler& s, const py::dict& events,
                             const std::vector<SpinCategory>& cats) {
        EventBatch b = dict_to_batch(events);
        std::vector<double> w = s.weights_for(b, cats);
        return move_array2(std::move(w),
                           static_cast<py::ssize_t>(cats.size()));
      }, py::arg("events"), py::arg("categories"),
         "Mode-W weight matrix w[i, k] from an events dict (uses cell, phi)");

  py::class_<PseudoExperiment>(m, "PseudoExperiment")
      .def_readonly("names", &PseudoExperiment::names)
      .def_readonly("lumi_pb", &PseudoExperiment::lumi_pb)
      .def("count", &PseudoExperiment::count, py::arg("name"))
      .def("lumi_of", &PseudoExperiment::lumi_of, py::arg("name"))
      .def("batch", [](const PseudoExperiment& p, const std::string& n) {
        EventBatch copy = p[n];
        return batch_to_dict(std::move(copy));
      }, py::arg("name"))
      .def("as_dict", [](const PseudoExperiment& p) {
        py::dict d;
        for (std::size_t i = 0; i < p.names.size(); ++i) {
          EventBatch copy = p.batches[i];
          d[py::str(p.names[i])] = batch_to_dict(std::move(copy));
        }
        return d;
      });
  m.def("run_pseudo_experiment", [](const InclusiveSampler& s,
                                    const RunPlan& plan, double lumi_pb,
                                    std::uint64_t seed, std::uint64_t run,
                                    bool poisson, unsigned nthreads) {
    py::gil_scoped_release unlock;
    return run_pseudo_experiment(s, plan, lumi_pb, seed, run, poisson, nthreads);
  }, py::arg("sampler"), py::arg("plan"), py::arg("total_lumi_pb"),
     py::arg("seed"), py::arg("run") = 1, py::arg("poisson") = true,
     py::arg("nthreads") = 1);

  py::class_<PhiHistogram>(m, "PhiHistogram")
      .def_readonly("counts", &PhiHistogram::counts)
      .def_readonly("edges", &PhiHistogram::edges);
  m.def("phi_histogram_pseudo", [](double n_expected, double a2, int nbins,
                                   Rng* rng, double a1, bool poisson) {
    return phi_histogram_pseudo(n_expected, a2, nbins, rng, a1, poisson);
  }, py::arg("n_expected"), py::arg("a2"), py::arg("nbins") = 36,
     py::arg("rng") = nullptr, py::arg("a1") = 0.0, py::arg("poisson") = true);

  py::module_ est = m.def_submodule("estimators",
      "Analysis-side counting / fit estimators (sampler.hpp)");
  est.def("yields", &estimators::yields, py::arg("counts"),
          py::arg("lumis") = std::vector<double>());
  est.def("apar_flip", &estimators::apar_flip, py::arg("n_plus"),
          py::arg("n_minus"), py::arg("pe"), py::arg("pz"),
          py::arg("l_plus") = 0.0, py::arg("l_minus") = 0.0);
  est.def("azz_thirds", &estimators::azz_thirds, py::arg("n_plus"),
          py::arg("n_minus"), py::arg("n_zero"), py::arg("pzz"),
          py::arg("lumis") = std::vector<double>());
  est.def("cos2phi_moment", &estimators::cos2phi_moment, py::arg("phi_prime"),
          py::arg("pzz"));
  est.def("cos2phi_fit_binned", &estimators::cos2phi_fit_binned,
          py::arg("counts"), py::arg("edges"), py::arg("pzz"),
          py::arg("acceptance") = nullptr, py::arg("nsub") = 9);
  est.def("cos2phi_fit", &estimators::cos2phi_fit, py::arg("phi_prime"),
          py::arg("pzz"), py::arg("nbins") = 36,
          py::arg("acceptance") = nullptr, py::arg("nsub") = 9);
  est.def("cos2phi_fit_err", &estimators::cos2phi_fit_err, py::arg("n"),
          py::arg("pzz"), py::arg("nbins") = 24);
  est.def("pull", &estimators::pull, py::arg("estimate"), py::arg("truth"),
          py::arg("expected_err"));

  m.def("nuclear_pdg", &nuclear_pdg, py::arg("z"), py::arg("a"));
}

// --------------------------------------------------------------- spectator

static void bind_spectator(py::module_& m) {
  m.def("cluster_mass", &cluster_mass, py::arg("name"));
  m.def("nuclear_mass", &nuclear_mass, py::arg("z"), py::arg("a"));
  m.def("nuclear_mass_known", &nuclear_mass_known, py::arg("z"), py::arg("a"));

  py::class_<ClusterChannel>(m, "ClusterChannel")
      .def(py::init<>())
      .def_readwrite("name", &ClusterChannel::name)
      .def_readwrite("beam_A", &ClusterChannel::beam_A)
      .def_readwrite("beam_Z", &ClusterChannel::beam_Z)
      .def_readwrite("spectator", &ClusterChannel::spectator)
      .def_readwrite("spectator_A", &ClusterChannel::spectator_A)
      .def_readwrite("spectator_Z", &ClusterChannel::spectator_Z)
      .def_readwrite("separation_energy", &ClusterChannel::separation_energy)
      .def_readwrite("l_wave", &ClusterChannel::l_wave)
      .def("m_spec", &ClusterChannel::m_spec)
      .def("m_beam", &ClusterChannel::m_beam)
      .def("m_partner", &ClusterChannel::m_partner)
      .def("partner_A", &ClusterChannel::partner_A)
      .def("partner_Z", &ClusterChannel::partner_Z)
      .def("kappa", &ClusterChannel::kappa)
      .def("r_at_k_zero", &ClusterChannel::r_at_k_zero);
  m.def("channel_by_name", [](const std::string& n) {
    return channel_by_name(n);
  }, py::arg("name"));

  m.def("momentum_density", &momentum_density, py::arg("k"), py::arg("kappa"),
        py::arg("beta"), py::arg("l_wave"));
  py::class_<MomentumSampler>(m, "MomentumSampler")
      .def(py::init<const ClusterChannel&, double, double, std::size_t>(),
           py::arg("channel"), py::arg("beta") = BETA_DEFAULT,
           py::arg("k_max") = 1.5, py::arg("n_grid") = 30000)
      .def("k_of_u", &MomentumSampler::k_of_u, py::arg("u"))
      .def("grid", [](const MomentumSampler& s) { return copy_array(s.grid()); })
      .def("cdf", [](const MomentumSampler& s) { return copy_array(s.cdf()); });

  py::class_<FragmentLab>(m, "FragmentLab")
      .def_readonly("pT", &FragmentLab::pT)
      .def_readonly("theta", &FragmentLab::theta)
      .def_readonly("phi", &FragmentLab::phi)
      .def_readonly("p_lab", &FragmentLab::p_lab)
      .def_readonly("R", &FragmentLab::R)
      .def_readonly("xL", &FragmentLab::xL)
      .def_readonly("k", &FragmentLab::k)
      .def_readonly("e_lab", &FragmentLab::e_lab)
      .def_readonly("pz_lab", &FragmentLab::pz_lab);
  m.def("boost_fragment", &boost_fragment, py::arg("channel"),
        py::arg("p_per_nucleon"), py::arg("kx"), py::arg("ky"), py::arg("kz"),
        py::arg("m"), py::arg("frag_Z"), py::arg("frag_A"));
  m.def("boost_spectator_fragment", &boost_spectator_fragment,
        py::arg("channel"), py::arg("p_per_nucleon"), py::arg("kx"),
        py::arg("ky"), py::arg("kz"));

  py::class_<Optics>(m, "Optics")
      .def(py::init<>())
      .def_readwrite("name", &Optics::name)
      .def_readwrite("sigma_theta", &Optics::sigma_theta)
      .def_readwrite("n_sigma", &Optics::n_sigma)
      .def_readwrite("sigma_theta_v", &Optics::sigma_theta_v)
      .def_readwrite("lumi_fraction", &Optics::lumi_fraction)
      .def("sigma_v", &Optics::sigma_v)
      .def("isotropic", &Optics::isotropic)
      .def("envelope_x", &Optics::envelope_x)
      .def("envelope_y", &Optics::envelope_y)
      .def("clears", &Optics::clears, py::arg("theta"), py::arg("phi"))
      .def("pt_cut_for", &Optics::pt_cut_for, py::arg("momentum"))
      .def("__repr__", [](const Optics& o) {
        char b[160];
        std::snprintf(b, sizeof b,
                      "Optics('%s', sigma_h=%.4g, sigma_v=%.4g, n=%.3g)",
                      o.name.c_str(), o.sigma_theta, o.sigma_v(), o.n_sigma);
        return std::string(b);
      });
  m.def("high_acceptance_optics", []() { return HIGH_ACCEPTANCE(); });
  m.def("high_divergence_optics", []() { return HIGH_DIVERGENCE(); });

  py::enum_<Route>(m, "Route")
      .value("Lost", kRouteLost)
      .value("RomanPots", kRouteRomanPots)
      .value("OMD", kRouteOMD)
      .value("B0", kRouteB0)
      .value("RPNearBeam", kRouteRPNearBeam)
      .value("ZDC", kRouteZDC)
      .value("RPInner", kRouteRPInner);
  m.def("rp_accepted", [](int route) { return rp_accepted(route); },
        py::arg("route"));

  py::class_<PotLevers>(m, "PotLevers")
      .def_readonly("r12", &PotLevers::r12)
      .def_readonly("r34", &PotLevers::r34)
      .def_readonly("dispersion", &PotLevers::dispersion)
      .def_readonly("dispersion2", &PotLevers::dispersion2)
      .def_readonly("blind_half_width", &PotLevers::blind_half_width);
  m.def("pot_levers", [](const std::string& c) { return pot_levers(c); },
        py::arg("config"));
  m.def("over_rigid_route", &over_rigid_route, py::arg("r"),
        py::arg("theta_x") = 0.0, py::arg("config") = "18x275");
  m.def("route_charged", &route_charged, py::arg("r"), py::arg("theta"),
        py::arg("pT"), py::arg("optics"), py::arg("phi") = std::nan(""),
        py::arg("theta_outer") = std::nan(""),
        py::arg("pot_config") = "18x275");
  m.def("route_neutral", &route_neutral, py::arg("theta"));
  m.def("yr_config_key", &yr_config_key, py::arg("ion_name"),
        py::arg("p_per_nucleon"));
  m.def("sigma_theta_for", [](const std::string& n, double p, bool ha) {
    double h = 0, v = 0;
    sigma_theta_for(n, p, ha, h, v);
    return py::make_tuple(h, v);
  }, py::arg("ion_name"), py::arg("p_per_nucleon"),
     py::arg("high_acceptance") = true);
  m.def("yr_optics", &yr_optics, py::arg("ion_name"), py::arg("p_per_nucleon"),
        py::arg("high_acceptance") = true, py::arg("n_sigma") = 10.0);
}

// ------------------------------------------------------------------ tagged

static void bind_tagged(py::module_& m) {
  // ---- the VMC tabulated backend (cluster.hpp) --------------------------
  py::enum_<ClusterWaveSource>(m, "ClusterWaveSource",
      "Which family of cluster radial forms a channel is built from.  "
      "`Hulthen` (the default everywhere) keeps every published number "
      "bit-for-bit; `VmcAV18` swaps in the ANL VMC tables.")
      .value("Hulthen", ClusterWaveSource::Hulthen)
      .value("VmcAV18", ClusterWaveSource::VmcAV18);

  m.def("data_dir", []() { return data_dir(); },
        "$LIPOLGEN_DATA_DIR, else the compiled-in ${CMAKE_SOURCE_DIR}/data.");
  m.def("data_path", &data_path, py::arg("relative"));

  py::class_<VmcRadial, std::shared_ptr<VmcRadial>>(m, "VmcRadial",
      "Tabulated, L-specific psi_L(k) [GeV], linearly interpolated and ZERO "
      "outside the tabulated range.  Signed: the S-D interference term of "
      "n_M(k, khat) goes as psi_0 psi_2.")
      .def(py::init([](std::vector<double> k, std::vector<double> psi, int l,
                       std::string prov) {
        return VmcRadial(std::move(k), std::move(psi), l, std::move(prov));
      }), py::arg("k_gev"), py::arg("psi"), py::arg("l"),
          py::arg("provenance") = std::string())
      .def("__call__", [](const VmcRadial& v, double k) { return v(k); },
           py::arg("k"))
      .def_property_readonly("l", &VmcRadial::l)
      .def_property_readonly("k", [](const VmcRadial& v) {
        return copy_array(v.k());
      })
      .def_property_readonly("dpsi", [](const VmcRadial& v) {
        return copy_array(v.dpsi());
      })
      .def_property_readonly("has_errors", &VmcRadial::has_errors)
      .def("shifted_by_sigma", &VmcRadial::shifted_by_sigma,
           py::arg("n_sigma"),
           "psi -> psi + n_sigma * dpsi, fully correlated across the table.")
      .def("norm2_error", [](const VmcRadial& v) {
        double cor = 0.0, quad = 0.0;
        v.norm2_error(&cor, &quad);
        return py::make_tuple(cor, quad);
      }, "(correlated, quadrature) 1-sigma MC error of norm2() -- the two "
         "limits the ANL errors admit, since the point-to-point correlation "
         "of one variational walk is unknown.")
      .def_property_readonly("psi", [](const VmcRadial& v) {
        return copy_array(v.psi());
      })
      .def_property_readonly("provenance", &VmcRadial::provenance)
      .def("norm2", &VmcRadial::norm2);

  m.def("vmc_from_overlap_k", &vmc_from_overlap_k, py::arg("path"),
        py::arg("column"), py::arg("l"));
  m.def("vmc_from_overlap_r", &vmc_from_overlap_r, py::arg("path"),
        py::arg("column"), py::arg("l"), py::arg("k_max_fm") = 5.0,
        py::arg("nk") = 251);
  m.def("vmc_from_momentum", [](const std::string& path, int block, int column,
                                int l, const VmcRadial* sign_from,
                                double node_max) {
    return vmc_from_momentum(path, block, column, l, sign_from, node_max);
  }, py::arg("path"), py::arg("block"), py::arg("column"), py::arg("l"),
     py::arg("sign_from") = nullptr, py::arg("node_search_max_fm") = 3.0);
  m.def("read_anl_momentum_norms", &read_anl_momentum_norms, py::arg("path"));

  py::class_<Wave>(m, "Wave")
      .def(py::init([](int l, double prob, double beta,
                       std::shared_ptr<const VmcRadial> vmc) {
        Wave w;
        w.l = l;
        w.prob = prob;
        w.beta = beta;
        w.vmc = std::move(vmc);
        return w;
      }), py::arg("l") = 0, py::arg("prob") = 1.0,
          py::arg("beta") = BETA_DEFAULT,
          py::arg("vmc") = std::shared_ptr<const VmcRadial>())
      .def_readwrite("l", &Wave::l)
      .def_readwrite("prob", &Wave::prob)
      .def_readwrite("beta", &Wave::beta)
      .def_readwrite("vmc", &Wave::vmc)
      .def("radial", [](const Wave& w, double k, double kappa) {
        return w.radial(k, kappa);
      }, py::arg("k"), py::arg("kappa"));
  m.def("theta_lm", &theta_lm, py::arg("l"), py::arg("m"), py::arg("c"));

  m.def("triton", []() { return TRITON(); });
  m.def("neutron_target", []() { return NEUTRON_TARGET(); });

  py::class_<TaggedChannel>(m, "TaggedChannel")
      .def(py::init<>())
      .def_readwrite("base", &TaggedChannel::base)
      .def_readwrite("j_ion", &TaggedChannel::j_ion)
      .def_readwrite("s_struck", &TaggedChannel::s_struck)
      .def_readwrite("s_spec", &TaggedChannel::s_spec)
      .def_readwrite("s_channel", &TaggedChannel::s_channel)
      .def_readwrite("waves", &TaggedChannel::waves)
      .def_readwrite("dis_target", &TaggedChannel::dis_target)
      .def_readwrite("label", &TaggedChannel::label)
      .def("validate", &TaggedChannel::validate);
  m.def("li6_alpha_channel", &li6_alpha_channel,
        py::arg("beta") = BETA_DEFAULT, py::arg("p_d") = P_D_LI6,
        py::arg("source") = ClusterWaveSource::Hulthen,
        py::arg("vmc_mc_sigma") = 0.0,
        "`vmc_mc_sigma` is the ANL tables' own 1-sigma Monte Carlo band, "
        "fully correlated across k and the same on both waves; 0 is today's "
        "table bit for bit and it THROWS on Hulthen.");
  m.def("li7_alpha_channel", &li7_alpha_channel, py::arg("beta") = BETA_DEFAULT,
        py::arg("source") = ClusterWaveSource::Hulthen,
        py::arg("vmc_mc_sigma") = 0.0);
  m.attr("VMC_N_ALPHA_D_LI6") = VMC_N_ALPHA_D_LI6;
  m.attr("VMC_P_D_LI6") = VMC_P_D_LI6;
  m.attr("VMC_S_ALPHA_D_LI6") = VMC_S_ALPHA_D_LI6;
  m.attr("VMC_S_ALPHA_T_LI7") = VMC_S_ALPHA_T_LI7;
  m.attr("LI6_CLUSTER_POLARIZATION_VMC") = LI6_CLUSTER_POLARIZATION_VMC;
  m.def("deuteron_channel", &deuteron_channel, py::arg("beta") = BETA_DEFAULT,
        py::arg("p_d") = P_D_DEUTERON,
        py::arg("source") = ClusterWaveSource::Hulthen,
        "`source = VmcAV18` swaps the analytic Hulthen pair at the scenario "
        "P_D_DEUTERON = 0.045 for the exact AV18 deuteron of fdeut.av18 at "
        "its own P_D (`deuteron_av18_p_d()`).");
  m.def("deuteron_av18", []() { return DEUTERON_AV18(); },
        "The embedded deuteron of ClusterWaveSource.VmcAV18: beams.deuteron() "
        "with eff_pol = vector_dilution_of(deuteron_av18_p_d()) = 0.913600 "
        "instead of the scenario 0.9325.  li6_alpha_channel() puts this on "
        "dis_target under VmcAV18, so the flag means ONE deuteron everywhere "
        "it is read (C5.5b).");
  m.def("deuteron_av18_p_d", &deuteron_av18_p_d,
        "P_D of the AV18 deuteron from fdeut.av18's own k-space block "
        "(derived, never typed; the file's r-space header prints 0.057599).");

  py::class_<TaggedModel>(m, "TaggedModel")
      .def(py::init<TaggedChannel, double, std::size_t, std::size_t>(),
           py::arg("channel"), py::arg("k_max") = 1.2, py::arg("nk") = 280,
           py::arg("nc") = 96)
      .def_property_readonly("channel", &TaggedModel::channel)
      .def_property_readonly("k", [](const TaggedModel& t) {
        return copy_array(t.k());
      })
      .def_property_readonly("c", [](const TaggedModel& t) {
        return copy_array(t.c());
      })
      .def("radial_table", [](const TaggedModel& t, int l) {
        return copy_array(t.radial_table(l));
      }, py::arg("l"),
         "The NORMALIZED radial table psihat_L(k) * sqrt(P_L) of wave L, on "
         "the model's own k grid.")
      .def_property_readonly("m_struck_values", [](const TaggedModel& t) {
        return copy_array(t.m_struck_values());
      })
      .def("n_of_kc", [](const TaggedModel& t, double m_ion, double k,
                         double c) { return t.n_of_kc(m_ion, k, c); },
           py::arg("m_ion"), py::arg("k"), py::arg("c"))
      .def("n_of_kc_table", [](const TaggedModel& t, double m_ion) {
        return move_array2(std::vector<double>(t.n_of_kc(m_ion)),
                           static_cast<py::ssize_t>(t.nc()));
      }, py::arg("m_ion"))
      .def("struck_populations", [](const TaggedModel& t, double m_ion) {
        // (n_mS, nk, nc), the same shape polligen's TaggedModel returns.
        // Exposed 2026-09-06 so validation/repin_tagged_from_lipolgen.py can
        // dump the reference through the INSTALLED module rather than
        // reimplementing build_amp2 in Python.
        const std::vector<std::vector<double>> p = t.struck_populations(m_ion);
        std::vector<double> flat;
        flat.reserve(p.size() * t.nk() * t.nc());
        for (const std::vector<double>& row : p)
          flat.insert(flat.end(), row.begin(), row.end());
        return move_array3(std::move(flat),
                           static_cast<py::ssize_t>(t.nk()),
                           static_cast<py::ssize_t>(t.nc()));
      }, py::arg("m_ion"),
         "p(m_S | M, k, cos theta_k) on the model grid, shape (n_mS, nk, nc).")
      .def("population_integrated", [](const TaggedModel& t, double m_ion) {
        return move_array(t.population_integrated(m_ion));
      }, py::arg("m_ion"))
      .def("norm", &TaggedModel::norm, py::arg("m_ion"))
      .def("vector_dilution", [](const TaggedModel& t) {
        return t.vector_dilution();
      })
      .def("tensor_dilution", [](const TaggedModel& t) {
        return t.tensor_dilution();
      })
      .def("p2_moment", &TaggedModel::p2_moment, py::arg("m_ion"))
      .def("p2_moment_mixture", &TaggedModel::p2_moment_mixture,
           py::arg("populations"));

  py::class_<SpectatorLab>(m, "SpectatorLab")
      .def_readonly("pT", &SpectatorLab::pT)
      .def_readonly("theta", &SpectatorLab::theta)
      .def_readonly("p_lab", &SpectatorLab::p_lab)
      .def_readonly("R", &SpectatorLab::R)
      .def_readonly("xL", &SpectatorLab::xL)
      .def_readonly("kx", &SpectatorLab::kx)
      .def_readonly("ky", &SpectatorLab::ky)
      .def_readonly("kz", &SpectatorLab::kz)
      .def_readonly("phi_spec", &SpectatorLab::phi_spec)
      .def_readonly("e_lab", &SpectatorLab::e_lab)
      .def_readonly("pz_lab", &SpectatorLab::pz_lab);
  m.def("boost_spectator", &boost_spectator, py::arg("channel"), py::arg("k"),
        py::arg("c"), py::arg("phi_k"), py::arg("p_per_nucleon"),
        py::arg("theta_s") = 0.0, py::arg("phi_s") = 0.0);

  py::class_<StruckCluster>(m, "StruckCluster")
      .def(py::init<>())
      .def_readonly("p_ion", &StruckCluster::p_ion)
      .def_readonly("p_spectator", &StruckCluster::p_spectator)
      .def_readonly("p", &StruckCluster::p)
      .def_readonly("m2", &StruckCluster::m2)
      .def_readonly("m_free", &StruckCluster::m_free)
      .def_readonly("virtuality", &StruckCluster::virtuality)
      .def_readonly("alpha_s", &StruckCluster::alpha_s)
      .def_readonly("alpha_x", &StruckCluster::alpha_x)
      .def_readonly("pt_s", &StruckCluster::pt_s)
      .def_readonly("p_per_nucleon_eff", &StruckCluster::p_per_nucleon_eff);
  m.def("struck_cluster", &struck_cluster, py::arg("channel"), py::arg("lab"),
        py::arg("p_per_nucleon"));

  m.def("azz_tensor_curve", [](const TaggedModel& t, std::size_t ic) {
    return move_array(azz_tensor_curve(t, ic));
  }, py::arg("model"), py::arg("ic"));
  m.def("azz_tensor_curve_weighted",
        [](const TaggedModel& t,
           py::array_t<double, py::array::c_style | py::array::forcecast> w) {
    const std::size_t want = t.nk() * t.nc();
    if (static_cast<std::size_t>(w.size()) != want) {
      throw std::runtime_error("azz_tensor_curve_weighted: weights must be an "
                               "(nk, nc) table of the model's own grid");
    }
    const double* p = static_cast<const double*>(w.data());
    return move_array(azz_tensor_curve_weighted(
        t, std::vector<double>(p, p + want)));
  }, py::arg("model"), py::arg("weights"),
     "The ACCEPTANCE-WEIGHTED wave-function tensor asymmetry vs k; "
     "`weights` is an (nk, nc) table such as `acceptance_weights`.");
  m.def("acceptance_weights", [](const TaggedModel& t, double p_u,
                                 const Optics& o, const std::string& pc,
                                 std::size_t n_phi, double theta_s,
                                 double phi_s) {
    return move_array2(acceptance_weights(t, p_u, o, pc, n_phi, theta_s, phi_s),
                       static_cast<py::ssize_t>(t.nc()));
  }, py::arg("model"), py::arg("p_per_nucleon"), py::arg("optics"),
     py::arg("pot_config"), py::arg("n_phi") = 64, py::arg("theta_s") = 0.0,
     py::arg("phi_s") = 0.0);

  py::class_<IonFill>(m, "IonFill")
      .def(py::init<>())
      .def_readwrite("name", &IonFill::name)
      .def_readwrite("j", &IonFill::j)
      .def_readwrite("populations", &IonFill::populations)
      .def_readwrite("lam_e", &IonFill::lam_e)
      .def_readwrite("pe", &IonFill::pe)
      .def_readwrite("theta_s", &IonFill::theta_s)
      .def_readwrite("phi_s", &IonFill::phi_s);

  py::class_<TaggedEvent>(m, "TaggedEvent")
      .def_readonly("x", &TaggedEvent::x)
      .def_readonly("q2", &TaggedEvent::q2)
      .def_readonly("y", &TaggedEvent::y)
      .def_readonly("phi", &TaggedEvent::phi)
      .def_readonly("m_ion", &TaggedEvent::m_ion)
      .def_readonly("m_struck", &TaggedEvent::m_struck)
      .def_readonly("k", &TaggedEvent::k)
      .def_readonly("cos_theta_k", &TaggedEvent::cos_theta_k)
      .def_readonly("phi_k", &TaggedEvent::phi_k)
      .def_readonly("lab", &TaggedEvent::lab)
      .def_readonly("route", &TaggedEvent::route)
      .def_readonly("weight", &TaggedEvent::weight);

  m.def("nuclide_pdg", &nuclide_pdg, py::arg("z"), py::arg("a"));
}

// -------------------------------------------------------------------- fsi

// -------------------------------------------------- cluster_config.hpp (G)
//
// OPT-IN AND INERT: nothing in the generator calls this module.  Every option
// field, every analytic predictor and every constant is bound, and
// python/tests/test_cluster_config.py exercises them.

static void bind_cluster_config(py::module_& m) {
  m.attr("CONFIG_STREAM") = py::int_(kConfigStream);
  m.attr("LI6_R_POINT_VMC_FM") = LI6_R_POINT_VMC_FM;
  m.attr("LI6_QUADRUPOLE_GFMC_FM2") = LI6_QUADRUPOLE_GFMC_FM2;
  m.attr("LI6_QUADRUPOLE_GFMC_ERR_FM2") = LI6_QUADRUPOLE_GFMC_ERR_FM2;
  m.attr("LI6_ETA_DS_GK") = LI6_ETA_DS_GK;
  m.attr("LI6_ETA_DS_GK_STAT") = LI6_ETA_DS_GK_STAT;
  m.attr("LI6_ETA_DS_GK_SYST") = LI6_ETA_DS_GK_SYST;
  m.attr("LI6_R2_POINT_FM2") = LI6_R2_POINT_FM2;
  m.attr("LI6_QUADRUPOLE_FM2") = LI6_QUADRUPOLE_FM2;
  m.attr("LI7_QUADRUPOLE_FM2") = LI7_QUADRUPOLE_FM2;
  m.attr("VMC_HE4_DENSITY") = VMC_HE4_DENSITY;
  m.attr("VMC_LI6_DENSITY") = VMC_LI6_DENSITY;
  m.attr("VMC_LI6_AD_FIT") = VMC_LI6_AD_FIT;
  m.attr("VMC_DEUTERON_WAVE") = VMC_DEUTERON_WAVE;

  py::enum_<AlphaCoreSource>(m, "AlphaCoreSource",
      "Where the alpha core's one-body point-nucleon density comes from.")
      .value("VmcHe4Density", AlphaCoreSource::VmcHe4Density)
      .value("Gaussian", AlphaCoreSource::Gaussian);

  py::enum_<AlphaDSource>(m, "AlphaDSource",
      "Where the alpha-d relative radial functions R_0, R_2 come from.")
      .value("FitRescaled", AlphaDSource::FitRescaled)
      .value("FitRaw", AlphaDSource::FitRaw)
      .value("OverlapRaw", AlphaDSource::OverlapRaw);

  py::enum_<SndConfigFormat>(m, "SndConfigFormat",
      "he3.dat-compatible (numeric rows only) or '#'-annotated.")
      .value("He3Compatible", SndConfigFormat::He3Compatible)
      .value("Annotated", SndConfigFormat::Annotated);

  py::class_<ClusterConfigOptions>(m, "ClusterConfigOptions",
      "All 13 knobs of the 6Li configuration sampler; every one of them is "
      "recorded in the writer's sidecar.")
      .def(py::init<>())
      .def_readwrite("alpha_source", &ClusterConfigOptions::alpha_source)
      .def_readwrite("alpha_d_source", &ClusterConfigOptions::alpha_d_source)
      .def_readwrite("theta_s", &ClusterConfigOptions::theta_s)
      .def_readwrite("phi_s", &ClusterConfigOptions::phi_s)
      .def_readwrite("alpha_cm_inflate", &ClusterConfigOptions::alpha_cm_inflate)
      .def_readwrite("min_nn_separation_fm",
                     &ClusterConfigOptions::min_nn_separation_fm)
      .def_readwrite("alpha_d_scale", &ClusterConfigOptions::alpha_d_scale)
      .def_readwrite("quadrupole_target_fm2",
                     &ClusterConfigOptions::quadrupole_target_fm2)
      .def_readwrite("exact_coherence", &ClusterConfigOptions::exact_coherence)
      .def_readwrite("n_r", &ClusterConfigOptions::n_r)
      .def_readwrite("n_c", &ClusterConfigOptions::n_c)
      .def_readwrite("r_max_fm", &ClusterConfigOptions::r_max_fm)
      .def_readwrite("rnp_max_fm", &ClusterConfigOptions::rnp_max_fm)
      .def("validate", &ClusterConfigOptions::validate);

  py::class_<NucleonPos>(m, "NucleonPos",
      "One nucleon: position [fm] in the ion rest frame, isospin (+1 p, -1 n) "
      "and cluster (0 alpha core, 1 deuteron).")
      .def_readonly("x", &NucleonPos::x)
      .def_readonly("y", &NucleonPos::y)
      .def_readonly("z", &NucleonPos::z)
      .def_readonly("isospin", &NucleonPos::isospin)
      .def_readonly("cluster", &NucleonPos::cluster);

  py::class_<ClusterConfig>(m, "ClusterConfig",
      "One 6Li configuration: six nucleons, c.m. at the origin.")
      .def_property_readonly("nucleon", [](const ClusterConfig& c) {
        std::vector<NucleonPos> v(c.nucleon.begin(), c.nucleon.end());
        return v;
      })
      .def_readonly("m_ion", &ClusterConfig::m_ion)
      .def_readonly("m_s", &ClusterConfig::m_s)
      .def_readonly("r_ad", &ClusterConfig::r_ad)
      .def_readonly("r_np", &ClusterConfig::r_np)
      .def_property_readonly("positions", [](const ClusterConfig& c) {
        std::vector<double> p;
        p.reserve(18);
        for (int i = 0; i < 6; ++i) {
          p.push_back(c.nucleon[i].x);
          p.push_back(c.nucleon[i].y);
          p.push_back(c.nucleon[i].z);
        }
        return move_array2(std::move(p), 3);
      });

  py::class_<ClusterConfigSet>(m, "ClusterConfigSet",
      "N configurations plus the moments and per-configuration standard "
      "deviations every Monte Carlo gate is written against.")
      .def_readonly("config", &ClusterConfigSet::config)
      .def_readonly("m_ion", &ClusterConfigSet::m_ion)
      .def_readonly("theta_s", &ClusterConfigSet::theta_s)
      .def_readonly("phi_s", &ClusterConfigSet::phi_s)
      .def_readonly("seed", &ClusterConfigSet::seed)
      .def_readonly("run", &ClusterConfigSet::run)
      .def_readonly("provenance", &ClusterConfigSet::provenance)
      .def_readonly("r2_mean_fm2", &ClusterConfigSet::r2_mean_fm2)
      .def_readonly("q_matter_fm2", &ClusterConfigSet::q_matter_fm2)
      .def_readonly("delta_perp_fm2", &ClusterConfigSet::delta_perp_fm2)
      .def_readonly("r2_sd_fm2", &ClusterConfigSet::r2_sd_fm2)
      .def_readonly("q_matter_sd_fm2", &ClusterConfigSet::q_matter_sd_fm2)
      .def_readonly("delta_perp_sd_fm2", &ClusterConfigSet::delta_perp_sd_fm2)
      .def_property_readonly("cm_max", [](const ClusterConfigSet& s) {
        return s.cm[0];
      })
      .def("__len__", [](const ClusterConfigSet& s) { return s.config.size(); })
      .def_property_readonly("positions", [](const ClusterConfigSet& s) {
        std::vector<double> p;
        p.reserve(s.config.size() * 18);
        for (const ClusterConfig& c : s.config) {
          for (int i = 0; i < 6; ++i) {
            p.push_back(c.nucleon[i].x);
            p.push_back(c.nucleon[i].y);
            p.push_back(c.nucleon[i].z);
          }
        }
        return move_array3(std::move(p), 6, 3);
      }, "(N, 6, 3) float64 array of nucleon positions [fm].")
      .def_property_readonly("isospin", [](const ClusterConfigSet& s) {
        std::vector<std::int64_t> t;
        t.reserve(s.config.size() * 6);
        for (const ClusterConfig& c : s.config) {
          for (int i = 0; i < 6; ++i) {
            t.push_back(static_cast<std::int64_t>(c.nucleon[i].isospin));
          }
        }
        return move_array2(std::move(t), 6);
      }, "(N, 6) int64 array of isospins (+1 proton, -1 neutron).")
      .def_property_readonly("m_ion_per_config", [](const ClusterConfigSet& s) {
        std::vector<std::int64_t> v;
        v.reserve(s.config.size());
        for (const ClusterConfig& c : s.config) {
          v.push_back(static_cast<std::int64_t>(c.m_ion));
        }
        return move_array(std::move(v));
      })
      .def_property_readonly("m_s", [](const ClusterConfigSet& s) {
        std::vector<double> v;
        v.reserve(s.config.size());
        for (const ClusterConfig& c : s.config) v.push_back(c.m_s);
        return move_array(std::move(v));
      });

  py::class_<RadialMoments>(m, "RadialMoments",
      "The (G4)/(G6) functionals of a tabulated radial pair.")
      .def_readonly("n0", &RadialMoments::n0)
      .def_readonly("n2", &RadialMoments::n2)
      .def_readonly("r2", &RadialMoments::r2)
      .def_readonly("q_int", &RadialMoments::q_int)
      .def_readonly("q_dd", &RadialMoments::q_dd)
      .def("norm", &RadialMoments::norm)
      .def("p_d", &RadialMoments::p_d)
      .def("quadrupole", &RadialMoments::quadrupole);

  m.def("cluster_amp2", &cluster_amp2, py::arg("f0"), py::arg("f2"),
        py::arg("m"), py::arg("m_s"), py::arg("c"),
        "(G1'): one m_S branch of the master angular density.");
  m.def("cluster_density", &cluster_density, py::arg("f0"), py::arg("f2"),
        py::arg("m"), py::arg("c"),
        "(G1): the m_S-summed master angular density.");
  m.def("radial_moments", &radial_moments, py::arg("x"), py::arg("f0"),
        py::arg("f2"));
  m.def("a2_from_quadrupole", &a2_from_quadrupole, py::arg("q_matter_fm2"),
        py::arg("a"), py::arg("t_abs"), py::arg("m"),
        "a_2(m) at |t| from the m = +-1 point-matter quadrupole; the m = 0 "
        "value follows internally from delta_0 = -2 delta_{+-1}.");
  m.def("quadrupole_from_a2_slope", &quadrupole_from_a2_slope,
        py::arg("a2_over_t"), py::arg("a"),
        "The EXACT inverse of a2_from_quadrupole at m = +-1: the point-matter "
        "quadrupole [fm^2] that a coefficient a_2(+-1)/|t| [GeV^-2] implies.  "
        "Feed it CoherentScenario.a2_m_state(1.0, 1) with a = 6 to read back "
        "what a scenario's eps_b0 assumes about Q(6Li).");

  py::class_<ClusterConfigSampler>(m, "ClusterConfigSampler",
      "Nucleon-position configurations of a polarized 6Li in the alpha + d "
      "cluster picture.  Immutable after construction, so it is thread safe "
      "with no lock.  Its tensor output carries a factor-of-several "
      "wave-function systematic -- always quote quadrupole_band_fm2().")
      .def(py::init<ClusterConfigOptions>(),
           py::arg("options") = ClusterConfigOptions{})
      .def_property_readonly("options", &ClusterConfigSampler::options)
      .def_property_readonly("provenance", &ClusterConfigSampler::provenance)
      .def("sample", &ClusterConfigSampler::sample, py::arg("rng"),
           py::arg("m_ion"))
      .def("sample_set", &ClusterConfigSampler::sample_set, py::arg("n"),
           py::arg("m_ion"), py::arg("seed"), py::arg("run") = 0)
      .def("sample_set_unpolarized",
           &ClusterConfigSampler::sample_set_unpolarized, py::arg("n"),
           py::arg("seed"), py::arg("run") = 0)
      .def("p_d_alpha_d", &ClusterConfigSampler::p_d_alpha_d)
      .def("s_alpha_d", &ClusterConfigSampler::s_alpha_d)
      .def("tensor_dilution", &ClusterConfigSampler::tensor_dilution)
      .def("q_matter_analytic_fm2",
           &ClusterConfigSampler::q_matter_analytic_fm2, py::arg("m"))
      .def("q_int_fm2", &ClusterConfigSampler::q_int_fm2)
      .def("q_dd_fm2", &ClusterConfigSampler::q_dd_fm2)
      .def("r2_analytic_fm2", &ClusterConfigSampler::r2_analytic_fm2)
      .def("delta_perp_analytic_fm2",
           &ClusterConfigSampler::delta_perp_analytic_fm2, py::arg("m"))
      .def("asymptotic_ds_ratio", &ClusterConfigSampler::asymptotic_ds_ratio)
      .def("quadrupole_for_eta", &ClusterConfigSampler::quadrupole_for_eta,
           py::arg("eta_target"),
           "The quadrupole_target_fm2 that puts asymptotic_ds_ratio() on "
           "`eta_target`.  A CONVERTER onto the one existing dial, not a "
           "second dial: eta is exactly linear in quadrupole_dial_s().")
      .def("a2_from_geometry", &ClusterConfigSampler::a2_from_geometry,
           py::arg("t_abs"), py::arg("m"))
      .def("eps_b0_equivalent", &ClusterConfigSampler::eps_b0_equivalent)
      .def("match_li6_radius", &ClusterConfigSampler::match_li6_radius,
           py::arg("target_rms_fm") = LI6_R_POINT_VMC_FM)
      .def("quadrupole_band_fm2", &ClusterConfigSampler::quadrupole_band_fm2,
           "{measured, GFMC AV18+IL7, THIS source's} Q_charge [fm^2].")
      .def("quadrupole_dial_s", &ClusterConfigSampler::quadrupole_dial_s)
      .def("r2_grid_fm2", &ClusterConfigSampler::r2_grid_fm2)
      .def("q_matter_grid_fm2", &ClusterConfigSampler::q_matter_grid_fm2,
           py::arg("m"))
      .def("p2_alpha_d_grid", &ClusterConfigSampler::p2_alpha_d_grid,
           py::arg("m"), py::arg("m_s"))
      .def("p_ms", &ClusterConfigSampler::p_ms, py::arg("m"), py::arg("m_s"))
      .def("rho_alpha", &ClusterConfigSampler::rho_alpha, py::arg("r_fm"))
      .def("rho_alpha_recentred", &ClusterConfigSampler::rho_alpha_recentred,
           py::arg("s_fm"))
      .def("rho_alpha_d", &ClusterConfigSampler::rho_alpha_d, py::arg("R_fm"),
           py::arg("c"), py::arg("m"), py::arg("m_s"))
      .def("rho_alpha_d_summed", &ClusterConfigSampler::rho_alpha_d_summed,
           py::arg("R_fm"), py::arg("c"), py::arg("m"))
      .def("rho_np", &ClusterConfigSampler::rho_np, py::arg("r_fm"),
           py::arg("c"), py::arg("m_s"))
      .def("r2_alpha_fm2", &ClusterConfigSampler::r2_alpha_fm2)
      .def("r2_np_fm2", &ClusterConfigSampler::r2_np_fm2)
      .def("r2_alpha_d_fm2", &ClusterConfigSampler::r2_alpha_d_fm2)
      .def_property_readonly("alpha_d_grid", [](const ClusterConfigSampler& s) {
        return copy_array(s.alpha_d_grid());
      })
      .def("alpha_d_wave", [](const ClusterConfigSampler& s, int l) {
        return copy_array(s.alpha_d_wave(l));
      }, py::arg("l"))
      .def_property_readonly("np_grid", [](const ClusterConfigSampler& s) {
        return copy_array(s.np_grid());
      })
      .def("np_wave", [](const ClusterConfigSampler& s, int l) {
        return copy_array(s.np_wave(l));
      }, py::arg("l"));

  m.def("write_snd_configs", &write_snd_configs, py::arg("set"),
        py::arg("sampler"), py::arg("path"),
        py::arg("fmt") = SndConfigFormat::He3Compatible,
        "Write the configuration table plus its mandatory `.meta.json` "
        "sidecar; returns the number of rows.  `md5` and `git` are emitted "
        "as null for lipolgen.configs to fill in.");
}

static void bind_fsi(py::module_& m) {
  m.attr("GEV2_TO_MB") = GEV2_TO_MB;
  m.def("cluster_point_a2_fm2", &cluster_point_a2_fm2, py::arg("z"),
        py::arg("a"),
        "a^2 [fm^2] of the Gaussian point-nucleon density of a spectator "
        "cluster (alpha: 0.700).");

  py::enum_<FsiVariant>(m, "FsiVariant",
      "Which X-cluster profile the FSI weight is built from (fsi.hpp).")
      .value("GlauberCluster", FsiVariant::GlauberCluster)
      .value("GlauberNucleon", FsiVariant::GlauberNucleon);

  py::enum_<PipelineFsi>(m, "PipelineFsi",
      "PipelineConfig.fsi: FSI weight model of a run.  Off = today's "
      "plane-wave impulse approximation, bit for bit.")
      .value("Off", PipelineFsi::Off)
      .value("GlauberCluster", PipelineFsi::GlauberCluster)
      .value("GlauberNucleon", PipelineFsi::GlauberNucleon);
  m.def("pipeline_fsi_name", &pipeline_fsi_name, py::arg("fsi"));

  py::class_<GlauberFsiOptions>(m, "GlauberFsiOptions")
      .def(py::init<>())
      .def_readwrite("variant", &GlauberFsiOptions::variant)
      .def_readwrite("sigma_xn_mb", &GlauberFsiOptions::sigma_xn_mb,
                     "sigma_XN [mb].  40 = free hadron; documented BAND "
                     "20-40, band it, never quote one row alone.")
      .def_readwrite("eps", &GlauberFsiOptions::eps)
      .def_readwrite("b_xn", &GlauberFsiOptions::b_xn)
      .def_readwrite("elastic_gain", &GlauberFsiOptions::elastic_gain)
      .def_readwrite("formation_ramp", &GlauberFsiOptions::formation_ramp)
      .def_readwrite("ramp_w_lo", &GlauberFsiOptions::ramp_w_lo)
      .def_readwrite("ramp_sigma_lo_mb", &GlauberFsiOptions::ramp_sigma_lo_mb)
      .def_readwrite("ramp_w_hi", &GlauberFsiOptions::ramp_w_hi)
      .def_readwrite("ramp_sigma_hi_mb", &GlauberFsiOptions::ramp_sigma_hi_mb)
      .def_readwrite("n_sigma_grid", &GlauberFsiOptions::n_sigma_grid)
      .def_readwrite("k_max", &GlauberFsiOptions::k_max)
      .def_readwrite("n_kz", &GlauberFsiOptions::n_kz)
      .def_readwrite("n_kt", &GlauberFsiOptions::n_kt)
      .def_readwrite("w_max", &GlauberFsiOptions::w_max);

  py::class_<FsiKinematics>(m, "FsiKinematics",
      "Everything the FSI weight may depend on: (k, cos_theta_k, phi_k) in "
      "the ion rest frame, (w, q2, x) for the sigma_XN(W) ramp, the "
      "spectator (Z, A) and the channel label the guard checks.")
      .def(py::init([](double k, double cos_theta_k, double phi_k, double w,
                       double q2, double x, int spectator_z, int spectator_a,
                       Channel channel) {
        FsiKinematics kin;
        kin.k = k;
        kin.cos_theta_k = cos_theta_k;
        kin.phi_k = phi_k;
        kin.w = w;
        kin.q2 = q2;
        kin.x = x;
        kin.spectator_z = spectator_z;
        kin.spectator_a = spectator_a;
        kin.channel = channel;
        return kin;
      }), py::arg("k") = 0.0, py::arg("cos_theta_k") = 0.0,
          py::arg("phi_k") = 0.0, py::arg("w") = 0.0, py::arg("q2") = 0.0,
          py::arg("x") = 0.0, py::arg("spectator_z") = 2,
          py::arg("spectator_a") = 4,
          py::arg("channel") = Channel::Inclusive)
      .def_readwrite("k", &FsiKinematics::k)
      .def_readwrite("cos_theta_k", &FsiKinematics::cos_theta_k)
      .def_readwrite("phi_k", &FsiKinematics::phi_k)
      .def_readwrite("w", &FsiKinematics::w)
      .def_readwrite("q2", &FsiKinematics::q2)
      .def_readwrite("x", &FsiKinematics::x)
      .def_readwrite("spectator_z", &FsiKinematics::spectator_z)
      .def_readwrite("spectator_a", &FsiKinematics::spectator_a)
      .def_readwrite("channel", &FsiKinematics::channel);

  py::class_<FsiWeight, std::shared_ptr<FsiWeight>>(m, "FsiWeight",
      "FSI as a per-event weight -- NEVER a shift of any four-vector "
      "(fsi.hpp).  Abstract; see GlauberFsiWeight.")
      .def("weight", &FsiWeight::weight, py::arg("kin"))
      .def("weight_normalised", &FsiWeight::weight_normalised, py::arg("kin"))
      .def("sigma_eff_mb", &FsiWeight::sigma_eff_mb, py::arg("w"));

  py::class_<GlauberFsiWeight, FsiWeight, std::shared_ptr<GlauberFsiWeight>>(
      m, "GlauberFsiWeight",
      "Cosyn-Weiss / Glauber FSI on the tagged spectator cluster, spin "
      "independent by construction.  Immutable after construction and "
      "thread-safe.  Quote it as an unpolarized-shape SYSTEMATIC, never as "
      "a correction to A_zz.")
      .def(py::init<TaggedChannel, GlauberFsiOptions>(), py::arg("channel"),
           py::arg("options") = GlauberFsiOptions())
      .def(py::init<const TaggedModel&, GlauberFsiOptions>(),
           py::arg("model"), py::arg("options") = GlauberFsiOptions())
      .def_property_readonly("options", &GlauberFsiWeight::options)
      .def_property_readonly("channel", &GlauberFsiWeight::channel)
      .def("sigma_cluster_mb", &GlauberFsiWeight::sigma_cluster_mb,
           py::arg("sigma_xn_mb"))
      .def("sigma_cluster_el_mb", &GlauberFsiWeight::sigma_cluster_el_mb,
           py::arg("sigma_xn_mb"))
      .def("slope_cluster_gev2", &GlauberFsiWeight::slope_cluster_gev2,
           py::arg("sigma_xn_mb"))
      .def("gtilde", &GlauberFsiWeight::gtilde, py::arg("q_gev"),
           py::arg("sigma_xn_mb"))
      .def("gamma_profile", &GlauberFsiWeight::gamma_profile, py::arg("b_fm"),
           py::arg("sigma_xn_mb"))
      .def("survival", &GlauberFsiWeight::survival, py::arg("w") = 0.0,
           "The tagged-cluster survival probability int w dGamma / int "
           "dGamma; what weight_normalised divides by.  LOG IT.")
      .def("clipped_grid_fraction", &GlauberFsiWeight::clipped_grid_fraction)
      .def("sigma_ladder", [](const GlauberFsiWeight& f) {
        return copy_array(f.sigma_ladder());
      })
      .def("grid", [](const GlauberFsiWeight& f, std::size_t i) {
        return copy_array(f.grid(i));
      }, py::arg("i") = 0,
         "The raw (n_kz * n_kt) weight table of ladder point i, row-major.");

  // A direct TaggedSampler, so an analysis can draw weighted tagged events
  // without a Pipeline.  Minimal surface: construct on a model (kept alive),
  // attach an FSI model, draw a category.
  py::class_<TaggedSampler>(m, "TaggedSampler",
      "Spin-correlated (e', spectator) draws for one tagged channel.  The "
      "pipeline builds its own; this direct binding is for wave-function / "
      "FSI studies (no DIS source, so x/q2/y stay 0).")
      .def(py::init<const TaggedModel&, double>(), py::arg("model"),
           py::arg("p_per_nucleon"), py::keep_alive<1, 2>())
      .def("set_fsi", &TaggedSampler::set_fsi, py::arg("fsi"),
           "Attach an FsiWeight: every subsequent draw carries "
           "TaggedEvent.weight.  None detaches (weights back to 1).")
      .def("sample_category", &TaggedSampler::sample_category, py::arg("fill"),
           py::arg("n"), py::arg("rng"))
      .def("sigma_tot_pb", &TaggedSampler::sigma_tot_pb, py::arg("fill"));
}

// ---------------------------------------------------------------------- rc

/// pybind11 trampoline so an analysis can supply its OWN 6Li elastic form
/// factors (a digitisation, a new calculation) without touching the library:
/// `RcOptions.ff` takes any `Spin1ElasticFF`, and this makes a Python subclass
/// one.  It is `shared_ptr`-held because `RcModel` keeps a
/// `shared_ptr<const Spin1ElasticFF>`.
class PySpin1ElasticFF : public Spin1ElasticFF {
 public:
  using Spin1ElasticFF::Spin1ElasticFF;
  double fc(double t) const override {
    PYBIND11_OVERRIDE_PURE(double, Spin1ElasticFF, fc, t);
  }
  double fm(double t) const override {
    PYBIND11_OVERRIDE_PURE(double, Spin1ElasticFF, fm, t);
  }
  double fq(double t) const override {
    PYBIND11_OVERRIDE_PURE(double, Spin1ElasticFF, fq, t);
  }
  std::string provenance() const override {
    PYBIND11_OVERRIDE_PURE(std::string, Spin1ElasticFF, provenance, );
  }
};

static void bind_rc(py::module_& m) {
  m.attr("RC_DELTA_HIGH_X") = RC_DELTA_HIGH_X;
  m.attr("RC_X_HIGH") = RC_X_HIGH;
  m.attr("RC_DELTA_LOW_X") = RC_DELTA_LOW_X;
  m.attr("RC_DELTA_LOW_X_OPTIMISTIC") = RC_DELTA_LOW_X_OPTIMISTIC;
  m.attr("RC_X_LOW") = RC_X_LOW;
  m.attr("RC_BAND_TAU_MAX") = RC_BAND_TAU_MAX;
  m.attr("RC_TAIL_Y_CEILING") = RC_TAIL_Y_CEILING;
  m.attr("RC_QE_KF_GEV") = RC_QE_KF_GEV;
  m.attr("kRcWeightCount") = kRcWeightCount;

  py::enum_<RcMode>(m, "RcMode",
      "PipelineConfig.rc: the RC weight family.  Off (the default) is today "
      "bit for bit; TensorBand adds rc_tensor_lo / rc_tensor_hi (the band on "
      "the tensor part of the rate) and rc_tail (the 6Li radiative tails).")
      .value("Off", RcMode::Off)
      .value("TensorBand", RcMode::TensorBand);

  py::enum_<RcTailModel>(m, "RcTailModel",
      "Which tail formulation.  TPeak is the DEFAULT and is bit for bit "
      "every published number (POLRAD Eqs. (37)-(39), (43) -- the t-peak "
      "alone, a LOWER BOUND on the dilution).  TPeakPlusLL adds the "
      "leading-log s- and p-peaks to the UNPOLARISED numerator and is a "
      "STATED MODEL of mixed approximation orders, not a controlled O(alpha) "
      "expansion: it has NO tensor s/p partner, so it LOWERS the tensor "
      "fraction of the tail where it matters -- run it BESIDE TPeak as a "
      "systematic, never instead of it.  PolradFull (2026-09-06) is POLRAD "
      "Eq. (18) + Appendix B + Eq. (A.4): ONE exact tau_A quadrature "
      "containing all three peaks WITH the s-/p-peaks' own Eq. (A.4) tensor "
      "content, so it needs no s/p tensor stand-in -- and it still is NOT "
      "checked against Mo-Tsai or any external exact tail (no such number is "
      "in this tree), only POLRAD-internally and against the leading-log "
      "fallback.  It does NOT bracket inside the TPeak/TPeakPlusLL band: "
      "measured, it sits between them on 56.5 % of the sampler's cells.  It "
      "needs n_eta >= 64 and ~10 s of table build.")
      .value("TPeak", RcTailModel::TPeak)
      .value("PolradFull", RcTailModel::PolradFull)
      .value("TPeakPlusLL", RcTailModel::TPeakPlusLL);

  py::enum_<C0Shape>(m, "C0Shape",
      "Which C0 (monopole) shape HoSpin1FF runs -- it is shared by F_c AND "
      "F_q, so it moves both together and cannot be rescaled out of one "
      "run.\n"
      "  Ho      the DEFAULT: the (a, alpha) harmonic oscillator, bit for bit "
      "what every published number was made with.  Its F_c and F_q both "
      "change sign at q_0 = 3.0999 fm^-1.\n"
      "  VmcFt   the OTHER BAND EDGE: the j0 transform of "
      "data/vmc/density/li6.density, r-rescaled to the same measured "
      "<r^2>_point = 6.0788 fm^2, exponentially continued above q = 3 fm^-1 "
      "and therefore WITHOUT a C0 zero at any q.  Run BOTH: the two differ by "
      "x2.5 in F_point at q = 2 fm^-1 and the tensor fraction of the elastic "
      "tail at x = 0.1 is +1.56e-4 on Ho and -4.07e-5 on VmcFt, i.e. its SIGN "
      "is not determined (phase_B_numbers.md sec. B1).")
      .value("Ho", C0Shape::Ho)
      .value("VmcFt", C0Shape::VmcFt);

  m.def("c0_shape_name", &c0_shape_name, py::arg("shape"),
        "\"ho\" | \"vmc-ft\" -- the meta key, the CLI spelling and the "
        "provenance string all read this one function.");

  py::enum_<RcScope>(m, "RcScope",
      "Which rank-2 terms the band rescales.  TensorRate (the default) is the "
      "b1..b4 sector POLRAD and Gakh-Shekhovtsova actually compute; TensorAll "
      "also takes the Delta cos 2phi term, for which NO RC calculation exists "
      "at all -- for PRICING the omission, never for correcting it.")
      .value("TensorRate", RcScope::TensorRate)
      .value("TensorAll", RcScope::TensorAll);

  m.def("rc_tail_model_name", &rc_tail_model_name, py::arg("model"),
        "'t-peak', 'polrad-full' or 't-peak+ll' -- the ONE spelling, the "
        "same one the npz meta['rc_tail_model'] key carries.");
  m.def("rc_mode_name", &rc_mode_name, py::arg("mode"),
        "\"off\" | \"tensor-band\".");
  m.def("pipeline_rc_name", &pipeline_rc_name, py::arg("rc"));
  m.def("rc_weight_name",
        [](std::size_t i, std::size_t slot) { return rc_weight_name(i, slot); },
        py::arg("i"), py::arg("slot") = 0,
        "The HepMC3 / npz name of RC weight i in slot `slot`: "
        "rc_weight_name(2, 0) == 'rc_tail', rc_weight_name(2, 3) == "
        "'rc_tail_3'.");
  m.def("rc_delta", &rc_delta, py::arg("x"),
        py::arg("delta_high") = RC_DELTA_HIGH_X,
        py::arg("delta_low") = RC_DELTA_LOW_X,
        py::arg("x_high") = RC_X_HIGH, py::arg("x_low") = RC_X_LOW,
        "The band half-width delta(x): log-linear between the two anchors, "
        "clamped outside them, monotone non-increasing.  delta_low = 0.30 is "
        "the SIZE of a correction this generator does not apply (Gakh-"
        "Shekhovtsova hep-ph/0403262, ZERO INSPIRE citations), taken as a "
        "1-sigma band; 0.19 is HERMES's measured fractional residual at its "
        "lowest-x bin.  BAND IT: run both, never quote one row alone.");

  py::class_<RcOptions>(m, "RcOptions",
      "The RC knobs.  NOTE there is no `mode` here: PipelineConfig.rc is the "
      "single source of truth and RcModel takes it as a constructor "
      "argument.")
      .def(py::init<>())
      .def_readwrite("scope", &RcOptions::scope)
      .def_readwrite("delta_high_x", &RcOptions::delta_high_x)
      .def_readwrite("x_high", &RcOptions::x_high)
      .def_readwrite("delta_low_x", &RcOptions::delta_low_x,
                     "0.30 default (the conservative end of the uncited "
                     "10-30 %); 0.19 = 'as good as HERMES actually achieved'.")
      .def_readwrite("x_low", &RcOptions::x_low)
      .def_readwrite("a_transfer_frac", &RcOptions::a_transfer_frac,
                     "The A = 2 -> A = 6 TRANSFER uncertainty on delta(x), as "
                     "a fraction of delta itself, added in QUADRATURE: "
                     "delta_eff = delta * sqrt(1 + a_transfer_frac^2).  "
                     "design Q8 -- every anchor the band interpolates between "
                     "is a DEUTERON number, and whether the deuteron's "
                     "fractional RC transfers to 6Li is the largest "
                     "unquantified assumption in the band.  0.0 (default) is "
                     "'it transfers exactly', which is what v0 shipped "
                     "unstated; 0.5 and 1.0 price the doubt.  PRICES the "
                     "assumption, never corrects it.  Must be >= 0.")
      .def_readwrite("with_tail", &RcOptions::with_tail)
      .def_readwrite("with_qe_tail", &RcOptions::with_qe_tail,
                     "The UNPOLARISED quasi-elastic tail (POLRAD Eq. (44)).  "
                     "Its A_zz is ~0, so it DILUTES A_zz just as the "
                     "unpolarised elastic tail does; HERMES subtracted both.  "
                     "Turning it off prices the elastic tail alone and MUST "
                     "be labelled so.")
      .def_readwrite("c0_shape", &RcOptions::c0_shape,
                     "C0Shape.Ho (default) or C0Shape.VmcFt -- which monopole "
                     "the DEFAULT form factor runs.  IGNORED when `ff` is "
                     "user-supplied.  A SHAPE, not a multiplier: run both "
                     "edges, never rescale one.")
      .def_readwrite("tail_model", &RcOptions::tail_model,
                     "RcTailModel.TPeak (default, bit for bit every "
                     "published number), RcTailModel.TPeakPlusLL (the "
                     "leading-log s+p peaks added to the UNPOLARISED "
                     "numerator; a STATED MODEL with no tensor s/p partner, "
                     "so it LOWERS the tensor fraction of the tail) or "
                     "RcTailModel.PolradFull (Eq. (18)'s exact tau_A "
                     "quadrature, all three peaks with their tensor "
                     "content).  Run the t-peak pair as a PRICE RANGE -- "
                     "measured, PolradFull is outside it on 43 % of the "
                     "sampler's cells, so it is not a confidence interval.")
      // NOT `def_readwrite`.  `RcOptions::ff` is a shared_ptr<const
      // Spin1ElasticFF>, so assigning a PYTHON SUBCLASS through a plain
      // readwrite stores only the C++ trampoline and drops the Python object
      // as soon as the caller's last reference goes: the next `fc()` call
      // then throws `Tried to call pure virtual function`.  Store an ALIASING
      // shared_ptr whose control block owns a py::object, so the Python half
      // of the trampoline lives exactly as long as the C++ half.
      .def_property("ff",
                    [](const RcOptions& o) { return o.ff; },
                    [](RcOptions& o, py::object obj) {
                      if (obj.is_none()) { o.ff.reset(); return; }
                      std::shared_ptr<Spin1ElasticFF> sp =
                          obj.cast<std::shared_ptr<Spin1ElasticFF>>();
                      std::shared_ptr<py::object> keep =
                          std::make_shared<py::object>(obj);
                      o.ff = std::shared_ptr<const Spin1ElasticFF>(keep,
                                                                   sp.get());
                    },
                    "A Spin1ElasticFF to use instead of the built-in "
                    "HoSpin1FF; None takes the default.  fq_scale and "
                    "tail_tensor_scale are then IGNORED (they are applied "
                    "when RcModel builds its own form factor).  A PYTHON "
                    "SUBCLASS is kept alive by this assignment -- you do not "
                    "have to hold your own reference.")
      .def_readwrite("fq_scale", &RcOptions::fq_scale,
                     "The +-100 % 6Li quadrupole band.  sigma^el_T is "
                     "QUADRATIC in it, so RUN it (0, 1, 2) -- never rescale "
                     "one run.")
      .def_readwrite("tail_tensor_scale", &RcOptions::tail_tensor_scale,
                     "Flat multiplier on F_m -- the eta F_m^2 tensor sector, "
                     "which fq_scale does NOT span.  Run 0.5, 1, 2.")
      .def_readwrite("qe_suppression", &RcOptions::qe_suppression,
                     "Flat multiplier on the UNPOLARISED quasi-elastic "
                     "radiative tail, applied ON TOP of the de Forest-"
                     "Walecka Pauli suppression S(q) that is ON BY DEFAULT "
                     "at k_F = RC_QE_KF_GEV (`qe_kf_gev` below; 0 disables "
                     "it).  POLRAD Eq. (44)'s S_E/S_M are that S(q), not "
                     "this: 1.0 = no EXTRA suppression, the conservative "
                     "direction for a DILUTION -- the band knob, not the "
                     "physics.  Run 0.0 / 0.5 / 1.0.")
      .def_readwrite("qe_tensor_scale", &RcOptions::qe_tensor_scale,
                     "The POLARISED quasi-elastic tail, priced by a STAND-IN "
                     "and not computed.  DEFAULT 0.0 = the shipped "
                     "tensor-blind quasi-elastic tail, bit for bit.  At 1.0 "
                     "the quasi-elastic tail is given the ELASTIC tail's own "
                     "sigma^el_T/sigma^el_U, which is a BORROWED magnitude "
                     "(a coherent nuclear quadrupole fraction on an "
                     "incoherent nucleon process), NOT a derived bound -- and "
                     "6Li's elastic tensor fraction is anomalously small for "
                     "a reason the quasi-elastic piece does not share, so 1.0 "
                     "may UNDERSTATE the omission by ~1e2 at x <= 0.1.  Read "
                     "the magnitude, never the sign.  It is LINEAR, so one "
                     "run rescales; it rides inside qe_suppression; and "
                     "with_qe_tail = False refuses a non-zero value.")
      .def_readwrite("sp_tensor_scale", &RcOptions::sp_tensor_scale,
                     "The TENSOR FRACTION OF THE LEADING-LOG s-/p-PEAKS, "
                     "priced by a STAND-IN and not computed.  DEFAULT 0.0 = "
                     "the shipped tensor-blind s+p, bit for bit, and "
                     "unreachable under the shipped t-peak tail anyway.  At "
                     "1.0 the ELASTIC s-/p-peaks are given the elastic "
                     "t-peak's own sigma^el_T/sigma^el_U, i.e. `the s/p "
                     "tensor fraction equals the elastic t-peak's'.  IT IS A "
                     "BOUND WITH NO DERIVATION: POLRAD's Eq. (38) supplies "
                     "no tensor s/p peak and none is derived here.  SAY WHICH "
                     "POLRAD -- Eq. (18) + Eq. (A.4) DOES carry it and "
                     "RcTailModel.PolradFull computes it, which is why this "
                     "scale is refused THERE for the opposite reason (the "
                     "term RAN); the unqualified sentence was true of "
                     "Eq. (38) and false of the paper.  It borrows less than "
                     "qe_tensor_scale does -- the SAME coherent 6Li vertex at "
                     "a different photon topology, not a coherent ratio on an "
                     "incoherent process -- but the s/p vertex sits at "
                     "Q'^2 ~ Q^2, where 6Li's coherent form factor is DEAD, "
                     "while the t-peak's sits at t ~ t_min where it is alive, "
                     "so the borrowed ratio may be the WRONG REFERENCE "
                     "ENTIRELY at high Q^2 -- and there the term it scales is "
                     "itself ~0, so the knob prices NOTHING exactly where the "
                     "s/p peaks dominate.  The QUASI-ELASTIC s/p column is "
                     "covered by neither this nor qe_tensor_scale and stays "
                     "exactly tensor-blind -- under PolradFull too, because "
                     "the quasi-elastic tail is a sum over SPIN-1/2 nucleons "
                     "whose Im_5..8 vanish identically.  Read the magnitude, "
                     "never the sign.  LINEAR, so one run rescales; it is "
                     "OUTSIDE qe_suppression; and any tail_model but "
                     "TPeakPlusLL refuses a non-zero value -- on TPeak "
                     "because the s/p column did not run, on PolradFull "
                     "because it DID (Eq. (18) carries the s-/p-peaks' own "
                     "Eq. (A.4) tensor content, so the stand-in would "
                     "double-count).")
      .def_readwrite("qe_kf_gev", &RcOptions::qe_kf_gev,
                     "POLRAD Eq. (44)'s S_E/S_M as `ffquas` codes them: the "
                     "de Forest-Walecka Fermi-gas factor "
                     "S(q) = (3/4)(q/k_F) - (q/k_F)^3/16 below q = 2 k_F.  "
                     "DEFAULT ON at 6Li's measured k_F = 0.169 GeV (Moniz et "
                     "al., PRL 26 (1971) 445).  It cuts the QRT to 0.47 at "
                     "x = 0.01 and 0.87 at x = 0.1, and the QRT is the "
                     "DOMINANT piece of rc_tail.  0 = the unsuppressed edge.")
      .def_readwrite("n_eta", &RcOptions::n_eta,
                     "The tail's one-dimensional quadrature resolution, and "
                     "it means two things: Gauss-Legendre nodes in ln(eta_A) "
                     "under TPeak/TPeakPlusLL, tanh-sinh nodes PER PANEL of "
                     "Eq. (18)'s tau_A integral (four panels) under "
                     "PolradFull, where the minimum is 64 and less is "
                     "refused -- sigma^el_U is 0.51 % low at 32.")
      .def_readwrite("m_lepton", &RcOptions::m_lepton,
                     "The LEPTON MASS, read by tail_model = PolradFull and by "
                     "nothing else: it is the m^2 of POLRAD Eq. (B.13)'s "
                     "C_1,2(tau), of F_IR = m^2 F_2+ - Q_m^2 F_d and of "
                     "lambda_s, i.e. what REGULATES the s- and p-peaks.  "
                     "UNUSED by the two t-peak tails, and PipelineConfig "
                     "REFUSES any value other than constants.hpp's "
                     "M_ELECTRON on them, so that a knob that did not run is "
                     "never recorded as if it had.  TPeakPlusLL's "
                     "leading-log radiator does carry ln(Q^2/m_e^2), but it "
                     "reads M_ELECTRON directly: a different lepton needs its "
                     "own elastic kinematics too, so honouring m_lepton in "
                     "the log alone would be a half-change dressed as a whole "
                     "one.  (It was refused on ALL THREE models, as "
                     "'RESERVED for PolradFull', until 2026-09-06.)")
      .def_readwrite("tail_max", &RcOptions::tail_max)
      .def_readwrite("band_tau_max", &RcOptions::band_tau_max,
                     "Ceiling on |tau| the BAND sees (default 1.0).  On the "
                     "TAGGED channels tau_tag = 1 - nbar/n_M is UNBOUNDED at "
                     "the nodes of the M-dependent spectator density and "
                     "without this the band edges go NEGATIVE.  With 1.0 "
                     "every edge stays in [1 - delta, 1 + delta].");

  py::class_<RcWeights>(m, "RcWeights",
      "One event's RC weight triple: lo = 1 - delta(x) tau, "
      "hi = 1 + delta(x) tau, tail = 1 + sigma_tail/sigma_Born.")
      .def(py::init<>())
      .def_readwrite("lo", &RcWeights::lo)
      .def_readwrite("hi", &RcWeights::hi)
      .def_readwrite("tail", &RcWeights::tail)
      .def_readonly("band_clipped", &RcWeights::band_clipped,
                    "|tau| hit RcOptions.band_tau_max on this event.")
      .def_readonly("tail_clipped", &RcWeights::tail_clipped,
                    "The tail ratio hit RcOptions.tail_max on this event.  "
                    "NOT the same quantity as RcModel.clipped_cell_fraction, "
                    "which counts TABLE NODES and is not event-weighted.")
      .def("__repr__", [](const RcWeights& w) {
        char b[192];
        std::snprintf(b, sizeof b,
                      "RcWeights(lo=%.8g, hi=%.8g, tail=%.8g, band_clipped=%d,"
                      " tail_clipped=%d)",
                      w.lo, w.hi, w.tail, static_cast<int>(w.band_clipped),
                      static_cast<int>(w.tail_clipped));
        return std::string(b);
      });

  py::class_<NucleonFF>(m, "NucleonFF")
      .def_readonly("ge_p", &NucleonFF::ge_p)
      .def_readonly("gm_p", &NucleonFF::gm_p)
      .def_readonly("ge_n", &NucleonFF::ge_n)
      .def_readonly("gm_n", &NucleonFF::gm_n);
  m.def("nucleon_ff", &nucleon_ff, py::arg("t_gev2"),
        "Dipole G_E^p, G_M^p, G_M^n plus Galster's G_E^n (5 %).");

  py::class_<Spin1ElasticFF, PySpin1ElasticFF,
             std::shared_ptr<Spin1ElasticFF>>(m, "Spin1ElasticFF",
      "POLRAD Eq. (A.4)'s (F_c, F_m, F_q) at the elastic-vertex t [GeV^2], in "
      "the Rosenbluth normalisation F_c(0) = Z, F_m(0) = (M_A/m_p) mu_A/mu_N, "
      "F_q(0) = M_A^2 Q_A.  Subclass it in Python to supply your own.")
      .def(py::init<>())
      .def("fc", &Spin1ElasticFF::fc, py::arg("t_gev2"))
      .def("fm", &Spin1ElasticFF::fm, py::arg("t_gev2"))
      .def("fq", &Spin1ElasticFF::fq, py::arg("t_gev2"))
      .def("provenance", &Spin1ElasticFF::provenance);

  py::class_<HoSpin1FFOptions>(m, "HoSpin1FFOptions",
      "Every 0.0 means 'take it from the Ion' -- HoSpin1FF.for_ion fills "
      "them.  There are NO ion-specific defaults, so a 7Li run cannot "
      "silently get 6Li form factors.")
      .def(py::init<>())
      .def_readwrite("a_fm", &HoSpin1FFOptions::a_fm)
      .def_readwrite("alpha", &HoSpin1FFOptions::alpha)
      .def_readwrite("z", &HoSpin1FFOptions::z)
      .def_readwrite("m_a_gev", &HoSpin1FFOptions::m_a_gev)
      .def_readwrite("mu_n", &HoSpin1FFOptions::mu_n)
      .def_readwrite("q_fm2", &HoSpin1FFOptions::q_fm2)
      .def_readwrite("fm_qz_fm", &HoSpin1FFOptions::fm_qz_fm)
      .def_readwrite("fm_b_fm", &HoSpin1FFOptions::fm_b_fm)
      .def_readwrite("fq_scale", &HoSpin1FFOptions::fq_scale)
      .def_readwrite("tail_tensor_scale",
                     &HoSpin1FFOptions::tail_tensor_scale)
      .def_readwrite("c0_shape", &HoSpin1FFOptions::c0_shape,
                     "Which monopole F_c and F_q share -- C0Shape.Ho or "
                     "C0Shape.VmcFt.")
      .def_readwrite("fold_nucleon", &HoSpin1FFOptions::fold_nucleon);

  py::class_<HoSpin1FF, Spin1ElasticFF, std::shared_ptr<HoSpin1FF>>(
      m, "HoSpin1FF",
      "The v0 6Li model: a harmonic-oscillator point-nucleon shape for C0/C2, "
      "its own two-parameter shape for the magnetic one, normalised on the "
      "MEASURED moments (NOT on VMC -- WS98's Q(6Li) = -0.23(9) fm^2 is 3x "
      "the measured -0.0818).  The shape parameters are UNFITTED STARTING "
      "VALUES; see the provenance string.")
      .def_static("for_ion", &HoSpin1FF::for_ion, py::arg("ion"),
                  py::arg("options") = HoSpin1FFOptions(),
                  "THROWS for any ion with no measured-moment block (today: "
                  "everything but 6Li).")
      .def_property_readonly("options", &HoSpin1FF::options);

  py::class_<TabulatedSpin1FF, Spin1ElasticFF,
             std::shared_ptr<TabulatedSpin1FF>>(m, "TabulatedSpin1FF",
      "A digitised (q, F_C0, F_C2, F_M1) table under data/ff/, log-linear in "
      "q and refusing to extrapolate.  No 6Li table ships with the library -- "
      "this is the hook for one.")
      .def_static("from_data_dir", &TabulatedSpin1FF::from_data_dir,
                  py::arg("relative"), py::arg("norm") = HoSpin1FFOptions());

  py::class_<RcModel, std::shared_ptr<RcModel>>(m, "RcModel",
      "The tensor-sector RC weight family.  IMMUTABLE after construction and "
      "safe to share between threads.  A Pipeline builds its own; read it "
      "back with Pipeline.rc_model.")
      .def_property_readonly("options", &RcModel::options)
      .def_property_readonly("mode", &RcModel::mode)
      .def_property_readonly("ff_provenance", &RcModel::ff_provenance,
                             "What form factor the run actually used -- PRINT "
                             "IT.")
      .def_property_readonly("applies", &RcModel::applies)
      .def_property_readonly("tail_applies", &RcModel::tail_applies,
                             "TRUE at every theta_S: the tail is emitted at "
                             "any axis through P_zz^eff = 3 Q_NN "
                             "P_2(cos theta_S) (POLRAD Eq. (43)).")
      .def_property_readonly("exclusion_reason", &RcModel::exclusion_reason)
      .def("delta", &RcModel::delta, py::arg("x"))
      .def("tensor_fraction",
           (double (RcModel::*)(const Event&) const) &RcModel::tensor_fraction,
           py::arg("event"))
      .def("tensor_fraction",
           (double (RcModel::*)(const Event&, std::size_t) const)
               &RcModel::tensor_fraction,
           py::arg("event"), py::arg("category"))
      .def("tail_ratio", &RcModel::tail_ratio, py::arg("cell"), py::arg("q_n"))
      .def("tail_ratio_at",
           [](const RcModel& r, double x, double q2, double q_n) {
             return r.tail_ratio_at(x, q2, q_n);
           },
           py::arg("x"), py::arg("q2"), py::arg("q_n"),
           "sigma_tail/sigma_Born at an arbitrary (x, Q2).  q_n is POLRAD's "
           "Q_N = P_zz^eff = 3 Q_NN P_2(cos theta_S); Eq. (37) carries it as "
           "Q_N/6.")
      .def("weights",
           (RcWeights (RcModel::*)(const Event&) const) &RcModel::weights,
           py::arg("event"), "Slot 0: the event's own PURE spin state.")
      .def("weights",
           (RcWeights (RcModel::*)(const Event&, std::size_t) const)
               &RcModel::weights,
           py::arg("event"), py::arg("category"),
           "Slot 1 + k: spin category k's population MIXTURE, mirroring "
           "InclusiveSampler.weights_for.  It equals slot 0 only when the "
           "category's population vector is pure.")
      .def_property_readonly("clipped_cell_fraction",
                             &RcModel::clipped_cell_fraction)
      .def_property_readonly("clipped_fraction_by_y",
                             [](const RcModel& r) {
                               const std::array<double, 3> v =
                                   r.clipped_fraction_by_y();
                               return std::vector<double>(v.begin(), v.end());
                             },
                             "Clipped node fraction in the three y bands "
                             "(y < 0.5, 0.5-0.9, > 0.9).  LOG ALL FOUR: the "
                             "tail grows like Y_+ ~ 1/(1-y), so a clipped "
                             "edge hides inside a small global number.")
      .def_property_readonly("table_x", [](const RcModel& r) {
        return copy_array(r.table_x());
      })
      .def_property_readonly("table_y", [](const RcModel& r) {
        return copy_array(r.table_y());
      })
      .def_property_readonly("sigma_tail_u", [](const RcModel& r) {
        return copy_array(r.sigma_tail_u());
      })
      .def_property_readonly("sigma_tail_t", [](const RcModel& r) {
        return copy_array(r.sigma_tail_t());
      })
      .def_property_readonly("sigma_tail_qe", [](const RcModel& r) {
        return copy_array(r.sigma_tail_qe());
      })
      .def_property_readonly("sigma_tail_u_sp", [](const RcModel& r) {
        return copy_array(r.sigma_tail_u_sp());
      })
      .def_property_readonly("sigma_tail_qe_sp", [](const RcModel& r) {
        return copy_array(r.sigma_tail_qe_sp());
      })
      .def("tail_sigma_at", [](const RcModel& r, double x, double q2) {
        const RcModel::TailTriple t = r.tail_sigma_at(x, q2);
        return py::make_tuple(t.u, t.t, t.qe, t.u_sp, t.qe_sp);
      }, py::arg("x"), py::arg("q2"),
         "(sigma^el_U, sigma^el_T, sigma^q_U, sigma^el_U(s+p), "
         "sigma^q_U(s+p)) per nucleon [GeV^-2], straight from the quadrature "
         "-- no table, no interpolation.  The last two are the leading-log "
         "s+p peaks and are IDENTICALLY ZERO unless tail_model = "
         "RcTailModel.TPeakPlusLL; there is deliberately no tensor s/p "
         "partner on that model (see RcTailModel.TPeakPlusLL).  Under "
         "RcTailModel.PolradFull the split does not exist at all -- Eq. (18) "
         "puts all three peaks in ONE tau_A integral, so u/t/qe already "
         "contain the s- and p-peaks and the last two columns stay zero for "
         "a DIFFERENT reason: there is no decomposition to read.");
}

// ---------------------------------------------------------------- coherent

static void bind_coherent(py::module_& m) {
  m.def("gaussian_slope", &gaussian_slope, py::arg("r_rms_fm"));
  m.attr("COHERENT_T_MAX_DEFAULT") = COHERENT_T_MAX_DEFAULT;
  m.attr("COHERENT_MX_MIN_DEFAULT") = COHERENT_MX_MIN_DEFAULT;
  py::class_<EstarlightLi6Row>(m, "EstarlightLi6Row")
      .def_readonly("vm", &EstarlightLi6Row::vm)
      .def_readonly("sigma_nb", &EstarlightLi6Row::sigma_nb)
      .def_readonly("sigma_q7_nb", &EstarlightLi6Row::sigma_q7_nb)
      .def_readonly("sigma_rmeas_nb", &EstarlightLi6Row::sigma_rmeas_nb)
      .def_readonly("b_fit", &EstarlightLi6Row::b_fit)
      .def_readonly("branching", &EstarlightLi6Row::branching)
      .def_readonly("b_rmeas", &EstarlightLi6Row::b_rmeas)
      .def_readonly("branching_all", &EstarlightLi6Row::branching_all)
      .def("__repr__", [](const EstarlightLi6Row& r) {
        return std::string("EstarlightLi6Row(") + r.vm + ", sigma_nb=" +
               std::to_string(r.sigma_nb) + ")";
      });
  m.def("estarlight_li6_coherent", &estarlight_li6_coherent,
        py::return_value_policy::reference,
        "The eSTARlight unpolarized coherent VM baseline for e 10 GeV x 6Li "
        "99.5 GeV/u (open item 11.1); see coherent.hpp for the caveats.");
  m.attr("COHERENT_JPSI_EFF_IR8_LI7") = COHERENT_JPSI_EFF_IR8_LI7;
  m.attr("COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV") =
      COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV;
  m.attr("COHERENT_JPSI_EFF_IR8_LI7_W_MEAN_GEV") =
      COHERENT_JPSI_EFF_IR8_LI7_W_MEAN_GEV;
  py::class_<Chang26EffEnergyRow>(m, "Chang26EffEnergyRow")
      .def_readonly("e_electron_gev", &Chang26EffEnergyRow::e_electron_gev)
      .def_readonly("e_ion_gev", &Chang26EffEnergyRow::e_ion_gev)
      .def_readonly("efficiency", &Chang26EffEnergyRow::efficiency)
      .def("__repr__", [](const Chang26EffEnergyRow& r) {
        return std::string("Chang26EffEnergyRow(") +
               std::to_string(r.e_electron_gev) + " x " +
               std::to_string(r.e_ion_gev) + ", eff=" +
               std::to_string(r.efficiency) + ")";
      });
  m.def("chang26_he3_energy_scan", &chang26_he3_energy_scan,
        py::return_value_policy::reference,
        "arXiv:2511.05638 sec. V.B: the SAME far-forward efficiency as "
        "COHERENT_JPSI_EFF_IR8_LI7, measured on e+3He at three EIC energy "
        "configurations.  The only lever on the beam-energy dependence "
        "this tree has, and open item O5 needs one: the 0.1775 is a TOP-"
        "energy number applied at 10 x 99.5 GeV/u.");
  py::class_<Chang26SpeciesEffRow>(m, "Chang26SpeciesEffRow")
      .def_readonly("nucleus", &Chang26SpeciesEffRow::nucleus)
      .def_readonly("a", &Chang26SpeciesEffRow::a)
      .def_readonly("z", &Chang26SpeciesEffRow::z)
      .def_readonly("e_ion_gev", &Chang26SpeciesEffRow::e_ion_gev)
      .def_readonly("efficiency", &Chang26SpeciesEffRow::efficiency)
      .def("__repr__", [](const Chang26SpeciesEffRow& r) {
        return std::string("Chang26SpeciesEffRow(") + r.nucleus +
               ", eff=" + std::to_string(r.efficiency) + ")";
      });
  m.def("chang26_species_efficiency", &chang26_species_efficiency,
        py::return_value_policy::reference,
        "arXiv:2511.05638 p. 4's species list, with Fig. 2's beam "
        "energies.  Every row is at the SAME rigidity A/Z x E = 275 "
        "GeV/e -- and that is ALL it fixes: at fixed rigidity E/u = R Z/A "
        "and p_z = Z R both still vary (E/u 118-183 GeV/u, Z 1-8), so the "
        "A-ordering is NOT a species lever on its own.  Read off the four "
        "rows that share a beam energy (2D, 4He, 12C, 16O at 137, "
        "bracketing 6Li's own 137.5), the 7Li -> 6Li substitution open "
        "item O5 makes is x0.99-1.33 -- it straddles 1.");
  py::class_<EstarlightLi6Q2Row>(m, "EstarlightLi6Q2Row")
      .def_readonly("vm", &EstarlightLi6Q2Row::vm)
      .def_readonly("q2_floor_gev2", &EstarlightLi6Q2Row::q2_floor_gev2)
      .def_readonly("sigma_nb", &EstarlightLi6Q2Row::sigma_nb)
      .def_readonly("b_fit", &EstarlightLi6Q2Row::b_fit)
      .def_readonly("sigma_rmeas_nb", &EstarlightLi6Q2Row::sigma_rmeas_nb)
      .def_readonly("b_rmeas", &EstarlightLi6Q2Row::b_rmeas)
      .def_readonly("w_mean_gev", &EstarlightLi6Q2Row::w_mean_gev)
      .def("__repr__", [](const EstarlightLi6Q2Row& r) {
        return std::string("EstarlightLi6Q2Row(") + r.vm + ", Q2 > " +
               std::to_string(r.q2_floor_gev2) + ", sigma_nb=" +
               std::to_string(r.sigma_nb) + ")";
      });
  m.def("estarlight_li6_q2_floors", &estarlight_li6_q2_floors,
        py::return_value_policy::reference,
        "The 2026-09-04 photoproduction scan: the same eSTARlight runs with "
        "MIN_GAMMA_Q2 lowered to 0.01 and to no floor at all.  Q^2 < 0.1 "
        "GeV^2 carries most of the coherent rate and was absent from the "
        "0.1 < Q^2 < 100 window the rest of the chain was priced with.");
  py::class_<CoherentXpomModel>(m, "CoherentXpomModel")
      .def(py::init<>())
      .def_readwrite("m_x_min", &CoherentXpomModel::m_x_min)
      .def_readwrite("x_pom_max", &CoherentXpomModel::x_pom_max)
      .def_static("x_pom_of", &CoherentXpomModel::x_pom_of, py::arg("m_x2"),
                  py::arg("q2"), py::arg("w2"))
      .def_static("m_x2_of", &CoherentXpomModel::m_x2_of, py::arg("x_pom"),
                  py::arg("q2"), py::arg("w2"))
      .def("x_pom_min", &CoherentXpomModel::x_pom_min, py::arg("q2"),
           py::arg("w2"))
      .def("draw", &CoherentXpomModel::draw, py::arg("q2"), py::arg("w2"),
           py::arg("rng"));
  py::class_<CoherentScenario>(m, "CoherentScenario")
      .def(py::init<>())
      .def_readwrite("f0", &CoherentScenario::f0)
      .def_readwrite("x_coh", &CoherentScenario::x_coh)
      .def_readwrite("slope_b", &CoherentScenario::slope_b)
      .def_readwrite("amp", &CoherentScenario::amp)
      .def_readwrite("eps_b0", &CoherentScenario::eps_b0)
      .def("coherent_fraction", &CoherentScenario::coherent_fraction,
           py::arg("x"))
      .def("dsigma_dt", &CoherentScenario::dsigma_dt, py::arg("t_abs"))
      .def("tag_acceptance", &CoherentScenario::tag_acceptance,
           py::arg("pt_cut"))
      .def("mean_t_tagged", &CoherentScenario::mean_t_tagged, py::arg("pt_cut"))
      .def("tag_acceptance_angular", &CoherentScenario::tag_acceptance_angular,
           py::arg("sigma_theta"), py::arg("p_per_nucleon"),
           py::arg("a_beam") = 6, py::arg("n_sigma") = 10.0)
      .def("delta_b_m", &CoherentScenario::delta_b_m, py::arg("m"),
           "Delta B_m [GeV^-2], the cos 2(Phi - Phi_S) coefficient of the m "
           "state's |F|^2 slope.  THE definition of Delta B in this library: "
           "eps_b0 = 2 Delta B_{+-1}/B = -Delta B_0/B.")
      .def("slope_at_azimuth", &CoherentScenario::slope_at_azimuth,
           py::arg("phi_rel"), py::arg("m"),
           "B + Delta B_m cos 2 phi_rel [GeV^-2]; every a_2 is its "
           "first-order expansion.")
      .def("a2_deformation", &CoherentScenario::a2_deformation,
           py::arg("t_abs"), py::arg("pzz"))
      .def("cos2phi_coefficient_deformation",
           &CoherentScenario::cos2phi_coefficient_deformation, py::arg("t_abs"),
           py::arg("pzz"))
      .def("a2_tagged", &CoherentScenario::a2_tagged, py::arg("pt_cut"),
           py::arg("pzz"))
      .def("a2_m_state", &CoherentScenario::a2_m_state, py::arg("t_abs"),
           py::arg("m"))
      .def("cos2phi_coefficient", &CoherentScenario::cos2phi_coefficient,
           py::arg("t_abs"), py::arg("pzz"))
      .def("positivity_margin", &CoherentScenario::positivity_margin,
           py::arg("t_max"), py::arg("pzz"),
           "1 - |c_2(t_max, pzz)|.  NEGATIVE means the azimuthal weight is "
           "not a density any more; CoherentSampler throws on it.")
      .def("t_positivity_edge", &CoherentScenario::t_positivity_edge,
           py::arg("pzz"),
           "|t| [GeV^2] where positivity_margin reaches zero, in closed "
           "form: (1 - sign(A) C)/|A| with c_2 = A|t| + C.  DERIVED, so it "
           "moves with eps_b0 -- 0.245 at the shipped -0.08 and P_zz = -2, "
           "2.80 at the measured 6Li quadrupole's -0.0070.  It is NOT the "
           "derivation of COHERENT_T_MAX_DEFAULT (the anchor range is); see "
           "coherent.hpp.");

  py::class_<CoherentRecoil>(m, "CoherentRecoil")
      .def_readonly("pT", &CoherentRecoil::pT)
      .def_readonly("theta", &CoherentRecoil::theta)
      .def_readonly("R", &CoherentRecoil::R)
      .def_readonly("xL", &CoherentRecoil::xL)
      .def_readonly("phi_t", &CoherentRecoil::phi_t);
  m.def("recoil_lab", &recoil_lab, py::arg("t_abs"), py::arg("phi_t"),
        py::arg("p_per_nucleon"), py::arg("x_pom") = 0.0,
        py::arg("a_beam") = 6);
  m.def("fragment_rigidity", &fragment_rigidity, py::arg("a"), py::arg("z"),
        py::arg("beam_a") = 6, py::arg("beam_z") = 3);
}

// ---------------------------------------------------------------- pipeline

static void bind_pipeline(py::module_& m) {
  py::enum_<PipelineChannel>(m, "PipelineChannel")
      .value("Inclusive", PipelineChannel::Inclusive)
      .value("TaggedLi6Alpha", PipelineChannel::TaggedLi6Alpha)
      .value("TaggedLi7Alpha", PipelineChannel::TaggedLi7Alpha)
      .value("TaggedDeuteronP", PipelineChannel::TaggedDeuteronP)
      .value("CoherentLi6", PipelineChannel::CoherentLi6);
  m.def("pipeline_channel_name", &pipeline_channel_name, py::arg("c"));
  m.def("is_tagged", &is_tagged, py::arg("c"));
  m.def("channel_isotope", &channel_isotope, py::arg("c"));

  // ---- tier T1: the cluster breakup (breakup.hpp) ------------------------
  py::enum_<Tier>(m, "Tier",
      "Fidelity tier of the final state: T0 keeps the struck cluster as one "
      "off-shell pseudo-particle, T1 resolves it into a struck nucleon plus "
      "on-shell partner spectators.")
      .value("T0", Tier::T0)
      .value("T1", Tier::T1);

  py::enum_<ClusterSpecies>(m, "ClusterSpecies")
      .value("Nucleon", ClusterSpecies::Nucleon)
      .value("Deuteron", ClusterSpecies::Deuteron)
      .value("Triton", ClusterSpecies::Triton);
  m.def("cluster_species", &cluster_species, py::arg("z"), py::arg("a"));
  m.attr("KAPPA_NN_VIRTUAL") = KAPPA_NN_VIRTUAL;
  m.attr("KAPPA_PN_SINGLET") = KAPPA_PN_SINGLET;

  // ---- the triton spectral function (triton_sf.hpp) ----------------------
  py::enum_<TritonChannel>(m, "TritonChannel",
      "Which of the triton's three breakup channels a draw came out in.")
      .value("NeutronD", TritonChannel::NeutronD)
      .value("NeutronPnCont", TritonChannel::NeutronPnCont)
      .value("ProtonNnCont", TritonChannel::ProtonNnCont);
  m.def("triton_channel_name", &triton_channel_name, py::arg("c"));

  py::enum_<TritonSfChoice>(m, "TritonSfChoice",
      "Spectral function of the 7Li alpha tag's T1 triton breakup: Hulthen "
      "(default, the sequential two-body decay, bit-compatible) or "
      "CiofiSimula (the three-channel Ciofi degli Atti-Simula model built by "
      "the Pipeline at the run's own cluster_beta).")
      .value("Hulthen", TritonSfChoice::Hulthen)
      .value("CiofiSimula", TritonSfChoice::CiofiSimula);

  py::enum_<B1Model>(m, "B1Model",
      "Which b1 backend `default_inclusive_kernel` fills the INCLUSIVE "
      "kernel's rank-2 slot with (CLI --b1-model).\n"
      "  Miller          the DEFAULT and bit-for-bit what every published "
      "inclusive tensor number was made with: Li6B1(MillerB1) through "
      "LI6_B1_RANK2_TRANSFER.  This IS the run PLAN's 'toy'.\n"
      "  Cdks            the same 6Li rank-2 transfer on the CDKS convolution "
      "camp; |b1| two orders of magnitude smaller below x ~ 0.1, comparable "
      "or larger above it (peak |x b1| 3.34e-4 at x = 0.766 against Miller's "
      "4.27e-4 at x = 0.084, both at Q2 = 2.5, and 7.7x larger at x = 0.3 "
      "with the opposite sign), and a different sign "
      "structure.  Miller and CDKS are different CAMPS -- say which one a "
      "plot used.  INCLUSIVE CHANNEL ONLY and 6Li ONLY, like Li6Convolution: "
      "it is Li6B1's 6Li transfer (PipelineConfig.validate refuses the "
      "rest).\n"
      "  Li6Convolution  b1_nuclear.hpp's four-term alpha-d convolution.  "
      "INCLUSIVE CHANNEL ONLY and 6Li ONLY (PipelineConfig.validate refuses "
      "the rest; on a tagged channel the alpha-d density is already in the "
      "event weight).  Its A = 2 magnitude gate PASSES since 2026-09-03 -- "
      "for the MSTW2008 LO unpolarised nucleon input at CDKS Eq. (21)'s "
      "delta-function (B1UnpolSource.Mstw, CLI --b1-unpol mstw), where G3b "
      "is 0.843243.  On the SHIPPED DEFAULT B1UnpolSource.Toy it is 0.440, "
      "OUTSIDE the gate's [0.5, 2] acceptance window, so quote numbers made "
      "with b1_unpol = Mstw; the toy backend's are not covered by the lift.  "
      "And it is an A = 2 gate either way -- band every number with "
      "--b1-band-scale 0/1/2.")
      .value("Miller", B1Model::Miller)
      .value("Cdks", B1Model::Cdks)
      .value("Li6Convolution", B1Model::Li6Convolution);
  m.def("b1_model_name", &b1_model_name, py::arg("model"));

  py::enum_<B1UnpolSource>(m, "B1UnpolSource",
      "Which UNPOLARISED backend b1_model = Li6Convolution folds its own F1 "
      "against (Li6ConvolutionOptions.unpol; CLI --b1-unpol).  It reaches the "
      "b1 and NOTHING else -- InclusiveKernel's own f2_source is set by a "
      "SEPARATE flag, --unpol-sf (UnpolSfSource), and at its shipped default "
      "'toy' it is ToyF2 on every --b1-unpol setting, so the SPIN-BLIND cell "
      "cross section (InclusiveSampler.cell_xsec_pb) is bit for bit under "
      "THIS flag and only the tensor shift moves.  The two are kept separate "
      "because they are different quantities (this one is the deuteron F1 "
      "inside CDKS Eq. (22), a CDKS-comparability choice the A = 2 gate "
      "verdict rests on); validate() refuses the one combination that would "
      "mislabel them, b1_unpol = Toy under a non-toy unpol_sf.\n"
      "  Toy      the DEFAULT: the kernel's own ToyF2, shared as ONE object.  "
      "Bit for bit what every published number was made with.\n"
      "  Mstw     MstwSF -- MSTW2008 LO over PYTHIA 8's pdfdata grid, the PDF "
      "CDKS computed their b1_d with, and the one the A = 2 gate passes on "
      "(G3b 0.843 against 0.440 on the toy).  MEASURED on "
      "Li6ConvolutionB1.b1(x, 2.5): x1.847766 / x1.275961 / x0.816971 at "
      "x = 0.10 / 0.30 / 0.50 -- not a normalisation.  Needs the PYTHIA "
      "tier.\n"
      "  Ct18Nlo  LhapdfSF('CT18NLO', 0), the phase-D stand-in, kept "
      "selectable so the systematic can be quoted.  x2.221302 / x1.238833 / "
      "x1.045870 at the same points.  Needs the LHAPDF tier.\n"
      "  Custom   whatever object is in PipelineConfig.b1_unpol_sf; the enum "
      "is the provenance meta['b1_unpol'] records.\n"
      "Set it with set_b1_unpol(config, source), which builds the backend and "
      "raises here if the tier or the grid is missing -- it is NEVER silently "
      "replaced by ToyF2.")
      .value("Toy", B1UnpolSource::Toy)
      .value("Mstw", B1UnpolSource::Mstw)
      .value("Ct18Nlo", B1UnpolSource::Ct18Nlo)
      .value("Custom", B1UnpolSource::Custom);
  m.def("b1_unpol_name", &b1_unpol_name, py::arg("source"));

  py::enum_<RSource>(m, "RSource",
      "Which R = sigma_L/sigma_T is threaded into BOTH halves of the 6Li "
      "tensor weight (CLI --r-source).  The shipped observable is a RATIO -- "
      "A_zz = -(2/3) K/D_phi with K = b1 + (1-y)/(x y^2) b2 and "
      "D_phi = F1 + (1-y)/(x y^2) F2 -- whose numerator's R is "
      "Li6ConvolutionOptions.r_func and whose denominator's is "
      "InclusiveKernel.Options.r_func.  This selector is the ONE hook "
      "default_inclusive_kernel puts in both, so whichever R is named the two "
      "halves are the same choice.  Read by b1_model = Li6Convolution ONLY "
      "(PipelineConfig.validate refuses the rest, the b1_unpol rule "
      "verbatim).\n"
      "  Unset     the DEFAULT: NOTHING is installed, both hooks stay null "
      "and each consumer resolves its own r_sigma_lt.  Today bit for bit, by "
      "construction -- the branch is not taken.\n"
      "  SigmaLt   ONE r_sigma_lt object in both.  MEASURED bit-identical to "
      "Unset (0 ulp) on b1, K/D_phi, A_zz and the cos 2phi amplitude at all "
      "six standard points; it is the wiring's own test, not a physics "
      "variation.\n"
      "  R1998     ONE r1998 object in both -- the A = 2 gate's R, on the "
      "numerator AND the denominator.  MEASURED against Unset at y = 0.5, "
      "Q2 = 2.5 on the 6Li li6-convolution tensor weight K/D_phi (and so on "
      "A_zz, which is -(2/3) of it): -3.8357 % / +26.3988 % / +3.3518 % at "
      "x = 0.05 / 0.10 / 0.30, against option (ii)'s (r1998 in the numerator "
      "alone) -5.5148 % / +24.6099 % / +2.6629 %; the cos 2phi amplitude, "
      "whose numerator carries no R at all, moves by -8.3268 % / -6.7268 % / "
      "-3.1441 % where (ii) leaves it EXACTLY unchanged -- the pure "
      "denominator effect.  The weight's shift is y-DEPENDENT here and "
      "y-INDEPENDENT under (ii): at x = 0.05, Q2 = 2.5 it runs -5.4705 % at "
      "y = 0.1 to +2.3678 % at y = 0.9, straddling zero, against (ii)'s flat "
      "-5.5148 %.  And unlike --b1-unpol this one moves the unpolarised rate "
      "too (sigma_pb -0.6847 % on a 2000-event 6Li run at x_max 0.95), "
      "because the kernel's r_func reaches F1, F_L, D(y) and ToyG1 "
      "(docs/open_items/run_2026-09-06/phase_A_numbers.md sec. A1).")
      .value("Unset", RSource::Unset)
      .value("SigmaLt", RSource::SigmaLt)
      .value("R1998", RSource::R1998);
  m.def("r_source_name", &r_source_name, py::arg("source"));

  // The `set_pythia_hadronizer` arrangement (bindings.cpp header note): the
  // OPTIONAL tier's object is built HERE, where the tier is visible, and put
  // into the config.  `lipolgen_core` links neither PYTHIA nor LHAPDF, so
  // `PipelineConfig::validate()` can refuse a named backend with an empty
  // slot but cannot fill it; this is the function that fills it.  A missing
  // tier throws with the tier named, and a missing grid file throws out of
  // MstwSF/LhapdfSF naming the path -- never a quiet fallback to ToyF2.
  m.def("set_b1_unpol", [](PipelineConfig& cfg, B1UnpolSource src) {
    switch (src) {
      case B1UnpolSource::Toy:
        cfg.b1_unpol_sf.reset();
        break;
      case B1UnpolSource::Mstw:
#if LIPOLGEN_HAVE_PYTHIA8
        cfg.b1_unpol_sf = std::make_shared<const MstwSF>();
#else
        throw std::runtime_error(
            "set_b1_unpol: b1_unpol = mstw needs the OPTIONAL PYTHIA 8 tier "
            "(MstwSF reads PYTHIA's own pdfdata/mstw2008lo.00.dat), and this "
            "build has none -- lipolgen.HAVE_PYTHIA8 is False.  Rebuild with "
            "-DLIPOLGEN_WITH_PYTHIA=ON.  It is never silently replaced by "
            "ToyF2, whose A = 2 gate ratio is 0.440 against MSTW's 0.843.");
#endif
        break;
      case B1UnpolSource::Ct18Nlo:
#if LIPOLGEN_HAVE_LHAPDF
        cfg.b1_unpol_sf = std::make_shared<const LhapdfSF>("CT18NLO", 0);
#else
        throw std::runtime_error(
            "set_b1_unpol: b1_unpol = ct18nlo needs the OPTIONAL LHAPDF tier, "
            "and this build has none -- lipolgen.HAVE_LHAPDF is False.  "
            "Rebuild with -DLIPOLGEN_WITH_LHAPDF=ON and install the CT18NLO "
            "set.  It is never silently replaced by ToyF2.");
#endif
        break;
      case B1UnpolSource::Custom:
        throw std::runtime_error(
            "set_b1_unpol: B1UnpolSource.Custom names an object this function "
            "cannot build -- assign config.b1_unpol_sf = <your UnpolSF> "
            "instead, which sets Custom for you");
    }
    cfg.b1_unpol = src;
  }, py::arg("config"), py::arg("source"),
     "config.b1_unpol = source, with the matching backend built and attached "
     "to config.b1_unpol_sf.  Raises RuntimeError naming the missing tier (or "
     "the missing grid file) rather than falling back to ToyF2.");

  py::enum_<UnpolSfSource>(m, "UnpolSfSource",
      "Which UNPOLARISED structure-function backend supplies F2 -- and "
      "through it F1, F_L and the whole unpolarised rate -- to EVERY kernel "
      "a Pipeline builds: the inclusive kernel, the coherent channel that "
      "rides its cell cross sections, and the tagged struck-cluster kernel "
      "(PipelineConfig.unpol_sf; CLI --unpol-sf).\n"
      "  Toy      the DEFAULT: the kernel's own ToyF2, bit for bit what "
      "every published number was made with.  ToyF2 is labelled TOY in "
      "sf.hpp and anchored BY EYE -- 'factor-1.5 rate estimates ONLY'.\n"
      "  Mstw     MstwSF (MSTW2008 LO over PYTHIA 8's pdfdata grid).  Needs "
      "the PYTHIA tier.\n"
      "  Ct18Nlo  LhapdfSF('CT18NLO', 0).  Needs the LHAPDF tier.  READ THE "
      "GRID CLAUSE: 36 % of a shipped 6Li run's own accepted cell cross "
      "section sits below CT18NLO's Q2 = 1.677 grid floor, where LHAPDF "
      "extrapolates downward rather than freezing.\n"
      "  Custom   whatever object is in PipelineConfig.unpol_sf_obj; the "
      "enum is the provenance meta['unpol_sf'] records.\n"
      "MEASURED on 6Li inclusive at config 1: the summed accepted cell cross "
      "section is 591846.2 pb (toy), 472571.9 (ct18nlo, x0.7985), 469556.5 "
      "(mstw, x0.7934) -- and it is not a normalisation, F2p at Q2 = 10 "
      "moves x0.92 / x1.11 / x1.37 / x1.26 (ct18nlo) at x = 0.01 / 0.10 / "
      "0.30 / 0.50.\n"
      "NOT ORTHOGONAL TO PolSfSource: InclusiveKernel's default ToyG1 is "
      "built on the kernel's OWN UnpolSF, so this flag alone moves g1 too "
      "(A1 = g1/F1 moves 5-8 %).  pol_sf = Toy does NOT mean 'g1 "
      "unchanged'.\n"
      "Set it with set_unpol_sf(config, source), which builds the backend "
      "and raises here if the tier or the grid is missing -- it is NEVER "
      "silently replaced by ToyF2.")
      .value("Toy", UnpolSfSource::Toy)
      .value("Mstw", UnpolSfSource::Mstw)
      .value("Ct18Nlo", UnpolSfSource::Ct18Nlo)
      .value("Custom", UnpolSfSource::Custom);
  m.def("unpol_sf_name", &unpol_sf_name, py::arg("source"));

  py::enum_<PolSfSource>(m, "PolSfSource",
      "Which POLARISED backend supplies g1 -- and through Wandzura-Wilczek "
      "g2 -- to the INCLUSIVE and TAGGED kernels a Pipeline builds "
      "(PipelineConfig.pol_sf; CLI --pol-sf).\n"
      "  Toy       the DEFAULT: the kernel's own ToyG1, built on the "
      "kernel's own UnpolSF and r_func.\n"
      "  NnpdfPol  LhapdfG1('NNPDFpol11_100', 0) -- the only POLARISED set "
      "installed in this tree's LHAPDF store, and already LhapdfG1's own "
      "default.  Needs the LHAPDF tier.\n"
      "  Custom    whatever object is in PipelineConfig.pol_sf_obj.\n"
      "MEASURED on 6Li per-nucleon g1A (A_par tracks it to better than "
      "0.1 % at y = 0.5): nnpdfpol/toy = 0.669 / 1.114 / 1.349 / 0.982 / "
      "0.636 at (x, Q2) = (0.01, 2.5) / (0.10, 10) / (0.20, 10) / "
      "(0.50, 25) / (0.70, 50) -- x0.64 to x1.35, and not monotone.\n"
      "THE NEUTRON IS A SIGN, not a factor: ToyG1's a1n crosses zero near "
      "x ~ 0.25 and is POSITIVE above it while NNPDFpol1.1's g1n stays "
      "negative to x ~ 0.6, so the shipped toy g1n has the WRONG SIGN over "
      "roughly 0.25 < x < 0.6.  That matters most on a neutron-tagged run "
      "(d + p tagging), least on isoscalar 6Li.\n"
      "SCOPE -- IT DOES NOT REACH EVERY KERNEL, and unlike UnpolSfSource it "
      "never could; this is the C++ header's own caveat "
      "(include/lipolgen/pipeline.hpp, PolSfSource), which this docstring "
      "dropped until 2026-09-05.  g1 enters through exactly one product, "
      "lam_e * P_e * (m/J) * cos(theta_S) * A_par, so the selector is read "
      "only where the CHANNEL evaluates g1 AND the FILL carries "
      "lam_e * P_e != 0.  It does NOT reach CoherentLi6, whose rate is "
      "spin-independent, and it does NOT reach any unpolarised-beam plan -- "
      "tensor_thirds_plan, transverse_tensor_plan and tensor_flip_plan build "
      "every category at lam_e = 0, and tensor-thirds is the CLI's own "
      "default, so at defaults --pol-sf is read on NO channel.  Neither is "
      "refused; both are LABELLED (pol_sf_is_read(config, plan), "
      "meta['pol_sf'] / meta['pol_sf_reach'], and the pol_sf row of "
      "Pipeline.knob_provenance).  'Every kernel the pipeline builds' is "
      "true of --unpol-sf alone.\n"
      "Set it with set_pol_sf(config, source).")
      .value("Toy", PolSfSource::Toy)
      .value("NnpdfPol", PolSfSource::NnpdfPol)
      .value("Custom", PolSfSource::Custom);
  m.def("pol_sf_name", &pol_sf_name, py::arg("source"));

  m.def("unpol_sf_grid_report", [](const Pipeline& p) {
    const std::pair<double, double> g = unpol_sf_grid_report(p);
    return py::make_tuple(g.first, g.second);
  }, py::arg("pipeline"),
     "(q2_min, below_frac) of the run's selected unpolarised backend: the "
     "loaded grid's OWN lower Q2 edge, and the fraction of THIS run's own "
     "accepted RATE that sits below it, where LHAPDF extrapolates downward "
     "rather than freezing.  The denominator is the CHANNEL's own per-cell "
     "rate (Pipeline.cell_rate_weights_pb), which on the coherent channel is "
     "sigma_cell * f_coh(x) and not the inclusive cells: measured on the "
     "shipped 6Li ct18nlo run, 0.447460 coherent against 0.361791 "
     "inclusive.  (0, 0) for the toy backend, which has no grid; (nan, nan) "
     "when the backend reports no floor (MstwSF keeps PYTHIA's qsqmin "
     "private; a Custom object or a caller-supplied kernel reports "
     "nothing).  ONE definition, read by both the run banner and "
     "meta['unpol_sf_below_grid_frac'].");

  // The polarised selector's REACH (pipeline.hpp).  Exposed in all three
  // pieces for the same reason the rank-2 trio is: the CLI needs the
  // predicate to decide whether to print, the sentence to print, and the
  // label to compare meta against; a pytest needs all three to check that
  // meta and the banner say the same true thing.
  m.def("pol_sf_is_read", &pol_sf_is_read, py::arg("config"),
        py::arg("plan"),
        "Does --pol-sf reach the RATE of THIS RUN?  TWO axes, and the second "
        "one is the shipped default.  False on the coherent channel -- the "
        "yield is f_coh(x) times the UNPOLARIZED cell cross sections and the "
        "tensor signal is the recoil azimuth's 1 + c2 cos 2(phi_t - phi_S), "
        "so nothing evaluates g1 -- and false under any UNPOLARISED-BEAM "
        "PLAN, because InclusiveKernel::amplitudes multiplies g1 by "
        "lam_e * P_e and the three tensor plans build every category at "
        "lam_e = 0, pe = 0 (helicity-flip at --pe 0 likewise).  Measured: "
        "--pol-sf nnpdfpol is bit-identical to toy on four channels under "
        "tensor-thirds, the CLI's own default plan.  The kernel still "
        "CARRIES the selected g1_model; the run never reads it.  --unpol-sf, "
        "by contrast, reaches every run.");
  m.def("pol_sf_reach_report", &pol_sf_reach_report, py::arg("config"),
        py::arg("plan"),
        "One sentence saying what --pol-sf reached on this run.  ONE "
        "definition, read by meta['pol_sf_reach'], by the pol_sf row of "
        "Pipeline.knob_provenance and printed by the CLI run banner.");
  m.def("pol_sf_unread_label", &pol_sf_unread_label, py::arg("config"),
        py::arg("plan"),
        "The label meta['pol_sf'] carries INSTEAD of a backend name where "
        "pol_sf_is_read is false -- the rank2_none_label rule for the "
        "polarised selector: a run whose g1 was never evaluated may not "
        "record 'nnpdfpol', and may not record 'toy' either.  Two shapes, "
        "one per axis: 'not read on channel coherent-6Li' and 'not read "
        "under this run's fill'.");
  m.def("plan_has_beam_helicity", &plan_has_beam_helicity, py::arg("plan"),
        "Does any category of this plan carry lam_e != 0 AND pe != 0 AND a "
        "non-zero population at some m != 0?  That product is the whole of "
        "what InclusiveKernel::amplitudes multiplies g1 by, so it is the "
        "exact condition under which a run evaluates a polarised structure "
        "function.  ONE definition, read by pol_sf_is_read and by the pe row "
        "of Pipeline.knob_provenance.");

  // The 7Li rank-2 zero, out loud (pipeline.hpp).  Both are exposed because
  // the CLI needs the PREDICATE to decide whether to print and the SENTENCE
  // to print, and a pytest needs both to check that meta and banner agree.
  m.def("inclusive_rank2_is_empty", &inclusive_rank2_is_empty,
        py::arg("config"),
        "True when this config's INCLUSIVE kernel has a rank-2 (tensor) "
        "sector that nothing fills -- i.e. isotope 7Li on the inclusive "
        "channel with no caller-supplied kernel.  default_inclusive_kernel "
        "fills b1_func/delta_func for spin 1 ONLY, so at spin 3/2 "
        "b1_32 = b2_32 = delta_32 = 0 and the tensor term of the rate, the "
        "cos 2phi amplitude and A_zz are IDENTICALLY zero -- not small.  "
        "False on a caller-supplied kernel, false off the inclusive channel "
        "(a tagged 7Li run DOES carry the alpha-t alignment, in the event "
        "weight), and false at spin 1/2 or 0, where there is no rank-2 "
        "sector at all.");
  m.def("rank2_none_label", &rank2_none_label,
        "The label meta['b1_model'] and meta['b1_unpol'] carry instead of a "
        "backend name when inclusive_rank2_is_empty(config) holds -- a "
        "'miller' there would name a backend that did not run.  The two "
        "SCALES b1_band_scale / b1_alpha_d_dwave_weight stay numeric and at "
        "their defaults, because validate()'s Miller branch already refuses "
        "any other value on that path.");
  m.def("rank2_input_report", &rank2_input_report, py::arg("config"),
        py::arg("plan"),
        "One sentence saying what filled -- or did not fill -- this run's "
        "rank-2 slots.  ONE definition, read by meta['rank2_input'] and "
        "printed by the CLI run banner, so the claim cannot drift between "
        "them.  The plan is read for one clause: a fill that carries a "
        "rank-2 moment (RunPlan.pzz_true, which is T and not P_zz at "
        "J = 3/2) is named, because that is the run that asked for the "
        "sector that is not there.");

  // The `set_b1_unpol` arrangement, for the two general selectors: the
  // OPTIONAL tier's object is built HERE, where the tier is visible.
  m.def("set_unpol_sf", [](PipelineConfig& cfg, UnpolSfSource src) {
    switch (src) {
      case UnpolSfSource::Toy:
        cfg.unpol_sf_obj.reset();
        break;
      case UnpolSfSource::Mstw:
#if LIPOLGEN_HAVE_PYTHIA8
        cfg.unpol_sf_obj = std::make_shared<const MstwSF>();
#else
        throw std::runtime_error(
            "set_unpol_sf: unpol_sf = mstw needs the OPTIONAL PYTHIA 8 tier "
            "(MstwSF reads PYTHIA's own pdfdata/mstw2008lo.00.dat), and this "
            "build has none -- lipolgen.HAVE_PYTHIA8 is False.  Rebuild with "
            "-DLIPOLGEN_WITH_PYTHIA=ON.  It is never silently replaced by "
            "ToyF2, which is 21 % away on the shipped 6Li rate.");
#endif
        break;
      case UnpolSfSource::Ct18Nlo:
#if LIPOLGEN_HAVE_LHAPDF
        cfg.unpol_sf_obj = std::make_shared<const LhapdfSF>("CT18NLO", 0);
#else
        throw std::runtime_error(
            "set_unpol_sf: unpol_sf = ct18nlo needs the OPTIONAL LHAPDF "
            "tier, and this build has none -- lipolgen.HAVE_LHAPDF is False. "
            " Rebuild with -DLIPOLGEN_WITH_LHAPDF=ON and install the CT18NLO "
            "set.  It is never silently replaced by ToyF2.");
#endif
        break;
      case UnpolSfSource::Custom:
        throw std::runtime_error(
            "set_unpol_sf: UnpolSfSource.Custom names an object this "
            "function cannot build -- assign config.unpol_sf_obj = <your "
            "UnpolSF> instead, which sets Custom for you");
    }
    cfg.unpol_sf = src;
  }, py::arg("config"), py::arg("source"),
     "config.unpol_sf = source, with the matching backend built and attached "
     "to config.unpol_sf_obj.  Raises RuntimeError naming the missing tier "
     "(or the missing grid file) rather than falling back to ToyF2.");

  m.def("set_pol_sf", [](PipelineConfig& cfg, PolSfSource src) {
    switch (src) {
      case PolSfSource::Toy:
        cfg.pol_sf_obj.reset();
        break;
      case PolSfSource::NnpdfPol:
#if LIPOLGEN_HAVE_LHAPDF
        cfg.pol_sf_obj = std::make_shared<const LhapdfG1>("NNPDFpol11_100", 0);
#else
        throw std::runtime_error(
            "set_pol_sf: pol_sf = nnpdfpol needs the OPTIONAL LHAPDF tier, "
            "and this build has none -- lipolgen.HAVE_LHAPDF is False.  "
            "Rebuild with -DLIPOLGEN_WITH_LHAPDF=ON and install the "
            "NNPDFpol11_100 set.  It is never silently replaced by ToyG1, "
            "whose g1n has the WRONG SIGN over 0.25 < x < 0.6.");
#endif
        break;
      case PolSfSource::Custom:
        throw std::runtime_error(
            "set_pol_sf: PolSfSource.Custom names an object this function "
            "cannot build -- assign config.pol_sf_obj = <your PolSF> "
            "instead, which sets Custom for you");
    }
    cfg.pol_sf = src;
  }, py::arg("config"), py::arg("source"),
     "config.pol_sf = source, with the matching backend built and attached "
     "to config.pol_sf_obj.  Raises RuntimeError naming the missing tier "
     "rather than falling back to ToyG1.");

  py::class_<CiofiSimulaOptions>(m, "CiofiSimulaOptions",
      "Configuration of CiofiSimulaTriton.  Everything here is a documented "
      "CHOICE; the transcribed CS coefficients are not options "
      "(triton_sf.hpp).")
      .def(py::init<>())
      .def_readwrite("k_max", &CiofiSimulaOptions::k_max)
      .def_readwrite("n_grid", &CiofiSimulaOptions::n_grid)
      .def_readwrite("q_max", &CiofiSimulaOptions::q_max)
      .def_readwrite("beta", &CiofiSimulaOptions::beta)
      .def_readwrite("kappa_nn", &CiofiSimulaOptions::kappa_nn)
      .def_readwrite("kappa_pn", &CiofiSimulaOptions::kappa_pn)
      .def_readwrite("n1_scale", &CiofiSimulaOptions::n1_scale)
      .def_readwrite("proton_n1_only", &CiofiSimulaOptions::proton_n1_only);

  py::class_<TritonSpectralFunction,
             std::shared_ptr<TritonSpectralFunction>>(m,
      "TritonSpectralFunction",
      "The triton spectral-function interface S_N(k, E) (triton_sf.hpp); "
      "CiofiSimulaTriton is the analytic implementation.")
      .def("n_of_k", &TritonSpectralFunction::n_of_k, py::arg("k"),
           py::arg("pdg"))
      .def("p_two_body", &TritonSpectralFunction::p_two_body, py::arg("k"),
           py::arg("pdg"))
      .def("e_rel_density", &TritonSpectralFunction::e_rel_density,
           py::arg("k"), py::arg("e_rel"), py::arg("pdg"))
      .def("e_max", &TritonSpectralFunction::e_max, py::arg("pdg"));

  py::class_<CiofiSimulaTriton, TritonSpectralFunction,
             std::shared_ptr<CiofiSimulaTriton>>(m, "CiofiSimulaTriton",
      "Ciofi degli Atti-Simula (PRC 53 (1996) 1689) triton spectral "
      "function: n_0/n_1 as transcribed, the 2-body/3-body split from their "
      "RATIO, the continuum pair split at its own virtual-state pole.")
      .def(py::init<CiofiSimulaOptions>(),
           py::arg("opt") = CiofiSimulaOptions())
      .def_property_readonly("options", &CiofiSimulaTriton::options)
      .def_property_readonly("s0", &CiofiSimulaTriton::s0)
      .def_property_readonly("s1", &CiofiSimulaTriton::s1)
      .def("n0_cs", &CiofiSimulaTriton::n0_cs, py::arg("k_gev"))
      .def("n1_cs", &CiofiSimulaTriton::n1_cs, py::arg("k_gev"));

  py::class_<BreakupOptions>(m, "BreakupOptions",
      "Configuration of the T1 cluster breakup.  `beta` and `f2` are "
      "overwritten by the Pipeline with the run's own values.")
      .def(py::init<>())
      .def_readwrite("beta", &BreakupOptions::beta)
      .def_readwrite("p_d", &BreakupOptions::p_d)
      .def_readwrite("source", &BreakupOptions::source,
                     "ClusterWaveSource of the DEUTERON the T1 breakup "
                     "resolves into.  Hulthen (default) = the analytic S + D "
                     "pair at p_d; VmcAV18 = the exact AV18 fdeut deuteron at "
                     "its own P_D = 0.057600.  Pipeline forwards "
                     "PipelineConfig.cluster_wave here so the spin the record "
                     "carries and the g1 the rate uses are one deuteron "
                     "(C5.5b).  The triton path is not affected.")
      .def_readwrite("kappa_nn", &BreakupOptions::kappa_nn)
      .def_readwrite("k_max", &BreakupOptions::k_max)
      .def_readwrite("nk", &BreakupOptions::nk)
      .def_readwrite("nc", &BreakupOptions::nc)
      .def_readwrite("n_grid", &BreakupOptions::n_grid)
      .def_property("f2",
          [](const BreakupOptions& o) {
            return std::const_pointer_cast<UnpolSF>(o.f2);
          },
          [](BreakupOptions& o, std::shared_ptr<UnpolSF> f) {
            o.f2 = std::move(f);
          })
      .def_property("triton_sf",
          [](const BreakupOptions& o) {
            return std::const_pointer_cast<TritonSpectralFunction>(o.triton_sf);
          },
          [](BreakupOptions& o, std::shared_ptr<TritonSpectralFunction> t) {
            o.triton_sf = std::move(t);
          });

  py::enum_<OpticsChoice>(m, "OpticsChoice")
      .value("YellowReportHighAcceptance",
             OpticsChoice::YellowReportHighAcceptance)
      .value("YellowReportHighDivergence",
             OpticsChoice::YellowReportHighDivergence)
      .value("Tagging", OpticsChoice::Tagging)
      .value("TaggingLegacyLevers", OpticsChoice::TaggingLegacyLevers)
      .value("Custom", OpticsChoice::Custom);
  m.def("tagging_optics", &tagging_optics, py::arg("ion_name"),
        py::arg("p_per_nucleon"), py::arg("legacy_levers") = false,
        py::arg("n_sigma") = 10.0);
  m.def("optics_for", &optics_for, py::arg("choice"), py::arg("ion_name"),
        py::arg("p_per_nucleon"), py::arg("n_sigma") = 10.0);

  py::class_<StruckClusterOptions>(m, "StruckClusterOptions")
      .def(py::init<>())
      .def_readwrite("inclusive_b1", &StruckClusterOptions::inclusive_b1)
      .def_readwrite("delta_func", &StruckClusterOptions::delta_func)
      .def_readwrite("scenario", &StruckClusterOptions::scenario)
      .def_readwrite("grid", &StruckClusterOptions::grid)
      .def_readwrite("with_perp", &StruckClusterOptions::with_perp)
      .def_property("f2_source",
          [](const StruckClusterOptions& o) {
            return std::const_pointer_cast<UnpolSF>(o.f2_source);
          },
          [](StruckClusterOptions& o, std::shared_ptr<UnpolSF> u) {
            o.f2_source = std::move(u);
          },
          "InclusiveKernel.Options.f2_source of the struck cluster's kernel "
          "-- the tagged channels' ONLY unpolarised injection point.  None "
          "(the default) is the kernel class's own ToyF2, bit for bit.  "
          "Pipeline fills it from PipelineConfig.unpol_sf_obj when it is "
          "left empty.")
      .def_property("g1_model",
          [](const StruckClusterOptions& o) {
            return std::const_pointer_cast<PolSF>(o.g1_model);
          },
          [](StruckClusterOptions& o, std::shared_ptr<PolSF> g) {
            o.g1_model = std::move(g);
          },
          "InclusiveKernel.Options.g1_model of the struck cluster's kernel.  "
          "None (the default) is the kernel's own ToyG1 on its own UnpolSF.  "
          "Pipeline fills it from PipelineConfig.pol_sf_obj when it is left "
          "empty.");

  py::class_<PipelineConfig>(m, "PipelineConfig")
      .def(py::init<>())
      .def_readwrite("isotope", &PipelineConfig::isotope)
      .def_readwrite("beam_config", &PipelineConfig::beam_config)
      .def_readwrite("channel", &PipelineConfig::channel)
      .def_readwrite("scenario", &PipelineConfig::scenario)
      .def_readwrite("grid", &PipelineConfig::grid)
      .def_property("kernel",
          [](const PipelineConfig& c) {
            return std::const_pointer_cast<InclusiveKernel>(c.kernel);
          },
          [](PipelineConfig& c, std::shared_ptr<InclusiveKernel> k) {
            c.kernel = std::move(k);
          })
      .def_readwrite("b1_model", &PipelineConfig::b1_model,
                     "Inclusive b1 backend (B1Model); Miller is the default "
                     "and is today bit for bit.  validate() refuses a "
                     "non-Miller model together with a caller-supplied "
                     "kernel, and refuses Li6Convolution off the inclusive "
                     "channel or off 6Li.")
      .def_readwrite("b1_band_scale", &PipelineConfig::b1_band_scale,
                     "The MANDATORY 100 % band: multiplies the whole b1 of "
                     "the Cdks and Li6Convolution backends.  Run 0/1/2 and "
                     "quote the envelope -- never a single row.")
      .def_readwrite("b1_alpha_d_dwave_weight",
                     &PipelineConfig::b1_alpha_d_dwave_weight,
                     "Li6ConvolutionOptions.w_alpha_d_dwave: the knob on "
                     "terms (2d) AND (2a) together, the SHAPE variant of the "
                     "band.")
      .def_readwrite("b1_unpol", &PipelineConfig::b1_unpol,
                     "B1UnpolSource: which UNPOLARISED backend "
                     "b1_model = Li6Convolution folds its F1 against.  Toy is "
                     "the default and is today bit for bit.  Assigning the "
                     "enum here does NOT build the backend -- use "
                     "set_b1_unpol(config, source), or validate() will refuse "
                     "the named-but-empty combination.")
      .def_property("b1_unpol_sf",
          [](const PipelineConfig& c) {
            return std::const_pointer_cast<UnpolSF>(c.b1_unpol_sf);
          },
          [](PipelineConfig& c, std::shared_ptr<UnpolSF> u) {
            // Assigning an object IS B1UnpolSource.Custom, and clearing it is
            // back to Toy: the two fields are the provenance and the
            // realisation of one choice (the optics_choice / optics
            // arrangement), so this setter keeps them from drifting apart.
            // `set_b1_unpol` is the other direction of the same invariant.
            c.b1_unpol = u ? B1UnpolSource::Custom : B1UnpolSource::Toy;
            c.b1_unpol_sf = std::move(u);
          },
          "The UnpolSF b1_unpol names.  Assigning one sets b1_unpol = Custom; "
          "assigning None clears it back to Toy.  For the named backends use "
          "set_b1_unpol(config, source), which builds them.")
      .def_readwrite("r_source", &PipelineConfig::r_source,
                     "RSource: the ONE R = sigma_L/sigma_T threaded into BOTH "
                     "Li6ConvolutionOptions.r_func (the tensor weight's "
                     "numerator) and InclusiveKernel.Options.r_func (its "
                     "denominator, and the unpolarised rate).  Unset is the "
                     "default and installs nothing, which is today bit for "
                     "bit.  Read by b1_model = Li6Convolution ONLY; "
                     "validate() refuses any other value elsewhere.")
      .def_readwrite("unpol_sf", &PipelineConfig::unpol_sf,
                     "UnpolSfSource: which UNPOLARISED backend EVERY kernel "
                     "this config builds takes its F2 from -- inclusive, "
                     "coherent and tagged alike.  Toy is the default and is "
                     "today bit for bit.  Assigning the enum here does NOT "
                     "build the backend -- use set_unpol_sf(config, source), "
                     "or validate() will refuse the named-but-empty "
                     "combination.")
      .def_property("unpol_sf_obj",
          [](const PipelineConfig& c) {
            return std::const_pointer_cast<UnpolSF>(c.unpol_sf_obj);
          },
          [](PipelineConfig& c, std::shared_ptr<UnpolSF> u) {
            // The provenance and the realisation are one choice; this setter
            // is what keeps them from drifting (the b1_unpol_sf arrangement).
            c.unpol_sf = u ? UnpolSfSource::Custom : UnpolSfSource::Toy;
            c.unpol_sf_obj = std::move(u);
          },
          "The UnpolSF unpol_sf names.  Assigning one sets unpol_sf = "
          "Custom; assigning None clears it back to Toy.  For the named "
          "backends use set_unpol_sf(config, source), which builds them.")
      .def_readwrite("pol_sf", &PipelineConfig::pol_sf,
                     "PolSfSource: which POLARISED backend the INCLUSIVE and "
                     "TAGGED kernels take their g1 (and, through "
                     "Wandzura-Wilczek, their g2) from -- NOT every kernel: "
                     "the coherent channel evaluates no g1, and neither does "
                     "any unpolarised-beam plan (tensor-thirds included, "
                     "which is the CLI default), so the run LABELS the "
                     "selector there rather than crediting it (PolSfSource's "
                     "own SCOPE note; pol_sf_is_read(config, plan)).  Toy is "
                     "the default and is today bit for bit.  Use "
                     "set_pol_sf(config, source) to build the backend.")
      .def_property("pol_sf_obj",
          [](const PipelineConfig& c) {
            return std::const_pointer_cast<PolSF>(c.pol_sf_obj);
          },
          [](PipelineConfig& c, std::shared_ptr<PolSF> g) {
            c.pol_sf = g ? PolSfSource::Custom : PolSfSource::Toy;
            c.pol_sf_obj = std::move(g);
          },
          "The PolSF pol_sf names.  Assigning one sets pol_sf = Custom; "
          "assigning None clears it back to Toy.")
      .def_readwrite("with_virtual_photon", &PipelineConfig::with_virtual_photon)
      .def_readwrite("seed", &PipelineConfig::seed)
      .def_readwrite("run", &PipelineConfig::run)
      .def_readwrite("lumi_pb", &PipelineConfig::lumi_pb)
      .def_readwrite("n_events", &PipelineConfig::n_events)
      .def_readwrite("poisson", &PipelineConfig::poisson)
      .def_readwrite("optics_choice", &PipelineConfig::optics_choice)
      .def_readwrite("optics", &PipelineConfig::optics)
      .def_readwrite("n_sigma", &PipelineConfig::n_sigma)
      .def_readwrite("pot_config", &PipelineConfig::pot_config)
      .def_readwrite("cluster_beta", &PipelineConfig::cluster_beta)
      .def_readwrite("p_d", &PipelineConfig::p_d)
      .def_readwrite("cluster_vmc_mc_sigma",
                     &PipelineConfig::cluster_vmc_mc_sigma,
                     "ANL VMC Monte Carlo band in units of the tables' own "
                     "printed 1 sigma, fully correlated across k.  0 = today "
                     "bit for bit; validate() refuses it off a lithium "
                     "alpha-tag channel with cluster_wave = VmcAV18.")
      .def_readwrite("cluster_wave", &PipelineConfig::cluster_wave,
                     "ClusterWaveSource for the lithium alpha-tag channels; "
                     "Hulthen by default (bit-compatible), VmcAV18 swaps in "
                     "the ANL VMC tables and ignores cluster_beta / p_d.")
      .def_readwrite("struck", &PipelineConfig::struck)
      .def_readwrite("tier", &PipelineConfig::tier,
                     "Fidelity tier of the tagged final state; Tier.T1 by "
                     "default (the struck cluster is resolved into a nucleon "
                     "plus partner spectators).")
      .def_readwrite("breakup", &PipelineConfig::breakup)
      .def_readwrite("triton_sf", &PipelineConfig::triton_sf,
                     "TritonSfChoice for the 7Li alpha tag's T1 triton "
                     "breakup; Hulthen by default (bit-compatible), "
                     "CiofiSimula swaps in the three-channel spectral "
                     "function built at the run's own cluster_beta.")
      .def_readwrite("fsi", &PipelineConfig::fsi,
                     "PipelineFsi: FSI of the DIS debris with the tagged "
                     "spectator, as a per-event WEIGHT on Event.weight (never "
                     "a shift).  Off by default = today's PWIA bit for bit; "
                     "tagged channels only.")
      .def_readwrite("fsi_sigma_mb", &PipelineConfig::fsi_sigma_mb,
                     "sigma_XN [mb] the FSI weight is built at.  40 = free "
                     "hadron; the documented band is 20-40 mb -- band it, "
                     "never quote one row alone.")
      .def_readwrite("rc", &PipelineConfig::rc,
                     "RcMode: tensor-sector radiative corrections as OPT-IN, "
                     "WEIGHT-ONLY families (rc.hpp).  They land on "
                     "Event.rc_weights and on NOTHING else -- never "
                     "Event.weight, never a four-vector, never a random "
                     "number -- so Off (the default) is today bit for bit.")
      .def_readwrite("rc_options", &PipelineConfig::rc_options,
                     "RcOptions: the RC knobs.  There is NO mode on it -- "
                     "PipelineConfig.rc is the single source of truth.")
      .def_readwrite("coherent", &PipelineConfig::coherent)
      .def_readwrite("coherent_t_max", &PipelineConfig::coherent_t_max)
      .def_readwrite("coherent_xpom", &PipelineConfig::coherent_xpom)
      .def_readwrite("apply_optics_lumi_fraction",
                     &PipelineConfig::apply_optics_lumi_fraction)
      .def_readwrite("coherent_weighted_azimuth",
                     &PipelineConfig::coherent_weighted_azimuth)
      .def_readwrite("hadronizer", &PipelineConfig::hadronizer)
      .def("validate", &PipelineConfig::validate);

  // The FULL form, defaulted so that `default_inclusive_kernel(ion)` is the
  // one-argument overload it has always been -- i.e. (Miller, 1, 1, null),
  // which `validation/reference/b1_default_li6.json` pins at rtol 1e-12.
  // `b1_unpol` is the object `PipelineConfig::b1_unpol_sf` carries; null is
  // the kernel's own ToyF2, shared as ONE object.
  m.def("default_inclusive_kernel", [](const Ion& ion, B1Model model,
                                       double band_scale, double w_alpha_d,
                                       std::shared_ptr<UnpolSF> b1_unpol,
                                       RSource r_source) {
    return std::const_pointer_cast<InclusiveKernel>(
        default_inclusive_kernel(ion, model, band_scale, w_alpha_d,
                                 std::move(b1_unpol), nullptr, nullptr,
                                 r_source));
  }, py::arg("ion"), py::arg("model") = B1Model::Miller,
     py::arg("band_scale") = 1.0, py::arg("w_alpha_d") = 1.0,
     py::arg("b1_unpol") = nullptr,
     py::arg("r_source") = RSource::Unset);

  m.def("momentum_residual", &momentum_residual, py::arg("ev"));
  m.def("momentum_scale", &momentum_scale, py::arg("ev"));
  m.def("charge_residual", &charge_residual, py::arg("ev"));
  m.def("route_of", &route_of, py::arg("ev"), py::arg("optics"),
        py::arg("pot_config") = "18x275");
  m.def("rp_tagged", &rp_tagged, py::arg("ev"), py::arg("optics"),
        py::arg("pot_config") = "18x275");
  m.def("struck_cluster_of", [](const Event& ev, const TaggedChannel& ch)
                                 -> py::object {
    StruckCluster sc;
    if (!struck_cluster_of(ev, ch, sc)) return py::none();
    return py::cast(sc);
  }, py::arg("ev"), py::arg("channel"));

  // ------------------------------------------- the knob-provenance table
  py::enum_<KnobStatus>(m, "KnobStatus",
      "What a run did with one user-settable knob (pipeline.hpp).")
      .value("Read", KnobStatus::Read,
             "the run consults it: some quantity it computes is a function "
             "of this knob")
      .value("NotRead", KnobStatus::NotRead,
             "nothing the run computes is a function of it; meta writes the "
             "reason in place of the value")
      .value("Refused", KnobStatus::Refused,
             "the AXIS is closed on this run -- validate() throws on any "
             "value but the one shown");
  m.def("knob_status_name", &knob_status_name, py::arg("status"),
        "'read', 'not-read' or 'refused' -- ONE spelling, shared by the npz "
        "meta block, the CLI banner and the tests.");
  py::class_<KnobProvenance>(m, "KnobProvenance",
      "One row of Pipeline.knob_provenance: what this run did with one "
      "user-settable knob, and why.")
      .def(py::init<>())
      .def_readwrite("name", &KnobProvenance::name)
      .def_readwrite("flag", &KnobProvenance::flag)
      .def_readwrite("value", &KnobProvenance::value)
      .def_readwrite("status", &KnobProvenance::status)
      .def_readwrite("reason", &KnobProvenance::reason)
      .def_readwrite("label", &KnobProvenance::label)
      .def_readwrite("at_default", &KnobProvenance::at_default)
      .def_property_readonly("meta_value", &KnobProvenance::meta_value,
        "What a meta key carries for this row: the bare value where the knob "
        "RAN or where its axis is refused (so the value is necessarily the "
        "default), and the reason sentence where it did not.")
      .def("__repr__", [](const KnobProvenance& r) {
        return "<KnobProvenance " + r.name + " = " + r.value + " (" +
               knob_status_name(r.status) + ")>";
      });
  py::class_<KnobRunContext>(m, "KnobRunContext",
      "The facts a Pipeline cannot see because they are resolved one tier "
      "up: the T2 bridge (which lives in the OPTIONAL PYTHIA tier) and what "
      "the caller asked the run-plan factory for (a RunPlan records its "
      "MOMENTS, not which flags produced them).  Rows the context cannot "
      "decide are omitted rather than guessed.")
      .def(py::init<>())
      .def_readwrite("plan_name", &KnobRunContext::plan_name)
      .def_readwrite("pz", &KnobRunContext::pz)
      .def_readwrite("pzz", &KnobRunContext::pzz)
      .def_readwrite("pe", &KnobRunContext::pe)
      .def_readwrite("rel_lumi_offset", &KnobRunContext::rel_lumi_offset)
      .def_readwrite("pzz_mode", &KnobRunContext::pzz_mode)
      .def_readwrite("t2_bound", &KnobRunContext::t2_bound)
      .def_readwrite("t2_pomeron", &KnobRunContext::t2_pomeron)
      .def_readwrite("pom_set", &KnobRunContext::pom_set)
      .def_readwrite("pom_rescale", &KnobRunContext::pom_rescale);

  py::class_<Pipeline>(m, "Pipeline",
      "One configured run: PipelineConfig + RunPlan -> Event records.")
      .def(py::init([](PipelineConfig cfg, RunPlan plan) {
        py::gil_scoped_release unlock;
        return std::make_unique<Pipeline>(std::move(cfg), std::move(plan));
      }), py::arg("config"), py::arg("plan"))
      .def_property_readonly("config", &Pipeline::config)
      .def_property_readonly("plan", &Pipeline::plan)
      .def_property_readonly("beam_config", &Pipeline::beam_config)
      .def_property_readonly("optics", &Pipeline::optics)
      .def_property_readonly("pot_config", &Pipeline::pot_config)
      .def_property_readonly("dis_sampler", &Pipeline::dis_sampler,
                             py::return_value_policy::reference_internal)
      .def_property_readonly("cell_rate_weights_pb", [](const Pipeline& p) {
        return copy_array(p.cell_rate_weights_pb());
      },
        "The per-cell weights of THIS run's own accepted rate [pb], in the "
        "cell order of dis_sampler.cell_xsec_pb / .x_cells / .q2_cells.  The "
        "same vector as dis_sampler.cell_xsec_pb on the inclusive and tagged "
        "channels; on the COHERENT one it is sigma_cell * f_coh(x) over the "
        "cells that admit a diffractive mass and 0 on the ones that do not, "
        "so it sums to sigma_pb() exactly.  It is the denominator of any "
        "'how much of this run sits in region R' fraction -- see "
        "unpol_sf_grid_report.  The spin modulation (1 + w_avg) is NOT in "
        "it: these are the spin-blind weights.")
      .def("knob_provenance", &Pipeline::knob_provenance,
           py::arg("context") = KnobRunContext(),
        "EVERY user-settable knob of this run, with what the run did with "
        "it: a list of KnobProvenance rows {name, flag, value, status, "
        "reason, at_default}.  THE mechanism the rule 'a knob that did not "
        "run may not be recorded as if it had' is enforced by -- the npz "
        "meta['knob_provenance'] block, the CLI banner's provenance block "
        "and python/tests/test_knob_provenance.py all read THIS table and "
        "nothing else.  `context` carries what the core cannot see "
        "(KnobRunContext).")
      .def_property_readonly("tagged_model", &Pipeline::tagged_model,
                             py::return_value_policy::reference_internal)
      .def_property_readonly("tagged_channel", &Pipeline::tagged_channel,
                             py::return_value_policy::reference_internal)
      .def_property_readonly("tier", &Pipeline::tier,
                             "The tier this run actually writes (Tier.T0 on "
                             "channels with no struck cluster).")
      .def_property_readonly("fsi_weight", &Pipeline::fsi_weight,
                             py::return_value_policy::reference_internal,
                             "The run's FSI weight model (GlauberFsiWeight), "
                             "or None when PipelineConfig.fsi is Off -- "
                             "print its sigma_eff_mb and survival().")
      .def_property_readonly("rc_model", &Pipeline::rc_model,
                             py::return_value_policy::reference_internal,
                             "The run's RC weight model (RcModel), or None "
                             "when PipelineConfig.rc is Off.  Print its "
                             "delta(x), ff_provenance and "
                             "clipped_fraction_by_y.")
      .def("sigma_per_category_pb", [](const Pipeline& p) {
        return copy_array(p.sigma_per_category_pb());
      })
      .def("sigma_pb", &Pipeline::sigma_pb)
      .def("optics_lumi_factor", &Pipeline::optics_lumi_factor)
      .def("lumi_per_category_pb", [](const Pipeline& p) {
        return copy_array(p.lumi_per_category_pb());
      })
      .def("counts", [](const Pipeline& p) {
        std::vector<std::int64_t> c;
        for (std::uint64_t v : p.counts()) c.push_back(static_cast<std::int64_t>(v));
        return move_array(std::move(c));
      })
      .def("size", &Pipeline::size)
      .def("__len__", [](const Pipeline& p) { return p.size(); })
      .def("category_of", &Pipeline::category_of, py::arg("index"))
      .def("event", [](const Pipeline& p, std::uint64_t i) {
        return p.event(i);
      }, py::arg("index"))
      .def("generate_range", [](const Pipeline& p, std::uint64_t a,
                                std::uint64_t b) {
        std::vector<Event> out;
        {
          py::gil_scoped_release unlock;
          p.generate_range(a, b, out);
        }
        return out;
      }, py::arg("first"), py::arg("last"),
         "list of Event records for the index range [first, last)")
      .def("generate", [](const Pipeline& p, std::uint64_t n, bool events,
                          unsigned nthreads, std::size_t chunk,
                          const KnobRunContext& ctx) {
        return pipeline_generate(p, n, events, nthreads, chunk, ctx);
      }, py::arg("n") = 0, py::arg("events") = false,
         py::arg("nthreads") = 1, py::arg("chunk") = 4096,
         py::arg("context") = KnobRunContext(),
         "Generate `n` events (0 = the whole run) and return a columnar dict\n"
         "of numpy arrays; `events=True` adds the Event records under\n"
         "'events'.  The GIL is released around the generation loop.\n"
         "`context` is the KnobRunContext the meta['knob_provenance'] block\n"
         "is written from -- what the core cannot see.  Its T2 fields are\n"
         "OVERWRITTEN from the bridge actually bound to the config, so only\n"
         "the run-plan half ever has to be supplied.")
      .def("generate_lumi", [](const Pipeline& p, double lumi_pb, bool poisson,
                               bool events, unsigned nthreads,
                               std::size_t chunk, const KnobRunContext& ctx) {
        PipelineConfig cfg = p.config();
        cfg.lumi_pb = lumi_pb;
        cfg.n_events = 0;
        cfg.poisson = poisson;
        std::unique_ptr<Pipeline> q;
        {
          py::gil_scoped_release unlock;
          q = std::make_unique<Pipeline>(cfg, p.plan());
        }
        return pipeline_generate(*q, 0, events, nthreads, chunk, ctx);
      }, py::arg("lumi_pb"), py::arg("poisson") = true,
         py::arg("events") = false, py::arg("nthreads") = 1,
         py::arg("chunk") = 4096, py::arg("context") = KnobRunContext(),
         "Rebuild this run at an integrated luminosity [pb^-1] and generate\n"
         "every event of it; same columnar dict as generate().")
      .def("write_hepmc", [](const Pipeline& p, const std::string& path,
                             std::uint64_t n) {
#if LIPOLGEN_HAVE_HEPMC3
        if (n == 0 || n > p.size()) n = p.size();
        py::gil_scoped_release unlock;
        HepMC3Writer w(path);
        std::vector<Event> buf;
        for (std::uint64_t i = 0; i < n; i += 2048) {
          const std::uint64_t j = std::min<std::uint64_t>(i + 2048, n);
          p.generate_range(i, j, buf);
          for (const Event& e : buf) w.write(e);
        }
        w.close();
        return n;
#else
        (void)p; (void)path; (void)n;
        throw std::runtime_error("LiPolGen was built without HepMC3");
#endif
      }, py::arg("path"), py::arg("n") = 0,
         "Write the first `n` events (0 = all) to a HepMC3 Asciiv3 file.");
}

// ---------------------------------------------------------------------- io

static void bind_io(py::module_& m) {
#if LIPOLGEN_HAVE_HEPMC3
  py::enum_<HepMC3Format>(m, "HepMC3Format")
      .value("Asciiv3", HepMC3Format::Asciiv3)
      .value("HepMC2Ascii", HepMC3Format::HepMC2Ascii);
  py::class_<HepMC3Writer>(m, "HepMC3Writer")
      .def(py::init<const std::string&, HepMC3Format, std::string,
                    std::string>(),
           py::arg("filename"), py::arg("format") = HepMC3Format::Asciiv3,
           py::arg("generator_name") = "LiPolGen",
           py::arg("generator_version") = LIPOLGEN_VERSION)
      .def("write", &HepMC3Writer::write, py::arg("event"))
      .def("close", &HepMC3Writer::close)
      .def("__enter__", [](HepMC3Writer& w) { return &w; },
           py::return_value_policy::reference)
      .def("__exit__", [](HepMC3Writer& w, py::object, py::object, py::object) {
        w.close();
        return false;
      });
#endif

#if LIPOLGEN_HAVE_PYTHIA8
  m.def("dot4", &dot4, py::arg("a"), py::arg("b"));
  m.def("dis_scattered_electron", [](double e_e, const Vec4& p_n, double x,
                                     double q2, double phi) -> py::object {
    Vec4 k;
    double y = 0;
    if (!dis_scattered_electron(e_e, p_n, x, q2, phi, &k, &y)) return py::none();
    return py::make_tuple(k, y);
  }, py::arg("e_e"), py::arg("p_n"), py::arg("x"), py::arg("q2"), py::arg("phi"));
  m.def("dis_parton_fraction", [](const Vec4& q, const Vec4& p_n,
                                  double xi_max, double m_out) -> py::object {
    double xi = 0;
    if (!dis_parton_fraction(q, p_n, &xi, xi_max, m_out)) return py::none();
    return py::cast(xi);
  }, py::arg("q"), py::arg("p_n"), py::arg("xi_max") = 1.0,
     py::arg("m_out") = 0.0);

  py::class_<HfsSummary>(m, "HfsSummary")
      .def_readonly("sigma_empz", &HfsSummary::sigma_empz)
      .def_readonly("px", &HfsSummary::px)
      .def_readonly("py", &HfsSummary::py)
      .def_readonly("pt", &HfsSummary::pt)
      .def_readonly("e", &HfsSummary::e)
      .def_readonly("pz", &HfsSummary::pz)
      .def_readonly("charge", &HfsSummary::charge)
      .def_readonly("n_total", &HfsSummary::n_total)
      .def_readonly("n_charged", &HfsSummary::n_charged)
      .def_readonly("n_neutral", &HfsSummary::n_neutral);
  m.def("hfs_summary", &hfs_summary, py::arg("ev"));
  m.def("hfs_sigma_empz_truth", &hfs_sigma_empz_truth, py::arg("e_e"),
        py::arg("y"), py::arg("p_n"));
  m.def("hfs_sigma_empz_exact", &hfs_sigma_empz_exact, py::arg("k"),
        py::arg("p_n"), py::arg("k_out"));

  py::enum_<CoherentT2>(m, "CoherentT2")
      .value("Pomeron", CoherentT2::Pomeron)
      .value("Off", CoherentT2::Off);

  py::enum_<NucleonChoice>(m, "NucleonChoice")
      .value("ByZN", NucleonChoice::ByZN)
      .value("Proton", NucleonChoice::Proton)
      .value("Neutron", NucleonChoice::Neutron);

  py::class_<PythiaBridgeOptions>(m, "PythiaBridgeOptions")
      .def(py::init<>())
      .def_readwrite("seed", &PythiaBridgeOptions::seed)
      .def_readwrite("settings", &PythiaBridgeOptions::settings)
      .def_readwrite("verbosity", &PythiaBridgeOptions::verbosity)
      .def_readwrite("max_retries", &PythiaBridgeOptions::max_retries)
      .def_readwrite("headroom", &PythiaBridgeOptions::headroom)
      .def_readwrite("include_strange", &PythiaBridgeOptions::include_strange)
      .def_readwrite("include_charm", &PythiaBridgeOptions::include_charm)
      .def_readwrite("include_bottom", &PythiaBridgeOptions::include_bottom)
      .def_readwrite("q2_pdf_min", &PythiaBridgeOptions::q2_pdf_min)
      .def_readwrite("with_neutron_instance",
                     &PythiaBridgeOptions::with_neutron_instance)
      .def_readwrite("coherent_t2", &PythiaBridgeOptions::coherent_t2,
                     "CoherentT2.Pomeron (default) builds the third, "
                     "Pomeron-beam (id 990) PYTHIA instance that hadronizes "
                     "the coherent channel; CoherentT2.Off skips it and "
                     "hadronize() then returns False on a coherent event.")
      .def_readwrite("pom_set", &PythiaBridgeOptions::pom_set,
                     "PYTHIA PDF:PomSet -- the Pomeron parton densities "
                     "(6 = H1 2006 Fit B LO, PYTHIA's own default).")
      .def_readwrite("pom_rescale", &PythiaBridgeOptions::pom_rescale,
                     "PYTHIA PDF:PomRescale; cancels out of the bridge's own "
                     "per-event flavour draw.")
      .def_readwrite("nucleon_choice", &PythiaBridgeOptions::nucleon_choice)
      .def_property("f2_source",
          [](const PythiaBridgeOptions& o) {
            return std::const_pointer_cast<UnpolSF>(o.f2_source);
          },
          [](PythiaBridgeOptions& o, std::shared_ptr<UnpolSF> u) {
            o.f2_source = std::move(u);
          },
          "Unpolarised structure functions of the T2 struck-nucleon SPECIES "
          "draw (nucleon_choice = ByStructureFunctions).  None = the "
          "library's ToyF2, which is also InclusiveKernel's default.  HAND "
          "IT THE SAME BACKEND THE KERNEL WAS BUILT WITH: the toy's "
          "F2n/F2p is a straight line and is up to 24 % away from CT18NLO's "
          "at x = 0.5, so a --unpol-sf ct18nlo run whose bridge is left on "
          "the toy draws its T2 species from a different nucleus than its "
          "rate.  python/lipolgen/cli.py does this wiring for --hadronize "
          "runs; a caller building a bridge directly must do it here.");

  py::class_<PythiaBridgeStats>(m, "PythiaBridgeStats")
      .def_readonly("n_called", &PythiaBridgeStats::n_called)
      .def_readonly("n_ok", &PythiaBridgeStats::n_ok)
      .def_readonly("n_failed", &PythiaBridgeStats::n_failed)
      .def_readonly("n_retries", &PythiaBridgeStats::n_retries)
      .def_readonly("n_no_surrogate", &PythiaBridgeStats::n_no_surrogate)
      .def_readonly("n_proton", &PythiaBridgeStats::n_proton)
      .def_readonly("n_neutron", &PythiaBridgeStats::n_neutron)
      .def_readonly("n_pomeron", &PythiaBridgeStats::n_pomeron,
                    "Coherent events hadronized off the Pomeron beam.")
      .def_readonly("n_pom_flavour_fallback",
                    &PythiaBridgeStats::n_pom_flavour_fallback,
                    "Coherent events whose Pomeron-PDF flavour weights were "
                    "all zero (the grid is pure gluon at small beta below "
                    "Q^2 ~ 1.5-1.75) and fell back to the bare e_q^2 charge "
                    "weights over the light flavours only; ~20% of a default "
                    "coherent run, so non-zero is routine -- and MEASURED "
                    "(D4.6), it costs exactly zero: raising q2_pdf_min to "
                    "1.75 drives the fallback to 0.00 % (sets 6 and 3) or a "
                    "few per cent (12, 13, 15) -- EXCEPT on pom_set = 4 "
                    "(H1 2006 Fit B NLO), where it only falls 35.52 % -> "
                    "29.32 % -- and leaves the final state bit-identical in "
                    "every one of those cases, set 4 INCLUDED, so that "
                    "residual costs nothing either, because every Pomeron "
                    "DPDF carries one light-quark singlet.  1.000 on "
                    "pom_set = 11, which is why that set is refused.  In the "
                    "npz meta since 2026-09-04, as n_pom_flavour_fallback "
                    "and pom_flavour_fallback_frac.")
      .def_readonly("n_flavour_dropped",
                    &PythiaBridgeStats::n_flavour_dropped)
      .def_readonly("n_cluster_fallback",
                    &PythiaBridgeStats::n_cluster_fallback,
                    "Events that took the DEPRECATED Role::StruckCluster "
                    "branch; non-zero means the run is at Tier.T0 and the "
                    "whole record does not conserve.")
      .def_readonly("max_rescale_dev", &PythiaBridgeStats::max_rescale_dev);

  py::class_<PythiaBridge, std::shared_ptr<PythiaBridge>>(m, "PythiaBridge")
      .def(py::init([](const BeamConfig& b, PythiaBridgeOptions o) {
        py::gil_scoped_release unlock;
        return std::make_shared<PythiaBridge>(b, std::move(o));
      }), py::arg("beams"), py::arg("options") = PythiaBridgeOptions())
      .def("hadronize", [](PythiaBridge& br, Event& ev, Rng& rng) {
        py::gil_scoped_release unlock;
        return br.hadronize(ev, rng);
      }, py::arg("event"), py::arg("rng"),
         "Shower and hadronize the event IN PLACE; False if PYTHIA vetoed.")
      .def_property_readonly("stats", &PythiaBridge::stats)
      .def_property_readonly("options", &PythiaBridge::options)
      .def("last_struck_nucleon", &PythiaBridge::last_struck_nucleon)
      .def("last_struck_nucleon_pdg", &PythiaBridge::last_struck_nucleon_pdg)
      .def("last_xi_pythia", &PythiaBridge::last_xi_pythia)
      .def("last_xi_physical", &PythiaBridge::last_xi_physical)
      .def("last_quark_id", &PythiaBridge::last_quark_id)
      .def("last_rescale", &PythiaBridge::last_rescale)
      .def("last_w_pythia", &PythiaBridge::last_w_pythia)
      .def("last_w_physical", &PythiaBridge::last_w_physical)
      .def("pythia_sigma_gen_mb", &PythiaBridge::pythia_sigma_gen_mb)
      .def("applied_settings", &PythiaBridge::applied_settings);

  // Bind the bridge straight into a PipelineConfig without a Python round
  // trip: the hook then runs with the GIL released.  The bridge is NOT
  // thread-safe, so a pipeline configured this way must generate with
  // nthreads = 1.
  m.def("set_pythia_hadronizer", [](PipelineConfig& cfg,
                                    std::shared_ptr<PythiaBridge> br) {
    // A NAMED callable, not a lambda: `columns_to_dict` recovers the bridge
    // from the type-erased hook with `std::function::target` and writes the
    // Pomeron / T2 block of the npz `meta` (D4.7).
    if (!br) throw std::invalid_argument(
        "set_pythia_hadronizer: bridge must not be None");
    // ONE UNPOLARISED BACKEND FOR THE WHOLE RUN, enforced where the two
    // objects finally meet.  `PipelineConfig::validate()` cannot do it -- the
    // core links no PYTHIA tier and has never heard of `PythiaBridge` -- and
    // the `Pipeline` cannot fill the slot the way it fills
    // `BreakupOptions::f2`, because the caller builds the bridge.  So this is
    // the `validate()` rule applied at the binding: a bridge whose T2
    // struck-nucleon SPECIES draw
    // (`NucleonChoice::ByStructureFunctions`, P(p) = Z F2p/(Z F2p + N F2n))
    // is not on the config's own `unpol_sf_obj` would draw from a backend
    // `meta["unpol_sf"]` does not name.  Measured price: the toy's F2n/F2p is
    // the straight line clip(1 - 0.75x, 0.25, 1), 0.9625 / 0.8500 / 0.6250 at
    // x = 0.05 / 0.20 / 0.50 against CT18NLO's 0.9218 / 0.7219 / 0.5035 --
    // up to 24 % away, on the species of every T2 event.
    //
    // IDENTITY, not equality: "the same object" is the rule everything else
    // in this run follows (the `Li6Convolution` b1 shares the kernel's own
    // `UnpolSF`; `BreakupOptions::f2` is the sampler kernel's own), and two
    // separately constructed `LhapdfSF("CT18NLO", 0)` are two grids whose
    // agreement nothing here checks.  Both null -- the shipped default -- is
    // the same object, so a default run is untouched.
    if (br->options().f2_source != cfg.unpol_sf_obj) {
      throw std::runtime_error(
          std::string("set_pythia_hadronizer: the bridge's species draw is "
                      "on a different unpolarised backend from this run's "
                      "kernel -- bridge f2_source is ") +
          (br->options().f2_source ? "an attached object"
                                   : "unset (the bridge's own ToyF2)") +
          " while the config's unpol_sf is " + unpol_sf_name(cfg.unpol_sf) +
          (cfg.unpol_sf_obj ? " with an object attached" : " with no object") +
          ".  The T2 struck-nucleon species draw is "
          "P(p) = Z F2p/(Z F2p + N F2n) on PythiaBridgeOptions::f2_source, "
          "so the run would draw its T2 species from one backend while "
          "meta[\"unpol_sf\"] recorded the other, and the toy's F2n/F2p is "
          "up to 24 % away from CT18NLO's.  Hand the bridge THE SAME OBJECT "
          "before constructing it: opts.f2_source = config.unpol_sf_obj "
          "(what lipolgen-run and lipolgen.run do), or leave both unset for "
          "the toy");
    }
    cfg.hadronizer = PythiaHadronizerHook{std::move(br)};
  }, py::arg("config"), py::arg("bridge"),
     "config.hadronizer = bridge.hadronize, as a pure C++ callable "
     "(no GIL per event; single-threaded generation only)");
#endif
}
