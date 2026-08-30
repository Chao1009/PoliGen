// pom_dis.cc -- COHERENT diffractive final state as gamma*-Pomeron DIS,
// driven through the SAME LHAup + surrogate + frame-map design LiPolGen's
// PythiaBridge already uses for gamma*-nucleon DIS.
//
// The identification that makes this work:
//
//     gamma* + N  ->  X            gamma* + IP  ->  X
//     W^2 = (q + p_N)^2            M_X^2 = (q + P_IP)^2
//     xi   = x_Bj                  zeta  = beta = Q^2/(M_X^2 + Q^2)   EXACT,
//                                          because P_IP^2 = 0 in PYTHIA
//
// so LiPolGen's existing bridge becomes the coherent bridge by swapping
//   Beams:idA  2212/2112 -> 990,
//   W^2        -> M_X^2  (the pipeline already computes it: CoherentEvent::m_x2),
//   the PDF     -> PYTHIA's own Pomeron PDF (PDF:PomSet), sampled at beta,
//   the target 4-vector -> P_IP = P_ion - P_recoil (already in the record as
//                          the Role::HadronicX pseudo-particle minus q).
//
// The beam remnant is then the POMERON remnant, so the string is
// struck-quark <-> Pomeron-remnant-antiquark: colour neutral, electrically
// neutral, and it carries the correct beta dependence (a low-beta event is
// gluon-initiated and gets a q qbar g topology from PYTHIA's own backward
// evolution, which is exactly what a hand-built q qbar string cannot do).
//
// Build:
//   source LiPolGen/env.sh
//   g++ -O2 -std=c++17 pom_dis.cc -o pom_dis $(pythia8-config --cxxflags --ldflags)
// Run:
//   ./pom_dis --nev 2000 --q2 5 --scan

#include "Pythia8/Pythia.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace Pythia8;

class LhaPomDis : public LHAup {
 public:
  LhaPomDis(double eA, double eB) : eA_(eA), eB_(eB) {
    setBeamA(990, eA);
    setBeamB(11, eB);
    addProcess(1, 1.0, 0.0, 1.0);
  }
  bool setInit() override { return true; }

  /// Arm the next event.  Returns false when the surrogate does not fit.
  bool arm(int idq, double q2, double mx, double xf) {
    idq_ = idq; q2_ = q2; mx_ = mx; xf_ = xf;
    zeta_ = q2 / (mx * mx + q2);                 // = beta, exactly
    const double B = 0.5 * (mx_ * mx_ + q2_);
    const double PAkB = 2.0 * eA_ * eB_;
    eprime_ = (PAkB - B + eA_ * q2_ / (2.0 * eB_)) / (2.0 * eA_);
    kz_ = q2_ / (2.0 * eB_) - eprime_;
    const double kt2 = eprime_ * eprime_ - kz_ * kz_;
    if (!(kt2 >= 0.0) || !(eprime_ > 0.0)) return false;
    kt_ = std::sqrt(kt2);
    return true;
  }
  bool setEvent(int = 0) override {
    const double pin = 0.5 * zeta_ * (2.0 * eA_);
    // Colour tag orientation MUST follow the sign of the initiator: an
    // incoming ANTIquark carries (0, 101), a quark (101, 0).  Getting this
    // wrong is a 50 % silent veto rate (measured).
    const int col  = (idq_ > 0) ? 101 : 0;
    const int acol = (idq_ > 0) ? 0 : 101;
    setProcess(1, 1.0, std::sqrt(q2_), 0.0072974, 0.118);
    addParticle(idq_, -1, 0, 0, col, acol, 0., 0., pin, pin, 0.);
    addParticle(11,   -1, 0, 0, 0,   0,    0., 0., -eB_, eB_, 0.);
    addParticle(11,    1, 1, 2, 0,   0,    kt_, 0., kz_, eprime_, 0.);
    addParticle(idq_,  1, 1, 2, col, acol, -kt_, 0.,
                -eB_ - kz_ + pin, eB_ - eprime_ + pin, 0.);
    setIdX(idq_, 11, zeta_, 1.0);
    setPdf(idq_, 11, zeta_, 1.0, std::sqrt(q2_), xf_ / zeta_, 1.0, true);
    return true;
  }
  double beta() const { return zeta_; }

 private:
  double eA_, eB_;
  int idq_ = 2;
  double zeta_ = 0.1, q2_ = 5.0, mx_ = 5.0, xf_ = 0.3;
  double eprime_ = 0.0, kz_ = 0.0, kt_ = 0.0;
};

int main(int argc, char** argv) {
  int nev = 2000, seed = 1;
  double q2 = 5.0, mxOne = 5.0;
  bool scan = false;
  for (int i = 1; i < argc; ++i) {
    std::string s = argv[i];
    if (s == "--nev") nev = std::atoi(argv[++i]);
    else if (s == "--q2") q2 = std::atof(argv[++i]);
    else if (s == "--mx") mxOne = std::atof(argv[++i]);
    else if (s == "--seed") seed = std::atoi(argv[++i]);
    else if (s == "--scan") scan = true;
  }

  Pythia pythia;
  auto lha = std::make_shared<LhaPomDis>(60.0, 20.0);
  pythia.setLHAupPtr(lha);
  pythia.readString("Beams:frameType = 5");
  pythia.readString("Check:beams = on");
  pythia.readString("SpaceShower:dipoleRecoil = on");
  pythia.readString("SpaceShower:pTmaxMatch = 2");
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");
  pythia.readString("SpaceShower:QEDshowerByQ = off");
  pythia.readString("LesHouches:setLeptonMass = 0");
  pythia.readString("Print:quiet = on");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  if (!pythia.init()) { printf("INIT FAILED\n"); return 1; }

  auto pdf = pythia.getPDFPtr(990);
  const int flav[6] = {1, 2, 3, -1, -2, -3};
  const double eq2[6] = {1./9, 4./9, 1./9, 1./9, 4./9, 1./9};

  std::vector<double> mxList =
      scan ? std::vector<double>{1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 40.0}
           : std::vector<double>{mxOne};

  printf("# Q^2 = %.2f GeV^2, nev = %d, PomSet = default (H1 2006 Fit B, "
         "LO)\n", q2, nev);
  printf("# %7s %8s %9s %9s %9s %9s %9s %9s\n", "M_X", "beta", "<n>", "<n_ch>",
         "<Q_tot>", "M_had/M_X", "veto", "P(gluon)");

  for (double mx : mxList) {
    const double beta = q2 / (mx * mx + q2);
    // Flavour draw P(q) ~ e_q^2 x f_q(beta, Q2) from PYTHIA's own Pomeron PDF,
    // and the gluon fraction of the Pomeron at this beta for reference.
    double w[6], wtot = 0.0;
    for (int i = 0; i < 6; ++i) {
      w[i] = eq2[i] * pdf->xf(flav[i], beta, q2);
      wtot += w[i];
    }
    const double xg = pdf->xf(21, beta, q2);
    double sumq = 0.0;
    for (int i = 0; i < 6; ++i) sumq += pdf->xf(flav[i], beta, q2);

    long nOk = 0, nVeto = 0;
    double sumN = 0, sumNch = 0, sumQ = 0, sumMratio = 0;
    for (int iev = 0; iev < nev; ++iev) {
      double r = pythia.rndm.flat() * wtot, acc = 0.0;
      int pick = 1;
      for (int i = 0; i < 6; ++i) { acc += w[i]; if (r < acc) { pick = i; break; } }
      if (!lha->arm(flav[pick], q2, mx, pdf->xf(flav[pick], beta, q2))) {
        ++nVeto; continue;
      }
      if (!pythia.next()) { ++nVeto; continue; }
      // Hadronic system = all final minus the scattered electron.
      Vec4 pH(0., 0., 0., 0.);
      double qtot = 0.0;
      int n = 0, nch = 0;
      bool foundE = false;
      for (int i = 0; i < pythia.event.size(); ++i) {
        if (!pythia.event[i].isFinal()) continue;
        if (!foundE && pythia.event[i].id() == 11 &&
            pythia.event[i].e() > 0.5 * lha->beta() * 0.0 + 1.0 &&
            pythia.event[i].status() > 60 && false) {}
        if (pythia.event[i].id() == 11 && !foundE &&
            pythia.event[i].mother1() > 0) { }
        pH += pythia.event[i].p();
        qtot += pythia.event[i].charge();
        ++n; if (pythia.event[i].isCharged()) ++nch;
      }
      // Subtract the scattered lepton: it is the bottom copy of entry 6 in
      // the LHA record (the outgoing e-).  Find the highest-|p| final e-.
      int ie = -1; double eBest = -1;
      for (int i = 0; i < pythia.event.size(); ++i)
        if (pythia.event[i].isFinal() && pythia.event[i].id() == 11 &&
            pythia.event[i].e() > eBest) { eBest = pythia.event[i].e(); ie = i; }
      if (ie >= 0) { pH -= pythia.event[ie].p();
                     qtot -= pythia.event[ie].charge();
                     --n; --nch; }
      sumN += n; sumNch += nch; sumQ += qtot;
      sumMratio += pH.mCalc() / mx;
      ++nOk;
    }
    if (nOk == 0) { printf("  %7.2f  all vetoed\n", mx); continue; }
    printf("  %7.2f %8.4f %9.3f %9.3f %9.4f %9.4f %9.4f %9.4f\n", mx, beta,
           sumN / nOk, sumNch / nOk, sumQ / nOk, sumMratio / nOk,
           double(nVeto) / nev, xg / (xg + sumq));
  }
  return 0;
}
