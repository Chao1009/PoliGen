// mx_string.cc -- prototype for LiPolGen T2 on the COHERENT channel.
//
// Hadronize the diffractive system X of mass M_X as a colour-singlet
// q-qbar (optionally q-qbar with FSR, or g-g) string, in the X rest frame,
// then boost the hadrons onto the physical P_X handed in by the pipeline.
//
// Build:
//   source LiPolGen/env.sh
//   g++ -O2 -std=c++17 mx_string.cc -o mx_string \
//       $(pythia8-config --cxxflags --ldflags)
//
// Usage:
//   ./mx_string --mx 3.0 --nev 2000 --seed 1 --mode qqbar
//   ./mx_string --scan --nev 2000
//
// What it demonstrates:
//   * PYTHIA 8.317 will hadronize a hand-filled colour-singlet string of any
//     mass above ~ 1 GeV with ProcessLevel:all = off + pythia.next().
//   * The final state can be boosted onto an arbitrary timelike P_X with a
//     PURE Lorentz transformation, so 4-momentum closes to the numerical
//     floor and the charge is exactly the string's (0 for q-qbar).
//   * <n_ch>(M_X) comes out on the ln(M_X^2) curve of e+e- / diffractive data.

#include "Pythia8/Pythia.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace Pythia8;

namespace {

struct Args {
  double mx = 3.0;
  double gamma = 0.0;      // boost of the X system along +z (0 = rest frame)
  double thetaX = 0.15;    // polar angle of P_X w.r.t. +z, to test rotations
  int nev = 2000;
  int seed = 1;
  std::string mode = "qqbar";   // qqbar | qqbarFSR | gg
  bool scan = false;
  bool labFill = false;    // fill partons already boosted (numerical test)
};

// Flavour of the q-qbar pair.  A DIS diffractive system at moderate beta is
// dominated by the gamma* -> q qbar dipole, so the flavour weights are the
// photon's e_q^2 (u : d : s : c = 4 : 1 : 1 : 4), NOT a Pomeron PDF; at low
// beta the gluon-initiated q qbar g configurations take over, which is what
// the "qqbarFSR" mode approximates.  Kept explicit so it can be replaced by
// a Pomeron-PDF draw.
int pickFlavour(Rndm& rndm, double mx) {
  const int id[4] = {2, 1, 3, 4};          // u, d, s, c
  double w[4] = {4.0, 1.0, 1.0, 0.0};      // e_q^2 weights
  if (mx > 4.0) w[3] = 4.0;                // open charm only above threshold
  double tot = w[0] + w[1] + w[2] + w[3];
  double r = rndm.flat() * tot, acc = 0.0;
  for (int i = 0; i < 4; ++i) { acc += w[i]; if (r < acc) return id[i]; }
  return 2;
}

}  // namespace

int main(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    std::string s = argv[i];
    auto nxt = [&]() { return std::string(argv[++i]); };
    if (s == "--mx") a.mx = std::stod(nxt());
    else if (s == "--gamma") a.gamma = std::stod(nxt());
    else if (s == "--theta") a.thetaX = std::stod(nxt());
    else if (s == "--nev") a.nev = std::stoi(nxt());
    else if (s == "--seed") a.seed = std::stoi(nxt());
    else if (s == "--mode") a.mode = nxt();
    else if (s == "--scan") a.scan = true;
    else if (s == "--labfill") a.labFill = true;
  }

  Pythia pythia;
  pythia.readString("ProcessLevel:all = off");
  pythia.readString("Check:event = off");
  pythia.readString("Print:quiet = on");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(a.seed));
  if (a.mode == "qqbarFSR") pythia.readString("PartonLevel:FSRinResonances = off");
  if (!pythia.init()) return 1;

  Event& event = pythia.event;
  ParticleData& pdt = pythia.particleData;

  std::vector<double> mxList;
  if (a.scan) mxList = {1.2, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 40.0};
  else mxList = {a.mx};

  printf("# mode = %s, nev = %d, seed = %d, gamma = %.3g, thetaX = %.3g, "
         "fill = %s\n",
         a.mode.c_str(), a.nev, a.seed, a.gamma, a.thetaX,
         a.labFill ? "lab" : "rest+boost");
  printf("# %8s %8s %8s %8s %10s %12s %12s %10s\n", "M_X", "<n>", "<n_ch>",
         "<n_neu>", "fail_frac", "max|dP|/E_X", "max|dQ|", "<pT_had>");

  for (double mx : mxList) {
    // The physical X four-vector the pipeline would hand us.
    double pXmag = (a.gamma > 0.0) ? mx * std::sqrt(a.gamma * a.gamma - 1.0) : 0.0;
    Vec4 pX(pXmag * std::sin(a.thetaX), 0.0, pXmag * std::cos(a.thetaX),
            std::sqrt(mx * mx + pXmag * pXmag));

    long nOk = 0, nFail = 0;
    double sumN = 0.0, sumNch = 0.0, sumPt = 0.0;
    double maxDP = 0.0, maxDQ = 0.0;

    for (int iev = 0; iev < a.nev; ++iev) {
      event.reset();

      if (a.mode == "gg") {
        double e = 0.5 * mx;
        // g g colour-singlet loop: the "gluonic Pomeron" configuration.
        event.append(21, 23, 101, 102, 0., 0.,  e, e);
        event.append(21, 23, 102, 101, 0., 0., -e, e);
      } else {
        int idq = pickFlavour(pythia.rndm, mx);
        double mq = pdt.m0(idq);
        double e = 0.5 * mx;
        double pp = sqrtpos(e * e - mq * mq);
        if (pp <= 0.0) { ++nFail; continue; }
        event.append( idq, 23, 101,   0, 0., 0.,  pp, e, mq);
        event.append(-idq, 23,   0, 101, 0., 0., -pp, e, mq);
      }

      if (a.mode == "qqbarFSR") {
        // Let the string radiate: this is what turns a 2-parton back-to-back
        // configuration into the q qbar g topology that dominates diffractive
        // DIS at low beta.  Max scale = M_X.
        event[1].scale(mx);
        event[2].scale(mx);
        pythia.forceTimeShower(1, 2, mx);
      }

      if (!pythia.next()) { ++nFail; continue; }

      // Boost the hadronized system onto the physical P_X.  This is the same
      // "pure Lorentz transformation" statement the DIS bridge relies on.
      Vec4 pSum(0., 0., 0., 0.);
      double qSum = 0.0;
      int nFin = 0, nCh = 0;
      double ptSum = 0.0;
      for (int i = 0; i < event.size(); ++i) {
        if (!event[i].isFinal()) continue;
        Vec4 p = event[i].p();
        if (pXmag > 0.0) p.bst(pX, mx);
        pSum += p;
        qSum += event[i].charge();
        ++nFin;
        if (event[i].isCharged()) ++nCh;
        ptSum += p.pT();
      }
      Vec4 d = pSum - pX;
      double dp = std::sqrt(d.e() * d.e() + d.px() * d.px() + d.py() * d.py()
                            + d.pz() * d.pz()) / pX.e();
      if (dp > maxDP) maxDP = dp;
      if (std::fabs(qSum) > maxDQ) maxDQ = std::fabs(qSum);
      sumN += nFin;
      sumNch += nCh;
      sumPt += (nFin > 0) ? ptSum / nFin : 0.0;
      ++nOk;
    }

    if (nOk == 0) { printf("  %8.3g  --- all failed\n", mx); continue; }
    printf("  %8.3g %8.3f %8.3f %8.3f %10.4f %12.3e %12.3g %10.4f\n", mx,
           sumN / nOk, sumNch / nOk, (sumN - sumNch) / nOk,
           double(nFail) / a.nev, maxDP, maxDQ, sumPt / nOk);
  }
  return 0;
}
