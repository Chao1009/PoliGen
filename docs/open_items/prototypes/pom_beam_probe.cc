// pom_beam_probe.cc -- can PYTHIA 8.317 be initialised with a POMERON beam
// (Beams:idA = 990) opposite an electron, in LHAup mode (frameType = 5)?
//
// If yes, the coherent channel needs NO new machinery: it is the existing
// LiPolGen DIS bridge with the target swapped from a nucleon to the Pomeron,
// W^2 -> M_X^2 and x_Bj -> beta.  The beam remnant is then the Pomeron
// remnant, so the string is q(struck) -- qbar(remnant), automatically colour
// neutral and electrically neutral.
//
// Build: g++ -O2 -std=c++17 pom_beam_probe.cc -o pom_beam_probe \
//        $(pythia8-config --cxxflags --ldflags)

#include "Pythia8/Pythia.h"
#include <cstdio>
#include <cmath>

using namespace Pythia8;

// Minimal LHAup: e + q(from Pomeron) -> e + q, exactly the shape of
// LiPolGen's src/pythia/lhaup_dis.hpp.
class LhaPomDis : public LHAup {
 public:
  LhaPomDis(double eA, double eB) : eA_(eA), eB_(eB) {
    setBeamA(990, eA);      // Pomeron along +z
    setBeamB(11, eB);       // electron along -z
    addProcess(1, 1.0, 0.0, 1.0);
  }
  bool setInit() override { return true; }
  // Set the next event: struck quark of flavour idq at LC fraction zeta,
  // photon virtuality Q2, gamma*-Pomeron invariant mass mx.
  void arm(int idq, double zeta, double q2, double mx, double xf) {
    idq_ = idq; zeta_ = zeta; q2_ = q2; mx_ = mx; xf_ = xf;
  }
  bool setEvent(int = 0) override {
    // Surrogate frame: Pomeron massless along +z at eA_, electron at -z.
    const double pAplus = 2.0 * eA_;
    // Solve for the scattered electron so that -q^2 = Q2 and (q + P_A)^2 = mx^2
    // with P_A^2 = 0 (Pomeron treated massless, as PYTHIA's ParticleData has it)
    const double B = 0.5 * (mx_ * mx_ + q2_);       // P_A . q
    const double PAkB = eA_ * eB_ * 2.0;            // P_A . k_B (massless, head-on)
    const double eprime = (PAkB - B + eA_ * q2_ / (2.0 * eB_)) / (2.0 * eA_);
    const double kz = q2_ / (2.0 * eB_) - eprime;
    const double kt2 = eprime * eprime - kz * kz;
    if (kt2 < 0.0) return false;
    const double kt = std::sqrt(kt2);
    const double pin = 0.5 * zeta_ * pAplus;

    setProcess(1, 1.0, std::sqrt(q2_), 0.0072974, 0.118);
    // 1: incoming quark out of beam A (must come first)
    addParticle(idq_, -1, 0, 0, 101, 0, 0., 0., pin, pin, 0.);
    // 2: incoming electron
    addParticle(11, -1, 0, 0, 0, 0, 0., 0., -eB_, eB_, 0.);
    // 3: outgoing electron
    addParticle(11, 1, 1, 2, 0, 0, kt, 0., kz, eprime, 0.);
    // 4: outgoing quark
    const double qx = -kt, qy = 0.0;
    const double qz = -eB_ - kz + pin, qe = eB_ - eprime + pin;
    addParticle(idq_, 1, 1, 2, 101, 0, qx, qy, qz, qe, 0.);
    setIdX(idq_, 11, zeta_, 1.0);
    setPdf(idq_, 11, zeta_, 1.0, std::sqrt(q2_), xf_ / zeta_, 1.0, true);
    return true;
  }
 private:
  double eA_, eB_;
  int idq_ = 2;
  double zeta_ = 0.1, q2_ = 5.0, mx_ = 5.0, xf_ = 0.3;
};

int main(int argc, char** argv) {
  int pomSet = (argc > 1) ? std::atoi(argv[1]) : -1;

  Pythia pythia;
  auto lha = std::make_shared<LhaPomDis>(50.0, 20.0);
  pythia.setLHAupPtr(lha);
  pythia.readString("Beams:frameType = 5");
  pythia.readString("Check:beams = on");
  pythia.readString("SpaceShower:dipoleRecoil = on");
  pythia.readString("SpaceShower:pTmaxMatch = 2");
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");
  pythia.readString("SpaceShower:QEDshowerByQ = off");
  pythia.readString("LesHouches:setLeptonMass = 0");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  if (pomSet >= 0)
    pythia.readString("PDF:PomSet = " + std::to_string(pomSet));

  printf("=== trying Beams:idA = 990 (Pomeron) vs idB = 11, frameType 5, "
         "PDF:PomSet = %d ===\n", pomSet);
  if (!pythia.init()) { printf("INIT FAILED\n"); return 1; }
  printf("INIT OK.  eCM = %.3f,  idA = %d, idB = %d\n",
         pythia.info.eCM(), pythia.info.idA(), pythia.info.idB());

  // What does the Pomeron PDF look like?
  auto pdf = pythia.getPDFPtr(990);
  if (pdf) {
    printf("\nPomeron PDF at Q2 = 10:  beta   xg      xu      xd      xs\n");
    for (double b : {0.05, 0.2, 0.5, 0.8, 0.95}) {
      printf("   %20.2f  %7.4f %7.4f %7.4f %7.4f\n", b,
             pdf->xf(21, b, 10.0), pdf->xf(2, b, 10.0),
             pdf->xf(1, b, 10.0), pdf->xf(3, b, 10.0));
    }
  } else printf("\nNO Pomeron PDF pointer\n");
  return 0;
}
