# 00 — Reference corpus inventory

**Purpose.** Before searching for new references for LiPolGen, this document
records what the project already has, so later searches target the gaps
and do not re-find what is already on disk. Two sources were read in full:

- `/home/cpeng/Projects/polli/PolarizedLithiumSim/refs/README.md` and
  `refs_dict.json` — the sibling repository's reference corpus: **54
  entries**, PDFs named by arXiv id or a short tag, committed under
  `refs/`. This is a *shared* corpus — most of its entries serve the
  sibling repo's own reconstruction-chain / detector-R&D work and are
  **not** cited anywhere in the LiPolGen tree; each row below says which.
- LiPolGen's own docs: `docs/PHYSICS_CHANNELS.md` (its `## 14. References`
  section, 58 `[KEY]`-labelled entries), `docs/open_items/physics_literature.md`,
  `docs/benchmarking/00_in_tree_checks.md` … `06_critic.md` (the
  201-resource survey, whose `T-1`…`T-46` table in `04_theory.md` is the
  single richest source of *candidate* references), `docs/CONVENTIONS.md`,
  and `docs/benchmarking/BENCHMARK_PLAN.md` (the `T1`…`T6` tier scheme
  used below).

**Method.** `grep`/`python` over the corpus JSON and every named doc;
arXiv ids cross-checked against the corpus's own id set (not re-fetched
over the network in this pass — this is an inventory, not a verification,
step). Every gap below is a reference *named* in the LiPolGen docs (by
author, journal, or arXiv id) that does **not** appear in `refs_dict.json`
and has no PDF under `refs/`. Nothing here has been downloaded; that is
the next stage's job.

---

## 1. The corpus already held (54 entries, `refs/refs_dict.json`)

`PDF on disk` = "yes" means a file is committed under
`PolarizedLithiumSim/refs/`; "no (no free copy)" means the dictionary
entry exists (title/authors/journal recorded) but no accessible copy
was ever found — these are the corpus's own pre-existing gaps, not new
ones found in this pass. `Tier` maps each entry onto
`BENCHMARK_PLAN.md`'s T1 (nucleon) / T2 (deuteron & b1) / T3 (nuclear
inputs) / T4 (generators) / T5 (MC closures, n/a here) / T6 ("chain":
detector, EIC machine, kinematic reconstruction) scheme, or
`theory-input` for structure-function formalism that doesn't fit T1–T6.

| # | Key | Citation | arXiv | PDF on disk | Tier | Where cited in LiPolGen docs |
|---|---|---|---|---|---|---|
| 1 | `abe1999-r1998-world-fit` | K. Abe et al. (E143 Collaboration), *Measurements of R = sigma_L/sigma_T for 0.03 < x < 0.1 and Fit to World Data*, Phys. Lett. B 452 (1999) 194 (SLAC-PUB-7927) | hep-ex/9808028 | yes | T1 | PHYSICS_CHANNELS.md [R1998]; benchmarking/00,02 |
| 2 | `aggarwal-caldwell2022-kinematic-fitting` | R. Aggarwal, A. Caldwell, *Kinematic fitting of neutral current events in deep inelastic ep collisions*, JINST 17 (2022) P09035 | 2206.04897 | yes | chain | not found by direct grep in LiPolGen/docs — sibling-repo kinematic-reconstruction work only |
| 3 | `arratia2022-dnn-dis-kinematics` | M. Arratia, D. Britzger, O. Long, B. Nachman, *Reconstructing the Kinematics of Deep Inelastic Scattering with Deep Learning*, Nucl. Instrum. Meth. A 1025 (2022) 166164 | 2110.05505 | yes | chain | not found by direct grep in LiPolGen/docs — sibling-repo kinematic-reconstruction work only |
| 4 | `athena2022-detector-proposal` | ATHENA Collaboration, J. Adam et al., *ATHENA detector proposal -- a totally hermetic electron nucleus apparatus proposed for IP6 at the Electron-Ion Collider*, JINST 17 (2022) P10019 | 2210.09048 | yes | chain | PHYSICS_CHANNELS.md |
| 5 | `bacchetta2007-sidis-small-pt` | A. Bacchetta, M. Diehl, K. Goeke, A. Metz, P. J. Mulders, M. Schlegel, *Semi-inclusive deep inelastic scattering at small transverse momentum*, JHEP 02 (2007) 093 | hep-ph/0611265 | yes | theory-input | PHYSICS_CHANNELS.md [Bacchetta07] (flagged there as 'not registered in the repository's literature list' even though cited in code); open_items/physics_literature.md |
| 6 | `ball2003-compass-6lid-target` | J. Ball et al., *First results of the large COMPASS 6LiD polarized target*, Nucl. Instrum. Meth. A 498 (2003) 101 | — | no (no free copy) | T3 | not found by direct grep in LiPolGen/docs — sibling-repo target-polarization sourcing only |
| 7 | `bultmann1999-lid-target-slac` | S. Bultmann et al., *The SLAC high density polarized 6LiD target*, Nucl. Instrum. Meth. A 425 (1999) 23 (SLAC-PUB-7904) | — | no (no free copy) | T3 | not found by direct grep in LiPolGen/docs — sibling-repo target-polarization sourcing only |
| 8 | `chang2026-ir8-light-nuclei` | W. Chang, E.-C. Aschenauer, A. Jentsch, A. Kumar, Z. Tu, Z. Yin, *Opportunities for imaging light nuclei with a second interaction region at the Electron-Ion Collider*, Phys. Rev. D 113 (2026) 032018 | 2511.05638 | yes | chain | PHYSICS_CHANNELS.md [Chang26]; benchmarking/00 (row E-16); open_items/physics_literature.md |
| 9 | `chaumette1989-6lid-7lih-saclay` | P. Chaumette, J. Deregel, G. Durand, J. Fabre, L. van Rossum, *Progress report on polarization of irradiated 6LiD and 7LiH*, AIP Conf. Proc. 187 (1989) 1275, doi:10.1063/1.38301 (body paywalled; the publisher abstract was read through the IAEA INIS record and OSTI) | — | no (no free copy) | T3 | not found by direct grep in LiPolGen/docs — sibling-repo target-polarization sourcing only |
| 10 | `cloet2006-polarized-emc` | I. C. Cloet, W. Bentz, A. W. Thomas, *Spin-dependent structure functions in nuclear matter and the polarized EMC effect*, Phys. Lett. B 642 (2006) 210 | nucl-th/0605061 | yes | theory-input | PHYSICS_CHANNELS.md [CBT06]; open_items/physics_literature.md |
| 11 | `compass2007-spectrometer-and-6lid-target` | P. Abbon et al. (COMPASS Collaboration), *The COMPASS experiment at CERN*, Nucl. Instrum. Meth. A 577 (2007) 455 | hep-ex/0703049 | yes | T3 | PHYSICS_CHANNELS.md (COMPASS 6LiD target dilution figures) |
| 12 | `cosyn-weiss-tagged-tensor2026` | W. Cosyn, C. Weiss, *Semi-inclusive deep-inelastic scattering on a polarized spin-1 target. I. Cross section and spin observables; II. Deuteron and spectator nucleon tagging*, Phys. Rev. C 114 (2026) 025201; 025202 | 2603.23699; 2603.23700 | yes | T2 | PHYSICS_CHANNELS.md [CWappD]/[CW26]; benchmarking/00,02,04,06; open_items/physics_literature.md |
| 13 | `cosyn2017-b1-convolution` | W. Cosyn, Yu-Bing Dong, S. Kumano, M. Sargsian, *Deuteron tensor structure function b1*, Phys. Rev. D 95 (2017) 074036 | 1702.05337 | yes | theory-input | PHYSICS_CHANNELS.md [CDKS17]; benchmarking/04,06; open_items/physics_literature.md |
| 14 | `cosyn2025-tensor-polarization-options` | W. Cosyn, B. Roldan Tomei, A. Sosa, A. Zec, *Polarization options in inclusive DIS off tensor polarized deuteron*, Eur. Phys. J. A 61 (2025) 83 | 2410.12764 | yes | T2 | PHYSICS_CHANNELS.md [Cosyn25]; benchmarking/00,02,04; open_items/physics_literature.md |
| 15 | `diefenthaler2022-deeply-learning-dis` | M. Diefenthaler, A. Farhat, A. Verbytskyi, Y. Xu, *Deeply Learning Deep Inelastic Scattering Kinematics*, Eur. Phys. J. C 82 (2022) 1064 | 2108.11638 | yes | chain | not found by direct grep in LiPolGen/docs — sibling-repo kinematic-reconstruction work only |
| 16 | `e155-1999-deuteron-g1-and-lid-target` | P. L. Anthony et al. (E155 Collaboration), *Measurement of the deuteron spin structure function g1d(x) for 1 (GeV/c)^2 < Q^2 < 40 (GeV/c)^2*, Phys. Lett. B 463 (1999) 339 | hep-ex/9904002 | yes | T3 | benchmarking/02_data_nucleon_deuteron.md (6LiD target dilution/polarization paragraph) |
| 17 | `eic-conceptual-design-report2021` | F. Willeke et al., J. Beebe-Wang (ed.), *Electron-Ion Collider Conceptual Design Report*, BNL-221006-2021-FORE (2021), doi:10.2172/1765663; 972 pp., not committed to refs/ (178 MB) | — | no (no free copy) | chain | not found by direct grep in LiPolGen/docs — sibling-repo beam-parameter sourcing only (ion cap, xi_p limit) |
| 18 | `eic-yellow-report2021` | R. Abdul Khalek et al., *Science Requirements and Detector Concepts for the Electron-Ion Collider: EIC Yellow Report*, Nucl. Phys. A 1026 (2022) 122447 | 2103.05419 | yes | chain | PHYSICS_CHANNELS.md [YR]; benchmarking/01_generators.md |
| 19 | `epic-inclusive-wg-resources` | P. Newman, S. Maple, B. Schmookler, C. Gwenlan (WG contacts), *ePIC Inclusive Physics WG wiki; BNL EIC 'DIS Kinematics' wiki; JeffersonLab/dis-reconstruction*, web resources (checked 2026-08-25) | — | no (no free copy) | chain | not found by direct grep in LiPolGen/docs — sibling-repo kinematic-reconstruction work only |
| 20 | `epic-preliminary-design-report2024` | ePIC Collaboration, *The ePIC Detector Preliminary Design Report (September 2024 draft, chapters 2 and 8)*, Zenodo record 13866213, doi:10.5281/zenodo.13866213, published 2024-10-01 (114 MB, 282 pp.), 'EIC preTDR - Chapters 2 and 8, DRAFT, Version0'; not committed to refs/ | — | no (no free copy) | chain | PHYSICS_CHANNELS.md |
| 21 | `epic-preliminary-tdr-v31-2026` | ePIC Collaboration, *The ePIC Detector Preliminary Technical Design Report, Version 3.1*, Zenodo, concept doi:10.5281/zenodo.18271601; record 18271602 (2026-01-16, 260 MB, 553 pp.) and record 19496158 (2026-04-10, the line-numbered build); CC-BY-4.0; not committed to refs/ | — | no (no free copy) | chain | PHYSICS_CHANNELS.md |
| 22 | `epios2026-polarized-ion-beams` | G. Atoian et al., *Realizing the scientific program with polarized ion beams at EIC (EPIOS white paper)*, Phys. Rev. C 113 (2026) 060501 | 2510.10794 | yes | chain | PHYSICS_CHANNELS.md [EPIOS]; benchmarking/00,06 |
| 23 | `gamage2021-second-interaction-region` | B. R. Gamage, V. Burkert, R. Ent, Y. Furletova, D. Higinbotham, A. Hutton, F. Lin, T. Michalski, V. S. Morozov, R. Rajput-Ghoshal, D. Romanov, T. Satogata, A. Seryi, A. Sy, C. Weiss, M. Wiseman, W. Wittmer, Y. Zhang, E.-C. Aschenauer, J. S. Berg, A. Jentsch, A. Kiselev, C. Montag, R. Palmer, B. Parker, V. Ptitsyn, F. Willeke, H. Witte, C. Hyde, P. Nadel-Turonski, *Design Concept for the Second Interaction Region for Electron-Ion Collider*, Proc. IPAC'21, TUPAB040, pp. 1435-1438, doi:10.18429/JACoW-IPAC2021-TUPAB040 | 2105.13564 | yes | chain | not found by direct grep in LiPolGen/docs (IR-2 dispersive-tagging threshold discussion lives in the sibling repo's Report 4) |
| 24 | `hamwi2026-li-resonances` | E. Hamwi, G. H. Hoffstaetter, *Polarization Transmission in the Electron-Ion Collider's Hadron Storage Ring*, Phys. Rev. Accel. Beams 29 (2026) 073501 | 2509.18558 | yes | chain | not found by direct grep in LiPolGen/docs — sibling-repo beam-polarization-transmission work only |
| 25 | `hoodbhoy-jaffe-manohar1989-and-jaffe-manohar1989` | P. Hoodbhoy, R. L. Jaffe, A. Manohar, R. L. Jaffe, A. Manohar, *Novel effects in deep inelastic scattering from spin-one hadrons; Nuclear gluonometry*, Nucl. Phys. B 312 (1989) 571; Phys. Lett. B 223 (1989) 218 | — | no (no free copy) | theory-input | CONVENTIONS.md; PHYSICS_CHANNELS.md [HJM89]; benchmarking/00,04,06,BENCHMARK_PLAN.md; open_items/physics_literature.md — NOTE: this corpus key pairs HJM89 (NPB312) with J&M's *Nuclear gluonometry* PLB223, NOT J&M's *Deep inelastic scattering from arbitrary spin targets* NPB321 that PHYSICS_CHANNELS.md's own [JM89] entry cites (see gap list) |
| 26 | `jacquet-blondel1979-and-bassler-bernardi1995` | F. Jacquet, A. Blondel, U. Bassler, G. Bernardi, *Jacquet-Blondel hadronic method (DESY 79/48, no free copy); U. Bassler, G. Bernardi, On the kinematic reconstruction of deep inelastic scattering at HERA: the Sigma method (DESY 94-231)*, DESY 79/48 (1979) 391; Nucl. Instrum. Meth. A 361 (1995) 197 | hep-ex/9412004 | yes | chain | not found by direct grep in LiPolGen/docs — sibling-repo kinematic-reconstruction work only |
| 27 | `jentsch-tu-weiss2021-ed-tagging` | A. Jentsch, Z. Tu, C. Weiss, *Deep-inelastic electron-deuteron scattering with spectator nucleon tagging at the future Electron-Ion Collider: Extracting free nucleon structure*, Phys. Rev. C 104 (2021) 065205 | 2108.08314 | yes | T2 | PHYSICS_CHANNELS.md (tagged e+d design comparator) |
| 28 | `jentsch2023-epic-far-forward-dis` | A. Jentsch (on behalf of the ePIC Collaboration), *Far-Forward Detectors and Physics with ePIC @ the EIC*, Talk, DIS 2023 (30th International Workshop on Deep-Inelastic Scattering), Michigan State University, 27-31 March 2023 | — | yes | chain | PHYSICS_CHANNELS.md (far-forward routing table) |
| 29 | `jlab-pr12-13-011-b1-and-2023-jeopardy` | K. Slifer et al., *The deuteron tensor structure function b1*, JLab proposal PR12-13-011 / E12-13-011 to PAC 40 (2013), and its jeopardy update to PAC 51 (2023) | — | no (no free copy) | T2 | not found by direct grep in LiPolGen/docs by that key name — but this IS the 'JLab b1 proposal' the task names; CONVENTIONS.md separately cites arXiv:2506.04506 p.8 for the same E12-13-011 kinematic edge, which may be a fetchable proxy — see gaps table |
| 30 | `jlab-pr12-14-001-polarized-emc` | W. K. Brooks, S. E. Kuhn, et al., *The EMC effect in spin structure functions*, JLab proposal PR12-14-001 / E12-14-001 to PAC 42 (22 May 2014) | — | no (no free copy) | T3 | not found by direct grep in LiPolGen/docs by that key name — but this IS the 'JLab b1/EMC proposal' the task names; see gaps table (proposal text itself has no free copy) |
| 31 | `koivuniemi2004-compass-polarization-buildup` | J. Koivuniemi et al., *Polarization build up in COMPASS 6LiD target*, Proc. SPIN 2004 (Trieste), p. 796 | — | no (no free copy) | T3 | not found by direct grep in LiPolGen/docs — sibling-repo target-polarization sourcing only |
| 32 | `lee2024-snspd-120gev-proton-beam-test` | S. Lee, T. Polakovic, W. Armstrong, A. Dibos, T. Draher, N. Pastika, Z.-E. Meziani, V. Novosad, *Beam Tests of SNSPDs with 120 GeV Protons*, Nucl. Instrum. Meth. A 1069 (2024) 169956 | 2312.13405 | yes | chain | not found by direct grep — sibling-repo near-beam nanowire-detector study only |
| 33 | `lee2026-cryogenic-test-station-alpha` | S. Lee, W. Armstrong, J. DiPreta, C. Dulya, V. Novosad, T. Polakovic, *Optimization of Cryogenic Detector Test Station by Rejecting Electromagnetic Interference*, Nucl. Instrum. Meth. A 1093 (2027) 171953 | 2601.03158 | yes | chain | not found by direct grep — sibling-repo near-beam nanowire-detector study only |
| 34 | `li-sick-whitney-yearian1971-6li-form-factor` | G. C. Li, I. Sick, R. R. Whitney, M. R. Yearian, *High-energy electron scattering from 6Li*, Nucl. Phys. A 162 (1971) 583 (also Suelzle et al., Phys. Rev. 162 (1967) 992) | — | no (no free copy) | T3 | not found by direct grep — but the SAME paper is the gap the task names as 'Suelzle/Li-Sick elastic tables' (see gaps table; benchmarking/03_data_nuclear.md discusses it at length without matching this key by author-name grep) |
| 35 | `li2023-eic-tracking-reference-design` | Xuan Li, *Latest vertex and tracking detector developments for the future Electron-Ion Collider*, arXiv:2305.15593 [physics.ins-det] (24 May 2023); JPS Conf. Proc. 42 (2024) 011032 | — | no (no free copy) | chain | not found by direct grep in LiPolGen/docs — sibling-repo tracking-resolution sourcing only |
| 36 | `mantysaari2024-polarized-deuteron-imaging` | H. Mantysaari, F. Salazar, B. Schenke, C. Shen, W. Zhao, *Spatial imaging of polarized deuterons at the Electron-Ion Collider*, Phys. Lett. B 858 (2024) 139053 | 2408.13213 | yes | T2 | PHYSICS_CHANNELS.md [Mant24]; benchmarking/04,06; open_items/physics_literature.md |
| 37 | `maple2024-epic-tracking-inclusive-dis-seminar` | S. Maple (ePIC / University of Birmingham), *Tracking and Inclusive DIS reconstruction with the ePIC detector at the EIC*, Seminar, University of Birmingham, 11 December 2024 (Indico Global event 4787; slides EIC_Seminar_SMaple.pdf, 67 pages) | — | yes | chain | not found by direct grep in LiPolGen/docs — sibling-repo kinematic-reconstruction work only |
| 38 | `miller2014-b1-pion-hidden-color` | G. A. Miller, *Pionic and hidden-color, six-quark contributions to the deuteron b1 structure function*, Phys. Rev. C 89 (2014) 045203 | 1311.4561 | yes | theory-input | PHYSICS_CHANNELS.md [Miller14]; benchmarking/04_theory.md |
| 39 | `nikolaev1998-azimuthal-asymmetry-ddis` | N. N. Nikolaev, A. V. Pronyaev, B. G. Zakharov, *Azimuthal asymmetry as a new handle on sigma_L/sigma_T in diffractive DIS*, (preprint, 1 Dec 1998; IKP Juelich / Landau Institute / Virginia Tech) | hep-ph/9812212 | yes | chain | benchmarking/06_critic.md |
| 40 | `pecar-vossen2022-sidis-ml-reconstruction` | C. Pecar, A. Vossen, *Reconstruction of event kinematics in semi-inclusive deep-inelastic scattering using the hadronic final state and machine learning*, Proceedings DIS2022 (Santiago de Compostela, May 2022) | 2209.14489 | yes | chain | not found by direct grep in LiPolGen/docs — sibling-repo kinematic-reconstruction work only |
| 41 | `pena2025-smspd-2mm-array` | C. Pena, C. Wang, S. Xie, et al. (Caltech / JPL / Fermilab), *Characterization of a Superconducting Microwire Single-Photon Detector Array for Charged-Particle Detection*, JINST 20 (2025) P03001 | 2410.00251 | yes | chain | not found by direct grep — sibling-repo near-beam nanowire-detector study only |
| 42 | `polakovic2019-snspd-high-magnetic-field` | T. Polakovic, W.R. Armstrong, V. Yefremenko, J.E. Pearson, K. Hafidi, G. Karapetrov, Z.-E. Meziani, V. Novosad, *Superconducting nanowires as high-rate photon detectors in strong magnetic fields*, Nucl. Instrum. Meth. A 959 (2020) 163543 | 1907.13059 | yes | chain | not found by direct grep — sibling-repo near-beam nanowire-detector study only |
| 43 | `pronyaev1998-forward-cone-lt` | A. V. Pronyaev, *The forward cone and L/T separation in diffractive DIS*, (conference proceedings, 27 Aug 1998) | hep-ph/9808432 | yes | chain | not found by direct grep in LiPolGen/docs — sibling-repo diffractive-DIS L/T work only |
| 44 | `ruspa2002-inclusive-diffraction-hera` | M. Ruspa, *Inclusive diffraction at HERA*, (proceedings, 14 Jun 2002; on behalf of H1 and ZEUS) | hep-ex/0206031 | yes | chain | not found by direct grep in LiPolGen/docs — sibling-repo diffractive-DIS work only |
| 45 | `sather-schmidt1990-double-helicity-flip-moment` | E. Sather, C. Schmidt, *Size and scaling of the double-helicity-flip hadronic structure function*, Phys. Rev. D 42 (1990) 1424 | — | no (no free copy) | theory-input | PHYSICS_CHANNELS.md [SS90]; benchmarking/04,06 (source of C_BAG) |
| 46 | `stone2019-magnetic-dipole-moments` | N. J. Stone, *Table of recommended nuclear magnetic dipole moments: Part I, long-lived states*, IAEA INDC(NDS)-0794 (November 2019) | — | no (no free copy) | T3 | PHYSICS_CHANNELS.md; benchmarking/03,06,BENCHMARK_PLAN.md |
| 47 | `stone2021-quadrupole-moments` | N. J. Stone, *Table of nuclear electric quadrupole moments*, IAEA INDC(NDS)-0833 (October 2021), doi:10.61092/iaea.a6te-dg7q; also At. Data Nucl. Data Tables 111-112 (2016) 1 | — | no (no free copy) | T3 | PHYSICS_CHANNELS.md; benchmarking/03,06,BENCHMARK_PLAN.md |
| 48 | `trento2004-and-diehl-sapeta2005` | A. Bacchetta, U. D'Alesio, M. Diehl, C. A. Miller, M. Diehl, S. Sapeta, *Single-spin asymmetries: the Trento conventions; On the analysis of lepton scattering on longitudinally or transversely polarized protons*, Phys. Rev. D 70 (2004) 117504; Eur. Phys. J. C 41 (2005) 515 | hep-ph/0410050; hep-ph/0503023 | yes | theory-input | CONVENTIONS.md; PHYSICS_CHANNELS.md; benchmarking/00,02,04,06,BENCHMARK_PLAN.md |
| 49 | `tronchin2018-polarized-emc` | S. Tronchin, H. H. Matevosyan, A. W. Thomas, *Polarized EMC effect in the QMC model*, Phys. Lett. B 783 (2018) 247 | 1806.00481 | yes | theory-input | PHYSICS_CHANNELS.md [TMT18]; benchmarking/00_in_tree_checks.md |
| 50 | `tunl-a6-evaluation2002` | D. R. Tilley et al. (TUNL), *Energy levels of light nuclei A = 5, 6, 7*, Nucl. Phys. A 708 (2002) 3 | — | yes | T3 | CONVENTIONS.md; PHYSICS_CHANNELS.md [TUNL6]; benchmarking/00,03,06,BENCHMARK_PLAN.md; open_items/physics_literature.md |
| 51 | `wang2021-gluon-spin-emc` | X. G. Wang, W. Bentz, I. C. Cloet, A. W. Thomas, *Polarized gluon EMC effect*, J. Phys. G 49 (2022) 03LT01 | 2109.03591 | yes | theory-input | PHYSICS_CHANNELS.md (gluon-spin EMC arm of plans/02 step 1.2.2) — not matched by author-surname grep because entry is 'Wang' (too common a token to grep reliably) |
| 52 | `wang2026-smspd-microwire-arrays` | C. Wang, C. Pena, A. Bornheim, S. Wu, A. Albert, T. Sievert, A. Apresyan, E. Knehr, B. Korzh, J. Luskin, L. Mori, S. Patel, G. Reales Gutierrez, M. Sahu, E. Schmidt, M. Shaw, E. Sledge, M. Spiropulu, T. Taher, S. Xie, *Towards High-Efficiency Particle Detection Using Superconducting Microwire Arrays*, submitted to JINST | 2510.11725 | yes | chain | not found by direct grep — sibling-repo near-beam nanowire-detector study only |
| 53 | `wiringa2014-momentum-distributions` | R. B. Wiringa, R. Schiavilla, S. C. Pieper, J. Carlson, *Nucleon and nucleon-pair momentum distributions in A <= 12 nuclei*, Phys. Rev. C 89 (2014) 024305 | 1309.3794 | yes | T3 | PHYSICS_CHANNELS.md [Wiringa14]; benchmarking/04_theory.md; open_items/physics_literature.md |
| 54 | `zeus2009-leading-protons-lrg` | ZEUS Collaboration, S. Chekanov et al., *Deep inelastic scattering with leading protons or large rapidity gaps at HERA*, Nucl. Phys. B 816 (2009) 1 (DESY-08-175) | 0812.2003 | yes | chain | benchmarking/00_in_tree_checks.md (diffractive azimuthal asymmetry gate) |

## 2. Gaps — the items named in the task prompt

These twelve are the references (or reference-families) the task
explicitly flagged as likely gaps. All twelve are confirmed absent from
`refs_dict.json` (checked by string search on author surname, arXiv id,
and DOI fragment) — except where noted, they are also absent as PDFs.

| # | Reference | Cited in LiPolGen as | Status | Priority |
|---|---|---|---|---|
| G1 | L. W. Mo, Y.-S. Tsai, *Radiative corrections to elastic and inelastic ep and mu-p scattering*, Rev. Mod. Phys. **41** (1969) 205 | `[MT69]`, PHYSICS_CHANNELS.md #50; `tests/test_rc.cpp` gate deliberately left unwritten; benchmarking/00 §9 item 5; 04_theory.md T-35 | **Paywalled at APS, no arXiv** (1088 INSPIRE citations). **A free substitute is already identified and verified in-tree**: Y.-S. Tsai, *Radiative corrections to electron scatterings*, SLAC-PUB-848 (1971), INSPIRE recid 67278 — `04_theory.md` T-36 records "I verified the full text is freely downloadable" at `https://inspirehep.net/files/3ce239706be17b0eefec145433700c64` (HTTP 200, `application/pdf`, 4.1 MB). This is what most of the field cites for the exact formulas anyway. | **HIGH** — closes a documented open item (`tests/test_rc.cpp`'s "deliberately not written" gate) |
| G2 | HERMES Collaboration (A. Airapetian et al.), *First Measurement of the Deuteron Spin Structure Function b1*, Phys. Rev. Lett. **95** (2005) 242001, arXiv:hep-ex/0506018 | `[HERMES05]`, PHYSICS_CHANNELS.md #5; discussed at length in `open_items/physics_literature.md` §3 (RC treatment) and §5 | **Open access on arXiv**, not fetched. The only tensor-DIS *datum* in the world (`BENCHMARK_PLAN.md` T2: "HERMES is the ONLY tensor-DIS datum in the world") | **HIGH** |
| G3 | (a) L. R. Suelzle, M. R. Yearian, H. Crannell, *Elastic Electron Scattering from Li⁶ and Li⁷*, Phys. Rev. **162** (1967) 992 (DOI `10.1103/PhysRev.162.992`); (b) G. C. Li, I. Sick, R. R. Whitney, M. R. Yearian, *High-energy electron scattering from ⁶Li*, Nucl. Phys. A**162** (1971) 583 | Already **one** corpus entry `li-sick-whitney-yearian1971-6li-form-factor` (`file: null`) that folds both citations together; `benchmarking/03_data_nuclear.md` discusses both individually | **No free copy located by the sibling repo**: APS/Elsevier both return 403 to non-browser clients from this environment; (b) may exist as DTIC report AD0714022, also 403'd. `03_data_nuclear.md` marks both "probably table-in-paper, unverified whether tabulated." Origin of the `HoSpin1FF` harmonic-oscillator parameters (⁶Li a=1.9069 fm, α=0.13822; ⁷Li a=1.77(2), α=0.327) currently shipped as **unfitted starting values** | **MEDIUM** |
| G4 | Ent/Mitchell (e,e′d) family — 6 papers on the ⁶Li α–d cluster-knockout channel (list in §3a below) | `benchmarking/03_data_nuclear.md` §7 (item 1–6), `06_critic.md` U-8, `BENCHMARK_PLAN.md` T3 | **None in corpus.** All pre-1998, mostly pre-arXiv; the two primary ones (Ent NPA578, Mitchell PRC44) are explicitly called "digitizable" | **HIGH** for Ent NPA578 (1994) 93 and Mitchell PRC44 (1991) 2002; **MEDIUM** for the other four (superseded/secondary) |
| G5 | R. E. Rand, R. F. Frosch, M. R. Yearian, Phys. Rev. **144** (1966) 859 (⁶Li/⁷Li magnetization) | `06_critic.md` N-6; `BENCHMARK_PLAN.md` T3 item 6 ("the only measured electromagnetic signature of α+d clustering") | **Article paywalled** (pre-arXiv); **a free compilation exists**: de Jager, ADNDT **14** Table V reproduces the same numbers — `06_critic.md` calls it "compilation free / article paywalled" | **HIGH** (½–1 day per `BENCHMARK_PLAN.md`; delivers the only external q→0 slope of F_M and the only ⁷Li magnetization number) |
| G6 | UVa NCD "Quasielastic Electron Nucleus Scattering Archive" Fourier–Bessel compilation, `FB_data.dat` (discovery.phys.virginia.edu) | `06_critic.md` N-5 ("no" on disk, `03` §8 says it does not exist in the tree); `BENCHMARK_PLAN.md` T3 | **Free, machine-readable, not yet fetched** — the sibling repo already downloaded the *sister* file from the same host (the QES archive's ⁶Li 133-point cross-section file, `06Li.dat`, IS committed via the survey though not as a corpus PDF entry), so this should be equally reachable | **HIGH** — `06_critic.md` itself rates this "high" priority, "hours" of effort |
| G7 | Z.-L. Zhou et al., Phys. Rev. Lett. **82** (1999) 687 | PHYSICS_CHANNELS.md's row on the polarised quasi-elastic tail: "named in `phase_C_numbers.md` and nowhere else"; `open_items/physics_literature.md` §3 | **Not in corpus at all** — a deuteron-only measurement/argument that both HERMES's RC treatment and this generator's `qe_tensor_scale=0` default lean on without the paper itself ever having been obtained | **HIGH** — load-bearing citation behind a default the code itself flags as possibly understating the omission by ~10² |
| G8 | R. Machleidt, *High-precision, charge-dependent Bonn nucleon-nucleon potential*, Phys. Rev. C **63** (2001) 024001, arXiv:nucl-th/0006014 | `[CDBonn01]`, PHYSICS_CHANNELS.md #58 (the *last* reference in the list); direct source of `CdBonnWave`'s twenty typed coefficients (Table XX) and published B_d/A_S/η/P_D/Q_d (Table XV) | **Open access on arXiv, simply never fetched** — this is the single most load-bearing reference missing from the corpus: the exact analytic CD-Bonn deuteron wave function transcribed into `cluster.hpp` is one arXiv download away | **HIGH** |
| G9 | E. A. George, L. D. Knutson, Phys. Rev. C **59** (1999) 598 (INSPIRE 522531) | `[GK99]`, PHYSICS_CHANNELS.md #55; the `C-1` consistency band in `benchmarking/00_in_tree_checks.md` | **No arXiv (1999 PRC), presumably APS-paywalled.** Note: `06_critic.md`'s own audit (`C-1`, `U-8`) already **refutes** this as a D/S-ratio benchmark — the measured η it supplies turns out to be an exact relabelling of the code's own quadrupole dial, not an independent wave-function test — so obtaining the PDF matters for correctly re-deriving/re-quoting η = −0.025±0.006±0.010, not for a fresh extraction | **MEDIUM** |
| G10 | N. J. Stone, IAEA INDC(NDS)-0833 (2021, quadrupole moments) and INDC(NDS)-0794 (2019, magnetic dipole moments) | Already corpus entries `stone2021-quadrupole-moments` / `stone2019-magnetic-dipole-moments` (both `file: null`) | IAEA INDC technical reports are normally freely downloadable from `nds.iaea.org` — this looks like a fetch nobody has attempted yet, not a genuine paywall | **MEDIUM** (low risk, likely quick) |
| G11 | K. Slifer et al., JLab PR12-13-011 / E12-13-011 (2013) + 2023 jeopardy update to PAC 51 | Already a corpus entry `jlab-pr12-13-011-b1-and-2023-jeopardy` (`file: null`) | JLab PAC proposals are not normally public documents. **Possible lead**: `CONVENTIONS.md:470` separately cites `arXiv:2506.04506` p.8 for "that experiment's own lower kinematic edge" — worth checking whether 2506.04506 *is* the jeopardy update (or a related published proposal) now posted to arXiv, which would close this gap directly | **MEDIUM** |
| G12 | "The ePIC/EIC documents the efficiencies come from" | far-forward efficiency numbers used throughout PHYSICS_CHANNELS.md and `benchmarking/00` row E-16 | The actual numbers trace to **already-in-corpus, PDF-on-disk** references — `chang2026-ir8-light-nuclei` (arXiv:2511.05638) and `eic-yellow-report2021` (arXiv:2103.05419, split 4 ways). The *documents* still missing are the EIC CDR (BNL-221006-2021-FORE, doi:10.2172/1765663, 972pp/178MB), ePIC PDR (Zenodo 13866213, 282pp/114MB) and ePIC TDR v3.1 (Zenodo 18271601/18271602, 553pp/260MB) — all open access, all deliberately **not** committed on size grounds per the sibling repo's own README (each already has a corpus entry with `file: null` and a DOI/Zenodo id) | **LOW** for the CDR/PDR/TDR themselves (huge; the numbers actually used are already sourced from Chang26 + YR); confirm those two remain the operative sources rather than re-deriving from the CDR |

---

## 3. Additional gaps found by the document survey

The task said the tree's docs cite "~50 more" references beyond the
54-entry corpus. The following subsections enumerate them by source
document. Every arXiv id below was transcribed from the LiPolGen docs'
own text (most of which record "abstract fetched" as part of the
201-resource survey); **none has been independently re-verified over
the network in this pass** — that is the next stage's job.

### 3a. `PHYSICS_CHANNELS.md` §14's own 58-entry reference list, not in the corpus

Cross-referencing every `[KEY]` in `## 14. References` against the 54
corpus arXiv ids found these entries **absent**:

| Key | Reference | arXiv / journal | Note |
|---|---|---|---|
| `[JM89]` | R. L. Jaffe, A. Manohar, *Deep inelastic scattering from arbitrary spin targets*, Nucl. Phys. B **321** (1989) 343 | No arXiv, INSPIRE recid 266865 | **Important finding**: the corpus's *own* entry keyed `hoodbhoy-jaffe-manohar1989-and-jaffe-manohar1989` pairs HJM89 (NPB 312) with a **different** Jaffe–Manohar paper — *Nuclear gluonometry*, PLB **223** (1989) 218 — not this one. NPB 321:343 (the complete arbitrary-spin-J basis, needed for the J=3/2 rank-2 normalisation, `open_items/physics_literature.md` §1) is genuinely missing and is called out repeatedly (`04_theory.md` T-11: "the one genuinely hard-to-get item in the survey") |
| `[POLRAD]` | I. Akushevich, A. Ilyichev, N. Shumeiko, A. Soroko, A. Tolkachev, *POLRAD 2.0*, Comput. Phys. Commun. **104** (1997) 201 | arXiv:hep-ph/9706516 | **Surprising gap** — POLRAD is the RC code this generator's tensor Born/elastic-tail sector is transcribed from (`04_theory.md` T-33 calls it "in-tree, D-3"), yet the paper itself has no PDF or corpus entry. Open access |
| `[WW]` | Wandzura–Wilczek relation (g₂ = −g₁ + WW integral) | No id recorded anywhere in the repo | Standard citation: R. L. Jaffe (or the original) S. Wandzura, F. Wilczek, Phys. Lett. B **72** (1977) 195 — pre-arXiv, should be freely available via INSPIRE scan |
| `[CK90]` | F. E. Close, S. Kumano, Phys. Rev. D **42** (1990) 2377 | No arXiv | The b₁ sum rule ∫b₁dx=0 |
| `[Bissey02]` | F. Bissey, Phys. Rev. C **65** (2002) 064317 | check for arXiv (likely nucl-th/0111034 or similar) | ³He per-nucleon polarizations |
| `[Schell93]` | N. C. Schellingerhout, Phys. Rev. C **48** (1993) 2714 | No arXiv likely | ⁶Li cluster product, named beside a scenario value |
| `[Forest96]` | J. L. Forest et al., Phys. Rev. C **54** (1996) 646 | check for arXiv | Cluster-overlap method |
| `[AV18]` | R. B. Wiringa, V. G. J. Stoks, R. Schiavilla, Phys. Rev. C **51** (1995) 38 | No arXiv (1995) | **Foundational** — the Hamiltonian behind every ANL VMC table the generator reads, and PHYSICS_CHANNELS.md itself flags it as "not cited with a specific reference in the repository" |
| `[Hulthen42]` | L. Hulthén, Ark. Mat. Astron. Fys. **28B** No. 5 (1942) | Pre-arXiv, historical | Analytic radial form, used by name only |
| `[CS96]` | C. Ciofi degli Atti, S. Simula, Phys. Rev. C **53** (1996) 1689 | check for arXiv (nucl-th/9508008-class) | Spectral function n₀/n₁ |
| `[Glauber59]` | R. J. Glauber, *High-Energy Collision Theory*, Lectures in Theoretical Physics Vol. 1 (1959) | Book chapter, no arXiv | Eikonal formalism, used by name |
| `[CT18]` | T.-J. Hou et al., Phys. Rev. D **103** (2021) 014013 | arXiv:1912.10053 | `CT18NLO`/`CT18ANLO` PDF set |
| `[NNPDFpol11]` | E. R. Nocera et al., Nucl. Phys. B **887** (2014) 276 | arXiv:1406.5539 | `NNPDFpol11_100` polarized PDF set |
| `[LHAPDF6]` | A. Buckley et al., Eur. Phys. J. C **75** (2015) 132 | arXiv:1412.7420 | The PDF-grid library itself |
| `[EPPS21]` | K. J. Eskola, P. Paakkinen, H. Paukkunen, C. A. Salgado, Eur. Phys. J. C **82** (2022) 413 | arXiv:2112.12462 | The nPDF set `EPPS21nlo_CT18Anlo_Li6` the EMC baseline uses |
| `[MSTW08]` | A. D. Martin, W. J. Stirling, R. S. Thorne, G. Watt, Eur. Phys. J. C **63** (2009) 189 | arXiv:0901.0002 | The PDF CDKS computed their b₁ with |
| `[PYTHIA83]` | C. Bierlich et al., SciPost Phys. Codebases **8** (2022) | arXiv:2203.11601 | **Core dependency manual** — PYTHIA 8.3, used everywhere, not fetched |
| `[HepMC3]` | A. Buckley et al., Comput. Phys. Commun. **260** (2021) 107310 | arXiv:1912.08005 | **Core dependency manual** — the event-record format the whole chain writes into ePIC |
| `[H1FitB]` | H1 Collaboration (A. Aktas et al.), Eur. Phys. J. C **48** (2006) 715 | arXiv:hep-ex/0606004 | The diffractive PDF PYTHIA's `PDF:PomSet=6` uses |
| `[Fu22]` | D. Fu, B.-D. Sun, Y. Dong, Phys. Rev. D **106** (2022) 116012 | arXiv:2209.12161 | Spin-3/2 GPDs/structure functions |
| `[E143]` | E143 Collaboration (K. Abe et al.), Phys. Rev. D **58** (1998) 112003 | Check for arXiv (distinct from the corpus's R1998/`hep-ex_9808028` entry, which is a *different* E143 paper) | The exact finite-γ target-mass formalism `a_parallel_exact` implements |
| `[Varsh]` | D. A. Varshalovich et al. (Wigner-d/Clebsch–Gordan conventions) | Book, no id | PHYSICS_CHANNELS.md itself: "no bibliographic identifier exists in the repository" |
| `[dFW66]` | T. de Forest, J. D. Walecka | No id recorded | Likely *Electron scattering and nuclear structure*, Adv. Phys. **15** (1966) 1 — needs confirmation |
| `[Moniz71]` | E. J. Moniz et al., Phys. Rev. Lett. **26** (1971) 445 | check for arXiv (unlikely, 1971) | Source of k_F(⁶Li)=0.169 GeV |
| `[eSTARlight]` | M. Lomnitz, S. R. Klein, Phys. Rev. C **99** (2019) 015203 | check for arXiv | The generator paper behind the unpolarized exclusive-VM rate scale (85–97 % of its rate at Q² < 0.1 GeV², below `coherent.hpp`'s q2_min = 0.7; not a coherent-channel baseline) |
| `[CW20]` | W. Cosyn, C. Weiss, Phys. Rev. C **102** (2020) 065204 | arXiv:2006.03033 | Polarized e-d DIS with spectator tagging — the LF polarized spectral function |
| `[BeAGLE22]` | W. Chang et al., Phys. Rev. D **106** (2022) 012007 | arXiv:2204.11998 | BeAGLE generator paper |
| `[CK11]` | C. Ciofi degli Atti, L. P. Kaptari, Phys. Rev. C **83** (2011) 044602 | arXiv:1011.5960 | Cluster-spectator Glauber survival product |
| `[CK04]` | C. Ciofi degli Atti, L. P. Kaptari | arXiv:nucl-th/0407024 | GEA vs. plain Glauber comparison |
| `[CS17]` | W. Cosyn, M. Sargsian, Int. J. Mod. Phys. E **26** (2017) 1730004 | arXiv:1704.06117 | σ_XN(W) Deeps fit |
| `[SW18]` | M. Strikman, C. Weiss, Phys. Rev. C **97** (2018) 035209 | arXiv:1706.02244 | S[IA]+S[FSI]+S[FSI²] and its unitarity sum rule |
| `[Blinov01]` | A. V. Blinov et al., Phys. At. Nucl. **64** (2001) 907 | arXiv:nucl-ex/9910012 | Measured α–p σ_tot, σ_el, B |
| `[Pastore13]` | S. Pastore, S. C. Pieper, R. Schiavilla, R. B. Wiringa, Phys. Rev. C **87** (2013) 035503 | arXiv:1212.3375 | GFMC Q(⁶Li) = −0.20(6) fm² |

(`[Bacchetta07]` — hep-ph/0611265 — is **already in corpus** as
`bacchetta2007-sidis-small-pt` and its PDF is on disk; PHYSICS_CHANNELS.md's
own entry flags only that it is "not registered in the repository's
literature list," a documentation-organization note, not a missing PDF.)

### 3b. The `04_theory.md` T-1…T-46 survey — candidates not covered above

This table (`benchmarking/04_theory.md`) is the single richest source of
candidate references in the tree — a systematic literature survey with
obtainability already checked for each row. Items not already listed in
§2 or §3a:

| Row | Reference | arXiv | What it would give |
|---|---|---|---|
| T-1 | Kumano & Kuroki, *Tensor-polarized PDFs of the deuteron by a convolution model* (2026) | arXiv:2607.09237 | Independent second convolution calculation of b₁ᵈ |
| T-2 | Kumano, PRD **82** (2010) 017501 | arXiv:1005.4524 | Parametrized b₁ᵈ(x,Q²) fitted directly to HERMES |
| T-3 | Khan & Hoodbhoy, PRC **44** (1991) 1219 | No arXiv | Few-parameter closed-form convolution b₁ |
| T-4 | Umnikov, PLB **391** (1997) 177 | arXiv:hep-ph/9605291 | Validity floor for non-relativistic convolution at small x |
| T-5 | Hirai, Kumano, Saito, Watanabe, PRC **83** (2011) 035202 | arXiv:1008.1313 | Closest A>2 unpolarized-SF clustering analogue |
| T-13 | Fu, Dong, Kumano, Xie | arXiv:2602.11587 | Positivity (Soffer-type) bounds for spin-3/2 PDFs |
| T-14 | Zhao, Zhang, Liang, Liu, Zhou | arXiv:2206.11742; arXiv:2401.10031 | Covariant spin-3/2 polarization parametrization (flagged UNVERIFIED in-tree already) |
| T-18 | Detmold & Shanahan, PRD **94** (2016) 014507 + Erratum PRD **95** (2017) 079902 | arXiv:1606.04505 | Lattice-QCD gluonic transversity (A₂≈0.23(2)(5) for φ meson) |
| T-19 | Kumano & Song | arXiv:1910.12523 | Gluon transversity in polarized p-d Drell–Yan |
| T-20 | Wiringa & Schiavilla ("WS98"), PRL **81** (1998) 4317 | arXiv:nucl-th/9807037 | VMC AV18+UIX ⁶Li elastic/transition form factors (figures only, no table) |
| T-22 | Piarulli, Pastore, Wiringa et al., PRC **107** (2023) 014314 | arXiv:2210.02421 | Newer chiral-EFT densities/momentum distributions, A≤12 |
| T-23 | Hebborn, Brune, Phillips (2025) | arXiv:2510.19067 | NCSMC ab-initio ⁶Li correlations ("read the paper first" — abstract does not confirm a D/S or Q result) |
| T-24 | Brida, Pieper, Wiringa, PRC **84** (2011) 024319 | arXiv:1106.3121 | QMC spectroscopic overlaps A≤7 (single-nucleon only, not α–d) |
| T-26 | Cosyn & Sargsian (2011); Cosyn, Melnitchouk, Sargsian (2014) | arXiv:1012.0293; arXiv:1311.3550 | GEA FSI vs. JLab data |
| T-28 | Kaptari, Del Dotto, Pace, Salmè, Scopetta, PRC **89** (2014) 035206 | arXiv:1307.2848 | **Polarized** cluster-spectator distorted spectral function, A=3 |
| T-30 | Sargsian & Strikman, PLB **639** (2006) 223 | arXiv:hep-ph/0511054 | No-loop / pole-extrapolation theorem |
| T-31 | Sargsian, IJMPE **10** (2001) 405 | arXiv:nucl-th/0110053 | GEA foundations |
| T-32 | `github.com/wcosyn/physics-code` | — (code repo) | Second independent FSI implementation |
| T-37 | Afanasev, Bernauer, Blunden, Blümlein et al. | arXiv:2306.14578 | Modern RC review (EPJ A topical review) |
| T-38 | Akushevich & Ilyichev | arXiv:1403.3421 (+related) | NLO RC codes for exclusive photon electroproduction |
| T-39 | Four-fit ⁶Li EMC cross-check: nCTEQ15, nNNPDF3.0, EPPS16 | (PDF sets, LHAPDF grids) | Fit-spread alternative to the single EPPS21 central curve |
| T-40 | ⁷Li nuclear PDFs: `nCTEQ15_7_3`, `TUJU19_*_7_3` | (PDF sets) | Correction: ⁷Li is not PDF-less as currently assumed |
| T-41 | Bissey, Guzey, Strikman et al. (³He g₁ convolution) | — | Nearest thing to a "polarized nuclear PDF" for a light nucleus |
| T-43 | Mäntysaari, Roch, Schenke, Shen, Zhao, PRD **114** (2026) 014068 | arXiv:2605.00454 | Coherent/incoherent J/ψ down to A=3,4 with α-clustering assessment |
| T-44 | Guzey, Rinaldi, Scopetta, Strikman, Viviani, PRL **129** (2022) 242503 | arXiv:2202.12200 | Coherent J/ψ on ⁴He, ³He at EIC |
| T-45 | `github.com/hejajama/subnucleondiffraction` | — (code repo) | The actual code behind Mäntysaari et al. 2408.13213 (Sartre itself is a dead end — see §3c) |

### 3c. Coherent-diffraction extras (`open_items/physics_literature.md` §7)

| Reference | arXiv | Note |
|---|---|---|
| Mäntysaari, Schenke, Shen, Zhao, PRL **131** (2023) 062301 | arXiv:2303.04866 | Multi-scale imaging of nuclear deformation at the EIC |
| Mäntysaari, Salazar, Schenke | arXiv:2207.03712 | Nuclear geometry at small-x |
| Toll & Ullrich (Sartre), CPC **185** (2014) 1835 | arXiv:1211.3048 | **Checked and rejected as a route**: Sartre's `Nucleus.cpp` is a hard-coded `switch` over A∈{1,2,16,27,40,63,90,110,197,208} with no ⁶Li branch and no user-defined nuclei |
| Toll, Ghosh, Srivastav (2026) | arXiv:2606.14633 | The paper whose abstract prompted (and whose source code then disproved) the Sartre-feasibility claim |
| Mondal, Kumar, Sarkar (2026) | arXiv:2608.23445 | Shell structure in coherent \|t\|-differential cross section, QMC densities |
| Magdy et al., EPJ A **60** (2024) | arXiv:2405.07844 | α-clustering EIC work — BeAGLE-based, not diffraction |

### 3d–3h. Data papers named individually in `benchmarking/03_data_nuclear.md`

| Reference | arXiv/DOI | What it gives | Corpus status |
|---|---|---|---|
| Arneodo et al. (NMC), Nucl. Phys. B**441** (1995) 12 | arXiv:hep-ex/9504002, HEPData `ins394050` | Direct, isoscalar F₂(⁶Li)/F₂(D), 24 points — the exact nucleus, currently replaced by an nPDF-fit proxy | **Not in corpus** |
| Lapikás, Wesseling, Wiringa, PRL **82** (1999) 4404 | arXiv:nucl-th/9904008 | ⁷Li(e,e′p)⁶He momentum distributions vs. VMC (spectroscopic factor 0.58±0.05 vs 0.60 VMC) | **Not in corpus** |
| Rodning & Knutson, PRC **41** (1990) 898 (+ PRL **57** (1986) 2248) | No arXiv (1990) | η_d = 0.0256(4), by the same group as George–Knutson, by a different technique (sub-Coulomb tensor analysing powers vs a restricted phase-shift analysis) | **Not in corpus** |
| Whitney, Sick, Ficenec, Kephart, Trower, Phys. Rev. C **9** (1974) 2230 | No arXiv | The full quasi-elastic (e,e′) paper behind the k_F fit (Moniz71 is only the letter) | **Not in corpus** |
| van Niftrik, Lapikás, de Vries, Box, Nucl. Phys. A**174** (1971) 173 | No arXiv | ⁷Li magnetization + quadrupole moment, 25–90 MeV | **Not in corpus** — but `06_critic.md` notes its Table V row is now readable via a compilation, without the paywall |

---

## 4. Summary counts

- Corpus already held: **54** entries (`refs_dict.json`), of which **38**
  have a PDF on disk and **16** are recorded with no free copy (already a
  known gap in the sibling repo, not new — these are the `stone2021`,
  `stone2019`, `jlab-pr12-14-001`, `jlab-pr12-13-011`, `hoodbhoy-jaffe-manohar1989...`,
  `sather-schmidt1990`, `li-sick-whitney-yearian1971`, `ball2003`,
  `koivuniemi2004`, `chaumette1989`, `bultmann1999`, `li2023-eic-tracking-reference-design`,
  `eic-conceptual-design-report2021`, `epic-preliminary-design-report2024`,
  `epic-preliminary-tdr-v31-2026`, and `epic-inclusive-wg-resources` rows
  in §1).
- Gaps confirmed in this pass: **12** explicitly named by the task (§2)
  + **~29** additional named-but-unheld references from
  `PHYSICS_CHANNELS.md` §14 (§3a) + **~26** from the `04_theory.md`
  T-survey (§3b) + **6** coherent-diffraction extras (§3c) + **5** data
  papers (§3d–h) ≈ **~78 distinct gap references**, most with an arXiv id
  already transcribed in the LiPolGen docs and ready to fetch without
  further searching; a minority (pre-1994 nuclear/atomic-physics papers:
  Ent, Mitchell, Rand–Frosch–Yearian, Suelzle, Li–Sick, Rodning–Knutson,
  Whitney PRC9, van Niftrik, AV18, Schellingerhout, Forest96, Hulthén,
  Glauber59) predate arXiv and will need INSPIRE/journal-DOI handling or
  a compilation-table workaround (as already found for Rand–Frosch–Yearian
  and George–Knutson-adjacent items).
