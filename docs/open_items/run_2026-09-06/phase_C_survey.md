# Phase C survey — a central-detector lepton reconstruction efficiency in `../PolarizedLithiumSim`?

Run 2026-09-06; this survey performed **2026-09-15**.  **Read-only survey of
the sibling repository.**  *(What LiPolGen then did with this result — task C1
took `PLAN.md`'s second branch, "say so and stop" — is recorded in
`phase_C_numbers.md`; the §C-S8 statement below that nothing in LiPolGen was
modified is scoped to this survey step, not to C1.)*  Nothing in
`PolarizedLithiumSim` was modified; nothing in LiPolGen was modified except
this file.

## C-S1 The question, stated exactly

`validation/o5_a2_reach.py` builds the O5 detection-efficiency chain on
`COHERENT_JPSI_EFF_IR8_LI7` and states, in its own assumption list
(`validation/o5_a2_reach.py:90-96`):

> `... it is a RECOIL-NUCLEUS number, so the chain sigma x BR(l+l-) x eps_det`
> `has NO central-detector acceptance and NO reconstruction efficiency for the`
> `DECAY LEPTONS.  ...  the GEOMETRIC part is bounded and is NOT the problem`
> `(0.9940 at this sample's <W> ...); the per-lepton reconstruction efficiency`
> `is UNBOUNDED HERE, and the 0.95/track used at the band's low end is the`
> `sibling ../PolarizedLithiumSim's own stand-in.`

The quantity wanted is therefore ONE number: **ε_reco per decay lepton
(e or μ) of a coherent J/ψ in the ePIC central detector**, given that the
lepton is already inside the acceptance.  The pair costs its square
(`EPS_TRACK_SIBLING_STANDIN ** 2` at `validation/o5_a2_reach.py:528`).

The brief asked for five shapes of candidate: a tracking efficiency per track;
a J/ψ → ℓℓ reconstruction efficiency; an `eic-smear` parameterisation; a
fullsim output with reconstructed J/ψ; an ePIC Yellow-Report number the project
already assumes.  Each is answered below.

## C-S2 What was searched

Every module of `fastsim/polli_fastsim/` (11 files, 3725 lines:
`asymmetries.py`, `beams.py`, `delta_models.py`, `farforward.py`, `fom.py`,
`__init__.py`, `inputs.py`, `kinematics.py`, `polarized.py`, `spectator.py`,
`structure.py`); all of `tools/fullsim/` (`README.md`, `ff_gun_scan.sh`,
`ff_gun_hits.py`, `ion_gun_hepmc.py`); all of `tools/analysis/`
(`dump_spectators.py`, `ed_control_analysis.py`);
`plans/03_phase2_full_simulation.md`; `plans/09_nearbeam_nanowire_far_forward.md`;
all of `refs/` (`README.md`, `refs_dict.json` — 54 entries — and the PDFs
`2210.09048.pdf` and `2511.05638.pdf`, read directly with `pdftotext`).

Because the brief's own pointer ("a 0.95/track labelled as the sibling's own
stand-in") did not resolve inside that set, the search was widened to the whole
sibling tree with repository-wide greps for `0.95`, `per.track|per-track|
track.eff|tracking effic`, `efficien`, `j/psi|jpsi|dilepton|muon|vector meson`,
`eic-smear|eic_smear|delphes|fun4all|craterlake|npsim|ddsim|eicrecon|juggler`,
and `eff_|_eff\b|eID|lepton`.  That is what found it, in `evgen/`, not in the
directories the brief named.

## C-S3 The candidate list

| # | where | the number | shape | class |
|---|---|---|---|---|
| 1 | `evgen/polligen/hfs.py:251` | `eff_track` = **0.95** plateau, per charged track | per-track tracking efficiency | **placeholder** — the file says "stand-in" |
| 2 | `evgen/polligen/reco.py:482-512` | `eps_eid(η)`, 9 anchors, **0.95** at η = −2.0 | scattered-electron ID efficiency | **placeholder** — "a CONSTRUCTED eta profile, not a published curve" |
| 3 | `fastsim/scripts/money_delta_20260729.py:798-835` | the same 9 anchors, duplicated | same | **placeholder** (second copy of #2) |
| 4 | `evgen/polligen/hfs.py:259-260` | `eff_nhad` = 0.9, neutral hadrons above threshold | calorimetric detection efficiency | **placeholder** — "stand-in" |
| 5 | ATHENA, JINST 17 (2022) P10019 Table 5, in `refs/2210.09048.pdf` | **95%** electron efficiency at > 99.8% π rejection, p ≥ 0.1 GeV/c, bECal | single-electron calorimetric ID efficiency | **measurement** (a published standalone-calorimeter simulation result) |
| 6 | ECCE, NIM A 1055 (2023) 168464 §3.5.2, cited at `evgen/polligen/reco.py:496-498` | EEMC ~95%, **BEMC ~70%**, FEMC ~90–95% | η-resolved electron ID efficiency | **unverifiable here** — no PDF in `refs/`, no `refs_dict.json` entry |
| 7 | `refs/README.md:73`, `refs/refs_dict.json:433-434`, `refs/2511.05638.pdf` | ⁷Li **17.75%** at IR-8 | intact-**recoil** far-forward number | **measurement of a different quantity** — and see C-S6 |
| 8 | ePIC preTDR v3.1 Fig. 3.56 / PDR Fig. 8.9 (`refs/refs_dict.json`, `epic-preliminary-tdr-v31-2026`) | dp/p = 3.00 / 1.08 / 0.38 / 1.13 / 2.55 % at p = 1 GeV/c | tracking **resolution**, not efficiency | **measurement of a different quantity** |
| 9 | `plans/09_nearbeam_nanowire_far_forward.md:223-226, 243` | 0.950 / 0.99 | far-forward Roman-Pot Z-ID working point | **not central-detector, not lepton** |
| 10 | EIC Yellow Report, `refs/refs_dict.json` (`eic-yellow-report2021`) | "10 fb⁻¹ = 30 weeks at 1e33 with **60% efficiency**" | machine uptime | **not a detector efficiency** |

### C-S3.1 Candidate 1 — the "0.95/track", found and read

`evgen/polligen/hfs.py`, class `HadronResponse`.  Its parameter block opens
(`hfs.py:241-243`):

> `Parameters.  Which are sourced and which are this programme's own`
> `stand-ins (labelled on 2026-08-28; the earlier docstring called them`
> `all "Yellow Report requirement values", which only the resolution tables are)`

and the row itself is `hfs.py:251`, in full:

> `      eff_track      plateau tracking efficiency (0.95, stand-in)`

The value and its label are on the same line, and the label is the sibling's
own.  The module docstring says the same thing again at `hfs.py:18-22`: the **resolution** tables are the Yellow Report
requirement values, while "the coverage, thresholds, efficiencies and the noise
are this programme's own stand-ins of the right magnitude ... to be replaced by
ePIC design values when available."

The default is set at `hfs.py:293` and applied at `hfs.py:360-362`:

```
eff = self.eff_track / (1.0 + np.exp(-(pt - self.pt_min_track) / self.pt_turn))
eff = np.where(pt < 0.5 * self.pt_min_track, 0.0, eff)
tracked = in_track & (rng.uniform(size=n) < eff)
```

So 0.95 is the **plateau** of a logistic p_T turn-on centred at
`pt_min_track` = 0.2 GeV with width 0.05 GeV, gating any charged particle
inside `|η| ≤ eta_track` = 3.5.  Three properties matter for how LiPolGen uses
it:

* it is applied to **hadronic-final-state particles** under a pion mass
  hypothesis (`hfs.py:372`), not to identified leptons;
* the J/ψ decay leptons carry p* = 1.548 GeV, far above the 0.2 GeV turn-on,
  so the plateau value is the value that would apply — LiPolGen's use of 0.95
  rather than the turned-on value is at least self-consistent;
* it carries **no PID, no material budget, no magnetic field and no fake or
  purity term**, and it is uniform in η across the entire |η| ≤ 3.5 coverage.

**Stated source: none.**  It is the sibling's own number of the right order of
magnitude.  As a bound it has no direction: nothing says a real ePIC per-lepton
efficiency is above it or below it.

### C-S3.2 Candidate 2 — `eps_eid(η)`, the only lepton-side efficiency in the tree

`evgen/polligen/reco.py:482-512`.  Its docstring is unusually explicit and is
worth quoting because it settles the "ePIC Yellow Report number the project
already assumes" question (`reco.py:485-487`):

> `Re-documented 2026-08-27.  This is a CONSTRUCTED eta profile, not a`
> `published curve, and it is anchored at exactly one point.  No ePIC`
> `electron-ID efficiency curve exists in any ePIC document.`

The anchors are `[0.85, 0.92, 0.95, 0.93, 0.90, 0.90, 0.85, 0.80, 0.70]` at
η = [−3.5, −3.0, −2.0, −1.0, 0, 1, 2, 3, 3.5] (`reco.py:509-510`), applied at
`evgen/polligen/recopseudo.py:278-281` to the **scattered DIS electron** and to
nothing else.  The docstring names its own defects (`reco.py:499-501`): the
0.95 at η = −2.0 is ECCE's EEMC value, "its barrel 0.90 matches neither source,
and its forward tail (0.85/0.80/0.70) runs OPPOSITE to ECCE's FEMC."  The
sibling's own consistency review flags the same thing as finding **F204**
(`docs/consistency_review_2026-09-02.md:346`).

This is the closest thing in the tree to a lepton reconstruction efficiency,
and it is a placeholder for a different lepton (the scattered electron, not a
decay lepton) with no tracking term in it at all.

### C-S3.3 Candidate 5 — ATHENA Table 5, the one published, verified number

Read directly from `refs/2210.09048.pdf` on 2026-09-15 (`pdftotext -f 20 -l 32`),
Table 5, "Expected bECal detector performance", row *e/π separation*, verbatim:

> `> 99.8% pion rejection with 95% electron efficiency at p ≥ 0.1 GeV/c`

with footnote *b*: "Based on simulation for a standalone bECal, see Fig. 9 for
detailed results", and the Fig. 9 caption, verbatim:

> `All the curves, including simulations and data, are obtained for the`
> `standalone calorimeter, i.e., no other materials are placed in front of the`
> `calorimeter and no magnetic field is involved.`

This is a **measurement** in the sense the brief means — a published simulation
result with a stated configuration — and it is the only one in the sibling's
whole reference corpus that is about identifying a lepton.  Five things stand
between it and the quantity the O5 chain needs, and every one of them has to be
said at any site that quotes it:

1. **It is ATHENA, not ePIC.**  ATHENA was a proposal for IP6 and was not
   built.  That ePIC's barrel imaging calorimeter inherits the bECal's
   performance is an inference; nothing in this tree or the sibling's
   establishes it.
2. **It is a standalone calorimeter** — no material in front, no magnetic
   field, by the paper's own caption.  Both omissions move a real efficiency
   **down**, so as a bound it is an **upper** bound, under the assumption that
   adding material and a 1.7 T field cannot raise a calorimetric ID efficiency.
3. **It is the barrel only.**  Table 5 prints no η range at all.  §2.4.2 gives
   the geometry — 12 staves, inner radius 103 cm, imaging layers 405 cm long —
   and adds that "the bECal not only functions as the barrel calorimeter, but
   also provides significant coverage in the electron-going direction".  Read
   as a symmetric barrel those two numbers give |η| ≲ 1.4; that reading is an
   **inference**, not a quoted attribute, and the asymmetry the same sentence
   describes is not in it.  What is certain is that Table 5 says nothing about
   the endcaps.
4. **It is an ID efficiency given a cluster**, not a reconstruction efficiency:
   there is no tracking term, no vertexing and no pair-finding in it.
5. **It is an electron number.**  For μ⁺μ⁻ the paper offers one qualitative
   sentence and no number (`2210.09048.pdf` §2.4.2: "the 3-D shower profiles
   measured by the bECal enable effective μ identification").  Nothing
   anywhere in the sibling quantifies a muon.

Point 3 is fatal for this sample in particular.  At the O5 sample's own
⟨W⟩ = 30.2 GeV the J/ψ sits at **y = −0.39**, and running
`validation/o5_a2_reach.py`'s own `lepton_pair_acceptance` at successively
tighter `eta_max` (SCHC transverse decay, `schc=True`, the pessimistic
weighting the function defaults to) gives the fraction of pairs with **both**
leptons inside:

| both leptons inside | fraction of pairs at ⟨W⟩ = 30.2 |
|---|---|
| \|η\| < 1.0 (barrel-ish) | **0.4477** |
| \|η\| < 2.0 | 0.8888 |
| \|η\| < 3.5 (`Scenario::eta_max`) | 0.9940 |

So an ATHENA-barrel number covers **45% of the pairs** and is silent on the
rest.  It cannot be extended to the other 55% without inventing the endcap
values — which is exactly what `eps_eid` (candidate 2) already does, and
labels as invention.

### C-S3.4 Candidate 6 — the number that would move the band, and cannot be checked

`evgen/polligen/reco.py:496-498` asserts, of ECCE (NIM A 1055 (2023) 168464
§3.5.2): "EEMC (−3.4 < η < −1.5) ~95%, BEMC (−1.72 < η < 1.31) **~70%**, FEMC
(1.3 < η < 3.5) ~90–95%."

`refs/2207.09437` **is not in `refs/`** and there is no `refs_dict.json` entry
for ECCE (checked: all 54 entries were enumerated; none is ECCE).  So the
one η-resolved electron-ID curve the sibling names cannot be verified here.

It is recorded anyway because of its direction.  **If** any proposal-era
electron-ID number were adopted as the per-lepton stand-in, the one that
applies to the barrel — where 45% of these pairs land — is **0.70, not 0.95**.
The 0.95/track LiPolGen currently uses to draw the band's low end is therefore
not even the pessimistic choice among the numbers the sibling itself names.
That does not make 0.70 a bound (it is unverifiable here, it is ECCE not ePIC,
and it is ID not reconstruction) — it makes the **openness below the band's low
rung larger than the current text implies**, which is a direction, not a
number.

## C-S4 The negative findings, each with its evidence

**There is no `eic-smear` parameterisation in the sibling.**  `eic-smear`
appears once, at `tools/beagle/README.md:109`, purely as a **format converter**
in the BeAGLE chain ("fort.92? → eic-smear `BuildTree`/`TreeToHepMC` (run
inside eic-shell) → HepMC3").  The only detector parameterisations reachable
through it are the Yellow-Report models in the external
`JeffersonLab/dis-reconstruction` repository, which `refs/refs_dict.json`
(`epic-inclusive-wg-resources`) records as a **URL only** — the repository is
not vendored, not mirrored, and not read by anything in the tree.

**There is no fullsim output with a reconstructed J/ψ — there is no fullsim
output with a reconstructed anything.**  Three independent statements, all the
sibling's own:

* `plans/03_phase2_full_simulation.md:62-63` — "The one clause never exercised
  is the literal `eicrecon` reconstruction command — **no EICrecon output
  exists anywhere in the tree**".
* `plans/03_phase2_full_simulation.md:224` — "**No central-detector
  npsim/EICrecon run exists anywhere in the tree.**"
* `plans/03_phase2_full_simulation.md:244-245` — "no reconstructed sample has
  ever come out of EICrecon (**`eicrecon` is invoked in zero scripts**)".

And `tools/fullsim/` matches that: it is an npsim **hit-level** far-forward gun
scan and nothing else.  `tools/fullsim/README.md:607` states its caveats —
"hit level only (no reconstruction)" — and `:92` says the same of the earlier
scan.  `ff_gun_hits.py` reads `B0Tracker*Hits`, `B0ECal*Hits` and their
siblings (`ff_gun_hits.py:27`); there is no track collection, no PID collection
and no central-detector collection anywhere in it.  A search of the whole
sibling for `*.root` and `*.hepmc*` returns **nothing** — no simulation output
is committed at all.

**There is no J/ψ reconstruction, and no decay leptons, anywhere.**  Every
J/ψ mention in the sibling is one of three things: (a) the Mäntysaari *et al.*
coherent-deuteron a₂ calculation, digitized into
`evgen/polligen/coherent.py:221` as `MANTYSAARI_A2_DEUTERON` — a theory curve
in the recoil azimuth, no final state; (b) the Chang *et al.* IR-8 recoil
tagging number (C-S6); (c) prose in `plans/06` and
`docs/note_cos2phi_coherent_6Li.md`.  `evgen/polligen/coherent.py` models the
coherent channel as `e + ⁶Li → e' + X + ⁶Li(g.s.)` — the **recoil** and the
scattered electron.  The vector meson's decay is not in the module.

**`fastsim/polli_fastsim/` contains no efficiency of any kind.**  A grep for
`efficien` across all 11 modules returns exactly four hits, every one of them
the substring inside the word *coefficient* (`delta_models.py:27` and `:84`,
`farforward.py:715`, `structure.py:81`) — not one is an efficiency.  The
package's only detector modelling is the geometric acceptance of `fom.py`
(`y_max = 0.95` at `fom.py:68`, |η_e| ≤ 3.5, E' ≥ 0.5 GeV) and the far-forward
angular windows of `farforward.py`.  `tools/analysis/` likewise: a grep for
`efficien|lepton|jpsi|track` across both of its scripts returns nothing.

**`plans/09` is far-forward only.**  A grep for `central detector|
central-detector|j/psi|lepton|reconstruction efficiency` over all 63 kB of it
returns **nothing**.  Its 0.950 figures (`plans/09:223-226`) are the matched
⁶Li **efficiency working point** of a Roman-Pot Z-ID fake-rate study, and its
0.99 (`plans/09:243`) is an assumed per-plane silicon hit efficiency in that
same study.  Neither is a central-detector number, and neither is about a
lepton.

## C-S5 Incidental finding (1): LiPolGen's citation names a class that does not exist

`validation/o5_a2_reach.py:466-468` (as of 7f68339) reads:

```
#: sibling ../PolarizedLithiumSim's own `HfsModel(eff_track=0.95)`, which that
#: file labels "stand-in" in as many words.  The pair costs its SQUARE.
EPS_TRACK_SIBLING_STANDIN = 0.95
```

and `:436-437` cites `HfsModel(pt_min_track=0.2)` the same way.

**There is no `HfsModel` in `PolarizedLithiumSim`.**  A grep for
`HfsModel|HadronResponse|class Hfs` over every `.py` in the sibling returns 26
hits, all of them `HadronResponse` (in `evgen/polligen/hfs.py`,
`evgen/scripts/hfs_*.py`, `evgen/scripts/money_cos2phi_reco.py:302`,
`evgen/tests/test_hfs.py`, `tools/pythia8/gen_dis_hfs.py:51`) and zero
`HfsModel`.  The correct citation is
**`polligen.hfs.HadronResponse(eff_track=0.95)`**, `evgen/polligen/hfs.py:251`
for the label and `:293` for the default.

The *substance* of LiPolGen's claim is exactly right — the file does label it
"stand-in" in as many words, at `hfs.py:251` — so this is a naming defect in
the citation, not a factual one.  But it is unresolvable as written: a reader
who greps the sibling for `HfsModel` finds nothing and cannot confirm the
label.  Two sites carry it (`o5_a2_reach.py:437` and `:466`);
`docs/OPEN_ITEMS_SOLUTIONS.md:2010-2012` repeats it a third time (re-pointed 2026-09-16, eight lines down, by the close-out's own edits to that file's state board).

## C-S6 Incidental finding (2): the 17.75% carries no detector efficiency either

This is the survey's one substantive result and it does not concern the
leptons.

`COHERENT_JPSI_EFF_IR8_LI7` = 0.1775 is the single largest multiplier in the
O5 chain.  The sibling records it twice — `refs/README.md:73` ("IR-8
intact-recoil efficiencies ... **17.75** ... for ... ⁷Li") and
`refs/refs_dict.json:433-434`, the latter as "Intact-recoil tagging
**efficiency x acceptance** at the IR-8 secondary focus".  LiPolGen's own
docstring (`validation/o5_a2_reach.py:38-39`) calls it "arXiv:2511.05638's ⁷Li
coherent-J/psi far-forward **efficiency**".

Read from `refs/2511.05638.pdf` directly on 2026-09-15, the last sentence of §IV
(*Event Generator*, p. 3→4), verbatim:

> `The current simulation only accounts for the acceptance effect and does not`
> `incorporate the efficiencies of the detector.  Additionally, we did not`
> `account for the efficiency and acceptance of the reconstructed distribution,`
> `which is to showcase our detection performance at its detector level form.`

The 17.75% is then quoted in §V (*Results*, p. 4) as the "global detection
efficiency" for ⁷Li.

So the number is a **pure geometric acceptance** of the intact recoil through
the IR-8 far-forward optics, after the EIC afterburner's crossing angle,
angular divergence and momentum spread — and it carries, by the paper's own
statement, **no Roman-Pot detector efficiency, no reconstruction efficiency and
no reconstruction acceptance**, on the recoil leg *or* anywhere else.

This sharpens rather than contradicts LiPolGen's assumption list, which says
the chain has no decay-lepton acceptance or efficiency.  True — and the recoil
leg has no detector efficiency either.  Both omissions point **down**.

**THIS IS A BOUND, and it is stated as one:**

> **BOUND.**  `COHERENT_JPSI_EFF_IR8_LI7` = 0.1775 is an **UPPER bound** on the
> far-forward intact-⁷Li tagging efficiency × acceptance at the IR-8 secondary
> focus.
> **Source:** W. Chang, E.-C. Aschenauer, A. Jentsch, A. Kumar, Z. Tu, Z. Yin,
> PRD 113 (2026) 032018 (`refs/2511.05638.pdf`), §IV last sentence for the
> exclusion, §V for the number.
> **The assumption that makes it a bound:** that a Roman-Pot detector
> efficiency, a reconstruction efficiency and a reconstruction acceptance are
> each ≤ 1, so restoring any of them can only lower the product.  Nothing here
> quantifies how far: the paper does not, and this tree does not.

It is a bound on the **recoil** leg.  It is **not** a bound on the decay
leptons, and it must not be quoted as one.

## C-S7 Verdict

**On the question asked — is there a defensible BOUND on the per-lepton
central-detector reconstruction efficiency at ePIC in `PolarizedLithiumSim`? —
the answer is NO.**

The honest deliverable is the sentence, and here it is with what was looked at:

> **Every module of `fastsim/polli_fastsim/`, all of `tools/fullsim/` and
> `tools/analysis/`, `plans/03` and `plans/09`, and all 54 reference entries
> plus two reference PDFs of `refs/` were read, and the whole sibling tree was
> grepped for the five shapes the brief named.  `PolarizedLithiumSim` contains
> no measurement, no parameterisation and no assumed value of a
> central-detector lepton reconstruction efficiency at ePIC.  It contains one
> per-track efficiency, `HadronResponse.eff_track` = 0.95 at
> `evgen/polligen/hfs.py:251`, which that file labels "stand-in" on the same
> line that states it; one constructed scattered-electron ID profile,
> `reco.eps_eid` at `evgen/polligen/reco.py:482`, whose own docstring says "No
> ePIC electron-ID efficiency curve exists in any ePIC document"; and no J/ψ
> reconstruction, no `eic-smear` parameterisation, and no EICrecon output of
> any kind — `eicrecon` is invoked in zero scripts (`plans/03:244`).**

Three qualifications, so that the "no" is not read as more than it is:

* **One conditional bound exists and is too narrow to use.**  ATHENA Table 5's
  95% (C-S3.3) is an **upper** bound on the calorimetric electron-ID efficiency
  of a barrel EMCal of the ATHENA bECal design, under the assumption that
  material in front and a magnetic field can only lower it.  It covers 45% of
  these pairs (both leptons within |η| < 1 at ⟨W⟩ = 30.2), contains no tracking
  term, and has no muon counterpart.  It does not bound ε_pair.
* **A bound was found on a neighbouring quantity** — C-S6, the 17.75% — and it
  points the same way (down) as the missing lepton factor.  It does **not**
  cancel the ×1.1224 beam-energy correction the way the lepton factor partly
  does; it compounds with it.
* **The direction below the band's low rung is now better known.**  0.95/track
  is the *optimistic* end of the numbers the sibling itself names; the barrel
  figure it attributes to ECCE is 0.70 (C-S3.4).  That is a direction, not a
  number, and its source is unverifiable in this tree.

**What this survey changes in LiPolGen: nothing numeric.**  The O5 chain's
statement that the per-lepton reconstruction efficiency is UNBOUNDED HERE
survives this survey intact and is **confirmed** by it — the sibling does not
supply one either, and the survey's job was to find out.  Two things are now
available to whoever acts next, and both are author decisions, not this
survey's to take:

1. the citation fix of C-S5 (`HfsModel` → `polligen.hfs.HadronResponse`,
   `evgen/polligen/hfs.py:251`), at `validation/o5_a2_reach.py:437` (as of 7f68339) and `:466` (as of 7f68339)
   and at `docs/OPEN_ITEMS_SOLUTIONS.md:1991-1993` (as of af3f415);
2. the sharpened reading of C-S6 — that 0.1775 is acceptance-only on the recoil
   leg too, by the paper's own §IV sentence — which, if adopted, would have to
   land at every site the ground rules name:
   `validation/o5_a2_reach.py` (the chain and its self-check),
   `docs/open_items/run_2026-09-03/phase_C_numbers.md` §C2,
   `docs/OPEN_ITEMS_SOLUTIONS.md` §11.3,
   `include/lipolgen/cluster_config.hpp` (`a2_from_geometry`),
   `python/lipolgen/configs.py`, the Mäntysaari draft, `SUMMARY.md`, and the
   registry §B25.

Neither was applied.  This file is the survey.

## C-S8 Baseline, measured before and after

The tree was clean at start and the only file added is this one, which lies
outside both gates' covered sets (`check_spdx_headers.py`'s `GROUPS` is
`include/`, `src/`, `tests/`, `python/bindings.cpp`, `python/lipolgen/*.py`,
`validation/*.py`; `check_physics_channels_links.py` covers
`docs/PHYSICS_CHANNELS.md`, `docs/theory/SPIN32_FINITE_GAMMA.md` and
`docs/PYTHIA_BRIDGE.md`).  Both were nevertheless re-run.

| | baseline | after |
|---|---|---|
| `build/lipolgen_tests` | 408/408 cases, 17 241 975 assertions, 3m17s | unchanged (no code touched) |
| `python -m pytest python/tests -q` | 1023 passed, 151 skipped, 261.8s | unchanged (no code touched) |
| `check_spdx_headers.py` | 102 files, all carry the header, rc 0 | 102 files, rc 0 |
| `check_physics_channels_links.py` | 1236 / 19 / 0 references strict, 96 / 115 / 0 ranges, **0 broken**, rc 0 | identical, rc 0 |
