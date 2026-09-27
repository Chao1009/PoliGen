<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Decisions of 2026-09-27

Taken on the author's delegation ("decide for me … make a brief note").
Each is reversible; the evidence is in `../run_2026-09-26/SUMMARY.md` and
`../../benchmarking/REPORT.md`.

**D1 — HepMC3 inclusive balance: scope the document, do not write a remnant.**
`docs/HEPMC3_CONVENTION.md` stated its conservation identity for every
channel, while inclusive records deliberately omit the (A−1) remnant
(`docs/USAGE.md` §2). The generator models the struck nucleon, not the
residual nucleus; writing a remnant would put an invented object (no
binding, no excitation, no fragmentation) on every inclusive file and into
everything downstream of it. The document now says an inclusive record
balances per nucleon, e + struck nucleon (its lines 140-147). No generator
output changed. `t6_hepmc3_convention`: FAIL 21/23 → PASS 23/23; tagged and
coherent files are still held to the full beams.

**D2 — the Hulthén deuteron stays the default; physics numbers use AV18.**
`t2_deuteron_static` finds the Hulthén pair's asymptotic D/S ratio
η = 0.0325, +17σ from the measured 0.0256(4); CD-Bonn (−0.07σ) and AV18
(−1.4σ) agree. The default is kept because every reference pin, the
polligen port gate and the knob-provenance tables are defined on it, and
moving it is a release-wide change of every tagged number. In exchange:
any tagged-channel number quoted as physics (the letter, a plot, a table)
is produced with `--cluster-wave vmc` (the AV18 tables, the configuration
the Cosyn–Weiss gate already runs on), and the row stays a recorded FAIL
so the gap stays visible. Revisit at M2 (numbers frozen): switching the
default then costs one re-pin, and should be done if any quoted number
still depends on the Hulthén D wave.

**D3 — a new benchmark-plan row: the CLAS Deeps spectator tail.**
`validation/benchmarks/data/clas_deeps_klimenko2006_f2P.json` was vendored
on 2026-09-23 but not wired, because it needed a new plan row. The row is
now added (`BENCHMARK_PLAN.md` §4, row 11) and wired as
`t2_deeps_spectator_tail`. The blocked rows (6, 7, 10) and the documents
they wait for are left to the author.
