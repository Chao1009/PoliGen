# Run 2026-09-02 — status

Cron job id: a07c52e5 (session-only, hourly at :11, first fire 2026-09-02 18:11 CDT,
auto-expires after 7 days) -- **deleted 2026-09-03 in the Phase H close-out**.
Started: 2026-09-02 17:05 CDT (planning); work ran from the first cron fire after
the usage-limit reset. **Closed out 2026-09-03.**

| phase | state | commit | notes |
|---|---|---|---|
| A physics-channels doc | DONE | 6f67d84 | docs/PHYSICS_CHANNELS.md 654 lines / 111 rows / 48 refs; gate 826 refs 0 broken (validation/check_physics_channels_links.py, --fix mode); README linked; refresh line refs + add C/D/G rows in Phase H |
| B VMC tag fractions + writer fix | DONE | 8c41899 | tag fractions at 5x41/10x100/18x275 (USAGE), beta band retired, generated-mass rule pinned (doctest + pytest); 276 doctest / 146 pytest green |
| C tensor RC band | DONE | fd339e2 + ce86a72 | rc.hpp/rc.cpp weight family (--rc tensor-band), POLRAD t-peak tails with the transcription check, 1/A^2 blocker caught in review; 314 doctest / 164 pytest green; numbers in phase_C_numbers.md and OPEN_ITEMS section 9 |
| D b1(6Li) convolution | DONE | b1071b1 | b1_nuclear.hpp four-term convolution, --b1-model li6-convolution (opt-in; miller default pinned at rtol 1e-12); A = 2 gate: shape passes, magnitude factor 2.27 low (open); 333 doctest / 183 pytest green |
| E spin-3/2 note | DONE | 472a7d4 + 0961c60 | docs/theory/SPIN32_FINITE_GAMMA.md (1195 lines, 2 review passes); g1_rank3/R_3 code comments; generate_full pure-fill pzz_true fixed |
| F packaging | DONE | 611e40a | pyproject + scikit-build-core; fresh-venv editable install 66 s, wheel 1.31 MB, in-tree flow unchanged; README snippet applied to README.md in Phase H |
| G coherent 6Li on-ramp | DONE | 6e393aa + 198c846 + a764383 | eSTARlight baseline (6Li J/psi 1.77 nb, slope_b band citable), reviewed sampler design, cluster_config.hpp sampler + lipolgen-configs + OPEN_ITEMS section 11 with the collaboration ask; 361 doctest / 208 pytest |
| H close-out | DONE | (this commit) | PHYSICS_CHANNELS.md refresh (1031 refs, 0 broken), status docs, whole-diff review + the 26 findings applied. Final: **363 doctest cases / 17 203 863 assertions / 1 skipped** (wall 1:17, down from 2:24 — the b1 A = 2 scan is threaded and T2 no longer rebuilds its densities three times) and **211 pytest** (32 s, down from 43 s). Cron job a07c52e5 deleted. |

States: TODO / IN PROGRESS (with what exists on disk) / DONE (with commit hash).

## Deliverables

Every file created in this run (`git diff --name-status 2e404b8 HEAD`, status `A`),
by phase.

Planning:

- `docs/open_items/run_2026-09-02/PLAN.md`
- `docs/open_items/run_2026-09-02/STATUS.md`

Phase A — physics-channels doc:

- `docs/PHYSICS_CHANNELS.md`
- `validation/check_physics_channels_links.py`

Phase C — tensor RC band:

- `include/lipolgen/rc.hpp`
- `src/core/rc.cpp`
- `tests/test_rc.cpp`
- `tests/test_rc_pipeline.cpp`
- `python/tests/test_rc.py`
- `docs/open_items/run_2026-09-02/design_C_tensor_rc.md`
- `docs/open_items/run_2026-09-02/phase_C_numbers.md`
- `docs/open_items/run_2026-09-02/polrad_transcription_check.md`

Phase D — b1(6Li) convolution:

- `include/lipolgen/b1_nuclear.hpp`
- `src/core/b1_nuclear.cpp`
- `tests/test_b1_nuclear.cpp`
- `python/tests/test_b1_model.py`
- `validation/b1_li6_table.py`
- `validation/dump_b1_default_li6.py`
- `validation/reference/b1_default_li6.json`
- `docs/open_items/run_2026-09-02/design_D_b1_li6.md`
- `docs/open_items/run_2026-09-02/phase_D_numbers.md`
- `docs/open_items/run_2026-09-02/phase_D_gate.md`

Phase E — spin-3/2 note:

- `docs/theory/SPIN32_FINITE_GAMMA.md`

Phase F — packaging:

- `pyproject.toml`
- `docs/open_items/run_2026-09-02/README_packaging_snippet.md`

Phase G — coherent 6Li on-ramp:

- `include/lipolgen/cluster_config.hpp`
- `src/core/cluster_config.cpp`
- `tests/test_cluster_config.cpp`
- `python/lipolgen/configs.py`
- `python/tests/test_cluster_config.py`
- `docs/open_items/run_2026-09-02/design_G_cluster_config.md`
- `docs/open_items/run_2026-09-02/phase_G_numbers.md`
- `docs/open_items/run_2026-09-02/estarlight_li6.md`
- `docs/open_items/run_2026-09-02/estarlight/li6_jpsi_10x99.5.in`
- `docs/open_items/run_2026-09-02/estarlight/li6_phi_10x99.5.in`
- `docs/open_items/run_2026-09-02/estarlight/li6_rho_10x99.5.in`
- `docs/open_items/run_2026-09-02/estarlight/li7_jpsi_10x99.5.in`
- `docs/open_items/run_2026-09-02/estarlight/li7_jpsi_18x117.9.in`
- `docs/open_items/run_2026-09-02/estarlight/li7_phi_18x117.9.in`
- `docs/open_items/run_2026-09-02/estarlight/li7_rho_18x117.9.in`
- `docs/open_items/run_2026-09-02/estarlight/nucleus_li_radii.patch`
- `docs/open_items/run_2026-09-02/estarlight/q7_li6_jpsi.in`
- `docs/open_items/run_2026-09-02/estarlight/rr_li6_jpsi_10x99.5.in`
- `docs/open_items/run_2026-09-02/example_li6_m1_configs.dat`
- `docs/open_items/run_2026-09-02/example_li6_m1_configs.dat.meta.json`
- `data/vmc/density/he4.density`
- `data/vmc/density/li6.density`
- `data/vmc/density/README.md`
