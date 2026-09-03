# Run 2026-09-02 — status

Cron job id: a07c52e5 (session-only, hourly at :11, first fire 2026-09-02 18:11 CDT, auto-expires after 7 days)
Started: 2026-09-02 17:05 CDT (planning). Work is deferred to the first cron fire
after the usage-limit reset.

| phase | state | commit | notes |
|---|---|---|---|
| A physics-channels doc | DONE | 6f67d84 | docs/PHYSICS_CHANNELS.md 654 lines / 111 rows / 48 refs; gate 826 refs 0 broken (validation/check_physics_channels_links.py, --fix mode); README linked; refresh line refs + add C/D/G rows in Phase H |
| B VMC tag fractions + writer fix | DONE | 8c41899 | tag fractions at 5x41/10x100/18x275 (USAGE), beta band retired, generated-mass rule pinned (doctest + pytest); 276 doctest / 146 pytest green |
| C tensor RC band | IN PROGRESS | | design docs/open_items/run_2026-09-02/design_C_tensor_rc.md (reviewed, revised); implementation workflow wf_e1bbc7bf-65d running: rc.hpp/rc.cpp, xsec/sampler refactor, wiring, tests, docs |
| D b1(6Li) convolution | IN PROGRESS | | design design_D_b1_li6.md (reviewed, revised); implementation script staged (scratchpad wf/phase_d.js), launches after C (shared build tree, same wiring files) |
| E spin-3/2 note | DONE | 472a7d4 + 0961c60 | docs/theory/SPIN32_FINITE_GAMMA.md (1195 lines, 2 review passes); g1_rank3/R_3 code comments; generate_full pure-fill pzz_true fixed |
| F packaging | DONE | 611e40a | pyproject + scikit-build-core; fresh-venv editable install 66 s, wheel 1.31 MB, in-tree flow unchanged; README snippet to apply after A |
| G coherent 6Li on-ramp | IN PROGRESS | 6e393aa (ii) | (ii) eSTARlight baseline DONE: estarlight_li6.md (6Li J/psi 1.77 nb, B 39-55 GeV^-2, slope_b band now citable, f0 stays scenario); (i) alpha+d configuration sampler: design pass wf_d557adab-23c running, implementation after D; (iii) OPEN_ITEMS section 11 write-up with (i) |
| H close-out | TODO | | |

States: TODO / IN PROGRESS (with what exists on disk) / DONE (with commit hash).
