# Run 2026-09-02 — status

Cron job id: a07c52e5 (session-only, hourly at :11, first fire 2026-09-02 18:11 CDT, auto-expires after 7 days)
Started: 2026-09-02 17:05 CDT (planning). Work is deferred to the first cron fire
after the usage-limit reset.

| phase | state | commit | notes |
|---|---|---|---|
| A physics-channels doc | DONE | 6f67d84 | docs/PHYSICS_CHANNELS.md 654 lines / 111 rows / 48 refs; gate 826 refs 0 broken (validation/check_physics_channels_links.py, --fix mode); README linked; refresh line refs + add C/D/G rows in Phase H |
| B VMC tag fractions + writer fix | DONE | 8c41899 | tag fractions at 5x41/10x100/18x275 (USAGE), beta band retired, generated-mass rule pinned (doctest + pytest); 276 doctest / 146 pytest green |
| C tensor RC band | IN PROGRESS | | design pass wf_38270d00-8c1 -> docs/open_items/run_2026-09-02/design_C_tensor_rc.md; implementation not started |
| D b1(6Li) convolution | IN PROGRESS | | design pass wf_38270d00-8c1 -> design_D_b1_li6.md; implementation not started |
| E spin-3/2 note | DONE (note) | 472a7d4 | docs/theory/SPIN32_FINITE_GAMMA.md (1195 lines, 2 review passes); follow-up after A commits: g1_rank3 code comments in spin.hpp/xsec.hpp and the generate_full.cpp pzz_true fix |
| F packaging | DONE | 611e40a | pyproject + scikit-build-core; fresh-venv editable install 66 s, wheel 1.31 MB, in-tree flow unchanged; README snippet to apply after A |
| G coherent 6Li on-ramp | TODO | | |
| H close-out | TODO | | |

States: TODO / IN PROGRESS (with what exists on disk) / DONE (with commit hash).
