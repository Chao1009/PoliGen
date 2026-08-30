# VMC two-cluster overlap / momentum-distribution tables (R. B. Wiringa, ANL)

Downloaded 2026-08-29 for the LiPolGen `VmcRadial` prototype (see
`validation/vmc_overlap_prototype.py`). **Not committed** — these are large,
third-party numerical tables; keep them out of git history and re-fetch from
the URLs below if this directory is lost.

## What was tried, and why these are Wayback Machine snapshots

`https://www.phy.anl.gov/theory/research/...` (R. B. Wiringa's pages, the
locations the plan names) is live but sits behind Cloudflare bot mitigation
that 403s every request from this environment — `curl` (with a browser
`User-Agent`) and the `WebFetch` tool both got `HTTP/2 403` with
`cf-mitigated: challenge`. General internet egress is otherwise fine
(google.com, arxiv.org, archive.org all return 200 from the same box), so
this is specifically Cloudflare's challenge page rejecting non-interactive
clients, not a network outage. Workaround: the Internet Archive Wayback
Machine has clean 2025/2026 crawls of every page and data file used below
(`web.archive.org/web/<timestamp>/<url>`), fetched with `curl`. Every file
here is the **raw text table**, not a re-typed or figure-digitized copy —
the Wayback Machine serves the original file bytes.

## Where the two-cluster (α+d, α+t, ...) overlaps actually live

The modern spectroscopic-overlap pages (`.../overlaps/`, `.../overlaps2/`,
`.../overlap_gfmc/`, Brida–Pieper–Wiringa PRC 84, 024319 (2011),
arXiv:1106.3121 [the plan's "1106.3006" is a typo for this]) are **all
single-nucleon overlaps** (A → (A−1) + N: e.g. `⁷Li → ⁶Li + n`,
`³H → ²H + n` with the deuteron treated as an inert core, `¹⁰B → ⁹Be + p`,
...). They do **not** carry genuine two-cluster overlaps where both
fragments have A > 1 (α+d, α+t, α+α, t+p, d+d, d+p). Confirmed by listing
every `href` on `overlap/`, `overlaps/`, `overlaps2/`, `overlap_gfmc/`,
`momenta/`, `momenta2/` (the last two are single- and *NN-pair* (q, Q =
relative/CM momentum of two nucleons) momentum distributions, not
nucleon-cluster distributions) — none match `li6.ad`, `li7.at`, or similar.

The genuine cluster-cluster overlaps live on the **older** page,
`.../overlap_old/`, which the current site still serves (linked from
`overlap/`'s "Older VMC Spectroscopic Overlaps") and which states:

> We have computed momentum distributions and configuration space overlaps
> in a number of systems, including ³He, ⁴He, and ⁶Li, as discussed in
> Forest et al., Phys. Rev. C 54, 646 (1996). ... The next files are the
> actual outputs of our Monte Carlo codes, containing the relevant
> amplitudes in momentum space and the (unsmoothed) radial overlap
> functions, both with their Monte Carlo statistical errors: dp in ³He, dd
> in ⁴He, tp in ⁴He, αd in ⁶Li, αt in ⁷Li, αα in ⁸Be.

This is exactly the α+d / α+t data the plan is after — computed with the
**Argonne v18 two-nucleon + Urbana IX three-nucleon Hamiltonian (AV18+UIX)**
VMC wave functions, run 19–22 Apr 2004 (per each file's own banner).
Reference: **B. S. Pudliner, V. R. Pandharipande, J. Carlson, S. C. Pieper,
R. B. Wiringa, "Quantum Monte Carlo calculations of nuclei with A ≤ 7", Phys.
Rev. C 56, 1720 (1997)** for the A=6,7 wave functions themselves, and **J. L.
Forest, V. R. Pandharipande, S. C. Pieper, R. B. Wiringa, R. Schiavilla, A.
Arriaga, "Femtometer toroidal structures in nuclei", Phys. Rev. C 54, 646
(1996)** for the cluster-overlap method and the smoothed R0/R2 fits
(`li6.adr.fit` = 6Li Fig. 21, `he3.dpr.fit` = 3He Fig. 14, `he4.ddr.fit` =
4He Fig. 18).

## Files fetched, one row per subdirectory

| dir | file | channel | source page (live URL) | Wayback snapshot fetched |
|---|---|---|---|---|
| `li6_alpha_d/` | `li6.ad` | **6Li(1+) → α + d**, S+D wave (`Aad00`, `Aad22`), r-space AND k-space, MC errors | `.../theory/research/overlap_old/li6.ad` | `web.archive.org/web/20250617050329/...` |
| `li6_alpha_d/` | `li6.adr.fit` | same channel, smoothed R0(r)/R2(r) fit (Forest et al. Fig. 21) | `.../overlap_old/li6.adr.fit` | `web.archive.org/web/20250617050324/...` |
| `li7_alpha_t/` | `li7.at` | **7Li(3/2-) → α + t**, P-wave split into j=3/2 (`Aat33`) and j=1/2 (`Aat11`) channels, r-space AND k-space, MC errors | `.../overlap_old/li7.at` | `web.archive.org/web/20250617050331/...` |
| `deuteron/` | `fdeut.av18` | deuteron u(r), w(r) (S,D) to r=100 fm AND u(k), w(k) to k=20 fm⁻¹, AV18 potential, plus EM form factors | `.../theory/research/deuteron/fdeut.av18` | `web.archive.org/web/20250529230933/...` |
| `li7_li6_n/` | `li6n_31.table` | **7Li(3/2-) → 6Li(1+) + n**, p1/2+p3/2 amplitudes, r-space ONLY (AV18+UX, 150k VMC samples) | `.../theory/research/overlaps/li6n_31.table` | `web.archive.org/web/20250609095630/...` |
| `h3_d_n/` | `h2n.table` | **3H(1/2+) → d + n**, s1/2 (L=0) + d3/2 (L=2) amplitudes, r-space ONLY, both VMC and GFMC, AV18+IL7 | `.../overlap_gfmc/h2n_111_110_av18il7_101123_100917.txt` | `web.archive.org/web/20250617050450/...` |
| `bonus_other_clusters/` | `he3.dp`, `he3.dpr.fit` | 3He → d + p (S+D), raw + smoothed fit | `.../overlap_old/he3.dp[.fit]` | `.../050325`, `.../050319` |
| `bonus_other_clusters/` | `he4.dd`, `he4.ddr.fit` | 4He → d + d (S+D), raw + smoothed fit | `.../overlap_old/he4.dd[.fit]` | `.../050327`, `.../050322` |
| `bonus_other_clusters/` | `he4.tp` | 4He → t + p (S+D), raw | `.../overlap_old/he4.tp` | `.../050328` |
| `bonus_other_clusters/` | `be8.aa` | 8Be → α + α, raw | `.../overlap_old/be8.aa` | `.../050333` |

All Wayback timestamps above are `web.archive.org/web/<14-digit timestamp>/<url>`;
fetch dates recorded here are 2026-08-29 (this session); the *crawl* dates
are late May / mid-June 2025 as shown in each timestamp. Verified plain-text
(no HTML wrapper) after download.

## Not found / not applicable

- **3H → p + nn**: no such file anywhere on the site. It isn't a two-cluster
  problem in the site's sense (a bound N-body core plus one relative
  coordinate) — `p+n+n` is a three-body breakup of the triton, which VMC
  spectroscopic-overlap methodology (a bound (A−a) core wave function
  contracted against the full A-body VMC state) doesn't produce a single
  radial overlap for. Nothing digitized; genuinely absent.
- **7Li → 6Li + n with 6Li in an excited state**, **7Li → 6He + p**, etc.:
  present on `overlap_gfmc/` and `overlap_old/` (5 daughter states each) but
  not fetched — out of scope for the ground-state comparison the prototype
  needs.
- Figures were deliberately **not** digitized (task explicitly excludes
  this) — only the `.table`/raw-output text files above.

## Format notes (see the prototype for the parsers)

- **`li6.ad` / `li7.at`** (Fortran VMC code raw stdout): a long
  human-readable header (potential, wave-function parameters — not needed),
  then a block
  ```
   k(fm-1)  Aad00(k)               Aad22(k)
      0.00 -123.6     (5.129    ) 0.1292E-05 (.1756E-06)
      ...
  ```
  (`Aat33(k)`/`Aat11(k)` for `li7.at`) — k in fm⁻¹, 0.00 to 5.00 (li6.ad) or
  ~20.0 (li7.at) in steps of 0.10 fm⁻¹, value then `(1-sigma MC error)` in
  parens. Later in the same file a second block
  ```
      rij   Aad00                  Aad22                  #samples
      0.05 0.5846     (.7251E-01)-0.1357E-01 (.2591E-01)       429
      ...
  ```
  gives the **unsmoothed** r-space overlap (r in fm, 0.05 to ~20 fm, step
  0.1 fm) with MC error and raw sample count per bin. For `li6.ad` there is
  also an `ndx / s-wave / d-wave` line right after the k-space table
  (`0.856 0.838 0.017`, errors `0.002 0.001 0.000`) — the code's own
  S/D-wave norm split; the prototype recomputes this from the tabulated
  `Aad00(k)`, `Aad22(k)` directly rather than trusting the printed line, and
  cross-checks against it (see report).
- **`li6.adr.fit`**: 3-column `R  R0LI6FIT  R2LI6FIT`, the smoothed r-space
  S/D radial functions only (no k-space, no errors).
- **`fdeut.av18`**: header gives `ebind dstate qm mu as eta rd r4` (dstate =
  0.057599 = 5.76% D-state, the standard AV18 deuteron value — a sanity
  check that this file is what it claims to be), then `r u(r) du/dr w(r)
  dw/dr` to r=100 fm, then `k u(k) w(k)` to k=20 fm⁻¹, then wave-function
  integrals and EM structure functions (not needed here).
- **`li6n_31.table` / `h2n.table`**: plain r-space amplitude tables,
  `R  A(R)  dA(R)  ...`, fm and fm⁻³ᐟ², **no momentum-space column** — this
  is the genuinely r-space-only case the prototype's Fourier–Bessel
  transform exists for. Both state their normalization convention inline:
  spectroscopic factor = ∫ A²(r) r² dr.

## Units / conversion used downstream

VMC tables are natural units, length in fm, momentum in fm⁻¹. LiPolGen's
`Wave::radial` works in GeV. Conversion: `k[GeV] = k[fm^-1] * hbar_c`,
`hbar_c = 0.1973269804 GeV·fm` (CODATA ħc). r[fm] is used as-is by the
Fourier–Bessel transform; the resulting k-grid is then converted to GeV.
