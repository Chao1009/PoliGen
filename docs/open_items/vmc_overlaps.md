## Summary

**Data found and provenance** (all in `/home/cpeng/Projects/polli/LiPolGen/data/vmc/`, uncommitted, with full URL/date/citation records in `data/vmc/README.md`):

`phy.anl.gov/theory/research/...` is live but sits behind Cloudflare bot-mitigation that 403s every request from this environment (confirmed via both `curl` and `WebFetch`, headers show `cf-mitigated: challenge`); general internet egress is otherwise fine. Fell back to the Internet Archive Wayback Machine, which has clean 2025/2026 crawls of every page and file used. All downloaded files are the **raw original text**, not digitized figures.

Key finding: the *current* ANL overlap pages (`overlaps/`, `overlaps2/`, `overlap_gfmc/` — the Brida–Pieper–Wiringa PRC 84,024319/arXiv:1106.3121 family) are **all single-nucleon overlaps** (A → (A−1)+N); they do **not** carry genuine two-cluster (α+d, α+t) overlaps. Those live on the **older** page `overlap_old/` (still linked from the site), which states it holds "the actual outputs of our Monte Carlo codes... dp in ³He, dd in ⁴He, tp in ⁴He, **αd in ⁶Li, αt in ⁷Li**, αα in ⁸Be" — AV18+UIX VMC, run 2004, per Forest et al., PRC 54, 646 (1996) and Pudliner et al., PRC 56, 1720 (1997). Downloaded:

- `li6_alpha_d/li6.ad` — **6Li→α+d**, S+D wave, **both r-space AND k-space**, with MC errors
- `li7_alpha_t/li7.at` — **7Li→α+t**, P-wave (j=3/2, j=1/2), **both r-space AND k-space**, with MC errors
- `deuteron/fdeut.av18` — deuteron u(r),w(r) to 100 fm and u(k),w(k) to 20 fm⁻¹ (AV18) — clean (non-MC) validation pair
- `li7_li6_n/li6n_31.table` (7Li→6Li+n) and `h3_d_n/h2n.table` (3H→d+n) — r-space only
- Bonus: `bonus_other_clusters/` (3He→d+p, 4He→d+d, 4He→t+p, 8Be→α+α)
- **Not found**: 3H→p+nn (not a two-cluster problem in this methodology; genuinely absent, not digitized from anywhere)

**Prototype**: `/home/cpeng/Projects/polli/LiPolGen/validation/vmc_overlap_prototype.py`, output plot `validation/vmc_overlap_comparison.png`. It parses these Fortran-output tables, implements a general Fourier–Bessel transform, and discovers (and documents) that the VMC cluster-overlap files use a **different normalization convention** (`A(k) = 4π∫A(r)j_L(kr)r²dr`) than the standard reduced-radial deuteron convention (`sqrt(2/π)∫`) — validated to ratio 1.000 against the clean deuteron pair, and to ~2-6% (MC noise level) against the independently-estimated VMC k-space columns.

**Headline physics result**: VMC k²n(k) has a genuine node/interference dip near k≈0.15–0.2 GeV that no two-parameter Hulthén form can reproduce, and its high-k tail is substantially **softer** than either β=0.30 or β=0.40 Hulthén, especially for 7Li α+t (e.g. P(k>0.45 GeV): VMC 0.0095 vs Hulthén β=0.40's 0.22 — over 20× harder in the current model). VMC 6Li α+d D-state probability = **0.020**, vs the Python scenario placeholder `P_D_LI6=0.0867` — about 4.3× too high.

**Interface change needed in `cluster.hpp`**: add a `VmcRadial` tabulated-ψ_L(k) class and a nullable member on `Wave` so `radial(k, kappa)` dispatches to either the existing analytic switch or interpolation into the table, with no signature change for callers:
```cpp
class VmcRadial {
 public:
  VmcRadial(std::vector<double> k_gev, std::vector<double> psi_k, int l);
  double operator()(double k) const;   // interpolated, unnormalized — 0 outside range
};
struct Wave {
  int l = 0; double prob = 1.0; double beta = BETA_DEFAULT;
  std::shared_ptr<const VmcRadial> vmc;   // null => current analytic Hulthen
  double radial(double k, double kappa) const;  // dispatches on vmc
};
```
`kappa` stays unused on the VMC branch (the table already encodes the physical falloff) but stays in the signature for compatibility. `VmcRadial` need not self-normalize — `TaggedModel` already normalizes any `Wave::radial` on its own grid, so this fits the existing "unnormalized radial forms" contract exactly.

Nothing was committed; `data/vmc/`, `validation/vmc_overlap_prototype.py`, and `validation/vmc_overlap_comparison.png` are untracked (`git status` confirms only `??` entries), and no C++ source was touched.