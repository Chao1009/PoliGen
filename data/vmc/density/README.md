# VMC one-body point-nucleon densities (R. B. Wiringa, ANL)

Fetched 2026-09-03 for `include/lipolgen/cluster_config.hpp`'s alpha core
(`AlphaCoreSource::VmcHe4Density`) and for its independent 6Li validation
target (design section 3.2,
`docs/open_items/run_2026-09-02/design_G_cluster_config.md`).  **Both files
are committed** (10 kB total): `ClusterConfigSampler` opens `he4.density` at
run time and `tests/test_cluster_config.cpp` reads `li6.density`, so the
build cannot depend on a directory that may or may not be there.  They are
found through `data_dir()` -- `$LIPOLGEN_DATA_DIR` if set, else the
compiled-in `${CMAKE_SOURCE_DIR}/data` -- and live inside the already
installed `data/vmc` tree, so packaging needs nothing.

Fetched exactly as `data/vmc/README.md` documents: the ANL site
(`https://www.phy.anl.gov/theory/research/density/`) is live but sits behind
Cloudflare bot mitigation that 403s every non-interactive client from this
environment, so the Internet Archive Wayback Machine's raw (`id_`) form is
used, which serves the original file bytes with no HTML wrapper:

    curl -o data/vmc/density/he4.density \
      "http://web.archive.org/web/20150905165459id_/http://www.phy.anl.gov/theory/research/density/he4.density"
    curl -o data/vmc/density/li6.density \
      "http://web.archive.org/web/20150905180021id_/http://www.phy.anl.gov/theory/research/density/li6.density"

## Files fetched, one row per file

| file | content | source page (live URL) | Wayback snapshot fetched | bytes | md5 |
|---|---|---|---|---|---|
| `he4.density` | **4He(0+) -- AV18+UX -- 18-Dec-12**, VMC 500k samples; columns `R RHORP DRHORP`, R = 0.05 .. 20.05 fm step 0.1 | `.../theory/research/density/he4.density` | `web.archive.org/web/20150905165459id_/...` | 4751 | `5d717cf0ed642dc7779aef3873e06317` |
| `li6.density` | **6Li(1+) -- AV18+UX -- 20-Mar-14**, VMC 200k samples; same columns and grid | `.../theory/research/density/li6.density` | `web.archive.org/web/20150905180021id_/...` | 5189 | `96c6b15f045785db0593c00e9a32f9b1` |

Fetch date 2026-09-03 (this session); the *crawl* date is 5-Sep-2015 for both,
as the 14-digit timestamps show.  Verified plain text after download, and both
md5s verified against the design's.

## Printed normalization and rms, and what LiPolGen recomputes

| file | file's own `4*PI*TOTINT(RHORP*R**2:R)` | recomputed | file's own rms | recomputed rms |
|---|---|---|---|---|
| `he4.density` | 1.9974 (= Z) | 1.99738 | 1.4404 = `SQRT(4*PI*TOTINT(RHORP*R**4:R)/2)` | 1.44131 |
| `li6.density` | 2.9991 (= Z) | 2.99909 | **2.4433** = `SQRT(.../3)` | 2.44331 |

The recomputed columns are `np.trapz` on the file's own printed grid; the file's
own numbers come from the Fortran code's finer internal quadrature, which is
why the rms differs in the 4th digit (0.07 %).  `LI6_R_POINT_VMC_FM = 2.4433`
(cluster_config.hpp) is the file's printed value, not the re-integration.

## Format

A title line, a samples line, a legend, a `Normalization` line, an `rms radius`
line, then a column header and a `*****  ********  ********` rule, then rows
`R RHORP DRHORP` (fm, fm^-3, fm^-3).  The `****` rule is exactly the block
introduction `read_anl_momentum` (`cluster.cpp`) already keys on, and each row
carries 3 numbers, so that reader lands `col[0]` = RHORP and `err[0]` = DRHORP
with no new parser -- see design section 3.3.

**These are point-nucleon densities** (ANL rho_1, no nucleon form factor folded
in) -- exactly what a dipole model wants.  For 4He, N = Z = 2, so the neutron
density equals the proton density by isospin and the per-nucleon probability
density is `RHORP / Z`.

## Not fetched

- `he4_v18.density` (AV18 without UX, snapshot `20170427210207`) -- the optional
  Hamiltonian systematic of design section 3.2.  Not needed by any gate.
