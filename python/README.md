# LiPolGen from Python

`lipolgen` is a thin pybind11 skin over the C++ core: the same objects, the
same names, the same conventions as `docs/USAGE.md`, with generation returning
**columnar numpy arrays** and exporters for the schemas
`PolarizedLithiumSim/evgen/polligen` consumes.

---

## 1. Install

There is **no `pip install` step**. The extension is a CMake target and the
package lives next to it in the build tree; `env.sh` already puts that on
`PYTHONPATH`:

```bash
cd LiPolGen
source env.sh                       # PYTHONPATH += build/python
cmake -S . -B build -DLIPOLGEN_WITH_PYTHON=ON
cmake --build build -j8             # builds _lipolgen and stages the package
python3 -c "import lipolgen; print(lipolgen.__version__)"
```

The build writes

```
build/python/lipolgen/_lipolgen.cpython-*.so   the extension
build/python/lipolgen/*.py                     __init__, export, cli
build/lipolgen-run                             the command-line entry point
```

`cmake --build build` re-stages the `.py` files on every build, so editing
`python/lipolgen/*.py` and rebuilding is enough — nothing is installed
anywhere.

Without `source env.sh`, set it by hand:

```bash
export PYTHONPATH=$PWD/build/python:$PYTHONPATH
export LD_LIBRARY_PATH=$PWD/../deps/install/lib:$LD_LIBRARY_PATH
```

Requirements: numpy (any 1.2x/2.x), pybind11 ≥ 2.10 at build time, and for the
test suite `pytest`, `pyhepmc`, optionally `PyYAML`.

Which optional tiers a build carries is visible at runtime:

```python
>>> import lipolgen as lg
>>> lg._lipolgen.HAVE_LHAPDF, lg._lipolgen.HAVE_HEPMC3, lg._lipolgen.HAVE_PYTHIA8
(True, True, True)
```

### Tests

```bash
source env.sh
python3 -m pytest python/tests -q          # 114 cases
python3 -m pytest python/tests/test_throughput.py -s     # prints ev/s
```

---

## 2. The shape of a run — same as the C++

A run is **one `PipelineConfig` + one `RunPlan`** (docs/USAGE.md §1).

```python
import lipolgen as lg

cfg = lg.PipelineConfig()
cfg.channel     = lg.PipelineChannel.TaggedLi6Alpha
cfg.isotope     = "6Li"
cfg.beam_config = 1                      # 10 GeV e x 99.5 GeV/u
cfg.n_events    = 100_000                # or cfg.lumi_pb = 2.0
cfg.seed        = 1
cfg.optics_choice = lg.OpticsChoice.YellowReportHighAcceptance

p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))

print(p.size(), p.sigma_pb(), p.counts())
ev = p.generate()                        # columnar dict of numpy arrays
```

`Pipeline` resolves everything in its constructor and then maps **event index
→ event** as a pure function, so `p.event(i)` is thread-safe and
`generate(..., nthreads=8)` is bit-identical to `nthreads=1`.

The shortcut, which is what `lipolgen-run` uses:

```python
ev = lg.run(isotope="6Li", channel="tagged-alpha", plan="tensor-thirds",
            events=100_000, seed=1)          # ev["pipeline"] is the Pipeline
```

### Generation entry points

| call | returns |
|---|---|
| `p.generate(n=0, events=False, nthreads=1)` | columnar dict; `n=0` is the whole run |
| `p.generate_lumi(lumi_pb, poisson=True, ...)` | the same run rebuilt at an integrated luminosity, all of it |
| `p.generate_range(first, last)` | `list[Event]` |
| `p.event(i)` | one `Event` |
| `p.write_hepmc(path, n=0)` | streams to HepMC3, stores nothing |

`events=True` adds the `Event` records under `ev["events"]` (needed by the HFS
exporter and by anything that walks particles).

### The columns

```
kinematics   x q2 y phi w2 nu s
spin         m (= m_ion) m_ion m_struck j pz pzz pe theta_s phi_s lam_e
electron     kp (n,4)  e_prime  theta_e  eta_e
tagging      k cos_theta_k phi_k alpha_s pt_s
spectator    pT theta p_lab R xL kx ky kz phi_spec route
coherent     t x_pom
bookkeeping  number weight category category_index category_names meta
```

`meta` is a dict carrying the beams, the seed, `sigma_pb` and `sigma_gen_mb`
(the polligen key), the optics name and the frame convention.

Every array **owns its buffer** (the C++ vector is moved into the array), so
it outlives the pipeline and costs no copy. Accessors like
`sampler.x_cells` return copies instead.

---

## 3. Inclusive, tagged, coherent

```python
# inclusive (docs/USAGE.md §2)
ev = lg.run(channel="inclusive", isotope="6Li", events=200_000, seed=1)
print(ev["x"].mean(), ev["q2"].mean(), ev["meta"]["sigma_pb"])

# tagged (§3) -- the Roman-Pot mask polligen applies
from lipolgen import export
ev = lg.run(channel="tagged-alpha", isotope="6Li", events=400_000, seed=1)
print("tag fraction", export.rp_accepted(ev).mean())     # 0.0248 at YR optics

# coherent (§4)
ev = lg.run(channel="coherent", isotope="6Li", events=300_000, seed=1,
            optics="tagging",
            coherent={"f0": 0.04, "slope_b": 50.0, "amp": 0.01})
print("<|t|> =", ev["t"].mean(), " ~ 1/B =", 1 / 50.0)
```

### Re-routing one sample at several optics (§5)

The route is recomputed from the record, never stored on it:

```python
p = ev["pipeline"]
p_u = p.beam_config.ion_momentum_per_nucleon
tag = lg.optics_for(lg.OpticsChoice.Tagging, "6Li", p_u)
seen = [lg.rp_tagged(e, tag, p.pot_config) for e in p.generate_range(0, 10_000)]
```

### Conservation checks (§7)

```python
e = p.event(0)
lg.momentum_residual(e), lg.momentum_scale(e), lg.charge_residual(e)
```

---

## 4. The pieces below the pipeline

Everything the C++ headers expose is bound, so a script can work at any level.

```python
# beams
cfg = lg.default_configs("6Li")[1]           # e(10) x 6Li(99.5/u)
cfg.sqrt_s_per_nucleon(), lg.li6().mass_per_nucleon()

# spin (rho, moments, max-entropy fills)
rho = lg.rho_from_populations(1.0, [0.5, 0.3, 0.2], theta=0.9, phi=0.4)
lg.vector_polarization(rho, 1.0), lg.tensor_polarization(rho, 1.0)
lg.populations_maxent(1.0, 8/13)             # (P_z, P_zz) = (8/13, 4/13)

# structure functions: a toy backend is always there, tables/LHAPDF optional
f2 = lg.ToyF2(); g1 = lg.ToyG1(f2)
b1 = lg.Li6B1(lg.MillerB1())                 # the rank-2 transfer to 6Li
if lg._lipolgen.HAVE_LHAPDF:
    f2 = lg.LhapdfSF("CT18NLO")
    epps = lg.Epps21Ratio()                  # EPPS21 nuclear ratio, per flavour

# the master formula
opt = lg.InclusiveKernel.Options()
opt.f2_source = f2
opt.g1_model  = g1
opt.b1_func   = b1.b1_func()
opt.delta_func = lambda x, q2, f1: lg.toy_delta_gluon(x, q2, f1, 1e-2)
opt.target_mass = False
kern = lg.InclusiveKernel(lg.li6(), opt)

st = lg.EventSpinState(lam_e=+1, pe=0.7, j=1.0, m=+1.0, theta_s=0.0)
t = kern.tables(0.1, 10.0, with_g2=True)
a = kern.amplitudes(t, 0.1, 10.0, cfg.s_per_nucleon(), st, True)
print(a.w_avg, a.a1, a.a2)
print(kern.dsigma(0.1, 10.0, 0.0, cfg.s_per_nucleon(), st))

# the sampler on its own -- polligen's `sample_category` dict, incl. `cell`
s = lg.InclusiveSampler(kern, cfg, lg.generator_scenario(lg.Scenario()))
cat = lg.tensor_thirds_plan(0.7, 0.6).categories[0]
batch = s.sample_n(cat, 100_000, seed=1, nthreads=8)
w = s.weights_for(batch, list(lg.tensor_thirds_plan(0.7, 0.6).categories))

# analysis-side estimators
lg.estimators.azz_thirds(n_plus, n_minus, n_zero, pzz=0.6)
lg.estimators.cos2phi_fit(phi_prime, pzz=0.6, nbins=36)
```

### T2: PYTHIA 8 hadronization

```python
beams  = lg.default_configs("6Li")[1]
bridge = lg.PythiaBridge(beams, lg.PythiaBridgeOptions())
lg.set_pythia_hadronizer(cfg_pipeline, bridge)   # a pure C++ hook, no GIL
p   = lg.Pipeline(cfg_pipeline, plan)
ev  = p.generate(0, events=True, nthreads=1)     # PythiaBridge is not re-entrant
print(bridge.stats.n_ok, bridge.stats.n_failed)
```

A hadronizer assigned as a **Python callable** (`cfg.hadronizer = f`) also
works, but pybind11 re-acquires the GIL for every event; `set_pythia_hadronizer`
binds the C++ method directly and stays GIL-free.

---

## 5. Exporting to polligen

`lipolgen.export` speaks the three schemas
`PolarizedLithiumSim/evgen/polligen` reads. Each function takes either the
columnar dict or a sequence of `Event` records.

```python
from lipolgen import export

ev = lg.run(channel="tagged-alpha", events=100_000, seed=1)

export.inclusive_dict(ev)      # polligen.sample -> x q2 y phi m category lam_e
export.tagged_dict(ev)         # polligen.tagged -> ... k cos_theta_k phi_k
                               #                    pT theta p_lab R xL kx ky kz
                               #                    phi_spec route
export.rp_accepted(ev).mean()  # polligen.tagged.rp_accepted
export.write_columns_npz(ev, "run.npz")
```

`polligen.hfs.HFSSample` (needs the T2 tier — there are no hadrons without it):

```python
export.write_hfs_npz(ev["events"], "hfs.npz", meta=dict(ev["meta"]))

# read back, with polligen
import sys; sys.path.insert(0, ".../PolarizedLithiumSim/evgen")
from polligen.hfs import HFSSample, hadronic_sums
s = HFSSample.load("hfs.npz")
sigma, ptx, pty = hadronic_sums(s.p4, s.offsets)
```

The file carries `offsets, pid, charge, p4 (E,px,py,pz), x, q2, y, kp, weight,
e_energy, p_per_nucleon` and a JSON `meta` with `sigma_gen_mb` — the key
`HFSSample.concatenate` reweights on. The **scattered electron is not in the
particle list** (it is `kp`), and spectator fragments are excluded by default,
both following `tools/pythia8/gen_dis_hfs.py`.

### Defaults that are a choice

Three physics defaults moved on 2026-08-29 and are reachable (and pinned by
`test_module.py`) rather than buried:

```python
opt = lg.InclusiveKernel.Options()
opt.target_mass          # True  -- the exact finite-gamma vector kernel
opt.g2_scale             # 1.0   -- multiplies g2; the twist-3 systematic
                         #          is re-run at 0.0 and 1.5
lg.EMC_BASELINE_DEFAULT  # EmcBaseline.Epps21 -- the unpolarized baseline the
                         # polarized-EMC transfer is referenced to
```

A hadronizer on the **coherent** channel is refused by `PipelineConfig.validate()`:
`PythiaBridge` v0 has no coherent-diffractive target, so it would invent a
nucleon that is not in the record's balance and the event would lose
four-momentum and charge. Set `hadronize_coherent=True` (or
`--hadronize-coherent`) to reproduce that known-broken behaviour on purpose.

### The one polligen key with no pipeline equivalent

`cell`, the flat accepted-(x, Q²)-cell index that
`polligen.sample.weights_for` needs for Mode-W reweighting, is an
`InclusiveSampler` internal and is **not** on the `Event` record. It is
present in `InclusiveSampler.sample_n()`'s dict and absent from
`Pipeline.generate()`'s. Reweight sampler batches, not pipeline output.

---

## 6. `lipolgen-run`

```bash
./build/lipolgen-run --isotope 6Li --config 1 --channel tagged-alpha \
                     --plan tensor-thirds --events 100000 --seed 1 \
                     --hepmc out.hepmc --npz out.npz
```

`python3 -m lipolgen ...` is the same command.

```
--channel     inclusive | tagged-alpha | tagged-6Li-alpha | tagged-7Li-alpha
              | tagged-d-p | coherent
--plan        tensor-thirds (azz) | helicity-flip (apar)
              | transverse-tensor (cos2phi) | tensor-flip (flip)
--optics      yr-high-acceptance | yr-high-divergence | tagging | tagging-legacy
--events N | --lumi PB       exclusive
--seed --run --pz --pzz --pe --rel-lumi-offset
--cluster-beta --p-d --inclusive-b1
--coherent-f0 --coherent-slope-b --coherent-amp
--nthreads N                 forced to 1 with --hadronize
--hadronize                  run the T2 (PYTHIA 8) tier
--hadronize-coherent         allow --hadronize on the coherent channel, which
                             the core refuses by default (C4)
--npz FILE                   the columnar sample
--hfs-npz FILE               polligen HFSSample (needs --hadronize)
--hepmc FILE                 HepMC3 Asciiv3
--config-file FILE           JSON, or YAML if PyYAML is installed
```

The config file is a flat mapping of the same names with dashes turned into
underscores; command-line switches win over it, and an unknown key is an
error rather than a silent no-op.

```json
{"isotope": "7Li", "config": 1, "channel": "tagged-alpha",
 "plan": "helicity-flip", "events": 200000, "seed": 5,
 "optics": "tagging", "pz": 0.7, "pe": 0.7}
```

```bash
./build/lipolgen-run --config-file run.json --events 20000
```

It prints σ per spin category, the generated event count, the throughput, and
(for the tagged and coherent channels) the tag fraction split into the main
Roman-Pot window and the near-beam tail.

---

## 7. Throughput

Measured on this machine, single core, from Python, columnar output only
(`pytest python/tests/test_throughput.py -s`):

| channel | ev/s |
|---|---|
| inclusive | 5.9 × 10⁵ |
| tagged 6Li α | 4.9 × 10⁵ |
| coherent | 5.4 × 10⁵ |
| inclusive, `Event` records kept | 1.8 × 10⁵ |
| inclusive, T2 (PYTHIA 8) | 1.9 × 10⁴ |

The bare C++ streaming loop does 3.4 × 10⁶ ev/s inclusive; the gap is the
column fill itself — ~40 scattered writes into 40 separate arrays per event —
which is the price of a structure-of-arrays result. The GIL is released around
the whole loop, so `nthreads=8` scales (modestly: at ~1.5 µs per event the
block scheduling is a visible fraction).

HepMC3 output is the bottleneck when it is on, ~18 k ev/s, exactly as in C++.
