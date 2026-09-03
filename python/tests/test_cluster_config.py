"""Phase G(i): the alpha+d configuration sampler through the bindings and CLI.

The C++ doctest (`tests/test_cluster_config.cpp`) owns the physics identities;
this file owns the binding surface -- array shapes and dtypes, every option
field, every analytic predictor, the constants -- and the `lipolgen-configs`
command line, including the two sidecar fields only Python can fill (`md5`,
`git`).
"""

import json
import math
import shutil
import subprocess
import sys

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import configs


@pytest.fixture(scope="module")
def sampler():
    return lg.ClusterConfigSampler()


# ------------------------------------------------------------- the constants

def test_constants_are_importable():
    assert lg.CONFIG_STREAM == 0x434F4E464947
    assert lg.LI6_R_POINT_VMC_FM == pytest.approx(2.4433)
    assert lg.LI6_QUADRUPOLE_GFMC_FM2 == pytest.approx(-0.20)
    assert lg.LI6_QUADRUPOLE_GFMC_ERR_FM2 == pytest.approx(0.06)
    assert lg.LI6_R2_POINT_FM2 == pytest.approx(6.0788)
    # r_point measured (2.4655) is LARGER than li6.density's VMC 2.4433
    assert math.sqrt(lg.LI6_R2_POINT_FM2) > lg.LI6_R_POINT_VMC_FM
    for p in (lg.VMC_HE4_DENSITY, lg.VMC_LI6_DENSITY, lg.VMC_LI6_AD_FIT,
              lg.VMC_DEUTERON_WAVE):
        assert p.startswith("vmc/")


def test_enums():
    assert lg.AlphaCoreSource.VmcHe4Density is not lg.AlphaCoreSource.Gaussian
    assert set(lg.AlphaDSource.__members__) == {"FitRescaled", "FitRaw",
                                               "OverlapRaw"}
    assert set(lg.SndConfigFormat.__members__) == {"He3Compatible", "Annotated"}


# ------------------------------------------------------------- the options

def test_every_option_field_is_bound():
    o = lg.ClusterConfigOptions()
    fields = ["alpha_source", "alpha_d_source", "theta_s", "phi_s",
              "alpha_cm_inflate", "min_nn_separation_fm", "alpha_d_scale",
              "quadrupole_target_fm2", "exact_coherence", "n_r", "n_c",
              "r_max_fm", "rnp_max_fm"]
    for f in fields:
        assert hasattr(o, f), f
    assert len(fields) == 13
    assert o.theta_s == pytest.approx(math.pi / 2)   # transverse by default
    assert o.n_r == 512 and o.n_c == 96
    assert o.alpha_cm_inflate is True
    assert o.exact_coherence is False
    o.validate()


def test_exact_coherence_is_reserved_and_throws():
    o = lg.ClusterConfigOptions()
    o.exact_coherence = True
    with pytest.raises(RuntimeError):
        o.validate()
    with pytest.raises(RuntimeError):
        lg.ClusterConfigSampler(o)


def test_quadrupole_target_out_of_range_throws():
    o = lg.ClusterConfigOptions()
    o.quadrupole_target_fm2 = 0.5           # above the s = 0 floor, +Q_d
    with pytest.raises(RuntimeError):
        lg.ClusterConfigSampler(o)


# ------------------------------------------------------- the analytic layer

def test_analytic_predictors_reproduce_the_design(sampler):
    assert sampler.p_d_alpha_d() == pytest.approx(0.020112, abs=2e-6)
    assert sampler.s_alpha_d() == pytest.approx(0.854232, abs=1e-5)
    assert sampler.tensor_dilution() == pytest.approx(
        1.0 - 0.9 * sampler.p_d_alpha_d(), rel=1e-12)
    assert sampler.r2_alpha_fm2() == pytest.approx(2.0774, abs=1e-4)
    assert sampler.r2_np_fm2() == pytest.approx(15.4817, abs=1e-3)
    assert sampler.r2_alpha_d_fm2() == pytest.approx(16.963, abs=1e-3)
    assert sampler.r2_analytic_fm2() == pytest.approx(6.4447, abs=1e-3)
    assert sampler.q_matter_analytic_fm2(1) == pytest.approx(-1.2309, abs=1e-3)
    assert sampler.q_matter_analytic_fm2(0) == pytest.approx(
        -2.0 * sampler.q_matter_analytic_fm2(1), rel=1e-14)
    assert sampler.delta_perp_analytic_fm2(1) == pytest.approx(-0.1026,
                                                               abs=1e-3)
    assert sampler.eps_b0_equivalent() == pytest.approx(-0.0506, abs=1e-3)
    assert sampler.a2_from_geometry(0.3, 1) == pytest.approx(0.1976, abs=1e-3)
    assert sampler.asymptotic_ds_ratio() == pytest.approx(-0.048, abs=0.01)


def test_the_quadrupole_band_is_mandatory_and_computed(sampler):
    band = sampler.quadrupole_band_fm2()
    assert len(band) == 3
    assert band[0] == pytest.approx(-0.0818)              # measured, TUNL
    assert band[1] == pytest.approx(lg.LI6_QUADRUPOLE_GFMC_FM2)
    assert band[2] == pytest.approx(0.5 * sampler.q_matter_analytic_fm2(1))
    # the model overshoots the measurement by a factor ~7.5
    assert 6.0 < abs(band[2] / band[0]) < 9.0
    # and it MOVES with the dial rather than being a literal
    o = lg.ClusterConfigOptions()
    o.quadrupole_target_fm2 = -0.0818
    dialed = lg.ClusterConfigSampler(o)
    assert dialed.quadrupole_band_fm2()[2] == pytest.approx(-0.0818, abs=1e-9)
    assert dialed.quadrupole_dial_s() == pytest.approx(0.4032, abs=1e-3)
    assert sampler.quadrupole_dial_s() == 1.0


def test_a2_from_quadrupole_reproduces_the_published_deuteron():
    q_d = 0.26967026805838806
    rows = {0.05: (-0.04, +0.08), 0.10: (-0.08, +0.15),
            0.20: (-0.17, +0.30), 0.30: (-0.28, +0.43)}
    for t, (m1, m0) in rows.items():
        p1 = lg.a2_from_quadrupole(2 * q_d, 2, t, 1)
        p0 = lg.a2_from_quadrupole(2 * q_d, 2, t, 0)
        assert p1 == pytest.approx(m1, rel=0.10)
        assert p0 == pytest.approx(m0, rel=0.25)
        assert p0 == pytest.approx(-2 * p1, rel=1e-14)


def test_cluster_density_is_the_published_closed_form():
    u, v = 0.37, -0.11
    for c in (-1.0, -0.6, 0.0, 0.35, 1.0):
        c2 = c * c
        g1 = (4 * u * u - 2 * math.sqrt(2) * (1 - 3 * c2) * u * v
              + (5 - 3 * c2) * v * v) / (16 * math.pi)
        g0 = (2 * u * u + 2 * math.sqrt(2) * (1 - 3 * c2) * u * v
              + (1 + 3 * c2) * v * v) / (8 * math.pi)
        assert lg.cluster_density(u, v, 1, c) == pytest.approx(g1, rel=1e-14)
        assert lg.cluster_density(u, v, 0, c) == pytest.approx(g0, rel=1e-14)
        # the m_S-summed density is the sum of the three branches
        s = sum(lg.cluster_amp2(u, v, 1, ms, c) for ms in (1, 0, -1))
        assert s == pytest.approx(lg.cluster_density(u, v, 1, c), rel=1e-15)


def test_radial_moments_on_the_deuteron_tables(sampler):
    m = lg.radial_moments(sampler.np_grid, sampler.np_wave(0),
                          sampler.np_wave(2))
    assert m.norm() == pytest.approx(1.0, rel=5e-5)
    assert m.p_d() == pytest.approx(0.057599, rel=5e-5)
    assert m.quadrupole() / 4.0 == pytest.approx(0.269673, rel=5e-5)
    assert math.sqrt(m.r2 / m.norm()) / 2.0 == pytest.approx(1.967364, rel=5e-5)


def test_the_grid_quadrature_tracks_the_analytic_layer(sampler):
    # They differ only by the trapezoid's own O(h^2) on the TABLES' 0.1 fm
    # nodes, which no n_r refinement removes.
    assert sampler.r2_grid_fm2() == pytest.approx(sampler.r2_analytic_fm2(),
                                                  rel=2e-3)
    assert sampler.q_matter_grid_fm2(1) == pytest.approx(
        sampler.q_matter_analytic_fm2(1), rel=5e-3)
    # the m_S conditioning of R-hat: exact -2/7 and +1/7 for M = +1
    assert sampler.p2_alpha_d_grid(1, -1.0) == pytest.approx(-2 / 7, abs=1e-6)
    assert sampler.p2_alpha_d_grid(1, 0.0) == pytest.approx(1 / 7, abs=1e-3)
    pd = sampler.p_d_alpha_d()
    assert sampler.p_ms(1, 1.0) == pytest.approx(1 - 0.9 * pd, rel=1e-12)
    assert sampler.p_ms(1, 0.0) == pytest.approx(0.3 * pd, rel=1e-12)
    assert sampler.p_ms(1, -1.0) == pytest.approx(0.6 * pd, rel=1e-12)


def test_match_li6_radius(sampler):
    o = lg.ClusterConfigOptions()
    o.alpha_d_scale = sampler.match_li6_radius()
    matched = lg.ClusterConfigSampler(o)
    assert math.sqrt(matched.r2_analytic_fm2()) == pytest.approx(
        lg.LI6_R_POINT_VMC_FM, rel=1e-12)
    assert 0.9 < o.alpha_d_scale < 1.0


def test_longitudinal_axis_gives_no_transverse_asymmetry():
    o = lg.ClusterConfigOptions()
    o.theta_s = 0.0
    s = lg.ClusterConfigSampler(o)
    assert s.delta_perp_analytic_fm2(1) == 0.0
    assert s.eps_b0_equivalent() == 0.0
    assert s.q_matter_analytic_fm2(1) < 0.0     # the geometry is unchanged


# ------------------------------------------------------------- the sampling

def test_sample_set_array_shapes_and_dtypes(sampler):
    s = sampler.sample_set(500, 1, 4242, 0)
    assert len(s) == 500
    assert s.positions.shape == (500, 6, 3)
    assert s.positions.dtype == np.float64
    assert s.isospin.shape == (500, 6)
    assert s.isospin.dtype == np.int64
    # 3 protons, 3 neutrons; nucleons 0-3 are the alpha core, 4 p and 5 n
    assert (s.isospin.sum(axis=1) == 0).all()
    assert (s.isospin[:, 4] == 1).all()
    assert (s.isospin[:, 5] == -1).all()
    # c.m. at the origin for EVERY configuration
    assert np.abs(s.positions.sum(axis=1)).max() < 1e-12
    assert s.cm_max < 1e-12
    assert s.m_ion == 1
    assert (s.m_ion_per_config == 1).all()
    assert set(np.unique(s.m_s)) <= {-1.0, 0.0, 1.0}
    assert s.seed == 4242 and s.run == 0
    assert "alpha-d R_0,R_2" in s.provenance


def test_unpolarized_splits_into_near_equal_thirds(sampler):
    s = sampler.sample_set_unpolarized(100000, 1, 0)
    counts = np.bincount(s.m_ion_per_config + 1, minlength=3)
    assert sorted(counts.tolist()) == [33333, 33333, 33334]
    assert s.m_ion == -2
    # the quadrupole of the equal-thirds mixture vanishes, within 5 sigma
    gate = 5.0 * s.q_matter_sd_fm2 / math.sqrt(len(s))
    assert abs(s.q_matter_fm2) < gate


def test_monte_carlo_closure_at_five_sigma(sampler):
    n = 100000
    s = sampler.sample_set(n, 1, 20260902, 0)
    assert s.r2_mean_fm2 == pytest.approx(
        sampler.r2_analytic_fm2(), abs=5.0 * s.r2_sd_fm2 / math.sqrt(n))
    assert s.q_matter_fm2 == pytest.approx(
        sampler.q_matter_analytic_fm2(1),
        abs=5.0 * s.q_matter_sd_fm2 / math.sqrt(n))
    assert s.delta_perp_fm2 == pytest.approx(
        sampler.delta_perp_analytic_fm2(1),
        abs=5.0 * s.delta_perp_sd_fm2 / math.sqrt(n))
    # the per-configuration spreads themselves
    assert s.r2_sd_fm2 == pytest.approx(3.66, rel=0.05)
    assert s.q_matter_sd_fm2 == pytest.approx(28.9, rel=0.05)
    assert s.delta_perp_sd_fm2 == pytest.approx(2.87, rel=0.05)


def test_counter_based_streams_are_reproducible(sampler):
    a = sampler.sample_set(64, 1, 777, 3)
    b = sampler.sample_set(128, 1, 777, 3)
    assert np.array_equal(a.positions, b.positions[:64])
    c = sampler.sample_set(64, 1, 777, 4)
    assert not np.array_equal(a.positions[0], c.positions[0])
    # config i on its own, from the header's stream slot
    one = sampler.sample(lg.Rng(777, 3, lg.CONFIG_STREAM, 5), 1)
    assert np.array_equal(one.positions, a.positions[5])


# ----------------------------------------------------------------- the CLI

def test_cli_defaults_and_switches():
    o = configs.resolve(["--n", "77", "--m", "-1"])
    assert o["n"] == 77
    assert o["m"] == "-1"
    assert o["seed"] == configs.DEFAULTS["seed"]
    assert o["theta_s"] == pytest.approx(math.pi / 2)
    assert o["alpha_d_source"] == "fit-rescaled"
    opts = configs.make_options(o)
    assert opts.alpha_d_source == lg.AlphaDSource.FitRescaled
    assert opts.alpha_source == lg.AlphaCoreSource.VmcHe4Density


def test_cli_writes_a_table_and_a_filled_sidecar(tmp_path):
    out = tmp_path / "li6_m1.dat"
    mom = tmp_path / "li6_m1.moments.json"
    rc = configs.main(["--n", "200", "--m", "+1", "--seed", "1",
                       "--out", str(out), "--moments-json", str(mom),
                       "--quiet"])
    assert rc == 0
    rows = out.read_text().splitlines()
    assert len(rows) == 200
    for line in rows:
        assert "#" not in line
        fields = line.split()
        assert len(fields) == 25          # 18 coordinates + 6 isospins + m
        assert int(fields[24]) == 1
    meta = json.loads((tmp_path / "li6_m1.dat.meta.json").read_text())
    # the 13 option fields, all recorded
    for f in ("alpha_source", "alpha_d_source", "theta_s", "phi_s",
              "alpha_cm_inflate", "min_nn_separation_fm", "alpha_d_scale",
              "quadrupole_target_fm2", "exact_coherence", "n_r", "n_c",
              "r_max_fm", "rnp_max_fm"):
        assert f in meta["options"], f
    assert meta["options"]["quadrupole_dial_s"] == 1.0
    assert meta["n_config"] == 200
    assert meta["rng"]["stream_value"] == lg.CONFIG_STREAM
    # the two fields only Python can fill
    assert len(meta["inputs"]) == 4
    for entry in meta["inputs"]:
        assert entry["md5"] is not None and len(entry["md5"]) == 32
    assert meta["git"] is None or len(meta["git"]) == 40
    assert meta["quadrupole_band_fm2"]["measured"] == pytest.approx(-0.0818)
    assert meta["caveats"]
    m = json.loads(mom.read_text())
    s = lg.ClusterConfigSampler()
    assert m["analytic"]["q_matter_fm2"]["1"] == pytest.approx(
        s.q_matter_analytic_fm2(1), rel=1e-12)
    assert m["analytic"]["eps_b0_equivalent"] == pytest.approx(
        s.eps_b0_equivalent(), rel=1e-12)
    assert m["sampled"]["n"] == 200
    assert m["timing"]["us_per_config"] > 0.0


def test_cli_unpolarized_thirds(tmp_path):
    out = tmp_path / "unpol.dat"
    configs.main(["--n", "1000", "--m", "unpolarized", "--seed", "2",
                  "--out", str(out), "--quiet"])
    ms = [int(line.split()[24]) for line in out.read_text().splitlines()]
    counts = [ms.count(v) for v in (1, 0, -1)]
    assert max(counts) - min(counts) <= 1
    meta = json.loads((tmp_path / "unpol.dat.meta.json").read_text())
    assert meta["m_ion"] == "unpolarized"


def test_cli_annotated_format_has_comments(tmp_path):
    out = tmp_path / "ann.dat"
    configs.main(["--n", "10", "--format", "annotated", "--out", str(out),
                  "--quiet"])
    lines = out.read_text().splitlines()
    assert any(line.startswith("#") for line in lines)
    assert len([l for l in lines if not l.startswith("#")]) == 10


def test_cli_warns_about_a_longitudinal_axis(capsys):
    configs.main(["--n", "10", "--theta-s", "0"])
    assert "longitudinal" in capsys.readouterr().out


def test_cli_match_li6_radius(tmp_path, capsys):
    out = tmp_path / "matched.dat"
    configs.main(["--n", "50", "--match-li6-radius", "--out", str(out)])
    assert "alpha_d_scale" in capsys.readouterr().out
    meta = json.loads((tmp_path / "matched.dat.meta.json").read_text())
    assert meta["options"]["alpha_d_scale"] < 1.0
    assert meta["moments"]["r_rms_fm"] == pytest.approx(lg.LI6_R_POINT_VMC_FM,
                                                        rel=1e-9)


def test_module_entry_point_runs():
    out = subprocess.run([sys.executable, "-m", "lipolgen.configs", "--help"],
                         capture_output=True, text=True)
    assert out.returncode == 0
    assert "lipolgen-configs" in out.stdout
    for flag in ("--moments-json", "--quadrupole-target", "--alpha-d-source",
                 "--match-li6-radius", "--theta-s"):
        assert flag in out.stdout


# ------------------------------------------------------ the provenance fields

def test_git_head_rejects_an_unrelated_repository(tmp_path):
    """A venv inside somebody ELSE's checkout must not be stamped as ours.

    `git rev-parse HEAD` succeeds in any work tree, so an installed copy
    sitting inside an unrelated repository would otherwise record that
    repository's SHA as the generator's provenance.
    """
    if shutil.which("git") is None:
        pytest.skip("git is not installed")
    repo = tmp_path / "other"
    (repo / "sub").mkdir(parents=True)
    (repo / "f.txt").write_text("not lipolgen\n")
    run = lambda *a: subprocess.run(["git", *a], cwd=repo, check=True,
                                    capture_output=True)
    run("init", "-q")
    run("add", "f.txt")
    run("-c", "user.name=t", "-c", "user.email=t@example.invalid",
        "commit", "-qm", "x")
    # the unrelated repository DOES have a HEAD ...
    assert subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo / "sub",
                          capture_output=True, text=True).returncode == 0
    # ... and git_head() still refuses it
    assert configs.git_head(str(repo / "sub")) is None
    assert configs.git_head(str(repo)) is None
    # our own tree, when it is one, still reports a full SHA
    head = configs.git_head()
    assert head is None or len(head) == 40
