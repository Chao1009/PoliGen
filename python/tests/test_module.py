"""The module imports, and the bound surface is the one docs/USAGE.md uses."""

import numpy as np
import pytest

import lipolgen as lg


def test_import_and_version():
    assert lg.__version__ == "0.1.0"
    assert lg._lipolgen.HAVE_HEPMC3 in (True, False)


def test_constants_have_one_definition():
    # constants.hpp is mirrored read-only; nothing here recomputes a number.
    assert lg.TENSOR_LL_SIGN == 1.0
    assert lg.ALPHA_EM == pytest.approx(1.0 / 137.036, rel=0, abs=0)
    assert lg.M_NUCLEON == 0.9383
    assert lg.GEV2_TO_PB == 0.3894e9


@pytest.mark.parametrize("name", [
    # beams / spin / sf / xsec / bookkeeping / sampler / pipeline / io
    "Ion", "BeamConfig", "default_configs", "nucleus_mass",
    "rho_from_populations", "moments_along_axis", "populations_maxent",
    "ToyF2", "ToyG1", "MillerB1", "Li6B1", "NuclearF2", "toy_b1",
    "InclusiveKernel", "EventSpinState", "Amplitudes", "SFTables",
    "azz", "a_cos2phi", "a_parallel", "err_azz",
    "SpinCategory", "RunPlan", "tensor_thirds_plan", "helicity_flip_plan",
    "transverse_tensor_plan", "tensor_flip_plan",
    "Scenario", "generator_scenario", "InclusiveSampler",
    "PipelineConfig", "Pipeline", "PipelineChannel", "OpticsChoice",
    "Event", "Particle", "Vec4", "Rng", "Optics", "Route",
    "route_of", "rp_tagged", "momentum_residual", "charge_residual",
])
def test_bound_names(name):
    assert hasattr(lg, name), name


def test_optional_tiers_follow_the_build():
    if lg._lipolgen.HAVE_HEPMC3:
        assert hasattr(lg, "HepMC3Writer")
    if lg._lipolgen.HAVE_PYTHIA8:
        assert hasattr(lg, "PythiaBridge")
        assert hasattr(lg, "hfs_summary")
    if lg._lipolgen.HAVE_LHAPDF:
        assert hasattr(lg, "LhapdfSF")
        assert hasattr(lg, "Epps21Ratio")


def test_beams_match_the_documented_menu():
    cfgs = lg.default_configs("6Li")
    assert [c.ion_momentum_per_nucleon for c in cfgs] == [40.8, 99.5, 137.5]
    cfgs7 = lg.default_configs("7Li")
    assert [c.ion_momentum_per_nucleon for c in cfgs7] == [40.8, 99.5, 117.9]


def test_spin_maxent_rational_anchors():
    # (P_z, P_zz) = (8/13, 4/13) at J=1 and (0.7, 0.4) at J=3/2 -- the t=3
    # spin-temperature anchors of the validation matrix.
    p1 = lg.populations_maxent(1.0, 8.0 / 13.0)
    m1 = lg.moments_along_axis(1.0, p1)
    assert m1.vector == pytest.approx(8.0 / 13.0, rel=1e-12)
    assert m1.tensor == pytest.approx(4.0 / 13.0, rel=1e-10)
    p32 = lg.populations_maxent(1.5, 0.7)
    m32 = lg.moments_along_axis(1.5, p32)
    assert m32.vector == pytest.approx(0.7, rel=1e-12)
    assert m32.tensor == pytest.approx(0.4, rel=1e-10)


def test_arrays_are_numpy_and_own_their_buffers():
    r = lg.Rng(1, 0, 0, 0)
    a = r.uniforms(10)
    assert isinstance(a, np.ndarray) and a.dtype == np.float64
    assert a.base is not None            # the capsule holding the vector
    assert ((a > 0) & (a < 1)).all()


def test_vec4_algebra():
    a = lg.Vec4(5.0, 1.0, 2.0, 3.0)
    b = lg.Vec4(1.0, 0.0, 0.0, 1.0)
    c = a - b
    assert (c.e, c.px, c.py, c.pz) == (4.0, 1.0, 2.0, 2.0)
    assert a.pt() == pytest.approx(np.hypot(1.0, 2.0))
    assert np.allclose(a.as_array(), [5.0, 1.0, 2.0, 3.0])
