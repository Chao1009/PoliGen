# SPDX-License-Identifier: GPL-3.0-or-later
"""A13 (2026-11-03): the VMC alpha-d port gate, from Python.

`validation/reference/tagged.json`'s `channels.li6_alpha_vmc` is polligen's
`li6_alpha_channel(wave="vmc")` (a REDUCED block: quantum numbers, the ANL
table's P_D, the two stretched-state dilutions), spliced into the committed
file byte for byte from the dump of 2026-11-06
(docs/open_items/run_2026-11-03/phase_P5_port_gate.md).  The bound module's
`TaggedModel(li6_alpha_channel(source=VmcAV18))` must meet it: the dilutions
at the DECLARED relative 1e-8 (tests/test_tagged.cpp `kVmcPortRtol`, written
before the comparison; measured 5.389e-9 vector / 3.196e-9 tensor), P_D and
the quantum numbers at rtol 1e-12.  The block covers the alpha-d motion only
-- not the embedded deuteron, not 7Li -- and the label is not compared.
"""

import json
import os

import pytest

import lipolgen as lg

VMC_PORT_RTOL = 1e-8   # A13 (2026-11-03): declared, never widened
RTOL = 1e-12


def _vmc_present():
    return os.path.isfile(os.path.join(lg._lipolgen.data_dir(),
                                       "vmc", "momenta", "li6_ad1.momentum"))


@pytest.fixture(scope="module")
def block(reference_dir):
    path = os.path.join(reference_dir, "tagged.json")
    if not os.path.isfile(path):
        pytest.skip("validation/reference/tagged.json absent")
    with open(path) as fh:
        doc = json.load(fh)
    assert "li6_alpha_vmc" in doc["channels"], "A13's block is missing"
    return doc["channels"]["li6_alpha_vmc"]


@pytest.fixture(scope="module")
def vmc_model():
    if not _vmc_present():
        pytest.skip("the ANL VMC tables are not installed")
    ch = lg.li6_alpha_channel(source=lg.ClusterWaveSource.VmcAV18)
    return ch, lg.TaggedModel(ch)


def test_the_block_declares_the_tolerance(block):
    assert block["port_gate_rtol"] == VMC_PORT_RTOL
    assert "alpha-d motion ONLY" in block["scope"]


def test_quantum_numbers_and_p_d(block, vmc_model):
    ch, _ = vmc_model
    for key in ("j_ion", "s_struck", "s_spec", "s_channel"):
        assert getattr(ch, key) == pytest.approx(block[key], rel=RTOL, abs=0.0)
    assert lg.VMC_P_D_LI6 == pytest.approx(block["P_D"], rel=RTOL, abs=0.0)


@pytest.mark.parametrize("name", ["vector_dilution", "tensor_dilution"])
def test_dilutions_at_the_declared_tolerance(block, vmc_model, name):
    _, model = vmc_model
    got = getattr(model, name)()
    want = block["model"][name]
    assert abs(got - want) <= VMC_PORT_RTOL * abs(want), (
        "%s: LiPolGen %.17g, polligen %.17g, rel %.3e > %g"
        % (name, got, want, abs(got / want - 1.0), VMC_PORT_RTOL))
