// SPDX-License-Identifier: GPL-3.0-or-later
// D7 (2026-10-27): the structure-function setting of the C++ core suite.
//
// Since D7 a default-constructed `PipelineConfig` names the SHIPPED grid
// backends -- `UNPOL_SF_DEFAULT` = ct18nlo (`LhapdfSF("CT18NLO", 0)`) and
// `POL_SF_DEFAULT` = nnpdfpol (`LhapdfG1("NNPDFpol11_100", 0)`) -- which the
// core library cannot build: they belong to the optional LHAPDF tier, which
// attaches them through its filler when a `Pipeline` is built
// (include/lipolgen/pipeline.hpp, "D7").  This suite installs no filler
// outside the LHAPDF-guarded D7 cases of tests/test_lhapdf.cpp, so every
// MACHINERY, closure and port-gate test of a `Pipeline` runs on the toy
// structure functions NAMED: `toy_config()` (pipeline.hpp), which is bit for
// bit the pre-D7 default and therefore keeps every number these tests pin.
//
// THE RULE, mechanical and greppable: `PipelineConfig cfg;` became
// `PipelineConfig cfg = toy_config();` wherever the test builds a Pipeline
// (or validates a config it then builds), with this header included in
// place of "lipolgen/pipeline.hpp", which it includes.  A test whose subject
// is the D7 wiring itself -- the default refusing without the tier, running
// with it, the caller-kernel rule, the struck-slot rule -- does NOT use it.
#ifndef LIPOLGEN_TESTS_SF_TEST_SUPPORT_HPP
#define LIPOLGEN_TESTS_SF_TEST_SUPPORT_HPP

#include "lipolgen/pipeline.hpp"

#endif  // LIPOLGEN_TESTS_SF_TEST_SUPPORT_HPP
