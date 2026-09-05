# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase E item E2: every first-party source file carries an SPDX header.

`validation/check_spdx_headers.py` is the gate; this file is what wires it
into the suite (mirroring how `test_doc_link_gate.py` wires in
`check_physics_channels_links.py`), loaded by path since it is a standalone
script with no package to import.
"""
import importlib.util
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
_SCRIPT = os.path.join(_ROOT, "validation", "check_spdx_headers.py")


def _load():
    spec = importlib.util.spec_from_file_location("_spdx_gate", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_no_covered_file_is_missing_its_spdx_header():
    """The real tree, under the real gate: exit 0, nothing missing."""
    mod = _load()
    assert mod.main([]) == 0, mod.missing_spdx()


def test_covered_set_is_not_accidentally_empty():
    """A glob that stopped matching would make the gate pass by checking
    nothing -- the same failure mode `test_doc_link_gate.py`'s reference-count
    test guards against for the citation gate."""
    mod = _load()
    files = mod.target_files()
    assert len(files) > 90, files
    by_dir = {"include": 0, "src": 0, "tests": 0, "python": 0, "validation": 0}
    for p in files:
        rel = str(p.relative_to(mod.ROOT))
        for d in by_dir:
            if rel.startswith(d + "/") or rel.startswith(d + os.sep):
                by_dir[d] += 1
                break
    assert all(n > 0 for n in by_dir.values()), by_dir


def test_header_is_the_right_token_for_the_language(tmp_path):
    mod = _load()
    cpp = tmp_path / "x.hpp"
    cpp.write_text("#pragma once\n")
    py = tmp_path / "x.py"
    py.write_text("import os\n")
    mod.insert_header(cpp)
    mod.insert_header(py)
    assert cpp.read_text().splitlines()[0] == "// SPDX-License-Identifier: GPL-3.0-or-later"
    assert py.read_text().splitlines()[0] == "# SPDX-License-Identifier: GPL-3.0-or-later"


def test_header_follows_a_shebang_rather_than_replacing_it(tmp_path):
    mod = _load()
    f = tmp_path / "x.py"
    f.write_text("#!/usr/bin/env python3\n\"\"\"doc\"\"\"\n")
    mod.insert_header(f)
    lines = f.read_text().splitlines()
    assert lines[0] == "#!/usr/bin/env python3"
    assert lines[1] == "# SPDX-License-Identifier: GPL-3.0-or-later"
    assert lines[2] == '"""doc"""'


def test_insert_header_is_idempotent(tmp_path):
    mod = _load()
    f = tmp_path / "x.hpp"
    f.write_text("#pragma once\n")
    mod.insert_header(f)
    once = f.read_text()
    mod.insert_header(f)
    assert f.read_text() == once


def test_missing_spdx_finds_a_planted_gap(tmp_path, monkeypatch):
    mod = _load()
    monkeypatch.setattr(mod, "ROOT", tmp_path)
    (tmp_path / "include").mkdir()
    (tmp_path / "include" / "a.hpp").write_text("// SPDX-License-Identifier: GPL-3.0-or-later\nint a;\n")
    (tmp_path / "include" / "b.hpp").write_text("int b;\n")
    missing = mod.missing_spdx()
    assert missing == ["include/b.hpp"], missing
