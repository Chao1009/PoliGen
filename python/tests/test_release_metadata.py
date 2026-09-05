# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase E item E2: release metadata (CITATION.cff, AUTHORS, the license
section of README.md) matches the tree instead of drifting from it silently.

No test here asserts a real author name: the placeholder
`<AUTHOR NAME — to be filled by the author>` is the tracked, correct state
until a human fills it in (`docs/open_items/run_2026-09-03/STATUS.md`'s
author-decision table), and every test below is written so it keeps passing
once that happens (it checks the placeholder is either still there OR has
been replaced everywhere consistently, never that it is any specific string).
"""
import os
import re
import subprocess
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))

PLACEHOLDER = "<AUTHOR NAME — to be filled by the author>"

CMAKE_VERSION_RE = re.compile(
    r"project\(LiPolGen\s+VERSION\s+(?P<value>[0-9]+\.[0-9]+\.[0-9]+)")
CFF_VERSION_RE = re.compile(r'^version:\s*"?([0-9]+\.[0-9]+\.[0-9]+)"?\s*$', re.M)


def _read(*parts):
    with open(os.path.join(_ROOT, *parts), encoding="utf-8") as f:
        return f.read()


def _cmake_version():
    m = CMAKE_VERSION_RE.search(_read("CMakeLists.txt"))
    assert m, "CMakeLists.txt's project() version regex did not match"
    return m.group("value")


def test_citation_cff_exists_and_is_gplv3():
    text = _read("CITATION.cff")
    assert "cff-version:" in text
    assert "license: GPL-3.0-or-later" in text
    assert "title: LiPolGen" in text


def test_citation_cff_version_matches_cmakelists():
    """Single source of truth: CMakeLists.txt's `project(... VERSION ...)`,
    exactly the way `pyproject.toml`'s own
    `[tool.scikit-build.metadata.version]` regex already reads it -- CFF has
    no build-time templating of its own, so this test is what catches a
    version bump that updated one and not the other."""
    cmake_version = _cmake_version()
    m = CFF_VERSION_RE.search(_read("CITATION.cff"))
    assert m, "CITATION.cff has no `version:` line"
    assert m.group(1) == cmake_version, (
        f"CITATION.cff version {m.group(1)!r} != CMakeLists.txt "
        f"project() version {cmake_version!r}")


def test_pyproject_and_cff_read_the_same_cmake_version():
    """Belt and suspenders: pyproject.toml's own regex (not re-implemented
    here -- read verbatim) against the same CMakeLists.txt, so a change to
    either the CMake version string's format or its value is caught from
    both directions at once."""
    pyproject = _read("pyproject.toml")
    m = re.search(
        r"regex\s*=\s*'([^']*project\\\(LiPolGen[^']*)'", pyproject)
    assert m, "pyproject.toml's scikit-build version regex was not found"
    pattern = re.compile(m.group(1))
    hit = pattern.search(_read("CMakeLists.txt"))
    assert hit and hit.group("value") == _cmake_version()


def test_citation_cff_has_an_author_entry():
    text = _read("CITATION.cff")
    assert "authors:" in text
    assert PLACEHOLDER in text or "name:" in text


def test_authors_file_exists_and_names_the_same_placeholder_as_citation_cff():
    authors = _read("AUTHORS")
    cff = _read("CITATION.cff")
    if PLACEHOLDER in cff:
        assert PLACEHOLDER in authors, (
            "CITATION.cff still carries the placeholder but AUTHORS does not "
            "-- the two must be filled in together, never one at a time")
    else:
        assert PLACEHOLDER not in authors, (
            "CITATION.cff's placeholder was replaced but AUTHORS still "
            "carries it")


def test_readme_license_section_points_at_the_release_files():
    readme = _read("README.md")
    assert "`LICENSE`" in readme
    assert "AUTHORS" in readme
    assert "CITATION.cff" in readme


def test_license_file_is_gplv3_or_later():
    text = _read("LICENSE")
    assert "GNU GENERAL PUBLIC LICENSE" in text
    assert "Version 3" in text


def test_citation_cff_repository_matches_the_configured_git_remote():
    """`git remote -v` is the source for `repository-code`; skipped (not
    failed) when this checkout has no remote configured at all -- a tarball
    export, or a clone with `origin` renamed -- since there is then nothing
    to check against. A FORK with a different origin is expected to update
    CITATION.cff to match its own remote, which is exactly what this test
    would then catch."""
    try:
        out = subprocess.run(
            ["git", "remote", "get-url", "origin"], cwd=_ROOT,
            capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        pytest.skip("git is not available")
    if out.returncode != 0 or not out.stdout.strip():
        pytest.skip("no 'origin' remote configured in this checkout")
    remote = out.stdout.strip()
    m = re.match(r"git@github\.com:(.+?)(?:\.git)?$", remote) or \
        re.match(r"https://github\.com/(.+?)(?:\.git)?$", remote)
    if not m:
        pytest.skip(f"origin remote {remote!r} is not a plain github.com URL")
    expected = f"https://github.com/{m.group(1)}"
    cff = _read("CITATION.cff")
    rm = re.search(r'repository-code:\s*"?([^"\n]+)"?', cff)
    assert rm, "CITATION.cff has no repository-code line"
    assert rm.group(1).rstrip("/") == expected, (
        rm.group(1), expected)
