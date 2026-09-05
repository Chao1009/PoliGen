#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Gate: every first-party source file carries an SPDX license header.

LiPolGen is GPL-3.0-or-later (`LICENSE`, `CITATION.cff`) but until 2026-09-05
no source file said so itself -- checked once, with a plain
`grep -rl SPDX-License-Identifier include src tests python validation`, which
returned nothing.  This script is the machine-checkable form of that same
check, run by `python/tests/test_spdx_headers.py`: every file under the
covered set (`GROUPS` below) must carry `SPDX-License-Identifier:
GPL-3.0-or-later` in its own first 5 lines (line 1 for a plain source file,
line 2 for one that opens with a `#!` shebang).

No copyright line is added anywhere -- the author's name is a separate,
still-open decision (`docs/open_items/run_2026-09-03/STATUS.md`'s
author-decision table); SPDX identifies the LICENSE, not the holder, so it is
uncontroversial to add on its own.  A generated or vendored file is not in
`GROUPS` and is not this script's concern; the set is deliberately the same
five directories (plus `python/bindings.cpp`) named in
`docs/open_items/run_2026-09-03/STATUS.md`'s Phase E item E2.

Run from the repository root; exit status 1 lists every file still missing
the header.  `--write` inserts it into every missing file (idempotent: a file
that already has it is left untouched) -- this is what did the one-time
insertion across all 100 covered files on 2026-09-05; a NEW file joining the
covered set afterward is expected to carry its own header from the day it is
written, exactly like any other project convention (an `#include` guard, a
copyright line elsewhere) that a linter checks but does not usually need to
author.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SPDX_ID = "GPL-3.0-or-later"
SPDX = f"SPDX-License-Identifier: {SPDX_ID}"

# (glob relative to ROOT) -- every source directory named in STATUS.md's E2,
# plus the one file (python/bindings.cpp) named explicitly there because it
# sits beside pure-Python siblings that take the other comment token.
GROUPS = [
    "include/**/*.hpp",
    "src/**/*.hpp",
    "src/**/*.cpp",
    "tests/**/*.hpp",
    "tests/**/*.cpp",
    "python/bindings.cpp",
    "python/lipolgen/*.py",
    "validation/*.py",
]

HEAD_WINDOW = 5  # lines scanned for the header: enough for a shebang + SPDX
                 # + a blank line, never so wide that a mid-file mention of
                 # the string (a docstring quoting this very rule) would count


def target_files() -> list:
    """Every file the covered globs name, in a stable (sorted, deduplicated)
    order.  `python/bindings.cpp` is also reachable via no other glob here,
    so there is nothing to deduplicate against it, but `**` glob results can
    overlap across two patterns in general, hence the `seen` guard."""
    seen: set = set()
    files: list = []
    for pattern in GROUPS:
        for p in sorted(ROOT.glob(pattern)):
            if p.is_file() and p not in seen:
                seen.add(p)
                files.append(p)
    return files


def has_spdx(path: Path) -> bool:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return any(SPDX in ln for ln in lines[:HEAD_WINDOW])


def missing_spdx() -> list:
    return sorted(str(p.relative_to(ROOT)) for p in target_files()
                 if not has_spdx(p))


def header_line(path: Path) -> str:
    """`//` for C++, `#` for Python -- decided by extension, not by
    directory, so `src/pythia/lhaup_dis.hpp` and `python/lipolgen/cli.py`
    each get the token their own language comments with."""
    token = "#" if path.suffix == ".py" else "//"
    return f"{token} {SPDX}\n"


def insert_header(path: Path) -> None:
    """Idempotent: a file that already has the header in its scan window is
    left untouched.  A `#!` shebang stays line 1 and the header becomes line
    2; every other file gets the header as its literal first line."""
    if has_spdx(path):
        return
    text = path.read_text(encoding="utf-8")
    line = header_line(path)
    if text.startswith("#!"):
        nl = text.find("\n") + 1
        path.write_text(text[:nl] + line + text[nl:], encoding="utf-8")
    else:
        path.write_text(line + text, encoding="utf-8")


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    write = "--write" in argv
    files = target_files()
    if write:
        for p in files:
            insert_header(p)
    missing = [str(p.relative_to(ROOT)) for p in files if not has_spdx(p)]
    missing.sort()
    if missing:
        print(f"{len(missing)} of {len(files)} covered file(s) missing "
              f"{SPDX!r} in their first {HEAD_WINDOW} lines:")
        for m in missing:
            print("  ", m)
        return 1
    print(f"{len(files)} files checked, all carry {SPDX!r}"
          + (" (written)" if write else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
