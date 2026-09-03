#!/usr/bin/env python3
"""Gate for docs/PHYSICS_CHANNELS.md: every `path:line` reference must point at a
line (within +-2) that contains the symbol named next to it.

A reference is recognised as `` `path:line` `` optionally followed by `` `name` ``
(possibly after a few characters) on the same table cell / line.  Exit status 1
lists every broken reference.  `--fix` rewrites a broken reference to the
nearest line that declares the named symbol (unambiguous nearest only) and
reports what it changed.  Run from the repository root.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "PHYSICS_CHANNELS.md"
REF = re.compile(r"`((?:include|src|python|data|validation|docs|tests|examples)/[^`:\s]+):(\d+)`(?:[^`\n]{0,12}`([^`\n]+)`)?")
WIN = 2


def resolve(lines: list[str], base: str, old: int) -> int | None:
    """Line (1-based) of the candidate declaration of `base` nearest to `old`,
    or None when there is no candidate or the nearest is ambiguous."""
    pat = re.compile(r"(?<![A-Za-z0-9_])" + re.escape(base) + r"(?![A-Za-z0-9_])")
    cands = [i + 1 for i, ln in enumerate(lines) if pat.search(ln)]
    if not cands:
        return None
    cands.sort(key=lambda i: abs(i - old))
    if len(cands) > 1 and abs(cands[0] - old) == abs(cands[1] - old):
        return None
    return cands[0]


def main() -> int:
    fix = "--fix" in sys.argv[1:]
    if not DOC.exists():
        print(f"missing {DOC}")
        return 1
    text = DOC.read_text()
    cache: dict[str, list[str]] = {}
    bad, n, fixed = [], 0, []

    def repl(m: re.Match) -> str:
        nonlocal n
        path, line, name = m.group(1), int(m.group(2)), m.group(3)
        n += 1
        f = ROOT / path
        if not f.exists():
            bad.append(f"{path}:{line}  (file missing)")
            return m.group(0)
        lines = cache.setdefault(path, f.read_text(errors="replace").splitlines())
        base = re.split(r"[(<]", name)[0].split("::")[-1].strip() if name else None
        ok_range = 1 <= line <= len(lines)
        if base and ok_range:
            lo, hi = max(1, line - WIN), min(len(lines), line + WIN)
            if any(base in lines[i - 1] for i in range(lo, hi + 1)):
                return m.group(0)
        elif ok_range:
            return m.group(0)
        if fix and base:
            new = resolve(lines, base, line)
            if new is not None:
                fixed.append(f"{path}:{line} -> :{new}  `{name}`")
                return m.group(0).replace(f"{path}:{line}`", f"{path}:{new}`", 1)
        bad.append(f"{path}:{line}  " + ("(only %d lines)" % len(lines) if not ok_range
                                        else f"`{name}` not within +-{WIN} lines"))
        return m.group(0)

    new_text = REF.sub(repl, text)
    if fix and fixed:
        DOC.write_text(new_text)
    print(f"{n} references checked, {len(bad)} broken" + (f", {len(fixed)} fixed" if fix else ""))
    for x in fixed:
        print("  fixed", x)
    for b in bad:
        print("  ", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
