"""The docs/PHYSICS_CHANNELS.md citation gate, exercised by the suite.

Two things are pinned here.

1. The document itself passes the STRICT gate, and every citation shape in it
   is actually reached by a rule.  Until 2026-09-04 the gate accepted any
   textual occurrence within +-2 lines of the referenced line, so an anchor
   could slide onto the doc comment above its declaration, onto a blank line,
   or onto the NEXT symbol's declaration and still be reported as "0 broken";
   160 anchors were re-anchored when that was fixed.  Later the same day three
   further blind spots were closed, each of which had let a real citation go
   stale: `path:first-last` RANGES matched no rule at all (93 of them, 20
   stale), citations with NO name were checked only for file existence and line
   range (one had slid onto a blank line), and a citation adjacent to another
   was swallowed as the first one's "name" and never checked in its own right
   (ten of them).  This test is what stops the fourth repeat.  A fourth and a
   fifth were closed on 2026-09-04: a use site in a file that declares its name
   nowhere was pinned to EVERY line of that file naming it as code (14 lines
   for one citation, 29 citations accepting more than one), and 5 point
   citations plus 1 range naming PYTHIA's own sources matched no rule at all.

2. The rules themselves, on fixtures.  A comment that merely names the symbol
   is not a declaration, a blank line is not, a neighbour's declaration is not,
   and a BASE CLASS is not -- while a use site in a file that does not declare
   the symbol still is a legitimate citation.  An unnamed citation must land on
   a line with something on it.  A range must be in bounds, must not have a
   blank edge, and must still hash to its recorded fingerprint -- a block that
   MOVED is named and fixable, a block that was EDITED is reported as needing a
   human.  And an allow-list exemption covers ONE citation and holds only while
   the line it names still carries its pin -- including the one whose name its
   own file never declares, where the escape was the carve-out beside the
   exemption rather than the exemption itself.  An ambiguous use site must
   carry a fingerprint of its line, a lone one need not, and an external
   PYTHIA citation is checked against the dependency tree or skipped BY NAME.
   And -- since 2026-09-06, Phase D item D1 -- `--record-ranges` REFUSES a
   block the citing sentence does not evidence: a wrong-block range is not
   fingerprinted at all, a right-block one is, and an exemption needs a pin
   that is inside the block and nowhere else in the file.

The gate is a standalone script with no package to import, so it is loaded by
path; nothing here needs the compiled extension.
"""
import importlib.util
import json
import os
import subprocess

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
_GATE = os.path.join(_ROOT, "validation", "check_physics_channels_links.py")


def _load():
    spec = importlib.util.spec_from_file_location("_doc_link_gate", _GATE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def gate():
    return _load()


def _capture(fn, *a):
    import contextlib
    import io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = fn(*a)
    return rc, buf.getvalue()


def _run(gate, tmp_path, doc_text, files, argv=(), ranges=None):
    """Point the gate at a throw-away tree and run it, returning (rc, stdout)."""
    for rel, text in files.items():
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    doc = tmp_path / "docs" / "PHYSICS_CHANNELS.md"
    doc.parent.mkdir(parents=True, exist_ok=True)
    doc.write_text(doc_text)
    gate.ROOT, gate.DOC = tmp_path, doc
    if ranges is not None:
        f = tmp_path / gate.RANGES_REL
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(json.dumps(ranges))
    rc, out = _capture(gate.main, list(argv))
    return rc, out, doc.read_text()


HDR = """\
namespace lipolgen {
/// Where the alpha core's one-body point-nucleon density comes from.
enum class AlphaCoreSource {
  Vmc = 0,
};

/// The coherent scenario knobs.
struct CoherentScenario {
  double f0 = 0.04;
  double x_coh = 0.01;
};
}  // namespace lipolgen
"""
# line numbers in HDR: 2 comment, 3 enum, 5 `};`, 6 blank, 7 comment, 8 struct,
# 9 f0, 10 x_coh

SRC = """\
#include "hdr.hpp"
void use() {
  // the AlphaCoreSource branch is chosen here
  AlphaCoreSource s = AlphaCoreSource::Vmc;
  (void)s;
}
"""
# line numbers in SRC: 3 comment mentioning it, 4 a genuine use site


# ----------------------------------------------------------------- the real document

def test_real_document_passes_strict():
    """The shipped docs/PHYSICS_CHANNELS.md, under the real gate."""
    gate = _load()
    assert gate.main([]) == 0


def test_real_document_reference_count():
    """Every citation is checked -- a regex that stopped matching would make
    the gate pass by checking nothing.  Both counts are asserted: the point
    citations and the ranges, which are a separate regex and a separate rule."""
    gate = _load()
    rc, out = _capture(gate.main, [])
    first = out.split()[0]
    assert int(first) > 1000, out
    n_ranges = int(out.split(",")[1].split()[0])
    assert n_ranges > 90, out
    assert rc == 0


def test_extra_documents_are_covered_and_clean():
    """Phase E item E1: `docs/PHYSICS_CHANNELS.md` is not the only document
    with more than a handful of `path:line` citations.  `docs/theory/
    SPIN32_FINITE_GAMMA.md` (134: 19 named points, 115 ranges) and
    `docs/PYTHIA_BRIDGE.md` (8 external PYTHIA upstream pointers) are both
    covered by the same gate run, each labelled by name so a broken citation
    in either is visible, not swallowed into the primary document's tally.
    `docs/USAGE.md`, `docs/CONVENTIONS.md` and `docs/T2_CHAIN.md` carry none;
    `docs/OPEN_ITEMS_SOLUTIONS.md` carries 3 (not "more than a handful") --
    all four stay uncovered, on purpose; see this file's module docstring and
    `docs/open_items/run_2026-09-03/phase_E_numbers.md` §E1 for the counts."""
    gate = _load()
    names = {p.rsplit("/", 1)[-1] for p, _ in gate.EXTRA_DOCS}
    assert names == {"SPIN32_FINITE_GAMMA.md", "PYTHIA_BRIDGE.md"}, gate.EXTRA_DOCS
    rc, out = _capture(gate.main, [])
    assert rc == 0, out
    assert "[SPIN32_FINITE_GAMMA.md]" in out, out
    assert "[PYTHIA_BRIDGE.md]" in out, out
    for label, min_n, min_ranges, min_ext in (
        ("SPIN32_FINITE_GAMMA.md", 15, 90, 0),
        ("PYTHIA_BRIDGE.md", 0, 0, 5),
    ):
        line = next(ln for ln in out.splitlines() if ln.endswith(f"[{label}]"))
        assert ", 0 broken" in line, line
        n = int(line.split()[0])
        n_ranges = int(line.split(",")[1].split()[0])
        n_ext = int(line.split(",")[2].split()[0])
        assert n >= min_n, line
        assert n_ranges >= min_ranges, line
        assert n_ext >= min_ext, line


def test_extra_document_missing_under_root_is_skipped_not_broken(gate, tmp_path):
    """A test fixture tree (or any `ROOT` that is not this repository) does
    not carry `docs/theory/SPIN32_FINITE_GAMMA.md` or `docs/PYTHIA_BRIDGE.md`
    -- only the PRIMARY document is required to exist.  This is what keeps
    every synthetic fixture in this file passing unchanged."""
    rc, out, _ = _run(gate, tmp_path, "no citations here\n", {})
    assert rc == 0, out
    assert "SPIN32_FINITE_GAMMA" not in out and "PYTHIA_BRIDGE" not in out, out


def test_real_document_range_fingerprints_are_all_live(gate):
    """Every entry in the sidecar is reached by a citation, and every citation
    has an entry.  A dead fingerprint hides the next drift exactly as a dead
    exemption does."""
    recorded = json.loads((gate.ranges_file()).read_text())
    assert recorded, "no range fingerprints recorded"
    rc, out = _capture(gate.main, [])
    assert rc == 0 and "fingerprint" not in out, out


def test_allow_list_entries_are_all_live(gate):
    """An exemption that no longer matches a citation is dead weight and hides
    the next one; every ALLOW key must be reached by the real run."""
    rc, out = _capture(gate.main, [])
    for path, line, name in gate.ALLOW:
        assert any(ln.startswith("  allowed %s:%d " % (path, line))
                   and "`%s`" % name in ln
                   for ln in out.splitlines()), (path, line, name, out)


def test_allow_list_pins_are_still_on_their_lines(gate):
    """The pin is the whole of an exemption's drift protection, so it is
    checked directly here as well as through the run."""
    for (path, line, _name), (pin, _why) in gate.ALLOW.items():
        text = (gate.ROOT / path).read_text().splitlines()
        assert pin in text[line - 1], (path, line, pin, text[line - 1])


# ----------------------------------------------------------------- named point citations

def test_declaration_line_passes(gate, tmp_path):
    doc = "`include/hdr.hpp:3` `AlphaCoreSource`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR})
    assert rc == 0, out


@pytest.mark.parametrize("line,what", [(2, "the doc comment above it"),
                                       (6, "a blank line"),
                                       (10, "another symbol's declaration")])
def test_strict_rejects_near_misses(gate, tmp_path, line, what):
    """The three shapes the +-2 window hid: a comment that names the symbol, a
    blank line, and the neighbouring declaration."""
    doc = "`include/hdr.hpp:%d` `AlphaCoreSource`\n" % line
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR})
    assert rc == 1, "%s was accepted (%s)" % (what, out)


def test_loose_accepts_what_strict_rejects(gate, tmp_path):
    """The old behaviour, kept behind --loose: this is the hole itself."""
    doc = "`include/hdr.hpp:2` `AlphaCoreSource`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR}, ["--loose"])
    assert rc == 0, out
    assert "loose" in out


def test_wrong_symbols_declaration_is_rejected(gate, tmp_path):
    """`x_coh` pointed at `f0`'s line reads as correct and is not."""
    doc = "`include/hdr.hpp:9` `x_coh`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR})
    assert rc == 1 and "x_coh" in out


def test_use_site_in_a_file_that_declares_it_is_rejected(gate, tmp_path):
    """A call is not a declaration, in the file that defines the callee."""
    hdr = HDR + "\nvoid f();\nvoid g() { f(); }\n"      # f declared :14, called :15
    doc = "`include/hdr.hpp:15` `f`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": hdr})
    assert rc == 1 and ":14" in out


def test_use_site_elsewhere_is_a_legitimate_citation(gate, tmp_path):
    """src/*.cpp cites the branch; the declaration lives in the header."""
    doc = "`src/use.cpp:4` `AlphaCoreSource`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR, "src/use.cpp": SRC})
    assert rc == 0, out


def test_comment_in_a_non_declaring_file_is_still_rejected(gate, tmp_path):
    """Line 3 of use.cpp names the symbol in a comment.  A comment is never
    the citation, whichever file it is in."""
    doc = "`src/use.cpp:3` `AlphaCoreSource`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR, "src/use.cpp": SRC})
    assert rc == 1, out


BASE = """\
class UnpolSF {
 public:
  virtual double f2(double x) const = 0;
};

class ToyF2 : public UnpolSF {
 public:
  double f2(double x) const override;
};
"""


def test_a_base_class_is_not_a_declaration(gate, tmp_path):
    """`class ToyF2 : public UnpolSF {` declares ToyF2.  Reading the base as a
    declaration of `UnpolSF` gave that name four extra "declaration" lines in
    sf.hpp alone -- lines --fix could have anchored onto."""
    files = {"include/sf.hpp": BASE}
    rc, out, _ = _run(gate, tmp_path, "`include/sf.hpp:6` `UnpolSF`\n", files)
    assert rc == 1, out
    rc, out, _ = _run(gate, tmp_path, "`include/sf.hpp:1` `UnpolSF`\n", files)
    assert rc == 0, out


def test_the_next_citation_is_not_this_ones_name(gate, tmp_path):
    """`` `a:1` and `b:2` `` is two citations, not one citation named "b:2".
    The consuming regex read the second as the first's name and then never
    matched it, so it was checked by nothing at all."""
    doc = "stated in `include/hdr.hpp:3` and `include/hdr.hpp:900`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR})
    assert rc == 1 and "only 12 lines" in out, out


# ----------------------------------------------------------------- unnamed citations

@pytest.mark.parametrize("line,what", [(6, "blank"), (5, "only punctuation")])
def test_unnamed_citation_must_land_on_something(gate, tmp_path, line, what):
    """Nothing names the target, so the one thing that can be asserted is that
    the line is a line somebody could have meant.  HDR:6 is blank, HDR:5 is
    `};`."""
    doc = "quoted at `include/hdr.hpp:%d`\n" % line
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR})
    assert rc == 1 and what in out, out


def test_unnamed_citation_on_a_real_line_needs_evidence_too(gate, tmp_path):
    """B2, 2026-09-15.  A line being non-blank asserts almost nothing -- 76.2 %
    of these citations survived a +1 drift under B1 alone -- so the citing
    sentence must quote something the cited LINE carries, R5's rule with the
    window closed to zero."""
    ok = ("the `alpha core's one-body point-nucleon density` note at "
          "`include/hdr.hpp:2`\n")
    rc, out, _ = _run(gate, tmp_path, ok, {"include/hdr.hpp": HDR})
    assert rc == 0, out
    rc, out, _ = _run(gate, tmp_path, "quoted at `include/hdr.hpp:2`\n",
                      {"include/hdr.hpp": HDR})
    assert rc == 1 and "unnamed citation" in out, out
    assert "quotes nothing" in out, out


def test_an_unnamed_citation_is_evidenced_by_a_symbol_too(gate, tmp_path):
    """The symbol branch of R5 (b), on one line: `x_coh` is used there."""
    doc = "the `x_coh` default is at `include/hdr.hpp:10`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR})
    assert rc == 0, out


def test_an_unnamed_citation_does_not_see_the_line_next_door(gate, tmp_path):
    """The whole point of B2 is that the window is ZERO: the phrase must be on
    the cited line, not on the one under it."""
    doc = "the `x_coh` default is at `include/hdr.hpp:9`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR})
    assert rc == 1 and "none of them is in the block" in out, out


def test_an_unnamed_citation_can_be_exempted_with_a_pin(gate, tmp_path):
    """A flagged citation is "this gate cannot tell", not "this is wrong", and
    the answer is a signed exemption -- keyed on the CITATION, pinned to the
    line's own text, printed on every run."""
    doc = "quoted at `include/hdr.hpp:2`\n"
    gate.UNNAMED_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:2"):
                          ("one-body point-nucleon density", "read 2026-09-15")}
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR})
    assert rc == 0 and "(unnamed) -- read 2026-09-15" in out, out
    # the pin has to be ON that line, and nowhere else in the file
    gate.UNNAMED_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:2"):
                          ("struct CoherentScenario", "read 2026-09-15")}
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR})
    assert rc == 1 and "is not inside the block" in out, out


def test_a_dead_unnamed_exemption_is_reported(gate, tmp_path):
    """The key carries the citing line, so a sentence that moves leaves an
    exemption nothing reaches -- reported like a dead fingerprint, because a
    dead exemption hides the next drift."""
    gate.UNNAMED_ALLOW = {("PHYSICS_CHANNELS.md", 9, "include/hdr.hpp:2"):
                          ("one-body point-nucleon density", "read 2026-09-15")}
    rc, out, _ = _run(gate, tmp_path, "the `x_coh` block at `include/hdr.hpp:7-10`\n",
                      {"include/hdr.hpp": HDR})
    assert rc == 1, out
    assert "no unnamed citation reaches it" in out, out


def test_an_unnamed_exemption_covers_one_citation_not_a_line(gate, tmp_path):
    """Two sentences can cite one line for different claims, so the key
    carries the citing line: the second citation is judged on its own."""
    doc = ("quoted at `include/hdr.hpp:2`\n"
           "and again, for another reason, at `include/hdr.hpp:2`\n")
    gate.UNNAMED_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:2"):
                          ("one-body point-nucleon density", "read 2026-09-15")}
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR})
    assert rc == 1, out
    assert "1 broken, 1 allow-listed" in out, out


def test_unnamed_citations_are_not_checked_in_loose_mode(gate, tmp_path):
    """--loose is for editing a file whose line numbers are in flux; it is not
    the gate."""
    rc, out, _ = _run(gate, tmp_path, "quoted at `include/hdr.hpp:2`\n",
                      {"include/hdr.hpp": HDR}, ["--loose"])
    assert rc == 0, out


# ----------------------------------------------------------------- ranges

# The citing sentence quotes `x_coh`, which the block declares: since R5 (Phase
# D item D1) `--record-ranges` REFUSES a block the citing sentence does not
# evidence, so every fixture that records a range has to say what it is citing,
# exactly as the real documents now do.
RANGED = "the `x_coh` block at `include/hdr.hpp:7-10`\n"


def _fingerprint(gate, tmp_path, files, doc=RANGED):
    """Record the sidecar the way a maintainer would, then hand it back."""
    _run(gate, tmp_path, doc, files, ["--record-ranges"])
    return json.loads((tmp_path / gate.RANGES_REL).read_text())


def test_recorded_range_passes(gate, tmp_path):
    fp = _fingerprint(gate, tmp_path, {"include/hdr.hpp": HDR})
    assert "include/hdr.hpp:7-10" in fp
    rc, out, _ = _run(gate, tmp_path, RANGED, {"include/hdr.hpp": HDR}, ranges=fp)
    assert rc == 0, out


def test_unrecorded_range_is_broken(gate, tmp_path):
    """A range with no fingerprint is a range checked by nothing, which is what
    every range in the document was until 2026-09-04."""
    rc, out, _ = _run(gate, tmp_path, RANGED, {"include/hdr.hpp": HDR}, ranges={})
    assert rc == 1 and "not fingerprinted" in out, out


@pytest.mark.parametrize("rng,why", [("7-99", "not a range inside"),
                                     ("10-7", "not a range inside"),
                                     ("6-10", "first line of the range is blank"),
                                     ("7-6", "not a range inside")])
def test_range_bounds_and_edges(gate, tmp_path, rng, why):
    doc = "block at `include/hdr.hpp:%s`\n" % rng
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR}, ranges={})
    assert rc == 1 and why in out, out


def test_range_blank_last_line_is_rejected(gate, tmp_path):
    doc = "block at `include/hdr.hpp:3-6`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR}, ranges={})
    assert rc == 1 and "last line of the range is blank" in out, out


def test_range_name_must_be_declared_inside_it(gate, tmp_path):
    """R3: a name given in the adjacent `` `path:a-b` `name` `` form."""
    files = {"include/hdr.hpp": HDR}
    fp = _fingerprint(gate, tmp_path, files, "`include/hdr.hpp:7-10` `x_coh`\n")
    rc, out, _ = _run(gate, tmp_path, "`include/hdr.hpp:7-10` `x_coh`\n", files,
                      ranges=fp)
    assert rc == 0, out
    rc, out, _ = _run(gate, tmp_path, "`include/hdr.hpp:2-3` `x_coh`\n", files,
                      ranges=fp)
    assert rc == 1 and "not declared anywhere inside the range" in out, out


def test_a_moved_block_is_caught_and_named(gate, tmp_path):
    """The shape that took 20 of the document's 93 ranges: the file grew above
    the block and the citation did not move."""
    files = {"include/hdr.hpp": HDR}
    fp = _fingerprint(gate, tmp_path, files)
    grown = "// two\n// new\n" + HDR                    # the block is now 9-12
    rc, out, _ = _run(gate, tmp_path, RANGED, {"include/hdr.hpp": grown}, ranges=fp)
    assert rc == 1 and "MOVED to :9-12" in out, out


def test_fix_relocates_a_moved_block_and_carries_the_fingerprint(gate, tmp_path):
    files = {"include/hdr.hpp": HDR}
    fp = _fingerprint(gate, tmp_path, files)
    grown = "// two\n// new\n" + HDR
    rc, out, after = _run(gate, tmp_path, RANGED, {"include/hdr.hpp": grown},
                          ["--fix"], ranges=fp)
    assert "`include/hdr.hpp:9-12`" in after, after
    assert rc == 0, out
    moved = json.loads((tmp_path / gate.RANGES_REL).read_text())
    assert "include/hdr.hpp:9-12" in moved and "include/hdr.hpp:7-10" not in moved
    assert moved["include/hdr.hpp:9-12"]["sha256"] == fp["include/hdr.hpp:7-10"]["sha256"]


def test_an_edited_block_asks_for_a_human(gate, tmp_path):
    """Editing a cited block is exactly when the row should be re-read, so it
    is reported rather than silently re-blessed."""
    files = {"include/hdr.hpp": HDR}
    fp = _fingerprint(gate, tmp_path, files)
    edited = HDR.replace("double f0 = 0.04;", "double f0 = 0.05;")
    rc, out, _ = _run(gate, tmp_path, RANGED, {"include/hdr.hpp": edited}, ranges=fp)
    assert rc == 1 and "EDITED in place" in out, out


def test_a_dead_fingerprint_is_reported(gate, tmp_path):
    fp = _fingerprint(gate, tmp_path, {"include/hdr.hpp": HDR})
    fp["include/hdr.hpp:1-2"] = {"sha256": "0" * 64, "first": "", "last": ""}
    rc, out, _ = _run(gate, tmp_path, RANGED, {"include/hdr.hpp": HDR}, ranges=fp)
    assert rc == 1 and "no longer cited" in out, out


def test_a_bare_range_inherits_the_path_before_it(gate, tmp_path):
    """`` `docs/X.md:20-40`, `:95-121` `` -- five of the document's ranges are
    written that way and matched no rule before 2026-09-04."""
    doc = ("the `x_coh` tables at `include/hdr.hpp:7-10` and the "
           "`AlphaCoreSource` one at `:2-3`\n")
    files = {"include/hdr.hpp": HDR}
    fp = _fingerprint(gate, tmp_path, files, doc)
    assert "include/hdr.hpp:2-3" in fp, fp
    rc, out, _ = _run(gate, tmp_path, doc, files, ranges=fp)
    assert rc == 0, out
    bad = dict(fp)
    bad["include/hdr.hpp:2-3"] = {"sha256": "0" * 64, "first": "x", "last": "y"}
    rc, out, _ = _run(gate, tmp_path, doc, files, ranges=bad)
    assert rc == 1 and "include/hdr.hpp:2-3" in out, out


# ----------------------------------------------------------------- --fix, allow list

def test_fix_re_anchors_onto_the_declaration(gate, tmp_path):
    """--fix must land on the declaration, not on the nearest text match --
    the doc comment on line 2 is nearer to line 1 than the enum on line 3."""
    doc = "`include/hdr.hpp:1` `AlphaCoreSource`\n"
    rc, out, after = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR}, ["--fix"])
    assert "`include/hdr.hpp:3`" in after, after
    assert rc == 0, out


# `f` declared on :1, named by a comment on :8, nothing else in between
GAP = "void f();\n" + "\n" * 6 + "/// `f` is documented here\n" + "\n" * 3 + "int q = 0;\n"


def test_fix_prefers_the_declaration_over_a_nearer_comment(gate, tmp_path):
    """The two --fix modes, on the same broken anchor.  Strict walks past the
    comment on :8 to the declaration on :1; loose stops at the comment, which
    is how an anchor lands on a doc comment in the first place."""
    doc = "`include/gap.hpp:11` `f`\n"
    rc, out, after = _run(gate, tmp_path, doc, {"include/gap.hpp": GAP}, ["--fix"])
    assert "`include/gap.hpp:1`" in after, after
    assert rc == 0, out
    rc, out, after = _run(gate, tmp_path, doc, {"include/gap.hpp": GAP},
                          ["--fix", "--loose"])
    assert "`include/gap.hpp:8`" in after, after


def test_allow_list_pin_catches_a_shift(gate, tmp_path):
    """An exemption is from the declaration rule, not from drift: the exact
    line it names must still carry the pinned text."""
    gate.ALLOW = {("src/use.cpp", 3, "AlphaCoreSource"):
                  ("branch is chosen here", "a use site, on purpose")}
    files = {"include/hdr.hpp": HDR, "src/use.cpp": SRC}
    rc, out, _ = _run(gate, tmp_path, "`src/use.cpp:3` `AlphaCoreSource`\n", files)
    assert rc == 0 and "allow-listed" in out, out
    shifted = "\n" + SRC                      # the same comment, now on line 4
    rc, out, _ = _run(gate, tmp_path, "`src/use.cpp:3` `AlphaCoreSource`\n",
                      {"include/hdr.hpp": HDR, "src/use.cpp": shifted})
    assert rc == 1 and "no longer contains its pin" in out, out


def test_an_exemption_covers_one_citation_not_a_name(gate, tmp_path):
    """Keyed on (path, name) and tested with a bare "is the name anywhere on
    the line", an exemption let its anchor drift across the whole file AND
    exempted every other citation of the same name in it.  Line 4 is a real use
    site of `AlphaCoreSource` in the same file and is NOT covered by line 3's
    exemption; it passes or fails on its own merits, and here `hdr.hpp` is
    absent, so use.cpp declares nothing and :4 stands as a use site while :5
    -- `(void)s;` -- does not."""
    gate.ALLOW = {("src/use.cpp", 3, "AlphaCoreSource"):
                  ("branch is chosen here", "a use site, on purpose")}
    files = {"src/use.cpp": SRC}
    rc, out, _ = _run(gate, tmp_path, "`src/use.cpp:4` `AlphaCoreSource`\n", files)
    assert rc == 0, out
    rc, out, _ = _run(gate, tmp_path, "`src/use.cpp:5` `AlphaCoreSource`\n", files)
    assert rc == 1 and "at another line only" in out, out


def test_missing_file_and_out_of_range_are_broken(gate, tmp_path):
    rc, out, _ = _run(gate, tmp_path, "`include/gone.hpp:3` `X`\n", {})
    assert rc == 1 and "file missing" in out
    rc, out, _ = _run(gate, tmp_path, "`include/hdr.hpp:900` `AlphaCoreSource`\n",
                      {"include/hdr.hpp": HDR})
    assert rc == 1 and "only 12 lines" in out


# ------------------------------------------------- S3: the use-site carve-out

# `helper` is declared nowhere in this file and used as code on three lines.
# Before 2026-09-04 a citation of any of the three passed, so a drift that
# moved the file under the citation landed it on a different use and was not
# reported.  The document's widest case was `python/lipolgen/__init__.py:249`
# `isotope` at 14 accepting lines, then `src/core/pipeline.cpp:531` `B1Model`
# at 11 -- 29 of the 42 use-site citations accepted more than one line.
USES = """\
#include "hdr.hpp"
void a() { helper(1); }
void b() { helper(2); }
void c() { helper(3); }
"""


def test_a_lone_use_site_needs_no_pin(gate, tmp_path):
    """One accepting line is already a pin: any drift moves it off."""
    doc = "`src/use.cpp:4` `AlphaCoreSource`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR, "src/use.cpp": SRC})
    assert rc == 0, out


def test_an_ambiguous_use_site_must_be_pinned(gate, tmp_path):
    """Three uses, so the citation alone says nothing about which one; the
    line is fingerprinted exactly as a range's block is."""
    files = {"src/uses.cpp": USES}
    doc = "`src/uses.cpp:3` `helper`\n"
    rc, out, _ = _run(gate, tmp_path, doc, files, ranges={})
    assert rc == 1 and "names it as code on 3 lines" in out, out
    fp = _fingerprint(gate, tmp_path, files, doc)
    assert fp["src/uses.cpp:3"]["name"] == "helper"
    rc, out, _ = _run(gate, tmp_path, doc, files, ranges=fp)
    assert rc == 0, out


def test_a_drifted_use_site_is_caught_and_fixable(gate, tmp_path):
    """The residual this rule exists to close: the file moves under the
    citation and it lands on a DIFFERENT use of the same name."""
    files = {"src/uses.cpp": USES}
    doc = "`src/uses.cpp:3` `helper`\n"
    fp = _fingerprint(gate, tmp_path, files, doc)
    grown = "// one\n" + USES                      # `helper(2)` is now on :4
    rc, out, _ = _run(gate, tmp_path, doc, {"src/uses.cpp": grown}, ranges=fp)
    assert rc == 1 and "MOVED to :4" in out, out
    rc, out, after = _run(gate, tmp_path, doc, {"src/uses.cpp": grown}, ["--fix"],
                          ranges=fp)
    assert "`src/uses.cpp:4`" in after and rc == 0, (after, out)
    moved = json.loads((tmp_path / gate.RANGES_REL).read_text())
    assert moved["src/uses.cpp:4"]["sha256"] == fp["src/uses.cpp:3"]["sha256"]


def test_an_edited_use_site_asks_for_a_human(gate, tmp_path):
    files = {"src/uses.cpp": USES}
    doc = "`src/uses.cpp:3` `helper`\n"
    fp = _fingerprint(gate, tmp_path, files, doc)
    edited = USES.replace("helper(2)", "helper(2, 0)")
    rc, out, _ = _run(gate, tmp_path, doc, {"src/uses.cpp": edited}, ranges=fp)
    assert rc == 1 and "EDITED in place" in out, out


def test_an_exempted_name_cannot_escape_through_the_carve_out(gate, tmp_path):
    """The one exemption of the six that was not per-citation in effect.
    `src/core/pipeline.cpp:567` `coherent` is pinned to its line, but
    `pipeline.cpp` declares no symbol `coherent`, so the carve-out accepted
    `:996` and `:1009` (`cfg_.coherent`) as "use sites" with nothing printed.
    They are ambiguous use sites, so S3 now demands a pin for them."""
    src = "// coherent is refused here\nvoid a() { cfg.coherent(1); }\n" \
          "void b() { cfg.coherent(2); }\n"
    gate.ALLOW = {("src/p.cpp", 1, "coherent"):
                  ("coherent is refused here", "the English word in a message")}
    files = {"src/p.cpp": src}
    rc, out, _ = _run(gate, tmp_path, "`src/p.cpp:1` `coherent`\n", files, ranges={})
    assert rc == 0 and "allowed" in out, out
    rc, out, _ = _run(gate, tmp_path, "`src/p.cpp:2` `coherent`\n", files, ranges={})
    assert rc == 1 and "names it as code on 2 lines" in out, out


# ------------------------------------------------- external (PYTHIA) citations

EXT_SRC = """\
// upstream
bool isHadronA = particleDataPtr->isHadron(idA) || idA == 990;
int  modeUnresolvedHadron = mode("BeamRemnants:unresolvedHadron");
"""


def _ext(gate, tmp_path, monkeypatch, doc, ranges=None, argv=(), content=EXT_SRC):
    up = tmp_path / "upstream"
    up.mkdir(exist_ok=True)
    (up / "BeamSetup.cc").write_text(content)
    monkeypatch.setenv("LIPOLGEN_PYTHIA_SRC", str(up))
    return _run(gate, tmp_path, doc, {}, argv, ranges=ranges) + (up,)


def test_external_citation_is_checked_against_the_deps_tree(gate, tmp_path, monkeypatch):
    """Five point citations and one range in the document name PYTHIA's own
    sources, which resolve outside the repository; before 2026-09-04 they
    matched no rule and were checked by nothing at all."""
    doc = "PYTHIA source pointers `BeamSetup.cc:2`\n"
    rc, out, _, up = _ext(gate, tmp_path, monkeypatch, doc, ranges={})
    assert rc == 1 and "not fingerprinted" in out and "1 external" in out, out
    _ext(gate, tmp_path, monkeypatch, doc, ranges={}, argv=["--record-ranges"])
    fp = json.loads((tmp_path / gate.RANGES_REL).read_text())
    assert "pythia8:BeamSetup.cc:2" in fp, fp
    rc, out, _, up = _ext(gate, tmp_path, monkeypatch, doc, ranges=fp)
    assert rc == 0, out
    rc, out, _, _ = _ext(gate, tmp_path, monkeypatch, doc, ranges=fp,
                         content=EXT_SRC.replace("990", "991"))
    assert rc == 1 and "no longer there" in out, out


def test_external_citation_takes_a_comma_list_and_a_range(gate, tmp_path, monkeypatch):
    doc = "pointers `BeamSetup.cc:2,3` and `BeamSetup.cc:2-3`\n"
    _ext(gate, tmp_path, monkeypatch, doc, ranges={}, argv=["--record-ranges"])
    fp = json.loads((tmp_path / gate.RANGES_REL).read_text())
    assert set(fp) == {"pythia8:BeamSetup.cc:2", "pythia8:BeamSetup.cc:3",
                       "pythia8:BeamSetup.cc:2-3"}, fp
    rc, out, _, _ = _ext(gate, tmp_path, monkeypatch, doc, ranges=fp)
    assert rc == 0 and "2 external" in out, out


def test_external_citation_is_skipped_by_name_when_pythia_is_absent(gate, tmp_path,
                                                                    monkeypatch):
    """A machine with no dependency tree must not fail the gate -- but an
    unchecked citation you cannot see is a hole, so each one is named.  All
    three names the resolver honours are cleared: with only two of them
    cleared this test failed on any machine that exported the third
    (`LIPOLGEN_DEPS_PREFIX`, the name README.md tells a pip user to set)."""
    monkeypatch.delenv("LIPOLGEN_PYTHIA_SRC", raising=False)
    monkeypatch.delenv("LIPOLGEN_DEPS", raising=False)
    monkeypatch.delenv("LIPOLGEN_DEPS_PREFIX", raising=False)
    doc = "PYTHIA source pointers `BeamSetup.cc:2`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {}, ranges={})
    assert rc == 0 and "1 skipped" in out, out
    assert "skipped BeamSetup.cc:2" in out, out


def test_skipping_does_not_drop_recorded_external_fingerprints(gate, tmp_path,
                                                               monkeypatch):
    """`--record-ranges` on a machine without PYTHIA must carry the external
    entries across, not silently delete the only record of them."""
    doc = "PYTHIA source pointers `BeamSetup.cc:2`\n"
    _ext(gate, tmp_path, monkeypatch, doc, ranges={}, argv=["--record-ranges"])
    fp = json.loads((tmp_path / gate.RANGES_REL).read_text())
    monkeypatch.delenv("LIPOLGEN_PYTHIA_SRC", raising=False)
    monkeypatch.delenv("LIPOLGEN_DEPS", raising=False)
    monkeypatch.delenv("LIPOLGEN_DEPS_PREFIX", raising=False)
    rc, out, _ = _run(gate, tmp_path, doc, {}, ["--record-ranges"], ranges=fp)
    after = json.loads((tmp_path / gate.RANGES_REL).read_text())
    assert after == fp, (fp, after)


def test_real_document_external_citations(gate):
    """On this machine (env.sh sourced) they are checked; without it they are
    skipped by name.  Either way the count is printed.

    7 since 2026-09-04 (D4): the Pomeron flavour row now cites
    `PartonDistributions.cc:2630`, the line at which PYTHIA's one `PomH1FitAB`
    class -- sets 3, 4 AND 6 -- assigns `xc = xcbar = 0.` unconditionally.
    The tree's old wording covered the LO set only.
    """
    rc, out = _capture(gate.main, [])
    assert rc == 0, out
    assert ", 7 external," in out, out
    assert out.count("  skipped ") in (0, 7), out


def test_external_citation_accepts_either_name_for_the_deps_prefix(gate, tmp_path,
                                                                   monkeypatch):
    """`env.sh` exports `$LIPOLGEN_DEPS`; CMakeLists.txt, README.md and
    docs/USAGE.md call the same directory `$LIPOLGEN_DEPS_PREFIX`.  Both are
    walked back to `<prefix>/../src/pythia8*/src`."""
    prefix = tmp_path / "deps" / "install"
    srcdir = tmp_path / "deps" / "src" / "pythia8317" / "src"
    srcdir.mkdir(parents=True)
    prefix.mkdir(parents=True)
    (srcdir / "BeamSetup.cc").write_text(EXT_SRC)
    doc = "PYTHIA source pointers `BeamSetup.cc:2`\n"
    for var in ("LIPOLGEN_DEPS", "LIPOLGEN_DEPS_PREFIX"):
        monkeypatch.delenv("LIPOLGEN_PYTHIA_SRC", raising=False)
        monkeypatch.delenv("LIPOLGEN_DEPS", raising=False)
        monkeypatch.delenv("LIPOLGEN_DEPS_PREFIX", raising=False)
        monkeypatch.setenv(var, str(prefix))
        rc, out, _ = _run(gate, tmp_path, doc, {}, ["--record-ranges"], ranges={})
        fp = json.loads((tmp_path / gate.RANGES_REL).read_text())
        assert "pythia8:BeamSetup.cc:2" in fp, (var, fp)
        rc, out, _ = _run(gate, tmp_path, doc, {}, ranges=fp)
        assert rc == 0 and "skipped" not in out, (var, out)


# ------------------------------------------------- R5: what --record-ranges will bless
#
# Phase D item D1 (2026-09-06).  R4 pins a block to the content it had WHEN IT
# WAS RECORDED and nothing checked that the content was the right content, so a
# range recorded onto the wrong block read as "0 broken" for ever: three
# `docs/CONVENTIONS.md` ranges on 2026-09-05 and two `docs/USAGE.md` ranges on
# 2026-09-06 were found that way, all of them by a human reading the row.
# `--record-ranges` now REFUSES a block the citing sentence does not evidence --
# an adjacent name that is declared or used in it, a quoted phrase that occurs
# in it, or an explicit `RANGE_ALLOW` exemption with a pin and a reason -- and a
# refused entry is simply not written, so the next ordinary run reports it as
# unfingerprinted.  Applied to the tree as committed at a3c9ecb the rule refused
# 89 of the 180 cited blocks, 20 of which were pointing at the wrong block.

def test_a_right_block_records(gate, tmp_path):
    """The citing sentence quotes `x_coh` and the block declares it."""
    rc, out, _ = _run(gate, tmp_path, RANGED, {"include/hdr.hpp": HDR},
                      ["--record-ranges"], ranges={})
    assert rc == 0, out
    assert ", 0 REFUSED" in out and "1 by quoted text" in out, out
    fp = json.loads((tmp_path / gate.RANGES_REL).read_text())
    assert "include/hdr.hpp:7-10" in fp, fp


def test_a_wrong_block_is_refused(gate, tmp_path):
    """The same sentence, pointed at the enum three lines up: `x_coh` is not
    there, so the fingerprint is not taken.  This is the whole of D1."""
    doc = "the `x_coh` block at `include/hdr.hpp:2-5`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR},
                      ["--record-ranges"], ranges={})
    assert rc == 1, out
    assert "REFUSED   PHYSICS_CHANNELS.md:1  include/hdr.hpp:2-5" in out, out
    assert "`x_coh`" in out, out
    fp = json.loads((tmp_path / gate.RANGES_REL).read_text())
    assert "include/hdr.hpp:2-5" not in fp, fp
    # ... and the ordinary run then says so, rather than passing
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR}, ranges=fp)
    assert rc == 1 and "not fingerprinted" in out, out


def test_a_sentence_that_evidences_nothing_is_refused_and_says_what_it_wanted(
        gate, tmp_path):
    doc = "block at `include/hdr.hpp:7-10`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR},
                      ["--record-ranges"], ranges={})
    assert rc == 1, out
    assert "carries no adjacent name and quotes nothing" in out, out


def test_a_refusal_drops_a_fingerprint_that_was_already_there(gate, tmp_path):
    """Recording is the only way into the sidecar, so a block that stops being
    evidenced does not stay blessed: it is dropped and the gate goes red."""
    files = {"include/hdr.hpp": HDR}
    fp = _fingerprint(gate, tmp_path, files)
    assert "include/hdr.hpp:7-10" in fp
    rc, out, _ = _run(gate, tmp_path, "block at `include/hdr.hpp:7-10`\n",
                      files, ["--record-ranges"], ranges=fp)
    assert rc == 1 and "REFUSED" in out, out
    after = json.loads((tmp_path / gate.RANGES_REL).read_text())
    assert "include/hdr.hpp:7-10" not in after, after


def test_the_quoted_phrase_may_be_prose_and_may_wrap(gate, tmp_path):
    """These documents are hard-wrapped, so a quotation that runs across a line
    break is still the quotation -- the shape that refused
    `docs/CONVENTIONS.md:72-73` while its quotation sat in front of it."""
    doc = ('the header calls them "The coherent scenario\n'
           'knobs." at `include/hdr.hpp:7-10`\n')
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR},
                      ["--record-ranges"], ranges={})
    assert rc == 0 and "1 by quoted text" in out, out


def test_rule_a_takes_a_use_as_well_as_a_declaration(gate, tmp_path):
    """An adjacent name only has to be IN the block: a range is a block
    citation and citing the lines that call a thing is legitimate.  (R3, which
    runs first, is the stricter test, so this branch is reached through the
    document only for a name R3 already accepted.)"""
    (tmp_path / "src").mkdir(parents=True, exist_ok=True)
    (tmp_path / "src" / "use.cpp").write_text(SRC)
    src = gate.Source(tmp_path / "src" / "use.cpp", "src/use.cpp")
    text = "`src/use.cpp:3-5` `AlphaCoreSource`\n"
    m = gate.RANGE.search(text)
    how, _ = gate.range_evidence(text, m, src, 4, 5, "AlphaCoreSource")
    assert how == "name"
    how, why = gate.range_evidence(text, m, src, 1, 2, "AlphaCoreSource")
    assert how is None and "neither declared nor used" in why, why


def test_an_exemption_records_what_the_sentence_cannot_say(gate, tmp_path):
    gate.RANGE_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:7-10"):
                        ("struct CoherentScenario {", "read 2026-09-06")}
    rc, out, _ = _run(gate, tmp_path, "block at `include/hdr.hpp:7-10`\n",
                      {"include/hdr.hpp": HDR}, ["--record-ranges"], ranges={})
    assert rc == 0 and "1 by RANGE_ALLOW exemption" in out, out
    fp = json.loads((tmp_path / gate.RANGES_REL).read_text())
    assert "include/hdr.hpp:7-10" in fp, fp


def test_an_exemption_pin_must_be_inside_the_block(gate, tmp_path):
    """The pin is the exemption's whole drift protection, exactly as it is for
    a point citation's ALLOW entry."""
    gate.RANGE_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:7-10"):
                        ("enum class AlphaCoreSource", "read 2026-09-06")}
    rc, out, _ = _run(gate, tmp_path, "block at `include/hdr.hpp:7-10`\n",
                      {"include/hdr.hpp": HDR}, ["--record-ranges"], ranges={})
    assert rc == 1 and "is not inside the block" in out, out


def test_an_exemption_pin_must_name_only_one_place(gate, tmp_path):
    """A pin that matches twice in the file names neither block."""
    dup = HDR + "\nstruct Other {\n  double f0 = 0.04;\n};\n"
    gate.RANGE_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:7-10"):
                        ("double f0 = 0.04;", "read 2026-09-06")}
    rc, out, _ = _run(gate, tmp_path, "block at `include/hdr.hpp:7-10`\n",
                      {"include/hdr.hpp": dup}, ["--record-ranges"], ranges={})
    assert rc == 1 and "occurs 2 times in the file" in out, out


def test_a_dead_range_exemption_is_reported(gate, tmp_path):
    """A dead exemption hides the next drift, so it is reported like a dead
    fingerprint."""
    gate.RANGE_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:1-2"):
                        ("namespace lipolgen {", "why")}
    rc, out, _ = _run(gate, tmp_path, RANGED, {"include/hdr.hpp": HDR},
                      ["--audit-ranges"], ranges={})
    assert rc == 1 and "no range citation reaches it" in out, out


def test_audit_ranges_writes_nothing(gate, tmp_path):
    """`_run(..., ranges=fp)` lays the sidecar down itself, so the file is
    untouched exactly when it still reads as that `json.dumps`."""
    files = {"include/hdr.hpp": HDR}
    fp = _fingerprint(gate, tmp_path, files)
    rc, out, _ = _run(gate, tmp_path, "block at `include/hdr.hpp:7-10`\n",
                      files, ["--audit-ranges"], ranges=fp)
    assert rc == 1 and "REFUSED" in out, out
    assert (tmp_path / gate.RANGES_REL).read_text() == json.dumps(fp)


def test_no_range_allow_suppresses_the_table(gate, tmp_path):
    """The measurement flag: it is what the published refusal counts are taken
    with, so it has to actually turn the exemptions off."""
    gate.RANGE_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:7-10"):
                        ("struct CoherentScenario {", "read 2026-09-06")}
    doc = "block at `include/hdr.hpp:7-10`\n"
    files = {"include/hdr.hpp": HDR}
    rc, out, _ = _run(gate, tmp_path, doc, files, ["--audit-ranges"], ranges={})
    assert rc == 0 and "0 REFUSED" in out, out
    rc, out, _ = _run(gate, tmp_path, doc, files,
                      ["--audit-ranges", "--no-range-allow"], ranges={})
    assert rc == 1 and "1 REFUSED" in out, out


def test_fix_says_which_exemption_key_moved(gate, tmp_path):
    """`--fix` re-points a moved block and carries its fingerprint over, but an
    exemption is keyed on the citation, so it names the key that has to move
    with it rather than leaving the next `--record-ranges` to refuse it."""
    files = {"include/hdr.hpp": HDR}
    fp = _fingerprint(gate, tmp_path, files)
    gate.RANGE_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:7-10"):
                        ("struct CoherentScenario {", "read 2026-09-06")}
    grown = "// two\n// new\n" + HDR                   # the block is now 9-12
    rc, out, _ = _run(gate, tmp_path, RANGED, {"include/hdr.hpp": grown},
                      ["--fix"], ranges=fp)
    assert rc == 0, out
    assert "should become" in out, out
    assert "'include/hdr.hpp:9-12'" in out, out


def test_evidence_is_per_citation_not_per_block(gate, tmp_path):
    """The leftover D1 left: the verdict was pooled on the block key, so a
    second citation of an already-evidenced block was accepted whatever its own
    sentence said.  Three of five deliberately wrong-block ranges passed that
    way."""
    doc = ("the `x_coh` block at `include/hdr.hpp:7-10`\n\n"
           "the gizmo calibration table at `include/hdr.hpp:7-10`\n")
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR},
                      ["--record-ranges"], ranges={})
    assert rc == 1, out
    assert "REFUSED   PHYSICS_CHANNELS.md:3  include/hdr.hpp:7-10" in out, out
    assert "2 range citation(s) checked" in out and "on 1 blocks" in out, out
    # one refused citation is enough to refuse the block: the fingerprint is
    # shared, and a block blessed on one sentence while another points at it
    # wrongly is exactly the hole this closed
    fp = json.loads((tmp_path / gate.RANGES_REL).read_text())
    assert "include/hdr.hpp:7-10" not in fp, fp


def test_an_exemption_covers_one_citation_not_one_range(gate, tmp_path):
    """`RANGE_ALLOW` is keyed on (document, citing line, range), so exempting
    one sentence does not bless the next sentence that points there."""
    gate.RANGE_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:7-10"):
                        ("struct CoherentScenario {", "read 2026-09-15")}
    doc = ("block at `include/hdr.hpp:7-10`\n\n"
           "and the gizmo table at `include/hdr.hpp:7-10`\n")
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR},
                      ["--audit-ranges"], ranges={})
    assert rc == 1, out
    assert "1 by RANGE_ALLOW exemption, 1 REFUSED" in out, out
    assert "PHYSICS_CHANNELS.md:3" in out, out


def test_a_moved_sentence_retires_its_exemption_loudly(gate, tmp_path):
    """The citing line is in the key on purpose: an exemption records that a
    human read THAT sentence, so a document edit that moves it is reported
    rather than quietly covering a sentence nobody read."""
    gate.RANGE_ALLOW = {("PHYSICS_CHANNELS.md", 1, "include/hdr.hpp:7-10"):
                        ("struct CoherentScenario {", "read 2026-09-15")}
    doc = "a new first line\n\nblock at `include/hdr.hpp:7-10`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR},
                      ["--audit-ranges"], ranges={})
    assert rc == 1, out
    assert "no range citation reaches it" in out, out
    assert "1 REFUSED" in out, out


# ---------------------------------------- R5's minimum evidence (2026-09-15)

WORDY = """\
/// the tier of the first thing
struct A { int a; };
/// the tier of the second thing
struct B { int b; };
/// the tier of the third thing
struct C { int c; };
/// the tier of the fourth thing
struct D { int d; };
/// the tier of the fifth thing
struct E { int e; };
/// the tier of the sixth thing
struct F { int f; };
"""


def test_a_common_word_does_not_evidence_a_block(gate, tmp_path):
    """`tier` is on 17 lines of `docs/USAGE.md` and `pzz` on 39 of
    `python/bindings.cpp`: a word that common cannot tell one block of that
    file from another, and twelve blocks rested on exactly that."""
    doc = "the `tier` rows at `include/hdr.hpp:1-2`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": WORDY},
                      ["--audit-ranges"], ranges={})
    assert rc == 1 and "1 REFUSED" in out, out
    # ... while the same word in a file that uses it twice still counts
    rc, out, _ = _run(gate, tmp_path, doc,
                      {"include/hdr.hpp": "/// the tier of the first thing\n"
                                          "struct A { int a; };\n"},
                      ["--audit-ranges"], ranges={})
    assert rc == 0 and "0 REFUSED" in out, out


def test_a_text_match_does_not_evidence_a_very_long_block(gate, tmp_path):
    """Measured by the decoy test, text evidence accepts a decoy 11.2 % of the
    time for a block of 1-20 lines and 50.0 % for one over 150."""
    long_file = "".join("/// line %d of the thing\n" % i for i in range(1, 401))
    doc = "the `line 7 of the thing` block at `include/hdr.hpp:1-200`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": long_file},
                      ["--audit-ranges"], ranges={})
    assert rc == 1 and "1 REFUSED" in out, out
    doc = "the `line 7 of the thing` block at `include/hdr.hpp:1-100`\n"
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": long_file},
                      ["--audit-ranges"], ranges={})
    assert rc == 0 and "0 REFUSED" in out, out


def test_a_list_item_is_its_own_citing_sentence(gate, tmp_path):
    """A bulleted list carries no blank lines, so the whole list was one
    "sentence" and a phrase quoted in the NEXT bullet evidenced this bullet's
    block -- which is how `docs/PHYSICS_CHANNELS.md:974` came to be evidenced
    by a word belonging to the bullet below it."""
    doc = ("* the enum at `include/hdr.hpp:3-4`\n"
           "* and the `AlphaCoreSource` it declares\n")
    rc, out, _ = _run(gate, tmp_path, doc, {"include/hdr.hpp": HDR},
                      ["--audit-ranges"], ranges={})
    assert rc == 1 and "1 REFUSED" in out, out
    same_bullet = "* the `AlphaCoreSource` enum at `include/hdr.hpp:3-4`\n"
    rc, out, _ = _run(gate, tmp_path, same_bullet, {"include/hdr.hpp": HDR},
                      ["--audit-ranges"], ranges={})
    assert rc == 0 and "0 REFUSED" in out, out


def test_recording_is_refused_in_loose_mode(gate, tmp_path):
    """--loose collects no fingerprints at all, so a recording pass under it
    would write an EMPTY sidecar over the real one."""
    files = {"include/hdr.hpp": HDR}
    fp = _fingerprint(gate, tmp_path, files)
    rc, out, _ = _run(gate, tmp_path, RANGED, files, ["--record-ranges", "--loose"],
                      ranges=fp)
    assert rc == 1 and "cannot be combined with --loose" in out, out
    assert json.loads((tmp_path / gate.RANGES_REL).read_text()) == fp


# ----------------------------------------------- R5 over the real documents

def test_every_real_range_is_evidenced(gate):
    """Every one of the cited blocks is either evidenced by its own sentence or
    exempted by name in RANGE_ALLOW -- and `--audit-ranges` touches nothing."""
    before = gate.ranges_file().read_bytes()
    rc, out = _capture(gate.main, ["--audit-ranges"])
    assert rc == 0, out
    assert ", 0 REFUSED" in out, out
    assert f"{len(gate.RANGE_ALLOW)} by RANGE_ALLOW exemption" in out, out
    assert gate.ranges_file().read_bytes() == before


def test_every_range_exemption_is_live_and_still_pinned(gate):
    """An exemption whose pin has slid off its block, or that no citation
    reaches any more, is dead weight that hides the next drift."""
    for key, (pin, why) in gate.RANGE_ALLOW.items():
        doc, line, ref = key
        assert isinstance(line, int) and doc.endswith(".md"), key
        path, span = ref.rsplit(":", 1)
        a, b = (int(v) for v in span.split("-"))
        src = gate.Source(gate.ROOT / path, path)
        assert gate.pin_state(src, a, b, pin) is None, (key, pin)
        assert len(why.split()) >= 5, (key, why)
    rc, out = _capture(gate.main, ["--audit-ranges"])
    assert "drop the exemption" not in out, out
    assert rc == 0, out


# ------------------------------------------------ R6: the dated run records
#
# Phase D item D2, 2026-09-15.  `docs/open_items/` -- run_2026-09-02/,
# run_2026-09-03/, run_2026-09-06/ and the older standing notes -- was gated by
# NOTHING until this date, and its 503 `path:line` citations into the live tree
# drifted with every line shift: measured against the committed a3c9ecb with
# this rule, 214 of the 448 citations the records carried there were broken.
# R6 is the relaxed rule for history: in bounds, non-blank, and one shared
# token of four characters or more with the citing sentence -- with an anchor
# that has stopped being true repaired by `(as of <commit>)` rather than moved,
# and one that lands right while the sentence names it in words the target
# never uses signed off by `(anchor read <date>)`.  Both annotations, and both
# section-marker forms, are counted and named on every run.


def _record(gate, tmp_path, records, files, argv=("--records-only",)):
    """A throw-away tree with dated run records in it."""
    for rel, text in {**files, **records}.items():
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    doc = tmp_path / "docs" / "PHYSICS_CHANNELS.md"
    doc.parent.mkdir(parents=True, exist_ok=True)
    if not doc.exists():
        doc.write_text("no citations here\n")
    gate.ROOT, gate.DOC = tmp_path, doc
    return _capture(gate.main, list(argv))


LIVE = """\
struct Widget {
  double flange_tolerance = 0.5;

  int spindle;
};
"""
# LIVE line numbers: 1 struct, 2 flange_tolerance, 3 blank, 4 spindle


def test_a_bare_point_citation_is_reached_at_all(gate, tmp_path):
    """`REF` requires a path, so a bare `` `:160` `` -- the point form of the
    inheritance these records use everywhere -- matched no rule and was checked
    by NOTHING: 152 of them in 16 records.  `RECORD_REF` is `REF` with the path
    made optional, and it is used by this pass only; widening `REF` itself
    would move the strict gate's own published counts."""
    text = ("| a | b |\n|---|---|\n"
            "| `src/w.cpp` | the spindle at `:4` |\n")
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert rc == 0 and "1 citations" in out, out
    assert gate.REF.search("`:4`") is None
    assert gate.RECORD_REF.search("`:4`") is not None


def test_records_pass_is_off_unless_asked_for(gate, tmp_path):
    """The default gate's output is what it always was: every count D1
    published is read off it with `tail -1`, and a second report printed after
    it would invalidate those recipes silently."""
    rec = {"docs/open_items/run_2026-09-02/PLAN.md": "`src/w.cpp:99` nonsense\n"}
    rc, out, _ = _run(gate, tmp_path, "no citations here\n",
                      {"src/w.cpp": LIVE, **rec})
    assert rc == 0, out
    assert "docs/open_items" not in out, out


def test_a_shared_token_passes_and_no_shared_token_is_broken(gate, tmp_path):
    ok = "the `flange_tolerance` default lives at `src/w.cpp:2`\n"
    bad = "the gizmo's calibration constant lives at `src/w.cpp:2`\n"
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": ok},
                      {"src/w.cpp": LIVE})
    assert rc == 0 and ", 0 broken" in out, out
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": bad},
                      {"src/w.cpp": LIVE})
    assert rc == 1 and "no token of 6+ characters" in out, out


def test_a_range_keeps_the_one_token_rule(gate, tmp_path):
    """The stronger POINT rule is not applied to a block: a range with one
    shared token of four characters passes, and one with none does not."""
    ok = "the `struct` and its `spindle` are at `src/w.cpp:1-4`\n"
    bad = "the gizmo calibration table is at `src/w.cpp:1-4`\n"
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": ok},
                      {"src/w.cpp": LIVE})
    assert rc == 0 and ", 0 broken" in out, out
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": bad},
                      {"src/w.cpp": LIVE})
    assert rc == 1 and "no token of 4+ characters" in out, out


def test_a_point_citation_needs_a_token_the_neighbours_do_not_carry(
        gate, tmp_path):
    """A rule the lines either side satisfy too cannot see the drift it exists
    to catch: over the records this phase leaves, the one-token rule accepts
    29.6 % of +1 drifts of a point citation and this one 12.2 %."""
    live = ("int compute() {\n"
            "  const double flange_tolerance = 0.5;   // flange_tolerance\n"
            "  return flange_tolerance;\n"
            "}\n")
    text = "the `flange_tolerance` default is at `src/w.cpp:2`\n"
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": live})
    assert rc == 1, out
    assert "absent from the lines either side of it" in out, out


def test_a_number_may_be_shorter_than_a_word(gate, tmp_path):
    """`1.848`, `0.001` and `0.05` are five characters that name one line of a
    file; refusing them cost four honest citations for nothing."""
    live = "int f() {\n  return 1;   // the ratio is 1.848 here\n}\n"
    text = "the ratio 1.848 is quoted at `src/w.cpp:2`\n"
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": live})
    assert rc == 0 and ", 0 broken" in out, out


def test_a_short_token_is_not_a_shared_token(gate, tmp_path):
    """Four characters, measured: at three the decoy rate over these records
    goes 35 % -> 60 % and the true-to-decoy ratio 2.86 -> 1.66."""
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md":
                       "the `int` it declares is at `src/w.cpp:4`\n"},
                      {"src/w.cpp": LIVE})
    assert rc == 1, out


def test_a_common_english_word_is_not_evidence(gate, tmp_path):
    """Both sides of the comparison are prose as often as they are code, so
    `which`, `that` and `because` would carry any citation anywhere."""
    assert "which" in gate.RELAX_STOP and "because" in gate.RELAX_STOP
    assert "line" not in gate.RELAX_STOP and "state" not in gate.RELAX_STOP
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md":
                       "the thing which is there, `src/w.cpp:2`\n"},
                      {"src/w.cpp": "double x;  // which one, though\n"})
    assert rc == 1, out


def test_blank_and_out_of_bounds_targets_are_broken(gate, tmp_path):
    for cite, why in (("`src/w.cpp:3`", "the target is blank"),
                      ("`src/w.cpp:99`", "out of bounds")):
        rc, out = _record(gate, tmp_path,
                          {"docs/open_items/run_2026-09-02/a.md":
                           "spindle flange_tolerance %s\n" % cite},
                          {"src/w.cpp": LIVE})
        assert rc == 1 and why in out, (cite, out)


def test_a_block_needs_one_non_blank_line_not_two(gate, tmp_path):
    """Deliberately weaker than R2, which wants a non-blank FIRST and LAST
    line: R2 protects a fingerprint and nothing here is fingerprinted."""
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md":
                       "the spindle block, `src/w.cpp:3-4`\n"},
                      {"src/w.cpp": LIVE})
    assert rc == 0, out


def test_a_citation_into_another_record_is_existence_only(gate, tmp_path):
    """One record quoting another is a pointer inside the history; the two
    move together or not at all."""
    recs = {"docs/open_items/run_2026-09-02/a.md":
            "nothing in common here, `docs/open_items/run_2026-09-03/b.md:1`\n",
            "docs/open_items/run_2026-09-03/b.md": "a line with words on it\n"}
    rc, out = _record(gate, tmp_path, recs, {})
    assert rc == 0 and "1 into another record" in out, out
    recs["docs/open_items/run_2026-09-03/b.md"] = "\n"
    rc, out = _record(gate, tmp_path, recs, {})
    assert rc == 1 and "the target is blank" in out, out


def test_a_target_outside_the_tree_is_skipped_by_name(gate, tmp_path):
    """Nine of them, none a LiPolGen path: eSTARlight's and PYTHIA's own
    sources, and the predecessor project's notes.  Named on every run, so
    deleting a live file shows up here rather than as silence."""
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md":
                       "eSTARlight's own `src/nucleons.cpp:118` reader\n"}, {})
    assert rc == 0, out
    assert "src/nucleons.cpp" in out and "no such file under ROOT" in out, out


def test_asof_covers_its_own_citation_and_not_its_neighbour(gate, tmp_path):
    """Keyed on the citing SENTENCE the annotation silently exempted every
    other citation within RELAX_WIN characters -- most of a table row."""
    text = ("gizmo calibration at `src/w.cpp:2` (as of a3c9ecb) and at "
            "`src/w.cpp:4`\n")
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert rc == 1, out
    assert "1 annotated historical" in out, out
    assert "src/w.cpp:4" in out and "src/w.cpp:2  (" not in out, out


def test_asof_sits_after_an_adjacent_name(gate, tmp_path):
    text = "gizmo at `src/w.cpp:2` `flange_tolerance` (as of a3c9ecb)\n"
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert rc == 0 and "1 annotated historical" in out, out


def test_asof_covers_a_target_that_has_since_gone_blank(gate, tmp_path):
    """Four of the real annotations land on a line that is now blank.  That is
    not a defect in the annotation, it is the thing the annotation is FOR --
    while `(anchor read <date>)`, which claims a human FOLLOWED the anchor,
    still has to land on something."""
    blank = "gizmo at `src/w.cpp:3` (as of a3c9ecb)\n"
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": blank},
                      {"src/w.cpp": LIVE})
    assert rc == 0 and "1 annotated historical" in out, out
    read = "gizmo at `src/w.cpp:3` (anchor read 2026-09-15)\n"
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": read},
                      {"src/w.cpp": LIVE})
    assert rc == 1 and "the target is blank" in out, out
    out_of_bounds = "gizmo at `src/w.cpp:99` (as of a3c9ecb)\n"
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": out_of_bounds},
                      {"src/w.cpp": LIVE})
    assert rc == 1 and "out of bounds" in out, out


def test_anchor_read_signs_off_a_citation_the_rule_cannot_check(gate, tmp_path):
    text = ("the gizmo's own constant, `src/w.cpp:2` "
            '(anchor read 2026-09-15 "flange_tolerance = 0.5")\n')
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert rc == 0 and "1 read by hand" in out, out


def test_an_anchor_read_point_without_a_pin_is_broken(gate, tmp_path):
    """Without the pin the annotation exempted the citation from (ii) and left
    only (i), so the anchor was accepted on ANY non-blank line of the file: 28
    point citations were in that state before 2026-09-15."""
    text = "the gizmo's own constant, `src/w.cpp:2` (anchor read 2026-09-15)\n"
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert rc == 1 and "must carry the text the reader saw" in out, out


def test_an_anchor_read_pin_that_has_slid_off_is_broken(gate, tmp_path):
    """The pin is what expires the reading when the line changes under it."""
    text = ("the gizmo's own constant, `src/w.cpp:4` "
            '(anchor read 2026-09-15 "flange_tolerance = 0.5")\n')
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert rc == 1 and "is no longer on that line" in out, out


def test_an_anchor_read_pin_may_be_written_in_backticks(gate, tmp_path):
    """A code line full of double quotes needs the other delimiter."""
    live = 'int f() {\n  return g("wide", "open");\n}\n'
    text = ("the call, `src/w.cpp:2` "
            "(anchor read 2026-09-15 `g(\"wide\"`)\n")
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": live})
    assert rc == 0 and "1 read by hand" in out, out


def test_a_section_marker_stops_at_the_next_heading(gate, tmp_path):
    text = ("## one\n\n**Anchors read 2026-09-15** because the rows quote "
            "keys.\n\ngizmo `src/w.cpp:1-2`\n\n## two\n\n"
            "gizmo `src/w.cpp:1-4`\n")
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert rc == 1, out
    assert "1 read by hand" in out and "src/w.cpp:1-4" in out, out
    assert "covering 1 citation(s)" in out, out


def test_a_section_marker_does_not_cover_a_point_citation(gate, tmp_path):
    """A marker signs for blocks; a point citation says which LINE, and the
    only thing that holds it there is its own pin."""
    text = ("## one\n\n**Anchors read 2026-09-15** because the rows quote "
            "keys.\n\ngizmo `src/w.cpp:2`\n")
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert rc == 1 and "a section marker covers ranges only" in out, out


def test_a_table_column_gives_a_bare_citation_its_path(gate, tmp_path):
    text = ("| # | site | form |\n|---|---|---|\n"
            "| 1 | `src/w.cpp` | the `flange_tolerance` default at `:2` |\n"
            "| 2 | | the spindle at `:4` |\n")
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert rc == 0 and "2 checked against a live file" in out, out


def test_a_row_naming_several_files_has_no_file_column(gate, tmp_path):
    """Guessing the last one is how `:133-136` landed on
    `vmc_tag_fractions.py` when the sentence meant the generated
    `vmc_reconciliation.md`, and reported a sound citation broken."""
    text = ("| a | b |\n|---|---|\n"
            "| `src/w.cpp` and `src/v.cpp` | the spindle at `:4` |\n")
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE, "src/v.cpp": LIVE})
    assert rc == 0, out
    assert "the path is carried by the prose" in out, out


def test_a_bare_filename_stops_the_same_line_rule(gate, tmp_path):
    """`` (`tagged.hpp:255-257` -> `:364-375`) ``: the main gate's rule reads
    past the un-prefixed name and inherits a path four citations to the left.
    Two sound citations in `run_2026-09-03/STATUS.md` were reported broken by
    exactly that on R6's first run."""
    text = ("`src/w.cpp:2` and then, re-pointed by hand "
            "(`tagged.hpp:255-257` -> `:364-375`)\n")
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert "src/w.cpp:364" not in out, out
    assert "the path is carried by the prose" in out, out


def test_a_bare_citation_inherits_past_another_bare_one(gate, tmp_path):
    text = ("| a | b |\n|---|---|\n"
            "| `src/w.cpp` | `:1-2` the flange_tolerance block, `:4` spindle |\n")
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md": text},
                      {"src/w.cpp": LIVE})
    assert rc == 0 and "2 checked against a live file" in out, out


def test_nothing_is_fingerprinted_by_the_records_pass(gate, tmp_path):
    """A fingerprint is a promise that a block will not change, and a record
    makes no such promise about anyone else's file."""
    fp = {"already": {"sha256": "x"}}
    rc, out = _record(gate, tmp_path,
                      {"docs/open_items/run_2026-09-02/a.md":
                       "flange_tolerance at `src/w.cpp:2`\n"},
                      {"src/w.cpp": LIVE, gate.RANGES_REL: json.dumps(fp)})
    assert rc == 0, out
    assert json.loads((tmp_path / gate.RANGES_REL).read_text()) == fp


# ------------------------------------------------ R6 over the real records

def test_the_real_records_pass_the_relaxed_gate(gate):
    """The whole of `docs/open_items/`, against the live tree."""
    rc, out = _capture(gate.main, ["--records-only"])
    line = out.splitlines()[0]
    assert rc == 0, out
    assert line.endswith("[docs/open_items]"), line
    assert ", 0 broken" in line, line
    n_cited = int(line.split()[0])
    n_records = int(line.split()[3])
    assert n_cited >= 640, line
    assert n_records >= 40 and n_records == len(gate.record_files()), line


def test_the_records_pass_reaches_every_record_and_most_citations(gate):
    """A rule that quietly stopped resolving paths would pass by checking
    nothing.  Both halves are pinned: how many citations are CHECKED against a
    live file, and how few are skipped for want of a path."""
    import re as _re
    rc, out = _capture(gate.main, ["--records-only"])
    line = out.splitlines()[0]
    checked = int(_re.search(r"(\d+) checked against a live file", line).group(1))
    skipped = int(_re.search(r"(\d+) skipped", line).group(1))
    assert checked >= 240, line
    assert skipped <= 200, line
    assert rc == 0


def test_every_as_of_commit_in_the_records_resolves(gate):
    """`(as of <commit>)` is a claim about a commit; a commit that is not in
    this repository makes it unverifiable, which is worse than no annotation."""
    import subprocess
    root = str(gate.ROOT)
    if subprocess.run(["git", "-C", root, "rev-parse", "--git-dir"],
                      capture_output=True).returncode != 0:
        pytest.skip("not a git checkout")
    seen = set()
    for doc in gate.record_files():
        for m in gate.ASOF_INLINE.finditer(doc.read_text()):
            seen.add(m.group(1))
        for _lo, _hi, sha in gate.marked_sections(doc.read_text(),
                                                  gate.ASOF_SECTION):
            seen.add(sha)
    assert seen, "no historical annotation in any record"
    for sha in sorted(seen):
        r = subprocess.run(["git", "-C", root, "cat-file", "-e", sha + "^{commit}"],
                           capture_output=True)
        assert r.returncode == 0, sha


def test_every_section_marker_in_the_records_covers_something(gate):
    """A dead marker hides the next drift exactly as a dead fingerprint and a
    dead exemption do."""
    rc, out = _capture(gate.main, ["--records-only"])
    assert rc == 0, out
    lines = [ln for ln in out.splitlines() if ln.startswith("  marker ")]
    assert lines, out
    for ln in lines:
        n = int(ln.split("covering")[1].split()[0])
        assert n > 0, ln


def test_every_unnamed_exemption_is_live_and_still_pinned(gate):
    """B2's exemptions carry the same contract R5's do: the pin is ON the cited
    line and nowhere else in the file, the reason is a sentence, and every
    entry is reached by a citation that is still there."""
    rc, out = _capture(gate.main, [])
    assert rc == 0, out
    for (doc, line, ref), (pin, why) in gate.UNNAMED_ALLOW.items():
        assert isinstance(line, int) and doc.endswith(".md"), (doc, line, ref)
        path, at = ref.rsplit(":", 1)
        src = gate.Source(gate.ROOT / path, path)
        assert gate.pin_state(src, int(at), int(at), pin) is None, (ref, pin)
        assert len(why.split()) >= 5, (ref, why)
        assert f"{ref}  (unnamed) -- {why}" in out, ref
    assert out.count("(unnamed) --") == len(gate.UNNAMED_ALLOW), out


def test_every_unnamed_citation_is_evidenced_or_exempted(gate):
    """The population B2 governs, counted the way D3.1 counts it."""
    cache = {}
    n = flagged = 0
    for doc in [gate.DOC] + [gate.ROOT / p for p, _ in gate.EXTRA_DOCS]:
        text = doc.read_text()
        for m in gate.REF.finditer(text):
            path, line, name = m.group(1), int(m.group(2)), m.group(3)
            if name and gate.CITATION.match(name):
                name = None
            if name:
                continue
            if path not in cache:
                cache[path] = gate.Source(gate.ROOT / path, path)
            src = cache[path]
            n += 1
            assert src.substantive(line), f"{path}:{line}"
            if gate.range_evidence(text, m, src, line, line, None)[0] is None:
                flagged += 1
                key = (doc.name, text.count("\n", 0, m.start()) + 1,
                       f"{path}:{line}")
                assert key in gate.UNNAMED_ALLOW, key
    assert n == 43 and flagged == len(gate.UNNAMED_ALLOW), (n, flagged)


def test_every_anchor_read_pin_in_the_records_is_still_on_its_line(gate):
    """The pin is what expires a reading when the line changes under it, so a
    stale pin is the one thing this annotation must not have."""
    pins = 0
    for rec in gate.record_files():
        text = rec.read_text()
        cites = [(m, int(m.group(2)), int(m.group(2)))
                 for m in gate.RECORD_REF.finditer(text)]
        cites += [(m, int(m.group(2)), int(m.group(3)))
                  for m in gate.RANGE.finditer(text)]
        for m, a, b in cites:
            got = gate.READ_INLINE.match(text, m.end())
            if got is None:
                continue
            pin = got.group(2) or got.group(3)
            if pin is None:
                assert a != b, f"{rec}:{a} point anchor with no pin"
                continue
            path = m.group(1)
            if path is None:
                path, online = gate.record_inherited_path(text, m.start())
                if path is None and not online:
                    path = gate.column_path(text, m.start())
            lines = (gate.ROOT / path).read_text(errors="replace").splitlines()
            assert pin in lines[a - 1], (rec.name, path, a, pin)
            pins += 1
    # 28 annotations that existed before 2026-09-15 and were pinned when the
    # pin became compulsory, plus the 10 D3.4's repairs wrote and the 10 in
    # this phase's own record, which gates itself like any other
    assert pins == 53  # 48 at the D3 count + 5 pinned anchor-read points added in phase_D_numbers.md D4 (2026-09-15), pins


# The six `(as of <commit>)` anchors that do NOT pass the gate's own bar at the
# commit they name, and were established by reading instead (D2.5's table).
# Every one was read at its commit and is right there by content; the list is
# here so that a SEVENTH cannot appear without someone signing for it.
ASOF_READ_BY_HAND = {
    ("run_2026-09-02/design_C_tensor_rc.md", "src/hepmc/hepmc_writer.cpp",
     135, 135, "0145885"),
    ("run_2026-09-03/phase_A_gate_mechanics.md", "docs/OPEN_ITEMS_SOLUTIONS.md",
     414, 414, "66dcda2"),
    ("run_2026-09-03/phase_A_miller_normalisation.md",
     "include/lipolgen/constants.hpp", 86, 86, "a98f0a0"),
    ("run_2026-09-03/phase_C_numbers.md", "python/lipolgen/cli.py",
     604, 604, "adec442"),
    ("run_2026-09-06/phase_C_survey.md", "validation/o5_a2_reach.py",
     466, 468, "7f68339"),
    # ... and this phase's own record, where D3.2's two tables quote the line
    # `include/lipolgen/constants.hpp:86` HELD at a3c9ecb -- the
    # b1-normalisation comment, which is the whole point of the row
    ("run_2026-09-06/phase_D_numbers.md", "include/lipolgen/constants.hpp",
     86, 86, "a3c9ecb"),
}


def test_every_as_of_anchor_passes_the_bar_at_its_commit(gate):
    """`(as of <commit>)` says the anchor was right AT THAT COMMIT, and for 93
    of the 107 that is machine-derived: the gate's own rule -- two shared tokens
    or one of six characters -- passes there.  The other nine are named above,
    with what each target reads at that commit, in D2.5."""
    checked = by_hand = 0
    for rec in gate.record_files():
        rel = str(rec.relative_to(gate.ROOT)).replace("docs/open_items/", "")
        text = rec.read_text()
        marks = gate.marked_sections(text, gate.ASOF_SECTION)
        cites = [(m, int(m.group(2)), int(m.group(2)))
                 for m in gate.RECORD_REF.finditer(text)]
        cites += [(m, int(m.group(2)), int(m.group(3)))
                  for m in gate.RANGE.finditer(text)]
        for m, a, b in cites:
            inline = gate.ASOF_INLINE.match(text, m.end())
            commit = inline.group(1) if inline else next(
                (n for lo, hi, n in marks if lo <= m.start() < hi), None)
            if commit is None:
                continue
            path = m.group(1)
            if path is None:
                path, online = gate.record_inherited_path(text, m.start())
                if path is None and not online:
                    path = gate.column_path(text, m.start())
            if path is None:
                continue
            checked += 1
            got = subprocess.run(["git", "show", f"{commit}:{path}"],
                                 cwd=gate.ROOT, capture_output=True, text=True)
            assert got.returncode == 0, (rel, commit, path)
            lines = got.stdout.splitlines()
            assert 1 <= a <= b <= len(lines), (rel, commit, path, a, b)
            sentence = gate.relax_sentence(text, m.start(), m.end())
            shared = (gate.relax_tokens(sentence)
                      & gate.relax_tokens(" ".join(lines[a - 1:b])))
            if len(shared) >= 2 or any(len(t) >= 6 for t in shared):
                continue
            by_hand += 1
            assert (rel, path, a, b, commit) in ASOF_READ_BY_HAND, (
                rel, path, a, b, commit, sorted(shared))
    assert checked == 107, checked
    assert by_hand == 9, by_hand      # the six of D2.5 + three in this record
