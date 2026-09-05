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

The gate is a standalone script with no package to import, so it is loaded by
path; nothing here needs the compiled extension.
"""
import importlib.util
import json
import os

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


def test_unnamed_citation_on_a_real_line_passes(gate, tmp_path):
    rc, out, _ = _run(gate, tmp_path, "quoted at `include/hdr.hpp:2`\n",
                      {"include/hdr.hpp": HDR})
    assert rc == 0, out


# ----------------------------------------------------------------- ranges

RANGED = "block at `include/hdr.hpp:7-10`\n"


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
    doc = "tables at `include/hdr.hpp:7-10` and `:2-3`\n"
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
