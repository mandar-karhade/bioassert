"""Layer-5 noise expansion + Bug 3a fix (Sub-phase 3.10).

Three concerns:

1. **Bug 3a fix** — ``whitespace`` substitutions (``tab``, ``newline_mid_sentence``,
   ``double_space``) must never land in prose-interior positions. A whitespace
   substitution is only permitted at a space character that is immediately
   adjacent to a labeled span (slot-boundary). The existing
   ``no_space_after_punct`` mode is a deletion — its scope is unchanged.

2. **OCR corruption** — new ``ocr_corruption`` sub-category. Length-preserving
   character substitutions (``l``↔``1``, ``O``↔``0``, ``S``↔``5``, ``B``↔``8``,
   ``Z``↔``2``, ``I``↔``l``) applied only to characters outside labeled spans.

3. **PDF hyphen-break artifact** — new ``pdf_artifact`` sub-category. Inserts
   ``-\\n`` at a random interior position of a non-labeled token, simulating a
   hyphenated line break at a PDF page edge. Requires span shifting like
   ``double_space``.
"""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.post_process import (
    OCR_CHAR_SWAPS,
    PDF_HYPHEN_LINEBREAK,
    _apply_ocr_corruption,
    _apply_pdf_artifact,
    _apply_whitespace,
    _slot_boundary_space_positions,
    apply_technical_noise,
)
from bioassert.generator.renderer import (
    AssertionFact,
    RenderedRecord,
    render_l1_record,
)


# ======================================================================
# Bug 3a — whitespace substitution restricted to slot-boundary spaces
# ======================================================================


def test_slot_boundary_positions_adjacent_to_spans_only() -> None:
    """Returned positions must have a span ending exactly at the space OR
    a span starting exactly after the space."""
    sentence = "EGFR was found to be positive by NGS."
    spans = {"gene": (0, 4), "status": (21, 29), "method": (33, 36)}
    positions = _slot_boundary_space_positions(sentence, spans)
    # Space positions in the sentence: 4, 8, 14, 17, 20, 29, 32.
    # Adjacent-to-span spaces: 4 (after gene), 20 (before status),
    # 29 (after status), 32 (before method).
    assert set(positions) == {4, 20, 29, 32}
    # Prose-interior spaces 8 ("was "), 14 ("found "), 17 ("to ")
    # must NOT be eligible.
    for prose_pos in (8, 14, 17):
        assert prose_pos not in positions


def test_slot_boundary_positions_handle_adjacent_slots() -> None:
    """Two slots separated by a single space: that space is both
    after-slot and before-slot. Return it once."""
    sentence = "EGFR positive."
    spans = {"gene": (0, 4), "status": (5, 13)}
    positions = _slot_boundary_space_positions(sentence, spans)
    assert positions == [4]


def test_whitespace_tab_only_at_slot_boundary(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Apply tab mode repeatedly to L1 records and check every tab that
    appears is adjacent to a labeled span — never mid-prose."""
    rng = random.Random(1001)
    tabs_seen = 0
    for _ in range(500):
        rec = render_l1_record("EGFR", biomarkers, common, rng)
        original = rec.sentence
        spans = dict(rec.assertions[0].spans)
        new_sentence, _new_spans = _apply_whitespace(
            original, spans, "tab", rng
        )
        if new_sentence == original:
            continue
        tab_pos = new_sentence.index("\t")
        tabs_seen += 1
        original_spans = rec.assertions[0].spans
        adjacent = any(
            tab_pos == end or tab_pos + 1 == start
            for start, end in original_spans.values()
        )
        assert adjacent, (
            f"tab landed at pos {tab_pos} in prose interior of {original!r}"
        )
    assert tabs_seen >= 100, (
        "tab mode should fire frequently enough to exercise boundary rule"
    )


def test_whitespace_newline_only_at_slot_boundary(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(1003)
    for _ in range(400):
        rec = render_l1_record("KRAS", biomarkers, common, rng)
        original = rec.sentence
        spans = dict(rec.assertions[0].spans)
        new_sentence, _ = _apply_whitespace(
            original, spans, "newline_mid_sentence", rng
        )
        if new_sentence == original:
            continue
        nl_pos = new_sentence.index("\n")
        original_spans = rec.assertions[0].spans
        assert any(
            nl_pos == end or nl_pos + 1 == start
            for start, end in original_spans.values()
        ), f"newline at {nl_pos} in prose interior of {original!r}"


def test_whitespace_double_space_only_at_slot_boundary(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """double_space insertion must also honor boundary rule."""
    rng = random.Random(1005)
    for _ in range(400):
        rec = render_l1_record("BRAF", biomarkers, common, rng)
        original = rec.sentence
        spans = dict(rec.assertions[0].spans)
        new_sentence, new_spans = _apply_whitespace(
            original, spans, "double_space", rng
        )
        if new_sentence == original:
            continue
        assert len(new_sentence) == len(original) + 1
        # ``pick`` was the index of the single space in ``original`` that got
        # doubled. In ``new_sentence`` that single space is now a "  " pair
        # starting at the same index. Recover pick by finding the first "  "
        # in new_sentence (L1 canonical renders never contain a double space).
        pick = new_sentence.index("  ")
        assert original[pick] == " "
        original_spans = rec.assertions[0].spans
        assert any(
            pick == end or pick + 1 == start
            for start, end in original_spans.values()
        ), (
            f"double_space inserted at prose-interior pos {pick} "
            f"in {original!r}"
        )


# ======================================================================
# OCR corruption — length-preserving substitutions outside labeled spans
# ======================================================================


def test_ocr_char_swap_pool_length_preserving() -> None:
    """Every swap in the pool maps a single char to a single char."""
    for src, dst in OCR_CHAR_SWAPS.items():
        assert len(src) == 1
        assert len(dst) == 1


def test_ocr_char_swap_pool_nonempty() -> None:
    assert len(OCR_CHAR_SWAPS) >= 6


def test_ocr_light_mode_swaps_at_most_two_chars() -> None:
    rng = random.Random(1007)
    sentence = "Result for EGFR: positive by NGS."
    spans = {"gene": (11, 15)}  # "EGFR" is labeled, must stay intact.
    new_sentence = _apply_ocr_corruption(sentence, spans, "light", rng)
    # Preserves length and preserves labeled span content.
    assert len(new_sentence) == len(sentence)
    s, e = spans["gene"]
    assert new_sentence[s:e] == sentence[s:e]
    # Mismatch count reflects 1-2 swapped chars.
    diffs = sum(1 for a, b in zip(sentence, new_sentence) if a != b)
    assert 0 <= diffs <= 2


def test_ocr_moderate_mode_swaps_up_to_five_chars() -> None:
    rng = random.Random(1009)
    sentence = "Result for EGFR: positive by NGS."
    spans = {"gene": (11, 15)}
    new_sentence = _apply_ocr_corruption(sentence, spans, "moderate", rng)
    assert len(new_sentence) == len(sentence)
    s, e = spans["gene"]
    assert new_sentence[s:e] == sentence[s:e]
    diffs = sum(1 for a, b in zip(sentence, new_sentence) if a != b)
    assert diffs <= 5


def test_ocr_canonical_is_identity() -> None:
    rng = random.Random(1011)
    sentence = "EGFR was positive."
    spans = {"gene": (0, 4), "status": (9, 17)}
    new_sentence = _apply_ocr_corruption(sentence, spans, "canonical", rng)
    assert new_sentence == sentence


def test_ocr_never_touches_labeled_chars(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Across many runs no swap ever lands inside a labeled span."""
    rng = random.Random(1013)
    for _ in range(200):
        rec = render_l1_record("EGFR", biomarkers, common, rng)
        original = rec.sentence
        spans = dict(rec.assertions[0].spans)
        new_sentence = _apply_ocr_corruption(
            original, spans, "moderate", rng
        )
        assert len(new_sentence) == len(original)
        for name, (s, e) in spans.items():
            assert new_sentence[s:e] == original[s:e], (
                f"OCR corrupted labeled span {name} in {original!r}"
            )


# ======================================================================
# PDF hyphen-break artifact — inserts "-\n" inside a non-labeled token
# ======================================================================


def test_pdf_hyphen_linebreak_constant() -> None:
    assert PDF_HYPHEN_LINEBREAK == "-\n"


def test_pdf_hyphen_linebreak_inserts_two_chars() -> None:
    rng = random.Random(1015)
    sentence = "Result for EGFR: positive by NGS."
    spans = {"gene": (11, 15), "status": (17, 25)}
    new_sentence, new_spans = _apply_pdf_artifact(
        sentence, spans, "hyphen_linebreak", rng
    )
    if new_sentence == sentence:
        # No eligible token-interior position picked on this roll; try more.
        for seed in range(1016, 1050):
            rng = random.Random(seed)
            new_sentence, new_spans = _apply_pdf_artifact(
                sentence, spans, "hyphen_linebreak", rng
            )
            if new_sentence != sentence:
                break
    assert new_sentence != sentence
    assert PDF_HYPHEN_LINEBREAK in new_sentence
    assert len(new_sentence) == len(sentence) + 2


def test_pdf_hyphen_linebreak_preserves_labeled_spans() -> None:
    """The inserted ``-\\n`` must land inside a non-labeled token; labeled
    spans must continue to point at the same surface after span shifting."""
    rng = random.Random(1017)
    sentence = "Result for EGFR: positive by NGS."
    spans = {"gene": (11, 15), "status": (17, 25), "method": (29, 32)}
    # Exercise enough seeds to land the artifact somewhere.
    for seed in range(1020, 1200):
        rng = random.Random(seed)
        new_sentence, new_spans = _apply_pdf_artifact(
            sentence, spans, "hyphen_linebreak", rng
        )
        if new_sentence == sentence:
            continue
        # Each labeled span's surface is preserved under the new spans.
        for name, (s, e) in new_spans.items():
            assert new_sentence[s:e] == sentence[spans[name][0]:spans[name][1]]
        return
    pytest.fail("no hyphen_linebreak fired across many seeds")


def test_pdf_canonical_is_identity() -> None:
    rng = random.Random(1019)
    sentence = "EGFR was positive."
    spans = {"gene": (0, 4), "status": (9, 17)}
    new_sentence, new_spans = _apply_pdf_artifact(
        sentence, spans, "canonical", rng
    )
    assert new_sentence == sentence
    assert new_spans == spans


# ======================================================================
# Corpus-level Bug 3a rate (exit gate) — stress check
# ======================================================================


def test_bug_3a_rate_is_zero_on_l1_prose_under_tab_mode(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Across many L1 records with tab whitespace mode applied, verify no
    tab lands in a prose-interior position. This is the sub-phase exit gate.
    """
    rng = random.Random(1021)
    for _ in range(1000):
        rec = render_l1_record("EGFR", biomarkers, common, rng)
        original = rec.sentence
        spans = dict(rec.assertions[0].spans)
        new_sentence, _ = _apply_whitespace(original, spans, "tab", rng)
        if "\t" not in new_sentence:
            continue
        tab_pos = new_sentence.index("\t")
        original_spans = rec.assertions[0].spans
        adjacent = any(
            tab_pos == end or tab_pos + 1 == start
            for start, end in original_spans.values()
        )
        assert adjacent, (
            f"Bug 3a regression: tab at pos {tab_pos} of prose in {original!r}"
        )
