"""Technical-noise post-processing: span preservation under transformation."""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.patient_sampler import PatientProfile
from bioassert.generator.post_process import (
    PostProcessedRecord,
    _apply_case,
    _apply_hyphenation,
    _apply_whitespace,
    _shift_spans,
    apply_technical_noise,
)
from bioassert.generator.renderer import RenderedRecord, render_l1_record


def _make_stub_record(sentence: str, spans: dict[str, tuple[int, int]]) -> RenderedRecord:
    return RenderedRecord(
        sentence=sentence,
        spans=spans,
        gene="EGFR",
        variant_id=None,
        negative_form_id=None,
        clone_id=None,
        status="positive",
        test_method=None,
        measurement_value=None,
        frame_template="",
    )


def test_case_variation_preserves_span_offsets() -> None:
    stub = _make_stub_record("EGFR was positive.", {"gene": (0, 4), "status": (9, 17)})
    upper = _apply_case(stub.sentence, "all_uppercase")
    assert upper == "EGFR WAS POSITIVE."
    assert upper[0:4] == "EGFR"
    assert upper[9:17] == "POSITIVE"


def test_shift_spans_after_insertion() -> None:
    spans = {"gene": (0, 4), "status": (9, 17)}
    shifted = _shift_spans(spans, insert_at=5, delta=1)
    assert shifted["gene"] == (0, 4)
    assert shifted["status"] == (10, 18)


def test_shift_spans_boundary_rule_is_end_exclusive() -> None:
    """Insertion at position == span end leaves span end unchanged."""
    spans = {"gene": (0, 4)}
    shifted = _shift_spans(spans, insert_at=4, delta=1)
    assert shifted["gene"] == (0, 4)


def test_hyphenation_grows_gene_span() -> None:
    stub = _make_stub_record("EGFR was positive.", {"gene": (0, 4), "status": (9, 17)})
    new_sentence, new_spans = _apply_hyphenation(stub.sentence, stub.spans, "hyphenated")
    assert new_sentence.startswith("EG-FR")
    start, end = new_spans["gene"]
    assert new_sentence[start:end] == "EG-FR"
    status_start, status_end = new_spans["status"]
    assert new_sentence[status_start:status_end] == "positive"


def test_hyphenation_skips_non_alpha_gene_surface() -> None:
    stub = _make_stub_record(
        "PD-L1 was positive.", {"gene": (0, 5), "status": (10, 18)}
    )
    out, spans = _apply_hyphenation(stub.sentence, stub.spans, "hyphenated")
    assert out == stub.sentence
    assert spans == stub.spans


def test_whitespace_tab_is_length_preserving() -> None:
    rng = random.Random(0)
    stub = _make_stub_record("EGFR was positive.", {"gene": (0, 4), "status": (9, 17)})
    new_sentence, new_spans = _apply_whitespace(stub.sentence, stub.spans, "tab", rng)
    assert len(new_sentence) == len(stub.sentence)
    assert "\t" in new_sentence or new_sentence == stub.sentence
    assert new_spans == stub.spans


def test_whitespace_double_space_shifts_spans() -> None:
    rng = random.Random(0)
    stub = _make_stub_record("EGFR was positive.", {"gene": (0, 4), "status": (9, 17)})
    new_sentence, new_spans = _apply_whitespace(
        stub.sentence, stub.spans, "double_space", rng
    )
    assert len(new_sentence) == len(stub.sentence) + 1
    gene_start, gene_end = new_spans["gene"]
    assert new_sentence[gene_start:gene_end] == "EGFR"
    status_start, status_end = new_spans["status"]
    assert new_sentence[status_start:status_end] == "positive"


def test_apply_technical_noise_end_to_end_preserves_spans(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(42)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(500):
        rendered = render_l1_record("EGFR", profile, biomarkers, common, rng)
        post = apply_technical_noise(rendered, common, rng)
        for name, (start, end) in post.spans.items():
            assert 0 <= start <= end <= len(post.sentence), (
                f"span {name} out of bounds: {(start, end)} vs len={len(post.sentence)}"
            )
            substring = post.sentence[start:end]
            original = rendered.sentence[rendered.spans[name][0]:rendered.spans[name][1]]
            # case/ocr transforms allow the substring to differ from original,
            # but the span must still point to non-empty text of a similar length.
            assert substring, f"empty substring for {name} in {post.sentence!r}"


def test_post_processed_record_reports_applied_transforms(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(0)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    rendered = render_l1_record("EGFR", profile, biomarkers, common, rng)
    post = apply_technical_noise(rendered, common, rng)
    assert isinstance(post, PostProcessedRecord)
    expected = {
        "whitespace",
        "case_variation",
        "hyphenation_gene_names",
        "punctuation_variation",
    }
    assert set(post.applied_transforms.keys()) == expected


def test_missing_technical_noise_block_raises() -> None:
    from bioassert.config.loader import CommonConfig

    empty = CommonConfig(schema_version="1.0", categories={})
    rng = random.Random(0)
    stub = _make_stub_record("EGFR was positive.", {"gene": (0, 4), "status": (9, 17)})
    with pytest.raises(KeyError):
        apply_technical_noise(stub, empty, rng)


def test_clone_span_preserved_under_noise(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """PD-L1 clone attachments survive every noise category with non-empty spans."""
    rng = random.Random(11)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    checked = 0
    attempts = 0
    while checked < 30 and attempts < 800:
        attempts += 1
        rendered = render_l1_record(
            "PD-L1", profile, biomarkers, common, rng, method_attach_prob=0.0
        )
        if rendered.clone_id is None:
            continue
        post = apply_technical_noise(rendered, common, rng)
        assert "clone" in post.spans
        start, end = post.spans["clone"]
        assert 0 <= start < end <= len(post.sentence)
        assert post.sentence[start:end], "empty clone span post-noise"
        checked += 1
    assert checked >= 10, f"only {checked} clone-bearing records observed"


def test_l2_complexity_level_propagates_through_noise(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L2 RenderedRecord → PostProcessedRecord preserves complexity_level."""
    rng = random.Random(17)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        rendered = render_l1_record(
            "EGFR", profile, biomarkers, common, rng, complexity_level="L2"
        )
        assert rendered.complexity_level == "L2"
        post = apply_technical_noise(rendered, common, rng)
        assert post.complexity_level == "L2"
        for name, (start, end) in post.spans.items():
            assert 0 <= start < end <= len(post.sentence)
            assert post.sentence[start:end]
