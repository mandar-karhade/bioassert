"""L6 temporal/certainty qualification (Sub-phase 3.8).

L6 records attach optional qualifier fields to each ``AssertionFact``:

- ``temporal``: when the observation applies (``previously``, ``currently``,
  ``at diagnosis``, ``at relapse``, ``post-TKI``)
- ``certainty``: hedging confidence (``suspected``, ``confirmed``,
  ``probable``, ``rule-out``, ``pending``)

Three structural shapes:

- ``"temporal"`` — 2 facts, same gene, 2 timepoints, divergent statuses
  (both ``temporal`` set, ``certainty`` None)
- ``"certainty"`` — 1 fact, hedged status (``certainty`` set,
  ``temporal`` None)
- ``"combined"`` — 2 facts, same gene, 2 timepoints, each carrying a
  ``certainty`` value (both ``temporal`` AND ``certainty`` set)

L6 records are single-gene for Sub-phase 3.8. Multi-gene L6 is deferred.
``compounding_tier`` is always ``"low"`` because the tier knob doesn't
apply to single-gene records.
"""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.renderer import (
    AssertionFact,
    L6_CERTAINTY_VOCAB,
    L6_SHAPES,
    L6_TEMPORAL_VOCAB,
    RenderError,
    RenderedRecord,
    render_l6_record,
)

GENES: tuple[str, ...] = (
    "EGFR",
    "KRAS",
    "BRAF",
    "ALK",
    "ROS1",
    "RET",
    "NTRK",
    "MET",
    "ERBB2",
)


# --- schema ------------------------------------------------------------


def test_assertion_fact_temporal_certainty_default_none() -> None:
    """Existing L1-L5 construction sites must remain valid without passing
    the new fields. Both default to None.
    """
    fact = AssertionFact(gene="EGFR", status="positive", spans={})
    assert fact.temporal is None
    assert fact.certainty is None


def test_vocab_constants_exact() -> None:
    """Vocabulary tuples match the spec in PHASE3_PLAN.md section 3.8."""
    assert L6_TEMPORAL_VOCAB == (
        "previously",
        "currently",
        "at diagnosis",
        "at relapse",
        "post-TKI",
    )
    assert L6_CERTAINTY_VOCAB == (
        "suspected",
        "confirmed",
        "probable",
        "rule-out",
        "pending",
    )


def test_shape_constants_exact() -> None:
    assert L6_SHAPES == ("temporal", "certainty", "combined")


# --- render_l6_record: temporal shape ----------------------------------


def test_l6_temporal_emits_two_facts_same_gene(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(801)
    for _ in range(100):
        rec = render_l6_record(
            list(GENES), biomarkers, common, rng, shape="temporal"
        )
        assert rec.complexity_level == "L6"
        assert rec.compounding_tier == "low"
        assert len(rec.assertions) == 2
        a, b = rec.assertions
        assert a.gene == b.gene, "temporal contrast should be same gene"
        assert a.temporal is not None and a.temporal in L6_TEMPORAL_VOCAB
        assert b.temporal is not None and b.temporal in L6_TEMPORAL_VOCAB
        assert a.temporal != b.temporal, "two timepoints should differ"
        assert a.certainty is None
        assert b.certainty is None
        # Temporal contrast is clinically interesting when statuses differ.
        assert a.status != b.status, "temporal contrast should differ in status"


# --- render_l6_record: certainty shape ---------------------------------


def test_l6_certainty_emits_one_fact_with_certainty(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(803)
    for _ in range(100):
        rec = render_l6_record(
            list(GENES), biomarkers, common, rng, shape="certainty"
        )
        assert rec.complexity_level == "L6"
        assert rec.compounding_tier == "low"
        assert len(rec.assertions) == 1
        f = rec.assertions[0]
        assert f.certainty is not None and f.certainty in L6_CERTAINTY_VOCAB
        assert f.temporal is None


# --- render_l6_record: combined shape ----------------------------------


def test_l6_combined_emits_two_facts_with_both_qualifiers(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(805)
    for _ in range(100):
        rec = render_l6_record(
            list(GENES), biomarkers, common, rng, shape="combined"
        )
        assert rec.complexity_level == "L6"
        assert len(rec.assertions) == 2
        a, b = rec.assertions
        assert a.gene == b.gene
        for f in (a, b):
            assert f.temporal is not None and f.temporal in L6_TEMPORAL_VOCAB
            assert f.certainty is not None and f.certainty in L6_CERTAINTY_VOCAB
        assert a.temporal != b.temporal


# --- shape dispatch ----------------------------------------------------


def test_l6_random_shape_covers_all_three(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """With ``shape=None`` the renderer picks one of the three shapes
    uniformly. Over enough draws every shape should surface.
    """
    rng = random.Random(807)
    shapes_seen: set[str] = set()
    for _ in range(300):
        rec = render_l6_record(list(GENES), biomarkers, common, rng)
        if len(rec.assertions) == 1:
            shapes_seen.add("certainty")
        elif any(f.certainty for f in rec.assertions):
            shapes_seen.add("combined")
        else:
            shapes_seen.add("temporal")
    assert shapes_seen == {"temporal", "certainty", "combined"}


def test_l6_rejects_unknown_shape(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(809)
    with pytest.raises(RenderError):
        render_l6_record(
            list(GENES), biomarkers, common, rng, shape="bogus"
        )


# --- span invariants ---------------------------------------------------


def test_l6_every_fact_has_gene_and_status_span(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(811)
    for _ in range(300):
        rec = render_l6_record(list(GENES), biomarkers, common, rng)
        for fact in rec.assertions:
            assert "gene" in fact.spans
            assert "status" in fact.spans
            gs, ge = fact.spans["gene"]
            ss, se = fact.spans["status"]
            assert 0 <= gs < ge <= len(rec.sentence)
            assert 0 <= ss < se <= len(rec.sentence)


def test_l6_gene_spans_are_literal_name_forms(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(813)
    for _ in range(300):
        rec = render_l6_record(list(GENES), biomarkers, common, rng)
        for fact in rec.assertions:
            gs, ge = fact.spans["gene"]
            surface = rec.sentence[gs:ge]
            known = set(
                biomarkers.get(fact.gene).name_forms.realizations.values()
            )
            assert surface in known, (
                f"{fact.gene} surface {surface!r} not in known name_forms"
            )


def test_l6_sentence_contains_temporal_markers(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """When ``temporal`` is set on a fact the corresponding marker must be a
    literal substring of the sentence — otherwise the label is unverifiable.
    """
    rng = random.Random(815)
    for _ in range(200):
        rec = render_l6_record(list(GENES), biomarkers, common, rng)
        for fact in rec.assertions:
            if fact.temporal is not None:
                assert fact.temporal in rec.sentence


def test_l6_sentence_contains_certainty_markers(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(817)
    for _ in range(200):
        rec = render_l6_record(list(GENES), biomarkers, common, rng)
        for fact in rec.assertions:
            if fact.certainty is not None:
                assert fact.certainty in rec.sentence


def test_l6_spans_non_overlapping_within_record(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Within a record, each distinct (start, end) position participates in
    at most one labeled token — but two facts may *share* the same surface
    (single gene mention carrying two temporal-contrast facts), in which
    case both point at the same span position. Dedupe by position before
    checking.
    """
    rng = random.Random(819)
    for _ in range(200):
        rec = render_l6_record(list(GENES), biomarkers, common, rng)
        unique_spans: set[tuple[int, int]] = set()
        for fact in rec.assertions:
            for s, e in fact.spans.values():
                unique_spans.add((s, e))
        sorted_spans = sorted(unique_spans)
        for i in range(1, len(sorted_spans)):
            prev_s, prev_e = sorted_spans[i - 1]
            cur_s, cur_e = sorted_spans[i]
            assert prev_e <= cur_s, (
                f"overlap in {rec.sentence!r}: "
                f"({prev_s},{prev_e}) vs ({cur_s},{cur_e})"
            )


# --- vocabulary coverage (smoke) --------------------------------------


def test_l6_temporal_vocab_is_exercised(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Over many draws every temporal term should surface."""
    rng = random.Random(821)
    seen: set[str] = set()
    for _ in range(800):
        rec = render_l6_record(
            list(GENES), biomarkers, common, rng, shape="temporal"
        )
        for fact in rec.assertions:
            if fact.temporal:
                seen.add(fact.temporal)
    assert seen == set(L6_TEMPORAL_VOCAB)


def test_l6_certainty_vocab_is_exercised(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(823)
    seen: set[str] = set()
    for _ in range(800):
        rec = render_l6_record(
            list(GENES), biomarkers, common, rng, shape="certainty"
        )
        for fact in rec.assertions:
            if fact.certainty:
                seen.add(fact.certainty)
    assert seen == set(L6_CERTAINTY_VOCAB)
