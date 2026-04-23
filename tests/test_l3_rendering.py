"""L3 shared-status rendering: N genes, one shared status, one frame.

Sub-phase 3.2 per docs/PHASE3_PLAN.md. Every listed gene becomes its own
AssertionFact with its own gene span; all facts share the same status span
and status value. L3 records carry no variant descriptors (deferred to L4/L5)
and keep the biomarker class homogeneous (mutation-class OR expression-class,
no mixing).
"""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.renderer import (
    RenderedRecord,
    _coordinate_gene_list,
    _L3_NAME_FORM_BLOCKLIST,
    render_l3_record,
)

MUTATION_BIOMARKERS: tuple[str, ...] = (
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

EXPRESSION_BIOMARKERS: tuple[str, ...] = ("PD-L1", "TMB")


def test_coordinate_gene_list_two_genes() -> None:
    rng = random.Random(0)
    joined, positions = _coordinate_gene_list(["EGFR", "KRAS"], rng)
    assert joined == "EGFR and KRAS"
    assert positions == [(0, 4), (9, 13)]
    assert joined[0:4] == "EGFR"
    assert joined[9:13] == "KRAS"


def test_coordinate_gene_list_three_genes_no_oxford() -> None:
    rng = random.Random(0)
    joined, positions = _coordinate_gene_list(
        ["EGFR", "KRAS", "BRAF"], rng, oxford_prob=0.0
    )
    assert joined == "EGFR, KRAS and BRAF"
    for (s, e), expected in zip(positions, ["EGFR", "KRAS", "BRAF"]):
        assert joined[s:e] == expected


def test_coordinate_gene_list_three_genes_with_oxford() -> None:
    rng = random.Random(0)
    joined, positions = _coordinate_gene_list(
        ["EGFR", "KRAS", "BRAF"], rng, oxford_prob=1.0
    )
    assert joined == "EGFR, KRAS, and BRAF"
    for (s, e), expected in zip(positions, ["EGFR", "KRAS", "BRAF"]):
        assert joined[s:e] == expected


def test_coordinate_gene_list_four_genes_with_oxford() -> None:
    rng = random.Random(1)
    joined, positions = _coordinate_gene_list(
        ["EGFR", "KRAS", "BRAF", "ALK"], rng, oxford_prob=1.0
    )
    assert joined == "EGFR, KRAS, BRAF, and ALK"
    for (s, e), expected in zip(positions, ["EGFR", "KRAS", "BRAF", "ALK"]):
        assert joined[s:e] == expected


def test_l3_record_has_n_assertions(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(42)
    for _ in range(200):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        assert isinstance(rec, RenderedRecord)
        n = len(rec.assertions)
        assert 2 <= n <= 4, f"expected N in [2,4], got {n}"
        assert rec.complexity_level == "L3"


def test_l3_all_assertions_share_status(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(7)
    for _ in range(200):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        statuses = {a.status for a in rec.assertions}
        assert len(statuses) == 1, f"status diverged across facts: {statuses}"


def test_l3_all_genes_are_literal_substrings(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(3)
    for _ in range(300):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        for fact in rec.assertions:
            start, end = fact.spans["gene"]
            assert 0 <= start < end <= len(rec.sentence)
            surface = rec.sentence[start:end]
            known = set(
                biomarkers.get(fact.gene).name_forms.realizations.values()
            )
            assert surface in known, (
                f"{fact.gene} surface {surface!r} not a configured name form"
            )


def test_l3_no_variants(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L3 is bare gene names — variant_id and negative_form_id always None."""
    rng = random.Random(9)
    for _ in range(200):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        for fact in rec.assertions:
            assert fact.variant_id is None
            assert fact.negative_form_id is None
            assert "variant" not in fact.spans


def test_l3_shared_status_span_identical_across_facts(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(11)
    for _ in range(200):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        status_spans = {fact.spans["status"] for fact in rec.assertions}
        assert len(status_spans) == 1, (
            f"status span diverged across facts: {status_spans}"
        )


def test_l3_spans_non_overlapping_within_record(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Gene spans must be disjoint; status span must not overlap any gene span."""
    rng = random.Random(13)
    for _ in range(200):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        gene_spans = sorted(fact.spans["gene"] for fact in rec.assertions)
        for i in range(1, len(gene_spans)):
            assert gene_spans[i - 1][1] <= gene_spans[i][0], (
                f"overlapping gene spans in {rec.sentence!r}: {gene_spans}"
            )
        status_span = rec.assertions[0].spans["status"]
        for gs, ge in gene_spans:
            assert status_span[1] <= gs or ge <= status_span[0], (
                f"status span {status_span} overlaps gene span ({gs},{ge}) "
                f"in {rec.sentence!r}"
            )


def test_l3_gene_order_in_sentence_matches_fact_order(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Fact ordering follows the sentence's left-to-right gene order."""
    rng = random.Random(17)
    for _ in range(100):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        starts = [fact.spans["gene"][0] for fact in rec.assertions]
        assert starts == sorted(starts), (
            f"fact order does not match sentence order: {starts}"
        )


def test_l3_genes_are_distinct_within_record(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(19)
    for _ in range(200):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        genes = [fact.gene for fact in rec.assertions]
        assert len(genes) == len(set(genes)), f"duplicate genes: {genes}"


def test_l3_expression_pool_only_supports_n2(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Expression class has only 2 biomarkers (PD-L1, TMB); N capped at 2."""
    rng = random.Random(23)
    for _ in range(50):
        rec = render_l3_record(
            list(EXPRESSION_BIOMARKERS), biomarkers, common, rng
        )
        assert len(rec.assertions) == 2
        assert {f.gene for f in rec.assertions} == set(EXPRESSION_BIOMARKERS)


def test_l3_frame_template_recorded(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Record stores the L3 frame template; must contain both placeholders."""
    rng = random.Random(29)
    for _ in range(100):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        assert "{gene_list}" in rec.frame_template
        assert "{status}" in rec.frame_template


def test_l3_rejects_pool_smaller_than_two(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(0)
    with pytest.raises(Exception):
        render_l3_record(["EGFR"], biomarkers, common, rng)


def test_l3_gene_surfaces_are_bare(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L3 gene surfaces must not embed alteration-type or status keywords.

    Forms like ``"ALK-positive"``, ``"KRAS mutation"``, ``"NTRK fusion"``
    would double-mark the assertion because L3 frames carry the status in
    its own slot. ``"tumor mutational burden"`` must still be permitted
    because ``mutational`` is a different word than ``mutation``.
    """
    rng = random.Random(31)
    for _ in range(500):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        for fact in rec.assertions:
            s, e = fact.spans["gene"]
            surface = rec.sentence[s:e]
            assert _L3_NAME_FORM_BLOCKLIST.search(surface) is None, (
                f"{fact.gene} surface {surface!r} embeds a blocklisted "
                f"keyword in sentence {rec.sentence!r}"
            )
    rng = random.Random(37)
    for _ in range(200):
        rec = render_l3_record(
            list(EXPRESSION_BIOMARKERS), biomarkers, common, rng
        )
        for fact in rec.assertions:
            s, e = fact.spans["gene"]
            surface = rec.sentence[s:e]
            assert _L3_NAME_FORM_BLOCKLIST.search(surface) is None, (
                f"{fact.gene} surface {surface!r} embeds a blocklisted "
                f"keyword in sentence {rec.sentence!r}"
            )
