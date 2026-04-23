"""L4 heterogeneous-prose rendering: N genes, each with its own status.

Sub-phase 3.4 per docs/PHASE3_PLAN.md. Same list structure as L3 but
every gene's status is drawn independently from its own population.
Each :class:`AssertionFact` carries its own gene span AND its own
status span (distributive — like L3 shorthand, unlike L3 prose). L4 is
"prose" meaning it uses the formal positive_phrases / negation_phrases
status vocabulary, not shorthand tokens.
"""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.config.schema import WeightedVariations
from bioassert.generator.renderer import (
    RenderedRecord,
    _L3_NAME_FORM_BLOCKLIST,
    _L4_FRAMES,
    render_l4_record,
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


def _formal_status_vocab(common: CommonConfig) -> set[str]:
    """All non-gene-embedding realizations across the four formal status
    categories. Used to sanity-check L4 status surfaces (surface text
    will often have prefix substrings due to expand_placeholders but each
    piece still comes from this vocab space).
    """
    out: set[str] = set()
    for name in (
        "positive_phrases",
        "negation_phrases",
        "equivocal_phrases",
        "not_tested_phrases",
    ):
        cat = common.categories[name]
        assert isinstance(cat, WeightedVariations)
        for vid, surf in cat.realizations.items():
            if "{gene}" in surf:
                continue
            out.add(surf)
    return out


def test_l4_frames_nonempty() -> None:
    assert len(_L4_FRAMES) >= 3
    for f in _L4_FRAMES:
        assert "inner" in f
        assert "sep" in f


def test_l4_record_has_n_assertions(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(201)
    for _ in range(200):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        assert isinstance(rec, RenderedRecord)
        n = len(rec.assertions)
        assert 2 <= n <= 4
        assert rec.complexity_level == "L4"


def test_l4_per_gene_statuses_can_diverge(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L4 statuses are drawn independently per gene. With realistic
    prevalences we must see at least some records where facts don't all
    share the same status value.
    """
    rng = random.Random(211)
    seen_divergent = False
    for _ in range(500):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        statuses = {a.status for a in rec.assertions}
        if len(statuses) > 1:
            seen_divergent = True
            break
    assert seen_divergent, "no L4 sample drew divergent per-gene statuses"


def test_l4_every_fact_has_gene_and_status_span(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(223)
    for _ in range(200):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        for fact in rec.assertions:
            assert "gene" in fact.spans
            assert "status" in fact.spans
            gs, ge = fact.spans["gene"]
            ss, se = fact.spans["status"]
            assert 0 <= gs < ge <= len(rec.sentence)
            assert 0 <= ss < se <= len(rec.sentence)


def test_l4_spans_literal_substrings(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(227)
    for _ in range(300):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        for fact in rec.assertions:
            gs, ge = fact.spans["gene"]
            known = set(
                biomarkers.get(fact.gene).name_forms.realizations.values()
            )
            assert rec.sentence[gs:ge] in known


def test_l4_gene_surfaces_are_bare(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(229)
    for _ in range(400):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        for fact in rec.assertions:
            s, e = fact.spans["gene"]
            surface = rec.sentence[s:e]
            assert _L3_NAME_FORM_BLOCKLIST.search(surface) is None, (
                f"{fact.gene} surface {surface!r} has blocklisted keyword"
            )


def test_l4_spans_non_overlapping_within_record(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(233)
    for _ in range(200):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        all_spans: list[tuple[int, int, str]] = []
        for fact in rec.assertions:
            for name, (s, e) in fact.spans.items():
                all_spans.append((s, e, f"{fact.gene}.{name}"))
        all_spans.sort()
        for i in range(1, len(all_spans)):
            prev = all_spans[i - 1]
            cur = all_spans[i]
            assert prev[1] <= cur[0], (
                f"overlap in {rec.sentence!r}: {prev} vs {cur}"
            )


def test_l4_gene_order_in_sentence_matches_fact_order(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(239)
    for _ in range(100):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        starts = [fact.spans["gene"][0] for fact in rec.assertions]
        assert starts == sorted(starts)


def test_l4_no_variants(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(241)
    for _ in range(200):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        for fact in rec.assertions:
            assert fact.variant_id is None
            assert fact.negative_form_id is None
            assert "variant" not in fact.spans


def test_l4_genes_distinct_within_record(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(251)
    for _ in range(200):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        genes = [fact.gene for fact in rec.assertions]
        assert len(genes) == len(set(genes))


def test_l4_expression_pool_only_supports_n2(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(257)
    for _ in range(50):
        rec = render_l4_record(
            list(EXPRESSION_BIOMARKERS), biomarkers, common, rng
        )
        assert len(rec.assertions) == 2
        assert {f.gene for f in rec.assertions} == set(EXPRESSION_BIOMARKERS)


def test_l4_rejects_pool_smaller_than_two(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(0)
    with pytest.raises(Exception):
        render_l4_record(["EGFR"], biomarkers, common, rng)


def test_l4_status_spans_are_distinct_per_fact(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L4 renders status surface per gene, so each fact gets its own
    status span — distinct start/end pairs, unlike L3 prose's shared
    span.
    """
    rng = random.Random(263)
    for _ in range(200):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        spans = [fact.spans["status"] for fact in rec.assertions]
        assert len(set(spans)) == len(spans), (
            f"status spans collapsed in {rec.sentence!r}: {spans}"
        )


def test_l4_frame_template_recorded(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Stores some non-empty frame template string for downstream analysis."""
    rng = random.Random(271)
    for _ in range(50):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS), biomarkers, common, rng
        )
        assert rec.frame_template
        assert "{gene}" in rec.frame_template
        assert "{status}" in rec.frame_template
