"""L4 shorthand/tabular rendering: N genes, each with its own status,
rendered in a compact `gene <inner> status<sep>gene <inner> status`
pair list using shorthand status vocab.

Sub-phase 3.5 per docs/PHASE3_PLAN.md. Parallels Sub-phase 3.3 L3S but
statuses are drawn independently per gene (like L4 prose, unlike L3S).
"""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.config.schema import WeightedVariations
from bioassert.generator.renderer import (
    RenderedRecord,
    _L3_NAME_FORM_BLOCKLIST,
    _L4_SHORTHAND_FRAMES,
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


def _shorthand_vocab(common: CommonConfig) -> tuple[set[str], set[str]]:
    pos = common.categories["positive_shorthand"]
    neg = common.categories["negation_shorthand"]
    assert isinstance(pos, WeightedVariations)
    assert isinstance(neg, WeightedVariations)
    return set(pos.realizations.values()), set(neg.realizations.values())


def test_l4_shorthand_frames_nonempty() -> None:
    assert len(_L4_SHORTHAND_FRAMES) >= 3
    for f in _L4_SHORTHAND_FRAMES:
        assert f["inner"]
        assert f["sep"]


def test_l4_shorthand_complexity_level(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(301)
    for _ in range(100):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L4S",
        )
        assert rec.complexity_level == "L4S"
        assert 2 <= len(rec.assertions) <= 4


def test_l4_shorthand_per_gene_statuses_can_diverge(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(307)
    seen = False
    for _ in range(600):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L4S",
        )
        statuses = {a.status for a in rec.assertions}
        if len(statuses) > 1:
            seen = True
            break
    assert seen, "no L4S sample drew divergent per-gene statuses"


def test_l4_shorthand_uses_shorthand_vocab_for_pos_neg(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Each fact whose status is positive or negative must render the
    status surface from the shorthand vocab. Equivocal / not_tested
    fall back to formal phrases (no shorthand exists).
    """
    pos_vocab, neg_vocab = _shorthand_vocab(common)
    rng = random.Random(311)
    checked = 0
    for _ in range(500):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L4S",
        )
        for fact in rec.assertions:
            if fact.status not in ("positive", "negative"):
                continue
            checked += 1
            expected = pos_vocab if fact.status == "positive" else neg_vocab
            s, e = fact.spans["status"]
            surface = rec.sentence[s:e]
            assert surface in expected, (
                f"{fact.status} surface {surface!r} not in shorthand vocab "
                f"for sentence {rec.sentence!r}"
            )
    assert checked > 0, "no positive/negative L4S fact sampled"


def test_l4_shorthand_gene_surfaces_are_bare(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(313)
    for _ in range(300):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L4S",
        )
        for fact in rec.assertions:
            s, e = fact.spans["gene"]
            surface = rec.sentence[s:e]
            assert _L3_NAME_FORM_BLOCKLIST.search(surface) is None, (
                f"{fact.gene} surface {surface!r} embeds blocklisted keyword"
            )


def test_l4_shorthand_spans_literal(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(317)
    for _ in range(300):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L4S",
        )
        for fact in rec.assertions:
            for name, (s, e) in fact.spans.items():
                assert rec.sentence[s:e], f"{name} span empty slice"


def test_l4_shorthand_spans_non_overlapping(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(319)
    for _ in range(200):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L4S",
        )
        all_spans: list[tuple[int, int, str]] = []
        for fact in rec.assertions:
            for name, (s, e) in fact.spans.items():
                all_spans.append((s, e, f"{fact.gene}.{name}"))
        all_spans.sort()
        for i in range(1, len(all_spans)):
            assert all_spans[i - 1][1] <= all_spans[i][0], (
                f"overlap in {rec.sentence!r}: "
                f"{all_spans[i - 1]} vs {all_spans[i]}"
            )


def test_l4_shorthand_no_variants(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(331)
    for _ in range(200):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L4S",
        )
        for fact in rec.assertions:
            assert fact.variant_id is None
            assert fact.negative_form_id is None
            assert "variant" not in fact.spans


def test_l4_shorthand_expression_pool(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(337)
    for _ in range(50):
        rec = render_l4_record(
            list(EXPRESSION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L4S",
        )
        assert len(rec.assertions) == 2
        assert {f.gene for f in rec.assertions} == set(EXPRESSION_BIOMARKERS)


def test_l4_shorthand_genes_distinct_within_record(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(347)
    for _ in range(200):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L4S",
        )
        genes = [fact.gene for fact in rec.assertions]
        assert len(genes) == len(set(genes))


def test_l4_shorthand_rejects_invalid_complexity_level(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(0)
    with pytest.raises(Exception):
        render_l4_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3",
        )
