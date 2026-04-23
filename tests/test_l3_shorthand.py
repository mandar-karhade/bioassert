"""L3 shorthand/tabular rendering: N genes, one shared status value,
rendered as a compact list of ``gene <inner> status`` pairs separated
by a delimiter (comma, newline, pipe, tab+newline, ...).

Sub-phase 3.3 per docs/PHASE3_PLAN.md. Parallels Sub-phase 3.2 L3 prose
but the status surface is drawn from ``positive_shorthand`` /
``negation_shorthand`` and repeated per gene. Each :class:`AssertionFact`
carries its own gene span AND its own status span (status surface is
distributed over the pair list, not shared like L3 prose).
"""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.config.schema import WeightedVariations
from bioassert.generator.renderer import (
    RenderedRecord,
    _L3_NAME_FORM_BLOCKLIST,
    _L3_SHORTHAND_FRAMES,
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


def _shorthand_vocab(common: CommonConfig) -> tuple[set[str], set[str]]:
    pos = common.categories["positive_shorthand"]
    neg = common.categories["negation_shorthand"]
    assert isinstance(pos, WeightedVariations)
    assert isinstance(neg, WeightedVariations)
    return set(pos.realizations.values()), set(neg.realizations.values())


def test_l3_shorthand_frames_nonempty() -> None:
    assert len(_L3_SHORTHAND_FRAMES) >= 3
    for f in _L3_SHORTHAND_FRAMES:
        assert "inner" in f and "sep" in f
        assert f["inner"]  # non-empty glue between gene and status
        assert f["sep"]  # non-empty glue between pairs


def test_l3_shorthand_complexity_level(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(101)
    for _ in range(100):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
        )
        assert rec.complexity_level == "L3S"
        assert isinstance(rec, RenderedRecord)
        assert 2 <= len(rec.assertions) <= 4


def test_l3_shorthand_all_share_status_value(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(103)
    for _ in range(200):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
        )
        statuses = {a.status for a in rec.assertions}
        assert len(statuses) == 1, f"status diverged: {statuses}"


def test_l3_shorthand_status_uses_shorthand_vocab_when_pos_or_neg(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Positive/negative L3S records must render the status from the
    shorthand vocabulary. Equivocal/not_tested have no shorthand, so
    they're allowed to draw from the formal phrase categories.
    """
    pos_vocab, neg_vocab = _shorthand_vocab(common)
    rng = random.Random(107)
    checked = 0
    for _ in range(500):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
        )
        status = rec.assertions[0].status
        if status not in ("positive", "negative"):
            continue
        checked += 1
        expected = pos_vocab if status == "positive" else neg_vocab
        for fact in rec.assertions:
            s, e = fact.spans["status"]
            surface = rec.sentence[s:e]
            assert surface in expected, (
                f"{status} shorthand surface {surface!r} not in vocab "
                f"{sorted(expected)} for sentence {rec.sentence!r}"
            )
    assert checked > 0, "no positive/negative L3S samples drawn"


def test_l3_shorthand_gene_spans_are_bare_surfaces(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(109)
    for _ in range(300):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
        )
        for fact in rec.assertions:
            s, e = fact.spans["gene"]
            surface = rec.sentence[s:e]
            assert _L3_NAME_FORM_BLOCKLIST.search(surface) is None, (
                f"{fact.gene} surface {surface!r} embeds blocklisted keyword"
            )
            known = set(
                biomarkers.get(fact.gene).name_forms.realizations.values()
            )
            assert surface in known


def test_l3_shorthand_spans_literal_substrings(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(113)
    for _ in range(300):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
        )
        for fact in rec.assertions:
            for name, (s, e) in fact.spans.items():
                assert 0 <= s < e <= len(rec.sentence), (
                    f"{name} span ({s},{e}) out of bounds"
                )
                assert rec.sentence[s:e], f"{name} span is empty slice"


def test_l3_shorthand_spans_non_overlapping(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(127)
    for _ in range(200):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
        )
        all_spans: list[tuple[int, int, str]] = []
        for fact in rec.assertions:
            for name, (s, e) in fact.spans.items():
                all_spans.append((s, e, f"{fact.gene}.{name}"))
        all_spans.sort()
        for i in range(1, len(all_spans)):
            prev_end = all_spans[i - 1][1]
            cur_start = all_spans[i][0]
            assert prev_end <= cur_start, (
                f"overlap in {rec.sentence!r}: "
                f"{all_spans[i - 1]} vs {all_spans[i]}"
            )


def test_l3_shorthand_gene_order_matches_fact_order(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(131)
    for _ in range(100):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
        )
        starts = [fact.spans["gene"][0] for fact in rec.assertions]
        assert starts == sorted(starts), f"fact order wrong: {starts}"


def test_l3_shorthand_no_variants(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(137)
    for _ in range(200):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
        )
        for fact in rec.assertions:
            assert fact.variant_id is None
            assert fact.negative_form_id is None
            assert "variant" not in fact.spans


def test_l3_shorthand_expression_pool(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(139)
    for _ in range(50):
        rec = render_l3_record(
            list(EXPRESSION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
        )
        assert len(rec.assertions) == 2
        assert {f.gene for f in rec.assertions} == set(EXPRESSION_BIOMARKERS)


def test_l3_shorthand_rejects_invalid_complexity_level(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(0)
    with pytest.raises(Exception):
        render_l3_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L2",
        )


def test_l3_shorthand_status_spans_distinct_per_fact(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Shorthand repeats the status surface per gene, so each fact has its
    own status span — distinct from the shared-span behavior of L3 prose.
    """
    rng = random.Random(149)
    seen_distinct = False
    for _ in range(100):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
        )
        status_spans = [fact.spans["status"] for fact in rec.assertions]
        if len(set(status_spans)) == len(status_spans):
            seen_distinct = True
            break
    assert seen_distinct, "expected at least one L3S record with distinct status spans"
