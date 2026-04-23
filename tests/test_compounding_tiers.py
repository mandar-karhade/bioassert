"""Compounding complexity tiers for L3+ records (Sub-phase 3.7).

The ``compounding_tier`` is orthogonal to ``complexity_level`` (L3/L3S/L4/
L4S/L5). It drives how many genes participate in one compound record:

- ``"low"``: < 3 genes (i.e., exactly 2 — the minimum compound size)
- ``"high"``: 3+ genes (3-8, clamped to pool size; allows cross-class
  pool mixing and multi-sentence emergence through ``sep=". "`` L4 frames)

Per user feedback 2026-04-23 the tier split is a simple 50/50 binary;
medium was collapsed into high. The corpus script takes ``--compound-low``
and ``--compound-high`` whose fractions must sum to 1.0.
"""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.patient_sampler import PatientProfile
from bioassert.generator.renderer import (
    COMPOUND_HIGH_MAX,
    COMPOUND_HIGH_MIN,
    COMPOUND_LOW_N,
    RenderedRecord,
    n_range_for_tier,
    render_l3_record,
    render_l4_record,
    render_l5_record,
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

MIXED_POOL: tuple[str, ...] = MUTATION_BIOMARKERS + EXPRESSION_BIOMARKERS


def test_compounding_tier_defaults_to_low() -> None:
    """RenderedRecord gains a compounding_tier field with default 'low' so
    existing construction sites (L1/L2 single-fact records, plus L3+ call
    sites that omit the argument) remain valid without adaptation.
    """
    rec = RenderedRecord(
        sentence="EGFR positive.",
        assertions=(),
        frame_template="{gene} {status}.",
    )
    assert rec.compounding_tier == "low"


def test_tier_constants_are_ordered() -> None:
    assert COMPOUND_LOW_N == 2
    assert COMPOUND_HIGH_MIN == 3
    assert COMPOUND_HIGH_MAX == 8
    assert COMPOUND_LOW_N < COMPOUND_HIGH_MIN
    assert COMPOUND_HIGH_MIN < COMPOUND_HIGH_MAX


def test_n_range_for_tier() -> None:
    assert n_range_for_tier("low") == (COMPOUND_LOW_N, COMPOUND_LOW_N)
    assert n_range_for_tier("high") == (COMPOUND_HIGH_MIN, COMPOUND_HIGH_MAX)


def test_n_range_for_tier_rejects_unknown() -> None:
    with pytest.raises(Exception):
        n_range_for_tier("bogus")


def test_n_range_for_tier_rejects_medium() -> None:
    """Sub-phase 3.7 per user feedback 2026-04-23: medium was removed."""
    with pytest.raises(Exception):
        n_range_for_tier("medium")


# --- L3 ----------------------------------------------------------------


def test_l3_low_tier_renders_two_genes(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(701)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(50):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
            compounding_tier="low",
        )
        assert rec.compounding_tier == "low"
        assert len(rec.assertions) == 2


def test_l3_high_tier_renders_three_to_eight_genes(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(703)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    observed_sizes: set[int] = set()
    for _ in range(300):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
            compounding_tier="high",
        )
        assert rec.compounding_tier == "high"
        assert COMPOUND_HIGH_MIN <= len(rec.assertions) <= COMPOUND_HIGH_MAX
        observed_sizes.add(len(rec.assertions))
    # Every size in [3, 8] should surface over 300 draws (pool=9 covers all).
    assert observed_sizes == set(range(COMPOUND_HIGH_MIN, COMPOUND_HIGH_MAX + 1))


def test_l3s_high_tier_renders_three_to_eight_genes(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(707)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
            compounding_tier="high",
        )
        assert rec.compounding_tier == "high"
        assert COMPOUND_HIGH_MIN <= len(rec.assertions) <= COMPOUND_HIGH_MAX


def test_l3_default_tier_is_low(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Calling render_l3_record without a tier defaults to 'low' (2 genes)."""
    rng = random.Random(709)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(50):
        rec = render_l3_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        assert rec.compounding_tier == "low"
        assert len(rec.assertions) == 2


# --- L4 ----------------------------------------------------------------


def test_l4_low_tier_renders_two_genes(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(711)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(50):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
            compounding_tier="low",
        )
        assert rec.compounding_tier == "low"
        assert len(rec.assertions) == 2


def test_l4_high_tier_renders_three_to_eight_genes(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(713)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    observed: set[int] = set()
    for _ in range(300):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
            compounding_tier="high",
        )
        assert COMPOUND_HIGH_MIN <= len(rec.assertions) <= COMPOUND_HIGH_MAX
        observed.add(len(rec.assertions))
    assert observed == set(range(COMPOUND_HIGH_MIN, COMPOUND_HIGH_MAX + 1))


def test_l4s_high_tier_renders_three_to_eight_genes(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(719)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l4_record(
            list(MUTATION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
            complexity_level="L4S",
            compounding_tier="high",
        )
        assert COMPOUND_HIGH_MIN <= len(rec.assertions) <= COMPOUND_HIGH_MAX


# --- L5 ----------------------------------------------------------------


def test_l5_low_tier_renders_small_wide_scope(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Low tier L5: enumerated wide scope has exactly 2 wide facts.
    Panel-wide frames always produce 1 fact regardless of tier.
    """
    rng = random.Random(721)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
            compounding_tier="low",
        )
        assert rec.compounding_tier == "low"
        wide = [a for a in rec.assertions if a.polarity_scope == "negation_wide"]
        if wide:
            assert len(wide) == 2


def test_l5_high_tier_renders_expanded_wide_scope(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """High tier L5: enumerated wide scope has >= 5 wide facts (pool size
    permitting). Panel-wide frames ignore tier (1 fact always).
    """
    rng = random.Random(723)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    seen_enumerated = 0
    for _ in range(200):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
            compounding_tier="high",
        )
        assert rec.compounding_tier == "high"
        wide = [a for a in rec.assertions if a.polarity_scope == "negation_wide"]
        if not wide:
            continue  # panel_wide
        seen_enumerated += 1
        assert COMPOUND_HIGH_MIN <= len(wide) <= COMPOUND_HIGH_MAX
    assert seen_enumerated > 0


# --- cross-class pool (high tier) --------------------------------------


def test_high_tier_cross_class_pool_l3(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """High tier accepts a mixed pool (mutation + expression) and emits a
    record whose assertions may include biomarkers from both classes.
    The renderer does not enforce pool homogeneity.
    """
    rng = random.Random(731)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    seen_mixed = False
    for _ in range(50):
        rec = render_l3_record(
            list(MIXED_POOL),
            profile,
            biomarkers,
            common,
            rng,
            complexity_level="L3S",
            compounding_tier="high",
        )
        genes = {a.gene for a in rec.assertions}
        has_mut = genes & set(MUTATION_BIOMARKERS)
        has_exp = genes & set(EXPRESSION_BIOMARKERS)
        if has_mut and has_exp:
            seen_mixed = True
            break
    assert seen_mixed, "no mixed-class high-tier record in 50 draws"


# --- tier validation ---------------------------------------------------


def test_renderer_rejects_unknown_tier(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(741)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    with pytest.raises(Exception):
        render_l3_record(
            list(MUTATION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
            compounding_tier="super-duper",
        )


def test_renderer_rejects_medium_tier(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Medium was collapsed into high per user feedback 2026-04-23."""
    rng = random.Random(742)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    with pytest.raises(Exception):
        render_l3_record(
            list(MUTATION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
            compounding_tier="medium",
        )


def test_high_tier_tolerates_small_pool(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """High tier against a 3-gene pool caps at the pool size rather than
    erroring — the caller is not forced to have >=5 biomarkers available.
    """
    rng = random.Random(743)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    small_pool = ["EGFR", "KRAS", "BRAF"]
    for _ in range(30):
        rec = render_l3_record(
            small_pool,
            profile,
            biomarkers,
            common,
            rng,
            compounding_tier="high",
        )
        assert 3 == len(rec.assertions), (
            f"expected all 3 genes, got {len(rec.assertions)}"
        )
