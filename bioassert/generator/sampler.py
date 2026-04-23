"""Deterministic weighted sampling over the Phase 2a config (§2.3–§8.4).

Wraps :class:`random.Random` (stdlib) with helpers that speak the shapes of
:mod:`bioassert.config.schema`. All entry points take an explicit ``rng``
so generation is reproducible from a single seed.
"""
from __future__ import annotations

import random
from typing import Optional

from bioassert.config.loader import BiomarkerConfig, CommonConfig
from bioassert.config.schema import (
    Biomarker,
    CloneAttribution,
    PopulationStatus,
    Variant,
    WeightedVariations,
)
from bioassert.generator.patient_sampler import (
    PatientProfile,
    resolve_status_distribution,
)

STATUS_NAMES: tuple[str, ...] = ("positive", "negative", "equivocal", "not_tested")


def _weighted_choice(
    rng: random.Random, weights: dict[str, float], subset: Optional[set[str]] = None
) -> str:
    if subset is not None:
        keys = [k for k in weights if k in subset]
        if not keys:
            raise ValueError(
                f"weighted_choice: no keys overlap subset {sorted(subset)}"
            )
    else:
        keys = list(weights)
    values = [weights[k] for k in keys]
    return rng.choices(keys, weights=values, k=1)[0]


def sample_variation(
    weighted: WeightedVariations | CloneAttribution,
    rng: random.Random,
    subset: Optional[set[str]] = None,
) -> str:
    """Sample a variation_id, optionally restricted to a subset."""
    return _weighted_choice(rng, weighted.variations, subset=subset)


def sample_status(
    distribution: PopulationStatus, rng: random.Random
) -> str:
    """Sample one of positive / negative / equivocal / not_tested."""
    weights = {
        "positive": distribution.positive,
        "negative": distribution.negative,
        "equivocal": distribution.equivocal,
        "not_tested": distribution.not_tested,
    }
    return _weighted_choice(rng, weights)


def sample_variant(
    biomarker: Biomarker, rng: random.Random
) -> tuple[str, Variant]:
    """Sample a variant_id from ``prevalence_within_biomarker`` weights."""
    prev = {
        vid: v.prevalence_within_biomarker for vid, v in biomarker.variants.items()
    }
    variant_id = _weighted_choice(rng, prev)
    return variant_id, biomarker.variants[variant_id]


def sample_biomarker_name_form(
    biomarker: Biomarker,
    rng: random.Random,
    variant: Optional[Variant] = None,
) -> str:
    """Sample a biomarker-level name_form ID.

    Honours ``variant.render_constraints.require_biomarker_name_forms`` when
    present — the subset restriction keeps ``unspecified_fusion`` variants from
    producing a sentence that loses the fusion signal
    (config_architecture.md §8.3).
    """
    subset: Optional[set[str]] = None
    if variant is not None and variant.render_constraints is not None:
        subset = set(variant.render_constraints.require_biomarker_name_forms)
    return sample_variation(biomarker.name_forms, rng, subset=subset)


def maybe_sample_clone(
    biomarker: Biomarker, rng: random.Random
) -> Optional[str]:
    """Roll the clone-attribution Bernoulli, sample a clone_id if attached.

    Independent of variant sampling — clones attach TO the TPS/CPS value, they
    aren't a sub-variant (§8.4). Returns ``None`` when no clone attaches or
    the biomarker has no ``clone_attribution`` block.
    """
    clone_block = biomarker.clone_attribution
    if clone_block is None:
        return None
    if rng.random() >= clone_block.attachment_probability:
        return None
    return sample_variation(clone_block, rng)


def sample_method(
    biomarker: Biomarker, common: CommonConfig, rng: random.Random
) -> str:
    """Sample a test_method variation_id.

    Uses ``biomarker.preferred_methods.variations`` when present, otherwise
    falls back to ``common.test_methods.variations``. Realization lookup is
    always via ``common.test_methods.realizations`` — biomarker-level configs
    supply weights only (§2.4).
    """
    if biomarker.preferred_methods is not None:
        weights = biomarker.preferred_methods.variations
        return _weighted_choice(rng, weights)
    methods = common.categories.get("test_methods")
    if not isinstance(methods, WeightedVariations):
        raise KeyError("common.test_methods is missing or wrong schema_type")
    return sample_variation(methods, rng)


def sample_measurement_value(
    measurement_range: tuple[float, float], rng: random.Random
) -> float:
    """Uniform sample in ``[lo, hi]`` for ``{value}`` placeholder expansion."""
    lo, hi = measurement_range
    return rng.uniform(lo, hi)


def resolve_population(
    biomarker: Biomarker, profile: PatientProfile
) -> PopulationStatus:
    dist, _key = resolve_status_distribution(biomarker, profile)
    return dist


__all__ = [
    "STATUS_NAMES",
    "maybe_sample_clone",
    "resolve_population",
    "sample_biomarker_name_form",
    "sample_measurement_value",
    "sample_method",
    "sample_status",
    "sample_variant",
    "sample_variation",
]
