"""Probability-weighted variation config.

Loads ``common_variations.json`` and ``biomarkers.json`` into typed Pydantic
models, enforces the cross-config invariants in ``validator.py``, and exposes
helpers the generator uses for weighted sampling.
"""
from bioassert.config.loader import (
    BiomarkerConfig,
    CommonConfig,
    load_biomarkers,
    load_common,
    load_configs,
)
from bioassert.config.schema import (
    Biomarker,
    CloneAttribution,
    NegativeForms,
    PostProcessTransformations,
    PreferredMethods,
    RenderConstraints,
    StatusDistribution,
    TransformationCategory,
    Variant,
    WeightedVariations,
)

__all__ = [
    "Biomarker",
    "BiomarkerConfig",
    "CloneAttribution",
    "CommonConfig",
    "NegativeForms",
    "PostProcessTransformations",
    "PreferredMethods",
    "RenderConstraints",
    "StatusDistribution",
    "TransformationCategory",
    "Variant",
    "WeightedVariations",
    "load_biomarkers",
    "load_common",
    "load_configs",
]
