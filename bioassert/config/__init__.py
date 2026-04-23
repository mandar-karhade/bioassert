"""Probability-weighted variation config (Phase 2a).

Loads `common_variations.json` and `biomarkers.json` into typed Pydantic
models, enforces the 12 validation rules from CONFIG_ARCHITECTURE.md Section 7,
and exposes helpers the generator uses for weighted sampling.
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
    PopulationStatus,
    PostProcessTransformations,
    PreferredMethods,
    RenderConstraints,
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
    "PopulationStatus",
    "PostProcessTransformations",
    "PreferredMethods",
    "RenderConstraints",
    "TransformationCategory",
    "Variant",
    "WeightedVariations",
    "load_biomarkers",
    "load_common",
    "load_configs",
]
