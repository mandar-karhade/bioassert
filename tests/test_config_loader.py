"""Loader + Pydantic-model invariants (config_architecture.md §7.1–7.4)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from bioassert.config import (
    BiomarkerConfig,
    CloneAttribution,
    CommonConfig,
    PostProcessTransformations,
    WeightedVariations,
    load_biomarkers,
    load_common,
    load_configs,
)
from bioassert.config.loader import _dispatch_common_category
from bioassert.config.schema import WEIGHT_TOLERANCE

ROOT = Path(__file__).resolve().parents[1]
PROJECT_DIR = ROOT / "projects" / "nsclc_adenocarcinoma"
COMMON_PATH = PROJECT_DIR / "configs" / "common_variations.json"
BIOMARKERS_PATH = PROJECT_DIR / "configs" / "biomarkers.json"


def test_configs_load_without_errors() -> None:
    common, biomarkers = load_configs(COMMON_PATH, BIOMARKERS_PATH)
    assert isinstance(common, CommonConfig)
    assert isinstance(biomarkers, BiomarkerConfig)
    assert biomarkers.biomarkers, "at least one biomarker must be loaded"


def test_common_categories_dispatched_by_schema_type(
    common: CommonConfig,
) -> None:
    types = {name: type(cat).__name__ for name, cat in common.categories.items()}
    assert types["assertion_verbs"] == "WeightedVariations"
    assert types["test_methods"] == "WeightedVariations"
    assert types["technical_noise"] == "PostProcessTransformations"


def test_every_weighted_distribution_sums_to_one(common: CommonConfig) -> None:
    """Structural invariant: sums within :data:`WEIGHT_TOLERANCE`."""
    for name, cat in common.categories.items():
        if isinstance(cat, (WeightedVariations, CloneAttribution)):
            total = sum(cat.variations.values())
            assert abs(total - 1.0) <= WEIGHT_TOLERANCE, (
                f"common.{name} weights sum to {total}"
            )


def test_post_process_sub_distributions_sum_to_one(common: CommonConfig) -> None:
    noise = common.categories["technical_noise"]
    assert isinstance(noise, PostProcessTransformations)
    for sub_name, category in noise.categories.items():
        total = sum(category.distribution.values())
        assert abs(total - 1.0) <= WEIGHT_TOLERANCE, (
            f"technical_noise.{sub_name} weights sum to {total}"
        )


def test_variant_prevalence_sums_to_one(
    biomarkers: BiomarkerConfig,
) -> None:
    for name, entry in biomarkers.biomarkers.items():
        total = sum(v.prevalence_within_biomarker for v in entry.variants.values())
        assert abs(total - 1.0) <= WEIGHT_TOLERANCE, (
            f"{name}: variant prevalence sums to {total}"
        )


def test_status_distribution_sums_to_one(
    biomarkers: BiomarkerConfig,
) -> None:
    for biomarker_name, entry in biomarkers.biomarkers.items():
        dist = entry.status_distribution
        total = dist.positive + dist.negative + dist.equivocal + dist.not_tested
        assert abs(total - 1.0) <= WEIGHT_TOLERANCE, (
            f"{biomarker_name} status_distribution sums to {total}"
        )


def test_variation_keys_match_realization_keys(
    common: CommonConfig,
) -> None:
    for name, cat in common.categories.items():
        if isinstance(cat, (WeightedVariations, CloneAttribution)):
            assert set(cat.variations.keys()) == set(cat.realizations.keys()), (
                f"common.{name} has mismatched variation/realization keys"
            )


def test_unknown_schema_type_raises() -> None:
    with pytest.raises(ValueError, match="unknown"):
        _dispatch_common_category(
            "bogus",
            {
                "$schema_type": "unsupported_type",
                "variations": {"a": 1.0},
                "realizations": {"a": "x"},
            },
        )


def test_missing_schema_type_defaults_to_weighted() -> None:
    cat = _dispatch_common_category(
        "plain",
        {
            "variations": {"a": 0.5, "b": 0.5},
            "realizations": {"a": "x", "b": "y"},
        },
    )
    assert isinstance(cat, WeightedVariations)


def test_weight_sum_failure_raises() -> None:
    with pytest.raises(ValidationError):
        WeightedVariations.model_validate(
            {
                "variations": {"a": 0.5, "b": 0.3},
                "realizations": {"a": "x", "b": "y"},
            }
        )


def test_key_mismatch_raises() -> None:
    with pytest.raises(ValidationError):
        WeightedVariations.model_validate(
            {
                "variations": {"a": 1.0},
                "realizations": {"a": "x", "extra": "y"},
            }
        )


def test_get_biomarker_case_insensitive_fallback(
    biomarkers: BiomarkerConfig,
) -> None:
    assert biomarkers.get("EGFR") is biomarkers.get("egfr")


def test_test_method_realization_lookup(common: CommonConfig) -> None:
    assert common.test_method_realization("NGS") == "NGS"
    assert common.test_method_realization("FISH") == "FISH"


def test_test_method_realization_raises_on_unknown(common: CommonConfig) -> None:
    with pytest.raises(KeyError):
        common.test_method_realization("NONEXISTENT_METHOD")


def test_biomarker_only_file_can_be_loaded_independently() -> None:
    bios = load_biomarkers(BIOMARKERS_PATH)
    assert bios.biomarkers, "biomarkers loaded without common is valid"


def test_common_file_can_be_loaded_independently() -> None:
    common = load_common(COMMON_PATH)
    assert "technical_noise" in common.categories


def test_loader_rejects_non_dict_payload(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"$schema_version": "1.0", "name_forms": 42}))
    with pytest.raises(ValueError, match="expected dict"):
        load_biomarkers(bad)
