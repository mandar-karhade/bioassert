"""Cross-config invariants (config_architecture.md §7.5, §7.6, §2.3)."""
from __future__ import annotations

from copy import deepcopy

import pytest
from pydantic import ValidationError

from bioassert.config import BiomarkerConfig, CommonConfig, Biomarker
from bioassert.config.validator import (
    SUPPORTED_PLACEHOLDERS,
    ConfigValidationError,
    _parse_population_key,
    describe_biomarker_shape,
    validate_configs,
)


def test_supported_placeholders_are_exactly_five() -> None:
    assert SUPPORTED_PLACEHOLDERS == {
        "gene",
        "method",
        "result",
        "value",
        "variant",
    }


def test_parse_population_key_accepts_canonical_order() -> None:
    assert _parse_population_key("adenocarcinoma") == [
        ("histology", "adenocarcinoma")
    ]
    assert _parse_population_key("adenocarcinoma_east_asian") == [
        ("histology", "adenocarcinoma"),
        ("ethnicity", "east_asian"),
    ]
    assert _parse_population_key(
        "adenocarcinoma_east_asian_nonsmoker_younger_female"
    ) == [
        ("histology", "adenocarcinoma"),
        ("ethnicity", "east_asian"),
        ("smoking", "nonsmoker"),
        ("age_group", "younger"),
        ("sex", "female"),
    ]


def test_parse_population_key_rejects_out_of_order() -> None:
    with pytest.raises(ConfigValidationError):
        _parse_population_key("adenocarcinoma_younger_nonsmoker")


def test_parse_population_key_rejects_non_histology_first() -> None:
    with pytest.raises(ConfigValidationError):
        _parse_population_key("east_asian_adenocarcinoma")


def test_validate_configs_passes_on_shipped_configs(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    validate_configs(common, biomarkers)


def test_validate_configs_detects_unknown_preferred_method(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    perturbed = _clone_biomarkers(biomarkers)
    egfr = perturbed.biomarkers["EGFR"]
    assert egfr.preferred_methods is not None
    bad_weights = dict(egfr.preferred_methods.variations)
    bad_weights["BOGUS_METHOD"] = 0.0
    object.__setattr__(
        egfr.preferred_methods, "variations", bad_weights
    )
    with pytest.raises(ConfigValidationError, match="BOGUS_METHOD"):
        validate_configs(common, perturbed)


def test_validate_configs_detects_unknown_placeholder(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    perturbed = _clone_common(common)
    verbs = perturbed.categories["assertion_verbs"]
    first_vid = next(iter(verbs.realizations))
    verbs.realizations[first_vid] = "{unsupported} verb"
    with pytest.raises(ConfigValidationError, match="unsupported"):
        validate_configs(perturbed, biomarkers)


def test_render_constraints_references_known_name_form_keys() -> None:
    """Biomarker model_validator rejects unknown name_form keys inside
    variant.render_constraints.require_biomarker_name_forms."""
    with pytest.raises(ValidationError):
        Biomarker.model_validate(
            {
                "canonical_name": "TEST",
                "alteration_type": "fusion",
                "name_forms": {
                    "variations": {"TEST": 1.0},
                    "realizations": {"TEST": "TEST"},
                },
                "variants": {
                    "unspecified_fusion": {
                        "actionable": True,
                        "prevalence_within_biomarker": 1.0,
                        "name_forms": {
                            "variations": {"blank": 1.0},
                            "realizations": {"blank": ""},
                        },
                        "render_constraints": {
                            "require_biomarker_name_forms": [
                                "not_a_real_name_form"
                            ]
                        },
                    }
                },
                "status_distribution_by_population": {
                    "adenocarcinoma": {
                        "positive": 0.1,
                        "negative": 0.8,
                        "equivocal": 0.05,
                        "not_tested": 0.05,
                    },
                    "squamous": {
                        "positive": 0.05,
                        "negative": 0.85,
                        "equivocal": 0.05,
                        "not_tested": 0.05,
                    },
                },
            }
        )


def test_describe_biomarker_shape_summarizes_features(
    biomarkers: BiomarkerConfig,
) -> None:
    egfr_shape = describe_biomarker_shape(biomarkers.get("EGFR"))
    assert egfr_shape["alteration_type"] == "mutation"
    assert egfr_shape["variant_count"] >= 1
    pdl1_shape = describe_biomarker_shape(biomarkers.get("PD-L1"))
    assert pdl1_shape["alteration_type"] == "expression"
    assert pdl1_shape["has_negative_forms"] is True


def test_value_placeholder_requires_measurement_range() -> None:
    """Loader rejects {value} realizations without a measurement_range."""
    with pytest.raises(ConfigValidationError, match="measurement_range"):
        common_bad = _make_minimal_common()
        biomarkers_bad = _make_minimal_biomarkers_with_value_no_range()
        validate_configs(common_bad, biomarkers_bad)


def _clone_biomarkers(b: BiomarkerConfig) -> BiomarkerConfig:
    return BiomarkerConfig.model_validate(b.model_dump())


def _clone_common(c: CommonConfig) -> CommonConfig:
    copied = deepcopy(c)
    return copied


def _make_minimal_common() -> CommonConfig:
    from bioassert.config.loader import _dispatch_common_category

    cat = _dispatch_common_category(
        "test_methods",
        {
            "$schema_type": "weighted_variations",
            "variations": {"unspecified": 1.0},
            "realizations": {"unspecified": ""},
        },
    )
    return CommonConfig(schema_version="1.0", categories={"test_methods": cat})


def _make_minimal_biomarkers_with_value_no_range() -> BiomarkerConfig:
    entry = Biomarker.model_validate(
        {
            "canonical_name": "TEST",
            "alteration_type": "mutation",
            "name_forms": {
                "variations": {"TEST": 1.0},
                "realizations": {"TEST": "TEST"},
            },
            "variants": {
                "only": {
                    "actionable": True,
                    "prevalence_within_biomarker": 1.0,
                    "name_forms": {
                        "variations": {"with_value": 1.0},
                        "realizations": {"with_value": "at {value}%"},
                    },
                }
            },
            "status_distribution_by_population": {
                "adenocarcinoma": {
                    "positive": 0.1,
                    "negative": 0.8,
                    "equivocal": 0.05,
                    "not_tested": 0.05,
                },
                "squamous": {
                    "positive": 0.05,
                    "negative": 0.85,
                    "equivocal": 0.05,
                    "not_tested": 0.05,
                },
            },
        }
    )
    return BiomarkerConfig(schema_version="1.0", biomarkers={"TEST": entry})
