"""Patient-profile axis sampling and population-key cascade (§2.3)."""
from __future__ import annotations

import logging
import random

import pytest

from bioassert.config import BiomarkerConfig
from bioassert.generator.patient_sampler import (
    DEFAULT_PROFILE_DISTRIBUTION,
    PatientProfile,
    PopulationCascadeMiss,
    resolve_status_distribution,
    sample_patient_profile,
)


def test_cascade_full_profile_emits_five_levels() -> None:
    profile = PatientProfile(
        patient_ref="p1",
        histology="adenocarcinoma",
        ethnicity="east_asian",
        smoking="nonsmoker",
        age_group="younger",
        sex="female",
    )
    cascade = profile.population_key_cascade()
    assert cascade == [
        "adenocarcinoma_east_asian_nonsmoker_younger_female",
        "adenocarcinoma_east_asian_nonsmoker_younger",
        "adenocarcinoma_east_asian_nonsmoker",
        "adenocarcinoma_east_asian",
        "adenocarcinoma",
    ]


def test_cascade_omits_optional_axes_that_are_none() -> None:
    profile = PatientProfile(
        patient_ref="p1", histology="squamous", smoking="smoker"
    )
    cascade = profile.population_key_cascade()
    assert "squamous_smoker" in cascade
    assert "squamous" in cascade
    for key in cascade:
        assert "east_asian" not in key
        assert "male" not in key


def test_cascade_deduplicates_when_all_optional_axes_are_none() -> None:
    profile = PatientProfile(patient_ref="p1", histology="adenocarcinoma")
    cascade = profile.population_key_cascade()
    assert cascade == ["adenocarcinoma"]


def test_resolve_status_distribution_matches_specific_key_when_present(
    biomarkers: BiomarkerConfig,
) -> None:
    """ALK has a specific 'adenocarcinoma_nonsmoker_younger' entry — the
    cascade should land on it before falling back to 'adenocarcinoma'.
    """
    profile = PatientProfile(
        patient_ref="p1",
        histology="adenocarcinoma",
        smoking="nonsmoker",
        age_group="younger",
    )
    alk = biomarkers.get("ALK")
    dist, matched = resolve_status_distribution(alk, profile)
    assert matched == "adenocarcinoma_nonsmoker_younger"
    assert dist.positive > 0


def test_resolve_status_distribution_falls_back_to_histology(
    biomarkers: BiomarkerConfig,
) -> None:
    profile = PatientProfile(
        patient_ref="p1",
        histology="squamous",
        ethnicity="east_asian",
        smoking="nonsmoker",
    )
    egfr = biomarkers.get("EGFR")
    dist, matched = resolve_status_distribution(egfr, profile)
    assert matched == "squamous"
    assert dist.positive >= 0


def test_resolve_status_distribution_emits_warning_on_fallback(
    caplog: pytest.LogCaptureFixture, biomarkers: BiomarkerConfig
) -> None:
    from bioassert.generator import patient_sampler

    patient_sampler._WARNED.clear()
    profile = PatientProfile(
        patient_ref="p1",
        histology="squamous",
        ethnicity="east_asian",
        smoking="nonsmoker",
    )
    with caplog.at_level(logging.WARNING):
        resolve_status_distribution(biomarkers.get("EGFR"), profile)
    assert any("population cascade fallback" in r.message for r in caplog.records)


def test_sample_patient_profile_honours_override_distribution() -> None:
    override = {
        "histology": {"adenocarcinoma": 1.0},
        "ethnicity": {"east_asian": 1.0},
        "smoking": {"nonsmoker": 1.0},
        "age_group": {"younger": 1.0},
        "sex": {"female": 1.0},
    }
    rng = random.Random(0)
    profile = sample_patient_profile("p", rng=rng, distribution=override)
    assert profile.histology == "adenocarcinoma"
    assert profile.ethnicity == "east_asian"
    assert profile.smoking == "nonsmoker"
    assert profile.age_group == "younger"
    assert profile.sex == "female"


def test_sample_patient_profile_default_distribution_sums_to_one() -> None:
    for axis, dist in DEFAULT_PROFILE_DISTRIBUTION.items():
        assert abs(sum(dist.values()) - 1.0) < 1e-6, axis


def test_cascade_miss_raises_population_cascade_miss() -> None:
    from bioassert.config.schema import Biomarker

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
                        "variations": {"only_form": 1.0},
                        "realizations": {"only_form": ""},
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
    object.__setattr__(
        entry, "status_distribution_by_population", {}
    )
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    with pytest.raises(PopulationCascadeMiss):
        resolve_status_distribution(entry, profile)
