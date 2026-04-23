"""Weighted sampling primitives (§2.2, §8.3, §8.4)."""
from __future__ import annotations

import random
from collections import Counter

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.patient_sampler import PatientProfile
from bioassert.generator.sampler import (
    maybe_sample_clone,
    resolve_population,
    sample_biomarker_name_form,
    sample_measurement_value,
    sample_method,
    sample_status,
    sample_variant,
    sample_variation,
)


def test_sample_variation_empirical_matches_weights(
    common: CommonConfig,
) -> None:
    verbs = common.categories["assertion_verbs"]
    rng = random.Random(0)
    n = 20000
    counts = Counter(sample_variation(verbs, rng) for _ in range(n))
    for vid, weight in verbs.variations.items():
        observed = counts[vid] / n
        assert abs(observed - weight) < 0.02, (
            f"assertion_verbs[{vid}]: observed={observed:.3f} weight={weight}"
        )


def test_sample_variation_subset_restricts_output(
    common: CommonConfig,
) -> None:
    verbs = common.categories["assertion_verbs"]
    rng = random.Random(0)
    subset = {"was", "is"}
    for _ in range(200):
        vid = sample_variation(verbs, rng, subset=subset)
        assert vid in subset


def test_sample_status_respects_distribution(
    biomarkers: BiomarkerConfig,
) -> None:
    egfr = biomarkers.get("EGFR")
    adeno_dist = egfr.status_distribution_by_population["adenocarcinoma"]
    rng = random.Random(0)
    n = 10000
    counts = Counter(sample_status(adeno_dist, rng) for _ in range(n))
    assert abs(counts["positive"] / n - adeno_dist.positive) < 0.02
    assert abs(counts["negative"] / n - adeno_dist.negative) < 0.02


def test_sample_variant_returns_known_variant(
    biomarkers: BiomarkerConfig,
) -> None:
    egfr = biomarkers.get("EGFR")
    rng = random.Random(0)
    for _ in range(100):
        vid, variant = sample_variant(egfr, rng)
        assert vid in egfr.variants
        assert variant is egfr.variants[vid]


def test_render_constraints_restrict_name_form_sampling(
    biomarkers: BiomarkerConfig,
) -> None:
    """§8.3: variant.render_constraints coerces biomarker name_form to a subset.

    ALK.unspecified_fusion has render_constraints; any sampled biomarker
    name_form must fall inside that subset.
    """
    alk = biomarkers.get("ALK")
    unspecified = alk.variants.get("unspecified_fusion")
    if unspecified is None or unspecified.render_constraints is None:
        pytest.skip("ALK.unspecified_fusion not configured with render_constraints")
    allowed = set(unspecified.render_constraints.require_biomarker_name_forms)
    rng = random.Random(0)
    for _ in range(200):
        vid = sample_biomarker_name_form(alk, rng, variant=unspecified)
        assert vid in allowed


def test_sample_biomarker_name_form_uses_full_distribution_without_constraints(
    biomarkers: BiomarkerConfig,
) -> None:
    egfr = biomarkers.get("EGFR")
    rng = random.Random(0)
    observed = set()
    for _ in range(500):
        vid = sample_biomarker_name_form(egfr, rng, variant=None)
        observed.add(vid)
    assert len(observed) >= 2


def test_maybe_sample_clone_never_fires_without_clone_block(
    biomarkers: BiomarkerConfig,
) -> None:
    egfr = biomarkers.get("EGFR")
    rng = random.Random(0)
    for _ in range(100):
        assert maybe_sample_clone(egfr, rng) is None


def test_maybe_sample_clone_bernoulli_matches_attachment_probability(
    biomarkers: BiomarkerConfig,
) -> None:
    pdl1 = biomarkers.get("PD-L1")
    if pdl1.clone_attribution is None:
        pytest.skip("PD-L1 has no clone_attribution block")
    attach_p = pdl1.clone_attribution.attachment_probability
    rng = random.Random(0)
    n = 10000
    fires = sum(1 for _ in range(n) if maybe_sample_clone(pdl1, rng) is not None)
    observed = fires / n
    assert abs(observed - attach_p) < 0.02


def test_sample_method_prefers_preferred_methods_when_present(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    egfr = biomarkers.get("EGFR")
    if egfr.preferred_methods is None:
        pytest.skip("EGFR has no preferred_methods override")
    rng = random.Random(0)
    observed: Counter[str] = Counter()
    for _ in range(5000):
        observed[sample_method(egfr, common, rng)] += 1
    observed_keys = set(observed)
    preferred_keys = set(egfr.preferred_methods.variations)
    assert observed_keys.issubset(preferred_keys)


def test_sample_method_falls_back_to_common_when_preferred_absent(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    candidate = None
    for name, entry in biomarkers.biomarkers.items():
        if entry.preferred_methods is None:
            candidate = entry
            break
    if candidate is None:
        pytest.skip("every biomarker has preferred_methods")
    rng = random.Random(0)
    vid = sample_method(candidate, common, rng)
    methods = common.categories["test_methods"]
    assert vid in methods.variations


def test_sample_measurement_value_stays_in_range() -> None:
    rng = random.Random(0)
    for _ in range(500):
        val = sample_measurement_value((0.0, 100.0), rng)
        assert 0.0 <= val <= 100.0


def test_resolve_population_returns_population_status(
    biomarkers: BiomarkerConfig,
) -> None:
    egfr = biomarkers.get("EGFR")
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    dist = resolve_population(egfr, profile)
    assert dist.positive >= 0
    assert dist.negative >= 0
