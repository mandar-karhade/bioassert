"""Sentence renderer: placeholder substitution, span tracking, frame picking."""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.patient_sampler import PatientProfile
from bioassert.generator.renderer import (
    RenderError,
    RenderedRecord,
    expand_placeholders,
    render_l1_record,
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

PANEL_BIOMARKERS: tuple[str, ...] = MUTATION_BIOMARKERS + EXPRESSION_BIOMARKERS


def test_expand_placeholders_fills_known_names() -> None:
    out = expand_placeholders("no {gene} alteration", {"gene": "EGFR"})
    assert out == "no EGFR alteration"


def test_expand_placeholders_raises_on_unknown_name() -> None:
    with pytest.raises(RenderError, match="unknown placeholder"):
        expand_placeholders("at {unsupported}", {"unsupported": "1"})


def test_expand_placeholders_raises_on_missing_context() -> None:
    with pytest.raises(KeyError, match="method"):
        expand_placeholders("by {method}", {})


def test_render_l1_record_returns_valid_sentence(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(42)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    rec = render_l1_record("EGFR", profile, biomarkers, common, rng)
    assert isinstance(rec, RenderedRecord)
    assert rec.sentence
    assert "gene" in rec.assertions[0].spans
    assert "status" in rec.assertions[0].spans


def test_gene_span_points_to_gene_surface(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(1)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        gene = rng.choice(MUTATION_BIOMARKERS)
        rec = render_l1_record(gene, profile, biomarkers, common, rng)
        start, end = rec.assertions[0].spans["gene"]
        assert 0 <= start < end <= len(rec.sentence)
        substring = rec.sentence[start:end]
        gene_forms = biomarkers.get(gene).name_forms.realizations.values()
        assert substring in set(gene_forms), (
            f"gene slice {substring!r} not in known name_forms for {gene}"
        )


def test_variant_span_when_present_matches_rendered_variant_surface(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(7)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    checked = 0
    attempts = 0
    while checked < 40 and attempts < 1000:
        attempts += 1
        rec = render_l1_record("EGFR", profile, biomarkers, common, rng)
        if "variant" not in rec.assertions[0].spans:
            continue
        start, end = rec.assertions[0].spans["variant"]
        assert rec.sentence[start:end], "variant substring must be non-empty"
        checked += 1
    assert checked >= 10, f"found only {checked} variant-bearing records"


def test_method_span_when_present_matches_common_realizations(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(3)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    realizations = {
        r
        for r in common.categories["test_methods"].realizations.values()
        if r
    }
    checked = 0
    attempts = 0
    while checked < 30 and attempts < 1000:
        attempts += 1
        rec = render_l1_record(
            "EGFR", profile, biomarkers, common, rng, method_attach_prob=0.9
        )
        if "method" not in rec.assertions[0].spans:
            continue
        start, end = rec.assertions[0].spans["method"]
        assert rec.sentence[start:end] in realizations
        checked += 1
    assert checked >= 10, f"only {checked} method-bearing records"


def test_expression_positive_renders_variant_with_value(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """PD-L1 positive: descriptor slot filled from variants, {value}% possible."""
    rng = random.Random(0)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    positive_seen = 0
    value_substituted = 0
    for _ in range(600):
        rec = render_l1_record(
            "PD-L1", profile, biomarkers, common, rng, method_attach_prob=0.0
        )
        if rec.assertions[0].status != "positive":
            continue
        positive_seen += 1
        assert rec.assertions[0].variant_id is not None, "positive expression needs variant_id"
        assert rec.assertions[0].negative_form_id is None
        if rec.assertions[0].measurement_value is not None:
            value_substituted += 1
            lo, hi = biomarkers.get("PD-L1").variants[rec.assertions[0].variant_id].measurement_range
            assert lo <= rec.assertions[0].measurement_value <= hi
    assert positive_seen >= 50, f"only {positive_seen} PD-L1 positives"
    assert value_substituted >= 10, (
        f"expected some {{value}} substitutions, saw {value_substituted}"
    )


def test_expression_negative_draws_from_negative_forms(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """TMB negative: descriptor filled from negative_forms when frame renders it."""
    rng = random.Random(1)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    tmb = biomarkers.get("TMB")
    neg = tmb.negative_forms
    assert neg is not None
    literal_forms = {r for r in neg.realizations.values() if "{value}" not in r}
    value_templates = [r for r in neg.realizations.values() if "{value}" in r]
    neg_seen = 0
    descriptor_rendered = 0
    for _ in range(600):
        rec = render_l1_record(
            "TMB", profile, biomarkers, common, rng, method_attach_prob=0.0
        )
        if rec.assertions[0].status != "negative":
            continue
        neg_seen += 1
        if rec.assertions[0].negative_form_id is None:
            continue
        assert rec.assertions[0].variant_id is None
        if "variant" not in rec.assertions[0].spans:
            continue
        descriptor_rendered += 1
        start, end = rec.assertions[0].spans["variant"]
        surface = rec.sentence[start:end]
        if surface in literal_forms:
            continue
        assert any(
            _matches_value_template(surface, t) for t in value_templates
        ), f"neg descriptor {surface!r} not a TMB negative_form realization"
    assert neg_seen >= 100
    assert descriptor_rendered >= 15, (
        f"expected some rendered neg descriptors, saw {descriptor_rendered}"
    )


def _matches_value_template(surface: str, template: str) -> bool:
    """True if ``surface`` matches ``template`` with ``{value}`` replaced by digits."""
    prefix, _, suffix = template.partition("{value}")
    if not surface.startswith(prefix) or not surface.endswith(suffix):
        return False
    middle = surface[len(prefix) : len(surface) - len(suffix) or None]
    return bool(middle) and all(c in "0123456789." for c in middle)


def test_clone_attribution_attaches_and_tracks_span(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """PD-L1 records occasionally get a clone suffix with its own span."""
    rng = random.Random(2)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    clone_attributions = biomarkers.get("PD-L1").clone_attribution
    assert clone_attributions is not None
    clone_realizations = set(clone_attributions.realizations.values())
    clone_seen = 0
    for _ in range(400):
        rec = render_l1_record(
            "PD-L1", profile, biomarkers, common, rng, method_attach_prob=0.0
        )
        if rec.assertions[0].clone_id is None:
            assert "clone" not in rec.assertions[0].spans
            continue
        clone_seen += 1
        assert "clone" in rec.assertions[0].spans
        start, end = rec.assertions[0].spans["clone"]
        assert rec.sentence[start:end] in clone_realizations
        other_spans = [s for name, s in rec.assertions[0].spans.items() if name != "clone"]
        for s, e in other_spans:
            assert e <= start, "clone should come after every other span"
    assert clone_seen >= 50, f"only {clone_seen} clone attachments in 400 tries"


def test_rendered_spans_are_non_overlapping(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(11)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l1_record("EGFR", profile, biomarkers, common, rng)
        spans = list(rec.assertions[0].spans.values())
        spans.sort()
        for i in range(1, len(spans)):
            assert spans[i - 1][1] <= spans[i][0], (
                f"overlapping spans in {rec.sentence!r}: {rec.assertions[0].spans}"
            )


def test_status_span_matches_status_phrase_realization(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(5)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    categories = {
        "positive": common.categories["positive_phrases"].realizations.values(),
        "negative": common.categories["negation_phrases"].realizations.values(),
        "equivocal": common.categories["equivocal_phrases"].realizations.values(),
        "not_tested": common.categories["not_tested_phrases"].realizations.values(),
    }
    for _ in range(200):
        rec = render_l1_record("EGFR", profile, biomarkers, common, rng)
        start, end = rec.assertions[0].spans["status"]
        surface = rec.sentence[start:end]
        allowed = {r for r in categories[rec.assertions[0].status]}
        expanded = set()
        for form in allowed:
            expanded.add(form.replace("{gene}", "EGFR"))
        assert surface in expanded or any(
            surface.startswith(a) for a in expanded
        ), f"status surface {surface!r} not from {rec.assertions[0].status} vocab"


def test_l2_positive_draws_from_positive_shorthand(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L2 positive records use the compact positive_shorthand realizations."""
    shorthand = set(common.categories["positive_shorthand"].realizations.values())
    rng = random.Random(0)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    seen = 0
    for _ in range(400):
        rec = render_l1_record(
            "EGFR", profile, biomarkers, common, rng, complexity_level="L2"
        )
        if rec.assertions[0].status != "positive":
            continue
        assert rec.complexity_level == "L2"
        start, end = rec.assertions[0].spans["status"]
        assert rec.sentence[start:end] in shorthand
        seen += 1
    assert seen >= 40, f"only {seen} L2 positives"


def test_l2_negative_draws_from_negation_shorthand(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L2 negative records use the compact negation_shorthand realizations."""
    shorthand = set(common.categories["negation_shorthand"].realizations.values())
    rng = random.Random(1)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    seen = 0
    for _ in range(400):
        rec = render_l1_record(
            "KRAS", profile, biomarkers, common, rng, complexity_level="L2"
        )
        if rec.assertions[0].status != "negative":
            continue
        start, end = rec.assertions[0].spans["status"]
        assert rec.sentence[start:end] in shorthand
        seen += 1
    assert seen >= 20, f"only {seen} L2 negatives"


def test_l2_suppresses_variant_and_method(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L2 frames never render variant descriptor or method slot."""
    rng = random.Random(2)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(300):
        rec = render_l1_record(
            "EGFR",
            profile,
            biomarkers,
            common,
            rng,
            complexity_level="L2",
            method_attach_prob=0.9,
        )
        assert "variant" not in rec.assertions[0].spans
        assert "method" not in rec.assertions[0].spans
        assert rec.assertions[0].test_method is None


def test_l2_frames_are_shorter_than_l1_on_average(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Sanity check: L2 sentences are meaningfully shorter than L1."""
    rng = random.Random(3)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    l1_len = 0
    l2_len = 0
    n = 300
    for _ in range(n):
        l1 = render_l1_record("EGFR", profile, biomarkers, common, rng)
        l2 = render_l1_record(
            "EGFR", profile, biomarkers, common, rng, complexity_level="L2"
        )
        l1_len += len(l1.sentence)
        l2_len += len(l2.sentence)
    assert l2_len < l1_len, "L2 corpus should be shorter than L1 on average"
    avg_l2 = l2_len / n
    assert avg_l2 < 25, f"L2 average length {avg_l2} too long"


def test_l2_spans_resolve_to_non_empty_substrings(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L2 records preserve the span-is-substring invariant."""
    rng = random.Random(4)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(500):
        gene = rng.choice(PANEL_BIOMARKERS)
        rec = render_l1_record(
            gene, profile, biomarkers, common, rng, complexity_level="L2"
        )
        assert rec.sentence.strip()
        for name, (start, end) in rec.assertions[0].spans.items():
            assert 0 <= start < end <= len(rec.sentence), (
                f"span {name} out of bounds: {(start, end)} in {rec.sentence!r}"
            )
            assert rec.sentence[start:end], (
                f"empty slice for {name} in {rec.sentence!r}"
            )


def test_l2_equivocal_falls_back_to_formal_vocabulary(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """No shorthand category exists for equivocal — L2 uses formal equivocal_phrases."""
    formal = set(common.categories["equivocal_phrases"].realizations.values())
    rng = random.Random(5)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    seen = 0
    for _ in range(1500):
        rec = render_l1_record(
            "EGFR", profile, biomarkers, common, rng, complexity_level="L2"
        )
        if rec.assertions[0].status != "equivocal":
            continue
        start, end = rec.assertions[0].spans["status"]
        assert rec.sentence[start:end] in formal
        seen += 1
        if seen >= 5:
            break
    assert seen >= 1, "no equivocal L2 records observed"


def test_render_rejects_unknown_complexity_level(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(0)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    with pytest.raises(RenderError, match="complexity_level"):
        render_l1_record(
            "EGFR", profile, biomarkers, common, rng, complexity_level="L3"
        )
