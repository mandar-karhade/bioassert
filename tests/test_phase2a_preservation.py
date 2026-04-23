"""End-to-end label preservation for the Phase 2a pipeline.

Every record emitted by the full pipeline (render → technical_noise) must
satisfy Non-Negotiable #1: the gene / variant / method / status spans
returned alongside the sentence are literal substrings at the recorded
positions, even after post-processing noise is applied.
"""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.patient_sampler import (
    PatientProfile,
    sample_patient_profile,
)
from bioassert.generator.post_process import (
    PostProcessedRecord,
    apply_technical_noise,
)
from bioassert.generator.renderer import render_l1_record

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

N_RECORDS = 1000


def _emit_records(
    n: int, common: CommonConfig, biomarkers: BiomarkerConfig, seed: int
) -> list[tuple[PatientProfile, PostProcessedRecord]]:
    rng = random.Random(seed)
    emitted: list[tuple[PatientProfile, PostProcessedRecord]] = []
    for i in range(n):
        profile = sample_patient_profile(f"p{i:05d}", rng)
        gene = rng.choice(PANEL_BIOMARKERS)
        rendered = render_l1_record(gene, profile, biomarkers, common, rng)
        post = apply_technical_noise(rendered, common, biomarkers, rng)
        emitted.append((profile, post))
    return emitted


def test_every_span_resolves_to_non_empty_substring(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    records = _emit_records(N_RECORDS, common, biomarkers, seed=2024)
    for _, rec in records:
        for name, (start, end) in rec.assertions[0].spans.items():
            assert 0 <= start < end <= len(rec.sentence), (
                f"span {name} out of bounds in {rec.sentence!r}: ({start},{end})"
            )
            assert rec.sentence[start:end], (
                f"empty slice for span {name} in {rec.sentence!r}"
            )


def test_spans_never_overlap_each_other(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    records = _emit_records(N_RECORDS, common, biomarkers, seed=99)
    for _, rec in records:
        sorted_spans = sorted(rec.assertions[0].spans.items(), key=lambda kv: kv[1])
        for i in range(1, len(sorted_spans)):
            prev_name, (_, prev_end) = sorted_spans[i - 1]
            cur_name, (cur_start, _) = sorted_spans[i]
            assert prev_end <= cur_start, (
                f"overlap: {prev_name}→{cur_name} in {rec.sentence!r}"
            )


def test_gene_surface_is_known_name_form_when_untransformed(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """When the hyphenation transform is canonical AND case is canonical,
    the gene substring must exactly match a configured name_form realization.
    """
    records = _emit_records(N_RECORDS, common, biomarkers, seed=17)
    checked = 0
    for _, rec in records:
        if rec.applied_transforms.get("hyphenation_gene_names") != "canonical":
            continue
        if rec.applied_transforms.get("case_variation") != "canonical":
            continue
        start, end = rec.assertions[0].spans["gene"]
        surface = rec.sentence[start:end]
        entry = biomarkers.get(rec.assertions[0].gene)
        known_forms = set(entry.name_forms.realizations.values())
        assert surface in known_forms, (
            f"{rec.assertions[0].gene}: gene surface {surface!r} not in configured forms"
        )
        checked += 1
    assert checked >= 500, f"only {checked} canonical-canonical records checked"


def test_status_surface_reflects_status_label(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """When case transform is canonical, the status substring (after gene
    substitution where applicable) is drawn from the status-phrase vocabulary
    for the sampled status.
    """
    vocab_for_status = {
        "positive": set(
            common.categories["positive_phrases"].realizations.values()
        ),
        "negative": set(
            common.categories["negation_phrases"].realizations.values()
        ),
        "equivocal": set(
            common.categories["equivocal_phrases"].realizations.values()
        ),
        "not_tested": set(
            common.categories["not_tested_phrases"].realizations.values()
        ),
    }

    records = _emit_records(N_RECORDS, common, biomarkers, seed=3)
    for _, rec in records:
        if rec.applied_transforms.get("case_variation") != "canonical":
            continue
        start, end = rec.assertions[0].spans["status"]
        surface = rec.sentence[start:end]
        allowed = vocab_for_status[rec.assertions[0].status]
        entry = biomarkers.get(rec.assertions[0].gene)
        known_forms = set(entry.name_forms.realizations.values())
        expanded = set()
        for form in allowed:
            for gene_form in known_forms:
                expanded.add(form.replace("{gene}", gene_form))
        assert surface in expanded, (
            f"status {rec.assertions[0].status!r} surface {surface!r} not in vocab; "
            f"sentence={rec.sentence!r}"
        )


def test_sentence_is_non_empty_and_printable(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    records = _emit_records(N_RECORDS, common, biomarkers, seed=55)
    for _, rec in records:
        assert rec.sentence.strip(), "empty rendered sentence"
        for ch in rec.sentence:
            assert ch.isprintable() or ch in "\t\n", (
                f"non-printable char {ch!r} in {rec.sentence!r}"
            )


def test_deterministic_with_same_seed(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Same seed ⇒ identical corpus."""
    a = _emit_records(200, common, biomarkers, seed=7)
    b = _emit_records(200, common, biomarkers, seed=7)
    for (p1, r1), (p2, r2) in zip(a, b):
        assert p1 == p2
        assert r1.sentence == r2.sentence
        assert r1.assertions[0].spans == r2.assertions[0].spans
        assert r1.assertions[0].status == r2.assertions[0].status


def test_applied_transforms_populated_for_every_record(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    records = _emit_records(N_RECORDS, common, biomarkers, seed=1)
    required = {
        "whitespace",
        "case_variation",
        "hyphenation_gene_names",
        "punctuation_variation",
        "ocr_corruption",
        "pdf_artifact",
        "abbreviation_inconsistency",
    }
    for _, rec in records:
        assert set(rec.applied_transforms.keys()) == required


def test_canonical_transforms_leave_sentence_identical_to_render(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(0)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    matched = 0
    for i in range(2000):
        rendered = render_l1_record("EGFR", profile, biomarkers, common, rng)
        post = apply_technical_noise(rendered, common, biomarkers, rng)
        if all(
            post.applied_transforms[k]
            in ("canonical", "single_space")
            for k in post.applied_transforms
        ):
            assert post.sentence == rendered.sentence
            assert post.assertions[0].spans == rendered.assertions[0].spans
            matched += 1
    assert matched >= 200, f"only {matched} pure-canonical records seen"
