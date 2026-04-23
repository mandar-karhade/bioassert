"""L5 negation scope: a scope marker ("No", "Absence of") negates a list
of genes wide; an optional exception clause ("but", "except") flips
polarity for a trailing gene that carries its own status surface.

Sub-phase 3.6 per docs/PHASE3_PLAN.md. Each :class:`AssertionFact` gains
a ``polarity_scope`` field: ``"direct"`` (default, L1-L4 compat),
``"negation_wide"`` (scope-marker driven, status always "negative",
status span is shared across the wide-scope facts), or ``"exception"``
(carries its own sampled status and status span).
"""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.patient_sampler import PatientProfile
from bioassert.generator.renderer import (
    AssertionFact,
    RenderedRecord,
    _L3_NAME_FORM_BLOCKLIST,
    _L5_FRAMES,
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


def test_polarity_scope_field_defaults_direct() -> None:
    """AssertionFact gains a polarity_scope field, default 'direct' so
    existing L1-L4 construction sites remain valid.
    """
    fact = AssertionFact(gene="EGFR", status="positive", spans={"gene": (0, 4)})
    assert fact.polarity_scope == "direct"


def test_l5_frames_nonempty() -> None:
    assert len(_L5_FRAMES) >= 3
    for f in _L5_FRAMES:
        assert f["kind"] in ("enumerated", "panel_wide")
        if f["kind"] == "enumerated":
            assert f["scope_marker"]
            assert f["neg_noun"]
        else:
            assert f["scope"]
            assert f["status_word"]
            assert f["exception_marker"]


def test_l5_record_structure(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(401)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        assert isinstance(rec, RenderedRecord)
        assert rec.complexity_level == "L5"
        # enumerated frames: 2..5 facts; panel_wide: exactly 1 fact (exception).
        assert 1 <= len(rec.assertions) <= 5


def test_l5_enumerated_has_at_least_two_negation_wide(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Enumerated L5 frames must carry >=2 negation_wide facts. Panel-wide
    frames emit 0 wide facts (implicit scope), so they are excluded here.
    """
    rng = random.Random(409)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    seen_enumerated = 0
    for _ in range(200):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        wide = [f for f in rec.assertions if f.polarity_scope == "negation_wide"]
        if not wide:
            continue  # panel_wide record
        seen_enumerated += 1
        assert len(wide) >= 2, (
            f"expected >=2 negation_wide facts, got "
            f"{[f.polarity_scope for f in rec.assertions]}"
        )
    assert seen_enumerated > 0, "no enumerated L5 records seen in 200 draws"


def test_l5_negation_wide_facts_are_all_negative(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(419)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        for fact in rec.assertions:
            if fact.polarity_scope == "negation_wide":
                assert fact.status == "negative", (
                    f"negation_wide fact had status={fact.status!r} in "
                    f"{rec.sentence!r}"
                )


def test_l5_negation_wide_facts_share_status_span(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(421)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        wide_status_spans = {
            f.spans["status"]
            for f in rec.assertions
            if f.polarity_scope == "negation_wide"
        }
        if not wide_status_spans:
            continue  # panel_wide record has no wide facts
        assert len(wide_status_spans) == 1, (
            f"negation_wide status spans diverged in {rec.sentence!r}: "
            f"{wide_status_spans}"
        )


def test_l5_exception_appears_sometimes(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Half the L5 frame family contains an exception clause. We should
    see at least one record in 200 samples with an exception fact.
    """
    rng = random.Random(431)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    seen_exception = False
    for _ in range(200):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        if any(f.polarity_scope == "exception" for f in rec.assertions):
            seen_exception = True
            break
    assert seen_exception, "no L5 record produced an exception fact"


def test_l5_exception_facts_carry_own_status_span(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Exception facts have their own status span — distinct from the
    shared negation_wide scope span.
    """
    rng = random.Random(433)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    checked = 0
    for _ in range(400):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        exceptions = [
            f for f in rec.assertions if f.polarity_scope == "exception"
        ]
        if not exceptions:
            continue
        wide_status_spans = {
            f.spans["status"]
            for f in rec.assertions
            if f.polarity_scope == "negation_wide"
        }
        for ex in exceptions:
            assert "status" in ex.spans
            assert ex.spans["status"] not in wide_status_spans, (
                f"exception fact reused wide-scope status span "
                f"in {rec.sentence!r}"
            )
            checked += 1
    assert checked > 0, "no exception facts observed to validate"


def test_l5_gene_surfaces_are_bare(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(439)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(300):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        for fact in rec.assertions:
            s, e = fact.spans["gene"]
            surface = rec.sentence[s:e]
            assert _L3_NAME_FORM_BLOCKLIST.search(surface) is None, (
                f"{fact.gene} surface {surface!r} embeds blocklisted keyword "
                f"in {rec.sentence!r}"
            )


def test_l5_spans_literal_substrings(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(443)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(300):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        for fact in rec.assertions:
            for name, (s, e) in fact.spans.items():
                assert 0 <= s < e <= len(rec.sentence), (
                    f"{name} span ({s},{e}) out of bounds in {rec.sentence!r}"
                )
                assert rec.sentence[s:e], f"{name} span empty slice"


def test_l5_gene_spans_disjoint_and_status_does_not_cross_genes(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Gene spans must be pairwise disjoint. Every status span (shared
    wide-scope or per-exception) must not overlap any gene span.
    """
    rng = random.Random(449)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        gene_spans = sorted(fact.spans["gene"] for fact in rec.assertions)
        for i in range(1, len(gene_spans)):
            assert gene_spans[i - 1][1] <= gene_spans[i][0], (
                f"overlapping gene spans in {rec.sentence!r}: {gene_spans}"
            )
        status_spans = {fact.spans["status"] for fact in rec.assertions}
        for ss, se in status_spans:
            for gs, ge in gene_spans:
                assert se <= gs or ge <= ss, (
                    f"status span ({ss},{se}) overlaps gene span "
                    f"({gs},{ge}) in {rec.sentence!r}"
                )


def test_l5_gene_order_matches_fact_order(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(457)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        starts = [fact.spans["gene"][0] for fact in rec.assertions]
        assert starts == sorted(starts)


def test_l5_genes_distinct_within_record(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(461)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        genes = [fact.gene for fact in rec.assertions]
        assert len(genes) == len(set(genes))


def test_l5_scope_marker_visible_in_sentence(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Every L5 sentence starts with a recognizable scope marker. For
    enumerated frames the marker ('No' / 'Absence of') must also precede
    every negation_wide gene. Panel-wide frames have no enumerated wide
    genes, so the precedence check is skipped for them.
    """
    rng = random.Random(463)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    enum_markers = {
        f["scope_marker"] for f in _L5_FRAMES if f["kind"] == "enumerated"
    }
    panel_prefixes = {
        f["scope"] for f in _L5_FRAMES if f["kind"] == "panel_wide"
    }
    for _ in range(100):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        starts_with_marker = any(
            rec.sentence.startswith(m) for m in enum_markers
        ) or any(rec.sentence.startswith(p) for p in panel_prefixes)
        assert starts_with_marker, f"no scope marker in {rec.sentence!r}"
        wide_spans = [
            f.spans["gene"] for f in rec.assertions
            if f.polarity_scope == "negation_wide"
        ]
        for (ws, _we) in wide_spans:
            head = rec.sentence[:ws]
            assert any(m in head for m in enum_markers), (
                f"scope marker must precede wide-scope gene at {ws} in "
                f"{rec.sentence!r}"
            )


def test_l5_exception_gene_follows_wide_genes(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """When an exception fact coexists with wide facts (enumerated frame
    with "and"/"except"), its gene appears after every wide gene. Panel-
    wide frames have no wide genes to compare against, so they're skipped
    when no wide facts exist in the record.
    """
    rng = random.Random(467)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    checked = 0
    for _ in range(400):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        wide_ends = [
            f.spans["gene"][1] for f in rec.assertions
            if f.polarity_scope == "negation_wide"
        ]
        ex_starts = [
            f.spans["gene"][0] for f in rec.assertions
            if f.polarity_scope == "exception"
        ]
        if not ex_starts or not wide_ends:
            continue  # panel_wide (no wide genes) or no-exception record
        for es in ex_starts:
            for we in wide_ends:
                assert we <= es, (
                    f"exception gene at {es} precedes wide gene ending "
                    f"at {we} in {rec.sentence!r}"
                )
            checked += 1
    assert checked > 0, "no enumerated exception facts with wide facts seen"


def test_l5_panel_wide_produces_single_exception_fact(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Panel-wide frames produce exactly one exception fact — no enumerated
    wide facts — because the wide scope is an implicit panel.
    """
    rng = random.Random(487)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    checked = 0
    for _ in range(500):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        has_panel_prefix = any(
            rec.sentence.startswith(f["scope"])
            for f in _L5_FRAMES
            if f["kind"] == "panel_wide"
        )
        if not has_panel_prefix:
            continue
        checked += 1
        assert len(rec.assertions) == 1, (
            f"panel_wide record must have 1 fact, got {len(rec.assertions)} "
            f"in {rec.sentence!r}"
        )
        fact = rec.assertions[0]
        assert fact.polarity_scope == "exception", (
            f"panel_wide sole fact must be polarity_scope='exception', got "
            f"{fact.polarity_scope!r} in {rec.sentence!r}"
        )
    assert checked > 0, "no panel_wide records seen in 500 draws"


def test_l5_panel_wide_surface_contains_other_than(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Panel-wide frames consistently use 'other than' as the exception
    marker in the rendered sentence.
    """
    rng = random.Random(491)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    checked = 0
    for _ in range(500):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        has_panel_prefix = any(
            rec.sentence.startswith(f["scope"])
            for f in _L5_FRAMES
            if f["kind"] == "panel_wide"
        )
        if not has_panel_prefix:
            continue
        checked += 1
        assert "other than" in rec.sentence, (
            f"panel_wide sentence missing 'other than': {rec.sentence!r}"
        )
    assert checked > 0, "no panel_wide records seen"


def test_l5_clinical_markers_only(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Clinical register: exception markers are restricted to 'and',
    'except', and 'other than'. 'but'/'however' must not appear as
    exception connectors (they're vanishingly rare in real pathology
    reports per user feedback 2026-04-23).
    """
    rng = random.Random(499)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(500):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        # These tokens should never appear as exception connectors in L5.
        # "but" as a whole word, "however" as a whole word.
        lower = rec.sentence.lower()
        assert " but " not in lower, (
            f"forbidden 'but' connector in {rec.sentence!r}"
        )
        assert " however" not in lower, (
            f"forbidden 'however' connector in {rec.sentence!r}"
        )


def test_l5_none_exception_is_dominant(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Clinical register: the majority of L5 records carry no exception
    clause at all. >= 40% no-exception is a conservative lower bound
    (frame mix currently targets ~50% None from enumerated frames alone).
    """
    rng = random.Random(503)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    no_exception = 0
    total = 500
    for _ in range(total):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        has_exception = any(
            f.polarity_scope == "exception" for f in rec.assertions
        )
        if not has_exception:
            no_exception += 1
    assert no_exception / total >= 0.40, (
        f"expected >=40% no-exception records, got {no_exception}/{total}"
    )


def test_l5_except_and_other_than_are_positive(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """'except' (enumerated polarity flip) and 'other than' (panel-wide)
    exception clauses always carry positive status on the trailing fact.
    """
    rng = random.Random(509)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    seen = 0
    for _ in range(1000):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        if " except " not in rec.sentence and "other than" not in rec.sentence:
            continue
        ex = [f for f in rec.assertions if f.polarity_scope == "exception"]
        assert ex, (
            f"'except'/'other than' sentence lacks exception fact: "
            f"{rec.sentence!r}"
        )
        seen += 1
        for fact in ex:
            assert fact.status == "positive", (
                f"'except'/'other than' exception should be positive, got "
                f"{fact.status!r} in {rec.sentence!r}"
            )
    assert seen > 0, "no 'except'/'other than' records observed"


def test_l5_and_exception_is_negative(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """'and' continuation exception carries negative polarity (same as
    wide scope). We identify 'and'-continuation frames by the presence
    of exactly one negative exception fact.
    """
    rng = random.Random(521)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    negative_exceptions = 0
    for _ in range(1000):
        rec = render_l5_record(
            list(MUTATION_BIOMARKERS), profile, biomarkers, common, rng
        )
        for fact in rec.assertions:
            if fact.polarity_scope == "exception" and fact.status == "negative":
                negative_exceptions += 1
    assert negative_exceptions > 0, (
        "no 'and'-continuation exception observed (negative-status exception)"
    )


def test_l5_rejects_pool_smaller_than_two(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(0)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    with pytest.raises(Exception):
        render_l5_record(["EGFR"], profile, biomarkers, common, rng)


def test_l5_expression_pool_supported(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(479)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(50):
        rec = render_l5_record(
            list(EXPRESSION_BIOMARKERS),
            profile,
            biomarkers,
            common,
            rng,
        )
        # Enumerated frames exhaust both expression genes (2 facts); panel_wide
        # emits a single exception fact (1 fact).
        assert 1 <= len(rec.assertions) <= 2
        assert {f.gene for f in rec.assertions} <= set(EXPRESSION_BIOMARKERS)
