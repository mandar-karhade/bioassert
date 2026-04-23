"""L7 cross-sentence coreference (Sub-phase 3.9).

L7 records are multi-sentence discourse units. Each ``AssertionFact`` carries
a ``sentence_index`` locating it in the record's ``sentences`` tuple; fact
spans index into the full ``sentence`` string (the joined document), and
must fall entirely within the character range of ``sentences[sentence_index]``
inside that joined string.

Three structural shapes:

- ``"setup_claim"`` — 2 sentences: generic setup then gene/status claim.
  The fact's ``sentence_index`` points at sentence 1.
- ``"claim_anaphora"`` — 2 sentences: gene/status claim then nominal
  anaphor sentence ("the mutation...", "this variant..."). Fact index 0.
  The anaphor sentence carries no labeled spans.
- ``"setup_claim_qualifier"`` — 3 sentences: setup, claim, qualifier.
  Fact index 1. Neither setup nor qualifier carries labeled spans.

L7 is single-fact (one assertion per record); multi-gene L7 is deferred.
``compounding_tier`` is always ``"low"`` because the tier knob doesn't apply
to single-fact records.
"""
from __future__ import annotations

import random

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.patient_sampler import PatientProfile
from bioassert.generator.renderer import (
    AssertionFact,
    L7_ANAPHOR_NOUNS,
    L7_SHAPES,
    RenderError,
    RenderedRecord,
    render_l7_record,
)

GENES: tuple[str, ...] = (
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


# --- schema ------------------------------------------------------------


def test_assertion_fact_sentence_index_default_zero() -> None:
    """L1-L6 construction sites must remain valid without passing
    sentence_index. It defaults to 0.
    """
    fact = AssertionFact(gene="EGFR", status="positive", spans={})
    assert fact.sentence_index == 0


def test_rendered_record_sentences_defaults_to_single_sentence_tuple() -> None:
    """L1-L6 constructors pass only ``sentence=``. The ``sentences`` tuple
    defaults to ``(sentence,)`` so existing records remain valid.
    """
    rec = RenderedRecord(
        sentence="EGFR was positive.",
        assertions=(),
        frame_template="{gene} was {status}.",
    )
    assert rec.sentences == ("EGFR was positive.",)


def test_shape_constants_exact() -> None:
    assert L7_SHAPES == ("setup_claim", "claim_anaphora", "setup_claim_qualifier")


def test_anaphor_vocab_nonempty() -> None:
    assert len(L7_ANAPHOR_NOUNS) >= 3
    for noun in L7_ANAPHOR_NOUNS:
        assert noun == noun.lower()


# --- render_l7_record: setup_claim shape -------------------------------


def test_l7_setup_claim_two_sentences_claim_has_span(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(901)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l7_record(
            list(GENES), profile, biomarkers, common, rng, shape="setup_claim"
        )
        assert rec.complexity_level == "L7"
        assert rec.compounding_tier == "low"
        assert len(rec.sentences) == 2
        assert len(rec.assertions) == 1
        fact = rec.assertions[0]
        assert fact.sentence_index == 1, "claim is the second sentence"


# --- render_l7_record: claim_anaphora shape ----------------------------


def test_l7_claim_anaphora_two_sentences_claim_first(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(903)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l7_record(
            list(GENES), profile, biomarkers, common, rng, shape="claim_anaphora"
        )
        assert rec.complexity_level == "L7"
        assert len(rec.sentences) == 2
        assert len(rec.assertions) == 1
        fact = rec.assertions[0]
        assert fact.sentence_index == 0, "claim is the first sentence"


def test_l7_claim_anaphora_second_sentence_carries_anaphor(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """The anaphor sentence must contain a known anaphor noun from
    ``L7_ANAPHOR_NOUNS``."""
    rng = random.Random(905)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l7_record(
            list(GENES), profile, biomarkers, common, rng, shape="claim_anaphora"
        )
        anaphor_sentence = rec.sentences[1]
        assert any(n in anaphor_sentence for n in L7_ANAPHOR_NOUNS), (
            f"anaphor sentence {anaphor_sentence!r} lacks a known anaphor noun"
        )


# --- render_l7_record: setup_claim_qualifier shape ---------------------


def test_l7_setup_claim_qualifier_three_sentences_claim_middle(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(907)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l7_record(
            list(GENES), profile, biomarkers, common, rng,
            shape="setup_claim_qualifier",
        )
        assert len(rec.sentences) == 3
        assert len(rec.assertions) == 1
        fact = rec.assertions[0]
        assert fact.sentence_index == 1, "claim is the middle sentence"


def test_l7_qualifier_sentence_carries_anaphor(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(909)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(100):
        rec = render_l7_record(
            list(GENES), profile, biomarkers, common, rng,
            shape="setup_claim_qualifier",
        )
        qualifier = rec.sentences[2]
        assert any(n in qualifier for n in L7_ANAPHOR_NOUNS)


# --- shape dispatch ----------------------------------------------------


def test_l7_random_shape_covers_all_three(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """With ``shape=None`` the renderer picks uniformly from the three shapes."""
    rng = random.Random(911)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    shapes_seen: set[str] = set()
    for _ in range(300):
        rec = render_l7_record(list(GENES), profile, biomarkers, common, rng)
        fact = rec.assertions[0]
        n = len(rec.sentences)
        if n == 3:
            shapes_seen.add("setup_claim_qualifier")
        elif n == 2 and fact.sentence_index == 1:
            shapes_seen.add("setup_claim")
        elif n == 2 and fact.sentence_index == 0:
            shapes_seen.add("claim_anaphora")
    assert shapes_seen == {"setup_claim", "claim_anaphora", "setup_claim_qualifier"}


def test_l7_rejects_unknown_shape(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(913)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    with pytest.raises(RenderError):
        render_l7_record(
            list(GENES), profile, biomarkers, common, rng, shape="bogus"
        )


# --- structural invariants --------------------------------------------


def test_l7_sentence_is_join_of_sentences(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """The record-level ``sentence`` text is the space-joined concatenation
    of ``sentences``. This gives downstream post-processing a flat surface
    while keeping sentence boundaries explicit.
    """
    rng = random.Random(915)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        rec = render_l7_record(list(GENES), profile, biomarkers, common, rng)
        assert rec.sentence == " ".join(rec.sentences)


def test_l7_sentence_index_within_bounds(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(917)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        rec = render_l7_record(list(GENES), profile, biomarkers, common, rng)
        for fact in rec.assertions:
            assert 0 <= fact.sentence_index < len(rec.sentences)


def test_l7_spans_lie_within_referenced_sentence(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Fact spans are global into ``sentence`` (the joined text) but must
    fall entirely inside the character range of ``sentences[sentence_index]``.
    """
    rng = random.Random(919)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(300):
        rec = render_l7_record(list(GENES), profile, biomarkers, common, rng)
        # Compute each sentence's byte range inside the joined surface. The
        # joiner is a single space, so sentence k begins at
        # sum(len(s_i)+1 for i<k) and ends at that start + len(s_k).
        offsets: list[tuple[int, int]] = []
        cursor = 0
        for i, s in enumerate(rec.sentences):
            offsets.append((cursor, cursor + len(s)))
            cursor += len(s) + (1 if i < len(rec.sentences) - 1 else 0)
        for fact in rec.assertions:
            s_start, s_end = offsets[fact.sentence_index]
            for span_name, (a, b) in fact.spans.items():
                assert s_start <= a < b <= s_end, (
                    f"{span_name} span ({a},{b}) escapes "
                    f"sentence {fact.sentence_index} range {offsets[fact.sentence_index]}"
                )


def test_l7_every_fact_has_gene_and_status_span(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(921)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        rec = render_l7_record(list(GENES), profile, biomarkers, common, rng)
        for fact in rec.assertions:
            assert "gene" in fact.spans
            assert "status" in fact.spans


def test_l7_gene_span_surface_is_known_name_form(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    rng = random.Random(923)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        rec = render_l7_record(list(GENES), profile, biomarkers, common, rng)
        for fact in rec.assertions:
            s, e = fact.spans["gene"]
            surface = rec.sentence[s:e]
            known = set(
                biomarkers.get(fact.gene).name_forms.realizations.values()
            )
            assert surface in known


# --- vocab coverage (smoke) -------------------------------------------


def test_l7_anaphor_vocab_is_exercised(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Over many draws of anaphor-bearing shapes, every anaphor noun should
    surface at least once.
    """
    rng = random.Random(925)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    seen: set[str] = set()
    for _ in range(1200):
        shape = rng.choice(("claim_anaphora", "setup_claim_qualifier"))
        rec = render_l7_record(
            list(GENES), profile, biomarkers, common, rng, shape=shape
        )
        anaphor_sentence = rec.sentences[-1]
        for n in L7_ANAPHOR_NOUNS:
            if n in anaphor_sentence:
                seen.add(n)
    assert seen == set(L7_ANAPHOR_NOUNS)


# --- discourse sentences carry no labeled spans ------------------------


def test_l7_non_claim_sentences_contain_no_fact_spans(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Setup and qualifier/anaphor sentences must not accidentally match
    labeled fact surfaces — they are pure discourse scaffolding.
    """
    rng = random.Random(927)
    profile = PatientProfile(patient_ref="p", histology="adenocarcinoma")
    for _ in range(200):
        rec = render_l7_record(list(GENES), profile, biomarkers, common, rng)
        offsets: list[tuple[int, int]] = []
        cursor = 0
        for i, s in enumerate(rec.sentences):
            offsets.append((cursor, cursor + len(s)))
            cursor += len(s) + (1 if i < len(rec.sentences) - 1 else 0)
        for i, _ in enumerate(rec.sentences):
            if i == rec.assertions[0].sentence_index:
                continue
            s_start, s_end = offsets[i]
            for fact in rec.assertions:
                for a, b in fact.spans.values():
                    assert not (s_start <= a < s_end), (
                        f"discourse sentence {i} overlaps fact span ({a},{b})"
                    )
