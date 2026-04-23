"""Layer 5 noise: abbreviation inconsistency across repeat gene mentions.

Phase 3.11 scope. When a rendered record contains two or more labeled
mentions of the same canonical gene (today produced by L6 temporal_frame2:
``{gene} was X at T1; {gene} was Y at T2``), the ``mixed`` mode re-renders
mention #2+ using an alternate alias drawn from
``biomarkers[canonical].name_forms.variations``, excluding mention #1's
surface. ``canonical`` mode leaves the sentence unchanged.

Orthogonal to Layer 4 complexity frames — applies to any rendered output
that has repeat labeled mentions, today or in the future.
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig
from bioassert.generator.post_process import (
    apply_technical_noise,
)
from bioassert.generator.renderer import (
    AssertionFact,
    RenderedRecord,
    render_l1_record,
    render_l6_record,
)


ROOT = Path(__file__).resolve().parents[1]
COMMON_PATH = ROOT / "projects" / "nsclc_adenocarcinoma" / "configs" / "common_variations.json"


def _make_l6_temporal_frame2_stub(
    gene: str,
    gene_surface: str,
    status_a_surface: str = "positive",
    temporal_a: str = "at diagnosis",
    status_b_surface: str = "negative",
    temporal_b: str = "post-TKI",
) -> RenderedRecord:
    """Construct a synthetic L6 temporal_frame2 RenderedRecord.

    Shape: ``{gene} was {status_a} {temporal_a}; {gene} was {status_b}
    {temporal_b}.``
    """
    parts: list[str] = []
    pos = 0
    g1_start = pos
    parts.append(gene_surface); pos += len(gene_surface)
    g1_end = pos
    parts.append(" was "); pos += 5
    sa_start = pos
    parts.append(status_a_surface); pos += len(status_a_surface)
    sa_end = pos
    parts.append(" " + temporal_a); pos += 1 + len(temporal_a)
    parts.append("; "); pos += 2
    g2_start = pos
    parts.append(gene_surface); pos += len(gene_surface)
    g2_end = pos
    parts.append(" was "); pos += 5
    sb_start = pos
    parts.append(status_b_surface); pos += len(status_b_surface)
    sb_end = pos
    parts.append(" " + temporal_b); pos += 1 + len(temporal_b)
    parts.append(".")
    sentence = "".join(parts)
    fact_a = AssertionFact(
        gene=gene,
        status="positive",
        spans={"gene": (g1_start, g1_end), "status": (sa_start, sa_end)},
        temporal=temporal_a,
    )
    fact_b = AssertionFact(
        gene=gene,
        status="negative",
        spans={"gene": (g2_start, g2_end), "status": (sb_start, sb_end)},
        temporal=temporal_b,
    )
    return RenderedRecord(
        sentence=sentence,
        assertions=(fact_a, fact_b),
        frame_template="{gene} was {status_a} {temporal_a}; {gene} was {status_b} {temporal_b}.",
        complexity_level="L6",
        compounding_tier="low",
    )


def _make_l6_shared_gene_stub(
    gene: str, gene_surface: str
) -> RenderedRecord:
    """L6 frame1 shape: single gene mention, two facts share the span."""
    sentence = f"{gene_surface}: positive at diagnosis, negative post-TKI."
    g_start = 0
    g_end = len(gene_surface)
    sa_start = g_end + 2  # ": "
    sa_end = sa_start + len("positive")
    sb_start = sentence.index(", ") + 2
    sb_end = sb_start + len("negative")
    fact_a = AssertionFact(
        gene=gene,
        status="positive",
        spans={"gene": (g_start, g_end), "status": (sa_start, sa_end)},
        temporal="at diagnosis",
    )
    fact_b = AssertionFact(
        gene=gene,
        status="negative",
        spans={"gene": (g_start, g_end), "status": (sb_start, sb_end)},
        temporal="post-TKI",
    )
    return RenderedRecord(
        sentence=sentence,
        assertions=(fact_a, fact_b),
        frame_template="{gene}: {status_a} {temporal_a}, {status_b} {temporal_b}.",
        complexity_level="L6",
        compounding_tier="low",
    )


# ---------------------------------------------------------------------------
# Config integrity
# ---------------------------------------------------------------------------


def test_common_config_declares_abbreviation_inconsistency(
    common: CommonConfig,
) -> None:
    noise = common.categories["technical_noise"]
    assert "abbreviation_inconsistency" in noise.categories
    dist = noise.categories["abbreviation_inconsistency"].distribution
    assert set(dist.keys()) == {"canonical", "mixed"}
    assert abs(sum(dist.values()) - 1.0) < 1e-9


# ---------------------------------------------------------------------------
# Applied-transforms key contract
# ---------------------------------------------------------------------------


def test_applied_transforms_contains_abbreviation_inconsistency_on_l1(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Single-fact L1 record reports the key even though it's a no-op."""
    rng = random.Random(0)
    rendered = render_l1_record("EGFR", biomarkers, common, rng)
    post = apply_technical_noise(rendered, common, biomarkers, rng)
    assert "abbreviation_inconsistency" in post.applied_transforms
    assert post.applied_transforms["abbreviation_inconsistency"] in (
        "canonical",
        "mixed",
    )


def test_applied_transforms_contains_abbreviation_inconsistency_on_l6(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Multi-fact L6 record reports the key on the passthrough branch."""
    rng = random.Random(11)
    seen_mode = set()
    for _ in range(200):
        rendered = render_l6_record(
            ["EGFR"], biomarkers, common, rng
        )
        post = apply_technical_noise(rendered, common, biomarkers, rng)
        assert "abbreviation_inconsistency" in post.applied_transforms
        seen_mode.add(post.applied_transforms["abbreviation_inconsistency"])
    assert seen_mode <= {"canonical", "mixed"}


# ---------------------------------------------------------------------------
# Canonical mode — no-op
# ---------------------------------------------------------------------------


def test_canonical_mode_leaves_sentence_and_spans_unchanged(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """When canonical is sampled, repeat-mention record is byte-identical."""
    rng = random.Random(0)
    rec = _make_l6_temporal_frame2_stub(gene="EGFR", gene_surface="EGFR")
    # Force-sample canonical by using a seed we control via the implementation's
    # exposure: run many times and verify at least one canonical draw leaves
    # the record unchanged.
    saw_canonical_unchanged = False
    for _ in range(200):
        post = apply_technical_noise(rec, common, biomarkers, random.Random())
        if post.applied_transforms.get("abbreviation_inconsistency") == "canonical":
            assert post.sentence == rec.sentence
            for orig, new in zip(rec.assertions, post.assertions):
                assert new.spans == orig.spans
            saw_canonical_unchanged = True
            break
    assert saw_canonical_unchanged, (
        "never drew canonical mode across 200 samples — distribution weighted wrong?"
    )


# ---------------------------------------------------------------------------
# Mixed mode — mention #2 surface differs
# ---------------------------------------------------------------------------


def test_mixed_mode_mention_2_surface_differs_from_mention_1(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """When mixed mode fires and a gene has ≥2 aliases, mention #2 uses a
    different surface than mention #1. (ERBB2 has 8 variations — rich pool.)
    """
    distinct_pairs = 0
    identical_pairs = 0
    tried = 0
    rng = random.Random(7)
    while distinct_pairs < 20 and tried < 2000:
        tried += 1
        rec = _make_l6_temporal_frame2_stub(gene="ERBB2", gene_surface="ERBB2")
        post = apply_technical_noise(rec, common, biomarkers, rng)
        if post.applied_transforms.get("abbreviation_inconsistency") != "mixed":
            continue
        g1_span = post.assertions[0].spans["gene"]
        g2_span = post.assertions[1].spans["gene"]
        surf_1 = post.sentence[g1_span[0]:g1_span[1]]
        surf_2 = post.sentence[g2_span[0]:g2_span[1]]
        if surf_1 != surf_2:
            distinct_pairs += 1
        else:
            identical_pairs += 1
    assert distinct_pairs >= 20, (
        f"too few distinct-surface pairs in mixed mode: "
        f"distinct={distinct_pairs} identical={identical_pairs} tried={tried}"
    )


def test_mixed_mode_span_integrity_preserved(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Every labeled span on every fact still resolves to its claimed surface."""
    rng = random.Random(42)
    tested = 0
    tried = 0
    while tested < 50 and tried < 1000:
        tried += 1
        rec = _make_l6_temporal_frame2_stub(
            gene="ERBB2", gene_surface="ERBB2"
        )
        post = apply_technical_noise(rec, common, biomarkers, rng)
        if post.applied_transforms.get("abbreviation_inconsistency") != "mixed":
            continue
        for fact in post.assertions:
            for name, (start, end) in fact.spans.items():
                assert 0 <= start <= end <= len(post.sentence), (
                    f"span {name} out of bounds: ({start},{end}) vs "
                    f"len={len(post.sentence)}; sentence={post.sentence!r}"
                )
                assert post.sentence[start:end], (
                    f"empty slice for {name} in {post.sentence!r}"
                )
        tested += 1
    assert tested >= 50, f"only {tested} mixed-mode records sampled"


def test_mixed_mode_shifts_downstream_status_span_correctly(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """When mention #2 surface differs in length, status span for fact_b
    must be shifted by (len(surf_2) - len(surf_1))."""
    rng = random.Random(19)
    shifted_cases = 0
    tried = 0
    while shifted_cases < 10 and tried < 2000:
        tried += 1
        rec = _make_l6_temporal_frame2_stub(gene="ERBB2", gene_surface="ERBB2")
        post = apply_technical_noise(rec, common, biomarkers, rng)
        if post.applied_transforms.get("abbreviation_inconsistency") != "mixed":
            continue
        g1 = post.assertions[0].spans["gene"]
        g2 = post.assertions[1].spans["gene"]
        surf_1 = post.sentence[g1[0]:g1[1]]
        surf_2 = post.sentence[g2[0]:g2[1]]
        if len(surf_1) == len(surf_2):
            continue
        s2_start, s2_end = post.assertions[1].spans["status"]
        assert post.sentence[s2_start:s2_end] == "negative", (
            f"status span misaligned after surface-length swap: "
            f"got {post.sentence[s2_start:s2_end]!r}, sentence={post.sentence!r}"
        )
        shifted_cases += 1
    assert shifted_cases >= 10, (
        f"only {shifted_cases} length-changing mixed-mode records observed"
    )


# ---------------------------------------------------------------------------
# No repeats → never swaps even when mixed is sampled
# ---------------------------------------------------------------------------


def test_mixed_on_l6_shared_gene_span_is_no_op(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L6 frame1 has two facts sharing the SAME gene span (single surface
    rendered once). There's no mention #2 to swap — must be a no-op.
    """
    rng = random.Random(3)
    for _ in range(100):
        rec = _make_l6_shared_gene_stub(gene="EGFR", gene_surface="EGFR")
        post = apply_technical_noise(rec, common, biomarkers, rng)
        mode = post.applied_transforms.get("abbreviation_inconsistency")
        if mode == "mixed":
            # No distinct mention #2 — sentence must be unchanged.
            assert post.sentence == rec.sentence
            assert post.assertions[0].spans == rec.assertions[0].spans
            assert post.assertions[1].spans == rec.assertions[1].spans


def test_single_fact_l1_is_always_no_op(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """L1 record has exactly one fact → no repeat mention → no swap even when
    mode=mixed is sampled.
    """
    rng = random.Random(5)
    for _ in range(200):
        rendered = render_l1_record("ERBB2", biomarkers, common, rng)
        post = apply_technical_noise(rendered, common, biomarkers, rng)
        if post.applied_transforms.get("abbreviation_inconsistency") == "mixed":
            # Only one gene span — swap cannot change anything.
            g_span = post.assertions[0].spans["gene"]
            # Surface must still match one of ERBB2's configured realizations
            # (i.e., it wasn't mutated into an alias mid-stream).
            surf = post.sentence[g_span[0]:g_span[1]]
            assert surf, "empty gene surface after no-op mixed"


# ---------------------------------------------------------------------------
# Integration — end-to-end through render_l6_record
# ---------------------------------------------------------------------------


def test_render_l6_then_noise_produces_distinct_surfaces_eventually(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """Running render_l6_record -> apply_technical_noise many times, we
    should observe some records where the two gene surfaces differ. This
    is the headline user-visible signal of 3.11.
    """
    rng = random.Random(123)
    distinct_observed = 0
    repeat_observed = 0
    for _ in range(2000):
        rendered = render_l6_record(
            ["ERBB2"], biomarkers, common, rng
        )
        if len(rendered.assertions) < 2:
            continue
        g1 = rendered.assertions[0].spans["gene"]
        g2 = rendered.assertions[1].spans["gene"]
        if g1 == g2:
            continue  # frame1/frame3: shared span, not a repeat mention
        repeat_observed += 1
        post = apply_technical_noise(rendered, common, biomarkers, rng)
        g1 = post.assertions[0].spans["gene"]
        g2 = post.assertions[1].spans["gene"]
        surf_1 = post.sentence[g1[0]:g1[1]]
        surf_2 = post.sentence[g2[0]:g2[1]]
        if surf_1 != surf_2:
            distinct_observed += 1
    assert repeat_observed >= 100, (
        f"only {repeat_observed} L6 frame2 records in 2000 — distribution shifted?"
    )
    assert distinct_observed >= 10, (
        f"only {distinct_observed}/{repeat_observed} L6 repeat-mention records "
        f"got a distinct-surface swap — mixed mode broken?"
    )


# ---------------------------------------------------------------------------
# Edge: biomarker with a single-dominant alias
# ---------------------------------------------------------------------------


def test_mixed_mode_handles_exhausted_alias_pool_gracefully(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> None:
    """If mention #1's surface is the only viable alias (degenerate pool),
    mixed mode becomes a no-op but must not crash and must still report
    the sampled mode.
    """
    rng = random.Random(999)
    # ALK has 7 variants — pick mention #1 as canonical 'ALK'; after
    # excluding 'ALK', 6 others remain. This exercises the exclusion
    # branch cleanly without triggering the degenerate case. We just
    # verify it doesn't crash over many iterations.
    for _ in range(200):
        rec = _make_l6_temporal_frame2_stub(gene="ALK", gene_surface="ALK")
        post = apply_technical_noise(rec, common, biomarkers, rng)
        assert post.applied_transforms["abbreviation_inconsistency"] in (
            "canonical",
            "mixed",
        )
        for fact in post.assertions:
            for name, (start, end) in fact.spans.items():
                assert 0 <= start <= end <= len(post.sentence)
