"""Technical-noise post-processing for rendered records (Phase 2a).

Applies the four ``technical_noise`` sub-categories
(``whitespace``, ``case_variation``, ``hyphenation_gene_names``,
``punctuation_variation``) AFTER semantic composition per
CONFIG_ARCHITECTURE.md §7.6. Every transformation updates the labeled
character spans in lockstep with the mutated sentence so Non-Negotiable #1
(spans resolve to literal substrings) continues to hold post-noise.

Design choices worth naming:
 * Case changes are length-preserving; spans untouched.
 * Same-length whitespace swaps (``tab``, ``newline_mid_sentence``) only
   replace a single gap that lies *outside* every labeled span. This avoids
   injecting whitespace inside a labeled token and losing readability.
 * Length-changing transforms (``double_space``, ``no_space_after_punct``,
   ``extra_comma``, ``missing_period``, hyphenation insertions) shift every
   downstream span using :func:`_shift_spans`.
 * Hyphenation only fires on pure-alphabetic gene surfaces with ≥4 chars;
   heterogeneous surfaces (``PD-L1``, ``HER2``) stay canonical because the
   naive middle-insert would produce nonsense.
"""
from __future__ import annotations

import random
import re
from dataclasses import dataclass
from typing import Optional

from bioassert.config.loader import CommonConfig
from bioassert.config.schema import PostProcessTransformations
from bioassert.generator.renderer import AssertionFact, RenderedRecord


@dataclass(frozen=True)
class PostProcessedRecord:
    """A :class:`RenderedRecord` plus the chosen noise decisions.

    ``sentence`` and every ``AssertionFact.spans`` reflect the POST-transformation
    state; the invariant ``sentence[start:end] == labeled_substring`` still holds
    for every entry in every fact's ``spans`` including ``clone`` when attached.
    """

    sentence: str
    assertions: tuple[AssertionFact, ...]
    frame_template: str
    applied_transforms: dict[str, str]
    complexity_level: str = "L1"


_ALPHA_ONLY = re.compile(r"^[A-Za-z]+$")
_PUNCT_CHARS = set(",.;:!?")


class PostProcessError(RuntimeError):
    """Raised when a noise transformation would violate span integrity."""


def apply_technical_noise(
    record: RenderedRecord,
    common: CommonConfig,
    rng: random.Random,
) -> PostProcessedRecord:
    """Run the four noise categories in order and return the new record."""
    noise = common.categories.get("technical_noise")
    if not isinstance(noise, PostProcessTransformations):
        raise KeyError(
            "common.technical_noise is missing or not a "
            "post_process_transformations block"
        )

    if len(record.assertions) != 1:
        return _passthrough_post_process(record)
    original_fact = record.assertions[0]
    sentence = record.sentence
    spans = dict(original_fact.spans)
    applied: dict[str, str] = {}

    if "hyphenation_gene_names" in noise.categories:
        mode = _sample_mode(
            noise.categories["hyphenation_gene_names"].distribution, rng
        )
        applied["hyphenation_gene_names"] = mode
        if mode != "canonical":
            sentence, spans = _apply_hyphenation(sentence, spans, mode)

    if "whitespace" in noise.categories:
        mode = _sample_mode(noise.categories["whitespace"].distribution, rng)
        applied["whitespace"] = mode
        if mode != "single_space":
            sentence, spans = _apply_whitespace(sentence, spans, mode, rng)

    if "case_variation" in noise.categories:
        mode = _sample_mode(noise.categories["case_variation"].distribution, rng)
        applied["case_variation"] = mode
        if mode != "canonical":
            sentence = _apply_case(sentence, mode)

    if "punctuation_variation" in noise.categories:
        mode = _sample_mode(
            noise.categories["punctuation_variation"].distribution, rng
        )
        applied["punctuation_variation"] = mode
        if mode != "canonical":
            sentence, spans = _apply_punctuation(sentence, spans, mode, rng)

    new_fact = AssertionFact(
        gene=original_fact.gene,
        status=original_fact.status,
        spans=spans,
        variant_id=original_fact.variant_id,
        negative_form_id=original_fact.negative_form_id,
        clone_id=original_fact.clone_id,
        test_method=original_fact.test_method,
        measurement_value=original_fact.measurement_value,
    )
    return PostProcessedRecord(
        sentence=sentence,
        assertions=(new_fact,),
        frame_template=record.frame_template,
        applied_transforms=applied,
        complexity_level=record.complexity_level,
    )


def _passthrough_post_process(record: RenderedRecord) -> PostProcessedRecord:
    """Return a PostProcessedRecord that copies ``record`` unchanged.

    Used for multi-assertion records (L3+) until Sub-phase 3.9 expands the
    noise system to coordinate span updates across N facts. Downstream
    consumers see ``applied_transforms`` marked ``"skipped"`` for every
    category so the schema stays consistent.
    """
    skipped = {
        "hyphenation_gene_names": "skipped",
        "whitespace": "skipped",
        "case_variation": "skipped",
        "punctuation_variation": "skipped",
    }
    return PostProcessedRecord(
        sentence=record.sentence,
        assertions=record.assertions,
        frame_template=record.frame_template,
        applied_transforms=skipped,
        complexity_level=record.complexity_level,
    )


def _sample_mode(distribution: dict[str, float], rng: random.Random) -> str:
    keys = list(distribution)
    weights = [distribution[k] for k in keys]
    return rng.choices(keys, weights=weights, k=1)[0]


def _apply_case(sentence: str, mode: str) -> str:
    """Apply casing; every mode preserves character count and spans."""
    if mode == "all_lowercase":
        return sentence.lower()
    if mode == "all_uppercase":
        return sentence.upper()
    if mode == "title_case":
        return sentence.title()
    raise PostProcessError(f"unknown case_variation mode {mode!r}")


def _shift_spans(
    spans: dict[str, tuple[int, int]], insert_at: int, delta: int
) -> dict[str, tuple[int, int]]:
    """Shift spans for an insertion/deletion at position ``insert_at``.

    Positions *strictly after* the cut move by ``delta``. A span whose end
    sits exactly at ``insert_at`` is NOT extended — the mutation happens
    outside it.
    """
    out: dict[str, tuple[int, int]] = {}
    for name, (start, end) in spans.items():
        new_start = start + delta if start >= insert_at else start
        new_end = end + delta if end > insert_at else end
        out[name] = (new_start, new_end)
    return out


def _grow_span(
    spans: dict[str, tuple[int, int]],
    span_name: str,
    insert_at: int,
    delta: int,
) -> dict[str, tuple[int, int]]:
    """Shift everything after ``insert_at`` by ``delta``, extending ``span_name``.

    Used when the insertion lands *inside* the named span (e.g., hyphenating
    ``EGFR`` to ``EGF-R`` grows the gene span by 1).
    """
    out: dict[str, tuple[int, int]] = {}
    target_start, target_end = spans[span_name]
    for name, (start, end) in spans.items():
        if name == span_name:
            out[name] = (target_start, target_end + delta)
            continue
        new_start = start + delta if start >= insert_at else start
        new_end = end + delta if end > insert_at else end
        out[name] = (new_start, new_end)
    return out


def _apply_hyphenation(
    sentence: str,
    spans: dict[str, tuple[int, int]],
    mode: str,
) -> tuple[str, dict[str, tuple[int, int]]]:
    """Hyphenate the gene surface if it's pure-alpha and long enough.

    Falls back to canonical silently when the gene surface can't be safely
    split — keeps noise injection conservative for multi-token surfaces like
    ``PD-L1`` where naive middle-insert would be nonsense.
    """
    if "gene" not in spans:
        return sentence, spans
    start, end = spans["gene"]
    gene_surface = sentence[start:end]
    if len(gene_surface) < 4 or not _ALPHA_ONLY.match(gene_surface):
        return sentence, spans

    split_at = len(gene_surface) // 2
    if mode == "hyphenated":
        insert = "-"
    elif mode == "spaced":
        insert = " "
    elif mode == "linebroken":
        insert = "-\n"
    else:
        raise PostProcessError(f"unknown hyphenation mode {mode!r}")

    new_gene = gene_surface[:split_at] + insert + gene_surface[split_at:]
    new_sentence = sentence[:start] + new_gene + sentence[end:]
    insert_pos = start + split_at
    delta = len(insert)
    new_spans = _grow_span(spans, "gene", insert_pos, delta)
    return new_sentence, new_spans


def _outside_any_span(pos: int, spans: dict[str, tuple[int, int]]) -> bool:
    """True iff ``pos`` is not inside any labeled span content.

    "Inside" means ``start <= pos < end`` (end is exclusive). A position
    equal to a span *end* is outside; a position equal to a span *start*
    is inside, since it refers to the span's first character.
    """
    for start, end in spans.values():
        if start <= pos < end:
            return False
    return True


def _space_positions_outside_spans(
    sentence: str, spans: dict[str, tuple[int, int]]
) -> list[int]:
    return [
        i
        for i, ch in enumerate(sentence)
        if ch == " " and _outside_any_span(i, spans)
    ]


def _apply_whitespace(
    sentence: str,
    spans: dict[str, tuple[int, int]],
    mode: str,
    rng: random.Random,
) -> tuple[str, dict[str, tuple[int, int]]]:
    positions = _space_positions_outside_spans(sentence, spans)
    if not positions:
        return sentence, spans
    pick = rng.choice(positions)

    if mode in ("tab", "newline_mid_sentence"):
        replacement = "\t" if mode == "tab" else "\n"
        new_sentence = sentence[:pick] + replacement + sentence[pick + 1 :]
        return new_sentence, spans

    if mode == "double_space":
        new_sentence = sentence[:pick] + "  " + sentence[pick + 1 :]
        new_spans = _shift_spans(spans, pick + 1, 1)
        return new_sentence, new_spans

    if mode == "no_space_after_punct":
        punct_space = _punct_space_positions(sentence, spans)
        if not punct_space:
            return sentence, spans
        pick = rng.choice(punct_space)
        new_sentence = sentence[:pick] + sentence[pick + 1 :]
        new_spans = _shift_spans(spans, pick, -1)
        return new_sentence, new_spans

    raise PostProcessError(f"unknown whitespace mode {mode!r}")


def _punct_space_positions(
    sentence: str, spans: dict[str, tuple[int, int]]
) -> list[int]:
    positions: list[int] = []
    for i in range(1, len(sentence)):
        if sentence[i] == " " and sentence[i - 1] in _PUNCT_CHARS:
            if _outside_any_span(i, spans):
                positions.append(i)
    return positions


def _apply_punctuation(
    sentence: str,
    spans: dict[str, tuple[int, int]],
    mode: str,
    rng: random.Random,
) -> tuple[str, dict[str, tuple[int, int]]]:
    if mode == "missing_period":
        if sentence.endswith("."):
            period_pos = len(sentence) - 1
            new_spans: dict[str, tuple[int, int]] = {}
            for name, (start, end) in spans.items():
                if end > period_pos:
                    end = period_pos
                new_spans[name] = (start, end)
            return sentence[:-1], new_spans
        return sentence, spans

    if mode == "extra_comma":
        candidates = _space_positions_outside_spans(sentence, spans)
        if not candidates:
            return sentence, spans
        pick = rng.choice(candidates)
        new_sentence = sentence[:pick] + "," + sentence[pick:]
        new_spans = _shift_spans(spans, pick, 1)
        return new_sentence, new_spans

    if mode == "ocr_artifact":
        return _ocr_artifact(sentence, spans, rng)

    raise PostProcessError(f"unknown punctuation mode {mode!r}")


_OCR_SUBS: dict[str, str] = {
    "o": "0",
    "O": "0",
    "l": "1",
    "I": "1",
    "S": "5",
    "B": "8",
}


def _ocr_artifact(
    sentence: str,
    spans: dict[str, tuple[int, int]],
    rng: random.Random,
) -> tuple[str, dict[str, tuple[int, int]]]:
    """Swap one character for its OCR-confusable twin, outside labeled spans.

    Length-preserving, spans untouched. No-op if no eligible character exists
    — the other noise categories are doing the variation work anyway.
    """
    candidates = [
        i
        for i, ch in enumerate(sentence)
        if ch in _OCR_SUBS and _outside_any_span(i, spans)
    ]
    if not candidates:
        return sentence, spans
    pick = rng.choice(candidates)
    swap = _OCR_SUBS[sentence[pick]]
    new_sentence = sentence[:pick] + swap + sentence[pick + 1 :]
    return new_sentence, spans


__all__ = [
    "PostProcessedRecord",
    "PostProcessError",
    "apply_technical_noise",
]
