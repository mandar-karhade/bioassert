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
    compounding_tier: str = "low"
    sentences: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.sentences:
            object.__setattr__(self, "sentences", (self.sentence,))


_ALPHA_ONLY = re.compile(r"^[A-Za-z]+$")
_PUNCT_CHARS = set(",.;:!?")

# OCR substitution pool — visually-confusable char pairs commonly seen in
# scanned documents. All swaps are length-preserving so labeled spans stay
# stable when corruption lands outside them.
OCR_CHAR_SWAPS: dict[str, str] = {
    "l": "1",
    "1": "l",
    "O": "0",
    "0": "O",
    "S": "5",
    "5": "S",
    "B": "8",
    "8": "B",
    "Z": "2",
    "2": "Z",
    "I": "l",
}

# PDF hyphen-linebreak artifact — a hyphen followed by a newline inserted
# inside a non-labeled token simulates a word broken across a PDF page edge.
PDF_HYPHEN_LINEBREAK: str = "-\n"


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
    if len(record.sentences) > 1:
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

    if "ocr_corruption" in noise.categories:
        mode = _sample_mode(
            noise.categories["ocr_corruption"].distribution, rng
        )
        applied["ocr_corruption"] = mode
        if mode != "canonical":
            sentence = _apply_ocr_corruption(sentence, spans, mode, rng)

    if "pdf_artifact" in noise.categories:
        mode = _sample_mode(
            noise.categories["pdf_artifact"].distribution, rng
        )
        applied["pdf_artifact"] = mode
        if mode != "canonical":
            sentence, spans = _apply_pdf_artifact(sentence, spans, mode, rng)

    new_fact = AssertionFact(
        gene=original_fact.gene,
        status=original_fact.status,
        spans=spans,
        variant_id=original_fact.variant_id,
        negative_form_id=original_fact.negative_form_id,
        clone_id=original_fact.clone_id,
        test_method=original_fact.test_method,
        measurement_value=original_fact.measurement_value,
        polarity_scope=original_fact.polarity_scope,
        temporal=original_fact.temporal,
        certainty=original_fact.certainty,
        sentence_index=original_fact.sentence_index,
    )
    return PostProcessedRecord(
        sentence=sentence,
        sentences=(sentence,),
        assertions=(new_fact,),
        frame_template=record.frame_template,
        applied_transforms=applied,
        complexity_level=record.complexity_level,
        compounding_tier=record.compounding_tier,
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
        "ocr_corruption": "skipped",
        "pdf_artifact": "skipped",
    }
    return PostProcessedRecord(
        sentence=record.sentence,
        sentences=record.sentences,
        assertions=record.assertions,
        frame_template=record.frame_template,
        applied_transforms=skipped,
        complexity_level=record.complexity_level,
        compounding_tier=record.compounding_tier,
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


def _slot_boundary_space_positions(
    sentence: str, spans: dict[str, tuple[int, int]]
) -> list[int]:
    """Return space positions that are immediately adjacent to a labeled span.

    A space at index ``i`` qualifies when either ``i == span_end`` for some
    labeled span (the space sits directly after a slot) or
    ``i + 1 == span_start`` (the space sits directly before a slot).
    Prose-interior spaces (between frame connective words like ``was`` and
    ``found to be``) are excluded so Bug 3a cannot reappear.
    """
    boundaries: set[int] = set()
    for start, end in spans.values():
        if start - 1 >= 0 and sentence[start - 1] == " ":
            boundaries.add(start - 1)
        if end < len(sentence) and sentence[end] == " ":
            boundaries.add(end)
    return sorted(boundaries)


def _apply_whitespace(
    sentence: str,
    spans: dict[str, tuple[int, int]],
    mode: str,
    rng: random.Random,
) -> tuple[str, dict[str, tuple[int, int]]]:
    if mode in ("tab", "newline_mid_sentence", "double_space"):
        positions = _slot_boundary_space_positions(sentence, spans)
    else:
        positions = _space_positions_outside_spans(sentence, spans)
    if mode == "no_space_after_punct":
        punct_space = _punct_space_positions(sentence, spans)
        if not punct_space:
            return sentence, spans
        pick = rng.choice(punct_space)
        new_sentence = sentence[:pick] + sentence[pick + 1 :]
        new_spans = _shift_spans(spans, pick, -1)
        return new_sentence, new_spans
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


_OCR_LIGHT_MAX = 2
_OCR_MODERATE_MAX = 5
_PDF_MIN_TOKEN_LEN = 4
_PDF_TOKEN_RE = re.compile(r"[A-Za-z]+")


def _apply_ocr_corruption(
    sentence: str,
    spans: dict[str, tuple[int, int]],
    mode: str,
    rng: random.Random,
) -> str:
    """Apply length-preserving OCR character swaps outside labeled spans.

    ``light`` swaps 1-2 eligible characters; ``moderate`` swaps up to 5.
    ``canonical`` is an identity. Labeled span characters are never touched
    so ``sentence[s:e]`` still resolves to the intended surface.
    """
    if mode == "canonical":
        return sentence
    if mode == "light":
        max_swaps = _OCR_LIGHT_MAX
    elif mode == "moderate":
        max_swaps = _OCR_MODERATE_MAX
    else:
        raise PostProcessError(f"unknown ocr_corruption mode {mode!r}")

    candidates = [
        i
        for i, ch in enumerate(sentence)
        if ch in OCR_CHAR_SWAPS and _outside_any_span(i, spans)
    ]
    if not candidates:
        return sentence
    k = rng.randint(1, min(max_swaps, len(candidates)))
    picks = rng.sample(candidates, k)
    chars = list(sentence)
    for pos in picks:
        chars[pos] = OCR_CHAR_SWAPS[chars[pos]]
    return "".join(chars)


def _apply_pdf_artifact(
    sentence: str,
    spans: dict[str, tuple[int, int]],
    mode: str,
    rng: random.Random,
) -> tuple[str, dict[str, tuple[int, int]]]:
    """Insert a PDF-style hyphenated linebreak inside a non-labeled token.

    ``hyphen_linebreak`` picks a purely-alphabetic word of ≥4 chars that
    does not overlap any labeled span, then splices ``-\\n`` at a random
    interior position. Spans strictly after the cut shift by +2.
    ``canonical`` is an identity.
    """
    if mode == "canonical":
        return sentence, spans
    if mode != "hyphen_linebreak":
        raise PostProcessError(f"unknown pdf_artifact mode {mode!r}")

    candidates: list[int] = []
    for match in _PDF_TOKEN_RE.finditer(sentence):
        start, end = match.start(), match.end()
        if (end - start) < _PDF_MIN_TOKEN_LEN:
            continue
        if any(
            not (end <= sstart or start >= send)
            for sstart, send in spans.values()
        ):
            continue
        candidates.extend(range(start + 1, end))

    if not candidates:
        return sentence, spans

    pick = rng.choice(candidates)
    new_sentence = sentence[:pick] + PDF_HYPHEN_LINEBREAK + sentence[pick:]
    new_spans = _shift_spans(spans, pick, len(PDF_HYPHEN_LINEBREAK))
    return new_sentence, new_spans


__all__ = [
    "OCR_CHAR_SWAPS",
    "PDF_HYPHEN_LINEBREAK",
    "PostProcessedRecord",
    "PostProcessError",
    "apply_technical_noise",
]
