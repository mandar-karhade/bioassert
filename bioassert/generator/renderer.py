"""Sentence renderer with placeholder substitution and span tracking.

Consumes the probability-weighted config and emits (sentence, spans,
metadata) tuples with char-level labels that satisfy Non-Negotiable #1: the
gene, variant, method, status, and clone spans returned alongside each
sentence are literal substrings at the recorded positions. Placeholders
(``{gene}``, ``{method}``, ``{result}``, ``{value}``) are expanded inside
realization strings before slot substitution per config_architecture.md §7.5.

Phase 2b extends Phase 2a L1 rendering to expression biomarkers (PD-L1, TMB):
positive status draws from ``variants`` with ``{value}`` substitution over the
variant's ``measurement_range``; negative status draws from ``negative_forms``
with ``{value}`` substitution over ``measurement_range_for_value_placeholder``
when present. IHC clone attribution attaches post-frame via an independent
Bernoulli (§8.4).
"""
from __future__ import annotations

import random
import re
from dataclasses import dataclass
from typing import Optional

from bioassert.config.loader import BiomarkerConfig, CommonConfig
from bioassert.config.schema import (
    Biomarker,
    NegativeForms,
    Variant,
    WeightedVariations,
)
from bioassert.config.validator import SUPPORTED_PLACEHOLDERS
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

_PLACEHOLDER_RE = re.compile(r"\{([^{}]+)\}")

# Name forms containing any of these whole-words are excluded from L3 rendering:
# they embed a status or variant-type descriptor ("EGFR mutation",
# "ALK rearrangement", "ALK-positive") which conflicts with L3's "bare gene
# names only" rule. "mutational burden" / "mutational load" survive because
# \bmutation\b does not match "mutational".
_L3_NAME_FORM_BLOCKLIST = re.compile(
    r"\b(mutation|mutations|fusion|fusions|alteration|alterations|"
    r"rearrangement|rearrangements|translocation|translocations|"
    r"amplification|amplifications|expression|positive|negative)\b",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class AssertionFact:
    """One biomarker claim extracted from a rendered sentence.

    Carries the structured truth (gene + status + optional variant/method/clone
    and measurement) plus the character spans that locate every labeled token in
    the rendering. Multiple facts can share the same status span (L3 shared
    status) or occupy disjoint regions of the same sentence (L4 heterogeneous).
    Invariant: ``sentence[s:e] == labeled_substring`` for every ``(s,e)`` in
    ``spans``.
    """

    gene: str
    status: str
    spans: dict[str, tuple[int, int]]
    variant_id: Optional[str] = None
    negative_form_id: Optional[str] = None
    clone_id: Optional[str] = None
    test_method: Optional[str] = None
    measurement_value: Optional[float] = None
    polarity_scope: str = "direct"
    temporal: Optional[str] = None
    certainty: Optional[str] = None
    sentence_index: int = 0


@dataclass(frozen=True)
class RenderedRecord:
    """One rendered sentence + the assertion facts it carries.

    L1/L2 records have exactly one fact; L3+ carries a tuple of facts
    whose order follows the gene order in the surface sentence.
    ``compounding_tier`` is orthogonal to ``complexity_level`` — it
    parameterizes how many genes participate in a compound record
    (``"low"`` = 2 genes, ``"high"`` = 3-8 genes clamped to pool size).
    L1/L2 single-gene records inherit the default ``"low"`` but the
    tier is only meaningful at L3+.

    L7 records span multiple sentences. ``sentences`` carries the
    per-sentence tuple and each ``AssertionFact.sentence_index`` points
    into it. For L1-L6 (single sentence), ``sentences`` defaults to
    ``(sentence,)``.
    """

    sentence: str
    assertions: tuple[AssertionFact, ...]
    frame_template: str
    complexity_level: str = "L1"
    compounding_tier: str = "low"
    sentences: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.sentences:
            object.__setattr__(self, "sentences", (self.sentence,))


_L1_FRAMES: list[dict] = [
    {"template": "{gene} was {status}.", "requires_variant": False, "has_method": False},
    {"template": "{gene} is {status}.", "requires_variant": False, "has_method": False},
    {"template": "{gene}: {status}.", "requires_variant": False, "has_method": False},
    {"template": "{gene} — {status}.", "requires_variant": False, "has_method": False},
    {"template": "{gene} testing: {status}.", "requires_variant": False, "has_method": False},
    {"template": "Result for {gene}: {status}.", "requires_variant": False, "has_method": False},
    {"template": "{gene} status: {status}.", "requires_variant": False, "has_method": False},
    {"template": "{gene} was found to be {status}.", "requires_variant": False, "has_method": False},
    {"template": "{gene} {variant} was {status}.", "requires_variant": True, "has_method": False},
    {"template": "{gene} {variant} {status}.", "requires_variant": True, "has_method": False},
    {"template": "{gene} was {status} by {method}.", "requires_variant": False, "has_method": True},
    {"template": "{gene} {status} ({method}).", "requires_variant": False, "has_method": True},
    {"template": "{method}: {gene} {status}.", "requires_variant": False, "has_method": True},
    {"template": "{gene} {variant} {status} by {method}.", "requires_variant": True, "has_method": True},
    {"template": "{gene} {variant} was {status} on {method}.", "requires_variant": True, "has_method": True},
]

_L2_FRAMES: list[dict] = [
    {"template": "{gene} {status}", "requires_variant": False, "has_method": False},
    {"template": "{gene} {status}.", "requires_variant": False, "has_method": False},
    {"template": "{gene}: {status}", "requires_variant": False, "has_method": False},
    {"template": "{gene}: {status}.", "requires_variant": False, "has_method": False},
    {"template": "{gene} - {status}", "requires_variant": False, "has_method": False},
    {"template": "{gene} [{status}]", "requires_variant": False, "has_method": False},
    {"template": "{gene} ({status})", "requires_variant": False, "has_method": False},
    {"template": "{gene}\t{status}", "requires_variant": False, "has_method": False},
]

_L3_FRAMES: list[dict] = [
    {"template": "{gene_list} were {status}."},
    {"template": "{gene_list} are {status}."},
    {"template": "{gene_list}: {status}."},
    {"template": "{gene_list} - {status}."},
    {"template": "Results for {gene_list}: {status}."},
    {"template": "{gene_list} showed {status}."},
    {"template": "{gene_list} {status}."},
]

# Shorthand/tabular L3 frames. Each gene repeats the status surface, glued
# by ``inner`` (between gene and its status) and ``sep`` (between pairs).
# Every gene-status pair gets its own gene span AND its own status span —
# unlike L3 prose where all facts share a single status span.
_L3_SHORTHAND_FRAMES: list[dict] = [
    {"inner": " ", "sep": ", "},
    {"inner": " ", "sep": "; "},
    {"inner": " ", "sep": " | "},
    {"inner": " ", "sep": "\n"},
    {"inner": ": ", "sep": ", "},
    {"inner": ": ", "sep": "\n"},
    {"inner": "\t", "sep": "\n"},
]

# Heterogeneous L4 frames. Same list structure as L3 but each gene's
# status is drawn independently from its own population, rendered with
# the formal positive_phrases / negation_phrases vocabulary. ``final``
# optionally inserts a coordinator (``"and"`` / ``"or"``) before the
# last pair; ``final_is_oxford`` controls whether a pre-conjunction
# comma is kept when the separator is ``", "``.
_L4_FRAMES: list[dict] = [
    {"inner": " ", "sep": ", ", "final": "and", "final_is_oxford": True},
    {"inner": " ", "sep": ", ", "final": "and", "final_is_oxford": False},
    {"inner": " ", "sep": ", ", "final": None, "final_is_oxford": False},
    {"inner": " ", "sep": "; ", "final": None, "final_is_oxford": False},
    {"inner": " was ", "sep": ", ", "final": "and", "final_is_oxford": True},
    {"inner": " was ", "sep": "; ", "final": None, "final_is_oxford": False},
    {"inner": ": ", "sep": ", ", "final": None, "final_is_oxford": False},
    {"inner": ": ", "sep": "; ", "final": None, "final_is_oxford": False},
    {"inner": " ", "sep": ". ", "final": None, "final_is_oxford": False},
]

# L4 shorthand/tabular frames: per-gene independent statuses rendered
# as compact pairs. Same separator/inner matrix as L3 shorthand.
_L4_SHORTHAND_FRAMES: list[dict] = [
    {"inner": " ", "sep": ", "},
    {"inner": " ", "sep": "; "},
    {"inner": " ", "sep": " | "},
    {"inner": " ", "sep": "\n"},
    {"inner": ": ", "sep": ", "},
    {"inner": ": ", "sep": "\n"},
    {"inner": "\t", "sep": "\n"},
]

# L5 negation-scope frames. Two structural ``kind`` values:
#
# - ``"enumerated"``: scope marker ("No" / "Absence of") negates an explicit
#   list of genes; optional exception clause introduces one trailing gene.
#   Clinical register drops "but"/"however" — the natural exception markers
#   are "and" (continuation, same negative polarity) and "except" (polarity
#   flip, positive). ``exception_marker=None`` is the dominant case.
#
# - ``"panel_wide"``: implicit "panel" as the scope (no enumerated genes);
#   one trailing gene surfaces through "other than" as the sole labeled
#   fact with ``polarity_scope="exception"``. No wide-scope facts emitted
#   because the wide scope is unnamed.
_L5_FRAMES: list[dict] = [
    # 6 enumerated None frames — dominant (~60%).
    {"kind": "enumerated", "scope_marker": "No", "neg_noun": "mutations", "exception_marker": None, "exception_status": None},
    {"kind": "enumerated", "scope_marker": "No", "neg_noun": "alterations", "exception_marker": None, "exception_status": None},
    {"kind": "enumerated", "scope_marker": "No", "neg_noun": "mutations detected", "exception_marker": None, "exception_status": None},
    {"kind": "enumerated", "scope_marker": "No", "neg_noun": "alterations detected", "exception_marker": None, "exception_status": None},
    {"kind": "enumerated", "scope_marker": "Absence of", "neg_noun": "mutations", "exception_marker": None, "exception_status": None},
    {"kind": "enumerated", "scope_marker": "Absence of", "neg_noun": "alterations", "exception_marker": None, "exception_status": None},
    # 2 enumerated "and" frames — continuation, same negative polarity.
    {"kind": "enumerated", "scope_marker": "No", "neg_noun": "mutations", "exception_marker": "and", "exception_status": "negative"},
    {"kind": "enumerated", "scope_marker": "No", "neg_noun": "alterations", "exception_marker": "and", "exception_status": "negative"},
    # 2 enumerated "except" frames — polarity flip, positive.
    {"kind": "enumerated", "scope_marker": "No", "neg_noun": "mutations", "exception_marker": "except", "exception_status": "positive"},
    {"kind": "enumerated", "scope_marker": "Absence of", "neg_noun": "alterations", "exception_marker": "except", "exception_status": "positive"},
    # 4 panel_wide frames — implicit panel, "other than" marks the sole
    # labeled fact. ``status_word`` is the literal positive/negative token
    # embedded in the scope phrase and owns the wide polarity reading.
    {"kind": "panel_wide", "scope": "No biomarker on the panel was", "status_word": "positive", "exception_marker": "other than"},
    {"kind": "panel_wide", "scope": "No gene on this panel was", "status_word": "positive", "exception_marker": "other than"},
    {"kind": "panel_wide", "scope": "No panel biomarker was", "status_word": "positive", "exception_marker": "other than"},
    {"kind": "panel_wide", "scope": "No panel gene tested", "status_word": "positive", "exception_marker": "other than"},
]

L5_MIN_WIDE = 2
L5_MAX_WIDE = 4

L3_MIN_N = 2
L3_MAX_N = 4
L3_OXFORD_PROB = 0.6

# Compounding tier bounds — controls how many genes participate in an
# L3+ compound record. Orthogonal to complexity_level (which picks the
# frame family). Per user direction 2026-04-23: only two tiers. "low"
# pins to the minimum compound size (< 3 genes, i.e., exactly 2);
# "high" covers 3+ genes (3-8, clamped to pool size), typically paired
# with a cross-class pool + multi-sentence emergence via L4 ". " frames.
COMPOUND_LOW_N = 2
COMPOUND_HIGH_MIN = 3
COMPOUND_HIGH_MAX = 8

_TIER_RANGES: dict[str, tuple[int, int]] = {
    "low": (COMPOUND_LOW_N, COMPOUND_LOW_N),
    "high": (COMPOUND_HIGH_MIN, COMPOUND_HIGH_MAX),
}


def n_range_for_tier(tier: str) -> tuple[int, int]:
    """Return the ``(min_n, max_n)`` gene-count range for ``tier``.

    Raises :class:`RenderError` when ``tier`` is not ``"low"`` or
    ``"high"``.
    """
    if tier not in _TIER_RANGES:
        raise RenderError(
            f"compounding_tier must be one of "
            f"{sorted(_TIER_RANGES)}, got {tier!r}"
        )
    return _TIER_RANGES[tier]

VARIANT_FRAME_PROB = 0.55
METHOD_FRAME_PROB = 0.35

_STATUS_CATEGORY: dict[str, str] = {
    "positive": "positive_phrases",
    "negative": "negation_phrases",
    "equivocal": "equivocal_phrases",
    "not_tested": "not_tested_phrases",
}

_SHORTHAND_STATUS_CATEGORY: dict[str, str] = {
    "positive": "positive_shorthand",
    "negative": "negation_shorthand",
}


class RenderError(RuntimeError):
    """Raised when renderer preconditions fail (frame/slot mismatch, missing
    placeholder context, no compatible status phrases). Loud failure is
    preferred over silent span corruption — see Non-Negotiable #1.
    """


def expand_placeholders(
    surface: str, context: dict[str, str]
) -> str:
    """Fill ``{name}`` placeholders in ``surface`` using ``context``.

    Raises :class:`KeyError` on any unfilled placeholder and
    :class:`RenderError` on unknown placeholder names. Unused context keys
    are silently allowed — callers supply the superset of supported names.
    """
    unknown = _extract_placeholders(surface) - SUPPORTED_PLACEHOLDERS
    if unknown:
        raise RenderError(
            f"unknown placeholder(s) {sorted(unknown)} in realization "
            f"{surface!r}"
        )
    try:
        return surface.format_map(context)
    except KeyError as exc:
        raise KeyError(
            f"placeholder {exc.args[0]!r} in {surface!r} has no context value "
            f"(available: {sorted(context)})"
        ) from exc


def _extract_placeholders(s: str) -> set[str]:
    return set(_PLACEHOLDER_RE.findall(s))


def _filter_gene_free(weighted: WeightedVariations) -> dict[str, float]:
    """Return a renormalized variations dict keeping only realizations that
    do not embed the ``{gene}`` placeholder. Phase 2a grammar frames always
    include an explicit ``{gene}`` slot, so status phrases that already embed
    the gene would produce redundant output like "EGFR: no EGFR alteration
    identified." Those realizations are available for future standalone
    frames in Phase 3+.
    """
    eligible: dict[str, float] = {}
    for vid, weight in weighted.variations.items():
        surface = weighted.realizations[vid]
        if "{gene}" in surface:
            continue
        eligible[vid] = weight
    if not eligible:
        raise RenderError(
            "no gene-free realizations available in this status-phrase category"
        )
    total = sum(eligible.values())
    return {k: v / total for k, v in eligible.items()}


def _sample_status_phrase(
    status: str,
    common: CommonConfig,
    rng: random.Random,
    complexity_level: str = "L1",
) -> tuple[str, str]:
    """Return ``(variation_id, surface_string)`` for the sampled status phrase.

    L2 draws positive/negative from the compact ``*_shorthand`` categories;
    equivocal and not_tested always draw from the formal phrase categories
    (no shorthand vocabulary exists for those statuses). L1 always uses the
    formal categories.
    """
    category_name: Optional[str] = None
    if complexity_level == "L2":
        category_name = _SHORTHAND_STATUS_CATEGORY.get(status)
    if category_name is None:
        category_name = _STATUS_CATEGORY[status]
    category = common.categories.get(category_name)
    if not isinstance(category, WeightedVariations):
        raise RenderError(
            f"common.{category_name} is missing or wrong schema_type"
        )
    eligible = _filter_gene_free(category)
    keys = list(eligible)
    weights = [eligible[k] for k in keys]
    vid = rng.choices(keys, weights=weights, k=1)[0]
    return vid, category.realizations[vid]


def _pick_frame(
    rng: random.Random,
    variant_sampled: bool,
    method_sampled: bool,
    complexity_level: str = "L1",
) -> dict:
    """Pick a template consistent with the sampled pieces.

    L2 frames are all no-variant/no-method compact forms; ``variant_sampled``
    and ``method_sampled`` are ignored for L2 because L2 targets tabular
    shorthand where descriptor/method decoration is suppressed.
    """
    if complexity_level == "L2":
        return rng.choice(_L2_FRAMES)
    candidates = _L1_FRAMES
    want_variant = variant_sampled and rng.random() < VARIANT_FRAME_PROB
    want_method = method_sampled and rng.random() < METHOD_FRAME_PROB

    def compatible(frame: dict) -> bool:
        if frame["requires_variant"] and not want_variant:
            return False
        if not frame["requires_variant"] and want_variant:
            return False
        if frame["has_method"] and not want_method:
            return False
        if not frame["has_method"] and want_method:
            return False
        return True

    matching = [f for f in candidates if compatible(f)]
    if not matching:
        want_variant = variant_sampled
        want_method = False
        matching = [
            f
            for f in candidates
            if f["requires_variant"] == want_variant and not f["has_method"]
        ]
    if not matching:
        raise RenderError(
            f"no compatible frames for variant_sampled={variant_sampled}, "
            f"method_sampled={method_sampled}"
        )
    return rng.choice(matching)


def _render_with_spans(
    template: str, slots: dict[str, str]
) -> tuple[str, dict[str, tuple[int, int]]]:
    out: list[str] = []
    spans: dict[str, tuple[int, int]] = {}
    pos = 0
    i = 0
    while i < len(template):
        ch = template[i]
        if ch == "{":
            end = template.index("}", i)
            slot_name = template[i + 1 : end]
            if slot_name not in slots:
                raise RenderError(
                    f"template slot {slot_name!r} missing from slots "
                    f"{sorted(slots)}"
                )
            value = slots[slot_name]
            spans[slot_name] = (pos, pos + len(value))
            out.append(value)
            pos += len(value)
            i = end + 1
        else:
            out.append(ch)
            pos += 1
            i += 1
    return "".join(out), spans


def _maybe_sample_method(
    biomarker: Biomarker,
    common: CommonConfig,
    rng: random.Random,
    attach_prob: float,
) -> Optional[tuple[str, str]]:
    """With probability ``attach_prob``, sample a (method_id, method_realized)
    pair. Returns ``None`` when no method should attach or when the sampled ID
    is the ``unspecified`` sentinel.
    """
    if rng.random() >= attach_prob:
        return None
    method_id = sample_method(biomarker, common, rng)
    if method_id == "unspecified":
        return None
    realization = common.test_method_realization(method_id)
    if not realization:
        return None
    return method_id, realization


def render_l1_record(
    biomarker_name: str,
    profile: PatientProfile,
    biomarkers: BiomarkerConfig,
    common: CommonConfig,
    rng: random.Random,
    method_attach_prob: float = 0.5,
    complexity_level: str = "L1",
) -> RenderedRecord:
    """Sample + render one record for ``biomarker_name`` under ``profile``.

    Flow mirrors config_architecture.md §3. Supports mutation, fusion,
    composite, and expression biomarkers. For expression biomarkers the
    descriptor slot draws from ``variants`` when status is positive and from
    ``negative_forms`` when status is negative; equivocal and not_tested leave
    the descriptor empty. Clone attribution (§8.4) is sampled independently and
    appended post-frame when attached.

    ``complexity_level`` controls surface form: ``"L1"`` produces formal prose
    frames ("EGFR was positive by NGS."); ``"L2"`` produces tabular shorthand
    ("EGFR: +"). L2 suppresses variant descriptor and method attachment; the
    shorthand vocabulary is drawn from ``positive_shorthand`` /
    ``negation_shorthand`` for positive/negative statuses.
    """
    if complexity_level not in ("L1", "L2"):
        raise RenderError(
            f"complexity_level must be 'L1' or 'L2', got {complexity_level!r}"
        )
    biomarker = biomarkers.get(biomarker_name)
    population = resolve_population(biomarker, profile)
    status = sample_status(population, rng)

    variant: Optional[Variant] = None
    variant_id: Optional[str] = None
    negative_form_id: Optional[str] = None
    measurement_value: Optional[float] = None

    if status == "positive":
        variant_id, variant = sample_variant(biomarker, rng)
        if variant.measurement_range is not None:
            measurement_value = sample_measurement_value(
                variant.measurement_range, rng
            )

    gene_name_form_id = sample_biomarker_name_form(biomarker, rng, variant)
    gene_surface = biomarker.name_forms.realizations[gene_name_form_id]

    descriptor_surface: Optional[str] = None
    if complexity_level == "L1" and variant is not None and variant_id is not None:
        vf_id = sample_variation(variant.name_forms, rng)
        raw_variant = variant.name_forms.realizations[vf_id]
        if raw_variant:
            ctx: dict[str, str] = {"gene": gene_surface}
            if measurement_value is not None:
                ctx["value"] = _format_measurement(measurement_value)
            descriptor_surface = expand_placeholders(raw_variant, ctx)
    elif (
        complexity_level == "L1"
        and status == "negative"
        and biomarker.alteration_type == "expression"
        and biomarker.negative_forms is not None
    ):
        nf_id, nf_surface, nf_value = _sample_negative_form(
            biomarker.negative_forms, gene_surface, rng
        )
        if nf_surface:
            negative_form_id = nf_id
            descriptor_surface = nf_surface
            if nf_value is not None:
                measurement_value = nf_value

    method_id: Optional[str] = None
    method_surface: Optional[str] = None
    if complexity_level == "L1":
        method_info = _maybe_sample_method(
            biomarker, common, rng, attach_prob=method_attach_prob
        )
        if method_info is not None:
            method_id, method_surface = method_info

    _, status_surface_raw = _sample_status_phrase(
        status, common, rng, complexity_level=complexity_level
    )
    status_surface = expand_placeholders(
        status_surface_raw, {"gene": gene_surface}
    )

    frame = _pick_frame(
        rng,
        variant_sampled=descriptor_surface is not None,
        method_sampled=method_surface is not None,
        complexity_level=complexity_level,
    )

    slots: dict[str, str] = {"gene": gene_surface, "status": status_surface}
    if frame["requires_variant"]:
        if descriptor_surface is None:
            raise RenderError(
                "frame requires descriptor but none was rendered"
            )
        slots["variant"] = descriptor_surface
    if frame["has_method"]:
        if method_surface is None:
            raise RenderError(
                "frame requires method but none was sampled"
            )
        slots["method"] = method_surface

    sentence, spans = _render_with_spans(frame["template"], slots)

    clone_id = maybe_sample_clone(biomarker, rng)
    clone_surface: Optional[str] = None
    if clone_id is not None and biomarker.clone_attribution is not None:
        clone_surface = biomarker.clone_attribution.realizations[clone_id]
        sentence, spans = _attach_clone(sentence, spans, clone_surface)

    fact = AssertionFact(
        gene=biomarker_name,
        status=status,
        spans=spans,
        variant_id=variant_id,
        negative_form_id=negative_form_id,
        clone_id=clone_id,
        test_method=method_id,
        measurement_value=measurement_value,
    )
    return RenderedRecord(
        sentence=sentence,
        assertions=(fact,),
        frame_template=frame["template"],
        complexity_level=complexity_level,
    )


def _sample_negative_form(
    neg: NegativeForms,
    gene_surface: str,
    rng: random.Random,
) -> tuple[str, Optional[str], Optional[float]]:
    """Sample one ``negative_forms`` entry and expand placeholders.

    Returns ``(variation_id, rendered_surface, value_used)``. ``rendered_surface``
    may be empty string — callers should treat empty as "no descriptor rendered".
    """
    keys = list(neg.variations)
    weights = [neg.variations[k] for k in keys]
    nf_id = rng.choices(keys, weights=weights, k=1)[0]
    raw = neg.realizations[nf_id]
    if not raw:
        return nf_id, "", None

    ctx: dict[str, str] = {"gene": gene_surface}
    value_used: Optional[float] = None
    if "{value}" in raw:
        if neg.measurement_range_for_value_placeholder is None:
            raise RenderError(
                f"negative_form {nf_id!r} uses {{value}} but no "
                "measurement_range_for_value_placeholder is configured"
            )
        value_used = sample_measurement_value(
            neg.measurement_range_for_value_placeholder, rng
        )
        ctx["value"] = _format_measurement(value_used)
    return nf_id, expand_placeholders(raw, ctx), value_used


def _attach_clone(
    sentence: str,
    spans: dict[str, tuple[int, int]],
    clone_surface: str,
) -> tuple[str, dict[str, tuple[int, int]]]:
    """Append ``clone_surface`` before terminal punctuation with a leading space.

    Registers a ``clone`` span pointing at the clone substring and leaves every
    other span untouched (the insertion lands at ``len(sentence)`` or just
    before a terminal period, both of which are after every existing span's
    ``end``).
    """
    if sentence.endswith("."):
        insert_at = len(sentence) - 1
    else:
        insert_at = len(sentence)
    new_sentence = sentence[:insert_at] + " " + clone_surface + sentence[insert_at:]
    clone_start = insert_at + 1
    clone_end = clone_start + len(clone_surface)
    new_spans = dict(spans)
    new_spans["clone"] = (clone_start, clone_end)
    return new_sentence, new_spans


def _format_measurement(value: float) -> str:
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.1f}"


def _bare_gene_name_form(
    biomarker: Biomarker, rng: random.Random
) -> str:
    """Sample a name_form whose surface is a bare gene identifier.

    Excludes realizations that embed alteration-type or status keywords
    (mutation, fusion, rearrangement, positive, ...). L3 frames carry the
    status in a separate slot, so forms like ``"ALK-positive"`` or
    ``"KRAS mutation"`` would double-mark the assertion. When every
    realization is blocklisted (defensive fallback), falls back to the
    full distribution.
    """
    eligible: dict[str, float] = {}
    for vid, weight in biomarker.name_forms.variations.items():
        surface = biomarker.name_forms.realizations[vid]
        if _L3_NAME_FORM_BLOCKLIST.search(surface):
            continue
        eligible[vid] = weight
    if not eligible:
        return sample_biomarker_name_form(biomarker, rng, None)
    keys = list(eligible)
    weights = [eligible[k] for k in keys]
    return rng.choices(keys, weights=weights, k=1)[0]


def _coordinate_gene_list(
    gene_surfaces: list[str],
    rng: random.Random,
    oxford_prob: float = L3_OXFORD_PROB,
) -> tuple[str, list[tuple[int, int]]]:
    """Join ``gene_surfaces`` into an English list and return per-gene spans.

    Two genes → ``"A and B"``. Three or more → ``"A, B, and C"`` (Oxford) or
    ``"A, B and C"`` (no Oxford), picked by Bernoulli on ``oxford_prob``. The
    returned positions are ``(start, end)`` tuples into the joined string for
    each input surface, in order.
    """
    if len(gene_surfaces) < 2:
        raise RenderError(
            f"gene list needs at least 2 surfaces, got {len(gene_surfaces)}"
        )
    if len(gene_surfaces) == 2:
        joined = f"{gene_surfaces[0]} and {gene_surfaces[1]}"
        first_end = len(gene_surfaces[0])
        second_start = first_end + len(" and ")
        return joined, [
            (0, first_end),
            (second_start, second_start + len(gene_surfaces[1])),
        ]

    use_oxford = rng.random() < oxford_prob
    prefix_parts = gene_surfaces[:-1]
    last = gene_surfaces[-1]
    prefix = ", ".join(prefix_parts)
    sep = ", and " if use_oxford else " and "
    joined = prefix + sep + last

    positions: list[tuple[int, int]] = []
    pos = 0
    for i, surface in enumerate(prefix_parts):
        positions.append((pos, pos + len(surface)))
        pos += len(surface)
        if i < len(prefix_parts) - 1:
            pos += len(", ")
    pos += len(sep)
    positions.append((pos, pos + len(last)))
    return joined, positions


def render_l3_record(
    biomarker_pool: list[str],
    profile: PatientProfile,
    biomarkers: BiomarkerConfig,
    common: CommonConfig,
    rng: random.Random,
    n: Optional[int] = None,
    complexity_level: str = "L3",
    compounding_tier: str = "low",
) -> RenderedRecord:
    """Render one L3 record: N genes share one status.

    Samples N distinct genes from ``biomarker_pool`` (N range controlled by
    ``compounding_tier`` via :func:`n_range_for_tier`, capped at pool size),
    draws one shared status from the first gene's resolved population, and
    renders either a prose frame (``complexity_level="L3"``) or a shorthand/
    tabular frame (``complexity_level="L3S"``).

    Prose: all facts share one status span (coordinator gene list + single
    formal status phrase). Shorthand: each fact gets its own status span
    because the shorthand surface is repeated per gene-status pair.

    L3 never emits variant descriptors or method slots — those land in L4/L5.
    The caller is responsible for keeping the pool homogeneous
    (mutation-class vs expression-class) for medium/low tiers; high tier is
    explicitly allowed to mix classes. The renderer does not enforce this.
    """
    if complexity_level not in ("L3", "L3S"):
        raise RenderError(
            f"complexity_level must be 'L3' or 'L3S', got {complexity_level!r}"
        )
    tier_min, tier_max = n_range_for_tier(compounding_tier)
    pool_size = len(biomarker_pool)
    if pool_size < L3_MIN_N:
        raise RenderError(
            f"L3 pool needs >= {L3_MIN_N} biomarkers, got {pool_size}"
        )
    min_n = min(tier_min, pool_size)
    max_n = min(tier_max, pool_size)
    if n is None:
        n = rng.randint(min_n, max_n)
    elif not min_n <= n <= max_n:
        raise RenderError(f"n must be in [{min_n}, {max_n}], got {n}")

    gene_names: list[str] = rng.sample(biomarker_pool, n)
    first_biomarker = biomarkers.get(gene_names[0])
    population = resolve_population(first_biomarker, profile)
    status = sample_status(population, rng)

    gene_surfaces: list[str] = []
    for name in gene_names:
        b = biomarkers.get(name)
        form_id = _bare_gene_name_form(b, rng)
        gene_surfaces.append(b.name_forms.realizations[form_id])

    if complexity_level == "L3S":
        return _render_l3_shorthand(
            gene_names, gene_surfaces, status, common, rng,
            compounding_tier=compounding_tier,
        )

    gene_list_surface, gene_positions_in_list = _coordinate_gene_list(
        gene_surfaces, rng
    )

    _, status_surface_raw = _sample_status_phrase(
        status, common, rng, complexity_level="L1"
    )
    status_surface = expand_placeholders(
        status_surface_raw, {"gene": gene_surfaces[0]}
    )

    frame = rng.choice(_L3_FRAMES)
    sentence, template_spans = _render_with_spans(
        frame["template"],
        {"gene_list": gene_list_surface, "status": status_surface},
    )
    gene_list_start = template_spans["gene_list"][0]
    status_span = template_spans["status"]

    assertions: list[AssertionFact] = []
    for name, (s, e) in zip(gene_names, gene_positions_in_list):
        fact_spans: dict[str, tuple[int, int]] = {
            "gene": (gene_list_start + s, gene_list_start + e),
            "status": status_span,
        }
        assertions.append(
            AssertionFact(
                gene=name,
                status=status,
                spans=fact_spans,
            )
        )

    return RenderedRecord(
        sentence=sentence,
        assertions=tuple(assertions),
        frame_template=frame["template"],
        complexity_level="L3",
        compounding_tier=compounding_tier,
    )


def _render_l3_shorthand(
    gene_names: list[str],
    gene_surfaces: list[str],
    status: str,
    common: CommonConfig,
    rng: random.Random,
    compounding_tier: str = "low",
) -> RenderedRecord:
    """Build the shorthand/tabular L3 surface with distributive status spans.

    The shorthand status surface (e.g. ``"-"``, ``"neg"``, ``"+"``) is drawn
    once and repeated per gene. Each returned :class:`AssertionFact` owns its
    own ``(gene, status)`` span pair. When the sampled ``status`` has no
    shorthand vocabulary (equivocal / not_tested), the formal phrase category
    is used — the resulting surface is longer but the pair structure is
    preserved.
    """
    _, status_surface_raw = _sample_status_phrase(
        status, common, rng, complexity_level="L2"
    )
    status_surface = expand_placeholders(
        status_surface_raw, {"gene": gene_surfaces[0]}
    )

    frame = rng.choice(_L3_SHORTHAND_FRAMES)
    inner = frame["inner"]
    sep = frame["sep"]

    pieces: list[str] = []
    assertions: list[AssertionFact] = []
    pos = 0
    for i, (name, gene_surface) in enumerate(zip(gene_names, gene_surfaces)):
        if i > 0:
            pieces.append(sep)
            pos += len(sep)
        gene_start = pos
        pieces.append(gene_surface)
        pos += len(gene_surface)
        gene_end = pos
        pieces.append(inner)
        pos += len(inner)
        status_start = pos
        pieces.append(status_surface)
        pos += len(status_surface)
        status_end = pos
        assertions.append(
            AssertionFact(
                gene=name,
                status=status,
                spans={
                    "gene": (gene_start, gene_end),
                    "status": (status_start, status_end),
                },
            )
        )

    sentence = "".join(pieces)
    frame_template = f"{{gene}}{inner}{{status}}{sep}"
    return RenderedRecord(
        sentence=sentence,
        assertions=tuple(assertions),
        frame_template=frame_template,
        complexity_level="L3S",
        compounding_tier=compounding_tier,
    )


def render_l4_record(
    biomarker_pool: list[str],
    profile: PatientProfile,
    biomarkers: BiomarkerConfig,
    common: CommonConfig,
    rng: random.Random,
    n: Optional[int] = None,
    complexity_level: str = "L4",
    compounding_tier: str = "low",
) -> RenderedRecord:
    """Render one L4 record: N genes, each with its own independently
    drawn status.

    Samples N distinct genes (N range from :func:`n_range_for_tier` using
    ``compounding_tier``, capped at pool size), draws one status per gene
    from the gene's resolved population, and glues the pairs into a single
    surface. Each :class:`AssertionFact` carries its own gene span AND its
    own status span — statuses may diverge.

    ``complexity_level="L4"`` renders formal prose frames (:data:`_L4_FRAMES`
    — ``positive_phrases``/``negation_phrases`` vocab + optional Oxford
    ``and``). ``complexity_level="L4S"`` renders shorthand/tabular frames
    (:data:`_L4_SHORTHAND_FRAMES` — ``positive_shorthand``/
    ``negation_shorthand`` vocab, compact separators, no coordinator).

    Variant descriptors and method slots are deferred to L5+. High tier
    ("high") is explicitly allowed to mix classes (mutation + expression)
    in the pool; medium/low tiers should be kept homogeneous by the caller.
    """
    if complexity_level not in ("L4", "L4S"):
        raise RenderError(
            f"complexity_level must be 'L4' or 'L4S', got {complexity_level!r}"
        )
    tier_min, tier_max = n_range_for_tier(compounding_tier)
    pool_size = len(biomarker_pool)
    if pool_size < L3_MIN_N:
        raise RenderError(
            f"L4 pool needs >= {L3_MIN_N} biomarkers, got {pool_size}"
        )
    min_n = min(tier_min, pool_size)
    max_n = min(tier_max, pool_size)
    if n is None:
        n = rng.randint(min_n, max_n)
    elif not min_n <= n <= max_n:
        raise RenderError(f"n must be in [{min_n}, {max_n}], got {n}")

    gene_names: list[str] = rng.sample(biomarker_pool, n)
    status_phrase_level = "L1" if complexity_level == "L4" else "L2"
    gene_surfaces: list[str] = []
    statuses: list[str] = []
    status_surfaces: list[str] = []

    for name in gene_names:
        b = biomarkers.get(name)
        population = resolve_population(b, profile)
        status = sample_status(population, rng)
        form_id = _bare_gene_name_form(b, rng)
        gene_surface = b.name_forms.realizations[form_id]
        _, raw = _sample_status_phrase(
            status, common, rng, complexity_level=status_phrase_level
        )
        status_surface = expand_placeholders(raw, {"gene": gene_surface})
        gene_surfaces.append(gene_surface)
        statuses.append(status)
        status_surfaces.append(status_surface)

    if complexity_level == "L4S":
        return _render_l4_shorthand(
            gene_names, gene_surfaces, statuses, status_surfaces, rng,
            compounding_tier=compounding_tier,
        )

    frame = rng.choice(_L4_FRAMES)
    inner: str = frame["inner"]
    sep: str = frame["sep"]
    final: Optional[str] = frame["final"]
    use_oxford: bool = frame["final_is_oxford"]

    pieces: list[str] = []
    assertions: list[AssertionFact] = []
    pos = 0
    last_idx = len(gene_names) - 1
    for i, (name, gene_surface, status, status_surface) in enumerate(
        zip(gene_names, gene_surfaces, statuses, status_surfaces)
    ):
        if i > 0:
            joiner = _l4_joiner(sep, final, use_oxford, i, last_idx)
            pieces.append(joiner)
            pos += len(joiner)
        gene_start = pos
        pieces.append(gene_surface)
        pos += len(gene_surface)
        gene_end = pos
        pieces.append(inner)
        pos += len(inner)
        status_start = pos
        pieces.append(status_surface)
        pos += len(status_surface)
        status_end = pos
        assertions.append(
            AssertionFact(
                gene=name,
                status=status,
                spans={
                    "gene": (gene_start, gene_end),
                    "status": (status_start, status_end),
                },
            )
        )

    pieces.append(".")
    sentence = "".join(pieces)
    frame_template = f"{{gene}}{inner}{{status}}{sep}"
    return RenderedRecord(
        sentence=sentence,
        assertions=tuple(assertions),
        frame_template=frame_template,
        complexity_level="L4",
        compounding_tier=compounding_tier,
    )


def _render_l4_shorthand(
    gene_names: list[str],
    gene_surfaces: list[str],
    statuses: list[str],
    status_surfaces: list[str],
    rng: random.Random,
    compounding_tier: str = "low",
) -> RenderedRecord:
    """Glue per-gene (gene, status) shorthand pairs into one surface.

    Parallels :func:`_render_l3_shorthand` but each gene has its own
    shorthand status surface (drawn per-gene from per-gene population).
    No terminal period — shorthand frames are tabular, not sentences.
    """
    frame = rng.choice(_L4_SHORTHAND_FRAMES)
    inner = frame["inner"]
    sep = frame["sep"]

    pieces: list[str] = []
    assertions: list[AssertionFact] = []
    pos = 0
    for i, (name, gene_surface, status, status_surface) in enumerate(
        zip(gene_names, gene_surfaces, statuses, status_surfaces)
    ):
        if i > 0:
            pieces.append(sep)
            pos += len(sep)
        gene_start = pos
        pieces.append(gene_surface)
        pos += len(gene_surface)
        gene_end = pos
        pieces.append(inner)
        pos += len(inner)
        status_start = pos
        pieces.append(status_surface)
        pos += len(status_surface)
        status_end = pos
        assertions.append(
            AssertionFact(
                gene=name,
                status=status,
                spans={
                    "gene": (gene_start, gene_end),
                    "status": (status_start, status_end),
                },
            )
        )

    sentence = "".join(pieces)
    frame_template = f"{{gene}}{inner}{{status}}{sep}"
    return RenderedRecord(
        sentence=sentence,
        assertions=tuple(assertions),
        frame_template=frame_template,
        complexity_level="L4S",
        compounding_tier=compounding_tier,
    )


def _l4_joiner(
    sep: str, final: Optional[str], use_oxford: bool, i: int, last_idx: int
) -> str:
    """Choose the glue that precedes pair ``i`` (1-indexed past the first).

    ``final`` inserts a coordinator before the last pair. For comma-separated
    frames, ``use_oxford=True`` keeps the pre-coordinator comma
    (``", and "``); otherwise the comma is dropped (``" and "``). Non-comma
    separators ignore the Oxford toggle.
    """
    is_last = i == last_idx
    if is_last and final is not None:
        if sep == ", ":
            return ", and " if use_oxford else " and "
        return f"{sep}{final} "
    return sep


def render_l5_record(
    biomarker_pool: list[str],
    profile: PatientProfile,
    biomarkers: BiomarkerConfig,
    common: CommonConfig,
    rng: random.Random,
    n_wide: Optional[int] = None,
    compounding_tier: str = "low",
) -> RenderedRecord:
    """Render one L5 negation-scope record.

    ``compounding_tier`` drives ``n_wide`` via :func:`n_range_for_tier`
    (capped at pool size). Panel-wide frames ignore the tier since they
    always yield one labeled fact. The caller is responsible for keeping
    the pool homogeneous (mutation vs expression) except in the high tier
    where cross-class mixing is explicitly allowed.

    Two structural frame kinds:

    - ``"enumerated"``: a scope marker ("No" / "Absence of") negates a list
      of N genes wide (``N_wide ∈ tier range``). All wide facts share one
      status span on the scope marker and carry
      ``polarity_scope="negation_wide"`` with status ``"negative"``. An
      optional exception clause ("and" continuation with same negative
      polarity, or "except" flip to positive) adds one trailing gene with
      ``polarity_scope="exception"`` and its own status span. The majority
      of frames carry no exception — this matches clinical register.

    - ``"panel_wide"``: the wide scope is an implicit panel ("No biomarker
      on the panel was positive"); no genes are enumerated, so no
      wide-scope facts are emitted. One trailing gene after "other than"
      becomes the sole fact with ``polarity_scope="exception"`` and status
      set to the frame's ``status_word``.
    """
    pool_size = len(biomarker_pool)
    if pool_size < L5_MIN_WIDE:
        raise RenderError(
            f"L5 pool needs >= {L5_MIN_WIDE} biomarkers, got {pool_size}"
        )

    frame = rng.choice(_L5_FRAMES)
    if frame["kind"] == "panel_wide":
        return _render_l5_panel_wide(
            frame,
            biomarker_pool,
            biomarkers,
            rng,
            compounding_tier=compounding_tier,
        )

    scope_marker: str = frame["scope_marker"]
    neg_noun: str = frame["neg_noun"]
    exception_marker: Optional[str] = frame["exception_marker"]
    exception_status: Optional[str] = frame["exception_status"]

    tier_min, tier_max = n_range_for_tier(compounding_tier)
    min_wide = min(tier_min, pool_size)
    max_wide = min(tier_max, pool_size)
    if n_wide is None:
        n_wide = rng.randint(min_wide, max_wide)
    elif not min_wide <= n_wide <= max_wide:
        raise RenderError(
            f"n_wide must be in [{min_wide}, {max_wide}] for tier "
            f"{compounding_tier!r} and pool size {pool_size}, got {n_wide}"
        )

    has_exception = exception_marker is not None and pool_size > n_wide
    total_n = n_wide + (1 if has_exception else 0)
    gene_names = rng.sample(biomarker_pool, total_n)
    wide_names = gene_names[:n_wide]
    ex_name = gene_names[n_wide] if has_exception else None

    wide_surfaces: list[str] = []
    for name in wide_names:
        b = biomarkers.get(name)
        form_id = _bare_gene_name_form(b, rng)
        wide_surfaces.append(b.name_forms.realizations[form_id])

    ex_surface: Optional[str] = None
    ex_status_surface: Optional[str] = None
    ex_status: Optional[str] = None
    if has_exception and ex_name is not None:
        b = biomarkers.get(ex_name)
        form_id = _bare_gene_name_form(b, rng)
        ex_surface = b.name_forms.realizations[form_id]
        ex_status = exception_status or "positive"
        _, raw = _sample_status_phrase(
            ex_status, common, rng, complexity_level="L1"
        )
        ex_status_surface = expand_placeholders(raw, {"gene": ex_surface})

    gene_list_surface, gene_positions_in_list = _coordinate_gene_list(
        wide_surfaces, rng
    )

    pieces: list[str] = []
    pos = 0
    scope_start = pos
    pieces.append(scope_marker)
    pos += len(scope_marker)
    scope_end = pos
    pieces.append(" ")
    pos += 1
    gene_list_start = pos
    pieces.append(gene_list_surface)
    pos += len(gene_list_surface)
    pieces.append(" ")
    pos += 1
    pieces.append(neg_noun)
    pos += len(neg_noun)

    ex_gene_span: Optional[tuple[int, int]] = None
    ex_status_span: Optional[tuple[int, int]] = None
    if has_exception and ex_surface is not None and ex_status_surface is not None:
        pieces.append(" ")
        pos += 1
        pieces.append(exception_marker)  # type: ignore[arg-type]
        pos += len(exception_marker)  # type: ignore[arg-type]
        pieces.append(" ")
        pos += 1
        ex_gene_start = pos
        pieces.append(ex_surface)
        pos += len(ex_surface)
        ex_gene_end = pos
        ex_gene_span = (ex_gene_start, ex_gene_end)
        pieces.append(" ")
        pos += 1
        ex_status_start = pos
        pieces.append(ex_status_surface)
        pos += len(ex_status_surface)
        ex_status_end = pos
        ex_status_span = (ex_status_start, ex_status_end)

    pieces.append(".")
    sentence = "".join(pieces)

    shared_status_span = (scope_start, scope_end)
    assertions: list[AssertionFact] = []
    for name, (s, e) in zip(wide_names, gene_positions_in_list):
        assertions.append(
            AssertionFact(
                gene=name,
                status="negative",
                spans={
                    "gene": (gene_list_start + s, gene_list_start + e),
                    "status": shared_status_span,
                },
                polarity_scope="negation_wide",
            )
        )
    if has_exception and ex_name is not None and ex_gene_span is not None and ex_status_span is not None:
        assertions.append(
            AssertionFact(
                gene=ex_name,
                status=ex_status or "positive",
                spans={
                    "gene": ex_gene_span,
                    "status": ex_status_span,
                },
                polarity_scope="exception",
            )
        )

    if has_exception:
        frame_template = (
            f"{scope_marker} {{gene_list}} {neg_noun} "
            f"{exception_marker} {{ex_gene}} {{ex_status}}."
        )
    else:
        frame_template = f"{scope_marker} {{gene_list}} {neg_noun}."
    return RenderedRecord(
        sentence=sentence,
        assertions=tuple(assertions),
        frame_template=frame_template,
        complexity_level="L5",
        compounding_tier=compounding_tier,
    )


def _render_l5_panel_wide(
    frame: dict,
    biomarker_pool: list[str],
    biomarkers: BiomarkerConfig,
    rng: random.Random,
    compounding_tier: str = "low",
) -> RenderedRecord:
    """Render a panel-wide L5 record: implicit panel scope, one exception.

    Surface: ``"{scope} {status_word} other than {gene}."`` — the scope
    phrase embeds the polarity ("positive") as a literal token that serves
    as the status span for the single labeled fact. No wide-scope facts
    are emitted because the panel genes are not named.
    """
    scope: str = frame["scope"]
    status_word: str = frame["status_word"]
    exception_marker: str = frame["exception_marker"]

    ex_name = rng.choice(biomarker_pool)
    b = biomarkers.get(ex_name)
    form_id = _bare_gene_name_form(b, rng)
    ex_surface = b.name_forms.realizations[form_id]

    pieces: list[str] = []
    pos = 0
    pieces.append(scope)
    pos += len(scope)
    pieces.append(" ")
    pos += 1
    status_start = pos
    pieces.append(status_word)
    pos += len(status_word)
    status_end = pos
    pieces.append(" ")
    pos += 1
    pieces.append(exception_marker)
    pos += len(exception_marker)
    pieces.append(" ")
    pos += 1
    gene_start = pos
    pieces.append(ex_surface)
    pos += len(ex_surface)
    gene_end = pos
    pieces.append(".")
    sentence = "".join(pieces)

    fact = AssertionFact(
        gene=ex_name,
        status=status_word,
        spans={
            "gene": (gene_start, gene_end),
            "status": (status_start, status_end),
        },
        polarity_scope="exception",
    )
    frame_template = f"{scope} {status_word} {exception_marker} {{ex_gene}}."
    return RenderedRecord(
        sentence=sentence,
        assertions=(fact,),
        frame_template=frame_template,
        complexity_level="L5",
        compounding_tier=compounding_tier,
    )


# --- L6 temporal / certainty qualification (Sub-phase 3.8) ---------------

L6_TEMPORAL_VOCAB: tuple[str, ...] = (
    "previously",
    "currently",
    "at diagnosis",
    "at relapse",
    "post-TKI",
)

L6_CERTAINTY_VOCAB: tuple[str, ...] = (
    "suspected",
    "confirmed",
    "probable",
    "rule-out",
    "pending",
)

L6_SHAPES: tuple[str, ...] = ("temporal", "certainty", "combined")

# Status realizations that cannot be used in L6 contexts:
# - "confirmed" overlaps with the certainty vocab; rendering it as a status
#   beside a certainty modifier produces ambiguous spans.
# - "no evidence of" is an incomplete phrase requiring a trailing noun and
#   doesn't fit the L6 frame shapes that place status mid-clause.
_L6_STATUS_EXCLUDE: frozenset[str] = frozenset({
    "confirmed",
    "no evidence of",
})


def _sample_l6_status_surface(
    status: str, common: CommonConfig, rng: random.Random
) -> str:
    """Draw a status surface suitable for L6 frames.

    Uses the formal phrase category (``positive_phrases`` /
    ``negation_phrases``) with two extra filters on top of
    :func:`_filter_gene_free`:

    - excludes realizations listed in :data:`_L6_STATUS_EXCLUDE` (which
      collide with certainty vocab or require trailing nouns)
    - renormalizes the remaining distribution before sampling
    """
    category_name = _STATUS_CATEGORY[status]
    category = common.categories.get(category_name)
    if not isinstance(category, WeightedVariations):
        raise RenderError(
            f"common.{category_name} is missing or wrong schema_type"
        )
    eligible: dict[str, float] = {}
    for vid, weight in category.variations.items():
        surface = category.realizations[vid]
        if "{gene}" in surface:
            continue
        if surface in _L6_STATUS_EXCLUDE:
            continue
        eligible[vid] = weight
    if not eligible:
        raise RenderError(
            f"no L6-eligible realizations for status {status!r}"
        )
    total = sum(eligible.values())
    norm = {k: v / total for k, v in eligible.items()}
    keys = list(norm)
    weights = [norm[k] for k in keys]
    vid = rng.choices(keys, weights=weights, k=1)[0]
    return category.realizations[vid]


def _capitalize_first(s: str) -> str:
    return s[0].upper() + s[1:] if s else s


def _pick_contrast_statuses(rng: random.Random) -> tuple[str, str]:
    """Pick a positive/negative pair for temporal contrast (order randomized)."""
    return rng.choice([("positive", "negative"), ("negative", "positive")])


def render_l6_record(
    biomarker_pool: list[str],
    profile: PatientProfile,
    biomarkers: BiomarkerConfig,
    common: CommonConfig,
    rng: random.Random,
    shape: Optional[str] = None,
) -> RenderedRecord:
    """Render one L6 qualifier record.

    Single-gene record in one of three structural shapes:

    - ``"temporal"``: 2 facts, same gene, contrasting statuses at two
      timepoints. Each fact carries a ``temporal`` marker; ``certainty``
      is None. Captures clinical change over time.
    - ``"certainty"``: 1 fact with a hedging certainty modifier; no
      temporal. Captures diagnostic confidence.
    - ``"combined"``: 2 facts like temporal contrast, each also carrying
      an independent certainty value.

    ``shape=None`` selects uniformly from :data:`L6_SHAPES`. Records stamp
    ``complexity_level="L6"`` and ``compounding_tier="low"`` (tier knob
    doesn't apply at single-gene level).

    ``profile`` is accepted for signature parity with other renderers but
    unused — L6 contrast statuses are drawn directly from a fixed
    positive/negative pair to guarantee divergence.
    """
    del profile  # reserved for parity; see docstring
    if shape is None:
        shape = rng.choice(L6_SHAPES)
    if shape not in L6_SHAPES:
        raise RenderError(
            f"L6 shape must be one of {L6_SHAPES}, got {shape!r}"
        )
    if shape == "temporal":
        return _render_l6_temporal(biomarker_pool, biomarkers, common, rng)
    if shape == "certainty":
        return _render_l6_certainty(biomarker_pool, biomarkers, common, rng)
    return _render_l6_combined(biomarker_pool, biomarkers, common, rng)


def _render_l6_temporal(
    biomarker_pool: list[str],
    biomarkers: BiomarkerConfig,
    common: CommonConfig,
    rng: random.Random,
) -> RenderedRecord:
    gene = rng.choice(biomarker_pool)
    b = biomarkers.get(gene)
    form_id = _bare_gene_name_form(b, rng)
    gene_surface = b.name_forms.realizations[form_id]

    status_a, status_b = _pick_contrast_statuses(rng)
    surf_a = _sample_l6_status_surface(status_a, common, rng)
    surf_b = _sample_l6_status_surface(status_b, common, rng)
    temporal_a, temporal_b = rng.sample(L6_TEMPORAL_VOCAB, 2)

    frame_kind = rng.randrange(3)
    if frame_kind == 0:
        sentence, gspans, sspans, frame_template = _l6_temporal_frame1(
            gene_surface, surf_a, surf_b, temporal_a, temporal_b
        )
    elif frame_kind == 1:
        sentence, gspans, sspans, frame_template = _l6_temporal_frame2(
            gene_surface, surf_a, surf_b, temporal_a, temporal_b
        )
    else:
        sentence, gspans, sspans, frame_template = _l6_temporal_frame3(
            gene_surface, surf_a, surf_b, temporal_a, temporal_b
        )

    fact_a = AssertionFact(
        gene=gene,
        status=status_a,
        spans={"gene": gspans[0], "status": sspans[0]},
        temporal=temporal_a,
    )
    fact_b = AssertionFact(
        gene=gene,
        status=status_b,
        spans={"gene": gspans[1], "status": sspans[1]},
        temporal=temporal_b,
    )
    return RenderedRecord(
        sentence=sentence,
        assertions=(fact_a, fact_b),
        frame_template=frame_template,
        complexity_level="L6",
        compounding_tier="low",
    )


def _l6_temporal_frame1(
    gene_surface: str,
    status_a: str,
    status_b: str,
    temporal_a: str,
    temporal_b: str,
) -> tuple[str, tuple[tuple[int, int], tuple[int, int]], tuple[tuple[int, int], tuple[int, int]], str]:
    """``{gene}: {status_a} {temporal_a}, {status_b} {temporal_b}.``

    Both facts share the gene span.
    """
    pieces: list[str] = []
    pos = 0
    g_start = pos
    pieces.append(gene_surface); pos += len(gene_surface)
    g_end = pos
    pieces.append(": "); pos += 2
    sa_start = pos
    pieces.append(status_a); pos += len(status_a)
    sa_end = pos
    pieces.append(" "); pos += 1
    pieces.append(temporal_a); pos += len(temporal_a)
    pieces.append(", "); pos += 2
    sb_start = pos
    pieces.append(status_b); pos += len(status_b)
    sb_end = pos
    pieces.append(" "); pos += 1
    pieces.append(temporal_b); pos += len(temporal_b)
    pieces.append(".")
    sentence = "".join(pieces)
    gene_span = (g_start, g_end)
    return (
        sentence,
        (gene_span, gene_span),
        ((sa_start, sa_end), (sb_start, sb_end)),
        "{gene}: {status_a} {temporal_a}, {status_b} {temporal_b}.",
    )


def _l6_temporal_frame2(
    gene_surface: str,
    status_a: str,
    status_b: str,
    temporal_a: str,
    temporal_b: str,
) -> tuple[str, tuple[tuple[int, int], tuple[int, int]], tuple[tuple[int, int], tuple[int, int]], str]:
    """``{gene} was {status_a} {temporal_a}; {gene} was {status_b} {temporal_b}.``

    Two gene mentions — each fact owns its own gene span. Gene leads to
    avoid sentence-start capitalization of the temporal marker (which
    would break the temporal-label-in-sentence invariant).
    """
    pieces: list[str] = []
    pos = 0
    g1_start = pos
    pieces.append(gene_surface); pos += len(gene_surface)
    g1_end = pos
    pieces.append(" was "); pos += 5
    sa_start = pos
    pieces.append(status_a); pos += len(status_a)
    sa_end = pos
    pieces.append(" "); pos += 1
    pieces.append(temporal_a); pos += len(temporal_a)
    pieces.append("; "); pos += 2
    g2_start = pos
    pieces.append(gene_surface); pos += len(gene_surface)
    g2_end = pos
    pieces.append(" was "); pos += 5
    sb_start = pos
    pieces.append(status_b); pos += len(status_b)
    sb_end = pos
    pieces.append(" "); pos += 1
    pieces.append(temporal_b); pos += len(temporal_b)
    pieces.append(".")
    sentence = "".join(pieces)
    return (
        sentence,
        ((g1_start, g1_end), (g2_start, g2_end)),
        ((sa_start, sa_end), (sb_start, sb_end)),
        "{gene} was {status_a} {temporal_a}; {gene} was {status_b} {temporal_b}.",
    )


def _l6_temporal_frame3(
    gene_surface: str,
    status_a: str,
    status_b: str,
    temporal_a: str,
    temporal_b: str,
) -> tuple[str, tuple[tuple[int, int], tuple[int, int]], tuple[tuple[int, int], tuple[int, int]], str]:
    """``{Status_a} {temporal_a}, {status_b} {temporal_b} for {gene}.``

    Trailing gene; both facts share the gene span.
    """
    status_a_cap = _capitalize_first(status_a)
    pieces: list[str] = []
    pos = 0
    sa_start = pos
    pieces.append(status_a_cap); pos += len(status_a_cap)
    sa_end = pos
    pieces.append(" "); pos += 1
    pieces.append(temporal_a); pos += len(temporal_a)
    pieces.append(", "); pos += 2
    sb_start = pos
    pieces.append(status_b); pos += len(status_b)
    sb_end = pos
    pieces.append(" "); pos += 1
    pieces.append(temporal_b); pos += len(temporal_b)
    pieces.append(" for "); pos += 5
    g_start = pos
    pieces.append(gene_surface); pos += len(gene_surface)
    g_end = pos
    pieces.append(".")
    sentence = "".join(pieces)
    gene_span = (g_start, g_end)
    return (
        sentence,
        (gene_span, gene_span),
        ((sa_start, sa_end), (sb_start, sb_end)),
        "{Status_a} {temporal_a}, {status_b} {temporal_b} for {gene}.",
    )


def _render_l6_certainty(
    biomarker_pool: list[str],
    biomarkers: BiomarkerConfig,
    common: CommonConfig,
    rng: random.Random,
) -> RenderedRecord:
    gene = rng.choice(biomarker_pool)
    b = biomarkers.get(gene)
    form_id = _bare_gene_name_form(b, rng)
    gene_surface = b.name_forms.realizations[form_id]

    status = rng.choice(("positive", "negative"))
    status_surface = _sample_l6_status_surface(status, common, rng)
    certainty = rng.choice(L6_CERTAINTY_VOCAB)

    frame_kind = rng.randrange(3)
    if frame_kind == 0:
        sentence, gene_span, status_span, frame_template = _l6_certainty_frame1(
            gene_surface, status_surface, certainty
        )
    elif frame_kind == 1:
        sentence, gene_span, status_span, frame_template = _l6_certainty_frame2(
            gene_surface, status_surface, certainty
        )
    else:
        sentence, gene_span, status_span, frame_template = _l6_certainty_frame3(
            gene_surface, status_surface, certainty
        )

    fact = AssertionFact(
        gene=gene,
        status=status,
        spans={"gene": gene_span, "status": status_span},
        certainty=certainty,
    )
    return RenderedRecord(
        sentence=sentence,
        assertions=(fact,),
        frame_template=frame_template,
        complexity_level="L6",
        compounding_tier="low",
    )


def _l6_certainty_frame1(
    gene_surface: str, status_surface: str, certainty: str
) -> tuple[str, tuple[int, int], tuple[int, int], str]:
    """``{gene} {status} ({certainty}).``"""
    pieces: list[str] = []
    pos = 0
    g_start = pos
    pieces.append(gene_surface); pos += len(gene_surface)
    g_end = pos
    pieces.append(" "); pos += 1
    s_start = pos
    pieces.append(status_surface); pos += len(status_surface)
    s_end = pos
    pieces.append(" ("); pos += 2
    pieces.append(certainty); pos += len(certainty)
    pieces.append(").")
    sentence = "".join(pieces)
    return sentence, (g_start, g_end), (s_start, s_end), "{gene} {status} ({certainty})."


def _l6_certainty_frame2(
    gene_surface: str, status_surface: str, certainty: str
) -> tuple[str, tuple[int, int], tuple[int, int], str]:
    """``{gene}: {status} ({certainty}).``"""
    pieces: list[str] = []
    pos = 0
    g_start = pos
    pieces.append(gene_surface); pos += len(gene_surface)
    g_end = pos
    pieces.append(": "); pos += 2
    s_start = pos
    pieces.append(status_surface); pos += len(status_surface)
    s_end = pos
    pieces.append(" ("); pos += 2
    pieces.append(certainty); pos += len(certainty)
    pieces.append(").")
    sentence = "".join(pieces)
    return sentence, (g_start, g_end), (s_start, s_end), "{gene}: {status} ({certainty})."


def _l6_certainty_frame3(
    gene_surface: str, status_surface: str, certainty: str
) -> tuple[str, tuple[int, int], tuple[int, int], str]:
    """``{gene} {status}, {certainty}.``"""
    pieces: list[str] = []
    pos = 0
    g_start = pos
    pieces.append(gene_surface); pos += len(gene_surface)
    g_end = pos
    pieces.append(" "); pos += 1
    s_start = pos
    pieces.append(status_surface); pos += len(status_surface)
    s_end = pos
    pieces.append(", "); pos += 2
    pieces.append(certainty); pos += len(certainty)
    pieces.append(".")
    sentence = "".join(pieces)
    return sentence, (g_start, g_end), (s_start, s_end), "{gene} {status}, {certainty}."


def _render_l6_combined(
    biomarker_pool: list[str],
    biomarkers: BiomarkerConfig,
    common: CommonConfig,
    rng: random.Random,
) -> RenderedRecord:
    gene = rng.choice(biomarker_pool)
    b = biomarkers.get(gene)
    form_id = _bare_gene_name_form(b, rng)
    gene_surface = b.name_forms.realizations[form_id]

    status_a, status_b = _pick_contrast_statuses(rng)
    surf_a = _sample_l6_status_surface(status_a, common, rng)
    surf_b = _sample_l6_status_surface(status_b, common, rng)
    temporal_a, temporal_b = rng.sample(L6_TEMPORAL_VOCAB, 2)
    certainty_a = rng.choice(L6_CERTAINTY_VOCAB)
    certainty_b = rng.choice(L6_CERTAINTY_VOCAB)

    frame_kind = rng.randrange(2)
    if frame_kind == 0:
        sentence, gspans, sspans, frame_template = _l6_combined_frame1(
            gene_surface, surf_a, surf_b, temporal_a, temporal_b, certainty_a, certainty_b
        )
    else:
        sentence, gspans, sspans, frame_template = _l6_combined_frame2(
            gene_surface, surf_a, surf_b, temporal_a, temporal_b, certainty_a, certainty_b
        )

    fact_a = AssertionFact(
        gene=gene,
        status=status_a,
        spans={"gene": gspans[0], "status": sspans[0]},
        temporal=temporal_a,
        certainty=certainty_a,
    )
    fact_b = AssertionFact(
        gene=gene,
        status=status_b,
        spans={"gene": gspans[1], "status": sspans[1]},
        temporal=temporal_b,
        certainty=certainty_b,
    )
    return RenderedRecord(
        sentence=sentence,
        assertions=(fact_a, fact_b),
        frame_template=frame_template,
        complexity_level="L6",
        compounding_tier="low",
    )


def _l6_combined_frame1(
    gene_surface: str,
    status_a: str,
    status_b: str,
    temporal_a: str,
    temporal_b: str,
    certainty_a: str,
    certainty_b: str,
) -> tuple[str, tuple[tuple[int, int], tuple[int, int]], tuple[tuple[int, int], tuple[int, int]], str]:
    """``{gene}: {status_a} {temporal_a} ({certainty_a}), {status_b} {temporal_b} ({certainty_b}).``"""
    pieces: list[str] = []
    pos = 0
    g_start = pos
    pieces.append(gene_surface); pos += len(gene_surface)
    g_end = pos
    pieces.append(": "); pos += 2
    sa_start = pos
    pieces.append(status_a); pos += len(status_a)
    sa_end = pos
    pieces.append(" "); pos += 1
    pieces.append(temporal_a); pos += len(temporal_a)
    pieces.append(" ("); pos += 2
    pieces.append(certainty_a); pos += len(certainty_a)
    pieces.append("), "); pos += 3
    sb_start = pos
    pieces.append(status_b); pos += len(status_b)
    sb_end = pos
    pieces.append(" "); pos += 1
    pieces.append(temporal_b); pos += len(temporal_b)
    pieces.append(" ("); pos += 2
    pieces.append(certainty_b); pos += len(certainty_b)
    pieces.append(").")
    sentence = "".join(pieces)
    gene_span = (g_start, g_end)
    return (
        sentence,
        (gene_span, gene_span),
        ((sa_start, sa_end), (sb_start, sb_end)),
        "{gene}: {status_a} {temporal_a} ({certainty_a}), {status_b} {temporal_b} ({certainty_b}).",
    )


def _l6_combined_frame2(
    gene_surface: str,
    status_a: str,
    status_b: str,
    temporal_a: str,
    temporal_b: str,
    certainty_a: str,
    certainty_b: str,
) -> tuple[str, tuple[tuple[int, int], tuple[int, int]], tuple[tuple[int, int], tuple[int, int]], str]:
    """``{gene} was {status_a} {temporal_a} ({certainty_a}); {gene} was {status_b} {temporal_b} ({certainty_b}).``

    Gene leads each clause to avoid sentence-start capitalization of a
    labeled temporal marker.
    """
    pieces: list[str] = []
    pos = 0
    g1_start = pos
    pieces.append(gene_surface); pos += len(gene_surface)
    g1_end = pos
    pieces.append(" was "); pos += 5
    sa_start = pos
    pieces.append(status_a); pos += len(status_a)
    sa_end = pos
    pieces.append(" "); pos += 1
    pieces.append(temporal_a); pos += len(temporal_a)
    pieces.append(" ("); pos += 2
    pieces.append(certainty_a); pos += len(certainty_a)
    pieces.append("); "); pos += 3
    g2_start = pos
    pieces.append(gene_surface); pos += len(gene_surface)
    g2_end = pos
    pieces.append(" was "); pos += 5
    sb_start = pos
    pieces.append(status_b); pos += len(status_b)
    sb_end = pos
    pieces.append(" "); pos += 1
    pieces.append(temporal_b); pos += len(temporal_b)
    pieces.append(" ("); pos += 2
    pieces.append(certainty_b); pos += len(certainty_b)
    pieces.append(").")
    sentence = "".join(pieces)
    return (
        sentence,
        ((g1_start, g1_end), (g2_start, g2_end)),
        ((sa_start, sa_end), (sb_start, sb_end)),
        "{gene} was {status_a} {temporal_a} ({certainty_a}); {gene} was {status_b} {temporal_b} ({certainty_b}).",
    )


# === L7 cross-sentence coreference (Sub-phase 3.9) =====================

L7_ANAPHOR_NOUNS: tuple[str, ...] = (
    "mutation",
    "variant",
    "alteration",
    "finding",
    "result",
)

L7_SHAPES: tuple[str, ...] = (
    "setup_claim",
    "claim_anaphora",
    "setup_claim_qualifier",
)

_L7_SETUP_TEMPLATES: tuple[str, ...] = (
    "Molecular profiling was performed on the specimen.",
    "Comprehensive NGS testing was completed for this case.",
    "The genomic panel report summarizes the following.",
    "Targeted sequencing analysis has been performed.",
    "An orthogonal validation was run on the specimen.",
)

_L7_ANAPHOR_TEMPLATES: tuple[str, ...] = (
    "The {anaphor} is associated with targeted therapy response.",
    "This {anaphor} has known clinical implications.",
    "The {anaphor} was confirmed by orthogonal methods.",
    "This {anaphor} was reviewed by the molecular tumor board.",
    "The {anaphor} appears in published treatment guidelines.",
)

_L7_QUALIFIER_TEMPLATES: tuple[str, ...] = (
    "This {anaphor} was classified as pathogenic.",
    "The {anaphor} was reported as clinically actionable.",
    "The {anaphor} warrants further clinical correlation.",
    "This {anaphor} is consistent with the reported histology.",
    "The {anaphor} prompted a tumor board review.",
)

_L7_CLAIM_TEMPLATES: tuple[str, ...] = (
    "{gene} was {status}.",
    "{gene}: {status}.",
)

_L7_STATUS_EXCLUDE: frozenset[str] = frozenset({"no evidence of"})


def _sample_l7_status_surface(
    status: str, common: CommonConfig, rng: random.Random
) -> str:
    """Draw a status surface suitable for L7 claim frames.

    Filters gene-dependent realizations (via ``{gene}`` marker) and the
    incomplete ``"no evidence of"`` phrase which requires a trailing noun.
    Remaining distribution is renormalized before sampling.
    """
    category_name = _STATUS_CATEGORY[status]
    category = common.categories.get(category_name)
    if not isinstance(category, WeightedVariations):
        raise RenderError(
            f"common.{category_name} is missing or wrong schema_type"
        )
    eligible: dict[str, float] = {}
    for vid, weight in category.variations.items():
        surface = category.realizations[vid]
        if "{gene}" in surface:
            continue
        if surface in _L7_STATUS_EXCLUDE:
            continue
        eligible[vid] = weight
    if not eligible:
        raise RenderError(
            f"no L7-eligible realizations for status {status!r}"
        )
    total = sum(eligible.values())
    norm = {k: v / total for k, v in eligible.items()}
    keys = list(norm)
    weights = [norm[k] for k in keys]
    vid = rng.choices(keys, weights=weights, k=1)[0]
    return category.realizations[vid]


def _build_l7_claim_sentence(
    gene_name: str, status_surface: str, rng: random.Random
) -> tuple[str, tuple[int, int], tuple[int, int]]:
    """Return (sentence, gene_span_local, status_span_local) for the claim."""
    template = rng.choice(_L7_CLAIM_TEMPLATES)
    if template == "{gene} was {status}.":
        mid = " was "
    elif template == "{gene}: {status}.":
        mid = ": "
    else:
        raise RenderError(f"unknown L7 claim template: {template!r}")
    g_start = 0
    g_end = len(gene_name)
    s_start = g_end + len(mid)
    s_end = s_start + len(status_surface)
    sentence = f"{gene_name}{mid}{status_surface}."
    return sentence, (g_start, g_end), (s_start, s_end)


def _assemble_l7_record(
    gene: str,
    status: str,
    sentences: list[str],
    claim_idx: int,
    claim_local_gene: tuple[int, int],
    claim_local_status: tuple[int, int],
    frame_template: str,
) -> RenderedRecord:
    joined = " ".join(sentences)
    offset = sum(len(sentences[i]) + 1 for i in range(claim_idx))
    g_global = (claim_local_gene[0] + offset, claim_local_gene[1] + offset)
    s_global = (claim_local_status[0] + offset, claim_local_status[1] + offset)
    fact = AssertionFact(
        gene=gene,
        status=status,
        spans={"gene": g_global, "status": s_global},
        sentence_index=claim_idx,
    )
    return RenderedRecord(
        sentence=joined,
        sentences=tuple(sentences),
        assertions=(fact,),
        frame_template=frame_template,
        complexity_level="L7",
        compounding_tier="low",
    )


def render_l7_record(
    biomarker_pool: list[str],
    profile: PatientProfile,
    biomarkers: BiomarkerConfig,
    common: CommonConfig,
    rng: random.Random,
    shape: Optional[str] = None,
) -> RenderedRecord:
    """Render an L7 multi-sentence coreference record.

    Three shapes:

    - ``"setup_claim"`` — generic setup sentence followed by a gene/status
      claim. Fact ``sentence_index=1``.
    - ``"claim_anaphora"`` — claim followed by a nominal anaphor sentence
      ("the mutation is associated with..."). Fact ``sentence_index=0``.
    - ``"setup_claim_qualifier"`` — setup, claim, qualifier (3 sentences).
      Fact ``sentence_index=1``.
    """
    if shape is None:
        shape = rng.choice(L7_SHAPES)
    if shape not in L7_SHAPES:
        raise RenderError(f"unknown L7 shape: {shape!r}")

    gene = rng.choice(biomarker_pool)
    status = rng.choice(("positive", "negative"))
    gene_cfg = biomarkers.get(gene)
    gene_name_forms = list(gene_cfg.name_forms.realizations.values())
    gene_name = rng.choice(gene_name_forms)
    status_surface = _sample_l7_status_surface(status, common, rng)

    claim, g_span, s_span = _build_l7_claim_sentence(
        gene_name, status_surface, rng
    )

    if shape == "setup_claim":
        setup = rng.choice(_L7_SETUP_TEMPLATES)
        return _assemble_l7_record(
            gene, status, [setup, claim], 1, g_span, s_span,
            frame_template="setup. claim.",
        )
    if shape == "claim_anaphora":
        anaphor_noun = rng.choice(L7_ANAPHOR_NOUNS)
        anaphor = rng.choice(_L7_ANAPHOR_TEMPLATES).format(anaphor=anaphor_noun)
        return _assemble_l7_record(
            gene, status, [claim, anaphor], 0, g_span, s_span,
            frame_template="claim. anaphor.",
        )
    # setup_claim_qualifier
    setup = rng.choice(_L7_SETUP_TEMPLATES)
    anaphor_noun = rng.choice(L7_ANAPHOR_NOUNS)
    qualifier = rng.choice(_L7_QUALIFIER_TEMPLATES).format(anaphor=anaphor_noun)
    return _assemble_l7_record(
        gene, status, [setup, claim, qualifier], 1, g_span, s_span,
        frame_template="setup. claim. qualifier.",
    )
