"""Sentence renderer with placeholder substitution and span tracking.

Consumes the probability-weighted config and emits (sentence, spans,
metadata) tuples with char-level labels that satisfy Non-Negotiable #1: the
gene, variant, method, status, and clone spans returned alongside each
sentence are literal substrings at the recorded positions. Placeholders
(``{gene}``, ``{method}``, ``{result}``, ``{value}``) are expanded inside
realization strings before slot substitution per CONFIG_ARCHITECTURE.md §7.5.

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


@dataclass(frozen=True)
class RenderedRecord:
    """One rendered sentence + the assertion facts it carries.

    L1/L2 records have exactly one fact; L3+ will carry a tuple of facts. The
    tuple ordering follows the order in which genes appear in the surface
    sentence. ``complexity_level`` is ``"L1"`` / ``"L2"`` (Phase 2b); L3+
    values land in later sub-phases.
    """

    sentence: str
    assertions: tuple[AssertionFact, ...]
    frame_template: str
    complexity_level: str = "L1"


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

    Flow mirrors CONFIG_ARCHITECTURE.md §3. Supports mutation, fusion,
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
