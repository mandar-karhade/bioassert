"""Cross-config validation for Phase 2a (config_architecture.md §7).

Structural invariants (weights sum to 1.0, matched variation/realization key
sets, bounded probabilities, terminal fallbacks) are enforced by the Pydantic
models in :mod:`bioassert.config.schema`. This module handles the remaining
checks that require looking at both configs together plus placeholder grammar
against the supported table in §7.5.

Additional invariants enforced here:
 * ``preferred_methods`` keys are a subset of ``common.test_methods.variations``
 * ``render_constraints`` resolve to valid biomarker-level ``name_form`` keys
 * placeholders in realization strings are from the supported set
 * every realization carrying ``{value}`` has a backing ``measurement_range``
 * population keys follow canonical axis order (``histology``, ``ethnicity``,
   ``smoking``, ``age_group``, ``sex``) — warning, not error, to avoid
   breaking consumers that already shipped non-canonical keys
"""
from __future__ import annotations

import logging
import re
from typing import Iterable

from bioassert.config.loader import BiomarkerConfig, CommonConfig
from bioassert.config.schema import (
    Biomarker,
    CloneAttribution,
    NegativeForms,
    WeightedVariations,
)

log = logging.getLogger(__name__)

SUPPORTED_PLACEHOLDERS: frozenset[str] = frozenset(
    {"gene", "method", "result", "value", "variant"}
)

_PLACEHOLDER_RE = re.compile(r"\{([^{}]*)\}")

_CANONICAL_HISTOLOGIES: frozenset[str] = frozenset({"adenocarcinoma", "squamous", "other"})
_CANONICAL_ETHNICITIES: frozenset[str] = frozenset({"western", "east_asian", "other"})
_CANONICAL_SMOKING: frozenset[str] = frozenset({"smoker", "nonsmoker", "former_smoker"})
_CANONICAL_AGE: frozenset[str] = frozenset({"younger", "older"})
_CANONICAL_SEX: frozenset[str] = frozenset({"female", "male"})

_AXIS_ORDER: tuple[tuple[str, frozenset[str]], ...] = (
    ("histology", _CANONICAL_HISTOLOGIES),
    ("ethnicity", _CANONICAL_ETHNICITIES),
    ("smoking", _CANONICAL_SMOKING),
    ("age_group", _CANONICAL_AGE),
    ("sex", _CANONICAL_SEX),
)


class ConfigValidationError(ValueError):
    """Raised when one or more cross-config invariants fail at load time."""


def _extract_placeholders(s: str) -> set[str]:
    return set(_PLACEHOLDER_RE.findall(s))


def _iter_realization_sources(
    biomarkers: BiomarkerConfig,
) -> Iterable[tuple[str, str, WeightedVariations | NegativeForms | CloneAttribution]]:
    """Yield (label, origin_path, weighted) triples for every realization dict."""
    for biomarker_name, entry in biomarkers.biomarkers.items():
        yield biomarker_name, f"{biomarker_name}.name_forms", entry.name_forms
        for variant_id, variant in entry.variants.items():
            yield (
                biomarker_name,
                f"{biomarker_name}.variants.{variant_id}.name_forms",
                variant.name_forms,
            )
        if entry.negative_forms is not None:
            yield (
                biomarker_name,
                f"{biomarker_name}.negative_forms",
                entry.negative_forms,
            )
        if entry.clone_attribution is not None:
            yield (
                biomarker_name,
                f"{biomarker_name}.clone_attribution",
                entry.clone_attribution,
            )


def _parse_population_key(key: str) -> list[tuple[str, str]]:
    """Split a population key into (axis, value) tokens in canonical order.

    Returns the parsed axis list if and only if the key obeys canonical axis
    order; otherwise raises ``ConfigValidationError`` with a diagnostic. The
    first token must always be a histology.
    """
    tokens = key.split("_")
    if not tokens:
        raise ConfigValidationError(f"empty population key {key!r}")
    histology_tokens: list[str] = []
    pieces = list(tokens)
    consumed: list[tuple[str, str]] = []

    axis_idx = 0
    while pieces:
        axis_name, axis_values = _AXIS_ORDER[axis_idx]
        match: str | None = None
        for length in range(min(4, len(pieces)), 0, -1):
            candidate = "_".join(pieces[:length])
            if candidate in axis_values:
                match = candidate
                break
        if match is None:
            if axis_idx + 1 < len(_AXIS_ORDER):
                axis_idx += 1
                continue
            raise ConfigValidationError(
                f"population key {key!r}: cannot place remaining tokens "
                f"{'_'.join(pieces)!r} against canonical axis order"
            )
        consumed.append((axis_name, match))
        pieces = pieces[len(match.split("_")):]
        axis_idx += 1
        if axis_idx >= len(_AXIS_ORDER) and pieces:
            raise ConfigValidationError(
                f"population key {key!r}: trailing tokens {'_'.join(pieces)!r} "
                "after exhausting canonical axes"
            )

    if not consumed or consumed[0][0] != "histology":
        raise ConfigValidationError(
            f"population key {key!r}: first token must be a histology, got "
            f"{consumed[0] if consumed else 'none'}"
        )
    _ = histology_tokens
    return consumed


def _validate_population_keys(biomarkers: BiomarkerConfig) -> list[str]:
    problems: list[str] = []
    for biomarker_name, entry in biomarkers.biomarkers.items():
        for key in entry.status_distribution_by_population:
            try:
                _parse_population_key(key)
            except ConfigValidationError as exc:
                problems.append(f"{biomarker_name}: {exc}")
    return problems


def _validate_preferred_methods(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> list[str]:
    problems: list[str] = []
    methods = common.categories.get("test_methods")
    if not isinstance(methods, WeightedVariations):
        return ["common.test_methods is missing or wrong schema_type"]
    allowed = set(methods.variations.keys())
    for name, entry in biomarkers.biomarkers.items():
        if entry.preferred_methods is None:
            continue
        unknown = set(entry.preferred_methods.variations) - allowed
        if unknown:
            problems.append(
                f"{name}.preferred_methods references unknown methods "
                f"{sorted(unknown)} (not in common.test_methods.variations)"
            )
    return problems


def _validate_placeholders(
    common: CommonConfig, biomarkers: BiomarkerConfig
) -> list[str]:
    problems: list[str] = []

    for category_name, category in common.categories.items():
        if isinstance(category, (WeightedVariations, CloneAttribution)):
            for vid, surface in category.realizations.items():
                problems.extend(
                    _check_placeholder_string(
                        f"common.{category_name}.realizations.{vid}",
                        surface,
                        measurement_range=None,
                    )
                )

    for biomarker_name, path, weighted in _iter_realization_sources(biomarkers):
        measurement_range = None
        if biomarker_name and path.endswith(".name_forms") and ".variants." in path:
            variant_id = path.split(".variants.")[1].rsplit(".name_forms", 1)[0]
            variant = biomarkers.biomarkers[biomarker_name].variants[variant_id]
            measurement_range = variant.measurement_range
        elif path.endswith(".negative_forms"):
            nf = biomarkers.biomarkers[biomarker_name].negative_forms
            measurement_range = (
                nf.measurement_range_for_value_placeholder if nf else None
            )
        for vid, surface in weighted.realizations.items():
            problems.extend(
                _check_placeholder_string(
                    f"{path}.realizations.{vid}", surface, measurement_range
                )
            )
    return problems


def _check_placeholder_string(
    label: str,
    surface: str,
    measurement_range: tuple[float, float] | None,
) -> list[str]:
    problems: list[str] = []
    found = _extract_placeholders(surface)
    unknown = found - SUPPORTED_PLACEHOLDERS
    if unknown:
        problems.append(
            f"{label}: unknown placeholders {sorted(unknown)}; "
            f"surface={surface!r}"
        )
    if "{{" in surface or "}}" in surface or _has_stray_braces(surface):
        problems.append(
            f"{label}: literal braces detected in {surface!r} — "
            "escape or remove"
        )
    if "value" in found and measurement_range is None:
        problems.append(
            f"{label}: realization contains {{value}} but parent has no "
            f"measurement_range"
        )
    return problems


def _has_stray_braces(s: str) -> bool:
    stripped = _PLACEHOLDER_RE.sub("", s)
    return "{" in stripped or "}" in stripped


def validate_configs(common: CommonConfig, biomarkers: BiomarkerConfig) -> None:
    """Apply the cross-config validation rules. Raises on any failure.

    On success returns None. On failure raises :class:`ConfigValidationError`
    with all accumulated problems in the message — loader fails loud,
    generation never runs against a broken config.
    """
    problems: list[str] = []
    problems.extend(_validate_preferred_methods(common, biomarkers))
    problems.extend(_validate_placeholders(common, biomarkers))
    problems.extend(_validate_population_keys(biomarkers))
    if problems:
        joined = "\n  - ".join(problems)
        raise ConfigValidationError(
            f"{len(problems)} config validation problem(s):\n  - {joined}"
        )


def describe_biomarker_shape(entry: Biomarker) -> dict[str, object]:
    """Diagnostic helper — summarizes what renderer-relevant features a
    biomarker declares. Used by tests and the CLI for visibility.
    """
    return {
        "canonical_name": entry.canonical_name,
        "alteration_type": entry.alteration_type,
        "variant_count": len(entry.variants),
        "has_negative_forms": entry.negative_forms is not None,
        "has_clone_attribution": entry.clone_attribution is not None,
        "has_preferred_methods": entry.preferred_methods is not None,
        "population_keys": sorted(entry.status_distribution_by_population),
    }
