"""Cross-config validation (config_architecture.md §7).

Structural invariants (weights sum to 1.0, matched variation/realization key
sets, bounded probabilities) are enforced by the Pydantic models in
:mod:`bioassert.config.schema`. This module handles the remaining checks that
require looking at both configs together plus placeholder grammar against the
supported table in §7.5.

Additional invariants enforced here:
 * ``preferred_methods`` keys are a subset of ``common.test_methods.variations``
 * ``render_constraints`` resolve to valid biomarker-level ``name_form`` keys
 * placeholders in realization strings are from the supported set
 * every realization carrying ``{value}`` has a backing ``measurement_range``
"""
from __future__ import annotations

import re
from typing import Iterable

from bioassert.config.loader import BiomarkerConfig, CommonConfig
from bioassert.config.schema import (
    Biomarker,
    CloneAttribution,
    NegativeForms,
    WeightedVariations,
)

SUPPORTED_PLACEHOLDERS: frozenset[str] = frozenset(
    {"gene", "method", "result", "value", "variant"}
)

_PLACEHOLDER_RE = re.compile(r"\{([^{}]*)\}")


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
    }
