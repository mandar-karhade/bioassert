"""JSON loaders for the probability-weighted config (Phase 2a).

Reads ``common_variations.json`` and ``biomarkers.json`` from disk, strips the
top-level metadata keys (``$schema_version``, ``description``, ``_changelog``,
``_schema_template``, ``notes``), dispatches each remaining entry to the
right Pydantic model by inspecting its ``$schema_type``, and wraps the
parsed models plus raw metadata in :class:`CommonConfig` / :class:`BiomarkerConfig`
containers. Cross-object invariants (preferred_methods ↔ common.test_methods,
placeholder grammar) are applied by ``validator.validate_configs``.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Union

from pydantic import BaseModel, ConfigDict

from bioassert.config.schema import (
    Biomarker,
    CloneAttribution,
    PostProcessTransformations,
    TransformationCategory,
    WeightedVariations,
)

_METADATA_KEYS: frozenset[str] = frozenset(
    {
        "$schema_version",
        "description",
        "_changelog",
        "_schema_template",
        "notes",
    }
)


CommonCategory = Union[WeightedVariations, CloneAttribution, PostProcessTransformations]


class CommonConfig(BaseModel):
    """Parsed ``common_variations.json`` — categories keyed by name."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    schema_version: str
    categories: dict[str, CommonCategory]

    def test_method_realization(self, variation_id: str) -> str:
        """Look up a test-method realization string via the canonical source.

        Biomarker ``preferred_methods`` only supplies weights — realizations
        are joined in from ``common.test_methods.realizations`` (see
        CONFIG_ARCHITECTURE.md §2.4).
        """
        methods = self.categories.get("test_methods")
        if not isinstance(methods, WeightedVariations):
            raise KeyError("common.test_methods is missing or wrong schema_type")
        if variation_id not in methods.realizations:
            raise KeyError(
                f"test method {variation_id!r} not in common.test_methods"
            )
        return methods.realizations[variation_id]


class BiomarkerConfig(BaseModel):
    """Parsed ``biomarkers.json`` — biomarkers keyed by canonical symbol."""

    schema_version: str
    biomarkers: dict[str, Biomarker]

    def get(self, name: str) -> Biomarker:
        """Case-preserving lookup with a case-insensitive fallback."""
        if name in self.biomarkers:
            return self.biomarkers[name]
        name_upper = name.upper()
        for key, entry in self.biomarkers.items():
            if key.upper() == name_upper:
                return entry
        raise KeyError(f"biomarker {name!r} not found")


def _strip_metadata(raw: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in raw.items() if k not in _METADATA_KEYS}


def _dispatch_common_category(name: str, payload: dict[str, Any]) -> CommonCategory:
    if not isinstance(payload, dict):
        raise ValueError(
            f"common.{name}: expected dict, got {type(payload).__name__}"
        )
    schema_type = payload.get("$schema_type")
    if schema_type in (None, "weighted_variations"):
        return WeightedVariations.model_validate(payload)
    if schema_type == "weighted_variations_with_attachment":
        return CloneAttribution.model_validate(payload)
    if schema_type == "post_process_transformations":
        return _parse_post_process(name, payload)
    raise ValueError(
        f"common.{name}: unknown $schema_type {schema_type!r}"
    )


def _parse_post_process(
    name: str, payload: dict[str, Any]
) -> PostProcessTransformations:
    reserved = {"$schema_type", "description"}
    categories: dict[str, TransformationCategory] = {}
    for sub_name, sub_payload in payload.items():
        if sub_name in reserved:
            continue
        if not isinstance(sub_payload, dict):
            raise ValueError(
                f"common.{name}.{sub_name}: expected dict, got "
                f"{type(sub_payload).__name__}"
            )
        try:
            categories[sub_name] = TransformationCategory.model_validate(sub_payload)
        except Exception as exc:
            raise ValueError(
                f"common.{name}.{sub_name}: {exc}"
            ) from exc
    if not categories:
        raise ValueError(f"common.{name}: no transformation sub-categories declared")
    return PostProcessTransformations(
        **{
            "$schema_type": payload.get("$schema_type"),
            "description": payload.get("description"),
            "categories": categories,
        }
    )


def load_common(path: Union[Path, str]) -> CommonConfig:
    """Load and parse ``common_variations.json``.

    Structural validation (weights sum to 1.0, key parity) is enforced by the
    Pydantic models. Cross-config invariants run in ``validator.validate_configs``.
    """
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    schema_version = str(raw.get("$schema_version", ""))
    entries = _strip_metadata(raw)
    categories: dict[str, CommonCategory] = {}
    for name, payload in entries.items():
        categories[name] = _dispatch_common_category(name, payload)
    return CommonConfig(schema_version=schema_version, categories=categories)


def load_biomarkers(path: Union[Path, str]) -> BiomarkerConfig:
    """Load and parse ``biomarkers.json``."""
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    schema_version = str(raw.get("$schema_version", ""))
    entries = _strip_metadata(raw)
    biomarkers: dict[str, Biomarker] = {}
    for name, payload in entries.items():
        if not isinstance(payload, dict):
            raise ValueError(
                f"biomarker {name!r}: expected dict, got "
                f"{type(payload).__name__}"
            )
        biomarkers[name] = Biomarker.model_validate(payload)
    return BiomarkerConfig(schema_version=schema_version, biomarkers=biomarkers)


def load_configs(
    common_path: Union[Path, str], biomarkers_path: Union[Path, str]
) -> tuple[CommonConfig, BiomarkerConfig]:
    """Convenience loader that runs cross-config validation."""
    from bioassert.config.validator import validate_configs

    common = load_common(common_path)
    biomarkers = load_biomarkers(biomarkers_path)
    validate_configs(common, biomarkers)
    return common, biomarkers
