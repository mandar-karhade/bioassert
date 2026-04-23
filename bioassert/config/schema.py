"""Pydantic models for the Phase 2a probability-weighted config.

Shape follows config_architecture.md Sections 2.1–2.4 and 8. Structural
invariants (sum-to-1 distributions, matched variations/realizations key sets,
positivity bounds) are enforced here at load time via model validators. Cross-
object invariants (preferred_methods ⊂ common.test_methods, render_constraints
resolvability, placeholder grammar) live in ``validator.py`` because they need
both configs.

Every weighted distribution sums to 1.0 within a tolerance of
:data:`WEIGHT_TOLERANCE`. Failures at this layer raise ``ValidationError``;
loader fails loudly, generation never sees a malformed config.
"""
from __future__ import annotations

from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator

WEIGHT_TOLERANCE: float = 0.001


AlterationType = Literal[
    "mutation", "fusion", "amplification", "expression", "composite"
]
StatusName = Literal["positive", "negative", "equivocal", "not_tested"]


def _assert_weights_sum_to_one(weights: dict[str, float], label: str) -> None:
    if not weights:
        raise ValueError(f"{label}: empty weighted distribution")
    total = sum(weights.values())
    if abs(total - 1.0) > WEIGHT_TOLERANCE:
        raise ValueError(
            f"{label}: weights sum to {total:.6f}, expected 1.0 "
            f"(tolerance {WEIGHT_TOLERANCE})"
        )
    for k, v in weights.items():
        if not 0.0 <= v <= 1.0:
            raise ValueError(
                f"{label}: weight for {k!r} = {v} is outside [0, 1]"
            )


def _assert_keys_match(
    a: dict[str, object], b: dict[str, object], a_label: str, b_label: str
) -> None:
    if set(a.keys()) != set(b.keys()):
        missing_in_b = set(a.keys()) - set(b.keys())
        missing_in_a = set(b.keys()) - set(a.keys())
        raise ValueError(
            f"key mismatch: {a_label} vs {b_label}; "
            f"missing in {b_label}: {sorted(missing_in_b)}; "
            f"missing in {a_label}: {sorted(missing_in_a)}"
        )


class WeightedVariations(BaseModel):
    """Probability-weighted surface variations.

    Schema type: ``weighted_variations``. The ``variations`` dict supplies the
    prior over variation IDs, and ``realizations`` maps each ID to its surface
    string. Realization strings may contain placeholder tokens (``{gene}``,
    ``{method}``, ``{result}``, ``{value}``) that the renderer fills at render
    time — see config_architecture.md §7.5.
    """

    model_config = ConfigDict(populate_by_name=True, extra="forbid")

    schema_type: Optional[str] = Field(default=None, alias="$schema_type")
    description: Optional[str] = None
    variations: dict[str, float]
    realizations: dict[str, str]

    @model_validator(mode="after")
    def _validate(self) -> "WeightedVariations":
        label = self.description or "WeightedVariations"
        _assert_weights_sum_to_one(self.variations, label)
        _assert_keys_match(
            self.variations,
            self.realizations,
            f"{label}.variations",
            f"{label}.realizations",
        )
        return self


class CloneAttribution(BaseModel):
    """Independent IHC clone attribution for PD-L1 and similar biomarkers.

    Schema type: ``weighted_variations_with_attachment``. Sampling is two-stage:
    first a Bernoulli roll with ``attachment_probability``, then (if attached)
    a weighted draw from ``variations``. See config_architecture.md §8.4.
    """

    model_config = ConfigDict(populate_by_name=True, extra="forbid")

    schema_type: Optional[str] = Field(default=None, alias="$schema_type")
    description: Optional[str] = None
    attachment_probability: float
    variations: dict[str, float]
    realizations: dict[str, str]

    @model_validator(mode="after")
    def _validate(self) -> "CloneAttribution":
        if not 0.0 <= self.attachment_probability <= 1.0:
            raise ValueError(
                f"clone_attribution.attachment_probability "
                f"{self.attachment_probability} outside [0, 1]"
            )
        label = self.description or "CloneAttribution"
        _assert_weights_sum_to_one(self.variations, label)
        _assert_keys_match(
            self.variations,
            self.realizations,
            f"{label}.variations",
            f"{label}.realizations",
        )
        return self


class TransformationCategory(BaseModel):
    """One sub-category inside ``technical_noise`` (whitespace, case, etc.)."""

    model_config = ConfigDict(extra="forbid")

    description: Optional[str] = None
    distribution: dict[str, float]

    @model_validator(mode="after")
    def _validate(self) -> "TransformationCategory":
        label = self.description or "TransformationCategory"
        _assert_weights_sum_to_one(self.distribution, label)
        return self


class PostProcessTransformations(BaseModel):
    """Technical-noise sub-categories keyed by name.

    Schema type: ``post_process_transformations``. Each sub-category
    (``whitespace``, ``case_variation``, …) is its own weighted distribution
    applied probabilistically AFTER semantic composition. See
    config_architecture.md §7.6.
    """

    model_config = ConfigDict(populate_by_name=True, extra="forbid")

    schema_type: Optional[str] = Field(default=None, alias="$schema_type")
    description: Optional[str] = None
    categories: dict[str, TransformationCategory]


class RenderConstraints(BaseModel):
    """Coupling hints that constrain sampling of sibling dimensions.

    Currently only ``require_biomarker_name_forms`` is defined: when a variant
    with this constraint is sampled, the biomarker-level name_form sampler is
    restricted to the listed subset (weights renormalized). Used to keep
    fusion ``unspecified_fusion`` variants legible — see
    config_architecture.md §8.3.
    """

    model_config = ConfigDict(extra="forbid")

    require_biomarker_name_forms: list[str] = Field(default_factory=list)


class Variant(BaseModel):
    """One sub-variant within a biomarker (e.g., L858R under EGFR)."""

    model_config = ConfigDict(extra="forbid")

    actionable: bool
    prevalence_within_biomarker: float
    category: Optional[str] = None
    notes: Optional[str] = None
    name_forms: WeightedVariations
    measurement_range: Optional[tuple[float, float]] = None
    render_constraints: Optional[RenderConstraints] = None

    @model_validator(mode="after")
    def _validate(self) -> "Variant":
        if not 0.0 <= self.prevalence_within_biomarker <= 1.0:
            raise ValueError(
                f"prevalence_within_biomarker "
                f"{self.prevalence_within_biomarker} outside [0, 1]"
            )
        if self.measurement_range is not None:
            lo, hi = self.measurement_range
            if lo > hi:
                raise ValueError(
                    f"measurement_range {self.measurement_range} has min > max"
                )
        return self


class NegativeForms(BaseModel):
    """Expression-biomarker negative renderings (PD-L1 TPS 0%, TMB-low, …).

    Treated as a ``weighted_variations`` block that is sampled from when
    ``status == "negative"``, separate from the positive-direction ``variants``
    list. See config_architecture.md §8.1.
    """

    model_config = ConfigDict(populate_by_name=True, extra="forbid")

    schema_type: Optional[str] = Field(default=None, alias="$schema_type")
    description: Optional[str] = None
    variations: dict[str, float]
    realizations: dict[str, str]
    measurement_range_for_value_placeholder: Optional[tuple[float, float]] = None

    @model_validator(mode="after")
    def _validate(self) -> "NegativeForms":
        label = self.description or "NegativeForms"
        _assert_weights_sum_to_one(self.variations, label)
        _assert_keys_match(
            self.variations,
            self.realizations,
            f"{label}.variations",
            f"{label}.realizations",
        )
        if self.measurement_range_for_value_placeholder is not None:
            lo, hi = self.measurement_range_for_value_placeholder
            if lo > hi:
                raise ValueError(
                    "measurement_range_for_value_placeholder min > max"
                )
        return self


class PopulationStatus(BaseModel):
    """Status-distribution row for one patient sub-population."""

    model_config = ConfigDict(extra="forbid")

    positive: float
    negative: float
    equivocal: float
    not_tested: float

    @model_validator(mode="after")
    def _validate(self) -> "PopulationStatus":
        total = self.positive + self.negative + self.equivocal + self.not_tested
        if abs(total - 1.0) > WEIGHT_TOLERANCE:
            raise ValueError(
                f"PopulationStatus probabilities sum to {total:.6f}, "
                f"expected 1.0 (tolerance {WEIGHT_TOLERANCE})"
            )
        for name, v in (
            ("positive", self.positive),
            ("negative", self.negative),
            ("equivocal", self.equivocal),
            ("not_tested", self.not_tested),
        ):
            if not 0.0 <= v <= 1.0:
                raise ValueError(f"PopulationStatus.{name} = {v} outside [0, 1]")
        return self


class PreferredMethods(BaseModel):
    """Biomarker-specific method-weights override.

    Weights only — realizations are canonical and come from
    ``common.test_methods.realizations``. See config_architecture.md §2.4.
    """

    model_config = ConfigDict(extra="forbid")

    variations: dict[str, float]

    @model_validator(mode="after")
    def _validate(self) -> "PreferredMethods":
        _assert_weights_sum_to_one(self.variations, "PreferredMethods.variations")
        return self


class Biomarker(BaseModel):
    """One biomarker entry: name forms, variants, population priors, methods."""

    model_config = ConfigDict(extra="forbid")

    canonical_name: str
    alteration_type: AlterationType
    gene_location: Optional[str] = None
    notes: Optional[str] = None
    name_forms: WeightedVariations
    variants: dict[str, Variant]
    negative_forms: Optional[NegativeForms] = None
    clone_attribution: Optional[CloneAttribution] = None
    status_distribution_by_population: dict[str, PopulationStatus]
    preferred_methods: Optional[PreferredMethods] = None

    @model_validator(mode="after")
    def _validate(self) -> "Biomarker":
        variant_prev = {vid: v.prevalence_within_biomarker for vid, v in self.variants.items()}
        _assert_weights_sum_to_one(variant_prev, f"{self.canonical_name}.variants")

        pops = self.status_distribution_by_population
        if "adenocarcinoma" not in pops:
            raise ValueError(
                f"{self.canonical_name}: missing terminal fallback "
                f"'adenocarcinoma' in status_distribution_by_population"
            )
        if "squamous" not in pops:
            raise ValueError(
                f"{self.canonical_name}: missing terminal fallback 'squamous' "
                f"in status_distribution_by_population"
            )

        for variant_id, variant in self.variants.items():
            if variant.render_constraints is not None:
                allowed = set(variant.render_constraints.require_biomarker_name_forms)
                if not allowed:
                    raise ValueError(
                        f"{self.canonical_name}.{variant_id}.render_constraints:"
                        " require_biomarker_name_forms is empty"
                    )
                unknown = allowed - set(self.name_forms.variations.keys())
                if unknown:
                    raise ValueError(
                        f"{self.canonical_name}.{variant_id}.render_constraints:"
                        f" unknown name_form keys {sorted(unknown)}"
                    )

        return self
