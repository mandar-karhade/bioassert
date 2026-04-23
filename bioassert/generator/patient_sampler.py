"""Patient profile sampling and population-key cascade (Phase 2a).

A :class:`PatientProfile` carries the five canonical population axes. The
cascade (:meth:`PatientProfile.population_key_cascade`) emits candidate keys
in decreasing specificity per config_architecture.md §2.3. The sampler pulls
status distributions for a biomarker by walking the cascade until it finds
a configured key; terminal fallback is always ``{histology}``.
"""
from __future__ import annotations

import logging
import random
from dataclasses import dataclass
from typing import Literal, Optional

from bioassert.config.schema import Biomarker, PopulationStatus

log = logging.getLogger(__name__)

Histology = Literal["adenocarcinoma", "squamous", "other"]
Ethnicity = Literal["western", "east_asian", "other"]
Smoking = Literal["smoker", "nonsmoker", "former_smoker"]
AgeGroup = Literal["younger", "older"]
Sex = Literal["female", "male"]


@dataclass(frozen=True)
class PatientProfile:
    """Composable five-axis patient profile.

    Canonical key order: histology → ethnicity → smoking → age_group → sex.
    Optional axes are omitted from the generated key when ``None``.
    """

    patient_ref: str
    histology: Histology
    ethnicity: Optional[Ethnicity] = None
    smoking: Optional[Smoking] = None
    age_group: Optional[AgeGroup] = None
    sex: Optional[Sex] = None

    def population_key_cascade(self) -> list[str]:
        """Return candidate population keys in decreasing specificity.

        Drop order follows §2.3: sex → age_group → smoking → ethnicity → (empty,
        leaving histology alone as the terminal fallback). Duplicate keys
        produced when an axis is ``None`` are removed while preserving order.
        """
        optional_axes: list[tuple[str, Optional[str]]] = [
            ("ethnicity", self.ethnicity),
            ("smoking", self.smoking),
            ("age_group", self.age_group),
            ("sex", self.sex),
        ]
        candidates: list[str] = []
        for drop_from_end in range(len(optional_axes) + 1):
            kept = optional_axes[: len(optional_axes) - drop_from_end]
            populated = [value for _, value in kept if value is not None]
            key = "_".join([self.histology, *populated])
            if key not in candidates:
                candidates.append(key)
        return candidates


class PopulationCascadeMiss(KeyError):
    """Raised when even the histology fallback is missing (loader bug)."""


_WARNED: set[tuple[str, str]] = set()


def resolve_status_distribution(
    biomarker: Biomarker, profile: PatientProfile
) -> tuple[PopulationStatus, str]:
    """Walk the cascade, return the first matching distribution and its key.

    Emits a one-time WARNING via :mod:`logging` the first time a
    ``(biomarker, profile-signature)`` combination has to fall back past the
    most specific candidate. When no histology-based key matches (typical
    for ``histology='other'`` patients since biomarker configs only guarantee
    ``adenocarcinoma``/``squamous``), the ultimate fallback is the configured
    ``adenocarcinoma`` row — the schema validator enforces that it exists.
    See config_architecture.md §2.3.
    """
    cascade = profile.population_key_cascade()
    pops = biomarker.status_distribution_by_population
    for idx, key in enumerate(cascade):
        if key in pops:
            if idx > 0:
                sig = (biomarker.canonical_name, cascade[0])
                if sig not in _WARNED:
                    _WARNED.add(sig)
                    log.warning(
                        "population cascade fallback: biomarker=%r "
                        "requested=%r matched=%r (depth=%d)",
                        biomarker.canonical_name,
                        cascade[0],
                        key,
                        idx,
                    )
            return pops[key], key
    if "adenocarcinoma" in pops:
        sig = (biomarker.canonical_name, cascade[0])
        if sig not in _WARNED:
            _WARNED.add(sig)
            log.warning(
                "population histology fallback: biomarker=%r requested=%r "
                "matched='adenocarcinoma' (no histology-cascade overlap)",
                biomarker.canonical_name,
                cascade[0],
            )
        return pops["adenocarcinoma"], "adenocarcinoma"
    raise PopulationCascadeMiss(
        f"{biomarker.canonical_name}: no population-key match for "
        f"cascade={cascade}"
    )


DEFAULT_PROFILE_DISTRIBUTION: dict[str, dict[str, float]] = {
    "histology": {"adenocarcinoma": 0.80, "squamous": 0.18, "other": 0.02},
    "ethnicity": {"western": 0.78, "east_asian": 0.18, "other": 0.04},
    "smoking": {"smoker": 0.55, "former_smoker": 0.25, "nonsmoker": 0.20},
    "age_group": {"older": 0.72, "younger": 0.28},
    "sex": {"male": 0.52, "female": 0.48},
}


def _sample_axis(
    rng: random.Random, distribution: dict[str, float], label: str
) -> str:
    keys = list(distribution)
    weights = [distribution[k] for k in keys]
    total = sum(weights)
    if abs(total - 1.0) > 1e-6:
        raise ValueError(
            f"{label} distribution sums to {total}, expected 1.0"
        )
    return rng.choices(keys, weights=weights, k=1)[0]


def sample_patient_profile(
    patient_ref: str,
    rng: Optional[random.Random] = None,
    distribution: Optional[dict[str, dict[str, float]]] = None,
) -> PatientProfile:
    """Sample a PatientProfile from NSCLC-plausible axis priors.

    Priors are rough clinical-epidemiology defaults usable for Phase 2a smoke
    testing; Phase 2b will swap in citation-backed distributions. Pass
    ``distribution`` to override per-axis priors for deterministic tests.
    """
    rng = rng or random.Random()
    dist = distribution or DEFAULT_PROFILE_DISTRIBUTION
    return PatientProfile(
        patient_ref=patient_ref,
        histology=_sample_axis(rng, dist["histology"], "histology"),
        ethnicity=_sample_axis(rng, dist["ethnicity"], "ethnicity"),
        smoking=_sample_axis(rng, dist["smoking"], "smoking"),
        age_group=_sample_axis(rng, dist["age_group"], "age_group"),
        sex=_sample_axis(rng, dist["sex"], "sex"),
    )
