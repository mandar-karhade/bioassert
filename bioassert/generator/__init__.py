"""Generator modules: assertion sampling, grammar rendering, adversarial records.

Phase 2a additions: :mod:`patient_sampler`, :mod:`sampler`, :mod:`renderer`,
:mod:`post_process` — the probability-weighted pipeline driven by the JSON
configs in :mod:`bioassert.config`.
"""
from bioassert.generator.patient_sampler import (
    PatientProfile,
    PopulationCascadeMiss,
    resolve_status_distribution,
    sample_patient_profile,
)
from bioassert.generator.post_process import (
    PostProcessError,
    PostProcessedRecord,
    apply_technical_noise,
)
from bioassert.generator.renderer import (
    RenderError,
    RenderedRecord,
    render_l1_record,
)
from bioassert.generator.sampler import (
    STATUS_NAMES,
    maybe_sample_clone,
    resolve_population,
    sample_biomarker_name_form,
    sample_measurement_value,
    sample_method,
    sample_status,
    sample_variant,
    sample_variation,
)

__all__ = [
    "PatientProfile",
    "PopulationCascadeMiss",
    "PostProcessError",
    "PostProcessedRecord",
    "RenderError",
    "RenderedRecord",
    "STATUS_NAMES",
    "apply_technical_noise",
    "maybe_sample_clone",
    "render_l1_record",
    "resolve_population",
    "resolve_status_distribution",
    "sample_biomarker_name_form",
    "sample_measurement_value",
    "sample_method",
    "sample_patient_profile",
    "sample_status",
    "sample_variant",
    "sample_variation",
]
