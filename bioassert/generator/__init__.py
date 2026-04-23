"""Generator modules: assertion sampling, grammar rendering, adversarial records.

The probability-weighted pipeline driven by the JSON configs in
:mod:`bioassert.config`. Status distributions are flat per biomarker for
the cohort described by ``biomarkers.json`` — no patient-profile
conditioning.
"""
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
    sample_biomarker_name_form,
    sample_measurement_value,
    sample_method,
    sample_status,
    sample_variant,
    sample_variation,
)

__all__ = [
    "PostProcessError",
    "PostProcessedRecord",
    "RenderError",
    "RenderedRecord",
    "STATUS_NAMES",
    "apply_technical_noise",
    "maybe_sample_clone",
    "render_l1_record",
    "sample_biomarker_name_form",
    "sample_measurement_value",
    "sample_method",
    "sample_status",
    "sample_variant",
    "sample_variation",
]
