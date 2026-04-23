# Architecture

**Status:** WIP — Phase 1 (stub). Populated as each layer lands.

## Purpose

Document the seven-layer architecture from [project_spec.md Section 3](project_spec.md) with concrete data-flow diagrams, integration points, and per-layer ownership as each phase completes.

## Outline

- Overview diagram: assertion → sentence → corpus → eval
- Layer 1 — Ontology: loader contract, schema validation, extension protocol
- Layer 2 — Canonical assertion schema: Pydantic contract, versioning rules
- Layer 3 — Compositional grammar: sub-layers 3a–3e, frame design principles
- Layer 4 — Complexity stratification: per-level generator orchestration
- Layer 5 — Surface noise: application order, probability budgeting
- Layer 6 — Clinical realism: prevalence resolution, profile conditioning
- Layer 7 — Document composition: distractor injection, section structuring
- Cross-cutting — LLM paraphrase layer (Section 4 of spec)
- Integration contract between `bioassert` generator and `bioassert_eval` harness

Cross-reference: [project_spec.md Section 3](project_spec.md)
