# Complexity Levels

**Status:** WIP — Phase 1 (stub). L1 implemented; L2–L7 deferred to Phase 2+.

## Purpose

Define the seven complexity strata (L1–L7) formally with canonical examples, grammar requirements, and per-level evaluation protocols. This is the file reviewers will read to understand the "stratified F1" claim.

## Outline

- Motivation: why aggregate F1 hides structural failure modes (link to [CONVERSATION_CONTEXT.md](CONVERSATION_CONTEXT.md))
- L1 — Single assertion, canonical name, explicit status
- L2 — Lexical variation, abbreviations, synonyms
- L3 — List distribution, one status shared across multiple genes
- L4 — Heterogeneous compound (mixed statuses)
- L5 — Negation scope complexity
- L6 — Temporal / certainty qualification
- L7 — Cross-sentence coreference
- Per level: canonical example, frame count target, grammar building blocks, metric focus, known failure modes for generic models
- Promotion rules: bordering cases (e.g., when is an example L3 vs L4)
- Sampling budget per level per biomarker combination
- Adversarial category placement (distractors sit orthogonal to L1–L7)

Cross-reference: [PROJECT_SPEC.md Section 3 / Layer 4](PROJECT_SPEC.md)
