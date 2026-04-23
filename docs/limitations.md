# Limitations

**Status:** WIP — Phase 1 (stub). Grows in honesty as scope grows. Non-Negotiable #6.

## Purpose

Honestly document what the corpus does NOT cover. Reviewers will read this before accepting the paper. Better to preempt criticism than to be surprised by it.

## Outline

- Surface-level artifacts NOT modeled
  - OCR-style corruption beyond simple character swaps
  - Handwriting and handwritten-scan artifacts
  - PDF multi-column interleaving and table-to-text extraction failures
- Linguistic scope
  - US-English clinical writing conventions only
  - No non-English pathology reports
  - Abbreviation conflicts across specialties not modeled (e.g., "MS" as multiple sclerosis vs mass spec)
- Domain scope
  - NSCLC only for v1; no breast, colorectal, prostate, hematologic
  - No germline / hereditary cancer syndromes
  - No pharmacogenomic assertions (CYP2D6, TPMT, etc.)
- Clinical nuance NOT captured
  - Pathologist idiosyncratic phrasing
  - Institution-specific report templates
  - Longitudinal progression narratives beyond L6 temporality
- Ontology gaps
  - Variant-level prevalence uniform in Phase 1 (skewed in reality)
  - Co-mutation dependencies not modeled until Phase 2
- LLM paraphrase layer caveats
  - Span preservation is necessary but not sufficient for semantic preservation
  - Potential to introduce subtle polarity drift; mitigated by regex verification but not eliminated
- Evaluation scope
  - Real-text validation limited to public datasets until collaborator agreement
