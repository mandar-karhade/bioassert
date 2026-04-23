# Claude Code Kickoff Prompt

Copy this into a new Claude Code session to begin development. Paste `PROJECT_SPEC.md` and `CONVERSATION_CONTEXT.md` as context files alongside.

---

## Prompt to paste

> I'm starting a new project: a synthetic corpus generator for biomarker assertion extraction in lung cancer (NSCLC). The full specification is in `PROJECT_SPEC.md` and the conversation backstory is in `CONVERSATION_CONTEXT.md`. Read both carefully before writing any code.
>
> The core architectural principle — never invert this — is: **ground truth is the structured assertion. Surface text is rendered from assertions. Labels are correct by construction.**
>
> For this first session, produce the minimum viable end-to-end slice defined in Section 13 of PROJECT_SPEC.md:
>
> 1. Repo scaffold with `pyproject.toml` (use `uv` or `hatch` for packaging), directory layout matching Section 9
> 2. `ontology/biomarkers.yaml` with 3 biomarkers fully profiled as the schema reference: **EGFR, ALK, KRAS**. Include all fields from the Layer 1 example (canonical, abbreviations, spelling_variants, gene_family, common_variants with aliases, test_methods, result_vocabulary, measurement_types, clinical_associations, prevalence with multiple subpopulations, prevalence_sources). Use published NSCLC molecular epidemiology as sources.
> 3. `bioassert/schema.py` with Pydantic `BiomarkerAssertion` and `Measurement` models exactly matching Section 3 Layer 2 of the spec
> 4. `bioassert/ontology.py` — loads YAML, validates against schema, exposes `get_biomarker(name)` and `sample_variant(gene)` helpers
> 5. `bioassert/generator/assertion_sampler.py` — prevalence-aware sampling. Given a simulated patient profile (histology, ethnicity), returns a list of `BiomarkerAssertion` instances with statuses sampled from the prevalence distribution for that subpopulation
> 6. `bioassert/generator/grammar.py` — 10–15 L1 sentence frames. Given an assertion, renders a surface sentence. Must also return the character spans for gene, status, (optionally) variant so labels are preserved
> 7. `scripts/generate_corpus.py` — emits 1,000 L1 assertions + rendered sentences as JSONL to `datasets/v0_smoke/corpus.jsonl`. Each record contains the full assertion object AND the rendered sentence AND the character spans
> 8. `tests/test_assertion_preservation.py` — pytest tests validating that for every generated record, the gene name, status term, and variant (if present) are actually found at the claimed character spans in the sentence. This test is the guardrail for Non-Negotiable #1
>
> Use Python 3.11+. Type-annotate everything. Pydantic v2. Pytest. No external LLM calls in this phase — pure rule-based generation.
>
> Do not build beyond this scope in this session. Once the smoke test passes with 1,000 records, we stop and validate the architecture end-to-end before adding L2+ complexity.
>
> Start by reading both context documents, then confirm your understanding of (a) the non-inversion principle and (b) the seven-layer architecture before writing any code.

---

## What to expect / how to validate the first session output

When Claude Code finishes, validate these before moving to Phase 2:

- [ ] `pytest tests/` passes with all assertion-preservation tests green
- [ ] `python scripts/generate_corpus.py` emits 1,000 records in `datasets/v0_smoke/corpus.jsonl`
- [ ] Each record has: full assertion object, rendered sentence, character spans
- [ ] Manual spot-check: read 20 random records. Do the sentences sound like plausible clinical text? Do the spans actually match?
- [ ] EGFR prevalence in the generated corpus matches the configured rate (~15% Western, ~45% East Asian) within sampling noise for a 1,000-record corpus
- [ ] No two identical sentences appear more than ~1% of the time (sanity check that the grammar has real variation)

If any of these fail, fix before proceeding. The whole project rests on this foundation.

---

## Phase 2 kickoff (for the next session)

After Phase 1 validates, the next Claude Code session should:
- Extend ontology to 15–20 biomarkers
- Add L2 (lexical variation) sentence frames
- Add L3 (list distribution) compound grammar
- Add L4 (heterogeneous compound) grammar
- Scale test corpus to 10K records stratified across L1–L4

Keep sessions tightly scoped. Don't let any single session try to build L1 through L7 at once — the grammar design for each level requires its own thought.
