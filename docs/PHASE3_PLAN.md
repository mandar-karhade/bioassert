# Phase 3 — Sub-phase Plan

**Goal:** deliver compositional complexity levels L3–L7 + Layer-5 noise expansion, one sub-phase at a time, each gated by passing tests + a sanity-check corpus run before the next starts.

**Driving principle:** every sub-phase touches the minimum surface needed to land one architectural idea. If it requires schema changes, the schema change ships alone first. Grammar/frame work ships separately from schema work. Tests for new behavior land with the behavior.

---

## Sub-phase 3.0 — Grammar compatibility pass (DEFERRED)

Bugs 1 and 2 from `known_issues.md` represent realistic noise in real-world medical records — pathology reports in the wild contain arbitrarily broken grammar far beyond what rule-based generation can anticipate. A deterministic rule-based grammar-compatibility engine now would over-specify the signal. Revisit with a non-deterministic local-LLM paraphrase pass in a later phase (likely Phase 4 paraphrase layer), where the LLM can both normalize and introduce diverse, realistic rewrites at the same time.

**Status:** skipped for now. Both bugs remain OPEN in `known_issues.md` as accepted background noise for Phase 3.

---

## Sub-phase 3.1 — Multi-assertion record schema  **[SHIPPED — PR #1]**

**Why:** L3 and beyond need the record to carry a list of `(gene, variant, status, method, …)` tuples instead of a single-valued record. Landing the schema change alone keeps the diff reviewable.

**Scope:**
- New `AssertionFact` frozen dataclass — the tuple fields currently flat on `RenderedRecord` move here.
- `RenderedRecord` gains `assertions: tuple[AssertionFact, ...]`. Existing L1/L2 renderer wraps its single fact in a 1-tuple. No behavior change for L1/L2.
- Corpus JSONL serializer emits `assertions: [...]` always; single-assertion records write a 1-element list.
- Preservation invariant generalizes: every `AssertionFact` has a matching gene span in the sentence.

**Exit gate:** 87 tests still pass after adjusting the 2–3 that read single-valued fields directly. Fresh corpus identical except the `assertion` key is now `assertions` (1-element list).

**Risk:** medium — schema change touches serializer, tests, docs. Keep strictly backward-compatible at the data level (same record content, wrapping-only).

---

## Sub-phase 3.2 — L3 shared-status prose frames  **[SHIPPED — PR #3]**

**Why:** the simplest compound — N genes share one status in one sentence. All machinery from 3.1 reused. No per-gene status divergence yet.

**Scope:**
- New `_L3_FRAMES` family: `{gene_list} were {status}.`, `No {gene_list} alterations detected.`, `{gene_list}: {status}.` with list-coordinator rendering (comma + Oxford-optional + "and" / "or").
- Sampler change: L3 record samples N genes from the panel (N ∈ {2,3,4}, configurable), one shared status, one shared frame.
- Status coherence: N gene-status tuples all get the same status; prevalence conditioning runs once per gene but status is drawn once and applied.
- Spans: one span per gene + one shared status span. Each `AssertionFact` points at its gene's span.

**Exit gate:** new L3 unit tests — every listed gene is a literal substring; shared status is consistent across assertions; L3 renders never include variants yet (deferred). Regen 50K corpus with `--l3-fraction 0.2` and verify 0 span violations.

**Risk:** medium — list coordinator is new; span tracking for N tokens needs careful `_render_with_spans` extension.

**Open design questions:** does L3 allow variant per gene? (Initial answer: no — L3 is "bare gene names only" shared status; variants land in L4/L5.) Can L3 mix biomarker types (mutation + expression in one list)? (Initial answer: no — homogeneous biomarker class per L3 record.)

---

## Sub-phase 3.3 — L3 shorthand / tabular  **[SHIPPED — PR #5]**

**Why:** the L2 analogue for L3. Tabular reports often list multiple genes with shorthand statuses per line.

**Scope:**
- New `_L3_SHORTHAND_FRAMES`: `EGFR -, ALK -, ROS1 -`, tab-separated variant, newline-separated variant.
- Uses `positive_shorthand` / `negation_shorthand` from Sub-phase 3.0's cleaned vocab.
- Shared status — same shorthand applied to each gene.

**Exit gate:** L3-shorthand unit tests; corpus shows L3-shorthand mixing with L1/L2/L3-prose per the fraction flags.

**Risk:** low — parallels L2 shorthand, schema already done.

---

## Sub-phase 3.4 — L4 heterogeneous compound (prose)  **[SHIPPED — PR #6]**

**Why:** same list structure as L3 but each gene carries independent status. This is where extraction models really start breaking.

**Scope:**
- `_L4_FRAMES`: `{gene1} {status1}, {gene2} {status2}, and {gene3} {status3}.` and semicolon-separated variants.
- Sampler: per-gene prevalence conditioning (one status per gene, independently drawn — reuses the existing per-biomarker distribution).
- Spans: N gene spans + N status spans. Each `AssertionFact` gets its own gene + status span.

**Exit gate:** L4 unit tests — every gene paired with its own status span, mixed-status record never conflates statuses. Regen corpus with L4 fraction.

**Risk:** medium — more spans to track, more room for off-by-one in separator positions.

---

## Sub-phase 3.5 — L4 shorthand / tabular  **[SHIPPED — PR #7]**

Parallel to 3.3 — shorthand form of heterogeneous multi-gene. `EGFR +, ALK -, ROS1 WT` style.

**Risk:** low.

---

## Sub-phase 3.6 — L5 negation scope  **[SHIPPED — PR #8]**

**Why:** a qualitatively different primitive — scope markers ("no", "absence of") apply to multiple entities; exception markers ("but", "except") flip polarity mid-sentence.

**Scope:**
- New frame family: `No {gene_list} mutations but {gene} {variant} detected.` / `Absence of {gene1}, {gene2} alterations; {gene3} {status}.`
- `AssertionFact` gains a `polarity_scope: Literal["direct", "negation_wide", "exception"]` field.
- Renderer tracks which genes fall inside the scope vs. after the exception.
- Status per assertion is computed from scope × base-status (negation-wide flips all to negative; exception segment carries its own status).

**Exit gate:** L5 tests verifying polarity is correctly attached per assertion across scope boundaries.

**Risk:** high — polarity bookkeeping is the trickiest part of Phase 3. Expect 1–2 iterations.

---

## Sub-phase 3.7 — Compounding tiers (low / high)  **[SHIPPED — PR #9]**

**Pivot note:** the original plan slotted L6 temporal/certainty here. Mid-Phase-3 user feedback (2026-04-23) redirected to compounding complexity first, because the tier dimension is orthogonal to every L3–L5 surface already shipped and unlocks more varied compound records without adding a new complexity level. L6 / L7 / noise expansion move down one slot each.

**Why:** L3/L3S/L4/L4S/L5 records previously used a fixed 2–4 gene range. Clinical reports vary from terse "KRAS and EGFR negative" to panel summaries listing 7–8 biomarkers. A tier knob lets the corpus mix both without inventing new frames.

**Scope (as shipped):**
- `RenderedRecord` / `PostProcessedRecord` gain `compounding_tier: str` (default `"low"`).
- Two tiers only — medium was collapsed into high per user feedback:
  - `low` = exactly 2 genes (minimum compound)
  - `high` = 3–8 genes, clamped to pool size; cross-class pool mixing (TMB/PD-L1 + mutation biomarkers) allowed
- `COMPOUND_LOW_N = 2`, `COMPOUND_HIGH_MIN = 3`, `COMPOUND_HIGH_MAX = 8`, `_TIER_RANGES`, `n_range_for_tier(tier)` helper.
- `render_l3_record` / `render_l4_record` / `render_l5_record` (and shorthand variants) accept `compounding_tier` and stamp it on the returned record. L5 panel-wide frames ignore tier (always 1 labeled fact).
- Corpus script CLI: `--compound-low` / `--compound-high` (default 0.5/0.5, required to sum to 1.0). Report emits `compound_*_fraction_requested` / `_observed` plus per-level tier breakdown.

**Exit gate (met):** 181 tests passing; 50K regen at `datasets/v1_phase3.7/` showing observed 0.5016 low / 0.4984 high, 107/112 within-2σ, 0 span violations; eyeball review approved.

**Risk:** low — additive parameter on existing renderers; medium pool-mixing decision landed cleanly.

---

## Sub-phase 3.8 — L6 temporal / certainty qualification

**Scope:**
- Optional `temporal: str | None` and `certainty: str | None` on each `AssertionFact`.
- Frame family: `Previously {status} for {gene}; current testing {certainty} {status}.`
- Certainty vocabulary: suspected, confirmed, probable, rule-out, pending.
- Temporal vocabulary: previously, currently, at diagnosis, at relapse, post-TKI.

**Exit gate:** L6 tests — qualifiers attached per assertion; confirmation-by-inspection sample.

**Risk:** medium — new vocab categories + optional qualifier slots in record schema.

---

## Sub-phase 3.9 — L7 cross-sentence coreference

**Scope:**
- Record becomes multi-sentence: `sentence` → `sentences: tuple[str, ...]`. Each `AssertionFact` carries a `sentence_index: int` into the tuple.
- Discourse frames: setup → claim → qualifier/distractor (3-sentence typical).
- Coreference devices: pronoun resolution ("Results showed…"), nominal anaphora ("the mutation", "this variant").

**Exit gate:** L7 tests ensure every assertion's spans are in the sentence referenced by its `sentence_index`; coreference frames exercise both pronoun and nominal forms.

**Risk:** high — multi-sentence serialization is a structural shift; rendering bugs get harder to diagnose across sentence boundaries.

---

## Sub-phase 3.10 — Layer 5 noise expansion + Bug 3a fix

**Why:** addresses `known_issues.md` Bug 3a (whitespace noise inside L1 prose) and extends Layer 5 per the spec's Section 3 Layer 5.

**Scope:**
- Fix Bug 3a: restrict whitespace substitution to slot-boundary positions (tracked via `_render_with_spans`), never interior spaces.
- Add OCR-style corruption (optional — configurable).
- Add PDF-extraction artifacts (hyphen-newline line breaks within tokens).
- Abbreviation inconsistency within a single record.

**Exit gate:** Bug 3a rate drops to 0% on L1 prose; noise fractions logged in the prevalence report.

**Risk:** low — purely post-process.

---

## What this does NOT cover

- **Paraphrase layer** (LLM surface rewrite) — Phase 4.
- **Document composition** (multi-paragraph pathology notes with distractors) — Phase 4.
- **Evaluation harness** (`bioassert_eval`, GLiNER / REBEL / ReLiK adapters, stratified F1 reporter) — Phase 5+.

---

## Execution rhythm

Each sub-phase follows the same loop:
1. Branch off `main` as `phase3.N/<slug>`.
2. Update `docs/PHASE3_PLAN.md` with any scope-shift decisions made during implementation.
3. TDD: write failing tests for the new behavior.
4. Implement to green.
5. Run full pytest.
6. Generate the 50K corpus + prevalence report.
7. **Surface a sample of the newly generated corpus to the user and wait for approval.** Tests verify code correctness, not data correctness — the user must eyeball representative records before merge.
8. After explicit user approval, merge the PR to `main` with `gh pr merge --merge` (no squash, no branch delete).
9. Commit prefix: `feat(phase3.N):` or `refactor(phase3.N):`.

Update the top-of-spec status line as each sub-phase ships.
