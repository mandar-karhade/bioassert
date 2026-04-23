# Phase 3 — Sub-phase Plan

**Goal:** deliver compositional complexity levels L3–L7 + Layer-5 noise expansion, one sub-phase at a time, each gated by passing tests + a sanity-check corpus run before the next starts.

**Driving principle:** every sub-phase touches the minimum surface needed to land one architectural idea. If it requires schema changes, the schema change ships alone first. Grammar/frame work ships separately from schema work. Tests for new behavior land with the behavior.

---

## Sub-phase 3.0 — Grammar compatibility pass (DEFERRED)

Bugs 1 and 2 from `known_issues.md` represent realistic noise in real-world medical records — pathology reports in the wild contain arbitrarily broken grammar far beyond what rule-based generation can anticipate. A deterministic rule-based grammar-compatibility engine now would over-specify the signal. Revisit with a non-deterministic local-LLM paraphrase pass in a later phase (likely Phase 4 paraphrase layer), where the LLM can both normalize and introduce diverse, realistic rewrites at the same time.

**Status:** skipped for now. Both bugs remain OPEN in `known_issues.md` as accepted background noise for Phase 3.

---

## Sub-phase 3.1 — Multi-assertion record schema

**Why:** L3 and beyond need the record to carry a list of `(gene, variant, status, method, …)` tuples instead of a single-valued record. Landing the schema change alone keeps the diff reviewable.

**Scope:**
- New `AssertionFact` frozen dataclass — the tuple fields currently flat on `RenderedRecord` move here.
- `RenderedRecord` gains `assertions: tuple[AssertionFact, ...]`. Existing L1/L2 renderer wraps its single fact in a 1-tuple. No behavior change for L1/L2.
- Corpus JSONL serializer emits `assertions: [...]` always; single-assertion records write a 1-element list.
- Preservation invariant generalizes: every `AssertionFact` has a matching gene span in the sentence.

**Exit gate:** 87 tests still pass after adjusting the 2–3 that read single-valued fields directly. Fresh corpus identical except the `assertion` key is now `assertions` (1-element list).

**Risk:** medium — schema change touches serializer, tests, docs. Keep strictly backward-compatible at the data level (same record content, wrapping-only).

---

## Sub-phase 3.2 — L3 shared-status prose frames

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

## Sub-phase 3.3 — L3 shorthand / tabular

**Why:** the L2 analogue for L3. Tabular reports often list multiple genes with shorthand statuses per line.

**Scope:**
- New `_L3_SHORTHAND_FRAMES`: `EGFR -, ALK -, ROS1 -`, tab-separated variant, newline-separated variant.
- Uses `positive_shorthand` / `negation_shorthand` from Sub-phase 3.0's cleaned vocab.
- Shared status — same shorthand applied to each gene.

**Exit gate:** L3-shorthand unit tests; corpus shows L3-shorthand mixing with L1/L2/L3-prose per the fraction flags.

**Risk:** low — parallels L2 shorthand, schema already done.

---

## Sub-phase 3.4 — L4 heterogeneous compound (prose)

**Why:** same list structure as L3 but each gene carries independent status. This is where extraction models really start breaking.

**Scope:**
- `_L4_FRAMES`: `{gene1} {status1}, {gene2} {status2}, and {gene3} {status3}.` and semicolon-separated variants.
- Sampler: per-gene prevalence conditioning (one status per gene, independently drawn — reuses the existing per-biomarker distribution).
- Spans: N gene spans + N status spans. Each `AssertionFact` gets its own gene + status span.

**Exit gate:** L4 unit tests — every gene paired with its own status span, mixed-status record never conflates statuses. Regen corpus with L4 fraction.

**Risk:** medium — more spans to track, more room for off-by-one in separator positions.

---

## Sub-phase 3.5 — L4 shorthand / tabular

Parallel to 3.3 — shorthand form of heterogeneous multi-gene. `EGFR +, ALK -, ROS1 WT` style.

**Risk:** low.

---

## Sub-phase 3.6 — L5 negation scope

**Why:** a qualitatively different primitive — scope markers ("no", "absence of") apply to multiple entities; exception markers ("but", "except") flip polarity mid-sentence.

**Scope:**
- New frame family: `No {gene_list} mutations but {gene} {variant} detected.` / `Absence of {gene1}, {gene2} alterations; {gene3} {status}.`
- `AssertionFact` gains a `polarity_scope: Literal["direct", "negation_wide", "exception"]` field.
- Renderer tracks which genes fall inside the scope vs. after the exception.
- Status per assertion is computed from scope × base-status (negation-wide flips all to negative; exception segment carries its own status).

**Exit gate:** L5 tests verifying polarity is correctly attached per assertion across scope boundaries.

**Risk:** high — polarity bookkeeping is the trickiest part of Phase 3. Expect 1–2 iterations.

---

## Sub-phase 3.7 — L6 temporal / certainty qualification

**Scope:**
- Optional `temporal: str | None` and `certainty: str | None` on each `AssertionFact`.
- Frame family: `Previously {status} for {gene}; current testing {certainty} {status}.`
- Certainty vocabulary: suspected, confirmed, probable, rule-out, pending.
- Temporal vocabulary: previously, currently, at diagnosis, at relapse, post-TKI.

**Exit gate:** L6 tests — qualifiers attached per assertion; confirmation-by-inspection sample.

**Risk:** medium — new vocab categories + optional qualifier slots in record schema.

---

## Sub-phase 3.8 — L7 cross-sentence coreference

**Scope:**
- Record becomes multi-sentence: `sentence` → `sentences: tuple[str, ...]`. Each `AssertionFact` carries a `sentence_index: int` into the tuple.
- Discourse frames: setup → claim → qualifier/distractor (3-sentence typical).
- Coreference devices: pronoun resolution ("Results showed…"), nominal anaphora ("the mutation", "this variant").

**Exit gate:** L7 tests ensure every assertion's spans are in the sentence referenced by its `sentence_index`; coreference frames exercise both pronoun and nominal forms.

**Risk:** high — multi-sentence serialization is a structural shift; rendering bugs get harder to diagnose across sentence boundaries.

---

## Sub-phase 3.9 — Layer 5 noise expansion + Bug 3a fix

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
1. Update `docs/PHASE3_PLAN.md` with any scope-shift decisions made during implementation.
2. TDD: write failing tests for the new behavior.
3. Implement to green.
4. Run full pytest + generate a sanity-size corpus (5K records) and inspect.
5. Run the full 50K corpus + prevalence report before advancing to the next sub-phase.
6. Commit with a clear `feat(phase3.N):` or `refactor(phase3.N):` prefix.

Update the top-of-spec status line as each sub-phase ships.
