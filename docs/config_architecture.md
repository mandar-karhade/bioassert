# Probability-Weighted Variation Config — Architecture Guide (v1.1)

**Purpose:** Design of `common_variations.json` and `biomarkers.json`, how they drive the generator, and how to extend them.

**Files:**
- `common_variations.json` — surface variations shared across all biomarkers
- `biomarkers.json` — per-biomarker surface variations

**Scope:** The config ships the full Tier 1 panel (mutation, fusion, composite, expression biomarkers). Status distributions are flat per biomarker for a single cohort described by `biomarkers.json` — to target a different cohort, author a separate config and mix outputs externally.

---

## Changelog

**v1.2** (flat-prevalence simplification):
- Dropped `PatientProfile` and the multi-axis population-key cascade — patient demographics don't affect how an NER/assertion extraction model reads text.
- Renamed `status_distribution_by_population` → `status_distribution` (flat per biomarker).
- Removed `bioassert/generator/patient_sampler.py`; cohort mixing is now an external concern (separate `biomarkers.json` files mixed by a downstream heuristic).

**v1.1** (post-review fixes):
- Unified `prevalence_within_biomarker` field name (removed `_given_positive` suffix)
- Pulled degenerate expression-biomarker "negative" variants into a separate `negative_forms` block
- Extracted PD-L1 IHC clones into independent `clone_attribution` dimension
- Documented fusion `unspecified_fusion` coupling via `render_constraints.require_biomarker_name_forms` (Section 8.3)
- Documented placeholder substitution (`{gene}`, `{method}`, `{result}`, `{value}`) as first-class feature (Section 7.5)
- Added `$schema_type` discriminator to `technical_noise`; documented loader dispatch (Section 7.6)
- Documented `preferred_methods` → `common.test_methods` realization join (Section 2.4)

---

## 1. Core Design Principles

**1.1 Separation of concerns**
- Each config entry controls ONE dimension of surface variation
- The grammar composes dimensions; the config supplies weighted vocabulary
- Adding a new biomarker = adding a JSON entry (no code change)
- Adding a new sentence frame = changing a grammar template (no config change)

**1.2 Probabilities, not templates**
- Every variation category is a probability distribution over surface forms
- Distributions sum to 1.0 within a 0.001 tolerance, validated at load time
- Sampling is `numpy.random.choice(keys, p=weights)` — deterministic given a seed

**1.3 Two-level structure: variation ID → realization string**
- The `variations` dict maps variation_id → probability
- The `realizations` dict maps variation_id → actual surface string
- Separation lets you tune weights without touching strings, and vice versa

**1.4 Prevalence is flat per biomarker for the cohort**
- `biomarkers.json` describes a single cohort (e.g., NSCLC adenocarcinoma)
- Status prior is positive / negative / equivocal / not_tested per biomarker
- Sub-variant prevalence is conditional on the biomarker having a meaningful result to report
- Multi-cohort generation is handled externally by mixing corpora produced from separate cohort configs

**1.5 Schema-type discriminated objects**
- Nested dicts declare their schema shape via `$schema_type`
- Loader dispatches validation based on the declared type
- Unknown or missing type on a nested dict is an error

---

## 2. File Structure

### 2.1 common_variations.json

Top-level keys are variation categories. Each category declares its schema type.

| Category | Schema type | What it controls |
|----------|------------|------------------|
| `assertion_verbs` | weighted_variations | Verbs linking biomarker to status |
| `negation_phrases` | weighted_variations | Ways of expressing negative |
| `positive_phrases` | weighted_variations | Ways of expressing positive |
| `equivocal_phrases` | weighted_variations | "equivocal", "indeterminate", "failed QC" |
| `not_tested_phrases` | weighted_variations | "not tested", "pending", "QNS" |
| `temporality_modifiers` | weighted_variations | "previously", "at diagnosis" |
| `certainty_modifiers` | weighted_variations | "suspected", "confirmed", "rule out" |
| `test_methods` | weighted_variations | NGS, IHC, FISH, PCR, etc. — canonical realizations source |
| `method_attribution_frames` | weighted_variations | How methods link to results (contains placeholders) |
| `list_connectors` | weighted_variations | ", and", ";", ", " for compound sentences |
| `technical_noise` | post_process_transformations | Whitespace, case, hyphenation, punctuation noise |

**Standard `weighted_variations` structure:**
```json
{
  "$schema_type": "weighted_variations",
  "description": "...",
  "variations": {"form_id_1": 0.x, "form_id_2": 0.y},
  "realizations": {"form_id_1": "actual string", "form_id_2": "actual string"}
}
```

**`post_process_transformations` structure** (see Section 7.6):
```json
{
  "$schema_type": "post_process_transformations",
  "description": "...",
  "<sub_category>": {
    "description": "...",
    "distribution": {"form_id_1": 0.x, "form_id_2": 0.y}
  }
}
```

### 2.2 biomarkers.json

Top-level keys are biomarker names. Each biomarker entry follows the `_schema_template` documented at the top of the file:

```json
{
  "canonical_name": "full formal name",
  "alteration_type": "mutation | fusion | amplification | expression | composite",
  "gene_location": "chromosomal band",

  "name_forms": {"variations": {...}, "realizations": {...}},
  "variants": {"<variant_id>": {...}, ...},
  "status_distribution": {"positive": x, "negative": y, "equivocal": z, "not_tested": w},
  "preferred_methods": {"variations": {...}},

  "negative_forms": { ... },        // ONLY for expression biomarkers
  "clone_attribution": { ... }      // ONLY for IHC-based biomarkers
}
```

**Variant structure:**
```json
"<variant_id>": {
  "actionable": true,
  "prevalence_within_biomarker": 0.46,
  "category": "classical_activating",
  "notes": "optional clinical context",
  "name_forms": {"variations": {...}, "realizations": {...}},
  "measurement_range": [min, max],                      // optional, for {value} placeholders
  "render_constraints": {                               // optional, for coupled sampling
    "require_biomarker_name_forms": ["form_id_1", "form_id_2"]
  }
}
```

**`prevalence_within_biomarker` semantics** (unified field name):
- For mutation/fusion biomarkers: conditional probability given status=positive
- For expression biomarkers: conditional probability given status is meaningful-positive (not negative). Negative renderings live in `negative_forms`.

### 2.3 Flat status distribution

Each biomarker carries a single `status_distribution` over the four status labels:

```json
"status_distribution": {
  "positive": 0.20,
  "negative": 0.78,
  "equivocal": 0.01,
  "not_tested": 0.01
}
```

Weights sum to 1.0 ± 0.001; values are bounded in `[0, 1]`. The distribution is the prior the sampler draws from every time it picks a status for this biomarker — no patient-profile conditioning, no cascade.

**Cohort scope.** One `biomarkers.json` = one cohort. The numbers shipped in the repo target lung adenocarcinoma (NSCLC). To target a different cohort (e.g., squamous, colorectal), author a new `biomarkers.json` with cohort-appropriate rates and mix the generated corpora externally at the ratio you want.

**Why flat.** Demographic stratification (histology × ethnicity × smoking × age × sex) was removed in v1.2. The NER/assertion extraction task is invariant to who the patient is — labels are correct by construction regardless of the prevalence used to draw the status. Demographic cascade added configuration complexity without serving the task.

### 2.4 `preferred_methods` — weights-only override

`preferred_methods` contains only a `variations` dict — NO `realizations` dict. Realizations are looked up in `common.test_methods.realizations` at render time.

Validation requirements:
- Every key in `preferred_methods.variations` MUST exist in `common.test_methods.variations`
- Weights in `preferred_methods.variations` MUST sum to 1.0 ± 0.001
- Loader fails loudly if a biomarker references an unknown method key

Rationale: method realizations are canonical and domain-wide. A biomarker shouldn't be able to say "NGS" renders as anything other than "NGS." Only the prior weight on picking NGS vs IHC vs FISH is biomarker-specific.

---

## 3. How Generation Uses the Config

Generation flow for one assertion:

```
1. For each biomarker in the panel:
   1a. Sample status from biomarker.status_distribution (Section 2.3):
       positive | negative | equivocal | not_tested
   1b. Branch on alteration_type and status:
       - expression + negative → sample from negative_forms
       - expression + positive → sample variant, then render with {value} if present
       - mutation/fusion + positive → sample variant from variants
       - any + equivocal → sample from common.equivocal_phrases
       - any + not_tested → sample from common.not_tested_phrases
   1c. If variant has render_constraints.require_biomarker_name_forms:
       - Restrict biomarker name_form sampling to that subset (renormalize weights)
     Else:
       - Sample biomarker name_form from full distribution
   1d. Sample assertion verb from common.assertion_verbs
   1e. Sample test method: if biomarker.preferred_methods exists, use it; 
       else common.test_methods. Look up realization in common.test_methods.realizations
   1f. Sample temporality modifier from common.temporality_modifiers
   1g. Sample certainty modifier from common.certainty_modifiers
   1h. For PD-L1 (or any IHC biomarker with clone_attribution):
       - Roll Bernoulli with clone_attribution.attachment_probability
       - If true, sample clone from clone_attribution.variations
     ↓
2. Pick a sentence frame from grammar templates matching complexity level
     ↓
3. Render frame with sampled slot values (Section 7.5 placeholder substitution)
     ↓
4. Apply technical_noise transformations (whitespace, case, hyphenation) with their sampled probabilities
     ↓
5. Emit (rendered_sentence, structured_assertion) pair with character spans
```

---

## 4. Extending the Config

### 4.1 Adding a new biomarker

1. Add entry to `biomarkers.json` with all required fields
2. Define `name_forms` (at least 3 variations summing to 1.0)
3. Define `variants` (at least one entry; can be `unspecified_*` if no real sub-variants)
4. Define `status_distribution` with positive/negative/equivocal/not_tested weights summing to 1.0
5. Define `preferred_methods` if different from common defaults (optional)
6. Add `negative_forms` if alteration_type is `expression`
7. Add `clone_attribution` if IHC-based and clone reporting is common in real data
8. Run the validator (Section 7) — catches sum errors, key mismatches, missing fallbacks

### 4.2 Adding a new sub-variant

1. Add to `biomarkers.<biomarker>.variants.<new_variant_id>`
2. Adjust other variants' `prevalence_within_biomarker` to keep sum = 1.0
3. Define `name_forms` for the new variant

### 4.3 Adding a new surface form to an existing category

1. Add to `variations` with a probability
2. Add matching `realizations` entry
3. Reduce one or more existing probabilities to keep sum = 1.0

### 4.4 Adding a new common variation category

1. Add top-level key to `common_variations.json` with `$schema_type: "weighted_variations"`
2. Reference it in the grammar templates as `{new_category}`
3. Update loader validation schema if the new category introduces any new placeholder syntax

---

## 5. Probability Calibration

**Weights in v1.1 are best-estimate priors, not measured frequencies.** They reflect:
- Published NSCLC molecular epidemiology (for prevalence layers)
- Clinical writing conventions inferred from pathology report samples
- Explicit design choices to ensure enough diversity

**Calibration process for v2:**
1. Collect 200–500 real pathology report sentences (from collaborator or n2c2 after DUA)
2. Annotate each for: biomarker name form, status phrase form, test method, etc.
3. Compute empirical frequencies per category
4. Update config weights to match empirical distribution (with some smoothing)
5. Re-run Phase 7 transfer evaluation

Weights are **not** meant to be exact. A 0.25 vs 0.30 difference is within calibration noise. A 0.05 vs 0.50 difference should trigger review.

---

## 6. Multi-Cohort Generation (External)

`biomarkers.json` is single-cohort. To target a different cohort (e.g., squamous NSCLC, colorectal, breast), author a separate `biomarkers.json` with that cohort's status priors. Mix generated corpora externally:

```bash
# Example: 70% adenocarcinoma, 30% squamous
python scripts/generate_corpus_v1.py \
    --config bioconfigs/biomarkers_adenocarcinoma.json --n 7000 --seed 1 \
    --output-dir /tmp/adeno
python scripts/generate_corpus_v1.py \
    --config bioconfigs/biomarkers_squamous.json --n 3000 --seed 2 \
    --output-dir /tmp/squamous
cat /tmp/adeno/corpus.jsonl /tmp/squamous/corpus.jsonl > mixed.jsonl
shuf mixed.jsonl > mixed_shuffled.jsonl
```

The mixing ratio is a heuristic choice, not a schema concern.

---

## 7. Validation (Required in Loader)

The Pydantic loader must validate at load time:

1. **Every `variations` dict has a matching `realizations` dict with identical keys**
2. **Every `variations` dict sums to 1.0 ± 0.001**
3. **Every biomarker's `variants` prevalences sum to 1.0 ± 0.001**
4. **Every `status_distribution` sums to 1.0 ± 0.001**
5. **Every `preferred_methods.variations` key exists in `common.test_methods.variations`**
6. **Every `render_constraints.require_biomarker_name_forms` entry is a valid biomarker-level name_form key**
7. **Every placeholder `{name}` in realization strings is in the supported table (Section 7.5)**
8. **Every variant containing `{value}` in any realization has a `measurement_range`**
9. **No duplicate variation_ids within a category**
10. **All prevalence values are in [0.0, 1.0]**
11. **`$schema_type` dispatch succeeds for every nested typed object**

Failure at load time (not generation time). Generation is deterministic only if configs are consistent.

### 7.5 Placeholder substitution

Realization strings may contain format-string placeholders filled at render time. The loader must enforce the supported set.

**Supported placeholders:**

| Placeholder | Filled with | Sampled when | Used in |
|-------------|------------|--------------|---------|
| `{gene}` | biomarker canonical or sampled name form realization | at render time | `common.negation_phrases.no_X_identified` |
| `{method}` | sampled test method realization | when method != unspecified | `common.method_attribution_frames.*` |
| `{result}` | sampled status phrase realization | at render time | `common.method_attribution_frames.METHOD_colon_result` |
| `{value}` | sampled numeric value from variant's `measurement_range` | at render time for expression biomarkers | PD-L1 `TPS_*`, CPS, TMB variants + `negative_forms` with `measurement_range_for_value_placeholder` |
| `{variant}` | reserved | not yet used | reserved |

**Rules:**
- Placeholders use single curly braces `{name}`
- Literal `{` or `}` in realization strings is an error (loader validates via regex scan)
- Renderer uses `str.format_map(context_dict)` at render time
- Any unfilled placeholder raises — loud failure, no silent empty substitution
- Loader validator warns if a placeholder appears in a realization string but the parent object lacks the machinery to fill it (e.g., `{value}` in a variant with no `measurement_range`)

### 7.6 Schema-type dispatch

Nested dicts declare their schema via `$schema_type`. The validator dispatches:

| `$schema_type` value | Validation rules |
|---------------------|------------------|
| `weighted_variations` | Must have `variations` + `realizations` with identical keys; weights sum to 1.0 |
| `weighted_variations_with_attachment` | As above, plus `attachment_probability` in [0, 1] |
| `post_process_transformations` | Each sub-key (except `description`, `$schema_type`) must be a dict with `distribution` summing to 1.0 |

**Objects without `$schema_type` are treated as `weighted_variations` (default)** for backward compatibility with simple cases — but the loader emits a deprecation warning recommending explicit annotation.

---

## 8. Special Cases

### 8.1 Expression biomarkers (PD-L1, TMB)

Expression biomarkers differ structurally from mutations:
- `variants` contains only **positive-direction** categories (TPS_high, TPS_intermediate, CPS; TMB_high)
- Negative renderings live in a separate `negative_forms` block (same `weighted_variations` schema)
- When status=negative is sampled, the generator renders from `negative_forms` instead of picking a variant
- `{value}` placeholders work in both `variants` (with variant-level `measurement_range`) and `negative_forms` (with `measurement_range_for_value_placeholder` at the block level)

### 8.2 Composite biomarkers (ERBB2, MET)

These have multiple fundamentally different alteration categories (mutation vs amplification vs overexpression). Each is encoded as a separate variant with its own `category` and `actionable` field. No special schema handling required — the variant list naturally accommodates this.

### 8.3 Fusion biomarkers and `unspecified_fusion` coupling

ALK, ROS1, RET, NTRK all have `unspecified_fusion` variants representing cases where the fusion partner isn't reported. When sampled, these variants produce an empty or near-empty string. The biomarker name form MUST then carry the "fusion"/"rearrangement" signal, or the sentence loses meaning.

**Mechanism:** `render_constraints.require_biomarker_name_forms` on the variant. When the sampler picks a variant with this field, it restricts biomarker-level name_form sampling to the listed subset (weights renormalized within the subset).

**Example (ALK):**
```json
"unspecified_fusion": {
  "name_forms": {"variations": {"unspecified_fusion_generic": 0.60, ...}, ...},
  "render_constraints": {
    "require_biomarker_name_forms": [
      "ALK_rearrangement",
      "ALK_fusion",
      "ALK_gene_rearrangement",
      "ALK_translocation"
    ]
  }
}
```

When `unspecified_fusion` is sampled, biomarker name MUST come from that 4-element subset — never "ALK" alone, which would lose the fusion signal.

### 8.4 IHC clone attribution (PD-L1)

Clones attach *to* TPS/CPS values. They are not mutually exclusive with variants — they're an independent dimension.

**Mechanism:** `clone_attribution` block at biomarker level (separate from `variants`).

```json
"clone_attribution": {
  "$schema_type": "weighted_variations_with_attachment",
  "attachment_probability": 0.40,
  "variations": {"clone_22C3": 0.45, ...},
  "realizations": {"clone_22C3": "(Clone 22C3)", ...}
}
```

**Sampler behavior:**
1. Sample TPS/CPS variant normally
2. Roll Bernoulli with `attachment_probability`
3. If true, sample clone from `clone_attribution.variations` and append to rendered sentence
4. If false, no clone rendered

Generalizes to any future IHC-based biomarker (HER2 IHC, etc.) — add `clone_attribution` to their entry, done.

### 8.5 Resistance mutations

EGFR T790M, C797S, ALK G1202R, ROS1 G2032R, etc. are encoded as variants with `actionable: false` and `category: resistance`. They have low base prevalence but appear frequently in rebiopsy reports after TKI progression. Document-level composition (Phase 4) should handle temporal context — this config layer doesn't.

---

## 9. Why This Architecture is Methodologically Defensible

Reviewers of a clinical NLP paper will look for:

- **Deterministic generation** — same seed → same output. ✓ (seeded stdlib random)
- **Clinically plausible prevalence** — not uniform random. ✓ (per-biomarker cohort priors)
- **Controlled diversity** — no single sentence frame dominates. ✓ (10+ variations per category)
- **Traceability** — every label ties back to a config entry. ✓ (structured assertions reference variation IDs)
- **Extensibility** — reproducible by other groups. ✓ (JSON is language-agnostic)
- **Version control** — weights can be tuned and diffed. ✓ (plain text JSON with changelog)

The alternative — LLM-generated synthetic data — fails on determinism, traceability, and (usually) clinical plausibility. This architecture wins on all three.

---

## 10. Module Layout

- `bioassert/config/schema.py` — Pydantic models for common + biomarkers configs (with `$schema_type` discriminated unions)
- `bioassert/config/loader.py` — JSON loading, model instantiation
- `bioassert/config/validator.py` — runs cross-config validations from Section 7 at load time
- `bioassert/generator/sampler.py` — weighted sampling with seeded `random.Random`; handles `render_constraints`, `clone_attribution` Bernoulli
- `bioassert/generator/renderer.py` — composes sampled values into final sentences; handles placeholder substitution (Section 7.5); L1–L7 complexity shapes
- `bioassert/generator/post_process.py` — applies `technical_noise` transformations
- Unit tests (`tests/`):
  - Config loads and validates
  - All probability distributions sum to 1.0
  - `render_constraints` actually restricts name_form sampling
  - `clone_attribution` Bernoulli attachment works as expected
  - Placeholder substitution fills all placeholders and raises on unfilled
  - Generated sentences contain labeled spans at claimed positions
  - Same seed → same output
  - 10,000 generations hit expected prevalence distribution within 2σ

---

## 11. What This Config Does NOT Cover (Yet)

Deferred to later phases or intentionally excluded:

- **Document-level composition** (Phase 4) — multi-paragraph reports with distractor content
- **Cross-sentence coreference** (L7 complexity) — "Patient underwent NGS. Results showed..."
- **Temporal evolution** — biomarker status changing across multiple reports for the same patient
- **OCR/PDF artifacts** — currently encoded in `technical_noise.punctuation_variation.ocr_artifact` but minimal
- **Pathologist writing idiosyncrasies** — requires real-data calibration (Phase 7)
- **Non-English terminology** — US-English-only in v1
- **Squamous-specific markers** (FGFR1, PIK3CA) — Tier 3 priority
- **Co-mutations** (TP53, STK11, KEAP1) — Tier 2, Phase 2b
- **Rare / research biomarkers** (NRG1, AXL, etc.) — Tier 2–3

These are listed in `docs/limitations.md`.

---

**The config files are the source of truth. The grammar is the combinator. The loader enforces every invariant. Together they produce deterministic, traceable, clinically plausible synthetic data at any scale — without a single LLM API call.**
