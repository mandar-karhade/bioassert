# Variation & Noise Profile

**Scope:** canonical reference for every layer of surface variation and technical noise currently shipped in `bioassert`. Lists every frame family, noise transform, vocabulary category, and distribution weight actually present in `main`, with pointers to the source.

**Audience:** anyone consuming the corpus, reviewing a PR that changes variation/noise, or writing the dataset card.

**Non-goal:** an execution log of how variation/noise was built over time. Phase/sub-phase numbering is intentionally not tracked here — it was creating confusion with the layer numbering (Phase 3.1 vs Layer 1 etc.) and has been retired. `git log` is authoritative for history.

---

## 1. Architecture recap

The generator is structured along the seven-layer spec in [project_spec.md §3](project_spec.md#3-seven-layer-architecture). For variation/noise purposes, four of the seven layers actively shape the surface text that lands in the corpus:

| Spec Layer | Role in variation | Mechanism | Source |
|---|---|---|---|
| L3 Sentence grammar | Frame templates (what structure the sentence has) | `_L?_FRAMES`, `_L?_SHORTHAND_FRAMES`, `_L5_FRAMES`, `_L6_*`, `_L7_*` | [`bioassert/generator/renderer.py`](../bioassert/generator/renderer.py) |
| L4 Complexity stratification | 7 complexity levels (L1–L7), orthogonal compounding tier (low/high) | `render_l{1..7}_record` entry points | [`bioassert/generator/renderer.py`](../bioassert/generator/renderer.py) |
| L5 Technical/surface noise | 7 post-process transforms applied after render | `apply_technical_noise` | [`bioassert/generator/post_process.py`](../bioassert/generator/post_process.py) |
| L6 Clinical realism | Status-prevalence sampling against each biomarker's flat per-cohort distribution | `sample_status`, `StatusDistribution` | [`bioassert/generator/sampler.py`](../bioassert/generator/sampler.py) |

Layer 7 (document-level composition — multi-paragraph pathology reports with distractor sections) is **deferred**. See [§7](#7-deferred-spec-layer-7-document-composition).

---

## 2. Layer 4 — complexity levels

Every corpus record carries a `complexity_level ∈ {L1, L2, L3, L3S, L4, L4S, L5, L6, L7}` and a `compounding_tier ∈ {"low", "high"}`. Level selection is a user-controlled fraction at corpus-gen time (`--l2-fraction`, `--l3-fraction`, etc.); tier is sampled per compound record from `--compound-low` / `--compound-high` (default 0.5/0.5, must sum to 1.0).

### 2.1 L1 — single assertion, canonical

**Frames:** 15 entries in `_L1_FRAMES` ([renderer.py:115](../bioassert/generator/renderer.py#L115)). Each frame has a `template` with `{gene}`, `{verb}`, `{status}`, optional `{method}`, `{clone}`, `{measurement}` placeholders.

**Representative templates:**
```
{gene} {verb} {status}.
{status} for {gene}.
{gene} mutation {status}.
{gene} testing: {status}.
Result: {gene} {status}.
```

**Assertions:** 1 fact per record. Spans: `{gene, status}` (+ `method`, `clone`, `measurement`, `variant` when slotted).

**Variation sources:**
- Gene name form — weighted draw from [`biomarkers.json`](../projects/nsclc_adenocarcinoma/configs/biomarkers.json) `name_forms.variations` (per-gene, see [§4](#4-per-gene-alias-vocabulary-name_forms)).
- Verb — weighted from `assertion_verbs` (14 entries, [§3.1](#31-vocabulary-categories)).
- Status surface — weighted from the matching `{positive,negation,equivocal,not_tested}_phrases` (6–12 entries each).
- Method attachment — Bernoulli per biomarker config (default ~30% for mutation biomarkers).
- Measurement — only for expression biomarkers (PD-L1 TPS %, TMB mut/Mb).

### 2.2 L2 — lexical shorthand / tabular single-fact

**Frames:** 8 entries in `_L2_FRAMES` ([renderer.py:133](../bioassert/generator/renderer.py#L133)). Compressed, punctuation-heavy.

**Representative templates:**
```
{gene} {shorthand_status}
{gene}: {shorthand_status}
{gene} [{shorthand_status}]
{gene}  {shorthand_status}       (tabular, multi-space)
```

**Status vocabulary:** drawn from `positive_shorthand` / `negation_shorthand` (7 entries each) — includes `+`, `-`, `(+)`, `(-)`, `pos`, `neg`, `WT`, `mut`.

**Assertions:** 1 fact. Same span model as L1.

### 2.3 L3 — list distribution, shared status (prose)

**Frames:** 7 entries in `_L3_FRAMES` ([renderer.py:144](../bioassert/generator/renderer.py#L144)).

**Representative templates:**
```
{gene_list} were {status}.
{gene_list} showed {status}.
No {gene_list} alterations detected.
{gene_list}: {status}.
```

**Assertions:** 2–8 facts (N genes share one status). Compounding tier: low=2 genes; high=3–8 (clamped to pool size). Cross-class pool mixing (TMB/PD-L1 mixed with mutation biomarkers) is allowed in the high tier.

**List coordinator:** comma + Oxford-optional + final `and`/`or`. `list_connectors` category (6 entries) controls separator style.

**Spans:** N gene spans + 1 shared status span. Every `AssertionFact.spans` points at its own gene span and the shared status span.

### 2.4 L3S — L3 shorthand / tabular

**Frames:** 9 entries in `_L3_SHORTHAND_FRAMES` ([renderer.py:158](../bioassert/generator/renderer.py#L158)).

**Representative templates:**
```
{gene1} -, {gene2} -, {gene3} -
{gene1}: WT; {gene2}: WT; {gene3}: WT
{gene1}\t{status} / {gene2}\t{status} / ...
```

Tab-separated and newline-separated variants present. Same shared-status invariant as L3.

### 2.5 L4 — heterogeneous compound (prose)

**Frames:** 21 entries in `_L4_FRAMES` ([renderer.py:174](../bioassert/generator/renderer.py#L174)). Includes comma, semicolon, and period separator variants.

**Representative templates:**
```
{gene1} {status1}, {gene2} {status2}, and {gene3} {status3}.
{gene1} was {status1}; {gene2} was {status2}.
{gene1} {status1}. {gene2} {status2}. {gene3} {status3}.
```

**Assertions:** 2–8 facts, each with independently sampled status drawn from that gene's flat `status_distribution` in `biomarkers.json`. Each fact gets its own gene span + own status span.

### 2.6 L4S — L4 shorthand / tabular

**Frames:** 9 entries in `_L4_SHORTHAND_FRAMES` ([renderer.py:188](../bioassert/generator/renderer.py#L188)).

**Representative templates:**
```
{gene1} +, {gene2} -, {gene3} WT
{gene1}: pos\t{gene2}: neg\t{gene3}: equivocal
```

### 2.7 L5 — negation scope

**Frames:** 46 entries in `_L5_FRAMES` ([renderer.py:210](../bioassert/generator/renderer.py#L210)) covering three polarity shapes.

**Representative templates:**
```
No evidence of {gene_list} alterations.                       (negation_wide)
Absence of {gene_list} mutations.                             (negation_wide)
No {gene_list} mutations but {gene} {variant} detected.       (exception)
No biomarker on the panel was positive other than {ex_gene}.  (panel-wide exception)
{gene_list} wild-type except {gene} {status}.                 (exception)
```

**Polarity scope field:** each `AssertionFact` carries `polarity_scope ∈ {"direct", "negation_wide", "exception"}`. Status is derived from scope × base-status: negation-wide flips all listed genes to `negative`; exception segment carries its own status.

**Panel-wide variants:** some L5 frames assert over the "panel" as a whole and surface a single labeled fact (the exception). Clinical register — the user-confirmed decision on 2026-04-23 is that clinical reports use `and`/`except`/`other than`, not `but`/`however`; `None` dominates the choice when stating the general case. See [`docs/known_issues.md`](known_issues.md) and the `feedback_l5_clinical_exception_markers.md` memory.

### 2.8 L6 — temporal / certainty qualification

**Shape dispatcher:** `render_l6_record` picks one of three shapes uniformly:

| Shape | Frame family | Facts | Qualifier fields populated |
|---|---|---|---|
| `temporal` | 3 frames in `_l6_temporal_frame{1,2,3}` | 2 | both `temporal` set, both `certainty=None` |
| `certainty` | Extended L1 frame with `{certainty}` slot | 1 | `certainty` set, `temporal=None` |
| `combined` | 2-fact frame with both qualifiers | 2 | both `temporal` and `certainty` set |

**Temporal vocabulary:** 5 entries in `L6_TEMPORAL_VOCAB` — `at diagnosis`, `at relapse`, `previously`, `currently`, `post-TKI`.

**Certainty vocabulary:** 7 entries in `certainty_modifiers` — `confirmed`, `suspected`, `probable`, `rule-out`, `pending`, plus hedged phrases.

**Temporal frame variants:**
- `temporal_frame1`: `{gene} was {status_a} {temporal_a}, {status_b} {temporal_b}.` (shared gene span, 2 status facts)
- `temporal_frame2`: `{gene} was {status_a} {temporal_a}; {gene} was {status_b} {temporal_b}.` (**two distinct gene mentions** — this is the sole site that triggers abbreviation_inconsistency mixed-mode swaps in the current frame library)
- `temporal_frame3`: `{status_a} {temporal_a}, {status_b} {temporal_b} for {gene}.` (shared gene, status-first)

### 2.9 L7 — cross-sentence coreference

**Shapes:** 3 entries in `L7_SHAPES`.

| Shape | Structure | Fact sentence_index |
|---|---|---|
| `setup_claim` | `setup. claim.` | 1 |
| `claim_anaphora` | `claim. anaphor.` | 0 |
| `setup_claim_qualifier` | `setup. claim. qualifier.` | 1 |

**Template pools:**
- 5 setup templates in `_L7_SETUP_TEMPLATES`
- 2 claim templates in `_L7_CLAIM_TEMPLATES`
- 5 anaphor templates in `_L7_ANAPHOR_TEMPLATES`
- 5 qualifier templates in `_L7_QUALIFIER_TEMPLATES`
- 5 anaphor nouns in `L7_ANAPHOR_NOUNS` (`this alteration`, `the mutation`, `this finding`, `the variant`, `this result`)

**Record shape:** multi-sentence — `RenderedRecord.sentences: tuple[str, ...]` with per-fact `sentence_index: int` pointing into it. `sentence` is the joined form for backward compatibility.

**No repeat gene mentions by design.** Coreference is pronoun/nominal, not repeat-alias.

### 2.10 Compounding tier (orthogonal)

Applies to L3, L3S, L4, L4S, L5. Not L1/L2 (always single-gene) or L6/L7 (single-gene semantics).

| Tier | N genes | Cross-class pool mix |
|---|---|---|
| `low` | 2 | no |
| `high` | 3–8 (clamped to pool size) | yes (TMB/PD-L1 can appear alongside mutation biomarkers) |

Constants: `COMPOUND_LOW_N = 2`, `COMPOUND_HIGH_MIN = 3`, `COMPOUND_HIGH_MAX = 8` in `renderer.py`.

---

## 3. Layer 3 — shared vocabulary

[`projects/nsclc_adenocarcinoma/configs/common_variations.json`](../projects/nsclc_adenocarcinoma/configs/common_variations.json) holds the shared surface-variation vocabularies. All are `$schema_type: weighted_variations` — a `{variation_id: weight}` map plus a `{variation_id: realization_string}` map. Weights sum to 1.0 (validated at load time, tolerance 0.001).

### 3.1 Vocabulary categories

| Category | Entries | Purpose |
|---|---|---|
| `assertion_verbs` | 14 | Linking verbs: `was`, `is`, `showed`, `revealed`, `detected as`, `:`, `—`, etc. |
| `positive_phrases` | 12 | `positive`, `detected`, `mutated`, `present`, `(+)`, `pos`, `mutation present`, ... |
| `negation_phrases` | 12 | `negative`, `not detected`, `wild-type`, `WT`, `absent`, `(-)`, `no {gene} alteration identified`, ... |
| `equivocal_phrases` | 6 | `equivocal`, `indeterminate`, `inconclusive`, `failed QC`, ... |
| `not_tested_phrases` | 7 | `not tested`, `pending`, `awaited`, `QNS` (quantity not sufficient), ... |
| `positive_shorthand` | 7 | `+`, `(+)`, `pos`, `mut`, `mut+`, `positive`, `detected` |
| `negation_shorthand` | 7 | `-`, `(-)`, `neg`, `WT`, `wt`, `negative`, `not detected` |
| `temporality_modifiers` | 8 | Used in L1/L2 optional qualifier slots (independent of L6 temporal qualifier) |
| `certainty_modifiers` | 7 | `suspected`, `probable`, `confirmed`, `rule-out`, `pending`, ... |
| `test_methods` | 11 | `NGS`, `PCR`, `FISH`, `IHC`, `Sanger`, `ddPCR`, `RT-PCR`, ... |
| `method_attribution_frames` | 7 | `by {method}`, `on {method}`, `via {method}`, ... |
| `list_connectors` | 6 | Comma-and, comma-only, semicolon, Oxford-comma variants |

### 3.2 Status-distribution resolution

Status (`positive` / `negative` / `equivocal` / `not_tested`) is **not** drawn from the vocab categories directly — it is drawn from the biomarker's flat `status_distribution` in `biomarkers.json` (see [§5](#5-layer-6--prevalence-calibration)), then **surfaced** via the matching phrase category.

---

## 4. Per-gene alias vocabulary (`name_forms`)

[`projects/nsclc_adenocarcinoma/configs/biomarkers.json`](../projects/nsclc_adenocarcinoma/configs/biomarkers.json) holds per-gene `name_forms` — the surface-string pool used by L1/L2/L3/... renderers when emitting `{gene}`. Also the swap pool for the `abbreviation_inconsistency` mixed-mode transform.

Tier 1 panel coverage (11 biomarkers):

| Gene | Alias pool size | Sample realizations |
|---|---|---|
| EGFR | 8 | `EGFR`, `EGFR mutation`, `Epidermal Growth Factor Receptor`, `EGF-R`, `HER1`, `ERBB1`, `EGFR gene`, `egfr` |
| KRAS | 6 | `KRAS`, `KRAS mutation`, `K-RAS`, `KRAS gene`, `K-ras`, `kras` |
| ALK | 7 | `ALK`, `ALK rearrangement`, `ALK fusion`, `ALK-positive`, `ALK gene rearrangement`, `ALK translocation`, `alk` |
| ROS1 | 6 | `ROS1`, `ROS1 rearrangement`, `ROS1 fusion`, `ROS-1`, `ROS1-positive`, `ros1` |
| BRAF | 6 | `BRAF`, `BRAF mutation`, `B-RAF`, `BRAF gene`, `b-raf`, `braf` |
| MET | 9 | `MET`, `MET alteration`, `c-MET`, `cMET`, `c-Met`, `MET gene`, `HGFR`, `met`, ... |
| RET | 6 | `RET`, `RET fusion`, `RET rearrangement`, `RET-positive`, `RET gene`, `ret` |
| ERBB2 | 8 | `ERBB2`, `HER2`, `HER-2`, `ERBB2 (HER2)`, `HER2/neu`, `ERBB2 gene`, ... |
| NTRK | 9 | `NTRK`, `NTRK fusion`, `NTRK1`, `NTRK3`, `TRK`, `pan-NTRK`, `NTRK2`, ... |
| PD-L1 | 8 | `PD-L1`, `PDL1`, `PD L1`, `PD-L1 expression`, `CD274`, `programmed death ligand 1`, ... |
| TMB | 5 | `TMB`, `tumor mutational burden`, `mutational burden`, `mutational load`, `tmb` |

**Status distribution:** each biomarker ships a single flat `status_distribution` (positive / negative / equivocal / not_tested, sums to 1.0) representing the target cohort. Multi-cohort corpora are produced by authoring separate `biomarkers.json` files and mixing generated output externally. See [§5](#5-layer-6--prevalence-calibration).

---

## 5. Layer 6 — prevalence calibration

Prevalence sampling is a covariate on status selection only. It does **not** drive surface-text variation — see [project_spec.md Layer 6](project_spec.md#layer-6--probabilistic-clinical-realism).

### 5.1 Flat per-biomarker distribution

Each biomarker in [`projects/nsclc_adenocarcinoma/configs/biomarkers.json`](../projects/nsclc_adenocarcinoma/configs/biomarkers.json) carries a single `status_distribution` field:

```json
"status_distribution": {
  "positive": 0.20,
  "negative": 0.70,
  "equivocal": 0.02,
  "not_tested": 0.08
}
```

The four probabilities sum to 1.0 (validated at load time, tolerance 0.001). `sampler.sample_status(entry, rng)` draws from this distribution directly — no patient profile, no cascade, no fallback.

**Design rationale:** the NER / assertion-extraction task is invariant to who the patient is. A model reading "EGFR was positive" processes the same tokens whether the patient is a 72-year-old East-Asian non-smoker or a 45-year-old Western smoker. Demographic cascade was removed in v1.2 — patient-level factors add complexity without changing the supervised signal the benchmark measures. See [config_architecture.md §1.4](config_architecture.md).

### 5.2 Cohort scope and multi-cohort corpora

A single `biomarkers.json` describes **one cohort**. The shipped config at [`projects/nsclc_adenocarcinoma/configs/biomarkers.json`](../projects/nsclc_adenocarcinoma/configs/biomarkers.json) targets lung adenocarcinoma, so prevalence values are adenocarcinoma-weighted averages across the PIONEER / GENIE / LCMC source studies.

For squamous cohorts, RNA-fusion-driven panels, or pan-tumor panels, create a sibling project directory (e.g., `projects/nsclc_squamous/`) with its own `configs/biomarkers.json`, then mix generated corpora externally at the target ratio:

```bash
bioassert generate --project projects/nsclc_adenocarcinoma --n 8000 --seed 1 --tag adeno
bioassert generate --project projects/nsclc_squamous       --n 2000 --seed 2 --tag squamous
cat projects/nsclc_adenocarcinoma/outputs/run_*_adeno_*/corpus.jsonl \
    projects/nsclc_squamous/outputs/run_*_squamous_*/corpus.jsonl \
    | shuf > corpus.jsonl
```

### 5.3 Source priors

Each distribution is sourced from PIONEER, AACR GENIE v13, LCMC, and published meta-analyses. Citation list will live in `projects/nsclc_adenocarcinoma/references/prevalence_sources.bib` (future).

Adenocarcinoma-weighted priors (shipped config):
- EGFR: ~20% positive (cohort-weighted average across Western + East Asian)
- KRAS: ~25% positive
- ALK: ~5% positive
- ROS1: ~1–2% positive
- BRAF V600E: ~1–2% positive
- PD-L1 ≥50%: ~30% positive
- TMB high (≥10 mut/Mb): ~15–20%

---

## 6. Layer 5 — technical / surface noise

Applied by `apply_technical_noise` in [`post_process.py`](../bioassert/generator/post_process.py) after render. Each category samples a mode from its distribution and may modify the sentence + shift spans. Single-fact records run through the full stack; multi-fact records (L3+) skip the six typographic transforms (whitespace/case/hyphenation/punctuation/OCR/PDF) but still get `abbreviation_inconsistency`. Multi-sentence records (L7) skip abbreviation_inconsistency too (no repeat-gene mentions by design).

### 6.1 Application order

```
1. hyphenation_gene_names
2. whitespace
3. case_variation
4. punctuation_variation
5. ocr_corruption
6. pdf_artifact
7. abbreviation_inconsistency
```

The order matters: hyphenation runs before whitespace so slot-boundary detection sees the post-hyphenation sentence; abbreviation_inconsistency runs last so it operates on the final span positions.

### 6.2 Noise categories and modes

Every record's `post_process` dict reports the mode applied for every category. A `canonical` mode means the transform sampled its no-op branch. A `skipped` mode means the record is in a passthrough path (multi-fact typographic skip; L7 multi-sentence abbreviation skip).

#### hyphenation_gene_names

Injects hyphens / spaces / line-breaks into gene surface forms when the surface is all-alpha. Skips gene names containing digits/hyphens (`PD-L1`, `ROS1`, `BRAF V600E`) — they already have intrinsic punctuation.

| Mode | Weight | Effect |
|---|---|---|
| `canonical` | 0.80 | no change |
| `hyphenated` | 0.12 | `EGFR` → `EG-FR` (insert `-` at random interior position); gene span grows by 1 |
| `spaced` | 0.05 | `EGFR` → `EG FR` |
| `linebroken` | 0.03 | `EGFR` → `EG\nFR` |

#### whitespace

Substitutes whitespace at slot-boundary positions only (Bug 3a fix — prose-interior subs no longer occur). Targets the space between the gene/status/method slot surfaces.

| Mode | Weight | Effect |
|---|---|---|
| `single_space` | 0.85 | no change (canonical) |
| `double_space` | 0.05 | one slot-boundary space → `"  "`; spans shift by +1 for downstream slots |
| `tab` | 0.03 | space → `\t`; length-preserving |
| `newline_mid_sentence` | 0.04 | space → `\n`; length-preserving |
| `no_space_after_punct` | 0.03 | drop a space after punctuation; spans shift by -1 |

#### case_variation

Applied to the full sentence. Labeled spans remain valid (character positions unchanged); the surface casing may drift.

| Mode | Weight | Effect |
|---|---|---|
| `canonical` | 0.92 | no change |
| `all_lowercase` | 0.04 | `sentence.lower()` |
| `all_uppercase` | 0.02 | `sentence.upper()` |
| `title_case` | 0.02 | `sentence.title()` |

#### punctuation_variation

| Mode | Weight | Effect |
|---|---|---|
| `canonical` | 0.90 | no change |
| `missing_period` | 0.05 | drop trailing `.`; length -1 |
| `extra_comma` | 0.03 | insert `,` at a clause boundary; length +1 |
| `ocr_artifact` | 0.02 | single punctuation glyph mutated (`.` → `,` or `,` → `.`); length-preserving |

#### ocr_corruption

Length-preserving char swaps on characters **outside labeled spans**. Swap pairs: `l↔1`, `O↔0`, `S↔5`, `B↔8`, `Z↔2`, `I↔l`.

| Mode | Weight | Effect |
|---|---|---|
| `canonical` | 0.93 | no change |
| `light` | 0.05 | ~1 char swap outside all labeled spans |
| `moderate` | 0.02 | ~3 char swaps outside all labeled spans |

#### pdf_artifact

| Mode | Weight | Effect |
|---|---|---|
| `canonical` | 0.96 | no change |
| `hyphen_linebreak` | 0.04 | insert `-\n` inside a non-labeled ≥4-char alphabetic token; downstream spans shift by +2 |

#### abbreviation_inconsistency (sub-phase 3.11)

When a record contains ≥2 labeled mentions of the same canonical gene, re-render mention #2+ with an alternate alias from `biomarkers.{gene}.name_forms.variations` (weighted, excluding mention #1's surface). Mention #1 stays verbatim; spans shift by length delta for all downstream labels.

| Mode | Weight | Effect |
|---|---|---|
| `canonical` | 0.70 | no change |
| `mixed` | 0.30 | swap mention #2+ to a drawn alternate |

**Trigger sites in the current frame library:** L6 `temporal_frame2` is the sole site today (~28% of L6 records). Future Layer 4 frames that introduce repeat mentions will be covered automatically.

Measured exit gate on 50K v1_phase3.11 corpus: **372/395 = 94.2%** of L6 mixed-mode repeat-mention records have distinct surfaces post-swap (target ≥80%).

### 6.3 Corpus-level prevalence of noise

Independent per record (products of per-category weights). On a 50K corpus, ~8% of records end up with **any** non-canonical typographic transform applied (complement of `0.85 × 0.92 × 0.80 × 0.90 × 0.93 × 0.96 ≈ 0.50` — but 6 of those 7 categories have a high canonical weight, so the composite floor stays high). 30% of eligible records end up with `abbreviation_inconsistency = mixed`, independent of the typographic stack.

---

## 7. Deferred — Spec Layer 7 (document composition)

Multi-paragraph pathology reports, consult notes, and tumor-board summaries with distractor content are **not planned for v1**. The planned approach when this layer ships (out-of-scope for current phases):

1. Hand-authored ~30–50 document-type scaffolds (pathology, molecular summary, consult, tumor board) with slot-filling for the biomarker section.
2. One-time frontier-model distillation for distractor-prose variety (~500 scaffolds, no ongoing API cost).
3. Biomarker-slot content is rendered by the existing L1–L7 generator — no new label surface area.
4. Optional local-LLM filler (Qwen 2.5 7B / Phi-3) for distractor prose with hard span-preservation verification (see [project_spec.md §4](project_spec.md#4-llm-paraphrase-layer-optional-post-hoc)).

**Explicitly out of scope for v1:**
- Multi-paragraph composition
- Distractor content corpora integration
- Document-length stratification metrics
- LLM paraphrase layer (Phase 4+)
- `bioassert/generator/document.py`, `bioassert/paraphrase/` (directories scaffolded in [project_spec.md §9](project_spec.md#9-repository-layout) but not populated)

---

## 8. What this document is not

- **Not a changelog.** `git log` is authoritative for implementation history. The previous `PHASE3_PLAN.md` execution log has been retired because its sub-phase numbers (3.1, 3.2, ...) were confusable with the layer numbers (Layer 1, Layer 2, ...).
- **Not a config reference.** For the JSON schema that drives variation, see [`config_architecture.md`](config_architecture.md).
- **Not an evaluation plan.** For stratified metrics and baselines, see [project_spec.md §5–§6](project_spec.md#5-evaluation-harness).
- **Not a bug list.** For known noise/grammar issues accepted as background noise, see [`known_issues.md`](known_issues.md).

When the shipped variation/noise stack changes — frame added, noise category added, distribution reweighted, new biomarker alias — this document must be updated in the same PR.
