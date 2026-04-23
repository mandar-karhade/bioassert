# Biomarker Assertion Corpus Generator — Project Spec

**Project codename:** `bioassert` (package name, PyPI)
**Dataset release name:** `BioAssert-NSCLC-v1` (human-readable), `bioassert_nsclc_v1` (identifier)
**Status:** Phase 2b shipped — full Tier 1 panel (11 biomarkers), L1 + L2 rendering, 50K mixed corpus; Phase 3 (L3+ compound sentences, document-level composition) planning
**Target:** Publishable synthetic corpus + benchmarking study for clinical biomarker extraction, starting with lung cancer (NSCLC)

---

## 1. Problem Statement

Generic NER/RE models (GLiNER, REBEL, etc.) are evaluated using micro-F1 across aggregated relation types. This evaluation methodology hides catastrophic failures on the relation structures that matter in clinical domains — specifically compound biomarker assertions.

Example of what breaks generic models:
> "EGFR L858R positive on NGS, ALK negative by FISH, ROS1 rearrangement not detected, PD-L1 TPS 80% (Clone 22C3), KRAS G12C detected."

This single sentence contains 5 n-ary assertions with heterogeneous structure (gene × variant × status × measurement × method), list distribution, and domain-specific measurements. Generic RE models extracting `(subject, predicate, object)` triples cannot represent this structure and will miss negation polarity, measurement semantics, and method attribution.

**Gap in the field:** No publicly available synthetic biomarker assertion corpus with:
- Controlled compositional complexity
- Clinically calibrated prevalence rates
- Perfect-by-construction ground truth labels
- Stratified evaluation splits
- Transfer validation to real clinical text

This project fills that gap.

---

## 2. Core Architectural Principle

> **Ground truth is the structured assertion, not the surface text. Surface text is generated from assertions with controlled variation. This guarantees label quality that no manual annotation process can match.**

Never invert this. The assertion graph comes first; the text is rendered from it.

---

## 3. Seven-Layer Architecture

### Layer 1 — Ontology / Biomarker Knowledge Base

Structured YAML/JSON entries for each biomarker. Each entry is a full profile:

```yaml
EGFR:
  canonical: "Epidermal Growth Factor Receptor"
  abbreviations: [EGFR, HER1, ERBB1]
  spelling_variants: [EGF-R, EGF R, EGFr, egfr]
  gene_family: ErbB
  common_variants:
    - id: L858R
      type: point_mutation
      exon: 21
      category: activating
      aliases: ["L858R", "Leu858Arg", "c.2573T>G"]
    - id: T790M
      type: point_mutation
      exon: 20
      category: resistance
      aliases: ["T790M", "Thr790Met"]
    - id: exon_19_deletion
      type: deletion
      aliases: ["exon 19 del", "ex19del", "del19", "E19del"]
    - id: exon_20_insertion
      type: insertion
      aliases: ["exon 20 ins", "ex20ins", "E20ins"]
  test_methods: [NGS, PCR, FISH, IHC, Sanger, ddPCR]
  result_vocabulary:
    positive: [positive, detected, mutated, present, "+", pos, identified]
    negative: [negative, "not detected", "wild-type", wt, absent, "-", neg, "no mutation"]
    equivocal: [equivocal, indeterminate, inconclusive, "failed QC"]
  measurement_types: [qualitative]
  clinical_associations: [NSCLC, adenocarcinoma]
  prevalence:
    NSCLC_western: 0.15
    NSCLC_east_asian: 0.45
    adenocarcinoma: 0.20
    squamous_cell: 0.03
  prevalence_sources:
    - "Shi et al. 2014 PIONEER study"
    - "AACR Project GENIE v13"
```

**Target scope for v1:** 15–20 biomarkers covering lung cancer completely.
- EGFR, ALK, ROS1, BRAF, KRAS, MET, RET, HER2/ERBB2
- NTRK1, NTRK2, NTRK3
- PD-L1, TMB, MSI
- NRG1
- TP53 (common co-mutation, relevant for context)
- STK11, KEAP1 (resistance context)

**Deliverable:** `biomarkers.yaml` — the asset that makes the rest defensible.

### Layer 2 — Canonical Assertion Schema

Every ground-truth label is a structured assertion object. Use Pydantic.

```python
from pydantic import BaseModel
from typing import Optional, Literal

class Measurement(BaseModel):
    value: float
    unit: str  # "percent", "mutations_per_mb", etc.
    clone: Optional[str] = None  # for IHC antibody clones like "22C3"

class BiomarkerAssertion(BaseModel):
    assertion_id: str
    patient_ref: str
    gene: str                                      # canonical from ontology
    variant: Optional[str] = None                  # canonical variant ID
    status: Literal["positive", "negative", "equivocal", "not_tested", "pending"]
    test_method: Optional[str] = None
    measurement: Optional[Measurement] = None
    temporality: Literal["current", "prior", "pending", "historical"] = "current"
    certainty: Literal["confirmed", "suspected", "rule_out"] = "confirmed"
    negation_scope: bool = False
    source_clause_id: str                          # trace back to surface text
    char_span: Optional[tuple[int, int]] = None    # for span-level eval
```

### Layer 3 — Compositional Sentence Grammar

Probabilistic context-free grammar or constraint-based generator. Four sub-layers:

**3a. Simple assertion realization** (~50+ frames):
```
[GENE] [is|was] [STATUS_TERM]
[GENE] mutation [STATUS_TERM]
[STATUS_TERM] for [GENE]
[GENE] testing returned [STATUS_TERM]
Result: [GENE] [STATUS_TERM]
```

**3b. Compound statement templates** (~50+ frames):
```
[GENE1], [GENE2], and [GENE3] were all [STATUS_TERM]
[GENE1] [STATUS_TERM]; [GENE2] [STATUS_TERM]
Molecular testing revealed [GENE1] [STATUS_TERM] and [GENE2] [STATUS_TERM]
NGS panel: [GENE1] [STATUS_TERM], [GENE2] [STATUS_TERM], [GENE3] [STATUS_TERM]
Comprehensive genomic profiling: [list]
```

**3c. Negation and assertion modifiers** (~30+ frames):
```
No evidence of [GENE] mutation
[GENE] mutation was not detected
Patient tested negative for [GENE]
Ruling out [GENE] alteration
Suspected [GENE] mutation, pending confirmation
Previously positive for [GENE], repeat testing [STATUS]
```

**3d. Measurement-bearing templates** (PD-L1, TMB, MSI):
```
PD-L1 TPS [VALUE]%
PD-L1 expression: [VALUE]% (Clone [CLONE])
PD-L1 CPS [VALUE]
TMB: [VALUE] mutations/Mb
MSI status: [STATUS]
```

**3e. Variant-bearing templates:**
```
[GENE] [VARIANT] detected
[GENE] mutation ([VARIANT]) identified by [METHOD]
[VARIANT] mutation in [GENE]
[GENE] c.[NOTATION] (p.[PROTEIN_CHANGE])
```

**Target:** 200–400 total sentence frames across all categories.

### Layer 4 — Complexity Stratification

The single methodological decision that differentiates this from every other synthetic data paper. Define formal complexity levels:

| Level | Description | Example |
|-------|-------------|---------|
| L1 | Single assertion, canonical name, explicit status | "EGFR was positive" |
| L2 | Lexical variation, abbreviations, synonyms | "EGF-R mutation detected" |
| L3 | List distribution, one status shared | "EGFR, ALK, and ROS1 were negative" |
| L4 | Heterogeneous compound (mixed statuses) | "EGFR positive, ALK negative, ROS1 negative" |
| L5 | Negation scope complexity | "No EGFR or ALK mutations but BRAF V600E detected" |
| L6 | Temporal/certainty qualification | "Previously tested negative for EGFR; current testing suspected positive" |
| L7 | Cross-sentence coreference | "Patient underwent NGS. Results showed EGFR L858R. No other actionable alterations were identified." |

**Generate N examples per level per biomarker combination. Report stratified F1.**

### Layer 5 — Technical / Surface Noise

Applied with controlled probability:
- Whitespace variants (single, double, tab, newline mid-sentence)
- Hyphenation (EGF-R vs EGFR vs EGF R vs EGF\nR)
- Punctuation (semicolons, commas, newlines as separators)
- Case variants (egfr, EGFR, Egfr, eGFR)
- OCR-style corruption (optional adversarial split)
- PDF-extraction line-break artifacts
- Abbreviation inconsistency within same document

Report performance with noise on vs off.

### Layer 6 — Probabilistic Clinical Realism

Pull positive/negative/equivocal/not-tested distributions from published NSCLC
molecular epidemiology and condition them on a simulated `PatientProfile`.

**Scope: prevalence sampling only.** Patient attributes influence *which*
status is sampled for a given biomarker — they do **not** drive surface-text
variation. We deliberately do not vary sentence frames, tone, or vocabulary
based on patient demographics. `PatientProfile` axes act exclusively as
covariates on the status distribution.

Implemented axes (see [`bioassert/generator/patient_sampler.py`](../bioassert/generator/patient_sampler.py)):

| Axis        | Values                                     |
|-------------|--------------------------------------------|
| histology   | adenocarcinoma, squamous, other            |
| ethnicity   | western, east_asian, other                 |
| smoking     | smoker, former_smoker, nonsmoker           |
| age_group   | older, younger                             |
| sex         | male, female                               |

Population-key cascade (§2.3): keys are joined as
`{histology}_{ethnicity}_{smoking}_{age_group}_{sex}` and progressively
dropped right-to-left until a biomarker config match is found. Terminal
fallback is `{histology}` alone.

**Not implemented, not planned for v1:** stage (I–IV), grade, comorbidities,
treatment line, site-of-biopsy, or any axis that would drive surface-text
variation beyond status prevalence.

Prevalence table (ground source: PIONEER, GENIE, LCMC, published meta-analyses):
- EGFR: 15% Western / 45% East Asian (adenocarcinoma)
- KRAS: 25–30% (smokers, adenocarcinoma)
- ALK: 5%
- ROS1: 1–2%
- BRAF V600E: 1–2%
- MET exon 14: 3–4%
- RET: 1–2%
- HER2: 2–3%
- PD-L1 ≥50%: 30%
- TMB high (≥10 mut/Mb): 15–20%
- Co-mutations: TP53 with EGFR ~55%, STK11 with KRAS ~20%

Cite every prevalence source in the dataset card. This is the credibility layer.

### Layer 7 — Document-Level Composition

Single sentences aren't enough. Compose:
- Pathology report fragments (1–3 paragraphs)
- Molecular testing summaries (structured-ish)
- Oncology consultation notes (biomarker section + irrelevant sections)
- Tumor board summaries
- Progress notes with biomarker mentions embedded

Each document has:
- Known biomarker-relevant spans (labeled)
- Known irrelevant spans (distractor content — bowel movements, social history, etc.)

This lets you evaluate whether models confuse signal and noise.

---

## 4. LLM Paraphrase Layer (Optional, Post-Hoc)

**Rule-based grammar first. LLM paraphrase as a diversity booster, never as the generator.**

Workflow:
1. Generate assertion + rule-based sentence
2. Pass sentence to local LLM (Qwen 2.5 7B or Phi-3) with instruction: "Paraphrase while preserving these exact spans: [GENE=EGFR, STATUS=positive, VARIANT=L858R]"
3. Verify paraphrase still contains labeled spans (regex check)
4. If labels preserved → keep paraphrase; if not → discard and retry

This gives lexical diversity without hallucinating assertions. Key methodological point: LLM never sees or invents the assertion structure — only rewords sentences around preserved anchors.

---

## 5. Evaluation Harness

Ship this as part of the release. Don't make other researchers build it.

```python
from bioassert_eval import evaluate_model, StratifiedReport

report = evaluate_model(
    model=my_model,
    dataset="bioassert_nsclc_v1",
    splits=["L1", "L2", "L3", "L4", "L5", "L6", "L7"],
    metrics=["span_f1", "assertion_f1", "negation_accuracy",
             "variant_accuracy", "method_attribution"],
)

report.print_stratified()
report.to_latex("results.tex")
```

Metrics to compute:
- **Span-level F1** per entity type (gene, variant, value, method)
- **Assertion-level F1** — true positive requires ALL n-tuple elements correct
- **Partial-credit F1** — for diagnostic analysis
- **Negation accuracy** — specifically score polarity
- **Variant binding accuracy** — variant correctly attached to correct gene
- **Method attribution accuracy** — test method correctly attached to correct assertion
- **List-distribution recall** — when one status applies to multiple genes, did all get extracted

Stratified reports by:
- Complexity level
- Biomarker type
- Noise on/off
- Document length

---

## 6. Models to Benchmark

Evaluate these on the corpus, both zero-shot and fine-tuned where applicable:

**Specialized small models:**
- GLiNER (urchade/gliner_large, gliner-bi-large-v2.0)
- GLiNER-biomed
- REBEL (Babelscape/rebel-large)
- ReLiK (sapienzanlp/relik-cie-small)
- GLiREL for relation extraction
- ClinicalBERT / BioClinicalBERT fine-tuned

**Rule-based baseline:**
- NegEx / ConText for assertion detection
- Regex + ontology lookup for the simplest baseline (must be included — may win on L1)

**LLM baselines:**
- Qwen 2.5 7B with constrained JSON output (Pydantic schema)
- Phi-3 mini with same schema
- Optionally a frontier LLM as ceiling reference (GPT-4o or Claude) — single eval run, not a cost-ongoing baseline

**The key story to tell:**
> Model X achieves 92% F1 on L1 but 34% F1 on L4. A fine-tuned small model achieves 78% on both. Stratified evaluation reveals what aggregate scores hide.

---

## 7. Transfer-to-Real-Text Validation

**Non-negotiable for publication.** A synthetic dataset is only useful if training on it improves performance on real clinical text.

Acquire a held-out real-text test set:
- n2c2 (i2b2) assertion challenge data (public after DUA)
- BioCreative corpora (BC5CDR, ChemProt)
- Collaborator dataset from a clinical partner
- MIMIC-III/IV discharge summaries with manual biomarker annotation

Show:
1. Model trained on synthetic only → performance on real
2. Model trained on real only → performance on real
3. Model trained on synthetic + fine-tuned on small real set → performance on real

The expected story: synthetic pretraining + small real fine-tuning > real-only training, at a fraction of the annotation cost.

Without this result, the paper is "we built an elaborate template matcher."

---

## 8. Publication Plan

**Two papers from one codebase.**

**Paper 1 — Dataset / Resource paper:**
- Venue: *Scientific Data*, ACL/EMNLP data track, BioNLP workshop, or AMIA Informatics Summit
- Contribution: the corpus, the generator, the ontology, the evaluation harness
- Key selling points: perfect labels by construction, complexity stratification, calibrated prevalence

**Paper 2 — Benchmarking study:**
- Venue: JAMIA, *Journal of Biomedical Informatics*, or ACL/EMNLP main
- Contribution: stratified evaluation of existing models, demonstration of synthetic-to-real transfer
- Key selling points: exposes where generic models fail, provides a rigorous baseline

Dataset citation accrual — if released well with good tooling — compounds over years.

---

## 9. Repository Layout

```
bioassert/
├── README.md
├── pyproject.toml                # uv + hatchling, PEP 621
├── uv.lock
├── LICENSE                       # Apache-2.0
├── config/
│   ├── common_variations.json    # shared surface variations (Phase 2a)
│   ├── biomarkers.json           # per-biomarker config (Phase 2a)
│   └── prevalence_sources.bib    # NSCLC epi citations
├── ontology/
│   ├── biomarkers.yaml           # Phase 1 legacy; deprecated after Phase 2a migration
│   └── sentence_frames.yaml      # grammar templates
├── bioassert/                    # Python package
│   ├── __init__.py
│   ├── schema.py                 # Pydantic assertion models
│   ├── config/
│   │   ├── __init__.py
│   │   ├── schema.py             # discriminated-union Pydantic models
│   │   ├── loader.py             # JSON loading
│   │   └── validator.py          # 12-check validation (CONFIG_ARCHITECTURE.md §7)
│   ├── generator/
│   │   ├── __init__.py
│   │   ├── patient_sampler.py    # PatientProfile + population key cascade
│   │   ├── sampler.py            # weighted sampling, render_constraints, clone Bernoulli
│   │   ├── renderer.py           # L1/L2 frames + placeholder substitution + span tracking
│   │   ├── post_process.py       # technical_noise transformations
│   │   ├── compound_builder.py   # L3+ list/heterogeneous (Phase 3, not yet landed)
│   │   ├── negation.py           # L5 negation scope (Phase 3, not yet landed)
│   │   ├── document.py           # Phase 4 document-level composition (not yet landed)
│   │   └── complexity.py         # L3–L7 stratification orchestrator (Phase 3+, not yet landed)
│   ├── paraphrase/                # optional, Phase 3+
│   │   ├── llm_paraphraser.py    # span-preserving
│   │   └── verifier.py
│   └── io/                        # Phase 3+
│       ├── exporters.py          # JSONL, HuggingFace datasets, CoNLL
│       └── splits.py             # train/val/test stratified
├── bioassert_eval/               # Evaluation harness package
│   ├── __init__.py
│   ├── metrics.py
│   ├── stratified.py
│   ├── adapters/
│   │   ├── gliner.py
│   │   ├── rebel.py
│   │   ├── relik.py
│   │   ├── regex_baseline.py
│   │   ├── negex.py
│   │   └── llm_structured.py
│   └── reports.py
├── experiments/
│   ├── 01_zero_shot_eval.py
│   ├── 02_finetune_gliner.py
│   ├── 03_finetune_rebel.py
│   ├── 04_transfer_to_n2c2.py
│   └── configs/
├── bioconfigs/                    # canonical JSON config (replaces Phase 1 YAML ontology)
│   ├── biomarkers.json
│   └── common_variations.json
├── datasets/
│   ├── v1_phase2a/               # 10K L1 mutation/fusion corpus
│   └── v1_phase2b/               # 50K mixed L1+L2 corpus (full Tier 1 panel)
├── docs/
│   ├── CONFIG_ARCHITECTURE.md    # Phase 2 config-driven architecture
│   ├── architecture.md           # WIP stub
│   ├── complexity_levels.md      # WIP stub
│   ├── limitations.md            # WIP stub
│   └── archived/                 # Phase 1 historical docs
├── tests/
│   ├── test_config_loader.py
│   ├── test_config_validator.py
│   ├── test_patient_sampler.py
│   ├── test_sampler.py
│   ├── test_renderer.py
│   ├── test_post_process.py
│   └── test_phase2a_preservation.py
└── scripts/
    └── generate_corpus_v1.py     # mixed L1+L2 corpus generator
```

---

## 10. Development Phases

**Phase 0 — Scoping & collaborators (Week 0)**
- Confirm clinical collaborator for ontology review
- Confirm access to real-text validation set
- Decide license (Apache-2.0 recommended for max adoption)

**Phase 1 — Architecture smoke test (Weeks 1–2)**
- Minimal YAML ontology: 3 biomarkers (EGFR, ALK, KRAS), minimal surface variation
- Pydantic `BiomarkerAssertion` + `Measurement` schema
- Prevalence-aware sampler (basic)
- L1 grammar with 10–15 sentence frames
- 1,000 L1 records emitted as JSONL with character spans
- Assertion-preservation tests passing green
- Adversarial layer stub (lexical_distractor, polarity_false_friend, irrelevant_medical)
- `docs/` stubs for architecture.md, complexity_levels.md, prevalence_calibration.md, limitations.md
- **Exit criterion:** end-to-end architecture validated. Do not expand complexity until this passes.

**Phase 2a — Probabilistic variation layer (Weeks 3–4)**
- Migrate the three Phase 1 biomarkers from YAML to the `biomarkers.json` + `common_variations.json` schema
- Pydantic discriminated-union models for `$schema_type` dispatch (weighted_variations, weighted_variations_with_attachment, post_process_transformations)
- Loader with all 12 validation checks from CONFIG_ARCHITECTURE.md Section 7
- Patient profile sampler + population key cascade lookup
- Sampler handling `render_constraints.require_biomarker_name_forms` coupling
- Sampler handling `clone_attribution` Bernoulli attachment
- Renderer handling placeholder substitution ({gene}, {method}, {result}, {value})
- Post-process pass applying `technical_noise` transformations
- Rewritten assertion-preservation tests against new schema
- **Exit criterion:** 10K records emitted, all preservation tests pass, prevalence hits population-stratified targets within 2σ

**Phase 2b — Full Tier 1 biomarker coverage (Weeks 5–6)**
- Extend `biomarkers.json` from 3 to all 11 Tier 1 biomarkers (ROS1, BRAF, MET, RET, ERBB2, NTRK, PD-L1, TMB added)
- Clinical collaborator review of ontology
- Source prevalence data with citations (bib file)
- Extend grammar to handle expression biomarkers (PD-L1, TMB with measurement rendering)
- L2 lexical variation fully exercised via the probability config
- 50K L1 + L2 records across the full Tier 1 panel

**Phase 3 — Compositional complexity (Weeks 7–9)**
- L3/L4 list distribution + heterogeneous compound
- L5 negation scope
- L6 temporal/certainty
- L7 cross-sentence coreference
- Noise layer (Layer 5)

**Phase 4 — Document composition (Weeks 10–11)**
- Multi-paragraph documents with distractor content
- Realistic pathology/consult note structure

**Phase 5 — Evaluation harness (Weeks 12–13)**
- Metrics library
- Model adapters for GLiNER, REBEL, ReLiK, NegEx, regex baseline
- Stratified reporting

**Phase 6 — Benchmarking (Weeks 14–16)**
- Zero-shot evaluation of all models
- Fine-tuning experiments on GLiNER, REBEL
- LLM baselines with constrained JSON

**Phase 7 — Transfer validation (Weeks 17–18)**
- Evaluate on real-text test set
- Synthetic-only vs synthetic+real vs real-only

**Phase 8 — Paper writing + release (Weeks 19–22)**
- Dataset card
- HuggingFace dataset upload
- arXiv preprint
- Submission to target venue

Total: ~5–6 months focused, assuming clinical collaborator available.

---

## 11. Non-Negotiables

These are the things that separate a publishable contribution from a curiosity:

1. **Ground truth is generated before surface text.** Never invert this.
2. **Complexity stratification is reported explicitly.** No aggregate-only F1.
3. **Prevalence is clinically calibrated with cited sources.**
4. **Transfer to real text is demonstrated.** Not just synthetic-on-synthetic.
5. **The generator and eval harness are released with the dataset.** Lower activation energy for reuse.
6. **Limitations are honestly documented.** No handwriting artifacts, no pathologist idiosyncrasies, US-English bias, etc.
7. **Adversarial examples are included.** "Patient is positive about the treatment plan" should be in the test set.
8. **Rule-based baselines are included.** Regex + ontology may win on L1, which is a feature of the eval, not a bug.

---

## 12. Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| Grammar collapses into obvious templates, transformer overfits trivially, transfer fails | 200+ sentence frames per level, probabilistic backoff, optional LLM paraphrase post-hoc |
| No access to real-text validation data | Start with public n2c2 after DUA; budget time for collaborator agreement |
| Clinical collaborator unavailable or slow | Start with published literature-only ontology; schedule clinical review as discrete milestone |
| Reviewers dismiss as "elaborate regex" | Transfer-to-real results are the answer. Also: stratified eval showing where existing SOTA fails |
| Scope creep into all of oncology | Stay lung-only for v1. Breast/colorectal/prostate are follow-up papers |
| LLM paraphrase introduces hallucinated assertions | Hard regex verification that labeled spans are preserved; discard-and-retry if not |

---

## 13. Initial Deliverable (Historical — Phase 1 Kickoff)

The original Phase 1 kickoff deliverable (now completed and superseded) targeted:

1. Repo scaffold with `pyproject.toml`, directory layout
2. `biomarkers.yaml` with 3 biomarkers fully profiled (EGFR, ALK, KRAS)
3. `schema.py` with Pydantic `BiomarkerAssertion` + `Measurement` models
4. `assertion_sampler.py` with prevalence-aware sampling for the 3 biomarkers
5. `grammar.py` with 10–15 L1 sentence frames and a working renderer
6. `scripts/generate_corpus.py` emitting 1,000 L1 assertions + sentences as JSONL
7. Unit tests validating that generated sentences always contain the labeled spans

Phase 1 validated the architecture end-to-end. The codebase has since moved to the JSON-driven config pipeline described in [CONFIG_ARCHITECTURE.md](CONFIG_ARCHITECTURE.md) — the YAML ontology, hand-rolled grammar module, and Pydantic schema have been replaced by `bioconfigs/` JSON files, `bioassert/config/schema.py`, and the data-driven renderer in `bioassert/generator/renderer.py`. Current status is summarized at the top of this document.
