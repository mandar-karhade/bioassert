# Conversation Context — Why This Project Exists

This document captures the reasoning chain that led to the project spec. Read before starting development. It's the "why" behind the "what."

---

## The Starting Question

**Are there small, specialized models for building knowledge graphs at scale without LLM API costs?**

Yes. KG construction decomposes into ~5 subtasks, each with dedicated small models that beat general-purpose LLMs on speed and often on accuracy. Total budget for a production setup: under 500M parameters, CPU-capable, tens of thousands of documents/hour on a single GPU.

### The Pipeline

1. **Coreference resolution** — `fastcoref`, `LingMess`, `maverick-coref`. Resolves "he", "the company", "it" to their referents. Skipping this is the #1 cause of orphan nodes.
2. **Named Entity Recognition** — GLiNER family. Zero-shot, encoder-only, runs on CPU. Natural language labels at inference time. `knowledgator/gliner-bi-*-v2.0` handles 1000+ entity types with near-constant inference speed using pre-computed label embeddings, up to 130× faster than uni-encoder at 1024 entity types.
3. **Relation Extraction** — REBEL (Babelscape) for Wikidata-style triples, GLiREL for zero-shot custom relations.
4. **Entity Linking / canonicalization** — ReLiK (Sapienza NLP) does EL + RE in one framework. Alternative: sentence-transformers bi-encoder against a custom knowledge base.
5. **Deduplication / clustering** — sentence-transformers (`bge-small`, `all-MiniLM-L6-v2`) + clustering.

All CPU or single consumer GPU. No API costs. Deterministic. Fine-tunable.

### Fine-tuning Evidence

GLiNER is fine-tunable on small domain datasets with published results:
- **Cybersecurity (AnnoCTR)**: fine-tuned GLiNER hit 89.7% precision, 74.3% recall, 80.5% F1, significantly outperforming ChatGPT's zero-shot NER.
- **Catalan**: fine-tuned on 9,242 manually annotated sentences with 13,732 named entities, using batch size 8 and focal loss.
- **PII/PHI (Gretel)**: fine-tuned bi-encoder variants on 50K synthetic documents covering 40+ PII types across 45 domains.
- **Biomedical**: GLiNER-biomed beats all general variants on 8 biomedical NER benchmarks.

Caveat: fine-tuning on 118-example datasets caused performance to drop sharply vs ~20k training rows. Floor is a few thousand examples. Below that you get catastrophic forgetting.

### Existing Pipelines

- **ReLiK + LlamaIndex** — closest to turnkey. NER + EL + RE in one framework.
- **spaCy + GLiNER + REBEL** — DIY composable stack, most common in production.
- **Strwythura (Derwen AI)** — explicit anti-monolithic, anti-LLM-API reference implementation integrating Senzing, Placekey, LanceDB, spaCy, GLiNER, RDFlib, NetworkX, designed to run air-gapped.

**Avoid:** Microsoft GraphRAG, LangChain graph builders, Zep's Graphiti — all LLM-API-native, defeat the cost-avoidance goal.

---

## The Insight That Redirected the Project

> "I am debating if a general completeness score matters or only the relationships that one cares about matters. For example for a task of annotating biomarker results for a patient, whether the patient had good bowel movement post surgery doesn't matter."

This is the critical insight. Generic NER/RE benchmarks report micro-F1 across aggregated relation types. This methodology hides catastrophic failure on the relations that actually matter:

- **Micro-F1 tyranny:** frequent easy classes (PERSON, ORG, "born in") dominate the score. Rare hard classes get lost in the average. A model with 90% micro-F1 can have 40% F1 on the class you actually care about.
- **Macro-F1 doesn't fix this** when the problem is structural complexity, not class frequency.

### Why Biomarker Extraction Breaks Generic RE

The example: "EGFR positive, ALK negative, BRAF V600E, PD-L1 50%"

This is not a triple extraction problem. It's structured attribute extraction disguised as text. Generic RE models like REBEL have:
- No concept of **negation** ("ALK negative" — REBEL extracts `(patient, has_mutation, ALK)` and misses polarity)
- No concept of **uncertainty** (suspected vs confirmed)
- No concept of **temporal scope** (prior vs current treatment)
- No concept of **n-ary relations** (biomarker result is a 4-tuple: patient × gene × status × measurement_method)
- No concept of **list distribution** ("EGFR, ALK, and BRAF were all tested negative")

### What Public Benchmarks Actually Test

- **CrossNER F1 61.5%** = span identification. Did it find "Cristiano Ronaldo" as PERSON. Binary match.
- **REBEL micro-F1 74** = given "X was born in Y", did it produce `(X, place_of_birth, Y)`. Maximally easy sentence structure.

Your multi-biomarker sentence is in a completely different complexity class. Nobody has benchmarked these models on compound clinical assertions because the public datasets don't exist at scale.

Closest existing clinical datasets:
- **n2c2 / i2b2** — medication, adverse events, temporal relations
- **BioRED** — biomedical relations, entity-level only
- **ChemProt / DrugProt** — narrow

None test compound statement distribution of real pathology reports.

### The Right Architecture for Biomarker Extraction

Not REBEL. Not generic NER + RE.

1. **NER**: GLiNER fine-tuned on biomarker entities (gene names, variants, values)
2. **Assertion/negation**: NegEx (2001 rule-based, still SOTA for clinical negation) or fine-tuned ClinicalBERT classifier. Don't rely on NER or RE for this.
3. **Structured extraction for compound statements**: small instruction-tuned model (Qwen 2.5 7B, Phi-3) with constrained JSON output against biomarker schema. This is the one place where generative models earn their cost — n-ary output with qualifiers.
4. **Your own eval set**, stratified by complexity, scored on your schema.

---

## The Synthetic Data Idea (The Project)

> "We need to create a library of biomarker names, abbreviations, variations, values, semantic variations, composite statement building from low to high complexity, technical variations. Start with ground truth and build up variation. This way the data acts as if it's annotated without manual annotation."

This is exactly right. Technically called a **grammar-based synthetic corpus with controlled variation and complexity stratification**.

### Why This Is Methodologically Defensible

Clinical NLP community has been burned by "we used GPT-4 to generate synthetic clinical notes" papers — the data inherits LLM biases, hallucinations, limited awareness of real clinical writing conventions.

A rule-based compositional generator with probabilistic control is defensible in a way that LLM-generated data is not.

### The Key Methodological Decisions

1. **Ground truth is the structured assertion, not the surface text.** Text is rendered from assertions. Labels correct by construction.
2. **Complexity stratification with explicit levels.** Report F1 per stratum, not aggregate. This single decision differentiates this from every other synthetic dataset paper.
3. **Clinically calibrated prevalence.** Pull from published NSCLC molecular epidemiology (PIONEER, GENIE, LCMC). Cite every source.
4. **Transfer to real text is non-negotiable.** Without demonstrating that training on synthetic improves performance on real clinical text, the paper is "an elaborate template matcher."

### Why the LLM Paraphrase Layer Is Safe (If Done Right)

Generate assertion + rule-based sentence first. Then optionally paraphrase with local LLM while preserving labeled spans. Regex verification that spans are still present; discard and retry if not.

LLM never sees or invents the assertion structure — only rewords text around preserved anchors. This gives lexical diversity without hallucinating assertions.

### Publication Positioning

Two separable contributions:
1. **Dataset/resource paper**: the corpus, generator, ontology, eval harness. Venue: *Scientific Data*, BioNLP, AMIA, ACL/EMNLP data track.
2. **Benchmarking study**: stratified evaluation of existing models + synthetic-to-real transfer demonstration. Venue: JAMIA, *Journal of Biomedical Informatics*, or main conference.

### What Separates Good From Mediocre Execution

- 200–400 sentence frames per level (not 50 frames swapped with gene names — transformer will trivially overfit)
- Prevalence calibration as first-class feature with citations
- Adversarial examples in the eval ("Patient is positive about the treatment plan")
- Release evaluation harness with the dataset (lower activation energy for reuse)
- Honest limitations documentation (no handwriting artifacts, US-English bias, etc.)
- Rule-based baselines included (regex + ontology may win L1 — feature of the eval, not bug)

### Realistic Scope for First Paper

- 15–20 biomarkers fully profiled (2–3 weeks with clinical collaborator)
- 200–400 sentence frames across 7 complexity levels (1–2 weeks)
- 50K–200K generated assertions with stratification (generator ~1 week; generation in hours)
- 4–6 models evaluated zero-shot + fine-tuned (2–3 weeks)
- Transfer eval on 200–500 real clinical sentences from n2c2 or collaborator (the bottleneck)
- Paper writing (3–4 weeks)

Total: ~3–6 months focused. Bottleneck is access to real-text validation data.

---

## Strategic Framing

The conventional move in clinical NLP is to scrape real notes and pay annotators. The unconventional move — and the stronger methodological contribution — is to **build a generator whose ground truth is correct by construction**.

That's a different and better story, if executed with enough rigor in the variation layers that reviewers don't dismiss it as templated.

The combination of:
- principled generator
- clinically calibrated prevalence
- complexity-stratified evaluation
- demonstrated transfer to real text

...is a genuinely stronger contribution than most clinical NLP dataset papers.

---

## Personal Context

This project leverages:
- Independent operator profile, content + platform work already running
- AnyWebAlert RSS infrastructure (could be adapted for corpus monitoring/updates)
- Medium blog for companion explanatory content
- Existing interest in AI/LLM space for benchmarking extensions

Publishing the dataset well creates a compounding citation asset. Follow-up papers: extend to breast/colorectal/prostate (each has own molecular landscape), extend to other cancer types, extend evaluation to pharmacogenomics.
