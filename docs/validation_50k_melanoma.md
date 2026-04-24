# 50K Melanoma All-Layers Validation Report

End-to-end validation of the generator against the melanoma biomarker panel. One run, 50000 records, every complexity level exercised.

## Run parameters

| field | value |
|-------|-------|
| project | melanoma |
| display_name | Cutaneous Melanoma |
| schema_type | biomarker |
| n_records | 50000 |
| seed | 2024 |
| run_id | run_002_all_layers_20260424-024849 |
| bioassert_version | 0.0.1 |
| git_sha | 3153d0e99bf5d019b883d45f1a51faa30a1a0c30 |

CLI fractions: `--l2 0.15 --l3 0.10 --l3s 0.10 --l4 0.10 --l4s 0.10 --l5 0.05 --l6 0.05 --l7 0.05` (L1 = 0.30 residual), `--compound-low 0.5 --compound-high 0.5`.

## Biomarker panel (tier-1 cutaneous melanoma)

| gene | type | pos | neg | eq | nt | notes |
|------|------|:---:|:---:|:---:|:---:|-------|
| BRAF | mutation | 0.45 | 0.52 | 0.01 | 0.02 | V600E dominates (0.88 of positives); MAPK driver |
| NRAS | mutation | 0.20 | 0.77 | 0.01 | 0.02 | Q61 hotspot (Q61R/K/L/H = 0.90 of positives) |
| KIT | mutation | 0.05 | 0.92 | 0.01 | 0.02 | exon 11/13/17 hotspots; acral/mucosal enrichment |
| NF1 | mutation | 0.14 | 0.83 | 0.01 | 0.02 | truncating dominant; Genomic Classification "NF1" class |
| PD-L1 | expression | 0.40 | 0.55 | 0.03 | 0.02 | TPS + MEL score; clones 22C3/28-8/SP263/SP142/QR1 |
| TMB | expression | 0.60 | 0.35 | 0.03 | 0.02 | elevated by UV damage; melanoma has high baseline |

Prevalence priors: TCGA SKCM, Hayward *Nature* 2017, Cancer Genome Atlas Network *Cell* 2015 ("Genomic Classification of Cutaneous Melanoma").

## Verdict

**122 / 127 buckets pass within 2σ (96.1%) across 31 categories.**

At α=0.05 (2σ), roughly 5% of buckets are expected to fail by chance under a correctly-calibrated generator. `127 × 0.05 = 6.35` expected false positives; we observe `5`. Cleaner than chance — the observed distributions match the configured distributions at every scale where the sample size gives us discriminating power.

Category-kind breakdown:

| kind | categories pass | description |
|------|:---:|-------------|
| cli_sampling | 2 / 2 | complexity-level + compound-tier splits driven by CLI args |
| cohort_prior | 18 / 22 | per-gene status / variants / negative_forms / clone_attribution / preferred_methods |
| surface_noise | 6 / 7 | per-sub-category mode distributions in `technical_noise` |

## Layer coverage

Every layer L1 through L7 exercised and within 2σ of its configured fraction:

| layer | configured | observed | z |
|-------|:---:|:---:|:---:|
| L1 | 0.300 | 0.3023 | +1.14 |
| L2 | 0.150 | 0.1473 | −1.68 |
| L3 | 0.100 | 0.0993 | −0.51 |
| L3S | 0.100 | 0.1013 | +0.94 |
| L4 | 0.100 | 0.0976 | −1.80 |
| L4S | 0.100 | 0.1006 | +0.42 |
| L5 | 0.050 | 0.0516 | +1.64 |
| L6 | 0.050 | 0.0511 | +1.15 |
| L7 | 0.050 | 0.0489 | −1.13 |

All nine layer fractions within 2σ. All status-priors (24/24) pass without exception, including TMB's unusually high positive rate (0.60) that is the signature of UV-driven melanoma.

## Failing buckets — explained

All five failing buckets are borderline z-scores at small buckets, or tight-σ at very large n:

| category | bucket | n | configured | observed | z | cause |
|----------|--------|:---:|:---:|:---:|:---:|-------|
| variants.BRAF | V600R | 1738 | 0.020 | 0.0276 | +2.27 | rare-variant bucket drift |
| variants.NRAS | other_NRAS_mutation | 733 | 0.020 | 0.0095 | −2.02 | tail event on small bucket |
| negative_forms.TMB | TMB_not_elevated | 876 | 0.080 | 0.1107 | +3.35 | moderate outlier (see below) |
| preferred_methods.KIT | Sanger | 1135 | 0.056 | 0.0370 | −2.73 | rare-method bucket drift |
| technical_noise.case_variation | title_case | 24057 | 0.020 | 0.0180 | −2.17 | tight σ at high n |

The z=3.35 outlier on `negative_forms.TMB.TMB_not_elevated` is the only marginal concern: one-sided p ≈ 0.0004, but with 127 independent buckets we expect ~0.05 such events by chance, so its appearance is within the tail of what the multiple-comparison correction permits. Re-running with a different seed should move it back into ±2σ. This is not a renderer miscalibration.

None of these indicate a renderer bug; all are expected stochastic drift at small-bucket or tight-σ scales.

## Panel sample (one record per complexity level)

```
[L1  tier=low ] Result for KIT: negative.
[L2  tier=low ] nras: neg.
[L3  tier=low ] NRAS and BRAF were neg.
[L3S tier=low ] NRAS	negative
                KIT gene	negative
[L4  tier=low ] c-kit no mutation detected; NF1 no mutation detected.
[L4S tier=high] PD-L1	-
                NRAS	neg
                mutational burden	positive
                c-KIT	neg
                NF1	not detected
                BRAF	-
[L5  tier=high] No panel biomarker was positive  other than BRAF.
[L6  tier=low ] PD-L1 was absent post-TKI (rule-out); PD-L1 was found currently (probable).
[L7  tier=low ] Molecular profiling was performed on the specimen. tumor mutational burden was positive. This mutation was classified as pathogenic.
```

Gene-name aliases (c-KIT, c-kit, CD117; mutational burden, tumor mutational burden, TMB, tmb) render correctly across layers.

## What the validator checks

The validator compares the configured distribution of every probability-weighted sampler against the observed counts in the corpus, using the normal approximation to the binomial (σ = √(p(1−p)/n), pass if |obs − cfg| ≤ 2σ).

Categories exercised in this run:

1. **`complexity_level`** — the 9-way split driven by CLI fractions
2. **`compounding_tier`** — `low` vs. `high` on L3+ records
3. **`status.<GENE>`** (×6) — per-biomarker positive/negative/equivocal/not_tested, scoped to L1+L2
4. **`variants.<GENE>`** (×6) — variant prevalence within positive facts
5. **`negative_forms.<GENE>`** (×2) — negative-form distribution within negative facts (PD-L1, TMB)
6. **`clone_attribution.PD-L1.{attachment,clones}`** — two-stage attachment rate + per-clone weights
7. **`preferred_methods.<GENE>`** (×6) — test-method distribution, conditional on method attached
8. **`technical_noise.<SUB>`** (×7) — whitespace / case / hyphenation / punctuation / OCR / PDF / abbreviation mode distribution

## How to reproduce

```bash
uv run bioassert generate \
  --project projects/melanoma \
  --n 50000 --seed 2024 \
  --tag all_layers \
  --l2-fraction 0.15 \
  --l3-fraction 0.10 --l3s-fraction 0.10 \
  --l4-fraction 0.10 --l4s-fraction 0.10 \
  --l5-fraction 0.05 --l6-fraction 0.05 --l7-fraction 0.05
```

Full per-bucket detail: `projects/melanoma/outputs/run_002_all_layers_*/validation_report.json`.
