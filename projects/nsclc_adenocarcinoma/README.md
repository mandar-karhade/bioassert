# NSCLC — Lung Adenocarcinoma

Tier-1 biomarker panel for synthetic corpus generation targeting lung adenocarcinoma.

## Panel

| Gene | Alteration type | Target prevalence (positive) |
|---|---|---|
| EGFR | mutation | ~20% (cohort-weighted Western + East-Asian avg) |
| KRAS | mutation | ~25% |
| ALK | fusion | ~5% |
| ROS1 | fusion | ~1–2% |
| BRAF | mutation | ~1–2% |
| MET | alteration | ~3–5% |
| RET | fusion | ~1–2% |
| ERBB2 (HER2) | mutation | ~2–3% |
| NTRK | fusion | <1% |
| PD-L1 | expression | ~30% (≥50% TPS) |
| TMB | expression | ~15–20% (≥10 mut/Mb) |

Prevalence priors are adenocarcinoma-weighted from PIONEER, AACR GENIE v13, LCMC, and published meta-analyses. For squamous cohorts or pan-tumor panels, author a separate project directory and mix output corpora externally.

## Layout

```
project.json              # project metadata + config paths
configs/
  common_variations.json  # shared surface vocab (assertion verbs, phrases, technical noise)
  biomarkers.json         # per-gene panel config
references/               # (future) prevalence_sources.bib, citation files
outputs/                  # versioned run dirs, gitignored
```

## Generate

```bash
bioassert generate --project projects/nsclc_adenocarcinoma --n 10000 --tag baseline
```

Each run writes into `outputs/run_NNN_<tag>_<UTC-timestamp>/` containing `corpus.jsonl`, `prevalence_report.json`, `manifest.json`, and a `snapshot/` of the configs used.
