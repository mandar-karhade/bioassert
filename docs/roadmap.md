# Roadmap

Planned future work, in priority order. Items here are **not scheduled** — they are the queue of intentionally-deferred ideas.

## Next schema: `lab_test` (non-biomarker domain)

The biomarker schema (gene / variant / status / method / clone / technical-noise) generalizes cleanly to lab chemistry. Queued as three progressively larger deliverables.

### A. Lipid panel — first lab-test project *(priority 1)*

- **Panel**: total cholesterol, LDL, HDL, triglycerides, VLDL, non-HDL
- **Why first**: smallest self-contained panel (6 analytes). Fastest path to a working second `schema_type`. Flushes out which biomarker-specific assumptions have leaked into the loader / renderer / validator.
- **Structure mapping**: biomarker → analyte; status_distribution → result_flag (normal / high / low / critical_H / critical_L); variants → value-bucket (measurement_range anchored to NCEP ATP III thresholds); preferred_methods → assay platform. Name aliases, negative_forms, technical_noise, and all 7 complexity levels carry over unchanged.
- **Exit criterion**: `bioassert generate --project projects/lipid_panel --n 50000` produces a corpus that passes the validator at ≥95% within-2σ, same as biomarker runs.

### B. CBC (Complete Blood Count) — mid-size workhorse *(priority 2)*

- **Panel**: WBC, RBC, HGB, HCT, MCV, MCH, MCHC, RDW, PLT, MPV, neutrophils %, lymphocytes %, monocytes %, eosinophils %, basophils % (≈12–15 analytes)
- **Why second**: stress-tests aliasing ("WBC" / "white count" / "leukocytes") and unit variation (10^3/µL vs. ×10⁹/L) at a scale comparable to the NSCLC/melanoma biomarker panels.
- **Depends on**: A (schema proven). May refine renderer hooks discovered during A.

### C. Cohort-contrast pair — CMP healthy vs. CMP hepatitis *(priority 3)*

- **Panels**: two `projects/cmp_*` directories sharing a common analyte catalog but differing in priors (normal liver enzymes vs. elevated ALT/AST/ALP/bilirubin).
- **Why third**: exercises the cohort-prior mechanic at the project-pair level — proves that the same schema with different priors produces clinically-distinct corpora. Also sets up future multi-cohort mixing.
- **Depends on**: A + B (schema and aliasing patterns stable).

## Not yet prioritized

- Pluggable `schema_type` registry (non-biomarker, non-lab domains: symptoms, signs/findings, imaging reports)
- `bioassert validate <project>` / `bioassert show <project>` CLI subcommands
- Cross-project corpus mixing (`bioassert mix p1 p2 --ratio 0.8 0.2`)
- Docker / container packaging for the CLI
