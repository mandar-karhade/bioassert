# bioassert

Synthetic corpus generator for biomarker assertion extraction in non-small-cell lung cancer (NSCLC).

- Package / repo: `bioassert`
- Dataset release (planned v1): `BioAssert-NSCLC-v1`
- Status: Phase 2b — full Tier 1 panel (11 biomarkers), L1 + L2 rendering

## Quickstart

```bash
uv sync
uv run pytest
uv pip install -e .
bioassert generate \
  --project projects/nsclc_adenocarcinoma \
  --n 50000 \
  --seed 42 \
  --l2-fraction 0.3 \
  --tag smoke
```

Each run lands in its own versioned folder under the project's `outputs/`:

```
projects/nsclc_adenocarcinoma/outputs/
  run_001_smoke_20260423-230756/
    corpus.jsonl               # generated records
    prevalence_report.json     # observed vs. configured status rates
    manifest.json              # CLI args, seed, UTC timestamp, git SHA, version
    snapshot/                  # byte-identical copy of configs + project.json
      configs/{common_variations,biomarkers}.json
      project.json
```

Every invocation is reproducible from the snapshot alone — the configs that produced a corpus travel with it.

## Project layout

A **project** is a self-contained directory describing one corpus-generation
target (configs, references, output history). The shipped example targets
NSCLC lung adenocarcinoma; additional cohorts or future diseases live as
sibling directories under `projects/`:

```
projects/
  nsclc_adenocarcinoma/
    project.json          # metadata + config paths (schema_type: "biomarker")
    configs/
      common_variations.json
      biomarkers.json
    references/           # citations, prevalence sources
    outputs/              # gitignored; versioned run dirs
    README.md
```

See [projects/nsclc_adenocarcinoma/README.md](projects/nsclc_adenocarcinoma/README.md) for the shipped cohort panel.

## Documents

- [docs/project_spec.md](docs/project_spec.md) — full project spec (seven-layer architecture, motivation, roadmap)
- [docs/config_architecture.md](docs/config_architecture.md) — live config-driven architecture (Phase 2)
- [docs/conversation_context.md](docs/conversation_context.md) — design rationale and decision log
- [docs/nsclc_biomarker_reference.md](docs/nsclc_biomarker_reference.md) — clinical reference for the Tier 1 panel
- [docs/variation_noise_profile.md](docs/variation_noise_profile.md) — canonical reference for every shipped complexity level, vocabulary category, noise transform, and prevalence axis
- [docs/known_issues.md](docs/known_issues.md) — tracked rendering variation bugs (with per-bug rates and deferred-fix notes)
- [docs/architecture.md](docs/architecture.md), [docs/complexity_levels.md](docs/complexity_levels.md), [docs/limitations.md](docs/limitations.md) — WIP stubs
- [docs/archived/](docs/archived/) — superseded Phase 1 artifacts (kickoff prompt, old package README)

## Core principle

Ground truth is the structured assertion. Surface text is rendered from assertions. Labels are correct by construction.

## License

[PolyForm Noncommercial 1.0.0](LICENSE) — free for noncommercial use (research, teaching, personal projects). Commercial use requires a separate license.
