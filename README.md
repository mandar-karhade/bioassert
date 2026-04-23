# bioassert

Synthetic corpus generator for biomarker assertion extraction in non-small-cell lung cancer (NSCLC).

- Package / repo: `bioassert`
- Dataset release (planned v1): `BioAssert-NSCLC-v1`
- Status: Phase 2b — full Tier 1 panel (11 biomarkers), L1 + L2 rendering

## Quickstart

```bash
uv sync
uv run pytest
PYTHONPATH=. uv run python scripts/generate_corpus_v1.py --n 50000 --l2-fraction 0.3
```

The corpus lands at `datasets/v1_phase2b/corpus.jsonl` with a `prevalence_report.json` summarizing observed vs configured status rates.

## Documents

- [docs/PROJECT_SPEC.md](docs/PROJECT_SPEC.md) — full project spec (seven-layer architecture, motivation, roadmap)
- [docs/CONFIG_ARCHITECTURE.md](docs/CONFIG_ARCHITECTURE.md) — live config-driven architecture (Phase 2)
- [docs/CONVERSATION_CONTEXT.md](docs/CONVERSATION_CONTEXT.md) — design rationale and decision log
- [docs/NSCLC_BIOMARKER_REFERENCE.md](docs/NSCLC_BIOMARKER_REFERENCE.md) — clinical reference for the Tier 1 panel
- [docs/PHASE3_PLAN.md](docs/PHASE3_PLAN.md) — Phase 3 sub-phase breakdown (L3–L7, noise expansion)
- [docs/known_issues.md](docs/known_issues.md) — tracked rendering variation bugs (with per-bug rates and deferred-fix notes)
- [docs/architecture.md](docs/architecture.md), [docs/complexity_levels.md](docs/complexity_levels.md), [docs/limitations.md](docs/limitations.md) — WIP stubs
- [docs/archived/](docs/archived/) — superseded Phase 1 artifacts (kickoff prompt, old package README)

## Core principle

Ground truth is the structured assertion. Surface text is rendered from assertions. Labels are correct by construction.

## License

[PolyForm Noncommercial 1.0.0](LICENSE) — free for noncommercial use (research, teaching, personal projects). Commercial use requires a separate license.
