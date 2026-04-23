"""Generate the v1 corpus: mixed L1 + L2 records from the JSON-driven pipeline.

Deterministic (seeded). Output: ``datasets/v1_phase2b/corpus.jsonl`` plus a
``prevalence_report.json`` summarizing observed vs configured status rates
per biomarker × matched population key.

Phase 2b scope covers the full Tier 1 panel: mutation/fusion/composite
biomarkers plus the two expression biomarkers (PD-L1, TMB) with
``negative_forms`` and IHC ``clone_attribution`` rendering. L1 records use
formal prose frames; L2 records use tabular shorthand (``positive_shorthand``
/ ``negation_shorthand`` vocabulary). The L2 mixing ratio is configurable
via ``--l2-fraction``.

Re-run with different ``--seed`` to regenerate.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Iterator

from bioassert.config import load_configs
from bioassert.generator.patient_sampler import (
    PatientProfile,
    resolve_status_distribution,
    sample_patient_profile,
)
from bioassert.generator.post_process import (
    PostProcessedRecord,
    apply_technical_noise,
)
from bioassert.generator.renderer import render_l1_record

ROOT = Path(__file__).resolve().parents[1]
COMMON_PATH = ROOT / "bioconfigs" / "common_variations.json"
BIOMARKERS_PATH = ROOT / "bioconfigs" / "biomarkers.json"
OUTPUT_DIR = ROOT / "datasets" / "v1_phase2b"

PANEL_BIOMARKERS: tuple[str, ...] = (
    "EGFR",
    "KRAS",
    "BRAF",
    "ALK",
    "ROS1",
    "RET",
    "NTRK",
    "MET",
    "ERBB2",
    "PD-L1",
    "TMB",
)

DEFAULT_N = 50_000
DEFAULT_SEED = 42
DEFAULT_L2_FRACTION = 0.3


def _iter_records(
    n: int, seed: int, l2_fraction: float
) -> Iterator[tuple[PatientProfile, str, str, PostProcessedRecord]]:
    common, biomarkers = load_configs(COMMON_PATH, BIOMARKERS_PATH)
    rng = random.Random(seed)
    for i in range(n):
        profile = sample_patient_profile(f"p_{i:06d}", rng)
        gene = rng.choice(PANEL_BIOMARKERS)
        _, matched_key = resolve_status_distribution(
            biomarkers.get(gene), profile
        )
        complexity_level = "L2" if rng.random() < l2_fraction else "L1"
        rendered = render_l1_record(
            gene, profile, biomarkers, common, rng,
            complexity_level=complexity_level,
        )
        post = apply_technical_noise(rendered, common, rng)
        yield profile, gene, matched_key, post


def _record_to_dict(
    profile: PatientProfile,
    gene: str,
    matched_key: str,
    record: PostProcessedRecord,
    idx: int,
) -> dict:
    labeled_spans = [
        {"span_type": name, "text": record.sentence[s:e], "char_span": [s, e]}
        for name, (s, e) in sorted(record.spans.items(), key=lambda kv: kv[1])
    ]
    return {
        "record_id": f"{record.complexity_level.lower()}_{idx:06d}",
        "record_type": record.complexity_level.lower(),
        "complexity_level": record.complexity_level,
        "sentence": record.sentence,
        "patient_profile": asdict(profile),
        "assertion": {
            "gene": gene,
            "variant_id": record.variant_id,
            "negative_form_id": record.negative_form_id,
            "clone_id": record.clone_id,
            "status": record.status,
            "test_method": record.test_method,
            "measurement_value": record.measurement_value,
            "frame_template": record.frame_template,
            "matched_population_key": matched_key,
        },
        "labeled_spans": labeled_spans,
        "post_process": record.applied_transforms,
    }


def _two_sigma_within(
    observed_rate: float, configured_rate: float, n: int
) -> bool:
    """Return True when observed is within 2σ of configured.

    σ is the binomial standard error for the configured rate. For tiny rates
    + small n the bar is loose; for 10K records and rates of 0.1+, ~±0.006.
    """
    if n == 0:
        return True
    sigma = math.sqrt(configured_rate * (1 - configured_rate) / n)
    if sigma == 0:
        return abs(observed_rate - configured_rate) < 1e-9
    return abs(observed_rate - configured_rate) <= 2 * sigma


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=DEFAULT_N)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument(
        "--l2-fraction",
        type=float,
        default=DEFAULT_L2_FRACTION,
        help="Fraction of records to render with L2 shorthand frames",
    )
    args = parser.parse_args()
    if not 0.0 <= args.l2_fraction <= 1.0:
        parser.error("--l2-fraction must be in [0, 1]")

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    corpus_path = output_dir / "corpus.jsonl"
    report_path = output_dir / "prevalence_report.json"

    status_counts: dict[tuple[str, str], Counter] = defaultdict(Counter)
    totals_by_pair: Counter = Counter()
    complexity_counts: Counter = Counter()

    sentences_seen: Counter = Counter()
    span_violations = 0

    with open(corpus_path, "w", encoding="utf-8") as f:
        for idx, (profile, gene, matched_key, record) in enumerate(
            _iter_records(args.n, args.seed, args.l2_fraction)
        ):
            f.write(
                json.dumps(_record_to_dict(profile, gene, matched_key, record, idx))
                + "\n"
            )

            for name, (s, e) in record.spans.items():
                if record.sentence[s:e] == "":
                    span_violations += 1

            pair_key = (gene, matched_key)
            status_counts[pair_key][record.status] += 1
            totals_by_pair[pair_key] += 1
            sentences_seen[record.sentence] += 1
            complexity_counts[record.complexity_level] += 1

    common, biomarkers = load_configs(COMMON_PATH, BIOMARKERS_PATH)

    per_pair_report: list[dict] = []
    within_2sigma = 0
    total_checks = 0
    for (gene, histology), total in sorted(totals_by_pair.items()):
        entry = biomarkers.get(gene)
        configured = entry.status_distribution_by_population.get(histology)
        if configured is None:
            continue
        counts = status_counts[(gene, histology)]
        for status_name in ("positive", "negative", "equivocal", "not_tested"):
            observed_count = counts[status_name]
            observed_rate = observed_count / total
            configured_rate = getattr(configured, status_name)
            in_band = _two_sigma_within(observed_rate, configured_rate, total)
            per_pair_report.append(
                {
                    "gene": gene,
                    "histology": histology,
                    "status": status_name,
                    "n": total,
                    "observed_count": observed_count,
                    "observed_rate": round(observed_rate, 4),
                    "configured_rate": round(configured_rate, 4),
                    "within_2sigma": in_band,
                }
            )
            total_checks += 1
            if in_band:
                within_2sigma += 1

    max_dupe_sentence, max_dupe_count = sentences_seen.most_common(1)[0]
    report = {
        "n_records": args.n,
        "seed": args.seed,
        "l2_fraction_requested": args.l2_fraction,
        "l2_fraction_observed": round(
            complexity_counts.get("L2", 0) / max(args.n, 1), 4
        ),
        "complexity_counts": dict(complexity_counts),
        "panel_biomarkers": list(PANEL_BIOMARKERS),
        "within_2sigma_count": within_2sigma,
        "total_checks": total_checks,
        "within_2sigma_ratio": round(within_2sigma / max(total_checks, 1), 4),
        "span_violations": span_violations,
        "max_duplicate_sentence": max_dupe_sentence,
        "max_duplicate_count": max_dupe_count,
        "per_pair": per_pair_report,
    }
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(f"wrote {args.n} records to {corpus_path}")
    print(f"wrote prevalence report to {report_path}")
    print(
        f"within-2σ: {within_2sigma}/{total_checks} "
        f"({100 * within_2sigma / max(total_checks, 1):.1f}%)"
    )
    print(f"span violations (empty slice): {span_violations}")
    print(
        f"max identical sentence count: {max_dupe_count} "
        f"({100 * max_dupe_count / args.n:.2f}% of corpus)"
    )


if __name__ == "__main__":
    main()
