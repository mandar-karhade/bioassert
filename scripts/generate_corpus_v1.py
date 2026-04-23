"""Generate the v1 corpus: mixed L1–L7 records from the JSON-driven pipeline.

Deterministic (seeded). Output: ``datasets/v1_phase3.10/corpus.jsonl`` plus a
``prevalence_report.json`` summarizing observed vs configured status rates
per biomarker.

Phase 3 scope covers the full Tier 1 panel: mutation/fusion/composite
biomarkers plus the two expression biomarkers (PD-L1, TMB) with
``negative_forms`` and IHC ``clone_attribution`` rendering. Status
distributions are flat per biomarker for the cohort described by
``biomarkers.json`` — to target a different cohort, author a separate
config and mix outputs externally.

Re-run with different ``--seed`` to regenerate.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator

from bioassert.config import load_configs
from bioassert.generator.post_process import (
    PostProcessedRecord,
    apply_technical_noise,
)
from bioassert.generator.renderer import (
    render_l1_record,
    render_l3_record,
    render_l4_record,
    render_l5_record,
    render_l6_record,
    render_l7_record,
)

ROOT = Path(__file__).resolve().parents[1]
COMMON_PATH = ROOT / "bioconfigs" / "common_variations.json"
BIOMARKERS_PATH = ROOT / "bioconfigs" / "biomarkers.json"
# Keep generated corpora separated by sub-phase so prior snapshots survive
# when we regen for a new sub-phase. Override via --output-dir for ad-hoc runs.
OUTPUT_DIR = ROOT / "datasets" / "v1_phase3.10"

MUTATION_BIOMARKERS: tuple[str, ...] = (
    "EGFR",
    "KRAS",
    "BRAF",
    "ALK",
    "ROS1",
    "RET",
    "NTRK",
    "MET",
    "ERBB2",
)

EXPRESSION_BIOMARKERS: tuple[str, ...] = ("PD-L1", "TMB")

PANEL_BIOMARKERS: tuple[str, ...] = MUTATION_BIOMARKERS + EXPRESSION_BIOMARKERS

# Probability that an L3 record draws from the expression class rather than
# the mutation class. Expression has only 2 biomarkers so L3 N is always 2;
# mutation has 9, so N ∈ {2,3,4}. Biased toward mutation to keep variety.
L3_EXPRESSION_CLASS_PROB = 0.15

DEFAULT_N = 50_000
DEFAULT_SEED = 42
DEFAULT_L2_FRACTION = 0.3
DEFAULT_L3_FRACTION = 0.0
DEFAULT_L3S_FRACTION = 0.0
DEFAULT_L4_FRACTION = 0.0
DEFAULT_L4S_FRACTION = 0.0
DEFAULT_L5_FRACTION = 0.0
DEFAULT_L6_FRACTION = 0.0
DEFAULT_L7_FRACTION = 0.0
# Compounding tier split (applies only to L3/L3S/L4/L4S/L5 records). The
# split is a binary low/high per user feedback 2026-04-23; medium was
# collapsed into high. Defaults to 0.5/0.5 and must sum to 1.0.
DEFAULT_COMPOUND_LOW = 0.5
DEFAULT_COMPOUND_HIGH = 0.5

_COMPOUND_LEVELS: frozenset[str] = frozenset({"L3", "L3S", "L4", "L4S", "L5"})


def _sample_tier(rng: random.Random, low: float) -> str:
    """Sample a compounding tier. Binary low/high split; low = exactly 2
    genes, high = 3+ (clamped to pool size).
    """
    return "low" if rng.random() < low else "high"


def _compound_pool(
    tier: str, rng: random.Random
) -> list[str]:
    """Pick the biomarker pool for a compound record.

    High tier may mix mutation and expression classes (per user direction
    "high use almost all extra"); low keeps the pool homogeneous so the
    surface stays clinically coherent for a 2-gene assertion.
    """
    if tier == "high":
        return list(PANEL_BIOMARKERS)
    pool_is_expression = rng.random() < L3_EXPRESSION_CLASS_PROB
    return list(EXPRESSION_BIOMARKERS if pool_is_expression else MUTATION_BIOMARKERS)


def _iter_records(
    n: int,
    seed: int,
    l2_fraction: float,
    l3_fraction: float,
    l3s_fraction: float,
    l4_fraction: float,
    l4s_fraction: float,
    l5_fraction: float,
    l6_fraction: float,
    l7_fraction: float,
    compound_low: float,
    compound_high: float,
) -> Iterator[tuple[str, PostProcessedRecord]]:
    common, biomarkers = load_configs(COMMON_PATH, BIOMARKERS_PATH)
    rng = random.Random(seed)
    for _ in range(n):
        roll = rng.random()
        l3_prose_bound = l3_fraction
        l3s_bound = l3_fraction + l3s_fraction
        l4_bound = l3s_bound + l4_fraction
        l4s_bound = l4_bound + l4s_fraction
        l5_bound = l4s_bound + l5_fraction
        l6_bound = l5_bound + l6_fraction
        l7_bound = l6_bound + l7_fraction
        if roll < l3s_bound:
            tier = _sample_tier(rng, compound_low)
            pool = _compound_pool(tier, rng)
            level = "L3" if roll < l3_prose_bound else "L3S"
            rendered = render_l3_record(
                pool, biomarkers, common, rng,
                complexity_level=level,
                compounding_tier=tier,
            )
            post = apply_technical_noise(rendered, common, biomarkers, rng)
            yield rendered.assertions[0].gene, post
            continue
        if roll < l4s_bound:
            tier = _sample_tier(rng, compound_low)
            pool = _compound_pool(tier, rng)
            level = "L4" if roll < l4_bound else "L4S"
            rendered = render_l4_record(
                pool, biomarkers, common, rng,
                complexity_level=level,
                compounding_tier=tier,
            )
            post = apply_technical_noise(rendered, common, biomarkers, rng)
            yield rendered.assertions[0].gene, post
            continue
        if roll < l5_bound:
            tier = _sample_tier(rng, compound_low)
            pool = _compound_pool(tier, rng)
            rendered = render_l5_record(
                pool, biomarkers, common, rng,
                compounding_tier=tier,
            )
            post = apply_technical_noise(rendered, common, biomarkers, rng)
            yield rendered.assertions[0].gene, post
            continue
        if roll < l6_bound:
            # L6 is single-gene regardless of shape; pool is the full panel
            # so any biomarker can surface with temporal/certainty qualifiers.
            rendered = render_l6_record(
                list(PANEL_BIOMARKERS), biomarkers, common, rng
            )
            post = apply_technical_noise(rendered, common, biomarkers, rng)
            yield rendered.assertions[0].gene, post
            continue
        if roll < l7_bound:
            # L7 is single-fact multi-sentence regardless of shape; pool is
            # the full panel so any biomarker can surface in the claim.
            rendered = render_l7_record(
                list(PANEL_BIOMARKERS), biomarkers, common, rng
            )
            post = apply_technical_noise(rendered, common, biomarkers, rng)
            yield rendered.assertions[0].gene, post
            continue
        gene = rng.choice(PANEL_BIOMARKERS)
        complexity_level = "L2" if roll < l7_bound + l2_fraction else "L1"
        rendered = render_l1_record(
            gene, biomarkers, common, rng,
            complexity_level=complexity_level,
        )
        post = apply_technical_noise(rendered, common, biomarkers, rng)
        yield gene, post


def _record_to_dict(
    record: PostProcessedRecord,
    idx: int,
) -> dict:
    assertions_out: list[dict] = []
    labeled_spans: list[dict] = []
    for fact in record.assertions:
        assertions_out.append(
            {
                "gene": fact.gene,
                "variant_id": fact.variant_id,
                "negative_form_id": fact.negative_form_id,
                "clone_id": fact.clone_id,
                "status": fact.status,
                "test_method": fact.test_method,
                "measurement_value": fact.measurement_value,
                "polarity_scope": fact.polarity_scope,
                "temporal": fact.temporal,
                "certainty": fact.certainty,
                "sentence_index": fact.sentence_index,
                "frame_template": record.frame_template,
            }
        )
        labeled_spans.extend(
            {"span_type": name, "text": record.sentence[s:e], "char_span": [s, e]}
            for name, (s, e) in sorted(fact.spans.items(), key=lambda kv: kv[1])
        )
    return {
        "record_id": f"{record.complexity_level.lower()}_{idx:06d}",
        "record_type": record.complexity_level.lower(),
        "complexity_level": record.complexity_level,
        "compounding_tier": record.compounding_tier,
        "sentence": record.sentence,
        "sentences": list(record.sentences),
        "assertions": assertions_out,
        "labeled_spans": labeled_spans,
        "post_process": record.applied_transforms,
    }


def _classify_l6_shape(assertions: tuple) -> str:
    """Derive the L6 structural shape from the qualifier fields on the facts.

    - ``temporal``: 2 facts, both ``temporal`` set, both ``certainty`` None
    - ``certainty``: 1 fact, ``certainty`` set
    - ``combined``: 2 facts, both qualifiers set on each fact
    """
    if len(assertions) == 1:
        return "certainty"
    if any(a.certainty is not None for a in assertions):
        return "combined"
    return "temporal"


def _classify_l7_shape(record) -> str:
    """Derive the L7 structural shape from ``sentences`` length and the
    claim's ``sentence_index``.

    - 3 sentences → ``setup_claim_qualifier``
    - 2 sentences, claim at index 1 → ``setup_claim``
    - 2 sentences, claim at index 0 → ``claim_anaphora``
    """
    n = len(record.sentences)
    fact = record.assertions[0]
    if n == 3:
        return "setup_claim_qualifier"
    if n == 2 and fact.sentence_index == 1:
        return "setup_claim"
    if n == 2 and fact.sentence_index == 0:
        return "claim_anaphora"
    return "unknown"


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
    parser.add_argument(
        "--l3-fraction",
        type=float,
        default=DEFAULT_L3_FRACTION,
        help="Fraction of records to render as L3 shared-status prose",
    )
    parser.add_argument(
        "--l3s-fraction",
        type=float,
        default=DEFAULT_L3S_FRACTION,
        help="Fraction of records to render as L3 shorthand/tabular",
    )
    parser.add_argument(
        "--l4-fraction",
        type=float,
        default=DEFAULT_L4_FRACTION,
        help="Fraction of records to render as L4 heterogeneous prose",
    )
    parser.add_argument(
        "--l4s-fraction",
        type=float,
        default=DEFAULT_L4S_FRACTION,
        help="Fraction of records to render as L4 shorthand/tabular",
    )
    parser.add_argument(
        "--l5-fraction",
        type=float,
        default=DEFAULT_L5_FRACTION,
        help="Fraction of records to render as L5 negation-scope",
    )
    parser.add_argument(
        "--l6-fraction",
        type=float,
        default=DEFAULT_L6_FRACTION,
        help=(
            "Fraction of records to render as L6 temporal/certainty "
            "(single-gene, shape drawn uniformly from temporal/certainty/"
            "combined)"
        ),
    )
    parser.add_argument(
        "--l7-fraction",
        type=float,
        default=DEFAULT_L7_FRACTION,
        help=(
            "Fraction of records to render as L7 multi-sentence "
            "coreference (shape drawn uniformly from setup_claim / "
            "claim_anaphora / setup_claim_qualifier)"
        ),
    )
    parser.add_argument(
        "--compound-low",
        type=float,
        default=DEFAULT_COMPOUND_LOW,
        help=(
            "Fraction of compound (L3/L3S/L4/L4S/L5) records rendered at the "
            "'low' compounding tier (< 3 genes, i.e., exactly 2)."
        ),
    )
    parser.add_argument(
        "--compound-high",
        type=float,
        default=DEFAULT_COMPOUND_HIGH,
        help=(
            "Fraction of compound records rendered at the 'high' tier "
            "(3+ genes, 3-8 clamped to pool size; cross-class pool mixing "
            "allowed). Must sum with --compound-low to 1.0."
        ),
    )
    args = parser.parse_args()
    for flag, val in (
        ("--l2-fraction", args.l2_fraction),
        ("--l3-fraction", args.l3_fraction),
        ("--l3s-fraction", args.l3s_fraction),
        ("--l4-fraction", args.l4_fraction),
        ("--l4s-fraction", args.l4s_fraction),
        ("--l5-fraction", args.l5_fraction),
        ("--l6-fraction", args.l6_fraction),
        ("--l7-fraction", args.l7_fraction),
        ("--compound-low", args.compound_low),
        ("--compound-high", args.compound_high),
    ):
        if not 0.0 <= val <= 1.0:
            parser.error(f"{flag} must be in [0, 1]")
    total = (
        args.l2_fraction
        + args.l3_fraction
        + args.l3s_fraction
        + args.l4_fraction
        + args.l4s_fraction
        + args.l5_fraction
        + args.l6_fraction
        + args.l7_fraction
    )
    if total > 1.0:
        parser.error(
            "--l2-fraction + --l3-fraction + --l3s-fraction + "
            "--l4-fraction + --l4s-fraction + --l5-fraction + "
            "--l6-fraction + --l7-fraction must be <= 1.0"
        )
    if abs(args.compound_low + args.compound_high - 1.0) > 1e-6:
        parser.error(
            "--compound-low + --compound-high must equal 1.0 "
            "(binary tier split; medium was removed in 3.7)"
        )

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    corpus_path = output_dir / "corpus.jsonl"
    report_path = output_dir / "prevalence_report.json"

    status_counts: dict[str, Counter] = defaultdict(Counter)
    totals_by_gene: Counter = Counter()
    complexity_counts: Counter = Counter()
    # Compound-only: tier counts overall and per complexity level.
    compound_tier_counts: Counter = Counter()
    tier_by_level: dict[str, Counter] = defaultdict(Counter)
    # L6-only: shape counts. Infer shape from the fact qualifiers since the
    # record doesn't carry an explicit shape label (temporal: 2 facts both
    # have ``temporal``, ``certainty`` all None; certainty: 1 fact has
    # ``certainty``; combined: 2 facts, both qualifiers set).
    l6_shape_counts: Counter = Counter()
    # L7-only: shape counts. Infer shape from ``sentences`` length and the
    # claim's ``sentence_index`` (the only fact's index).
    l7_shape_counts: Counter = Counter()

    sentences_seen: Counter = Counter()
    span_violations = 0

    with open(corpus_path, "w", encoding="utf-8") as f:
        for idx, (gene, record) in enumerate(
            _iter_records(
                args.n,
                args.seed,
                args.l2_fraction,
                args.l3_fraction,
                args.l3s_fraction,
                args.l4_fraction,
                args.l4s_fraction,
                args.l5_fraction,
                args.l6_fraction,
                args.l7_fraction,
                args.compound_low,
                args.compound_high,
            )
        ):
            f.write(json.dumps(_record_to_dict(record, idx)) + "\n")

            for fact in record.assertions:
                for name, (s, e) in fact.spans.items():
                    if record.sentence[s:e] == "":
                        span_violations += 1

            # Prevalence bookkeeping: only L1/L2 records contribute to the
            # per-gene status counts. L3/L3S draw status once from the first
            # gene's distribution and apply it to every listed gene, so the
            # per-gene empirical rate would be misleading. L4/L4S draw
            # per-gene independently but span multiple genes. L5 forces
            # negative on all wide-scope facts and positive on the exception,
            # bypassing sampling entirely.
            if record.complexity_level in ("L1", "L2"):
                for fact in record.assertions:
                    status_counts[gene][fact.status] += 1
                    totals_by_gene[gene] += 1
            sentences_seen[record.sentence] += 1
            complexity_counts[record.complexity_level] += 1
            if record.complexity_level in _COMPOUND_LEVELS:
                compound_tier_counts[record.compounding_tier] += 1
                tier_by_level[record.complexity_level][record.compounding_tier] += 1
            if record.complexity_level == "L6":
                l6_shape_counts[_classify_l6_shape(record.assertions)] += 1
            if record.complexity_level == "L7":
                l7_shape_counts[_classify_l7_shape(record)] += 1

    common, biomarkers = load_configs(COMMON_PATH, BIOMARKERS_PATH)

    per_gene_report: list[dict] = []
    within_2sigma = 0
    total_checks = 0
    for gene, total in sorted(totals_by_gene.items()):
        configured = biomarkers.get(gene).status_distribution
        counts = status_counts[gene]
        for status_name in ("positive", "negative", "equivocal", "not_tested"):
            observed_count = counts[status_name]
            observed_rate = observed_count / total
            configured_rate = getattr(configured, status_name)
            in_band = _two_sigma_within(observed_rate, configured_rate, total)
            per_gene_report.append(
                {
                    "gene": gene,
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
    compound_total = sum(compound_tier_counts.values())

    def _tier_fraction(tier: str) -> float:
        if compound_total == 0:
            return 0.0
        return round(compound_tier_counts.get(tier, 0) / compound_total, 4)

    tier_breakdown_by_level: dict[str, dict[str, int]] = {
        level: dict(counts) for level, counts in sorted(tier_by_level.items())
    }
    report = {
        "n_records": args.n,
        "seed": args.seed,
        "l2_fraction_requested": args.l2_fraction,
        "l2_fraction_observed": round(
            complexity_counts.get("L2", 0) / max(args.n, 1), 4
        ),
        "l3_fraction_requested": args.l3_fraction,
        "l3_fraction_observed": round(
            complexity_counts.get("L3", 0) / max(args.n, 1), 4
        ),
        "l3s_fraction_requested": args.l3s_fraction,
        "l3s_fraction_observed": round(
            complexity_counts.get("L3S", 0) / max(args.n, 1), 4
        ),
        "l4_fraction_requested": args.l4_fraction,
        "l4_fraction_observed": round(
            complexity_counts.get("L4", 0) / max(args.n, 1), 4
        ),
        "l4s_fraction_requested": args.l4s_fraction,
        "l4s_fraction_observed": round(
            complexity_counts.get("L4S", 0) / max(args.n, 1), 4
        ),
        "l5_fraction_requested": args.l5_fraction,
        "l5_fraction_observed": round(
            complexity_counts.get("L5", 0) / max(args.n, 1), 4
        ),
        "l6_fraction_requested": args.l6_fraction,
        "l6_fraction_observed": round(
            complexity_counts.get("L6", 0) / max(args.n, 1), 4
        ),
        "l6_shape_counts": dict(l6_shape_counts),
        "l7_fraction_requested": args.l7_fraction,
        "l7_fraction_observed": round(
            complexity_counts.get("L7", 0) / max(args.n, 1), 4
        ),
        "l7_shape_counts": dict(l7_shape_counts),
        "complexity_counts": dict(complexity_counts),
        "compound_low_fraction_requested": args.compound_low,
        "compound_low_fraction_observed": _tier_fraction("low"),
        "compound_high_fraction_requested": args.compound_high,
        "compound_high_fraction_observed": _tier_fraction("high"),
        "compound_total": compound_total,
        "compound_tier_counts": dict(compound_tier_counts),
        "compound_tier_by_level": tier_breakdown_by_level,
        "panel_biomarkers": list(PANEL_BIOMARKERS),
        "within_2sigma_count": within_2sigma,
        "total_checks": total_checks,
        "within_2sigma_ratio": round(within_2sigma / max(total_checks, 1), 4),
        "span_violations": span_violations,
        "max_duplicate_sentence": max_dupe_sentence,
        "max_duplicate_count": max_dupe_count,
        "per_gene": per_gene_report,
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
