"""``bioassert`` command-line entry point.

Installed via ``pyproject.toml`` ``[project.scripts]`` as ``bioassert``. After
``uv pip install -e .`` (or ``pip install -e .``) the command is available
from anywhere, and a project directory can live outside the source tree::

    bioassert generate --project /path/to/my_project --n 10000 --tag smoke

Each invocation writes into ``<project>/outputs/run_NNN_<tag?>_<UTC>/`` with
a snapshotted copy of the configs and a ``manifest.json`` recording the CLI
args, seed, timestamp, package version, and (best-effort) git SHA.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Iterator

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
from bioassert.project import Project, ProjectError
from bioassert.validator import ValidationError, validate_run, write_validation_report

L3_EXPRESSION_CLASS_PROB = 0.15
_COMPOUND_LEVELS: frozenset[str] = frozenset({"L3", "L3S", "L4", "L4S", "L5"})


def _split_panel(project: Project) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Derive the mutation/expression panel lists from the project's biomarkers.

    Panels are ordered by the biomarker insertion order in ``biomarkers.json``,
    so each project controls its own canonical gene list.
    """
    mutations: list[str] = []
    expressions: list[str] = []
    for name, biomarker in project.biomarkers.biomarkers.items():
        if biomarker.alteration_type == "mutation":
            mutations.append(name)
        elif biomarker.alteration_type == "expression":
            expressions.append(name)
    if not mutations and not expressions:
        raise ProjectError(
            f"project {project.name!r} has no biomarkers"
        )
    return tuple(mutations), tuple(expressions)

DEFAULT_N = 50_000
DEFAULT_SEED = 42
DEFAULT_L2_FRACTION = 0.3
DEFAULT_COMPOUND_LOW = 0.5
DEFAULT_COMPOUND_HIGH = 0.5


def _sample_tier(rng: random.Random, low: float) -> str:
    return "low" if rng.random() < low else "high"


def _compound_pool(
    tier: str,
    rng: random.Random,
    mutation_panel: tuple[str, ...],
    expression_panel: tuple[str, ...],
    full_panel: tuple[str, ...],
) -> list[str]:
    if tier == "high":
        return list(full_panel)
    if not expression_panel:
        return list(mutation_panel)
    if not mutation_panel:
        return list(expression_panel)
    pool_is_expression = rng.random() < L3_EXPRESSION_CLASS_PROB
    return list(expression_panel if pool_is_expression else mutation_panel)


def _iter_records(
    project: Project,
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
) -> Iterator[tuple[str, PostProcessedRecord]]:
    common = project.common
    biomarkers = project.biomarkers
    mutation_panel, expression_panel = _split_panel(project)
    full_panel = mutation_panel + expression_panel
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
            pool = _compound_pool(tier, rng, mutation_panel, expression_panel, full_panel)
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
            pool = _compound_pool(tier, rng, mutation_panel, expression_panel, full_panel)
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
            pool = _compound_pool(tier, rng, mutation_panel, expression_panel, full_panel)
            rendered = render_l5_record(
                pool, biomarkers, common, rng,
                compounding_tier=tier,
            )
            post = apply_technical_noise(rendered, common, biomarkers, rng)
            yield rendered.assertions[0].gene, post
            continue
        if roll < l6_bound:
            rendered = render_l6_record(
                list(full_panel), biomarkers, common, rng
            )
            post = apply_technical_noise(rendered, common, biomarkers, rng)
            yield rendered.assertions[0].gene, post
            continue
        if roll < l7_bound:
            rendered = render_l7_record(
                list(full_panel), biomarkers, common, rng
            )
            post = apply_technical_noise(rendered, common, biomarkers, rng)
            yield rendered.assertions[0].gene, post
            continue
        gene = rng.choice(full_panel)
        complexity_level = "L2" if roll < l7_bound + l2_fraction else "L1"
        rendered = render_l1_record(
            gene, biomarkers, common, rng,
            complexity_level=complexity_level,
        )
        post = apply_technical_noise(rendered, common, biomarkers, rng)
        yield gene, post


def _record_to_dict(record: PostProcessedRecord, idx: int) -> dict:
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
    if len(assertions) == 1:
        return "certainty"
    if any(a.certainty is not None for a in assertions):
        return "combined"
    return "temporal"


def _classify_l7_shape(record) -> str:
    n = len(record.sentences)
    fact = record.assertions[0]
    if n == 3:
        return "setup_claim_qualifier"
    if n == 2 and fact.sentence_index == 1:
        return "setup_claim"
    if n == 2 and fact.sentence_index == 0:
        return "claim_anaphora"
    return "unknown"


def _two_sigma_within(observed_rate: float, configured_rate: float, n: int) -> bool:
    if n == 0:
        return True
    sigma = math.sqrt(configured_rate * (1 - configured_rate) / n)
    if sigma == 0:
        return abs(observed_rate - configured_rate) < 1e-9
    return abs(observed_rate - configured_rate) <= 2 * sigma


def _bioassert_version() -> str:
    try:
        return version("bioassert")
    except PackageNotFoundError:
        return "0.0.0+unknown"


def _git_sha(project_root: Path) -> str | None:
    """Best-effort ``git rev-parse HEAD`` scoped to the project directory.

    Returns ``None`` if git isn't available, the directory isn't a repo, or
    the command fails for any reason. Never raises.
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=project_root,
            capture_output=True,
            text=True,
            check=False,
            timeout=2,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0:
        return None
    sha = result.stdout.strip()
    return sha or None


def _snapshot_configs(project: Project, run_dir: Path) -> dict[str, str]:
    """Copy configs + project.json into ``run_dir/snapshot/`` byte-for-byte.

    Returns a dict mapping logical name → relative path (for the manifest).
    """
    snapshot_root = run_dir / "snapshot"
    configs_dst = snapshot_root / "configs"
    configs_dst.mkdir(parents=True, exist_ok=True)
    common_dst = configs_dst / project.common_path.name
    biomarkers_dst = configs_dst / project.biomarkers_path.name
    shutil.copy2(project.common_path, common_dst)
    shutil.copy2(project.biomarkers_path, biomarkers_dst)
    project_json_dst = snapshot_root / "project.json"
    shutil.copy2(project.project_json_path, project_json_dst)
    return {
        "common": str(common_dst.relative_to(run_dir)),
        "biomarkers": str(biomarkers_dst.relative_to(run_dir)),
        "project_json": str(project_json_dst.relative_to(run_dir)),
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bioassert")
    subparsers = parser.add_subparsers(dest="command", required=True)

    gen = subparsers.add_parser(
        "generate",
        help="Generate a corpus for a project.",
    )
    gen.add_argument(
        "--project",
        type=Path,
        required=True,
        help="Path to the project directory (containing project.json).",
    )
    gen.add_argument("--n", type=int, default=DEFAULT_N)
    gen.add_argument("--seed", type=int, default=DEFAULT_SEED)
    gen.add_argument(
        "--tag",
        type=str,
        default=None,
        help="Short label for the run folder; must match [A-Za-z0-9._-]+.",
    )
    for flag, default in (
        ("--l2-fraction", DEFAULT_L2_FRACTION),
        ("--l3-fraction", 0.0),
        ("--l3s-fraction", 0.0),
        ("--l4-fraction", 0.0),
        ("--l4s-fraction", 0.0),
        ("--l5-fraction", 0.0),
        ("--l6-fraction", 0.0),
        ("--l7-fraction", 0.0),
    ):
        gen.add_argument(flag, type=float, default=default)
    gen.add_argument(
        "--compound-low", type=float, default=DEFAULT_COMPOUND_LOW,
        help="Fraction of compound records rendered at the 'low' tier (exactly 2 genes).",
    )
    gen.add_argument(
        "--compound-high", type=float, default=DEFAULT_COMPOUND_HIGH,
        help="Fraction of compound records at the 'high' tier (3+ genes).",
    )
    gen.add_argument(
        "--skip-validation",
        action="store_true",
        help="Skip the post-generation validation_report.json step.",
    )

    val = subparsers.add_parser(
        "validate",
        help="Validate a completed run directory against its snapshotted configs.",
    )
    val.add_argument(
        "run_dir",
        type=Path,
        help="Path to a run directory "
             "(e.g., projects/<name>/outputs/run_001_<tag>_<UTC>/).",
    )
    return parser


def _generate(args: argparse.Namespace) -> int:
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
            print(f"error: {flag} must be in [0, 1] (got {val})", file=sys.stderr)
            return 2
    total = (
        args.l2_fraction + args.l3_fraction + args.l3s_fraction
        + args.l4_fraction + args.l4s_fraction + args.l5_fraction
        + args.l6_fraction + args.l7_fraction
    )
    if total > 1.0:
        print(
            "error: sum of --l2..l7 fractions must be <= 1.0",
            file=sys.stderr,
        )
        return 2
    if abs(args.compound_low + args.compound_high - 1.0) > 1e-6:
        print(
            "error: --compound-low + --compound-high must equal 1.0",
            file=sys.stderr,
        )
        return 2

    try:
        project = Project.load(args.project)
    except ProjectError as exc:
        print(f"error loading project: {exc}", file=sys.stderr)
        return 2

    started_at = datetime.now(timezone.utc)
    run_dir = project.next_run_dir(tag=args.tag, now=started_at)
    run_dir.mkdir(parents=True, exist_ok=True)
    corpus_path = run_dir / "corpus.jsonl"
    report_path = run_dir / "prevalence_report.json"
    manifest_path = run_dir / "manifest.json"

    snapshot_paths = _snapshot_configs(project, run_dir)

    status_counts: dict[str, Counter] = defaultdict(Counter)
    totals_by_gene: Counter = Counter()
    complexity_counts: Counter = Counter()
    compound_tier_counts: Counter = Counter()
    tier_by_level: dict[str, Counter] = defaultdict(Counter)
    l6_shape_counts: Counter = Counter()
    l7_shape_counts: Counter = Counter()
    sentences_seen: Counter = Counter()
    span_violations = 0

    with open(corpus_path, "w", encoding="utf-8") as f:
        for idx, (gene, record) in enumerate(
            _iter_records(
                project, args.n, args.seed,
                args.l2_fraction, args.l3_fraction, args.l3s_fraction,
                args.l4_fraction, args.l4s_fraction, args.l5_fraction,
                args.l6_fraction, args.l7_fraction,
                args.compound_low,
            )
        ):
            f.write(json.dumps(_record_to_dict(record, idx)) + "\n")
            for fact in record.assertions:
                for name, (s, e) in fact.spans.items():
                    if record.sentence[s:e] == "":
                        span_violations += 1
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

    per_gene_report: list[dict] = []
    within_2sigma = 0
    total_checks = 0
    for gene, total_g in sorted(totals_by_gene.items()):
        configured = project.biomarkers.get(gene).status_distribution
        counts = status_counts[gene]
        for status_name in ("positive", "negative", "equivocal", "not_tested"):
            observed_count = counts[status_name]
            observed_rate = observed_count / total_g
            configured_rate = getattr(configured, status_name)
            in_band = _two_sigma_within(observed_rate, configured_rate, total_g)
            per_gene_report.append(
                {
                    "gene": gene,
                    "status": status_name,
                    "n": total_g,
                    "observed_count": observed_count,
                    "observed_rate": round(observed_rate, 4),
                    "configured_rate": round(configured_rate, 4),
                    "within_2sigma": in_band,
                }
            )
            total_checks += 1
            if in_band:
                within_2sigma += 1

    max_dupe_sentence, max_dupe_count = (
        sentences_seen.most_common(1)[0] if sentences_seen else ("", 0)
    )
    compound_total = sum(compound_tier_counts.values())

    def _tier_fraction(tier: str) -> float:
        if compound_total == 0:
            return 0.0
        return round(compound_tier_counts.get(tier, 0) / compound_total, 4)

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
        "compound_tier_by_level": {
            level: dict(counts) for level, counts in sorted(tier_by_level.items())
        },
        "panel_biomarkers": list(project.biomarkers.biomarkers.keys()),
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

    cli_args = {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}
    cli_args.pop("command", None)
    manifest = {
        "run_id": run_dir.name,
        "project": project.name,
        "schema_type": project.schema_type,
        "timestamp_utc": started_at.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "tag": args.tag,
        "seed": args.seed,
        "n_records": args.n,
        "cli_args": cli_args,
        "bioassert_version": _bioassert_version(),
        "git_sha": _git_sha(project.root),
        "config_snapshot": snapshot_paths,
    }
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"project: {project.name} ({project.display_name})")
    print(f"run dir: {run_dir}")
    print(f"corpus : {corpus_path} ({args.n} records)")
    print(f"report : {report_path}")
    print(f"within-2σ (status only): {within_2sigma}/{total_checks} "
          f"({100 * within_2sigma / max(total_checks, 1):.1f}%)")
    print(f"span violations: {span_violations}")
    print(f"max identical sentence count: {max_dupe_count} "
          f"({100 * max_dupe_count / max(args.n, 1):.2f}%)")

    if not args.skip_validation:
        vreport_path = write_validation_report(run_dir)
        vreport = json.loads(vreport_path.read_text(encoding="utf-8"))
        summary = vreport["summary"]
        print(
            f"validation: {summary['within_2sigma_pass']}/{summary['buckets_checked']} "
            f"buckets pass ({100 * summary['pass_rate']:.1f}%) across "
            f"{summary['categories_checked']} categories "
            f"({vreport_path.name})"
        )
        if summary["failing_categories"]:
            print("  failing categories:")
            for fail in summary["failing_categories"]:
                print(
                    f"    - {fail['category_id']} "
                    f"(n={fail['n']}, failed={fail['failed_buckets']})"
                )
    return 0


def _validate(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.expanduser().resolve()
    try:
        report_path = write_validation_report(run_dir)
    except ValidationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    report = json.loads(report_path.read_text(encoding="utf-8"))
    summary = report["summary"]
    print(f"run: {report['run_id']}  n_records={report['n_records']}")
    print(
        f"validation: {summary['within_2sigma_pass']}/{summary['buckets_checked']} "
        f"buckets pass ({100 * summary['pass_rate']:.1f}%) across "
        f"{summary['categories_checked']} categories"
    )
    if summary["categories_skipped"]:
        print(f"  skipped (n=0): {len(summary['categories_skipped'])} categories")
    if summary["failing_categories"]:
        print("  failing categories:")
        for fail in summary["failing_categories"]:
            print(
                f"    - {fail['category_id']} "
                f"(n={fail['n']}, failed={fail['failed_buckets']})"
            )
    print(f"wrote {report_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command == "generate":
        return _generate(args)
    if args.command == "validate":
        return _validate(args)
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
