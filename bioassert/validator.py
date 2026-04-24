"""Distributional validation for a completed generation run.

Takes a run directory, reads the snapshotted configs alongside ``corpus.jsonl``,
and checks every bucket of every probability-weighted distribution in the
config against what was actually observed. Emits ``validation_report.json``
with per-category ``within_2σ`` pass/fail flags and z-scores.

Categories checked (each against its configured distribution):

* **complexity_level** — L1/L2/L3/L3S/L4/L4S/L5/L6/L7 fractions (CLI args).
* **compounding_tier** — low vs. high, on compound-eligible levels.
* **status.<GENE>** — positive/negative/equivocal/not_tested per biomarker
  (scoped to L1/L2 records, matching the sampler's decision surface).
* **variants.<GENE>** — ``variant_id`` prevalence within each biomarker.
* **negative_forms.<GENE>** — ``negative_form_id`` among negative records of
  biomarkers that have a ``negative_forms`` block.
* **clone_attribution.<GENE>** — attachment rate + per-clone weights for
  biomarkers with a ``clone_attribution`` block.
* **preferred_methods.<GENE>** — ``test_method`` for biomarkers that pin a
  weighted methods override.
* **technical_noise.<SUB>** — whitespace / case / hyphenation / punctuation /
  OCR / PDF / abbreviation_inconsistency mode usage, filtered to records that
  actually applied the transform (``"skipped"`` excluded).

Per-bucket pass uses the normal approximation to the binomial:
``σ = √(p(1-p)/n)``; a bucket passes when
``|observed_rate - configured_rate| ≤ 2σ``. Pass-by-construction cases
(``σ == 0`` with matching rates) are counted as passing.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from bioassert.config import BiomarkerConfig, CommonConfig, load_configs
from bioassert.config.schema import (
    CloneAttribution,
    PostProcessTransformations,
)

L1_L2_LEVELS: frozenset[str] = frozenset({"L1", "L2"})
COMPOUND_LEVELS: frozenset[str] = frozenset({"L3", "L3S", "L4", "L4S", "L5"})
NOISE_SENTINELS: frozenset[str] = frozenset({"skipped"})


class ValidationError(RuntimeError):
    """Raised when the validator cannot locate the inputs it needs."""


def within_2sigma(
    observed_count: int, configured_rate: float, n: int
) -> tuple[bool, float]:
    """Return ``(pass, z)`` for one bucket against its configured rate.

    ``pass`` is the binomial 2σ check. ``z`` is the standardized residual
    (0.0 when ``n == 0`` or when the distribution is degenerate — caller
    can treat 0.0 as "not discriminating").
    """
    if n == 0:
        return True, 0.0
    observed_rate = observed_count / n
    sigma = math.sqrt(configured_rate * (1 - configured_rate) / n)
    if sigma == 0.0:
        return abs(observed_rate - configured_rate) < 1e-9, 0.0
    z = (observed_rate - configured_rate) / sigma
    return abs(z) <= 2.0, z


def _check_distribution(
    category_id: str,
    kind: str,
    description: str,
    configured: dict[str, float],
    observed: Counter,
    n: int,
    *,
    extra_counts_are_failures: bool = True,
) -> dict[str, Any]:
    """Compare ``observed`` counts against ``configured`` rates bucket-by-bucket.

    ``extra_counts_are_failures``: when True, any observed bucket not declared
    in ``configured`` is recorded as a failing bucket with a ``note`` set to
    ``"unexpected_bucket"``. Useful for catching renderer bugs that emit an
    unknown variant_id or method.
    """
    buckets: list[dict[str, Any]] = []
    all_pass = True
    for bucket_name in sorted(configured):
        rate = configured[bucket_name]
        count = observed.get(bucket_name, 0)
        ok, z = within_2sigma(count, rate, n)
        buckets.append(
            {
                "bucket": bucket_name,
                "configured_rate": round(rate, 6),
                "observed_count": count,
                "observed_rate": round(count / n, 6) if n else 0.0,
                "within_2sigma": ok,
                "z": round(z, 3),
            }
        )
        if not ok:
            all_pass = False

    extras: list[str] = [k for k in observed if k not in configured]
    for extra in sorted(extras):
        count = observed[extra]
        if extra_counts_are_failures and count > 0:
            all_pass = False
        buckets.append(
            {
                "bucket": extra,
                "configured_rate": None,
                "observed_count": count,
                "observed_rate": round(count / n, 6) if n else 0.0,
                "within_2sigma": False if (extra_counts_are_failures and count > 0) else True,
                "z": None,
                "note": "unexpected_bucket",
            }
        )

    return {
        "category_id": category_id,
        "kind": kind,
        "description": description,
        "n": n,
        "buckets": buckets,
        "overall_within_2sigma": all_pass and n > 0,
        "skipped": n == 0,
    }


def _observed_complexity_levels(records: list[dict]) -> Counter:
    return Counter(r["complexity_level"] for r in records)


def _observed_compound_tiers(records: list[dict]) -> Counter:
    return Counter(
        r["compounding_tier"] for r in records if r["complexity_level"] in COMPOUND_LEVELS
    )


def _configured_complexity_levels(manifest: dict, n: int) -> dict[str, float]:
    """Reconstruct per-level configured rate from the manifest CLI args.

    Missing fractions (L1) are the residual so the distribution sums to 1.
    """
    cli = manifest.get("cli_args", {})
    configured: dict[str, float] = {
        "L2": float(cli.get("l2_fraction", 0.0)),
        "L3": float(cli.get("l3_fraction", 0.0)),
        "L3S": float(cli.get("l3s_fraction", 0.0)),
        "L4": float(cli.get("l4_fraction", 0.0)),
        "L4S": float(cli.get("l4s_fraction", 0.0)),
        "L5": float(cli.get("l5_fraction", 0.0)),
        "L6": float(cli.get("l6_fraction", 0.0)),
        "L7": float(cli.get("l7_fraction", 0.0)),
    }
    residual = 1.0 - sum(configured.values())
    configured["L1"] = max(residual, 0.0)
    return configured


def _configured_compound_tiers(manifest: dict) -> dict[str, float]:
    cli = manifest.get("cli_args", {})
    return {
        "low": float(cli.get("compound_low", 0.5)),
        "high": float(cli.get("compound_high", 0.5)),
    }


def _check_status_per_gene(
    records: list[dict], biomarkers: BiomarkerConfig
) -> list[dict[str, Any]]:
    """Per-gene status distribution, scoped to L1/L2 records only.

    L3+ records mix genes and condition on compound status choices, so they'd
    pollute the signal we actually care about (per-gene cohort prior).
    """
    status_by_gene: dict[str, Counter] = {}
    n_by_gene: Counter = Counter()
    for r in records:
        if r["complexity_level"] not in L1_L2_LEVELS:
            continue
        for fact in r["assertions"]:
            gene = fact["gene"]
            status_by_gene.setdefault(gene, Counter())[fact["status"]] += 1
            n_by_gene[gene] += 1

    out: list[dict[str, Any]] = []
    for gene in sorted(n_by_gene):
        bio = biomarkers.get(gene)
        dist = bio.status_distribution
        configured = {
            "positive": dist.positive,
            "negative": dist.negative,
            "equivocal": dist.equivocal,
            "not_tested": dist.not_tested,
        }
        out.append(
            _check_distribution(
                category_id=f"status.{gene}",
                kind="cohort_prior",
                description=f"Status distribution for {gene} (L1/L2 only)",
                configured=configured,
                observed=status_by_gene[gene],
                n=n_by_gene[gene],
            )
        )
    return out


def _check_variants_per_gene(
    records: list[dict], biomarkers: BiomarkerConfig
) -> list[dict[str, Any]]:
    variants_by_gene: dict[str, Counter] = {}
    for r in records:
        for fact in r["assertions"]:
            vid = fact.get("variant_id")
            if vid is None:
                continue
            variants_by_gene.setdefault(fact["gene"], Counter())[vid] += 1

    out: list[dict[str, Any]] = []
    for gene in sorted(variants_by_gene):
        bio = biomarkers.get(gene)
        configured = {
            vid: v.prevalence_within_biomarker for vid, v in bio.variants.items()
        }
        n = sum(variants_by_gene[gene].values())
        out.append(
            _check_distribution(
                category_id=f"variants.{gene}",
                kind="cohort_prior",
                description=f"Variant prevalence within {gene}",
                configured=configured,
                observed=variants_by_gene[gene],
                n=n,
            )
        )
    return out


def _check_negative_forms_per_gene(
    records: list[dict], biomarkers: BiomarkerConfig
) -> list[dict[str, Any]]:
    neg_by_gene: dict[str, Counter] = {}
    for r in records:
        for fact in r["assertions"]:
            nfid = fact.get("negative_form_id")
            if nfid is None:
                continue
            neg_by_gene.setdefault(fact["gene"], Counter())[nfid] += 1

    out: list[dict[str, Any]] = []
    for gene in sorted(neg_by_gene):
        bio = biomarkers.get(gene)
        if bio.negative_forms is None:
            continue
        configured = dict(bio.negative_forms.variations)
        n = sum(neg_by_gene[gene].values())
        out.append(
            _check_distribution(
                category_id=f"negative_forms.{gene}",
                kind="cohort_prior",
                description=f"Negative-form prevalence within {gene}",
                configured=configured,
                observed=neg_by_gene[gene],
                n=n,
            )
        )
    return out


def _check_clone_attribution_per_gene(
    records: list[dict], biomarkers: BiomarkerConfig
) -> list[dict[str, Any]]:
    """Two-stage check: attachment rate, then per-clone weights given attached.

    Scoped to L1/L2 records: the renderer only calls ``maybe_sample_clone``
    from ``render_l1_record`` (which handles both L1 and L2). L3+ facts never
    enter clone sampling, so counting them as unattached would be wrong.
    """
    attached_by_gene: dict[str, Counter] = {}
    eligible_by_gene: Counter = Counter()
    for r in records:
        if r["complexity_level"] not in L1_L2_LEVELS:
            continue
        for fact in r["assertions"]:
            gene = fact["gene"]
            try:
                bio = biomarkers.get(gene)
            except KeyError:
                continue
            if bio.clone_attribution is None:
                continue
            eligible_by_gene[gene] += 1
            cid = fact.get("clone_id")
            if cid is not None:
                attached_by_gene.setdefault(gene, Counter())[cid] += 1

    out: list[dict[str, Any]] = []
    for gene in sorted(eligible_by_gene):
        bio = biomarkers.get(gene)
        clone: CloneAttribution = bio.clone_attribution  # type: ignore[assignment]
        attached_counts = attached_by_gene.get(gene, Counter())
        n_eligible = eligible_by_gene[gene]
        n_attached = sum(attached_counts.values())

        attach_configured = {
            "attached": clone.attachment_probability,
            "unattached": 1.0 - clone.attachment_probability,
        }
        attach_observed = Counter(
            {
                "attached": n_attached,
                "unattached": n_eligible - n_attached,
            }
        )
        out.append(
            _check_distribution(
                category_id=f"clone_attribution.{gene}.attachment",
                kind="cohort_prior",
                description=(
                    f"Clone attachment rate for {gene} "
                    f"(Bernoulli p={clone.attachment_probability})"
                ),
                configured=attach_configured,
                observed=attach_observed,
                n=n_eligible,
            )
        )

        if n_attached > 0:
            out.append(
                _check_distribution(
                    category_id=f"clone_attribution.{gene}.clones",
                    kind="cohort_prior",
                    description=f"Per-clone weights within {gene}'s attached records",
                    configured=dict(clone.variations),
                    observed=attached_counts,
                    n=n_attached,
                )
            )
    return out


def _check_preferred_methods_per_gene(
    records: list[dict], biomarkers: BiomarkerConfig
) -> list[dict[str, Any]]:
    """Per-gene test-method distribution, conditional on a method being attached.

    The renderer ( ``renderer._maybe_sample_method`` ) collapses two independent
    decisions into ``fact.test_method=None``:

    1. The attachment gate (``rng.random() >= method_attach_prob``) — fully
       drops the method slot for this record.
    2. Sampling the ``unspecified`` bucket from ``preferred_methods`` — also
       stored as ``None``.

    We can't disentangle these from the corpus alone, so we check the
    distribution *conditional on method being attached*: we drop ``None``
    observations, drop the ``unspecified`` bucket from the configured
    distribution, and renormalize the remaining weights to sum to 1. That's
    the distribution we'd expect among facts that actually name a method.
    """
    methods_by_gene: dict[str, Counter] = {}
    for r in records:
        for fact in r["assertions"]:
            method = fact.get("test_method")
            if method is None:
                continue
            gene = fact["gene"]
            try:
                bio = biomarkers.get(gene)
            except KeyError:
                continue
            if bio.preferred_methods is None:
                continue
            methods_by_gene.setdefault(gene, Counter())[method] += 1

    out: list[dict[str, Any]] = []
    for gene in sorted(methods_by_gene):
        bio = biomarkers.get(gene)
        variations = dict(bio.preferred_methods.variations)  # type: ignore[union-attr]
        non_unspecified = {k: v for k, v in variations.items() if k != "unspecified"}
        total = sum(non_unspecified.values())
        if total <= 0:
            continue
        configured = {k: v / total for k, v in non_unspecified.items()}
        n = sum(methods_by_gene[gene].values())
        out.append(
            _check_distribution(
                category_id=f"preferred_methods.{gene}",
                kind="cohort_prior",
                description=(
                    f"Test-method distribution for {gene} "
                    "(conditional on method attached; 'unspecified' bucket "
                    "collapses to None in output and is renormalized out)"
                ),
                configured=configured,
                observed=methods_by_gene[gene],
                n=n,
            )
        )
    return out


def _check_technical_noise(
    records: list[dict], common: CommonConfig
) -> list[dict[str, Any]]:
    noise = common.categories.get("technical_noise")
    if not isinstance(noise, PostProcessTransformations):
        return []

    per_subcat: dict[str, Counter] = {sub: Counter() for sub in noise.categories}
    for r in records:
        applied = r.get("post_process", {}) or {}
        for sub, mode in applied.items():
            if sub not in per_subcat:
                continue
            if mode in NOISE_SENTINELS:
                continue
            per_subcat[sub][mode] += 1

    out: list[dict[str, Any]] = []
    for sub in sorted(per_subcat):
        configured = dict(noise.categories[sub].distribution)
        observed = per_subcat[sub]
        n = sum(observed.values())
        out.append(
            _check_distribution(
                category_id=f"technical_noise.{sub}",
                kind="surface_noise",
                description=f"technical_noise.{sub} mode distribution (skipped records excluded)",
                configured=configured,
                observed=observed,
                n=n,
            )
        )
    return out


def _load_snapshot(run_dir: Path) -> tuple[CommonConfig, BiomarkerConfig]:
    snapshot_common = run_dir / "snapshot" / "configs" / "common_variations.json"
    snapshot_biomarkers = run_dir / "snapshot" / "configs" / "biomarkers.json"
    if not snapshot_common.is_file() or not snapshot_biomarkers.is_file():
        raise ValidationError(
            f"snapshot configs not found under {run_dir}/snapshot/configs/"
        )
    return load_configs(snapshot_common, snapshot_biomarkers)


def _load_records(run_dir: Path) -> list[dict]:
    corpus = run_dir / "corpus.jsonl"
    if not corpus.is_file():
        raise ValidationError(f"corpus.jsonl not found at {corpus}")
    records: list[dict] = []
    with corpus.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def _load_manifest(run_dir: Path) -> dict:
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.is_file():
        raise ValidationError(f"manifest.json not found at {manifest_path}")
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def _summarize(categories: list[dict]) -> dict[str, Any]:
    buckets_total = 0
    buckets_pass = 0
    failures: list[dict[str, Any]] = []
    skipped: list[str] = []
    for cat in categories:
        if cat["skipped"]:
            skipped.append(cat["category_id"])
            continue
        failed_buckets = [
            b["bucket"] for b in cat["buckets"] if not b["within_2sigma"]
        ]
        for b in cat["buckets"]:
            buckets_total += 1
            if b["within_2sigma"]:
                buckets_pass += 1
        if failed_buckets:
            failures.append(
                {
                    "category_id": cat["category_id"],
                    "n": cat["n"],
                    "failed_buckets": failed_buckets,
                }
            )
    return {
        "categories_checked": len([c for c in categories if not c["skipped"]]),
        "categories_skipped": skipped,
        "buckets_checked": buckets_total,
        "within_2sigma_pass": buckets_pass,
        "within_2sigma_fail": buckets_total - buckets_pass,
        "pass_rate": round(buckets_pass / buckets_total, 4) if buckets_total else 1.0,
        "failing_categories": failures,
    }


def validate_run(run_dir: Path) -> dict[str, Any]:
    """Build and return the full validation report for ``run_dir``.

    Does not write to disk — caller is responsible for persisting the return
    value. Reads the run's snapshotted configs so the check is authoritative
    against the exact configs that produced the corpus (not whatever's on
    disk today).
    """
    run_dir = Path(run_dir).resolve()
    if not run_dir.is_dir():
        raise ValidationError(f"run directory does not exist: {run_dir}")

    manifest = _load_manifest(run_dir)
    common, biomarkers = _load_snapshot(run_dir)
    records = _load_records(run_dir)
    n = len(records)

    categories: list[dict[str, Any]] = []

    categories.append(
        _check_distribution(
            category_id="complexity_level",
            kind="cli_sampling",
            description="Complexity level fractions (L1 = residual of --l{2..7}-fraction)",
            configured=_configured_complexity_levels(manifest, n),
            observed=_observed_complexity_levels(records),
            n=n,
        )
    )

    compound_total = sum(
        1 for r in records if r["complexity_level"] in COMPOUND_LEVELS
    )
    categories.append(
        _check_distribution(
            category_id="compounding_tier",
            kind="cli_sampling",
            description="Compound-tier split on L3/L3S/L4/L4S/L5 records only",
            configured=_configured_compound_tiers(manifest),
            observed=_observed_compound_tiers(records),
            n=compound_total,
        )
    )

    categories.extend(_check_status_per_gene(records, biomarkers))
    categories.extend(_check_variants_per_gene(records, biomarkers))
    categories.extend(_check_negative_forms_per_gene(records, biomarkers))
    categories.extend(_check_clone_attribution_per_gene(records, biomarkers))
    categories.extend(_check_preferred_methods_per_gene(records, biomarkers))
    categories.extend(_check_technical_noise(records, common))

    return {
        "run_id": manifest.get("run_id", run_dir.name),
        "project": manifest.get("project"),
        "schema_type": manifest.get("schema_type"),
        "n_records": n,
        "categories": categories,
        "summary": _summarize(categories),
    }


def write_validation_report(run_dir: Path) -> Path:
    """Validate ``run_dir`` and write ``validation_report.json`` into it."""
    report = validate_run(run_dir)
    out = run_dir / "validation_report.json"
    with out.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    return out


__all__ = [
    "ValidationError",
    "validate_run",
    "within_2sigma",
    "write_validation_report",
]
