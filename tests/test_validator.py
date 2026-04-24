"""Validator tests: bucket sigma math, category shapes, end-to-end via CLI."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from bioassert.cli import main
from bioassert.validator import (
    ValidationError,
    validate_run,
    within_2sigma,
    write_validation_report,
)

ROOT = Path(__file__).resolve().parents[1]
SHIPPED_PROJECT = ROOT / "projects" / "nsclc_adenocarcinoma"


def _clone_project(tmp_path: Path) -> Path:
    dest = tmp_path / "nsclc"
    shutil.copytree(SHIPPED_PROJECT, dest)
    outputs = dest / "outputs"
    if outputs.exists():
        shutil.rmtree(outputs)
    return dest


def _run(argv: list[str]) -> int:
    return main(argv)


def _latest_run(project_dir: Path) -> Path:
    runs = sorted((project_dir / "outputs").glob("run_*"))
    assert runs
    return runs[-1]


# --- within_2sigma unit tests ---------------------------------------------


def test_within_2sigma_zero_n_passes_by_default() -> None:
    ok, z = within_2sigma(0, 0.5, 0)
    assert ok is True
    assert z == 0.0


def test_within_2sigma_exact_match() -> None:
    ok, z = within_2sigma(50, 0.5, 100)
    assert ok is True
    assert abs(z) < 1e-9


def test_within_2sigma_small_deviation_passes() -> None:
    # p=0.5, n=100 → σ=0.05. Observed 53 (rate 0.53) → z=0.6, passes.
    ok, z = within_2sigma(53, 0.5, 100)
    assert ok is True
    assert 0 < z < 2.0


def test_within_2sigma_large_deviation_fails() -> None:
    # p=0.5, n=100 → σ=0.05. Observed 65 (rate 0.65) → z=3.0, fails.
    ok, z = within_2sigma(65, 0.5, 100)
    assert ok is False
    assert z > 2.0


def test_within_2sigma_degenerate_distribution() -> None:
    # p=0.0, any observation of 1+ fails (sigma=0, not exact match)
    ok, z = within_2sigma(1, 0.0, 10)
    assert ok is False
    assert z == 0.0
    # p=0.0, exact-zero observations pass
    ok, z = within_2sigma(0, 0.0, 10)
    assert ok is True


# --- end-to-end: generate then inspect validation_report.json -------------


@pytest.fixture
def smoke_run(tmp_path: Path) -> Path:
    """Generate a small corpus, return path to the run dir."""
    dest = _clone_project(tmp_path)
    rc = _run(
        [
            "generate",
            "--project", str(dest),
            "--n", "500",
            "--seed", "2024",
            "--tag", "vtest",
        ]
    )
    assert rc == 0
    return _latest_run(dest)


def test_generate_writes_validation_report(smoke_run: Path) -> None:
    report_path = smoke_run / "validation_report.json"
    assert report_path.is_file()
    report = json.loads(report_path.read_text())
    assert report["n_records"] == 500
    assert report["project"] == "nsclc_adenocarcinoma"
    assert report["schema_type"] == "biomarker"
    assert report["run_id"] == smoke_run.name


def test_validation_report_has_complexity_level_category(smoke_run: Path) -> None:
    report = json.loads((smoke_run / "validation_report.json").read_text())
    by_id = {c["category_id"]: c for c in report["categories"]}
    cx = by_id["complexity_level"]
    assert cx["kind"] == "cli_sampling"
    bucket_names = {b["bucket"] for b in cx["buckets"]}
    # Default run has --l2-fraction=0.3; expect L1 + L2 at minimum
    assert "L1" in bucket_names
    assert "L2" in bucket_names


def test_validation_report_has_per_gene_status(smoke_run: Path) -> None:
    report = json.loads((smoke_run / "validation_report.json").read_text())
    ids = {c["category_id"] for c in report["categories"]}
    # Shipped panel must have at least these
    assert "status.EGFR" in ids
    assert "status.KRAS" in ids
    assert "status.PD-L1" in ids


def test_validation_report_has_variants_per_gene(smoke_run: Path) -> None:
    report = json.loads((smoke_run / "validation_report.json").read_text())
    ids = {c["category_id"] for c in report["categories"]}
    assert "variants.EGFR" in ids


def test_validation_report_has_technical_noise_categories(smoke_run: Path) -> None:
    report = json.loads((smoke_run / "validation_report.json").read_text())
    ids = {c["category_id"] for c in report["categories"]}
    # At minimum the shipped config includes these sub-categories
    expected = {
        "technical_noise.whitespace",
        "technical_noise.case_variation",
        "technical_noise.hyphenation_gene_names",
        "technical_noise.punctuation_variation",
    }
    assert expected.issubset(ids), f"missing: {expected - ids}"


def test_validation_report_summary_shape(smoke_run: Path) -> None:
    report = json.loads((smoke_run / "validation_report.json").read_text())
    s = report["summary"]
    assert s["buckets_checked"] > 0
    assert s["within_2sigma_pass"] + s["within_2sigma_fail"] == s["buckets_checked"]
    assert 0.0 <= s["pass_rate"] <= 1.0
    assert s["categories_checked"] >= 10  # shipped panel + noise + sampling


def test_validate_run_matches_written_report(smoke_run: Path) -> None:
    written = json.loads((smoke_run / "validation_report.json").read_text())
    in_memory = validate_run(smoke_run)
    # run_id + n_records + summary must match exactly
    assert written["run_id"] == in_memory["run_id"]
    assert written["n_records"] == in_memory["n_records"]
    assert written["summary"] == in_memory["summary"]


# --- CLI: validate subcommand --------------------------------------------


def test_validate_subcommand_rewrites_report(smoke_run: Path) -> None:
    report_path = smoke_run / "validation_report.json"
    # Delete the auto-written one
    report_path.unlink()
    rc = _run(["validate", str(smoke_run)])
    assert rc == 0
    assert report_path.is_file()


def test_validate_subcommand_errors_on_missing_run(tmp_path: Path) -> None:
    rc = _run(["validate", str(tmp_path / "nope")])
    assert rc == 2


def test_validate_errors_on_missing_snapshot(tmp_path: Path) -> None:
    fake = tmp_path / "run"
    fake.mkdir()
    (fake / "corpus.jsonl").write_text("")
    (fake / "manifest.json").write_text("{}")
    with pytest.raises(ValidationError, match="snapshot configs not found"):
        validate_run(fake)


# --- skip flag -----------------------------------------------------------


def test_skip_validation_flag(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    rc = _run(
        [
            "generate",
            "--project", str(dest),
            "--n", "10",
            "--seed", "1",
            "--skip-validation",
        ]
    )
    assert rc == 0
    run_dir = _latest_run(dest)
    assert not (run_dir / "validation_report.json").exists()


# --- write_validation_report returns path -------------------------------


def test_write_validation_report_returns_output_path(smoke_run: Path) -> None:
    out = write_validation_report(smoke_run)
    assert out == smoke_run / "validation_report.json"
    assert out.is_file()
