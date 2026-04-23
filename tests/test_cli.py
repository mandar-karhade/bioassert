"""CLI smoke test: argparse, manifest, snapshot structure via tmp_path project."""
from __future__ import annotations

import filecmp
import json
import shutil
from pathlib import Path

import pytest

from bioassert.cli import main

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
    assert runs, f"expected at least one run dir in {project_dir}/outputs"
    return runs[-1]


def test_generate_smoke_writes_corpus_report_manifest(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    rc = _run(
        [
            "generate",
            "--project", str(dest),
            "--n", "25",
            "--seed", "2024",
            "--tag", "smoke",
        ]
    )
    assert rc == 0

    run_dir = _latest_run(dest)
    assert run_dir.name.startswith("run_001_smoke_")

    corpus = run_dir / "corpus.jsonl"
    assert corpus.is_file()
    lines = corpus.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 25

    report = json.loads((run_dir / "prevalence_report.json").read_text())
    assert report["n_records"] == 25
    assert report["seed"] == 2024

    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["project"] == "nsclc_adenocarcinoma"
    assert manifest["schema_type"] == "biomarker"
    assert manifest["tag"] == "smoke"
    assert manifest["seed"] == 2024
    assert manifest["n_records"] == 25
    assert manifest["run_id"] == run_dir.name
    assert "cli_args" in manifest
    assert "bioassert_version" in manifest
    assert "git_sha" in manifest  # may be None, but key must exist
    assert manifest["config_snapshot"]["common"].startswith("snapshot/")


def test_generate_snapshot_is_byte_identical(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    rc = _run(
        [
            "generate",
            "--project", str(dest),
            "--n", "5",
            "--seed", "1",
        ]
    )
    assert rc == 0
    run_dir = _latest_run(dest)

    snapshot = run_dir / "snapshot"
    assert (snapshot / "project.json").is_file()
    assert filecmp.cmp(
        snapshot / "configs" / "common_variations.json",
        dest / "configs" / "common_variations.json",
        shallow=False,
    )
    assert filecmp.cmp(
        snapshot / "configs" / "biomarkers.json",
        dest / "configs" / "biomarkers.json",
        shallow=False,
    )
    assert filecmp.cmp(
        snapshot / "project.json",
        dest / "project.json",
        shallow=False,
    )


def test_generate_without_tag_omits_tag_segment(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    rc = _run(
        ["generate", "--project", str(dest), "--n", "3", "--seed", "7"]
    )
    assert rc == 0
    run_dir = _latest_run(dest)
    # run_NNN_YYYYMMDD-HHMMSS — exactly 3 underscore-joined segments
    assert run_dir.name.count("_") == 2


def test_generate_second_run_increments(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    assert _run(["generate", "--project", str(dest), "--n", "2", "--seed", "1"]) == 0
    assert _run(["generate", "--project", str(dest), "--n", "2", "--seed", "2"]) == 0
    runs = sorted((dest / "outputs").glob("run_*"))
    assert len(runs) == 2
    assert runs[0].name.startswith("run_001_")
    assert runs[1].name.startswith("run_002_")


def test_generate_rejects_fraction_out_of_range(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    rc = _run(
        [
            "generate",
            "--project", str(dest),
            "--n", "1",
            "--l2-fraction", "1.5",
        ]
    )
    assert rc == 2


def test_generate_rejects_fraction_sum_above_one(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    rc = _run(
        [
            "generate",
            "--project", str(dest),
            "--n", "1",
            "--l2-fraction", "0.6",
            "--l3-fraction", "0.6",
        ]
    )
    assert rc == 2


def test_generate_rejects_compound_tiers_that_dont_sum_to_one(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    rc = _run(
        [
            "generate",
            "--project", str(dest),
            "--n", "1",
            "--compound-low", "0.3",
            "--compound-high", "0.3",
        ]
    )
    assert rc == 2


def test_generate_reports_missing_project(tmp_path: Path) -> None:
    rc = _run(
        [
            "generate",
            "--project", str(tmp_path / "does_not_exist"),
            "--n", "1",
        ]
    )
    assert rc == 2


def test_main_without_command_prints_help_and_fails() -> None:
    with pytest.raises(SystemExit):
        main([])  # argparse exits on missing subcommand
