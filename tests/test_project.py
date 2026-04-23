"""Project.load validation + next_run_dir naming + schema_type dispatch."""
from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pytest

from bioassert.project import Project, ProjectError

ROOT = Path(__file__).resolve().parents[1]
SHIPPED_PROJECT = ROOT / "projects" / "nsclc_adenocarcinoma"


def _clone_project(tmp_path: Path) -> Path:
    """Copy the shipped project into ``tmp_path`` so tests can mutate it."""
    dest = tmp_path / "nsclc"
    shutil.copytree(SHIPPED_PROJECT, dest)
    # ensure outputs/ is empty — shipped dir may contain a .gitkeep or nothing
    outputs = dest / "outputs"
    if outputs.exists():
        shutil.rmtree(outputs)
    return dest


def test_load_shipped_project_succeeds() -> None:
    project = Project.load(SHIPPED_PROJECT)
    assert project.name == "nsclc_adenocarcinoma"
    assert project.schema_type == "biomarker"
    assert project.common.categories, "common config should have categories"
    assert project.biomarkers.biomarkers, "biomarkers config should be populated"
    assert project.common_path.is_file()
    assert project.biomarkers_path.is_file()
    assert project.project_json_path.is_file()


def test_load_accepts_string_path(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    project = Project.load(str(dest))
    assert project.root == dest.resolve()


def test_load_rejects_missing_directory(tmp_path: Path) -> None:
    with pytest.raises(ProjectError, match="does not exist"):
        Project.load(tmp_path / "nope")


def test_load_rejects_missing_project_json(tmp_path: Path) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ProjectError, match="missing project.json"):
        Project.load(empty)


def test_load_rejects_invalid_json(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    (dest / "project.json").write_text("not { valid json", encoding="utf-8")
    with pytest.raises(ProjectError, match="not valid JSON"):
        Project.load(dest)


def test_load_rejects_non_object_json(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    (dest / "project.json").write_text("[]", encoding="utf-8")
    with pytest.raises(ProjectError, match="must contain a JSON object"):
        Project.load(dest)


def test_load_rejects_missing_required_keys(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    (dest / "project.json").write_text(
        json.dumps({"name": "incomplete"}), encoding="utf-8"
    )
    with pytest.raises(ProjectError, match="missing required keys"):
        Project.load(dest)


def test_load_rejects_unsupported_schema_type(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    meta = json.loads((dest / "project.json").read_text())
    meta["schema_type"] = "symptom"
    (dest / "project.json").write_text(json.dumps(meta))
    with pytest.raises(ProjectError, match="unsupported schema_type"):
        Project.load(dest)


def test_load_rejects_configs_missing_common(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    meta = json.loads((dest / "project.json").read_text())
    meta["configs"] = {"biomarkers": "configs/biomarkers.json"}
    (dest / "project.json").write_text(json.dumps(meta))
    with pytest.raises(ProjectError, match="missing required key 'common'"):
        Project.load(dest)


def test_load_rejects_configs_pointing_at_missing_file(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    meta = json.loads((dest / "project.json").read_text())
    meta["configs"]["common"] = "configs/does_not_exist.json"
    (dest / "project.json").write_text(json.dumps(meta))
    with pytest.raises(ProjectError, match="common config not found"):
        Project.load(dest)


def test_load_rejects_configs_not_an_object(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    meta = json.loads((dest / "project.json").read_text())
    meta["configs"] = ["configs/common_variations.json"]
    (dest / "project.json").write_text(json.dumps(meta))
    with pytest.raises(ProjectError, match="'configs' must be an object"):
        Project.load(dest)


def test_outputs_dir_derived_from_root() -> None:
    project = Project.load(SHIPPED_PROJECT)
    assert project.outputs_dir == project.root / "outputs"
    assert project.references_dir == project.root / "references"


def test_next_run_dir_starts_at_001_when_outputs_empty(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    project = Project.load(dest)
    frozen = datetime(2026, 4, 23, 23, 7, 56, tzinfo=timezone.utc)
    run_dir = project.next_run_dir(tag="smoke", now=frozen)
    assert run_dir.name == "run_001_smoke_20260423-230756"
    assert run_dir.parent == project.outputs_dir


def test_next_run_dir_increments_past_existing_runs(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    project = Project.load(dest)
    project.outputs_dir.mkdir()
    (project.outputs_dir / "run_001_20260101-000000").mkdir()
    (project.outputs_dir / "run_007_smoke_20260101-000000").mkdir()
    # non-matching names should be ignored
    (project.outputs_dir / "notes.txt").write_text("ignored")
    (project.outputs_dir / "run_broken").mkdir()
    frozen = datetime(2026, 4, 23, 23, 7, 56, tzinfo=timezone.utc)
    run_dir = project.next_run_dir(now=frozen)
    assert run_dir.name == "run_008_20260423-230756"


def test_next_run_dir_omits_tag_segment_when_none(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    project = Project.load(dest)
    frozen = datetime(2026, 1, 1, tzinfo=timezone.utc)
    run_dir = project.next_run_dir(now=frozen)
    assert run_dir.name == "run_001_20260101-000000"


def test_next_run_dir_rejects_invalid_tag(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    project = Project.load(dest)
    with pytest.raises(ProjectError, match=r"\[A-Za-z0-9"):
        project.next_run_dir(tag="bad tag/slash")
    with pytest.raises(ProjectError):
        project.next_run_dir(tag="")


def test_next_run_dir_does_not_create_directory(tmp_path: Path) -> None:
    dest = _clone_project(tmp_path)
    project = Project.load(dest)
    run_dir = project.next_run_dir(tag="smoke")
    assert not run_dir.exists(), "next_run_dir is pure — caller does mkdir"
