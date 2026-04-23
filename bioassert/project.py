"""Project directory loader.

A **project** is a self-contained directory describing one corpus-generation
target: its configs, its reference materials, and its output history. The
shipped example is ``projects/nsclc_adenocarcinoma/``; additional projects
(other cohorts, future disease domains) get their own sibling directories.

Layout::

    <project_root>/
      project.json        # metadata + config paths (this file's schema)
      configs/
        common_variations.json
        biomarkers.json
      references/         # citations, prevalence sources, etc.
      outputs/            # versioned run directories

``Project.load`` validates ``project.json``, loads the configs via
:mod:`bioassert.config.loader`, and returns an immutable handle. Only
``schema_type: "biomarker"`` is recognised today; future domains will
extend :data:`_SUPPORTED_SCHEMA_TYPES`.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Union

from bioassert.config import BiomarkerConfig, CommonConfig, load_configs
from bioassert.config.validator import validate_configs

_SUPPORTED_SCHEMA_TYPES: frozenset[str] = frozenset({"biomarker"})
_REQUIRED_KEYS: frozenset[str] = frozenset(
    {"name", "schema_type", "configs"}
)
_RUN_DIR_PATTERN: re.Pattern[str] = re.compile(r"^run_(\d{3,})(?:_|$)")


class ProjectError(ValueError):
    """Raised when a project directory is malformed."""


@dataclass(frozen=True)
class Project:
    """Loaded project handle.

    Immutable. Holds resolved paths plus fully-parsed configs ready to hand
    to the generator.
    """

    root: Path
    name: str
    display_name: str
    description: str
    schema_type: str
    common: CommonConfig
    biomarkers: BiomarkerConfig
    common_path: Path
    biomarkers_path: Path
    project_json_path: Path

    @classmethod
    def load(cls, project_dir: Union[Path, str]) -> "Project":
        """Load a project directory.

        Parameters
        ----------
        project_dir:
            Path to the project root (the directory containing ``project.json``).

        Raises
        ------
        ProjectError
            If ``project.json`` is missing, malformed, references unknown
            paths, or declares an unsupported ``schema_type``.
        """
        root = Path(project_dir).expanduser().resolve()
        if not root.is_dir():
            raise ProjectError(f"project directory does not exist: {root}")

        project_json = root / "project.json"
        if not project_json.is_file():
            raise ProjectError(f"missing project.json at {project_json}")

        meta = _read_project_json(project_json)
        missing = _REQUIRED_KEYS - meta.keys()
        if missing:
            raise ProjectError(
                f"project.json is missing required keys: {sorted(missing)}"
            )

        schema_type = meta["schema_type"]
        if schema_type not in _SUPPORTED_SCHEMA_TYPES:
            raise ProjectError(
                f"unsupported schema_type {schema_type!r}; "
                f"supported: {sorted(_SUPPORTED_SCHEMA_TYPES)}"
            )

        configs_meta = meta["configs"]
        if not isinstance(configs_meta, dict):
            raise ProjectError("project.json 'configs' must be an object")
        for required_key in ("common", "biomarkers"):
            if required_key not in configs_meta:
                raise ProjectError(
                    f"project.json 'configs' missing required key {required_key!r}"
                )

        common_path = (root / configs_meta["common"]).resolve()
        biomarkers_path = (root / configs_meta["biomarkers"]).resolve()
        for label, p in (("common", common_path), ("biomarkers", biomarkers_path)):
            if not p.is_file():
                raise ProjectError(f"{label} config not found at {p}")

        common, biomarkers = load_configs(common_path, biomarkers_path)
        validate_configs(common, biomarkers)

        return cls(
            root=root,
            name=meta["name"],
            display_name=meta.get("display_name", meta["name"]),
            description=meta.get("description", ""),
            schema_type=schema_type,
            common=common,
            biomarkers=biomarkers,
            common_path=common_path,
            biomarkers_path=biomarkers_path,
            project_json_path=project_json,
        )

    @property
    def outputs_dir(self) -> Path:
        return self.root / "outputs"

    @property
    def references_dir(self) -> Path:
        return self.root / "references"

    def next_run_dir(
        self,
        tag: str | None = None,
        *,
        now: datetime | None = None,
    ) -> Path:
        """Compute the path for the next run directory.

        Scans ``outputs/`` for existing ``run_NNN_...`` folders and returns a
        ``Path`` for a sibling ``run_{NNN+1}_{tag?}_{UTC-timestamp}``. Does
        **not** create the directory — the caller is responsible for
        ``mkdir(parents=True)``.

        Parameters
        ----------
        tag:
            Optional short label inserted after the run number. Validated
            against ``[A-Za-z0-9._-]+``; pass ``None`` to omit.
        now:
            Override the UTC timestamp (useful for tests).
        """
        if tag is not None:
            if not tag or not re.fullmatch(r"[A-Za-z0-9._-]+", tag):
                raise ProjectError(
                    f"--tag must match [A-Za-z0-9._-]+ (got {tag!r})"
                )
        stamp = (now or datetime.now(timezone.utc)).strftime("%Y%m%d-%H%M%S")
        next_n = _next_run_number(self.outputs_dir)
        parts = [f"run_{next_n:03d}"]
        if tag:
            parts.append(tag)
        parts.append(stamp)
        return self.outputs_dir / "_".join(parts)


def _read_project_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ProjectError(f"{path} is not valid JSON: {exc}") from exc
    if not isinstance(data, dict):
        raise ProjectError(f"{path} must contain a JSON object")
    return data


def _next_run_number(outputs_dir: Path) -> int:
    if not outputs_dir.is_dir():
        return 1
    highest = 0
    for child in outputs_dir.iterdir():
        if not child.is_dir():
            continue
        match = _RUN_DIR_PATTERN.match(child.name)
        if not match:
            continue
        highest = max(highest, int(match.group(1)))
    return highest + 1


__all__ = ["Project", "ProjectError"]
