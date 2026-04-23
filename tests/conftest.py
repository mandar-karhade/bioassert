"""Shared fixtures for Phase 2a tests."""
from __future__ import annotations

import random
from pathlib import Path

import pytest

from bioassert.config import BiomarkerConfig, CommonConfig, load_configs

ROOT = Path(__file__).resolve().parents[1]
COMMON_PATH = ROOT / "bioconfigs" / "common_variations.json"
BIOMARKERS_PATH = ROOT / "bioconfigs" / "biomarkers.json"

MUTATION_BIOMARKERS: tuple[str, ...] = (
    "EGFR",
    "KRAS",
    "BRAF",
    "ERBB2",
    "MET",
    "ALK",
    "ROS1",
    "RET",
    "NTRK",
)

EXPRESSION_BIOMARKERS: tuple[str, ...] = ("PD-L1", "TMB")

PANEL_BIOMARKERS: tuple[str, ...] = MUTATION_BIOMARKERS + EXPRESSION_BIOMARKERS


@pytest.fixture(scope="session")
def configs() -> tuple[CommonConfig, BiomarkerConfig]:
    return load_configs(COMMON_PATH, BIOMARKERS_PATH)


@pytest.fixture(scope="session")
def common(configs) -> CommonConfig:
    return configs[0]


@pytest.fixture(scope="session")
def biomarkers(configs) -> BiomarkerConfig:
    return configs[1]


@pytest.fixture
def rng() -> random.Random:
    return random.Random(12345)
