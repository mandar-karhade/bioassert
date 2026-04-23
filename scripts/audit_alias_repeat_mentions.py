"""Instrument the v1_phase3.10 corpus for same-gene repeat mentions.

Question for Sub-phase 3.11 scoping: do existing L3/L4/L5 frames already
emit ≥2 mentions of the same canonical gene within one record (possibly
using different alias surfaces)? If yes, Option B (plumb alias-mixing
through existing sites) is justified. If no, Option A (add a dedicated
L5 alias-inconsistency frame family) is the only path to produce the
phenomenon at all.

Metric per record:
- mentions_per_gene[canonical]: total alias occurrences in the sentence
- distinct_surfaces_per_gene[canonical]: number of unique alias surfaces
- has_repeat_mention: any gene has mentions_per_gene >= 2
- has_alias_variance: any gene has distinct_surfaces >= 2

Aggregated by complexity_level and record_type.
"""
from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import sys

CORPUS = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    ROOT / "datasets" / "v1_phase3.10" / "corpus.jsonl"
)
BIOMARKERS = ROOT / "bioconfigs" / "biomarkers.json"


def _load_aliases() -> dict[str, list[str]]:
    cfg = json.loads(BIOMARKERS.read_text())
    aliases: dict[str, list[str]] = {}
    for key, b in cfg.items():
        if key.startswith(("$", "_")) or not isinstance(b, dict):
            continue
        if "name_forms" not in b:
            continue
        surfaces = list(b["name_forms"]["realizations"].values())
        aliases[key] = sorted(set(surfaces), key=len, reverse=True)
    return aliases


def _normalize(sentence: str) -> str:
    # Undo PDF hyphen-linebreak artifact for matching purposes
    return sentence.replace("-\n", "")


def _count_mentions(
    sentence: str, alias_list: list[str]
) -> tuple[int, set[str]]:
    """Greedy longest-first match. Returns (total_count, distinct_surfaces_matched).

    Case-insensitive. Uses word-boundary where alias is alphanumeric-only;
    for aliases with hyphens/spaces, uses literal substring match.
    """
    text = _normalize(sentence)
    lower = text.lower()
    consumed = [False] * len(lower)
    total = 0
    matched: set[str] = set()
    for alias in alias_list:
        alow = alias.lower()
        start = 0
        while True:
            idx = lower.find(alow, start)
            if idx < 0:
                break
            end = idx + len(alow)
            if any(consumed[idx:end]):
                start = idx + 1
                continue
            # word-boundary check (prevent "RET" inside "RETroactive", etc.)
            before_ok = idx == 0 or not (
                lower[idx - 1].isalnum() or lower[idx - 1] == "-"
            )
            after_ok = end == len(lower) or not (
                lower[end].isalnum() or lower[end] == "-"
            )
            if not (before_ok and after_ok):
                start = idx + 1
                continue
            for k in range(idx, end):
                consumed[k] = True
            total += 1
            matched.add(alias)
            start = end
    return total, matched


def main() -> None:
    aliases = _load_aliases()
    by_level: dict[str, Counter] = defaultdict(Counter)
    by_type: dict[str, Counter] = defaultdict(Counter)

    total_records = 0
    repeat_examples: dict[str, list[dict]] = defaultdict(list)
    variance_examples: list[dict] = []

    with CORPUS.open() as f:
        for line in f:
            rec = json.loads(line)
            total_records += 1
            level = rec["complexity_level"]
            rtype = rec["record_type"]
            sentence = rec["sentence"]

            canonical_genes = {a["gene"] for a in rec["assertions"]}

            has_repeat = False
            has_variance = False
            for gene in canonical_genes:
                if gene not in aliases:
                    continue
                count, matched = _count_mentions(sentence, aliases[gene])
                if count >= 2:
                    has_repeat = True
                if len(matched) >= 2:
                    has_variance = True
                    if len(variance_examples) < 15:
                        variance_examples.append(
                            {
                                "record_id": rec["record_id"],
                                "level": level,
                                "gene": gene,
                                "surfaces": sorted(matched),
                                "sentence": sentence,
                            }
                        )

            by_level[level]["total"] += 1
            by_type[rtype]["total"] += 1
            if has_repeat:
                by_level[level]["repeat_mention"] += 1
                by_type[rtype]["repeat_mention"] += 1
                if len(repeat_examples[level]) < 5:
                    repeat_examples[level].append(
                        {"record_id": rec["record_id"], "sentence": sentence}
                    )
            if has_variance:
                by_level[level]["alias_variance"] += 1
                by_type[rtype]["alias_variance"] += 1

    print(f"Total records: {total_records}\n")
    print("=" * 80)
    print(
        f"{'level':<8} {'total':>8} {'repeat':>8} {'repeat%':>8} "
        f"{'variance':>10} {'variance%':>10}"
    )
    print("=" * 80)
    for level in sorted(by_level):
        c = by_level[level]
        total = c["total"]
        rep = c["repeat_mention"]
        var = c["alias_variance"]
        rp = 100 * rep / total if total else 0
        vp = 100 * var / total if total else 0
        print(
            f"{level:<8} {total:>8} {rep:>8} {rp:>7.2f}% "
            f"{var:>10} {vp:>9.2f}%"
        )

    print("\nBy record_type:")
    for rtype in sorted(by_type):
        c = by_type[rtype]
        total = c["total"]
        rep = c["repeat_mention"]
        var = c["alias_variance"]
        rp = 100 * rep / total if total else 0
        vp = 100 * var / total if total else 0
        print(
            f"  {rtype:<6} total={total:>6} repeat={rep:>5} ({rp:>5.2f}%) "
            f"variance={var:>5} ({vp:>5.2f}%)"
        )

    print("\nExamples of repeat-mention by level (5 per level):")
    for level in sorted(repeat_examples):
        print(f"\n  {level}:")
        for ex in repeat_examples[level]:
            print(f"    {ex['record_id']}: {ex['sentence']!r}")

    print("\nExamples of alias-variance (distinct surfaces for same gene):")
    for ex in variance_examples[:15]:
        print(
            f"  {ex['record_id']} [{ex['level']}] gene={ex['gene']} "
            f"surfaces={ex['surfaces']}"
        )
        print(f"    {ex['sentence']!r}")


if __name__ == "__main__":
    main()
