"""Ad-hoc Bug 3a audit on the sub-phase 3.10 corpus.

Bug 3a: whitespace substitutions (tab, newline_mid_sentence, double_space)
landing in prose-interior positions. Sub-phase 3.10 restricts them to
slot-boundary positions only. This script flags any residual prose-interior
whitespace substitution on L1 records.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

CORPUS = Path(__file__).resolve().parents[1] / "datasets" / "v1_phase3.10" / "corpus.jsonl"


def main() -> None:
    records = [json.loads(line) for line in CORPUS.read_text().splitlines()]
    print(f"total records: {len(records)}")

    l1 = [r for r in records if r["complexity_level"] == "L1"]
    print(f"L1 records: {len(l1)}")

    ws_modes = Counter(r["post_process"].get("whitespace") for r in l1)
    ocr_modes = Counter(r["post_process"].get("ocr_corruption") for r in l1)
    pdf_modes = Counter(r["post_process"].get("pdf_artifact") for r in l1)
    print(f"L1 whitespace modes: {dict(ws_modes)}")
    print(f"L1 ocr_corruption modes: {dict(ocr_modes)}")
    print(f"L1 pdf_artifact modes: {dict(pdf_modes)}")

    violations: list[tuple[str, str, int, str]] = []
    by_mode: Counter[str] = Counter()

    for r in l1:
        mode = r["post_process"].get("whitespace")
        if mode not in ("tab", "newline_mid_sentence", "double_space"):
            continue
        sentence: str = r["sentence"]
        spans = [
            (s["char_span"][0], s["char_span"][1])
            for s in r.get("labeled_spans", [])
        ]

        def adjacent(pos: int) -> bool:
            return any(pos == end or pos + 1 == start for start, end in spans)

        def inside(pos: int) -> bool:
            return any(start <= pos < end for start, end in spans)

        if mode == "tab":
            for i, ch in enumerate(sentence):
                if ch == "\t" and not adjacent(i) and not inside(i):
                    violations.append((r["record_id"], mode, i, sentence))
                    by_mode[mode] += 1
                    break
        elif mode == "newline_mid_sentence":
            # A standalone newline inserted by this mode; exclude the "-\n"
            # pairs owned by the pdf_artifact / hyphenation_gene_names modes.
            for i, ch in enumerate(sentence):
                if ch != "\n":
                    continue
                preceded_by_hyphen = i > 0 and sentence[i - 1] == "-"
                if preceded_by_hyphen:
                    continue
                if not adjacent(i) and not inside(i):
                    violations.append((r["record_id"], mode, i, sentence))
                    by_mode[mode] += 1
                    break
        elif mode == "double_space":
            # After insertion, the "  " pair occupies pos and pos+1. Adjacency
            # means either a span ends exactly where the pair begins (pos ==
            # end), a span starts exactly after the pair (pos + 2 == start),
            # OR a post-whitespace punctuation_variation (extra_comma) has
            # pushed the pair one position away from a span boundary
            # (sentence[pos-1] in ",.;" and pos-1 == end).
            idx = sentence.find("  ")
            if idx >= 0:
                preceding_punct = idx >= 1 and sentence[idx - 1] in ",.;:"
                adjacent_ds = any(
                    idx == end
                    or idx + 2 == start
                    or (preceding_punct and idx - 1 == end)
                    for start, end in spans
                )
                if not adjacent_ds and not inside(idx):
                    violations.append((r["record_id"], mode, idx, sentence))
                    by_mode[mode] += 1

    print(f"\nL1 Bug 3a violations: {len(violations)}")
    print(f"by mode: {dict(by_mode)}")
    if violations:
        print("\nfirst 10 violations:")
        for rid, mode, pos, sent in violations[:10]:
            print(f"  {rid} [{mode}] pos={pos} sent={sent!r}")


if __name__ == "__main__":
    main()
