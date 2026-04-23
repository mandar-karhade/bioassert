"""Print 3.10-specific samples (OCR corrupted, PDF hyphenated, boundary whitespace) for eyeball review."""
from __future__ import annotations

import json
from pathlib import Path

CORPUS = Path(__file__).resolve().parents[1] / "datasets" / "v1_phase3.10" / "corpus.jsonl"


def show_bucket(title: str, records: list[dict]) -> None:
    print("=" * 72)
    print(title)
    print("=" * 72)
    for r in records:
        spans = {s["span_type"]: tuple(s["char_span"]) for s in r["labeled_spans"]}
        post = r["post_process"]
        print(f"\n{r['record_id']}  L={r['complexity_level']}  post={post}")
        print(f"  sentence: {r['sentence']!r}")
        for span_type, (s, e) in spans.items():
            print(f"    {span_type}=({s},{e})  -> {r['sentence'][s:e]!r}")


def main() -> None:
    records = [json.loads(line) for line in CORPUS.read_text().splitlines()]

    ocr_light = [r for r in records if r["post_process"].get("ocr_corruption") == "light"][:6]
    ocr_moderate = [r for r in records if r["post_process"].get("ocr_corruption") == "moderate"][:6]
    pdf_hyphen = [r for r in records if r["post_process"].get("pdf_artifact") == "hyphen_linebreak"][:6]
    ws_tab = [r for r in records if r["post_process"].get("whitespace") == "tab"][:4]
    ws_newline = [r for r in records if r["post_process"].get("whitespace") == "newline_mid_sentence"][:4]
    ws_double = [r for r in records if r["post_process"].get("whitespace") == "double_space"][:4]

    show_bucket("OCR corruption = light (1-2 swaps)", ocr_light)
    show_bucket("OCR corruption = moderate (1-5 swaps)", ocr_moderate)
    show_bucket("PDF artifact = hyphen_linebreak", pdf_hyphen)
    show_bucket("Whitespace = tab (boundary)", ws_tab)
    show_bucket("Whitespace = newline_mid_sentence (boundary)", ws_newline)
    show_bucket("Whitespace = double_space (boundary)", ws_double)


if __name__ == "__main__":
    main()
