"""Span-integrity audit across the 3.10 corpus.

For every record, every labeled span must resolve to its claimed surface
text when indexed into the sentence. Verifies that OCR corruption and PDF
hyphen-break transforms never mutate labeled text.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

CORPUS = Path(__file__).resolve().parents[1] / "datasets" / "v1_phase3.10" / "corpus.jsonl"


def main() -> None:
    records = [json.loads(line) for line in CORPUS.read_text().splitlines()]
    print(f"total records: {len(records)}")

    violations = 0
    by_level: Counter[str] = Counter()
    by_mode_ocr: Counter[str] = Counter()
    by_mode_pdf: Counter[str] = Counter()
    sample_violations: list[tuple[str, str, str, str]] = []

    for r in records:
        sentence = r["sentence"]
        for span in r.get("labeled_spans", []):
            s, e = span["char_span"]
            expected = span["text"]
            actual = sentence[s:e]
            if actual != expected:
                violations += 1
                by_level[r["complexity_level"]] += 1
                post = r.get("post_process", {})
                by_mode_ocr[post.get("ocr_corruption", "?")] += 1
                by_mode_pdf[post.get("pdf_artifact", "?")] += 1
                if len(sample_violations) < 10:
                    sample_violations.append((
                        r["record_id"],
                        span["span_type"],
                        expected,
                        actual,
                    ))

    print(f"\nspan-integrity violations: {violations}")
    if violations:
        print(f"by level: {dict(by_level)}")
        print(f"by ocr_corruption mode: {dict(by_mode_ocr)}")
        print(f"by pdf_artifact mode: {dict(by_mode_pdf)}")
        print("samples:")
        for rid, span_type, expected, actual in sample_violations:
            print(f"  {rid} [{span_type}] expected={expected!r} actual={actual!r}")
    else:
        print("ALL LABELED SPANS RESOLVE TO CLAIMED SURFACES.")


if __name__ == "__main__":
    main()
