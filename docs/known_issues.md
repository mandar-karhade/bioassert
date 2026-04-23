# Known Rendering Issues

**Status:** Tracking variation bugs observed in generated corpora. Fixes deferred except where noted.

Every bug here was measured on the 50K `v1_phase2b` corpus (seed=42, `--l2-fraction 0.3`). Rates are record-level unless otherwise noted.

---

## Bug 1 — Frame / negative-form collision (OPEN)

**Rate:** ~5.1% (2,563 / 50K records).

**Symptom:** Adjectival frames like `{gene} was {status}.` get paired with sentence-form status phrases (e.g. "no mutation detected"), producing ungrammatical output.

Examples:
- `ROS1 was no mutation detected.`
- `PDL1 TPS 0% was no mutation detected on IHC (Clone 22C3).`
- `ALK rearrangement is no mutation detected.`
- `NT-RK was found to be no mutation detected.`
- `ALK was no evidence of.`

**Root cause:** The negation vocabulary in `common_variations.json` (`negation_phrases`, `negative_forms`) contains sentence-form phrases ("no mutation detected", "no evidence of", "none detected"). These are inserted into the `{status}` slot of adjectival frames without checking frame/phrase compatibility. `{gene} was <sentence>` is grammatically broken.

**Proposed fix (deferred):** Tag each status phrase with a `phrase_type` attribute (`adjectival` vs `sentential`). Frame selector filters phrase pool by compatible type. Adjectival frames (`was/is/were/found to be`) → adjectival-only status phrases. Fragment frames (`{gene}: {status}.`, `Result for {gene}: {status}.`) → either type.

---

## Bug 2 — Truncated `no evidence of` (OPEN)

**Rate:** ~5.5% (2,760 / 50K records).

**Symptom:** The phrase "no evidence of" renders without its required tail noun (mutation, amplification, fusion), producing incomplete sentences.

Examples:
- `NT-RK status: no evidence of.`
- `ROS1-positive — no evidence of.`
- `ALK was no evidence of.`
- `BRAF mutation was found to be no evidence of.`
- `mutational load testing: no evidence of.`

**Root cause:** The `negation_phrases` category includes a vocabulary entry that is just "no evidence of" (a prefix expecting a noun complement). When sampled standalone (no variant descriptor), the phrase is rendered truncated.

**Proposed fix (deferred):** Either (a) remove bare "no evidence of" from the vocabulary and require callers to supply the noun tail via frame context, or (b) require this entry to carry a `{noun}` placeholder that the renderer fills from a biomarker-specific noun list (mutation, fusion, amplification, rearrangement) at render time.

---

## Bug 3a — Whitespace noise corrupts L1 prose (OPEN)

**Rate:** ~7.6% of L1 records (2,634 / 34,901 L1).

**Symptom:** The post-process `whitespace` transform replaces spaces with tabs or newlines anywhere in the sentence, including mid-prose. Real NSCLC reports have tab/newline separators between *fields*, not inside running sentences.

Examples:
- `ROS1\twas found to be negative.`
- `Result\nfor c-MET: negative.`
- `RET\t— wild-type.`
- `ALK rearrangement:\twild-type.`

**Root cause:** `apply_technical_noise` whitespace-mode swap operates on any space character without respecting field boundaries. Spans are still preserved, so labels remain correct, but surface realism is broken.

**Proposed fix (deferred:)** Restrict whitespace substitution to the *first* space character that separates a frame slot from its neighbor (the `{gene}`→`{status}` or `{gene}`→`{variant}` boundary), not interior spaces within slot content or frame connective tissue (`was`, `found to be`, `Result for`, etc.). Requires tracking slot-boundary indices through `_render_with_spans`.

---

## Bug 3b — Duplicated gene name (FIXED, 2026-04-23)

**Rate before fix:** ~2.0% (1,021 / 50K records).

**Symptom:** Frames that combine `{gene}` and `{variant}` slots produced output with the gene symbol twice (once from `{gene}`, once embedded in the variant realization).

Examples (pre-fix):
- `KRAS KRAS G12D alteration detected.`
- `TMB TMB 53.2 mut/Mb was pos.`
- `c-MET MET exon 14 skipping was identified.`
- `ERBB2 ERBB2 amplification.`

**Root cause:** 44 variant realizations across the Tier 1 panel started with the gene symbol or a gene alias (e.g., `"KRAS G12C"`, `"HER2 amplification"`, `"MET exon 14 skipping"`). When a frame like `{gene} {variant}` selected an independent gene surface form, the result duplicated the gene token.

**Fix applied:** Stripped the leading gene-name prefix from every variant and `negative_forms` realization in `bioconfigs/biomarkers.json` — both space-boundary (`"KRAS G12C"`) and hyphen-boundary (`"TMB-high"`). Minimal change: 49 realizations updated in place (44 space-boundary + 5 TMB hyphen-boundary/negative_forms), all IDs and weight budgets preserved. Post-fix verification on a fresh 50K corpus: 0/50,000 duplicated-gene records.

Trade-offs:
- NTRK1/NTRK2/NTRK3 subtype realizations collapse to bare `"fusion"`. Subtype info still carried by `variant_id` in assertion metadata (labels correct by construction); only the surface form loses the subtype token.
- TMB-H/TMB-L collapse to `"H"`/`"L"` and TMB-high/TMB-low to `"high"`/`"low"`. Semantics recoverable from the `TMB` gene slot + status, but these surface forms are less canonical than the clinical `TMB-high` / `TMB-low` terms. Revisit during grammar-compatibility pass.

**See:** `_changelog` entries in `bioconfigs/biomarkers.json` for the per-realization diff.

---

## Combined impact

Pre-fix union of the four bugs: ~15% of the 50K corpus had at least one visible rendering issue. After 3b fix, the union drops by ~2pp to ~13%. Bugs 1 and 2 are the largest remaining contributors and should be addressed together as a grammar-compatibility pass on the status-phrase vocabulary.
