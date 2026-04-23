# NSCLC Biomarker Reference — Ground Truth for `bioassert` Ontology

**Purpose:** This is the canonical biomarker reference for the `bioassert` project's Layer 1 ontology. Sourced from NCCN Clinical Practice Guidelines v1.2026, published NSCLC molecular epidemiology (peer-reviewed, 2021–2026), and FDA-approved companion diagnostics.

**How to use:** Each biomarker section below maps directly to one entry in `ontology/biomarkers.yaml`. Copy structure verbatim when building the YAML.

**Last updated:** Based on sources accessed through 2026-Q1. NCCN updates quarterly; re-verify before each release.

---

## 1. The Current NCCN-Recommended Biomarker Panel (v1.2026)

Per NCCN Clinical Practice Guidelines in Oncology: Non-Small Cell Lung Cancer, Version 1.2026, the **minimum required biomarkers** for all patients with advanced or metastatic NSCLC are:

| Biomarker | Alteration type | Current targeted therapy status |
|-----------|-----------------|--------------------------------|
| **EGFR** | Mutations (activating + ex20ins) | Multiple TKIs approved |
| **ALK** | Gene fusions/rearrangements | Multiple TKIs approved |
| **ROS1** | Gene fusions/rearrangements | Multiple TKIs approved |
| **BRAF** | V600E mutation (primarily) | Dabrafenib + trametinib approved |
| **KRAS** | G12C mutation (actionable); others | Sotorasib, adagrasib approved for G12C |
| **MET** | Exon 14 skipping + amplification | Capmatinib, tepotinib approved |
| **RET** | Gene fusions/rearrangements | Selpercatinib, pralsetinib approved |
| **ERBB2 (HER2)** | Mutations (primarily ex20ins) | T-DXd, sevabertinib (2026) approved |
| **NTRK1/2/3** | Gene fusions | Larotrectinib, entrectinib, repotrectinib |
| **PD-L1** | Expression (IHC, TPS/CPS) | Multiple IO agents approved |

**Also increasingly tested:**
- **NRG1** fusions (rare but actionable with zenocutuzumab)
- **TMB** (tumor mutational burden) — pembrolizumab approval for TMB-high
- **MSI / dMMR** — rare in NSCLC but screened
- **TP53, STK11, KEAP1** — co-mutations with prognostic/predictive value

**Testing method:** NCCN strongly advises broad molecular profiling via NGS. ctDNA (liquid biopsy) is acceptable and can be done concurrently with tissue testing. PD-L1 must be done on tissue (IHC); ctDNA not acceptable for PD-L1.

---

## 2. Prevalence Summary — Adenocarcinoma (Primary Histology)

Prevalence in advanced/metastatic lung adenocarcinoma. Ranges reflect variation across studies and populations. Sources cited per biomarker below.

| Biomarker | Western/Caucasian | East Asian | Notes |
|-----------|-------------------|------------|-------|
| KRAS | 25–32% | 10–15% | Most common in Western; lower in Asian |
| EGFR | 15–20% | 40–50% | Strongly ethnicity-associated |
| ALK fusions | 3–7% | 3–7% | Similar across populations |
| MET ex14 skipping | 2–4% | 2–4% | |
| MET amplification | 1–5% | 1–5% | De novo amplification |
| BRAF | 1.5–4% | 1–2% | V600E ~50% of BRAF; non-V600 ~50% |
| ERBB2/HER2 mutations | 1–3% | 2–3% | Mostly ex20ins |
| ROS1 fusions | 1–2% | 1–2% | |
| RET fusions | 1–2% | 1–2% | |
| NTRK fusions | <0.5% | <0.5% | Very rare |
| NRG1 fusions | <0.5% | <0.5% | Enriched in invasive mucinous adenocarcinoma |
| PD-L1 ≥50% (TPS) | ~30% | ~30% | |
| PD-L1 1–49% | ~30% | ~30% | |
| PD-L1 <1% | ~40% | ~40% | |
| TMB-high (≥10 mut/Mb) | 15–20% | 15–20% | |

## 3. Prevalence Summary — Squamous Cell Carcinoma

Substantially different landscape. Fewer actionable drivers.

| Biomarker | Prevalence |
|-----------|-----------|
| KRAS G12C | 2–10% |
| EGFR mutations | 1–5% |
| FGFR1 amplification | 10–20% |
| PIK3CA mutations | 5–15% |
| DDR2 mutations | 3–4% |
| All other major drivers (ALK, ROS1, RET, NTRK, BRAF V600E, MET ex14, ERBB2) | <0.5% each |
| PD-L1 ≥50% | ~30% |

---

## 4. Detailed Biomarker Profiles

### 4.1 EGFR — Epidermal Growth Factor Receptor

**Canonical:** Epidermal Growth Factor Receptor
**Aliases / abbreviations:** EGFR, HER1, ERBB1, ERBB, c-erbB-1
**Spelling variants:** EGF-R, EGF R, EGFr, egfr
**Gene location:** 7p11.2

**Prevalence:**
- Overall NSCLC: 10–15% (Western), 40–50% (East Asian)
- Adenocarcinoma (Western): 15–20%
- Adenocarcinoma (East Asian): 40–50%
- Enriched in: women, never-smokers, adenocarcinoma histology
- Squamous cell: 1–5%

**Subvariant distribution (within EGFR-mutated tumors):**
| Variant | % of EGFR mutations | Clinical category |
|---------|---------------------|-------------------|
| Exon 19 deletion (ex19del) | ~46% | Classical activating |
| Exon 21 L858R | ~38% | Classical activating |
| Exon 20 insertion (ex20ins) | 4–12% | Distinct entity, TKI-resistant |
| Uncommon (G719X, L861Q, S768I) | 5–10% | Atypical |
| T790M | Mostly resistance (post-TKI) | Resistance |
| C797S | Post-osimertinib resistance | Resistance |

**Specific variants to include in ontology:**

*Classical activating:*
- `L858R` — exon 21 point mutation (Leu858Arg; c.2573T>G)
- `exon_19_deletion` — includes E746_A750del, L747_P753delinsS, L747_T751del, others
  - Aliases: ex19del, E19del, del19, exon 19 del

*Exon 20 insertions (60+ unique variants known):*
- `V769_D770insASV` — A767_V769dup (ASV insertion); most common, PCR-detectable
- `S768_D770dup` — SVD insertion; second most common, PCR-detectable
- `N771_H773dup` — NPH insertion; PCR undetectable (requires NGS)
- `D770delinsGY` — NGS-required
- `H773_V774insAH` — NGS-required
- Generic: `exon_20_insertion` — catch-all for heterogeneous variants

*Uncommon/atypical:*
- `G719X` — includes G719A, G719C, G719S
- `L861Q`
- `S768I`

*Resistance mutations (secondary):*
- `T790M` — most common resistance to 1st/2nd gen TKIs
- `C797S` — acquired resistance to osimertinib
- `L718Q` — rare osimertinib resistance

**Test methods:** NGS (preferred), PCR (may miss >50% of ex20ins variants), Sanger, ddPCR, liquid biopsy NGS

**Approved therapies (as of 2026):** osimertinib, amivantamab + lazertinib, gefitinib, erlotinib, afatinib, dacomitinib, mobocertinib (withdrawn 2023), sunvozertinib, datopotamab deruxtecan (Dato-DXd; post-osimertinib)

**Clinical associations:** NSCLC (adenocarcinoma predominantly), never-smokers, female, Asian ethnicity

**Prevalence sources:**
- NCCN Guidelines NSCLC v1.2026
- AACR Project GENIE v13+
- Shi et al. PIONEER study (Asian NSCLC)
- Flatiron database (US real-world, 2011–2020)

---

### 4.2 KRAS — Kirsten Rat Sarcoma Viral Oncogene Homolog

**Canonical:** Kirsten Rat Sarcoma Viral Oncogene Homolog
**Aliases:** KRAS, K-RAS, KRAS2, C-K-RAS, c-Ki-ras2
**Spelling variants:** K-Ras, K-ras, kras
**Gene location:** 12p12.1

**Prevalence:**
- Adenocarcinoma (Western): 25–32%
- Adenocarcinoma (East Asian): 10–15%
- Enriched in: smokers, men, adenocarcinoma
- Squamous cell: 2–10% (mostly G12C)

**Subvariant distribution (within KRAS-mutated NSCLC):**
| Variant | % of KRAS mutations |
|---------|---------------------|
| G12C | ~40% |
| G12V | 18–21% |
| G12D | 17–18% |
| G12A | 7% |
| G13C | 5–7% |
| Q61H | 3–5% |
| Other G12/G13/Q61 | remainder |

**Specific variants to include in ontology:**
- `G12C` — glycine to cysteine at codon 12; **only FDA-approved actionable KRAS mutation** (sotorasib, adagrasib)
- `G12V` — glycine to valine
- `G12D` — glycine to aspartate (MRTX1133 in development)
- `G12A` — glycine to alanine
- `G12R` — glycine to arginine
- `G13C` — glycine to cysteine at codon 13
- `G13D`
- `Q61H` — glutamine to histidine at codon 61
- `Q61L`
- `Q61R`
- `A146T` — atypical

**Test methods:** NGS, PCR (good coverage for G12C), liquid biopsy NGS

**Approved therapies (G12C only):** sotorasib (AMG 510), adagrasib (MRTX849)

**Clinical associations:** NSCLC (adenocarcinoma), smokers, mutually exclusive with EGFR/ALK/ROS1 in most cases

**Common co-mutations:**
- TP53 (~40–50%)
- STK11 (~20%, negative predictive for IO)
- KEAP1 (~15%, negative predictive for IO)
- SMARCA4 (associated with smoking)

**Prevalence sources:**
- Nature Signal Transduction and Targeted Therapy (2025) — KRAS variant distribution
- Project GENIE
- Multiple published cohort studies

---

### 4.3 ALK — Anaplastic Lymphoma Kinase

**Canonical:** Anaplastic Lymphoma Kinase
**Aliases:** ALK, CD246
**Spelling variants:** ALK, alk
**Gene location:** 2p23.2-p23.1

**Alteration type:** Gene fusions/rearrangements (NOT point mutations in primary disease)

**Prevalence:**
- NSCLC overall: 3–7%
- Adenocarcinoma: 3–7%
- Enriched in: younger patients (median ~50), never or light smokers

**Fusion partners (EML4 dominant, many others exist):**
- `EML4-ALK` — dominant partner (~85% of ALK fusions)
  - Multiple variants: V1 (E13:A20), V2 (E20:A20), V3 (E6:A20), V5 (E2:A20)
  - V3a/V3b associated with worse outcomes
- `KIF5B-ALK`
- `KLC1-ALK`
- `TFG-ALK`
- `HIP1-ALK`
- Generic: `ALK_fusion` or `ALK_rearrangement`

**Resistance mutations (post-TKI, secondary point mutations in kinase domain):**
- `G1202R` — solvent-front mutation, resistant to multiple TKIs
- `L1196M` — gatekeeper mutation
- `F1174L`
- `C1156Y`
- `1151Tins`
- `G1269A`
- `I1171N/T/S`

**Test methods:**
- FISH (gold standard historically, Vysis ALK Break Apart)
- IHC (Ventana ALK D5F3) — highly sensitive, often used for screening
- NGS (DNA and RNA-based; RNA more sensitive for fusions)
- PCR (fusion-specific, limited partner coverage)

**Approved therapies:** crizotinib, ceritinib, alectinib (1L preferred), brigatinib, lorlatinib, ensartinib

**Clinical associations:** NSCLC (adenocarcinoma), signet-ring cell features, younger, never-smokers

**Prevalence sources:**
- Multiple NSCLC registries
- NCCN Guidelines

---

### 4.4 ROS1 — ROS Proto-Oncogene 1, Receptor Tyrosine Kinase

**Canonical:** ROS Proto-Oncogene 1, Receptor Tyrosine Kinase
**Aliases:** ROS1, ROS, c-ros
**Spelling variants:** ROS-1, ros1
**Gene location:** 6q22.1

**Alteration type:** Gene fusions/rearrangements

**Prevalence:**
- NSCLC: 1–2%
- Adenocarcinoma: 1–2%
- Enriched in: younger, never-smokers

**Fusion partners:**
- `CD74-ROS1` — most common
- `EZR-ROS1`
- `SLC34A2-ROS1`
- `TPM3-ROS1`
- `SDC4-ROS1`
- `GOPC-ROS1`
- Generic: `ROS1_fusion` / `ROS1_rearrangement`

**Resistance mutations (post-TKI):**
- `G2032R` — most common, solvent-front mutation, crizotinib-resistant
- `D2033N`
- `L2026M` — gatekeeper
- `S1986F/Y`

**Test methods:** FISH, IHC (less established than ALK), NGS (RNA-based preferred), PCR

**Approved therapies:** crizotinib, entrectinib, lorlatinib, repotrectinib

**Clinical associations:** NSCLC (adenocarcinoma), never-smokers, younger patients

---

### 4.5 BRAF — B-Raf Proto-Oncogene, Serine/Threonine Kinase

**Canonical:** B-Raf Proto-Oncogene, Serine/Threonine Kinase
**Aliases:** BRAF, BRAF1, B-raf, v-Raf murine sarcoma viral oncogene homolog B
**Spelling variants:** B-RAF, b-raf, braf
**Gene location:** 7q34

**Prevalence:**
- NSCLC: 1.5–4%
- V600E: ~50% of BRAF-mutated NSCLC
- Non-V600: ~50%

**Subvariant categories:**

*Class I (V600 mutations — highest activity, monomer-dependent):*
- `V600E` — valine to glutamate at codon 600; **most clinically actionable**
- `V600K`
- `V600D`
- `V600R`
- `V600M`

*Class II (non-V600, RAS-independent dimer activators):*
- `G469A/V/R`
- `K601E`
- `L597Q/R/V`

*Class III (kinase-impaired, require RAS activation):*
- `G466V/A/E`
- `N581S/I`
- `D594G/N/V`
- `G596R`

*Fusions (rare):*
- `TRIM24-BRAF`
- Other rare BRAF fusions

**Test methods:** NGS, PCR (V600E-specific kits widely available), IHC (anti-V600E antibody)

**Approved therapies:** dabrafenib + trametinib (BRAF V600E), encorafenib + binimetinib (BRAF V600E)

**Clinical associations:** NSCLC, smokers (more common than in EGFR mutants), mutually exclusive with EGFR/ALK/ROS1/MET ex14/RET in most cases

---

### 4.6 MET — MET Proto-Oncogene, Receptor Tyrosine Kinase

**Canonical:** MET Proto-Oncogene, Receptor Tyrosine Kinase (hepatocyte growth factor receptor)
**Aliases:** MET, HGFR, c-Met, c-MET, RCCP2
**Spelling variants:** c-Met, cMET, c-met, met
**Gene location:** 7q31.2

**Alteration types (three distinct actionable forms):**

*1. MET exon 14 skipping mutations:*
- Prevalence: 2–4% of NSCLC
- Splice site mutations that cause skipping of exon 14 → loss of Y1003 ubiquitination site → MET stabilization
- Many distinct splice site variants; ontology should have `MET_ex14_skipping` as a category
- Test methods: RNA-seq (most sensitive), DNA NGS (captures some), liquid biopsy

*2. MET amplification (de novo):*
- Prevalence: 1–5% de novo
- Measured as MET/CEP7 ratio (FISH) or copy number (NGS)
- High-level amplification (≥5 copies or ratio ≥2.2) more likely actionable
- Often co-occurs with EGFR resistance
- Variable definition thresholds across assays

*3. MET amplification (acquired, post-EGFR TKI):*
- Most common resistance mechanism to osimertinib (~15–25%)
- Different therapeutic implications than de novo

**Approved therapies:** capmatinib, tepotinib (both for MET ex14 skipping)

**Test methods:** NGS (DNA + RNA), FISH, IHC (screening only), liquid biopsy

---

### 4.7 RET — Ret Proto-Oncogene

**Canonical:** Ret Proto-Oncogene (Rearranged during Transfection)
**Aliases:** RET, CDHF12, PTC, RET51
**Spelling variants:** RET, ret
**Gene location:** 10q11.21

**Alteration type:** Gene fusions (primary); rarely point mutations

**Prevalence:** 1–2% of NSCLC (adenocarcinoma)

**Fusion partners:**
- `KIF5B-RET` — most common (~70% of RET fusions in NSCLC)
- `CCDC6-RET`
- `NCOA4-RET`
- `TRIM33-RET`
- Generic: `RET_fusion`

**Resistance mutations:**
- `G810R/S/C` — solvent-front mutation

**Test methods:** NGS (RNA-based preferred), FISH, RT-PCR

**Approved therapies:** selpercatinib, pralsetinib

---

### 4.8 ERBB2 (HER2)

**Canonical:** Erb-B2 Receptor Tyrosine Kinase 2
**Aliases:** ERBB2, HER2, HER-2, NEU, CD340
**Spelling variants:** HER2, HER-2, Her-2, her2
**Gene location:** 17q12

**Alteration types (three distinct, differing therapeutic implications):**

*1. ERBB2 mutations (primarily exon 20 insertions):*
- Prevalence: 1–3% of NSCLC
- `YVMA_insertion` (A775_G776insYVMA) — most common variant
- Other ex20 insertions: `G776delinsVC`, `P780_Y781insGSP`
- Point mutations: `L755S`, `V842I`, `G660D`
- Generic: `ERBB2_exon_20_insertion`

*2. ERBB2 amplification:*
- Distinct from mutation
- Less clearly actionable in NSCLC than in breast/gastric

*3. ERBB2 overexpression (IHC):*
- NCCN now recommends IHC testing
- Different from amplification by ISH

**Test methods:** NGS (mutations), IHC (overexpression), ISH (amplification)

**Approved therapies:** trastuzumab deruxtecan (T-DXd; ERBB2-mutant), sevabertinib (2026 approval for ERBB2-mutant NSCLC)

---

### 4.9 NTRK1, NTRK2, NTRK3

**Canonical:** Neurotrophic Receptor Tyrosine Kinase 1/2/3
**Aliases:**
- NTRK1 = TRKA
- NTRK2 = TRKB
- NTRK3 = TRKC

**Alteration type:** Gene fusions (primary)

**Prevalence:** <0.5% of NSCLC (very rare)

**Fusion partners (highly variable, many partners):**
- NTRK1: `MPRIP-NTRK1`, `CD74-NTRK1`, `TPR-NTRK1`
- NTRK2: `TRIM24-NTRK2` (rare)
- NTRK3: `ETV6-NTRK3` (most common NTRK3 partner)
- Generic: `NTRK_fusion` (pan-NTRK)

**Resistance mutations:**
- `G595R` (NTRK1) — solvent-front
- `G623R` (NTRK3) — solvent-front

**Test methods:** NGS (RNA-based highly preferred due to large intronic regions in DNA), pan-TRK IHC (screening), FISH (limited)

**Approved therapies:** larotrectinib (pan-NTRK), entrectinib (NTRK + ROS1 + ALK), repotrectinib

---

### 4.10 NRG1 — Neuregulin 1

**Canonical:** Neuregulin 1
**Aliases:** NRG1, HRG, HRG1, NDF, GGF, GGF2, SMDF
**Spelling variants:** NRG1, nrg1

**Alteration type:** Gene fusions

**Prevalence:** <0.5% of NSCLC (enriched in invasive mucinous adenocarcinoma, ~7–27%)

**Fusion partners:**
- `CD74-NRG1` — most common
- `SLC3A2-NRG1`
- `ATP1B1-NRG1`
- `SDC4-NRG1`
- Generic: `NRG1_fusion`

**Test methods:** RNA NGS (strongly preferred), DNA NGS (misses many), RT-PCR

**Approved therapies (2026):** zenocutuzumab (FDA-approved for NRG1-positive NSCLC)

---

### 4.11 PD-L1 — Programmed Death Ligand 1

**Canonical:** Programmed Death Ligand 1 (CD274)
**Aliases:** PD-L1, PDL1, CD274, B7-H1, B7H1
**Spelling variants:** PDL-1, PD L1, pd-l1, pdl1

**Alteration type:** Protein expression level (NOT a mutation)

**Measurement methods:**
- **TPS** (Tumor Proportion Score) — % of tumor cells with membrane staining
- **CPS** (Combined Positive Score) — tumor + immune cells / tumor cells × 100
- **IC** (Immune Cell score) — % of tumor area with immune cell staining

**Clinical cutoffs (TPS):**
- `PD-L1 ≥50%` — high expression; pembrolizumab monotherapy eligible
- `PD-L1 1–49%` — low-intermediate; IO-chemo combinations
- `PD-L1 <1%` — negative

**IHC antibody clones (each has different cutoffs and companion diagnostics):**
- `22C3` (Dako) — pembrolizumab CDx
- `28-8` (Dako) — nivolumab
- `SP263` (Ventana) — durvalumab, atezolizumab
- `SP142` (Ventana) — atezolizumab

**Prevalence:**
- PD-L1 ≥50%: ~30%
- PD-L1 1–49%: ~30%
- PD-L1 <1%: ~40%

**Test methods:** IHC only (ctDNA NOT acceptable)

**Approved therapies:** pembrolizumab, nivolumab + ipilimumab, atezolizumab, cemiplimab, durvalumab, tislelizumab, and combinations

---

### 4.12 TMB — Tumor Mutational Burden

**Canonical:** Tumor Mutational Burden
**Aliases:** TMB, mutational load
**Unit:** mutations per megabase (mut/Mb)

**Alteration type:** Quantitative summary metric across whole exome or large panel

**Clinical cutoffs:**
- TMB-high: ≥10 mut/Mb (FDA-approved threshold for pembrolizumab tissue-agnostic)
- TMB-low: <10 mut/Mb

**Prevalence:** 15–20% of NSCLC is TMB-high

**Test methods:** Whole exome sequencing (gold standard), large targeted NGS panels (FoundationOne, MSK-IMPACT, etc. — validated against WES)

**Approved therapies:** pembrolizumab (tissue-agnostic TMB-high approval)

**Caveats:** Assay-dependent; cutoffs not universally reproducible; correlation with IO response less strong than PD-L1 in NSCLC specifically

---

### 4.13 MSI / dMMR

**Canonical:** Microsatellite Instability / Mismatch Repair Deficiency
**Aliases:** MSI, MSI-H, dMMR, MMR-D

**Alteration type:** Genomic instability phenotype

**Categories:**
- MSI-H (MSI-High)
- MSI-L (MSI-Low)
- MSS (Microsatellite Stable)
- dMMR (deficient mismatch repair — IHC measure)

**Prevalence:** <1% of NSCLC (rare)

**Test methods:** PCR (Bethesda markers), NGS (MSIsensor, mSINGS), IHC (MLH1, MSH2, MSH6, PMS2 proteins)

**Approved therapies:** pembrolizumab (tissue-agnostic for MSI-H/dMMR)

---

## 5. Co-Mutations (Important for Prognostic/Predictive Context)

These are NOT primary drivers but appear frequently in NSCLC reports and are clinically meaningful. Ontology should include them with a `biomarker_category: co_mutation` field.

### 5.1 TP53 — Tumor Protein P53
- Most common co-mutation across NSCLC subtypes (~50% in EGFR-mutated, 25–50% in ALK/ROS1/RET fusions, 40–50% in KRAS-mutant)
- Associated with worse outcomes on EGFR TKIs
- Aliases: TP53, p53

### 5.2 STK11 — Serine/Threonine Kinase 11
- ~20% co-mutation with KRAS
- Predictor of IO resistance
- Aliases: STK11, LKB1

### 5.3 KEAP1 — Kelch-like ECH-Associated Protein 1
- ~15% co-mutation with KRAS
- Predictor of IO resistance
- Aliases: KEAP1, INRF2

### 5.4 SMARCA4
- Co-mutates with KRAS in smokers
- Associated with worse outcomes
- Aliases: SMARCA4, BRG1

### 5.5 CDKN2A/B
- Tumor suppressor loss
- Affects response to various therapies

### 5.6 RB1, PTEN, PIK3CA, MDM2
- Co-alterations with clinical significance in specific contexts

---

## 6. Squamous-Specific Markers (Lower Priority for Lung Adenocarcinoma Focus)

If the corpus extends beyond adenocarcinoma:
- **FGFR1 amplification** (10–20% squamous)
- **PIK3CA mutations** (5–15% squamous)
- **DDR2 mutations** (3–4% squamous)

---

## 7. Emerging / Research Biomarkers (Not Standard-of-Care Yet)

Mention in ontology as lower priority, but include because they appear in research-oriented pathology reports:
- FGFR2/3 alterations
- AXL alterations
- AKT1 mutations
- MAP2K1 (MEK1) mutations
- NRAS mutations
- HRAS mutations
- EGFR fusions (rare, distinct from point mutations)

---

## 8. Complete Biomarker List for v1 Ontology (Prioritized)

### Priority Tier 1 (MUST include in Phase 1 ontology — NCCN minimum panel)
1. EGFR
2. ALK
3. ROS1
4. BRAF
5. KRAS
6. MET
7. RET
8. ERBB2 (HER2)
9. NTRK1
10. NTRK2
11. NTRK3
12. PD-L1

### Priority Tier 2 (include by end of Phase 1 or Phase 2)
13. NRG1
14. TMB
15. MSI / dMMR
16. TP53 (co-mutation)
17. STK11 (co-mutation)
18. KEAP1 (co-mutation)

### Priority Tier 3 (Phase 3+ or follow-up corpora)
19. SMARCA4
20. CDKN2A/B
21. RB1
22. PTEN
23. PIK3CA
24. FGFR1/2/3
25. NRAS
26. HRAS

---

## 9. Recommendations for Phase 1 Implementation

Start the Phase 1 ontology with the **three biomarkers originally specified** (EGFR, ALK, KRAS) from the project_spec.md, but use this document to fill in their full profiles with clinical accuracy. These three give you:

1. **EGFR** — point mutations with sub-variants, including the compositionally complex ex20ins category. Tests the generator on heterogeneous variant distribution.
2. **ALK** — gene fusions with fusion partner structure. Tests the generator on a fundamentally different alteration type from mutations.
3. **KRAS** — point mutations with clear single-locus actionable variant (G12C). Tests the generator on actionable vs non-actionable subvariant classification.

Expanding to the full Tier 1 panel is the Phase 2 deliverable.

---

## 10. Critical Testing Methodology Notes

**For the NER generator, keep these in mind because they appear in real reports:**

- **IHC vs NGS vs FISH vs PCR** — each method has vocabulary conventions. IHC reports use "positive/negative/equivocal" with quantitative stains (TPS, H-score). FISH reports break-apart signals. NGS reports variant allele frequencies (VAF).
- **Concordance between ctDNA and tissue** — reports often mention both. "Tissue negative but ctDNA positive for..." is a real pattern.
- **"Wild-type" language** — very common for negative results. Don't just train on "negative" — include "wild-type", "wt", "WT", "no mutation detected."
- **VAF and allele frequency** — reports often include "EGFR L858R detected at 42% VAF". Your ontology should support a `measurement` field for quantitative mutation data.
- **Tumor proportion vs cell fraction** — different metrics for PD-L1 vs other biomarkers. Don't conflate.
- **"Not tested" vs "failed QC" vs "insufficient sample"** — all distinct statuses; important for completeness tracking in EHR extraction.

---

## 11. Prevalence Sources (Citation List for `prevalence_sources.bib`)

Primary sources that should be cited in the dataset card:

1. **NCCN Guidelines NSCLC v1.2026** — current standard-of-care panel
2. **AACR Project GENIE** — real-world cancer genomics consortium, v13+
3. **Shi et al., PIONEER study** — East Asian NSCLC EGFR prevalence
4. **Campbell et al. 2016** — distinct patterns of somatic genome alterations in lung adenocarcinoma and squamous cell (TCGA)
5. **Jordan et al. 2017 Cancer Discovery** — MSK-IMPACT cohort, co-mutation landscape
6. **Skoulidis & Heymach 2019 Nat Rev Cancer** — co-mutations and therapy response
7. **Robichaux et al. 2021 Nature** — structure-function EGFR classification
8. **Flatiron Health database** — US real-world NSCLC cohorts
9. **LC-SCRUM-Japan / Asia** — Asian cohort data
10. **Riely et al. 2024 JNCCN** — NCCN guidelines v4.2024 updates
11. **Gainor et al. 2016** — ALK resistance mutations
12. **Advances in molecular pathology and therapy of NSCLC — Signal Transduction and Targeted Therapy, 2025**

---

**This reference document should be the authoritative source that Claude Code pulls from when building `ontology/biomarkers.yaml`. Every prevalence number and variant name in this document is traceable to peer-reviewed or guideline sources. Do not introduce biomarkers or variants that aren't in this document without checking a source first.**
