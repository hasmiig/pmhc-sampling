# Fine-tuning ProteinMPNN on pMHC Structures — Data Methodology

## Overview

ProteinMPNN is the leading model for structure-based protein sequence design. However, its training data contains relatively few stable **peptide–MHC (pMHC) structures**, limiting its performance on this class of complexes. This document describes the data curation methodology used to prepare a large, high-quality dataset of predicted pMHC class I structures for fine-tuning ProteinMPNN.

This work is part of the **PMGen / FinetuneMPNN** project.

---

## Data Preparation

### Data Source

We downloaded eluted ligand and binding affinity data from the **Immune Epitope Database (IEDB)**. We restrict to **MHC class I** data exclusively.

We retain only **confirmed binders**:
- Eluted ligands (`ELlabel = 1`)
- High-affinity binding measurements (BA ≤ 500 nM)

**Non-binders are excluded.** For a sequence design model, training on structures of non-binding peptides would teach the model conformations that are not stabilized at the pMHC interface — the opposite of what we want.

### MHC-Allele Balancing

Raw pMHC datasets are heavily imbalanced — some MHC alleles have thousands of peptides while others have only a handful. This imbalance would cause the fine-tuned model to over-specialise on high-frequency alleles.

#### Iterative Median Sampling Strategy

We applied an **iterative median sampling** approach to balance allele representation:

1. Count the number of pMHC pairs observed for each MHC allele
2. Calculate the **median** of these counts across all alleles
3. Randomly sample that number of pMHC pairs for each allele
4. Remove alleles with counts below the median from the next iteration
5. Continue with a **sampling cap of 1,000 samples per allele**
   - Once an allele reaches 1,000 samples, it is excluded from further iterations
   - This ensures highly abundant alleles are capped while rarer alleles are fully represented

**Result:** A well-distributed dataset with 426 MHC-I alleles and balanced representation across the allele frequency spectrum.

### Anchor Residue Analysis

We analyzed anchor residue diversity at position **P2** (anchor 1) and the **C-terminal** position (anchor 2) to assess whether further allele-specific balancing was needed.

For each allele independently, we calculated amino acid frequencies at both anchor positions and examined 20×20 anchor combination heatmaps (P2 × C-terminal).

**Finding:** Per-allele anchor profiles show strong allele-specific dominance (e.g. HLA-A*02:01 strongly prefers leucine/methionine at P2). This is consistent with known MHC-I binding preferences.

**Decision:** Anchor residue balancing (Phase 2) was **not applied** to the final dataset. Anchor preferences in binders are a direct reflection of the MHC groove chemistry — preserving this variation is important for the design model to learn allele-specific constraints. Artificial homogenisation would destroy this signal.

---

## Structure Prediction with PMGen

### Input Preparation

Allele-balanced binder peptides were processed using the **PMGen pipeline** in initial guess mode.

For each peptide–MHC pair, anchor combinations were:
- **Enumerated without relying on NetMHCpan** (to reduce external dependencies)
- **Restricted to known MHC-I constraints:**
  - First anchor at position P1 or P2
  - Last anchor at the C-terminal residue

### Structure Generation

For each peptide–MHC pair:
1. Generated **two structures** per pair
2. Retained the structure with **higher mean peptide pLDDT**
3. Excluded samples with **mean peptide pLDDT < 80** (confidence threshold)

The pLDDT threshold ensures we only train on structures where the peptide conformation is reliably predicted — low-confidence structures would introduce noise into the design model's training signal.

**Initial yield:** 87,187 high-confidence structures spanning 426 MHC-I alleles

### Re-balancing After Prediction

Structure prediction introduces residual allele bias (some alleles fail more often due to sequence length, unusual anchor positions, or computational artifacts). We re-apply iterative median sampling to the predicted structures.

**Final dataset:** 63,817 structures (73% retention)

---

## Train/Val/Test Splitting

We use **two complementary splitting strategies** to evaluate how well the fine-tuned model generalises:

### Strategy 1: Allele-based Splitting (HLA split)

- Rank MHC alleles by frequency
- Hold out the **rarest 20%** (85 alleles out of 426) as an independent test set
- Divide remaining alleles into **5 folds** for cross-validation

**Use case:** Evaluate fine-tuning generalization to MHC alleles unseen during training — a realistic scenario for rare or newly characterized alleles.

### Strategy 2: Anchor-Combination-based Splitting (Anchor split)

- Divide folds based on **unique anchor combinations** (P2 × C-terminal pairs)
- Each observed anchor combination is held out exactly once across 5 rotating folds:
  - **Test set:** one chunk of combinations
  - **Validation set:** next chunk
  - **Training set:** remaining three chunks

**Use case:** Evaluate generalization to unseen anchor residue contexts within known alleles.

---

## ProteinMPNN Input Format

After splitting, each parquet split is converted to ProteinMPNN-ready format by `prepare_pmhc_data.py`:

- **`data.jsonl`** — one record per structure in ProteinMPNN's expected format
- **`fixed_positions.json`** — chain A (MHC heavy chain) is **fixed**; chain P (peptide) is **designable**

This setup trains ProteinMPNN to generate peptide sequences conditioned on the pMHC complex backbone — teaching it which peptide sequences are compatible with a given MHC groove geometry.

---

## Summary Statistics

| Metric | Value |
|--------|-------|
| Initial IEDB MHC-I binders | ~250,000 |
| After allele balancing | ~100,000 |
| After structure prediction (pLDDT ≥ 80) | 87,187 |
| After final re-balancing | **63,817** |
| Unique MHC-I alleles | 426 |
| Sampling cap per allele | 1,000 |
| Train/Val/Test split | 5-fold (60% / 20% / 20%) |

---

## Quality Control

### pLDDT Confidence

Mean peptide pLDDT is the primary quality metric:
- **Threshold:** ≥ 80 — high-confidence structure
- **Rationale:** Low-confidence predictions have poorly defined peptide backbones. Training ProteinMPNN on these would degrade its ability to generate sequences consistent with stable binding conformations.

### Allele-Specific Anchor Profiles

We confirmed that the final training set preserves known allele-specific anchor preferences through per-allele anchor analysis. Balancing did not over-homogenize the data.

---

## References

- IEDB: Vita, R., et al. "The immune epitope database (IEDB): 2018 update." *Nucleic Acids Research* (2019)
- ProteinMPNN: Dauparas, J., et al. "Robust deep learning–based protein sequence design using ProteinMPNN." *Science* (2022)
- AlphaFold: Jumper, J., et al. "Highly accurate protein structure prediction with AlphaFold2." *Nature* (2021)
- PMGen: structure prediction pipeline for pMHC complexes
