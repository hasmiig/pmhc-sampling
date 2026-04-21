# ProteinMPNN Fine-tuning — pMHC Data Preparation

End-to-end pipeline for curating a high-quality dataset of pMHC class I binder structures for fine-tuning ProteinMPNN. Starting from raw IEDB exports, the pipeline produces balanced, structure-predicted, QC-filtered datasets in ProteinMPNN-ready JSONL format with MHC chain fixed and peptide chain designable.

This project is part of the **PMGen / FinetuneMPNN** effort. For full methodology and design decisions, see [docs/METHODOLOGY.md](docs/METHODOLOGY.md).

---

## Overview

```
IEDB export (binders only)
    │
    ▼
01_data_preparation/      ← filter to binders, deduplicate, balance alleles
    │
    ▼
02_structure_prediction/  ← prepare PMGen input, run structure prediction
    │
    ▼
03_filtering_analysis/    ← pLDDT QC, re-balance, produce ProteinMPNN JSONL splits
    │
    ▼
data.jsonl + fixed_positions.json  (per train/val/test fold)
```

---

## Repository Structure

```
proteinmpnn-finetune-data-prep/
├── scripts/
│   ├── 01_data_preparation/
│   │   ├── pmhc_sampling.py         # Iterative median sampling (allele balancing)
│   ├── 02_structure_prediction/
│   │   ├── prepare_pmgen_input.py   # Build TSV input for PMGen
│   │   └── chunk_tsv.py             # Split large files for parallel processing
│   │   └── run_pmgen.sh             # SLURM job script for PMGen structure prediction
│   ├── 03_filtering_analysis/
│   │   ├── compute_plddt_means.py   # Extract mean pLDDT per structure
│   │   ├── analyse_plddt.py         # Visualize pLDDT distributions, apply threshold
│   │   ├── filter_map_train_prep.py # Filter, re-sample, produce train/val/test splits
│   │   └── prepare_pmhc_data.py     # Convert split parquets → ProteinMPNN JSONL
│   └── utils/
│       ├── setup_env.sh             # Environment setup for cluster jobs
│       └── submit_prepare_pmhc.sh   # Batch submission helper
├── docs/
│   ├── METHODOLOGY.md               # Full design rationale and statistics
│   └── <dated run logs>/            # Per-run exploration outputs
└── README.md
```

---

## Requirements

```
pandas
numpy
matplotlib
scipy
pyarrow     # for .parquet files
```

```bash
pip install pandas numpy matplotlib scipy pyarrow
```

Structure prediction requires PMGen and a GPU node (see `scripts/utils/setup_env.sh`).

---

## Stage 1 — Data Preparation

**Directory:** `scripts/01_data_preparation/`

**Methodology:** See [docs/2026-03-25_subsampling/README.md](docs/2026-03-25_subsampling/README.md) for detailed rationale on two-phase sampling strategy.

### 1a. Inspect raw IEDB export

Verify data structure and columns before processing:

```bash
python pmhc_sampling.py --input iedb_export.csv --inspect
```

Prints schema and first rows. Verify column names match the constants at the top of `pmhc_sampling.py` (see [Input column requirements](#input-column-requirements) below).

### 1b. Explore raw binder distributions (optional)

Understand allele and anchor composition before sampling:

```bash
python pmhc_sampling.py \
    --input iedb_export.csv \
    --mode binder \
    --explore \
    --plots exploration_raw/
```

### 1c. Filter, deduplicate, balance, and generate plots (**REQUIRED**)

After inspection/exploration, perform the full data preparation:

```bash
python pmhc_sampling.py \
    --input iedb_export.csv \
    --mode binder \
    --phases only_phase1 \
    --output binders_sampled.parquet \
    --plots results/binders/
```

This single command performs:

1. **Filter** — Keep MHC class I binders only (`assigned_label == 1.0`, `mhc_class == 1.0`)
2. **Deduplicate** — Remove exact duplicates by (peptide, allele) pair
3. **Validate** — Remove peptides with non-standard amino acids or length ≥ 15 residues
4. **Phase 1 Balancing** — Iterative median sampling to balance allele representation (cap: 1000 peptides per allele)
5. **Generate Plots** — Summary statistics and allele/anchor distributions to `results/binders/`

**Phase 1 only (not Phase 2):** Anchor residue preferences in binders reflect biological MHC-peptide constraints and should not be down-sampled. Phase 2 anchor balancing is only applied to non-binder datasets.

**Output:** `binders_sampled.parquet` — ready for Stage 2 (structure prediction).

### Input column requirements

| Column | Description |
|---|---|
| `long_mer` | Peptide amino acid sequence |
| `allele` | MHC allele name (e.g. `HLA-A*02:01`) |
| `assigned_label` | Binding label: `1.0` = binder, `0.0` = non-binder |
| `mhc_class` | MHC class: `1.0` = class I, `2.0` = class II |

If your file uses different column names, update the constants at the top of `pmhc_sampling.py`:

```python
COL_PEPTIDE = "long_mer"
COL_MHC     = "allele"
COL_LABEL   = "assigned_label"
MHC_CLASS   = "mhc_class"
SAMPLE_CAP  = 1000
```

---

## Stage 2 — Structure Prediction

**Directory:** `scripts/02_structure_prediction/`

### 2a. Prepare PMGen input

```bash
python prepare_pmgen_input.py
```

Merges sampled peptides with MHC sequences from `mhc1_encodings.csv` and writes a TSV for PMGen. Update the hardcoded paths at the top of the script to point to your parquet and MHC encoding files.

### 2b. Chunk TSV for parallel processing

Split the large TSV into smaller chunks for parallel structure prediction:

```bash
python chunk_tsv.py --input pmgen_input.tsv --size 4000 --output chunks/
```

Creates multiple `chunk_*.tsv` files in the `chunks/` directory (adjust `--size` based on available GPU memory and job duration preferences).

### 2c. Run PMGen on the cluster

```bash
for chunk in chunks/*.tsv; do
    sbatch run_pmgen.sh $chunk
done
```

Each job runs PMGen in `--initial_guess --multiple_anchors` mode, generating two structures per pMHC pair. PDB outputs are written to `outputs/<chunk>/`.

---

## Stage 3 — Filtering, Splitting, and ProteinMPNN Format Conversion

**Directory:** `scripts/03_filtering_analysis/`

### 3a. Extract pLDDT means

```bash
python compute_plddt_means.py
```

Reads PMGen PDB outputs and writes `plddt_means_binder.csv` with per-structure mean peptide pLDDT and anchor pLDDT.

### 3b. Analyse pLDDT distributions

```bash
python analyse_plddt.py
```

Generates pLDDT histograms, applies the 80-threshold filter, and saves the best structure per (allele, peptide) pair to `plddt_means_binder_best.csv`.

### 3c. Filter, re-sample, and split

```bash
python filter_map_train_prep.py \
    --plddt_csv      outputs/binder/plddt_means_binder.csv \
    --parquet        iedb_mhc1_binders.parquet \
    --mhc_encodings  data/mhc1_encodings.csv \
    --pdb_base_dir   outputs/binder \
    --output_dir     trainprep/binder/ \
    --mode           binder \
    --plddt_threshold 80 \
    --k              5 \
    --split_mode     hla
```

Four sequential stages:

1. **Filter** — select best structure per (allele, peptide); apply pLDDT threshold
2. **Map** — join back to original parquet to recover metadata and MHC sequences
3. **Resample** — re-run iterative median sampling to correct any allele bias from prediction
4. **Split** — train/val/test splits in one of two modes:
   - `hla` — rare alleles (bottom 20% by frequency) held out as test; k-fold CV on remainder
   - `anchor` — k-fold CV stratified by anchor residue combinations (P2 × C-terminal)

### 3d. Convert to ProteinMPNN JSONL format

```bash
python prepare_pmhc_data.py \
    --splits_dir trainprep/binder/splits/hla \
    --output_dir proteinmpnn_input/binder_hla \
    --split_mode hla
```

Reads each split parquet, finds the corresponding PDB files, and writes:
- `data.jsonl` — ProteinMPNN-format structure records
- `fixed_positions.json` — chain A (MHC) fixed, chain P (peptide) designable

Output structure (HLA mode):
```
proteinmpnn_input/binder_hla/
  test/
    data.jsonl
    fixed_positions.json
  fold_1/
    train/  val/
  fold_2/ ... fold_5/
```

---

## Key Numbers

| Stage | Count |
|---|---|
| Initial IEDB MHC-I binders | ~250,000 |
| After allele balancing | ~100,000 |
| After structure prediction (pLDDT ≥ 80) | 87,187 |
| After final re-balancing | **63,817** |
| Unique MHC-I alleles | 426 |

---

## Documentation

- [docs/METHODOLOGY.md](docs/METHODOLOGY.md) — rationale for binders-only approach, sampling strategy, pLDDT thresholds, splitting design
- `docs/<date>/` — per-run exploration plots and stats
