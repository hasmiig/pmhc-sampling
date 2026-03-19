# pMHC Peptide Sampling Pipeline

A two-phase iterative sampling pipeline for balancing peptide-MHC (pMHC) class I datasets — for both binder and non-binder peptides — prior to downstream tasks such as structure prediction and model training.

---

## Background

Raw pMHC datasets are heavily imbalanced: some MHC alleles have thousands of peptides while others have only a handful, and within each allele certain amino acids dominate the anchor positions (P2 and last position). This imbalance can introduce inductive bias into any model trained on the data.

This pipeline addresses this in two sequential phases:

- **Phase 1** — balances peptide counts *across* MHC alleles using iterative median sampling, capped at `SAMPLE_CAP` (default: 1000) peptides per allele
- **Phase 2** — within each allele independently, balances anchor residue diversity at position P2 (anchor 1) and the last peptide position (anchor 2) using the same iterative median sampling strategy

The pipeline can be run in `--mode nonbinder` or `--mode binder` to process each class separately.

---

## Requirements

```
pandas
numpy
matplotlib
scipy
pyarrow        # for .parquet files
```

Install with:

```bash
pip install pandas numpy matplotlib scipy pyarrow
```

---

## Input Format

The pipeline accepts `.parquet`, `.tsv`, `.txt`, or `.csv` files. The input file must contain at minimum the following columns:

| Column | Description |
|---|---|
| `long_mer` | Peptide amino acid sequence |
| `allele` | MHC allele name (e.g. `HLA-A*02:01`) |
| `assigned_label` | Binding label: `1.0` = binder, `0.0` = non-binder |
| `mhc_class` | MHC class: `1.0` = class I, `2.0` = class II |

> **Important:** these column names are defined as constants at the top of the script (`COL_PEPTIDE`, `COL_MHC`, `COL_LABEL`, `MHC_CLASS`). If your input file uses different column names, update these constants before running — see [Adapting to a New Dataset](#adapting-to-a-new-dataset) below.

---

## Usage

### Step 0 — Inspect the file

Before running the pipeline on a new dataset, always inspect it first:

```bash
python pmhc_sampling.py --input data.parquet --inspect
```

This prints the file schema, column names, data types, and the first 5 rows without loading the full dataset. Use this to verify that:
- The column names in your file match the constants at the top of the script (`COL_PEPTIDE`, `COL_MHC`, `COL_LABEL`, `MHC_CLASS`)
- The label values are as expected (`0.0` / `1.0`)
- The MHC class column is correctly populated

If the column names differ from the defaults, update the constants at the top of `pmhc_sampling.py` to match your dataset before proceeding.

---

### Step 1 — Explore the raw data (optional but recommended)

Run in exploration mode to generate diagnostic plots of the raw data without doing any sampling:

```bash
# non-binders
python pmhc_sampling.py --input data.parquet --mode nonbinder --explore --plots plots/nonbinders/

# binders
python pmhc_sampling.py --input data.parquet --mode binder --explore --plots plots/binders/
```

This saves plots showing:
- Allele distribution (sorted by count)
- Top 20 most and least frequent alleles
- Peptide length distribution
- Anchor residue (P2 and last position) frequencies
- Data source distribution (if a `source` column is present)

Use these to understand the shape of the raw data before committing to a full run.

---

### Step 2 — Run the full pipeline

```bash
# non-binders
python pmhc_sampling.py \
    --input data.parquet \
    --mode nonbinder \
    --output nonbinders_sampled.parquet \
    --plots plots/nonbinders/

# binders
python pmhc_sampling.py \
    --input data.parquet \
    --mode binder \
    --output binders_sampled.parquet \
    --plots plots/binders/
```

The pipeline will:
1. Filter to MHC class I peptides of the specified type (binder or non-binder)
2. Remove duplicates (same peptide + allele pair)
3. Remove peptides with non-standard amino acid characters (e.g. placeholders, rare AAs)
4. Remove peptides of length ≥ 15 (too few survive sampling to be useful)
5. Run Phase 1: iterative median sampling across alleles
6. Run Phase 2: iterative median sampling of anchor residues within each allele
7. Save all plots, statistics, and the final sampled dataset

---

## All CLI Arguments

| Argument | Required | Description |
|---|---|---|
| `--input` | Yes | Path to input file (`.parquet`, `.tsv`, `.csv`) |
| `--mode` | Yes (for sampling) | `binder` or `nonbinder` |
| `--output` | Yes (for sampling) | Path to save the final sampled dataset |
| `--plots` | No | Directory for plots and stats (default: `plots/`) |
| `--inspect` | No | Print schema and first rows, then exit. `--mode` not required. |
| `--explore` | No | Generate raw data plots only, no sampling |

---

## Outputs

All plots and stats are saved to the `--plots` directory. The final sampled dataset is saved to `--output`.

### Statistics
| File | Description |
|---|---|
| `stats_report.txt` | Full text summary of row counts, allele counts, peptide lengths, anchor frequencies at each pipeline stage |
| `stats_summary.csv` | Numeric summary table, one row per stage (raw / phase 1 / phase 2) |

### Comparison Plots
| File | Description |
|---|---|
| `comparison_allele_distribution.png` | Peptide count per allele at each stage |
| `comparison_peptide_lengths.png` | Peptide length distribution at each stage |
| `comparison_anchor1_residues.png` | Global P2 anchor residue frequencies at each stage |
| `comparison_anchor2_residues.png` | Global last-position anchor residue frequencies at each stage |

### Per-allele Anchor Plots
| File | Description |
|---|---|
| `per_allele_anchor1_postphase1.png` | Per-allele P2 anchor AA fractions after Phase 1 (50 alleles per panel) |
| `per_allele_anchor2_postphase1.png` | Per-allele last-position anchor AA fractions after Phase 1 |
| `per_allele_anchor1_postphase2.png` | Per-allele P2 anchor AA fractions after Phase 2 |
| `per_allele_anchor2_postphase2.png` | Per-allele last-position anchor AA fractions after Phase 2 |

### KL Divergence
| File | Description |
|---|---|
| `anchor_kl_divergence_boxplot.png` | Per-allele KL divergence vs Uniform(1/20) for both anchor positions, post-phase 1 and post-phase 2. Includes reference lines for canonical dominance baselines (top-1, 2, 4, 8, 16 AA). Broken y-axis separates the bulk distribution from outliers. |
| `high_kl_alleles.csv` | Alleles where anchor dominance persists after Phase 2 — flagged when KL is above the top-16 AA baseline AND reduced by less than 5% between phases. Includes sample counts at each stage. |

### Anchor Combination
| File | Description |
|---|---|
| `anchor_combo_heatmap.png` | 20×20 heatmap of anchor combination counts (P2 × last position) at each stage |
| `anchor_combo_stats.csv` | Per-allele anchor combination diversity statistics after Phase 2 |

### Sampled Dataset
| File | Description |
|---|---|
| `<output>` | Final sampled dataset after both phases, ready for downstream use |

---

## Adapting to a New Dataset

If your input file uses different column names, update the constants at the top of `pmhc_sampling.py`:

```python
COL_PEPTIDE = "long_mer"       # column containing the peptide sequence
COL_MHC     = "allele"         # column containing the MHC allele name
COL_LABEL   = "assigned_label" # column containing the binding label (0.0 / 1.0)
MHC_CLASS   = "mhc_class"      # column containing the MHC class (1.0 / 2.0)
```

You can also adjust the following constants:

```python
SAMPLE_CAP      = 1000   # maximum peptides per allele after Phase 1
MAX_PEPTIDE_LEN = 15     # peptides of this length or longer are removed
```

---
