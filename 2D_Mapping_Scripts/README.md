# 2D Visualization: MASS-PRF Output → CSV, Site Lists, and Plots

**Updated by Yide Jin — March 2026**

This script reads plain-text output files from MASS-PRF, verifies that each run completed successfully, parses the per-site results table, expands scaled positions to per-nucleotide resolution, classifies sites by selection intensity, and produces clean CSV files, site lists, and publication-ready 2D plots of the scaled selection coefficient (γ) across each gene.

These outputs feed directly into the 3D structure-mapping step (see `batchMASSPRF_To_Chimera.R`).

---

## Table of Contents

1. [Requirements](#1-requirements)
2. [Quick Start](#2-quick-start)
3. [What the Script Does](#3-what-the-script-does)
4. [Outputs](#4-outputs)
5. [Understanding the Plots](#5-understanding-the-plots)
6. [Understanding the Site Lists](#6-understanding-the-site-lists)
7. [Hand-off to 3D Visualization](#7-hand-off-to-3d-visualization)
8. [Troubleshooting](#8-troubleshooting)
9. [Acknowledgments](#9-acknowledgments)

---

## 1. Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| Python | 3.6+ | Run the script |
| pandas | any recent | Table parsing |
| matplotlib | any recent | Plot generation |
| numpy | any recent | Numeric operations |

Install dependencies:

```bash
pip install pandas matplotlib numpy
```

No other dependencies are required. The script uses only the Python standard library beyond these three packages.

---

## 2. Quick Start

### Step 1: Collect your MASS-PRF output files

Place all MASS-PRF result `.txt` files in the **same folder** as the script. Each file should be the direct output of one MASS-PRF run (one gene per file).

```
your_working_folder/
├── 2D_process_massprf_res.py
├── GeneA_results.txt
├── GeneB_results.txt
└── GeneC_results.txt
```

### Step 2: Run the script

```bash
python3 2D_process_massprf_res.py
```

That is all. No flags or arguments are needed.

### Step 3: Check the outputs

All outputs are written to `./processed_output/`:

```
processed_output/
├── GeneA.csv
├── GeneA_pos_sites.txt
├── GeneB.csv
├── Positive_genes.txt
├── Negative_genes.txt
├── Boring_genes.txt
├── Failed_genes.txt
└── plots/
    ├── GeneA.pdf
    └── GeneB.pdf
```

---

## 3. What the Script Does

For each `.txt` file in the working directory, the script:

1. **Checks for completion** — looks for the `Mission accomplished.` marker. Files without it are logged to `Failed_genes.txt` and skipped.

2. **Parses the results table** — locates the per-site results block using robust anchor detection. The table contains one row per site with columns including `Position`, `Gamma`, `Lower_CI_Gamma`, and `Upper_CI_Gamma`.

3. **Detects and applies scaling** — if the MASS-PRF run used the `-SCL` flag (sequence compression), the script reads the scaling factor from the output and expands each row back to per-nucleotide resolution. This ensures that position numbers in the CSV and plots always correspond to actual nucleotide positions. If no scaling was used, positions remain 1:1.

4. **Classifies sites** — applies fixed γ thresholds to label sites as positively selected, strongly positively selected, or negatively selected (see [Site Lists](#6-understanding-the-site-lists)).

5. **Writes outputs** — per-gene CSV, site list text files, a 2D γ plot (PDF), and summary gene lists.

---

## 4. Outputs

### Per-gene CSV (`processed_output/<GENE>.csv`)

One row per nucleotide position with the following columns (when present in the MASS-PRF output):

| Column | Description |
|--------|-------------|
| `Position` | Nucleotide position (1-indexed, continuous) |
| `DivergentTime` | Estimated divergence time at this site |
| `Gamma` | Model-averaged scaled selection coefficient |
| `Lower_CI_Gamma` | Lower bound of 95% confidence interval |
| `Upper_CI_Gamma` | Upper bound of 95% confidence interval |

This CSV is the primary input for the 3D coloring step.

### Site list files

| File | Sites included |
|------|---------------|
| `<GENE>_pos_sites.txt` | Positively selected: γ > 4 and lower CI > 0 |
| `<GENE>_str_pos_sites.txt` | Strongly positively selected: γ > 4 and lower CI > 4 |
| `<GENE>_neg_sites.txt` | Negatively selected: γ < −1 and upper CI < 0 |

Each file lists one nucleotide position per line. Only created when at least one qualifying site exists.

### Summary gene lists (`processed_output/`)

| File | Contents |
|------|----------|
| `Positive_genes.txt` | Genes with at least one positively selected site |
| `Strongly_positive_genes.txt` | Genes with at least one strongly positively selected site |
| `Negative_genes.txt` | Genes with at least one negatively selected site |
| `Boring_genes.txt` | Genes with no significantly selected sites |
| `Failed_genes.txt` | Files that could not be processed (missing completion marker or unparsable table) |

### Plots (`processed_output/plots/<GENE>.pdf`)

One PDF per gene. See [Understanding the Plots](#5-understanding-the-plots).

---

## 5. Understanding the Plots

Each plot shows γ (scaled selection coefficient) along the gene:

- **X-axis:** nucleotide position (1 to gene length)
- **Y-axis:** γ value (model-averaged)
- **Shaded band:** 95% confidence interval around γ
- **Horizontal line at γ = 0:** neutral expectation

**How to read it:**

| γ value | Interpretation |
|---------|---------------|
| γ >> 0 | Positive selection (site evolves faster than neutral) |
| γ ≈ 0 | Neutral evolution |
| γ << 0 | Purifying (negative) selection (site is conserved) |

Sites where the entire CI band lies above or below zero are the most statistically confident calls.

---

## 6. Understanding the Site Lists

The classification thresholds are:

| Category | Gamma threshold | CI threshold | Biological meaning |
|----------|----------------|--------------|-------------------|
| Positively selected | γ > 4 | lower CI > 0 | Confidently above neutral |
| Strongly positively selected | γ > 4 | lower CI > 4 | CI entirely above neutral |
| Negatively selected | γ < −1 | upper CI < 0 | Confidently below neutral |

These thresholds are fixed in the script. Sites that do not meet either a positive or negative criterion are counted as "boring" at the gene level.

---

## 7. Hand-off to 3D Visualization

The 3D coloring step (`batchMASSPRF_To_Chimera.R`) reads the per-gene CSV files produced here.

**Checklist before running the 3D step:**

- Confirm `<GENE>.csv` exists in `processed_output/` for each gene of interest.
- If your 3D workflow colors **amino acid residues**, you will need to map nucleotide positions to codon/residue numbers (divide by 3, accounting for reading frame).
- Use `Positive_genes.txt` or `Strongly_positive_genes.txt` to prioritize genes with interesting selection signatures.
- The `_str_pos_sites.txt` files are the most conservative input for highlighting residues in Chimera.

---

## 8. Troubleshooting

### A gene appears in `Failed_genes.txt`

The most common reasons:

1. **No `Mission accomplished.` line** — the MASS-PRF run did not finish. Re-run MASS-PRF on that gene.
2. **Table not found** — the output format may be unusual. Open the `.txt` file and confirm it contains a results table between the `//Results based on model averaging` header and the `Abbreviation:` footer.
3. **Encoding issue** — the script reads files with `errors='ignore'`. If a file is binary or severely corrupted, it will still fail.

### The plot shows a flat line or all zeros

This usually means the results table was parsed but γ values are all zero or missing. Check that the MASS-PRF run completed the selection coefficient estimation step (look for the `r_stochastic` section in the output).

### Positions look wrong (e.g., jumps in the x-axis)

This can happen if MASS-PRF was run with `-SCL` (sequence compression) and the scaling factor is not being detected. Open the `.txt` file and search for `Scaling by supplied factor of` or `No Scaling requested`. If neither line is present, scaling detection will default to 1 (no expansion).

### The script runs but produces no output

Make sure you are running the script from the folder that contains the `.txt` files, not from a parent or child directory. The script scans only the current working directory.

---

## 9. Acknowledgments

We thank **Dr. Hayley Hassler** for foundational 2D visualization scripts that inspired this tool.
