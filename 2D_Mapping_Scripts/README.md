**Update by Yide Jin at 2025-10-15 **  

## 2D Processing: MASS-PRF Output → CSV, Site Lists, and Plots

## 2D Usage 

### 1) Prepare inputs
Place your MASS-PRF `*.txt` result files in the **same folder** as the 2D script.

The scripts look for:
- A completion line like **`Mission accomplished`**.  
- A results table **between** the header line and the abbreviation block.  
  The **original** script expects the exact header & footer lines (string match).    
  The **new-fixed** script locates the table using **anchor strings** with fallback logic (more tolerant).  

### 2) Run
```bash

python3 2D_process_massprf_res.py
MASS-PRF 2D processing is the bridge between the **core MASS-PRF** run and the **3D structure-mapping** step. It reads per-gene plain‑text outputs from MASS-PRF, verifies completion, parses the results table, optionally expands scaled sites to **per‑nucleotide resolution**, classifies selection at the site and gene level, and emits clean, machine‑readable files (CSV and site lists) plus publication‑ready **2D Gamma plots**. These outputs are then used by the 3D mapping scripts to color protein structures in UCSF Chimera/ChimeraX. Use this README as a practical guide to run the 2D script, understand its outputs, and hand them off to the 3D stage.


```

### 3) Outputs (auto-created under `./processed_output/`)

**Per-gene CSV:** `<GENE>.csv` with at least  
`Position, DivergentTime, Gamma, Lower_CI_Gamma, Upper_CI_Gamma`  (new-fixed guarantees per-NT Position from 1..N)  

**Site lists:**  
- `<GENE>_pos_sites.txt` — `Gamma > 4` **and** `Lower_CI_Gamma > 0`  
- `<GENE>_str_pos_sites.txt` — `Gamma > 4` **and** `Lower_CI_Gamma > 4`  
- `<GENE>_neg_sites.txt` — `Gamma < -1` **and** `Upper_CI_Gamma < 0`  
(Thresholds are identical in both versions.)  

**Summary lists:** `Positive_genes.txt`, `Strongly_positive_genes.txt`, `Negative_genes.txt`, `Boring_genes.txt`, `Failed_genes.txt` (both versions).  

**Plots:** `./processed_output/plots/<GENE>.pdf`  
- X-axis: position (original labels "Nucleotide position"; new-fixed uses generic "position")  
- Y-axis: scaled selection coefficient (gamma); shaded CI band; 0-reference line  

### 4) Hand-off to 3D
- The **CSV** is the numeric backbone for coloring (continuous ramp by `Gamma` or CI-aware bins).  
- The **site lists** act as masks to select/highlight residues (e.g., strongly positive).  
- If your 3D workflow needs **amino-acid numbering**, map nucleotide positions to residues (NT→codon→AA) before coloring.

### Practical tips & troubleshooting

- **File naming**: keep each result as one `*.txt` per gene/region. Avoid duplicate stems; the script uses the file stem as the gene ID.
- **Where are my outputs?** If a gene is missing its CSV or plot, check `Failed_genes.txt` for the reason (no completion flag; table not found).
- **Scaling messages**: if MASS-PRF prints `Scaling by supplied factor of X`, the script expands each row **X-fold**; if it prints `No Scaling requested`, positions remain 1:1.
- **CI-based filters**: “positive/strongly positive/negative” rely on both `Gamma` **and** its CI bound. Borderline sites may be excluded by design.
- **Reproducibility**: record upstream flags (`gap_policy`, `gap_threshold`, `-n`) with your outputs. They affect which codons/sites contribute to estimates.
- **Large batches**: processing is per‑file; you can parallelize by splitting input files into subfolders and running multiple instances.
- **Plot sanity checks**: expect a horizontal zero line and a CI band; sudden position jumps usually indicate missing scaling expansion upstream.

### Hand-off checklist for 3D

- Have `<GENE>.csv` and any site lists (`_str_pos_sites.txt`, `_pos_sites.txt`, `_neg_sites.txt`) ready.
- If your 3D workflow colors **amino‑acid residues**, prepare a nucleotide→codon→residue map so that `Position` aligns with residue numbering.
- For batch jobs, first filter genes by `Positive_genes.txt` / `Strongly_positive_genes.txt` to prioritize interesting cases.


---


## What’s updated

- The script is now **more robust to upstream format changes**: results table detection uses tolerant anchors and fallbacks; missing scaling defaults to 1 instead of stopping.
- **Per‑nucleotide expansion is standardized**: after parsing, positions are rebuilt as a continuous `1..N` grid to keep downstream 3D mapping consistent.
- **Failure handling is cleaner**: files without a completion flag or with unparsable tables are logged to `Failed_genes.txt`, while the rest continue to process.
- **Output consistency**: per‑gene CSV columns are normalized (`Position, DivergentTime, Gamma, Lower_CI_Gamma, Upper_CI_Gamma` when present); plot labels are generic enough for NT or AA contexts.
- **Filename handling is safer**: gene identifiers are derived from the full basename (before `.txt`) to avoid accidental truncation when underscores are present.


## Minimal example layout

```
./
├── 2D_process_massprf_res.py        
├── *.txt                            # MASS-PRF outputs (inputs)
└── processed_output/
    ├── <GENE>.csv
    ├── <GENE>_pos_sites.txt
    ├── <GENE>_str_pos_sites.txt
    ├── <GENE>_neg_sites.txt
    ├── Positive_genes.txt
    ├── Strongly_positive_genes.txt
    ├── Negative_genes.txt
    ├── Boring_genes.txt
    ├── Failed_genes.txt
    └── plots/
        └── <GENE>.pdf
```

---

## Requirements

- Python 3  
- Packages: `pandas`, `matplotlib`, `numpy` (plus stdlib `os`, `re`)

Install:
```bash
pip install pandas matplotlib numpy
```

---

## Notes & tips

- If you use **scaling** upstream, the new-fixed script will **always** expand to per-NT; this makes 3D stage simpler and standardized.   
- Both versions share the **same thresholds** for site classification; differences in gene lists are usually due to **parsing robustness** and **position expansion** differences rather than logic changes. 


## Acknowledgment

We thank **Dr. Hayley Hassler** for foundational 2D visualization scripts that inspired this tool.
