## 2D Processing: MASS-PRF Output → CSV, Site Lists, and Plots

## 2D Usage 

### 1) Prepare inputs
Place your MASS-PRF `*.txt` result files in the **same folder** as the 2D script.

The scripts look for:
- A completion line like **`Mission accomplished`**.  
- A results table **between** the header line and the abbreviation block.  
  The **original** script expects the exact header & footer lines (string match).  fileciteturn2file1  
  The **new-fixed** script locates the table using **anchor strings** with fallback logic (more tolerant).  fileciteturn2file0

### 2) Run
```bash
# Original (baseline)
python3 2D_process_massprf_res.py

# New (fixed)
python3 "2D_new_fixed (3).py"
```

### 3) Outputs (auto-created under `./processed_output/`)

**Per-gene CSV:** `<GENE>.csv` with at least  
`Position, DivergentTime, Gamma, Lower_CI_Gamma, Upper_CI_Gamma`  (new-fixed guarantees per-NT Position from 1..N)  fileciteturn2file0

**Site lists:**  
- `<GENE>_pos_sites.txt` — `Gamma > 4` **and** `Lower_CI_Gamma > 0`  
- `<GENE>_str_pos_sites.txt` — `Gamma > 4` **and** `Lower_CI_Gamma > 4`  
- `<GENE>_neg_sites.txt` — `Gamma < -1` **and** `Upper_CI_Gamma < 0`  
(Thresholds are identical in both versions.)  fileciteturn2file1

**Summary lists:** `Positive_genes.txt`, `Strongly_positive_genes.txt`, `Negative_genes.txt`, `Boring_genes.txt`, `Failed_genes.txt` (both versions).  fileciteturn2file1

**Plots:** `./processed_output/plots/<GENE>.pdf`  
- X-axis: position (original labels "Nucleotide position"; new-fixed uses generic "position")  
- Y-axis: scaled selection coefficient (gamma); shaded CI band; 0-reference line  fileciteturn2file1turn2file0

### 4) Hand-off to 3D
- The **CSV** is the numeric backbone for coloring (continuous ramp by `Gamma` or CI-aware bins).  
- The **site lists** act as masks to select/highlight residues (e.g., strongly positive).  
- If your 3D workflow needs **amino-acid numbering**, map nucleotide positions to residues (NT→codon→AA) before coloring.

---

## What’s updated (new-fixed vs original)

| Area | Original: `2D_process_massprf_res.py` | New-fixed: `2D_new_fixed (3).py` |
|---|---|---|
| **Robustness of parsing** | Relies on **exact string match** for the header line and a very long **exact** abbreviation line; may fail if upstream print format changes. fileciteturn2file1 | Uses **anchor substrings** and a **fallback** (first blank line) to find the table; more tolerant to format/whitespace differences. fileciteturn2file0 |
| **Scaling detection** | Requires seeing “Scaling by supplied factor of X” or “No Scaling requested”; otherwise **raises error** and stops. fileciteturn2file1 | Extracts scaling when present; if missing, **defaults to 1** to avoid hard failure. fileciteturn2file0 |
| **Nucleotide expansion** | Expands to NT **only when needed**; position rebuilt 1..N in code blocks that handle pos/neg cases separately. fileciteturn2file1 | **Always** uses a single `descale_to_nt()` to expand to per-NT grid, then **uniformly** rebuilds `Position = 1..N` for downstream consistency. fileciteturn2file0 |
| **Type coercion & columns** | Coerces `Gamma`/CIs; table parsing via simple `split()`; column set depends on header line tokens. fileciteturn2file1 | Centralized `parse_table()` coerces numerics for `Position, DivergentTime, Gamma, Lower_CI_Gamma, Upper_CI_Gamma` when present; trims/normalizes row tokens. fileciteturn2file0 |
| **File name → gene id** | `gene_name = name[0] + '_' + name[1]` (first two underscore tokens only). fileciteturn2file1 | `gene_name = basename(entry).split('.')[0]` (keeps full stem; avoids accidental truncation). fileciteturn2file0 |
| **Plotting** | Plots in a **second pass** over CSVs; X label is `"Nucleotide position"`. fileciteturn2file1 | Plots **immediately** per gene after CSV write; X label is generic `"position"` to be NT/AA-agnostic. fileciteturn2file0 |
| **Failure handling** | Prints “Failed: …” when missing “Mission accomplished”; may stop on scaling not found. fileciteturn2file1 | Logs to `Failed_genes.txt` and **continues**; scaling missing falls back to 1; parsing has defensive fallbacks. fileciteturn2file0 |

**Bottom line:** The new version is more robust, more tolerant of variations in upstream output, and produces more consistent results — reducing the risk of crashes or false “failed run” detections caused by minor formatting differences.

---

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

- If you use **scaling** upstream, the new-fixed script will **always** expand to per-NT; this makes 3D stage simpler and standardized. fileciteturn2file0  
- Both versions share the **same thresholds** for site classification; differences in gene lists are usually due to **parsing robustness** and **position expansion** differences rather than logic changes. fileciteturn2file1turn2file0


## Acknowledgment

We thank **Dr. Hayley Hassler** for foundational 2D visualization scripts that inspired this tool.
