# MASS-PRF

**MASS-PRF** (McDonald and Kreitman with Simultaneous Synonymous-to-Replacement rate estimation) infers and quantifies natural selection across regions within a coding gene. It jointly analyzes within-species polymorphism and between-species divergence data to estimate site-specific selection coefficients using a maximum likelihood clustering framework.

**Version:** 2.0 | **Last Updated:** March 2026
**License:** Creative Commons CC BY-NC
**Reference:** Zi-Ming Zhao, Ning Li, Zhang Zhang and Jeffrey P. Townsend. (2016) Regions within coding gene sequences experience diverse intensities of natural selection inferred from polymorphism and divergence. *G3: Genes, Genomes, Genetics.*

**Tested on:**
- AWS EC2 (Ubuntu 24.04 LTS, ARM64 / x86\_64, t4g/t3 series) — functional testing of all gap policy modes and example datasets
- Yale University McCleary HPC cluster — large-scale dataset validation and runtime benchmarking

---

## Table of Contents

1. [System Requirements](#1-system-requirements)
2. [Installation](#2-installation)
3. [Quick Start Tutorial: First Example](#3-quick-start-tutorial-first-example)
4. [Input File Format](#4-input-file-format)
5. [All Options](#5-all-options)
6. [Gap Handling (–gap_policy)](#6-gap-handling)
7. [Understanding the Output](#7-understanding-the-output)
8. [Downstream Visualization: 2D and 3D Tools](#8-downstream-visualization)
9. [Common Situations and Warnings](#9-common-situations-and-warnings)
10. [Pipeline for Batch Processing](#10-pipeline-for-batch-processing)
11. [Contact](#11-contact)

---

## 1. System Requirements

### Minimum (for typical genes ≤700 bp, ≤50 sequences)

| Component | Requirement |
|-----------|-------------|
| CPU | Any modern dual-core processor |
| RAM | 2 GB |
| Disk | 100 MB (source + lookup tables) |
| OS | Linux (x86\_64 or ARM64), macOS 10.9+, or Windows via WSL |
| Compiler | g++ with C++11 support (GCC 4.8 or later) |

### Recommended (genes >700 bp or large datasets)

| Component | Requirement |
|-----------|-------------|
| CPU | 4+ cores (MASS-PRF uses multi-threading automatically) |
| RAM | 8 GB or more |
| Environment | Linux server, HPC, or cloud instance (e.g., AWS t3.medium or larger) |

### Threading behavior

MASS-PRF automatically detects the number of available CPU cores at startup and uses all of them for multi-threaded steps (confidence interval calculation and selection coefficient estimation). There is no command-line flag to limit thread count. On a shared HPC node, be aware that MASS-PRF will attempt to use all cores visible to the process.

### Scaling: when do you need more resources?

| Gene length | Sequences | Peak RAM (observed) | Typical runtime (ARM64) |
|-------------|-----------|---------------------|-------------------------|
| very short (≤200 bp) | any | < 1 GB | seconds |
| ~700 bp | ~12 | ~9 GB | ~30 min |
| ~900 bp | ~21 | ~15 GB | ~17 min |
| ~620 bp | ~34 | >15 GB | ~90 min |
| >1000 bp | >50 | 32+ GB | hours |

> Memory scales primarily with the **number of sequences** (the clustering step builds models for every unique haplotype combination). The values above are from empirical testing of the four provided example datasets on AWS EC2 (ARM64, Ubuntu 24.04 LTS). See `examples/` for expected output files (`_v200_ARM64` suffix).

---

## 2. Installation

### Required software

| Software | Purpose | Required? |
|----------|---------|-----------|
| g++ (GCC 4.8+) | Compile MASS-PRF | **Required** |
| make | Optional build automation | Optional |
| R + bio3d + Biostrings + muscle | 3D protein visualization only | Optional |
| UCSF Chimera | 3D protein visualization only | Optional |
| Python 3 + pandas + matplotlib + numpy | 2D plot generation only | Optional |

MASS-PRF itself has **no external library dependencies**. It uses only the C++ standard library (`<thread>`, `<iostream>`, `<vector>`, etc.), which ship with any standard g++ installation.

### Install the compiler

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get update && sudo apt-get install -y g++
```

**Linux (CentOS/RHEL):**
```bash
sudo yum install -y gcc-c++
```

**macOS:**
```bash
xcode-select --install
```

**Windows:** Use WSL (Windows Subsystem for Linux) and follow the Linux instructions above. Native Windows builds are not officially supported.

### Clone and compile

> **Important — which repository to clone:**
> This repository (`MASSPRF`) contains the core analysis program source code.
> A separate repository (`massprf-pipeline`) provides batch-processing scripts for running MASS-PRF across many genes. They are independent. You do not need the pipeline to run MASS-PRF.

```bash
# Clone the MASS-PRF source code
git clone https://github.com/Townsend-Lab-Yale/MASSPRF.git
cd MASSPRF

# Create output directory and compile
mkdir -p bin
g++ -std=c++0x -O3 -pthread MASSprf.cpp PRFCluster.cpp base.cpp kfunc.cpp -o bin/massprf
```

### Verify the installation

```bash
./bin/massprf -h
```

You should see the help message with version information and all available options. If you see `command not found`, confirm you are in the `MASSPRF` directory and the binary compiled successfully.

### Files in this repository

```
MASSPRF/
├── MASSprf.cpp                          # Main program entry point
├── PRFCluster.cpp / PRFCluster.h        # Core clustering and ML engine
├── base.cpp / base.h                    # Sequence utilities
├── kfunc.cpp / kfunc.h                  # Mathematical functions (gamma, beta)
├── LookupTable_*.dat                    # Pre-computed lookup tables (required at runtime)
├── Attacin-C_DmDs_pol.fas               # Example: polymorphism data (Drosophila)
├── Attacin-C_DmDs_div.fas               # Example: divergence data (Drosophila)
├── FZF1_yeast_pol.fas / _div.fas        # Example: yeast FZF1 gene
├── PanI_cod_pol.fas / _div.fas          # Example: Atlantic cod PanI gene
├── Pol_all_YIR024C.fas / Div_Spar_*.fas # Example: yeast YIR024C gene
└── bin/                                 # Compiled binary goes here
```

> The four `LookupTable_*.dat` files must remain in the same directory as the binary, or be accessible from the working directory when you run MASS-PRF.

---

## 3. Quick Start Tutorial: First Example

This tutorial walks you through running MASS-PRF from scratch using the Attacin-C example data included in the repository. No prior experience is required.

### What is this example?

- **Gene:** Attacin-C, an immunity gene in *Drosophila*
- **Polymorphism data:** 12 sequences from *Drosophila melanogaster* (`Attacin-C_DmDs_pol.fas`)
- **Divergence data:** 1 sequence from *Drosophila simulans* (`Attacin-C_DmDs_div.fas`)
- **Gene length:** 723 bp (241 codons)

### Step 1: Confirm you are in the right directory

```bash
ls examples/*.fas
```

You should see files including `Attacin-C_DmDs_pol.fas` and `Attacin-C_DmDs_div.fas`.

### Step 2: Run the analysis

```bash
./bin/massprf -p examples/Attacin-C_DmDs_pol.fas -d examples/Attacin-C_DmDs_div.fas
```

**What each flag does:**
- `-p` — the polymorphism file: multiple aligned sequences from the same species or population
- `-d` — the divergence file: one or more outgroup sequences from a related species

### Step 3: Read the output

The program prints results to the terminal. You will see lines like:

```
[GAP] policy=strict threshold=0.5

The gene length: 723 bp
PS: 18
DS: 20
PR: 4
DR: 9
```

**What do these numbers mean?**

| Symbol | Full name | Meaning |
|--------|-----------|---------|
| PS | Polymorphic Synonymous | Silent (synonymous) nucleotide differences among the ingroup sequences |
| DS | Divergent Synonymous | Silent differences between ingroup and outgroup |
| PR | Polymorphic Replacement | Amino-acid-changing differences among the ingroup sequences |
| DR | Divergent Replacement | Amino-acid-changing differences between ingroup and outgroup |

The ratio DR/DS relative to PR/PS is the core signal. When DR/DS > PR/PS, the gene (or a region of it) has experienced positive selection. MASS-PRF then clusters the gene into regions with distinct selection intensities and estimates a selection coefficient (γ) for each site.

### Step 4: Save the output to a file

```bash
./bin/massprf -p examples/Attacin-C_DmDs_pol.fas -d examples/Attacin-C_DmDs_div.fas > output_Attacin-C.txt
```

Open `output_Attacin-C.txt` in any text editor to review the full results, including per-site γ estimates, clustering boundaries, and confidence intervals.

### Step 5: Try the other example datasets

```bash
# Yeast FZF1 gene (~17 min, requires ~15 GB RAM)
./bin/massprf -p examples/FZF1_yeast_pol.fas -d examples/FZF1_yeast_div.fas > output_FZF1.txt

# Atlantic cod PanI gene (~90 min, requires >15 GB RAM)
./bin/massprf -p examples/PanI_cod_pol.fas -d examples/PanI_cod_div.fas > output_PanI.txt

# Yeast YIR024C gene (~45 seconds)
# This gene has no synonymous sites (PS=DS=0), so divergence time must be provided with -t
./bin/massprf -p examples/Pol_all_YIR024C.fas -d examples/Div_Spar_YIR024C.fas -t 5 > output_YIR024C.txt
```

### Expected runtime (v2.0, ARM64)

| Dataset | Sequences | RAM required | Runtime (ARM64) |
|---------|-----------|--------------|-----------------|
| YIR024C | 19 | < 5 GB | ~45 seconds |
| Attacin-C | 12 | ~9 GB | ~30 min |
| FZF1 | 21 | ~15 GB | ~17 min |
| PanI | 34 | >15 GB | ~90 min |

Tested on AWS EC2 ARM64 (Ubuntu 24.04 LTS). Runtime depends on number of sequences and gene length. The installation procedure is identical on x86\_64 and ARM64.

Expected output files for all four datasets are provided in `examples/` (suffixed `_v200_ARM64` for v2.0 ARM64).

---

## 4. Input File Format

Both input files must be in **standard FASTA format** with sequences pre-aligned to the same length.

```
>Sequence_Name_1
ATGAGCAAAACTGTTCTCCTAATT...
>Sequence_Name_2
ATGAGCAAAATTGTTCTCCTAATT...
```

**Rules:**
- All sequences in each file must be the **same length**
- Sequences must be **in-frame coding sequences** starting at the first codon position
- MASS-PRF does **not** perform alignment — use a tool such as MAFFT, MUSCLE, or ClustalW first
- Gap characters (`-`) are accepted and handled according to `-gap_policy`
- Ambiguous bases (non-A/T/G/C) are handled according to `-n`
- Each file can have one blank line at the end without causing an error (you will see a notice, which is safe to ignore)

**Polymorphism file (`-p`):** multiple sequences from the same species or population (ingroup).
**Divergence file (`-d`):** one outgroup sequence, or a consensus of outgroup sequences.

---

## 5. All Options

```
./bin/massprf -p <pol_file> -d <div_file> [options]
```

### Input / Output

| Flag | Description | Default |
|------|-------------|---------|
| `-p` | Polymorphism FASTA file | **required** |
| `-d` | Divergence FASTA file | **required** |
| `-ic` | Input format: `0`=raw DNA sequences, `1`=pre-built consensus | `0` |
| `-sn` | Number of polymorphism sequences (only needed when `-ic 1`) | — |
| `-o` | Output level: `0`=amino acid, `1`=nucleotide | `0` |
| `-v` | Verbose output: `0`=minimal, `1`=full | `1` |

### Clustering

| Flag | Description | Default |
|------|-------------|---------|
| `-c` | Criterion: `0`=BIC, `1`=AIC, `2`=AICc, `3`=LRT | `0` (BIC) |
| `-s` | Include synonymous site clustering: `0`=no, `1`=yes | `1` |
| `-m` | `0`=model selection + averaging, `1`=selection only | `0` |
| `-ci_m` | Compute 95% CI for model averaging: `0`/`1` | `0` |
| `-a` | Minimum gene length to trigger scaling | auto |
| `-SCL` | Compress sequences by 3, 6, or 9 nucleotides | off |

> **Recommendation:** Use BIC (`-c 0`) for datasets with many sites or sequences. Use AIC or AICc for smaller sample sizes. Keep `-s 1` (default) for improved accuracy.

### Selection Coefficient Estimation

| Flag | Description | Default |
|------|-------------|---------|
| `-r` | Estimate selection coefficient per site: `0`/`1` | `1` |
| `-ci_r` | Compute 95% CI for selection coefficient: `0`/`1` | `1` |
| `-exact` | CI algorithm: `0`=stochastic, `1`=exact | `0` |
| `-mn` | Number of models for stochastic algorithm | `10000` |
| `-NI` | Estimate Neutrality Index per site: `0`/`1` | `0` |
| `-rMAp` | Output gamma from model-averaged PR and DR: `0`/`1` | `0` |

### Divergence Time

| Flag | Description | Default |
|------|-------------|---------|
| `-t` | Species divergence time in million years (MY) | auto-estimated |
| `-ssd` | Use site-specific divergence time from silent-site clustering | off |

> `-t` and `-ssd` cannot be used together.

### Sequence Handling

| Flag | Description | Default |
|------|-------------|---------|
| `-g` | Genetic code (NCBI table number) | `1` (Standard) |
| `-n` | Non-A/T/G/C base treatment: `0`=treat as gap, `1`=replace with most frequent base | `0` |
| `-ng` | Non-genic mode (every base treated as nonsynonymous): `0`/`1` | `0` |
| `-gap_policy` | Gap handling: `strict`, `pairwise`, or `threshold` | `strict` |
| `-gap_threshold` | Gap fraction cutoff for threshold mode (0.0–1.0) | `0.5` |

### Other

| Flag | Description |
|------|-------------|
| `-h` | Show help and exit |

---

## 6. Gap Handling

Aligned sequences often contain gap characters (`-`) where a sequence has an insertion or deletion relative to the alignment. MASS-PRF offers three strategies for how to treat codons that contain gaps.

| Policy | Rule | When to use |
|--------|------|-------------|
| `strict` (default) | Exclude any codon where **any** sequence has a gap | Most conservative. Use when gap-containing sequences are unreliable or rare. |
| `pairwise` | Keep a codon if at least **2 sequences** have no gap at that position | Retains more data. Appropriate for most standard alignment quality levels. |
| `threshold` | `pairwise` rule, plus additionally exclude if the fraction of gapped sequences exceeds `-gap_threshold` | Fine-grained control. Useful when some positions are highly gapped but others are not. |

The active gap policy is always printed at the start of every run:
```
[GAP] policy=pairwise threshold=0.5
```

This line confirms which policy was applied and allows results to be reproduced exactly.

### Examples

```bash
# Default: strict
./bin/massprf -p pol.fas -d div.fas

# Pairwise: keep codon if 2+ sequences are gap-free
./bin/massprf -p pol.fas -d div.fas -gap_policy pairwise

# Threshold: pairwise + reject if more than 30% of sequences have a gap
./bin/massprf -p pol.fas -d div.fas -gap_policy threshold -gap_threshold 0.3
```

### Interaction with `-N` flag

The legacy `-N` flag (imputation mode) implicitly sets gap policy. If `-gap_policy` is also specified, the explicit `-gap_policy` takes priority.

---

## 7. Understanding the Output

### Key summary lines

```
The gene length: 723 bp
PS: 18    # polymorphic synonymous sites
DS: 20    # divergent synonymous sites
PR: 4     # polymorphic replacement (nonsynonymous) sites
DR: 9     # divergent replacement (nonsynonymous) sites
```

### Selection coefficient (γ) per site

After clustering, MASS-PRF reports per-site γ values. Positive γ indicates positive selection; negative γ indicates purifying selection; γ near 0 indicates neutrality.

### Clustering boundaries

The output shows which regions of the gene were assigned to the same selection class, with log-likelihood scores for model selection.

### "Mission accomplished"

The final line of a successful run is:
```
Mission accomplished.
```

If you do not see this line, the run did not complete successfully.

---

## 8. Downstream Visualization

MASS-PRF output can be visualized in two ways using companion tools.

### 2D Plot (per-gene summary)

The `2D_process_massprf_res.py` Python script processes MASS-PRF result files and generates 2D selection intensity plots.

**Requirements:**
```bash
pip install pandas matplotlib numpy
```

**Usage:** Place the script in the directory containing your MASS-PRF output files and run:
```bash
python3 2D_process_massprf_res.py
```

Plots and processed results are saved to `./processed_output/`.

### 3D Protein Structure (Chimera coloring)

The `batchMASSPRF_To_Chimera.R` R script maps MASS-PRF selection coefficients onto a protein structure and generates coloring commands for UCSF Chimera.

**Requirements:**
- R with packages: `bio3d`, `Biostrings`, `muscle`
- UCSF Chimera (free for academic use)

**Usage in R:**
```r
source("batchMASSPRF_To_Chimera.R")

# Single gene
batchMASSPRF_Chimera(
  designFile = "example_inputs/S-design.tsv",
  hasHeader  = TRUE,
  onlySig    = FALSE,
  bins       = 10
)
```

The design file (TSV format) specifies the PDB structure file, the nucleotide FASTA, the MASS-PRF output table, and the scaling factor. Example design files are provided in `example_inputs/`.

### Workflow summary

```
Raw sequences (FASTA)
        │
        ▼
   [Alignment tool: MAFFT / MUSCLE / ClustalW]
        │
        ▼
   MASS-PRF  (this tool)
        │
        ├──► 2D plot via 2D_process_massprf_res.py
        │
        └──► 3D structure coloring via batchMASSPRF_To_Chimera.R + UCSF Chimera
```

---

## 9. Common Situations and Warnings

### Using non-genic mode (`-ng 1`)

Non-genic mode treats every nucleotide position as a replacement site, bypassing synonymous/nonsynonymous classification. This is intended for non-coding sequences (e.g., introns, UTRs) where the standard codon-based distinction does not apply.

**Required flags when using `-ng 1`:**

```bash
./bin/massprf -p pol.fas -d div.fas -ng 1 -s 0 -t <divergence_time>
```

- `-s 0` — synonymous clustering must be disabled (there are no synonymous sites)
- `-t <value>` — divergence time in million years must be provided, since it cannot be auto-estimated without synonymous sites

**Example:**
```bash
./bin/massprf -p ncds_pol.fas -d ncds_div.fas -ng 1 -o 1 -s 0 -n 0 -t 5
```

### "Warning: There is not enough information from synonymous sites to estimate divergence time!"

This warning appears when PS=0 and DS=0 (the gene has no synonymous polymorphic or divergent sites). The program cannot auto-estimate the divergence time `r` from silent sites.

**What to do:** Provide the divergence time manually with `-t <value>` (in million years). For example, for *S. cerevisiae* vs *S. paradoxus*, use `-t 5`. Without `-t`, the program will crash after this warning.

### "Warning: PR is too low for clustering, skipping PR clustering."

This means the dataset has fewer than 2 nonsynonymous polymorphic sites (PR < 2). The clustering step requires at least 2 events to statistically distinguish regions. This is not a program error — it is expected for highly conserved genes or small datasets. The program continues and outputs the PS/DS/PR/DR counts it did compute.

**What to do:** This is normal for genes under strong purifying selection. If you expected more variation, check that your polymorphism file contains enough sequences and that the sequences are correctly aligned in-frame.

### "Warning: DR is too low for clustering, skipping DR clustering."

Same situation as above, but for divergent replacement sites. The program continues.

### "Error clustering (PS is too low!)" or "Error clustering (DS is too low!)"

These errors indicate PS or DS < 2. They occur inside the clustering engine when the clustering was already attempted. Increasing your dataset size or using a gene with more variation is the remedy.

### "Notice: End-of-file reached."

This is a harmless notice that appears if your FASTA file has a trailing blank line. You can ignore it.

### The program is killed (exit code 137)

This is an out-of-memory (OOM) condition. The system ran out of RAM.

**Solutions:**
- Use a machine or instance with more RAM (see System Requirements)
- Add swap space on Linux: `sudo fallocate -l 4G /swapfile && sudo chmod 600 /swapfile && sudo mkswap /swapfile && sudo swapon /swapfile`
- For genes >1000 bp, use an HPC node or a cloud instance with ≥8 GB RAM

### Results are identical across gap policy modes

This is expected when gap-containing positions do not overlap with polymorphic sites, or when the gap fraction is low in all gap-containing positions. The gap policy difference only affects results when gapped codons coincide with non-synonymous polymorphic sites in sufficient numbers.

### Multi-threading on shared HPC nodes

MASS-PRF uses all CPU cores it can detect (`std::thread::hardware_concurrency()`). On a shared HPC node, this may conflict with other jobs. Submit your job with a CPU allocation that matches your node's available cores, or consult your HPC system administrator about core binding.

---

## 10. Pipeline for Batch Processing

For analyzing many genes automatically, a separate pipeline repository is available:

```bash
git clone https://github.com/Townsend-Lab-Yale/massprf-pipeline.git
```

> This is a **distinct repository** from the MASS-PRF source code. It provides shell/Python scripts to run MASS-PRF across directories of gene files. Refer to its own README for installation and usage instructions. The pipeline requires MASS-PRF to already be compiled and accessible.

---

## 11. Contact

For bug reports or suggestions:

- Yide Jin（Maintainer） — jinyide0202@gmail.com
- Zi-Ming Zhao（Creator） — ziming.gt@gmail.com
- Michael C. Campbell — mc44680@usc.edu
- Jeffrey Townsend — Jeffrey.Townsend@Yale.edu

More information: https://medicine.yale.edu/lab/townsend/sand/
