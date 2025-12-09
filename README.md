# MASS-PRF: Model Averaged Site Selection via Poisson Random Field 
Updated in 2025/12/9 by Yide Jin

## Overview
MASS-PRF is a computational tool designed to detect regional variation in selection intensity within protein-coding genes using DNA sequence polymorphism and divergence data. This repository includes the program, preprocessing scripts, and a pipeline for genome-wide analysis.

---
## Updates in Version 2.0 (December 9, 2025)

- **Non-genic support — `-ng 1` (default 0)**
  - **What it does**
    - Turns on per-nucleotide mode for non-coding regions (intergenic / intronic / UTR).
    - Disables codon / synonymous logic.
  - **Output changes (by design)**
    - Every variable nucleotide site is counted as replacement (R).
    - PS/DS (synonymous polymorphism/divergence) are not computed → effectively PS = 0 and DS = 0.
    - MKT and NI are skipped (reported as NA) because S = 0 makes them undefined.
    - PR/DR and Gamma are still computed as usual.
  - **Backwards compatibility**
    - With `-ng 0` (the default), behavior is unchanged.
  - **Flag interactions**
    - `-s 1` and `-ssd` have no effect in `-ng 1` (no synonymous sites / clusters); they are ignored safely.
    - Provide a species divergence time from coding runs (or an external estimate) with `-t <div_time>` (same as `Div_time`).
      - This matches the two-step workflow Jeff requested: estimate `dt` in coding → reuse in non-genic.
  - **Ambiguous bases & gaps (same knobs as before)**
    - `-n 0/1/2` (ambiguous base handling):
      - `0` — treat non-ATGC as missing at that site (skip). *(default)*
      - `1` — replace ambiguous bases by the site majority (A/T/G/C).
      - `2` — replace both gaps and ambiguous bases by the site majority  
        (equivalent to enabling a majority rule for gaps; internally aligns with `gap_policy = 2`).
  - **Gap handling (unchanged behavior)**
    - Default behavior matches previous releases.
    - You can optionally override via environment variables:
      ```bash
      export GAP_POLICY=0|1|2    # 0=skip site if any gap; 1=keep site with ≥2 gap-free seqs; 2=majority rule with threshold
      export GAP_THRESHOLD=0.5   # used only when GAP_POLICY=2
      ```
    - On startup, MASS-PRF prints the active settings, e.g. `[GAP] policy=2 threshold=0.7`.
  - **Recommended for non-genic**
    - `-o 1` → nucleotide-level output.
    - `-t <dt>` → pass the divergence time estimated from coding CDS.
  - **Examples**
    ```bash
    # Coding run (unchanged)
    massprf -p pol.fasta -d div.fasta -s 1 -m 0 -o 1

    # Non-genic run (per-nucleotide, reuse dt from coding)
    massprf -p pol_ng.fasta -d div_ng.fasta -ng 1 -n 1 -t 3.2 -o 1
    ```

- **Stability fix (no change to results)**
  - **What was fixed**
    - Rare runs could crash with:
      ```text
      blew fast model num ... 0.999193
      ```
    - This came from a floating-point edge case when drawing models by the cumulative weights (inverse-CDF sampling).
  - **What we changed (implementation detail)**
    - Normalize the cumulative weights so the last entry is exactly `1.0` (valid CDF).
    - Draw the random number `u` in `[0,1)` and select with a standard `lower_bound` search.
  - **Why this is safe**
    - This is the textbook inverse-CDF method.
    - It does **not** change model probabilities or downstream estimates.
    - It only removes the boundary crash; behavior is otherwise identical and default pipelines are unaffected.
  - **2D script notes (compatibility)**
    - 2D now accepts the scaling factor from MASS-PRF header lines (e.g. `Scaling by supplied factor of 3 / Scaled size: …`) **or** via CLI `--scale <k>`.
    - Added a `--non_genic` switch so PS/DS-based panels are skipped when input comes from `-ng 1`.
    - If you don’t use `-ng 1`, you don’t need any new flags—plots match the old behavior.


## Updates in Version 1.32 (September 19, 2025)

- **New gap-handling controls**
  - `gap_policy` (controls how alignment gaps “–” are treated at each **site**):
    - `0` — **Skip gapped sites** (default): if any sequence has a gap at the site, the site is excluded from analysis.  
    - `1` — **Keep site, drop gapped sequences**: analyze the site using only sequences without a gap at that position (**requires ≥ 2 sequences with valid nucleotides**).  
    - `2` — **Majority-rule replacement**: replace gaps using the consensus nucleotide **only if** the non-gap majority at the site is ≥ `gap_threshold`; otherwise exclude the site.
  - `gap_threshold` (float, default `0.5`): used **only** when `gap_policy = 2`; sets the required non-gap majority fraction for replacement.
  - Both settings can be provided via environment variables:
    ```bash
    export GAP_POLICY=2 GAP_THRESHOLD=0.7
    ```
  - On startup the program now reports the active settings, e.g. `[GAP] policy=2 threshold=0.7`.

- **Revised default for `-n` (ambiguous bases)**
  - **Old default:** `-n 1` (replace ambiguous, non-ATGC nucleotides).  
  - **New default:** `-n 0` (treat ambiguous nucleotides as missing; exclude the site).  
  - If ambiguous bases are common and you wish to impute them, specify `-n 1` (majority replacement) or `-n 2` (apply the same majority rule to **both** ambiguous bases **and** gaps, equivalent to `gap_policy = 2` for gaps).

- **Performance**
  - Automatic CPU core detection (`NUM_CPU`) enables faster parallel execution across platforms.


## Pipeline Overview

- **This stage (MASS-PRF core):**  
  Run `massprf` per gene/region to produce plain-text results (`*.txt`).  
  Core options such as `gap_policy`, `gap_threshold`, and the default `-n 0` affect which sites/codons contribute to the results table.

- **Downstream (2D processing):**  
  - Read `*.txt` and confirm the run finished successfully  
  - Parse the results table  
  - Expand to per-nucleotide positions when a scaling factor was used  
  - Classify selection and save **tidy CSV** + **site lists**  
  - Generate **2D Gamma + CI** plots (PDF)

- **Downstream (3D mapping):**  
  Use the per-gene CSVs and site lists to color protein structures (UCSF Chimera/ChimeraX).  
  Typical usage: color by `Gamma` and highlight sites from the positive/negative lists.  
  See `3D_Mapping_Scripts/` for examples.

> Note: The current 2D script name reflects the **new fixed** implementation, while keeping the original filename (`2D_process_massprf_res.py`) for backward compatibility.
---

## Installation

### Prerequisites
The pipeline assumes a Unix environment with bash shell. Advanced users can adjust instructions for other environments.

---

### Steps

#### 0) Install MASS-PRF and MASS-PRF Preprocess
You can clone this repository and read within for build instructions:  
[https://github.com/zimingz/MASSPRF_10July2016](https://github.com/zimingz/MASSPRF_10July2016)

It is important to note that you may need to build both `massprf` and `massprf_preprocess` independently:
- The final executable for `massprf` will be in `./bin`.
- The `massprf_preprocess` executable will be in `./MASSPRF_preprocessing_08July2016`.

---

#### 0.1) Build Symlinks to MASS-PRF and MASS-PRF Preprocess
Get the absolute path of your compiled `massprf` & `MASSPRF_preprocess`. These will be something like:

```plaintext
~/PATH/TO/MASSPRF/MASSPRF_10July2016/bin/massprf
~/PATH/TO/MASSPRF/MASSPRF_10July2016/MASSPRF_preprocessing_08July2016/MASSPRF_preprocess
```
Createa custom symlinks folder in your home directory (if you haven't already):
```
mkdir ~/symbolics
cd ~/symbolics
```
Create links to massprf and massprf_preprocess in that directory:
```
ln -s ~/PATH/TO/MASSPRF/MASSPRF_10July2016/bin/massprf massprf
ln -s ~/PATH/TO/MASSPRF/MASSPRF_10July2016/MASSPRF_preprocessing_08July2016/MASSPRF_preprocess massprf_preprocess
```
Add the symbolic link folder to your $PATH in ~/.bash_profile:
```
vim ~/.bash_profile
```
Scroll to the line that says something like:
```
PATH=$PATH:OTHERPATHS
```
Append a colon and add the path to your symbolic links, such that it looks like this:
```
PATH=$PATH:OTHERPATHS:$HOME/symbolics
Save and reload your profile:
source ~/.bash_profile
```
Save the file and reload your profile:
```
source ~/.bash_profile
```
#### 1) Install Conda Package Manager
   ```
   pip install conda
   ```
If you are working on a cluster, you may need to install Miniconda directly from the package due to file permission issues or if pip is unavailable.
In this case, download the file from https://conda.io/miniconda.html, run it on the cluster, and make sure the installation directory is added to your path via ~/.bashrc.
After installing Miniconda, source the path and verify Python version:
```
source ~/.bashrc
python
```
#### 2) Update Conda
Make sure your Conda is up to date by running:
```
conda update conda
```
#### 3) Update Python
Update Python to version 3.5, or optionally create a Python 3.5 virtual environment:
```
conda update python
```
#### 4) Install Package Dependencies
Add the Bioconda channel to Conda:
```
conda config --add channels bioconda
```
Install required packages:
```
conda install pandas
conda install biopython
conda install gffutils
conda install pyvcf
```
#### 5) Clone This Repository
```
git clone https://github.com/Townsend-Lab-Yale/massprf-pipeline.git
```
Usage
To run the example pipeline, see jobs.list in the massprf-pipeline folder. Example commands include:

Prepare Input Files

Ensure the input files are in the required format (e.g., FASTA or consensus FASTA). See the examples/ folder for sample files.
Run the Program

Example command:
```
./massprf -p examples/input_pol.fasta -d examples/input_div.fasta -o 1 -s 1 -exact 0
```

***Files and Folders***<br>
massprf-pipeline/: Scripts and documentation for genome-wide analysis.<br>
examples/: Sample input and output files.<br>
Source Code: All necessary .cpp and .h files for MASS-PRF.<br>
2D_Mapping_Scripts/: Scripts for generating 2D mapping visualizations.<br>
3D_Mapping_Scripts/: Scripts for generating 3D mapping visualizations.

***Citation***<br>
If you use MASS-PRF in your research, please cite:

Z.-M. Zhao, M. C. Campbell, N. Li, Z. Zhang, and J. P. Townsend. Detection of regional variation in selection intensity within protein-coding genes using DNA sequence polymorphism and divergence. Molecular Biology and Evolution, 2017. 34 (11), 3006-3022.<br> 
https://doi.org/10.1093/molbev/msx213

***Contact***<br>
For questions or support, contact:<br> 
Jeffrey Townsend<br>
Elihu Professor of Biostatistics and Ecology & Evolutionary Biology<br>
Email:Jeffrey.Townsend@Yale.edu
or <br>
Yide Jin <br>
Email: yide.jin@yale.edu/jinyide0202@gmail.com
 
