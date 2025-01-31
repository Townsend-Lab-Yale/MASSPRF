# MASS-PRF: Model Averaged Site Selection via Poisson Random Field 
Updated in 2025/1/30

## Overview
MASS-PRF is a computational tool designed to detect regional variation in selection intensity within protein-coding genes using DNA sequence polymorphism and divergence data. This repository includes the program, preprocessing scripts, and a pipeline for genome-wide analysis.

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
 
