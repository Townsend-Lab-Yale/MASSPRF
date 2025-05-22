# massprf-protein-coloring
Maps the output of MASS-PRF into a file that can be used in Chimera to color proteins.

## Installation
Source the main script in R:
```r
source("batchMASSPRF_To_Chimera.R")
```

## Dependencies
- R packages:
  - `bio3d`
  - `Biostrings`
  - `msa`
- UCSF Chimera (for visualization)

## Example Calls
```r
# Single gene (S gene)
batchMASSPRF_Chimera(
  designFile = "example_inputs/S-design.tsv",
  hasHeader  = TRUE,
  onlySig    = FALSE,
  bins       = 10,
  midColor   = c(240,245,240),
  logT       = FALSE
)

# Batch genes
batchMASSPRF_Chimera(
  designFile = "example_inputs/batch-design.tsv",
  hasHeader  = TRUE,
  onlySig    = FALSE,
  bins       = 15,
  midColor   = c(240,245,240),
  logT       = 2
)
```

## Inputs
Your `designFile` is a TSV (tab-separated) file with columns:

1. **pdbFile**: Path to PDB structure file.  
2. **MASSPRF_Nuc_Fasta**: Nucleotide FASTA from MASS-PRF (in-frame).  
3. **MASSPRF_Table**: MASS-PRF output CSV with γ and CI values.  
4. **scaling**: Scaling factor (1 or multiple of 3).  
5. **outfile**: Path to write Chimera coloring script.

Additional parameters:
- **doOnly**: Vector of row indices to process. Default: all.
- **hasHeader**: `TRUE`/`FALSE` (does design file include header?). Default: `FALSE`.
- **sigSetting**: Significance rule (`average`, `any`, `majority`, `strict`). Default: `average`.
- **onlySig**: `TRUE` to color only significant sites. Default: `TRUE`.
- **rgb1**, **rgb2**: RGB triples for positive/negative γ. Default: red (`c(250,30,30)`), blue (`c(30,30,250)`).
- **midColor**: Optional RGB triple for midpoint gradient. Default: `NULL`.
- **bins**: Number of color bins. Default: `10` in examples.
- **ehColor**: RGB triple for gaps/non-significant. Default: grey (`c(128,128,128)`).

## Usage in Chimera
After running the script, load the generated file in Chimera’s command line:
```
read <outfile>;
```
