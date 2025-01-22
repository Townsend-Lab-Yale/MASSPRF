# massprf-protein-coloring
maps the output of massprf into a file that can be used in Chimera to color proteins
execute 

    source("batchMASSPRF_To_Chimera.R") 

from R, which creates the batchMASSPRF_Chimera command in your local enviornment.

## dependencies
  require("bio3d")
  
  require("Biostrings")
  
  require("muscle")
  
The coloring happens in Chimera, so you'll need that program too.


## example calls
    #do one gene (S gene)
    batchMASSPRF_Chimera(designFile = "example_inputs/S-design.tsv",hasHeader = T,onlySig = F,bins=10,midColor = c(240,245,240),logT = F)

    #do multiple genes jointly
    batchMASSPRF_Chimera(designFile = "example_inputs/batch-design.tsv",hasHeader = T,onlySig = F,bins=15,midColor = c(240,245,240),logT = 2)


## Inputs
designFile: a tsv file with each of the following columns (column names ignored)

    pdbFile: The protein structure you are trying to map onto, in pdb format. REQUIRED INPUT


    MASSPRF_Nuc_Fasta: the nucleotide sequence representing your MASSPRF input, in open reading frame and in fasta format. REQUIRED INPUT


    MASSPRF_Table: The output table from MASSPRF containing the gamma values and CI. REQUIRED INPUT
    

    scaling: the scale factor from MASSPRF--look in the output of MASSPRF to find this. 1 or a multiple of 3. REQUIRED INPUT


    outfile: where the resultant commands for Chimera to color will be stored. Defaults chimeraColoring.txt


    
    

doOnly: a vector of rows from the design file to use. Default is use all. (ex. c(1,3,4))
    

hasHeader: logical. Does your design file contain a header? Defaults FALSE
    

sigSetting: what are you using to determine significance? only affects scale factor 1; otherwise sitewise significance is used. Defaults average


onlySig: T or F do you want to only color significant sites. Defaults T


rgb1 and rgb2: vectors of size 3, corresponding to a rgb value for coloring the selection intensity. Defaults red and blue

    
midColor: vector of size 3, optional, corresponding to an rgb value. A midpoint color for a gradient between rgb1 and 2. Defaults NULL (i.e. 2 color gradient)


bins: how many equally spaced apart color categories you want to use for your data. Defaults 10


ehColor: what do you want to color 'eh' residues (missing data, nonsignificant) as a vector size three of rgb. Defaults grey


## To Use
run "read \<outfile\>;" on the Chimera command line 

