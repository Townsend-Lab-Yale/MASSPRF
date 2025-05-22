
# 1. Read the original TSV (assuming it is called test.tsv)
library(readr)
library(dplyr)
design <- read_tsv("temp.tsv", col_types = cols())

# 2. Change the MASSPRF_Table column to the full path of the corresponding file in the working directory
# Here we assume that all the files are in the C:/Users/bjyid/Documents directory
wd <- normalizePath("C:/Users/bjyid/Documents", winslash = "/")
design <- design %>%
mutate(
MASSPRF_Table = file.path(wd, basename(MASSPRF_Table), fsep = "/")
)

# 3. Write it back (you can write it to a new file to avoid overwriting the original file)
write_tsv(design, "test_fixed.tsv")

library(msa)

# Main wrapper: reads a design file, parses inputs,
# computes a global color lookup for all γ values,
# then processes each protein to emit a Chimera coloring script.
batchMASSPRF_Chimera <-
  function(designFile,
           doOnly = NULL,            # optionally subset rows of the design file
           hasHeader = FALSE,        # does the design file include a header?
           sigSetting = "average",   # how to aggregate significance calls: "average", "any", "majority", "strict"
           onlySig = TRUE,           # if TRUE, color only significant sites (others get ehColor)
           rgb1 = c(250, 30, 30),     # color for positive γ
           rgb2 = c(30, 30, 250),     # color for negative γ
           bins = 510,               # number of bins in the color ramp
           ehColor = c(128, 128, 128), # “eh, default” color for gaps / non-sig
           midColor = NULL,          # optional midpoint color
           logT = FALSE,             # if numeric, base for log-scaling; if FALSE, linear
           verbose = FALSE) {        # verbose progress printing
  
    # load dependencies
    require("bio3d")
    require("Biostrings")
    require("msa")
  
    # sanity-check the significance‐aggregation parameter
    if (!sigSetting %in% c("average", "any", "majority", "strict")) {
      stop("sigSetting must be average, any, majority, strict")
    }
    if(verbose) print("Processing Design File")
  
    # read and optionally subset the design file
    myDesignFile <- read.csv(designFile, sep = "\t", header = hasHeader)
    if (!is.null(doOnly)) myDesignFile <- myDesignFile[doOnly, ]
    colnames(myDesignFile) <- c("pdbList",
                                "MASSPRF_Nuc_Fasta_List",
                                "MASSPRF_Table_List",
                                "scalingList",
                                "outList")
  
    # unpack lists of inputs/outputs
    mypdbList                <- myDesignFile$pdbList
    myMASSPRF_Nuc_Fasta_List <- myDesignFile$MASSPRF_Nuc_Fasta_List
    myMASSPRF_Table_List     <- myDesignFile$MASSPRF_Table_List
    myScalingList            <- as.numeric(myDesignFile$scalingList)
    myOutList                <- myDesignFile$outList
  
    # verify equal lengths
    if (length(unique(lengths(list(
      mypdbList,
      myMASSPRF_Nuc_Fasta_List,
      myMASSPRF_Table_List,
      myScalingList,
      myOutList
    )))) != 1) {
      stop("Input Lists of different length")
    }
    numToDo <- length(mypdbList)
    if(verbose) print("Design File Looks Good!")
  
    # check scaling values: must be 1 or multiple of 3
    for (i in myScalingList) {
      if (i == 1) next()
      if (i %% 3 == 0) next()
      stop("Scaling must be 1 or a multiple of 3 (check MASSPRF file for details)")
    }
    if(verbose) print("Scalings are Valid")
    if(verbose) print("Discerning RGB Scaling for all Samples")
  
    #
    # STEP 1: Build a global γ → color lookup table
    #
  
    # convert input rgb triples into hex
    redHex  <- rgb(rgb1[1], rgb1[2], rgb1[3], maxColorValue = 255)
    blueHex <- rgb(rgb2[1], rgb2[2], rgb2[3], maxColorValue = 255)
    midHex  <- if (!is.null(midColor)) {
      rgb(midColor[1], midColor[2], midColor[3], maxColorValue = 255)
    } else NULL
  
    numberOfColors <- bins
    ehHex <- rgb(ehColor[1], ehColor[2], ehColor[3], maxColorValue = 255)
  
    # gather all γ values across every table
    allGammas <- NULL
    for (i in seq_len(numToDo)) {
      tbl <- read.csv(myMASSPRF_Table_List[i],
                      sep = ",", header = TRUE,
                      stringsAsFactors = FALSE)
      # assign the expected 15 columns
      colnames(tbl) <- c(
        "Position","MS_PolSys","MA_PolSys","MS_PolRep","MA_PolRep",
        "MS_DivSys","MA_DivSys","MS_DivRep","MA_DivRep",
        "DivergentTime","Gamma","Lower_CI_Gamma","Upper_CI_Gamma",
        "PolymorphismMutationStatus","DivergenceMutationStatus"
      )
      # cast γ and CI to numeric
      tbl$Gamma          <- as.numeric(tbl$Gamma)
      tbl$Lower_CI_Gamma <- as.numeric(tbl$Lower_CI_Gamma)
      tbl$Upper_CI_Gamma <- as.numeric(tbl$Upper_CI_Gamma)
  
      # collect for global palette
      gamVec <- tbl$Gamma
      allGammas <- if (is.null(allGammas)) gamVec else c(allGammas, gamVec)
    }
  
    # unique sorted γ
    allGammas <- sort(unique(allGammas))
    G <- length(allGammas)
  
    # prepare data.frame to hold linear or log γ and hex colors
    gamma2Color <- data.frame(
      gamma = allGammas,
      logGamma = NA_real_,
      color    = "",
      logColor = "",
      stringsAsFactors = FALSE
    )
  
    # compute log-Gamma if requested
    if (!is.logical(logT)) {
      negIdx <- which(allGammas < 0)
      gamma2Color$logGamma <- log(abs(allGammas) + 1, base = logT)
      gamma2Color$logGamma[negIdx] <- -gamma2Color$logGamma[negIdx]
    } else {
      gamma2Color$logGamma <- allGammas
    }
  
    # build linear palette → assign gamma2Color$color
    cuts <- cut(
      gamma2Color$gamma,
      breaks = seq(min(allGammas), max(allGammas), length.out = numberOfColors),
      include.lowest = TRUE
    )
    pal1 <- if (is.null(midHex)) {
      colorRampPalette(c(blueHex, redHex))(numberOfColors - 1)
    } else {
      colorRampPalette(c(blueHex, midHex, redHex))(numberOfColors - 1)
    }
    gamma2Color$color <- pal1[as.integer(cuts)]
  
    # build log palette → assign gamma2Color$logColor
    cutsLog <- cut(
      gamma2Color$logGamma,
      breaks = seq(min(gamma2Color$logGamma),
                   max(gamma2Color$logGamma),
                   length.out = numberOfColors),
      include.lowest = TRUE
    )
    pal2 <- if (is.null(midHex)) {
      colorRampPalette(c(blueHex, redHex))(numberOfColors - 1)
    } else {
      colorRampPalette(c(blueHex, midHex, redHex))(numberOfColors - 1)
    }
    gamma2Color$logColor <- pal2[as.integer(cutsLog)]
  
    if (verbose) {
      plot(seq_len(G), gamma2Color$logGamma, col = gamma2Color$logColor,
           main = "LogGamma Color Ramp")
    }
    if(verbose) print("RGB Bins for all Data Defined!")
  
    # helper: test if a CI pair is significant (pulls out zeros/NA)
    isSignificant <- function(LCI, UCI) {
      if (is.na(LCI) || is.na(UCI)) return(FALSE)
      if (LCI > 0 || UCI < 0) return(TRUE)
      FALSE
    }
  
    #
    # STEP 2: Loop over each gene, align sequences, map colors, write Chimera script
    #
    if(verbose) print("Begin Processing Individual Proteins")
  
    for (fileNum in seq_len(numToDo)) {
      if(verbose) cat("Processing", myMASSPRF_Table_List[fileNum], "\n")
  
      # re-read table, cast numeric
      tbl <- read.csv(myMASSPRF_Table_List[fileNum],
                      sep = ",", header = TRUE,
                      stringsAsFactors = FALSE)
      colnames(tbl) <- c(
        "Position","MS_PolSys","MA_PolSys","MS_PolRep","MA_PolRep",
        "MS_DivSys","MA_DivSys","MS_DivRep","MA_DivRep",
        "DivergentTime","Gamma","Lower_CI_Gamma","Upper_CI_Gamma",
        "PolymorphismMutationStatus","DivergenceMutationStatus"
      )
      tbl$Gamma          <- as.numeric(tbl$Gamma)
      tbl$Lower_CI_Gamma <- as.numeric(tbl$Lower_CI_Gamma)
      tbl$Upper_CI_Gamma <- as.numeric(tbl$Upper_CI_Gamma)
  
      # get PDB‐derived sequence
      structUnaligned <- paste0(pdbseq(read.pdb(mypdbList[fileNum])), collapse = "")
      # get original translated sequence
      origUnaligned <- as.character(
        translate(readDNAStringSet(myMASSPRF_Nuc_Fasta_List[fileNum])[[1]])
      )
      origUnaligned <- gsub("\\*", "", origUnaligned)
  
      # perform pairwise MSA for orig vs struct
      writeLines(c(">orig", origUnaligned, ">struct", structUnaligned),
                 "_tmpAln.fasta")
      aln <- msa(readAAStringSet("_tmpAln.fasta"), method = "Muscle")
      alnSeqs <- as.character(aln)
      origSeq   <- alnSeqs[1]
      structSeq <- alnSeqs[2]
  
      # prepare tmpMASS with Position/Gamma/LCI/UCI
      tmpMASS <- tbl[, c("Position","Gamma","Lower_CI_Gamma","Upper_CI_Gamma")]
      colnames(tmpMASS) <- c("Position","Gamma","LCI","UCI")
  
      if(verbose) print("Determining Coloring...")
  
      # initialize default color (ehHex) for every row
      tmpMASS$color <- ehHex
      # only recolor where Gamma is non-NA
      valid <- which(!is.na(tmpMASS$Gamma))
      for (it in valid) {
        idx <- which.min(abs(gamma2Color$gamma - tmpMASS$Gamma[it]))
        if (is.logical(logT)) {
          tmpMASS$color[it] <- gamma2Color$color[idx]
        } else {
          tmpMASS$color[it] <- gamma2Color$logColor[idx]
        }
      }
  
      # build the per‐residue Raw data.frame
      Raw <- data.frame(
        AA   = strsplit(origUnaligned, "")[[1]],
        Gamma = NA_real_,
        LCI   = NA_real_,
        UCI   = NA_real_,
        color = ehHex,
        Significant = FALSE,
        stringsAsFactors = FALSE
      )
  
      # fill in values for scaling = 1 vs >1
      if (myScalingList[fileNum] == 1) {
        # group every 3 rows → one AA
        massPRFCOPY <- tmpMASS
        pos <- 1
        while (nrow(massPRFCOPY) >= 3) {
          group <- massPRFCOPY[1:3, ]
          group <- group[complete.cases(group), ]
          Raw$Gamma[pos] <- mean(group$Gamma)
          Raw$LCI[pos]   <- mean(group$LCI)
          Raw$UCI[pos]   <- mean(group$UCI)
          Raw$Significant[pos] <- {
            sigs <- mapply(isSignificant, LCI = group$LCI, UCI = group$UCI)
            switch(sigSetting,
                   any     = any(sigs),
                   majority= sum(sigs) >= 2,
                   strict  = all(sigs),
                   average = isSignificant(mean(group$LCI), mean(group$UCI)))
          }
          # color by nearest γ
          idx <- which.min(abs(gamma2Color$gamma - Raw$Gamma[pos]))
          if (is.logical(logT)) {
            Raw$color[pos] <- gamma2Color$color[idx]
          } else {
            Raw$color[pos] <- gamma2Color$logColor[idx]
          }
          # drop processed
          massPRFCOPY <- massPRFCOPY[-(1:3), ]
          pos <- pos + 1
        }
      } else {
        # scaled >1: each row in tmpMASS expands to aaScale residues
        aaScale <- myScalingList[fileNum] / 3
        pos <- 1
        for (i in seq_len(nrow(tmpMASS))) {
          block <- pos:(pos + aaScale - 1)
          Raw$Gamma[block] <- tmpMASS$Gamma[i]
          Raw$LCI[block]   <- tmpMASS$LCI[i]
          Raw$UCI[block]   <- tmpMASS$UCI[i]
          Raw$color[block] <- tmpMASS$color[i]
          Raw$Significant[block] <- isSignificant(tmpMASS$LCI[i],
                                                  tmpMASS$UCI[i])
          pos <- pos + aaScale
        }
      }
  
      # align Raw → struct positions
      Align <- data.frame(
        orig = strsplit(origSeq, "")[[1]],
        struct= strsplit(structSeq, "")[[1]],
        gamma = NA_real_,
        color = ehHex,
        Significant = FALSE,
        stringsAsFactors = FALSE
      )
      aIndex <- 1
      rIndex <- 1
      while (rIndex <= nrow(Raw) && aIndex <= nrow(Align)) {
        if (Align$orig[aIndex] == Raw$AA[rIndex]) {
          Align$gamma[aIndex] <- Raw$Gamma[rIndex]
          Align$color[aIndex] <- Raw$color[rIndex]
          Align$Significant[aIndex] <- Raw$Significant[rIndex]
          rIndex <- rIndex + 1
        } else if (Align$orig[aIndex] == "-") {
          Align$color[aIndex] <- ehHex
        }
        aIndex <- aIndex + 1
      }
  
      # filter down to only struct residues
      fStruct <- subset(Align, struct != "-")
      # optionally mask non-significant
      if (onlySig) {
        fStruct$color[!fStruct$Significant] <- ehHex
      }
  
      # now build the Chimera command script
      if(verbose) print("Begin Constructing Chimera Script")
      holdAllCommands <- c(paste0("color ", ehHex, ";"))
  
      uniqueCols <- unique(fStruct$color)
      # define custom colors
      for (k in seq_along(uniqueCols)) {
        holdAllCommands <- c(
          holdAllCommands,
          paste0("colordef Custom", k, " ", uniqueCols[k], ";")
        )
      }
  
      # assign ranges
      k <- 1
      startPos <- 1
      curCol   <- fStruct$colorName <- NA_character_
      # tag each residue with its CustomN
      for (i in seq_len(nrow(fStruct))) {
        ix <- which(uniqueCols == fStruct$color[i])
        fStruct$colorName[i] <- paste0("Custom", ix)
      }
  
      # collapse contiguous runs
      lastName <- fStruct$colorName[1]
      for (i in 2:nrow(fStruct)) {
        if (fStruct$colorName[i] != lastName) {
          holdAllCommands <- c(
            holdAllCommands,
            paste0("color ", lastName,
                   " #0:", startPos, "-", i - 1, ";")
          )
          startPos <- i
          lastName <- fStruct$colorName[i]
        }
      }
      # final segment
      holdAllCommands <- c(
        holdAllCommands,
        paste0("color ", lastName,
               " #0:", startPos, "-", nrow(fStruct), ";")
      )
  
      # write to output file
      writeLines(holdAllCommands, con = myOutList[fileNum])
    }
  
    if(verbose) print("Done!")
  }

# example invocation
setwd("C:/Users/bjyid/Documents")
batchMASSPRF_Chimera(
  designFile = "test_fixed.tsv",
  hasHeader  = TRUE,
  onlySig    = FALSE,
  bins       = 510,
  verbose    = TRUE
)
#example call
#batchMASSPRF_Chimera(designFile = "~/../Desktop/5-design.tsv",
#                     hasHeader = TRUE,
#                     onlySig   = FALSE,
#                     bins      = 10,
#                     midColor  = c(240,245,240),
#                     logT      = 2)
#
# designFile is the path to a TSV where the columns are:
#   1) path to the PDB file
#   2) path to the nucleotide FASTA from MASS-PRF
#   3) path to the table of results from MASS-PRF
#   4) the scaling factor from MASS-PRF
#   5) the output file name (Chimera script)
# hasHeader: does your design file include a header line? Default FALSE.
# doOnly: optionally specify which rows of your design file to process.
# logT: FALSE (linear) or a numeric base for log-transforming γ values.
#       Negative γ are handled by log(abs(γ)+1) and restoring sign.
#       (Hint: after Chimera-script generation, in Chimera CLI: read <outfile>;)
# sigSetting: how to determine significance when scaling == 1 only.
# onlySig: TRUE/FALSE to color only significant sites (others get ehColor).
# rgb1 / rgb2: vectors of 3 RGB values—positive vs negative selection colors.
# bins: number of equally spaced color categories (max ~510 per colorRamp).
# ehColor: RGB triple for “eh” (missing data or nonsignificant) residues.
