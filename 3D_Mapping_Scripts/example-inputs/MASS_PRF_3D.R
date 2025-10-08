library(readr)
library(dplyr)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")
install.packages("bio3d")

library(msa)
library(bio3d)

# Main wrapper
batchMASSPRF_Chimera <-
  function(designFile,
           doOnly = NULL,            # optionally subset rows of the design file
           hasHeader = FALSE,        # does the design file include a header?
           sigSetting = "average",   # "average", "any", "majority", "strict"
           onlySig = FALSE,           # if TRUE, color only significant sites
           rgb1 = c(250, 30, 30),     # color for positive γ
           rgb2 = c(30, 30, 250),     # color for negative γ
           bins = 510,               # number of bins
           ehColor = c(180, 180, 180), # default color
           midColor = c(240, 240, 240), # optional midpoint color
           logT = FALSE,             # if numeric, base for log-scaling; if FALSE, linear
           verbose = FALSE) {        # verbose printing
    
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
                                "formatList",
                                "outList")
    
    # unpack lists of inputs/outputs
    mypdbList                <- myDesignFile$pdbList
    myMASSPRF_Nuc_Fasta_List <- myDesignFile$MASSPRF_Nuc_Fasta_List
    myMASSPRF_Table_List     <- myDesignFile$MASSPRF_Table_List
    myFormatList             <- myDesignFile$formatList
    myOutList                <- myDesignFile$outList
    
    # verify equal lengths
    if (length(unique(lengths(list(
      mypdbList,
      myMASSPRF_Nuc_Fasta_List,
      myMASSPRF_Table_List,
      myFormatList,
      myOutList
    )))) != 1) {
      stop("Input Lists of different length")
    }
    numToDo <- length(mypdbList)
    if(verbose) print("Design File Looks Good!")
    
    # check format values
    for (f in myFormatList) {
      if (!f %in% c("NT", "AA")) {
        stop("Format column must contain only 'NT' or 'AA'")
      }
    }
    if(verbose) print("Formats are Valid")
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
      colnames(tbl) <- c(
        "Position","MS_PolSys","MA_PolSys","MS_PolRep","MA_PolRep",
        "MS_DivSys","MA_DivSys","MS_DivRep","MA_DivRep",
        "DivergentTime","Gamma","Lower_CI_Gamma","Upper_CI_Gamma",
        "PolymorphismMutationStatus","DivergenceMutationStatus"
      )
      tbl$Gamma          <- as.numeric(tbl$Gamma)
      tbl$Lower_CI_Gamma <- as.numeric(tbl$Lower_CI_Gamma)
      tbl$Upper_CI_Gamma <- as.numeric(tbl$Upper_CI_Gamma)
      allGammas <- c(allGammas, tbl$Gamma)
    }
    
    allGammas <- sort(unique(allGammas))
    G <- length(allGammas)
    
    gamma2Color <- data.frame(
      gamma = allGammas,
      logGamma = NA_real_,
      color    = "",
      logColor = "",
      stringsAsFactors = FALSE
    )
    
    if (!is.logical(logT)) {
      negIdx <- which(allGammas < 0)
      gamma2Color$logGamma <- log(abs(allGammas) + 1, base = logT)
      gamma2Color$logGamma[negIdx] <- -gamma2Color$logGamma[negIdx]
    } else {
      gamma2Color$logGamma <- allGammas
    }
    
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
    
    isSignificant <- function(LCI, UCI) {
      if (is.na(LCI) || is.na(UCI)) return(FALSE)
      if (LCI > 0 || UCI < 0) return(TRUE)
      FALSE
    }
    
    #
    # STEP 2: Loop over each gene
    #
    if(verbose) print("Begin Processing Individual Proteins")
    
    for (fileNum in seq_len(numToDo)) {
      thisFormat <- myFormatList[fileNum]
      if(verbose) cat("Processing", myMASSPRF_Table_List[fileNum], "with format", thisFormat, "\n")
      
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
      
      structUnaligned <- paste0(pdbseq(read.pdb(mypdbList[fileNum])), collapse = "")
      origUnaligned <- as.character(
        translate(readDNAStringSet(myMASSPRF_Nuc_Fasta_List[fileNum])[[1]])
      )
      origUnaligned <- gsub("\\*", "", origUnaligned)
      
      writeLines(c(">orig", origUnaligned, ">struct", structUnaligned),
                 "_tmpAln.fasta")
      aln <- msa(readAAStringSet("_tmpAln.fasta"), method = "Muscle")
      alnSeqs <- as.character(aln)
      origSeq   <- alnSeqs[1]
      structSeq <- alnSeqs[2]
      
      tmpMASS <- tbl[, c("Position","Gamma","Lower_CI_Gamma","Upper_CI_Gamma")]
      colnames(tmpMASS) <- c("Position","Gamma","LCI","UCI")
      
      if(verbose) print("Determining Coloring...")
      
      tmpMASS$color <- ehHex
      valid <- which(!is.na(tmpMASS$Gamma))
      for (it in valid) {
        idx <- which.min(abs(gamma2Color$gamma - tmpMASS$Gamma[it]))
        if (is.logical(logT)) {
          tmpMASS$color[it] <- gamma2Color$color[idx]
        } else {
          tmpMASS$color[it] <- gamma2Color$logColor[idx]
        }
      }
      
      Raw <- data.frame(
        AA   = strsplit(origUnaligned, "")[[1]],
        Gamma = NA_real_,
        LCI   = NA_real_,
        UCI   = NA_real_,
        color = ehHex,
        Significant = FALSE,
        stringsAsFactors = FALSE
      )
      
      massPRFCOPY <- tmpMASS
      pos <- 1
      while (nrow(massPRFCOPY) >= ifelse(thisFormat == "NT", 3, 1) &
             pos <= nrow(Raw)) {
        if(thisFormat == "NT"){
          group <- massPRFCOPY[1:3, ]
        } else if(thisFormat == "AA"){
          group <- massPRFCOPY[1, , drop = F]
        }
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
        idx <- which.min(abs(gamma2Color$gamma - Raw$Gamma[pos]))
        if (is.logical(logT)) {
          Raw$color[pos] <- gamma2Color$color[idx]
        } else {
          Raw$color[pos] <- gamma2Color$logColor[idx]
        }
        if(thisFormat == "NT"){
          massPRFCOPY <- massPRFCOPY[-(1:3), , drop = F]
        } else if(thisFormat == "AA"){
          massPRFCOPY <- massPRFCOPY[-1, , drop = F]
        }
        pos <- pos + 1
      }
      
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
      
      fStruct <- subset(Align, struct != "-")
      if (onlySig) {
        fStruct$color[!fStruct$Significant] <- ehHex
      }
      
      if(verbose) print("Begin Constructing Chimera Script")
      holdAllCommands <- c(paste0("color ", ehHex, ";"))
      
      uniqueCols <- unique(fStruct$color)
      for (k in seq_along(uniqueCols)) {
        holdAllCommands <- c(
          holdAllCommands,
          paste0("colordef Custom", k, " ", uniqueCols[k], ";")
        )
      }
      
      for (i in seq_len(nrow(fStruct))) {
        ix <- which(uniqueCols == fStruct$color[i])
        fStruct$colorName[i] <- paste0("Custom", ix)
      }
      
      lastName <- fStruct$colorName[1]
      startPos <- 1
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
      holdAllCommands <- c(
        holdAllCommands,
        paste0("color ", lastName,
               " #0:", startPos, "-", nrow(fStruct), ";")
      )
      
      writeLines(holdAllCommands, con = myOutList[fileNum])
    }
    
    if(verbose) print("Done!")
  }

# Example run
batchMASSPRF_Chimera(
  designFile = "design.tsv",
  hasHeader  = TRUE,
  onlySig    = FALSE,
  bins       = 15,
  midColor   = c(240, 240, 240),
  ehColor    = c(0, 180, 0),
  logT       = 2
)
