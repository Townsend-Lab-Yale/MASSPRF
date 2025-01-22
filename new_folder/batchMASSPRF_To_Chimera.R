#example call
#batchMASSPRF_Chimera(designFile = "~/../Desktop/5-design.tsv",hasHeader = T,onlySig = F,bins=10,midColor = c(240,245,240),logT = 2)


#designFile is the path to a tsv where the columns are the path to the pdb, path to the nucleotide fasta from massprf,
## path to the table of results from massprf, the scaling factor from massprf, and the output file name. 
#hasHeader just specifies if your design file has a header sequence or not. Default F.
#doOnly allows users to specify what rows of their overall design file they want to do the analysis on
#logT defaults false, otherwise takes in the log to transform the data by for coloring (2,10, etc).
## for negative values, you take the log of the absolute value and then return the sign. 
### (hint: execute the following on the Chimera command line): read <outfile>;
#sigSetting: what are you using to determine significance? only affects scale factor 1; otherwise sitewise significance is used
#onlySig: T or F do you want to only color significant sites
#rgb1 and rgb2: vectors of size 3, corresponding to a rgb value for coloring the selection intensity
#bins: how many equally spaced apart color categories you want to use for your data. Defaults 510, rbg precision limit for 2 cols
#ehColor: what do you want to color 'eh' residues (missing data, nonsignificant) as a vector size three of rgb


batchMASSPRF_Chimera <-
  function(designFile,
           doOnly=NULL,
           hasHeader=F,
           sigSetting = "average",
           onlySig = T,
           rgb1 = c(250, 30, 30),
           rgb2 = c(30, 30, 250),
           bins = 510,
           ehColor = c(128, 128, 128),
           midColor = NULL,
           logT = F,
           verbose=F) {
    require("bio3d")
    require("Biostrings")
    require("muscle")
    if (!sigSetting %in% c("average", "any", "majority", "strict")) {
      stop("sigSetting must be average, any, majority, strict")
    }
    if(verbose){print("Processing Design File")}
    myDesignFile<-read.csv(file=designFile,sep="\t",header=hasHeader)
    if(!is.null(doOnly)){myDesignFile<-myDesignFile[doOnly,]}
    colnames(myDesignFile)<-c("pdbList","MASSPRF_Nuc_Fasta_List","MASSPRF_Table_List","scalingList","outList")
    #verify lengths of lists all the same
    mypdbList <- myDesignFile$pdbList
    myMASSPRF_Nuc_Fasta_List <- myDesignFile$MASSPRF_Nuc_Fasta_List
    myMASSPRF_Table_List <- myDesignFile$MASSPRF_Table_List
    myScalingList <- as.numeric(myDesignFile$scalingList)
    myOutList <- myDesignFile$outList
    if (length(unique(lengths(
      list(
        mypdbList,
        myMASSPRF_Nuc_Fasta_List,
        myMASSPRF_Table_List,
        myScalingList,
        myOutList
      )
    ))) != 1) {
      stop("Input Lists of different length")
    }
    numToDo <- length(mypdbList)
    if(verbose){print("Design File Looks Good!")}
    for (i in myScalingList) {
      if (i == 1) {
        next()
      }
      if (i  %%  3  ==  0) {
        next()
      }
      stop("Scaling must be 1 or a multiple of 3 (check MASSPRF file for details)")
    }
    if(verbose){print("Scalings are Valid")}
    if(verbose){print("Discerning RGB Scaling for all Samples")}
    #lets figure out the best scale factor for everything
    redHex <- rgb(rgb1[1], rgb1[2], rgb1[3], maxColorValue = 255)
    blueHex <- rgb(rgb2[1], rgb2[2], rgb2[3], maxColorValue = 255)
    midHex <- NULL
    if (!is.null(midColor)) {
      midHex <-
        rgb(midColor[1], midColor[2], midColor[3], maxColorValue = 255)
    }
    numberOfColors <- bins
    ehHex <- rgb(ehColor[1], ehColor[2], ehColor[1], maxColorValue = 255)
    allGammas <- NULL
    for (i in 1:numToDo) {
      MASSPRF_FILE <- read.csv(myMASSPRF_Table_List[i], sep = "\t")
      tmpMASS <-
        MASSPRF_FILE[, c("Position", "Gamma", "Lower_CI_Gamma", "Upper_CI_Gamma")]
      colnames(tmpMASS) <- c("Position", "Gamma", "LCI", "UCI")
      if (is.null(allGammas)) {
        allGammas <- tmpMASS$Gamma
      } else{
        allGammas <- c(allGammas, tmpMASS$Gamma)
      }
    }
    allGammas <- unique(allGammas)
    allGammas <- allGammas[order(allGammas)]
    
    gamma2Color <- data.frame(matrix(ncol = 4, nrow = length(allGammas)))
    colnames(gamma2Color) <- c("gamma", "logGamma", "color", "logColor")
    gamma2Color$gamma <- allGammas
    if (!is.logical(logT)) {
      
      isNegative<-(which(allGammas<0))
      gamma2Color$logGamma <-
        log((abs(allGammas)+1), base = logT)
      gamma2Color$logGamma[isNegative]<-(gamma2Color$logGamma[isNegative]*(-1))
      
      
    } else{
      gamma2Color$logGamma <- allGammas
    }
    
    #get color, as normal
    values <- gamma2Color$gamma
    ii <-
      cut(
        values,
        breaks = seq(min(values), max(values), len = numberOfColors),
        include.lowest = TRUE
      )
    ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
    if (is.null(midHex)) {
      myColors <-
        colorRampPalette(c(blueHex, redHex))((numberOfColors - 1))[ii]
    }
    if (!is.null(midHex)) {
      myColors <-
        colorRampPalette(c(blueHex, midHex, redHex))((numberOfColors - 1))[ii]
    }
    gamma2Color$color <- myColors
    
    #get color for log values
    values <- gamma2Color$logGamma
    ii <-
      cut(
        values,
        breaks = seq(min(values), max(values), len = numberOfColors),
        include.lowest = TRUE
      )
    ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
    if (is.null(midHex)) {
      myColors <-
        colorRampPalette(c(blueHex, redHex))((numberOfColors - 1))[ii]
    }
    if (!is.null(midHex)) {
      myColors <-
        colorRampPalette(c(blueHex, midHex, redHex))((numberOfColors - 1))[ii]
    }
    gamma2Color$logColor <- myColors
    

    if(verbose){plot(1:length(gamma2Color$logColor),gamma2Color$logGamma,col=gamma2Color$logColor)}

    
    if(verbose){print("RGB Bins for all Data Defined!")}
    
    isSignificant <- function(LCI, UCI) {
      if (is.na(LCI)) {
        return(FALSE)
      }
      if (is.na(UCI)) {
        return(FALSE)
      }
      if (LCI > 0) {
        return(TRUE)
      }
      if (UCI < 0) {
        return(TRUE)
      }
      return(FALSE)
    }
    if(verbose){print("Begin Processing Individual Proteins")}
    for (fileNum in 1:numToDo) {
      print(paste0("Processing ",myMASSPRF_Table_List[fileNum]))
      MASSPRF_FILE <- read.csv(myMASSPRF_Table_List[fileNum], sep = "\t")
      structUnaligned <-
        paste0(pdbseq(read.pdb(mypdbList[fileNum])), collapse = "")
      origUnaligned <-
        gsub("\\*", "", as.character(Biostrings::translate(
          Biostrings::readDNAStringSet(myMASSPRF_Nuc_Fasta_List[fileNum])
        )[[1]]))
      tmpAln <- c(">orig", origUnaligned, ">struct", structUnaligned)
      writeLines(tmpAln, con = "_tmpAln.fasta")
      if(verbose){print("Aligning...")}
      res <- muscle(Biostrings::readAAStringSet("_tmpAln.fasta"))
      #file.remove("_tmpAln.fasta")
      res <- as.character(res)
      origSeq <- res[1]
      structSeq <- res[2]
      scaling <- myScalingList[fileNum]
      onlySig <- onlySig
      sigSetting <- sigSetting
      tmpMASS <-
        MASSPRF_FILE[, c("Position", "Gamma", "Lower_CI_Gamma", "Upper_CI_Gamma")]
      colnames(tmpMASS) <- c("Position", "Gamma", "LCI", "UCI")
      if(verbose){print("Determining Coloring...")}
      for (it in 1:nrow(tmpMASS)) {
        if (is.logical(logT)) {
          tmpMASS[it, "color"] <-
            gamma2Color[which(gamma2Color$gamma == tmpMASS$Gamma[it]), "color"]
        } else{
          tmpMASS[it, "color"] <-
            gamma2Color[which(gamma2Color$gamma == tmpMASS$Gamma[it]), "logColor"]
        }
      }

      Raw <- data.frame(matrix(nrow = nchar(origUnaligned), ncol = 6))
      colnames(Raw) <- c("AA", "Gamma", "LCI", "UCI", "color", "Significant")
      Raw$AA <- unlist(strsplit(x = origUnaligned, split = ""))
      i <- 1
      if (scaling == 1) {
        meanAAPosGamma <- data.frame(matrix(nrow = nchar(origUnaligned), ncol=6))
        colnames(meanAAPosGamma) <-
          c("Position", "Gamma", "LCI", "UCI", "color", "Significant")
        meanAAPosGamma$Position <- ((1:nchar(origUnaligned)) - 1)
        massPRFCOPY <- tmpMASS[, c("Gamma", "LCI", "UCI", "color")]
        for (j in 1:nrow(meanAAPosGamma)) {
          theseThree <- massPRFCOPY[1:3, ]
          theseThree <- theseThree[complete.cases(theseThree), ]
          meanGam <- mean(theseThree$Gamma)
          meanLCI <- mean(theseThree$LCI)
          meanUCI <- mean(theseThree$UCI)
          meanSig <- NA
          sigs <-
            mapply(theseThree$LCI, FUN = isSignificant, UCI = theseThree$UCI)
          if (sigSetting == "any") {
            if (sum(sigs) >= 1) {
              meanSig <- TRUE
            } else{
              meanSig <- FALSE
            }
          }
          if (sigSetting == "majority") {
            if (sum(sigs) >= 2) {
              meanSig <- TRUE
            } else{
              meanSig <- FALSE
            }
          }
          if (sigSetting == "strict") {
            if (FALSE %in% sigs) {
              meanSig <- FALSE
            } else{
              meanSig <- TRUE
            }
          }
          if (sigSetting == "average") {
            meanSig <- isSignificant(meanLCI, meanUCI)
          }
          meanAAPosGamma[j, ] <-
            c(meanAAPosGamma$Position[j],
              meanGam,
              meanLCI,
              meanUCI,
              NA,
              meanSig)
          massPRFCOPY <- massPRFCOPY[-c(1, 2, 3), ]
        }
        Raw$Gamma <- meanAAPosGamma$Gamma
        Raw$LCI <- meanAAPosGamma$LCI
        Raw$UCI <- meanAAPosGamma$UCI
        Raw$Significant <- meanAAPosGamma$Significant
        for (tc in 1:nrow(Raw)) {
          if (is.logical(logT)) {
            Raw[tc, "color"] <-
              gamma2Color$color[which.min(abs(gamma2Color$gamma-Raw$Gamma[tc]))]
          } else {
            Raw[tc, "color"] <-
              gamma2Color$logColor[which.min(abs(gamma2Color$gamma-Raw$Gamma[tc]))]
          }
        }
      } else{
        if (scaling %% 3 != 0) {
          stop("Non 1 scaling needs to be divisible by 3!!!!")
        }
        aaScale <- scaling / 3
        while (i < nrow(tmpMASS)) {
          thisG <- tmpMASS$Gamma[i]
          thisC <- tmpMASS$color[i]
          thisL <- tmpMASS$LCI[i]
          thisU <- tmpMASS$UCI[i]
          thisSig <- isSignificant(thisL, thisU)
          
          relevantIndecies <- ((((i * aaScale) - (aaScale - 1))):(i * aaScale))
          if(max(relevantIndecies)>length(Raw$Gamma)){i<-i+1;next()} #if sizes don't match because terminal gap!
          Raw$Gamma[relevantIndecies] <- thisG
          Raw$color[relevantIndecies] <- thisC
          Raw$LCI[relevantIndecies] <- thisL
          Raw$UCI[relevantIndecies] <- thisU
          Raw$Significant[relevantIndecies] <- thisSig
          i <- i + 1
        }
      }
      Align <- data.frame(matrix(ncol = 5, nrow = nchar(structSeq)))
      colnames(Align) <- c("orig", "struct", "gamma", "color", "Significant")
      
      Align$orig <- unlist(strsplit(origSeq, ""))
      Align$struct <- unlist(strsplit(structSeq, ""))
      
      
      aIndex <- 1
      rIndex <- 1
      while (rIndex < nrow(Raw)) {
        if (Align$orig[aIndex] == Raw$AA[rIndex]) {
          Align$gamma[aIndex] <- Raw$Gamma[rIndex]
          Align$color[aIndex] <- Raw$color[rIndex]
          Align$Significant[aIndex] <- Raw$Significant[rIndex]
          rIndex <- rIndex + 1
        } else if (Align$orig[aIndex] == "-") {
          Align$gamma[aIndex] <- NA
          Align$color[aIndex] <- ehHex #grey for missing
        }
        aIndex <- aIndex + 1
      }
      
      fStruct <- Align[(Align$struct != "-"), ]
      fStruct[is.na(fStruct$gamma), "color"] <- ehHex
      fStruct[is.na(fStruct$Significant), "Significant"] <- FALSE
      if (onlySig) {
        fStruct[which(fStruct$Significant == FALSE), "color"] <- ehHex
      }
      
      
      setmyColor <- function(HEX, COLORNAME, POS) {
        myString <- "color COLORNAME #0:POS;"
        myString <- gsub("HEX", HEX, myString)
        myString <- gsub("COLORNAME", COLORNAME, myString)
        myString <- gsub("POS", POS, myString)
        return(myString)
      }
      if(verbose){print("Begin Constructing Chimera Script")}
      if(verbose){print("Defining all Colors")}
      holdAllCommands <- ""
      uniqueColors<-unique(fStruct$color)
      for(i in 1:length(uniqueColors)){
        holdAllCommands<-c(holdAllCommands,paste0("colordef Custom",i," ",uniqueColors[i],";"))
      }
      
      fStruct$colorName<-""
      for (i in 1:nrow(fStruct)) {
        myHex <- fStruct$color[i]
        myHexPos<-which(grepl(pattern = myHex,uniqueColors))
        #myPOS<-(i-1) #convert to python indexing
        myPOS <- (i) #convert to python indexing
        fStruct$colorName[i]<-paste0("Custom",(myHexPos))
        #tmp <- setmyColor(myHex, paste0("Custom", (myHexPos)), myPOS)
        #holdAllCommands <- c(holdAllCommands, tmp)
      }
      if(verbose){print("Assigning colors and Optimizing Chimera script...")}
      mcInd<-1
      comTemp<-"color COLORNAME #0:POS1-POS2;"
      lastCol<-fStruct$colorName[mcInd]
      spos<-1
      cpos<-0
      while(mcInd<=nrow(fStruct)){
        currCol<-fStruct$colorName[mcInd]
        if(currCol!=lastCol){
          ctmp<-gsub("COLORNAME",lastCol,comTemp)
          ctmp<-gsub("POS1",as.character(spos),ctmp)
          ctmp<-gsub("POS2",as.character(mcInd-1),ctmp)
          #ctmp<-gsub("POS1",as.character(spos-1),ctmp) #convert to python indexing
          #ctmp<-gsub("POS2",as.character(mcInd-2),ctmp)
          holdAllCommands<-c(holdAllCommands,ctmp)
          spos<-mcInd
          lastCol<-currCol
        }
        mcInd<-mcInd+1
      }
      ctmp<-gsub("COLORNAME",lastCol,comTemp)
      ctmp<-gsub("POS1",as.character(spos),ctmp)
      ctmp<-gsub("POS2",as.character(mcInd-1),ctmp)
      #ctmp<-gsub("POS1",as.character(spos-1),ctmp) #convert to python indexing
      #ctmp<-gsub("POS2",as.character(mcInd-2),ctmp)
      holdAllCommands<-c(holdAllCommands,ctmp)
      
      holdAllCommands <- holdAllCommands[nchar(holdAllCommands) > 1]
      holdAllCommands<-c(paste0("color ",ehHex,";"),holdAllCommands)
      writeLines(holdAllCommands, con = myOutList[fileNum])
    }
    if(verbose){print("Done!")}
  }

