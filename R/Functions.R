
utils::globalVariables(c("coverage", "HighestPeakReadCoverage",
    "start","outputFrame"))


#' inputCheck
#'
#' Checks user inputs to ensure PIPETS can run with given parameters
#' @importFrom utils read.delim file_test
#' @param inputBedFile See PIPETS_Run for full explanation
#' @param readLength See PIPETS_Run for full explanation
#' @param slidingWindowSize See PIPETS_Run for full explanation
#' @param slidingWindowMovementDistance See PIPETS_Run for full explanation
#' @param threshAdjust See PIPETS_Run for full explanation
#' @param pValue See PIPETS_Run for full explanation
#' @param topEndPercentage See PIPETS_Run for full explanation
#' @param adjacentPeakDistance See PIPETS_Run for full explanation
#' @param peakCondensingDistance See PIPETS_Run for full explanation
#' @return Returns kicker variable that will stop PIPETS if error is detected
#' @noRd
#'
inputCheck <- function(inputBedFile,readLength,slidingWindowSize,
                       slidingWindowMovementDistance,threshAdjust,
                       pValue,topEndPercentage,
                       adjacentPeakDistance, peakCondensingDistance){
    kicker <- 0
    if(file_test("-f", as.character(inputBedFile))){
        test <- try(read.delim(file = as.character(inputBedFile),
                               header = FALSE, stringsAsFactors = FALSE))
    } else if(!file_test("-f", as.character(inputBedFile))){
        kicker <- 1
        message("Input File Not Found")
        return(kicker)
    }
    while(kicker ==0 ){
        if(is.na(test)){
            message("Input file path not found")
            kicker <- 1
        }
        else if(nrow(test) < 100){
            message("Input file is too short")
            kicker <- 1
        }
        else if(ncol(test) < 4){
            message("Not enough columns in input file")
            kicker <- 1
        }
        if(slidingWindowSize ==0 | slidingWindowMovementDistance ==0 |
           slidingWindowMovementDistance == 0 | threshAdjust ==0 |
           pValue == 0 | adjacentPeakDistance ==0 | peakCondensingDistance == 0){
            message("One or more parameters is 0 and PIPETS cannot run")
            kicker <- 1
        }
    }
    return(kicker)
}





#' thresCalc
#'
#' Used to calculate a cutoff threshold for the data
#' @param rf Dataframe containing strand specific reads
#' @param threshAdjust This parameter is used to establish a global cutoff threshold informed by the data. PIPETS sorts the genomic positions of each strand from highest to lowest, and starts with the highest read coverage position and subtracts that value from the total read coverage for that strand. By default, this continues until 75% of the total read coverage has been accounted for. Increasing the percentage (e.x. 0.9) will lower the strictness of the cutoff, thus increasing the total number of significant results.
#' @param topEndPercentage This parameter is used along with threshAdjust to trim off the influence exerted by high read coverage outliers. By default, it removes the top 0.01 percent of the highest read coverage positions from the calculation of the global threshold (e.x. if there are 200 positions that make up 75% of the total reads, then this parameter will take the top 2 read coverage positions and remove them from the calculation of the global threshold). This parameter can be tuned to account for datasets with outliers that would otherwise severely skew the global threshold.
#' @return Outputs threshold used for cutoff
#' @noRd
#'
thresCalc <- function(rf, threshAdjust,topEndPercentage){
    threshCalc <- rf$coverage[order(rf$coverage, decreasing = TRUE)]
    tempMax <- sum(threshCalc) * threshAdjust
    posCount <- 1
    for(x in seq_along(threshCalc)){
        tempMax <- tempMax - threshCalc[[x]]
        if(tempMax > 0){
            posCount <- posCount + 1
        }
        else if (tempMax <= 0 ){
            posCount <- posCount - 1
            break()
        }
    }
    rmTopEnd <- round(posCount * topEndPercentage)
    threshold <- sum(threshCalc[rmTopEnd:posCount])/(posCount - rmTopEnd)
    return(threshold)
}


#' consecutiveCheck
#' Identifies significant positions that are consecutive
#' @param OF Dataframe with significant positions (after initial Poisson test)
#' @param OMF Dataframe that will be the output, specified in function this is nested in
#' @param aPD Adjacent Peak distance speicified by user in full run
#' @return Returns list of merged termination peaks
#' @noRd
#'
consecutiveCheck <- function(OF, OMF, aPD){
    x <- 1
    peakFrameCoord <- 1
    for(x in seq_along(outputFrame[,1])){
        tempSubset <- ""
        if(x < nrow(OF) & (OF$stop[x+1] - OF$stop[x]) <= aPD ){
            TPH <- rbind(TPH,OF[x,])
        }
        else if(x == nrow(OF)){
            if((OF$stop[x] - OF$stop[(x-1)]) <= aPD) {
                TPH <- rbind(TPH,OF[x,])
                OMF$chrom[peakFrameCoord] <- OF$chrom[1]
                OMF$strand[peakFrameCoord] <- OF$strand[1]
                tS <- subset(TPH, coverage == max(TPH$coverage))
                tS <- subset(tS, subset = !duplicated(tS[c("coverage")]),
                    select = c("chrom", "start", "stop","coverage","strand"))
                OMF$HighestPeak[peakFrameCoord] <- tS$stop[1]
                OMF$HighestPeakReadCoverage[peakFrameCoord] <-tS$coverage[1]
                OMF$LowestPeakCoord[peakFrameCoord] <- min(TPH$stop)
                OMF$HighestPeakCoord[peakFrameCoord] <- max(TPH$stop)
                peakFrameCoord <- peakFrameCoord + 1
                OMF[peakFrameCoord,] <- NA
                TPH <- as.data.frame(matrix(nrow = 0, ncol = ncol(OF)))
                colnames(TPH) <- colnames(OF)
            }
            else if(!(OF$stop[x] - OF$stop[x-1]) <= aPD){
                OMF[peakFrameCoord,] <- c(OF$chrom[1],
                    OF$strand[x], OF$stop[x],OF$coverage[x],
                    OF$stop[x],OF$stop[x])
            }
        }
        else {
            TPH <- rbind(TPH,OF[x,])
            OMF$chrom[peakFrameCoord] <- OF$chrom[1]
            OMF$strand[peakFrameCoord] <- OF$strand[1]
            tS <- subset(TPH, coverage == max(TPH$coverage))
            tS <- subset(tS,subset = !duplicated(tempSubset[c("coverage")]),
                    select = c("chrom", "start", "stop","coverage","strand"))
            OMF$HighestPeak[peakFrameCoord] <- tS$stop[1]
            OMF$HighestPeakReadCoverage[peakFrameCoord] <-tS$coverage[1]
            OMF$LowestPeakCoord[peakFrameCoord] <- min(TPH$stop)
            OMF$HighestPeakCoord[peakFrameCoord] <- max(TPH$stop)
            peakFrameCoord <- peakFrameCoord + 1
            OMF[peakFrameCoord,] <- NA
            TPH <- as.data.frame(matrix(nrow = 0, ncol = ncol(OF)))
            colnames(TPH) <- colnames(OF)
        }
    }
    return(OMF)
}


#' consecutivePeakCheck
#' Combines Proximal Peaks
#' @param OMF OutputMergeFrame from function this is nested in
#' @param SWR Output of file that is created in nested function
#' @param pCD peak condensing distance established in full run
#' @return Returns significant peaks to be fixed and output after
#' @noRd
#'
consecutivePeakCheck <- function(OMF, SWR, pCD){
    PFC <- 1
    for(x in seq_along(OMF[,1])){
        tempSubset <- ""
        if(x < nrow(OMF) & (OMF[x+1,5] - OMF[x,6]) <= pCD){
            TWH <- rbind(TWH,OMF[x,])}
        else if(x == nrow(OMF)){
            if((OMF[x,6] - OMF[(x-1),5]) <= pCD){
                TWH <- rbind(TWH,OMF[x,])
                SWR$chrom[PFC] <- OMF$chrom[1]
                SWR$strand[PFC] <- OMF$strand[1]
                tS <- subset(TWH,HighestPeakReadCoverage == max(TWH[,4]))
                tS <- subset(tS,subset =
                    !duplicated(tS[c("HighestPeakReadCoverage")]),select =
                    c("chrom","strand","HighestPeak","HighestPeakReadCoverage",
                    "LowestPeakCoord","HighestPeakCoord"))
                SWR$HighestPeak[PFC] <- tS$HighestPeak[1]
                SWR$HighestPeakReadCoverage[PFC] <-tS$HighestPeakReadCoverage[1]
                SWR$LowestPeakCoord[PFC] <- min(TWH$LowestPeakCoord)
                SWR$HighestPeakCoord[PFC] <- max(TWH$HighestPeakCoord)
                PFC <- PFC + 1
                SWR[PFC,] <- NA
                TWH <- as.data.frame(matrix(nrow = 0, ncol = ncol(OMF)))
                colnames(TWH) <- colnames(OMF)}
            else if(!(OMF$LowestPeakCoord[x]-OMF$HighestPeakCoord[x-1])<= pCD){
                SWR[PFC,] <- c(OMF$chrom[1],
                    OMF$strand[x],OMF$HighestPeak[x],
                    OMF$HighestPeakReadCoverage[x],OMF$HighestPeak[x],
                    OMF$HighestPeak[x])}
        }
        else {
            TWH <- rbind(TWH,OMF[x,])
            SWR$chrom[PFC] <- OMF$chrom[1]
            SWR$strand[PFC] <- OMF$strand[1]
            tS <- subset(TWH,HighestPeakReadCoverage == max(TWH[,4]))
            tS <- subset(tS,
                subset = !duplicated(tS[c("HighestPeakReadCoverage")]),
                select = c("chrom","strand","HighestPeak",
                "HighestPeakReadCoverage","LowestPeakCoord","HighestPeakCoord"))
            SWR$HighestPeak[PFC] <- tS$HighestPeak[1]
            SWR$HighestPeakReadCoverage[PFC] <- tS$HighestPeakReadCoverage[1]
            SWR$LowestPeakCoord[PFC] <- min(TWH$LowestPeakCoord)
            SWR$HighestPeakCoord[PFC] <- max(TWH$HighestPeakCoord)
            PFC <- PFC + 1
            SWR[PFC,] <- NA
            TWH <- as.data.frame(matrix(nrow = 0, ncol = ncol(OMF)))
            colnames(TWH) <- colnames(OMF)}
    }
    return(SWR)
}




#' Bed_Split
#' First step of PIPETS. Takes input Bed files and splits them by strand while also assigning read coverage to each genomic position
#' @title Split Input Bed Data By Strand
#' @importFrom dplyr arrange distinct
#' @importFrom stats aggregate ppois complete.cases
#' @importFrom utils write.csv write.table read.table read.delim
#' @param inputBedFile Input BED file that is not strand split. For PIPETS, the first column must be the chromosome name, the second column must be the start coordinate, the third column must be the stop coordinate, and the 6th column must be the strand. Columns 4 and 5 must be present but their information will not be used.
#' @param readLength The user must input the expected length of each “proper” read from the 3’-seq protocol. PIPETS gives a +/- 2 range for detecting reads to be used in the input BED files based on distributions of read coverage in parameter testing.
#' @return Returns a list containing the Plus Strand Reads, the Minus Strand Reads, and the user defined name for the files. Also writes out the strand split bed files to the project directory.
#' @examples
#' ## Split input bed file into stranded files without running PIPETS
#' ## The user will be prompted to enter a string to name output files
#' Bed_Split(inputBedFile="Test_Data.bed", readLength=58)
#' @noRd
#'
Bed_Split <- function(inputBedFile,readLength){
    message("+-----------------------------------+")
    OutputFileName <- readline(prompt="Please Enter the Sample Name: ")
    message("Splitting Input Bed File By Strand")
    rB <- read.delim(file = as.character(inputBedFile),
        header = FALSE, stringsAsFactors = FALSE)
    startBed <- rB[,c(1,2,3)]
    if(length(names(rB[grepl("-", rB[1,], fixed = TRUE)])) >0){
        startBed$V4 <- rB[,names(rB[grepl('-', rB[1,], fixed = TRUE)])]
    } else if (length(names(rB[grepl("+", rB[1,], fixed = TRUE)])) >0){
        startBed$V4 <- rB[,names(rB[grepl('+', rB[1,], fixed = TRUE)])]
  }
    colnames(startBed) <- c("chrom","start","stop","strand")
    startBed$coverage <- 0
    startBed <- startBed[ , c("chrom", "start", "stop", "coverage", "strand")]
    PSR <- startBed[startBed$strand %in% "+",]
    PSR <- PSR[(PSR$stop-PSR$start) %in% (readLength-2):(readLength+2),]
    PSC <- as.data.frame(table(PSR$stop))
    PSR <- distinct(PSR, stop, .keep_all = TRUE)
    PSR <- arrange(PSR, stop)
    PSR$coverage <- PSC$Freq[match(PSR$stop,PSC$Var1)]
    tempChrom <- PSR$chrom[1]
    tempStrand <- PSR$strand[1]
    PSR <- aggregate(coverage~start,PSR,sum)
    PSR$chrom <- tempChrom
    PSR$stop <- PSR$start + readLength
    PSR$strand <- tempStrand
    PSR <- PSR[,c("chrom","start","stop","coverage","strand")]
    MSR <- startBed[startBed$strand %in% "-",]
    MSR <- MSR[(MSR$stop-MSR$start) %in% (readLength-2):(readLength+2),]
    MSC <- as.data.frame(table(MSR$start))
    MSR <- distinct(MSR, start, .keep_all = TRUE)
    MSR <- arrange(MSR, start)
    MSR$coverage <- MSC$Freq[match(MSR$start,MSC$Var1)]
    tempChrom <- MSR$chrom[1]
    tempStrand <- MSR$strand[1]
    MSR <- aggregate(coverage~stop,MSR,sum)
    MSR$chrom <- tempChrom
    MSR$start <- MSR$stop - readLength
    MSR$strand <- tempStrand
    MSR <- MSR[,c("chrom","start","stop","coverage","strand")]
    write.table(PSR,
    file = paste(as.character(OutputFileName),"PlusStrandCounts.bed", sep = "_")
    ,quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(MSR,
    file = paste(as.character(OutputFileName),
        "MinusStrandCounts.bed", sep = "_"),
    quote = FALSE, row.names = FALSE, col.names = FALSE)
    return(list(OutputFileName,PSR,MSR))
}


#' TopStrand_InitialPoisson
#' Poisson Significant Peak Identification Test for the Top Strand Data
#' @importFrom dplyr arrange distinct
#' @importFrom stats aggregate ppois complete.cases
#' @importFrom utils write.csv write.table read.table read.delim
#' @param MinusStrandReads Minus Strand Read DataFrame from the Bed_Split method. The minus strand reads inform the Top Strand termination signal
#' @param slidingWindowSize This parameter establishes the distance up and down stream of each position that a sliding window will be created around. The default value is 25, and this will result in a sliding window of total size 51 (25 upstream + position (1) + 25 downstream).
#' @param slidingWindowMovementDistance This parameter sets the distance that the sliding window will be moved. By default, it is set to move by half of the sliding window size in order to ensure that almost every position in the data is tested twice.
#' @param threshAdjust This parameter is used to establish a global cutoff threshold informed by the data. PIPETS sorts the genomic positions of each strand from highest to lowest, and starts with the highest read coverage position and subtracts that value from the total read coverage for that strand. By default, this continues until 75% of the total read coverage has been accounted for. Increasing the percentage (e.x. 0.9) will lower the strictness of the cutoff, thus increasing the total number of significant results.
#' @param pValue Choose the minimum pValue that the Poisson distribution test must pass in order to be considered significant
#' @param topEndPercentage This parameter is used along with threshAdjust to trim off the influence exerted by high read coverage outliers. By default, it removes the top 0.01 percent of the highest read coverage positions from the calculation of the global threshold (e.x. if there are 200 positions that make up 75% of the total reads, then this parameter will take the top 2 read coverage positions and remove them from the calculation of the global threshold). This parameter can be tuned to account for datasets with outliers that would otherwise severely skew the global threshold.
#' @return Returns a dataframe with all genomic positions that were identified as having significant read coverage.
#' @noRd
#'
TopStrand_InitialPoisson <- function(MinusStrandReads,slidingWindowSize = 25,
    slidingWindowMovementDistance = 25,threshAdjust = 0.75,pValue = 0.005,
    topEndPercentage= 0.01){
    MSR <- MinusStrandReads
    SWMD <- slidingWindowMovementDistance
    SWS <- slidingWindowSize
    outputFrame <- as.data.frame(matrix(nrow = 0, ncol = 5))
    colnames(outputFrame) <- colnames(MSR)
    SWF <- as.data.frame(matrix(nrow =
        MSR[nrow(MSR),"stop"], ncol = 2))
    colnames(SWF) <- c("position","coverage")
    SWF$coverage <- 0
    SWF$position <- seq_along(SWF[,1])
    posThreshold <- thresCalc(rf = MSR, threshAdjust = threshAdjust,
        topEndPercentage = topEndPercentage)
    message(paste("Top Strand Cutoff", posThreshold, sep = " "))
    nm <- "coverage"
    SWF[nm] <- lapply(nm, function(x) MSR[[x]][match(SWF$position, MSR$stop)])
    SWF$coverage[is.na(SWF$coverage)] <- 0
    z <- 1+ SWS
    breakCond <- 0
    numericBreakCond <- 0
    while(breakCond == 0){
        sumOfWindowCoverage <- sum(SWF$coverage[(z-SWS):(z+SWS)])
        AOW <- sumOfWindowCoverage/((2*SWS)+1)
        for(y in (z-SWS):(z+SWS)){
            if(SWF$coverage[y] >= 10){
                probabilityAtY <- (1- ppois(q = SWF$coverage[y], lambda = AOW))
                if(SWF$coverage[y] >= posThreshold & probabilityAtY <= pValue){
                    outputFrame <- rbind(MSR[MSR$stop == SWF$position[y],],
                    outputFrame)
                }
            }
        }
        if(numericBreakCond == 1){
            breakCond <- 1
        }
        if(numericBreakCond == 0 & (z+SWMD) >= (nrow(SWF)-SWS)){
            z <- (nrow(SWF)-SWS)
            numericBreakCond <- 1
        }
        if(numericBreakCond == 0 & (z+SWMD) < (nrow(SWF)-SWS)){
            z <- (z + SWMD)
        }
    }
    outputFrame <- outputFrame[!duplicated(outputFrame),]
    outputFrame <- outputFrame[order(outputFrame$stop),]
    rownames(outputFrame) <- seq_along(outputFrame[,1])
    return(outputFrame)
}

#' TopStrand_InitialCondense
#' Takes the significant top strand positions from TopStrand_InitialPoisson and condenses all proximal positions into termination "peaks"
#' @importFrom dplyr arrange distinct
#' @importFrom stats aggregate ppois complete.cases
#' @importFrom utils write.csv write.table read.table read.delim
#' @param TopInititalPoisson Uses the output of TopStrand_InitialPoisson as the input
#' @param adjacentPeakDistance adjacentPeakDistance During the peak condensing step, this parameter is used to define “adjacent” for significant genomic positions. This is used to identify initial peak structures in the data. By default this value is set to 2 to ensure that single instances of loss of signal are not sufficient to prevent otherwise contiguous peak signatures from being combined.
#' @return Returns a dataframe that contains the list of termination "peaks" for the top strand
#' @noRd
#'
TopStrand_InitialCondense <- function(TopInititalPoisson,
    adjacentPeakDistance = 2){
    outputFrame <- TopInititalPoisson
    TPH <- as.data.frame(matrix(nrow = 0, ncol = ncol(outputFrame)))
    colnames(TPH) <- colnames(outputFrame)
    OMF <- as.data.frame(matrix(nrow = 1, ncol = 6))
    colnames(OMF) <- c("chrom","strand","HighestPeak",
    "HighestPeakReadCoverage","LowestPeakCoord",
    "HighestPeakCoord")
    OMF <- consecutiveCheck(OF = outputFrame, OMF = OMF,
        aPD = adjacentPeakDistance)
    OMF_Top <- OMF[complete.cases(OMF),]
    rownames(OMF_Top) <- seq_along(OMF_Top[,1])
    OMF_Top$LowestPeakCoord <- as.numeric(OMF_Top$LowestPeakCoord)
    OMF_Top$HighestPeakCoord <- as.numeric(OMF_Top$HighestPeakCoord)
    OMF_Top$HighestPeakReadCoverage <-
        as.numeric(OMF_Top$HighestPeakReadCoverage)
    OMF_Top$HighestPeak <- as.numeric(OMF_Top$HighestPeak)
    return(OMF_Top)
}

#' TopStrand_SecondaryCondense
#' Takes the list of termination "peaks" from TopStrand_InitialCondense and condenses peaks that are proximal to each other
#' @importFrom dplyr arrange distinct
#' @importFrom stats aggregate ppois complete.cases
#' @importFrom utils write.csv write.table read.table read.delim
#' @param TopInititalCondense Uses the output of TopStrand_InitialCondense as the input
#' @param peakCondensingDistance peakCondensingDistance Following the initial peak condensing step, this parameter is used to identify peak structures in the data that are close enough to be considered part of the same termination signal. In testing, we have not identified cases in which two distinct termination signals so proximal that the default parameters incorrectly combine the signals together.
#' @param OutputFileName A string that will be used to identify printed results files. When run with PIPETS_FullRun, this string will be input by the user in the beginning
#' @return The method writes a csv file in the project directory that contains the results for the Top Strand termination peaks
#' @noRd
#'
TopStrand_SecondaryCondense <- function(TopInititalCondense,
    peakCondensingDistance = 20,
    OutputFileName){
    OFN <- OutputFileName
    OMF_Top <- TopInititalCondense
    TWH <- as.data.frame(matrix(nrow = 0, ncol = ncol(OMF_Top)))
    colnames(TWH) <- colnames(OMF_Top)
    SWR_Top <- as.data.frame(matrix(nrow = 1, ncol = 6))
    colnames(SWR_Top) <- c("chrom","strand","HighestPeak",
        "HighestPeakReadCoverage","LowestPeakCoord","HighestPeakCoord")
    SWR_Top <- consecutivePeakCheck(OMF_Top, SWR_Top, peakCondensingDistance)
    SWR_Top <- SWR_Top[complete.cases(SWR_Top),]
    SWR_Top <- SWR_Top[,-2]
    write.csv(SWR_Top,
        file = paste(as.character(OFN),"topStrandResults.csv", sep = "_"),
        row.names = FALSE)
}

#' CompStrand_InitialPoisson
#' Poisson Significant Peak Identification Test for the Complement Strand Data
#' @importFrom dplyr arrange distinct
#' @importFrom stats aggregate ppois complete.cases
#' @importFrom utils write.csv write.table read.table read.delim
#' @param PlusStrandReads Plus Strand Read DataFrame from the Bed_Split method. The Plus strand reads inform the Complement Strand termination signal
#' @param slidingWindowSize This parameter establishes the distance up and down stream of each position that a sliding window will be created around. The default value is 25, and this will result in a sliding window of total size 51 (25 upstream + position (1) + 25 downstream).
#' @param slidingWindowMovementDistance This parameter sets the distance that the sliding window will be moved. By default, it is set to move by half of the sliding window size in order to ensure that almost every position in the data is tested twice.
#' @param threshAdjust This parameter is used to establish a global cutoff threshold informed by the data. PIPETS sorts the genomic positions of each strand from highest to lowest, and starts with the highest read coverage position and subtracts that value from the total read coverage for that strand. By default, this continues until 75% of the total read coverage has been accounted for. Increasing the percentage (e.x. 0.9) will lower the strictness of the cutoff, thus increasing the total number of significant results.
#' @param pValue Choose the minimum pValue that the Poisson distribution test must pass in order to be considered significant
#' @param topEndPercentage This parameter is used along with threshAdjust to trim off the influence exerted by high read coverage outliers. By default, it removes the top 0.01 percent of the highest read coverage positions from the calculation of the global threshold (e.x. if there are 200 positions that make up 75% of the total reads, then this parameter will take the top 2 read coverage positions and remove them from the calculation of the global threshold). This parameter can be tuned to account for datasets with outliers that would otherwise severely skew the global threshold.
#' @return Returns a dataframe with all genomic positions that were identified as having significant read coverage.
#' @noRd
#'
CompStrand_InitialPoisson <- function(PlusStrandReads,slidingWindowSize = 25,
    slidingWindowMovementDistance = 25,threshAdjust = 0.75,pValue = 0.005,
    topEndPercentage= 0.01){
    SWMD <- slidingWindowMovementDistance
    SWS <- slidingWindowSize
    PSR <- PlusStrandReads
    OF <- as.data.frame(matrix(nrow = 0, ncol = 5))
    colnames(OF) <- colnames(PSR)
    SWF <- as.data.frame(matrix(nrow = PSR[nrow(PSR),"start"], ncol = 2))
    colnames(SWF) <- c("position","coverage")
    SWF$coverage <- 0
    SWF$position <- seq_along(SWF[,1])
    compThreshold <- thresCalc(rf = PSR,threshAdjust = threshAdjust,
        topEndPercentage = topEndPercentage)
    message(paste("Complement Strand Cutoff", compThreshold, sep = " "))
    nm <- "coverage"
    SWF[nm] <- lapply(nm, function(x) PSR[[x]][match(SWF$position, PSR$start)])
    SWF$coverage[is.na(SWF$coverage)] <- 0
    z <- 1+ slidingWindowSize
    breakCond <- 0
    numericBreakCond <- 0
    while(breakCond == 0){
        sumOfWindowCoverage <- sum(SWF$coverage[(z-SWS):(z+SWS)])
        AOW <- sumOfWindowCoverage/((2*SWS)+1)
        for(y in (z-SWS):(z+SWS)){
            if(SWF$coverage[y] >= 10){
                probabilityAtY <- (1- ppois(q = SWF$coverage[y], lambda = AOW))
                if(SWF$coverage[y] >= compThreshold & probabilityAtY <= pValue){
                    OF <- rbind(PSR[PSR$start == SWF$position[y],], OF)
                }
            }
        }
        if(numericBreakCond == 1){
            breakCond <- 1
        }
        if(numericBreakCond == 0 & (z+SWMD) >= (nrow(SWF)-SWS)){
            z <- (nrow(SWF)-SWS)
            numericBreakCond <- 1
        }
        if(numericBreakCond == 0 & (z+SWMD) < (nrow(SWF)-SWS)){
            z <- (z + SWMD)
        }
    }
    OF <- OF[!duplicated(OF),]
    OF <- OF[order(OF$stop),]
    rownames(OF) <- seq_along(OF[,1])
    return(OF)
}

#' CompStrand_InitialCondense
#' Takes the significant top strand positions from CompStrand_InitialPoisson and condenses all proximal positions into termination "peaks"
#' @importFrom dplyr arrange distinct
#' @importFrom stats aggregate ppois complete.cases
#' @importFrom utils write.csv write.table read.table read.delim
#' @param CompInitialPoisson Uses the output of CompStrand_InitialPoisson as the input
#' @param adjacentPeakDistance adjacentPeakDistance During the peak condensing step, this parameter is used to define “adjacent” for significant genomic positions. This is used to identify initial peak structures in the data. By default this value is set to 2 to ensure that single instances of loss of signal are not sufficient to prevent otherwise contiguous peak signatures from being combined.
#' @return Returns a dataframe that contains the list of termination "peaks" for the complement strand
#' @noRd
#'
CompStrand_InitialCondense <- function(CompInitialPoisson,
    adjacentPeakDistance = 2){
    outputFrame <- CompInitialPoisson
    tempPeakHold <- as.data.frame(matrix(nrow = 0, ncol = ncol(outputFrame)))
    colnames(tempPeakHold) <- colnames(outputFrame)
    OMF <- as.data.frame(matrix(nrow = 1, ncol = 6))
    colnames(OMF) <- c("chrom","strand","HighestPeak",
        "HighestPeakReadCoverage",
        "LowestPeakCoord","HighestPeakCoord")
    OMF <- consecutiveCheck(OF = outputFrame, OMF = OMF,
        aPD = adjacentPeakDistance)
    OMF_Comp <- OMF[complete.cases(OMF),]
    rownames(OMF_Comp) <- seq_along(OMF_Comp[,1])
    OMF_Comp$LowestPeakCoord <- as.numeric(OMF_Comp$LowestPeakCoord)
    OMF_Comp$HighestPeakCoord <- as.numeric(OMF_Comp$HighestPeakCoord)
    OMF_Comp$HighestPeakReadCoverage <-
        as.numeric(OMF_Comp$HighestPeakReadCoverage)
    OMF_Comp$HighestPeak <- as.numeric(OMF_Comp$HighestPeak)
    return(OMF_Comp)
}

#' CompStrand_SecondaryCondense
#' Takes the list of termination "peaks" from CompStrand_InitialCondense and condenses peaks that are proximal to each other
#' @importFrom dplyr arrange distinct
#' @importFrom stats aggregate ppois complete.cases
#' @importFrom utils write.csv write.table read.table read.delim
#' @param CompInitialCondense Uses the output of CompStrand_InitialCondense as the input
#' @param peakCondensingDistance peakCondensingDistance Following the initial peak condensing step, this parameter is used to identify peak structures in the data that are close enough to be considered part of the same termination signal. In testing, we have not identified cases in which two distinct termination signals so proximal that the default parameters incorrectly combine the signals together.
#' @param OutputFileName A string that will be used to identify printed results files. When run with PIPETS_FullRun, this string will be input by the user in the beginning
#' @return The method writes a csv file in the project directory that contains the results for the Complement Strand termination peaks
#' @noRd
#'
CompStrand_SecondaryCondense <- function(CompInitialCondense,
    peakCondensingDistance = 20,
    OutputFileName){
    OFN <- OutputFileName
    OMF_Comp <- CompInitialCondense
    TWH <- as.data.frame(matrix(nrow = 0, ncol = ncol(OMF_Comp)))
    colnames(TWH) <- colnames(OMF_Comp)
    SWR_Comp <- as.data.frame(matrix(nrow = 1, ncol = 6))
    colnames(SWR_Comp) <- c("chrom","strand","HighestPeak",
        "HighestPeakReadCoverage","LowestPeakCoord","HighestPeakCoord")
    SWR_Comp <- consecutivePeakCheck(OMF_Comp, SWR_Comp, peakCondensingDistance)
    SWR_Comp <- SWR_Comp[complete.cases(SWR_Comp),]
    SWR_Comp <- SWR_Comp[,-2]
    write.csv(SWR_Comp, file = paste(as.character(OFN),
        "CompStrandResults.csv", sep = "_"),  row.names = FALSE)
}


#' PIPETS_FullRun
#' @title Analyze 3'-seq Data with PIPETS
#' Poisson Identification of PEaks from Term-Seq data. This is the full run method that begins with input Bed file and returns the strand split results
#' @importFrom dplyr arrange distinct
#' @importFrom stats aggregate ppois complete.cases
#' @importFrom utils write.csv write.table read.table read.delim
#' @param inputBedFile Input BED file that is not strand split. For PIPETS, the first column must be the chromosome name, the second column must be the start coordinate, the third column must be the stop coordinate, and the 6th column must be the strand. Columns 4 and 5 must be present but their information will not be used.
#' @param readLength The user must input the expected length of each “proper” read from the 3’-seq protocol. PIPETS gives a +/- 2 range for detecting reads to be used in the input BED files based on distributions of read coverage in parameter testing.
#' @param slidingWindowSize This parameter establishes the distance up and down stream of each position that a sliding window will be created around. The default value is 25, and this will result in a sliding window of total size 51 (25 upstream + position (1) + 25 downstream).
#' @param slidingWindowMovementDistance This parameter sets the distance that the sliding window will be moved. By default, it is set to move by half of the sliding window size in order to ensure that almost every position in the data is tested twice.
#' @param adjacentPeakDistance During the peak condensing step, this parameter is used to define “adjacent” for significant genomic positions. This is used to identify initial peak structures in the data. By default this value is set to 2 to ensure that single instances of loss of signal are not sufficient to prevent otherwise contiguous peak signatures from being combined.
#' @param peakCondensingDistance Following the initial peak condensing step, this parameter is used to identify peak structures in the data that are close enough to be considered part of the same termination signal. In testing, we have not identified cases in which two distinct termination signals so proximal that the default parameters incorrectly combine the signals together.
#' @param threshAdjust This parameter is used to establish a global cutoff threshold informed by the data. PIPETS sorts the genomic positions of each strand from highest to lowest, and starts with the highest read coverage position and subtracts that value from the total read coverage for that strand. By default, this continues until 75% of the total read coverage has been accounted for. Increasing the percentage (e.x. 0.9) will lower the strictness of the cutoff, thus increasing the total number of significant results.
#' @param pValue Choose the minimum pValue that the Poisson distribution test must pass in order to be considered significant
#' @param topEndPercentage This parameter is used along with threshAdjust to trim off the influence exerted by high read coverage outliers. By default, it removes the top 0.01 percent of the highest read coverage positions from the calculation of the global threshold (e.x. if there are 200 positions that make up 75% of the total reads, then this parameter will take the top 2 read coverage positions and remove them from the calculation of the global threshold). This parameter can be tuned to account for datasets with outliers that would otherwise severely skew the global threshold.
#' @examples
#' ## When run, the user will be prompted to provide a string for file names
#' ## During the run, PIPETS will output the minumum read coverage cutoff for each strand
#' ## After completion, the output files will be created in the R project directory
#'
#' ## For run with defualt strictness of analysis
#' PIPETS_FullRun(inputBedFile = "TestData.bed", readLength = 58)
#'
#' ## For a more strict run (can be run for files with high total read depth)
#' PIPETS_FullRun(inputBedFile = "TestData.bed", readLength = 58, threshAdjust = 0.6)
#'
#' ## For a less strict run (for data with low total read depth)
#' PIPETS_FullRun(inputBedFile = "TestData.bed", readLength = 58, threshAdjust = 0.9)
#'
#' @return PIPETS outputs strand specific results files as well as strand specific bed files to the directory that the R project is in.
#' @export
#'
PIPETS_FullRun <- function(inputBedFile,readLength,slidingWindowSize = 25,
    slidingWindowMovementDistance = 25,threshAdjust = 0.75,
    pValue = 0.005,topEndPercentage= 0.01,
    adjacentPeakDistance = 2, peakCondensingDistance = 20){
    kicker <- inputCheck(inputBedFile,readLength,slidingWindowSize,
                         slidingWindowMovementDistance,threshAdjust,
                         pValue,topEndPercentage,
                         adjacentPeakDistance, peakCondensingDistance)
    if(kicker ==1){
        message("PIPETS cannot run, see above error")
        return()
    }
    AllReads <- Bed_Split(inputBedFile, readLength)
    message("+-----------------------------------+")
    message("Performing Top Strand Analysis")
    TopInititalPoisson <- TopStrand_InitialPoisson(
        MinusStrandReads <- AllReads[[3]],slidingWindowSize = slidingWindowSize,
        slidingWindowMovementDistance = slidingWindowMovementDistance,
        threshAdjust = threshAdjust, pValue = pValue,
        topEndPercentage= topEndPercentage)
    TopInititalCondense <- TopStrand_InitialCondense(
        TopInititalPoisson = TopInititalPoisson,
        adjacentPeakDistance = adjacentPeakDistance)
    TopStrand_SecondaryCondense(TopInititalCondense = TopInititalCondense,
        peakCondensingDistance = peakCondensingDistance,
        OutputFileName = AllReads[[1]])
    message("+-----------------------------------+")
    message("Performing Complement Strand Analysis")
    CompInitialPoisson <- CompStrand_InitialPoisson(
        PlusStrandReads = AllReads[[2]],slidingWindowSize = slidingWindowSize,
        slidingWindowMovementDistance = slidingWindowMovementDistance,
        threshAdjust = threshAdjust, pValue = pValue,
        topEndPercentage= topEndPercentage)
    CompInitialCondense <- CompStrand_InitialCondense(
        CompInitialPoisson = CompInitialPoisson,
        adjacentPeakDistance = adjacentPeakDistance)
    CompStrand_SecondaryCondense(CompInitialCondense = CompInitialCondense,
        peakCondensingDistance = peakCondensingDistance,
        OutputFileName = AllReads[[1]])
}



