#' Combines data for probeset summarization
#' @param fsetid fsetid
#' @return a Indices whhich are used for probeset summarization 
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
getSingleProbeSetSize <- function(fsetid) {
    uniqueProbeSets <- unique(fsetid)
    endIndex <- length(uniqueProbeSets)
    orderedProbesets <-  order(uniqueProbeSets)
    NumOfPInPS <- as.numeric(table(fsetid))
    NumOfPInPS[orderedProbesets] <- NumOfPInPS
    
    index.mat <- matrix(NA, endIndex, 2)
    index.mat[,1] <- cumsum(c(1,NumOfPInPS))[1:endIndex] 
    index.mat[,2] <- cumsum(NumOfPInPS)
    colnames(index.mat) <- c("start", "end")
    PS.info <- data.frame(index.mat)
    return(PS.info)
}


#' Finds SNPs which belong to one fragment
#' @param fragLength fragLength
#' @return windows for fragments
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
getFragmentSet <- function(fragLength) {
    ##FIXME: improve implementation
    tmpIdx <- which(is.na(fragLength))
    if (length(tmpIdx) != 0) {
        fragLength[tmpIdx] <- 1200
    }
    
    y <- diff(fragLength, lag = 1)
    y[y != 0] <- 1
    z <- diff(y, lag = 1)
    
    x03 <- cbind(fragLength, c(1, y), c(-1, z, 1))
    
    a01 <- which(x03[, 3] == -1)
    a02 <- which(x03[, 3] == 1)
    
    x <- cbind(a01, a02)
    colnames(x) <- c("start", "end")
    x <- data.frame(x)
    return(x)    
}


#' Combines neighbouring locations to windows
#' @param phInf The locations on the chromosomes.
#' @param windowSize Size of how many Locations should be combined.
#' @param overlap States if the windows should overlap.
#' @return Indices for summarization 
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @examples
#' ## create toy physical data
#' sizeTmp <- 30
#' phInf <- data.frame(
#'         chrom = rep("15", sizeTmp),
#'         start = seq(from = 1, by = 300, length.out = sizeTmp), 
#'         end = seq(from = 3600, by = 300, length.out = sizeTmp),
#'         man_fsetid = paste("SNP_A-", seq(sizeTmp)+1000, sep = ""))
#' summarizeWindowStd(phInf)
summarizeWindowStd <- function(phInf, windowSize = 3, overlap = TRUE) {
    endIndex <- length(unique(phInf$"man_fsetid"))
    
    if (windowSize == 1) {
        print("used different function")
    } else {
        if(overlap){
            start <- 1:endIndex
            end <- (1:endIndex) + (windowSize - 1)
            max.index <- which(end==endIndex)
        } else {
            start <- seq(1, endIndex, by=(windowSize))
            end <- seq(windowSize,endIndex, by=(windowSize))
            max.index <- min(length(start), length(end))
        }
        startTmp <- as.numeric(start[1:max.index])
        endTmp <- as.numeric(end[1:max.index])
        PS.info <- data.frame(
                start = startTmp, 
                end = endTmp)
    }
    return(PS.info)
}

#' Combines neighbouring locations to windows
#' @param phInf The locations on the chromosomes.
#' @param fixedBps Size of the window in basepairs.
#' @param upperLimit Maximal number of neigbouring locations to combine.
#' @return Indices for summarization 
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @examples
#' ## create toy physical data
#' sizeTmp <- 30
#' phInf <- data.frame(
#'         chrom = rep("15", sizeTmp),
#'         start = seq(from = 1, by = 300, length.out = sizeTmp), 
#'         end = seq(from = 3600, by = 300, length.out = sizeTmp),
#'         man_fsetid = paste("SNP_A-", seq(sizeTmp)+1000, sep = ""))
#' summarizeWindowBps(phInf)
summarizeWindowBps <- function(phInf, fixedBps = 10000, upperLimit = 6) {
    startIdx <- c()
    endIdx <- c()
    probesInChrom <- as.numeric(table(phInf$chrom))  # occurancies per chromosome 
    numberOfProbesInChrom <- cumsum(probesInChrom) # vector!!!
    chr <- unique(phInf$chrom)
    
    winStart <- list()
    winEnd <- list()
    for(iCnt in 1:length(chr)){
        if(iCnt == 1){
            startIndexOfChr <- 1
        } else {
            startIndexOfChr <- numberOfProbesInChrom[(iCnt - 1)] + 1
            
        } 
        
        phInfoSelectedChromosome <- 
                phInf[startIndexOfChr:numberOfProbesInChrom[iCnt], ]    

        phPos <- as.numeric(as.character(phInfoSelectedChromosome[, 2]))
        diffPhPositions <- c(diff(phPos, lag = 1))
        
        startIdx <- 1:length(phPos)
        endIdx <- vector(length = length(phPos))
        for (startPos in 1:length(phPos - 1)) {
            endPos <- startPos
            check <- phPos[startPos] + fixedBps
            while (phPos[endPos + 1] <= check & endPos < length(phPos)) {    
                endPos <- endPos + 1
            }
            endIdx[startPos] <- endPos    
        }
        windows <- cbind(startIdx, endIdx)
        windows <- windows[-which(apply(windows, 1, diff) < 2), ]
        windows <- windows + (startIndexOfChr - 1)
        winStart <- append(winStart, windows[, 1])
        winEnd <- append(winEnd, windows[, 2])
    }
    winStart <- unlist(winStart)
    winEnd <- unlist(winEnd)
    
    ## truncate 
    truncIndices <- which(winEnd - winStart >= upperLimit)
    winEnd[truncIndices] <- winStart[truncIndices] + upperLimit
    return(data.frame(
                    start = as.numeric(winStart), 
                    end = as.numeric(winEnd)))
}
