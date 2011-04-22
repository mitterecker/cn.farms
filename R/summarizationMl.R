#' Does summarization
#' @param object an instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' @param windowMethod Method for combination of neighbouring SNPs. 
#' Possible values are Std and Bps.
#' @param windowParam further parameters as the window size
#' @param summaryMethod allowed versions for the summarization step are: 
#' Gaussian, Variational, Exact. Default is Variational.
#' @param summaryParam summaryParam
#' @param callParam callParam
#' @param returnValues List with return values. 
#' For possible values see summaryMethod.
#' @return Some data
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @examples 
#' load(system.file("exampleData/slData.RData", package="cn.farms"))
#' windowMethod <- "std"
#' windowParam <- list()
#' windowParam$windowSize <- 5
#' windowParam$overlap <- TRUE
#' summaryMethod <- "Variational"
#' summaryParam <- list()
#' summaryParam$cyc <- c(20)
#' mlData <- mlSummarization(slData, windowMethod, windowParam, 
#'         summaryMethod, summaryParam)
#' assayData(mlData)
mlSummarization <- function(object, windowMethod, windowParam, 
        summaryMethod, summaryParam, callParam=list(runtype="ff"), 
        returnValues) {
    t00 <- Sys.time()
    summaryWindowName <- paste("summarizeWindow", paste(
                    toupper(substring(windowMethod, 1,1)), 
                    substring(windowMethod, 2),
                    sep="", collapse=" "), sep="")
    
    if (!exists(summaryWindowName)) {
        stop(paste("Unknown method (can't find function", summaryWindowName, 
                        ")"))
    } else if (callParam$runtype=="bm" & file.exists("mlData.RData")) {
        message("Multi-loci summarization has already been done")
        message("Trying to load data ...")
        load("mlData.RData", envir=globalenv())
        return(invisible())    
    }
    
    phInf <- featureData(object)[, c("chrom", "start", "man_fsetid")]
    runIdx <- do.call(summaryWindowName, c(alist(phInf@data), windowParam))
    
    phInfTmp <- data.frame(
            chrom=phInf$"chrom"[runIdx[, 1]],
            start=phInf$"start"[runIdx[, 1]],
            end=phInf$"start"[runIdx[, 2]],
            man_fsetid=phInf$"man_fsetid"[runIdx[, 2]], 
            stringsAsFactors=F)
    idxDel <- which(phInfTmp$end - phInfTmp$start < 0)
    
    if (length(idxDel) != 0) {
        phInfTmp <- phInfTmp[-idxDel, ]
        runIdx <- runIdx[-idxDel, ]
        #message("Warning: lines have been deleted due to data inconsistencies")
    }
    
    summaryMethodName <- paste("summarizeFarms", paste(
                    toupper(substring(summaryMethod, 1,1)), 
                    substring(summaryMethod, 2),
                    sep="", collapse=" "), sep="")
    
    if (!exists(summaryMethodName)) {
        stop(paste("Unknown method (can't find function", summaryMethodName, 
                        ")"))
    }
    
    myData <- do.call("callSummarize", c(alist(
                            object = assayData(object)$intensity, 
                            psInfo = runIdx,
                            batchList = phenoData(object)$batch,
                            summaryMethod = summaryMethodName, 
                            summaryParam = summaryParam, 
                            returnValues = returnValues), 
                    callParam))
    
    eSet <- new("ExpressionSet")
    assayData(eSet) <- myData
    phenoData(eSet) <- phenoData(object)
    protocolData(eSet)<- protocolData(object)
    featureData(eSet) <- new("AnnotatedDataFrame", data = phInfTmp)    
    experimentData(eSet) <- experimentData(object)
    experimentData(eSet)@other$summaryMethod <- summaryMethod
    experimentData(eSet)@other$summaryParam <- summaryParam
    experimentData(eSet)@other$type <- "mlData"
    annotation(eSet) <- annotation(object)
    
    if (callParam$runtype=="bm") {
        cat(paste(Sys.time(), "|   Saving data \n"))
        mlData <- eSet
        save(mlData, file="mlData.RData")
    }
    
    return(eSet)
}
