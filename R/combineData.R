#' Combine two ExpressionSet objects
#' 
#' Suitable for SNP or non-polymorphic data which were already processed with 
#' single locus FARMS
#' 
#' @param object01 An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' either with SNP or non-polymorphic data
#' @param object02 An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' either with SNP or non-polymorphic data
#' @param obj01Var States the variable which should be combined from the 
#' assayData slot. Default is intensity.
#' @param obj02Var States the variable which should be combined from the 
#' assayData slot. Default is intensity.
#' @param runtype Mode how the results are saved. Possible values are ff or bm. 
#' If ff is chosen the data will not be saved automatically. 
#' With bm the results will be saved permanently. 
#' @return An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @examples 
#' load(system.file("exampleData/normData.RData", package="cn.farms"))
#' experimentData(normData)@@other$annotDir <- 
#'         system.file("exampleData/annotation/pd.genomewidesnp.6/1.1.0",
#'                 package="cn.farms")
#' summaryMethod <- "Variational"
#' summaryParam <- list()
#' summaryParam$cyc <- c(10)
#' slData <- slSummarization(normData, 
#'         summaryMethod = summaryMethod, 
#'         summaryParam = summaryParam)
#' assayData(slData)$L_z[1:10, ] 
#' combData <- combineData(slData, slData)
#' combData
#' 
#' @export
combineData <- function (object01, object02, 
        obj01Var="intensity", obj02Var="intensity", runtype="ff") {

    if (runtype=="bm" & file.exists("combData.RData")) {
        message("Combining the data has already been done")
        message("Trying to load combined data ...")
        load("combData.RData", envir=globalenv())
        return(invisible())
    }
    
    dataSnpType <- experimentData(object01)@other$type
    dataNpType <- experimentData(object02)@other$type
    
    colNamesSnp <- c("chrom", "start", "end", "man_fsetid")
    colNamesNp <-  c("chrom", "start", "end", "man_fsetid")
    
    if (dataSnpType == "slData" & dataNpType == "npData") {
        a <- featureData(object01)@data[, colNamesSnp]
        b <- featureData(object02)@data[, colNamesNp]
        colnames(b) <- colNamesSnp
    } else if (dataSnpType == "npData" & dataNpType == "slData") {
        a <- featureData(object01)@data[, colNamesNp]
        b <- featureData(object02)@data[, colNamesSnp]
        colnames(a) <- colNamesSnp
    } else if (dataSnpType == "slData" & dataNpType == "slData") {
        a <- featureData(object01)@data[, colNamesSnp]
        b <- featureData(object02)@data[, colNamesSnp]
    } else {
        stop("Wrong input objects")
    }
    
    phInfTmp <- rbind(a, b) 
    idx <- order(phInfTmp[,"chrom"], phInfTmp[, "start"])
    phInf <- phInfTmp[idx, ]
    rm(phInfTmp)
    cutoff <- nrow(a)
    nbrOfProbes <- nrow(phInf)
    nbrOfSamples <- ncol(assayData(object01)[[obj01Var]])
    
#    if (runtype == "bm") {
#        intensity <- initializeBigMatrix("", nbrOfProbes, nbrOfSamples, "double")    
#    } else if (runtype == "ff") {
#        intensity <- ff(vmode="double", dim=c(nbrOfProbes, nbrOfSamples))    
#    }
    
    intensity <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, type="double",
            bmName="comb_")
    
    
    intensity[which(idx <= cutoff), ] <- assayData(object01)[[obj01Var]][]
    intensity[which(idx > cutoff), ] <- assayData(object02)[[obj02Var]][]
    
    eSet <- new("ExpressionSet")
    
    ## assay data    
    assayData(eSet) <- list(intensity=intensity)
    
    ## protocol data
    protocolData(eSet) <- protocolData(object01)
    
    ## pheno data
    phenoData(eSet) <- phenoData(object01)
    
    ## feature data
    featureData(eSet) <- new("AnnotatedDataFrame", 
            data = phInf)
    
    ## experiment data
    experimentData(eSet) <- experimentData(object01)
    experimentData(eSet)@other$type <- "combData"
    
    ## annotation
    annotation(eSet) <- annotation(object01)
    
    if (runtype=="bm") {
        cat(paste(Sys.time(), "|   Saving normalized data \n"))
        combData <- eSet
        save(combData, file="combData.RData")
    }
    
    return(eSet)
}