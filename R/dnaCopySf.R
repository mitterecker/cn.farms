#' Runs DNAcopy in parallel mode
#' 
#' This function even works very well with ff matrices,
#' 
#' @param x A matrix with data of the copy number experiments
#' @param chrom The chromosomes (or other group identifier) from which the markers came
#' @param maploc  The locations of marker on the genome 
#' @param cores Number of cores to use
#' @param smoothing States if smoothing of the data should be done 
#' @param ... Further parameter for the function segment of DNAcopy
#' @return An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' containing the segments.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @importFrom DNAcopy CNA
#' @importFrom DNAcopy smooth.CNA
#' @importFrom DNAcopy segment
#' @export
#' @examples
#' load(system.file("exampleData/mlData.RData", package="cn.farms"))
#' mlData <- mlData[, 1:3]
#' colnames(assayData(mlData)$L_z) <- sampleNames(mlData)
#' segments <- dnaCopySf(
#'         x         = assayData(mlData)$L_z, 
#'         chrom     = featureData(mlData)@@data$chrom, 
#'         maploc    = featureData(mlData)@@data$start, 
#'         cores     = 1, 
#'         smoothing = FALSE)
#' featureData(segments)@@data
dnaCopySf <- function (x, chrom, maploc, cores=1, smoothing, ...) {
    t00 <- Sys.time()
       
    if (is.null(colnames(x))) {
        stop("Colnames of x must not be empty")    
    }
    
    if (nrow(x) != length(chrom)) {
        stop("Rownames of x must have correct dimension")    
    }
    
    if (cores == 1) {
        sfInit(parallel=FALSE)
    } else {
        sfInit(parallel=TRUE, cpus=cores, type="SOCK")        
    }
    
    sfLibrary("ff", character.only=TRUE, verbose=FALSE)
    sfLibrary("DNAcopy", character.only=TRUE, verbose=FALSE)
    suppressWarnings(sfExport("x"))
    suppressWarnings(sfExport("chrom"))
    suppressWarnings(sfExport("maploc"))
    suppressWarnings(sfExport("smoothing"))
    
    res <- suppressWarnings(sfLapply(1:ncol(x), dnaCopySfH01, ...))
    
    ## assign ID
    idSamples <- colnames(x)
    for (i in 1:length(res)) {
        res[[i]]$ID[] <- idSamples[i]
    }
    sfStop()
    
    eSet <- new("ExpressionSet")
    phInf <- do.call("rbind", res)
    phInf <- cbind(phInf[, -1], individual=phInf[, 1])
    colnames(phInf)[1:3] <- c("chrom", "start", "end")
    featureData(eSet) <- new("AnnotatedDataFrame", 
            data = phInf)
    nbrOfSamples <- length(idSamples)
    nbrOfCnvrs <- nrow(featureData(eSet))
    assayData(eSet) <- list(cnv=matrix(rep(2, nbrOfSamples * nbrOfCnvrs), 
                    ncol=nbrOfSamples)) 
    experimentData(eSet)@other$cnvLabels <- list(
            color=c("red", "black", "green"), 
            desc=c("deletion", "normal", "duplication"))
    print(difftime(Sys.time(), t00))
    return(eSet)
}

#' Helper function
#' @param i i
#' @param ... ...
#' @return Some data
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
dnaCopySfH01 <- function (i, ...) {
    if (!exists("min.width")) min.width <- 3
    if (!exists("undo.splits")) undo.splits <- "sdundo"
    if (!exists("undo.SD")) undo.SD <- 1
    
    ## non-visible bindings
    x <- x
    chrom <- chrom
    maploc <- maploc
    smoothing <- smoothing
    
    cnaObj <- DNAcopy::CNA(
            x[, i], 
            chrom, 
            maploc,
            data.type="logratio")

    if(smoothing == T) {
        cnaObj <- DNAcopy::smooth.CNA(cnaObj)
    } 
    segs <- DNAcopy::segment(cnaObj, min.width=min.width, 
            undo.splits=undo.splits, undo.SD=undo.SD, ...)     
    return(segs$output)
}