#' Wrapper for the normalization functions
#' 
#' This functions provides different normalization methods for microarray data.
#' At the moment only SOR and quantile normalization are implemented.
#' 
#' @param filenames The absolute path of the CEL files as a list.
#' @param method The normalization method. Possible methods so far: 
#' SOR, quantiles
#' @param cores Number of cores for used for parallelization.
#' @param alleles States if information for allele A and B should be given back.
#' @param runtype Mode how the results are saved. Possible values are ff or bm. 
#' If ff is chosen the data will not be saved automatically. 
#' With bm the results will be saved permanently. 
#' @param annotDir An optional annotation directory.
#' @param ... Further parameters for the normalization method. 
#' @return An ExpressionSet object with the normalized data.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @examples 
#' \dontrun{
#' library("hapmapsnp6") 
#' celDir <- system.file("celFiles", package="hapmapsnp6")
#' filenames <- dir(path=celDir, full.names=TRUE)
#' createAnnotation(filenames=filenames)
#' normData <- normalizeCels(filenames, method="SOR")
#' }
normalizeCels <- function (
        filenames, 
        method = c("SOR", "quantiles"), 
        cores=1, 
        alleles=FALSE, 
        runtype="bm", 
        annotDir=NULL, ...) {

    mapping <- affxparser::readCelHeader(filenames[1])$chiptype
    pkgname <- oligo::cleanPlatformName(mapping)
    
    if (pkgname == "pd.genomewideex.6") {
        pkgname <- "pd.genomewidesnp.6"
    }
    
    normAdd <- normAdd(pkgname)
    
    if (normAdd %in% c("Nsp", "Sty", "Hind240", "Xba240")) {
        loadFile <- paste("normData", normAdd, ".RData", sep="")
    } else {
        loadFile <- paste("normData.RData")
    }

    method <- match.arg(method)
    normMethods <- c("SOR", "quantiles")
       
    if (!method %in% normMethods) {
        stop("Normalization method not found!")
    } else if (runtype=="bm" & 
            file.exists(loadFile)) {
        message("Normalization has already been done")
        message("Trying to load normalized data ...")
        load(loadFile)
        return(normData)
    }
    normData <- switch(method, 
            SOR = normalizeSor(filenames=filenames, cores=cores, 
                    alleles=alleles, runtype=runtype, annotDir=annotDir, 
                    pkgname=pkgname, ...), 
            quantiles = normalizeQuantiles(filenames=filenames, cores=cores, 
                    runtype=runtype, annotDir=annotDir, pkgname=pkgname, ...))
    
    if (runtype=="bm") {
        cat(paste(Sys.time(), "|   Saving normalized data \n"))
        save(normData, file=loadFile)
    }
    return(normData)
}
