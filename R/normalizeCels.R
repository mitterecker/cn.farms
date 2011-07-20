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
#' @param saveFile Name of the file to save.
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
        cores = 1, 
        alleles = FALSE, 
        runtype = "bm", 
        annotDir = NULL, 
        saveFile = "normData", 
        ...) {

    ## assure correct file extension
    saveFile <- gsub("\\.RData", "", saveFile)
    saveFile <- gsub("\\.rda", "", saveFile)
    saveFile <- paste(saveFile, ".RData", sep="")    
    
    mapping <- affxparser::readCelHeader(filenames[1])$chiptype
    pkgname <- oligo::cleanPlatformName(mapping)
    pkgname <- correctPkgname(pkgname)
    
    normAdd <- normAdd(pkgname)
    
    if (normAdd %in% c("Nsp", "Sty", "Hind240", "Xba240")) {
        saveFile <- paste(gsub("\\.RData", "", saveFile), 
                normAdd, ".RData", sep="")
    } 
    
    method <- match.arg(method)
    normMethods <- c("SOR", "quantiles")
       
    if (!method %in% normMethods) {
        stop("Normalization method not found!")
    } else if (runtype == "bm" & file.exists(saveFile)) {
        message("Normalization has already been done")
        message("Trying to load normalized data ...")
        load(saveFile)
        return(normData)
    }
    normData <- switch(method, 
            SOR = normalizeSor(filenames = filenames, cores = cores, 
                    alleles = alleles, runtype = runtype, annotDir = annotDir, 
                    pkgname = pkgname, saveFile = saveFile, ...), 
            quantiles = normalizeQuantiles(filenames = filenames, cores = cores, 
                    runtype = runtype, annotDir = annotDir, pkgname = pkgname, 
                    saveFile = saveFile, ...))
    
    if (runtype == "bm") {
        cat(paste(Sys.time(), "|   Saving normalized data \n"))
        save(normData, file=saveFile)
    }
    return(normData)
}
