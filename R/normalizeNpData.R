#' Processes the non-polymorphic data
#' 
#' Normalization for non-polymorphic data for Affymetrix SNP5 and SNP6
#' 
#' @param filenames the absolute path of the CEL files as a list
#' @param cores number of cores for used for parallelization
#' @param annotDir Optional annotation directory.
#' @param runtype Mode how the results are saved. Possible values are ff or bm. 
#' If ff is chosen the data will not be saved automatically. 
#' With bm the results will be saved permanently. 
#' @return An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}} 
#' containing the non-polymorphic data of the microarray.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @examples 
#' \dontrun{
#' library("hapmapsnp6") 
#' celDir <- system.file("celFiles", package="hapmapsnp6")
#' filenames <- dir(path=celDir, full.names=TRUE)
#' createAnnotation(filenames=filenames)
#' npData <- normalizeNpData(filenames)
#' }
#' @export
normalizeNpData <- function(filenames, cores=1, 
        annotDir=NULL, runtype="ff"){
    
    mapping <- affxparser::readCelHeader(filenames[1])$chiptype
    pkgname <- oligo::cleanPlatformName(mapping)
    
    if (pkgname == "pd.genomewideex.6") {
        pkgname <- "pd.genomewidesnp.6"
    }
    
    if (is.null(annotDir)) {
        vers <- dir(file.path("annotation", pkgname))[1]
        annotDir <- normalizePath(file.path("annotation", pkgname, vers))
    }
    
    if (runtype=="bm" & file.exists("npData.RData")) {
        message("Normalization has already been done")
        message("Trying to load normalized data ...")
        load("npData.RData")
        return(npData)
    }
    
    cat(paste(Sys.time(), "|   Annotation directory: ", 
                    annotDir, " \n"))
    
    featureFile <- file.path(annotDir, "pmfeatureCNV.RData")
    
    if(file.exists(featureFile)) {
        pmfeatureCNV <- data.frame()
        featureSetCNV <- data.frame()
        load(featureFile)    
    } else {
        stop("Something went wrong with the annotation")
    }    
    
    nbrOfSamples <- length(filenames)
    nbrOfProbes <- nrow(pmfeatureCNV)
    
    intensity <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, type="double",
            bmName="np_")
    
    if (cores == 1) {
        sfInit(parallel=FALSE)
    } else {
        sfInit(parallel=TRUE, cpus=cores, type="SOCK")        
    }
    
    #sfLibrary("affxparser", character.only=TRUE)
    suppressWarnings(sfExport("normalizeAverage", namespace="cn.farms"))
    suppressWarnings(sfExport("readCelIntensities", namespace="affxparser"))
    suppressWarnings(sfExport("pmfeatureCNV"))
    suppressWarnings(sfExport("intensity"))
    suppressWarnings(sfExport("filenames"))
    
    cat(paste(Sys.time(), "|   Starting processing \n"))
    res <- suppressWarnings(sfLapply(1:nbrOfSamples, normalizeNpDataH01))
    cat(paste(Sys.time(), "|   Non-polymorphic data done \n"))
    
    gc()
    
    sfStop()
    
    featureSetFile <- file.path(annotDir, "featureSetCNV.RData")
    
    if(file.exists(featureSetFile)) {
        load(featureSetFile)    
    } else {
        stop("Something went wrong with the annotation")
    }    
    
    phInf <- featureSetCNV[match(pmfeatureCNV$fsetid, featureSetCNV$fsetid), ]    
    tmpIdx <- which(colnames(phInf) %in% c("chrom_start", "chrom_stop"))
    colnames(phInf)[tmpIdx] <- c("start", "end")
    
    eSet <- new("ExpressionSet")
    assayData(eSet) <- list(intensity=intensity)
    phenoData(eSet) <- new("AnnotatedDataFrame", data = data.frame(filenames))
    featureData(eSet) <- new("AnnotatedDataFrame", data = phInf)
    experimentData(eSet) <- new("MIAME", 
            other=list(
                    annotDir=annotDir, 
                    normalization="SOR",
                    type="npData"))    
    annotation(eSet) <- pkgname
    
    sampleNames(eSet) <- gsub("\\.[Cc][Ee][Ll]", "", basename(filenames))    
    
    if (runtype=="bm") {
        cat(paste(Sys.time(), "|   Saving normalized data \n"))
        npData <- eSet
        save(npData, file="npData.RData")
    }
    
    return(eSet)
}


#' Helper function 
#' @param i i
#' @return data
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @importFrom affxparser readCelIntensities
normalizeNpDataH01 <- function(i){
    
    ## non-visible bindings
    filenames <- filenames
    pmfeatureCNV <- pmfeatureCNV
    
    LZExprs <- affxparser::readCelIntensities(filenames[i], indices=pmfeatureCNV$fid)
    LZExprs <- normalizeAverage(LZExprs)
    LZExprs <- log2(LZExprs)
    intensity[, i] <- LZExprs
}
