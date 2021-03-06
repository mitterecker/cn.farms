#' Runs the SOR normalization on microarray data
#' @param filenames an absolute path of the CEL files 
#' @param cores cores
#' @param annotDir annotDir
#' @param alleles alleles
#' @param cyc states the number of cycles for the EM algorithm.
#' @param runtype Mode how the results are saved. Possible values are ff or bm. 
#' If ff is chosen the data will not be saved automatically. 
#' With bm the results will be saved permanently. 
#' @param pkgname Optional parameter for the CEL mapping.
#' @param saveFile Name of the file to save.
#' @return An instance of \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @importMethodsFrom oligoClasses db
#' @importMethodsFrom DBI dbGetQuery
#' @importFrom affxparser readCelHeader
#' @importFrom oligo cleanPlatformName
normalizeNone <- function (filenames, cores = 1, annotDir = NULL, alleles = FALSE, 
        runtype = "ff", cyc = 5, pkgname = NULL, saveFile = "Sor") {
    
    ## assure correct file extension
    saveFile <- gsub("\\.RData", "", saveFile)
    saveFile <- gsub("\\.rda", "", saveFile)
    saveFile <- paste(saveFile, ".RData", sep = "")
    
    if(is.null(pkgname)) {
        mapping <- affxparser::readCelHeader(filenames[1])$chiptype
        pkgname <- oligo::cleanPlatformName(mapping)
    }
    
    if (is.null(annotDir)) {
        vers <- dir(file.path("annotation", pkgname))[1]
        annotDir <- normalizePath(file.path("annotation", pkgname, vers))
    }
    
    cat(paste(Sys.time(), "|   Annotation directory: ", 
                    annotDir, " \n"))
    
    normFile <- file.path(annotDir, "annotNormalization.RData")
    if(file.exists(normFile)) {
        pmfeature <- data.frame
        load(normFile)    
    } else {
        stop("Something went wrong with the annotation")
    }    
    
    if (cores == 1) {
        sfInit(parallel = FALSE)
    } else {
        sfInit(parallel = TRUE, cpus = cores, type = "SOCK")        
    }
    
    nbrOfSamples <- length(filenames)
    nbrOfProbes <- length(which(pmfeature[, "allele"] == 1))
    
    intensity <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, 
            type = "double", bmName = gsub("\\.rda", "", saveFile))
    
    if (alleles) {
        intensityA <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, 
                type = "double", bmName = gsub("\\.rda", "", saveFile))
        intensityB <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, 
                type = "double", bmName = gsub("\\.rda", "", saveFile))        
    }
    
    cnLibrary("cn.farms", character.only = TRUE)
    cnLibrary("affxparser", character.only = TRUE)
    cnLibrary("oligo", character.only = TRUE)
    
    suppressWarnings(
            sfExport(list = c(
                            "pmfeature", "uniquePairs", 
                            "idxOfAlleleA", "idxOfStrandA",
                            "idxOfStrandB", "idxOfAlleleB", 
                            "pairs", "alleles", "intensity")))
    
    gc()
    sfClusterEval(open(intensity))
    if (alleles) {
        suppressWarnings(sfExport("intensityA"))
        sfClusterEval(open(intensityA))
        suppressWarnings(sfExport("intensityB"))
        sfClusterEval(open(intensityB))
    }
    cat(paste(Sys.time(), "|   Starting normalization \n"))
    res <- sfLapply(1:nbrOfSamples, normalizeNoneH01, filenames, cyc)
    cat(paste(Sys.time(), "|   Normalization done \n"))
    sfStop()
    
    eSet <- new("ExpressionSet")
    
    ## assay data    
    if (alleles) {
        assayData(eSet) <- list(
                intensity  = intensity, 
                intensityA = intensityA, 
                intensityB = intensityB)
    } else {
        assayData(eSet) <- list(intensity = intensity)
    }
    
    ## pheno data
    phenoData(eSet) <- new("AnnotatedDataFrame", data = data.frame(filenames))
    
    ## feature data
    featureData(eSet) <- new("AnnotatedDataFrame", 
            data = pmfeature[pmfeature$allele == 1, c("fid", "fsetid")])
    
    ## experiment data
    experimentData(eSet) <- new("MIAME", 
            other = list(
                    annotDir = annotDir, 
                    normalization = "SOR", 
                    type = "normData"))    
    
    ## annotation
    annotation(eSet) <- pkgname
    
    ## sample names
    sampleNames(eSet) <- gsub("\\.[Cc][Ee][Ll]", "", basename(filenames))    
    
    return(eSet)
}


#' Helper function
#' @param ii ii
#' @param filenames filenames
#' @return Some data
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @noRd
normalizeNoneH01 <- function (ii, filenames, cyc) {
    ## non-visible bindings
    pmfeature    <- pmfeature
    uniquePairs  <- uniquePairs
    idxOfStrandA <- idxOfStrandA
    idxOfStrandB <- idxOfStrandB
    idxOfAlleleA <- idxOfAlleleA
    idxOfAlleleB <- idxOfAlleleB
    alleles <- alleles
    
    tmpExprs <- affxparser::readCelIntensities(filenames[ii], 
            indices = pmfeature$fid)
    alleleA <- tmpExprs[idxOfAlleleA]
    alleleB <- tmpExprs[idxOfAlleleB]
    
    if (alleles) {
        intensityA[, ii] <- alleleA
        intensityB[, ii] <- alleleB
    }
    
    intensity[, ii] <- log2(alleleA + alleleB)
    
    gc()
}