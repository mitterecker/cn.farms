#' Normalization Quantiles
#' @param filenames filenames
#' @param cores cores
#' @param batch batch
#' @param annotDir annotDir
#' @param runtype Mode how the results are saved. Possible values are ff or bm. 
#' If ff is chosen the data will not be saved automatically. 
#' With bm the results will be saved permanently. 
#' @param pkgname Optional parameter for the CEL mapping.
#' @param saveFile Name of the file to save.
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @return The normalized data.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
normalizeQuantiles <- function (
        filenames, 
        cores = 1, 
        batch = NULL, 
        annotDir = NULL, 
        runtype = "ff", 
        pkgname = NULL, 
        saveFile = "normDataQuant") { 
    
    ## assure correct file extension
    saveFile <- gsub("\\.RData", "", saveFile)
    saveFile <- gsub("\\.rda", "", saveFile)
    saveFile <- paste(saveFile, ".RData", sep="")
    
    if(is.null(batch)) {
        batch <- rep(1, length(filenames))
    } 
    
    if (length(filenames) != length(batch)) {
        stop("Filenames and batch need to have the same dimension")
    }
    
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
        pmfeature <- data.frame()
        idxOfAlleleA <- vector()
        idxOfAlleleB <- vector()
        load(normFile)    
    } else {
        stop("Something went wrong with the annotation")
    }    
    
    nbrOfSamples <- length(filenames)
    nbrOfProbes <- length(which(pmfeature[, "allele"] == 1))
    
    intensity <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, 
            type = "double", bmName = gsub("\\.RData", "", saveFile))
    
       
    if (cores == 1) {
        tmpExprs <- affxparser::readCelIntensities(filenames, 
                indices = pmfeature$fid)
        tmpExprs <- (tmpExprs[idxOfAlleleA, ] + tmpExprs[idxOfAlleleB, ]) / 2
        tmpExprs <- log(tmpExprs)
        x <- preprocessCore::normalize.quantiles(tmpExprs)
        intensity[] <- x[]
        rm(x, tmpExprs)
    } else {
        
        ## determine the best array 
        baseline <- determineBaselineArray(
                filenames, nbrOfProbes = 10000, runtype = "ff", pmfeature$fid, 
                cores = cores)
        target <- as.vector(
                affxparser::readCelIntensities(filenames[baseline], 
                        indices = pmfeature$fid))
        target <- target[idxOfAlleleA] + target[idxOfAlleleB]
        
        fidTmp <- pmfeature$fid
        
        gc()
        sfInit(parallel = TRUE, cpus = cores, type = "SOCK")   
        sfLibrary("cn.farms", character.only = TRUE, verbose = FALSE)
        sfLibrary("affxparser", character.only = TRUE, verbose = FALSE)
        sfLibrary("oligo", character.only = TRUE, verbose = FALSE)
        sfLibrary("preprocessCore", character.only = TRUE, verbose = FALSE)
        exportList <- c("fidTmp", "target", "idxOfAlleleA", 
                "idxOfAlleleB", "intensity")
        suppressWarnings(sfExport(list = exportList))
        
        cat(paste(Sys.time(), "|   Starting normalization \n"))
        res <- suppressWarnings(sfLapply(1:nbrOfSamples, normalizeQuantilesH01, 
                        filenames))
    }
    
    sfStop()    
    cat(paste(Sys.time(), "|   Normalization done \n"))
    
    
    eSet <- new("ExpressionSet")
    
    ## assay data    
    assayData(eSet) <- list(intensity = intensity)
    
    ## pheno data
    phenoData(eSet) <- new("AnnotatedDataFrame", data = data.frame(filenames,
                    batch=batch))
    
    ## feature data
    featureData(eSet) <- new("AnnotatedDataFrame", 
            data = pmfeature[pmfeature$allele == 1, c("fid", "fsetid")])
    
    ## experiment data
    experimentData(eSet) <- new("MIAME", 
            other=list(
                    annotDir = annotDir, 
                    normalization = "quantiles", 
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
normalizeQuantilesH01 <- function (ii, filenames) {
    
    ## non-visible bindings
    fidTmp <- fidTmp
    idxOfAlleleA <- idxOfAlleleA
    idxOfAlleleB <- idxOfAlleleB
    target <- target
    
    
    tmpExprs <- affxparser::readCelIntensities(filenames[ii], indices = fidTmp)
    tmpExprs <- (tmpExprs[idxOfAlleleA] + tmpExprs[idxOfAlleleB]) / 2
#    intensity[, ii] <- log2(preprocessCore::normalize.quantiles.use.target(
#                    as.matrix(tmpExprs), target))
    intensity[, ii] <- log(target[rank(tmpExprs)])
}



