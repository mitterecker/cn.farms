#' Selects the best array (index of filenames) for some kind of baseline 
#' normalization 
#' @param filenames 
#' @param nbrOfProbes 
#' @param runtype 
#' @param pmfeatureFid 
#' @param cores 
#' @param saveFile 
#' @return Some data
#' @author Andreas Mitterecke
#' @noRd
determineBaselineArray <- function(
        filenames, 
        nbrOfProbes = 10000, 
        runtype = "ff", 
        pmfeatureFid, 
        cores = 1, 
        saveFile = "tmp") {

    nbrOfSamples <- length(filenames)
    intensityTmp <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, 
            type = "double", bmName = gsub("\\.RData", "", saveFile))

    sampleFidIdx <- sample(pmfeatureFid, nbrOfProbes)

    if (cores == 1) {
        sfInit(parallel = FALSE)
    } else {
        sfInit(parallel = TRUE, cpus = cores, type = "SOCK")        
    }
    
    sfLibrary("ff", character.only = TRUE, verbose = FALSE, keep.source = FALSE)
    suppressWarnings(sfExport("readCelIntensities", namespace = "affxparser"))
    suppressWarnings(sfExport("sampleFidIdx"))
    suppressWarnings(sfExport("intensityTmp"))
    suppressWarnings(sfExport("filenames"))

    res <- suppressWarnings(
            sfLapply(1:nbrOfSamples, determineBaselineArrayH01)) 
    tmp <- apply(intensityTmp[], 2, median)
    sfStop()
    return(which(rank(tmp) == round(nbrOfSamples / 2, 0)))    
}
    
#' Helper function
#' @param ii ii
#' @param filenames filenames
#' @return Some data
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @noRd
determineBaselineArrayH01 <- function (i) {
    ## non-visible bindings
    filenames <- filenames
    sampleFidIdx <- sampleFidIdx
    
    tmpExprs <- affxparser::readCelIntensities(filenames[i], 
            indices = sampleFidIdx)
    intensityTmp[, i] <- log(tmpExprs)
}
