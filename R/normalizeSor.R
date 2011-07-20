#' Runs the SOR normalization on microarray data
#' @param filenames an absolute path of the CEL files 
#' @param cores cores
#' @param annotDir annotDir
#' @param alleles alleles
#' @param runtype Mode how the results are saved. Possible values are ff or bm. 
#' If ff is chosen the data will not be saved automatically. 
#' With bm the results will be saved permanently. 
#' @param pkgname Optional parameter for the CEL mapping.
#' @param saveFile Name of the file to save.
#' @return An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @importMethodsFrom oligoClasses db
#' @importMethodsFrom DBI dbGetQuery
#' @importFrom affxparser readCelHeader
#' @importFrom oligo cleanPlatformName
#' @importFrom snowfall sfInit
#' @importFrom snowfall sfExport
#' @importFrom snowfall sfLibrary
#' @importFrom snowfall sfLapply
#' @importFrom snowfall sfClusterEval
#' @importFrom snowfall sfStop
normalizeSor <- function (filenames, cores=1, annotDir = NULL, alleles = FALSE, 
        runtype = "ff", cyc = 5, pkgname = NULL, saveFile = "Sor", ...) {
    
    ## assure correct file extension
    saveFile <- gsub("\\.RData", "", saveFile)
    saveFile <- gsub("\\.rda", "", saveFile)
    saveFile <- paste(saveFile, ".RData", sep="")
    
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
        sfInit(parallel=FALSE)
    } else {
        sfInit(parallel=TRUE, cpus=cores, type="SOCK")        
    }
    
    nbrOfSamples <- length(filenames)
    nbrOfProbes <- length(which(pmfeature[, "allele"]==1))
    
    intensity <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, type="double",
            bmName=gsub("\\.rda", "", saveFile))
    
    if (alleles) {
        intensityA <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, 
                type="double", bmName=gsub("\\.rda", "", saveFile))
        intensityB <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, 
                type="double", bmName=gsub("\\.rda", "", saveFile))        
    }

    sfLibrary("cn.farms", character.only=TRUE)
    sfLibrary("affxparser", character.only=TRUE)
    sfLibrary("oligo", character.only=TRUE)
    
    suppressWarnings(
            sfExport(list=c(
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
    res <- sfLapply(1:nbrOfSamples, normalizeSorH01, filenames, cyc)
    cat(paste(Sys.time(), "|   Normalization done \n"))
    sfStop()
    
    eSet <- new("ExpressionSet")
    
    ## assay data    
    if (alleles) {
        assayData(eSet) <- list(
                intensity=intensity, 
                intensityA=intensityA, 
                intensityB=intensityB)
    } else {
        assayData(eSet) <- list(intensity=intensity)
    }
    
    ## pheno data
    phenoData(eSet) <- new("AnnotatedDataFrame", data = data.frame(filenames))
    
    ## feature data
    featureData(eSet) <- new("AnnotatedDataFrame", 
            data = pmfeature[pmfeature$allele==1, c("fid", "fsetid")])
    
    ## experiment data
    experimentData(eSet) <- new("MIAME", 
            other=list(
                    annotDir=annotDir, 
                    normalization="SOR", 
                    type="normData"))    
    
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
normalizeSorH01 <- function (ii, filenames, cyc) {
    
    ## non-visible bindings
    pmfeature <- pmfeature
    uniquePairs <- uniquePairs
    idxOfStrandA <- idxOfStrandA
    idxOfStrandB <- idxOfStrandB
    idxOfAlleleA <- idxOfAlleleA
    idxOfAlleleB <- idxOfAlleleB
    alleles <- alleles
    
    tmpExprs <- affxparser::readCelIntensities(filenames[ii], 
            indices=pmfeature$fid)
    
    LZExprs <- matrix(NA, dim(tmpExprs)[1], dim(tmpExprs)[2])
    for (jj in 1:length(uniquePairs)) {
        for (kk in 1:2) { # strand
            tmp_indices <- which(pairs==uniquePairs[jj])
            if (kk == 1) {
                tmp_indices <- intersect(idxOfStrandA, tmp_indices)    
            } else {
                tmp_indices <- intersect(idxOfStrandB, tmp_indices)    
            }
            
            tmp_exprs_A <- tmpExprs[idxOfAlleleA[tmp_indices]]
            tmp_exprs_B <- tmpExprs[idxOfAlleleB[tmp_indices]]
            foo1 <- matrix(c(tmp_exprs_A, tmp_exprs_B), length(tmp_exprs_A), 2, 
                    byrow=FALSE)
            foo1[, 1] <- foo1[, 1] - mean(sort(foo1[, 1])[1:100])
            foo1[, 2] <- foo1[, 2] - mean(sort(foo1[, 2])[1:100])
            res <- sparseFarmsC(foo1, cyc)
            LzId1 <- t(res$Lz)
            LzId1 <- normalizeAverage(LzId1)
            LZExprs[idxOfAlleleA[tmp_indices]] <- LzId1[, 1]
            LZExprs[idxOfAlleleB[tmp_indices]] <- LzId1[, 2]
            
            tmp2 <- log2(LZExprs + 175)
            if (alleles) {
                intensityA[, ii] <- tmp2[idxOfAlleleA]
                intensityB[, ii] <- tmp2[idxOfAlleleB]
            }
            intensity[, ii] <- (tmp2[idxOfAlleleA] + tmp2[idxOfAlleleB]) / 2
        }
    }
    gc()
}