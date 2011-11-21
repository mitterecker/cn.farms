#' Does a fragment length correction
#' @param object An instance of 
#' \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' @param runtype Mode how the results are saved. Possible values are ff or bm. 
#' @param ... Further parameters passed to the correction method.
#' @param saveFile Name of the file to save.
#' @return An instance of 
#' \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @examples
#' load(system.file("exampleData/slData.RData", package = "cn.farms"))
#' slDataFlc <- fragLengCorr(slData)
#' @export
fragLengCorr <- function (
        object, 
        runtype = "ff", 
        saveFile = "slDataFlc", 
        ...) {
    
    ## assure correct file extension
    saveFile <- gsub("\\.RData", "", saveFile)
    saveFile <- gsub("\\.rda", "", saveFile)
    saveFile <- paste(saveFile, ".RData", sep = "")
    
    normAdd <- normAdd(object@annotation)
    if (normAdd %in% c("Nsp", "Sty", "Hind240", "Xba240")) {
        saveFile <- paste(gsub("\\.RData", "", saveFile), 
                normAdd, ".RData", sep = "")
    } 
    
    if (runtype == "bm" & file.exists(saveFile)) {
        message("Fragment length normalization has already been done")
        message("Trying to load  data ...")
        load(saveFile)
        return(slData)
    }
    
    cat(paste(Sys.time(), " |   Starting fragment length correction \n", 
                    sep = ""))
    
    if (length(experimentData(object)@other$flc) != 0) {
        
        stop("Fragment length correction already done")
        
    } else if (annotation(object) %in% 
            c("pd.genomewidesnp.5", "pd.genomewidesnp.6")) {
        
        fl <- featureData(object)@data[, 
                c("fragment_length", "fragment_length2")]
        x <- assayData(object)$intensity
        flcD <- flcSnp6Std(x, fl, saveFile = saveFile, runtype = runtype, ...)
        
        slData <- new("ExpressionSet")
        assayData(slData) <- list(intensity = flcD)
        phenoData(slData) <- phenoData(object)
        featureData(slData) <- featureData(object)
        experimentData(slData) <- experimentData(object)
        annotation(slData) <- annotation(object)
        sampleNames(slData) <- sampleNames(object)  
        experimentData(slData)@other$flc <- 1
        cat(paste(Sys.time(), " |   Fragment length correction done \n", 
                        sep = ""))
        
    } else if (annotation(object) != "pd.genomewidesnp.6") {
        
        fl <- featureData(object)@data[, c("fragment_length")]
        x <- assayData(object)$intensity
        flcD <- flcStd(x[], fl, saveFile = saveFile, runtype = runtype, ...)
               
        slData <- new("ExpressionSet")
        assayData(slData) <- list(intensity = flcD)
        phenoData(slData) <- phenoData(object)
        featureData(slData) <- featureData(object)
        experimentData(slData) <- experimentData(object)
        annotation(slData) <- annotation(object)
        sampleNames(slData) <- sampleNames(object)  
        experimentData(slData)@other$flc <- 1
        cat(paste(Sys.time(), " |   Fragment length correction done \n", 
                        sep = ""))
        
    } else {
        stop("We have a problem")
    }
    
    ##FIXME: implement also properly for npData
    if (runtype == "bm") {
        cat(paste(Sys.time(), "|   Saving normalized data \n"))
        save(slData, file = saveFile)
    }
    
    return(slData)
}

#' Does a fragment length correction on intensities 
#' @param y y
#' @param fragmentLengths fragmentLengths
#' @param targetFcn targetFcn
#' @param subsetToFit subsetToFit
#' @param runtype Mode how the results are saved. Possible values are ff or bm. 
#' If ff is chosen the data will not be saved automatically. 
#' With bm the results will be saved permanently. 
#' @param cores cores
#' @param saveFile Name of the file to save.
#' @param ... ...
#' @return data frame
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
flcStd <- function(y, fragmentLengths, targetFcn = NULL, subsetToFit = NULL,
        runtype = "ff", cores = 1, saveFile = "flc", ...) {
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Estimate normalization function
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit smooth curve
    
    
    ##FIXME: flc has a problem with bm!
    
    nbrOfSamples <- ncol(y)
    nbrOfProbes <- nrow(y)
    yN <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, type = "double", 
            bmName = gsub("\\.RData", "", saveFile))
    ok <- (is.finite(fragmentLengths) & is.finite(y[,i]))
    
    if (cores == 1) {
        sfInit(parallel = FALSE)
    } else {
        sfInit(parallel = TRUE, cpus = cores, type = "SOCK")        
    }
    
    sfLibrary("stats", character.only = TRUE, verbose = FALSE, keep.source = FALSE)
    sfLibrary("ff", character.only = TRUE, verbose = FALSE, keep.source = FALSE)
    
    suppressWarnings(sfExport(list = c(
                            "nbrOfSamples", "nbrOfProbes", 
                            "yN", "ok", "y", "subsetToFit", 
                            "fragmentLengths")))
    
    res <- suppressWarnings(sfLapply(1:nbrOfSamples, flcStdH01))
    
    sfStop()
    
    fit <- lowess(fragmentLengths[ok], Biobase::rowMedians(y[ok,]))
    
    yN_average <- approx(fit, xout = fragmentLengths, ties = mean)$y
    
    y_hat <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, 
            type = "double", bmName = gsub("\\.RData", "", saveFile))
    
    for (i in 1:nbrOfSamples) {
        
        y_hat[, i] <- y[, i] - (yN[, i] - yN_average)
        
    }
    
    return(y_hat)
}


#' Helper function
#' @param i i
#' @param ... ...
#' @return Some data
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @noRd
flcStdH01 <- function (i, ...) {
    
    ## non-visible bindings
    subsetToFit <- subsetToFit
    fragmentLengths <- fragmentLengths
    y <- y
    
    if (!is.null(subsetToFit)) {
        ok[-subsetToFit] <- FALSE
    }
    
    suppressWarnings({
                
                fit <- lowess(fragmentLengths[ok], y[ok,i], ...)
                
            })
    
    yN[,i] <- approx(fit, xout = fragmentLengths, ties = mean)$y
    
}






#' Does a fragment length correction on intensities 
#' @param y y
#' @param fragmentLengths fragmentLengths
#' @param targetFcn targetFcn
#' @param subsetToFit subsetToFit
#' @param runtype runtype
#' @param cores cores
#' @param saveFile Name of the file to save.
#' @param ... ...
#' @return data frame 
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
flcSnp6Std <- function(y, fragmentLengths, targetFcn = NULL, 
        subsetToFit = NULL, runtype = "ff", cores = 1, saveFile = "flc", ...) {
    
    ## adapted from the aroma.affymetrix package (www.aroma.project.org)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Estimate normalization function
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit smooth curve
    nbrOfSamples <- ncol(y)
    nbrOfProbes <- nrow(y)
    yN <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, type = "double", 
            bmName = gsub("\\.RData", "", saveFile))
    
    ## FIXME: hack if fragment length is NA
    idxTmp <- is.na(fragmentLengths[, 1])
    if (length(idxTmp) != 0) {
        fragmentLengths[idxTmp, 1] <- 1200
    }
    
    idxTmp <- is.na(fragmentLengths[, 2])
    if (length(idxTmp) != 0) {
        fragmentLengths[idxTmp, 2] <- 1200
    }
    
    fl_1_ok <- is.finite(fragmentLengths[, 1]) 
    fl_2_ok <- is.finite(fragmentLengths[, 2]) 
    
    if (cores == 1) {
        sfInit(parallel = FALSE)
    } else {
        sfInit(parallel = TRUE, cpus = cores, type = "SOCK")        
    }
    
    sfLibrary("stats", character.only = TRUE, verbose = FALSE, keep.source = FALSE)
    sfLibrary("ff", character.only = TRUE, verbose = FALSE, keep.source = FALSE)
    
    suppressWarnings(sfExport(list = c(
                            "nbrOfSamples", "nbrOfProbes", 
                            "yN", "fl_1_ok", "fl_2_ok", "y", 
                            "subsetToFit", "fragmentLengths")))
    
    res <- suppressWarnings(sfLapply(1:nbrOfSamples, flcSnp6StdH01))
    sfStop()
    
    y_median <- ffapply(X = y, MARGIN = 1, AFUN = "median", RETURN = TRUE)
    yN_average <- rep(0, nbrOfProbes)
    y_ok <- is.finite(y_median[])
    ok_1 <- ok <- (fl_1_ok & y_ok)
    if (!is.null(subsetToFit)) {
        ok[-subsetToFit] <- FALSE
    }
    
    fit <- lowess(fragmentLengths[ok,1], y_median[ok])
    yN_average[ok_1] <- approx(fit, xout = fragmentLengths[ok_1, 1], ties = mean)$y
    ok_2 <- ok <- (fl_2_ok & y_ok)
    
    if (!is.null(subsetToFit)) {
        ok[-subsetToFit] <- FALSE
    }
    
    fit <- lowess(fragmentLengths[ok, 2], y_median[ok])
    yN_average[ok_2] <- yN_average[ok_2] + 
            approx(fit, xout = fragmentLengths[ok_2,2], ties = mean)$y
    yN_average[ok_1 & ok_2] <- yN_average[ok_1 & ok_2] / 2
    
    y_hat <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, type = "double", 
            bmName = gsub("\\.RData", "", saveFile))
    
    for (i in 1:nbrOfSamples) {
        y_hat[, i] <- y[, i] - (yN[, i] - yN_average)
    }
    
    return(y_hat)
}


#' Helper function
#' @param i i
#' @param ... ...
#' @return Data
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @noRd
flcSnp6StdH01 <- function (i, ...) {
    
    ## non-visible bindings
    y <- y
    fl_1_ok <- fl_1_ok
    subsetToFit <- subsetToFit
    fragmentLengths <- fragmentLengths
    fl_2_ok <- fl_2_ok
    
    
    ## adapted from the aroma.affymetrix package (www.aroma.project.org)
    
    y_ok <- is.finite(y[, i])
    ok_1 <- ok <- (fl_1_ok & y_ok)
    if (!is.null(subsetToFit)) {
        ok[-subsetToFit] <- FALSE
    }
    
    ## fit lowess 
    suppressWarnings({
                fit <- lowess(fragmentLengths[ok, 1], y[ok, i], ...)
            })
    
    ## predict fl effect
    yN[ok_1, i] <- approx(fit, xout = fragmentLengths[ok_1, 1], ties = mean)$y
    ok_2 <- ok <- (fl_2_ok & y_ok)
    
    if (!is.null(subsetToFit)) {
        ok[-subsetToFit] <- FALSE
    }
    
    ## fit lowess 
    suppressWarnings({
                fit <- lowess(fragmentLengths[ok, 2], y[ok, i], ...)
            })
    
    ## predict fl effect
    yN[ok_2, i] <- yN[ok_2, i] + 
            approx(fit, xout = fragmentLengths[ok_2, 2], ties = mean)$y
    yN[(ok_1 & ok_2), i] <- yN[(ok_1 & ok_2), i] / 2
}
