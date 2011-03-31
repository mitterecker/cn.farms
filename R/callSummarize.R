#' Defines which variables should be written back
#' @param object Normalized intensity values
#' @param psInfo Physical position
#' @param summaryMethod summaryMethod
#' @param summaryParam summaryParam
#' @param batchList batchList
#' @param cores cores
#' @param runtype Mode how the results are saved. Possible values are ff or bm. 
#' If ff is chosen the data will not be saved automatically. 
#' With bm the results will be saved permanently. 
#' @return Results of FARMS run with specified parameters - exact FARMS version
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
callSummarize <- function(object, psInfo, summaryMethod, summaryParam,
        batchList=NULL, cores=1, runtype="ff", returnValues){
    cat(paste(Sys.time(), " |   Starting summarization \n", sep = ""))
    if (is.null(batchList)) {
        batchList <- rep(1, ncol(object))
    }

    cat(paste(Sys.time(), " |   Computations will take some time, please be patient \n", sep = ""))
    if(runtype=="ram") cores=1
    nbrOfSamples <- ncol(object)
    maxNbrOfProbes <- max((psInfo$end - psInfo$start) + 1)
    nbrOfProbes <- length(psInfo$start)
    
    varNames <- rbind(
            c("slot01", "intensity",    nbrOfSamples, nbrOfProbes),
            c("slot02", "maxZ",         nbrOfSamples, nbrOfProbes),
            c("slot03", "ICtransform",  nbrOfSamples, nbrOfProbes),
            c("slot04", "KL",           nbrOfSamples, nbrOfProbes),
            c("slot05", "L_z",          nbrOfSamples, nbrOfProbes),
            c("slot06", "SNR",          nbrOfSamples, nbrOfProbes),
            c("slot07", "IC",           nbrOfSamples, nbrOfProbes),
            c("slot08", "cns",          nbrOfSamples, nbrOfProbes),
            c("slot09", "medianSample", nbrOfSamples, nbrOfProbes),
            c("slot10", "z",            nbrOfSamples, nbrOfProbes),
            c("slot11", "INICall",      1,            nbrOfProbes),
            c("slot12", "probesize",    1,            nbrOfProbes))
    
    
    if (missing(returnValues)) {
        if (summaryMethod == "summarizeFarmsVariational") {
            idxNames <- c(1, 5, 6)    
        } else if (summaryMethod == "summarizeFarmsExact") {
            idxNames <- c(1, 5, 6)    
        } else if (summaryMethod == "summarizeFarmsGaussian"){
            idxNames <- c(1, 5, 11)
        } else {
            stop("Problem with summaryMethod")
        }
    } else {
        idxNames <- which(varNames[, 2] %in% returnValues)
    }
    
    for (i in idxNames) {
        tmp <- paste(
                varNames[i, 1], 
                " <- createMatrix(\'", 
                runtype, "\'",  
                ", nrow=", varNames[i, 4],
                ", ncol=", varNames[i, 3],
                ", bmName=\'summData_\'",
                ")", 
                sep="")   
        eval(parse(text=tmp))
    }
    
    if (cores == 1) {
        sfInit(parallel=FALSE)
    } else {
        sfInit(parallel=TRUE, cpus=cores, type="SOCK")        
    }
    snp_index <- 1:nbrOfProbes
    
    sfLibrary("cn.farms", character.only=TRUE)
    sfLibrary("ff", character.only=TRUE)
    suppressWarnings(sfExport("snp_index"))
    suppressWarnings(sfExport("psInfo"))
    suppressWarnings(sfExport("object"))
    suppressWarnings(sfExport("summaryMethod"))
    suppressWarnings(sfExport("summaryParam"))
    suppressWarnings(sfExport("idxNames"))
    suppressWarnings(sfExport("varNames"))
    for (i in idxNames) {
        tmp <- paste("sfExport(\"", varNames[i, 1],"\")", sep="")
        suppressWarnings(eval(parse(text=tmp)))
    }
    
    batches <- as.character(sort(unique(batchList)))
    for(i in 1:length(batches)){
        cat(paste(Sys.time(), " |   Summarizing ", 
                        batches[i], " ... \n", sep = ""))
        sampleIndices <- which(batchList == batches[i])
        res <- sfLapply(snp_index, farmsExactH01, sampleIndices)    
    }
    
    sfStop()
    result <- list()
    tmp <- paste("result <- list(", 
            paste(paste(varNames[idxNames, 2], varNames[idxNames, 1], sep="="), 
                    collapse=", "), ")")
    eval(parse(text=tmp))    
    cat(paste(Sys.time(), " |   Summarization done \n", sep = ""))
    return(result)
}


#' Helper function
#' @param i i
#' @param sampleIndices sampleIndices
#' @return Data
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
farmsExactH01 <- function(i, sampleIndices) {
    
    ## non-visible bindings
    snp_index <- snp_index
    psInfo <- psInfo
    object <- object
    summaryMethod <- summaryMethod
    summaryParam <- summaryParam
    idxNames <- idxNames
    varNames <- varNames
    
    
    mc_index <- which(snp_index==i)
    for(k in 1:length(mc_index)){
        tmpIdx <- mc_index[k]
        probeIndices_Start <- psInfo$start[tmpIdx]
        probeIndices_End <- psInfo$end[tmpIdx]
        tmpIndices <- probeIndices_Start:probeIndices_End
        pps <- object[tmpIndices, sampleIndices]
        tmp <- apply(pps, 1, var) > 1e-20
        if (!all(tmp)) {
            for (m in which(!tmp)) {
                pps[m, ] <- pps[m, ] + rnorm(length(sampleIndices), 0, 0.00005)
            }
        }
        
        probesize <- length(tmpIndices)
        exprs <- do.call(summaryMethod, c(alist(pps), summaryParam))
        for (i in idxNames) {
            if (as.numeric(varNames[i, 3]) != 1) {
                tmp <- paste(varNames[i, 1], 
                        "[tmpIdx, sampleIndices] <- exprs$", 
                        varNames[i, 2], sep="")
            } else {
                tmp <- paste(varNames[i, 1], 
                        "[tmpIdx] <- exprs$", 
                        varNames[i, 2], sep="")
            }
            eval(parse(text=tmp))
        }
    }
    
    invisible()
}
