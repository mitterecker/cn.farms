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
#' @param returnValues List with return values. 
#' For possible values see summaryMethod.
#' @return Results of FARMS run with specified parameters - exact FARMS version
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
callSummarize <- function(object, psInfo, summaryMethod, summaryParam,
        batchList=NULL, cores=1, runtype="ff", returnValues){

    cat(paste(Sys.time(), " |   Starting summarization \n", sep = ""))
    
    if (is.null(batchList)) {
        batchList <- rep(1, ncol(object))
    }
    
    cat(paste(Sys.time(), " |   Computations will take some time,", 
                    " please be patient \n", sep = ""))
    
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
            c("slot11", "lapla",        nbrOfSamples, nbrOfProbes),
            c("slot12", "INICall",      1,            nbrOfProbes),
            c("slot13", "probesize",    1,            nbrOfProbes), 
            c("slot14", "INI_sigVar",   1,            nbrOfProbes),
            c("slot15", "ICtransform",  nbrOfSamples, nbrOfProbes))
    varNames <- data.frame(varNames, stringsAsFactors=FALSE)
    varNames[, 3] <- as.numeric(varNames[, 3])
    varNames[, 4] <- as.numeric(varNames[, 4])
    
    if (missing(returnValues)) {
        if (summaryMethod == "summarizeFarmsVariational") {
            idxNames <- c(1, 5, 7, 11)   # 12 
        } else if (summaryMethod == "summarizeFarmsExact") {
            idxNames <- c(1, 2, 5, 7, 15)    
        } else if (summaryMethod == "summarizeFarmsGaussian"){
            idxNames <- c(1, 5, 12) # 14
        } else {
            idxNames <- c(1, 5)
        }
    } else {
        idxNames <- which(varNames[, 2] %in% returnValues)
    }
    
    ## create objects
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
    
    ## create assignment calls
    resCommands <- vector(length=length(idxNames))
    for (i in seq(length(idxNames))) {
        if (varNames[idxNames[i], 3] != 1) {
            resCommands[i] <- paste(varNames[idxNames[i], 1], 
                    "[i, sampleIndices] <- exprs$", 
                    varNames[idxNames[i], 2], sep="")
        } else {
            resCommands[i] <- paste(varNames[idxNames[i], 1], 
                    "[i] <- exprs$", 
                    varNames[idxNames[i], 2], sep="")
        }
    }
    
    if (cores == 1) {
        sfInit(parallel=FALSE)
    } else {
        sfInit(parallel=TRUE, cpus=cores, type="SOCK")        
    }
    
    sfLibrary("cn.farms", character.only=TRUE)
    #suppressWarnings(sfExport("summarizeFarmsGaussian", namespace="cn.farms"))
    #suppressWarnings(sfExport("summarizeFarmsVariational", namespace="cn.farms"))
    #suppressWarnings(sfExport("summarizeFarmsExact", namespace="cn.farms"))
    
    suppressWarnings(sfExport("psInfo"))
    suppressWarnings(sfExport("object"))
    suppressWarnings(sfExport("summaryParam"))
    suppressWarnings(sfExport("resCommands"))
    for (i in idxNames) {
        tmp <- paste("sfExport(\"", varNames[i, 1],"\")", sep="")
        suppressWarnings(eval(parse(text=tmp)))
    }
    
    batches <- as.character(sort(unique(batchList)))
    for(i in 1:length(batches)){
        cat(paste(Sys.time(), " |   Summarizing batch ", 
                        batches[i], " ... \n", sep = ""))
        sampleIndices <- which(batchList == batches[i])
        res <- sfLapply(seq(nbrOfProbes), callSummarizeH01, sampleIndices, 
                summaryMethod)    
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
callSummarizeH01 <- function(i, sampleIndices, summaryMethod) {
    
    ## non-visible bindings
    psInfo <- psInfo
    object <- object
    summaryParam <- summaryParam
    resCommands <- resCommands
    
    tmpIndices <- psInfo$start[i]:psInfo$end[i]
    pps <- object[tmpIndices, sampleIndices]
    
    ## special case if there is no variation in the data
    tmp <- apply(pps, 1, var) > 1e-20
    if (!all(tmp)) {
        for (m in which(!tmp)) {
            pps[m, ] <- pps[m, ] + rnorm(length(sampleIndices), 0, 0.00005)
        }
    }
    
    exprs <- do.call(summaryMethod, c(alist(pps), summaryParam))
    sapply(resCommands, function (x)  eval(parse(text=x)))
    
    invisible()
}
