#' Method for computation of the single-locus summarization
#'
#' The different probes of the SNPs of the array are summarized to a probeset.
#' 
#' @name slSummarization
#' @param object An instance of 
#' \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
#' @param summaryMethod allowed versions for the summarization step are: 
#' Gaussian,Variational, Exact. Default is Variational.
#' @param summaryParam The parameters for the summaryMethod. Further information
#' can be obtained via the according functions: 
#' \code{\link[cn.farms:summarizeFarmsGaussian]{cn.farms}},
#' \code{\link[cn.farms:summarizeFarmsVariational]{cn.farms}} or 
#' \code{\link[cn.farms:summarizeFarmsExact]{cn.farms}}
#' @param callParam Additional parameters for runtype (ff or bm) as well as 
#' cores for parallelization.
#' @param summaryWindow Method for combination of the SNPs. Possible values are 
#' sl and fragment.
#' @param returnValues List with return values. 
#' @param saveFile Name of the file to save.
#' @seealso \code{\link{summarizeFarmsExact}}
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @return Single-locus summarized data of an instance of 
#' \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
#' @examples 
#' load(system.file("exampleData/normData.RData", package = "cn.farms"))
#' notes(experimentData(normData))$annotDir <- 
#'         system.file("exampleData/annotation/pd.genomewidesnp.6/1.1.0",
#'                 package = "cn.farms")
#' summaryMethod <- "Variational"
#' summaryParam <- list()
#' summaryParam$cyc <- c(10)
#' slData <- slSummarization(normData, 
#'         summaryMethod = summaryMethod, 
#'         summaryParam = summaryParam)
#' assayData(slData)$L_z[1:10, 1:10]
#' 
#' summaryMethod <- "Gaussian"
#' summaryParam <- list()
#' summaryParam$cyc <- c(10)
#' slData <- slSummarization(normData, 
#'         summaryMethod = summaryMethod, 
#'         summaryParam = summaryParam)
#' assayData(slData)$L_z[1:10, 1:10]
#' 
#' summaryMethod <- "Exact"
#' summaryParam <- list()
#' summaryParam$cyc <- c(10, 20)
#' slData <- slSummarization(normData, 
#'         summaryMethod = summaryMethod, 
#'         summaryParam = summaryParam)
#' assayData(slData)$L_z[1:10, 1:10]
slSummarization <- function(
        object, 
        summaryMethod = "Exact", 
        summaryParam = list(), 
        callParam = list(runtype = "ff", cores = 1), 
        summaryWindow = c("std", "fragment"), 
        returnValues, 
        saveFile = "slData") {
    
    if ("cores" %in% names(callParam) == FALSE) {
        callParam$cores <- 1
    } 
    if ("runtype" %in% names(callParam) == FALSE) {
        callParam$runtype <- "ff"
    } 
    
    ## assure correct file extension
    saveFile <- gsub("\\.RData", "", saveFile)
    saveFile <- gsub("\\.rda", "", saveFile)
    saveFile <- paste(saveFile, ".RData", sep = "")
    
    normAdd <- normAdd(object@annotation)
    if (normAdd %in% c("Nsp", "Sty", "Hind240", "Xba240")) {
        saveFile <- paste(gsub("\\.RData", "", saveFile), 
                normAdd, ".RData", sep = "")
    }
    
    if (callParam$runtype == "bm" & file.exists(saveFile)) {
        message("Single-locus summarization has already been done.")
        message("Trying to load data ...")
        load(saveFile)
        return(slData)    
    }
    
    summaryWindow <- match.arg(summaryWindow)
    if (summaryWindow == "fragment") {
        featureSet <- data.frame()
        load(file.path(notes(experimentData(object))$annotDir, 
                        "featureSet.RData"))
        tmp <- match(featureData(object)$fsetid, featureSet$fsetid)
        runIdx <- getFragmentSet(featureSet$fragment_length[tmp])
        rm(tmp)
    } else {
        runIdx <- getSingleProbeSetSize(featureData(object)$fsetid)        
    }
    
    t00 <- Sys.time()
    
    if (tolower(summaryMethod) == "median") {
        nbrOfSamples <- ncol(object)
        slRes <- createMatrix(callParam$runtype, nrow(runIdx), nbrOfSamples, 
                type = "double", bmName = gsub("\\.RData", "", saveFile))
        if (callParam$cores == 1) {
            x <- assayData(object)$intensity
            for (i in 1:nrow(runIdx)) {
                slRes[i, ] <- rowMedians(
                        t(x[runIdx[i, 1]:runIdx[i, 2],]))
            }
            rm(x)
        } else {
            sfInit(parallel = TRUE, cpus = callParam$cores, type = "SOCK")
            suppressWarnings(cnLibrary("Biobase", character.only = TRUE))
            suppressWarnings(cnLibrary("ff", character.only = TRUE))
            suppressWarnings(sfExport("slRes"))
            sfClusterEval(open(slRes))
            suppressWarnings(sfExport("runIdx"))
            intensity <- assayData(object)$intensity
            suppressWarnings(sfExport("intensity"))
            
            #splitIdx <- clusterSplit(seq_len(cores), 1:nrow(runIdx))
            dummy <- sfLapply(seq_len(nrow(runIdx)), function(i) {
                        slRes[i, ] <- rowMedians(
                                t(intensity[runIdx[i, 1]:runIdx[i, 2], ]))
                    })
            sfClusterEval(close(slRes))
            sfStop()
        }
        
        myData <- list(intensity = slRes)
    } else {
        summaryMethodName <- paste("summarizeFarms", paste(
                        toupper(substring(summaryMethod, 1,1)), 
                        substring(summaryMethod, 2),
                        sep = "", collapse = " "), sep = "")
        
        if (!exists(summaryMethodName)) {
            stop(paste("Unknown method (can't find function", summaryMethodName, 
                            ")"))
        }
        
        myData <- do.call("callSummarize", c(alist(
                                object = assayData(object)$intensity, 
                                psInfo = runIdx,
                                batchList = phenoData(object)$batch,
                                summaryMethod = summaryMethodName, 
                                summaryParam = summaryParam,
                                returnValues = returnValues, 
                                saveFile = saveFile), 
                        callParam))
    }
    
    
    eSet <- new("ExpressionSet")
    
    ## assay data    
    assayData(eSet) <- myData
    
    ## pheno data
    phenoData(eSet) <- phenoData(object)
    
    ## feature data
    load(file.path(notes(experimentData(object))$annotDir, "featureSet.RData"))
    if (object@annotation == "pd.mapping250k.nsp" | 
            object@annotation == "pd.mapping250k.sty" |
            object@annotation == "pd.mapping50k.hind240" |
            object@annotation == "pd.mapping50k.xba240" ) {
        tmp <- featureSet[, c("chrom", "physical_pos", "physical_pos", "fsetid", "man_fsetid")]
    } else {
        tmp <- featureSet[, c("chrom", "physical_pos", "physical_pos", "fsetid", "man_fsetid")]
    }
    
    colnames(tmp)[2:3] <- c("start", "end") 
    featureData(eSet) <- new("AnnotatedDataFrame", 
            data = tmp[match(featureData(object)$fsetid[runIdx[, 1]], 
                            tmp$fsetid), ])
    
    ## experiment data
    experimentData(eSet) <- experimentData(object)
    experimentData(eSet)@other$summaryMethod <- summaryMethod
    experimentData(eSet)@other$summaryParam <- summaryParam
    experimentData(eSet)@other$type <- "slData"
    
    ## annotation
    annotation(eSet) <- annotation(object)
    
    t01 <- Sys.time()
    print(difftime(t01, t00))
    
    if (callParam$runtype == "bm") {
        cat(paste(Sys.time(), "|   Saving data \n"))
        slData <- eSet
        save(slData, file = saveFile)
    }
    return(eSet)
}


