## Code adopted and improved by Chris Cowing-Zitron
## 2011-01-03

## at the moment excluded due to a very unstable version of VanillaICE
## don't forget to include VanillaICE into Imports of DESCRIPTION files
## and NAMESPACE: 
## @importFrom VanillaICE hmm.setup
## @importFrom VanillaICE cnEmission
## @export


#' Postprocessing of summarized data with VanillaICE
#' @param mlData data from cn.farms
#' @param runmed Integer. Size of window for running median. A running median 
#' of the total copy number is used to estimate the probability that a copy 
#' number estimate is an outlier. 
#' @param tau Numeric. Factor for scaling the distance weighted transition probability. 
#' @param copynumberstates The mean values for each of the hidden states. 
#' @param states Character vector indicating the names of the hidden states
#' @param normalindex Integer. Indicates which element of the states vector 
#' corresponds to the 'normal' state.
#' @param initiallogprobs Numeric vector of the same length as the number of 
#' states. Specifies the initial state probabilities on the log scale.
#' @return Some data
#' @author Andreas Mitterecker
#' @noRd
#mlDataToHmm <- function(
#        mlData, 
#        runmed = 5, 
#        tau = 1e+08, 
#        copynumberstates = c(0.1, 1, 2, 3, 4), 
#        states = c("homDel", "hemDel", "normal", "3copy", "4copy"), 
#        normalindex = 3, 
#        initiallogprobs = log(rep(1 / length(states), length(states)))) {
#    
#    #library(VanillaICE)
#   
#    
#    lzdata <- assayData(mlData)$L_z[,1]
#    for (sample in 2:length(sampleNames(mlData))) { 
#        lzdata <- cbind(lzdata,assayData(mlData)$L_z[,sample]) 
#    }
#    rownames(lzdata) <- NULL
#    colnames(lzdata) <- NULL
#    mycnset <- new("CopyNumberSet", copyNumber = lzdata, annotation = "pd.mapping250k.nsp")
#    phInf <- fData(mlData)
#    phInf$man_fsetid[length(phInf$man_fsetid)] <- "None"
#    featureNames(mycnset) <- phInf$man_fsetid
#    sampleNames(mycnset) <- sampleNames(mlData)
#    featureData(mycnset)[["chromosome"]] <- phInf$chrom
#    featureData(mycnset)[["position"]] <- phInf$start
#    featureData(mycnset)[["isSnp"]] <- rep(TRUE, length(featureNames(mycnset)))
#    hmmOpts <- hmm.setup(mycnset, is.log = TRUE, tau = tau, 
#            copynumberStates = copynumberstates, states = states, 
#            normalIndex = normalindex, log.initialPr = initiallogprobs)
#    mylogbeta <- cnEmission(mycnset, hmmOpts, k = runmed, 
#            cnStates = hmmOpts[["copynumberStates"]], 
#            is.log = TRUE, normalIndex = normalindex, verbose = TRUE)
#    dimnames(mylogbeta) <- list(featureNames(mycnset), sampleNames(mycnset), hmmOpts$states)
#    myviterbi <- VanillaICE:::viterbi(mycnset, hmmOpts, log.E = mylogbeta)
#    namedstate <- vector(length = length(myviterbi[["state"]]))
#    for (state in 1:length(states)) { 
#        namedstate[myviterbi[["state"]] == state] <- states[state] 
#    }
#    myviterbi[["namedstate"]] <- namedstate
#    
#    return(myviterbi)
#    
#}
