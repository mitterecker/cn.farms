#' Summarization Laplacian approach with exact computation
#'
#' This function implements an exact Laplace FARMS algorithm. Users should be 
#' aware, that a change of weight in comparison to the default parameter
#' might also entail a need to change of eps1 and eps2. Unexperienced users 
#' should not change weightZ, since a change in weightZ is also connected to 
#' weight, eps1 and eps2. 
#' 
#' @param probes A matrix with numeric values.
#' @param mu Hyperparameter value which allows to quantify different aspects of
#' potential prior knowledge. Values near zero assumes that most positions do 
#' not contain a signal, and introduces a bias for loading matrix elements near 
#' zero. Default value is 0.
#' @param weight Hyperparameter value which determines the influence of the 
#' Gaussian prior of the loadings
#' @param weightZ Hyperparameter value which determines how strong the Laplace 
#' prior of the factor should be at 0.
#' @param eps1 Epsilon parameter that determines at which values the algorithm 
#' should do an exception handling for low values of the loadings
#' @param eps2 Epsilon parameter that determines at which values the algorithm 
#' should do an exception handling for low values of the data point likelihood 
#' or posterior functions, that tend to be Gaussian
#' @param cyc Number of cycles. If the length is two, it is assumed, that a 
#' minimum and a maximum number of cycles is given. If the length is one, the 
#' value is interpreted as the exact number of cycles to be executed 
#' (minimum == maximum).
#' @param tol States the termination tolerance if cyc[1]!=cyc[2]. 
#' Default is 0.00001.
#' @param weightType Flag, that is used to summarize the loading matrix. 
#' @param centering States how the data is centered. Default is median.
#' @param rescale Rescales the Moments.
#' @param maxIntensity Use of the mode values for building expression values, 
#' if set to TRUE. 
#' @param refIdx index or indices which are used for computation of the 
#' centering
#' @param ... Further parameters for expert users.
#' @return A list including:
#' the found parameters: lambda0, lambda1, Psi
#' 
#' the estimated factors: z (expectation), maxZ (maximum)
#' 
#' p: log-likelihood of the data given the found lambda0, lambda1, 
#' Psi (not the posterior likelihood that is optimized)
#' 
#' varzx: variances of the hidden variables given the data
#' 
#' KL: Kullback Leibler divergences between between posterior and prior 
#' distribution of the hidden variables
#' 
#' IC: Information Content considering the hidden variables and data
#' 
#' ICtransform: transformed Information Content
#' 
#' Case: Case for computation of a sample point (non-exception, special exception)
#' 
#' L1median: Median of the lambda vector components
#' 
#' intensity: back-computed summarized probeset values with mean correction
#' 
#' L_z: back-computed summarized probeset values without mean correction
#' 
#' rawCN: transformed values of L_z
#' 
#' SNR: some additional signal to noise ratio value
#' 
#' @export  
#' @author Andreas Mayr \email{mayr@@bioinf.jku.at} and 
#' Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @examples 
#' x <- matrix(rnorm(100, 11), 20, 5)
#' summarizeFarmsExact(x)

summarizeFarmsExact <- function(
        probes, 
        mu = 0, 
        weight = 0.5, 
        weightZ = 1.0, 
        eps1 = 0.01, 
        eps2 = 10^-10, 
        cyc = c(20, 20), 
        tol = 0.00001, 
        weightType = "mean", 
        centering = "median",
        rescale = FALSE,
        maxIntensity = FALSE,
        refIdx,
        ...) {
    
    probes <- as.matrix(probes)
    sigmaZ <- 1.0 / weightZ
    if(length(cyc) == 1)
        cyc <- c(cyc, cyc)
    additional <- list(...)
    
    if("initPsi" %in% names(additional)) 
        initPsi <- additional$initPsi 
    else 
        initPsi <- 0.5
    
    if("algorithm" %in% names(additional)) 
        algorithm <- additional$algorithm
    else
        algorithm <- "exact"
    
    if(algorithm == "v") {
        if("boundedLapla" %in% names(additional))
            boundedLapla <- additional$boundedLapla
        else
            boundedLapla <- TRUE
    }
    
    ## Initialize variables
    Data <- probes
    DataT <- t(probes)
    n <- ncol(Data)
    dimension <- nrow(Data)
    
    
    
    if (missing(refIdx)) { refIdx <- 1:n }
    
    if (centering == "median") {
        if (length(refIdx) == 1) {
            mean.Data <- median(probes[, refIdx])
        } else {
            mean.Data <- Biobase::rowMedians(probes[, refIdx])
        }
        
    } else if (centering == "mean") {
        if (length(refIdx) == 1) {
            mean.Data <- mean(probes[, refIdx])  
        } else {
            mean.Data <- rowMeans(probes[, refIdx])
        } 
    }
    
    NData <- Data - mean.Data
    
    sd.Data <- sqrt(rowSums(NData^2) / n)
    NData <- NData / (sd.Data + 10^-3)
    NDataT <- t(NData)
    DataCov <- NData %*% NDataT / n
    DataCov[DataCov <= 0] <- 10^-3
    DataCov <- (DataCov + t(DataCov)) / 2
    
    myLambda <- rep(0, dimension)
    PsiLambda <- rep(mean(diag(DataCov)) / (weight * dimension), dimension)
    
    Psi <- initPsi * diag(DataCov)
    lambda <- sqrt(diag(DataCov) - Psi)
    lapla <- rep(1, n)
    
    NData2 <- NData^2
    
    ## EM algorithm
    nrCyc <- cyc[2]
    PsiOld <- Psi
    
    results <- list()
    
    log_p_ges <- 0
    lambda_old <- lambda
    Psi_old <- Psi
    
    for (i in 1:cyc[2]) {
        
        ## E-Step
        if (algorithm == "exact") {
            InvPsi <- 1 / Psi
            av <- rep(-0.5 * (t(lambda) %*% (InvPsi * lambda))[1], n)
            bv <- as.vector((lambda * InvPsi) %*% NData)
            cv <- -0.5 * colSums(NData2 * InvPsi)
            nv <- rep(1 / ((2 * pi)^(dimension / 2) * 
                                prod(Psi)^(1/2) * 2 * sigmaZ), n)
            moments <- .Call("momentsGauss", i, eps1, eps2, av, bv, cv, sigmaZ, 
                    nv, 1, 0, PACKAGE="cn.farms")
            
            if(rescale) {
                sdmom <- sqrt(1 / n * sum(moments$moment2)) / sigmaZ + 10^-3
                moments$moment1 <- moments$moment1 / sdmom 
                moments$moment2 <- moments$moment2 / sdmom^2
            }
            
            avg_xEz <- as.vector(NData %*% moments$moment1 / n)
            avg_Ez2 <- mean(moments$moment2)
            
        } else if (algorithm == "v") {
            InvPsi <- 1 / Psi
            sigma2 <- 1 / (lapla + (t(lambda) %*% (InvPsi * lambda))[1])
            Ez <- as.vector((lambda * InvPsi) %*% NData) * sigma2
            Ez2 <- Ez^2 + sigma2
            avg_xEz <- as.vector(NData %*% Ez / n)
            avg_Ez2 <- mean(Ez2)
            
            ## M-Step only for v
            lapla <- 1 / sqrt(Ez2)
            if(boundedLapla) { lapla[lapla < 1] <- 1 } 
        } else if (algorithm == "g") {
            PsiM1_Lambda <- lambda / Psi
            sigmaZ2 <- 1 / (1 + sum(lambda * PsiM1_Lambda))
            avg_xEz <- (as.vector(DataCov %*% PsiM1_Lambda) * sigmaZ2)
            avg_Ez2 <- sum(PsiM1_Lambda * avg_xEz) * sigmaZ2 + sigmaZ2
        }
        
        
        ## M-Step
        Psi_PsiLambdaM1 <- Psi / PsiLambda
        lambda <- (avg_xEz + Psi_PsiLambdaM1 * myLambda) / 
                (avg_Ez2 * rep(1, dimension) + Psi_PsiLambdaM1)
        lambda <- ifelse(lambda < 0.0, rep(0, length(lambda)), lambda)
        Psi <- diag(DataCov) - avg_xEz * lambda + 
                Psi_PsiLambdaM1 * (myLambda - lambda) * lambda
        Psi[Psi < 10^-3] <- 10^-3
        
        if(i > cyc[1] && max(abs(Psi - PsiOld)) / max(abs(PsiOld)) < tol) {
            nrCyc <- i + 1
            break
        }
        
        PsiOld <- Psi
        
        if(algorithm == "exact") {
            log_p_ges <- sum(log(moments$normConst))
        }    
    }
    
    
    InvPsi <- 1 / Psi
    av <- rep(-0.5 * (t(lambda) %*% (InvPsi * lambda))[1], n)
    bv <- as.vector((lambda * InvPsi) %*% NData)
    cv <- -0.5 * colSums(NData2 * InvPsi)
    nv <- rep(1 / ((2 * pi)^(dimension / 2) * prod(Psi)^(1 / 2) * 2 * sigmaZ), 
            n)
    moments <- .Call("momentsGauss", i, eps1, eps2, av, bv, cv, sigmaZ, nv, 
            1, 0, PACKAGE="cn.farms")
    log_p_ges <- sum(log(moments$normConst))      
    
    
    if(algorithm == "exact") {
        z <- moments$moment1
        varzx <- moments$moment2 - moments$moment1^2
        KL <- moments$CrossEntropy - moments$Entropy
        IC <- log(2 * exp(1.0) * sigmaZ) - moments$Entropy    
    } else if (algorithm == "v") {
        sigma2 <- 1 / (lapla + (t(lambda) %*% (InvPsi * lambda))[1])
        z <- as.vector((lambda * InvPsi) %*% NData) * sigma2
        varzx <- sigma2
        KL <- (0.5 * log(2 * pi) - log(lapla) + 
                    ((z^2 + sigma2) * lapla) / 2.0) - 
                (log(sigma2 * sqrt(2 * pi * exp(1))))
        IC <- 0.5 * log(1.0 + (t(lambda) %*% (InvPsi * lambda))[1] / lapla)    
    } else if (algorithm == "g") {
        PsiM1_Lambda <- lambda / Psi
        sigmaZ2 <- 1 / (1 + sum(lambda * PsiM1_Lambda))
        z <- as.vector(NDataT %*% PsiM1_Lambda) * sigmaZ2
        varzx <- rep(sigmaZ2, ncol(Data))    
        KL <- (0.5 * log(2 * pi) + (z^2 + sigmaZ2 / 2.0)) - 
                (log(sigmaZ2 * sqrt(2 * pi * exp(1))))
        IC <- 0.5 * log(1.0 + (t(lambda) %*% (InvPsi * lambda))[1])        
    }
    
    ICtransform <- 1 / exp(IC * 2.0)
    
	sdz <- sqrt(1 / n * sum(moments$moment2)) / sigmaZ
	if(sdz == 0.0) {
        sdz <- 1
    }
    
    z <- z / sdz
    lambda <- lambda * sdz
    lambda0 <- mean.Data
    lambda1 <- lambda * sd.Data
    Psi <- Psi * sd.Data^2
    
    names(lambda0) <- rownames(probes)
    names(lambda1) <- rownames(probes)
    names(Psi) <- rownames(probes)
    names(z) <- colnames(probes)
    
    if(maxIntensity) {
        zint <- moments$max
    } else {
        zint <- z    
    }
    
    
    if(weightType == "square") {
        PsiLL <- (lambda^2 / Psi)^2
        sumPsiLL <- sum(PsiLL)
        if (sumPsiLL == 0) { sumPsiLL <- 1 }
        propPsiLL <- PsiLL / sumPsiLL
        L_c <- as.vector(crossprod(lambda1, propPsiLL)) * zint
        mean_int <- mean(lambda0)
        express <- L_c + mean_int
        median_int <- median(lambda0)
        rawCN <- (2^(L_c + mean_int) / 2^median_int)
    } else if (weightType == "linear") {
        PsiLL <- (lambda^2 / Psi)
        sumPsiLL <- sum(PsiLL)
        if (sumPsiLL == 0) { sumPsiLL <- 1 }
        propPsiLL <- PsiLL / sumPsiLL
        L_c <- as.vector(crossprod(lambda1, propPsiLL)) * zint
        mean_int <- mean(lambda0)
        express <- L_c + mean_int
        median_int <- median(lambda0)
        rawCN <- (2^(L_c + mean_int) / 2^median_int)
    } else if (weightType == "median") {
        L_c <- median(lambda1) * zint
        mean_int <- median(lambda0)
        express <- L_c + mean_int
        median_int <- median(lambda0)
        rawCN <- (2^(L_c + mean_int) / 2^median_int)
    } else if (weightType == "mean") {
        L_c <- mean(lambda1) * zint
        mean_int <- mean(lambda0)
        express <- L_c + mean_int
        median_int <- median(lambda0)
        rawCN <- (2^(L_c + mean_int) / 2^median_int)
    } else if (weightType == "softmax") {
        PsiLL <- exp(lambda1)
        sumPsiLL <- sum(PsiLL)
        if (sumPsiLL == 0) { sumPsiLL <- 1 }
        propPsiLL <- PsiLL / sumPsiLL
        L_c <- as.vector(crossprod(lambda1, propPsiLL)) * zint
        mean_int <- mean(lambda0)
        express <- L_c + mean_int
        median_int <- median(lambda0)
        rawCN <- (2^(L_c + mean_int) / 2^median_int)
    }
    
    L1median <- median(lambda1)
    
    SNR <- 1 / (1 + (lambda1 %*% (1 / Psi * lambda1)))
    
    return(list(    lambda0     = lambda0, 
                    lambda1     = lambda1, 
                    Psi         = Psi, 
                    z           = z, 
                    maxZ        = moments$max, 
                    p           = log_p_ges, 
                    varzx       = varzx, 
                    KL          = KL, 
                    IC          = IC, 
                    ICtransform = ICtransform, 
                    INICall     = min(ICtransform),
                    INI         = 1 / (1 + sum(lambda1^2 / Psi)),
                    Case        = moments$Case, 
                    L1median    = L1median, 
                    intensity   = as.numeric(express), 
                    L_z         = as.numeric(L_c), 
                    rawCN       = as.numeric(rawCN), 
                    SNR         = SNR))
}
