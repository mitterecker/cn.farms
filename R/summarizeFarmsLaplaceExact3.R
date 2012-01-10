#' Summarization Laplacian approach with exact computation
#'
#' This function implements an exact Laplace FARMS algorithm. 
#' 
#' @param probes A matrix with numeric values.
#' @param mu Hyperparameter value which allows to quantify different aspects of potential prior knowledge. Values near zero assumes that most positions do not contain a signal, and introduces a bias for loading matrix elements near zero. Default value is 0 and it's recommended not to change it.
#' @param weight Hyperparameter value which determines the influence of the Gaussian prior of the loadings
#' @param weightSignal Hyperparameter value on the signal.
#' @param weightZ Hyperparameter value which determines how strong the Laplace prior of the factor should be at 0. Users should be aware, that a change of weightZ in comparison to the default parameter might also entail a need to change other parameters. Unexperienced users should not change weightZ.
#' @param weightProbes Parameter (TRUE/FALSE), that determines, if the number of probes should additionally be considered in weight. If TRUE, weight will be modified.
#' @param updateSignal updateSignal.
#' @param cyc Number of cycles. If the length is two, it is assumed, that a minimum and a maximum number of cycles is given. If the length is one, the value is interpreted as the exact number of cycles to be executed (minimum == maximum).
#' @param tol States the termination tolerance if cyc[1]!=cyc[2]. Default is 0.00001.
#' @param weightType Flag, that is used to summarize the probes of a sample. 
#' @param centering States how the data should be centered ("mean", "median"). Default is median.
#' @param rescale Parameter (TRUE/FALSE), that determines, if moments in exact Laplace FARMS are rescaled in each iteration. Default is FALSE.
#' @param backscaleComputation Parameter (TRUE/FALSE), that determines if the moments of hidden variables should be reestimated after rescaling the parameters.
#' @param maxIntensity Parameter (TRUE/FALSE), that determines if the expectation value (=FALSE) or the maximum value (=TRUE) of p(z|x_i) should be used for an estimation of the hidden varaible. 
#' @param refIdx index or indices which are used for computation of the centering
#' @param ... Further parameters for expert users.
#' 
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
#' 
#' @author Andreas Mayr \email{mayr@@bioinf.jku.at} and Djork-Arne Clevert \email{okko@@clevert.de} and Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @examples 
#' x <- matrix(rnorm(100, 11), 20, 5)
#' summarizeFarmsExact(x)
#' @export

summarizeFarmsExact3 <- function(
        probes, 
        mu=1, 
        weight=100, 
        weightSignal=1.0,
        weightZ=30.0, 
        weightProbes=TRUE,
        updateSignal=FALSE,
        cyc=c(10, 10), 
        tol=0.00001, 
        weightType="mean", 
        centering="median",
        rescale=FALSE,
        backscaleComputation=FALSE,
        maxIntensity=TRUE,
        refIdx,
        ...) {
    eps1 <- 0.01
    eps2 <- 10^-10
    updateWeights <- FALSE
    
    probes <- as.matrix(probes)
    sigmaZ <- 1/sqrt(2*weightZ)
    if(length(cyc) == 1)
        cyc <- c(cyc, cyc)
    additional <- list(...)
    
    if("initPsi" %in% names(additional)) 
        initPsi <- additional$initPsi 
    else 
        initPsi <- 0.1
    
    if("init" %in% names(additional)) 
        init <- additional$init
    else
        init <- "mean"
    
    if("algorithm" %in% names(additional)) 
        algorithm <- additional$algorithm
    else
        algorithm <- "exact"
    
    if(algorithm=="v"||algorithm=="ard") {
        if("boundedLapla" %in% names(additional))
            boundedLapla <- additional$boundedLapla
        else
            boundedLapla <- TRUE
        
        if("spuriousCorrelation" %in% names(additional))
            spuriousCorrelation <- additional$spuriousCorrelation
        else
            spuriousCorrelation <- 1.0
    }
    
    if(algorithm=="ard") {
        if("ard2Max" %in% names(additional))
            ard2Max <- additional$ard2Max
        else
            ard2Max <- 0
    }
    if(algorithm=="g")
        sigmaZ <- 1.0 / sqrt(weightZ)
    
    ## Initialize variables
    Data <- probes
    DataT <- t(probes)
    n <- ncol(Data)
    dimension <- nrow(Data)
    
    if(algorithm=="mean") {
        allMedian <- median(Data) 
        twoMedian <- median(apply(Data, 1, median))
        allMean <- mean(Data) 
        twoMean <- mean(apply(Data, 1, mean))
        if (centering == "median") { mean.Data <- apply(Data, 1, median) }
        if (centering == "mean") { mean.Data <- rowMeans(Data) } 
        NData <- Data - mean.Data
        sd.Data <- sqrt(rowSums(NData^2) / n) + 10^-3
        NData <- NData / (sd.Data)
        return(list(colMean=colMeans(NData), colMedian=apply(NData, 2, median), allMedian=allMedian, twoMedian=twoMedian, allMean=allMean, twoMean=twoMean))
    }
    
    if (missing(refIdx)) {
        refIdx <- 1:n 
    }
    if (centering == "median") { 
        if (length(refIdx) == 1) {
            mean.Data <- Data[, refIdx]
        }
        else {
            mean.Data <- Biobase::rowMedians(Data[, refIdx])
        }
    }
    if (centering == "mean") { 
        if (length(refIdx) == 1) {
            mean.Data <- Data[, refIdx]
        }
        else {
            mean.Data <- rowMeans(Data[, refIdx])
        }
    }
    NData <- Data - mean.Data
    
    sd.Data <- sqrt(rowSums(NData^2) / n) + 10^-3
    NData <- NData / (sd.Data)
    NDataT <- t(NData)
    DataCov <- NData %*% NDataT / n
    DataCov[DataCov <= 0] <- 0.0
    DataCov <- (DataCov + t(DataCov)) / 2
    
    myLambda <- rep(mu, dimension)
    PsiLambda <- rep(1.0/weight, dimension)
    if(weightProbes==TRUE)
        PsiLambda <- rep(1.0/(weight*dimension), dimension)
    sigmaS <- 1/sqrt(weightSignal)
    
    s <- 1.0
    Psi <- initPsi * diag(DataCov)
    lambda <- sqrt(diag(DataCov) - Psi)
    lambdaOld <- lambda
    lapla <- rep(1, n)
    if(init=="maxAbsSample") {
        sampleIdx <- which.max(abs(apply(NData, 2, mean)))
        lambda <- (abs(NData[, sampleIdx])/max(abs(NData[, sampleIdx])))*sqrt(1-initPsi)
        Psi <- diag(DataCov)-lambda^2
    }
    else if(init=="maxAbs") {
        maxAcrossSamples <- apply(NData, 1, function(x) max(abs(x)))
        lambda <- (maxAcrossSamples/max(maxAcrossSamples))*sqrt(1-initPsi)
        Psi <- diag(DataCov)-lambda^2
    }
    else if(init=="absEqual") {
        sumSamples <- abs(rowSums(NData))
        lambda <- (sumSamples/max(sumSamples))*sqrt(1-initPsi)
        Psi <- diag(DataCov)-lambda^2
    }
    
    
    NData2 <- NData^2
    
    #EM algorithm
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
            nv <- rep(1 / ((2 * pi)^(dimension / 2) * prod(Psi)^(1/2) * 2 * sigmaZ), n)
            moments <- .Call("momentsGauss", i, eps1, eps2, av, bv, cv, sigmaZ, nv, 1, 0, PACKAGE="cn.farms")
            if(rescale) {
                sdmom <- sqrt(1/n*sum(moments$moment2))/(sqrt(2.0)*sigmaZ)+10^-3
                moments$moment1 <- moments$moment1/sdmom 
                moments$moment2 <- moments$moment2/sdmom^2
            }
            avg_xEz <- as.vector(NData %*% moments$moment1 / n)
            avg_Ez2 <- mean(moments$moment2)
            
            if(updateWeights) {
                avg_absEz <- mean(moments$absMoment1)
                sigmaZ <- avg_absEz
            }
        } else if (algorithm == "v") {
            InvPsi <- 1 / Psi
            sigma2 <- 1 / (lapla + (t(lambda) %*% (InvPsi * lambda))[1])
            Ez <- as.vector((lambda * InvPsi) %*% NData) * sigma2
            Ez2 <- Ez^2 + sigma2
            avg_xEz <- as.vector(NData %*% Ez / n)
            avg_Ez2 <- mean(Ez2)
            
            ## M-Step only for v
            lapla <- 1 / sqrt(Ez2)
            if(boundedLapla)
                lapla[lapla < spuriousCorrelation] <- spuriousCorrelation 
        } else if (algorithm == "ard") {
            InvPsi <- 1 / Psi
            sigma2 <- 1 / (lapla + (t(lambda) %*% (InvPsi * lambda))[1])
            Ez <- as.vector((lambda * InvPsi) %*% NData) * sigma2
            Ez2 <- Ez^2 + sigma2
            avg_xEz <- as.vector(NData %*% Ez / n)
            avg_Ez2 <- mean(Ez2)
            
            ## M-Step only for ard
            lapla <- 1 / Ez2
            if(boundedLapla)
                lapla[lapla < spuriousCorrelation] <- spuriousCorrelation 
        } else if (algorithm == "g") {
            InvPsi <- 1 / Psi
            sigma2 <- 1 / (1/sigmaZ^2 + (t(lambda) %*% (InvPsi * lambda))[1])
            avg_xEz <- (as.vector(DataCov %*% (InvPsi * lambda)) * sigma2)
            avg_Ez2 <- sum((InvPsi * lambda) * avg_xEz) * sigma2 + sigma2
            if(updateWeights) {
                sigmaZ <- sqrt(avg_Ez2)
            }
        }
        
        lambda <- lambdaOld
        
        ## M-Step
        if(updateSignal) {
            Psi_PsiLambdaM1 <- Psi / PsiLambda
            lambda <- (avg_xEz * s + Psi_PsiLambdaM1 * myLambda) / (avg_Ez2 * s^2 * rep(1, dimension) + Psi_PsiLambdaM1)
            lambda <- ifelse(lambda < 0.0, rep(0, length(lambda)), lambda)
            lambda <- (lambda/sqrt(sum(lambda^2)))*sqrt(sum(myLambda^2))
            s <- sum(lambda*avg_xEz/Psi)/(1/sigmaS^2+sum(lambda^2/Psi)*avg_Ez2)
            Psi <- diag(DataCov) - avg_xEz * s * lambda + Psi_PsiLambdaM1 * (myLambda - lambda) * lambda
            Psi[Psi < 10^-3] <- 10^-3
        }
        else {
            s <- 1
            Psi_PsiLambdaM1 <- Psi / PsiLambda
            lambda <- (avg_xEz * s + Psi_PsiLambdaM1 * myLambda) / (avg_Ez2 * s^2 * rep(1, dimension) + Psi_PsiLambdaM1)
            lambda <- ifelse(lambda < 0.0, rep(0, length(lambda)), lambda)
            Psi <- diag(DataCov) - avg_xEz * s * lambda + Psi_PsiLambdaM1 * (myLambda - lambda) * lambda
            Psi[Psi < 10^-3] <- 10^-3
        }
        
        lambdaOld <- lambda
        lambda <- lambda*s
        
        
        if(i > cyc[1] && max(abs(Psi - PsiOld)) / max(abs(PsiOld)) < tol) {
            nrCyc <- i + 1
            break
        }
        
        PsiOld <- Psi
    }
    
    InvPsi <- 1 / Psi
    av <- rep(-0.5 * (t(lambda) %*% (InvPsi * lambda))[1], n)
    bv <- as.vector((lambda * InvPsi) %*% NData)
    cv <- -0.5 * colSums(NData2 * InvPsi)
    nv <- rep(1 / ((2 * pi)^(dimension/2) * prod(Psi)^(1 / 2) * 2 * sigmaZ), n)
    moments <- .Call("momentsGauss", i, eps1, eps2, av, bv, cv, sigmaZ, nv, 1, 0, PACKAGE="cn.farms")
    
    if(algorithm=="exact") {
        z <- moments$moment1
        maxZ <- moments$max
        varzx <- moments$moment2 - moments$moment1^2
        KL <- moments$CrossEntropy - moments$Entropy
        IC <- log(2 * exp(1.0) * sigmaZ) - moments$Entropy  
        sdz <- sqrt(1/n*sum(moments$moment2))/(sqrt(2.0)*sigmaZ)
    } else if ((algorithm == "v")||(algorithm == "ard")) {
        sigma2 <- 1 / (lapla + (t(lambda) %*% (InvPsi * lambda))[1])
        z <- as.vector((lambda * InvPsi) %*% NData) * sigma2
        maxZ <- z
        if(algorithm=="ard"&&(ard2Max>0)) {
            disc <- 8*av+bv^2*ard2Max
            maxZ[] <- 0.0
            maxZP <- maxZ
            maxZM <- maxZ
            d0 <- disc>0.0
            maxZP[d0] <- 2.0/(bv[d0]*ard2Max+sqrt(ard2Max*disc[d0]))
            maxZM[d0] <- 2.0/(bv[d0]*ard2Max-sqrt(ard2Max*disc[d0]))
            d0P <- (disc>0.0)&(maxZP>0.0)
            maxZ[d0P] <- maxZM[d0P]
            d0M <- (disc>0.0)&(maxZM<0.0)
            maxZ[d0M] <- maxZP[d0M]
        }
        varzx <- sigma2
        KL <- (0.5 * log(2 * pi) - 0.5*log(lapla) + ((z^2 + sigma2) * lapla) / 2.0) - 0.5*(log(2 * sigma2 * pi * exp(1)))
        IC <- 0.5 * log(1.0 + (t(lambda) %*% (InvPsi * lambda))[1] / lapla) 
        sdz <- sqrt(1/n*sum(z^2+varzx))
    } else if (algorithm == "g") {
        InvPsi <- 1 / Psi
        sigma2 <- 1 / (1/sigmaZ^2 + (t(lambda) %*% (InvPsi * lambda))[1])
        z <- as.vector(NDataT %*% (InvPsi * lambda)) * sigma2
        maxZ <- z
        varzx <- rep(sigma2, ncol(Data))
        KL <- (0.5 * log(2 * pi) - 0.5*log(1/sigmaZ) + ((z^2 + sigma2) * (1/sigmaZ)) / 2.0) - 0.5*(log(2 * sigma2 * pi * exp(1)))
        IC <- 0.5 * log(1.0 + (t(lambda) %*% (InvPsi * lambda))[1] / (1/sigmaZ)) 
        sdz <- sqrt(1/n*sum(z^2+varzx))/sigmaZ
    }
    ICtransform <- 1 / exp(IC * 2.0)
    
    if(sdz == 0.0)
        sdz <- 1
    
    z <- z / sdz
    maxZ <- maxZ / sdz
    lambda <- lambda * sdz
    if(algorithm=="v")
        lapla <- lapla/sdz
    else if (algorithm=="ard")
        lapla <- lapla/sdz^2
    
    avRET <- rep(-0.5 * (t(lambda) %*% (InvPsi * lambda))[1], n)
    bvRET <- as.vector((lambda * InvPsi) %*% NData)
    cvRET <- -0.5 * colSums(NData2 * InvPsi)
    nvRET <- rep(1 / ((2 * pi)^(dimension/2) * prod(Psi)^(1 / 2) * 2 * sigmaZ), n)
    momentsRET <- .Call("momentsGauss", i, eps1, eps2, avRET, bvRET, cvRET, sigmaZ, nvRET, 1, 0, PACKAGE="cn.farms")
    cvCOMP <- -0.5 * colSums(NData2 * rep(1.0, dimension))
    nvCOMP <- rep(1 / ((2 * pi)^(dimension/2)), n)
    MLQ <- pchisq(2*(sum(log(momentsRET$normConst)-(cvCOMP+log(nvCOMP)))), dimension)
    log_p_ges <- sum(log(momentsRET$normConst))
    
    if(algorithm=="exact") {
        zScaled <- momentsRET$moment1
        maxZScaled <- momentsRET$max
        varzxScaled <- momentsRET$moment2 - momentsRET$moment1^2
        KLScaled <- momentsRET$CrossEntropy - momentsRET$Entropy
        ICScaled <- log(2 * exp(1.0) * sigmaZ) - momentsRET$Entropy    
    } else if ((algorithm == "v")||(algorithm == "ard")) {
        sigma2Scaled <- 1 / (lapla + (t(lambda) %*% (InvPsi * lambda))[1])
        zScaled <- as.vector((lambda * InvPsi) %*% NData) * sigma2Scaled
        maxZScaled <- zScaled
        if(algorithm=="ard"&&(ard2Max>0)) {
            disc <- 8*avRET+bvRET^2*ard2Max
            maxZScaled[] <- 0.0
            maxZScaledP <- maxZScaled
            maxZScaledM <- maxZScaled
            d0 <- disc>0.0
            maxZScaledP[d0] <- 2.0/(bvRET[d0]*ard2Max+sqrt(ard2Max*disc[d0]))
            maxZScaledM[d0] <- 2.0/(bvRET[d0]*ard2Max-sqrt(ard2Max*disc[d0]))
            d0P <- (disc>0.0)&(maxZScaledP>0.0)
            maxZScaled[d0P] <- maxZScaledM[d0P]
            d0M <- (disc>0.0)&(maxZScaledM<0.0)
            maxZScaled[d0M] <- maxZScaledP[d0M]
        }
        varzxScaled <- sigma2Scaled
        KLScaled <- (0.5 * log(2 * pi) - 0.5*log(lapla) + ((zScaled^2 + sigma2Scaled) * lapla) / 2.0) - 0.5*(log(2 * sigma2Scaled * pi * exp(1)))
        ICScaled <- 0.5 * log(1.0 + (t(lambda) %*% (InvPsi * lambda))[1] / lapla)    
    } else if (algorithm == "g") {
        sigma2Scaled <- 1 / (1/sigmaZ^2 + (t(lambda) %*% (InvPsi * lambda))[1])
        zScaled <- as.vector(NDataT %*% (InvPsi * lambda)) * sigma2Scaled
        maxZScaled <- zScaled
        varzxScaled <- rep(sigma2Scaled, ncol(Data))
        KLScaled <- (0.5 * log(2 * pi) - 0.5*log(1/sigmaZ) + ((zScaled^2 + sigma2Scaled) * (1/sigmaZ)) / 2.0) - 0.5*(log(2 * sigma2Scaled * pi * exp(1)))
        ICScaled <- 0.5 * log(1.0 + (t(lambda) %*% (InvPsi * lambda))[1] / (1/sigmaZ)) 
    }
    ICtransformScaled <- 1 / exp(ICScaled * 2.0)
    
    if(backscaleComputation==TRUE) {
        z <- zScaled
        maxZ <- maxZScaled
    }
    
    if(maxIntensity)
        zint <- maxZ
    else
        zint <- z
    
    lambdaNormalized <- lambda
    PsiNormalized <- Psi
    lambda0 <- mean.Data
    lambda1 <- lambda * sd.Data
    Psi <- Psi * sd.Data^2
    rm("lambda")
    
    names(lambdaNormalized) <- rownames(probes)
    names(PsiNormalized) <- rownames(probes)
    names(z) <- colnames(probes)
    names(maxZ) <- colnames(probes)
    names(lambda0) <- rownames(probes)
    names(lambda1) <- rownames(probes)
    names(Psi) <- rownames(probes)
    
    lambdaMeanNormalized <- mean(lambdaNormalized)
    lambdaMedianNormalized <- median(lambdaNormalized)
    lambda1Mean <- mean(lambda1)
    lambda1Median <- median(lambda1)
    SNR <- 1/(1+sum(lambda1^2/Psi))
    
    #linear
    PsiLL <- (lambda1^2 / Psi)
    sumPsiLL <- sum(PsiLL)
    if (sumPsiLL == 0) { 
        sumPsiLL <- 1
    }
    propPsiLL <- PsiLL / sumPsiLL
    L_c_Linear <- as.vector(crossprod(lambda1, propPsiLL)) * zint
    L_c_Linear_Normalized <- as.vector(crossprod(lambdaNormalized, propPsiLL)) * zint
    expressLinear <- mean(lambda0) + L_c_Linear
    
    #square
    PsiLL <- (lambda1^2 / Psi)^2
    sumPsiLL <- sum(PsiLL)
    if (sumPsiLL == 0) { 
        sumPsiLL <- 1
    }
    propPsiLL <- PsiLL / sumPsiLL
    L_c_Square <- as.vector(crossprod(lambda1, propPsiLL)) * zint
    L_c_Square_Normalized <- as.vector(crossprod(lambdaNormalized, propPsiLL)) * zint
    expressSquare <- mean(lambda0) + L_c_Square
    
    #softmax
    PsiLL <- exp(lambda1)
    sumPsiLL <- sum(PsiLL)
    if (sumPsiLL == 0) { 
        sumPsiLL <- 1
    }
    propPsiLL <- PsiLL / sumPsiLL
    L_c_Softmax <- as.vector(crossprod(lambda1, propPsiLL)) * zint
    L_c_Softmax_Normalized <- as.vector(crossprod(lambdaNormalized, propPsiLL)) * zint
    expressSoftmax <- mean(lambda0) + L_c_Softmax
    
    #mean
    L_c_Mean <- mean(lambda1)*zint
    L_c_Mean_Normalized <- mean(lambdaNormalized)*zint
    expressMeanBC <- mean(lambda0)+mean(lambda1)*zint
    expressBCMean <- colMeans(lambda0+lambda1%o%zint)
    
    #median
    L_c_Median <- median(lambda1)*zint
    L_c_Median_Normalized <- median(lambdaNormalized)*zint
    expressMedianBC <- median(lambda0)+median(lambda1)*zint
    expressBCMedian <- apply(lambda0+lambda1%o%zint, 2, median)
    
    if (weightType == "linear") {
        L_c <- L_c_Linear
        L_c_Normalized <- L_c_Linear_Normalized
        express <- expressLinear
    } else if(weightType == "square") {
        L_c <- L_c_Square
        L_c_Normalized <- L_c_Square_Normalized
        express <- expressSquare
    } else if (weightType == "softmax") {
        L_c <- L_c_Softmax
        L_c_Normalized <- L_c_Softmax_Normalized
        express <- expressSoftmax
    } else if ((weightType == "medianBC")||(weightType == "median")) {
        L_c <- L_c_Median
        L_c_Normalized <- L_c_Median_Normalized
        express <- expressMedianBC
    } else if ((weightType == "meanBC")||(weightType == "mean")) {
        L_c <- L_c_Mean
        L_c_Normalized <- L_c_Mean_Normalized
        express <- expressMeanBC
    } else if (weightType == "BCMedian") { 
        L_c <- L_c_Median
        L_c_Normalized <- L_c_Median_Normalized
        express <- expressBCMedian
    } else if (weightType == "BCMean") {
        L_c <- L_c_Mean
        L_c_Normalized <- L_c_Mean_Normalized
        express <- expressBCMean
    }
    
    cnINIMean <- mean(lambda1)*zint
    cnINIMedian <- median(lambda1)*zint
    cnINIMeanSignal <- cnINIMean*((sum(lambda1^2/Psi))/(1+sum(lambda1^2/Psi)))
    cnINIMedianSignal <- cnINIMedian*((sum(lambda1^2/Psi))/(1+sum(lambda1^2/Psi)))
    
    return(list(
                    sdData=sd.Data,
                    lambdaNormalized=lambdaNormalized,
                    PsiNormalized=PsiNormalized,
                    lambda0=lambda0, 
                    lambda1=lambda1, 
                    Psi=Psi,
                    z=z, 
                    maxZ=maxZ, 
                    
                    ExactMaxZ=moments$max,
                    ExactMaxZScaled=momentsRET$max,
                    ExactP=log_p_ges, 
                    ExactMLQ=MLQ,
                    ExactCase=momentsRET$Case,
                    IndividualExactKL=momentsRET$CrossEntropy - momentsRET$Entropy,
                    SummaryExactKL=mean(momentsRET$CrossEntropy - momentsRET$Entropy),
                    IndividualExactKLScaled=moments$CrossEntropy - moments$Entropy,
                    SummaryExactKLScaled=mean(moments$CrossEntropy - moments$Entropy),
                    
                    IndividualKL=KL, 
                    SummaryKL=mean(KL),
                    IndividualKLScaled=KLScaled,
                    SummaryKLScaled=mean(KLScaled),
                    
                    IndividualIC=IC, 
                    SummaryIC=mean(IC),
                    IndividualICScaled=ICScaled,
                    SummaryICScaled=mean(ICScaled),
                    
                    ICtransform=ICtransform, 
                    ICtransformScaled=ICtransformScaled,
                    IndividualINICall=ICtransformScaled,
                    SummaryINICall=min(ICtransformScaled),
                    
                    varzx=varzx, 
                    varzxScaled=varzxScaled,
                    
                    lambdaMeanNormalized=lambdaMeanNormalized,
                    lambdaMedianNormalized=lambdaMedianNormalized,
                    lambda1Mean=lambda1Mean,
                    lambda1Median=lambda1Median,
                    
                    SNR=SNR,
                    cnINIMean=cnINIMean,
                    cnINIMedian=cnINIMedian,
                    cnINIMeanSignal=cnINIMeanSignal,
                    cnINIMedianSignal=cnINIMedianSignal,
                    
                    
                    L_z=as.numeric(L_c),
                    L_z_Normalized=as.numeric(L_c_Normalized),
                    intensity=as.numeric(express), 
                    
                    L_c_Linear=L_c_Linear,
                    L_c_Linear_Normalized=L_c_Linear_Normalized,
                    expressLinear=expressLinear,
                    L_c_Square=L_c_Square,
                    L_c_Square_Normalized=L_c_Square_Normalized,
                    expressSquare=expressSquare,
                    L_c_Softmax=L_c_Softmax,
                    L_c_Softmax_Normalized=L_c_Softmax_Normalized,
                    expressSoftmax=expressSoftmax,
                    L_c_Mean=L_c_Mean,
                    L_c_Mean_Normalized=L_c_Mean_Normalized,
                    expressMeanBC=expressMeanBC,
                    expressMedianBC=expressMedianBC,
                    L_c_Median=L_c_Median,
                    L_c_Median_Normalized=L_c_Median_Normalized,
                    expressBCMean=expressBCMean,
                    expressBCMedian=expressBCMedian
            ))
}
