#' Summarization Gaussian approach
#' 
#' This function runs the FARMS algorithm.
#' 
#' @param probes A matrix with numeric values.
#' @param weight Hyperparameter value in the range of [0,1] which determines 
#' the influence of the prior.
#' @param mu Hyperparameter value which allows to quantify different aspects of
#' potential prior knowledge. Values near zero assumes that most genes do not
#' contain a signal, and introduces a bias for loading matrix elements near 
#' zero. Default value is 0.
#' @param cyc Number of cycles for the EM algorithm.
#' @param tol States the termination tolerance. Default is 0.00001.
#' @param weightType Flag, that is used to summarize the loading matrix. 
#' The default value is set to mean.
#' @param init Parameter for estimation.
#' @param correction Value that indicates whether the covariance matrix should 
#' be corrected for negative eigenvalues which might emerge from the 
#' non-negative correlation constraints or not. 
#' Default = O  (means that no correction is done), 
#' 1 (minimal noise (0.0001) is added to the diagonal elements of the 
#' covariance matrix to force positive definiteness), 
#' 2 (Maximum Likelihood solution to compute the nearest positive definite 
#' matrix under the given non-negative correlation constraints of the covariance
#' matrix)
#' @param minNoise States the minimal noise. Default is 0.35.
#' @param centering States how the data is centered. Default is median.
#' @param refIdx index or indices which are used for computation of the 
#' centering
#' @return A list containing the results of the run.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @examples
#' x <- matrix(rnorm(100, 11), 20, 5)
#' summarizeFarmsGaussian(x)
summarizeFarmsGaussian <- function(probes, 
        weight = 0.15, 
        mu = 0, 
        cyc = 10, 
        tol = 0.0001, 
        weightType = "mean", 
        init = 0.6, 
        correction = 0, 
        minNoise = 0.35, 
        centering = "median", 
        refIdx) {
    
    a_old <- 0.5
    
    n_array <-  ncol(probes)
    
    n_probes <- nrow(probes)
    
    if (missing(refIdx)) { refIdx <- 1:n_array }
    
    if (centering == "median") {
        if (length(refIdx) == 1) {
            mean.probes <- median(probes[, refIdx])
        } else {
            mean.probes <- Biobase::rowMedians(probes[, refIdx])
        }
        
    } else if (centering == "mean") {
        if (length(refIdx) == 1) {
            mean.probes <- mean(probes[, refIdx])  
        } else {
            mean.probes <- rowMeans(probes[, refIdx])
        } 
    }
    
    centered.probes <- probes - mean.probes
    
    sd.probes <- sqrt(diag(crossprod(t(centered.probes))) / n_array) 
    
    if (0 %in% sd.probes) {
        
        index <- which(sd.probes == 0)
        
        sd.probes[index] <- 1
        
        probes <- probes / sd.probes 
        
        x <- t(probes)
        
        if (centering == "median") {y_v <- apply(x, 2, median)}
        
        if (centering == "mean") {y_v <- colMeans(x)}
        
        xmean <- matrix(y_v, n_array, n_probes, byrow = TRUE)
        
        X <- x - xmean 
        
        XX <- crossprod(X,X) / n_array
        
        diag(XX)[index] <- 1
        
    } else {
        
        probes <- probes / sd.probes
        
        x <- t(probes)
        
        if (centering == "median") {y_v <- apply(x, 2, median)}
        
        if (centering == "mean") {y_v <- colMeans(x)}
        
        xmean <- matrix(y_v, n_array, n_probes, byrow = TRUE)
        
        X <- x - xmean
        
        XX <- crossprod(X, X) / n_array
        
    }
    
    XX <- (XX + t(XX)) / 2 
    
    XX[which(XX < 0.1)] <- 0
    
    minEigenValues <- -1
    
    if (correction >= 1) {
        
        while(minEigenValues < 0){
            
            eigen_XX <- eigen(XX)
            
            eigenValues_XX <- eigen_XX$values
            
            eigenVectors_XX <- eigen_XX$vectors
            
            minEigenValues <- min(eigenValues_XX)
            
            if (correction < 2) {
                
                if(minEigenValues<minNoise){
                    
                    diag(XX)<-diag(XX)+(minNoise - minEigenValues)
                    
                }
                
            } else {
                
                if(minEigenValues < minNoise){
                    
                    eigenValues_XX[which(eigenValues_XX<minNoise)] <- minNoise
                    
                    XX <- eigenVectors_XX%*%diag(eigenValues_XX) %*% t(eigenVectors_XX)
                    
                }
                
            }
            
        }
        
    }
    
    diagXX <- diag(XX)
    
    XX_tmp <- XX
    
    diag(XX_tmp) <- 0
    
    L <- init * sqrt(apply(XX_tmp,1,max))
    
    Ph <- 1 - L^2
    
    alpha <- weight 
    
    bbeta <- mu * alpha
    
    
    
    for (i in 1:cyc) {
        
        # E Step
        
        PsiL <- (1 / Ph) * L
        
        a <- 1 / as.vector(1 + crossprod(L, PsiL))
        
        bar <- PsiL * a
        
        beta <- t(bar)
        
        XXbeta <- XX %*% bar
        
        EZZ <- a + beta %*% XXbeta
        
        t_XXbeta <- XXbeta + Ph * bbeta
        
        t_EZZ <- as.vector(EZZ) + Ph * alpha
        
        ## M Step
        
        L <- t_XXbeta / t_EZZ
        
        Ph <- diagXX - XXbeta * L + Ph * alpha * L * (mu - L) 
        
        
        if (sqrt((a_old - a)^2) < tol) {
            
            break
            
        }
        
        a_old <- a
        
    }
    
    c <- X %*% bar ## hidden variable c - factor
    
    if (EZZ == 0) {
        
        var_z_scale <- 1 ## avoiding division by zero
        
    } else {
        
        var_z_scale <- sqrt(EZZ)
        
    }
    
    c <- c / as.vector(var_z_scale)
    
    L <- L * as.vector(var_z_scale)
    
    PsiL <- (1 / Ph) * L
    
    a <- as.vector(1 + crossprod(L,PsiL))
    
    SNR <- 1 / a 
    
    SIG <- as.vector(crossprod(L, diag(as.vector(1/Ph)))) %*% XX %*% 
            diag(as.vector(1/Ph)) %*% L * a^-2
    
    signal_info <- numeric(length = 4)
    
    signal_info[1] <- SNR
    
    signal_info[2] <- SIG
    
    signal_info[3] <- EZZ
    
    signal_info[4] <- i
    
    if(weightType == "square") {
        
        PsiLL <- ((1 / Ph) * L^2 )^2
        
        sumPsiLL <- sum(PsiLL)
        
        if(sumPsiLL == 0) { sumPsiLL <- 1 }
        
        propPsiLL <- PsiLL / sumPsiLL
        
        L_c <- as.vector(crossprod(L * sd.probes, propPsiLL)) * c
        
        mean_int <- mean(y_v * sd.probes)
        
        express <- L_c + mean_int
        
        median_int <- median(y_v * sd.probes)
        
        rawCN <- (2^(L_c + mean_int) / 2^median_int)
        
    } else if (weightType == "linear") {
        
        PsiLL <- (1 / Ph) * L^2
        
        sumPsiLL <- sum(PsiLL)
        
        if (sumPsiLL == 0) { sumPsiLL <- 1 }
        
        propPsiLL <- PsiLL / sumPsiLL
        
        L_c <- as.vector(crossprod(L * sd.probes, propPsiLL)) * c
        
        mean_int <- mean(y_v * sd.probes)
        
        express <- L_c + mean_int
        
        median_int <- median(y_v * sd.probes)
        
        rawCN <- (2^(L_c + mean_int) / 2^median_int)
        
    } else if (weightType == "median") {
        
        L_c <- median(L * sd.probes) * c
        
        mean_int <- median(y_v * sd.probes)
        
        express <- L_c + mean_int
        
        median_int <- median(y_v * sd.probes)
        
        rawCN <- (2^(L_c + mean_int) / 2^median_int)
        
    } else if (weightType == "mean") {
        
        L_c <- mean(L * sd.probes) * c
        
        mean_int <- mean(y_v * sd.probes)
        
        express <- L_c + mean_int
        
        median_int <- median(y_v * sd.probes)
        
        rawCN <- (2^(L_c + mean_int) / 2^median_int)
        
    } else if (weightType == "softmax") {
        
        PsiLL <- exp( L * sd.probes)
        
        sumPsiLL <- sum(PsiLL)
        
        if (sumPsiLL == 0) { sumPsiLL <- 1 }
        
        propPsiLL <- PsiLL / sumPsiLL
        
        L_c <- as.vector(crossprod(L * sd.probes, propPsiLL)) * c
        
        mean_int <- mean(y_v * sd.probes)
        
        express <- L_c + mean_int
        
        median_int <- median(y_v * sd.probes)
        
        rawCN <- (2^(L_c + mean_int) / 2^median_int)
        
    }
    
    medianSample <- median(express)
    
    return(list(intensity = as.numeric(express), 
                    INICall = as.numeric(signal_info)[1], 
                    L_z = as.numeric(L_c), 
                    INI_sigVar = var(mean(L * sd.probes) * c),
                    rawCN = as.numeric(rawCN), 
                    z = as.numeric(c)))
    
}

