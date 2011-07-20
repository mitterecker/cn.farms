## FIXME: delete tol?

#' Summarization variational Laplacian approach
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
#' @param spuriousCorrelation Numeric value for suppression of spurious 
#' correlation.
#' @param minNoise States the minimal noise. Default is 0.35.
#' @param centering States how the data is centered. Default is median.
#' @return A list containing the results of the run. 
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @examples 
#' x <- matrix(rnorm(100, 11), 20, 5)
#' summarizeFarmsVariational(x)
summarizeFarmsVariational <- function(
        probes, 
        weight = 0.15, 
        mu = 0, 
        cyc = 10, 
        weightType = "median", 
        init = 0.6, 
        correction = 0, 
        minNoise = 0.35, 
        spuriousCorrelation = 0.3,
        centering = "median") {
    
    a_old <- 0.5
    
    n_array <- ncol(probes)
    
    n_probes <- nrow(probes)
    
    eps <- 1E-5
    
    if (centering == "median") { mean.probes <- Biobase::rowMedians(probes) }
    
    if (centering == "mean") { mean.probes <- rowMeans(probes) } 
    
    centered.probes <- probes - mean.probes
    
    sd.probes <- sqrt(diag(crossprod(t(centered.probes))) / n_array)
    
    if(0 %in% sd.probes) {
        
        index <- which(sd.probes == 0)
        
        sd.probes[index] <- eps
        
        probes <- probes / sd.probes
        
        x <- t(probes)
        
        if (centering=="median") { y_v <- apply(x, 2, median) }
        
        if (centering=="mean") { y_v <- colMeans(x) }
        
        xmean <- matrix(y_v, n_array, n_probes, byrow = TRUE)
        
        X <- x - xmean
        
        XX <- crossprod(X, X) / n_array
        
        diag(XX)[index] <- 1
        
    } else {
        
        probes <- probes / sd.probes
        
        x <- t(probes)
        
        if (centering == "median") { y_v <- apply(x, 2, median) }
        
        if (centering == "mean") { y_v <- colMeans(x) }
        
        xmean <- matrix(y_v, n_array, n_probes, byrow = TRUE)
        
        X <- x - xmean
        
        XX <- crossprod(X, X) / n_array
        
    }
    
    XX <- (XX + t(XX)) / 2 
    
    XX[which(XX < 0.1)] <- 0
    
    minEigenValues <- -1
    
    if(correction >= 1) {
        
        while(minEigenValues < 0) {
            
            eigen_XX <- eigen(XX)
            
            eigenValues_XX <- eigen_XX$values
            
            eigenVectors_XX <- eigen_XX$vectors
            
            minEigenValues <- min(eigenValues_XX)
            
            if(correction < 2){
                
                if(minEigenValues < minNoise){
                    
                    diag(XX) <- diag(XX) + (minNoise - minEigenValues)
                    
                }
                
            } else {
                
                if(minEigenValues < minNoise) {
                    
                    eigenValues_XX[which(eigenValues_XX < minNoise)] <- minNoise
                    
                    XX <- eigenVectors_XX %*% diag(eigenValues_XX) %*% t(eigenVectors_XX)
                    
                }
                
            }
            
        }
        
    }
    
    diagXX <- diag(XX)
    
    XX_tmp <- XX
    
    diag(XX_tmp) <- 0
    
#    L <- init * sqrt(apply(XX_tmp, 1, max))
#
#    Ph <- 1 - L^2
#
#    Ph[which(Ph < eps)] <- eps
    
    L <- sqrt(init * diagXX)
    
    Ph <- diagXX - L^2
    
    alpha <- weight
    
    bbeta <- mu * alpha
    
    PsiL <- (1/Ph) * L
    
    a <- as.vector(1 + crossprod(L, PsiL))
    
    bar <- PsiL / a
    
    mu_ZX <- X %*% bar
    
    lapla <- 1 / sqrt(mu_ZX^2)
    
    for (i in 1:cyc) {
        
        ## E Step
        
        PsiL <- (1 / Ph) * L
        
        a <- 1 / as.vector(as.vector(lapla) + crossprod(L, PsiL))
        
        mu_ZX <- X %*% PsiL * a
        
        EZZ <- mu_ZX^2 + a
        
        ## M Step
        
        sumXMU <- 1 / n_array * crossprod(X, mu_ZX)
        
        L <- (sumXMU + Ph * bbeta) / (mean(EZZ) + Ph * alpha)
        
        L[which(L < eps)] <- 0
        
        Ph <- diagXX - L * sumXMU + Ph * alpha * L * (mu - L) 
        
        Ph[which(Ph < eps)] <- eps
        
        lapla <- 1 / sqrt(EZZ)
        
        if (spuriousCorrelation != 0) { 
            
            lapla[which(lapla < spuriousCorrelation)] <- spuriousCorrelation
            
        }
               
        a_old <- a
        
    }
    
    c <- mu_ZX  
    
    var_scale <- sd(c) * (1 - 1 / n_array)
    
    if(var_scale == 0){
        
        var_z_scale <- eps
        
    } else {
        
        var_z_scale <- var_scale
        
    }
    
    c <- c / as.vector(var_z_scale)
    
    L <- L * as.vector(var_z_scale)
    
    PsiL <- (1 / Ph) * L
    
    a <- as.vector(1 + crossprod(L, PsiL))    
    
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
        
        PsiLL <- (1 / Ph) * L
        
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
    
    results <- list()
    
    ## recomputed intensity values
    results$intensity <- as.numeric(express)
    
    ## L_z
    results$L_z <- as.numeric(L_c)
    
    ## L
    results$L <- L
    
    ## z
    results$z <- c
    
    ## rawCN
    results$rawCN <- as.numeric(rawCN) 
    
    ## original I/NI call
    results$SNR <- 1 / a
    
    ## original I/NI call
    results$INICall <- 1 / a
       
    ## laplacian IC from farms package
    results$IC <- 1 / as.vector(1 + crossprod(L, PsiL) / as.vector(lapla))
    
    ## laINI
    results$laINI <- log(1 + as.vector(crossprod(L, PsiL) / as.vector(lapla)))
    
    ## lapla
    results$lapla <- lapla
    
    ## other IC 
    results$SNR_a <- log(1 + (1 - 1 / a) * median(sd.probes))
    results$SNR_b <- log(1 + (1 - 1 / a) * mean(sd.probes))
    results$SNR_c <- var(L_c)
    results$a_lapla <- as.vector(as.vector(lapla) + crossprod(L, PsiL))
    results$SNR_d <- log(1 + median(sd.probes) * c^2 / (1 / results$a_lapla + c^2))
    results$a_lapla1 <- as.vector(1 + crossprod(L, PsiL) / as.vector(lapla))
    results$IC1   <- log(lapla / results$SNR)
    results$IC3   <- 0.5 * log(1 + lapla / results$SNR)
        
    return(results)
    
}