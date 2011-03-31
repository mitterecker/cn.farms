#' Normalize array with non-polymorphic data
#' @param x Data matrix
#' @param baseline Choose the baseline channel.
#' @param avg The function for averaging.
#' @param targetAvg Value to which the array should be averaged.
#' @param ... Further optional parameters.
#' @return Normalized non-polymorphic data.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export 
#' @examples
#' x <- matrix(rnorm(100, 11), 20, 5)
#' normalizeAverage(x)
normalizeAverage <- function(x, baseline=1, avg=median, targetAvg=2200, ...) {
    
    # Estimate the scale for each channel
    scale <- apply(x, MARGIN=2, FUN=avg, ...)
    
    scale1 <- scale[baseline]
    
    # Standardize so that the 'baseline' channel has scale one.
    scale <- scale / scale1
    
    # Rescale to target averages?
    if (!is.null(targetAvg)) {
        
        rho <- (scale1 / targetAvg)
        
        scale <- rho * scale
        
    }
    
    # Rescale so that all channels have the same scale
    for (cc in 1:ncol(x)) {
        
        x[,cc] <- x[,cc] / scale[cc]
        
    }
    
    return(x)
}
