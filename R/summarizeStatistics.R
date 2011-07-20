#' Mean or median instead of the FARMS model
#' @param probes A matrix with numeric values.
#' @param type The statistic which you want to apply.
#' @param ... Further parameters
#' @return Some data
#' @author Andreas Mitterecker
#' @export
summarizeFarmsStatistics <- function(probes, type = "median", ...) {
    if (type == "median") {
        intensity <- apply(probes, 2, median)
        L_z <- intensity - median(intensity)
    } else {
        intensity <- apply(probes, 2, mean)
        L_z <- intensity - mean(intensity)
    }
    return(list(intensity = intensity, L_z = L_z))
}