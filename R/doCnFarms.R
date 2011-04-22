
#' Does the whole cn.farms process in one call
#' 
#' Works for all kind of Affymetrix SNP arrays
#' 
#' @param celfiles The celfiles which you want to process with the whole path. 
#' Either a vector or a matrix with two columns for combined analysis e.g. 500K 
#' Array.  
#' @param samplenames An optional vector with the same dimension as the number
#' of cel files 
#' @param normalization The normalization method you want to use. 
#' @return The ready cn.FARMS results.
#' @author Andreas Mitterecker
doCnFarmsSingle <- function (celfiles, samplenames, normalization) {
    
}