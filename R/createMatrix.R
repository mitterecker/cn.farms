#' Creates the needed matrix
#' @param runtype Mode how the results are saved. Possible values are ff or bm. 
#' If ff is chosen the data will not be saved automatically. 
#' With bm the results will be saved permanently.  
#' @param nrow nrow
#' @param ncol ncol
#' @param type type
#' @param bmName Identifier for ff name
#' @return a matrix 
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @importFrom oligoClasses initializeBigMatrix
#' @export 
createMatrix <- function(runtype, nrow, ncol, type = "double", bmName = "NA") {
    bmName <- paste(bmName, "_", sep = "")
    if (runtype == "ff") {
        x <- ff(vmode = type, dim = c(nrow, ncol))    
    } else if (runtype == "bm") {
        x <- initializeBigMatrix(name = bmName, nrow, ncol, vmode = type)
    } else {
        x <- matrix(NA, nrow, ncol)
    }
    return(x)
}

