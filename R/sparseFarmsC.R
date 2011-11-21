#' Normalizes the data with SOR
#'  
#' @param probes The intensity matrix.
#' @param cyc Number of cycles.
#' @return Normalized Data.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @examples 
#' x <- matrix(rnorm(100, 11), 20, 5)
#' sparseFarmsC(x, 50)
sparseFarmsC <- function(probes, cyc = 5){
    
    ## probes - data matrix
    ## cyc - maximum number of cycles of EM (default 100)
    ## L - factor loadings
    
    x <- probes
    
    n <- length(x[,1])
    
    n_probes <- length(x[1,])
    
    nn <- length(x[,1])
    
    XX <- crossprod(x) / n
    
    pointer <- .Call("sparseFarmsC", x, as.integer(cyc) , XX, nn, 
            PACKAGE = "cn.farms")
    
    L <- matrix(.Call("getL", pointer, PACKAGE = "cn.farms"), 2, 3)
    
    E_SX_n <- matrix(.Call("getEss", pointer, PACKAGE = "cn.farms"), nn, 3)
    
    lapla <- matrix(.Call("getLap", pointer, PACKAGE = "cn.farms"), nn, 3)
    
    .Call("deinit", pointer, PACKAGE = "cn.farms")
    
    Lz <- L%*%t(E_SX_n)
    
    return(list(Lz = Lz))
}

