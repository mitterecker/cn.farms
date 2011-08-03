#' Corrects improper annotation names 
#' @param pkgname The pkgname.
#' @return Some data
#' @author Andreas Mitterecker
#' @noRd
correctPkgname <- function (pkgname) {
    if (pkgname == "pd.genomewideex.6") {
        pkgname <- "pd.genomewidesnp.6"
    }
    return(pkgname)
}