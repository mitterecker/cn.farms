#' Extracts info from the package name 
#' @param pkgname The package name according to the bioconductor annotation 
#' names.  
#' @return Additional info for save files.
#' @author Andreas Mitterecker
normAdd <- function(pkgname) {
    tmp <- unlist(strsplit(pkgname, "\\."))
    .simpleCap <- function(x) {
        s <- strsplit(x, " ")[[1]]
        paste(toupper(substring(s, 1, 1)), substring(s, 2),
                sep = "", collapse = " ")
    }
    return(.simpleCap(tmp[length(tmp)]))
}
