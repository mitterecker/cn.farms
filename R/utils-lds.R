## adopted from the oligoClasses package

## utilities for Large Dataset Support
##
## Summary:
##   - is.ffmatrix: test specifically for ff_matrix
##   - ldStatus: TRUE if Large Dataset Support is available
##   - ldPath: where to save ff files
##   - createFF: creates an ff object setting path appropriately
##               (leaving out of fftempdir b/c parallel processes
##               access the object very easily)

initializeBigMatrix <- function(name=basename(tempfile()), nr=0L, nc=0L, vmode="integer", initdata=NA){
    if(isPackageLoaded("ff")){
        if(prod(nr, nc) > 2^31){
            ##Need multiple matrices
            ## -- use ffdf
            ## How many samples per ff object
            S <- floor(2^31/nr - 1)
            ## How many ff objects
            L <- ceiling(nc/S)
            name <- paste(name, 1:L, sep="_")
            resultsff <- vector("list", L)
            for(i in 1:(L-1)){  ## the Lth object may have fewer than nc columns
                resultsff[[i]] <- createFF(name=name[i],
                        dim=c(nr, S),
                        vmode=vmode, initdata=initdata)
            }
            ##the Lth element
            leftOver <- nc - ((L-1)*S)
            resultsff[[L]] <- createFF(name=name[L],
                    dim=c(nr, leftOver),
                    vmode=vmode, initdata=initdata)
            results <- do.call(ffdf, resultsff)
            rm(resultsff); gc()
        } else {
            results <- createFF(name=name,
                    dim=c(nr, nc),
                    vmode=vmode, initdata=initdata)
        }
    }  else {
        init <- switch(vmode,
                integer=as.integer(initdata),
                double=as.double(initdata),
                character=as.character(initdata),
                stop("Mode ", vmode, " not implemented for regular matrices"))
        results <- matrix(init, nr, nc)
    }
    return(results)
}

createFF <- function(name, dim, vmode="double", initdata=NULL)
    ff(initdata=initdata, vmode=vmode, dim=dim, dimorder=2:1, pattern=file.path(ldPath(), basename(name)))
