#' Creates a smooth scatter plot
#' @param object An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}.
#' @param variable States which variable of the assayData should be plotted.
#' @param chrom The chromosome you want to plot.
#' @param start The physical start position.
#' @param end The physical end position.
#' @param ylim The limits for the y axis.
#' @param pdfname The name of the pdf file.
#' @param ... Further arguments passed to smoothScatter function.
#' @return A graph.
#' @author Andreas Mitterecker
#' @export 
#' @examples
#' load(system.file("exampleData/slData.RData", package="cn.farms"))
#' plotSmoothScatter(slData[, 1:3], chrom="23")
plotSmoothScatter <- function(
        object, 
        variable,
        chrom, 
        start, 
        end, 
        ylim, 
        pdfname, 
        ...) {

    if (missing(variable)) variable <- "L_z"
    if (missing(ylim)) ylim <- c(-1, 1)
    if (missing(pdfname)) pdfname <- "envSumHapMap.pdf"
    if (missing(chrom)) {
        warning("No chromosome provided: assuming chrom 23 (X)")
        chrom <- 23
    }
    
    chrom <- oligoClasses::chromosome2integer(chrom)
    
    if (missing(start) | missing(end)) {
        tmpIdx <- which(featureData(object)@data$chrom == chrom)
    } else {
        tmpIdx <- which(featureData(object)$chrom == chrom &
                        ((featureData(object)$start <= start & 
                                featureData(object)$end >= start) | 
                            (featureData(object)$start >= start & 
                                featureData(object)$end <= end) |
                            (featureData(object)$start <= start & 
                                featureData(object)$end >= end) |
                            (featureData(object)$start <= end & 
                                featureData(object)$end >= end)))
    }
    
    if (length(tmpIdx) == 0) {
        stop("No data to plot!")
    }
    
    cat(paste(Sys.time(), " |   Writing graphics to ", getwd(), "/", pdfname, 
                    "\n", sep = ""))
    
    pdf(pdfname)
    x <- object[tmpIdx, ]

    xSl <- featureData(x)$start
    nbrOfSamples <- ncol(assayData(x)[[variable]])
    for (i in seq(nbrOfSamples)) {
        cat(paste(Sys.time(), " |   Plotting graph ", i, " of ", 
                        nbrOfSamples, "\n", sep = ""))
        
        ySl <- assayData(x)[[variable]][, i]
        smoothScatter(x=xSl, y = ySl, ylim=ylim, ylab = "L_z", xlab = "bp", 
                nrpoints = 1000, ...) 
        abline(h=0)
        loess.out <- loess(ySl ~ xSl)
        points(loess.out$x, loess.out$fitted, col = "red", pch = ".")      
    }
    dev.off()
}


