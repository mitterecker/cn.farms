#' Function to create a density plot
#' 
#' Simple density plot. Adapted from the aroma.affymetrix package (www.aroma-project.org)
#' 
#' @param x Matrix with numeric values.  
#' @param xlim The limits for the x axis.
#' @param ylim The limits for the y axis.
#' @param col Vector with colors corresponding to the columns of the matrix.
#' @param lty The line type (see \code{\link[graphics:par]{graphics}}). 
#' @param lwd The line width, a positive number, defaulting to 1 
#' (see \code{\link[graphics:par]{graphics}}).
#' @param add If FALSE (the default) then a new plot is produced. 
#' If TRUE, density lines are added to the open graphics device.
#' @param xlab The labeling of the x axis.
#' @param ylab The labeling of the y axis.
#' @param log Logical values which states if the log2 should be taken from the data.
#' @param ... Further arguments of the plot function
#' ' @return A plot written to the graphics device.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @examples 
#' load(system.file("exampleData/slData.RData", package = "cn.farms"))
#' plotDensity(assayData(slData)$intensity)
#' @export
plotDensity <- function(x, xlim = c(0, 16), ylim, col, lty, 
        lwd, add = FALSE, xlab, ylab, log = TRUE, ...) {
    
    if (missing(ylab)) {
        ylab <- "density (integrates to one)"
    } 
    nbrOfSamples <- ncol(x)
    
    if (missing(xlab)) {
        if (log) {
            xlab <- expression(log[2](y))
        } else {
            xlab <- expression(y)
        }
    }
    
    if (missing(col)) {
        col <- seq(length = nbrOfSamples)
    } else {
        col <- rep(col, length.out = nbrOfSamples)
    }

    if (missing(lty)) {
        lty <- NULL
    } else {
        lty <- rep(lty, length.out = nbrOfSamples)
    }
    
    if (missing(lwd)) {
        lwd <- NULL
    } else {
        lwd <- rep(lwd, length.out = nbrOfSamples)
    }
    
    ds <- list()
    xlimDef <- c(NA, NA)
    ylimDef <- c(0, NA)
    for(kk in 1:nbrOfSamples) {
        xx <- x[, kk]
        xx <- xx[is.finite(xx)]
        suppressWarnings({
                    d <- density(xx, ...)
                })
        ds[[kk]] <- d
        xlimDef <- range(c(xlimDef, range(d$x, na.rm = TRUE)), na.rm = TRUE)
        ylimDef <- range(c(ylimDef, range(d$y, na.rm = TRUE)), na.rm = TRUE)
    }
    if (missing(xlim))
        xlim <- xlimDef
    if (missing(ylim))
        ylim <- ylimDef
    if (add == FALSE) {
        suppressWarnings({
                    plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, 
                            ylab = ylab, ...)
                })
    }
    for(kk in 1:nbrOfSamples) {
        suppressWarnings({
                    lines(ds[[kk]], col = col[kk], lty = lty[kk], 
                            lwd = lwd[kk], ...)
                })
    }
}
