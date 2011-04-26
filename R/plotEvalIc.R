#' Creates a plot with known regions and a numeric vector 
#' @param object  an instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' @param segments A data.frame with known regions.
#' @param chrom the chromosome.
#' @param variable The numeric vector which should be plotted.
#' @param ylim the limits of the y axis. 
#' @param stripCol color of points.
#' @param regionCol color of regions. 
#' @param pointSize size of the points. 
#' @param pointType type of the points.
#' @param bandwidth for the color of the plot.
#' @param nbin number of bins for the coloring. 
#' @return Some data
#' @author Andreas Mitterecker
#' @export
#' @examples 
#' load(system.file("exampleData/slData.RData", package="cn.farms"))
#' load(system.file("exampleData/testSegments.RData", package="cn.farms"))
#' plotEvalIc(slData, featureData(testSegments)@@data, 
#'      variable=assayData(slData)$L_z[, 1],  23)
plotEvalIc <- function(
        object,
        segments, 
        chrom,
        variable, 
        ylim,
        stripCol="lightgray", 
        regionCol=rgb(130, 0, 139, max=255),
        pointSize=0.75, 
        pointType=4, 
        bandwidth=c(0.01, 1000), 
        nbin=100) { 

    if (missing(chrom)) {
        stop("Please state for which chromosome you want to create the plot!")
    }
    
    if (missing(variable)) {
        
        if ("L_z" %in% names(assayData(object))) {
            message("Assuming I/NI call as plot variable.")
            variable <- assayData(object)$INICall
        } else {
            stop("Please provide function with a variable (vector) to plot!")    
        }
        variable <- assayData(object)$L_z[, 3]
        #variable <- abs(assayData(object)$L_z[, 3])
        #variable <- variable / max(variable)
    }

    if (missing(ylim)) {
        ylim <- c(min(variable), max(variable))        
    }
    
    chrom <- oligoClasses::chromosome2integer(chrom)
    phInf <- featureData(object)@data
    chrIdx <- which(phInf$chrom == chrom)
    
    xlim <- c(
            min(phInf[chrIdx, 2] / (1E+6)),
            max(phInf[chrIdx, 2] / (1E+6)))
    
    cnvrPosChr <- segments[which(segments[, 1] == chrom), c(2, 3)]
    Chr1_regions_conrad <- unlist(apply(cnvrPosChr, 1, 
                    function(x){ x[1]:x[2] }))
    tmp01 <- phInf[chrIdx, 2] %in% Chr1_regions_conrad
    tmp02 <- phInf[chrIdx, 3] %in% Chr1_regions_conrad
    cnRegions <- tmp01 | tmp02
    labels <- chrIdx
    labels[] <- 0
    labels[cnRegions] <- 1
    
    yOutCn <- variable[which(labels == 0)]
    yInCn <- variable[which(labels == 1)]
    
    ## create plot
    plot(-1, -1, 
            xlab=paste("Position on chromosome", 
                    chrom, "[mb]", sep=" "), 
            ylab="I/NI call", 
            ylim=ylim,
            cex.lab=1.5,
            cex.axis=1.5, 
            xlim=xlim)
    
    ## positions of true CNVRs
    rect(cnvrPosChr[, 1] / (1E+6), 
            ybottom=0, 
            cnvrPosChr[, 2] / (1E+6), 
            ytop=1, 
            col=stripCol,
            border=stripCol,
            lwd=0.5)
    
    if (length(which(labels == 1)) == 0) {
        x_all <- phInf[chrIdx, 2] / (1E+6)
        points( y=yOutCn, 
                x=x_all,
                cex=.1, 
                pch=8,
                col=densCols(
                        y=yOutCn,
                        x=x_all, 
                        nbin=nbin, 
                        bandwidth=bandwidth))
    } else {
        x_regions <- phInf[chrIdx[which(labels == 1)], 2] / (1E+6)
        x_all <- phInf[chrIdx[-which(labels == 1)], 2] / (1E+6)
        
        points( y=yOutCn, 
                x=x_all,
                cex=.1, 
                pch=8,
                col=densCols(
                        y=yOutCn,
                        x=x_all, 
                        nbin=nbin, 
                        bandwidth=bandwidth))
        
        points( y=yInCn,
                x=x_regions,
                col=regionCol,
                cex=pointSize, 
                pch=pointType)
    }
    
    mtext("cn.FARMS", side=3, adj=c(0, 0),  at=1, cex=1)
}