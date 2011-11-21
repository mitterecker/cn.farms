#' Plots given regions by segments 
#' 
#' A pdf in the working directory is produced.
#' 
#' @param object An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' @param segments An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' with the segments to plot
#' @param addInd States how many indices should be plotted besides the region
#' @param ylim The limits for the y axis.
#' @param variable States which variable of the assayData should be plotted.
#' @param colorVersion States different color versions.
#' @param plotLegend If a legend should be plotted or not. 
#' @param pdfname The name of the pdf file.
#' @return A graph. Normally a pdf in the current work directory. 
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export 
#' @examples
#' load(system.file("exampleData/slData.RData", package = "cn.farms"))
#' load(system.file("exampleData/testSegments.RData", package = "cn.farms"))
#' plotRegions(slData, testSegments, addInd = 10, ylim = c(-2, 2), 
#'         variable = "L_z", colorVersion = 1, plotLegend = TRUE, 
#'         pdfname = "slData.pdf")
plotRegions <- function(
        object, 
        segments, 
        addInd = NULL, 
        ylim, 
        variable, 
        colorVersion = 0, 
        plotLegend = TRUE, 
        pdfname) {
    if (missing(pdfname)) pdfname <- "plotRegions.pdf"
    if (is.null(addInd)) addInd <- 10
    if (missing(variable)) variable <- "intensity"
    if (missing(ylim)) ylim <- c(7, 15)
    
    if (!all(phenoData(object)$samples %in% colnames(assayData(segments)$cn))) {
        warning("Sample names do not fit!")
    }
    
    cat(paste(Sys.time(), " |   Writing graphics to ", getwd(), "/", pdfname, 
                    "\n", sep = ""))
    
    pdf(pdfname)
    nbrOfSegments <- nrow(featureData(segments))
    for (i in 1:nbrOfSegments) {
        cat(paste(Sys.time(), " |   Plotting graph ", i, " of ", 
                        nbrOfSegments, "\n", sep = ""))
        
        chr <- featureData(segments)$chrom[i]
        start <- featureData(segments)$start[i]
        end <- featureData(segments)$end[i]
        tmpIdx <- which(featureData(object)$chrom == chr &
                        ((featureData(object)$start <= start & 
                                featureData(object)$end >= start) | 
                            (featureData(object)$start >= start & 
                                featureData(object)$end <= end) |
                            (featureData(object)$start <= start & 
                                featureData(object)$end >= end) |
                            (featureData(object)$start <= end & 
                                featureData(object)$end >= end)))
        
        myLength <- length(tmpIdx)
        
        if (myLength == 0) {
            plot(c(1), ylab = "", xaxt = "n", col = "white", 
                    main = paste("chr", chr, ":", start, "-", end, sep=""), 
                    yaxt = "n")
            text(1, 1, "no probe in region")
        } else if (myLength > 500) {
            plot(c(1), ylab = "", xaxt = "n", col = "white", 
                    main = paste("chr", chr, ":", start, "-", end, sep = ""), 
                    yaxt = "n")
            text(1, 1, "region too large to plot")
        } else {
            indexSlExt <- max(1, (min(tmpIdx) - addInd)):
                    min((max(tmpIdx) + addInd), nrow(featureData(object)))
            ySl <- assayData(object)[[variable]][indexSlExt, ]
            xSl <- featureData(object)$start[indexSlExt]
            par(mar = c(7, 4, 4, 2) + 0.1)
            
            if (colorVersion != 0) {
                sampSegIdx <- match(sampleNames(object), 
                        sampleNames(segments))
                
                sampCn <- unlist(assayData(segments)$cn[i, ])[sampSegIdx]
                sampCn[is.na(sampCn)] <- 99
                
                ## FIXME: assure that these guys exist!!!
                myColors <- experimentData(segments)@other$cnvLabels$color
                myDescr <- experimentData(segments)@other$cnvLabels$desc      
                myNum <- experimentData(segments)@other$cnvLabels$numeric
                
                ## transform sampCn in color vector
                myColPlot <- myColors[match(sampCn, myNum)]
                
                if (colorVersion == 1) {
                    
                    ## plot cn == 2
                    myPlotInstant <- c(which(sampCn == 2), which(sampCn == 99))
                    matplot(ySl[, myPlotInstant], 
                            type = "l", 
                            pch  = "", 
                            xaxt = "n", 
                            col  = myColPlot[myPlotInstant], 
                            ylab = variable, 
                            main = paste("chr", chr, ":", start, "-", end, sep = ""),  
                            ylim = ylim)
                    abline(h = 0, col = "grey")
                    
                    ## plot others
                    myPlotInstant <- setdiff(seq(length(sampCn)), myPlotInstant)
                    if (length(myPlotInstant) != 0) {
                        matlines(ySl[, myPlotInstant], 
                                col = myColPlot[myPlotInstant])    
                    }
                    if(plotLegend) {
                        legend("top", 
                                legend = myDescr,
                                fill = myColors)
                    }
                } else if (colorVersion == 2) {
                    
                    matplot(ySl[], 
                            type = "l", 
                            pch  = "", 
                            xaxt = "n", 
                            col  = myColPlot,  
                            ylab = variable, 
                            main = paste("chr", chr, ":", start, "-", end, sep=""),  
                            ylim = ylim)
                    abline(h = 0, col = "grey")
                    if(plotLegend) {
                        legend("topright", myDescr,
                                fill = myColors)
                    }
                }
            } else {
                matplot(ySl, 
                        type = "l", 
                        pch  = "", 
                        xaxt = "n",  
                        ylab = variable, 
                        main = paste("chr", chr, ":", start, "-", end, sep=""),  
                        ylim = ylim)
                abline(h = 0, col = "grey")
            }
            
            axis(1, at = seq(dim(ySl)[1]), labels = FALSE)
            struges <- round(seq(1, dim(ySl)[1], 
                            length.out = nclass.Sturges(xSl)), 0)
            text(struges, par("usr")[3] - diff(ylim) / 20, srt = 45, adj = 1,
                    labels = xSl[struges], xpd = TRUE)
            
            
            abline(v = (addInd + 1))
            abline(v = (nrow(ySl) - addInd))
        }
    }
    dev.off()    
}
