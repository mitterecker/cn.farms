#' Plots a dendrogram
#' @param DivMetric The input data (see example).
#' @param colorLabels A color label with the dimension of the columns.
#' @return A dendrogram.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @examples
#' load(system.file("exampleData/normData.RData", package = "cn.farms"))
#' x <- assayData(normData)$intensity[, 1:3]
#' y <- distributionDistance(x)
#' attr(y, "Labels") <- substr(sampleNames(normData), 1, 7)
#' plotDendrogram(y)
plotDendrogram <- function(DivMetric, colorLabels) {
    
    if (missing(colorLabels)) colorLabels <- seq(DivMetric)
    
    colLab <- function(n, colors) {
        if(is.leaf(n)) {
            a <- attributes(n)
            i <<- i + 1
            attr(n, "nodePar") <-
                    c(a$nodePar, list(lab.col = colors[i], 
                                    lab.font= i%%3))
        }
        n
    }  
    i <- 0
    
    
    dhc <- as.dendrogram(hclust(DivMetric))
    colorLabels <- colorLabels[order.dendrogram(dhc)]
    dL <- dendrapply(dhc, colLab, colorLabels)
    plot(dL)
    
}
