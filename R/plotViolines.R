#' Create a violine plot
#' 
#' This function creates a violine plot on intensity values
#' 
#' @param object An instance of 
#' \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' @param variable states which variable of assayData should be plotted.
#' @param groups Vector with the dimension of the samples for coloring.
#' @param ... Further arguments passed to the lattice graph.
#' @return Creates a violine plot.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @importFrom lattice bwplot
#' @importFrom lattice panel.violin
#' @importFrom lattice panel.bwplot
#' @examples
#' load(system.file("exampleData/normData.RData", package="cn.farms"))
#' normData <- normData[, 1:10]
#' groups <- seq(sampleNames(normData))
#' plotViolines(normData, variable="intensity", groups, xlab="Intensity values")
#' @export 
plotViolines <- function(object, variable="intensity", groups, ...) {
    x <- assayData(object)[[variable]]
    if (missing(groups)) {
        groups <- rep(1, length(sampleNames(object)))  
    } 
    myData <- data.frame(
            x=as.vector(x), 
            y=rep(sampleNames(object), each = nrow(x)))
    
    lattice::bwplot(y~x, myData, ...,
            names=sampleNames(object),
            panel = function(..., box.ratio) {
                lattice::panel.violin(
                        ..., 
                        col = "transparent",    
                        varwidth = FALSE, 
                        box.ratio = box.ratio)
                lattice::panel.bwplot(
                        ..., 
                        fill = groups, 
                        box.ratio = .1)
            })
}
