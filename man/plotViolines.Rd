\name{plotViolines}
\alias{plotViolines}
\title{Create a violine plot...}
\usage{plotViolines(object, variable="intensity", groups, ...)
}
\description{Create a violine plot}
\details{This function creates a violine plot on intensity values}
\value{Creates a violine plot.}
\author{Djork-Arne Clevert \email{okko@clevert.de} and 
Andreas Mitterecker \email{mitterecker@bioinf.jku.at}}
\arguments{\item{object}{An instance of 
\code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}}
\item{variable}{states which variable of assayData should be plotted.}
\item{groups}{Vector with the dimension of the samples for coloring.}
\item{...}{Further arguments passed to the lattice graph.}
}
\examples{load(system.file("exampleData/normData.RData", package="cn.farms"))
normData <- normData[, 1:10]
groups <- seq(sampleNames(normData))
plotViolines(normData, variable="intensity", groups, xlab="Intensity values")}

