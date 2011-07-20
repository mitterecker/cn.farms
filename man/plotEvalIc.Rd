\name{plotEvalIc}
\alias{plotEvalIc}
\title{Creates a plot with known regions and a numeric vector...}
\usage{plotEvalIc(object, segments, chrom, variable, ylim, ylab="CN indicator",
    stripCol="lightgray", regionCol=rgb(130, 0, 139, max = 255),
    pointSize=0.75, pointType=4, bandwidth=c(0.01, 1000), nbin=100)
}
\description{Creates a plot with known regions and a numeric vector}
\value{Some data}
\author{Andreas Mitterecker}
\arguments{\item{object}{an instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}}
\item{segments}{A data.frame with known regions.}
\item{chrom}{the chromosome.}
\item{variable}{The numeric vector which should be plotted.}
\item{ylim}{the limits of the y axis.}
\item{ylab}{the ylab from function par.}
\item{stripCol}{color of points.}
\item{regionCol}{color of regions.}
\item{pointSize}{size of the points.}
\item{pointType}{type of the points.}
\item{bandwidth}{for the color of the plot.}
\item{nbin}{number of bins for the coloring.}
}
\examples{load(system.file("exampleData/slData.RData", package="cn.farms"))
load(system.file("exampleData/testSegments.RData", package="cn.farms"))
plotEvalIc(slData, featureData(testSegments)@data, 
variable=assayData(slData)$L_z[, 1], 23)}

