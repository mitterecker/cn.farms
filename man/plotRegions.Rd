\name{plotRegions}
\alias{plotRegions}
\title{Plots given regions by segments}
\usage{
  plotRegions(object, segments, addInd = NULL, ylim,
    variable, colorVersion = 0, plotLegend = TRUE, pdfname)
}
\arguments{
  \item{object}{An instance of
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}}

  \item{segments}{An instance of
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
  with the segments to plot}

  \item{addInd}{States how many indices should be plotted
  besides the region}

  \item{ylim}{The limits for the y axis.}

  \item{variable}{States which variable of the assayData
  should be plotted.}

  \item{colorVersion}{States different color versions.}

  \item{plotLegend}{If a legend should be plotted or not.}

  \item{pdfname}{The name of the pdf file.}
}
\value{
  A graph. Normally a pdf in the current work directory.
}
\description{
  A pdf in the working directory is produced.
}
\examples{
load(system.file("exampleData/slData.RData", package = "cn.farms"))
load(system.file("exampleData/testSegments.RData", package = "cn.farms"))
plotRegions(slData, testSegments, addInd = 10, ylim = c(-2, 2),
        variable = "L_z", colorVersion = 1, plotLegend = TRUE,
        pdfname = "slData.pdf")
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

