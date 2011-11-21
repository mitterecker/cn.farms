\name{plotDensity}
\alias{plotDensity}
\title{Function to create a density plot}
\usage{
  plotDensity(x, xlim = c(0, 16), ylim, col, lty, lwd,
    add = FALSE, xlab, ylab, log = TRUE, ...)
}
\arguments{
  \item{x}{Matrix with numeric values.}

  \item{xlim}{The limits for the x axis.}

  \item{ylim}{The limits for the y axis.}

  \item{col}{Vector with colors corresponding to the
  columns of the matrix.}

  \item{lty}{The line type (see
  \code{\link[graphics:par]{graphics}}).}

  \item{lwd}{The line width, a positive number, defaulting
  to 1 (see \code{\link[graphics:par]{graphics}}).}

  \item{add}{If FALSE (the default) then a new plot is
  produced. If TRUE, density lines are added to the open
  graphics device.}

  \item{xlab}{The labeling of the x axis.}

  \item{ylab}{The labeling of the y axis.}

  \item{log}{Logical values which states if the log2 should
  be taken from the data.}

  \item{...}{Further arguments of the plot function '}
}
\value{
  A plot written to the graphics device.
}
\description{
  Simple density plot. Adapted from the aroma.affymetrix
  package (www.aroma-project.org)
}
\examples{
load(system.file("exampleData/slData.RData", package = "cn.farms"))
plotDensity(assayData(slData)$intensity)
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

