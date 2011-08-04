\name{plotDendrogram}
\alias{plotDendrogram}
\title{Plots a dendrogram}
\usage{
  plotDendrogram(DivMetric, colorLabels)
}
\arguments{
  \item{DivMetric}{The input data (see example).}

  \item{colorLabels}{A color label with the dimension of
  the columns.}
}
\value{
  A dendrogram.
}
\description{
  Plots a dendrogram
}
\examples{
load(system.file("exampleData/normData.RData", package="cn.farms"))
x <- assayData(normData)$intensity[, 1:3]
y <- distributionDistance(x)
attr(y, "Labels") <- substr(sampleNames(normData), 1, 7)
plotDendrogram(y)
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

