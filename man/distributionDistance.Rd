\name{distributionDistance}
\alias{distributionDistance}
\title{Computes the distribution distance...}
\usage{distributionDistance(intensityData, method=c("JSDiv", "KLDiv", "KLInf"), useSubset=T,
    subsetFraction=0.25, useQuantileReference=FALSE)
}
\description{Computes the distribution distance}
\details{Be aware that this function is implemented quite slow.}
\value{Computes the distribution distance}
\author{Djork-Arne Clevert \email{okko@clevert.de} and 
Andreas Mitterecker \email{mitterecker@bioinf.jku.at}}
\arguments{\item{intensityData}{A matrix or an AffyBatch object.}
\item{method}{The method you want to use.}
\item{useSubset}{Logical. States if only a subset should be used.}
\item{subsetFraction}{The fraction of the subset.}
\item{useQuantileReference}{Logical for a quantile reference.}
}
\examples{load(system.file("exampleData/normData.RData", package="cn.farms"))
x <- assayData(normData)$intensity[, 1:3]
y <- distributionDistance(x)
attr(y, "Labels") <- substr(sampleNames(normData), 1, 7)
plotDendrogram(y)}

