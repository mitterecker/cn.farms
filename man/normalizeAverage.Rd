\name{normalizeAverage}
\alias{normalizeAverage}
\title{Scales the range of the non-polymorphic data to the range of a given array.}
\usage{normalizeAverage(x, baselineArray, avg=median, targetAvg=2200, ...)
}
\description{Scales the range of the non-polymorphic data to the range of a given array.}
\value{Normalized non-polymorphic data.}
\author{Djork-Arne Clevert \email{okko@clevert.de} and 
Andreas Mitterecker \email{mitterecker@bioinf.jku.at}}
\arguments{\item{x}{Data matrix}
\item{baselineArray}{Choose the baseline channel array.}
\item{avg}{The function for averaging.}
\item{targetAvg}{Value to which the array should be averaged.}
\item{...}{Further optional parameters.}
}
\examples{x <- matrix(rnorm(100, 11), 20, 5)
normalizeAverage(x, x[, 1])}

