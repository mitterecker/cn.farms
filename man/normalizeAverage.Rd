\name{normalizeAverage}
\alias{normalizeAverage}
\title{Normalize array with non-polymorphic data...}
\usage{normalizeAverage(x, baseline=1, avg=median, targetAvg=2200, ...)
}
\description{Normalize array with non-polymorphic data}
\value{Normalized non-polymorphic data.}
\author{Djork-Arne Clevert \email{okko@clevert.de} and 
Andreas Mitterecker \email{mitterecker@bioinf.jku.at}}
\arguments{\item{x}{Data matrix}
\item{baseline}{Choose the baseline channel.}
\item{avg}{The function for averaging.}
\item{targetAvg}{Value to which the array should be averaged.}
\item{...}{Further optional parameters.}
}
\examples{x <- matrix(rnorm(100, 11), 20, 5)
normalizeAverage(x)}

