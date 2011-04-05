\name{sparseFarmsC}
\alias{sparseFarmsC}
\title{Normalizes the data with SOR...}
\usage{sparseFarmsC(probes, cyc=5)
}
\description{Normalizes the data with SOR}
\value{Normalized Data.}
\author{Djork-Arne Clevert \email{okko@clevert.de} and 
Andreas Mitterecker \email{mitterecker@bioinf.jku.at}}
\arguments{\item{probes}{The intensity matrix.}
\item{cyc}{Number of cycles.}
}
\examples{x <- matrix(rnorm(100, 11), 20, 5)
sparseFarmsC(x, 50)}

