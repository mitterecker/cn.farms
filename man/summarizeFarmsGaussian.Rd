\name{summarizeFarmsGaussian}
\alias{summarizeFarmsGaussian}
\title{Summarization Gaussian approach}
\usage{
  summarizeFarmsGaussian(probes, weight = 0.15, mu = 0, cyc
  = 10, tol = 1e-04, weightType = "mean", init = 0.6,
  correction = 0, minNoise = 0.35, centering = "median",
  refIdx)
}
\arguments{
  \item{probes}{A matrix with numeric values.}

  \item{weight}{Hyperparameter value in the range of [0,1]
  which determines the influence of the prior.}

  \item{mu}{Hyperparameter value which allows to quantify
  different aspects of potential prior knowledge. Values
  near zero assumes that most genes do not contain a
  signal, and introduces a bias for loading matrix elements
  near zero. Default value is 0.}

  \item{cyc}{Number of cycles for the EM algorithm.}

  \item{tol}{States the termination tolerance. Default is
  0.00001.}

  \item{weightType}{Flag, that is used to summarize the
  loading matrix. The default value is set to mean.}

  \item{init}{Parameter for estimation.}

  \item{correction}{Value that indicates whether the
  covariance matrix should be corrected for negative
  eigenvalues which might emerge from the non-negative
  correlation constraints or not. Default = O (means that
  no correction is done), 1 (minimal noise (0.0001) is
  added to the diagonal elements of the covariance matrix
  to force positive definiteness), 2 (Maximum Likelihood
  solution to compute the nearest positive definite matrix
  under the given non-negative correlation constraints of
  the covariance matrix)}

  \item{minNoise}{States the minimal noise. Default is
  0.35.}

  \item{centering}{States how the data is centered. Default
  is median.}

  \item{refIdx}{index or indices which are used for
  computation of the centering}
}
\value{
  A list containing the results of the run.
}
\description{
  This function runs the FARMS algorithm.
}
\examples{
x <- matrix(rnorm(100, 11), 20, 5)
summarizeFarmsGaussian(x)
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

