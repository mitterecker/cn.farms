\name{summarizeFarmsExact}
\alias{summarizeFarmsExact}
\title{Summarization Laplacian approach with exact computation}
\usage{
  summarizeFarmsExact(probes, mu = 0, weight = 0.5, weightZ
  = 1, weightProbes = TRUE, cyc = c(10, 10), tol = 1e-05,
  weightType = "mean", centering = "median", rescale =
  FALSE, backscaleComputation = FALSE, maxIntensity = TRUE,
  refIdx, ...)
}
\arguments{
  \item{probes}{A matrix with numeric values.}

  \item{mu}{Hyperparameter value which allows to quantify
  different aspects of potential prior knowledge. Values
  near zero assumes that most positions do not contain a
  signal, and introduces a bias for loading matrix elements
  near zero. Default value is 0.}

  \item{weight}{Hyperparameter value which determines the
  influence of the Gaussian prior of the loadings}

  \item{weightZ}{Hyperparameter value which determines how
  strong the Laplace prior of the factor should be at 0.}

  \item{weightProbes}{States if the probes should be
  weighted.}

  \item{cyc}{Number of cycles. If the length is two, it is
  assumed, that a minimum and a maximum number of cycles is
  given. If the length is one, the value is interpreted as
  the exact number of cycles to be executed (minimum ==
  maximum).}

  \item{tol}{States the termination tolerance if
  cyc[1]!=cyc[2]. Default is 0.00001.}

  \item{weightType}{Flag, that is used to summarize the
  loading matrix.}

  \item{centering}{States how the data is centered. Default
  is median.}

  \item{rescale}{Rescales the Moments.}

  \item{backscaleComputation}{New estimation of z values
  after backscaling.}

  \item{maxIntensity}{Use of the mode values for building
  expression values, if set to TRUE.}

  \item{refIdx}{index or indices which are used for
  computation of the centering}

  \item{...}{Further parameters for expert users.}
}
\value{
  A list including: the found parameters: lambda0, lambda1,
  Psi

  the estimated factors: z (expectation), maxZ (maximum)

  p: log-likelihood of the data given the found lambda0,
  lambda1, Psi (not the posterior likelihood that is
  optimized)

  varzx: variances of the hidden variables given the data

  KL: Kullback Leibler divergences between between
  posterior and prior distribution of the hidden variables

  IC: Information Content considering the hidden variables
  and data

  ICtransform: transformed Information Content

  Case: Case for computation of a sample point
  (non-exception, special exception)

  L1median: Median of the lambda vector components

  intensity: back-computed summarized probeset values with
  mean correction

  L_z: back-computed summarized probeset values without
  mean correction

  rawCN: transformed values of L_z

  SNR: some additional signal to noise ratio value
}
\description{
  This function implements an exact Laplace FARMS
  algorithm. Users should be aware, that a change of weight
  in comparison to the default parameter might also entail
  a need to change of eps1 and eps2. Unexperienced users
  should not change weightZ, since a change in weightZ is
  also connected to weight, eps1 and eps2.
}
\examples{
x <- matrix(rnorm(100, 11), 20, 5)
summarizeFarmsExact(x)
}
\author{
  Andreas Mayr \email{mayr@bioinf.jku.at} and Djork-Arne
  Clevert \email{okko@clevert.de} and Andreas Mitterecker
  \email{mitterecker@bioinf.jku.at}
}

