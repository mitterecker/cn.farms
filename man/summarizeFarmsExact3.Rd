\name{summarizeFarmsExact3}
\alias{summarizeFarmsExact3}
\title{Summarization Laplacian approach with exact computation}
\usage{
  summarizeFarmsExact3(probes, mu = 1, weight = 100,
    weightSignal = 1, weightZ = 30, weightProbes = TRUE,
    updateSignal = FALSE, cyc = c(10, 10), tol = 1e-05,
    weightType = "mean", centering = "median",
    rescale = FALSE, backscaleComputation = FALSE,
    maxIntensity = TRUE, refIdx, ...)
}
\arguments{
  \item{probes}{A matrix with numeric values.}

  \item{mu}{Hyperparameter value which allows to quantify
  different aspects of potential prior knowledge. Values
  near zero assumes that most positions do not contain a
  signal, and introduces a bias for loading matrix elements
  near zero. Default value is 0 and it's recommended not to
  change it.}

  \item{weight}{Hyperparameter value which determines the
  influence of the Gaussian prior of the loadings}

  \item{weightSignal}{Hyperparameter value on the signal.}

  \item{weightZ}{Hyperparameter value which determines how
  strong the Laplace prior of the factor should be at 0.
  Users should be aware, that a change of weightZ in
  comparison to the default parameter might also entail a
  need to change other parameters. Unexperienced users
  should not change weightZ.}

  \item{weightProbes}{Parameter (TRUE/FALSE), that
  determines, if the number of probes should additionally
  be considered in weight. If TRUE, weight will be
  modified.}

  \item{updateSignal}{updateSignal.}

  \item{cyc}{Number of cycles. If the length is two, it is
  assumed, that a minimum and a maximum number of cycles is
  given. If the length is one, the value is interpreted as
  the exact number of cycles to be executed (minimum ==
  maximum).}

  \item{tol}{States the termination tolerance if
  cyc[1]!=cyc[2]. Default is 0.00001.}

  \item{weightType}{Flag, that is used to summarize the
  probes of a sample.}

  \item{centering}{States how the data should be centered
  ("mean", "median"). Default is median.}

  \item{rescale}{Parameter (TRUE/FALSE), that determines,
  if moments in exact Laplace FARMS are rescaled in each
  iteration. Default is FALSE.}

  \item{backscaleComputation}{Parameter (TRUE/FALSE), that
  determines if the moments of hidden variables should be
  reestimated after rescaling the parameters.}

  \item{maxIntensity}{Parameter (TRUE/FALSE), that
  determines if the expectation value (=FALSE) or the
  maximum value (=TRUE) of p(z|x_i) should be used for an
  estimation of the hidden varaible.}

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
  algorithm.
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

