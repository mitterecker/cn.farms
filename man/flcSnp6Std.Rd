\name{flcSnp6Std}
\alias{flcSnp6Std}
\title{Does a fragment length correction on intensities}
\usage{
  flcSnp6Std(y, fragmentLengths, targetFcn = NULL,
  subsetToFit = NULL, runtype = "ff", cores = 1, saveFile =
  "flc", ...)
}
\arguments{
  \item{y}{y}

  \item{fragmentLengths}{fragmentLengths}

  \item{targetFcn}{targetFcn}

  \item{subsetToFit}{subsetToFit}

  \item{runtype}{runtype}

  \item{cores}{cores}

  \item{saveFile}{Name of the file to save.}

  \item{...}{...}
}
\value{
  data frame
}
\description{
  Does a fragment length correction on intensities
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

