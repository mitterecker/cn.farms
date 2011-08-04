\name{flcStd}
\alias{flcStd}
\title{Does a fragment length correction on intensities}
\usage{
  flcStd(y, fragmentLengths, targetFcn = NULL, subsetToFit
  = NULL, runtype = "ff", cores = 1, saveFile = "flc", ...)
}
\arguments{
  \item{y}{y}

  \item{fragmentLengths}{fragmentLengths}

  \item{targetFcn}{targetFcn}

  \item{subsetToFit}{subsetToFit}

  \item{runtype}{Mode how the results are saved. Possible
  values are ff or bm. If ff is chosen the data will not be
  saved automatically. With bm the results will be saved
  permanently.}

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

