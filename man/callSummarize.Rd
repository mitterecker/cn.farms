\name{callSummarize}
\alias{callSummarize}
\title{Defines which variables should be written back}
\usage{
  callSummarize(object, psInfo, summaryMethod,
  summaryParam, batchList = NULL, cores = 1, runtype =
  "ff", returnValues, saveFile = "summData")
}
\arguments{
  \item{object}{Normalized intensity values}

  \item{psInfo}{Physical position}

  \item{summaryMethod}{summaryMethod}

  \item{summaryParam}{summaryParam}

  \item{batchList}{batchList}

  \item{cores}{cores}

  \item{runtype}{Mode how the results are saved. Possible
  values are ff or bm. If ff is chosen the data will not be
  saved automatically. With bm the results will be saved
  permanently.}

  \item{returnValues}{List with return values. For possible
  values see summaryMethod.}

  \item{saveFile}{Name of the file to save.}
}
\value{
  Results of FARMS run with specified parameters - exact
  FARMS version
}
\description{
  Defines which variables should be written back
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

