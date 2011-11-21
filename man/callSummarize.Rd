\name{callSummarize}
\alias{callSummarize}
\title{Defines which variables should be written back when calling a cn.farms run}
\usage{
  callSummarize(object, psInfo, summaryMethod,
    summaryParam, batchList = NULL, cores = 1,
    runtype = "ff", returnValues, saveFile = "summData")
}
\arguments{
  \item{object}{an matrix with normalized intensity
  values.}

  \item{psInfo}{a data frame stating the physical
  position.}

  \item{summaryMethod}{the summarization method.}

  \item{summaryParam}{a list with the parameters of the
  summarization method.}

  \item{batchList}{batchList}

  \item{cores}{cores}

  \item{runtype}{mode how the results are saved. Possible
  values are ff or bm. If ff is chosen the data will not be
  saved automatically. With bm the results will be saved
  permanently.}

  \item{returnValues}{list with return values. For possible
  values see summaryMethod.}

  \item{saveFile}{name of the file to save.}
}
\value{
  Results of FARMS run with specified parameters - exact
  FARMS version
}
\description{
  Defines which variables should be written back when
  calling a cn.farms run
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

