\name{mlSummarization}
\alias{mlSummarization}
\title{Does summarization}
\usage{
  mlSummarization(object, windowMethod, windowParam,
    summaryMethod, summaryParam,
    callParam = list(runtype = "ff"), returnValues,
    saveFile = "mlData")
}
\arguments{
  \item{object}{an instance of
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}}

  \item{windowMethod}{Method for combination of
  neighbouring SNPs. Possible values are Std and Bps.}

  \item{windowParam}{further parameters as the window size}

  \item{summaryMethod}{allowed versions for the
  summarization step are: Gaussian, Variational, Exact.
  Default is Variational.}

  \item{summaryParam}{summaryParam}

  \item{callParam}{callParam}

  \item{returnValues}{List with return values.}

  \item{saveFile}{Name of the file to save. For possible
  values see summaryMethod.}
}
\value{
  Some data
}
\description{
  Does summarization
}
\examples{
load(system.file("exampleData/slData.RData", package = "cn.farms"))
windowMethod <- "std"
windowParam <- list()
windowParam$windowSize <- 5
windowParam$overlap <- TRUE
summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(20)
mlData <- mlSummarization(slData, windowMethod, windowParam,
        summaryMethod, summaryParam)
assayData(mlData)
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

