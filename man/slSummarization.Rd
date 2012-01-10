\name{slSummarization}
\alias{slSummarization}
\title{Method for computation of the single-locus summarization}
\usage{
  slSummarization(object, summaryMethod = "Exact",
    summaryParam = list(),
    callParam = list(runtype = "ff", cores = 1),
    summaryWindow = c("std", "fragment"), returnValues,
    saveFile = "slData")
}
\arguments{
  \item{object}{An instance of
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}}

  \item{summaryMethod}{allowed versions for the
  summarization step are: Gaussian,Variational, Exact.
  Default is Variational.}

  \item{summaryParam}{The parameters for the summaryMethod.
  Further information can be obtained via the according
  functions:
  \code{\link[cn.farms:summarizeFarmsGaussian]{cn.farms}},
  \code{\link[cn.farms:summarizeFarmsVariational]{cn.farms}}
  or \code{\link[cn.farms:summarizeFarmsExact]{cn.farms}}}

  \item{callParam}{Additional parameters for runtype (ff or
  bm) as well as cores for parallelization.}

  \item{summaryWindow}{Method for combination of the SNPs.
  Possible values are sl and fragment.}

  \item{returnValues}{List with return values. For possible
  values see summaryMethod.}

  \item{saveFile}{Name of the file to save.}
}
\value{
  Single-locus summarized data of an instance of
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
}
\description{
  The different probes of the SNPs of the array are
  summarized to a probeset.
}
\examples{
load(system.file("exampleData/normData.RData", package = "cn.farms"))
notes(experimentData(normData))$annotDir <-
        system.file("exampleData/annotation/pd.genomewidesnp.6/1.1.0",
                package = "cn.farms")
summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(10)
slData <- slSummarization(normData,
        summaryMethod = summaryMethod,
        summaryParam = summaryParam)
assayData(slData)$L_z[1:10, 1:10]

summaryMethod <- "Gaussian"
summaryParam <- list()
summaryParam$cyc <- c(10)
slData <- slSummarization(normData,
        summaryMethod = summaryMethod,
        summaryParam = summaryParam)
assayData(slData)$L_z[1:10, 1:10]

summaryMethod <- "Exact"
summaryParam <- list()
summaryParam$cyc <- c(10, 20)
slData <- slSummarization(normData,
        summaryMethod = summaryMethod,
        summaryParam = summaryParam)
assayData(slData)$L_z[1:10, 1:10]
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}
\seealso{
  \code{\link{summarizeFarmsExact}}
}

