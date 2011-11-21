\name{fragLengCorr}
\alias{fragLengCorr}
\title{Does a fragment length correction}
\usage{
  fragLengCorr(object, runtype = "ff",
    saveFile = "slDataFlc", ...)
}
\arguments{
  \item{object}{An instance of
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}}

  \item{runtype}{Mode how the results are saved. Possible
  values are ff or bm.}

  \item{...}{Further parameters passed to the correction
  method.}

  \item{saveFile}{Name of the file to save.}
}
\value{
  An instance of
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}.
}
\description{
  Does a fragment length correction
}
\examples{
load(system.file("exampleData/slData.RData", package = "cn.farms"))
slDataFlc <- fragLengCorr(slData)
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

