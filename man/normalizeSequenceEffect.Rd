\name{normalizeSequenceEffect}
\alias{normalizeSequenceEffect}
\title{Correction for probe sequence effects}
\usage{
  normalizeSequenceEffect(object, annotDir = NULL,
    runtype = "ff", saveFile = "seqNorm")
}
\arguments{
  \item{object}{an instance of
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}}

  \item{annotDir}{the directory where the anntotation can
  be found}

  \item{runtype}{mode how the results are saved. Possible
  values are ff or bm. If ff is chosen the data will not be
  saved automatically.}

  \item{saveFile}{name of the file to save.}
}
\value{
  Some data
}
\description{
  Correction for probe sequence effects
}
\author{
  Andreas Mitterecker
}

