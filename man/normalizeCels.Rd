\name{normalizeCels}
\alias{normalizeCels}
\title{Wrapper for the normalization functions}
\usage{
  normalizeCels(filenames,
    method = c("SOR", "quantiles", "none"), cores = 1,
    alleles = FALSE, runtype = "bm", annotDir = NULL,
    saveFile = "normData", ...)
}
\arguments{
  \item{filenames}{The absolute path of the CEL files as a
  list.}

  \item{method}{The normalization method. Possible methods
  so far: SOR, quantiles}

  \item{cores}{Number of cores for used for
  parallelization.}

  \item{alleles}{States if information for allele A and B
  should be given back.}

  \item{runtype}{Mode how the results are saved. Possible
  values are ff or bm. If ff is chosen the data will not be
  saved automatically. With bm the results will be saved
  permanently.}

  \item{annotDir}{An optional annotation directory.}

  \item{saveFile}{Name of the file to save.}

  \item{...}{Further parameters for the normalization
  method.}
}
\value{
  An ExpressionSet object with the normalized data.
}
\description{
  This functions provides different normalization methods
  for microarray data. At the moment only SOR and quantile
  normalization are implemented.
}
\examples{
\dontrun{
library("hapmapsnp6")
celDir <- system.file("celFiles", package = "hapmapsnp6")
filenames <- dir(path = celDir, full.names = TRUE)
createAnnotation(filenames = filenames)
normData <- normalizeCels(filenames, method = "SOR")
}
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

