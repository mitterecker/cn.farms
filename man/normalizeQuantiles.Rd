\name{normalizeQuantiles}
\alias{normalizeQuantiles}
\title{Normalization Quantiles}
\usage{
  normalizeQuantiles(filenames, cores = 1, batch = NULL,
  annotDir = NULL, runtype = "ff", pkgname = NULL, saveFile
  = "normDataQuant")
}
\arguments{
  \item{filenames}{filenames}

  \item{cores}{cores}

  \item{batch}{batch}

  \item{annotDir}{annotDir}

  \item{runtype}{Mode how the results are saved. Possible
  values are ff or bm. If ff is chosen the data will not be
  saved automatically. With bm the results will be saved
  permanently.}

  \item{pkgname}{Optional parameter for the CEL mapping.}

  \item{saveFile}{Name of the file to save.}
}
\value{
  The normalized data.
}
\description{
  Normalization Quantiles
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

