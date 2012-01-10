\name{normalizeNone}
\alias{normalizeNone}
\title{Runs the SOR normalization on microarray data}
\usage{
  normalizeNone(filenames, cores = 1, annotDir = NULL,
    alleles = FALSE, runtype = "ff", cyc = 5,
    pkgname = NULL, saveFile = "Sor")
}
\arguments{
  \item{filenames}{an absolute path of the CEL files}

  \item{cores}{cores}

  \item{annotDir}{annotDir}

  \item{alleles}{alleles}

  \item{cyc}{states the number of cycles for the EM
  algorithm.}

  \item{runtype}{Mode how the results are saved. Possible
  values are ff or bm. If ff is chosen the data will not be
  saved automatically. With bm the results will be saved
  permanently.}

  \item{pkgname}{Optional parameter for the CEL mapping.}

  \item{saveFile}{Name of the file to save.}
}
\value{
  An instance of
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
}
\description{
  Runs the SOR normalization on microarray data
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

