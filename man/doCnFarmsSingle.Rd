\name{doCnFarmsSingle}
\alias{doCnFarmsSingle}
\title{Does the whole cn.farms process in one call}
\usage{
  doCnFarmsSingle(celfiles, samplenames, normalization)
}
\arguments{
  \item{celfiles}{The celfiles which you want to process
  with the whole path. Either a vector or a matrix with two
  columns for combined analysis e.g. 500K Array.}

  \item{samplenames}{An optional vector with the same
  dimension as the number of cel files}

  \item{normalization}{The normalization method you want to
  use.}
}
\value{
  The ready cn.FARMS results.
}
\description{
  Works for all kind of Affymetrix SNP arrays
}
\author{
  Andreas Mitterecker
}

