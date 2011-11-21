\name{dnaCopySf}
\alias{dnaCopySf}
\title{Runs DNAcopy in parallel mode}
\usage{
  dnaCopySf(x, chrom, maploc, cores = 1, smoothing, ...)
}
\arguments{
  \item{x}{A matrix with data of the copy number
  experiments}

  \item{chrom}{The chromosomes (or other group identifier)
  from which the markers came}

  \item{maploc}{The locations of marker on the genome}

  \item{cores}{Number of cores to use}

  \item{smoothing}{States if smoothing of the data should
  be done}

  \item{...}{Further parameter for the function segment of
  DNAcopy}
}
\value{
  An instance of
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
  containing the segments.
}
\description{
  This function even works very well with ff matrices,
}
\examples{
load(system.file("exampleData/mlData.RData", package = "cn.farms"))
mlData <- mlData[, 1:3]
colnames(assayData(mlData)$L_z) <- sampleNames(mlData)
segments <- dnaCopySf(
        x         = assayData(mlData)$L_z,
        chrom     = fData(mlData)$chrom,
        maploc    = fData(mlData)$start,
        cores     = 1,
        smoothing = FALSE)
fData(segments)$data
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

