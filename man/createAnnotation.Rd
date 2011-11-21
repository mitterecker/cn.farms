\name{createAnnotation}
\alias{createAnnotation}
\title{Creation of annotation files}
\usage{
  createAnnotation(filenames = NULL, annotation = NULL,
    annotDir = NULL, checks = TRUE)
}
\arguments{
  \item{filenames}{An absolute path of the CEL files to
  process.}

  \item{annotation}{Optional parameter stating the
  annotation from a pd-mapping.}

  \item{annotDir}{Optional parameter stating where the
  annotation should go.}

  \item{checks}{States if sanity checks should be done.}
}
\value{
  \code{NULL}
}
\description{
  Annotation files for cn.farms are created
}
\note{
  The annotation files used for cn.farms will be placed in
  the current work directory under annotations.
}
\examples{
\dontrun{
library("hapmapsnp6")
celDir <- system.file("celFiles", package = "hapmapsnp6")
filenames <- dir(path = celDir, full.names = TRUE)
createAnnotation(filenames = filenames)
}
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

