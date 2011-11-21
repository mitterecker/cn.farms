\name{createMatrix}
\alias{createMatrix}
\title{Creates the needed matrix}
\usage{
  createMatrix(runtype, nrow, ncol, type = "double",
    bmName = "NA")
}
\arguments{
  \item{runtype}{Mode how the results are saved. Possible
  values are ff or bm. If ff is chosen the data will not be
  saved automatically. With bm the results will be saved
  permanently.}

  \item{nrow}{nrow}

  \item{ncol}{ncol}

  \item{type}{type}

  \item{bmName}{Identifier for ff name}
}
\value{
  a matrix
}
\description{
  Creates the needed matrix
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

