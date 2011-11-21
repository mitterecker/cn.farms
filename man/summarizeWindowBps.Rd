\name{summarizeWindowBps}
\alias{summarizeWindowBps}
\title{Combines neighbouring locations to windows}
\usage{
  summarizeWindowBps(phInf, fixedBps = 10000,
    upperLimit = 6)
}
\arguments{
  \item{phInf}{The locations on the chromosomes.}

  \item{fixedBps}{Size of the window in basepairs.}

  \item{upperLimit}{Maximal number of neigbouring locations
  to combine.}
}
\value{
  Indices for summarization
}
\description{
  Combines neighbouring locations to windows
}
\examples{
## create toy physical data
sizeTmp <- 30
phInf <- data.frame(
        chrom = rep("15", sizeTmp),
        start = seq(from = 1, by = 300, length.out = sizeTmp),
        end = seq(from = 3600, by = 300, length.out = sizeTmp),
        man_fsetid = paste("SNP_A-", seq(sizeTmp)+1000, sep = ""))
summarizeWindowBps(phInf)
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

