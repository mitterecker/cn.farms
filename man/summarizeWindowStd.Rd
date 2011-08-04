\name{summarizeWindowStd}
\alias{summarizeWindowStd}
\title{Combines neighbouring locations to windows}
\usage{
  summarizeWindowStd(phInf, windowSize = 3, overlap = TRUE)
}
\arguments{
  \item{phInf}{The locations on the chromosomes.}

  \item{windowSize}{Size of how many Locations should be
  combined.}

  \item{overlap}{States if the windows should overlap.}
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
chrom=rep("15", sizeTmp),
start=seq(from=1, by=300, length.out=sizeTmp),
end=seq(from=3600, by=300, length.out=sizeTmp),
man_fsetid=paste("SNP_A-", seq(sizeTmp)+1000, sep=""))
summarizeWindowStd(phInf)
}
\author{
  Djork-Arne Clevert \email{okko@clevert.de} and Andreas
  Mitterecker \email{mitterecker@bioinf.jku.at}
}

