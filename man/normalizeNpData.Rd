\name{normalizeNpData}
\alias{normalizeNpData}
\title{Processes the non-polymorphic data...}
\usage{normalizeNpData(filenames, cores=1, annotDir, runtype="ff", saveFile="npData",
    method=c("baseline", "quantiles"))
}
\description{Processes the non-polymorphic data}
\details{Normalization for non-polymorphic data for Affymetrix SNP5 and SNP6}
\value{An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}} 
containing the non-polymorphic data of the microarray.}
\author{Djork-Arne Clevert \email{okko@clevert.de} and 
Andreas Mitterecker \email{mitterecker@bioinf.jku.at}}
\arguments{\item{filenames}{the absolute path of the CEL files as a list}
\item{cores}{number of cores for used for parallelization}
\item{annotDir}{Optional annotation directory.}
\item{runtype}{Mode how the results are saved. Possible values are ff or bm. 
If ff is chosen the data will not be saved automatically. 
With bm the results will be saved permanently.}
\item{saveFile}{Name of the file to save.}
\item{method}{The method for the normalization.}
}
\examples{\dontrun{
library("hapmapsnp6") 
celDir <- system.file("celFiles", package="hapmapsnp6")
filenames <- dir(path=celDir, full.names=TRUE)
createAnnotation(filenames=filenames)
npData <- normalizeNpData(filenames)
}}

