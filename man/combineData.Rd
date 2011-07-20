\name{combineData}
\alias{combineData}
\title{Combine two ExpressionSet objects...}
\usage{combineData(object01, object02, obj01Var="intensity", obj02Var="intensity",
    runtype="ff", saveFile="combData")
}
\description{Combine two ExpressionSet objects}
\details{Suitable for SNP or non-polymorphic data which were already processed with 
single locus FARMS}
\value{An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}.}
\author{Djork-Arne Clevert \email{okko@clevert.de} and 
Andreas Mitterecker \email{mitterecker@bioinf.jku.at}}
\arguments{\item{object01}{An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
either with SNP or non-polymorphic data}
\item{object02}{An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
either with SNP or non-polymorphic data}
\item{obj01Var}{States the variable which should be combined from the 
assayData slot. Default is intensity.}
\item{obj02Var}{States the variable which should be combined from the 
assayData slot. Default is intensity.}
\item{runtype}{Mode how the results are saved. Possible values are ff or bm. 
If ff is chosen the data will not be saved automatically. 
With bm the results will be saved permanently.}
\item{saveFile}{Name of the file to save.}
}
\examples{load(system.file("exampleData/normData.RData", package="cn.farms"))
experimentData(normData)@other$annotDir <- 
system.file("exampleData/annotation/pd.genomewidesnp.6/1.1.0",
package="cn.farms")
summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(10)
slData <- slSummarization(normData, 
summaryMethod = summaryMethod, 
summaryParam = summaryParam)
assayData(slData)$L_z[1:10, ] 
combData <- combineData(slData, slData)
combData}

