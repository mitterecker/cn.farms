\name{plotSmoothScatter}
\alias{plotSmoothScatter}
\title{Creates a smooth scatter plot...}
\usage{plotSmoothScatter(object, variable, chrom, start, end, ylim, pdfname, ...)
}
\description{Creates a smooth scatter plot}
\value{A graph.}
\author{Andreas Mitterecker}
\arguments{\item{object}{An instance of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}.}
\item{variable}{States which variable of the assayData should be plotted.}
\item{chrom}{The chromosome you want to plot.}
\item{start}{The physical start position.}
\item{end}{The physical end position.}
\item{ylim}{The limits for the y axis.}
\item{pdfname}{The name of the pdf file.}
\item{...}{Further arguments passed to smoothScatter function.}
}
\examples{load(system.file("exampleData/slData.RData", package="cn.farms"))
plotSmoothScatter(slData[, 1:3], chrom="23")}

