## The examples from the man pages
library(cn.farms)
dontrun <- T



## combineData
load(system.file("exampleData/normData.RData", package="cn.farms"))
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
combData


## createAnnotation
if (!dontrun) {
    library("hapmapsnp6") 
    celDir <- system.file("celFiles", package="hapmapsnp6")
    filenames <- dir(path=celDir, full.names=TRUE)
    createAnnotation(filenames=filenames)
}


## distributionDistance (see plotDendrogram)
load(system.file("exampleData/normData.RData", package="cn.farms"))
x <- assayData(normData)$intensity[, 1:3]
y <- distributionDistance(x)
attr(y, "Labels") <- substr(sampleNames(normData), 1, 7)
plotDendrogram(y)


## dnaCopySf
load(system.file("exampleData/mlData.RData", package="cn.farms"))
mlData <- mlData[, 1:3]
colnames(assayData(mlData)$L_z) <- sampleNames(mlData)
segments <- dnaCopySf(
        x         = assayData(mlData)$L_z, 
        chrom     = featureData(mlData)@data$chrom, 
        maploc    = featureData(mlData)@data$start, 
        cores     = 1, 
        smoothing = FALSE)
featureData(segments)@data


## fragLengCorr
load(system.file("exampleData/slData.RData", package="cn.farms"))
slDataFlc <- fragLengCorr(slData)


## mlSummarization
load(system.file("exampleData/slData.RData", package="cn.farms"))
windowMethod <- "std"
windowParam <- list()
windowParam$windowSize <- 5
windowParam$overlap <- TRUE
summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(20)
mlData <- mlSummarization(slData, windowMethod, windowParam, 
        summaryMethod, summaryParam)
assayData(mlData)


## normalizeCels
if (!dontrun) {
    library("hapmapsnp6") 
    celDir <- system.file("celFiles", package="hapmapsnp6")
    filenames <- dir(path=celDir, full.names=TRUE)
    createAnnotation(filenames=filenames)
    normData <- normalizeCels(filenames, method="SOR")
}


## normalizeNpData
if (!dontrun) {
    library("hapmapsnp6") 
    celDir <- system.file("celFiles", package="hapmapsnp6")
    filenames <- dir(path=celDir, full.names=TRUE)
    createAnnotation(filenames=filenames)
    npData <- normalizeNpData(filenames)
}


## plotDendrogram
load(system.file("exampleData/normData.RData", package="cn.farms"))
x <- assayData(normData)$intensity[, 1:3]
y <- distributionDistance(x)
attr(y, "Labels") <- substr(sampleNames(normData), 1, 7)
plotDendrogram(y)


## plotDensity
load(system.file("exampleData/slData.RData", package="cn.farms"))
plotDensity(assayData(slData)$intensity)


## plotEvalIc
load(system.file("exampleData/slData.RData", package="cn.farms"))
load(system.file("exampleData/testSegments.RData", package="cn.farms"))
plotEvalIc(slData, featureData(testSegments)@data, 
        variable=assayData(slData)$L_z[, 1],  chrom=23)


## plotRegions
load(system.file("exampleData/slData.RData", package="cn.farms"))
load(system.file("exampleData/testSegments.RData", package="cn.farms"))
plotRegions(slData, testSegments, addInd=10, ylim=c(-2, 2), 
        variable="L_z", colorVersion=1, plotLegend=TRUE, pdfname="slData.pdf")

## plotSmoothScatter
load(system.file("exampleData/slData.RData", package="cn.farms"))
plotSmoothScatter(slData[, 1:3], chrom="23")


## plotVioline
load(system.file("exampleData/normData.RData", package="cn.farms"))
normData <- normData[, 1:10]
groups <- seq(sampleNames(normData))
plotViolines(normData, variable="intensity", groups, xlab="Intensity values")


## slSummarization
load(system.file("exampleData/normData.RData", package="cn.farms"))
experimentData(normData)@other$annotDir <- 
        system.file("exampleData/annotation/pd.genomewidesnp.6/1.1.0",
                package="cn.farms")

summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(10)
slData <- slSummarization(normData, 
        summaryMethod = summaryMethod, 
        summaryParam = summaryParam)
assayData(slData)$L_z[1:10, 1:10]

load(system.file("exampleData/normData.RData", package="cn.farms"))
summaryMethod <- "Gaussian"
summaryParam <- list()
summaryParam$cyc <- c(10)
slData <- slSummarization(normData, 
        summaryMethod = summaryMethod, 
        summaryParam = summaryParam)
assayData(slData)$L_z[1:10, 1:10]

load(system.file("exampleData/normData.RData", package="cn.farms"))
summaryMethod <- "Exact"
summaryParam <- list()
summaryParam$cyc <- c(10, 20)
slData <- slSummarization(normData, 
        summaryMethod = summaryMethod, 
        summaryParam = summaryParam)
assayData(slData)$L_z[1:10, 1:10]


## sparseFarmsC
x <- matrix(rnorm(100, 11), 20, 5)
sparseFarmsC(x, 50)


## summarizeFarmsExact
x <- matrix(rnorm(100, 11), 20, 5)
summarizeFarmsExact(x)



## summarizeFarmsGaussian
x <- matrix(rnorm(100, 11), 20, 5)
summarizeFarmsGaussian(x)

## summarizeFarmsMethods
summarizeFarmsMethods()

## summarizeFarmsVariational
x <- matrix(rnorm(100, 11), 20, 5)
summarizeFarmsVariational(x)


## summarizeWindowBps
## create toy physical data
sizeTmp <- 30
phInf <- data.frame(
        chrom=rep("15", sizeTmp),
        start=seq(from=1, by=300, length.out=sizeTmp), 
        end=seq(from=3600, by=300, length.out=sizeTmp),
        man_fsetid=paste("SNP_A-", seq(sizeTmp)+1000, sep=""))
summarizeWindowBps(phInf)


## summarizeWindowMethods
summarizeWindowMethods()

## summarizeWindowStd
sizeTmp <- 30
phInf <- data.frame(
        chrom=rep("15", sizeTmp),
        start=seq(from=1, by=300, length.out=sizeTmp), 
        end=seq(from=3600, by=300, length.out=sizeTmp),
        man_fsetid=paste("SNP_A-", seq(sizeTmp)+1000, sep=""))
summarizeWindowStd(phInf)



