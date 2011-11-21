## HapMap analysis on 250K
library(hapmap500knsp)  
library(cn.farms)

###############################################################################
## general settings
###############################################################################

workDir <- tempdir()
dir.create(workDir, showWarnings=F, recursive=T)
setwd(workDir)

cores <- 3
runtype <- "ff"
#runtype <- "bm"

## settings for test run
testing <- T
myChr <- "16"

## settings for ff
dir.create("ffObjects/ff", showWarnings=F, recursive=T)
oligoClasses::ldPath(file.path(getwd(), "ffObjects"))
options(fftempdir = file.path(oligoClasses::ldPath(), "ff"))

## CEL files
celDir <- system.file("celFiles", package="hapmap500knsp")
filenames <- dir(path=celDir, full.names=TRUE)

###############################################################################
## process annotation
###############################################################################

if(exists("annotDir")) {
	createAnnotation(filenames=filenames, annotDir=annotDir)	
} else {
	createAnnotation(filenames=filenames)
}


###############################################################################
## process SNP data
###############################################################################


normMethod <- "SOR"

## normalization of SNP data
if(exists("annotDir")) {
	normData <- normalizeCels(filenames, method=normMethod, cores, alleles=T, 
			annotDir=annotDir, runtype=runtype)
} else {
	normData <- normalizeCels(filenames, method=normMethod, cores, alleles=T, 
			runtype=runtype)
}

assayData(normData)$intensity[1:10, ]
head(featureData(normData)@data)

## include more phenoData e.g. gender
phenoData(normData)$gender <- rep(1, length(phenoData(normData)$filenames))
phenoData(normData)
sampleNames(normData)
experimentData(normData)@title <- "HapMap data"


if (testing) {
	## select one chromosome for further analysis
	if (!exists("normDataBak")) normDataBak <- normData
	load(file.path(notes(experimentData(normData))$annotDir, "featureSet.RData"))
	tmp <- featureSet$chrom[match(featureData(normData)@data$fsetid, 
					featureSet$fsetid)]
	normData <- normData[which(tmp == myChr), ]
}


## sl FARMS
summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(10, 10)

callParam <- list(cores=cores, runtype=runtype)

slData <- slSummarization(normData, 
		summaryMethod = summaryMethod, 
		summaryParam = summaryParam, 
		callParam = callParam, 
		summaryWindow = "std")
assayData(slData)$L_z[1:10, ]


###############################################################################
## run multi loci FARMS
###############################################################################

windowMethod <- "std"
windowParam <- list()
windowParam$windowSize <- 5
windowParam$overlap <- F

summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(20, 20)

callParam <- list()
callParam = list(cores=cores, runtype=runtype)

mlData <- mlSummarization(slData, windowMethod, windowParam, 
		summaryMethod, summaryParam, callParam = callParam)

assayData(mlData)$intensity
