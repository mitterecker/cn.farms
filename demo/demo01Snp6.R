## HapMap analysis on SNP6
library(hapmapsnp6)  
library(oligoClasses)
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
testing <- TRUE
myChr <- "16"

## settings for ff
dir.create("ffObjects/ff", showWarnings=F, recursive=T)
ldPath(file.path(getwd(), "ffObjects"))
options(fftempdir = file.path(ldPath(), "ff"))

## CEL files
celDir <- system.file("celFiles", package="hapmapsnp6")
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
	load(file.path(experimentData(normData)@other$annotDir, "featureSet.RData"))
	tmp <- featureSet$chrom[match(featureData(normData)@data$fsetid, 
					featureSet$fsetid)]
	normData <- normData[which(tmp == myChr), ]
}


## sl FARMS
summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(10)

callParam <- list(cores=cores, runtype=runtype)

slData <- slSummarization(normData, 
		summaryMethod = summaryMethod, 
		summaryParam = summaryParam, 
		callParam = callParam, 
		summaryWindow = "std")
assayData(slData)$L_z[1:10, ]

## fragment length correction
#slData <- fragLengCorr(slData, cores=cores)
#experimentData(slData)@other


###############################################################################
## process non polymorphic data
###############################################################################

if (exists("annotDir")) {
	npData <- normalizeNpData(filenames, cores, annotDir=annotDir)	
} else {
	npData <- normalizeNpData(filenames, cores, runtype=runtype)
}

if (testing) {
	## select one chromosome for further analysis
	npDataBak <- npData
	load(file.path(experimentData(npData)@other$annotDir, 
					"featureSetCNV.RData"))
	tmp <- featureSetCNV$chrom[match(
					featureData(npData)@data$fsetid, featureSetCNV$fsetid)]
	npData <- npData[which(tmp == myChr), ]
}

### fragment length correction
#npData <- fragLengCorr(npData, cores=cores)
#experimentData(npData)@other
#assayData(npData)$intensity

###############################################################################
## combine SNP and CN data
###############################################################################

combData <- combineData(slData, npData, runtype=runtype)


###############################################################################
## run multi loci FARMS
###############################################################################

windowMethod <- "std"
windowParam <- list()
windowParam$windowSize <- 5
windowParam$overlap <- F

summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(20)

callParam <- list()
callParam <- list(cores=cores, runtype=runtype)

mlData <- mlSummarization(combData, windowMethod, windowParam, 
		summaryMethod, summaryParam, callParam = callParam)

assayData(mlData)$intensity

