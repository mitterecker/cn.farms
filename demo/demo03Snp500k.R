## HapMap analysis on 500K
library(hapmap500knsp)
library(hapmap500ksty)  
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
filenames01 <- dir(path=celDir, full.names=TRUE)

celDir <- system.file("celFiles", package="hapmap500ksty")
filenames02 <- dir(path=celDir, full.names=TRUE)


###############################################################################
## process annotation
###############################################################################

if(exists("annotDir")) {
	createAnnotation(filenames=filenames01, annotDir=annotDir)
	createAnnotation(filenames=filenames02, annotDir=annotDir)	
} else {
	createAnnotation(filenames=filenames01)
	createAnnotation(filenames=filenames02)
}


###############################################################################
## process SNP data
###############################################################################


normMethod <- "SOR"

## normalization of SNP data
if(exists("annotDir")) {
	normData01 <- normalizeCels(filenames01, method=normMethod, cores, alleles=T, 
			annotDir=annotDir, runtype=runtype)
	normData02 <- normalizeCels(filenames02, method=normMethod, cores, alleles=T, 
			annotDir=annotDir, runtype=runtype)
} else {
	normData01 <- normalizeCels(filenames01, method=normMethod, cores, alleles=T, 
			runtype=runtype)
	normData02 <- normalizeCels(filenames02, method=normMethod, cores, alleles=T, 
			runtype=runtype)
}

if (testing) {
	## select one chromosome for further analysis
	load(file.path(notes(experimentData(normData))$annotDir, "featureSet.RData"))
	tmp <- featureSet$chrom[match(featureData(normData01)@data$fsetid, 
					featureSet$fsetid)]
	normData01 <- normData01[which(tmp == myChr), ]
	
	load(file.path(notes(experimentData(normData))$annotDir, "featureSet.RData"))
	tmp <- featureSet$chrom[match(featureData(normData02)@data$fsetid, 
					featureSet$fsetid)]
	normData02 <- normData02[which(tmp == myChr), ]
	
}


## sl FARMS
summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(10, 10)

callParam <- list(cores=cores, runtype=runtype)

slData01 <- slSummarization(normData01, 
		summaryMethod = summaryMethod, 
		summaryParam = summaryParam, 
		callParam = callParam, 
		summaryWindow = "std")
assayData(slData01)$L_z[1:10, ]

slData02 <- slSummarization(normData02, 
        summaryMethod = summaryMethod, 
        summaryParam = summaryParam, 
        callParam = callParam, 
        summaryWindow = "std")
assayData(slData02)$L_z[1:10, ]


###############################################################################
## combine SNP and CN data
###############################################################################

combData <- combineData(slData01, slData02, runtype=runtype)


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

mlData <- mlSummarization(combData, windowMethod, windowParam, 
		summaryMethod, summaryParam, callParam = callParam)

assayData(mlData)$intensity
