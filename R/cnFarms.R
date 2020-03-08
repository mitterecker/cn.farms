#' cn.farms
#' 
#' Wrapper for the cn.farms algorithm
#' 
#' @param filenames the absolute filepaths of the CEL files. 
#' @param cores number of parallel instances.
#' @param runtype either ff (results will be lost after closing the R session) or bm (results will be saved after closing the R session). 
#' @return An instance of \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
#' containing the results of the analysis.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @examples
#' \dontrun{
#' require('hapmapsnp6') 
#' celDir <- system.file('celFiles', package = 'hapmapsnp6')
#' filenames <- dir(path = celDir, full.names = TRUE)
#' cn.farms(filenames = filenames)
#' }
cn.farms <- function(filenames, cores = 1, runtype = "bm") {
  
  message("Be aware that results will be stored in the current working directory!")
  message("Also assure that all CEL files are from the same type!")
  
  ## settings for ff
  dir.create("ffObjects/ff", showWarnings = FALSE, recursive = TRUE)
  ldPath(file.path(getwd(), "ffObjects"))
  options(fftempdir = file.path(ldPath(), "ff"))
  
  
  
  ## check if filenames are correct
  if (!all(file.exists(filenames))) {
    stop("The location of the CEL files doesn't exist! Please provide the full path!")
  }
  
  ## create annotation
  message("Creating annotation")
  createAnnotation(filenames = filenames, checks = FALSE)
  
  
  ## normalize SNP data
  normMethod <- "quantiles"
  normData <- normalizeCels(filenames, method = normMethod, cores = cores, alleles = TRUE, 
      runtype = runtype)
  
  ## run single locus FARMS on SNP data
  summaryMethod <- "Variational"
  summaryParam <- list()
  summaryParam$cyc <- c(10)
  callParam <- list(cores = cores, runtype = runtype)
  slData <- slSummarization(normData, summaryMethod = summaryMethod, summaryParam = summaryParam, 
      callParam = callParam, summaryWindow = "std")
  
  ## process non polymorphic data
  npData <- normalizeNpData(filenames, cores, runtype = runtype)
  
  ## combine SNP and CN data
  combData <- combineData(slData, npData, runtype = runtype)
  
  ## run multi loci FARMS
  windowMethod <- "std"
  windowParam <- list()
  windowParam$windowSize <- 9
  windowParam$overlap <- TRUE
  summaryMethod <- "Variational"
  summaryParam <- list()
  summaryParam$cyc <- c(20)
  callParam <- list()
  callParam <- list(cores = cores, runtype = runtype)
  
  mlData <- mlSummarization(combData, windowMethod, windowParam, summaryMethod, 
      summaryParam, callParam = callParam)
  
  separateAnalysis <- FALSE
  if (separateAnalysis) {
    mlDataSnp <- mlSummarization(slData, windowMethod, windowParam, summaryMethod, 
        summaryParam, callParam = callParam, saveFile = "mlSnp")
    
    mlDataCn <- mlSummarization(npData, windowMethod, windowParam, summaryMethod, 
        summaryParam, callParam = callParam, saveFile = "mlCn")
  }
  
  return(mlData)
  
}
