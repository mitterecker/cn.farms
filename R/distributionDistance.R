## Some distance measures as mentioned in Singh (2006), "Study Of Some
## Distance Measures For Language And Encoding Identification",
## http://clair.si.umich.edu/clair/anthology/query.cgi?type=Paper&id=W06-1109.


#' Computes the distribution distance 
#' 
#' Be aware that this function is implemented quite slow.
#' 
#' @param intensityData A matrix or an AffyBatch object.
#' @param method The method you want to use.
#' @param useSubset Logical. States if only a subset should be used.
#' @param subsetFraction The fraction of the subset.
#' @param useQuantileReference Logical for a quantile reference.
#' @return Computes the distribution distance
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export
#' @examples
#' load(system.file("exampleData/normData.RData", package="cn.farms"))
#' x <- assayData(normData)$intensity[, 1:3]
#' y <- distributionDistance(x)
#' attr(y, "Labels") <- substr(sampleNames(normData), 1, 7)
#' plotDendrogram(y)
distributionDistance <- function(
        intensityData, 
        method = c("JSDiv", "KLDiv", "KLInf"),
	    useSubset = T, 
	    subsetFraction = 0.25, 
	    useQuantileReference = FALSE) {
    
	## Parse arguments, stored in variables
	## If no value is specified, the first option in the list is taken as 
    ## the default value  
	method <- match.arg(method)
	if(!is.matrix(intensityData)){
        stop("A matrix as input is needed!")  
	}
	
	nbrOfRows <- nrow(intensityData)
	nbrOfCols <- ncol(intensityData)
	
	
	## Computing the CDFs is computational expensive. Costs can be reduced by 
    ## data subsetting.
	## The distance measure will be computed on a subset if useSubset is set 
    ## to TRUE (default)
	## The subset size can be controlled by subsetFraction (25% of the data 
    ## is used by default)
	if (useSubset){ 
		### use fix seed for reproducibility 
		set.seed(123)    
		### get data subset
		intensityData <- intensityData[sample.int(nbrOfRows, 
                        floor(nbrOfRows * subsetFraction)),]
	}
	
	
	## Scale data to remove offset and scale issues
	intensityData <- scale(intensityData, center=T, scale=T)
	
  	## use a fix reference distribution 
	if(useQuantileReference){
		### compute quantile distribution as reference 
      	referenceDistribution <- Biobase::rowMedians(apply(intensityData, 2, sort))
       	res <- rep(NA, nbrOfCols)
       	for(i in 1:nbrOfCols){
   			### compute distance against the reference
			res[i] <- calcDistance(
                    x=intensityData[, i], 
                    y=referenceDistribution, 
                    method=method)
		}
		
		attributes(res) <- list(Size = nbrOfCols, Labels = colnames(intensityData), 
		        methods=method, class = "vector")
	} else {
       	res <- rep(NA, nbrOfCols * (nbrOfCols - 1) / 2)
   		count <- 1
   		## compute pairwise distances
		for(i in 1:(nbrOfCols - 1)){
   			for(j in (i + 1):nbrOfCols){
      			res[count] <- calcDistance(
                        x=intensityData[, i], 
                        y=intensityData[, j], 
                        method=method)
        		count <- count + 1
           	}
        }
	    
		attributes(res) <- list(
                Size=nbrOfCols, 
                Labels=colnames(intensityData), 
		        methods=method, 
                class="dist")
        
    }
    
    return(res)
    
}


calcDistance <- function(x, y, method=c("JSDiv", "KLDiv", "KLInf")){ 
    
	## eps is added to avoid numerical issues
    eps <- .Machine$double.eps 
    
    ## compute kernel density estimates for both distributions
	P <- density(x, n=10000, from=min(c(x, y)), to=max(c(x, y)))$y + eps  
	Q <- density(y, n=10000, from=min(c(x, y)), to=max(c(x, y)))$y + eps
    
    ### calc distance measure
	distance <- switch(method, 
			JSDiv = JSDivergence(P, Q),
            KLDiv = KLDivergence(P, Q),
            KLInf = KLInformation(P, Q)
    )
    
	return(distance) 	
	
}

## Kullback-Leibler divergences
## http://en.wikipedia.org/wiki/Kullback–Leibler_divergence:
KLDivergence <- function(p, q){
    sum(p * log(p / q))
}    

KLInformation <- function(p, q){
    sum((p - q) * log(p / q))
}


## Jensen-Shannon divergence
## http://en.wikipedia.org/wiki/Jensen–Shannon_divergence
JSDivergence <- function(p, q){
	m <- 0.5 * (p + q)
	0.5 * (KLDivergence(p, m) +  KLDivergence(q, m))
}
