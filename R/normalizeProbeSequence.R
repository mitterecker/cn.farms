#' Correction for probe sequence effects
#' @param object an instance of 
#' \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
#' @param annotDir the directory where the anntotation can be found 
#' @param runtype mode how the results are saved. Possible values are ff or bm. 
#' If ff is chosen the data will not be saved automatically. 
#' @param saveFile name of the file to save.
#' @return Some data
#' @author Andreas Mitterecker
#' @export
#' @
normalizeSequenceEffect <- function(
        object, 
        annotDir = NULL, 
        runtype = "ff", 
        saveFile = "seqNorm") {

    pkgname <- object@annotation
    nbrOfProbes <- nrow(object)
    nbrOfSamples <- ncol(object)
    
   
    if (pkgname != "pd.genomewidesnp.6") {
        stop("Sequence effect normalization is not implemented for this array type!")
    }
    
    ## assure correct file extension
    saveFile <- gsub("\\.RData", "", saveFile)
    saveFile <- gsub("\\.rda", "", saveFile)
    saveFile <- paste(saveFile, ".RData", sep = "")
       
    if (runtype == "bm" & file.exists(saveFile)) {
        message("Sequence normalization has already been done")
        message("Trying to load  data ...")
        load(saveFile)
        return(slData)
    }
    
    ## load and assign sequences
    if (is.null(annotDir)) {
        vers <- dir(file.path("annotation", pkgname))[1]
        annotDir <- normalizePath(file.path("annotation", pkgname, vers))
    }
    
    load(file.path(annotDir, "sequence.RData"))
    load(file.path(annotDir, "pmfeature.RData"))
    load(file.path(annotDir, "featureSet.RData"))
    
    featureSet <- featureSet
    pmfeature <- pmfeature
    sequence <- sequence
    
    if (pkgname == "pd.genomewidesnp.6") {
        load(file.path(annotDir, "sequenceCNV.RData"))
        load(file.path(annotDir, "pmfeatureCNV.RData"))
        load(file.path(annotDir, "featureSetCNV.RData"))
        
        sequenceCNV <- sequenceCNV
        pmfeatureCNV <- pmfeatureCNV
        featureSetCNV <- featureSetCNV
        
        seq <- rbind(sequenceCNV[, c("fid", "seq")], sequence[, c("fid", "seq")])
        feat <- rbind(featureSetCNV[, c("fsetid", "man_fsetid")], 
                featureSet[, c("fsetid", "man_fsetid")])
        pmfeat <- rbind(pmfeatureCNV[, c("fid", "fsetid")], 
                        pmfeature[, c("fid", "fsetid")])
        
    } else {
        seq <- sequence[, c("fid", "seq")]
        feat <-  featureSet[, c("fsetid", "man_fsetid")]
        pmfeat <- pmfeature[, c("fid", "fsetid")]
    }

    tmp01 <- feat[match(featureData(object)$man_fsetid, feat$man_fsetid), "fsetid"]
    tmp02 <- pmfeat[match(tmp01, pmfeat$fsetid), "fid"]
    seqs <- seq[match(tmp02, seq$fid), "seq"]
       
    rm(feat, pmfeat, tmp01, tmp02)
    
    ## sequence normalization
    trainIdx <- which(featureData(object)$chrom %in% 1:22)
    trainSeqs <- seqs[trainIdx]
    designMat <- getProbePositionEffectDesignMatrix(trainSeqs, verbose = FALSE)
    
    out <- createMatrix(runtype, nbrOfProbes, nbrOfSamples, 
            type = "double", bmName = gsub("\\.RData", "", saveFile))
    
    for (i in 1:ncol(object)){
        
        fit <- fit_model(
                designMat, 
                assayData(object)$intensity[trainIdx, i])
        
        out[, i] <- suppressWarnings(predict_model(
                fit, 
                seqs, 
                assayData(object)$intensity[, i]))
        
    }

    slData <- new("ExpressionSet")
    assayData(slData) <- list(intensity = out)
    phenoData(slData) <- phenoData(object)
    featureData(slData) <- featureData(object)
    experimentData(slData) <- experimentData(object)
    annotation(slData) <- annotation(object)
    sampleNames(slData) <- sampleNames(object)  
    experimentData(slData)@other$seqNorm <- 1
    cat(paste(Sys.time(), " |   Sequence correction done \n", 
                    sep = ""))
    
    if (runtype == "bm") {
        cat(paste(Sys.time(), "|   Saving normalized data \n"))
        save(slData, file = saveFile)
    }
    
    return(slData)
    
}

## taken from CRMA v2 (www.aroma-project.org)

getEffects <- function(fit, intercept=FALSE, ...) {
	params <- fit$params;
	B <- fit$B;
	map <- fit$map;
	
	factors <- names(params);
	factors <- setdiff(factors, "intercept");
	F <- length(factors);
	rho <- matrix(0, nrow=nrow(B), ncol=F+1);
	for (kk in 1:F) {
		key <- factors[kk];
		rho[,kk] <- B %*% params[[key]];
	}
	
	if (intercept) {
		rho <- rho[,1:F];
		colnames(rho) <- factors;
	} else {
		rho <- rho - rowSums(rho[,1:F])/(F+1);
		colnames(rho) <- names(map)[-1];
	}
	
	rho;
} # getEffects()



predict_effect <- function(object, seqs, ..., verbose=FALSE) {
    
	fit <- object;
	
	if (is.character(seqs)) {
		
		K <- length(seqs);
		P <- nchar(seqs[1]);
		
		seqs <- paste(seqs, collapse="");
		seqs <- strsplit(seqs, split="", fixed=TRUE)[[1]];
		map <- c("NA"=0, A=1, C=2, G=3, T=4);
		names <- names(map);
		map <- as.raw(map);
		names(map) <- names;
		values <- map[-1];
		seqs <- match(seqs, names(values));
		seqs <- as.raw(seqs);
		seqs <- matrix(seqs, nrow=K, ncol=P, byrow=TRUE);
		attr(seqs, "map") <- map;
		
	}
	
	
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Get probe-position effects
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	rho <- getEffects(fit);
	
	
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Calculate predicted values
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	K <- nrow(seqs);
	P <- ncol(seqs);
	
	# Allocate probe-affinity vector
	phi <- double(K);
	
	map <- attr(seqs, "map");
	values <- map[-1];
	factors <- names(values);
	#  values <- map[c("A", "C", "G", "T")];
	
	
	# Is it safe to use the "quick" approach for prediction?
	# The quick approach is 6-7 times faster. /HB 2008-12-03
	safeValues <- as.raw(1:4);
	names(safeValues) <- c("A", "C", "G", "T");
	safe <- identical(values, safeValues);
	
	
	if (safe) {
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# This approach assumes that the 'values' are A=01, C=02, G=03, T=04
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# Identify for which cells the sequences are known
		known <- whichVector(seqs[,1] != as.raw(0));
		K2 <- length(known);
		phi2 <- double(K2);
		
		# For each position
		for (pp in seq(length=P)) {
			
			# Get the nucleotides at this position for all sequences
			seqsPP <- seqs[known,pp];
			
			seqsPP <- as.integer(seqsPP);
			rhoPP <- rho[pp,];
			names(rhoPP) <- NULL;
			phi2 <- phi2 + rhoPP[seqsPP];
			
			rm(seqsPP);
			
		} # for (pp ...)
		phi[known] <- phi2;
		rm(phi2, known, K2);
	} else {
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# This approach assumes nothing about the 'values'
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# For each position
		for (pp in seq(length=P)) {
			
			
			# Get the nucleotides at this position for all sequences
			seqsPP <- seqs[,pp];
			
			allIdxs <- 1:length(seqsPP);
			for (bb in 1:ncol(rho)) {
				##      verbose && enter(verbose, sprintf("Factor #%d ('%s') of %d", bb, factors[bb], ncol(rho)));
				
				# Identify sequences with nucleotide 'bb' at position 'pp'.
				##      verbose && enter(verbose, "Identifying subset");
				subset <- whichVector(seqsPP == values[bb]);
				##      verbose && exit(verbose);
				
				# Add the nucleotide effect rho(pp,bb) to the probe-affinity
				idxs <- allIdxs[subset];
				phi[idxs] <- phi[idxs] + rho[pp,bb];
				
				# Skip already found cells
				allIdxs <- allIdxs[-subset];
				seqsPP <- seqsPP[-subset];
				
				##      verbose && exit(verbose);
			} # for (bb ...)
			
			rm(seqsPP);
			
		} # for (pp ...)
	}
	
	phi;
}




getProbePositionEffectDesignMatrix <-  function(seqs, verbose) {
	# Author: Hendrik Bengtsson; addapted by Djork Clevert
	# Argument 'seqs':
	P <- nchar(seqs);
	P <- unique(P);
	if (length(P) != 1) {
		stop("Argument 'seqs' contains sequences of different lengths: ",
				paste(head(sort(P)), collapse=", "));
	}
	
	K <- length(seqs);
	P <- nchar(seqs[1]);
	
	seqs <- paste(seqs, collapse="");
	seqs <- charToRaw(seqs);
	map <- as.raw(0:4);
	names(map) <- c(NA, "A", "C", "G", "T");
	values <- map[2:5];
	from <- charToRaw(paste(names(values), collapse=""));
	for (kk in seq(along=values)) {
		idxs <- whichVector(seqs == from[kk]);
		seqs[idxs] <- values[kk];
	}
	seqs[seqs > map[length(map)]] <- as.raw(0);
	seqs <- matrix(seqs, nrow=K, ncol=P, byrow=TRUE);
	gc <- gc();
	
	intercept=TRUE
	
	df=5
	
	P <- ncol(seqs);  # Number of positions in sequences
	
	B <- splines::ns(1:P, df=df);
	
	K <- nrow(seqs);  # Number of sequences
	
	df <- ncol(B);    # Number of basis vectors
	
	# Exclude NA factors and last factor
	map <- as.raw(0:4);
	names(map) <- c(NA, "A", "C", "G", "T"); 
	
	# Sanity check
	
	factors <- map[seq(from=2, to=length(map)-1)];
	
	L <- df*length(factors);
	if (intercept)
		L <- L + 1;
	
	# Allocate an KxL prediction matrix in the model
	
	X <- matrix(0, nrow=K, ncol=L);
	
	
	# Intercept?
	if (intercept) {
		X[,1] <- 1;
		gc <- gc();
	}
	
	for (bb in seq(along=factors)) {
		
		# For every position in the sequences
		for (pp in 1:P) {
			
			# Identify sequences with factor 'bb' in position 'pp'
			idxs <- whichVector(seqs[,pp] == factors[bb]);
			
			# For every dimension in the base vector
			for (jj in 1:df) {
				cc <- 1 + df*(bb-1) + jj;
				X[idxs,cc] <- X[idxs,cc] + B[pp,jj];
			}
			rm(idxs);
			
		} # for (pp ...)
		
		gc <- gc();
		
	} # for (bb ...)
	
	rm(seqs);
	gc <- gc();  
	
	
	res <- list(X=X, map=map, factors=factors, B=B);
	class(res) <- "ProbePositionEffectDesignMatrix";
	
	res;
	
}

fit_model <- function(X,y){
	# Author: Hendrik Bengtsson; addapted by Djork Clevert
	B <- X$B
	map <- X$map
	factors <- X$factors
	mX <- X$X;
	xtx <- crossprod(mX); 
	xty <- crossprod(mX, y);
	coefs <- solve(xtx, xty);
	coefs <- as.vector(coefs);
	params <- list();
	intercept <- TRUE;
	if (intercept) {
		params$intercept <- coefs[1];
		coefs <- coefs[-1];
	}
	df <- length(coefs)/length(factors);
	idxs <- 1:df;
	for (kk in seq(along=factors)) {
		key <- names(factors)[kk];
		if (is.null(key)) {
			key <- sprintf("factor%02d", kk);
		}
		params[[key]] <- coefs[idxs];
		coefs <- coefs[-idxs];
	}
	fit <- list(params=params, map=map, B=B, algorithm="solve");
	class(fit) <- "ProbePositionEffects";
	rm(mX, B)
	#gc()
	return(fit)
}


predict_model <- function(fit, seqs, y){
	x <- -predict_effect(fit,seqs) + y;
	return(x)
}




whichVector <- function(x, na.rm=TRUE, use.names=TRUE, ...) {                                                                                                                             
	if (!is.vector(x)) {                                                                                                                                                                                           
		stop("Argument 'x' is not a vector: ", class(x)[1]);                                                                                                                                                         
	}                                                                                                                                                                                                              
	
	idxs <- seq_along(x);                                                                                                                                                                                          
	
	# Identify TRUE and NA elements                                                                                                                                                                                
	idxs <- idxs[x];                                                                                                                                                                                               
	
	# Remove missing values?                                                                                                                                                                                       
	if (na.rm) {                                                                                                                                                                                                   
		idxs <- idxs[!is.na(idxs)];                                                                                                                                                                                  
	}                                                                                                                                                                                                              
	
	# Use names                                                                                                                                                                                                    
	if (use.names) {                                                                                                                                                                                               
		names(idxs) <- names(x)[idxs];                                                                                                                                                                               
	}                                                                                                                                                                                                              
	
	idxs;                                                                                                                                                                                                          
} 
