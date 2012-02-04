#' Creation of annotation files
#'
#' Annotation files for cn.farms are created
#' 
#' @aliases createAnnotation
#' @note The annotation files used for cn.farms will be placed in the current 
#' work directory under annotations. 
#' @param filenames An absolute path of the CEL files to process.
#' @param annotation Optional parameter stating the annotation from a pd-mapping.
#' @param annotDir Optional parameter stating where the annotation should go.
#' @param checks States if sanity checks should be done.
#' @return \code{NULL}
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @useDynLib cn.farms
#' @importMethodsFrom oligoClasses db
#' @importMethodsFrom DBI dbGetQuery
#' @importFrom affxparser readCelHeader
#' @importFrom oligo cleanPlatformName
#' @importFrom oligoClasses ldPath
#' @importFrom oligoClasses chromosome2integer
#' @export
#' @examples
#' \dontrun{
#' library("hapmapsnp6") 
#' celDir <- system.file("celFiles", package = "hapmapsnp6")
#' filenames <- dir(path = celDir, full.names = TRUE)
#' createAnnotation(filenames = filenames)
#' }
createAnnotation <- function(filenames = NULL, annotation = NULL, 
        annotDir = NULL, checks = TRUE) {
    
    if (checks) {
            ## check for correct CEL-files
            celfiles <- filenames[grep("\\.[cC][eE][lL]$", filenames)]
            if (length(celfiles) != length(filenames)) {
            message("It looks like that not all filenames are CEL-files!")
            message("Especially: ")
            message(paste("   ", setdiff(filenames, celfiles), sep = " ", 
                            collapse = "\n"))
            stop("Check CEL-files!")
            }
            
            ## check for equal mappings
            mappings <- sapply(filenames, 
                            function (x) affxparser::readCelHeader(x)$chiptype)
            if (length(unique(mappings)) != 1) {
                    print(mappings)
                    stop("Different mappings found!")
            }
    }
    
    ## check parameters
    if (is.null(filenames) & is.null(annotation) & is.null(annotDir)) {
        stop("Either provide celfile names or annotation string")
    } else if (!is.null(filenames) & is.null(annotation)) {
        mapping <- affxparser::readCelHeader(filenames[1])$chiptype
        pkgname <- oligo::cleanPlatformName(mapping)
        pkgname <- correctPkgname(pkgname)
        
    } else if (!is.null(annotation)){
        mapping <- annotation
        pkgname <- oligo::cleanPlatformName(mapping)
    }
    
    if (!is.element(pkgname, installed.packages()[,1])) {
        stop(paste("Package", pkgname, 
                        "from bioconductor.org not installed but needed!"))
    }
    
    require(pkgname, character.only = TRUE, quietly = TRUE)
    version <- installed.packages()[pkgname, "Version"]

    ## check annotation directory
    if (is.null(annotDir)) {
        annotDir <- file.path(getwd(), "annotation", pkgname, version)
        dir.create(annotDir, showWarnings = FALSE, recursive = TRUE)
    }
    if (length(dir(annotDir)) != 0) {
        cat(paste(Sys.time(), "|   Directory", annotDir, "not empty \n",
                        "                   |  ",
                        "Annotation probably already done \n"), sep = "")
        return(invisible())
    }
    
    knownPackages <- c(
            "pd.mapping50k.hind240",
            "pd.mapping50k.xba240",
            "pd.genomewidesnp.5",
            "pd.genomewidesnp.6", 
            "pd.mapping250k.nsp", 
            "pd.mapping250k.sty", 
            "pd.cytogenetics.array")
    
    if (!(pkgname %in% knownPackages)) {
        stop("Unknown annotation")
    } else {
        cat(paste(Sys.time(), "|   Reading annotation from package", 
                        pkgname, version, " \n"))
        cat(paste(Sys.time(), "|   Annotation will be saved in", 
                        annotDir, " \n"))
        
        ## featureSet
        sql <- "SELECT * FROM featureSet"
        tmp <- DBI::dbGetQuery(oligoClasses::db(get(pkgname)), sql)
        tmp$chrom <- oligoClasses::chromosome2integer(tmp$chrom)
        idx <- !is.na(tmp$chrom)
        tmp <- tmp[idx, ]
        
        if (pkgname == "pd.cytogenetics.array") {
            tmp$physical_pos <- tmp$position
            tmp$position <- NULL
        }
        
        featureSetFull <- tmp[order(tmp$chrom, tmp$physical_pos, tmp$fsetid), ]
        save(featureSetFull, file = file.path(annotDir, "featureSetFull.RData"))
        
        if (pkgname %in% c("pd.genomewidesnp.5", "pd.genomewidesnp.6")) {
            featureSet <- featureSetFull[, c("fsetid", "man_fsetid",  
                            "chrom", "physical_pos", "allele_a",
                            "allele_b")] ## "dbsnp_rs_id" missing
        } else if (pkgname == "pd.cytogenetics.array") { 
            featureSet <- featureSetFull
        } else {
            featureSet <- featureSetFull[, c("fsetid", "man_fsetid",  
                            "chrom", "physical_pos", "allele_a",
                            "allele_b")] ## "dbsnp_rs_id"
        }
        save(featureSet, file = file.path(annotDir, "featureSet.RData"))
        gc()
        rm(featureSetFull)
        
        ## pmfeature
        sql <- "SELECT * FROM pmfeature"
        pmfeatureTmp <- DBI::dbGetQuery(db(get(pkgname)), sql)
        idxTmp <- match(pmfeatureTmp$fsetid, featureSet$fsetid)
        idx <- order(idxTmp)
        pmfeature <- pmfeatureTmp[idx, ]
        save(pmfeature, file = file.path(annotDir, "pmfeature.RData"))
        rm(idxTmp, idx, pmfeatureTmp)
        gc()
        
        ## sequence
        if (pkgname != "pd.cytogenetics.array") {
            sql <- "SELECT * FROM sequence"    
            sequence <- DBI::dbGetQuery(db(get(pkgname)), sql)
            save(sequence, file = file.path(annotDir, "sequence.RData"))
            rm(sequence)
            gc()
        }
        
        cat(paste(Sys.time(), "|   SNP information done \n"))
        
        ## available only for newer Affymetrix arrays
        if (pkgname %in% c("pd.genomewidesnp.5", "pd.genomewidesnp.6", 
                "pd.cytogenetics.array")) {
            
            ## featureSetCNV
            sql <- "SELECT * FROM featureSetCNV"
            tmp <- DBI::dbGetQuery(db(get(pkgname)), sql)
            tmp$chrom <- oligoClasses::chromosome2integer(tmp$chrom)
            idx <- !is.na(tmp$chrom)
            tmp <- tmp[idx, ]
            
            if (pkgname == "pd.cytogenetics.array") {
                tmp$chrom_start <- tmp$position
                tmp$chrom_stop <- tmp$position
                tmp$position <- NULL
            }
            
            featureSetCNVFull <- tmp[order(tmp$chrom, 
                            tmp$chrom_start, tmp$fsetid), ]
            
            if (pkgname == "pd.cytogenetics.array") { 
                featureSetCNV <- featureSetCNVFull
            } else {
                featureSetCNV <- featureSetCNVFull[, c("fsetid", "man_fsetid", 
                                "chrom", "chrom_start", "chrom_stop")]
            }
            
            save(featureSetCNV, file = file.path(annotDir, "featureSetCNV.RData"))
            save(featureSetCNVFull, 
                    file = file.path(annotDir, "featureSetCNVFull.RData"))
            gc()
            
            sql <- "SELECT * FROM pmfeatureCNV"
            pmfeatureCNVTmp <- DBI::dbGetQuery(db(get(pkgname)), sql)
            idxTmp <- match(pmfeatureCNVTmp$fsetid, featureSetCNV$fsetid)
            idx <- order(idxTmp)
            pmfeatureCNV <- pmfeatureCNVTmp[idx, ]
            save(pmfeatureCNV, file = file.path(annotDir, "pmfeatureCNV.RData"))
            rm(idxTmp, idx, pmfeatureCNVTmp)
            gc()
            
            if (pkgname != "pd.cytogenetics.array") {
                sql <- "SELECT * FROM sequenceCNV"    
                sequenceCNV <- DBI::dbGetQuery(db(get(pkgname)), sql)
                save(sequenceCNV, file = file.path(annotDir, "sequenceCNV.RData"))
                gc()
            }
        }
        
        cat(paste(Sys.time(), "|   Non polymorphic information done \n"))
        
        ## for normalization
        tmp <- match(pmfeature$fsetid, featureSet$fsetid)
        pmfeatureAllele <- featureSet[tmp, ]
        rm(featureSet)
        
        idxOfAlleleA <- which(pmfeature[, "allele"] == 0)
        idxOfAlleleB <- which(pmfeature[, "allele"] == 1)
        alleleA <- pmfeatureAllele[idxOfAlleleA, "allele_a"]
        alleleB <- pmfeatureAllele[idxOfAlleleB, "allele_b"]
        
        idxOfStrandA <- which(pmfeature[, "strand"] == 0)
        idxOfStrandB <- which(pmfeature[, "strand"] == 1)
        
        pairs <- paste(alleleA, alleleB, sep = "")
        uniquePairs <- unique(pairs)
        
        save(pmfeature, uniquePairs, idxOfAlleleA, idxOfAlleleB, 
                idxOfStrandA, idxOfStrandB, pairs, 
                file = file.path(annotDir, "annotNormalization.RData"))
    }
    gc()
    cat(paste(Sys.time(), "|   Annotation processed \n"))
    invisible()
}
