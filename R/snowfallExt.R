##*****************************************************************************
## Functions which extend Snow usage or implement some higher level usage.
## Wrappers for Snow functions are in snowWrappers.R
##
## Functions:
##    sfLoadLib   - Load library in cluster (path conversion, flags...)
##    sfSource    - Load source (path conversion...)
##
##    sfExport    - Export local and global objects to cluster
##    sfExportAll - Export all global objects (with given exception list)
##    sfRemove    - Remove objects from nodes
##    sfRemoveAll - Remove all objects from nodes (with given excpetion list)
##
##    sfCat       - Cat something on cluster
##    sfPrint     - Print something on cluster
##*****************************************************************************

##*****************************************************************************
## Load a library depending on sequential or parallel execution.
##
## Should behave most likely the build-in library() function.
##
## Running in sequential mode: normal "library()" command.
## Running in parallel mode: the library is loaded on ALL nodes.
##
## PARAMETERS: 
## RETURN:     Logical TRUE/FALSE on success. On failure, if noStopOnError is not
##             set, stop immidiately.
##*****************************************************************************
#' This function was taken from snowfall and edited due to some deprecated function calls. 
#' 
#' @param package name of the package. Check 'library' for details.
#' @param pos position in search path to load library.
#' @param lib.loc a character vector describing the location of the R library  trees to search through, or 'NULL'. Check 'library' for details.
#' @param character.only a logical indicating package can be assumed to be a character string. Check 'library' for details.
#' @param warn.conflicts warn on conflicts (see "library").
#' @param keep.source DEPRECATED (see "library").
#' @param verbose enable verbose messages.
#' @param version version of library to load (see "library").
#' @param stopOnError logical. 
#' @author snowfall
#' @return for more information see "library".
#' @export
#' @importFrom snowfall sfCat
#' @importFrom snowfall sfExport
cnLibrary <- function( package,
    pos = 2,
    lib.loc = NULL,
    character.only = FALSE,
    warn.conflicts = TRUE,
    keep.source = getOption("keep.source.pkgs"),
    verbose = getOption("verbose"),
    version,
    stopOnError = TRUE ) {
  snowfall:::sfCheck();
  
  ## Generate (global) names list with all parameters.
#  setVar( ".sfPars", list() )
  sfPars <- list()
  
  absFilePath <- function( file ) {
    ## If not starting with separator, path is most likely relative.
    ## Make it absolute then.
    ## On Windows absolute path can contain drive chars.
    if( .Platform$OS.type == "windows" ) {
      if( ( substr( file, 1, 1 ) != .Platform$file.sep ) &&
          ( substr( file, 2, 2 ) != ":" ) )
        file <- file.path( getwd(), file )
    }
    else
    if( substr( file, 1, 1 ) != .Platform$file.sep )
      file <- file.path( getwd(), file )
    
    return( file )
  }
  
  ## Help does not make sense.
  ##  if( !missing( help ) )
  ##    stop( "Help is not allowed in cnLibrary. Use 'library' instead." )
  
  if( !missing( package ) ) {
    if( character.only ) {
      if( is.character( package ) )
        sfPars$package <- package
      else
        stop( paste( "Package", package, "is no character string." ) )
    }
    else
      sfPars$package <- deparse( substitute( package ) )
  }
  
  ## package is now a string in any case.
  sfPars$character.only <- TRUE
  
  sfPars$pos            <- pos
  sfPars$lib.loc        <- lib.loc
  sfPars$warn.conflicts <- warn.conflicts
  #sfPars$keep.source    <- keep.source
  sfPars$verbose        <- verbose
  
  ## All libraries are loaded internally with logical.return.
  sfPars$logical.return <- TRUE
  
  if( !missing( version ) )
    sfPars$version <- version
  
  if( sfParallel() ) {
    ## On Nodes load location with absolute path.
    if( !is.null( sfPars$lib.loc ) )
      sfPars$lib.loc <- absFilePath( sfPars$lib.loc )
        
    ## Export to namespace.
    snowfall:::setVar( ".sfPars", sfPars )
    .sfPars <- sfPars
    
    ## Weird enough ".sfPars" need to be exported (else it would not be found
    ## on slave, although it is a parameter)
    sfExport( ".sfPars", local=FALSE, namespace="snowfall" )
    
    ## Load libs using require as Exception on nodes doesn't help us here.
    ## @todo Check on correct execution via logical.return
    ## @todo Exporting of .sfPars needed?
    result <- try( sfClusterEval( do.call( "library", .sfPars ) ) )
    ##    result <- try( sfClusterEval( library( .sfPars ) ) )
    
    if( inherits( result, "try-error" ) ||
        ( length( result ) != sfCpus() ) ||
        !all( snowfall:::checkTryErrorAny( result ) ) ||
        !all( unlist( result ) ) ) {
      if( stopOnError )
        stop( paste( "Stop: error loading library on slave(s):",
                .sfPars$package ) )
      else {
        warning( paste( "Error loading library on slave(s):", package ) )
        return( invisible( FALSE ) )
      }
    }
    else {
      ## Load message in slave logs.
      sfCat( paste( "Library", .sfPars$package, "loaded.\n" ) )
      
      ## Message in masterlog.
      message( paste( "Library", .sfPars$package, "loaded in cluster.\n" ) )
    }
  }
  
  result <- try( do.call( "library", sfPars ) )
  
  ## Remove global var from cluster (and local).
  ## Do before exception checks, as there might by a stop.
  sfRemove( ".sfPars" )
  
  if( inherits( result, "try-error" ) || !result ) {
    if( stopOnError ) {
      warning( paste( "Unable to load library:", package ) )
      return( invisible( FALSE ) )
    }
    else
      stop( paste( "Unable to load library:", package ) )
  }
  else {
    if( verbose )
      message( paste( "Library", package, "loaded.\n" ) )
    
    ## If logical return is requested here it comes.
    ## In clustermode the programm immidiately stops it a library couldn't be
    ## load on a slave. In sequentially mode it behaves like library().
    return( invisible( TRUE ) )
  }
}
