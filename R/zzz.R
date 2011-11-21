.onLoad <- function(libname, pkgname) {
#    suppressMessages(library(Biobase, quietly = TRUE))
#    suppressMessages(library(ff, quietly = TRUE))
#    suppressMessages(library(oligo, quietly = TRUE))
#    suppressMessages(library(oligoClasses, quietly = TRUE))
#    suppressMessages(library(snowfall, quietly = TRUE))
#    suppressMessages(library(snow, quietly = TRUE))
#    suppressMessages(library(DBI, quietly = TRUE))
#    suppressMessages(library(affxparser, quietly = TRUE))
#    suppressMessages(library(DNAcopy, quietly = TRUE))
#    suppressMessages(library(preprocessCore, quietly = TRUE))
#    suppressMessages(library(lattice, quietly = TRUE))
#    suppressMessages(library(methods, quietly = TRUE))

	
	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui")
		addVigs2WinMenu("cn.farms")
	
	packageStartupMessage(
			"               _                                " ,"\n",
			"              | |                               " ,"\n",
			" __   _  _    | |  __,   ,_    _  _  _    ,     " ,"\n",
			"/    / |/ |   |/  /  |  /  |  / |/ |/ |  / \\_  " ,"\n",
			"\\___/  |  |_/o|__/\\_/|_/   |_/  |  |  |_/ \\/ " ,"\n",
			"              |\\                               " ,"\n",
			"              |/                                " ,"\n")
	
	packageStartupMessage("Citation: D.-A. Clevert et al.," ,"\n",
			"cn.FARMS: a latent variable model to detect copy number ", 
			"variations in microarray data with a low false discovery rate","\n",
			"Nucleic Acids Research, 2011.","\n",
			"BibTex: enter 'toBibtex(citation(\"cn.farms\"))'","\n\n",
			"Homepage: http://www.bioinf.jku.at/software/cnfarms/cnfarms.html","\n\n",
			"cn.farms Package Version ", 
            utils::packageDescription("cn.farms")$Version, "\n")
	#suppressMessages()
    
    pi <- utils::packageDescription(pkgname);
      packageStartupMessage(pkgname, " v", pi$Version, " (", 
                pi$Date, ") successfully loaded. See ?", pkgname, " for help."); 
    
}

.onUnload <- function(libpath) {
	library.dynam.unload("cn.farms", libpath)
}
