.onLoad <- function(libname, pkgname) {
  
  pi <- utils::packageDescription(pkgname)
  
  packageStartupMessage(
      "               _                                " ,"\n",
      "              | |                               " ,"\n",
      " __   _  _    | |  __,   ,_    _  _  _    ,     " ,"\n",
      "/    / |/ |   |/  /  |  /  |  / |/ |/ |  / \\_  " ,"\n",
      "\\___/  |  |_/o|__/\\_/|_/   |_/  |  |  |_/ \\/ " ,"\n",
      "              |\\                               " ,"\n",
      "              |/                                " ,"\n",
      "Citation: D.-A. Clevert et al.," ,"\n",
      "cn.FARMS: a latent variable model to detect copy number ", 
      "variations in microarray data with a low false discovery rate","\n",
      "Nucleic Acids Research, 2011.","\n",
      "BibTex: enter 'toBibtex(citation(\"cn.farms\"))'","\n\n",
      "Homepage: http://www.bioinf.jku.at/software/cnfarms/cnfarms.html","\n\n",
      "cn.farms Package Version ", 
      utils::packageDescription("cn.farms")$Version, "\n",
      pkgname, " v", pi$Version, " (", 
      pi$Date, ") successfully loaded.") 
  
}

.onUnload <- function(libpath) {
  library.dynam.unload("cn.farms", libpath)
}
