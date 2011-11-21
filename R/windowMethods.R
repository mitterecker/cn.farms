#' Lists methods for possible window methods
#' 
#' Function to list how neighbouring positions can be combined.
#' 
#' @return Returns a data frame with all possible methods. 
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @export 
#' @examples
#' summarizeWindowMethods()
summarizeWindowMethods <- function() {
    pattern <- "summarizeWindow"
    func <- ls("package:cn.farms")
    funcSel <- func[grep(pattern, func)]
    funcSel <- funcSel[-which(funcSel == "summarizeWindowMethods")]
    funcCall <- gsub(pattern, "", funcSel)
    funcSel <- paste(funcSel, "()", sep = "")
    return(data.frame("call"=funcCall, "assignedFunction"=funcSel))
}
