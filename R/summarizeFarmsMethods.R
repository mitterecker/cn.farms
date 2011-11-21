#' Lists methods for possible FARMS summarization
#' 
#' Possible FARMS summarization
#' 
#' @return Returns a data frame with all possible FARMS calls.
#' @author Djork-Arne Clevert \email{okko@@clevert.de} and 
#' Andreas Mitterecker \email{mitterecker@@bioinf.jku.at}
#' @examples
#' summarizeFarmsMethods()
#' @export
summarizeFarmsMethods <- function() {
    pattern <- "summarizeFarms"
    func <- ls("package:cn.farms")
    funcSel <- func[grep(pattern, func)]
    funcSel <- funcSel[-which(funcSel == "summarizeFarmsMethods")]
    funcCall <- gsub(pattern, "", funcSel)
    funcSel <- paste(funcSel, "()", sep = "")
    return(data.frame("call"=funcCall, "assignedFunction"=funcSel))
}
