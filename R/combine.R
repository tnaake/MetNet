load("/home/thomas/Documents/University/Master/MScArbeit/Metabolic_profiling/WOS_MeJA_0h_72h/CAMERA_complete.RData")

peaklist <- peaklist[sort(sample(x = 1:dim(peaklist)[1], size = 100)),]

#' @name combine_structural_statistical
#' @title Combine structural and statistical network
#' @description The function \code{combine_structural_statistical} takes as 
#' input the structural and statistical network, that were create in former 
#' steps, adds them together and will report a connection between metabolites
#' when the edge is reported in networks. \code{combine_structural_statistical}
#' returns this consensus matrix of the structural and statistical network.
#' @usage combine_structural_statistical(structure, statistical, treshold = 1)
#' @param structure matrix containing structural network
#' @param statistical matrix containing statistical netowrk
#' @details \code{combine_structural_statistical} will call internally the 
#' function \code{consensus_network} to create the consensus network of the 
#' statistical networks. 
#' @return a matrix containing the consensus network as described above 
#' harbouring connections reported by the structual network and statistcal 
#' networks. 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples combine_structural_statistical(structure, statistical, threshold)
#' @export
combine_structural_statistical <- function(structure, statistical, 
                    threshold = 1) {
    
    if (!is.matrix(structure) & !is.numeric(structure)) 
        stop("structure is not a numeric matrix")
    if (!is.matrix(statistical) & !is.numeric(statistical)) 
        stop("statistical is not a numeric matrix")
    if (!all(rownames(structure) == rownames(statistical))) 
        stop("rownames are not identical")
    if (!all(colnames(structure) == colnames(statistical))) 
        stop("colnames are not identical")
    if (!is.numeric(threshold)) stop("threshold is not numeric")
    ##l <- list(structure, statistical)
    ##consensus_mat <- consensus_network(statistical, ...)
    consensus_mat <- structure + statistical
    consensus_mat <- ifelse(consensus_mat > threshold, 1, 0)
    return(consensus_mat)
}


combine_structural_statistical(mat_fg[[1]], l, method = "central.graph")


combine_structural_statistical(mat_fg[[1]], l, method = "central.graph")