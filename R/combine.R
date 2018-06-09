#' @name combine_structural_statistical
#' @aliases combine_structural_statistical
#' @title Combine structural and statistical network
#' @description The function \code{combine\_structural\_statistical} takes as 
#' input the structural and statistical network, that were create in former 
#' steps, adds them together and will report a connection between metabolites
#' when the edge is reported in networks. 
#' \code{combine\_structural\_statistical} returns this consensus matrix of 
#' the structural and statistical network.
#' @usage combine_structural_statistical(structure, statistical, threshold = 1)
#' @param structure matrix containing structural network
#' @param statistical matrix containing statistical network
#' @param threshold numeric, threshold value to be applied to define a 
#' connection as present 
#' @details \code{combine\_structural\_statistical} will call internally the 
#' function \code{consensus\_network} to create the consensus network of the 
#' statistical networks. 
#' @return a matrix containing the consensus network as described above 
#' harbouring connections reported by the structual network and statistcal 
#' networks. 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x_test <- as.matrix(x_test)
#' functional_groups <- rbind(
#'     c("Hydroxylation (-H)", "O", "15.9949146221"),
#'     c("Malonyl group (-H2O)", "C3H2O3", "86.0003939305"),
#'     c("C6H10O6", "C6H10O6", "178.0477380536"),
#'     c("D-ribose (-H2O) (ribosylation)", "C5H8O4", "132.0422587452"),
#'     c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851"),
#'     c("Glucuronic acid (-H2O)", "C6H8O6", "176.0320879894"),
#'     c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"),
#'     c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945"))
#' functional_groups <- data.frame(group = functional_groups[,1],
#'                                 formula = functional_groups[,2],
#'                                 mass = as.numeric(functional_groups[,3]))
#' struct_net <- create_structural_network(x_test, functional_groups, ppm = 5)
#' stat_net <- create_statistical_network(x_test, 
#'     model = c("pearson", "spearman","bayes"), 
#'     correlation_adjust = "bonferroni")
#' combine_structural_statistical(struct_net[[1]], stat_net)
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