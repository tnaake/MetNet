#' @name combineStructuralStatistical
#' @aliases combineStructuralStatistical
#' @title Combine structural and statistical adjacency matrix
#' @description The function \code{combineStructuralStatistical} takes as 
#' input the structural and statistical adjacency matrix, created in former 
#' steps, adds them together and will report a connection between metabolites
#' in the returned when the sum exceeds the \code{threshold} . 
#' \code{combineStructuralStatistical} returns this consensus matrix supported
#' by the structural and statistical adjacency matrices.
#' @usage combineStructuralStatistical(structure, statistical, threshold = 1)
#' @param structure matrix containing structural adjacency matrix
#' @param statistical matrix containing statistical adjacency matrix
#' @param threshold numeric, threshold value to be applied to define a 
#' connection as present 
#' @details The matrices will be added and a unweighted connection will 
#' be reported when the value exceeds a certain value. 
#' @return a matrix containing the consensus adjacency matrix as described 
#' above harbouring connections reported by the structual and 
#' statistcal adjacency matrices. 
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
#' struct_adj <- createStructuralAdjacency(x_test, functional_groups, ppm = 5)
#' stat_adj <- createStatisticalAdjacency(x_test, 
#'     model = c("pearson", "spearman", "bayes"), 
#'     correlation_adjust = "bonferroni")
#' combineStructuralStatistical(struct_adj[[1]], stat_adj)
#' @export
combineStructuralStatistical <- function(structure, statistical, 
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
    
    ## sum the matrices
    consensus_mat <- structure + statistical
    consensus_mat <- ifelse(consensus_mat > threshold, 1, 0)
    return(consensus_mat)
}
