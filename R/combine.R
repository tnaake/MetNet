#' @name combine
#'
#' @aliases combine
#'
#' @title Combine structural and statistical adjacency matrix
#'
#' @description
#' The function `combine` takes as
#' input the structural and statistical adjacency matrix, created in former
#' steps, adds them together and will report a connection between metabolites
#' in the returned when the sum exceeds the `threshold`.
#' \code{combine} returns this consensus matrix supported
#' by the structural and statistical adjacency matrices.
#'
#' @param structural list containing `numeric` structural adjacency matrix in
#' the first entry and `character` structural adjanceny matrix in the second
#' entry
#'
#' @param statistical matrix containing `numeric` statistical adjacency matrix
#'
#' @param threshold numeric, threshold value to be applied to define a
#' connection as present
#'
#' @details The matrices will be added and a unweighted connection will
#' be reported when the value exceeds a certain value.
#'
#' @return `list`, in the first entry `matrix` of type `numeric`containing the
#' consensus adjacency matrix as described
#' above harbouring connections reported by the structual and
#' statistcal adjacency matrices. In the second entry a `matrix` of type
#' `character` the corresonding type/putative link at this position.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x_test <- as.matrix(x_test)
#' functional_groups <- rbind(
#'     c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"),
#'     c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851"),
#'     c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945"))
#' functional_groups <- data.frame(group = functional_groups[, 1],
#'      formula = functional_groups[, 2],
#'      mass = as.numeric(functional_groups[, 3]))
#' struct_adj <- structural(x_test, functional_groups, ppm = 5)
#' stat_adj_l <- statistical(x_test,
#'     model = c("pearson", "spearman"),
#'     correlation_adjust = "bonferroni")
#' stat_adj <- threshold(stat_adj_l, type = "top2", args = list(n = 10))
#' combine(struct_adj, stat_adj)
#'
#' @export
combine <- function(structural, statistical, threshold = 1) {

    if (!is.list(structural) | length(structural) != 2)
        stop("structural is not a list of length 2")

    if (!is.matrix(structural[[1]]) | !is.numeric(structural[[1]]))
        stop("structural[[1]] is not a numeric matrix")

    if (!is.matrix(structural[[2]]) | !is.character(structural[[2]]))
        stop("structural[[2]] is not a character matrix")

    if (!is.matrix(statistical) | !is.numeric(statistical))
        stop("statistical is not a numeric matrix")

    if (!all(rownames(structural[[1]]) == rownames(structural[[2]])))
        stop(c("rownames of structural[[1]] are not identical to rownames of ",
             "structural[[2]]"))

    if (!all(colnames(structural[[1]]) == colnames(structural[[2]])))
        stop(c("colnames of structural[[1]] are not identical to colnames of ",
                 "structural[[2]]"))

    if (!all(rownames(structural[[1]]) == rownames(statistical)))
        stop("rownames are not identical")

    if (!all(colnames(structural[[1]]) == colnames(statistical)))
        stop("colnames are not identical")

    if (!is.numeric(threshold)) stop("threshold is not numeric")

    ## create list to store results
    res <- list()

    ## create the first entry of the list
    ## sum the matrices structural and statistical, if the value is above
    ## threshold then assign 1, otherwise 0
    cons_num <- structural[[1]] + statistical
    cons_num <- ifelse(cons_num > threshold, 1, 0)

    ## create the second entry of the list
    ## if element in cons_num is equal to 1, take the element in structural[[2]]
    ## (the type of link), otherwise ""
    cons_char <- ifelse(cons_num == 1, structural[[2]], "")

    ## assign to list
    res[[1]] <- cons_num
    res[[2]] <- cons_char

    return(res)
}