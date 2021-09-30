#' @name combine
#'
#' @aliases combine
#'
#' @title Combine structural and statistical `AdjacencyMatrix` objects
#'
#' @description
#' The function `combine` takes as nput the structural and statistical 
#' `AdjacencyMatrix` objects, created in former
#' steps. It will access the assays `binary` and `consensus`, adds them 
#' together and will report a connection between metabolites
#' if the edge is present in both matrices. 
#' 
#' `combine` returns an `AdjacencyMatrix` containing this consensus matrix 
#' supported by the structural and statistical adjacency matrices (assay
#' `combine_binary`). The `AdjacencyMatrix` object furthermore contains the 
#' assays from the statistical `AdjacencyMatrix` and the 
#' combined assays from the structural `AdjacencyMatrix`, e.g. if the 
#' structural `AdjacencyMatrix` has the assays `group` and `mass`, the 
#' combine `AdjacencyMatrix` object will contain the assays `combine_group` and
#' `combine_mass` that have support from the structural and statistical 
#' `AdjacencyMatrix` object.
#'
#' @param am_structural `AdjacencyMatrix` containing the `numeric` structural 
#' adjacency matrix (assay `binary`) and other `character` structural adjacency 
#' matrices (e.g. `group` and  `mass`).
#'
#' @param am_statistical `AdjacencyMatrix` containing the assay `consensus` and
#' other `numeric` adjacency matrices depending on the chosen statistical 
#' models
#'
#' @details The matrices from the assays `binary` and `consensus` will be added 
#' and an unweighted connection will
#' be reported when the edges are respectively present in both `binary` and 
#' `consensus`.
#'
#' @return 
#' `AdjacencyMatrix` object containing the assays
#' `combine_binary` (`numeric` adjacency matrix), and the combined matrices
#' derived from the structural `AdjacencyMatrix` (`character` adjacency 
#' matrices).
#' 
#' The `AdjacencyMatrix` object will also contain all other assays contained
#' in `am_structural` and `am_statistical`.
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x_test <- as.matrix(x_test)
#' transformation <- rbind(
#'     c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"),
#'     c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851"),
#'     c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945"))
#' transformation <- data.frame(group = transformation[, 1],
#'     formula = transformation[, 2],
#'     mass = as.numeric(transformation[, 3]))
#'      
#' ## create AdjacencyMatrix object of type structural
#' am_struct <- structural(x_test, transformation, var = c("group", "mass"), 
#'     ppm = 10)
#' 
#' ## create AdjacencyMatrix object of type statistical
#' x_test_cut <- as.matrix(x_test[, -c(1:2)])
#' am_stat <- statistical(x_test_cut, model = c("pearson", "spearman"),
#'     correlation_adjust = "bonferroni")
#' am_stat <- threshold(am_stat, type = "top2", args = list(n = 10))
#' 
#' ## combine
#' combine(am_structural = am_struct, am_statistical = am_stat)
#'
#' @importFrom SummarizedExperiment assay assayNames rowData
#' @export
combine <- function(am_structural, am_statistical) {

    if (!is(am_structural, "AdjacencyMatrix")) 
        stop("'am_structural' is not an 'AdjacencyMatrix' object")
    
    if (!validObject(am_structural)) 
        stop("'am_structural' must be a valid 'AdjacencyMatrix' object")
    
    if (!is(am_statistical, "AdjacencyMatrix")) 
        stop("'am_statistical' is not an 'AdjacencyMatrix' object")
    
    if (!validObject(am_statistical)) 
        stop("'am_statistical' must be a valid 'AdjacencyMatrix' object")
    
    if (!("consensus" %in% SummarizedExperiment::assayNames(am_statistical))) 
        stop("'am_statistical' must contain assay 'consensus'")
    
    if (!all(names(am_structural) == names(am_statistical)))
        stop("names of 'am_structural' do not match names of 'am_statistical'")

    ## create the first entry of the list
    ## sum the matrices structural and statistical, if the value is above
    ## threshold then assign 1, otherwise 0
    l_comb_binary <- SummarizedExperiment::assay(am_structural, "binary") + 
        SummarizedExperiment::assay(am_statistical, "consensus")
    l_comb_binary <- ifelse(l_comb_binary > 1, 1, 0)
    
    ## create the entries of the list
    ## if element in mat_bin is equal to 1, take the element in 
    ## the assay .nms_i, otherwise take ""
    .nms <- SummarizedExperiment::assayNames(am_structural)
    .nms_cut <- .nms[!.nms %in% "binary"]
    
    l_comb <- lapply(.nms_cut, function(.nms_cut_i) {
        l_comb_i <- ifelse(l_comb_binary == 1, 
            yes = SummarizedExperiment::assay(am_structural, .nms_cut_i), 
            no = "")
        l_comb_i
    })
    names(l_comb) <- .nms_cut
    
    ## add the l_comb_binary to l_comb
    l_comb <- c(list(binary = l_comb_binary), l_comb)
    names(l_comb) <- paste("combine_", names(l_comb), sep = "")
    
   
    ## create the AdjacencyMatrix object, start with the structural 
    ## AdjacencyMatrix, then the statistical AdjacencyMatrix and add 
    ## the l_comb
    .nms <- SummarizedExperiment::assayNames(am_structural)
    l_structural <- lapply(.nms, function(.nms_i) 
        SummarizedExperiment::assay(am_structural, .nms_i))
    names(l_structural) <- .nms
    
    .nms <- SummarizedExperiment::assayNames(am_statistical)
    l_statistical <- lapply(.nms, function(.nms_i) 
        SummarizedExperiment::assay(am_statistical, .nms_i))
    names(l_statistical) <- .nms
    
    ## concatenate the lists from structural, statistical and combine
    l <- c(l_structural, l_statistical, l_comb)
    
    ## directed slot
    if (directed(am_structural) | directed(am_statistical)) {
        directed <- TRUE 
    } else {
        directed <- FALSE
    }
    
    ## thresholded slot
    if (am_structural@thresholded | am_statistical@thresholded) {
        thresholded <- TRUE 
    } else {
        thresholded <- FALSE
    }
    
    ## create the rowData
    rD <- SummarizedExperiment::rowData(am_statistical)
    
    ## finally, combine the information and create the AdjacencyMatrix object
    am <- AdjacencyMatrix(l, rowData = rD, type = "combine", 
        directed = directed, thresholded = thresholded)
    
    return(am)
}

