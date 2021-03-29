#' @name combine
#'
#' @aliases combine
#'
#' @title Combine structural and statistical `AdjacencyMatrix` objects
#'
#' @description
#' The function `combine` takes as
#' input the structural and statistical `AdjacencyMatrix` objects, created in former
#' steps. It will access the assays `binary` and `consensus`, adds them 
#' together and will report a connection between metabolites
#' if the edge is present in both matrices. 
#' 
#' `combine` returns an `AdjacencyMatrix` containing this consensus matrix 
#' supported by the structural and statistical adjacency matrices (assays
#' `combine_binary`, `combine_transformation`, and `combine_mass_difference`.
#'
#' @param am_structural `AdjacencyMatrix` containing `numeric` structural 
#' adjacency matrices (assays `binary`, `transformation`, and 
#' `mass_difference`).
#'
#' @param am_statistical `AdjacencyMatrix` containing the assay `consensus`
#'
#' @details The matrices from the assays `binary` and `consensus` will be added 
#' and an unweighted connection will
#' be reported when the edges are respectively present in both `binary` and 
#' `consensus`.
#'
#' @return 
#' `AdjacencyMatrix` object containing the assays
#' `combine_binary` (`numeric` adjacency matrix), 
#' `combine_transformation` (`character` adjacency matrix), and 
#' `combine_mass_difference` (`character` adjacency matrix). 
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
#'      formula = transformation[, 2],
#'      mass = as.numeric(transformation[, 3]))
#'      
#' am_struct <- structural(x_test, transformation, ppm = 5)
#' am_stat <- statistical(x_test, model = c("pearson", "spearman"),
#'     correlation_adjust = "bonferroni")
#' 
#' am_stat <- threshold(am_stat, type = "top2", args = list(n = 10))
#' combine(am_structural = am_struct, am_statistical = am_stat)
#'
#' @export
combine <- function(am_structural, am_statistical) {

    if (!is(am_structural, "AdjacencyMatrix")) 
        stop("'am_structural' is not an 'AdjacencyMatrix' object")
    
    if (!validObject(am_structural)) 
        stop("'am_structural' must be a valid 'AdjacencyMatrix' object")
    
    if (!is(am_statistical, "AdjacencyMatrix")) 
        stop("'am_structural' is not an 'AdjacencyMatrix' object")
    
    if (!validObject(am_statistical)) 
        stop("'am_structural' must be a valid 'AdjacencyMatrix' object")
    
    if (!("consensus" %in% assayNames(am_statistical))) 
        stop("'am_statistical' must contain assay 'consensus'")
    
    ## create list to store results
    res <- list()

    ## create the first entry of the list
    ## sum the matrices structural and statistical, if the value is above
    ## threshold then assign 1, otherwise 0
    mat_bin <- assay(am_structural, "binary") + assay(am_statistical, "consensus")
    mat_bin <- ifelse(mat_bin > 1, 1, 0)
    
    ## create the second entry of the list
    ## if element in mat_bin is equal to 1, take the element in 
    ## the assay "transformation" (the type of link), otherwise ""
    mat_type <- ifelse(mat_bin == 1, assay(am_structural, "transformation"), "")
    
    ## create the third entry of the list
    ## if element in mat_bin is equal to 1, take the element in 
    ## the assay "mass_difference" (the mass difference of link), otherwise ""
    mat_mass <- ifelse(mat_bin == 1, assay(am_structural, "mass_difference"), "")
    
    ## create the AdjacencyMatrix object
    l <- list(combine_binary = mat_bin, combine_transformation = mat_type, 
        combine_mass_difference = mat_mass)
    rD <- DataFrame(names = rownames(mat_bin), row.names = rownames(mat_bin))
    
    directed <- if (directed(am_structural) | directed(am_statistical)) TRUE else FALSE
    thresholded <- if (thresholded(am_structural) | thresholded(am_statistical)) TRUE else FALSE
    
    am <- AdjacencyMatrix(l, rowData = rD, type = "combine", 
        directed = directed, thresholded = thresholded)
    
    return(am)
    
}
