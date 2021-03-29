#' @name .AdjacencyMatrix
#'
#' @title Create S4 class AdjacencyMatrix
#'
#' @description
#' The class `AdjacencyMatrix` extends the `SummarizedExperiment` class. It 
#' will add the slots `type`, `directed`, and `thresholded`.
#' 
#' @details
#' The slot `type` is of type `"character"`, storing the type of the 
#' `"AdjacencyMatrix"`, i.e. `"structural"`, `"statistical"`, or `"combined"`.
#' The slot `directed` is of type `"logical"`, storing if the adjacency matrix
#' is directed or not.
#' The slot `thresholded` is of type `"logical"`, storing if the adjacency 
#' matrix was thresholded, e.g. if the functions `rtCorrection` or `threshold`
#' were applied on the `structural` or `statistical` `AdjacencyMatrix` objects.
#' 
#' If any of the `AdjacencyMatrix` objects passed to the `combine` function
#' was `directed = TRUE` or `thresholded = TRUEs` the `combine` 
#' `AdjacencyMatrix` object will be `directed = TRUE` or `thresholded = TRUE`.
#'
#' @return
#' class generator function for class `AdjacencyMatrix`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
.AdjacencyMatrix <- setClass("AdjacencyMatrix",
    contains = c("SummarizedExperiment"),
    slots = c(
        type = "character", ## structural, statistical, combined,
        directed = "logical",
        thresholded = "logical" ## thresholded, e.g. by rtCorrection or by threshold
    )
)

#' @name AdjacencyMatrix
#'
#' @title Wrapper to create an instance of S4 class AdjacencyMatrix
#'
#' @description
#' The function `AdjacencyMatrix` will create an object of type 
#' `AdjacencyMatrix`.
#'
#' @details
#' `adj_l` is a list of adjacency matrices. The adjacency matrices have
#' identical dimensions and `dimnames` and each adjacency matrix has the 
#' same number of columns and rows and identical `rownames` and `colnames`.
#' `rowData` will be also used for the `colData` slot (since the `rownames` 
#' and `colnames` are identical).
#'
#' @param adj_l `list` of adjacency matrices
#' 
#' @param rowData information on the features
#' 
#' @param type `character`, either `"structural"`, `"statistical"`, or 
#' `"combined"`
#' 
#' @param directed `logical`, if the adjacency matrix underlying the graph is
#' directed or undirected 
#' 
#' @param thresholded `logical`, if the functions `rtCorrection` or `threshold`
#' were applied on the `structural` or `statistical` `AdjacencyMatrix` objects
#'
#' @return
#' object of S4 class `AdjacencyMatrix`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' binary <- matrix(0, ncol = 10, nrow = 10)
#' type <- matrix("", ncol = 10, nrow = 10)
#' mass_differene("", ncol = 10, nrow = 10)
#' 
#' rownames(binary) <- rownames(type) <- rownames(mass_difference) <- 
#'     paste("feature", 1:10)
#' colnames(binary) <- rownames(type) <- rownames(mass_difference) <- 
#'     paste("feature", 1:10)
#' 
#' binary[5, 4] <- 1
#' type[5, 4] <- "glucose addition"
#' mass_difference[5, 4] <- "162"
#' 
#' ## create adj_l and rowData
#' adj_l <- list(binary = binary, type = type, 
#'     mass_difference = mass_difference)
#' rowData <- DataFrame(features = rownames(binary), 
#'     row.names = rownames(binary))
#'     
#' AdjacencyMatrix(adj_l = adj_l, rowData = rowData, type = structural", 
#'     directed = TRUE, thresholded = FALSE)
#'
#' @export
AdjacencyMatrix <- function(adj_l, rowData, 
    type = c("structural", "statistical", "combined"), 
    directed = c(TRUE, FALSE), thresholded = c(TRUE, FALSE)) {
    
    type <- match.arg(type)
    directed <- match.arg(directed)
    thresholded <- match.arg(thresholded)
    
    se <- SummarizedExperiment(adj_l, rowData = rowData, colData = rowData)
    .AdjacencyMatrix(se, type = type, directed = directed, 
        thresholded = thresholded)
}


### 
### Validity
###
setValidity2("AdjacencyMatrix", function(object) {
    msg <- NULL
    
    type_obj <- type(object)
    
    a <- assays(object)[[1]]
    if (nrow(a) != ncol(a)) {
        msg <- c(msg, "nrow not equal ncol for assays")
    }
    
    ## valid assay names for structural, statistical, and combine
    valid_names_struct <- c("binary", "transformation", "mass_difference")
    valid_names_stat <- c(
        "lasso_coef", "randomForest_coef", "clr_coef", "aracne_coef",
        "pearson_coef", "pearson_pvalue", "spearman_coef", "spearman_pvalue",
        "pearson_partial_coef", "pearson_partial_pvalue",
        "spearman_partial_coef", "spearman_partial_pvalue",
        "pearson_semipartial_coef", "pearson_semipartial_pvalue",
        "spearman_semipartial_coef", "spearman_semipartial_pvalue",
        "bayes_coef", "consensus")
    valid_names_combine <- c("combine_binary", "combine_transformation", 
        "combine_mass_difference")
    
    ## check if colnames, rownames for one assay are the same
    .assays_have_identical_colnames_rownames(object)    

    ## check if colnames and rownames across assays are the same
    .assays_have_identical_dimnames(object)
    
    ############################################################################
    ############################################################################
    
    if (type_obj == "structural" | type_obj == "combine") {
        
        if (!all(valid_names_struct %in% assayNames(object)))
            msg <- c(msg, "assay names must be 'binary', 'transformation', 'mass_difference'")
        
        .obj <- assay(object, "binary")
        if (!is.numeric(.obj))
            msg <- c(msg, "slot 'binary' must be numeric")
        
        .obj <- assay(object, "transformation")
        if (!is.character(.obj))
            msg <- c(msg, "slot 'transformation' must be character")
        
        .obj <- assay(object, "mass_difference")
        if (!is.character(.obj))
            msg <- c(msg, "slot 'mass_difference' must be character")
    }
    
    if (type_obj == "statistical" | type_obj == "combine") {
        
        if (any(!(assayNames(object) %in% valid_names_stat)))
            msg <- c(msg, "assay names contain invalid entries")
        
        .nms <- assayNames(object)
        
        if (type_obj == "combine")
            .nms <- .nms[!(.nms %in% c("binary", "type", "mass_difference"))]
        
        assay_num <- lapply(seq_along(.nms),
                                    function(i) is.numeric(assay(object, i)))
        if (!all(unlist(assay_num)))
            msg <- c(msg, "assays must be all numeric")
    }
    
    if (type_obj == "combine") {
        
        if (!all(valid_names_combine %in% assayNames(object)))
            msg <- c(msg, 
                "assay names must be 'comine_binary', 'combine_transformation', 'combine_mass_difference'")
        
        .obj <- assay(object, "combine_binary")
        if (!is.numeric(.obj))
            msg <- c(msg, "slot 'combine_binary' must be numeric")

        .obj <- assay(object, "combine_transformation")
        if (!is.character(.obj))
            msg <- c(msg, "slot 'combine_transformation' must be character")

        .obj <- assay(object, "combine_mass_difference")
        if (!is.character(.obj))
            msg <- c(msg, "slot 'combine_mass_difference' must be character")
    }

    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }

})

#' @name .assays_have_identical_dimnames
#'
#' @title Check if the assays in the `AdjacencyMatrix` object have
#' identical dimnames
#'
#' @description
#' The function will check if the assays in the `AdjacencyMatrix` object
#' have identical dimnames.
#'
#' @details
#' Helper function for validity check of `AdjacencyMatrix` objects.
#' 
#' @param object `AdjacencyMatrix` object
#'
#' @return `logical` of length 1
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
.assays_have_identical_dimnames <- function(object) {
    .assays <- assays(object)
    a_1 <- getListElement(.assays, 1)
    ident <- vapply(seq_along(.assays),
        function(i) {
            a <- getListElement(.assays, i)
            a_dimnames <- dimnames(a)
            identical(rownames(a_dimnames), dimnames(a_1))
            },
        logical(1)
    )
    all(ident)
}

#' @name .assays_have_identical_colnames_rownames
#'
#' @title Check if all the assays in the `AdjacencyMatrix` object have
#' identical colnames and rownames
#'
#' @description
#' The function will check if all the assays in the `AdjacencyMatrix` object
#' have identical colnames and rownames.
#'
#' @details
#' Helper function for validity check of `AdjacencyMatrix` objects.
#' 
#' @param object `AdjacencyMatrix` object
#'
#' @return `logical` of length 1
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
.assays_have_identical_colnames_rownames <- function(object) {
    .assays <- assays(object)
    ident <- vapply(seq_along(.assays),
        function(i) {
            a <- getListElement(.assays, i)
            identical(colnames(a), rownames(a))
        }, 
        logical(1)
    )
    all(ident) 
}


