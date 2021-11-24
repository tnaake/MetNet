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
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay `assay<-` assays assayNames
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
#' @section Accessors:
#' 
#' - The `AdjacencyMatrix` class extends the
#'   [SummarizedExperiment::SummarizedExperiment] class and inherits
#'   all its accessors and replacement methods.
#'
#' - The `type` accessor returns the `type` (`"structural"`, `"statistical"`, 
#'   `"combine"`) slot.
#'   
#' - The `directed` accessor returns the `directed` (`logical` of length 1)
#'   slot.
#'   
#' - The `thresholded` accessor returns the `thresholded` 
#' (`logical` of length 1) slot.
##'
#' @details
#' `adj_l` is a list of adjacency matrices. The adjacency matrices have
#' identical dimensions and `dimnames` and each adjacency matrix has the 
#' same number of columns and rows and identical `rownames` and `colnames`.
#' `rowData` will be also used for the `colData` slot (since the `rownames` 
#' and `colnames` are identical).
#'
#' @param adj_l `list` of adjacency matrices
#' 
#' @param rowData `data.frame`, containing information on the features
#' 
#' @param type `character`, either `"structural"`, `"statistical"`, or 
#' `"combine"`
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
#' transformation <- matrix("", ncol = 10, nrow = 10)
#' mass_difference <- matrix("", ncol = 10, nrow = 10)
#' 
#' rownames(binary) <- rownames(transformation) <- rownames(mass_difference) <- paste("feature", 1:10)
#' colnames(binary) <- rownames(transformation) <- rownames(mass_difference) <- paste("feature", 1:10)
#' 
#' binary[5, 4] <- 1
#' transformation[5, 4] <- "glucose addition"
#' mass_difference[5, 4] <- "162"
#' 
#' ## create adj_l and rowData
#' adj_l <- list(binary = binary, transformation = transformation, 
#'     mass_difference = mass_difference)
#' rowData <- DataFrame(features = rownames(binary),
#'     row.names = rownames(binary))
#' 
#' AdjacencyMatrix(adj_l = adj_l, rowData = rowData, type = "structural", 
#'     directed = TRUE, thresholded = FALSE)
#'
#' @export
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assays 
#' @importFrom methods new is
AdjacencyMatrix <- function(adj_l, rowData,
    type = c("structural", "statistical", "combine"),
    directed = c(TRUE, FALSE), thresholded = c(TRUE, FALSE)) {

    type <- match.arg(type)

    se <- SummarizedExperiment::SummarizedExperiment(adj_l, 
        rowData = rowData, colData = rowData)
    .AdjacencyMatrix(se, type = type, directed = directed, 
        thresholded = thresholded)
}

### 
### Validity
###
#' @importFrom S4Vectors setValidity2
#' @importFrom SummarizedExperiment assay assays assayNames
#' @importFrom methods validObject
setValidity2("AdjacencyMatrix", function(object) {
    msg <- NULL
    
    type_obj <- object@type
    
    a <- SummarizedExperiment::assays(object)[[1]]
    if (nrow(a) != ncol(a)) {
        msg <- c(msg, "nrow not equal ncol for assays")
    }
    
    ## valid assay names for structural, statistical, and combine
    valid_names_struct <- "binary"
    valid_names_stat <- c(
        "lasso_coef", "randomForest_coef", "clr_coef", "aracne_coef",
        "pearson_coef", "pearson_pvalue", "spearman_coef", "spearman_pvalue",
        "pearson_partial_coef", "pearson_partial_pvalue",
        "spearman_partial_coef", "spearman_partial_pvalue",
        "pearson_semipartial_coef", "pearson_semipartial_pvalue",
        "spearman_semipartial_coef", "spearman_semipartial_pvalue",
        "ggm_coef", "ggm_pvalue",
        "bayes_coef", "consensus")
    valid_names_combine <- "combine_binary"

    ## check if colnames, rownames for one assay are the same
    .assays_have_identical_colnames_rownames(object)

    ## check if colnames and rownames across assays are the same
    .assays_have_identical_dimnames(object)
    
    ############################################################################
    ############################################################################
    
    if (type_obj == "structural") {
        
        .nms <- SummarizedExperiment::assayNames(object)
        
        if (!all(valid_names_struct %in% .nms))
            msg <- c(msg, "assay names must contain 'binary'")

        .obj <- SummarizedExperiment::assay(object, "binary")
        if (!is.numeric(.obj))
            msg <- c(msg, "slot 'binary' must be numeric")
        
        .nms_cut <- .nms[!.nms %in% "binary"]
        
        msg_l <- lapply(.nms_cut, function(.nms_i) {
            .obj <- SummarizedExperiment::assay(object, .nms_i)
            if (!is.character(.obj))
                sprintf("slot '%s' must be character", .nms_i)
        })
        msg <- c(msg, unlist(msg_l))
    }

    if (type_obj == "statistical") {
        
        .nms <- SummarizedExperiment::assayNames(object)
        
        # if (type_obj == "combine")
        #     .nms <- .nms[!(.nms %in% c("binary", "transformation", "mass_difference",
        #         "combine_binary", "combine_transformation", "combine_mass_difference"))]
        
        if (any(!(.nms %in% valid_names_stat)))
            msg <- c(msg, "assay names contain invalid entries")
        
        assay_num <- lapply(seq_along(.nms),
                                    function(i) is.numeric(assay(object, i)))
        if (!all(unlist(assay_num)))
            msg <- c(msg, "assays must be all numeric")
    }
    
    if (type_obj == "combine") {
        
        .nms <- SummarizedExperiment::assayNames(object)
        
        ## check for columns binary and combine binary 
        ## (these should be numeric)
        if (!all(valid_names_struct %in% .nms))
            msg <- c(msg, "assay names must contain 'binary'")
        
        .obj <- SummarizedExperiment::assay(object, "binary")
        if (!is.numeric(.obj))
            msg <- c(msg, "slot 'binary' must be numeric")

        if (!all(valid_names_combine %in% .nms))
            msg <- c(msg, "assay names must contain 'combine_binary'")

        .obj <- SummarizedExperiment::assay(object, "combine_binary")
        if (!is.numeric(.obj))
            msg <- c(msg, "slot 'combine_binary' must be numeric")

        ## check for columns that originate from structural (expect for *binary)
        ## these should be character
        .nms_cut <- .nms[!.nms %in% 
            c("binary", "combine_binary", valid_names_stat)]
        
        msg_l <- lapply(.nms_cut, function(.nms_cut_i) {
            .obj <- SummarizedExperiment::assay(object, .nms_cut_i)
            if (!is.character(.obj))
                sprintf("slot '%s' must be character", .nms_cut_i)
        })
        msg <- c(msg, unlist(msg_l))
        
        ## check for columns that originate from statistical 
        ## (these should be numeric)
        .nms_cut <- .nms[.nms %in% valid_names_stat]
        
        msg_l <- lapply(.nms_cut, function(.nms_cut_i) {
            .obj <- SummarizedExperiment::assay(object, .nms_cut_i)
            if (!is.numeric(.obj))
                sprintf("slot '%s' must be numeric", .nms_cut_i)
        })
        msg <- c(msg, unlist(msg_l))
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
#' 
#' @importFrom S4Vectors getListElement
#' @importFrom SummarizedExperiment SummarizedExperiment
.assays_have_identical_dimnames <- function(object) {
    .assays <- SummarizedExperiment::assays(object)
    a_1 <- S4Vectors::getListElement(.assays, 1)
    ident <- vapply(seq_along(.assays),
        function(i) {
            a <- S4Vectors::getListElement(.assays, i)
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
#' 
#' @importFrom S4Vectors getListElement
#' @importFrom SummarizedExperiment SummarizedExperiment
.assays_have_identical_colnames_rownames <- function(object) {
    .assays <- SummarizedExperiment::assays(object)
    ident <- vapply(seq_along(.assays),
        function(i) {
            a <- S4Vectors::getListElement(.assays, i)
            identical(colnames(a), rownames(a))
        }, 
        logical(1)
    )
    all(ident) 
}


