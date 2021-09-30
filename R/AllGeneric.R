##' @title Placeholder for generics functions documentation
##'
##' @name AllGenerics
##' 
##' @rdname AllGenerics
NULL


###
### Accessors
###

#' @name AdjacencyMatrix-class
#' 
#' @rdname AdjacencyMatrix-class
#' 
#' @title Methods for `AdjacencyMatrix` objects
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL

#' @rdname AdjacencyMatrix-class
#' 
#' @description 
#' `length` returns the length of an `AdjacencyMatrix` object
#' (number of rows of an assay). `length` returns a `numeric` of length 1.
#' 
#' @param x an instance of class `AdjacencyMatrix`
#' 
#' @exportMethod length
setMethod("length", "AdjacencyMatrix",
    function(x) nrow(assay(x, 1))
)

 
#' @rdname AdjacencyMatrix-class
#' 
#' @aliases dim
#' 
#' @description 
#' `dim` returns the length of an `AdjacencyMatrix` object
#' (number of rows of an assay, number of cols of an assay). `dim` returns
#' a `numeric` of length 2.
#' 
#' @param x an instance of class `AdjacencyMatrix`
#' 
#' @exportMethod dim
setMethod("dim", "AdjacencyMatrix",
    function(x) c(nrow(assay(x, 1)), ncol(assay(x, 1)))
)

setGeneric("type",
    function(x) standardGeneric("type"))

#' @rdname AdjacencyMatrix-class
#' 
#' @aliases type
#' 
#' @description 
#' `type` will return the type of an `AdjacencyMatrix` (`statistical`, 
#' `structural` or `combine`). `type` returns a `character` of length 1
#' 
#' @param x instance of class `AdjacencyMatrix`
#' 
#' @exportMethod type
setMethod("type", "AdjacencyMatrix",
    function(x) x@type)

setGeneric("directed",
    function(object) standardGeneric("directed"))

#' @rdname AdjacencyMatrix-class
#' 
#' @aliases directed
#' 
#' @description 
#' `directed` returns the information on directed of an `AdjacencyMatrix`, i.e.
#' if the underlying graph is directed or undirected.
#' `directed` returns `logical` of length 1.
#' 
#' @param object instance of class `AdjacencyMatrix`
#' 
#' @exportMethod directed
setMethod("directed", "AdjacencyMatrix",
    function(object) object@directed)

setGeneric("thresholded",
    function(object) standardGeneric("thresholded"))

#' @rdname AdjacencyMatrix-class
#' 
#' @aliases thresholded
#' 
#' @description 
#' `thresholded` returns the information if the adjacency matrix is 
#' thresholded, i.e. if the function `rtCorrection` or `threshold` was applied
#' to the `AdjacencyMatrix` object. `thresholded` returns a `logical` of length 
#' 1.
#' 
#' @param object instance of class `AdjacencyMatrix`
#' 
#' @exportMethod thresholded
setMethod("thresholded", "AdjacencyMatrix",
    function(object) object@thresholded)



###
### Show
###
#' @rdname AdjacencyMatrix-class
#'
#' @aliases show
#' 
#' @description 
#' `show` prints summary information on an object of class 
#' `AdjacencyMatrix`.
#' 
#' @param object instance of class `AdjacencyMatrix`
#' 
#' @exportMethod show
#' 
#' @importFrom S4Vectors coolcat
#' @importFrom methods show
setMethod("show", "AdjacencyMatrix",
    function(object) {
          
        cat("class:", class(object), "\n")
        cat("dim:", dim(object), "\n")
              
        ## type()
        .type <- object@type
        if (is.null(.type))
            .type <- character(length(object@type))
        coolcat("type(%d): %s\n", .type)
        
        ## directed()
        .directed <- object@directed
        if (is.null(.directed))
            .directed <- character(length(object@directed))
        coolcat("directed(%d): %s\n", .directed)
              
        ## threshold()
        .threshold <- object@thresholded
        if (is.null(.threshold))
            .threshold <- character(length(object@thresholded))
        coolcat("thresholded(%d): %s\n", .threshold)
              
        ## assay()
        .nms <- assayNames(object)
        if (is.null(.nms))
            .nms <- character(length(assay(object, withDimnames=FALSE)))
        coolcat("assay(%d): %s\n", .nms)
              
        ## rownames()
        rownames <- rownames(object)
        if (!is.null(rownames)) coolcat("rownames(%d): %s\n", rownames)
        else cat("rownames: NULL\n")
              
        ## colnames()
        colnames <- colnames(object)
        if (!is.null(colnames)) coolcat("colnames(%d): %s\n", colnames)
        else cat("colnames: NULL\n")
    }
)

###
### as.data.frame
###

#' @rdname AdjacencyMatrix-class
#'
#' @aliases as.data.frame
#'
#' @description 
#' `as.data.frame` returns the adjacency matrices (stored in the `assays` slot)
#' and returns information on the nodes and the associated information on 
#' edges as a data frame. `as.data.frame` returns a `data.frame`.
#' 
#' @param x instance of class `AdjacencyMatrix` 
#' 
#' @export
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select `%>%`
setMethod("as.data.frame", "AdjacencyMatrix",
    function(x) {
              
        .nms <- assayNames(x)
                  
        l <- lapply(seq_along(.nms),
            function(.nms_i) assay(x, .nms_i))
        if (!directed(x)) {
            l <- lapply(l, function(l_i) {
                l_i[upper.tri(l_i)] <- NA
                l_i   
            })
        }
        l <- lapply(seq_along(l), function(i) l[[i]] %>%
            as.data.frame() %>%
            rownames_to_column(var = "Row") %>%
            pivot_longer(-"Row", names_to = "Col", values_to = .nms[i]))
        
        tbl <- dplyr::select(l[[1]], "Row", "Col")
        l <- lapply(l, function(x) x[, 3])
        df <- do.call("cbind", l)
        df <- cbind(tbl, df)
        
        ## remove rows with NA (directed(x) == FALSE)
        df <- df[!apply(df[, .nms, drop = FALSE], 1, 
            function(x) all(is.na(x))), ]
                
        return(df)
    }
)
