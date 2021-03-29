##' @title Placeholder for generics functions documentation
##'
##' @name AllGenerics
##' @rdname AllGenerics
NULL


###
### Accessors
###

#' @name length.AdjacencyMatrix
#' 
#' @aliases length
#' 
#' @rdname AdjacencyMatrix-class
#' 
#' @title Return length of `AdjacencyMatrix` object
#' 
#' @description 
#' `length` will return the length of an `AdjacencyMatrix` object
#' (number of rows of an assay).
#' 
#' @param x an instance of class `AdjacencyMatrix`
#' 
#' @return `numeric` of length ``
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @exportMethod length
setMethod("length", "AdjacencyMatrix",
    function(x) nrow(assay(x, 1))
)

#' @name dim.AdjacencyMatrix
#' 
#' @aliases dim
#' 
#' @rdname AdjacencyMatrix-class
#' 
#' @title Return dimension of `AdjacencyMatrix` object
#' 
#' @description 
#' `dim` will return the length of an `AdjacencyMatrix` object
#' (number of rows of an assay, number of cols of an assay).
#' 
#' @param x an instance of class `AdjacencyMatrix`
#' 
#' @return `numeric` of length 2
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @exportMethod dim
setMethod("dim", "AdjacencyMatrix",
    function(x) c(nrow(assay(x, 1)), ncol(assay(x, 1)))
)

#' @name type.AdjacencyMatrix
#' 
#' @aliases type
#' 
#' @rdname AdjacencyMatrix-class
#' 
#' @title Return type of `AdjacencyMatrix` object
#' 
#' @description 
#' `type` will return the type of an `AdjacencyMatrix` (`statistical`, 
#' `structural` or `combine`).
#' 
#' @param x instance of class `AdjacencyMatrix`
#' 
#' @return `character` of length 1
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @exportMethod type
setMethod("type", "AdjacencyMatrix",
    function(x) x@type)

setGeneric("directed",
    function(object) standardGeneric("directed"))

#' @name directed.AdjacencyMatrix
#' 
#' @aliases directed
#' 
#' @title Return the information on directed of an `AdjacencyMatrix`
#' 
#' @rdname AdjacencyMatrix-class
#' 
#' @description 
#' `directed` will return the information on directed of an `AdjacencyMatrix`
#' 
#' @param object instance of class `AdjacencyMatrix`
#' 
#' @return `character` of length 1
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @exportMethod directed
setMethod("directed", "AdjacencyMatrix",
    function(object) object@directed)

setGeneric("thresholded",
    function(object) standardGeneric("thresholded"))

#' @name thresholded.AdjacencyMatrix
#' 
#' @aliases thresholded
#' 
#' @title Return the information on thresholded of an `AdjacencyMatrix`
#' 
#' @rdname AdjacencyMatrix-class
#' 
#' @description 
#' `thresholded` will return the information on directed of an `AdjacencyMatrix`
#' 
#' @param object instance of class `AdjacencyMatrix`
#' 
#' @return `logical` of length 1
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @exportMethod thresholded
setMethod("thresholded", "AdjacencyMatrix",
    function(object) object@thresholded)



###
### Show
###

#' @name show.AdjacencyMatrix
#' 
#' @aliases show
#' 
#' @rdname AdjacencyMatrix-class
#' 
#' @title show method for object of class `AdjacencyMatrix`
#' 
#' @description 
#' `show` will print summary information on an object of class `AdjacencyMatrix`
#' 
#' @param object instance of class `AdjacencyMatrix`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @exportMethod show
setMethod("show", "AdjacencyMatrix",
    function(object) {
          
        cat("class:", class(object), "\n")
        cat("dim:", dim(object), "\n")
              
        ## type()
        .type <- type(object)
        if (is.null(.type))
            .type <- character(length(type(object)))
        coolcat("type(%d): %s\n", .type)
        
        ## directed()
        .directed <- directed(object)
        if (is.null(.directed))
            .directed <- character(length(directed(object)))
        coolcat("directed(%d): %s\n", .directed)
              
        ## threshold()
        .threshold <- thresholded(object)
        if (is.null(.threshold))
            .threshold <- character(length(thresholded(object)))
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
#' @name as.data.frame.AdjacencyMatrix
#' 
#' @aliases as.data.frame
#' 
#' @rdname AdjacencyMatrix-class
#' 
#' @title Return `data.frame` of adjacency matrices
#' 
#' @description 
#' `as.data.frame` returns the adjacency matrices (stored in the `assays` slot)
#' and returns information on the nodes and the associated information on 
#' edges as a data frame.
#' 
#' @param x instance of class `AdjacencyMatrix` 
#' 
#' @return `data.frame`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select
setMethod("as.data.frame", "AdjacencyMatrix",
    function(x) {
              
        .nms <- assayNames(x)
                  
        l <- lapply(seq_along(.nms),
            function(i) assay(x, i))
        if (!directed(x)) {
            l <- lapply(l, function(i) {
                i[upper.tri(i)] <- NA
                i   
            })
        }
        l <- lapply(seq_along(l), function(i) l[[i]] %>%
            as.data.frame() %>%
            rownames_to_column(var = "Row") %>%
            pivot_longer(-Row, names_to = "Col", values_to = .nms[i]) %>%
            ## remove rows with NA (directed(x) == FALSE)
            filter(., !is.na(get(.nms[i]))))
        
        tbl <- dplyr::select(l[[1]], "Row", "Col")
        l <- lapply(l, function(x) x[, 3])
        df <- do.call("cbind", l)
        df <- cbind(tbl, df)
                
        return(df)
    }
)
