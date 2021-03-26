
###
### Accessors
###

setMethod("length", "AdjacencyMatrix",
    function(object) nrow(assay(object, 1))
)

setMethod("dim", "AdjacencyMatrix",
    function(object) c(length(object), ncol(assay(object, 1)))
)

setMethod("type", "AdjacencyMatrix",
    function(object) object@type)

setGeneric("directed",
    function(object) standardGeneric("directed"))

setMethod("directed", "AdjacencyMatrix",
    function(object) object@directed)

setGeneric("thresholded",
    function(object) standardGeneric("thresholded"))

setMethod("thresholded", "AdjacencyMatrix",
    function(object) object@thresholded)



###
### Show
###

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
