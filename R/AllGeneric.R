
###
### Accessors
###

setMethod("length", "AdjacencyMatrix",
    function(x) nrow(x@elementMetadata)
)

setMethod("dim", "AdjacencyMatrix",
    function(x) c(length(x), length(x))
)

setMethod("type", "AdjacencyMatrix",
    function(x) x@type)

setGeneric("thresholded",
    function(object) standardGeneric("thresholded"))

setMethod("thresholded", "AdjacencyMatrix",
    function(x) x@thresholded)



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
              
        ## threshold()
        .threshold <- thresholded(object)
        if (is.null(.threshold))
            .threshold <- character(length(thresholded(object)))
        coolcat("type(%d): %s\n", .threshold)
              
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

#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select
setMethod("as.data.frame", "AdjacencyMatrix",
    function(x) {
              
    .nms <- assayNames(x)
              
    l <- lapply(seq_along(.nms),
        function(i) assay(x, i) %>%
                as.data.frame() %>%
                rownames_to_column(var = "Row") %>%
                pivot_longer(-Row, names_to = "Col", values_to = .nms[i]))
            tbl <- select(l[[1]], "Row", "Col")
            l <- lapply(l, function(x) x[, 3])
            df <- do.call("cbind", l)
            tbl <- cbind(tbl, df) %>% as.data.frame()
            
            return(tbl)
    }
)
