
#' @name threshold
#'
#' @aliases threshold
#'
#' @title  Threshold the statistical adjacency matrices
#'
#' @description
#' The function `threshold` takes as input an `AdjacencyMatrix` object 
#' containing adjacency matrices as returned from the function `statistical` or 
#' `spectral`. Depending on the `type` argument, `threshold` will identify the 
#' strongest link that are lower or higher a certain threshold 
#' (`type = "threshold"`) or identify the top `n` links (`type` either `"top1`, 
#' `"top2` or `"mean"`).
#' It will return this kind of information as a binary matrix in the form
#' of an `AdjacencyMatrix` object. 
#'
#' @param am `AdjacencyMatrix` object of `type` `"statistical"` or `"spectral"` 
#' as created from the function `statistical`. The object will contain the
#' adjacency matrices in the `assay` slot.
#'
#' @param type `character`, either `"threshold"`, `"top1`, `"top2` or
#' `"mean"`
#'
#' @param args `list`. Depending on the `type` arguments the list element
#' will be different. 
#' 
#' In the case of `type == "threshold"`, `args` has the 
#' entry `filter` (`character` of length 1). The character vector will specify the 
#' kind of filtering applied to the adjacency matrices. Elements in `filter`
#' will refer to the `assayNames`, e.g. 
#' `list(filter = "pearson_coef > 0.8")` will retain all edges with Pearson
#' correlation coefficients > 0.8. 
#' `list(filter = "pearson_coef > 0.8 & spearman_coef > 0.5")` will retain all 
#' edges with Pearson correlation coefficients > 0.8 AND Spearman correlation
#' coefficients > 0.5. 
#' `list(filter = "abs(pearson_coef) > 0.8 & spearman_coef > 0.5")` will retain all 
#' edges with Pearson correlation coefficients > 0.8 and < -0.8.
#' 
#' In the case of `type == "top1"`, `type == "top2"`, or `type == "mean"`, 
#' `args` has the entry `n` (`numeric` of length 1), that 
#' denotes the number of top ranks written to the
#' consensus matrix. Optionally, `args` has the entry `abs` which will take 
#' absolute values of the coefficients (will default to `FALSE` if `args$abs` 
#' is not specified).
#' 
#' @param values `character`, take from the adjacency matrix all values ("all"),
#' the minimum of the pairs ("min") or the maximum ("max")
#' a^*_{ij} = min(a_ij, a_ji)
#' a^*_{ij} = max(a_ij, a_ji)
#' 
#' @param na.rm `logical`, if set to `TRUE`, the `NA`s in the assay slots will 
#' not be taken into account when creating the `"consensus"` assay. If set 
#' to `FALSE`, the `NA`s will be taken into account and might be passed to the
#' `"consensus"` assay. If `FALSE` the user can set the filter e.g. to 
#' `(ggm_coef > 0 | is.na(ggm_coef))`, when there are `NA`s in 
#' `ggm_coef` to disregard `NA`s. 
#'
#' @details
#' `values` has to be set carefully depending on if the `AdjacencyMatrix` object
#' `am` is `directed` or not.
#' 
#' In the case of `type == "threshold"`, `args` has the 
#' entry `filter` (`character` of length 1). The character vector will specify the 
#' kind of filtering applied to the adjacency matrices. Elements in `filter`
#' will refer to the `assayNames`, e.g. 
#' `list(filter = "pearson_coef > 0.8")` will retain all edges with Pearson
#' correlation coefficients > 0.8. 
#' `list(filter = "pearson_coef > 0.8 & spearman_coef > 0.5")` will retain all 
#' edges with Pearson correlation coefficients > 0.8 AND Spearman correlation
#' coefficients > 0.5. 
#' `list(filter = "abs(pearson_coef) > 0.8 & spearman_coef > 0.5")` will retain all 
#' edges with Pearson correlation coefficients > 0.8 and < -0.8.
#'
#' If  `type` is equal to `"top1"`, `"top2"` or
#' `"mean"`, then `args` has to contain a numeric vector of length 1 that
#' gives the number of top ranks included in the returned adjacency matrix.
#' In this case
#' values that are 0 for the models `lasso`, `randomForest` and `bayes` are
#' set to `NaN`; values from correlation (Pearson and Spearman, including
#' for partial and semipartial correlation) and `clr` and `aracne` are
#' taken as they are.
#'
#' For `type = "top1"`, the best (i.e. lowest) rank in `am` is taken.
#' For `type = "top2"`, the second best (i.e. second lowest) rank in
#' `am` is taken.
#' For `type = "mean"`, the average rank in `am` is taken.
#' Subsequently the first `n` unique ranks are returned.
#'
#' @return `AdjacencyMatrix` object containing a binary adjacency matrix
#' given the links supported by the `type` and the `args` (in the slot 
#' `"consensus"`. The object will furthermore contain the supplied data input,
#' i.e. all assays from `am`. The slot `threshold` is set to `TRUE`.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' model <- c("pearson", "spearman")
#' args <- list()
#' am_stat <- statistical(x, model = model)
#'
#' ## type = "threshold"
#' args <- list(filter = "pearson_coef > 0.95 & spearman_coef > 0.95")
#' threshold(am = am_stat, type = "threshold", args = args)
#'
#' ## type = "top1"
#' args <- list(n = 10)
#' threshold(am = am_stat, type = "top1", args = args)
#'
#' ## type = "top2"
#' threshold(am = am_stat, type = "top2", args = args)
#'
#' ## type = "mean"
#' threshold(am = am_stat, type = "mean", args = args)
#'
#' @export
#' 
#' @importFrom rlang parse_expr
threshold <- function(am, 
    type = c("threshold", "top1", "top2", "mean"), 
    args, values = c("all", "min", "max"), na.rm = TRUE) {

    ## check match.arg for values
    type <- match.arg(type)
    values <- match.arg(values)
    if (!is.list(args)) stop("'args' is not a list")
    
    if (!is(am, "AdjacencyMatrix")) {
        stop("'am' is not an 'AdjacencyMatrix' object")
    }
    
    if (!validObject(am)) {
        stop("'am' must be a valid 'AdjacencyMatrix' object")
    }
    
    if (!type(am) %in% c("spectral", "statistical")) {
        stop("'am' must be an 'AdjacencyMatrix' either type 'spectral' or 'statistical'")
    }

    if (am@thresholded) {
        stop("'am' has been already thresholded")
    }
    
    ## args, either N for tops
    ## or a list of threshold
    if (any(duplicated(names(args)))) {
        stop("names(args) contain duplicated entries")
    }

    ## check args
    if (type == "threshold" & !is.character(args$filter)) {
        stop("'filter' must be character")
    }
    
    if (type != "threshold" && length(args$n) != 1 && !is.numeric(args$n)) {
        stop("args does not contain the numeric entry `n` of length 1")
    }

    ## create the cons matrix to store the consensus information
    a_1 <- assay(am, 1)
    cons <- matrix(0, nrow = ncol(a_1), ncol = ncol(a_1))
    rownames(cons) <- colnames(cons) <- colnames(a_1)
    diag(cons) <- NaN
    
    ## 
    df <- as.data.frame(am)
    n_nas <- apply(df, 2, function(x) sum(is.na(x)))
    n_nas <- n_nas[!names(n_nas) %in% c("Row", "Col")]
    sprintf("There are %s NAs in %s", as.numeric(n_nas), names(n_nas))
    
    if (type == "threshold") {
        
        if (na.rm) {
            ## get the Row and Col of the elements that match our filter 
            ## condition
            df_filter <- df |>
                dplyr::filter(!!rlang::parse_expr(args$filter))    
        } else {
            ## get the Row and Col of the elements where there is at least one
            ## NA
            df_na <- df[apply(df, 1, function(x) any(is.na(x))), ]
            
            ## get the Row and Col of the elements that match our filter 
            ## condition
            df_filter <- df |>
                dplyr::filter(!!rlang::parse_expr(args$filter))
            
            ## find the rows in df_na, that are not present in df_filter
            df_na <- df_na |>
                dplyr::filter(!("Row" %in% df_filter[, "Row"] &
                                                "Col" %in% df_filter[, "Col"]))
                
            ## set the elements in cons to NaN that are defined by df_na
            inds_row <- match(df_na[, "Row"], rownames(cons))
            inds_col <- match(df_na[, "Col"], colnames(cons))
            cons[cbind(inds_row, inds_col)] <- NaN
        }
        
        inds_row <- match(df_filter[, "Row"], rownames(cons))
        inds_col <- match(df_filter[, "Col"], colnames(cons))
        
        ## write the remaining elements to the consensus matrix
        cons[cbind(inds_row, inds_col)] <- 1
        if (!am@directed) 
            cons[cbind(inds_col, inds_row)] <- 1

    } else { ## if type is in "top1", "top2" or "mean"
        
        if (type(am) == "spectral") {
          ind_coef <- assayNames(am)
        }
      else {## if am is "statistical"
        ind_coef <- grep(pattern = "_coef", x = assayNames(am))
      }
        
        l <- as.list(assays(am)[ind_coef])

        l_df <- lapply(seq_along(l), function(x) {

            ## find corresponding model in l
            name_x <- names(l)[x]

            ## get corresponding adjacency matrix in l
            l_x <- l[[name_x]]
            
            if (is.null(args$abs)) args$abs <- FALSE
            if (args$abs) 
                l_x <- abs(l_x)

            ## take the respective minimum or maximum depending on `values`,
            ## do not do anything if `values` is equal to `all`
            if (values %in% c("min", "max")) {

                ## get values from the lower triangle
                lower_tri <- l_x[lower.tri(l_x)]

                ## get values from the upper triangle (requires transposing)
                l_x_t <- t(l_x)
                upper_tri <- l_x_t[lower.tri(l_x_t)]

                ## get min'max (argument values) of lower_tri and upper_tri
                values <- apply(rbind(lower_tri, upper_tri), 2, get(values))

                ## write back to the matrix
                l_x[lower.tri(l_x)] <- values
                l_x <- t(l_x)
                l_x[lower.tri(l_x)] <- values
            }

            ## for pearson/spearman correlation (incl. partial and
            ## semi-partial), lasso, randomForest, ggm, clr, aracne and bayes
            ## higher values correspond to higher confidence
            if (grepl(name_x, 
                pattern = "lasso_coef|randomForest_coef|bayes_coef")) {
                
                ## set values that are equal to 0 to NaN (values that are 0)
                ## do not explain the variability
                res <- getLinks(l_x, exclude = "== 0", decreasing = TRUE)
            }
            else {
                
                res <- getLinks(l_x, exclude = NULL, decreasing = TRUE)
            }
            res
        })

        ## bind together the ranks of the models, stored in l_df
        ranks <- lapply(l_df, function(l_i) l_i$rank)
        ranks <- do.call("cbind", ranks)
        colnames(ranks) <- names(l)

        ## calculate the consensus information, i.e. either get the first or
        ## second top rank per row or calculate the average across rows
        ## depending on the type argument
        cons_val <- topKnet(ranks, type, na.rm = na.rm)

        ## bind row and col information with cons information
        row_col <- l_df[[1]][, c("row", "col")]
        ranks <- cbind(row_col, cons_val)

        ## get the top N features
        top_n <- sort(unique(cons_val))[1:args$n]
        ranks_top <- ranks[cons_val %in% top_n, ]

        ## write links in ranks_top to binary adjacency matrix cons
        cons[as.numeric(rownames(ranks_top))] <- 1
        
        ## write NA to the ones where ranks are NA if na.rm = FALSE
        if (!na.rm) {
            ranks_NA <- ranks[is.na(cons_val), ]
            cons[as.numeric(rownames(ranks_NA))] <- NA    
        }
    }

    ## assign the consensus matrix to a new slot
    assay(am, "consensus") <- cons
    
    if (type %in% c("top1", "top2", "mean") & values %in% c("min", "max")) 
        am@directed <- FALSE
    am@thresholded <- TRUE
    
    return(am)
}


#' @name getLinks
#'
#' @aliases getLinks
#'
#' @title Write an adjacency matrix to a `data.frame`
#'
#' @description
#' `getLinks` vectorizes a numerical square `matrix` and writes the values
#' and their corresponding ranks to a `data.frame`.
#'
#' @param mat matrix containing the values of confidence for a link
#'
#' @param
#' exclude `character`, logical statement as `character` to set `TRUE`
#' values to NaN in `mat`, will be omitted if `exclude = NULL`
#'
#' @param decreasing `logical`, if `TRUE`, the highest confidence value will 
#' get the first rank, if `FALSE`, the lowest confidence value will get the 
#' first rank
#'
#' @details
#' `getLinks` is a helper function used in the function `threshold`.
#'
#' @return `data.frame` with entries `row`, `col`, `confidence` and `rank`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' mat <- matrix(0:8, ncol = 3, nrow = 3)
#' MetNet:::getLinks(mat, exclude = "== 0", decreasing = TRUE)
#'
#' @export
getLinks <- function(mat, exclude = "== 1", decreasing = TRUE) {

    if (ncol(mat) != nrow(mat)) {
        stop("`mat` is not a square matrix")
    }

    ## replace values that match exclude by NaN,
    ## will be omitted if exclude = NULL
    if (!is.null(exclude)) {
        exclude <- paste0("mat", exclude)
        mat[which(eval(parse(text = exclude)))] <- NaN
    }

    ## vectorize mat and write values of mat to confidence
    df <- data.frame(row = c(row(mat)), col = c(col(mat)), confidence = c(mat))

    ## if decreasing = TRUE, the highest confidence value will get the first 
    ## rank, recalculate the confidence values that the values with highest
    ## support have low values
    ## if decreasing = FALSE, the lowest confidence value will get the first 
    ## rank
    if (decreasing)
        conf <- max(df$confidence, na.rm = TRUE) - df$confidence
    else
        conf <- df$confidence

    ## calculate rank and add to data.frame
    df <- data.frame(df, rank = NaN)
    df$rank[!is.na(df$confidence)] <- rank(conf, na.last = NA)

    ## return
    return(df)
}


#' @name topKnet
#' @aliases topKnet
#' @title Return consensus ranks from a matrix containing ranks
#'
#' @description
#' `topKnet` returns consensus ranks depending on the `type` argument from
#' `ranks`, a matrix containing the ranks per statistical `model`.
#'
#' @param ranks `matrix` containing the ranks per statistical model (in columns)
#' and per feature pair (in rows)
#'
#' @param type `character`, either `"top1"`, `"top2"` or `"mean"`
#' 
#' @param na.rm `logical`, if set to `TRUE`, the `NA`s in the assay slots will 
#' not be taken into account when creating the `"top1"`, `"top2"` or `"mean"` 
#' of ranks. If set to `FALSE`, the `NA`s will be taken into account 
#' when creating the `"top1"`, `"top2"` or `"mean"` ranks. If `FALSE` the 
#' resulting aggregations will be `NA` if an `NA` is present in the coeffients
#' of one feature pair.
#'
#' @details
#' See Hase et al. (2014) for further details.
#'
#' @references
#' Hase et al. (2014):  Harnessing Diversity towards the Reconstructing of
#' Large Scale Gene Regulatory Networks. PLoS Computational Biology, 2013,
#' e1003361, doi:
#' [10.1371/journal.pcbi.1003361](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003361)
#'
#' @return `numeric` `vector`` with consensus ranks
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' ## na.rm == TRUE
#' ranks <- matrix(c(c(1, 2, 3), c(2, 1, 3)), ncol = 2)
#'
#' ## type = "top1"
#' MetNet:::topKnet(ranks = ranks, type = "top1", na.rm = TRUE)
#'
#' ## type = "top2"
#' MetNet:::topKnet(ranks = ranks, type = "top2", na.rm = TRUE)
#'
#' ## type = "mean"
#' MetNet:::topKnet(ranks = ranks, type = "mean", na.rm = TRUE)
#' 
#' ## na.rm == FALSE
#' ranks <- matrix(c(c(1, 2, 3), c(2, 1, 3)), ncol = 2)
#'
#' ## type = "top1"
#' MetNet:::topKnet(ranks = ranks, type = "top1", na.rm = FALSE)
#'
#' ## type = "top2"
#' MetNet:::topKnet(ranks = ranks, type = "top2", na.rm = FALSE)
#'
#' ## type = "mean"
#' MetNet:::topKnet(ranks = ranks, type = "mean", na.rm = FALSE)
#' 
#' ## na.rm == FALSE
#' ranks <- matrix(c(c(1, 2, NA), c(2, 1, 3)), ncol = 2)
#'
#' ## type = "top1"
#' MetNet:::topKnet(ranks = ranks, type = "top1", na.rm = FALSE)
#'
#' ## type = "top2"
#' MetNet:::topKnet(ranks = ranks, type = "top2", na.rm = FALSE)
#'
#' ## type = "mean"
#' MetNet:::topKnet(ranks = ranks, type = "mean", na.rm = FALSE)
topKnet <- function(ranks, type, na.rm = TRUE) {

    if ((!is.matrix(ranks)) || !is.numeric(ranks)) {
        stop("ranks is not a numerical matrix")
    }

    if (! type %in% c("top1", "top2", "mean")) {
        stop("type neither 'top1', 'top2' or 'mean'")
    }

    ## calculate the consensus information, i.e. either get the first or
    ## second top rank per row or calculate the average across rows
    ## depending on the type argument
    if (type == "top1") {
        ## get the lowest rank
        cons_val <- apply(ranks, 1, function(x) {
            if (!all(is.na(x))) {
                min(x, na.rm = na.rm)
            } else {
                NaN
            }
        })
    }

    if (type == "top2") {

        ## "top2" will work only if there are at least two `model`s
        ## return error if ncol(ranks) == 1
        if (ncol(ranks) == 1) {
            stop("ncol(ranks) has to be > 1")
        }

        ## get the second lowest rank (only if there are at least two or more
        ## non-NA values, otherwise return NA --> sort(x)[2] will return
        ## NA if there are less elements)
        cons_val <- apply(ranks, 1, function(x) {
           
            if (!all(is.na(x))) {
                if (na.rm)
                    sort(x)[2]
                else
                    sort(x, na.last = TRUE)[2]
            } else {
                NaN
            }    
            
        })
    }

    if (type == "mean") {
        ## get the average of all ranks
        cons_val <- apply(ranks, 1, mean, na.rm = na.rm)
    }

    return(cons_val)
}
