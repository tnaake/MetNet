#' @importFrom stabs stabsel.matrix glmnet.lasso
NULL
#' @importFrom GENIE3 GENIE3
NULL
#' @importFrom mpmi cmi
NULL
#' @importFrom parmigene clr aracne.a
NULL
#' @importFrom bnlearn fast.iamb arcs
NULL
#' @importFrom sna consensus
NULL
#' @importFrom BiocParallel bplapply
NULL
#' @importFrom stats formula p.adjust sd cor
NULL
#' @importFrom methods formalArgs
NULL
#' @importFrom ppcor pcor spcor
NULL

#' @name lasso
#' @aliases lasso
#' @title Create an adjacency matrix based on LASSO
#'
#' @description
#' `lasso` infers a adjacency matrix using
#' LASSO using the `stabsel.matrix` function from the
#' `stabs` package. `lasso` extracts the  predictors from the
#' function `stabsel.matrix` and writes the coefficients
#' to an adjacency matrix.
#'
#' @param
#' x matrix, where columns are the samples and the rows are features
#' (metabolites), cell entries are intensity values
#'
#' @param
#' parallel logical, should computation be parallelized? If
#' `parallel = TRUE` the `bplapply` will be applied if
#' `parallel = FALSE` the `lapply` function will be applied.
#'
#' @param ... parameters passed to `stabsel.matrix`
#'
#' @details For use of the parameters used in the `stabsel.matrix` function,
#' refer to `?stabs::stabsel.matrix`.
#'
#' @return 
#' matrix, matrix with edges inferred from LASSO algorithm
#' `stabsel.matrix`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' x_z <- t(apply(x, 1, function(y) (y - mean(y)) / sd(y)))
#' \dontrun{lasso(x_z, PFER = 0.95, cutoff = 0.95)}
#'
#' @export
lasso <- function(x, parallel = FALSE, ...) {

    ## x should be z-scaled
    if (parallel) {
        l1 <- bplapply(seq_len(nrow(x)), function(i) {
            x_l1 <- t(x[-i, ]); y_l1 <- x[i, ]

            ## lasso: alpha set to 1
            ## allow for compatibility of arguments 
            l1 <- threeDotsCall("stabsel.matrix", x = as.matrix(x_l1),
                    y = y_l1, fitfun = glmnet.lasso,
                    args.fitfun = list("alpha" = 1), ...)

            ## return selection probabilities of features that are not 0
            return(l1$max[l1$max != 0])
        }) 
        
    } else {
        l1 <- lapply(seq_len(nrow(x)), function(i) {
            x_l1 <- t(x[-i, ]); y_l1 <- x[i, ]
            
            ## lasso: alpha set to 1
            ## allow for compatibility of arguments
            l1 <- threeDotsCall("stabsel.matrix", x = as.matrix(x_l1),
                    y = y_l1, fitfun = glmnet.lasso,
                    args.fitfun = list("alpha" = 1), ...)
            
            ## return selection probabilities of features that are not 0
            return(l1$max[l1$max != 0])
        })
    }

    ## create a matrix to store the values
    l1_mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    colnames(l1_mat) <- rownames(l1_mat) <- rownames(x)

    ## write the selection probabilitiy values in l1 to l1_mat
    for (i in seq_len(length(l1))) {
        l1_mat[names(l1[[i]]), i] <- l1[[i]]
    }

    return(l1_mat)
}

#' @name randomForest
#' 
#' @aliases randomForest
#' 
#' @title Create an adjacency matrix based on random forest
#' 
#' @description
#' `randomForest` infers an adjacency matrix using
#' random forest using the `GENIE3` function from the 
#' `GENIE3` package. `randomForest` returns the importance of the link
#' between features in the form of an adjacency matrix.
#' 
#' @param
#' x matrix, where columns are the samples and the rows are features
#' (metabolites), cell entries are intensity values
#' 
#' @param ... parameters passed to `GENIE3`
#' 
#' @details For use of the parameters used in the `GENIE3` function,
#' refer to `?GENIE3::GENIE3`. The arguments `regulators` and `targets` are
#' set to `NULL`. Element \eqn{w_{i,j}} (row i, column j) gives the importance
#' of the link from i to j.
#' 
#' @return matrix, matrix with the importance of the links inferred from
#' random forest algorithm implemented by `GENIE3`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' randomForest(x)
#' 
#' @export
randomForest <- function(x, ...) {
    
    ## GENIE3 returns the importance of the link from "regulator gene" i to
    ## target gene "j" in the form of a weighted adjacency matrix
    ## set regulators and targets to NULL that they cannot be changed
    rf <- threeDotsCall(GENIE3::GENIE3, exprMatrix = x, regulators = NULL,
        targets = NULL, ...)
    
    return(rf)
}

#' @name clr
#' 
#' @aliases clr
#' 
#' @title Create an adjacency matrix based on context likelihood or
#' relatedness network
#' 
#' @description  `clr` infers an adjacency matrix using
#' context likelihood/relatedness network using the `clr` function from
#' the `parmigene` package. `clr` will
#' return the adjacency matrix containing the Context Likelihood of
#' Relatedness Network-adjusted scores of Mutual
#' Information values.
#' 
#' @param mi matrix, where columns and the rows are features
#' (metabolites), cell entries are mutual information values between the
#' features. As input, the mutual information (e.g. raw MI estimates or
#' Jackknife bias corrected MI estimates) from the `cmi` function of the
#' `mpmi` package can be used.
#'
#' @details For more details on the `clr` function,
#' refer to `?parmigene::clr`. CLR computes the score
#' \eqn{sqrt(z_i ^2 + z_j ^2)} for each pair of variables i, j, where
#' \eqn{z_i = max(0, ( I(X_i, X_j) - mean(X_i) ) / sd(X_i) )}.
#' \eqn{mean(X_i)} and \eqn{sd(X_i)} are the mean and standard deviation of the
#' mutual information values \eqn{I(X_i, X_k)} for all \eqn{k = 1, ..., n}.
#' For more information on the CLR algorithm
#' see Faith et al. (2007).
#'
#' @references
#' Faith et al. (2007): Large-Scale Mapping and Validation of Escherichia coli
#' Transcriptional Regulation from a Compendium of Expression Profiles.
#' PLoS Biology, e8, doi:
#' [10.1371/journal.pbio.0050008]( https://doi.org/10.1371/journal.pbio.0050008)
#'
#' @return matrix, matrix with edges inferred from Context Likelihood of
#' Relatedness Network algorithm `clr`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' x_z <- t(apply(x, 1, function(y) (y - mean(y)) / sd(y)))
#' mi_x_z <- mpmi::cmi(x_z)$bcmi
#' clr(mi_x_z)
#'
#' @export
clr <- function(mi) {

    ## call the clr function from the parmigene package
    clr_mat <- parmigene::clr(mi)

    ## assign col- and rownames to clr_mat
    colnames(clr_mat) <- rownames(clr_mat) <- rownames(mi)
    return(clr_mat)
}

#' @name aracne
#'
#' @aliases aracne
#'
#' @title Create an adjacency matrix based on algorithm for the reconstruction
#' of accurate cellular networks
#'
#' @description  `aracne` infers an adjacency matrix using
#' the algorithm for the reconstruction of accurate cellular networks
#' using the `aracne.a` function from the
#' `parmigene` package. The function `aracne` will return the weighted
#' adjacency matrix of the inferred network after applying `aracne.a`.
#'
#' @param mi matrix, where columns and the rows are features
#' (metabolites), cell entries are mutual information values between the
#' features. As input, the mutual information (e.g. raw MI estimates or
#' Jackknife bias corrected MI estimates) from the `cmi` function of the
#' `mpmi` package can be used.
#'
#' @param eps numeric, used to remove the weakest edge of each triple of nodes
#'
#' @details For more details on the `aracne.a` function,
#' refer to `?parmigene::aracne.a`. `aracne.a` considers each triple of
#' edges independently and removes the weakest one if
#' \eqn{MI(i, j) < MI(j, k) - eps} and \eqn{MI(i, j) < MI(i, k) - eps}. See
#' Margolin et al. (2006) for further information.
#'
#' @references
#' Margolin et al. (2006): ARACNE : An algorithm for the reconstruction of
#' gene regulatory networks in a mammalian cellular context. BMC Bioinformatics,
#' S7, doi:
#' [10.1186/1471-2105-7-S1-S7](https://doi.org/10.1186/1471-2105-7-S1-S7)
#'
#' @return
#' matrix, matrix with edges inferred from Reconstruction of accurate
#' cellular networks algorithm `aracne.a`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' x_z <- t(apply(x, 1, function(y) (y - mean(y)) / sd(y)))
#' mi_x_z <- mpmi::cmi(x_z)$bcmi
#' aracne(mi_x_z, eps = 0.05)
#'
#' @export
aracne <- function(mi, eps = 0.05) {

    ## call the aracne.a function from the parmigene package
    aracne_mat <- parmigene::aracne.a(mi, eps = eps)

    ## assign col- and rownames to aracne_mat
    colnames(aracne_mat) <- rownames(aracne_mat) <- rownames(mi)
    return(aracne_mat)
}

#' @name correlation
#'
#' @aliases correlation
#'
#' @title Create an adjacency matrix based on correlation
#'
#' @description
#' `correlation` infers an adjacency matrix using
#' correlation using the `cor` function (from the
#' `stats` package), `pcor` (from `ppcor`) or
#' `spcor` (from `ppcor`). `correlation` extracts the reported pair-wise
#' correlation coefficients from the function `corAndPvalue`, `pcor` or `spcor`
#' and will return
#' the weighted adjacency matrix of the absolute correlation values.
#'
#' @param
#' x matrix, where columns are the samples and the rows are features
#' (metabolites), cell entries are intensity values
#'
#' @param
#' type `character`, either "pearson", "spearman", "pearson_partial",
#' "spearman_partial", "pearson_semipartial" or "spearman_semipartial".
#'
#' @param
#' use `character` string giving a method for computing covariance in the
#' presence of missing values, Only for `type = "pearson"` or
#' `type = "spearman"`. For further information see `?stats::cor`
#'
#' @details
#' If `"pearson"` or `"spearman"` is used as a `method`, the function
#' `corAndPvalue` from `stats` will be employed.
#'
#' If `"pearson_partial"` or `"spearman_partial"?` is used as a `method` the
#' function `pcor` from `spcor` will be employed.
#'
#' If `"pearson_semipartial"` or `"spearman_semipartial"` is used as a
#' `method` the function `spcor` from `spcor` will be employed.
#'
#' `type` will be passed to argument `method` in `cor`
#' (in the case of `"pearson"` or `"spearman"`) or to `method` in `pcor`
#' (`"pearson"` and `"spearman"` for `"pearson_partial"` and
#' `"spearman_partial"`, respectively) or to `method` in `spcor`
#' (`"pearson"` or `"spearman"` for `"pearson_semipartial"` and
#' `"spearman_semipartial"`, respectively).
#'
#' @return
#' matrix, matrix with edges inferred from correlation algorithm
#' `corAndPvalue`, `pcor` or `spcor` (depending on the chosen `method`)
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' correlation(x, type = "pearson")
#'
#' @export
correlation <- function(x, type = "pearson", use = "pairwise.complete.obs") {

    ## pearson/spearman
    if (type %in% c("pearson", "spearman")) {
        cor_mat <- cor(x = t(x), method = type, use = use)
    }

    ## partial pearson/spearman
    if (type %in% c("pearson_partial", "spearman_partial")) {
        if (type == "pearson_partial") method <- "pearson"
        if (type == "spearman_partial") method <- "spearman"
        cor_mat <- ppcor::pcor(t(x), method = method)$estimate
    }

    ## semipartial pearson/spearman
    if (type %in% c("pearson_semipartial", "spearman_semipartial")) {
        if (type == "pearson_semipartial") method <- "pearson"
        if (type == "spearman_semipartial") method <- "spearman"
        cor_mat <- ppcor::spcor(t(x), method = method)$estimate
    }

    ## assign col- and rownames to cor_mat
    colnames(cor_mat) <- rownames(cor_mat) <- rownames(x)

    ## get absolute values
    cor_mat <- abs(cor_mat)

    return(cor_mat)
}

#' @name bayes
#'
#' @aliases bayes
#'
#' @title Create an adjacency matrix based on score-based structure learning
#' algorithm
#'
#' @description  
#' `bayes` infers an adjacency matrix using score-based structure learning
#' algorithm `boot.strength` from the
#' `bnlearn` package. `bayes` extracts then the reported
#' connections from running the `boot.strength` function and assigns the
#' strengths of the arcs of the Bayesian connections to an adjacency matrix.
#' `bayes` returns this weighted adjacency matrix.
#'
#' @param
#' x `matrix` where columns are the samples and the rows are features
#' (metabolites), cell entries are intensity values
#'
#' @param algorithm `character`,
#' structure learning to be applied to the
#' bootstrap replicates (default is `"tabu"`)
#'
#' @param R `numeric`, number of bootstrap replicates
#'
#' @param ... parameters passed to `boot.strength`
#'
#' @details
#' `boot.strength` measures the strength of the
#' probabilistic relationships by the arcs of a Bayesian network, as learned
#' from bootstrapped data. By default `bayes` uses the
#' Tabu greedy search.
#'
#' For use of the parameters used in the `boot.strength` function,
#' refer to `?bnlearn::boot.strength`. For further information see also
#' Friedman et al. (1999) and Scutari and Nagarajan (2001).
#'
#' @references
#' Friedman et al. (1999): Data Analysis with Bayesian Networks: A Bootstrap
#' Approach. Proceedings of the 15th Annual Conference on Uncertainty in
#' Artificial Intelligence, 196-201.
#'
#' Scutari and Nagarajan (2011): On Identifying Significant Edges in Graphical
#' Models. Proceedings of the Workshop Probabilistic Problem Solving in
#' Biomedicine of the 13th Artificial Intelligence in Medicine Conference,
#' 15-27.
#'
#' @return
#' `matrix` with edges inferred from score-based structure
#' learning algorithm `boot.strength`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' bayes(x, algorithm = "tabu", R = 100)
#'
#' @export
bayes <- function(x, algorithm = "tabu", R = 100, ...) {

    x_df <- data.frame(t(x))

    ## allow for compatibility of arguments
    strength <- threeDotsCall(bnlearn::boot.strength, data = x_df,
        algorithm = algorithm, R = R, ...)

    ## create empty bs_mat to be filled with connections
    bs_mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    colnames(bs_mat) <- rownames(bs_mat) <- rownames(x)

    ## write to bs_mat
    for(i in seq_len(nrow(strength))) {
        tmp <- as.character(strength[i, ])
        names(tmp) <- names(strength[i, ])
        bs_mat[tmp["from"], tmp["to"] ] <- tmp["strength"]
    }

    mode(bs_mat) <- "numeric"

    return(bs_mat)
}

#' @name addToList
#'
#' @aliases addToList
#'
#' @title Add adjacency matrix to list
#'
#' @description 
#' This helper function used in the function
#' `statistical` adds an adjacency matrix to a `list` of
#' adjacency matrices.
#'
#' @param l `list` of adjacency matrices
#'
#' @param name `character`, name of added entry
#'
#' @param object `matrix` that will be added
#'
#' @details
#' The function `addToList` is a helper function used internally in
#' `statistical`.
#'
#' @return
#' `list` containing the existing matrices and the added matrix
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' cor_pearson <- correlation(x, type = "pearson")
#' cor_spearman <- correlation(x, type = "spearman")
#' l <- list(pearson = cor_pearson)
#' MetNet:::addToList(l, "spearman", cor_spearman)
addToList <- function(l, name, object) {

    ## test validity of objects
    if (!is.list(l)) {
        stop("l is not a list")
    }

    if (!is.character(name)) {
        stop("name is not a character")
    }

    if (!is.matrix(object)) {
        stop("object is not a matrix")
    }

    ## add object to l
    new_index <- length(l) + 1
    l[[new_index]] <- object

    ## assign the name to the newly added entry
    names(l)[new_index] <- name

    return(l)
}

#' @name statistical
#'
#' @aliases statistical
#'
#' @title Create a list of adjacency matrices from statistical methods
#'
#' @description
#' The function `statitical` infers adjacency matrix topologies from
#' statistical methods and returns matrices of these networks in a `list`. The
#' function includes functionality to calculate adjacency matrices based on
#' LASSO (L1 norm)-regression, random forests, context likelihood of
#' relatedness (CLR), the algorithm for the reconstruction of accurate
#' cellular networks (ARACNE), Pearson correlation (also partial and
#' semipartial), Spearman correlation (also partial and semipartial)
#' and score-based structure learning (Bayes). The function returns a
#' list of adjacency matrices that are defined by `model`.
#'
#' @param
#' x `matrix` that contains intensity values of
#' features/metabolites (rows) per sample (columns).
#'
#' @param
#' model `character` vector containing the methods that will be used
#' (`"lasso"`, `"randomForest"`, `"clr"`, `"aracne"`, `"pearson"`,
#' `"pearson_partial"`, `"pearson_semipartial"`, `"spearman"`,
#' `"spearman_partial"`, `"spearman_semipartial"`, `"bayes"`)
#'
#' @param
#' ... parameters passed to the functions  `lasso`, `randomForest`,
#' `clr`, `aracne`, `correlation` and/or `bayes`
#'
#' @details
#' The function `statistical` includes functionality to calculate adjacency
#' matrices based on
#' LASSO (L1 norm)-regression, random forests, context likelihood of
#' relatedness (CLR), the algorithm for the reconstruction of accurate
#' cellular networks (ARACNE), Pearson correlation (also partial and
#' semipartial), Spearman correlation (also partial and semipartial)
#' and Constraint-based structure learning (Bayes).
#'
#' `statistical` calls the function
#' `lasso`, `randomForest`, `clr`, `aracne`,
#' `correlation` (for `"pearson"`, `"pearson_partial"`, `"pearson_semipartial"`,
#' `"spearman"`, `"spearman_partial"`, `"spearman_semipartial"`) and/or `bayes`
#' as specified by `model`. It will create adjacency matrices using the
#' specified methods and will return a `list` containing the weighted
#' adjacency matrices.
#'
#' Internally `x` will be z-scaled and the z-scaled object
#' will be used in `lasso`, `clr` and/or `aracne`.
#'
#' @return `list` containing the respective adjacency matrices specified by
#' `model`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' statistical(x = x, model = c("pearson", "spearman"))
#'
#' @export
statistical <- function(x, model, ...) {

    ## check if model complies with the implemented model and return error
    ## if not so
    if (!(all(model %in% c("lasso", "randomForest", "clr", "aracne",
            "pearson", "pearson_partial", "pearson_semipartial",
            "spearman", "spearman_partial", "spearman_semipartial", "bayes"))))
        stop("'model' not implemented in statistical")
  
    ## check if x is numeric matrix and return error if not so
    if (mode(x) != "numeric") stop("x is not a numerical matrix")

    ## z-scale x and transpose
    x_z <- apply(x, 1, function(y) {
        (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
    })
    x_z <- t(x_z)

    l <- list()

    ## add entry for lasso if "lasso" is in model
    if ("lasso" %in% model) {
        lasso <- lasso(x = x_z, ...)
        diag(lasso) <- NaN
        l <- addToList(l, "lasso", lasso)
        print("lasso finished")
    }

    ## add entry for randomForest if "randomForest" is in model
    if ("randomForest" %in% model) {
        randomForest <- randomForest(x = x, ...)
        diag(randomForest) <- NaN
        l <- addToList(l, "randomForest", randomForest)
        print("randomForest finished.")
    }

    ## calculate mutual information if "clr" or "aracne" is in model
    if (any(c("clr", "aracne") %in% model)) {
        mi_x_z <- mpmi::cmi(t(x_z))$bcmi
        rownames(mi_x_z) <- colnames(mi_x_z) <- rownames(x)
    }

    ## add entry for clr if "clr" is in model
    if ("clr" %in% model) {
        clr <- threeDotsCall("clr", mi = mi_x_z, ...)
        diag(clr) <- NaN
        l <- addToList(l, "clr", clr)
        print("clr finished.")
    }

    ## add entry for aracne if "aracne" is in model
    if ("aracne" %in% model) {
        aracne <- threeDotsCall("aracne", mi = mi_x_z, ...)
        diag(aracne) <- NaN
        l <- addToList(l, "aracne", aracne)
        print("aracne finished.")
    }

    ## add entry for pearson if "pearson" is in model
    if ("pearson" %in% model) {
        pearson <- threeDotsCall("correlation", x = x, type = "pearson", ...)
        diag(pearson) <- NaN
        l <- addToList(l, "pearson", pearson)
        print("pearson finished.")
    }

    ## add entry for pearson_partial if "pearson_partial" is in model
    if ("pearson_partial" %in% model) {
        pearson_partial <- threeDotsCall("correlation", x = x, 
            type = "pearson_partial", ...) 
        diag(pearson_partial) <- NaN
        l <- addToList(l, "pearson_partial", pearson_partial)
        print("pearson_partial finished.")
    }

    ## add entry for pearson_semipartial if "pearson_semipartial" is in model
    if ("pearson_semipartial" %in% model) {
        pearson_sp <- threeDotsCall("correlation", x = x, 
            type = "pearson_semipartial", ...)
        diag(pearson_sp) <- NaN
        l <- addToList(l, "pearson_semipartial", pearson_sp)
        print("pearson_semipartial finished.")
    }

    ## add entry for spearman if "spearman" is in model
    if ("spearman" %in% model) {
        spearman <- threeDotsCall("correlation", x = x, type = "spearman", ...)
        diag(spearman) <- NaN
        l <- addToList(l, "spearman", spearman)
        print("spearman finished.")
    }

    ## add entry for spearman_partial if "spearman_partial" is in model
    if ("spearman_partial" %in% model) {
        spearman_partial <- threeDotsCall("correlation", x = x,
            type = "spearman_partial", ...)
        diag(spearman_partial) <- NaN
        l <- addToList(l, "spearman_partial", spearman_partial)
        print("spearman_partial finished.")
    }

    ## add entry for spearman_semipartial if "spearman_semipartial" is in model
    if ("spearman_semipartial" %in% model) {
        spearman_sp <- threeDotsCall("correlation", x = x,
            type = "spearman_semipartial", ...)
        diag(spearman_sp) <- NaN
        l <- addToList(l, "spearman_semipartial", spearman_sp)
        print("spearman_semipartial finished.")
    }

    ## add entry for bayes if "bayes" is in model
    if ("bayes" %in% model) {
        bayes <- threeDotsCall("bayes", x = x, ...)
        diag(bayes) <- NaN
        l <- addToList(l, "bayes", bayes)
        print("bayes finished.")
    }

    return(l)
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
#' @details
#' `getLinks` is a helper function used in the function `threshold`.
#'
#' @return `data.frame` with entries `row`, `col`, `confidence` and `rank`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' mat <- matrix(0:8, ncol = 3, nrow = 3)
#' MetNet:::getLinks(mat, exclude = "== 0")
#'
#' @export
getLinks <- function(mat, exclude = "== 1") {

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

    ## the highest confidence value should get the first rank
    ## recalculate the confidence values that the values with highest
    ## support have low values
    conf <- max(df$confidence, na.rm = TRUE) - df$confidence

    ## calculate rank and add to data.frame
    df <- data.frame(df, rank = NaN) 
    df$rank[!is.na(df$confidence)] <- rank(conf, na.last = NA)

    ## return
    return(df)
}

#' @name threshold
#'
#' @aliases threshold
#'
#' @title  Threshold the statistical adjacency matrices
#'
#' @description 
#' The function `threshold` takes as input a list of adjacency matrices
#' as returned from the function `statistical`. Depending on the `type`
#' argument, 'threshold` will identify the strongest link that are 
#' lower or higher a certain threshold (`type = "threshold"`) or
#' identify the top `n` links (`type` either `"top1`, `"top2` or `"mean"`).
#'
#' @param statistical `list` containing adjacency matrices
#'
#' @param type `character`, either `"threshold"`, `"top1`, `"top2` or
#' `"mean"`
#'
#' @param args `list` of arguments, has to contain thresholds for weighted
#' adjacency matrices depending on the statistical model
#' (a named list, where names are identical to `model`s in `statistical`)
#' or a numerical
#' vector of length 1 that denotes the number of top ranks written to the
#' consensus matrix (a named list with entry `n`)
#'
#' @param ... parameters passed to the function `consensus` in the
#' `sna` package (only for `type = "threshold"`)
#'
#' @details
#' The entries of `args` differ depending on the argument `type`.
#' If `type = "treshhold"`, then `args` has to contain numeric vector of
#' length 1 with names equal to `names(statistical)` for each `model`
#' (`names(statistical)`) and the entry `threshold`, a numerical `vector(1)`
#' to threshold the consensus
#' matrix after using the `consensus` function from the `sna` package.
#' Depending on the chosen `method` in `consensus`, the `threshold` value of the
#' consensus adjacency matrix should be chosen accordingly to report a
#' connection by different statistical modelss.
#'
#' When combining the adjacency matrices the
#' `threshold` value defines if an edge is reported or not. For
#' `method = "central.graph"` threshold should be set to 1 by default. For other
#' values of `method`, the value should be carefully defined by the user.
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
#' For `type = "top1"`, the best (i.e. lowest) rank in `statistical` is taken.
#' For `type = "top2"`, the second best (i.e. second lowest) rank in
#' `statistical` is taken.
#' For `type = "mean"`, the average rank in `statistical` is taken.
#' Subsequently the first `n` unique ranks are returned.
#'
#' @return `matrix`, binary adjacency matrix given the links supported by the
#' `type` and the `args`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' model <- c("pearson", "spearman")
#' args <- list("pearson" = 0.95, "spearman" = 0.95, n = 10)
#' l <- statistical(x, model = model)
#'
#' ## type = "threshold" 
#' args <- list("pearson" = 0.95, "spearman" = 0.95, threshold = 1)
#' threshold(statistical = l, type = "threshold", args = args)
#'
#' ## type = "top1" 
#' args <- list(n = 10)
#' threshold(statistical = l, type = "top1", args = args)
#'
#' ## type = "top2"
#' threshold(statistical = l, type = "top2", args = args)
#'
#' ## type = "mean"
#' threshold(statistical = l, type = "mean", args = args)
#' @export
threshold <- function(statistical, type, args, ...) {

    l <- statistical
    ## args, either N for tops
    ## or a list of threshold
    if (any(duplicated(names(args)))) {
        stop("names(args) contain duplicated entries")
    }

    if (!type %in% c("top1", "top2", "mean", "threshold"))
        stop("type not in 'top1', 'top2', 'mean', 'threshold'")

    ## check args
    if (type %in% c("threshold")) {
        if (!(all(names(l) %in% names(args)))) {
            stop("'args' does not contain entries for all 'model's in ", 
                "'statistical'")
        }
  
        if (!"threshold" %in% names(args) && length(args$threshold) != 1) {
            stop("'args' does not contain entry 'threshold' of length 1")
        }
    }

    if (type %in% c("top1", "top2", "mean")) {
        if (! ("n"  %in% names(args) && length(args$n) == 1 && 
            is.numeric(args$n)) )
            stop("args does not contain the numeric entry `n` of length 1")
    }

    if (type == "threshold") {
        ## iterate through the list and remove the links below or above the 
        ## threshold and write to list
        l <- lapply(seq_along(l), function(x) {

            ## find corresponding model in l 
            name_x <- names(l)[x]

            ## get corresponding threshold in args
            threshold_x <- args[[names(l)[x]]]

            ## get corresponding adjacency matrix in l
            l_x <- l[[name_x]]

            ## for pearson/spearman correlation models (incl. partial and
            ## semi-partial), lasso, randomForest, clr, aracne and bayes higher
            ## values corresond to higher confidence
            ## only assign 1 to values that are above the threshold
            ifelse(l_x > threshold_x, 1, 0)
        })

        ## allow for compatibility of arguments
        ## calculate consenses from the binary matrices
        cons <- threeDotsCall(sna::consensus, dat = l, ...)

        ## threshold consensus that it is a binary matrix
        cons <- ifelse(cons >= args$threshold, 1, 0)

        rownames(cons) <- colnames(cons) <- colnames(l[[1]])

    } else { ## if type is in "top1", "top2" or "mean"
        l_df <- lapply(seq_along(l), function(x) {

            ## find corresponding model in l
            name_x <- names(l)[x]

            ## get corresponding adjacency matrix in l
            l_x <- l[[name_x]]

            ## for pearson/spearman correlation (incl. partial and
            ## semi-partial), lasso, randomForest, clr, aracne and bayes
            ## higher values corresond to higher confidence
            if (grepl(name_x, pattern = "lasso|randomForest|bayes")) {
                ## set values that are equal to 0 to NaN (values that are 0)
                ## do not explain the variability
                res <- getLinks(l_x, exclude = "== 0")
            } 
            if (grepl(name_x, pattern = "pearson|spearman|clr|aracne")) {
                res <- getLinks(l_x, exclude = NULL)
            }

            res
        })

        names(l_df) <- names(l)

        ## bind together the ranks of the models, stored in l_df
        ranks <- lapply(l_df, function(x) x$rank)
        ranks <- do.call("cbind", ranks)
        colnames(ranks) <- names(l_df)

        ## calculate the consensus information, i.e. either get the first or
        ## second top rank per row or calculate the average across rows
        ## depending on the type argument
        cons_val <- topKnet(ranks, type)

        ## bind row and col information with cons information
        row_col <- l_df[[1]][, c("row", "col")]
        ranks <- cbind(row_col, cons_val)

        ## get the top N features
        n <- args$n
        top_n <- sort(unique(cons_val))[1:n]
        ranks_top <- ranks[cons_val %in% top_n, ]

        ## write links in ranks_top to binary adjacency matrix cons
        cons <- matrix(0, nrow = ncol(l[[1]]), ncol = ncol(l[[1]]))
        rownames(cons) <- colnames(cons) <- colnames(l[[1]])
        cons[as.numeric(rownames(ranks_top))] <- 1

    }

    return(cons)

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
#' ranks <- matrix(c(c(1, 2, 3), c(2, 1, 3)), ncol = 2)
#'
#' ## type = "top1"
#' MetNet:::topKnet(ranks = ranks, type = "top1")
#'
#' ## type = "top2"
#' MetNet:::topKnet(ranks = ranks, type = "top2")
#'
#' ## type = "mean"
#' MetNet:::topKnet(ranks = ranks, type = "mean")
topKnet <- function(ranks, type) {

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
                min(x, na.rm = TRUE)
            } else {NaN}
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
                sort(x)[2] 
            } else {NaN}
        })
    }

    if (type == "mean") {
        ## get the average of all ranks
        cons_val <- apply(ranks, 1, mean, na.rm = TRUE)
    }

    return(cons_val)
}

#' @name threeDotsCall
#'
#' @aliases threeDotsCall
#'
#' @title Check if passed arguments match the function's formal arguments and
#' call the function with the checked arguments
#'
#' @description 
#' The function `threeDotsCall` gets the formal arguments
#' of a function `fun` and checks if the passed arguments `...`
#' matches the formal arguments. `threeDotsCall` will call the function
#' `fun` with the filtered arguments and will return the result of the function
#' call and the given arguments.
#'
#' @param fun `function` to check for arguments and to call
#'
#' @param ... arguments to be tested to be passed to `fun`
#'
#' @details
#' Used internally in `lasso`, `randomForest`, `bayes`,
#' `statistical` and `threshold`.
#'
#' `threeDotsCall` will not remove duplicated arguments and throw an error.
#'
#' @return Returned object given the function call with passed arguments
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' MetNet:::threeDotsCall(stats::sd, x = 1:10, y = 1:10)
#' ## in contrast to the above example, the following example will result in an
#' ## error
#' \dontrun{stats::sd(x = 1:10, y = 1:10)}
threeDotsCall <- function(fun, ...) {

    formal_args <- formalArgs(fun)
    args <- list(...)
    if (any(duplicated(names(args)))) stop("duplicated args in ...")

    input <- args[names(args) %in% formal_args]

    ## call the function
    res <- do.call(fun, input)
    return(res)
}