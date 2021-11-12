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
#' @importFrom BiocParallel bplapply
#' @importFrom stabs stabsel.matrix glmnet.lasso
#' 
#' @export
lasso <- function(x, parallel = FALSE, ...) {

    ## x should be z-scaled
    if (parallel) {
        l1 <- BiocParallel::bplapply(seq_len(nrow(x)), function(i) {
            x_l1 <- t(x[-i, ]); y_l1 <- x[i, ]

            ## lasso: alpha set to 1
            ## allow for compatibility of arguments
            l1 <- threeDotsCall("stabsel.matrix", x = as.matrix(x_l1),
                    y = y_l1, fitfun = stabs::glmnet.lasso,
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
                    y = y_l1, fitfun = stabs::glmnet.lasso,
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
#' @importFrom GENIE3 GENIE3
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
#' @importFrom parmigene clr
#' @importFrom stats sd
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
#' @importFrom parmigene aracne.a
#' @importFrom stats sd
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
#' correlation using the `corr.test` function (from the
#' `psych` package), `pcor` (from `ppcor`) or
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
#' method `character`, either "pearson", "spearman", "pearson_partial",
#' "spearman_partial", "pearson_semipartial" or "spearman_semipartial".
#' 
#' @param 
#' p.adjust `character`, method of p-value adjustment passed to `p.adjust`
#'
#' @details
#' If `"pearson"` or `"spearman"` is used as a `method`, the function
#' `corr.test` from `psych` will be employed.
#'
#' If `"pearson_partial"` or `"spearman_partial"` is used as a `method` the
#' function `pcor` from `spcor` will be employed.
#'
#' If `"pearson_semipartial"` or `"spearman_semipartial"` is used as a
#' `method` the function `spcor` from `spcor` will be employed.
#'
#' `method` will be passed to argument `method` in `corr.test`
#' (in the case of `"pearson"` or `"spearman"`) or to `method` in `pcor`
#' (`"pearson"` and `"spearman"` for `"pearson_partial"` and
#' `"spearman_partial"`, respectively) or to `method` in `spcor`
#' (`"pearson"` or `"spearman"` for `"pearson_semipartial"` and
#' `"spearman_semipartial"`, respectively).
#'
#' @return
#' `list` containing two matrices, 
#' the first matrix contains correlation coefficients and 
#' the second matrix contains the corresponding p-values as obtained from the 
#' correlation algorithms `corr.test`, `pcor` or `spcor` (depending on the 
#' chosen `method`) and optionally the adjusted p.values (argument
#' `p.adjust`)
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' correlation(x, method = "pearson")
#'
#' @export
#' 
#' @importFrom psych corr.test
#' @importFrom ppcor pcor spcor
#' @importFrom stats p.adjust
correlation <- function(x, method = "pearson", p.adjust = "none") {

    ## for pearson/spearman correlation
    if (method %in% c("pearson", "spearman")) {
        cor_mat <- psych::corr.test(x = t(x), method = method, adjust = "none")
        cor_mat$p <- matrix(
            stats::p.adjust(as.vector(cor_mat$p), method = p.adjust),
            ncol = ncol(cor_mat$p), nrow = nrow(cor_mat$p), byrow = TRUE)
        cor_mat <- list(r = cor_mat[["r"]], p = cor_mat[["p"]])
    }

    ## for partial pearson/spearman correlation
    if (method %in% c("pearson_partial", "spearman_partial")) {
        if (method == "pearson_partial") {
            method <- "pearson"
        } else {
            method <- "spearman"
        }
        cor_mat <- ppcor::pcor(t(x), method = method)
        cor_mat$p.value <- matrix(
            stats::p.adjust(as.vector(cor_mat$p.value), method = p.adjust),
            ncol = ncol(cor_mat$p.value), nrow = nrow(cor_mat$p.value), 
            byrow = TRUE)
    }

    ## for semipartial pearson/spearman corelation
    if (method %in% c("pearson_semipartial", "spearman_semipartial")) {
        if (method == "pearson_semipartial") {
            method <- "pearson"
        } else {
            method <- "spearman"
        }
        cor_mat <- ppcor::spcor(t(x), method = method)
        cor_mat$p.value <- matrix(
            stats::p.adjust(as.vector(cor_mat$p.value), method = p.adjust),
            ncol = ncol(cor_mat$p.value), nrow = nrow(cor_mat$p.value), 
            byrow = TRUE)
    }
    
    ## for correlation based on graphical Gaussian models (ggm)
    if (method == "ggm") {
        # remove NaN from dataset
        x <- na.omit(x)
        
        cor_mat <- GeneNet::ggm.estimate.pcor(t(x), method = "static")
        cor_mat <- cor_mat[seq_len(nrow(x)), seq_len(nrow(x))] 
        
        # calculate p-values
        
        # n2kappa converts sample size to the corresponding degree of freedom
        kappa <- GeneNet::n2kappa(n = length(x), p = nrow(x))
        p <- GeneNet::cor0.test(r = cor_mat, kappa = kappa, method = "student")
        cor_mat <- list("estimate" = cor_mat, "p" = p)
        
        cor_mat$p <- matrix(
            stats::p.adjust(as.vector(cor_mat$p), method = p.adjust),
            ncol = ncol(cor_mat$p), nrow = nrow(cor_mat$p), 
            byrow = TRUE)
        
    }

    ## assign col- and rownames to cor_mat
    colnames(cor_mat[[1]]) <- rownames(cor_mat[[1]]) <- rownames(x)
    colnames(cor_mat[[2]]) <- rownames(cor_mat[[2]]) <- rownames(x)

    # ## get absolute values
    # cor_mat <- abs(cor_mat)

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
#' @importFrom bnlearn boot.strength
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
    for (i in seq_len(nrow(strength))) {
        tmp <- as.character(strength[i, ])
        names(tmp) <- names(strength[i, ])
        bs_mat[tmp["from"], tmp["to"]] <- tmp["strength"]
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
#' cor_pearson <- correlation(x, method = "pearson")
#' cor_spearman <- correlation(x, method = "spearman")
#' l <- list(pearson = cor_pearson)
#' MetNet:::addToList(l, "spearman_coef", cor_spearman$r)
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
#' @title Create an `AdjacencyMatrix` object containing assays of adjacency 
#' matrices from statistical methods
#'
#' @description
#' The function `statitical` infers adjacency matrix topologies from
#' statistical methods and returns matrices of these networks in an
#' `AdjacencyMatrix` object. The
#' function includes functionality to calculate adjacency matrices based on
#' LASSO (L1 norm)-regression, random forests, context likelihood of
#' relatedness (CLR), the algorithm for the reconstruction of accurate
#' cellular networks (ARACNE), Pearson correlation (also partial and
#' semipartial), Spearman correlation (also partial and semipartial)
#' and score-based structure learning (Bayes). The function returns an
#' `AdjacencyMatrix` object of adjacency matrices that are defined by `model`.
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
#' specified methods and will return an `AdjacencyMatrix` containing the weighted
#' adjacency matrices in the `assays` slot.
#'
#' Internally `x` will be z-scaled and the z-scaled object
#' will be used in `lasso`, `clr` and/or `aracne`.
#' 
#' The slot `type` is set to `statistical`. The slot `directed` is set to
#' `TRUE` if the methods `"lasso"`, `"randomForest"`, or `"bayes"` were used, 
#' otherwise `directed` is set to `FALSE`.
#' The slot `threshold` is set to `FALSE`.
#'
#' @return `AdjacencyMatrix` containing the respective adjacency matrices in the
#' `assay` slot as specified by `model`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' statistical(x = x, model = c("pearson", "spearman"))
#' statistical(x = x, model = c("pearson", "spearman"), p.adjust = "BH")
#'
#' @export
#' 
#' @importFrom S4Vectors DataFrame
#' @importFrom mpmi cmi
#' @importFrom stats sd
statistical <- function(x, model, ...) {

    ## check if model complies with the implemented model and return error
    ## if not so
    if (!(all(model %in% c("lasso", "randomForest", "clr", "aracne",
            "pearson", "pearson_partial", "pearson_semipartial",
            "spearman", "spearman_partial", "spearman_semipartial", "bayes", "ggm"))))
        stop("'model' not implemented in statistical")

    ## check if x is numeric matrix and return error if not so
    if (mode(x) != "numeric") stop("x is not a numerical matrix")

    ## z-scale x and transpose
    x_z <- apply(x, 1, function(y) {
        (y - mean(y, na.rm = TRUE)) / stats::sd(y, na.rm = TRUE)
    })
    x_z <- t(x_z)

    l <- list()

    ## add entry for lasso if "lasso" is in model
    if ("lasso" %in% model) {
        res <- lasso(x = x_z, ...)
        diag(res) <- NaN
        l <- addToList(l, "lasso_coef", res)
        print("lasso finished")
    }

    ## add entry for randomForest if "randomForest" is in model
    if ("randomForest" %in% model) {
        res <- randomForest(x = x, ...)
        diag(res) <- NaN
        res <- res[rownames(x), rownames(x)]
        l <- addToList(l, "randomForest_coef", res)
        print("randomForest finished.")
    }

    ## calculate mutual information if "clr" or "aracne" is in model
    if (any(c("clr", "aracne") %in% model)) {
        mi_x_z <- mpmi::cmi(t(x_z))$bcmi
        rownames(mi_x_z) <- colnames(mi_x_z) <- rownames(x)
    }

    ## add entry for clr if "clr" is in model
    if ("clr" %in% model) {
        res <- threeDotsCall("clr", mi = mi_x_z, ...)
        diag(res) <- NaN
        l <- addToList(l, "clr_coef", res)
        print("clr finished.")
    }

    ## add entry for aracne if "aracne" is in model
    if ("aracne" %in% model) {
        res <- threeDotsCall("aracne", mi = mi_x_z, ...)
        diag(res) <- NaN
        l <- addToList(l, "aracne_coef", res)
        print("aracne finished.")
    }

    ## add entry for pearson if "pearson" is in model
    if ("pearson" %in% model) {
        res <- threeDotsCall("correlation", x = x, method = "pearson", ...)
        pearson_coef <- res[["r"]]
        diag(pearson_coef) <- NaN
        pearson_pvalue <- res[["p"]]
        diag(pearson_pvalue) <- NaN
        l <- addToList(l, "pearson_coef", pearson_coef)
        l <- addToList(l, "pearson_pvalue", pearson_pvalue)
        print("pearson finished.")
    }

    ## add entry for pearson_partial if "pearson_partial" is in model
    if ("pearson_partial" %in% model) {
        res <- threeDotsCall("correlation", x = x,
            method = "pearson_partial", ...)
        pearson_p_coef <- res[["estimate"]]
        diag(pearson_p_coef) <- NaN
        pearson_p_pvalue <- res[["p.value"]]
        diag(pearson_p_pvalue) <- NaN
        l <- addToList(l, "pearson_partial_coef", pearson_p_coef)
        l <- addToList(l, "pearson_partial_pvalue", pearson_p_pvalue)
        print("pearson_partial finished.")
    } ## estimate, p.value

    ## add entry for pearson_semipartial if "pearson_semipartial" is in model
    if ("pearson_semipartial" %in% model) {
        res <- threeDotsCall("correlation", x = x,
            method = "pearson_semipartial", ...)
        pearson_sp_coef <- res[["estimate"]]
        diag(pearson_sp_coef) <- NaN
        pearson_sp_pvalue <- res[["p.value"]]
        diag(pearson_sp_pvalue) <- NaN
        l <- addToList(l, "pearson_semipartial_coef", pearson_sp_coef)
        l <- addToList(l, "pearson_semipartial_pvalue", pearson_sp_pvalue)
        print("pearson_semipartial finished.")
    }

    ## add entry for spearman if "spearman" is in model
    if ("spearman" %in% model) {
        res <- threeDotsCall("correlation", x = x, method = "spearman", ...)
        spearman_coef <- res[["r"]]
        diag(spearman_coef) <- NaN
        spearman_pvalue <- res[["p"]]
        diag(spearman_pvalue) <- NaN
        l <- addToList(l, "spearman_coef", spearman_coef)
        l <- addToList(l, "spearman_pvalue", spearman_pvalue)
        print("spearman finished.")
    }

    ## add entry for spearman_partial if "spearman_partial" is in model
    if ("spearman_partial" %in% model) {
        res <- threeDotsCall("correlation", x = x, 
            method = "spearman_partial", ...)
        spearman_p_coef <- res[["estimate"]]
        diag(spearman_p_coef) <- NaN
        spearman_p_pvalue <- res[["p.value"]]
        diag(spearman_p_pvalue) <- NaN
        l <- addToList(l, "spearman_partial_coef", spearman_p_coef)
        l <- addToList(l, "spearman_partial_pvalue", spearman_p_pvalue)
        print("spearman_partial finished.")
    }

    ## add entry for spearman_semipartial if "spearman_semipartial" is in model
    if ("spearman_semipartial" %in% model) {
        res <- threeDotsCall("correlation", x = x,
            method = "spearman_semipartial", ...)
        spearman_sp_coef <- res[["estimate"]]
        diag(spearman_sp_coef) <- NaN
        spearman_sp_pvalue <- res[["p.value"]]
        diag(spearman_sp_pvalue) <- NaN
        l <- addToList(l, "spearman_semipartial_coef", spearman_sp_coef)
        l <- addToList(l, "spearman_semipartial_pvalue", spearman_sp_pvalue)
        print("spearman_semipartial finished.")
    }
    
    ## add entry for ggm if "ggm" is in model
    if ("ggm" %in% model) {
        res <- threeDotsCall("correlation", x = x, method = "ggm", ...)
        ggm_coef <- res[["estimate"]]
        diag(ggm_coef) <- NaN
        ggm_pvalue <- res[["p"]]
        diag(ggm_pvalue) <- NaN
        l <- addToList(l, "ggm_coef", ggm_coef)
        l <- addToList(l, "ggm_pvalue", ggm_pvalue)
        print("ggm finished.")
    }

    ## add entry for bayes if "bayes" is in model
    if ("bayes" %in% model) {
        res <- threeDotsCall("bayes", x = x, ...)
        diag(res) <- NaN
        l <- addToList(l, "bayes_coef", res)
        print("bayes finished.")
    }

    rD <- DataFrame(names = rownames(l[[1]]))
    rownames(rD) <- rownames(l[[1]])
    
    directed <- if (any(c("lasso", "randomForest", "bayes") %in% model)) {
        TRUE
    } else {
        FALSE
    }

    adj <- AdjacencyMatrix(l, rowData = rD,
        type = "statistical", directed = directed, thresholded = FALSE)
    
    return(adj)
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
#' The function `threshold` takes as input an `AdjacencyMatrix` object 
#' containing adjacency matrices
#' as returned from the function `statistical`. Depending on the `type`
#' argument, `threshold` will identify the strongest link that are
#' lower or higher a certain threshold (`type = "threshold"`) or
#' identify the top `n` links (`type` either `"top1`, `"top2` or `"mean"`).
#' It will return this kind of information as a binary matrix in the form
#' of an `AdjacencyMatrix` object. 
#'
#' @param am `AdjacencyMatrix` object of `type` `"statistical"` as 
#' created from the function `statistical`. The object will contain the
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
#' consensus matrix.
#' 
#' @param values `character`, take from the adjacency matrix all values ("all"),
#' the minimum of the pairs ("min") or the maximum ("max")
#' a^*_{ij} = min(a_ij, a_ji)
#' a^*_{ij} = max(a_ij, a_ji)
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
    args, values = c("all", "min", "max")) {

    ## check match.arg for values
    type <- match.arg(type)
    values <- match.arg(values)
    
    if (!is(am, "AdjacencyMatrix")) {
        stop("'am' is not an 'AdjacencyMatrix' object")
    }
    
    if (!validObject(am)) {
        stop("'am' must be a valid 'AdjacencyMatrix' object")
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
    
    if (type == "threshold") {
        
        df_filter <- as.data.frame(am) %>% 
            dplyr::filter(!!rlang::parse_expr(args$filter))
        
        inds_row <- match(df_filter[, "Row"], rownames(cons))
        inds_col <- match(df_filter[, "Col"], colnames(cons))
        
        ## write the remaining elements to the consensus matrix
        cons[cbind(inds_row, inds_col)] <- 1

    } else { ## if type is in "top1", "top2" or "mean"
        
        ind_coef <- grep(pattern = "_coef", x = assayNames(am))
        l <- as.list(assays(am)[ind_coef])

        l_df <- lapply(seq_along(l), function(x) {

            ## find corresponding model in l
            name_x <- names(l)[x]

            ## get corresponding adjacency matrix in l
            l_x <- l[[name_x]]

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
            ## semi-partial), lasso, randomForest, clr, aracne and bayes
            ## higher values corresond to higher confidence
            if (grepl(name_x, 
                pattern = "lasso_coef|randomForest_coef|bayes_coef")) {
                
                ## set values that are equal to 0 to NaN (values that are 0)
                ## do not explain the variability
                res <- getLinks(l_x, exclude = "== 0")
            }
            if (grepl(name_x, 
                pattern = "pearson_coef|pearson_partial_coef|pearson_semipartial_coef|spearman_coef|spearman_partial_coef|spearman_semipartial_coef|clr_coef|ggm_coef|aracne_coef")) {
                
                res <- getLinks(l_x, exclude = NULL)
            }

            res
        })

        ## bind together the ranks of the models, stored in l_df
        ranks <- lapply(l_df, function(x) x$rank)
        ranks <- do.call("cbind", ranks)
        colnames(ranks) <- names(l)

        ## calculate the consensus information, i.e. either get the first or
        ## second top rank per row or calculate the average across rows
        ## depending on the type argument
        cons_val <- topKnet(ranks, type)

        ## bind row and col information with cons information
        row_col <- l_df[[1]][, c("row", "col")]
        ranks <- cbind(row_col, cons_val)

        ## get the top N features
        top_n <- sort(unique(cons_val))[1:args$n]
        ranks_top <- ranks[cons_val %in% top_n, ]

        ## write links in ranks_top to binary adjacency matrix cons
        cons[as.numeric(rownames(ranks_top))] <- 1
    }

    ## assign the consensus matrix to a new slot
    assay(am, "consensus") <- cons
    
    if (type %in% c("top1", "top2", "mean") & values %in% c("min", "max")) 
        am@directed <- FALSE
    am@thresholded <- TRUE
    
    return(am)
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
                sort(x)[2]
            } else {
                NaN
            }
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
#' @importFrom methods formalArgs
#' @examples
#' MetNet:::threeDotsCall(stats::sd, x = 1:10, y = 1:10)
#' ## in contrast to the above example, the following example will result in an
#' ## error
#' \dontrun{stats::sd(x = 1:10, y = 1:10)}
threeDotsCall <- function(fun, ...) {

    formal_args <- methods::formalArgs(fun)
    args <- list(...)
    if (any(duplicated(names(args)))) stop("duplicated args in ...")

    input <- args[names(args) %in% formal_args]

    ## call the function
    res <- do.call(fun, input)
    return(res)
}
