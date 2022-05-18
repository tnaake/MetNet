#' @name lasso
#' 
#' @aliases lasso
#' 
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
#' \dontrun{lasso(x = x_z, PFER = 0.95, cutoff = 0.95)}
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom stabs stabsel.matrix glmnet.lasso
#' 
#' @export
lasso <- function(x, parallel = FALSE, ...) {

    ## define list of arguments
    args_l <- list(...)
    args_l[["fitfun"]] <- stabs::glmnet.lasso
    args_l[["args.fitfun"]] <- list("alpha" = 1)
    
    ## x should be z-scaled
    if (parallel) {
        l1 <- BiocParallel::bplapply(seq_len(nrow(x)), function(i) {
            x_l1 <- t(x[-i, ])
            y_l1 <- x[i, ]

            ## lasso: alpha set to 1
            ## allow for compatibility of arguments
            args_l[["x"]] <- as.matrix(x_l1)
            args_l[["y"]] <- y_l1
            l1 <- do.call("stabsel.matrix", args_l)

            ## return selection probabilities of features that are not 0
            return(l1$max[l1$max != 0])
        })

    } else {
        l1 <- lapply(seq_len(nrow(x)), function(i) {
            x_l1 <- t(x[-i, ])
            y_l1 <- x[i, ]

            ## lasso: alpha set to 1
            ## allow for compatibility of arguments
            args_l[["x"]] <- as.matrix(x_l1)
            args_l[["y"]] <- y_l1
            l1 <- do.call("stabsel.matrix", args_l)

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
#' @param x matrix, where columns are the samples and the rows are features
#' (metabolites), cell entries are intensity values
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

    ## define a list with arguments
    args_l <- list(...)
    args_l[["exprMatrix"]] <- x
    args_l[["regulators"]] <- NULL
    args_l[["targets"]] <- NULL
    args_genie3 <- names(formals("GENIE3"))
    args_l <- args_l[names(args_l) %in% args_genie3]
    
    ## GENIE3 returns the importance of the link from "regulator gene" i to
    ## target gene "j" in the form of a weighted adjacency matrix
    ## set regulators and targets to NULL that they cannot be changed
    rf <- do.call("GENIE3", args_l)

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
#' @param mi matrix, where columns are samples and the rows are features
#' (metabolites), cell entries are mutual information values between the
#' features. As input, the mutual information (e.g. raw MI estimates) from the 
#' `knnmi.all` function of the `parmigene` package can be used.
#' @param ... not used here
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
#' x_z <- apply(x, 1, function(y) (y - mean(y)) / sd(y))
#' mi_x_z <- parmigene::knnmi.all(x_z)
#' clr(mi_x_z)
#'
#' @importFrom parmigene clr
#' @importFrom stats sd
#' 
#' @export
clr <- function(mi, ...) {

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
#' @param mi matrix, where columns are the samples and the rows are features
#' (metabolites), cell entries are mutual information values between the
#' features. As input, the mutual information (e.g. raw MI estimates) from the 
#' `knnmi.all` function of the `parmigene` package can be used.
#'
#' @param eps numeric, used to remove the weakest edge of each triple of nodes
#' @param ... not used here
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
#' x_z <- apply(x, 1, function(y) (y - mean(y)) / sd(y))
#' mi_x_z <- parmigene::knnmi.all(x_z)
#' aracne(mi_x_z, eps = 0.05)
#'
#' @importFrom parmigene aracne.a
#' @importFrom stats sd
#'
#' @export
aracne <- function(mi, eps = 0.05, ...) {

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
#' `psych` package) or partialCorrelation. `correlation` extracts the 
#' reported pair-wise correlation coefficients from the function 
#' `corr.test` and `partialCorrelation` and will return
#' the weighted adjacency matrix of the correlation coefficients, together 
#' with the associated p-values.
#'
#' @param
#' x `matrix`, where columns are the samples and the rows are features
#' (metabolites), cell entries are intensity values
#'
#' @param
#' method `character`, either "pearson", "spearman", "pearson_partial",
#' "spearman_partial", or "ggm".
#' 
#' @param 
#' p.adjust `character`, method of p-value adjustment passed to `p.adjust`
#' 
#' @param
#' ... additional arguments passed to `corr.test` or `partialCorrelation`
#'
#' @details
#' If `"pearson"` or `"spearman"` is used as a `method`, the function
#' `corr.test` from `psych` will be employed.
#' 
#' If `"ggm"` is used as a `method`, the function `ggm.estimate.pcor` from
#' `GeneNet` will be employed.
#'
#' If `"pearson_partial"` or `"spearman_partial"` is used as a `method` the
#' function `partialCorrelation` will be employed.
#'
#' `method` will be passed to argument `method` in `corr.test`
#' (in the case of `"pearson"` or `"spearman"`) or to `method` in 
#' `partialCorrelation` (`"pearson"` and `"spearman"` for `"pearson_partial"` 
#' and `"spearman_partial"`, respectively).
#'
#' @return
#' `list` containing two matrices, 
#' the first matrix contains correlation coefficients and 
#' the second matrix contains the corresponding p-values as obtained from the 
#' correlation algorithms `corr.test` or `partialCorrelation` (depending on the 
#' chosen `method`) and optionally the adjusted p.values (argument
#' `p.adjust`)
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com},
#' Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
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
#' @importFrom stats p.adjust
#' @importFrom GeneNet ggm.estimate.pcor n2kappa cor0.test
correlation <- function(x, method = "pearson", p.adjust = "none", ...) {

    ## define list with arguments
    args_l <- list(...)
    args_l[["method"]] <- method
    
    ## for pearson/spearman correlation
    if (method %in% c("pearson", "spearman")) {
        args_l[["x"]] <- t(x)
        args_l[["adjust"]] <- "none"
        args_l[["ci"]] <- FALSE
        args_corrtest <- names(formals("corr.test"))
        args_l <- args_l[names(args_l) %in% args_corrtest]
        cor_mat <- do.call("corr.test", args_l)
        
        ## calculate adjusted p-values
        cor_mat$p <- matrix(
            stats::p.adjust(as.vector(cor_mat$p), method = p.adjust),
            ncol = ncol(cor_mat$p), nrow = nrow(cor_mat$p), byrow = TRUE)
        cor_mat <- list(r = cor_mat[["r"]], p = cor_mat[["p"]])
        
        ## assign col- and rownames to cor_mat
        colnames(cor_mat[[1]]) <- rownames(cor_mat[[1]]) <- rownames(x)
        colnames(cor_mat[[2]]) <- rownames(cor_mat[[2]]) <- rownames(x)
    }

    ## for partial pearson/spearman correlation
    if (method %in% c("pearson_partial", "spearman_partial", "ggm")) {
        cor_mat <- list(r = matrix(NA, nrow = nrow(x), ncol = nrow(x)),
                        p = matrix(NA, nrow = nrow(x), ncol = nrow(x)))
        
        ## assign col- and rownames to cor_mat
        colnames(cor_mat[[1]]) <- rownames(cor_mat[[1]]) <- rownames(x)
        colnames(cor_mat[[2]]) <- rownames(cor_mat[[2]]) <- rownames(x)
        
        ## remove the rows that have NA values
        x_na <- na.omit(x)
        
        if (method %in% c("pearson_partial", "spearman_partial")) {
            if (method == "pearson_partial") {
                method <- "pearson"
            } else {
                method <- "spearman"
            }
            cor_mat_na <- partialCorrelation(x = t(x_na), method = method, ...)
            cor_mat_na$p <- matrix(
                stats::p.adjust(as.vector(cor_mat_na$p), method = p.adjust), 
                ncol = ncol(cor_mat_na$p), nrow = nrow(cor_mat_na$p), 
                byrow = TRUE)
        }
        
        ## for correlation based on graphical Gaussian models (ggm)
        if (method == "ggm") {
            
            cor_mat_na <- GeneNet::ggm.estimate.pcor(x = t(x_na), 
                method = "static")
            cor_mat_na <- cor_mat_na[seq_len(nrow(x_na)), seq_len(nrow(x_na))] 
            
            # calculate p-values
            # n2kappa converts sample size to the corresponding degree of freedom
            kappa <- GeneNet::n2kappa(n = length(x_na), p = nrow(x_na))
            p <- GeneNet::cor0.test(r = cor_mat_na, kappa = kappa, 
                method = "student")
            cor_mat_na <- list("r" = cor_mat_na, "p" = p)
            
            cor_mat_na$p <- matrix(
                stats::p.adjust(as.vector(cor_mat_na$p), method = p.adjust),
                ncol = ncol(cor_mat_na$p), nrow = nrow(cor_mat_na$p), 
                byrow = TRUE)
            
        }
        ## assign col- and rownames to cor_mat
        colnames(cor_mat_na[[1]]) <- rownames(cor_mat_na[[1]]) <- rownames(x_na)
        colnames(cor_mat_na[[2]]) <- rownames(cor_mat_na[[2]]) <- rownames(x_na)
        
        ## overwrite values from cor_mat_na into cor_mat
        cor_mat$r[rownames(cor_mat_na$r), colnames(cor_mat_na$r)] <- cor_mat_na$r
        cor_mat$p[rownames(cor_mat_na$p), colnames(cor_mat_na$p)] <- cor_mat_na$p
        
    }

    return(cor_mat)
}

#' @name partialCorrelation
#'
#' @aliases partialCorrelation
#'
#' @title Calculate the partial correlation and p-values
#'
#' @description
#' `partialCorrelation` infers an adjacency matrix of partial correlation 
#' values and associated p-values using using the `partial.r` function (from the
#' `psych` package). `partialCorrelation` calculates the p-values from the 
#' number of samples (`n`) and the number of controlling variables (`g`). 
#' The function will return a list containing the 
#' weighted adjacency matrix of the correlation values, together with the 
#' associated p-values.
#'
#' @param
#' x `matrix`, where columns are the features (metabolites) and the rows are 
#' samples, cell entries are intensity values
#'
#' @param
#' method `character`, either "pearson", "spearman"
#' 
#' @param 
#' ... further arguments passed to `partial.r` from `psych`
#' 
#' @details
#' The correlation coefficients $r_{ij|S}$ are obtained from `partial.r`
#' (`psych` package).
#' 
#' The t-values are calculated via
#' 
#' \eqn{t_{ij|S} = r_{ij|S} \cdot \sqrt{\frac{n-2-g}{1-r_{ij|S}^2}}},
#' where $n$ are the number of samples and $g$ the number of controlling
#' variables (number of features - 2).
#' 
#' The p-values are calculated as follows
#' \eqn{p_{ij|S} = 2 \cdot pt(-abs(t_{ij|S}), df = n - 2 - g)}
#' 
#' @return
#' `list` containing two matrices, 
#' the first matrix contains correlation coefficients and 
#' the second matrix contains the corresponding p-values
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' x <- t(x)
#' partialCorrelation(x, use = "pairwise", method = "pearson")
#'
#' @export
#' 
#' @importFrom psych partial.r
#' @importFrom stats pt
partialCorrelation <- function(x, method = "pearson", ...) {
    
    ## create a list with arguments
    args_l <- list(...)
    args_l[["data"]] <- x
    args_l[["method"]] <- method
    args_partialr <- names(formals("partial.r"))
    args_l <- args_l[names(args_l) %in% args_partialr]
    
    ## calculate the partial correlation coefficients
    r <- do.call("partial.r", args_l)
    
    ## obtain n, the sample size
    n <- nrow(x)
    
    ## obtain g, the number of controlling variables
    g <- ncol(x) - 2
    
    ## calculate the t-values and p-values
    df <- n - 2 - g
    t <- r * sqrt(df / (1 - r^2))
    diag(t) <- 0
    p <- 2 * stats::pt(-abs(t), df)
    
    return(list(r = r, p = p))
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

    ## create list of arguments that is compatible with arguments of 
    ## boot.strenth
    args_l <- list(...)
    args_l[["data"]] <- x_df
    args_l[["algorithm"]] <- algorithm
    args_l[["R"]] <- R
    args_bootstrength <- names(formals("boot.strength"))
    args_l <- args_l[names(args_l) %in% args_bootstrength]
    
    ## allow for compatibility of arguments by using do.call
    strength <- do.call("boot.strength", args_l)

    ## create empty bs_mat to be filled with connections
    bs_mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    colnames(bs_mat) <- rownames(bs_mat) <- make.names(rownames(x))

    ## write to bs_mat
    for (i in seq_len(nrow(strength))) {
        bs_mat[strength[i, "from"], strength[i, "to"]] <- strength[i, "strength"]
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
#' cellular networks (ARACNE), Pearson correlation (also partial), 
#' Spearman correlation (also partial)
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
#' `"pearson_partial"`, `"spearman"`, `"spearman_partial"`, `ggm`, `"bayes"`)
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
#' cellular networks (ARACNE), Pearson correlation (also partial), 
#' Spearman correlation (also partial) and Constraint-based structure learning 
#' (Bayes).
#'
#' `statistical` calls the function
#' `lasso`, `randomForest`, `clr`, `aracne`,
#' `correlation` (for `"pearson"`, `"pearson_partial"`, `"spearman"`, 
#' `"spearman_partial"`, `"ggm"`) and/or `bayes`
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
#' @importFrom stats sd
#' @importFrom parmigene knnmi.all
statistical <- function(x, model, ...) {

    ## check if model complies with the implemented model and return error
    ## if not so
    if (!(all(model %in% c("lasso", "randomForest", "clr", "aracne",
            "pearson", "pearson_partial", "spearman", "spearman_partial", 
            "ggm", "bayes"))))
        stop("'model' not implemented in statistical")

    ## check if x is numeric matrix and return error if not so
    if (mode(x) != "numeric") stop("x is not a numerical matrix")

    ## z-scale x and transpose
    x_z <- apply(x, 1, function(y) {
        (y - mean(y, na.rm = TRUE)) / stats::sd(y, na.rm = TRUE)
    })
    x_z <- t(x_z)

    ## create list to store results of running the models
    l <- list()
    
    ## create list to store the arguments in ...
    args_l <- list(...)

    ## add entry for lasso if "lasso" is in model
    if ("lasso" %in% model) {
        args_l[["x"]] <- x_z
        res <- do.call("lasso", args_l)
        diag(res) <- NaN
        l <- addToList(l, "lasso_coef", res)
        print("lasso finished")
    }

    ## add entry for randomForest if "randomForest" is in model
    if ("randomForest" %in% model) {
        args_l[["x"]] <- x
        res <- do.call("randomForest", args_l)
        diag(res) <- NaN
        res <- res[rownames(x), rownames(x)]
        l <- addToList(l, "randomForest_coef", res)
        print("randomForest finished.")
    }

    ## calculate mutual information if "clr" or "aracne" is in model
    if (any(c("clr", "aracne") %in% model)) {
        mi_x_z <- parmigene::knnmi.all(x_z)
    }

    ## add entry for clr if "clr" is in model
    if ("clr" %in% model) {
        args_l[["mi"]] <- mi_x_z
        res <- do.call("clr", args_l)
        diag(res) <- NaN
        l <- addToList(l, "clr_coef", res)
        print("clr finished.")
    }

    ## add entry for aracne if "aracne" is in model
    if ("aracne" %in% model) {
        args_l[["mi"]] <- mi_x_z
        res <- do.call("aracne", args_l)
        diag(res) <- NaN
        l <- addToList(l, "aracne_coef", res)
        print("aracne finished.")
    }

    ## add entry for pearson if "pearson" is in model
    if ("pearson" %in% model) {
        args_l[["x"]] <- x
        args_l[["method"]] <- "pearson"
        res <- do.call("correlation", args_l)
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
        args_l[["x"]] <- x
        args_l[["method"]] <- "pearson_partial"
        res <- do.call("correlation", args_l)
        pearson_p_coef <- res[["r"]]
        diag(pearson_p_coef) <- NaN
        pearson_p_pvalue <- res[["p"]]
        diag(pearson_p_pvalue) <- NaN
        l <- addToList(l, "pearson_partial_coef", pearson_p_coef)
        l <- addToList(l, "pearson_partial_pvalue", pearson_p_pvalue)
        print("pearson_partial finished.")
    } ## estimate, p.value


    ## add entry for spearman if "spearman" is in model
    if ("spearman" %in% model) {
        args_l[["x"]] <- x
        args_l["method"] <- "spearman"
        res <- do.call("correlation", args_l)
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
        args_l[["x"]] <- x
        args_l[["method"]] <- "spearman_partial"
        res <- do.call("correlation", args_l)
        spearman_p_coef <- res[["r"]]
        diag(spearman_p_coef) <- NaN
        spearman_p_pvalue <- res[["p"]]
        diag(spearman_p_pvalue) <- NaN
        l <- addToList(l, "spearman_partial_coef", spearman_p_coef)
        l <- addToList(l, "spearman_partial_pvalue", spearman_p_pvalue)
        print("spearman_partial finished.")
    }

    
    ## add entry for ggm if "ggm" is in model
    if ("ggm" %in% model) {
        args_l[["x"]] <- x
        args_l[["method"]] <- "ggm"
        res <- do.call("correlation", args_l)
        ggm_coef <- res[["r"]]
        diag(ggm_coef) <- NaN
        ggm_pvalue <- res[["p"]]
        diag(ggm_pvalue) <- NaN
        l <- addToList(l, "ggm_coef", ggm_coef)
        l <- addToList(l, "ggm_pvalue", ggm_pvalue)
        print("ggm finished.")
    }

    ## add entry for bayes if "bayes" is in model
    if ("bayes" %in% model) {
        args_l[["x"]] <- x
        res <- do.call("bayes", args_l)
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
