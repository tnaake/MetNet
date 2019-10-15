#' @importFrom stabs stabsel.matrix glmnet.lasso
NULL
#' @importFrom GENIE3 GENIE3 ##rfPermute rfPermute rfPermute.formula rp.importance
NULL
#' @importFrom mpmi cmi
NULL
#' @importFrom parmigene clr aracne.a
NULL
#' @importFrom WGCNA corAndPvalue
NULL
#' @importFrom bnlearn fast.iamb arcs
NULL
#' @importFrom sna consensus
NULL
#' @importFrom BiocParallel bplapply
NULL
#' @importFrom stats formula p.adjust sd
NULL
#' @importFrom methods formalArgs
NULL
#' @importFrom ppcor pcor spcor
NULL

#' @name lasso
#' @aliases lasso
#' @title Create a adjacency matrix based on LASSO
#' @description  \code{lasso} infers a adjacency matrix using 
#' LASSO using the \code{stabsel.matrix} function from the 
#' \code{stabs} package. \code{lasso} extracts the  predictors from the 
#' function \code{stabsel.matrix} and writes the presence/absence 
#' of this connection to a matrix that is returned. 
#' @usage lasso(x, parallel = FALSE, ...)
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values 
#' @param parallel logical, should computation be parallelized? If 
#' \code{parallel = TRUE} the \code{bplapply} will be applied if 
#' \code{parallel = FALSE} the \code{lapply} function will be applied. 
#' @param ... parameters passed to \code{stabsel.matrix} 
#' @details For use of the parameters used in the \code{stabsel.matrix} function, 
#' refer to ?stabs::stabsel.matrix.
#' @return matrix, matrix with edges inferred from LASSO algorithm 
#' \code{stabsel.matrix}
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' x_z <- t(apply(x, 1, function(y) (y - mean(y)) / sd(y)))
#' \dontrun{lasso(x_z, PFER = 0.95, cutoff = 0.95)}
#' @export
lasso <- function(x, parallel = FALSE, ...) {
    ## x should be z-scaled
    if (parallel) {
        l1 <- bplapply(seq_len(nrow(x)), function(i) {
            x_l1 <- t(x[-i, ]); y_l1 <- x[i, ]
            
            ## lasso: alpha = 1
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
            
            ## lasso: alpha = 1
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
    for (i in seq_len(length(l1))) {l1_mat[names(l1[[i]]), i] <- l1[[i]]}
    
    return(l1_mat)
}

#' @name randomForest
#' @aliases randomForest
#' @title Create a adjacency matrix based on random forest
#' @description  \code{randomForest} infers an adjacency matrix using 
#' random forest using the \code{rfPermute} function from the 
#' \code{rfPermute} package. \code{randomForest} extracts the p-values  
#' by the function \code{rp.importance} and writes the presence/absence based 
#' on the significance value (\eqn{\alpha \leq 0.05}) of this 
#' connection to a matrix. The adjacency matrix is returned. 
#' @usage randomForest(x, parallel = FALSE, randomForest_adjust = "none", ...)
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values 
#' @param parallel logical, should computation be parallelized? If 
#' \code{parallel = TRUE} the \code{bplapply} will be applied if 
#' \code{parallel = FALSE} the \code{lapply} function will be applied. 
#' @param randomForest_adjust character, correction method for p-values from
#' \code{rp.importance}, \code{randomForest_adjust} will be passed to the 
#' \code{p.adjust} function and should be one of "holm", "hochberg", "hommel", 
#' "bonferroni", "BH", "BY", "fdr", "none"
#' @param ... parameters passed to \code{rfPermute.default} 
#' @details For use of the parameters used in the \code{rfPermute} function, 
#' refer to ?rfPermute::rfPermute.default.
#' @return matrix, matrix with edges inferred from random forest algorithm 
#' \code{rfPermute} and \code{rp.importance}
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' randomForest(x)
#' @export
randomForest <- function(x, ...) {#parallel = FALSE, 
                         #randomForest_adjust = "none", ...) {
    
    ## GENIE3 returns the importance of the link from "regulator gene" i to 
    ## target gene "j" in the form of a weighted adjacency matrix
    ## set regulators and targets to NULL that they cannot be changed
    rf <- threeDotsCall(GENIE3, exprMatrix = x, regulators = NULL, 
        targets = NULL, ...)
    
    return(rf)
}
#     if (parallel) {
#         rf <- bplapply(seq_len(nrow(x)), function(i) {
#             x_rf <- df_x[, -i]
#             y_rf <- df_x[, i]
#             ##formula_rf <- paste(rownames(x)[i], "~", ".")
#             ## allow for compatibility of arguments
#             rf <- threeDotsCall(rfPermute::rfPermute.default,
#                                 x = x_rf, y = y_rf, ...)
#             GENIE3(x)
#             rf_p <- rp.importance(rf)[,"IncNodePurity.pval"]
#             return(rf_p)
#         })
#     } else {
#         rf <- lapply(seq_len(nrow(x)), function(i) {
#             x_rf <- df_x[, -i]
#             y_rf <- df_x[, i]
#             ##formula_rf <- paste(rownames(x)[i], "~", ".")
#             ## allow for compatibility of arguments 
#             rf <- threeDotsCall(rfPermute::rfPermute.default,
#                                 x = x_rf, y = y_rf, ...)
#             rf_p <- rp.importance(rf)[,"IncNodePurity.pval"]
#             return(rf_p)
#         })
#     }
#     rf_mat <- matrix(1, nrow = nrow(x), ncol = nrow(x))
#     colnames(rf_mat) <- rownames(rf_mat) <- rownames(x)
#     
#     for (i in seq_len(length(rf))) {
#         rf_mat[names(rf[[i]]), rownames(x)[i]] <- rf[[i]]
#     }
#     rf_mat <- stats::p.adjust(rf_mat, method = randomForest_adjust)
#     rf_mat <- matrix(rf_mat, ncol = nrow(x), nrow = nrow(x), byrow = FALSE)
#     ##rf_mat <- ifelse(rf_mat > 0.05, 0, 1)
#     colnames(rf_mat) <- rownames(rf_mat) <- rownames(x)
#     
#     return(rf_mat)
# }
#' @name clr
#' @aliases clr
#' @title Create an adjacency matrix based on context likelihood or 
#' relatedness network
#' @description  \code{clr} infers an adjacency matrix using 
#' context likelihood/relatedness network using the \code{clr} function from 
#' the \code{parmigene} package. The presence/absence is based on if the 
#' returned value exceeds a user-defined threshold value. \code{clr} will 
#' return the adjacency matrix containing the presence/absence value.
#' @usage clr(mi, clr_threshold = 0)
#' @param mi matrix, where columns and the rows are features 
#' (metabolites), cell entries are mutual information values between the 
#' features. As input, the mutual information (e.g. raw MI estimates or 
#' Jackknife bias corrected MI estimates) from the \code{cmi} function of the
#' \code{mpmi} package can be used.
#' @param clr_threshold numeric, if the clr value exceeds the threshold 
#' (clr$_{i,j}$ > threshold, where clr$_{i, j}$ is the clr value of the ith row 
#' feature and of the jth column feature), the connection is defined as 
#' present, if the clr value is lower than the threshold value 
#' (clr$_{i,j}$ \eqn{\leq} threshold)
#' there is no statistical connection reported. 
#' @details For more details on the \code{clr} function, 
#' refer to ?parmigene::clr.
#' @return matrix, matrix with edges inferred from Context Likelihood or 
#' Relatedness Network algorithm \code{clr}
#' @author Thomas Naake, \email{thomasnaake @googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' x_z <- t(apply(x, 1, function(y) (y - mean(y)) / sd(y)))
#' mi_x_z <- mpmi::cmi(x_z)$bcmi
#' clr(mi_x_z, clr_threshold = 0)
#' @export
clr <- function(mi) {#, clr_threshold = 0) {
    ##if (!is.numeric(clr_threshold)) stop("clr_threshold is not numeric")
    clr_mat <- parmigene::clr(mi)
    ##clr_mat <- ifelse(clr_mat > clr_threshold, 1, 0)
    colnames(clr_mat) <- rownames(clr_mat) <- rownames(mi)
    return(clr_mat)
}

#' @name aracne
#' @aliases aracne
#' @title Create an adjacency matrix based on algorithm for the reconstruction
#' of accurate cellular networks 
#' @description  \code{.information} infers an adjacency matrix using 
#' the algorithm for the reconstruction of accurate cellular networks 
#' using the \code{aracne.a} function from the 
#' \code{parmigene} package. The presence/absence is based on if the 
#' returned value exceeds a user-defined threshold value. \code{aracne} will 
#' return the adjacency matrix containing the presence/absence value.
#' @usage aracne(mi, eps = 0.05, aracne_threshold = 0)
#' @param mi matrix, where columns and the rows are features 
#' (metabolites), cell entries are mutual information values between the 
#' features. As input, the mutual information (e.g. raw MI estimates or 
#' Jackknife bias corrected MI estimates) from the \code{cmi} function of the
#' \code{mpmi} package can be used.
#' @param eps numeric, used to remove the weakest edge of each triple of nodes
#' @param aracne_threshold numeric, if the aracne value exceeds the threshold 
#' (aracne$_{i,j}$ > threshold, where aracne$_{i, j}$ is the aracne value of 
#' the ith row feature and of the jth column feature), the connection is 
#' defined as present, if the aracne value is lower than the threshold value 
#' (aracne$_{i,j}$ <=  threshold) there is no statistical connection reported. 
#' @details For more details on the \code{aracne.a} function, 
#' refer to ?parmigene::aracne.a.
#' @return matrix, matrix with edges inferred from Reconstruction of accurate
#' cellular networks algorithm \code{aracne}
#' @author Thomas Naake, \email{thomasnaake @googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' x_z <- t(apply(x, 1, function(y) (y - mean(y)) / sd(y)))
#' mi_x_z <- mpmi::cmi(x_z)$bcmi
#' aracne(mi_x_z, eps = 0.05, aracne_threshold = 0)
#' @export
aracne <- function(mi, eps = 0.05) {#, aracne_threshold = 0) {
    #if (!is.numeric(aracne_threshold)) stop("aracne_threshold is not numeric")
    aracne_mat <- parmigene::aracne.a(mi, eps = eps)
    #aracne_mat <- ifelse(aracne_mat > aracne_threshold, 1, 0)
    colnames(aracne_mat) <- rownames(aracne_mat) <- rownames(mi)
    return(aracne_mat)
}

#' @name correlation 
#' @aliases correlation
#' @title Create an adjacency matrix based on correlation 
#' @description  \code{correlation} infers an adjacency matrix using 
#' correlation using the \code{corAndPvalue} function (from the 
#' \code{WGCNA} package), \code{pcor} (from \code{ppcor}) or 
#' \code{spcor} (from \code{ppcor}). \code{correlation} extracts the reported 
#' p-values from the function \code{corAndPvalue}, \code{pcor} or \code{spcor} 
#' that can be adjusted for 
#' multiple testing (\code{correlation_adjust} parameter) and will return 
#' an unweighted adjacency matrix containing edges if the (adjusted) p-value
#' is below the value defined by \code{correlation_threshold}. 
#' @usage correlation(x, correlation_adjust = "none", type = "pearson", 
#'                                         correlation_threshold = 0.05, ...) 
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values 
#' @param type character, either "pearson", "spearman", "pearson_partial",
#' "spearman_partial", "pearson_semipartial" or "spearman_semipartial". 
#' \code{type} will be passed to argument \code{method} in \code{corAndPvalue} 
#' (in the case of "pearson" or "spearman") or to \code{method} in \code{pcor} 
#' ("pearson" and "spearman" for "pearson_partial" and "spearman_partial", 
#' respectively) or to \code{method} in \code{spcor} ("pearson" or "spearman"
#' for "pearson_semipartial" and "spearman_semipartial", respectively)
#' @param correlation_adjust character 
#' @param correlation_threshold numeric, significance level \eqn{\alpha}
#' (default: 0.05), if the (adjusted) p-values exceed this value, there 
#' is no statistical connection between features 
#' @param ... parameters passed to \code{corAndPvalue} (argument \code{adjust} 
#' will be ignored)
#' @details If "pearson" or "spearman" is used as a \code{method} the function 
#' \code{corAndPvalue} from \code{WGCNA} will be employed. 
#' If "pearson_partial" or "spearman_partial" is used as a \code{method} the 
#' function \code{pcor} from \code{spcor} will be employed. 
#' If "pearson_semipartial" or "spearman_semipartial" is used as a 
#' \code{method} the function \code{spcor} from \code{spcor} will be employed. 
#' For use of the parameters used in the \code{corAndPvalue} function, 
#' refer to ?WGCNA::corAndPvalue. 
#' @return matrix, matrix with edges inferred from correlation algorithm 
#' \code{corAndPvalue}, \code{pcor} or \code{spcor} (depending on the chosen 
#' method)
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' correlation(x, correlation_adjust = "bonferroni", type = "pearson")
#' @export
correlation <- function(x, correlation_adjust = "none", type = "pearson", 
    correlation_threshold = 0.05, ...) {
    if (!is.numeric(correlation_threshold))
        stop("correlation_threshold is not numeric")
    ## get character vector for p-value adjustment
    adjust <- correlation_adjust
    
    ## allow for compatibility of arguments
    if (type %in% c("pearson", "spearman")) {
        cor_mat <- threeDotsCall(WGCNA::corAndPvalue, x = t(x),
            method = type, ...)$p
        cor_mat <- stats::p.adjust(cor_mat, method = adjust)
        cor_mat <- matrix(cor_mat, ncol = nrow(x), nrow = nrow(x),
            byrow = FALSE)
    }
    if (type %in% c("pearson_partial", "spearman_partial")) {
        if (type == "pearson_partial") method <- "pearson"
        if (type == "spearman_partial") method <- "spearman"
        cor_mat <- ppcor::pcor(t(x), method = method)$p.value
        cor_mat <- stats::p.adjust(cor_mat, method = adjust)
        cor_mat <- matrix(cor_mat, ncol = nrow(x), nrow = nrow(x),
            byrow = FALSE)
    }
    if (type %in% c("pearson_semipartial", "spearman_semipartial")) {
        if (type == "pearson_semipartial") method <- "pearson"
        if (type == "spearman_semipartial") method <- "spearman"
        cor_mat <- ppcor::spcor(t(x), method = method)$p.value
        cor_mat <- stats::p.adjust(cor_mat, method = adjust)
        cor_mat <- matrix(cor_mat, ncol = nrow(x), nrow = nrow(x), 
            byrow = FALSE)
    }
    
    #cor_mat <- ifelse(cor_mat_p > correlation_threshold, 0, 1)
    colnames(cor_mat) <- rownames(cor_mat) <- rownames(x)
    return(cor_mat)
}

#' @name bayes
#' @aliases bayes
#' @title Create an adjacency matrix based on constraint-based structure 
#' learning algorithm
#' @description  \code{bayes} infers an adjacency matrix using 
#' constraint-based structure learning algorithm \code{fast.iamb} from the 
#' \code{bnlearn} package. \code{bayes} extracts then the reported 
#' connections from running the \code{fast.iamb} function and assigns the
#' arcs of the discrete Bayesian connections to binary values. The 
#' adjacency matrix is returned by \code{bayes}. 
#' @usage bayes(x, ...)
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values 
#' @param ... parameters passed to \code{fast.iamb}
#' @details For use of the parameters used in the \code{fast.iamb} function, 
#' refer to ?bnlearn::fast.iamb.
#' @return matrix, matrix with edges inferred from constraint-based structure 
#' learning algorithm \code{fast.iamb}
#' @author Thomas Naake, \email{thomasnaake @googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' bayes(x)
#' @export
bayes <- function(x, algorithm = "tabu", R = 100, ...) {
    x_df <- data.frame(t(x))
    
    ## allow for compatibility of arguments 
    strength <- threeDotsCall(bnlearn::boot.strength, data = x_df, ...)
    
    ## create empty bs_mat to be filled with connections
    bs_mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    colnames(bs_mat) <- rownames(bs_mat) <- rownames(x)
    
    ## write to bs_mat
    for(i in seq_len(nrow(x_strength))) {
        tmp <- as.character(strength[i, ])
        names(tmp) <- names(strength[i, ])
        bs_mat[tmp["from"], tmp["to"] ] <- tmp["strength"]
    }
    mode(bs_mat) <- "numeric"
    
    return(bs_mat)
}

#' @name addToList
#' @aliases addToList
#' @title Add adjacency matrix to list
#' @description This helper function used in the function
#' \code{createStatisticalAdjacencyList} adds a adjacency matrix to a list of 
#' adjacency matrices.
#' @usage addToList(l, name, object)
#' @param l list of adjacency matrices
#' @param name character, name of newly created entry
#' @param object matrix containing the adjacency matrix to be added
#' @details Used internally in \code{createStatisticalAdjacencyList}
#' @return list containing the existing adjacency matrices and the added 
#' adjacency matrix
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' cor_pearson <- correlation(x, type = "pearson")
#' cor_spearman <- correlation(x, type = "spearman")
#' l <- list(pearson = cor_pearson)
#' MetNet:::addToList(l, "spearman", cor_spearman)
addToList <- function(l, name, object) {
    if (!is.list(l)) stop("l is not a list")
    if (!is.character(name)) stop("name is not a character")
    if (!is.matrix(object)) stop("object is not a matrix")
    new_index <- length(l) + 1
    l[[new_index]] <- object
    names(l)[new_index] <- name
    return(l)
}

#' @name createStatisticalAdjacencyList
#' @aliases createStatisticalAdjacencyList
#' @title Create a list of statistical adjacency matrices
#' @description The function infers adjacency matrix topologies from 
#' statistical methods and returns matrices of these networks in a list. The 
#' function includes functionality to caluclate adjacency matrices based on 
#' LASSO (L1 norm)-regression, random forests, context likelihood of 
#' relatedness (CLR), the algorithm for the reconstruction of accurate 
#' cellular networks (ARACNE), Pearson correlation (also partial and 
#' semipartial), Spearman correlation (also partial and semipartial) 
#' and Constraint-based structure learning (Bayes). The function returns a 
#' list of adjacency matrices that are defined by \code{model}. 
#' @usage createStatisticalAdjacencyList(x, model, ...)
#' @param x matrix that contains intensity values of features/metabolites (rows)
#' per sample (columns). 
#' @param model, character vector containing the methods that will be used 
#' ("lasso", "randomForest", "clr", "aracne", "pearson", "pearson_partial", 
#' "pearson_semipartial","spearman", "spearman_partial", "spearman_semipartial",
#' "bayes")
#' @param ... parameters passed to the functions  \code{lasso}, 
#' \code{randomForest}, \code{clr}, \code{aracne}, \code{correlation} and/or
#' \code{bayes}
#' @details \code{createStatisticalAdjacencyList} calls the function
#' \code{lasso}, \code{randomForest}, \code{clr}, \code{aracne}, 
#' \code{correlation} (for "pearson", "pearson_partial", "pearson_semipartial",
#' "spearman", "spearman_partial", "spearman_semipartial") and/or \code{bayes} 
#' as specified by \code{model}. It will create adjacency matrices using the 
#' specified methods and will return a list containing the unweighted
#' adjacency matrix (if \code{model} is of length 1) or append these 
#' unweighted adjacency matrices to a list (if \code{model} is of 
#' length > 1). Internally \code{x} will be z-scaled and the z-scaled object 
#' will be used in \code{lasso}, \code{clr} and/or \code{aracne}.
#' @return list containing the respective adjacency matrices specified by 
#' \code{model}
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' createStatisticalAdjacencyList(x, c("pearson", "spearman"))
#' @export
statistical <- function(x, model, ...) {
    
    ## check if model complies with the implemented model and return error 
    ## if not so
    if (!(all(model %in% c("lasso", "randomForest", "clr", "aracne", 
            "pearson", "pearson_partial", "pearson_semipartial", 
            "spearman", "spearman_partial", "spearman_semipartial", "bayes"))))
        stop("method not implemented in createModelList")
        
    ## check if x is numeric matrix and return error if not so
    if (mode(x) != "numeric") stop("x is not a numerical matrix")
    
    ## z-scale
    x_z <- apply(x, 1, function(x) {
        (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    })
    x_z <- t(x_z)
    
    l <- list()
    
    ## add entry for lasso if "lasso" is in model
    if ("lasso" %in% model) {
        lasso <- threeDotsCall("lasso", x_z, ...)
        l <- addToList(l, "lasso", lasso)
        print("lasso finished")
    }
    ## add entry for randomForest if "randomForest" is in model
    if ("randomForest" %in% model) {
        randomForest <- threeDotsCall("randomForest", x, ...)
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
        l <- addToList(l, "clr", clr)
        print("clr finished.")
    }
    ## add entry for aracne if "aracne" is in model
    if ("aracne" %in% model) {
        aracne <- threeDotsCall("aracne", mi = mi_x_z, ...)
        l <- addToList(l, "aracne", aracne)
        print("aracne finished.")
    }
    ## add entry for pearson if "pearson" is in model
    if ("pearson" %in% model) {
        pearson <- correlation(x, type = "pearson", ...)
        l <- addToList(l, "pearson", pearson)
        print("pearson finished.")
    }
    ## add entry for pearson_partial if "pearson_partial" is in model
    if ("pearson_partial" %in% model) {
        pearson_partial <- correlation(x, type = "pearson_partial", ...)
        l <- addToList(l, "pearson_partial", pearson_partial)
        print("pearson_partial finished.")
    }
    ## add entry for pearson_semipartial if "pearson_semipartial" is in model
    if ("pearson_semipartial" %in% model) {
        pearson_sp <- correlation(x, type = "pearson_semipartial", ...)
        l <- addToList(l, "pearson_semipartial", pearson_sp)
        print("pearson_semipartial finished.")
    }
    ## add entry for spearman if "spearman" is in model
    if ("spearman" %in% model) {
        spearman <- correlation(x, type = "spearman", ...)
        l <- addToList(l, "spearman", spearman)
        print("spearman finished.")
    }
    ## add entry for spearman_partial if "spearman_partial" is in model
    if ("spearman_partial" %in% model) {
        spearman_partial <- correlation(x, type = "spearman_partial", ...)
        l <- addToList(l, "spearman_partial", spearman_partial)
        print("spearman_partial finished.")
    }
    ## add entry for spearman_semipartial if "spearman_semipartial" is in model
    if ("spearman_semipartial" %in% model) {
        spearman_sp <- correlation(x, type = "spearman_semipartial", ...)
        l <- addToList(l, "spearman_semipartial", spearman_sp)
        print("spearman_semipartial finished.")
    }
    ## add entry for bayes if "bayes" is in model
    if ("bayes" %in% model) {
        bayes <- bayes(x, ...)
        l <- addToList(l, "bayes", bayes)
        print("bayes finished.")
    }
    return(l)
}

#' @name matrixToDataFrame
#' @aliases matrixToList
#' @title Write an adjacency matrix to a data.frame
#' @description 
#' @param mat
#' @details 
#' @return 
#' @author
#' @examples 
#' @export
getLinks <- function(mat, decreasing = TRUE) {
    
    ## replace 0 by NaN
    mat[which(mat == 0)] <- NaN
    df <- data.frame(row = c(row(mat)), col = c(col(mat)), confidence = c(mat))
    
    ## treat confidence values depending on decreasing parameter
    ## if TRUE, then the highest confidence value should get the first rank
    ## if TRUE, recalculate the confidence values that the values with highest
    ## support have low values
    if (decreasing) {
        conf <- max(df$confidence, na.rm = TRUE) - df$confidence
    } else {
        conf <- df$confidence - min(df$confidence, na.rm = TRUE)
    }
    
    ## calculate rank and add to data.frame
    df <- data.frame(df, rank = rank(conf, na.last = TRUE))
    
    ## return
    return(df)
}

#' 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' model <- c("pearson", "spearman")
#' args <- list("pearson" = 0.05, "spearman" = 0.05)

l <- statistical(x, model = model)

threshold(statistical, method, args) {
    
    l <- statistical
    ## args, either N for tops
    ## or a list of threshold
    if (!method %in% c("top1", "top2", "mean", "threshold"))
        stop("method not in 'top1', 'top2', 'mean', 'threshold'")
    
    ## check args
    if (method %in% c("top1", "top2", "mean")) {
        if (!all(model %in% names(args))) {
            stop()
        }
    }
    
    if (method %in% threshold) {
        if (! "n"  %in% names(args))
            stop()
    }
    
    if (method == "threshold") {
        ## iterate through the list and remove the links below or above the threshold
        ## write to list
        l <- lapply(seq_along(l), function(x) {
            
            ## find corresponding model in l 
            name_x <- names(l)[x]
            
            ## get corresponding threshold in args
            threshold_x <- args[[names(l)[x]]]
            
            ## get corresponding adjacency matrix in l
            l_x <- l[[name_x]]
            
            ## for pearson/spearman correlation methods (incl. partial and 
            ## semi-partial), low values correspond to higher confidence,
            ## only assign 1 to values that are below the threshold
            if (grepl(name_x, pattern = "pearson|spearman")) {
                ifelse(l_x < threshold_x, 1, 0)    
                ## for lasso, randomForest, clr, aracne and bayes higher values 
                ## corresond to higher confidence 
                ## only assign 1 to values that are above the threshold
            } else {
                ifelse(l_x > threshold_x, 1, 0)
            }
        })

        ## allow for compatibility of arguments 
        ## calculate consenses from the binary matrices
        cons <- threeDotsCall(sna::consensus, dat = l, ...)
        
        ##if (method == "central.graph") threshold <- 1
        ## threshold consensus that it is a binary matrix
        cons <- ifelse(cons >= threshold, 1, 0)
        
        rownames(cons) <- colnames(cons) <- colnames(l[[1]])
        
    } else { ## if method %in% c("top1", "top2", "mean")
        l_df <- lapply(seq_along(l), function(x) {
            
            ## find corresponding model in l
            name_x <- names(l)[x]
            
            ## get corresponding adjacency matrix in l
            l_x <- l[[name_x]]
            
            ## for pearson/spearman correlation methods (incl. partial and 
            ## semi-partial), low values correspond to higher confidence
            if (grepl(name_x, pattern = "pearson|spearman")) {
                getLinks(l_x, decreasing = FALSE)   
            ## for lasso, randomForest, clr, aracne and bayes higher values 
            ## corresond to higher confidence 
              
            } else {
                getLinks(l_x, decreasing = TRUE)
            }
        })
        
        names(l_df) <- names(l)
        
        ## bind together the ranks of the models, stored in l_df
        ranks <- lapply(l_df, function(x) x$rank)
        ranks <- do.call("cbind", ranks)
        colnames(ranks) <- names(l_df)
        
        ## calculate the consensus information, i.e. either get the first or 
        ## second top rank per row or calculate the average across rows
        ## depending on the method argument
        cons_val <- topKnet(ranks, method)
        
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
    
#' statistical()
#' l_df <- getLinks()
#' getRanks
getRanks <- function(l_df) {
    
}

#' 
topKnet <- function(ranks, method) {
    
    if (!is.matrix(ranks) && is.numeric(ranks)) {
        stop("ranks is not a numerical matrix")
    }
    if (! method %in% c("top1", "top2", "mean")) {
        stop("method neither 'top1', 'top2' or 'mean'")
    }
    ## calculate the consensus information, i.e. either get the first or 
    ## second top rank per row or calculate the average across rows
    ## depending on the method argument
    if (method == "top1") {
        ## get the lowest rank
        cons_val <- apply(ranks, 1, min)
    }
    if (method == "top2") {
        ## get the second lowest rank
        cons_val <- apply(ranks, 1, function(x) sort(x)[2])
    }
    if (method == "mean") {
        ## get the average of all ranks
        cons_val <- apply(ranks, 1, mean, na.rm = TRUE)
    }
    
    return(cons_val)
}





#' @name consensusAdjacency
#' @aliases consensusAdjacency
#' @title 
#' Create a consensus adjacency matrix of statistical adjacency matrices
#' 
#' @description 
#' The function takes a list of parameters (\code{l}) as input and
#' creates a consensus adjacency matrix from these adjacency matrices by 
#' calling the function \code{consensus} from the \code{sna} package. Depending 
#' on the chosen \code{method} in \code{consensus}, the threshold of the 
#' consensus adjacency matrix should be chosen accordingly to report a 
#' connection by different statistical methods. 
#' 
#' @param l list, each entry of the list contains an adjacency matrix
#' @param threshold numeric, when combining the adjacency matrices the 
#' \code{threshold} parameter defines if an edge is reported or not. For 
#' \code{method = "central.graph"} threshold is set to 1 by default. For other
#' values of method, the value should be carefully defined by the user. If 
#' threshold is set to \code{NULL} (default), it will be set to 1 internally. 
#' 
#' @param ... parameters passed to the function \code{consensus} in the
#' \code{sna} package 
#' 
#' @details \code{consensusAdjacency} is a wrapper function of the 
#' \code{consensus} function of the \code{sna} package. For use of the 
#' parameters used in the \code{consensus} function, refer to ?sna::consensus.
#' 
#' @return matrix, consensus matrix from adjacency matrices
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' stat_adj_l <- createStatisticalAdjacencyList(x, c("pearson", "spearman"))
#' consensusAdjacency(stat_adj_l)
#' 
#' @export
consensusAdjacency <- function(l, threshold = 1, ...) {
    
    if (!is.numeric(threshold)) stop("threshold is not numeric")
    if (!is.list(l)) stop("l is not a list")
    
    ## check compatibility of matrices
    ncol_1 <- ncol(l[[1]])
    nrow_1 <- nrow(l[[1]])
    rownames_1 <- rownames(l[[1]])
    colnames_1 <- colnames(l[[1]])
    
    if (!all(rownames_1 == colnames_1))
        stop("colnames and rownames are not identical")
    
    if (length(l) > 1) {
        for (i in 2:length(l)) {
            ## check if number of columns and rows are identical for 
            ## all matrices
            if (ncol(l[[i]]) != ncol_1)
                stop("ncol of matrices are not identical")
            if (nrow(l[[i]]) != nrow_1)
                stop("nrow of matrices are not identical")
            ## check if colnames and rownames are identical for all matrices
            if (!all(colnames(l[[i]]) == colnames_1))
                stop("matrices have different colnames")
            if (!all(rownames(l[[i]]) == rownames_1))
                stop("matrices have different rownames")
            if (!is.null(colnames_1) & is.null(colnames(l[[i]])))
                stop("matrices have different colnames")
            if (!is.null(rownames_1) & is.null(rownames(l[[i]])))
                stop("matrices have different rownames")
        }
    }
    ## end check compatibility
    
    
    ## allow for compatibility of arguments
    consensus_mat <- threeDotsCall(sna::consensus, dat = l, ...)
    ## was sna::consensus(dat = l, ...)
    
    ##if (method == "central.graph") threshold <- 1
    consensus_mat <- ifelse(consensus_mat >= threshold, 1, 0)
    
    rownames(consensus_mat) <- colnames(consensus_mat) <- colnames(l[[1]])
    return(consensus_mat)
}

#' @name threeDotsCall
#' @aliases threeDotsCall
#' @title Check if passed arguments match the function's formal arguments and
#' call the function with the checked arguments
#' @description The function \code{threeDotsCall} gets the formal arguments
#' of a function \code{fun} and checks if the passed arguments \code{...} 
#' matches the formal arguments. \code{threeDotsCall} will remove 
#' duplicated arguments. \code{threeDotsCall} will call the function 
#' \code{fun} with the filtered arguments and will return the result.
#' @usage threeDotsCall(fun, ...)
#' @param fun function to check for arguments and to call
#' @param ... arguments to be tested to be passed to fun
#' @details Used internally in \code{lasso}, \code{randomForest},
#' \code{correlation}, \code{bayes}, \code{consensusAdjacency}
#' @return Function call with passed arguments 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
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
    input <- input[!duplicated(names(input))]
    res <- do.call(fun, input)
    return(res)
}

#' @name createStatisticalAdjacency
#' @aliases createStatisticalAdjacency
#' @title Create statistical adjacency matrix
#' @description \code{createStatisticalAdjacency} creates a consensus 
#' adjacency matrix given the models to use. 
#' @usage createStatisticalAdjacency(x, model, threshold = 1, ...)
#' @param x matrix that contains intensity values of features/metabolites 
#' (rows) per sample (columns). 
#' @param model, character, vector containing the model that will be used 
#' ("lasso", "randomForest", "clr", "aracne", "pearson", "pearson_partial", 
#' "pearson_semipartial","spearman", "spearman_partial", 
#' "spearman_semipartial", "bayes")
#' @param threshold numeric, when combining the adjacency matrices the 
#' threshold parameter defines if an edge is reported or not. For 
#' \code{method = "central.graph"} threshold is set to 1 by default. For other
#' values of method, the value should be carefully defined by the user. If 
#' threshold is set to NULL (default), it will be set to 1 internally. 
#' @param ... parameters passed to the functions  \code{lasso}, 
#' \code{randomForest}, \code{clr}, \code{aracne}, \code{correlation},
#' \code{bayes} and/or \code{consensusAdjacency}
#' @return matrix, containing binary values if a connection is present or not
#' @details \code{createStatisticalAdjacency} is a wrapper function for the 
#' functions \code{createStatisticalAdjacencyList} and 
#' \code{consensusAdjacency}. See \code{?createStatisticalAdjacencyList} and 
#' \code{?consensusAdjacency} for further details. The function 
#' \code{createStatisticalAdjacencyList}
#' includes functionality to caluclate adjacency matrices based on 
#' LASSO (L1 norm)-regression, random forests, context likelihood of 
#' relatedness (CLR), the algorithm for the reconstruction of accurate 
#' cellular networks (ARACNE), Pearson correlation (also partial and 
#' semipartial), Spearman correlation (also partial and semipartial)
#' and Constraint-based structure learning (Bayes). 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' createStatisticalAdjacency(x, c("pearson", "spearman"))
#' @export
createStatisticalAdjacency <- function(x, model, threshold = 1, ...) {
    ## first use function createStatisticalAdjacency_list
    l <- createStatisticalAdjacencyList(x = x, model = model, ...)
    ## combine statistical adjacency matrices by the function
    ## consensusAdjacency
    consensus_mat <- consensusAdjacency(l = l, threshold = threshold, ...)
    return(consensus_mat)
}


