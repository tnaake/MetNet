#' @importFrom stabs stabsel stabsel.matrix glmnet.lasso
#' @importFrom rfPermute rfPermute rfPermute.formula rp.importance
#' @importFrom mpmi cmi
#' @importFrom parmigene clr aracne.a
#' @importFrom psych corr.test
#' @importFrom bnlearn fast.iamb arcs
#' @importFrom sna consensus
#' @importFrom parallel mclapply
#' @importFrom stats formula p.adjust sd
#' @importFrom methods formalArgs
#' @importFrom ppcor pcor spcor


#' @name lasso
#' @aliases lasso
#' @title Create a statistical network based on LASSO
#' @description  \code{lasso} infers a statistical network using 
#' LASSO using the \code{stabsel} function from the 
#' \code{stabs} package. \code{lasso} extracts the  
#' predictors from the function \code{stabsel} and the presence/absence of this 
#' connection to a matrix that is returned. 
#' @usage lasso(x, parallel = FALSE, ...)
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values 
#' @param parallel logical, should computation be parallelized? If 
#' \code{parallel = TRUE} the \code{mclapply} will be applied if 
#' \code{parallel = FALSE} the \code{lapply} function will be applied. 
#' @param ... parameters passed to \code{corr.test} and \code{mclapply} (if 
#' \code{parallel = TRUE})
#' @details For use of the parameters used in the \code{corr.test} function, 
#' refer to ?psych::corr.test.
#' @return matrix, matrix with edges inferred from LASSO algorithm 
#' \code{stabsel}
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' x_z <- t(apply(x, 1, function(y) (y - mean(y)) / sd(y)))
#' \dontrun{lasso(x_z, PFER = 0.75, cutoff = 0.95)}
#' @export
lasso <- function(x, parallel = FALSE, ...) {
    ## x should be z-scaled
    if (parallel) {
        l1 <- mclapply(seq_len(nrow(x)), function(i) {
            x_l1 <- t(x[-i, ]); y_l1 <- x[i, ]
            ## lasso: alpha = 1
            ## allow for compatibility of arguments 
            l1 <- threeDots_call("stabsel.matrix", x = as.matrix(x_l1), 
                    y = y_l1, fitfun = glmnet.lasso, 
                    args.fitfun = list("alpha" = 1), ...)
            return(l1$selected)
        }, ...)    
    } else {
        l1 <- lapply(seq_len(nrow(x)), function(i) {
            x_l1 <- t(x[-i, ]); y_l1 <- x[i, ]
            ## lasso: alpha = 1
            ## allow for compatibility of arguments 
            l1 <- threeDots_call("stabsel.matrix", x = as.matrix(x_l1),
                    y = y_l1, fitfun = glmnet.lasso, 
                    args.fitfun = list("alpha" = 1), ...)
            return(l1$selected)
        })   
    }
    
    l1_mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    colnames(l1_mat) <- rownames(l1_mat) <- rownames(x)
    for (i in seq_len((l1))) {l1_mat[names(l1[[i]]), i] <- 1} 
    ## ; l1_mat[i, l1[[i]]] <- 1}
    return(l1_mat)
}

#' @name randomForest
#' @aliases randomForest
#' @title Create a statistical network based on random forest
#' @description  \code{randomForest} infers a statistical network using 
#' random forest using the \code{rfPermute} function from the 
#' \code{rfPermute} package. \code{randomForest} extracts the p-values  
#' by the function \code{rp.importance} and the presence/absence based on the 
#' significance value (\eqn{\alpha \leq 0.05}) of this 
#' connection to a matrix that is returned. 
#' @usage randomForest(x, parallel = FALSE, randomForest_adjust = "none", ...)
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values 
#' @param parallel logical, should computation be parallelized? If 
#' \code{parallel = TRUE} the \code{mclapply} will be applied if 
#' \code{parallel = FALSE} the \code{lapply} function will be applied. 
#' @param randomForest_adjust character, correction method for p-values from
#' \code{rp.importance}, \code{randomForest_adjust} will be passed to the 
#' \code{p.adjust} function and should be one of "holm", "hochberg", "hommel", 
#' "bonferroni", "BH", "BY", "fdr", "none"
#' @param ... parameters passed to \code{corr.test} and \code{mclapply} (if 
#' \code{parallel = TRUE})
#' @details For use of the parameters used in the \code{rfPermute} function, 
#' refer to ?rfPermute::rfPermute.
#' @return matrix, matrix with edges inferred from random forest algorithm 
#' \code{rfPermute} and \code{rp.importance}
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' \dontrun{randomForest(x)}
#' @export
randomForest <- function(x, parallel = FALSE, randomForest_adjust = "none", 
                                                                        ...) {
    df_x <- data.matrix(t(x))
    
    if (parallel) {
        rf <- parallel::mclapply(seq_len(nrow(x)), function(i) {
            formula_rf <- paste(rownames(x)[i], "~", ".")    
            ## allow for compatibility of arguments 
            rf <- threeDots_call(rfPermute::rfPermute.formula, 
                    formula = stats::formula(formula_rf), data = df_x, ...)
            rf_p <- rp.importance(rf)[,"IncNodePurity.pval"]
            return(rf_p)
        }, mc.cores = 4)
    } else {
        rf <- lapply(seq_len(nrow(x)), function(i) {
            formula_rf <- paste(rownames(x)[i], "~", ".")    
            ## allow for compatibility of arguments 
            rf <- threeDots_call(rfPermute::rfPermute.formula, 
                    formula = stats::formula(formula_rf), data = df_x, ...) 
            rf_p <- rp.importance(rf)[,"IncNodePurity.pval"]
            return(rf_p)
        })
    }
    rf_mat <- matrix(1, nrow = nrow(x), ncol = nrow(x))    
    colnames(rf_mat) <- rownames(rf_mat) <- rownames(x)
    
    for (i in seq_len(length(rf))) {
        rf_mat[names(rf[[i]]), rownames(x)[i]] <- rf[[i]]
    }
    rf_mat <- stats::p.adjust(rf_mat, method = randomForest_adjust)     
    rf_mat <- matrix(rf_mat, ncol = nrow(x), nrow = nrow(x), byrow = FALSE)
    rf_mat <- ifelse(rf_mat > 0.05, 0, 1)
    colnames(rf_mat) <- rownames(rf_mat) <- rownames(x)
    
    return(rf_mat)
}

#' @name clr
#' @aliases clr
#' @title Create a statistical network based on context likelihood or 
#' relatedness network
#' @description  \code{clr} infers a statistical network using 
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
clr <- function(mi, clr_threshold = 0) {
    if (!is.numeric(clr_threshold)) stop("clr_threshold is not numeric")
    clr_mat <- parmigene::clr(mi)
    clr_mat <- ifelse(clr_mat > clr_threshold, 1, 0)
    colnames(clr_mat) <- rownames(clr_mat) <- rownames(mi)
    return(clr_mat)
}

#' @name aracne
#' @aliases aracne
#' @title Create a statistical network based on algorithm for the reconstruction
#' of accurate cellular networks 
#' @description  \code{.information} infers a statistical network using 
#' the algortihm for the reconstruction of accurate cellular networks 
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
#' (aracne$_{i,j}$ <= threshold) there is no statistical connection reported. 
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
aracne <- function(mi, eps = 0.05, aracne_threshold = 0) {
    if (!is.numeric(aracne_threshold)) stop("aracne_threshold is not numeric")
    aracne_mat <- parmigene::aracne.a(mi, eps = eps)  
    aracne_mat <- ifelse(aracne_mat > aracne_threshold, 1, 0)
    colnames(aracne_mat) <- rownames(aracne_mat) <- rownames(mi)
    return(aracne_mat)
}

#' @name correlation 
#' @aliases correlation
#' @title Create a statistical network based on correlation 
#' @description  \code{correlation} infers a statistical network using 
#' correlation using the \code{corr.test} function (from the 
#' \code{psych} package), \code{pcor} (from \code{ppcor}) or 
#' \code{spcor} (from \code{ppcor}). \code{correlation} extracts the reported 
#' p-values from the function \code{corr.test}, \code{pcor} or \code{spcor} 
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
#' \code{type} will be passed to argument \code{method} in \code{corr.test} 
#' (in the case of "pearson" or "spearman") or to \code{method} in \code{pcor} 
#' ("pearson" and "spearman" for "pearson_partial" and "spearman_partial", 
#' respectively) or to \code{method} in \code{spcor} ("pearson" or "spearman"
#' for "pearson_semipartial" and "spearman_semipartial", respectively)
#' @param correlation_adjust character 
#' @param correlation_threshold numeric, significance level \eqn{\alpha}
#' (default: 0.05), if the (adjusted) p-values exceed this value, there 
#' is no statistical connection between features 
#' @param ... parameters passed to \code{corr.test} (argument \code{adjust} 
#' will be ignored)
#' @details If "pearson" or "spearman" is used as a \code{method} the function 
#' \code{corr.test} from \code{psych} will be employed. 
#' If "pearson_partial" or "spearman_partial" is used as a \code{method} the 
#' function \code{pcor} from \code{spcor} will be employed. 
#' If "pearson_semipartial" or "spearman_semipartial" is used as a \code{method}
#' the function \code{spcor} from \code{spcor} will be employed. 
#' For use of the parameters used in the \code{corr.test} function, 
#' refer to ?psych::corr.test. 
#' @return matrix, matrix with edges inferred from correlation algorithm 
#' \code{corr.test}, \code{pcor} or \code{spcor} (depending on the chosen 
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
        cor_mat_p <- threeDots_call(psych::corr.test, x = t(x), 
                                    adjust = adjust, method = type, ...)$p    
    }
    if (type %in% c("pearson_partial", "spearman_partial")) {
        if (type == "pearson_partial") method <- "pearson"
        if (type == "spearman_partial") method <- "spearman"
        cor_mat_p <- ppcor::pcor(t(x), method = method)$p.value
        cor_mat_p <- stats::p.adjust(cor_mat_p, 
                            method = correlation_adjust)     
        cor_mat_p <- matrix(cor_mat_p, ncol = nrow(x), 
                            nrow = nrow(x), byrow = FALSE)
    }
    if (type %in% c("pearson_semipartial", "spearman_semipartial")) {
        if (type == "pearson_semipartial") method <- "pearson"
        if (type == "spearman_semipartial") method <- "spearman"
        cor_mat_p <- ppcor::spcor(t(x), method = method)$p.value
        cor_mat_p <- stats::p.adjust(cor_mat_p, 
                            method = correlation_adjust)     
        cor_mat_p <- matrix(cor_mat_p, ncol = nrow(x), 
                            nrow = nrow(x), byrow = FALSE)
    }
    
    ## was cor_mat_p <- corr.test(t(x), adjust = adjust, ...)$p
    cor_mat <- ifelse(cor_mat_p > correlation_threshold, 0, 1)
    colnames(cor_mat) <- rownames(cor_mat) <- rownames(x)
    return(cor_mat)
}

#' @name bayes
#' @aliases bayes
#' @title Create a statistical network based on constraint-based structure 
#' learning algorithm
#' @description  \code{bayes} infers a statistical network using 
#' constraint-based structure learning algorithm \code{fast.iamb} from the 
#' \code{bnlearn} package. \code{bayes} extracts then the reported 
#' connections from running the \code{fast.iamb} function and assigns the
#' arcs of the discrete Bayesian network to binary values into a matrix that is 
#' returned by \code{bayes}. 
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
bayes <- function(x, ...) {
    x_df <- data.frame(t(x))
    ## allow for compatibility of arguments 
    x_fast.iamb <- threeDots_call(bnlearn::fast.iamb, x = x_df, ...)
    ## was x_fast.iamb <- fast.iamb(x_df, ...) 
    bs_mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    colnames(bs_mat) <- rownames(bs_mat) <- rownames(x)
    arcs_fast.iamb <- bnlearn::arcs(x_fast.iamb)
    for(i in seq_len(nrow(arcs_fast.iamb))) {
        bs_mat[arcs_fast.iamb[i, "from"], arcs_fast.iamb[i, "to"] ] <- 1} 
    bs_mat <- as.matrix(bs_mat)
    return(bs_mat)
}

#' @name add_to_list
#' @aliases add_to_list
#' @title Add network to list
#' @description This helper function used in the function
#' \code{create_statistical_networks_list} adds a network to a list of networks.
#' @usage add_to_list(l, name, object)
#' @param l list of networks
#' @param name character, name of newly created entry
#' @param object matrix containing the network to be added
#' @details Used internally in \code{create_statistical_networks_list}
#' @return list containing the existing networks and the added network 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' cor_pearson <- correlation(x, type = "pearson")
#' cor_spearman <- correlation(x, type = "spearman")
#' l <- list(pearson = cor_pearson)
#' add_to_list(l, "spearman", cor_spearman)
#' @export
add_to_list <- function(l, name, object) {
    if (!is.list(l)) stop("l is not a list")
    if (!is.character(name)) stop("name is not a character")
    if (!is.matrix(object)) stop("object is not a matrix")
    new_index <- length(l) + 1
    l[[new_index]] <- object
    names(l)[new_index] <- name
    return(l)
}



#' @name create_statistical_networks_list
#' @aliases create_statistical_networks_list
#' @title Create a list of statistical networks
#' @description The function infers network topologies from statistical 
#' methods and returns matrices of these networks in a list. The function
#' includes functionality to caluclate statistical networks based on 
#' LASSO (L1 norm)-regression, random forests, context likelihood of 
#' relatedness (CLR), the algorithm for the reconstruction of accurate 
#' cellular networks (ARACNE), Pearson correlation (also partial and 
#' semipartial), Spearman correlation (also partial and semipartial) 
#' and Constraint-based structure learning (Bayes). 
#' @usage create_statistical_networks_list(x, model, ...)
#' @param x matrix that contains intensity values of features/metabolites (rows)
#' per sample (columns). 
#' @param model, character vector containing the methods that will be used 
#' ("lasso", "randomForest", "clr", "aracne", "pearson", "pearson_partial", 
#' "pearson_semipartial","spearman", "spearman_partial", "spearman_semipartial",
#' "bayes")
#' @param ... parameters passed to the functions  \code{lasso}, 
#' \code{randomForest}, \code{clr}, \code{aracne}, \code{correlation} and/or
#' \code{bayes}
#' @details \code{create_statistical_networks_list} calls the function
#' \code{lasso}, \code{randomForest}, \code{clr}, \code{aracne}, 
#' \code{correlation} (for "pearson", "pearson_partial", "pearson_semipartial",
#' "spearman", "spearman_partial", "spearman_semipartial") and/or \code{bayes} 
#' as specified by \code{model}. It will create network(s) using the specified 
#' methods and will return a list containing the binary network (if model is 
#' of length 1) or append these binary networks to a list (if model is of 
#' length > 1). Internally x will be z-scaled and the z-scaled object will 
#' be used in \code{lasso}, \code{clr} and/or \code{aracne}.
#' @return list containing the respective statistical networks specified by 
#' \code{model}
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' create_statistical_networks_list(x, c("pearson", "spearman"))
#' @export
create_statistical_networks_list <- function(x, model, ...) {
    
    ## check if model complies with the implemented model and return error 
    ## if not so
    if (!(all(model %in% c("lasso", "randomForest", "clr", "aracne", 
            "pearson", "pearson_partial", "pearson_semipartial", 
            "spearman", "spearman_partial", "spearman_semipartial", "bayes"))))
        stop("method not implemented in create_statistical_networks_list")
        
    ## check if x is numeric matrix and return error if not so
    if (mode(x) != "numeric") stop("x is not a numerical matrix")
    
    x_z <- apply(x, 1, function(x) (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))
    x_z <- t(x_z)
    
    l <- list()
    
    ## add entry for lasso if "lasso" is in model
    if ("lasso" %in% model) {
        lasso <- lasso(x_z, ...)
        l <- add_to_list(l, "lasso", lasso)
    }
    ## add entry for randomForest if "randomForest" is in model
    if ("randomForest" %in% model) {
        randomForest <- randomForest(x, ...)
        l <- add_to_list(l, "randomForest", randomForest)
    }
    
    ## calculate mutual information if "clr" or "aracne" is in model
    if (any(c("clr", "aracne") %in% model)) {
        mi_x_z <- mpmi::cmi(t(x_z))$bcmi
        rownames(mi_x_z) <- colnames(mi_x_z) <- rownames(x)
    }
    
    ## add entry for clr if "clr" is in model
    if ("clr" %in% model) {
        clr <- threeDots_call("clr", mi = mi_x_z, ...)
        l <- add_to_list(l, "clr", clr)
    }
    ## add entry for aracne if "aracne" is in model
    if ("aracne" %in% model) {
        aracne <- threeDots_call("aracne", mi = mi_x_z, ...)
        l <- add_to_list(l, "aracne", aracne)
    }
    ## add entry for pearson if "pearson" is in model
    if ("pearson" %in% model) {
        pearson <- correlation(x, type = "pearson", ...)
        l <- add_to_list(l, "pearson", pearson)
    }
    ## add entry for pearson_partial if "pearson_partial" is in model
    if ("pearson_partial" %in% model) {
        pearson <- correlation(x, type = "pearson_partial", ...)
        l <- add_to_list(l, "pearson_partial", pearson)
    }
    ## add entry for pearson_semipartial if "pearson_semipartial" is in model
    if ("pearson_semipartial" %in% model) {
        pearson <- correlation(x, type = "pearson_semipartial", ...)
        l <- add_to_list(l, "pearson_semipartial", pearson)
    }
    ## add entry for spearman if "spearman" is in model
    if ("spearman" %in% model) {
        spearman <- correlation(x, type = "spearman", ...)
        l <- add_to_list(l, "spearman", spearman)
    }
    ## add entry for spearman_partial if "spearman_partial" is in model
    if ("spearman_partial" %in% model) {
        pearson <- correlation(x, type = "spearman_partial", ...)
        l <- add_to_list(l, "spearman_partial", pearson)
    }
    ## add entry for spearman_semipartial if "spearman_semipartial" is in model
    if ("spearman_semipartial" %in% model) {
        pearson <- correlation(x, type = "spearman_semipartial", ...)
        l <- add_to_list(l, "spearman_semipartial", pearson)
    }
    ## add entry for bayes if "bayes" is in model
    if ("bayes" %in% model) {
        bayes <- bayes(x, ...)
        l <- add_to_list(l, "bayes", bayes)
    }
    return(l)
}

#' @name consensus_network
#' @aliases consensus_network
#' @title Create a consensus network of statistical networks
#' @description The function takes a list of parameters (\code{l}) as input and
#' creates a consensus network from these networks by calling the function 
#' \code{consensus} from the \code{sna} package. Depending on the chosen 
#' \code{method } in \code{consensus}, the threshold of the consensus network
#' should be chosen accordingly to report a connection by different 
#' statistical methods. 
#' @usage consensus_network(l, threshold = 1, ...)
#' @param l list, each entry of the list contains a network
#' @param threshold numeric, when combining the networks the threshold 
#' parameter defines if an edge is reported or not. For 
#' \code{method = "central.graph"} threshold is set to 1 by default. For other
#' values of method, the value should be carefully defined by the user. If 
#' threshold is set to NULL (default), it will be set to 1 internally. 
#' @param ... parameters passed to the function \code{consensus} in the
#' \code{sna} package 
#' @details \code{consensus_network} is a wrapper function of the 
#' \code{consensus} function of the \code{sna} package. For use of the 
#' parameters used in the \code{consensus} function, refer to ?sna::consensus.
#' @return matrix, consensus matrix from statistical networks
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- x_test[, 3:dim(x_test)[2]]
#' x <- as.matrix(x)
#' stat_net_l <- create_statistical_networks_list(x, c("pearson", "spearman"))
#' consensus_network(stat_net_l)
#' @export
consensus_network <- function(l, threshold = 1, ...) {
    
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
    consensus_mat <- threeDots_call(sna::consensus, dat = l, ...)
    ## was sna::consensus(dat = l, ...)
    
    ##if (method == "central.graph") threshold <- 1
    consensus_mat <- ifelse(consensus_mat >= threshold, 1, 0)
    
    rownames(consensus_mat) <- colnames(consensus_mat) <- colnames(l[[1]])
    return(consensus_mat)
}

#' @name threeDots_call
#' @aliases threeDots_call
#' @title Check if passed arguments match the function's formal arguments and
#' call the function with the checked arguments
#' @description The function \code{threeDots_call} gets the formal arguments
#' of a function \code{fun} and checks if the passed arguments \code{...} 
#' matches the formal arguments. \code{threeDots_call} will remove 
#' duplicated arguments. \code{threeDots_call} will call the function 
#' \code{fun} with the filtered arguments and will return the result. 
#' @usage threeDots_call(fun, ...)
#' @param fun function to check for arguments and to call
#' @param ... arguments to be tested to be passed to fun
#' @details Used internally in \code{lasso}, \code{randomForest}, 
#' \code{correlation}, \code{bayes}, \code{consensus_network}
#' @return Function call with passed arguments 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' threeDots_call(stats::sd, x = 1:10, y = 1:10)
#' ## in contrast to the above example, the following example will result in an 
#' ## error
#' \dontrun{stats::sd(x = 1:10, y = 1:10)}
#' @export
threeDots_call <- function(fun, ...) {
    formal_args <- formalArgs(fun)
    args <- list(...)
    if (any(duplicated(names(args)))) stop("duplicated args in ...")
    
    input <- args[names(args) %in% formal_args]
    input <- input[!duplicated(names(input))]
    res <- do.call(fun, input)
    return(res)
}

#' @name create_statistical_network
#' @aliases create_statistical_network
#' @title Create statistical network
#' @description \code{create_statistical_network} creates a consensus network
#' given the models to use. 
#' @usage create_statistical_network(x, model, threshold = 1, ...)
#' @param x matrix that contains intensity values of features/metabolites (rows)
#' per sample (columns). 
#' @param model, character, vector containing the model that will be used 
#' ("lasso", "randomForest", "clr", "aracne", "pearson", "pearson_partial", 
#' "pearson_semipartial","spearman", "spearman_partial", "spearman_semipartial",
#' "bayes")
#' @param threshold numeric, when combining the networks the threshold 
#' parameter defines if an edge is reported or not. For 
#' \code{method = "central.graph"} threshold is set to 1 by default. For other
#' values of method, the value should be carefully defined by the user. If 
#' threshold is set to NULL (default), it will be set to 1 internally. 
#' @param ... parameters passed to the functions  \code{lasso}, 
#' \code{randomForest}, \code{clr}, \code{aracne}, \code{correlation},
#' \code{bayes} and/or \code{consensus_network}
#' @return matrix, containing binary values if a connection is present or not
#' @details \code{create_statistical_network} is a wrapper function for the 
#' functions \code{create_statistical_networks_list} and 
#' \code{consensus_network}. See \code{?create_statistical_networks_list} and 
#' \code{?consensus_network} for further details. The function 
#' \code{create_statistical_networks_list}
#' includes functionality to caluclate statistical networks based on 
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
#' create_statistical_network(x, c("pearson", "spearman"))
#' @export
create_statistical_network <- function(x, model, threshold = 1, ...) {
    ## first use function create_statistical_network_list
    l <- create_statistical_networks_list(x = x, model = model, ...)
    ## combine statistical networks by the function consensus_network
    consensus_mat <- consensus_network(l = l, threshold = threshold, ...)
    return(consensus_mat)
}


