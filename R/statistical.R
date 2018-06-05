#' @import glmnet
#' @import stabs
#' @import randomForest
#' @import rfPermute
#' @import mpmi
#' @import parmigene
#' @import psych
#' @import bnlearn


#' @name .lasso
#' @title Create a statistical network based on LASSO
#' @description  \code{.lasso} infers a statistical network using 
#' LASSO using the \code{stabsel} function from the 
#' \code{stabs} package. \code{.lasso} extracts the  
#' predictors from the function \code{stabsel} and the presence/absence of this 
#' connection to a matrix that is returned. 
#' @usage .lasso(x, parallel = FALSE, ...)
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
#' @examples .lasso(x, parallel = FALSE, ...)
.lasso <- function(x, parallel = FALSE, ...) {
    ## x should be z-scaled
    if (parallel) {
        l1 <- mclapply(1:dim(x)[1], function(i) {
            x_l1 <- t(x[-i, ]); y_l1 <- x[i, ]
            ## lasso: alpha = 1
            ## allow for compatibility of arguments 
            l1 <- .threeDots_call("stabsel.matrix", x = as.matrix(x_l1), y = y_l1, 
                                  fitfun = glmnet.lasso, 
                                  args.fitfun = list("alpha" = 1), ...)
            ## was l1 <- stabsel(as.matrix(x_l1), y_l1, fitfun = glmnet.lasso, 
            ##          args.fitfun = list("alpha" = 1), ...)
            return(l1$selected)
        }, ...)    
    } else {
        l1 <- lapply(1:dim(x)[1], function(i) {
            x_l1 <- t(x[-i, ]); y_l1 <- x[i, ]
            ## lasso: alpha = 1
            ## allow for compatibility of arguments 
            l1 <- .threeDots_call("stabsel.matrix", x = as.matrix(x_l1), y = y_l1, 
                                  fitfun = glmnet.lasso, 
                                  args.fitfun = list("alpha" = 1), ...)
            ## was l1 <- stabsel(as.matrix(x_l1), y_l1, fitfun = glmnet.lasso, 
            ##          args.fitfun = list("alpha" = 1), ...)
            return(l1$selected)
        })   
    }
    
    l1_mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    colnames(l1_mat) <- rownames(l1_mat) <- rownames(x)
    for (i in 1:length(l1)) {l1_mat[names(l1[[i]]), i] <- 1} ## ; l1_mat[i, l1[[i]]] <- 1}
    
    return(l1_mat)
}

#' @name .randomForest
#' @title Create a statistical network based on random forest
#' @description  \code{.randomForest} infers a statistical network using 
#' random forest using the \code{rfPermute} function from the 
#' \code{rfPermute} package. \code{.randomForest} extracts the p-values  
#' by the function \code{rp.importance} and the presence/absence based on the 
#' significance value ($\alpha$ <= 0.05) of this 
#' connection to a matrix that is returned. 
#' @usage .randomForest(x, parallel = FALSE, randomforest_adjust = NULL, ...)
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values 
#' @param parallel logical, should computation be parallelized? If 
#' \code{parallel = TRUE} the \code{mclapply} will be applied if 
#' \code{parallel = FALSE} the \code{lapply} function will be applied. 
#' @param ... parameters passed to \code{corr.test} and \code{mclapply} (if 
#' \code{parallel = TRUE})
#' @details For use of the parameters used in the \code{rfPermute} function, 
#' refer to ?rfPermute::rfPermute.
#' @return matrix, matrix with edges inferred from random forest algorithm 
#' \code{rfPermute} and \code{rp.importance}
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples .randomForest(x, parallel = FALSE, randomForest_adjust = "none", ...)
.randomForest <- function(x, parallel = FALSE, randomforest_adjust = "none", ...) {
    df_x <- data.matrix(t(x))
    
    if (parallel) {
        rf <- mclapply(1:dim(x)[1], function(i) {
            formula_rf <- paste(rownames(x)[i], "~", ".")    
            ## allow for compatibility of arguments 
            rf <- .threeDots_call("rfPermute.formula", formula = formula(formula_rf), data = df_x, ...)
            ## was rf <- rfPermute(formula = formula(formula_rf), data = df_x, ...) 
            rf_p <- rp.importance(rf)[,"IncNodePurity.pval"]
            return(rf_p)
        }, mc.cores = 4)
    } else {
        rf <- lapply(1:dim(x)[1], function(i) {
            formula_rf <- paste(rownames(x)[i], "~", ".")    
            ## allow for compatibility of arguments 
            rf <- .threeDots_call("rfPermute.formula", formula = formula(formula_rf), data = df_x, ...)
            ## was rf <- rfPermute(formula = formula(formula_rf), data = df_x, ...)   
            rf_p <- rp.importance(rf)[,"IncNodePurity.pval"]
            return(rf_p)
        })
    }
    rf_mat <- matrix(1, nrow = nrow(x), ncol = nrow(x))    
    colnames(rf_mat) <- rownames(rf_mat) <- rownames(x)
    
    for (i in 1:length(rf)) {rf_mat[names(rf[[i]]), rownames(x)[i]] <- rf[[i]]}
    rf_mat <- p.adjust(rf_mat, method = randomforest_adjust)     
    rf_mat <- matrix(rf_mat, ncol = nrow(x), nrow = nrow(x), byrow = FALSE)
    rf_mat <- ifelse(rf_mat > 0.05, 0, 1)
    colnames(rf_mat) <- rownames(rf_mat) <- rownames(x)
    
    return(rf_mat)
}

#' @name .clr
#' @title Create a statistical network based on context likelihood or 
#' relatedness network
#' @description  \code{.clr} infers a statistical network using 
#' context likelihood/relatedness network using the \code{clr} function from the 
#' \code{parmigene} package. The presence/absence is based on if the 
#' returned value exceeds a user-defined threshold value. \code{.clr} will 
#' return the adjacency matrix containing the presence/absence value.
#' @usage .clr(mi, threshold_clr = 0)
#' @param mi matrix, where columns and the rows are features 
#' (metabolites), cell entries are mutual information values between the 
#' features. As input, the mutual information (e.g. raw MI estimates or 
#' Jackknife bias corrected MI estimates) from the \code{cmi} function of the
#' \code{mpmi} package can be used.
#' @param threshold_clr numeric, if the clr value exceeds the threshold 
#' (clr$_{i,j}$ > threshold, where clr$_{i, j}$ is the clr value of the ith row 
#' feature and of the jth column feature), the connection is defined as present, if 
#' the clr value is lower than the threshold value (clr$_{i,j}$ <= threshold)
#' there is no statistical connection reported. 
#' @details For more details on the \code{clr} function, 
#' refer to ?parmigene::clr.
#' @return matrix, matrix with edges inferred from Context Likelihood or 
#' Relatedness Network algorithm \code{clr}
#' @author Thomas Naake, \email{thomasnaake @googlemail.com}
#' @examples .clr(mi_x_z, threshold_clr = 0)
.clr <- function(mi, threshold_clr = 0) {
    if (!is.numeric(threshold_clr)) stop("threshold_clr is not numeric")
    clr_mat <- clr(mi)
    clr_mat <- ifelse(clr_mat > threshold_clr, 1, 0)
    colnames(clr_mat) <- rownames(clr_mat) <- rownames(mi)
    return(clr_mat)
}

#' @name .aracne
#' @title Create a statistical network based on algorithm for the reconstruction
#' of accurate cellular networks 
#' @description  \code{.information} infers a statistical network using 
#' the algortihm for the reconstruction of accurate cellular networks 
#' using the \code{aracne.a} function from the 
#' \code{parmigene} package. The presence/absence is based on if the 
#' returned value exceeds a user-defined threshold value. \code{.aracne} will 
#' return the adjacency matrix containing the presence/absence value.
#' @usage .aracne(mi = mi_x_z, eps = 0.05, threshold_aracne = 0)
#' @param mi matrix, where columns and the rows are features 
#' (metabolites), cell entries are mutual information values between the 
#' features. As input, the mutual information (e.g. raw MI estimates or 
#' Jackknife bias corrected MI estimates) from the \code{cmi} function of the
#' \code{mpmi} package can be used.
#' @param threshold_aracne numeric, if the aracne value exceeds the threshold 
#' (aracne$_{i,j}$ > threshold, where aracne$_{i, j}$ is the aracne value of 
#' the ith row feature and of the jth column feature), the connection is 
#' defined as present, if the aracne value is lower than the threshold value 
#' (aracne$_{i,j}$ <= threshold) there is no statistical connection reported. 
#' @details For more details on the \code{aracne.a} function, 
#' refer to ?parmigene::aracne.a.
#' @return matrix, matrix with edges inferred from Reconstruction of accurate
#' cellular networks algorithm \code{aracne}
#' @author Thomas Naake, \email{thomasnaake @googlemail.com}
#' @examples .aracne(mi_x_z, eps = 0.05, threshold_aracne = 0)
.aracne <- function(mi, eps = 0.05, threshold_aracne = 0) {
    if (!is.numeric(threshold_aracne)) stop("threshold_aracne is not numeric")
    aracne_mat <- aracne.a(mi, eps = eps)  
    aracne_mat <- ifelse(aracne_mat > threshold_aracne, 1, 0)
    colnames(aracne_mat) <- rownames(aracne_mat) <- rownames(mi)
    return(aracne_mat)
}

#' @name .correlation 
#' @title Create a statistical network based on correlation 
#' @description  \code{.correlation} infers a statistical network using 
#' correlation using the \code{corr.test} function from the 
#' \code{psych} package. \code{.correlation} extracts the reported 
#' p-values from the function \code{corr.test} that can be adjusted for 
#' multiple testing (\code{adjust_correlation} parameter). Significant 
#' correlations (p <= 0.05) will be assigned as 
#' connections from running the \code{fast.iamb} function and assigns the 
#' arcs of the discrete Bayesian network to binary values into a matrix that is 
#' returned by \code{.bayes}. In the \code{create_statistical_networks} 
#' @usage .correlation(x, adjust_correlation = "none", threshold_correlation = 0.05, ...) 
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values 
#' @param adjust_correlation character 
#' @param threshold_correlation numeric, significance level $\alpha$ 
#' (default: 0.05), if the (adjusted) p-values exceed this value, there 
#' is no statistical connection between features 
#' @param ... parameters passed to \code{corr.test} (argument \code{adjust} will be ignored) 
#' @details For use of the parameters used in the \code{corr.test} function, 
#' refer to ?psych::corr.test. 
#' @return matrix, matrix with edges inferred from correlation algorithm 
#' \code{corr.test}
#' @author Thomas Naake, \email{thomasnaake @googlemail.com}
#' @examples .correlation(x, adjust_correlation = "none", threshold_correlation = 0.05, ...)
.correlation <- function(x, adjust_correlation = "none", type = "pearson", 
                         threshold_correlation = 0.05, ...) {
    if (!is.numeric(threshold_correlation)) 
        stop("threshold_correlation is not numeric")
    ## get character vector for p-value adjustment
    adjust <- adjust_correlation
    ## allow for compatibility of arguments 
    cor_mat_p <- .threeDots_call("corr.test", x = t(x), adjust = adjust, 
                                 method = type, ...)$p
    ## was cor_mat_p <- corr.test(t(x), adjust = adjust, ...)$p
    cor_mat <- ifelse(cor_mat_p > threshold_correlation, 0, 1)
    colnames(cor_mat) <- rownames(cor_mat) <- rownames(x)
    return(cor_mat)
}

#' @name .bayes
#' @title Create a statistical network based on constraint-based structure 
#' learning algorithm
#' @description  \code{.bayes} infers a statistical network using 
#' constraint-based structure learning algorithm \code{fast.iamb} from the 
#' \code{bnlearn} package. \code{.bayes} extracts then the reported 
#' connections from running the \code{fast.iamb} function and assigns the
#' arcs of the discrete Bayesian network to binary values into a matrix that is 
#' returned by \code{.bayes}. 
#' @usage .bayes(x, ...)
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values 
#' @param 
#' @param  
#' @details For use of the parameters used in the \code{fast.iamb} function, 
#' refer to ?bnlearn::fast.iamb.
#' @return matrix, matrix with edges inferred from constraint-based structure 
#' learning algorithm \code{fast.iamb}
#' @author Thomas Naake, \email{thomasnaake @googlemail.com}
#' @examples .bayes(x, ...)
.bayes <- function(x, ...) {
    x_df <- data.frame(t(x))
    ## allow for compatibility of arguments 
    x_fast.iamb <- .threeDots_call("fast.iamb", x = x_df, ...)
    ## was x_fast.iamb <- fast.iamb(x_df, ...) 
    bs_mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    colnames(bs_mat) <- rownames(bs_mat) <- rownames(x)
    arcs_fast.iamb <- arcs(x_fast.iamb)
    for(i in 1:dim(arcs_fast.iamb)[1]) {
        bs_mat[arcs_fast.iamb[i, "from"], arcs_fast.iamb[i, "to"] ] <- 1} 
    bs_mat <- as.matrix(bs_mat)
    return(bs_mat)
}

#' @name .add_to_list
#' @title Add network to list
#' @description This helper function used in the function
#' \code{create_statistical_networks_list} adds a network to a list of networks.
#' @usage .add_to_list(l, name, object)
#' @param l list of networks
#' @param name character, name of newly created entry
#' @param object matrix containing the network to be added
#' @details Used internally in \code{create_statistical_networks_list}
#' @return list containing the existing networks and the added network 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples .add_to_list(l, name, object) 
#' @export
.add_to_list <- function(l, name, object) {
    if (!is.list(l)) stop("l is not a list")
    if (!is.character(name)) stop("name is not a character")
    if (!is.matrix(object)) stop("object is not a matrix")
    new_index <- length(l) + 1
    l[[new_index]] <- object
    names(l)[new_index] <- name
    return(l)
}



#' @name create_statistical_networks_list
#' @title Create a list of statistical networks
#' @description The function infers network topologies from statistical 
#' methods and returns matrices of these networks in a list. The function
#' includes functionality to caluclate statistical networks based on 
#' LASSO (L1 norm)-regression, random forests, context likelihood of 
#' relatedness (CLR), the algorithm for the reconstruction of accurate 
## cellular networks (ARACNE), Pearson correlation,
#' Spearman correlation and Constraint-based structure learning (Bayes). 
#' @usage create_statistical_networks_list(x, model, ...)
#' @param x matrix that contains intensity values of features/metabolites (rows)
#' per sample (columns). 
#' @param model, character vector containing the methods that will be used 
#' ("lasso", "randomForest", "clr", "aracne", "pearson", "spearman", "bayes")
#' @param ... parameters passed to the functions  \code{.lasso}, 
#' \code{.randomForest}, \code{.clr}, \code{.aracne}, \code{.correlation} and/or
#' \code{.bayes}
#' @details \code{create_statistical_networks_list} calls the function
#' \code{.lasso}, \code{.randomForest}, \code{.clr}, \code{.aracne}, 
#' \code{.pearson}, \code{.spearman} and/or \code{.bayes} as specified by 
#' \code{model}. It will create network(s) using the specified methods and will 
#' return a list containing the binary network (if model is of length 1) or
#' append these binary networks to a list (if model is of length > 1). 
#' Internally x will be z-scaled and the z-scaled object will be used in 
#' \code{.lasso}, \code{.clr} and/or \code{.aracne}.
#' @return list containing the respective statistical networks specified by 
#' \code{model}
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples create_statistical_networks_list(x, model, ...) 
#' @export
create_statistical_networks_list <- function(x, model, ...) {
    
    ## check if model complies with the implemented model and return error 
    ## if not so
    if (!(all(model %in% c("lasso", "randomForest", "clr", "aracne", 
                           "pearson", "spearman", "bayes"))))
        stop("method not implemented in create_statistical_networks_list")
        
    ## check if x is numeric matrix and return error if not so
    if (mode(x) != "numeric") stop("x is not a numerical matrix")
    
    x_z <- apply(x, 1, function(x) (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))
    x_z <- t(x_z)
    
    l <- list()
    
    ## add entry for lasso if "lasso" is in model
    if ("lasso" %in% model) {
        lasso <- .lasso(x_z, ...)
        l <- .add_to_list(l, "lasso", lasso)
    }
    ## add entry for randomForest if "randomForest" is in model
    if ("randomForest" %in% model) {
        randomForest <- .randomForest(x, ...)
        l <- .add_to_list(l, "randomForest", randomForest)
    }
    if (any(c("clr", "aracne") %in% model)) {
        mi_x_z <- cmi(t(x_z))$bcmi
        rownames(mi_x_z) <- colnames(mi_x_z) <- rownames(x)
    }
    
    ## add entry for clr if "clr" is in model
    if ("clr" %in% model) {
        clr <- .threeDots_call(".clr", mi = mi_x_z, ...)
                               ##.clr(mi = mi_x_z, ...)
        l <- .add_to_list(l, "clr", clr)
    }
    ## add entry for aracne if "aracne" is in model
    if ("aracne" %in% model) {
        aracne <- .threeDots_call(".aracne", mi = mi_x_z, ...)
        l <- .add_to_list(l, "aracne", aracne)
    }
    ## add entry for pearson if "pearson" is in model
    if ("pearson" %in% model) {
        pearson <- .correlation(x, type = "pearson", ...)
        l <- .add_to_list(l, "pearson", pearson)
    }
    ## add entry for spearman if "spearman" is in model
    if ("spearman" %in% model) {
        spearman <- .correlation(x, type = "spearman", ...)
        l <- .add_to_list(l, "spearman", spearman)
    }
    ## add entry for bayes if "bayes" is in model
    if ("bayes" %in% model) {
        bayes <- .bayes(x, ...)
        l <- .add_to_list(l, "bayes", bayes)
    }
    return(l)
}


#' @name consensus_network
#' @title Create a consensus network of statistical networks
#' @description The function takes a list of parameters (\code{l}) as input and
#' creates a consensus network from these networks by calling the function 
#' \code{consensus} from the \code{sna} package. Depending on the chosen 
#' \code{method } in \code{consensus}, the threshold of the consensus network
#' should be chosen accordingly to report a connection by different 
#' statistical methods. 
#' @usage consensus_network(l, threshold = 1, ...)
#' @param l list, each entry of the list contains a network
#' @param threshold numerical value, when combining the networks the threshold 
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
#' @examples consensus_network(l, threshold = 1, ... )
#' @export
consensus_network <- function(l, threshold = 1, ...) {
    
    if (!is.numeric(threshold)) stop("threshold is not numeric")
    if (!is.list(l)) stop("l is not a list")
    
    ## check compatibility of matrices
    ncol_1 <- ncol(l[[1]])
    nrow_1 <- nrow(l[[1]])
    rownames_1 <- rownames(l[[1]])
    colnames_1 <- colnames(l[[1]])
    
    if (!all(rownames_1 == colnames_1)) stop("colnames and rownames are not identical")
    
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
    consensus_mat <- .threeDots_call(consensus, dat = l, ...)
    ## was sna::consensus(dat = l, ...)
    
    ##if (method == "central.graph") threshold <- 1
    consensus_mat <- ifelse(consensus_mat >= threshold, 1, 0)
    
    rownames(consensus_mat) <- colnames(consensus_mat) <- colnames(l[[1]])
    return(consensus_mat)
}

#' @name .threeDots_call
#' @title Check if passed arguments match the function's formal arguments and
#' call the function with the checked arguments
#' @description The function \code{.threeDots_call} gets the formal arguments
#' of a function \code{fun} and checks if the passed arguments \code{...} 
#' matches the formal arguments. \code{.threeDots_call} will remove 
#' duplicated arguments. \code{.threeDots_call} will call the function 
#' \code{fun} with the filtered arguments and will return the result. 
#' @usage .threeDots_call(fun, ...)
#' @param fun function to check for arguments and to call
#' @param ... arguments to be tested to be passed to fun
#' @details Used internally in \code{.lasso}, \code{.randomForest}, 
#' \code{.correlation}, \code{.bayes}, \code{consensus_network}
#' @return Function call with passed arguments 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples .threeDots_call("sd", x = 1:10, y = 1:10)
#' ## in contrast to the above example, the following example will result in an 
#' ## error
#' \donotrun{sd(x = 1:10, y = 1:10)}
.threeDots_call <- function(fun, ...) {
    formal_args <- formalArgs(fun)
    args <- list(...)
    if (any(duplicated(names(args)))) stop("duplicated args in ...")
    
    input <- args[names(args) %in% formal_args]
    input <- input[!duplicated(names(input))]
    res <- do.call(fun, input)
    return(res)
}

stat_net <- create_statistical_network(x, model = c("pearson", "spearman"), adjust_correlation = "BH", cutoff = 0.7, PFER = 0.7, method = "PCA.reweight")

l <- create_statistical_networks_list(x, model = c("pearson", "spearman"),  adjust_correlation = "BH", cutoff = 0.7, PFER = 0.7)
stat_net2 <- consensus_network(l = l, threshold = 1, method = "PCA.reweight")

#' @name create_statistical_network
#' @title Create statistical network
#' @description 
#' @usage create_statistical_network(x, model, threshold = 1, ...)
#' @param x matrix that contains intensity values of features/metabolites (rows)
#' per sample (columns). 
#' @param model, character vector containing the model that will be used 
#' ("lasso", "randomForest", "clr", "aracne", "pearson", "spearman", "bayes")
#' @param ... parameters passed to the functions  \code{.lasso}, 
#' \code{.randomForest}, \code{.clr}, \code{.aracne}, \code{.correlation},
#' \code{.bayes} and/or \code{consensus_network}
#' @return matrix, containing binary values if a connection is present or not
#' @details create_statistical_network is a wrapper function for the 
#' functions \code{create_statistical_networks_list} and 
#' \code{consensus_network}. See \code{?create_statistical_networks_list} and 
#' \code{?consensus_network} for further details. The function 
#' \code{create_statistical_networks_list}
#' includes functionality to caluclate statistical networks based on 
#' LASSO (L1 norm)-regression, random forests, context likelihood of 
#' relatedness (CLR), the algorithm for the reconstruction of accurate 
## cellular networks (ARACNE), Pearson correlation,
#' Spearman correlation and Constraint-based structure learning (Bayes). 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples create_statistical_network(x, model, threshold = 1, ...)
#' @export
create_statistical_network <- function(x, model, threshold = 1, ...) {
    ##l <- .threeDots_call(create_statistical_networks_list, x = x, model = model, ...)
    ##consensus_mat <- .threeDots_call(consensus_network, l = l, threshold = threshold, ...)
    l <- create_statistical_networks_list(x = x, model = model, ...)
    consensus_mat <- consensus_network(l = l, threshold = threshold, ...)
    return(consensus_mat)
}


