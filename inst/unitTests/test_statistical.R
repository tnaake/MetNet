## load RUnit
library(RUnit)

## create toy example data set
data("mat_test", package = "MetNet")

## START unit test lasso ##
lasso_mat <- lasso(t(mat_test_z), parallel = FALSE, PFER = 0.75, cutoff = 0.95)
lasso_mat_parallel <- lasso(t(mat_test_z)[1:5,], parallel = TRUE, 
    PFER = 0.75, cutoff = 0.95)
test_lasso <- function() {
    checkException(lasso(NULL, PFER = 0.75, cutoff = 0.95), 
        msg = "must be coercible to non-negative integer")
    
    checkException(lasso(mat_test))
    checkEquals(rownames(lasso_mat), colnames(lasso_mat))
    checkEquals(rownames(lasso_mat), rownames(mat_test))
    checkEquals(ncol(lasso_mat), nrow(lasso_mat))
    checkEquals(nrow(lasso_mat), nrow(mat_test))
    checkTrue(is.numeric(lasso_mat))
    checkTrue(is.matrix(lasso_mat))
    checkTrue(max(lasso_mat) <= 1)
    checkTrue(min(lasso_mat) >= 0)
}
## END unit test lasso ##

## START unit test randomForest ##
rf_mat <- randomForest(mat_test)
test_randomForest <- function() {
    checkException(randomForest(NULL), msg = "find an inherited method for ")
    checkEquals(sum(rf_mat), 7, tolerance = 1e-06)
    checkEquals(rownames(rf_mat), colnames(rf_mat))
    checkEquals(rownames(rf_mat), rownames(mat_test))
    checkEquals(ncol(rf_mat), nrow(rf_mat))
    checkEquals(nrow(rf_mat), nrow(mat_test))
    checkTrue(is.numeric(rf_mat))
    checkTrue(is.matrix(rf_mat))
    checkTrue(max(rf_mat) <= 1)
    checkTrue(min(rf_mat) >= 0)
}
## END unit test randomForest ##

## START unit test clr ##
mi_mat_test_z <- mpmi::cmi(mat_test_z)$bcmi
rownames(mi_mat_test_z) <- colnames(mi_mat_test_z) <- colnames(mat_test_z)
clr_mat <- clr(mi_mat_test_z)
test_clr <- function() {
    checkException(clr(NULL), msg = "must be a matrix")
    checkEquals(sum(clr_mat), 9.042674, tolerance = 1e-06)
    checkEquals(rownames(clr_mat), colnames(clr_mat))
    checkEquals(rownames(clr_mat), rownames(mat_test))
    checkEquals(ncol(clr_mat), nrow(clr_mat))
    checkEquals(nrow(clr_mat), nrow(mat_test))
    checkTrue(is.numeric(clr_mat))
    checkTrue(is.matrix(clr_mat))
    checkTrue(max(clr_mat) <= 1)
    checkTrue(min(clr_mat) >= 0)
}
## END unit test clr ##

## START unit test aracne ##
aracne_mat <- aracne(mi_mat_test_z)
test_aracne <- function() {
    checkException(aracne(NULL, eps = 0.05), msg = "must be a matrix")
    checkException(aracne(mi_mat_test_z, eps = "a"), msg = "in foreign")
    checkEquals(sum(aracne_mat), 22.92893, tolerance = 1e-06)
    checkEquals(rownames(aracne_mat), colnames(aracne_mat))
    checkEquals(rownames(aracne_mat), rownames(mat_test))
    checkEquals(ncol(aracne_mat), nrow(aracne_mat))
    checkEquals(nrow(aracne_mat), nrow(mat_test))
    checkTrue(is.numeric(aracne_mat))
    checkTrue(is.matrix(aracne_mat))
    checkTrue(max(aracne_mat) <= 1.3)
    checkTrue(min(aracne_mat) >= 0)
}
## END unit test aracne ##

## START unit test correlation ##
correlation_p_mat <- correlation(mat_test, type = "pearson")
correlation_p_p_mat <- correlation(mat_test, type = "pearson_partial")
correlation_p_sp_mat <- correlation(mat_test, type = "pearson_semipartial")
correlation_s_mat <- correlation(mat_test, type = "spearman")
correlation_s_p_mat <- correlation(mat_test, type = "spearman_partial")
correlation_s_sp_mat <- correlation(mat_test, type = "spearman_semipartial")

test_correlation <- function() {
    
    checkException(correlation(NULL, type = "pearson"), msg = "not a matrix")
    checkException(correlation(mat_test, type = "a"), msg = "not found")
    checkException(correlation(mat_test, type = "pearson", 
        correlation_adjust = "abc"), msg = "should be one of ")
    
    ## pearson
    checkIdentical(correlation_p_mat, 
        WGCNA::corAndPvalue(t(mat_test), method = "pearson")$p)
    checkIdentical(correlation(mat_test, alternative = "greater"), 
        WGCNA::corAndPvalue(t(mat_test), method = "pearson", 
            alternative = "greater")$p)
    checkIdentical(correlation(mat_test, alternative = "less"), 
        WGCNA::corAndPvalue(t(mat_test), method = "pearson", 
            alternative = "less")$p)
    checkIdentical(correlation(mat_test, alternative = "two.sided"), 
        WGCNA::corAndPvalue(t(mat_test), method = "pearson", 
            alternative = "two.sided")$p)
    checkEquals(sum(correlation_p_mat), 0.3561751, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "bonferroni")), 
        11.28729, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "none")), 
        0.3561751, tolerance = 1e-06)
    checkEquals(rownames(correlation_p_mat), colnames(correlation_p_mat))
    checkEquals(rownames(correlation_p_mat), rownames(mat_test))
    checkEquals(ncol(correlation_p_mat), nrow(correlation_p_mat))
    checkEquals(nrow(correlation_p_mat), nrow(mat_test))
    checkTrue(is.numeric(correlation_p_mat))
    checkTrue(is.matrix(correlation_p_mat))
    checkTrue(max(correlation_p_mat) <= 1)
    checkTrue(min(correlation_p_mat) >= 0)
    
    ## spearman
    checkIdentical(correlation_s_mat, 
        WGCNA::corAndPvalue(t(mat_test), method = "spearman")$p)
    checkEquals(sum(correlation_s_mat), 0.4854916, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "spearman")), 12, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "none", 
        type = "spearman")), 0.4854916, tolerance = 1e-06)
    checkEquals(rownames(correlation_s_mat), colnames(correlation_s_mat))
    checkEquals(rownames(correlation_s_mat), rownames(mat_test))
    checkEquals(ncol(correlation_s_mat), nrow(correlation_s_mat))
    checkEquals(nrow(correlation_s_mat), nrow(mat_test))
    checkTrue(is.numeric(correlation_s_mat))
    checkTrue(is.matrix(correlation_s_mat))
    checkTrue(max(correlation_s_mat) <= 1)
    checkTrue(min(correlation_s_mat) >= 0)
    
    ## partial pearson
    checkIdentical(correlation_p_p_mat, 
        ppcor::pcor(t(mat_test), method = "pearson")$p.value)
    checkEquals(sum(correlation_p_p_mat), 19.78391, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "pearson_partial")), 42, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "none", 
        type = "pearson_partial")), 19.78391, tolerance = 1e-06)
    checkEquals(rownames(correlation_p_p_mat), colnames(correlation_p_p_mat))
    checkEquals(rownames(correlation_p_p_mat), rownames(mat_test))
    checkEquals(ncol(correlation_p_p_mat), nrow(correlation_p_p_mat))
    checkEquals(nrow(correlation_p_p_mat), nrow(mat_test))
    checkTrue(is.numeric(correlation_p_p_mat))
    checkTrue(is.matrix(correlation_p_p_mat))
    checkTrue(max(correlation_p_p_mat) <= 1)
    checkTrue(min(correlation_p_p_mat) >= 0)
    
    ## semi-partial pearson
    checkIdentical(correlation_p_sp_mat, 
        ppcor::spcor(t(mat_test), method = "pearson")$p.value)
    checkEquals(sum(correlation_p_sp_mat), 35.74922, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "pearson_semipartial")), 42, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "none", 
        type = "pearson_semipartial")), 35.74922, tolerance = 1e-06)
    checkEquals(rownames(correlation_p_sp_mat), colnames(correlation_p_sp_mat))
    checkEquals(rownames(correlation_p_sp_mat), rownames(mat_test))
    checkEquals(ncol(correlation_p_sp_mat), nrow(correlation_p_sp_mat))
    checkEquals(nrow(correlation_p_sp_mat), nrow(mat_test))
    checkTrue(is.numeric(correlation_p_sp_mat))
    checkTrue(is.matrix(correlation_p_sp_mat))
    checkTrue(max(correlation_p_sp_mat) <= 1)
    checkTrue(min(correlation_p_sp_mat) >= 0)
    
    ## partial spearman
    checkTrue(all(correlation_s_p_mat == 
        ppcor::pcor(t(mat_test), method = "spearman")$p.value, na.rm = TRUE))
    checkEquals(sum(correlation_s_p_mat, na.rm = TRUE), 0.9986576, 
        tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "spearman_partial"), na.rm = TRUE), 12, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "none", 
        type = "spearman_partial"), na.rm = TRUE), 0.9986576, tolerance = 1e-06)
    checkEquals(rownames(correlation_s_p_mat), colnames(correlation_s_p_mat))
    checkEquals(rownames(correlation_s_p_mat), rownames(mat_test))
    checkEquals(ncol(correlation_s_p_mat), nrow(correlation_s_p_mat))
    checkEquals(nrow(correlation_s_p_mat), nrow(mat_test))
    checkTrue(is.numeric(correlation_s_p_mat))
    checkTrue(is.matrix(correlation_s_p_mat))
    checkTrue(max(correlation_s_p_mat, na.rm = TRUE) <= 1)
    checkTrue(min(correlation_s_p_mat, na.rm = TRUE) >= 0)
    
    ## semi-partial spearman
    checkTrue(all(correlation_s_sp_mat ==
        ppcor::spcor(t(mat_test), method = "spearman")$p.value, na.rm = TRUE))
    checkEquals(sum(correlation_s_sp_mat, na.rm = TRUE), 
        0.4993288, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "spearman_semipartial"), na.rm = TRUE), 6, tolerance = 1e-06)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "none", 
        type = "spearman_semipartial"), na.rm = TRUE), 0.4993288, 
        tolerance = 1e-06)
    checkEquals(rownames(correlation_s_sp_mat), colnames(correlation_s_sp_mat))
    checkEquals(rownames(correlation_s_sp_mat), rownames(mat_test))
    checkEquals(ncol(correlation_s_sp_mat), nrow(correlation_s_sp_mat))
    checkEquals(nrow(correlation_s_sp_mat), nrow(mat_test))
    checkTrue(is.numeric(correlation_s_sp_mat))
    checkTrue(is.matrix(correlation_s_sp_mat))
    checkTrue(max(correlation_s_sp_mat, na.rm = TRUE) <= 1)
    checkTrue(min(correlation_s_sp_mat, na.rm = TRUE) >= 0)
}
## END unit test correlation ##

## START unit test bayes ##
bayes_mat <- bayes(mat_test)
test_bayes <- function() {
    checkTrue(sum(bayes_mat) < 21 & sum(bayes_mat) > 17)
    checkEquals(rownames(bayes_mat), colnames(bayes_mat))
    checkEquals(rownames(bayes_mat), rownames(mat_test))
    checkEquals(ncol(bayes_mat), nrow(bayes_mat))
    checkEquals(nrow(bayes_mat), nrow(mat_test))
    checkTrue(is.numeric(bayes_mat))
    checkTrue(is.matrix(bayes_mat))
    checkTrue(max(bayes_mat) <= 1)
    checkTrue(min(bayes_mat) >= 0)
}
## END unit test bayes ##

## START unit test addToList ##
l <- list()
l <- MetNet:::addToList(l, "newEntry", matrix())
test_addToList <- function() {
    checkException(MetNet:::addToList(l, "newEntry", NULL), 
        msg = "not a matrix")
    checkException(MetNet:::addToList(l, NULL, matrix()), 
        msg = "not a character")
    checkException(MetNet:::addToList(NULL, "newEntry", matrix()),
        msg = "not a list")
    checkEquals(length(l), 1)
    checkTrue(is.matrix(l[[1]]))
    checkEquals(l[[1]], matrix())
}
## END unit test addToList ##


## START unit test statistical ##
stat_adj_l <- statistical(mat_test,
    model = c("lasso", "randomForest", "clr", "aracne", "pearson", 
        "pearson_partial", "pearson_semipartial", "spearman", 
        "spearman_partial", "spearman_semipartial", "bayes"), 
    correlation_adjust = "bonferroni", PFER = 0.75, cutoff = 0.95)
test_statistical <- function() {
    checkException(statistical(NULL, model = "lasso"), msg = "not a numerical")
    checkException(statistical(mat_test, model = "foo"), 
        msg = "model not implemented in ")
    checkException(statistical(mat_test, model = c("lasso")), 
        msg = "Two of the three arguments 'PFER', 'cutoff' and 'q'" )
    
    ## take a high tolerance value for LASSO, randomForest and bayes 
    ## since these models are probabilistic
    tmp <- lasso_mat; diag(tmp) <- NaN
    checkEquals(stat_adj_l[["lasso"]], tmp, tolerance = 5e-01)
    tmp <- rf_mat; diag(tmp) <- NaN
    checkEquals(stat_adj_l[["randomForest"]], tmp, tolerance = 5e-01)
    tmp <- bayes_mat; diag(tmp) <- NaN
    checkEquals(stat_adj_l[["bayes"]], tmp, tolerance = 5e-01)
    tmp <- clr_mat; diag(tmp) <- NaN
    checkIdentical(stat_adj_l[["clr"]], tmp)
    tmp <- aracne_mat; diag(tmp) <- NaN
    checkIdentical(stat_adj_l[["aracne"]], tmp)
    tmp <-  correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "pearson")
    diag(tmp) <- NaN
    checkIdentical(stat_adj_l[["pearson"]], tmp)
    tmp <- correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "pearson_partial")
    diag(tmp) <- NaN
    checkIdentical(stat_adj_l[["pearson_partial"]], tmp)
    tmp <- correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "pearson_semipartial")
    diag(tmp) <- NaN
    checkIdentical(stat_adj_l[["pearson_semipartial"]], tmp)
    tmp <- correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "spearman")
    diag(tmp) <- NaN
    checkIdentical(stat_adj_l[["spearman"]], tmp)
    tmp <- correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "spearman_partial")
    diag(tmp) <- NaN
    checkIdentical(stat_adj_l[["spearman_partial"]], tmp)
    tmp <- correlation(mat_test, correlation_adjust = "bonferroni", 
        type = "spearman_semipartial")
    diag(tmp) <- NaN
    checkIdentical(stat_adj_l[["spearman_semipartial"]], tmp)
    checkEquals(length(stat_adj_l), 11)
    checkEquals(as.numeric(lapply(stat_adj_l, nrow)), rep(7, 11))
    checkEquals(as.numeric(lapply(stat_adj_l, ncol)), rep(7, 11))
    checkEquals(as.character(unlist((lapply(stat_adj_l, rownames)))),
                        rep(c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), 11))
    checkEquals(as.character(unlist((lapply(stat_adj_l, colnames)))),
                rep(c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), 11))
    checkEquals(names(stat_adj_l), 
        c("lasso", "randomForest", "clr", "aracne", "pearson", 
          "pearson_partial", "pearson_semipartial", "spearman",         
          "spearman_partial", "spearman_semipartial", "bayes"))
    checkTrue(all(unlist(lapply(stat_adj_l, function(x) is.numeric(x)))))
}
## END unit test statistical ##


## START unit test getLinks ##
mat <- matrix(0:8, ncol = 3, nrow = 3)
getLinks_df <- MetNet:::getLinks(mat, decreasing = TRUE, exclude = "== 0")

test_getLinks <- function() {
    
    checkException(MetNet:::getLinks(NULL), msg = "argument is of length zero")
    checkException(MetNet:::getLinks(mat[, 1:2]), msg = "not a square matrix")
    checkException(MetNet:::getLinks(mat, exclude = "foo"), 
        msg = "object 'matfoo' not found")
    
    ## exclude = "== 0"
    checkTrue(is.data.frame(getLinks_df))
    checkEquals(getLinks_df$row, rep(c(1, 2, 3), 3))
    checkEquals(getLinks_df$col, rep(c(1, 2, 3), each = 3))
    checkEquals(getLinks_df$confidence, c(NaN, 1:8))
    checkEquals(getLinks_df$rank, c(NaN, 8:1))
    
    ## exclude = NULL
    getLinks_df <- MetNet:::getLinks(mat, decreasing = TRUE, exclude = NULL)
    checkEquals(getLinks_df$row, rep(c(1, 2, 3), 3))
    checkEquals(getLinks_df$col, rep(c(1, 2, 3), each = 3))
    checkEquals(getLinks_df$confidence, 0:8)
    checkEquals(getLinks_df$rank, 9:1)
    
    ## decreasing = FALSE
    getLinks_df <- MetNet:::getLinks(mat, decreasing = FALSE, exclude = NULL)
    checkEquals(getLinks_df$row, rep(c(1, 2, 3), 3))
    checkEquals(getLinks_df$col, rep(c(1, 2, 3), each = 3))
    checkEquals(getLinks_df$confidence, 0:8)
    checkEquals(getLinks_df$rank, 1:9)
}
## END unit test getLinks ##


## START unit test threshold  ##
## remove partial/semipartial correlation from stat_adj_l
stat_adj_l_cut <- stat_adj_l[!names(stat_adj_l) %in% 
    c("pearson_partial", "pearson_semipartial", "spearman_partial", 
        "spearman_semipartial", "bayes", "randomForest", "lasso")]
args_thr <- list(clr = 0.5, aracne = 0.8, pearson = 0.05, spearman = 0.05, 
    threshold = 1)
thr_thr <- threshold(stat_adj_l_cut, type = "threshold", args = args_thr)

args_top <- list(n = 5)
thr_top1 <- threshold(stat_adj_l_cut, type = "top1", args = args_top)
thr_top2 <- threshold(stat_adj_l_cut, type = "top2", args = args_top)
thr_mean <- threshold(stat_adj_l_cut, type = "mean", args = args_top)

test_threshold <- function() {
    ## test arguments
    checkException(threshold(NULL, type = "threshold", args = args_thr), 
        msg = "consensus requires graphs of identical order")
    checkException(threshold(1:3, type = "threshold", args = args_thr),
        msg = "attempt to select less than one element in")
    args_foo <- list(clr = 0.5, aracne = 0.8, threshold = 1)
    checkException(threshold(
        data.frame(clr = 1:3, aracne = 1:3), 
        type = "threshold", args = args_foo), 
        msg = "input must be an adjacency matrix/array, network, or list")
    checkException(threshold(cbind(clr = c("a", "b", "c"), 
        aracne = c("a", "b", "c")), type = "threshold", args = args_foo),
        msg = "attempt to select less than one element in")
    args_thr_double <- list(lasso = 0.8, randomForest = 0.2, 
        clr = 0.5, clr = 0.5, aracne = 0.8, pearson = 0.05, spearman = 0.05, 
        threshold = 1)
    checkException(
        threshold(stat_adj_l_cut, type = "threshold", args = args_thr_double),
        msg = "contain duplicated entries")
    checkException(
        threshold(stat_adj_l_cut, type = "foo", args = args_thr),
        msg = "type not in")
    
    ## check that args contains all models
    checkException(threshold(stat_adj_l_cut, type = "threshold", args_thr[1:3]),
        msg = "does not contain entries for all 'model's in 'statistical'")
    
    ## check that args contains threshold
    checkException(threshold(stat_adj_l_cut, type = "threshold", args_thr[1:4]),
        msg = "'args' does not contain entry 'threshold' of length 1")
    
    ## check args for top1, top2, mean
    checkException(threshold(stat_adj_l, type = "top1", args = list(x = 1)),
        msg = "does not contain the numeric entry `n` of length 1")
    checkException(threshold(stat_adj_l, type = "top2", args = list(x = 1)),
        msg = "does not contain the numeric entry `n` of length 1")
    checkException(threshold(stat_adj_l, type = "mean", args = list(x = 1)),
        msg = "does not contain the numeric entry `n` of length 1")
    checkException(threshold(stat_adj_l, type = "top1", args = list(n = 1:2)),
        msg = "does not contain the numeric entry `n` of length 1")
    checkException(threshold(stat_adj_l, type = "top2", args = list(n = 1:2)),
        msg = "does not contain the numeric entry `n` of length 1")
    checkException(threshold(stat_adj_l, type = "mean", args = list(n = 1:2)),
        msg = "does not contain the numeric entry `n` of length 1")
    checkException(threshold(stat_adj_l, type = "top1", args = list(n = "a")),
        msg = "does not contain the numeric entry `n` of length 1")
    checkException(threshold(stat_adj_l, type = "top2", args = list(n = "a")),
        msg = "does not contain the numeric entry `n` of length 1")
    checkException(threshold(stat_adj_l, type = "mean", args = list(n = "a")),
        msg = "does not contain the numeric entry `n` of length 1")
    
    ## check output
    checkTrue(is.matrix(thr_thr))
    checkTrue(is.matrix(thr_top1))
    checkTrue(is.matrix(thr_top2))
    checkTrue(is.matrix(thr_mean))
    checkTrue(is.numeric(thr_thr))
    checkTrue(is.numeric(thr_top1))
    checkTrue(is.numeric(thr_top2))
    checkTrue(is.numeric(thr_mean))
    checkEquals(rownames(thr_thr), c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkEquals(colnames(thr_thr), c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkEquals(rownames(thr_top1), c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkEquals(colnames(thr_top1), c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkEquals(rownames(thr_top2), c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkEquals(colnames(thr_top2), c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkEquals(rownames(thr_mean), c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkEquals(colnames(thr_mean), c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkEquals(sum(thr_thr), 30)
    checkEquals(sum(thr_top1), 20)
    checkEquals(sum(thr_top2), 14)
    checkEquals(sum(thr_mean), 10)
    checkTrue(all(thr_thr %in% c(0, 1)))
    checkTrue(all(thr_top1 %in% c(0, 1)))
    checkTrue(all(thr_top2 %in% c(0, 1)))
    checkTrue(all(thr_mean %in% c(0, 1)))
}
## END unit test threshold  ##


## START unit test topKnet ##
ranks <- matrix(c(c(1, 2, 3), c(2, 1, 3)), ncol = 2)
 
test_topKnet <- function() {
    
    ## type of ranks
    checkException(MetNet:::topKnet(ranks = 1:3))
    checkException(MetNet:::topKnet(ranks = data.frame(x = 1:3, y = 3:1), 
        type = "top1"))
    checkException(MetNet:::topKnet(ranks = matrix(c("a", "b", "c")), 
        type = "top1"))
    
    ## type argument
    checkException(MetNet:::topKnet(ranks = ranks, type = "foo"))
    
    ## check results
    checkEquals(MetNet:::topKnet(ranks = ranks, type = "top1"), c(1, 1, 3))
    checkEquals(MetNet:::topKnet(ranks = ranks, type = "top2"), c(2, 2, 3))
    checkEquals(MetNet:::topKnet(ranks = ranks, type = "mean"), c(1.5, 1.5, 3))
    
    ## matrix with ncol 1
    checkEquals(MetNet:::topKnet(ranks = matrix(1:3, nrow = 3),
        type = "top1"), 1:3)
    checkException(MetNet:::topKnet(ranks = matrix(1:3, nrow = 3), 
        type = "top2"))
    checkEquals(MetNet:::topKnet(ranks = matrix(1:3, nrow = 3), 
        type = "mean"), 1:3)
}
## END unit test topKnet ##


## START unit test threeDotsCall ##
test_threeDots_call <- function() {
    checkException(MetNet:::threeDotsCall("mean", x = 1:10, x = 1:10), 
                   msg = "duplicated args in ...")
    checkEquals(MetNet:::threeDotsCall("mean", x = 1:10, foo = 1), mean(1:10))
}
## END unit test threeDotsCall ##

