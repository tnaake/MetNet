## load RUnit
library(RUnit)

## create toy example data set
data("mat_test", package = "MetNet")

## START unit test lasso ##
lasso_mat <- lasso(t(mat_test_z), parallel = FALSE, PFER = 0.75, cutoff = 0.95)
test_lasso <- function() {
    checkException(lasso(mat_test))
    checkEquals(rownames(lasso_mat), colnames(lasso_mat))
    checkEquals(rownames(lasso_mat), rownames(mat_test))
    checkEquals(ncol(lasso_mat), nrow(lasso_mat))
    checkEquals(nrow(lasso_mat), nrow(mat_test))
    checkTrue(is.numeric(lasso_mat))
    checkTrue(is.matrix(lasso_mat))
    checkTrue(all(lasso_mat %in% c(0, 1)))
}
## END unit test lasso ##

## START unit test randomForest ##
rf_mat <- randomForest(mat_test, parallel = FALSE, num.cores = 1)
test_randomForest <- function() {
    checkEquals(rownames(rf_mat), colnames(rf_mat))
    checkEquals(rownames(rf_mat), rownames(mat_test))
    checkEquals(ncol(rf_mat), nrow(rf_mat))
    checkEquals(nrow(rf_mat), nrow(mat_test))
    checkTrue(is.numeric(rf_mat))
    checkTrue(is.matrix(rf_mat))
    checkTrue(all(rf_mat %in% c(0, 1)))
}
## END unit test randomForest ##

## START unit test clr ##
mi_mat_test_z <- mpmi::cmi(mat_test_z)$bcmi
rownames(mi_mat_test_z) <- colnames(mi_mat_test_z) <- colnames(mat_test_z)
clr_mat <- clr(mi_mat_test_z)
test_clr <- function() {
    checkException(clr(mi_mat_test_z, clr_threshold = "a"))
    checkEquals(sum(clr_mat), 30)
    checkEquals(rownames(clr_mat), colnames(clr_mat))
    checkEquals(rownames(clr_mat), rownames(mat_test))
    checkEquals(ncol(clr_mat), nrow(clr_mat))
    checkEquals(nrow(clr_mat), nrow(mat_test))
    checkTrue(is.numeric(clr_mat))
    checkTrue(is.matrix(clr_mat))
    checkTrue(all(clr_mat %in% c(0, 1)))
}
## END unit test clr ##

## START unit test aracne ##
aracne_mat <- aracne(mi_mat_test_z)
test_aracne <- function() {
    checkException(aracne(mi_mat_test_z, aracne_threshold = "a"))
    checkException(aracne(mi_mat_test_z, eps = "a"))
    checkEquals(sum(aracne_mat), 31)
    checkEquals(rownames(aracne_mat), colnames(aracne_mat))
    checkEquals(rownames(aracne_mat), rownames(mat_test))
    checkEquals(ncol(aracne_mat), nrow(aracne_mat))
    checkEquals(nrow(aracne_mat), nrow(mat_test))
    checkTrue(is.numeric(aracne_mat))
    checkTrue(is.matrix(aracne_mat))
    checkTrue(all(aracne_mat %in% c(0, 1)))
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
    checkException(correlation(mat_test, correlation_threshold = "a"))
    checkEquals(sum(correlation_p_mat), 47)
    checkEquals(sum(correlation(mat_test,
        correlation_adjust = "bonferroni", correlation_threshold = 1e-11)), 21)
    checkEquals(sum(correlation(mat_test,
        correlation_adjust = "none", correlation_threshold = 1e-11)), 33)
    checkEquals(sum(correlation_s_mat), 47)
    checkEquals(rownames(correlation_p_mat), colnames(correlation_p_mat))
    checkEquals(rownames(correlation_p_p_mat), colnames(correlation_p_p_mat))
    checkEquals(rownames(correlation_p_sp_mat), colnames(correlation_p_sp_mat))
    checkEquals(rownames(correlation_s_mat), colnames(correlation_s_mat))
    checkEquals(rownames(correlation_s_p_mat), colnames(correlation_s_p_mat))
    checkEquals(rownames(correlation_s_sp_mat), colnames(correlation_s_sp_mat))
    checkEquals(rownames(correlation_p_mat), rownames(mat_test))
    checkEquals(rownames(correlation_p_p_mat), rownames(mat_test))
    checkEquals(rownames(correlation_p_sp_mat), rownames(mat_test))
    checkEquals(rownames(correlation_s_mat), rownames(mat_test))
    checkEquals(rownames(correlation_s_p_mat), rownames(mat_test))
    checkEquals(rownames(correlation_s_sp_mat), rownames(mat_test))
    checkEquals(ncol(correlation_p_mat), nrow(correlation_p_mat))
    checkEquals(ncol(correlation_p_p_mat), nrow(correlation_p_p_mat))
    checkEquals(ncol(correlation_p_sp_mat), nrow(correlation_p_sp_mat))
    checkEquals(ncol(correlation_s_mat), nrow(correlation_s_mat))
    checkEquals(ncol(correlation_s_p_mat), nrow(correlation_s_p_mat))
    checkEquals(ncol(correlation_s_sp_mat), nrow(correlation_s_sp_mat))
    checkEquals(nrow(correlation_p_mat), nrow(mat_test))
    checkEquals(nrow(correlation_p_p_mat), nrow(mat_test))
    checkEquals(nrow(correlation_p_sp_mat), nrow(mat_test))
    checkEquals(nrow(correlation_s_mat), nrow(mat_test))
    checkEquals(nrow(correlation_s_p_mat), nrow(mat_test))
    checkEquals(nrow(correlation_s_sp_mat), nrow(mat_test))
    checkTrue(is.numeric(correlation_p_mat))
    checkTrue(is.numeric(correlation_p_p_mat))
    checkTrue(is.numeric(correlation_p_sp_mat))
    checkTrue(is.numeric(correlation_s_mat))
    checkTrue(is.numeric(correlation_s_p_mat))
    checkTrue(is.numeric(correlation_s_sp_mat))
    checkTrue(is.matrix(correlation_p_mat))
    checkTrue(is.matrix(correlation_p_p_mat))
    checkTrue(is.matrix(correlation_p_sp_mat))
    checkTrue(is.matrix(correlation_s_mat))
    checkTrue(is.matrix(correlation_s_p_mat))
    checkTrue(is.matrix(correlation_s_sp_mat))
    checkTrue(all(correlation_p_mat %in% c(0, 1)))
    checkTrue(all(correlation_p_p_mat %in% c(0, 1)))
    checkTrue(all(correlation_p_sp_mat %in% c(0, 1)))
    checkTrue(all(correlation_s_mat %in% c(0, 1)))
    checkTrue(all(correlation_s_p_mat %in% c(0, 1, NA)))
    checkTrue(all(correlation_s_sp_mat %in% c(0, 1, NA)))
}
## END unit test correlation ##

## START unit test bayes ##
bayes_mat <- bayes(mat_test)
test_bayes <- function() {
    checkEquals(sum(bayes_mat), 6)
    checkEquals(rownames(bayes_mat), colnames(bayes_mat))
    checkEquals(rownames(bayes_mat), rownames(mat_test))
    checkEquals(ncol(bayes_mat), nrow(bayes_mat))
    checkEquals(nrow(bayes_mat), nrow(mat_test))
    checkTrue(is.numeric(bayes_mat))
    checkTrue(is.matrix(bayes_mat))
    checkTrue(all(bayes_mat %in% c(0, 1)))
}
## END unit test bayes ##

## START unit test addToList ##
l <- list()
l <- MetNet:::addToList(l, "newEntry", matrix())
test_addToList <- function() {
    checkException(MetNet:::addToList(l, "newEntry", NULL))
    checkException(MetNet:::addToList(l, NULL, matrix()))
    checkException(MetNet:::addToList(NULL, "newEntry", matrix()))
    checkEquals(length(l), 1)
    checkTrue(is.matrix(l[[1]]))
    checkEquals(l[[1]], matrix())
}
## END unit test addToList ##

## START unit test threeDotsCall ##
test_threeDots_call <- function() {
    checkException(MetNet:::threeDotsCall("mean", x = 1:10, x = 1:10))
    checkEquals(MetNet:::threeDotsCall("mean", x = 1:10, foo = 1), mean(1:10))
}
## END unit test threeDotsCall ##

## START unit test createStatisticalAdjacencyList ##
stat_adj_l <- createStatisticalAdjacencyList(mat_test,
    model = c("clr", "aracne", "pearson", "spearman", "bayes"))
test_createStatisticalAdjacencyList <- function() {
    checkException(createStatisticalAdjacencyList(NULL, model = "lasso"))
    checkException(createStatisticalAdjacencyList(mat_test, model = "foo"))
    checkException(createStatisticalAdjacencyList(mat_test, model = c("lasso")))
    checkEquals(length(stat_adj_l), 5)
    checkEquals(as.numeric(lapply(stat_adj_l, nrow)), rep(7, 5))
    checkEquals(as.numeric(lapply(stat_adj_l, ncol)), rep(7, 5))
    checkEquals(as.character(unlist((lapply(stat_adj_l, rownames)))),
                        rep(c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), 5))
    checkEquals(as.character(unlist((lapply(stat_adj_l, colnames)))),
                rep(c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), 5))
    checkEquals(names(stat_adj_l), 
        c("clr", "aracne", "pearson", "spearman", "bayes"))
    checkTrue(all(unlist(lapply(stat_adj_l, function(x) is.numeric(x)))))
}
## END unit test createStatisticalAdjacencyList ##

## START unit test consensusAdjacency ##
stat_adj_cons <- consensusAdjacency(stat_adj_l)
list_foo <- list(matrix(1:49, ncol = 7), matrix(1:49, ncol = 7))
rownames(list_foo[[1]]) <- rownames(stat_adj_l[[1]][])
colnames(list_foo[[1]]) <- rownames(stat_adj_l[[1]][])
test_consensusAdjacency <- function() {
    checkException(consensusAdjacency(NULL))
    checkException(consensusAdjacency(stat_adj_l, threshold = "a"))
    checkException(consensusAdjacency(stat_adj_l, method = "foo"))
    checkException(consensusAdjacency(
        list(matrix(1:25, ncol = 5), matrix(1:36, ncol = 6))))
    checkException(consensusAdjacency(list_foo))
    checkEquals(ncol(stat_adj_cons), 7)
    checkEquals(nrow(stat_adj_cons), 7)
    checkEquals(rownames(stat_adj_cons), 
        c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkEquals(colnames(stat_adj_cons), 
        c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkTrue(is.numeric(stat_adj_cons))
    checkTrue(is.matrix(stat_adj_cons))
}
## END unit test consensusAdjacency ##

## START unit test createStatisticalAdjacency ##
stat_adj_l <- createStatisticalAdjacencyList(mat_test, 
        model = c("clr", "aracne","pearson", "spearman", "bayes"))
stat_adj_cons <- consensusAdjacency(stat_adj_l)
stat_adj <- createStatisticalAdjacency(mat_test, 
        model = c("clr", "aracne","pearson", "spearman", "bayes"))
test_createStatisticalAdjacency <- function() {
    checkEquals(stat_adj, stat_adj_cons)
}
## END unit test createStatisticalAdjacency
