## create toy example data set
set.seed(1)
mat_test <- matrix(rnorm(140, mean = 10, sd = 2), nrow = 7)
mat_test[1:3, ] <- t(apply(mat_test[1:3, ], 1, sort))
mat_test[5:7, ] <- t(apply(mat_test[5:7, ], 1, sort, decreasing = TRUE))
rownames(mat_test) <- paste("x", 1:7, sep = "")
mat_test_z <- apply(mat_test, 1, function(x) (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))

## START unit test lasso ##
lasso_mat <- lasso(t(mat_test_z), parallel = FALSE, PFER = 0.75, cutoff = 0.95)
test_lasso <- function() {
    checkException(lasso(mat_test))
    checkEquals(sum(lasso_mat), 6)
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
randomForest_mat <- randomForest(mat_test)
test_randomForest <- function() {
    checkEquals(sum(randomForest_mat), 6)
    checkEquals(rownames(randomForest_mat), colnames(randomForest_mat))
    checkEquals(rownames(randomForest_mat), rownames(mat_test))
    checkEquals(ncol(randomForest_mat), nrow(randomForest_mat))
    checkEquals(nrow(randomForest_mat), nrow(mat_test))
    checkTrue(is.numeric(randomForest_mat))
    checkTrue(is.matrix(randomForest_mat))
    checkTrue(all(randomForest_mat %in% c(0, 1)))
}
## END unit test randomForest ## 

## START unit test clr ## 
mi_mat_test_z <- mpmi::cmi(mat_test_z)$bcmi
rownames(mi_mat_test_z) <- colnames(mi_mat_test_z) <- colnames(mat_test_z)
clr_mat <- clr(mi_mat_test_z)
test_clr <- function() {
    checkException(clr(mi_mat_test_z, clr_threshold = "a"))
    checkEquals(sum(clr_mat), 10)
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
    checkEquals(sum(aracne_mat), 27)
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
correlation_s_mat <- correlation(mat_test, type = "spearman")
test_correlation <- function() {
    checkException(correlation(mat_test, correlation_threshold = "a"))
    checkEquals(sum(correlation_p_mat), 37)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "bonferroni", correlation_threshold = 1e-11)), 21)
    checkEquals(sum(correlation(mat_test, correlation_adjust = "none", correlation_threshold = 1e-11)), 23)
    checkEquals(sum(correlation_s_mat), 37)
    checkEquals(rownames(correlation_p_mat), colnames(correlation_p_mat))
    checkEquals(rownames(correlation_s_mat), colnames(correlation_s_mat))
    checkEquals(rownames(correlation_p_mat), rownames(mat_test))
    checkEquals(rownames(correlation_s_mat), rownames(mat_test))
    checkEquals(ncol(correlation_p_mat), nrow(correlation_p_mat))
    checkEquals(ncol(correlation_s_mat), nrow(correlation_s_mat))
    checkEquals(nrow(correlation_p_mat), nrow(mat_test))
    checkEquals(nrow(correlation_s_mat), nrow(mat_test))
    checkTrue(is.numeric(correlation_p_mat))
    checkTrue(is.numeric(correlation_s_mat))
    checkTrue(is.matrix(correlation_p_mat))
    checkTrue(is.matrix(correlation_s_mat))
    checkTrue(all(correlation_p_mat %in% c(0, 1)))
    checkTrue(all(correlation_s_mat %in% c(0, 1)))
}
## END unit test correlation ##

## START unit test bayes ##
bayes_mat <- bayes(mat_test)
test_bayes <- function() {
    checkEquals(sum(bayes_mat), 8)
    checkEquals(rownames(bayes_mat), colnames(bayes_mat))
    checkEquals(rownames(bayes_mat), rownames(mat_test))
    checkEquals(ncol(bayes_mat), nrow(bayes_mat))
    checkEquals(nrow(bayes_mat), nrow(mat_test))
    checkTrue(is.numeric(bayes_mat))
    checkTrue(is.matrix(bayes_mat))
    checkTrue(all(bayes_mat %in% c(0, 1)))
}
## END unit test bayes ## 

## START unit test add_to_list ## 
l <- list()
l <- add_to_list(l, "newEntry", matrix())
test_add_to_list <- function() {
    checkException(add_to_list(l, "newEntry", NULL))
    checkException(add_to_list(l, NULL, matrix()))
    checkException(add_to_list(NULL, "newEntry", matrix()))
    checkEquals(length(l), 1)
    checkTrue(is.matrix(l[[1]]))
    checkEquals(l[[1]], matrix())
}
## END unit test add_to_list ##

## START unit test threeDots_call ##
test_threeDots_call <- function() {
    checkException(threeDots_call("mean", x = 1:10, x = 1:10))
    checkEquals(threeDots_call("mean", x = 1:10, foo = 1), mean(1:10))
}
## END unit test threeDots_call ## 

## START unit test create_statistical_networks_list ##
stat_net_l <- create_statistical_networks_list(mat_test, 
    model = c("lasso", "randomForest", "clr", "aracne","pearson", 
              "spearman", "bayes"), PFER = 0.75, cutoff = 0.95)
test_create_statistical_networks_list <- function() {
    checkException(create_statistical_networks_list(NULL, model = "lasso"))
    checkException(create_statistical_networks_list(mat_test, model = "foo"))
    checkException(create_statistical_networks_list(mat_test, model = c("lasso")))
    checkEquals(length(stat_net_l), 7)
    checkEquals(as.numeric(lapply(stat_net_l, nrow)), rep(7, 7))
    checkEquals(as.numeric(lapply(stat_net_l, ncol)), rep(7, 7))
    checkEquals(as.character(unlist((lapply(stat_net_l, rownames)))),
                        rep(c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), 7))
    checkEquals(as.character(unlist((lapply(stat_net_l, colnames)))),
                rep(c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), 7))
    checkEquals(names(stat_net_l), c("lasso", "randomForest", "clr", "aracne", 
                         "pearson", "spearman", "bayes"))
    checkTrue(all(unlist(lapply(stat_net_l, function(x) is.numeric(x)))))
}
## END unit test create_statistical_networks_list ## 

## START unit test consensus_network ##
stat_net_cons <- consensus_network(stat_net_l)
list_foo <- list(matrix(1:49, ncol = 7), matrix(1:49, ncol = 7))
rownames(list_foo[[1]]) <- rownames(stat_net_l[[1]][])
colnames(list_foo[[1]]) <- rownames(stat_net_l[[1]][])
test_consensus_network <- function() {
    checkException(consensus_network(NULL))
    checkException(consensus_network(stat_net_l, threshold = "a"))
    checkException(consensus_network(stat_net_l, method = "foo"))
    checkException(consensus_network(
        list(matrix(1:25, ncol = 5), matrix(1:36, ncol = 6))))
    checkException(consensus_network(list_foo))
    checkEquals(ncol(stat_net_cons), 7)
    checkEquals(nrow(stat_net_cons), 7)
    checkEquals(rownames(stat_net_cons), 
        c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkEquals(colnames(stat_net_cons), 
        c("x1", "x2", "x3", "x4", "x5", "x6", "x7"))
    checkTrue(is.numeric(stat_net_cons))
    checkTrue(is.matrix(stat_net_cons))
}
## END unit test consensus_network ## 

## START unit test create_statistical_network ## 
stat_net_l <- create_statistical_networks_list(mat_test, 
        model = c("clr", "aracne","pearson", "spearman", "bayes"))
stat_net_cons <- consensus_network(stat_net_l)
stat_net <- create_statistical_network(mat_test, 
        model = c("clr", "aracne","pearson", "spearman", "bayes"))
test_create_statistical_network <- function() {
    checkEquals(stat_net, stat_net_cons)
}
## END unit test create_statistical_network
