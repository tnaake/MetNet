## create toy example data set
data("mat_test", package = "MetNet")

## START unit test lasso ##
lasso_mat <- lasso(x = t(mat_test_z[, 1:5]), parallel = TRUE, PFER = 0.75, cutoff = 0.95)

test_that("lasso", {
    expect_error(lasso(NULL, PFER = 0.75, cutoff = 0.95),
        "must be coercible to non-negative integer")

    expect_error(lasso(mat_test), "Two of the three argumnets")
    expect_equal(rownames(lasso_mat), colnames(lasso_mat))
    expect_equal(rownames(lasso_mat), rownames(mat_test)[1:5])
    expect_equal(ncol(lasso_mat), nrow(lasso_mat))
    expect_equal(nrow(lasso_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(lasso_mat))
    expect_true(is.matrix(lasso_mat))
    expect_true(max(lasso_mat) <= 1)
    expect_true(min(lasso_mat) >= 0)
})
## END unit test lasso ##

## START unit test randomForest ##
rf_mat <- randomForest(mat_test[1:5, ])

test_that("randomForest", {
    expect_error(randomForest(NULL), "find an inherited method for")
    expect_equal(sum(rf_mat), 5, tolerance = 1e-06)
    expect_equal(rownames(rf_mat), colnames(rf_mat))
    expect_equal(rownames(rf_mat), rownames(mat_test[1:5, ]))
    expect_equal(ncol(rf_mat), nrow(rf_mat))
    expect_equal(nrow(rf_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(rf_mat))
    expect_true(is.matrix(rf_mat))
    expect_true(max(rf_mat) <= 1)
    expect_true(min(rf_mat) >= 0)
})
## END unit test randomForest ##

## START unit test clr ##
mi_mat_test_z <- mpmi::cmi(mat_test_z[, 1:5])$bcmi
rownames(mi_mat_test_z) <- colnames(mi_mat_test_z) <- colnames(mat_test_z)[1:5]
clr_mat <- clr(mi_mat_test_z)

test_that("clr", {
    expect_error(clr(NULL), msg = "mi must be a matrix")
    expect_equal(sum(clr_mat), 3.586808, tolerance = 1e-03)
    expect_equal(rownames(clr_mat), colnames(clr_mat))
    expect_equal(rownames(clr_mat), rownames(mat_test)[1:5])
    expect_equal(ncol(clr_mat), nrow(clr_mat))
    expect_equal(nrow(clr_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(clr_mat))
    expect_true(is.matrix(clr_mat))
    expect_true(max(clr_mat) <= 1)
    expect_true(min(clr_mat) >= 0)
})
## END unit test clr ##

## START unit test aracne ##
aracne_mat <- aracne(mi_mat_test_z)

test_that("aracne", {
    expect_error(aracne(NULL, eps = 0.05), "mi must be a matrix")
    suppressWarnings(
        expect_error(aracne(mi_mat_test_z, eps = "a"), 
            "in foreign function call"))
    expect_equal(sum(aracne_mat), 12.69407, tolerance = 1e-02)
    expect_equal(rownames(aracne_mat), colnames(aracne_mat))
    expect_equal(rownames(aracne_mat), rownames(mat_test)[1:5])
    expect_equal(ncol(aracne_mat), nrow(aracne_mat))
    expect_equal(nrow(aracne_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(aracne_mat))
    expect_true(is.matrix(aracne_mat))
    expect_true(max(aracne_mat) <= 1.3)
    expect_true(min(aracne_mat) >= 0)
})
## END unit test aracne ##

## START unit test partialCorrelation/correlation ##
correlation_p_mat <- correlation(mat_test[1:5, ], method = "pearson")
correlation_p_p_mat <- correlation(mat_test[1:5, ], method = "pearson_partial")
correlation_s_mat <- correlation(mat_test[1:5, ], method = "spearman")
suppressWarnings(
    correlation_s_p_mat <- correlation(mat_test[1:5, ], 
                                                method = "spearman_partial"))
correlation_g_mat <- correlation(mat_test[1:5, ], method = "ggm")

test_that("partialCorrelation", {
    pcor_p <- partialCorrelation(t(mat_test), method = "pearson")
    expect_equal(names(pcor_p), c("r", "p"))
    expect_equal(pcor_p$r, psych::partial.r(t(mat_test), method = "pearson"))
    expect_equal(sum(pcor_p$r), 9.60788, tolerance = 1e-06)
    expect_equal(sum(pcor_p$p, na.rm = TRUE), 26.78391, tolerance = 1e-06)
    suppressWarnings(
        pcor_s <- partialCorrelation(t(mat_test), method = "spearman"))
    expect_equal(names(pcor_s), c("r", "p"))
    suppressWarnings(
        expect_equal(pcor_s$r, psych::partial.r(t(mat_test), 
                                                method = "spearman")))
    expect_equal(sum(pcor_s$r), 7.211808, tolerance = 1e-06)
    expect_equal(sum(pcor_s$p, na.rm = TRUE), 43.12084, tolerance = 1e-06)
})

test_that("correlation", {
    expect_error(correlation(NULL, method = "pearson"),
        "argument is not a matrix")
    expect_error(correlation(mat_test, method = "a"),
        msg = "object [']cor_mat['] not found")

    ## pearson
    expect_equal(names(correlation_p_mat), c("r", "p"))
    expect_equal(correlation_p_mat$r,
        cor(t(mat_test[1:5, ]), method = "pearson"), tolerance = 1e-06)
    expect_equal(sum(correlation_p_mat$r), 3.27161, tolerance = 1e-06)
    expect_equal(sum(correlation_p_mat$p, na.rm = TRUE), 0.2890875, tolerance = 1e-06)
    expect_equal(rownames(correlation_p_mat$r), colnames(correlation_p_mat$r))
    expect_equal(rownames(correlation_p_mat$r), rownames(mat_test)[1:5])
    expect_equal(ncol(correlation_p_mat$r), nrow(correlation_p_mat$r))
    expect_equal(nrow(correlation_p_mat$r), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_p_mat$r))
    expect_true(is.matrix(correlation_p_mat$r))
    expect_true(max(correlation_p_mat$r) <= 1)
    expect_true(min(correlation_p_mat$r) >= -1)

    ## spearman
    expect_equal(names(correlation_s_mat), c("r", "p"))
    expect_equal(correlation_s_mat$r,
        cor(t(mat_test[1:5, ]), method = "spearman"), tolerance = 1e-06)
    expect_equal(sum(correlation_s_mat$r), 3.153383, tolerance = 1e-06)
    expect_equal(sum(correlation_s_mat$p, na.rm = TRUE), 0.3236611, tolerance = 1e-06)
    expect_equal(rownames(correlation_s_mat$r), colnames(correlation_s_mat$r))
    expect_equal(rownames(correlation_s_mat$r), rownames(mat_test)[1:5])
    expect_equal(ncol(correlation_s_mat$r), nrow(correlation_s_mat$r))
    expect_equal(nrow(correlation_s_mat$r), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_s_mat$r))
    expect_true(is.matrix(correlation_s_mat$r))
    expect_true(max(correlation_s_mat$r) <= 1)
    expect_true(min(correlation_s_mat$r) >= -1)

    ## partial pearson
    expect_equal(names(correlation_p_p_mat), c("r", "p"))
    expect_true(all(correlation_p_p_mat$r -
        ppcor::pcor(t(mat_test[1:5, ]), method = "pearson")$r == 0))
    expect_equal(sum(correlation_p_p_mat$r), 7.053181, tolerance = 1e-06)
    expect_equal(sum(correlation_p_p_mat$p), 12.44441, tolerance = 1e-06)
    expect_equal(rownames(correlation_p_p_mat$r), 
        colnames(correlation_p_p_mat$r))
    expect_equal(rownames(correlation_p_p_mat$r), rownames(mat_test[1:5, ]))
    expect_equal(ncol(correlation_p_p_mat$r), 
        nrow(correlation_p_p_mat$r))
    expect_equal(nrow(correlation_p_p_mat$r), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_p_p_mat$r))
    expect_true(is.matrix(correlation_p_p_mat$r))
    expect_true(max(correlation_p_p_mat$r) <= 1)
    expect_true(min(correlation_p_p_mat$r) >= -1)

    ## partial spearman
    expect_equal(names(correlation_s_p_mat), c("r", "p"))
    suppressWarnings(expect_true(all(correlation_s_p_mat$r - 
        ppcor::pcor(t(mat_test[1:5, ]), method = "spearman")$r == 0)))
    expect_equal(sum(correlation_s_p_mat$r), 4.479568, tolerance = 1e-06)
    expect_equal(sum(correlation_s_p_mat$p, na.rm = TRUE), 19.09195, 
        tolerance = 1e-06)
    expect_equal(rownames(correlation_s_p_mat$r), 
        colnames(correlation_s_p_mat$r))
    expect_equal(rownames(correlation_s_p_mat$r), 
        rownames(mat_test)[1:5])
    expect_equal(ncol(correlation_s_p_mat$r), 
        nrow(correlation_s_p_mat$r))
    expect_equal(nrow(correlation_s_p_mat$r), 
        nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_s_p_mat$r))
    expect_true(is.matrix(correlation_s_p_mat$r))
    expect_true(max(correlation_s_p_mat$r) <= 1.0000001)
    expect_true(min(correlation_s_p_mat$r) >= -1.0000001)
    
    ## ggm
    expect_equal(names(correlation_g_mat), c("r", "p"))
    expect_true(all(correlation_g_mat$r -
        GeneNet::ggm.estimate.pcor(t(mat_test[1:5, ]), 
            method = "static")[seq_len(dim(mat_test[1:5, ])[1]), 
                seq_len(dim(mat_test[1:5, ])[1])]== 0))
    expect_equal(sum(correlation_g_mat$r), 5.421633, tolerance = 1e-06)
    expect_equal(sum(correlation_g_mat$p), 2.889295, tolerance = 1e-06)
    expect_equal(rownames(correlation_g_mat$r), 
                 colnames(correlation_g_mat$r))
    expect_equal(rownames(correlation_g_mat$r), rownames(mat_test[1:5, ]))
    expect_equal(ncol(correlation_g_mat$r), 
                 nrow(correlation_g_mat$r))
    expect_equal(nrow(correlation_g_mat$r), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_g_mat$r))
    expect_true(is.matrix(correlation_g_mat$r))
    expect_true(max(correlation_g_mat$r) <= 1)
    expect_true(min(correlation_g_mat$r) >= -1)
})
## END unit test correlation ##

## START unit test bayes ##
bayes_mat <- bayes(mat_test[1:5, ])

test_that("bayes", {
    expect_true(sum(bayes_mat) > 9 & sum(bayes_mat) < 13)
    expect_equal(rownames(bayes_mat), colnames(bayes_mat))
    expect_equal(rownames(bayes_mat), rownames(mat_test)[1:5])
    expect_equal(ncol(bayes_mat), nrow(bayes_mat))
    expect_equal(nrow(bayes_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(bayes_mat))
    expect_true(is.matrix(bayes_mat))
    expect_true(max(bayes_mat) <= 1)
    expect_true(min(bayes_mat) >= 0)
})
## END unit test bayes ##

## START unit test addToList ##
l <- list()
l <- MetNet:::addToList(l, "newEntry", matrix())

test_that("addToList", {
    expect_error(MetNet:::addToList(l, "newEntry", NULL),
       "not a matrix")
    expect_error(MetNet:::addToList(l, NULL, matrix()),
        "not a character")
    expect_error(MetNet:::addToList(NULL, "newEntry", matrix()),
        "not a list")
    expect_equal(length(l), 1)
    expect_true(is.matrix(l[[1]]))
    expect_equal(l[[1]], matrix())
})
## END unit test addToList ##


## START unit test statistical ##
stat_adj <- statistical(mat_test[1:5, ],
    model = c("randomForest", "clr", "aracne", "pearson",
        "pearson_partial", "pearson_semipartial", "spearman",
        "spearman_partial", "spearman_semipartial", "ggm", "bayes"),
    PFER = 0.75, cutoff = 0.95)

test_that("statistical", {
    expect_error(statistical(NULL, model = "lasso"), "not a numerical matrix")
    expect_error(statistical(mat_test, model = "foo"),
        "'model' not implemented in statistical")
    expect_error(statistical(mat_test, model = c("lasso")),
        "Two of the three argumnets")

    ## take a high tolerance value for randomForest and bayes
    ## since these models are probabilistic
    tmp <- rf_mat; diag(tmp) <- NaN
    expect_equal(assay(stat_adj, "randomForest_coef"), tmp, tolerance = 5e-01)
    tmp <- bayes_mat; diag(tmp) <- NaN
    expect_equal(assay(stat_adj, "bayes_coef"), tmp, tolerance = 5e-01)
    tmp <- clr_mat; diag(tmp) <- NaN
    expect_true(all(assay(stat_adj, "clr_coef") == tmp, na.rm = TRUE))
    tmp <- aracne_mat; diag(tmp) <- NaN
    expect_true(all(assay(stat_adj, "aracne_coef") == tmp, na.rm = TRUE))
    tmp <-  correlation(mat_test[1:5, ], method = "pearson")
    diag(tmp$r) <- NaN; diag(tmp$p) <- NaN
    suppressWarnings(
        expect_true(all(assay(stat_adj, "pearson_coef"), tmp$r, na.rm = TRUE)))
    suppressWarnings(
        expect_true(all(assay(stat_adj, "pearson_pvalue"), tmp$p, na.rm = TRUE)))
    tmp <- correlation(mat_test[1:5, ], method = "pearson_partial")
    diag(tmp$r) <- NaN; diag(tmp$p) <- NaN
    expect_true(
        all(assay(stat_adj, "pearson_partial_coef") == tmp$r, 
        na.rm = TRUE))
    expect_true(
        all(assay(stat_adj, "pearson_partial_pvalue") == tmp$p,
        na.rm = TRUE))
    tmp <- correlation(mat_test[1:5, ], method = "spearman")
    diag(tmp$r) <- NaN; diag(tmp$p) <- NaN
    expect_true(all(stat_adj[["spearman_coef"]] == tmp$r, na.rm = TRUE))
    expect_true(all(stat_adj[["spearman_pvalue"]] == tmp$p, na.rm = TRUE))
    suppressWarnings(tmp <- correlation(mat_test[1:5, ], method = "spearman_partial"))
    diag(tmp$r) <- NaN; diag(tmp$p) <- NaN
    expect_true(all(stat_adj[["spearman_partial_coef"]] == tmp$estimate, na.rm = TRUE))
    expect_true(all(stat_adj[["spearman_partial_pvalue"]] == tmp$p.value, na.rm = TRUE))
    tmp <- correlation(mat_test[1:5, ], method = "ggm")
    diag(tmp$r) <- NaN; diag(tmp$p) <- NaN
    expect_true(all(stat_adj[["ggm_coef"]] == tmp$estimate, na.rm = TRUE))
    expect_true(all(stat_adj[["ggm_pvalue"]] == tmp$p, na.rm = TRUE))
    expect_equal(length(assays(stat_adj)), 14) 
    expect_equal(as.numeric(lapply(assays(stat_adj), nrow)), rep(5, 14))
    expect_equal(as.numeric(lapply(assays(stat_adj), ncol)), rep(5, 14))
    expect_equal(as.character(unlist((lapply(assays(stat_adj), rownames)))),
        rep(c("x1", "x2", "x3", "x4", "x5"), 14))
    expect_equal(as.character(unlist((lapply(assays(stat_adj), colnames)))),
        rep(c("x1", "x2", "x3", "x4", "x5"), 14))
    expect_equal(assayNames(stat_adj),
        c("randomForest_coef", "clr_coef", "aracne_coef", "pearson_coef",
            "pearson_pvalue", "pearson_partial_coef", "pearson_partial_pvalue", 
            "spearman_coef", "spearman_pvalue", "spearman_partial_coef", 
            "spearman_partial_pvalue", "ggm_coef", "ggm_pvalue","bayes_coef"))
    expect_true(all(unlist(lapply(assays(stat_adj), function(x) is.numeric(x)))))
    expect_equal(length(stat_adj), 5)
    expect_equal(dim(stat_adj), c(5, 5))
    expect_equal(directed(stat_adj), TRUE)
    expect_equal(type(stat_adj), "statistical")
    expect_equal(thresholded(stat_adj), FALSE)
    expect_equal(directed(statistical(mat_test, model = "pearson")), FALSE)
    expect_equal(directed(statistical(mat_test, model = "spearman")), FALSE)
    expect_equal(directed(statistical(mat_test, model = "ggm")), FALSE)
    expect_equal(directed(statistical(mat_test, model = "randomForest")), TRUE)
})
## END unit test statistical ##


## START unit test getLinks ##
mat <- matrix(0:8, ncol = 3, nrow = 3)
getLinks_df <- MetNet:::getLinks(mat, exclude = "== 0")

test_that("getLinks", {
    expect_error(MetNet:::getLinks(NULL), "argument is of length zero")
    expect_error(MetNet:::getLinks(mat[, 1:2]), "not a square matrix")
    expect_error(MetNet:::getLinks(mat, exclude = "foo"),
        "object 'matfoo' not found")

    ## checks for exclude = "== 0"
    expect_true(is.data.frame(getLinks_df))
    expect_equal(getLinks_df$row, rep(c(1, 2, 3), 3))
    expect_equal(getLinks_df$col, rep(c(1, 2, 3), each = 3))
    expect_equal(getLinks_df$confidence, c(NaN, 1:8))
    expect_equal(getLinks_df$rank, c(NaN, 8:1))

    ## checks for exclude = NULL
    getLinks_df <- MetNet:::getLinks(mat, exclude = NULL)
    expect_equal(getLinks_df$row, rep(c(1, 2, 3), 3))
    expect_equal(getLinks_df$col, rep(c(1, 2, 3), each = 3))
    expect_equal(getLinks_df$confidence, 0:8)
    expect_equal(getLinks_df$rank, 9:1)

})
## END unit test getLinks ##


## START unit test threshold  ##
## remove partial/semipartial correlation from stat_adj_l
stat_adj_cut <- statistical(mat_test[1:5, ], 
    model = c("clr", "aracne", "pearson", "spearman", "ggm"))

args_thr <- list(filter = "clr_coef > 0.5 & aracne_coef > 0.8 & abs(pearson_coef) > 0.95 & abs(spearman_coef) > 0.95 & abs(ggm_coef) > 0.2")
thr_thr <- threshold(stat_adj_cut, type = "threshold", args = args_thr)

args_top <- list(n = 2)
thr_top1 <- threshold(stat_adj_cut, type = "top1", args = args_top)
thr_top2 <- threshold(stat_adj_cut, type = "top2", args = args_top)
thr_mean <- threshold(stat_adj_cut, type = "mean", args = args_top)

test_that("threshold", {

    ## test arguments
    expect_error(threshold(NULL, type = "threshold", args = args_thr),
        "is not an 'AdjacencyMatrix' object")
    expect_error(threshold(1:3, type = "threshold", args = args_thr),
        "is not an 'AdjacencyMatrix' object")
    args_thr_double <- list(
        filter = "randomForest_coef > 0.2 & clr_coef > 0.5 & foo > 0.5")
    expect_error(
        threshold(stat_adj, type = "threshold", args = args_thr_double),
        "not found")
    expect_error(
        threshold(stat_adj_cut, type = "foo", args = args_thr),
        "should be one of ")
    expect_error(threshold(stat_adj_cut, type = "top1", args = args_top, 
        values = "foo"), "should be one of ")

    ## check args for top1, top2, mean
    expect_error(threshold(stat_adj, type = "top1", args = list(x = 1)),
        "does not contain the numeric entry `n` of length 1")
    expect_error(threshold(stat_adj, type = "top2", args = list(x = 1)),
        "does not contain the numeric entry `n` of length 1")
    expect_error(threshold(stat_adj, type = "mean", args = list(x = 1)),
        "does not contain the numeric entry `n` of length 1")
    expect_warning(threshold(stat_adj, type = "top1", args = list(n = 1:2)),
        "numerical expression has 2 elements")
    expect_warning(threshold(stat_adj, type = "top2", args = list(n = 1:2)),
        "numerical expression has 2 elements")
    expect_warning(threshold(stat_adj, type = "mean", args = list(n = 1:2)),
        "numerical expression has 2 elements")
    suppressWarnings(expect_error(threshold(stat_adj, type = "top1", 
        args = list(n = "a")), "NA/NaN argument"))
    suppressWarnings(expect_error(threshold(stat_adj, type = "top2", 
        args = list(n = "a")), "NA/NaN argument"))
    suppressWarnings(expect_error(threshold(stat_adj, type = "mean", 
        args = list(n = "a")), "NA/NaN argument"))

    ## check output
    expect_true(is(thr_thr, "AdjacencyMatrix"))
    expect_true(is(thr_top1, "AdjacencyMatrix"))
    expect_true(is(thr_top2, "AdjacencyMatrix"))
    expect_true(is(thr_mean, "AdjacencyMatrix"))
    expect_true(is(thr_thr, "AdjacencyMatrix"))
    assay_names <- c("clr_coef", "aracne_coef", "pearson_coef",
        "pearson_pvalue", "spearman_coef", "spearman_pvalue",
        "ggm_coef", "ggm_pvalue", "consensus")
    expect_equal(assayNames(thr_top1), assay_names)
    expect_equal(assayNames(thr_top2), assay_names)
    expect_equal(assayNames(thr_mean), assay_names)
    expect_true(is.numeric(assay(thr_top1, "consensus")))
    expect_true(is.numeric(assay(thr_top2, "consensus")))
    expect_true(is.numeric(assay(thr_mean, "consensus")))
    expect_equal(rownames(assay(thr_thr, "consensus")), 
        c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(colnames(assay(thr_thr, "consensus")), 
        c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(rownames(assay(thr_top1, "consensus")), 
        c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(colnames(assay(thr_top1, "consensus")), 
        c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(rownames(assay(thr_top2, "consensus")), 
        c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(colnames(assay(thr_top2, "consensus")), 
        c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(rownames(assay(thr_mean, "consensus")), 
        c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(colnames(assay(thr_mean, "consensus")), 
        c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(sum(assay(thr_thr, "consensus"), na.rm = TRUE), 1)
    expect_equal(sum(assay(thr_top1, "consensus"), na.rm = TRUE), 4)
    expect_equal(sum(assay(thr_top2, "consensus"), na.rm = TRUE), 4)
    expect_equal(sum(assay(thr_mean, "consensus"), na.rm = TRUE), 2)
    expect_true(all(assay(thr_thr, "consensus") %in% c(0, 1, NaN)))
    expect_true(all(assay(thr_top1, "consensus") %in% c(0, 1, NaN)))
    expect_true(all(assay(thr_top2, "consensus") %in% c(0, 1, NaN)))
    expect_true(all(assay(thr_mean, "consensus") %in% c(0, 1, NaN)))
    
    stat_adj_cut <- statistical(mat_test[1:5, ], model = "bayes")
    am_all <- threshold(stat_adj_cut, type = "top1", args = args_top, 
        values = "all")
    am_min <- threshold(stat_adj_cut, type = "top1", args = args_top, 
        values = "min")
    am_max <- threshold(stat_adj_cut, type = "top1", args = args_top, 
        values = "max")
    expect_equal(directed(am_all), TRUE)
    expect_equal(directed(am_min), FALSE)
    expect_equal(directed(am_max), FALSE)
    expect_equal(as.vector(assay(am_all, "consensus")[, 2]), c(1, NaN, 0, 0, 0))
    expect_equal(as.vector(assay(am_min, "consensus")[, 2]), c(1, NaN, 0, 0, 0))
    expect_equal(as.vector(assay(am_max, "consensus")[, 2]), c(1, NaN, 0, 0, 0))
    expect_equal(as.vector(assay(am_all, "consensus")[2, ]), c(1, NaN, 0, 0, 0))
    expect_equal(as.vector(assay(am_min, "consensus")[2, ]), c(1, NaN, 0, 0, 0))
    expect_equal(as.vector(assay(am_max, "consensus")[2, ]), c(1, NaN, 0, 0, 0))
    expect_true(thresholded(am_all))
    expect_true(thresholded(am_min))
    expect_true(thresholded(am_max))
    expect_equal(type(am_all), "statistical")
    expect_equal(type(am_min), "statistical")
    expect_equal(type(am_max), "statistical")
    
    ## test NA values
    mat_test_NA <- mat_test[1:5, ]
    mat_test_NA[1, 5] <- NA
    stat_adj_cut_NA <- statistical(mat_test_NA, model = c("pearson", "ggm"))
    thr_NA_thr <- threshold(stat_adj_cut_NA, type = "threshold", 
        args = list(filter = "ggm_coef > 0.01"), na.rm = TRUE)
    expect_equal(as.vector(assay(thr_NA_thr, "consensus")[, 1]), 
        c(NaN, 0, 0, 0, 0))
    thr_NA_thr <- threshold(stat_adj_cut_NA, type = "threshold", 
        args = list(filter = "ggm_coef > 0.01"), na.rm = FALSE)
    expect_equal(as.vector(assay(thr_NA_thr, "consensus")[, 1]), 
        c(NaN, NaN, NaN, NaN, NaN))
    thr_NA_thr <- threshold(stat_adj_cut_NA, type = "threshold", 
        args = list(filter = "ggm_coef > 0.01 | is.na(ggm_coef)"), 
        na.rm = FALSE)
    expect_equal(as.vector(assay(thr_NA_thr, "consensus")[, 1]), 
        c(NaN, 1, 1, 1, 1))
    thr_NA_top1 <- threshold(stat_adj_cut_NA, type = "top1", 
        args = list(n = 5), na.rm = TRUE)
    expect_equal(as.vector(assay(thr_NA_top1, "consensus")[, 1]), 
        c(NaN, 1, 1, 0, 0))
    thr_NA_top1 <- threshold(stat_adj_cut_NA, type = "top1", 
        args = list(n = 5), na.rm = FALSE)
    expect_equal(as.vector(assay(thr_NA_top1, "consensus")[, 1]), ##########
        c(NaN, NaN, NaN, NaN, NaN))
    thr_NA_top2 <- threshold(stat_adj_cut_NA, type = "top2", 
        args = list(n = 5), na.rm = TRUE)
    expect_equal(as.vector(assay(thr_NA_top2, "consensus")[, 1]), 
        c(NaN, 0, 0, 0, 0))
    thr_NA_top2 <- threshold(stat_adj_cut_NA, type = "top2", 
        args = list(n = 5), na.rm = FALSE)
    expect_equal(as.vector(assay(thr_NA_top2, "consensus")[, 1]), ###########
        c(NaN, NaN, NaN, NaN, NaN))
    thr_NA_mean <- threshold(stat_adj_cut_NA, type = "top2", 
        args = list(n = 5), na.rm = TRUE)
    expect_equal(as.vector(assay(thr_NA_mean, "consensus")[, 1]), 
        c(NaN, 0, 0, 0, 0))
    thr_NA_mean <- threshold(stat_adj_cut_NA, type = "top2", 
        args = list(n = 5), na.rm = FALSE)
    expect_equal(as.vector(assay(thr_NA_mean, "consensus")[, 1]), #########
        c(NaN, NaN, NaN, NaN, NaN))
})
## END unit test threshold  ##

## START unit test topKnet ##
ranks <- matrix(c(c(1, 2, 3), c(2, 1, 3)), ncol = 2)

test_that("topKnet", {

    ## type of ranks
    ## na.rm = TRUE
    ranks <- matrix(c(c(1, 2, 3), c(2, 1, 3)), ncol = 2)
    expect_error(MetNet:::topKnet(ranks = 1:3, na.rm = TRUE), "ranks is not a numerical")
    expect_error(MetNet:::topKnet(ranks = data.frame(x = 1:3, y = 3:1),
        type = "top1", na.rm = TRUE), "ranks is not a numeric")
    expect_error(MetNet:::topKnet(ranks = matrix(c("a", "b", "c")),
        type = "top1", na.rm = TRUE), "ranks is not a numeric")

    ## type argument
    expect_error(MetNet:::topKnet(ranks = ranks, type = "foo", na.rm = TRUE), 
        "type neither")

    ## check results
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "top1", na.rm = TRUE), 
        c(1, 1, 3))
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "top2", na.rm = TRUE), 
        c(2, 2, 3))
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "mean", na.rm = TRUE), 
        c(1.5, 1.5, 3))

    ## matrix with ncol 1
    expect_equal(MetNet:::topKnet(ranks = matrix(1:3, nrow = 3),
        type = "top1", na.rm = TRUE), 1:3)
    expect_error(MetNet:::topKnet(ranks = matrix(1:3, nrow = 3),
        type = "top2", na.rm = TRUE), "ncol[(]ranks[])] has to be")
    expect_equal(MetNet:::topKnet(ranks = matrix(1:3, nrow = 3),
        type = "mean", na.rm = TRUE), 1:3)
    
    ## na.rm = FALSE
    ranks <- matrix(c(c(1, 2, 3), c(2, 1, 3)), ncol = 2)
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "top1", na.rm = FALSE), 
        c(1, 1, 3))
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "top2", na.rm = FALSE), 
        c(2, 2, 3))
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "mean", na.rm = FALSE), 
        c(1.5, 1.5, 3))
    ranks <- matrix(c(c(1, 2, 3), c(2, 1, NA)), ncol = 2)
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "top1", na.rm = FALSE), 
        c(1, 1, NA))
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "top2", na.rm = FALSE), 
        c(2, 2, NA))
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "mean", na.rm = FALSE), 
        c(1.5, 1.5, NA))
})
## END unit test topKnet ##


## START unit test threeDotsCall ##

test_that("threeDotsCall", {
    expect_error(MetNet:::threeDotsCall("mean", x = 1:10, x = 1:10),
        "duplicated args in ...")
    expect_equal(MetNet:::threeDotsCall("mean", x = 1:10, foo = 1), mean(1:10))
})
## END unit test threeDotsCall ##

