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
    expect_equal(sum(clr_mat), 3.586808, tolerance = 1e-06)
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
    expect_error(aracne(mi_mat_test_z, eps = "a"), "in foreign function call")
    expect_equal(sum(aracne_mat), 12.69407, tolerance = 1e-06)
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

## START unit test correlation ##
correlation_p_mat <- correlation(mat_test[1:5, ], type = "pearson")
correlation_p_p_mat <- correlation(mat_test[1:5, ], type = "pearson_partial")
correlation_p_sp_mat <- correlation(mat_test[1:5, ],
    type = "pearson_semipartial")
correlation_s_mat <- correlation(mat_test[1:5, ], type = "spearman")
correlation_s_p_mat <- correlation(mat_test[1:5, ], type = "spearman_partial")
correlation_s_sp_mat <- correlation(mat_test[1:5, ],
    type = "spearman_semipartial")

test_that("correlation", {

    expect_error(correlation(NULL, type = "pearson"),
        "argument is not a matrix")
    expect_error(correlation(mat_test, type = "a"),
        msg = "object [']cor_mat['] not found")

    ## pearson
    expect_equal(correlation_p_mat,
        abs(cor(t(mat_test[1:5, ]), method = "pearson")), tolerance = 1e-06)
    expect_equal(sum(correlation_p_mat), 20.44286, tolerance = 1e-06)
    expect_equal(rownames(correlation_p_mat), colnames(correlation_p_mat))
    expect_equal(rownames(correlation_p_mat), rownames(mat_test)[1:5])
    expect_equal(ncol(correlation_p_mat), nrow(correlation_p_mat))
    expect_equal(nrow(correlation_p_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_p_mat))
    expect_true(is.matrix(correlation_p_mat))
    expect_true(max(correlation_p_mat) <= 1)
    expect_true(min(correlation_p_mat) >= 0)

    ## spearman
    expect_equal(correlation_s_mat,
        abs(cor(t(mat_test[1:5, ]), method = "spearman")), tolerance = 1e-06)
    expect_equal(sum(correlation_s_mat), 20.69323, tolerance = 1e-06)
    expect_equal(rownames(correlation_s_mat), colnames(correlation_s_mat))
    expect_equal(rownames(correlation_s_mat), rownames(mat_test)[1:5])
    expect_equal(ncol(correlation_s_mat), nrow(correlation_s_mat))
    expect_equal(nrow(correlation_s_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_s_mat))
    expect_true(is.matrix(correlation_s_mat))
    expect_true(max(correlation_s_mat) <= 1)
    expect_true(min(correlation_s_mat) >= 0)

    ## partial pearson
    expect_true(all(correlation_p_p_mat -
        abs(ppcor::pcor(t(mat_test[1:5, ]), method = "pearson")$estimate) == 0))
    expect_equal(sum(correlation_p_p_mat), 11.02719, tolerance = 1e-06)
    expect_equal(rownames(correlation_p_p_mat), colnames(correlation_p_p_mat))
    expect_equal(rownames(correlation_p_p_mat), rownames(mat_test[1:5, ]))
    expect_equal(ncol(correlation_p_p_mat), nrow(correlation_p_p_mat))
    expect_equal(nrow(correlation_p_p_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_p_p_mat))
    expect_true(is.matrix(correlation_p_p_mat))
    expect_true(max(correlation_p_p_mat) <= 1)
    expect_true(min(correlation_p_p_mat) >= 0)

    ## semi-partial pearson
    expect_true(all(correlation_p_sp_mat - abs(ppcor::spcor(t(mat_test[1:5, ]),
        method = "pearson")$estimate) == 0))
    expect_equal(sum(correlation_p_sp_mat), 6.690744, tolerance = 1e-06)
    expect_equal(rownames(correlation_p_sp_mat), colnames(correlation_p_sp_mat))
    expect_equal(rownames(correlation_p_sp_mat), rownames(mat_test)[1:5])
    expect_equal(ncol(correlation_p_sp_mat), nrow(correlation_p_sp_mat))
    expect_equal(nrow(correlation_p_sp_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_p_sp_mat))
    expect_true(is.matrix(correlation_p_sp_mat))
    expect_true(max(correlation_p_sp_mat) <= 1)
    expect_true(min(correlation_p_sp_mat) >= 0)

    ## partial spearman
    expect_true(all(correlation_s_p_mat - abs(ppcor::pcor(t(mat_test[1:5, ]),
        method = "spearman")$estimate) == 0))
    expect_equal(sum(correlation_s_p_mat), 20.69323, tolerance = 1e-06)
    expect_equal(rownames(correlation_s_p_mat), colnames(correlation_s_p_mat))
    expect_equal(rownames(correlation_s_p_mat), rownames(mat_test)[1:5])
    expect_equal(ncol(correlation_s_p_mat), nrow(correlation_s_p_mat))
    expect_equal(nrow(correlation_s_p_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_s_p_mat))
    expect_true(is.matrix(correlation_s_p_mat))
    expect_true(max(correlation_s_p_mat) <= 1.0000001)
    expect_true(min(correlation_s_p_mat) >= 0)

    ## semi-partial spearman
    expect_true(all(correlation_s_sp_mat - abs(ppcor::spcor(t(mat_test[1:5, ]),
        method = "spearman")$estimate) == 0, na.rm = TRUE))
    expect_equal(rownames(correlation_s_sp_mat), colnames(correlation_s_sp_mat))
    expect_equal(rownames(correlation_s_sp_mat), rownames(mat_test)[1:5])
    expect_equal(ncol(correlation_s_sp_mat), nrow(correlation_s_sp_mat))
    expect_equal(nrow(correlation_s_sp_mat), nrow(mat_test[1:5, ]))
    expect_true(is.numeric(correlation_s_sp_mat))
    expect_true(is.matrix(correlation_s_sp_mat))
    expect_true(min(correlation_s_sp_mat, na.rm = TRUE) >= 0)
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
stat_adj_l <- statistical(mat_test[1:5, ],
    model = c("randomForest", "clr", "aracne", "pearson",
        "pearson_partial", "pearson_semipartial", "spearman",
        "spearman_partial", "spearman_semipartial", "bayes"),
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
    expect_equal(stat_adj_l[["randomForest"]], tmp, tolerance = 5e-01)
    tmp <- bayes_mat; diag(tmp) <- NaN
    expect_equal(stat_adj_l[["bayes"]], tmp, tolerance = 5e-01)
    tmp <- clr_mat; diag(tmp) <- NaN
    expect_true(all(stat_adj_l[["clr"]] == tmp, na.rm = TRUE))
    tmp <- aracne_mat; diag(tmp) <- NaN
    expect_true(all(stat_adj_l[["aracne"]] == tmp, na.rm = TRUE))
    tmp <-  correlation(mat_test[1:5, ], type = "pearson")
    diag(tmp) <- NaN
    expect_true(all(stat_adj_l[["pearson"]], tmp, na.rm = TRUE))
    tmp <- correlation(mat_test[1:5, ], type = "pearson_partial")
    diag(tmp) <- NaN
    expect_true(all(stat_adj_l[["pearson_partial"]] == tmp, na.rm = TRUE))
    tmp <- correlation(mat_test[1:5, ], type = "pearson_semipartial")
    diag(tmp) <- NaN
    expect_true(all(stat_adj_l[["pearson_semipartial"]] == tmp, na.rm = TRUE))
    tmp <- correlation(mat_test[1:5, ], type = "spearman")
    diag(tmp) <- NaN
    expect_true(all(stat_adj_l[["spearman"]] == tmp, na.rm = TRUE))
    tmp <- correlation(mat_test[1:5, ], type = "spearman_partial")
    diag(tmp) <- NaN
    expect_true(all(stat_adj_l[["spearman_partial"]] == tmp, na.rm = TRUE))
    tmp <- correlation(mat_test[1:5, ], type = "spearman_semipartial")
    diag(tmp) <- NaN
    expect_true(all(stat_adj_l[["spearman_semipartial"]] == tmp, na.rm = TRUE))
    expect_equal(length(stat_adj_l), 10)
    expect_equal(as.numeric(lapply(stat_adj_l, nrow)), rep(5, 10))
    expect_equal(as.numeric(lapply(stat_adj_l, ncol)), rep(5, 10))
    expect_equal(as.character(unlist((lapply(stat_adj_l, rownames)))),
        rep(c("x1", "x2", "x3", "x4", "x5"), 10))
    expect_equal(as.character(unlist((lapply(stat_adj_l, colnames)))),
        rep(c("x1", "x2", "x3", "x4", "x5"), 10))
    expect_equal(names(stat_adj_l),
        c("randomForest", "clr", "aracne", "pearson",
            "pearson_partial", "pearson_semipartial", "spearman",
            "spearman_partial", "spearman_semipartial", "bayes"))
    expect_true(all(unlist(lapply(stat_adj_l, function(x) is.numeric(x)))))
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
stat_adj_l_cut <- stat_adj_l[!names(stat_adj_l) %in%
    c("pearson_partial", "pearson_semipartial", "spearman_partial",
        "spearman_semipartial", "bayes", "randomForest")]
args_thr <- list(clr = 0.5, aracne = 0.8, pearson = 0.95, spearman = 0.95,
    threshold = 1)
thr_thr <- threshold(stat_adj_l_cut, type = "threshold", args = args_thr)

args_top <- list(n = 5)
thr_top1 <- threshold(stat_adj_l_cut, type = "top1", args = args_top)
thr_top2 <- threshold(stat_adj_l_cut, type = "top2", args = args_top)
thr_mean <- threshold(stat_adj_l_cut, type = "mean", args = args_top)

mat <- matrix(1:25, nrow = 5, ncol = 5)
diag(mat) <- NaN
rownames(mat) <- colnames(mat) <- paste("x", 1:5, sep = "")
mat_l <- list(pearson = mat)
thr_top1_all <- threshold(mat_l, type = "top1", args = args_top, values = "all")
thr_top1_min <- threshold(mat_l, type = "top1", args = args_top, values = "min")
thr_top1_max <- threshold(mat_l, type = "top1", args = args_top, values = "max")

test_that("threshold", {

    ## test arguments
    expect_error(threshold(NULL, type = "threshold", args = args_thr),
        "consensus requires graphs of identical order")
    expect_error(threshold(1:3, type = "threshold", args = args_thr),
        "attempt to select less than one element in")
    args_foo <- list(clr = 0.5, aracne = 0.8, threshold = 1)
    expect_error(threshold(
        data.frame(clr = 1:3, aracne = 1:3),
        type = "threshold", args = args_foo),
        "input must be an adjacency matrix/array, network, or list")
    expect_error(threshold(cbind(clr = c("a", "b", "c"),
        aracne = c("a", "b", "c")), type = "threshold", args = args_foo),
        "attempt to select less than one element in")
    args_thr_double <- list(randomForest = 0.2,
        clr = 0.5, clr = 0.5, aracne = 0.8, pearson = 0.05, spearman = 0.05,
        threshold = 1)
    expect_error(
        threshold(stat_adj_l_cut, type = "threshold", args = args_thr_double),
        "contain duplicated entries")
    expect_error(
        threshold(stat_adj_l_cut, type = "foo", args = args_thr),
        "type not in")
    expect_error(threshold(stat_adj_l_cut, type = "top1", args = args_top, 
        values = "foo"), "should be one of ")

    ## check that args contains all models
    expect_error(threshold(stat_adj_l_cut, type = "threshold", args_thr[1:3]),
        "does not contain entries for all 'model's in 'statistical'")

    ## check that args contains threshold
    expect_error(threshold(stat_adj_l_cut, type = "threshold", args_thr[1:4]),
        "'args' does not contain entry 'threshold' of length 1")

    ## check args for top1, top2, mean
    expect_error(threshold(stat_adj_l, type = "top1", args = list(x = 1)),
        "does not contain the numeric entry `n` of length 1")
    expect_error(threshold(stat_adj_l, type = "top2", args = list(x = 1)),
        "does not contain the numeric entry `n` of length 1")
    expect_error(threshold(stat_adj_l, type = "mean", args = list(x = 1)),
        "does not contain the numeric entry `n` of length 1")
    expect_error(threshold(stat_adj_l, type = "top1", args = list(n = 1:2)),
        "does not contain the numeric entry `n` of length 1")
    expect_error(threshold(stat_adj_l, type = "top2", args = list(n = 1:2)),
        "does not contain the numeric entry `n` of length 1")
    expect_error(threshold(stat_adj_l, type = "mean", args = list(n = 1:2)),
        "does not contain the numeric entry `n` of length 1")
    expect_error(threshold(stat_adj_l, type = "top1", args = list(n = "a")),
        "does not contain the numeric entry `n` of length 1")
    expect_error(threshold(stat_adj_l, type = "top2", args = list(n = "a")),
        "does not contain the numeric entry `n` of length 1")
    expect_error(threshold(stat_adj_l, type = "mean", args = list(n = "a")),
        "does not contain the numeric entry `n` of length 1")

    ## check output
    expect_true(is.matrix(thr_thr))
    expect_true(is.matrix(thr_top1))
    expect_true(is.matrix(thr_top2))
    expect_true(is.matrix(thr_mean))
    expect_true(is.numeric(thr_thr))
    expect_true(is.numeric(thr_top1))
    expect_true(is.numeric(thr_top2))
    expect_true(is.numeric(thr_mean))
    expect_equal(rownames(thr_thr), c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(colnames(thr_thr), c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(rownames(thr_top1), c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(colnames(thr_top1), c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(rownames(thr_top2), c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(colnames(thr_top2), c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(rownames(thr_mean), c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(colnames(thr_mean), c("x1", "x2", "x3", "x4", "x5"))
    expect_equal(sum(thr_thr), 12)
    expect_equal(sum(thr_top1), 16)
    expect_equal(sum(thr_top2), 14)
    expect_equal(sum(thr_mean), 10)
    expect_equal(sum(thr_top1_all), 5)
    expect_equal(sum(thr_top1_min), 10)
    expect_equal(sum(thr_top1_max), 10)
    expect_equal(c(thr_top1_all), c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0))
    expect_equal(c(thr_top1_min), c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 
        1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0))
    expect_equal(c(thr_top1_max), c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
        1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0))
    expect_true(all(thr_thr %in% c(0, 1)))
    expect_true(all(thr_top1 %in% c(0, 1)))
    expect_true(all(thr_top2 %in% c(0, 1)))
    expect_true(all(thr_mean %in% c(0, 1)))
    
    threshold(stat_adj_l_cut, type = "top1", args = args_top, values = "all")
    threshold(stat_adj_l_cut, type = "top1", args = args_top, values = "min")
    threshold(stat_adj_l_cut, type = "top1", args = args_top, values = "max")
})
## END unit test threshold  ##


## START unit test topKnet ##
ranks <- matrix(c(c(1, 2, 3), c(2, 1, 3)), ncol = 2)

test_that("topKnet", {

    ## type of ranks
    expect_error(MetNet:::topKnet(ranks = 1:3), "ranks is not a numerical")
    expect_error(MetNet:::topKnet(ranks = data.frame(x = 1:3, y = 3:1),
        type = "top1"), "ranks is not a numeric")
    expect_error(MetNet:::topKnet(ranks = matrix(c("a", "b", "c")),
        type = "top1"), "ranks is not a numeric")

    ## type argument
    expect_error(MetNet:::topKnet(ranks = ranks, type = "foo"), "type neither")

    ## check results
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "top1"), c(1, 1, 3))
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "top2"), c(2, 2, 3))
    expect_equal(MetNet:::topKnet(ranks = ranks, type = "mean"), c(1.5, 1.5, 3))

    ## matrix with ncol 1
    expect_equal(MetNet:::topKnet(ranks = matrix(1:3, nrow = 3),
        type = "top1"), 1:3)
    expect_error(MetNet:::topKnet(ranks = matrix(1:3, nrow = 3),
        type = "top2"), "ncol[(]ranks[])] has to be")
    expect_equal(MetNet:::topKnet(ranks = matrix(1:3, nrow = 3),
        type = "mean"), 1:3)
})
## END unit test topKnet ##


## START unit test threeDotsCall ##

test_that("threeDotsCall", {
    expect_error(MetNet:::threeDotsCall("mean", x = 1:10, x = 1:10),
        "duplicated args in ...")
    expect_equal(MetNet:::threeDotsCall("mean", x = 1:10, foo = 1), mean(1:10))
})
## END unit test threeDotsCall ##