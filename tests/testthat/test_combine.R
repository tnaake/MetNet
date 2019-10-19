## create toy example data set
data("mat_test", package = "MetNet")
colnames(mat_test) <- paste0("intensity_", 1:20)
mz <- c(100, 150, 200, 200, 262.0528, 348.0532, 448.0532)
mat_test <- cbind(mz = mz, mat_test)

## transformations object for structual calculation
transformations <- rbind(
    c("Malonyl group (–H2O)", "C3H2O3", "86.0003939305"),
    c("Monosaccharide (–H2O)", "C6H10O5", "162.0528234315"))
transformations <- data.frame(group = as.character(transformations[, 1]),
                                formula = as.character(transformations[, 2]),
                                mass = as.numeric(transformations[, 3]))

## create structural network
struct_adj <- structural(mat_test,
        transformation = transformations, ppm = 5)

## create statistical network
stat_adj_l <- statistical(mat_test,
    model = c("clr", "aracne", "pearson", "spearman", "bayes"))
stat_adj <- threshold(stat_adj_l, type = "top2", args = list(n = 10))

## START unit test combine ##
cons_adj <- combine(struct_adj, stat_adj)

test_that("combine", {
    expect_error(combine(NULL, stat_adj), "not a list of length 2")
    expect_error(combine(struct_adj, NULL), "not a numeric matrix")
    expect_error(combine(struct_adj, stat_adj, threshold = "a"),
        "is not numeric")
    expect_equal(sum(cons_adj[[1]]), 4)
    expect_equal(dim(cons_adj[[1]]), c(7, 7))
    expect_equal(rownames(cons_adj[[1]]), paste0("x", 1:7))
    expect_equal(rownames(cons_adj[[2]]), paste0("x", 1:7))
    expect_equal(colnames(cons_adj[[1]]), paste0("x", 1:7))
    expect_equal(colnames(cons_adj[[2]]), paste0("x", 1:7))
    expect_equal(rownames(cons_adj[[1]]), colnames(cons_adj[[1]]))
    expect_equal(rownames(cons_adj[[2]]), colnames(cons_adj[[2]]))
    expect_true(is.matrix(cons_adj[[1]]))
    expect_true(is.matrix(cons_adj[[2]]))
    expect_true(is.numeric(cons_adj[[1]]))
    expect_true(is.character(cons_adj[[2]]))

    ## check rownames/colnames
    mock <- struct_adj
    mock_num <- mock[[1]]
    colnames(mock_num)[1] <- "foo"
    mock[[1]] <- mock_num
    expect_error(combine(mock, stat_adj), "are not identical to ")

    mock <- struct_adj
    mock_num <- mock[[1]]
    rownames(mock_num)[1] <- "foo"
    mock[[1]] <- mock_num
    expect_error(combine(mock, stat_adj), "are not identical to")

    mock <- stat_adj
    colnames(mock)[1] <- "foo"
    expect_error(combine(struct_adj, mock), "are not identical")

    mock <- stat_adj
    rownames(mock)[1] <- "foo"
    expect_error(combine(struct_adj, mock), "are not identical")

    ## check for structural
    l <- list(1, cons_adj[[2]])
    expect_error(combine(structural = l, statistical = stat_adj),
        "not a numeric matrix")
    l <- list(cons_adj[[1]], 1)
    expect_error(combine(structural = l, statistical = stat_adj),
        "not a character ")
})
## END unit test combine ##
