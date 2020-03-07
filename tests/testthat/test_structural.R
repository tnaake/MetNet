## create toy example data set
data("mat_test", package = "MetNet")
colnames(mat_test) <- paste0("intensity_", 1:20)
mz <- c(100, 150, 262.0528, 262.0528, 262.0528, 348.0532, 448.0532)
rt <- c(100, 100, 50, 150, 150, 150, 150)
mat_test <- cbind(mz = mz, rt = rt, mat_test)

## transformations object for structual calculation
transformations <- rbind(
    c("Malonyl group (–H2O)", "C3H2O3", 86.0003939305, "+"),
    c("Monosaccharide (–H2O)", "C6H10O5", 162.0528234315, "-"))
transformations_neg <- transformations <- data.frame(
    group = transformations[, 1],
    formula = transformations[, 2],
    mass = as.numeric(transformations[, 3]),
    rt = transformations[, 4])
transformations_neg[, 3] <- -1 * transformations_neg[, 3]


## START unit test structural ##
struct_adj <- structural(mat_test,
        transformation = transformations, ppm = 5, directed = FALSE)
struct_adj_neg <- structural(mat_test,
        transformation = transformations_neg, ppm = 5, directed = FALSE)
struct_adj_dir <- structural(mat_test,
        transformation = transformations, ppm = 5, directed = TRUE)
struct_adj_dir_neg <- structural(mat_test,
        transformation = transformations_neg, ppm = 5, directed = TRUE)

test_that("structural", {
    expect_error(structural(mat_test[, -1], transformations),
        "does not contain the column mz")
    expect_error(structural(NULL, transformations),
        "does not contain the column mz")
    expect_error(structural(mat_test, transformations[, -1]),
        "does not contain the column group")
    expect_error(structural(mat_test, transformations[, -3]),
        "does not contain the column mass")
    expect_error(structural(mat_test, matrix()),
        "is not a data.frame")
    expect_error(structural(mat_test, transformations, ppm = "a"),
        "is not numeric")
    expect_equal(length(struct_adj), 2)
    expect_equal(dim(struct_adj[[1]]), c(7, 7))
    expect_equal(dim(struct_adj[[2]]), c(7, 7))
    expect_equal(rownames(struct_adj[[1]]), colnames(struct_adj[[1]]))
    expect_equal(rownames(struct_adj[[2]]), colnames(struct_adj[[2]]))
    expect_equal(rownames(struct_adj[[1]]), rownames(struct_adj[[2]]))
    expect_equal(rownames(struct_adj[[1]]), paste0("x", 1:7))
    expect_equal(sum(struct_adj[[1]]), 12)
    expect_equal(sum(struct_adj_neg[[1]]), 0)
    expect_equal(sum(struct_adj_dir[[1]]), 6)
    expect_equal(sum(struct_adj_dir_neg[[1]]), 6)
    expect_equal(unique(as.vector(struct_adj[[2]])),
                c("", "Monosaccharide (–H2O)", "Malonyl group (–H2O)"))
    expect_true(is.matrix(struct_adj[[1]]))
    expect_true(is.matrix(struct_adj[[2]]))
    expect_true(is.numeric(struct_adj[[1]]))
    expect_true(is.character(struct_adj[[2]]))
})
## END unit test structural ##

## START unit test rtCorrection ##
struct_adj_rt <- rtCorrection(struct_adj, mat_test, transformations)
struct_adj_rt_dir <- rtCorrection(struct_adj_dir, mat_test, transformations)

test_that("rtCorrection", {
    expect_error(rtCorrection(struct_adj[[1]], mat_test, transformations),
        "is not a list")
    expect_error(
        rtCorrection(list(struct_adj[[1]]), mat_test, transformations),
        "is not a list of length 2")
    expect_error(rtCorrection(list(struct_adj[[2]], struct_adj[[2]]),
        mat_test, transformations))
    expect_error(rtCorrection(list(struct_adj[[1]], struct_adj[[1]]),
        mat_test, transformations), "is not a character matrix")
    expect_error(rtCorrection(struct_adj, NULL, transformations),
        "does not contain the column rt")
    expect_error(rtCorrection(struct_adj, mat_test[, -1], transformations),
        "subscript out of bounds")
    expect_error(rtCorrection(struct_adj, mat_test[, -2], transformations),
        "does not contain the column rt")
    expect_error(rtCorrection(struct_adj, mat_test, NULL),
        "does not contain the column group")
    expect_error(rtCorrection(struct_adj, mat_test, transformations[, -1]),
        "does not contain the column group")
    expect_error(rtCorrection(struct_adj, mat_test, transformations[, -3]),
        "does not contain the column mz")
    expect_error(rtCorrection(struct_adj, mat_test, transformations[, -4]),
        "does not contain the column rt")
    expect_error(rtCorrection(struct_adj, mat_test,
        transformation = cbind(transformations[, -4], rt = rep("a", 4))))#,
    #    "does contain other")
    expect_true(is.matrix(struct_adj_rt[[1]]))
    expect_true(is.numeric(struct_adj_rt[[1]]))
    expect_true(is.matrix(struct_adj_rt[[2]]))
    expect_true(is.character(struct_adj_rt[[2]]))
    expect_equal(colnames(struct_adj_rt[[1]]), paste0("x", 1:7))
    expect_equal(colnames(struct_adj_rt[[1]]), rownames(struct_adj_rt[[1]]))
    expect_equal(colnames(struct_adj_rt[[1]]), colnames(struct_adj_rt[[2]]))
    expect_equal(colnames(struct_adj_rt[[1]]), rownames(struct_adj_rt[[2]]))
    expect_true(table(struct_adj_rt[[1]])[1] == 41)
    expect_true(table(struct_adj_rt_dir[[1]])[1] == 45)
    expect_equal(sum(struct_adj_rt[[1]]), 8)
    expect_equal(sum(struct_adj_rt_dir[[1]]), 4)

    ## dims of struct_adj[[1]] and struct[[2]], rownames/colnames
    foo_1 <- struct_adj[[1]]
    foo_2 <- struct_adj[[2]]
    foo <- list(foo_1, foo_2[, 1:6])
    expect_error(rtCorrection(foo, mat_test, transformations),
        " is not equal to dim")

    rownames(foo_1)[1] <- "foo"
    foo <- list(foo_1, foo_2)
    expect_error(rtCorrection(foo, mat_test, transformations),
        "are not identical to rownames of ")

    foo_1 <- struct_adj[[1]]
    rownames(foo_2)[1] <- "foo"
    foo <- list(foo_1, foo_2)
    expect_error(rtCorrection(foo, mat_test, transformations),
        "are not identical to ")

    foo_2 <- struct_adj[[2]]
    foo <- list(foo_1, foo_2)
    expect_error(rtCorrection(foo, mat_test[1:6, ], transformations),
        "do not fit rownames[(]x[])]")

    foo_1 <- as.data.frame(struct_adj[[1]])
    foo <- list(foo_1, foo_2)
    expect_error(rtCorrection(foo, mat_test, transformations),
        "is not a numeric matrix")

    foo_1 <- struct_adj[[1]]
    foo_2 <- as.data.frame(struct_adj[[2]])
    foo <- list(foo_1, foo_2)
    expect_error(rtCorrection(foo, mat_test, transformations),
        "is not a character matrix")
})
## END unit test rtCorrection ##
