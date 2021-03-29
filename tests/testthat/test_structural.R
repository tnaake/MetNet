## create toy example data set
data("mat_test", package = "MetNet")
colnames(mat_test) <- paste0("intensity_", 1:20)
mz <- c(100, 150, 262.0528, 262.0528, 262.0528, 348.0532, 448.0532)
rt <- c(100, 100, 50, 150, 150, 150, 150)
mat_test <- cbind(mz = mz, rt = rt, mat_test)

## transformations object for structual calculation
transformations <- rbind(
    c("Malonyl group (-H2O)", "C3H2O3", 86.0003939305, "+"),
    c("Monosaccharide (-H2O)", "C6H10O5", 162.0528234315, "-"))
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
    expect_true(validObject(struct_adj))
    expect_equal(assayNames(struct_adj), 
        c("binary", "transformation", "mass_difference"))
    expect_equal(length(struct_adj), 7)
    expect_equal(dim(struct_adj), c(7, 7))
    expect_equal(assay(struct_adj, "binary"), c(7, 7))
    expect_equal(assay(struct_adj, "transformation"), c(7, 7))
    expect_equal(assay(struct_adj, "mass_difference"), c(7, 7))
    expect_equal(rownames(assay(struct_adj, 1)), colnames(assay(struct_adj, 1)))
    expect_equal(rownames(assay(struct_adj, 2)), colnames(assay(struct_adj, 2)))
    expect_equal(rownames(assay(struct_adj, 3)), rownames(assay(struct_adj, 3)))
    expect_equal(rownames(assay(struct_adj, 1)), rownames(assay(struct_adj, 2)))
    expect_equal(rownames(assay(struct_adj, 1)), rownames(assay(struct_adj, 3)))
    expect_equal(rownames(assay(struct_adj, 1)), paste0("x", 1:7))
    expect_equal(sum(assay(struct_adj, "binary")), 12)
    expect_equal(sum(assay(struct_adj_neg, "binary")), 0)
    expect_equal(sum(assay(struct_adj_dir, "binary")), 6)
    expect_equal(sum(assay(struct_adj_dir_neg, "binary")), 6)
    expect_equal(unique(as.vector(assay(struct_adj, "transformation"))),
        c("", "Monosaccharide (-H2O)", "Malonyl group (-H2O)"))
    expect_equal(unique(as.vector(assay(struct_adj, "mass_difference"))),
                 c("", "162.0528234315", "86.0003939305"))
    expect_true(is.matrix(assay(struct_adj, "binary")))
    expect_true(is.matrix(assay(struct_adj, "transformation")))
    expect_true(is.matrix(assay(struct_adj, "mass_difference")))
    expect_true(is.numeric(assay(struct_adj, "binary")))
    expect_true(is.character(assay(struct_adj, "transformation")))
    expect_true(is.character(assay(struct_adj, "mass_difference")))
    expect_equal(directed(struct_adj), FALSE)
    expect_equal(directed(struct_adj_neg), FALSE)
    expect_equal(directed(struct_adj_dir), TRUE)
    expect_equal(directed(struct_adj_dir_neg), TRUE)
    expect_equal(type(struct_adj), "structural")
    expect_equal(type(struct_adj_neg), "structural")
    expect_equal(type(struct_adj_dir), "structural")
    expect_equal(type(struct_adj_dir_neg), "structural")
    expect_equal(thresholded(struct_adj), FALSE)
    expect_equal(thresholded(struct_adj_neg), FALSE)
    expect_equal(thresholded(struct_adj_dir), FALSE)
    expect_equal(thresholded(struct_adj_dir_neg), FALSE)
})
## END unit test structural ##

## START unit test rtCorrection ##
struct_adj_rt <- rtCorrection(struct_adj, mat_test, transformations)
struct_adj_rt_dir <- rtCorrection(struct_adj_dir, mat_test, transformations)

test_that("rtCorrection", {
    expect_error(rtCorrection(struct_adj[[1]], mat_test, transformations),
        "is not an 'AdjacencyMatrix'")
    tmp <- struct_adj
    assay(tmp, "transformation") <- assay(tmp, "binary")
    expect_error(rtCorrection(tmp, mat_test, transformations), 
        "must be character")
    tmp <- struct_adj
    assay(tmp, "mass_difference") <- assay(tmp, "binary")
    expect_error(rtCorrection(tmp, mat_test, transformations), 
        "must be character")
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
        transformation = cbind(transformations[, -4], rt = rep("a", 2))),
        "does contain other")
    expect_true(is.matrix(assay(struct_adj_rt, "binary")))
    expect_true(is.numeric(assay(struct_adj_rt, "binary")))
    expect_true(is.matrix(assay(struct_adj_rt, "transformation")))
    expect_true(is.character(assay(struct_adj_rt, "transformation")))
    expect_true(is.matrix(assay(struct_adj_rt, "mass_difference")))
    expect_true(is.character(assay(struct_adj_rt, "mass_difference")))
    expect_equal(colnames(assay(struct_adj_rt, "binary")), paste0("x", 1:7))
    expect_equal(colnames(assay(struct_adj_rt, 1)),
        rownames(assay(struct_adj_rt, 1)))
    expect_equal(colnames(assay(struct_adj_rt, 2)),
        rownames(assay(struct_adj_rt, 2)))
    expect_equal(colnames(assay(struct_adj_rt, 3)),
        rownames(assay(struct_adj_rt, 3)))
    expect_equal(colnames(assay(struct_adj_rt, 1)),
        colnames(assay(struct_adj_rt, 2)))
    expect_equal(colnames(assay(struct_adj_rt, 1)),
        colnames(assay(struct_adj_rt, 3)))
    expect_true(table(assay(struct_adj_rt, "binary"))[1] == 41)
    expect_true(table(assay(struct_adj_rt_dir, "binary"))[1] == 45)
    expect_equal(sum(assay(struct_adj_rt, "binary")), 8)
    expect_equal(sum(assay(struct_adj_rt_dir, "binary")), 4)
    expect_true(all(
        table(assay(struct_adj_rt, "transformation")) == c(41, 4, 4)))
    expect_true(all(
        table(assay(struct_adj_rt, "mass_difference")) == c(41, 4, 4)))
    expect_true(all(
        table(assay(struct_adj_rt_dir, "transformation")) == c(45, 2, 2)))
    expect_true(all(
        table(assay(struct_adj_rt_dir, "mass_difference")) == c(45, 2, 2)))

    expect_equal(directed(struct_adj_rt), FALSE)
    expect_equal(directed(struct_adj_rt_dir), TRUE)
    expect_equal(type(struct_adj_rt), "structural")
    expect_equal(type(struct_adj_rt_dir), "structural")
    expect_equal(thresholded(struct_adj_rt), TRUE)
    expect_equal(thresholded(struct_adj_rt_dir), TRUE)
})
## END unit test rtCorrection ##
