## create toy example data set
data("mat_test", package = "MetNet")
colnames(mat_test) <- paste0("intensity_", 1:20)
mz <- c(100, 150, 262.0528, 262.0528, 262.0528, 348.0532, 448.0532)
rt <- c(100, 100, 50, 150, 150, 150, 150)
mat_test <- cbind(mz = mz, rt = rt, mat_test)
rownames(mat_test) <- paste(mz, rt, sep = "_")

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
struct_adj <- structural(mat_test, transformation = transformations, 
    var = c("group", "formula", "mass"), ppm = 5, directed = FALSE)
struct_adj_neg <- structural(mat_test, transformation = transformations_neg, 
    var = c("group", "formula", "mass"), ppm = 5, directed = FALSE)
struct_adj_dir <- structural(mat_test, transformation = transformations, 
    var = c("group", "formula", "mass"), ppm = 5, directed = TRUE)
struct_adj_dir_neg <- structural(mat_test, transformation = transformations_neg, 
    var = c("group", "formula", "mass"), ppm = 5, directed = TRUE)

g_undir <- igraph::graph_from_adjacency_matrix(
    assay(struct_adj, "binary", mode = "directed", weighted = NULL))
g_undir_neg <- igraph::graph_from_adjacency_matrix(
    assay(struct_adj_neg, "binary", mode = "directed", weighted = NULL))
g_dir <- igraph::graph_from_adjacency_matrix(
    assay(struct_adj_dir, "binary"), mode = "directed", weighted = NULL)
g_dir_neg <- igraph::graph_from_adjacency_matrix(
    assay(struct_adj_dir_neg, "binary"), mode = "directed", weighted = NULL)

plot(g_undir, edge.width = 1, edge.arrow.size = 0.5, 
    vertex.label.cex = 0.8, edge.color = "grey")
plot(g_undir_neg, edge.width = 1, edge.arrow.size = 0.5, 
     vertex.label.cex = 0.8, edge.color = "grey")
plot(g_dir, edge.width = 1, edge.arrow.size = 0.5, 
     vertex.label.cex = 0.8, edge.color = "grey")
plot(g_dir_neg, edge.width = 1, edge.arrow.size = 0.5, 
     vertex.label.cex = 0.8, edge.color = "grey")

test_that("structural", {
    expect_error(structural(mat_test[, -1], transformations),
        "does not contain the column mz")
    expect_error(structural(NULL, transformations),
        "'x' has to be a matrix or data.frame")
    expect_error(structural(mat_test, transformations[, -3]),
        "does not contain the column mass")
    expect_error(structural(mat_test, matrix()),
        "is not a data.frame")
    expect_error(structural(mat_test, transformations, ppm = "a"),
        "'ppm' has to be a numeric of length 1")
    expect_error(structural(mat_test, transformations, var = c("group", "foo")),
        "'transformation' does not contain the column 'foo'")
    expect_error(structural(mat_test, transformations, var = "foo"),
        "'transformation' does not contain the column 'foo'")
    expect_error(structural(mat_test, transformations, 
        var = c("group", "foo", "foo2")),
        "'transformation' does not contain the column 'foo', 'foo2'")
    expect_error(structural(mat_test, transformations, 
        var = c("foo", "foo2")),
        "'transformation' does not contain the column 'foo', 'foo2'")
    expect_error(structural(mat_test, transformations, var = ""),
        "'transformation' does not contain the column ''")
    expect_error(structural(mat_test, transformations, var = NULL), 
        "'var' is not a character vector")
    expect_error(structural(mat_test, transformations, var = numeric()), 
        "'var' is not a character vector")
    expect_error(structural(mat_test, transformations, var = logical()), 
        "'var' is not a character vector")
    expect_equal(assayNames(structural(mat_test, transformations, 
        var = character())), "binary")
    expect_true(validObject(struct_adj))
    expect_equal(assayNames(struct_adj), 
        c("binary", "group", "formula", "mass"))
    expect_equal(length(struct_adj), 7)
    expect_equal(dim(struct_adj), c(7, 7))
    expect_equal(dim(assay(struct_adj, "binary")), c(7, 7))
    expect_equal(dim(assay(struct_adj, "group")), c(7, 7))
    expect_equal(dim(assay(struct_adj, "formula")), c(7, 7))
    expect_equal(dim(assay(struct_adj, "mass")), c(7, 7))
    expect_equal(rownames(assay(struct_adj, 1)), colnames(assay(struct_adj, 1)))
    expect_equal(rownames(assay(struct_adj, 2)), colnames(assay(struct_adj, 2)))
    expect_equal(rownames(assay(struct_adj, 3)), rownames(assay(struct_adj, 3)))
    expect_equal(rownames(assay(struct_adj, 4)), rownames(assay(struct_adj, 4)))
    expect_equal(rownames(assay(struct_adj, 1)), rownames(assay(struct_adj, 2)))
    expect_equal(rownames(assay(struct_adj, 1)), rownames(assay(struct_adj, 3)))
    expect_equal(rownames(assay(struct_adj, 1)), rownames(assay(struct_adj, 4)))
    expect_equal(rownames(assay(struct_adj, 1)), paste(mz, rt, sep = "_"))
    expect_equal(sum(assay(struct_adj, "binary")), 12)
    expect_equal(sum(assay(struct_adj_neg, "binary")), 0)
    expect_equal(sum(assay(struct_adj_dir, "binary")), 6)
    expect_equal(sum(assay(struct_adj_dir_neg, "binary")), 6)
    expect_equal(as.vector(assay(struct_adj, "group")[, 1]),
        c("", "", "Monosaccharide (-H2O)", "Monosaccharide (-H2O)", 
            "Monosaccharide (-H2O)", "", ""))
    expect_equal(as.vector(assay(struct_adj, "group")[, 2]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj, "group")[, 3]),
        c("Monosaccharide (-H2O)", "", "", "", "", "Malonyl group (-H2O)", ""))
    expect_equal(as.vector(assay(struct_adj, "group")[, 4]),
        c("Monosaccharide (-H2O)", "", "", "", "", "Malonyl group (-H2O)", ""))
    expect_equal(as.vector(assay(struct_adj, "group")[, 5]),
        c("Monosaccharide (-H2O)", "", "", "", "", "Malonyl group (-H2O)", ""))
    expect_equal(as.vector(assay(struct_adj, "group")[, 6]),
        c("", "", "Malonyl group (-H2O)", "Malonyl group (-H2O)", 
            "Malonyl group (-H2O)", "", ""))
    expect_equal(as.vector(assay(struct_adj, "group")[, 7]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_neg, "group")[, 1]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_neg, "group")[, 2]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_neg, "group")[, 3]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_neg, "group")[, 4]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_neg, "group")[, 5]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_neg, "group")[, 6]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_neg, "group")[, 7]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir, "group")[, 1]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir, "group")[, 2]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir, "group")[, 3]),
        c("Monosaccharide (-H2O)", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir, "group")[, 4]),
        c("Monosaccharide (-H2O)", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir, "group")[, 5]),
        c("Monosaccharide (-H2O)", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir, "group")[, 6]),
        c("", "", "Malonyl group (-H2O)", "Malonyl group (-H2O)", 
            "Malonyl group (-H2O)", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir, "group")[, 7]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir_neg, "group")[, 1]),
        c("", "", "Monosaccharide (-H2O)", "Monosaccharide (-H2O)", 
          "Monosaccharide (-H2O)", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir_neg, "group")[, 2]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir_neg, "group")[, 3]),
        c("", "", "", "", "", "Malonyl group (-H2O)", ""))
    expect_equal(as.vector(assay(struct_adj_dir_neg, "group")[, 4]),
        c("", "", "", "", "", "Malonyl group (-H2O)", ""))
    expect_equal(as.vector(assay(struct_adj_dir_neg, "group")[, 5]),
        c("", "", "", "", "", "Malonyl group (-H2O)", ""))
    expect_equal(as.vector(assay(struct_adj_dir_neg, "group")[, 6]),
        c("", "", "", "", "", "", ""))
    expect_equal(as.vector(assay(struct_adj_dir_neg, "group")[, 7]),
        c("", "", "", "", "", "", ""))
    expect_equal(unique(as.vector(assay(struct_adj, "group"))),
        c("", "Monosaccharide (-H2O)", "Malonyl group (-H2O)"))
    expect_equal(unique(as.vector(assay(struct_adj, "mass"))),
        c("", "162.0528234315", "86.0003939305"))
    expect_equal(unique(as.vector(assay(struct_adj, "formula"))),
        c("", "C6H10O5", "C3H2O3"))
    expect_true(is.matrix(assay(struct_adj, "binary")))
    expect_true(is.matrix(assay(struct_adj, "group")))
    expect_true(is.matrix(assay(struct_adj, "mass")))
    expect_true(is.matrix(assay(struct_adj, "formula")))
    expect_true(is.numeric(assay(struct_adj, "binary")))
    expect_true(is.character(assay(struct_adj, "group")))
    expect_true(is.character(assay(struct_adj, "mass")))
    expect_true(is.character(assay(struct_adj, "formula")))
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
struct_adj_rt <- rtCorrection(struct_adj, mat_test, transformations, 
    var = "group")
struct_adj_rt_dir <- rtCorrection(struct_adj_dir, mat_test, transformations,
    var = "group")

g_undir <- igraph::graph_from_adjacency_matrix(
    assay(struct_adj, "binary", mode = "directed", weighted = NULL))
g_undir_rt <- igraph::graph_from_adjacency_matrix(
    assay(struct_adj_rt, "binary", mode = "directed", weighted = NULL))
g_dir <- igraph::graph_from_adjacency_matrix(
    assay(struct_adj_dir, "binary"), mode = "directed", weighted = NULL)
g_dir_rt <- igraph::graph_from_adjacency_matrix(
    assay(struct_adj_rt_dir, "binary"), mode = "directed", weighted = NULL)

plot(g_undir, edge.width = 1, edge.arrow.size = 0.5, 
     vertex.label.cex = 0.8, edge.color = "grey")
plot(g_undir_rt, edge.width = 1, edge.arrow.size = 0.5, 
     vertex.label.cex = 0.8, edge.color = "grey")
plot(g_dir, edge.width = 1, edge.arrow.size = 0.5, 
     vertex.label.cex = 0.8, edge.color = "grey")
plot(g_dir_rt, edge.width = 1, edge.arrow.size = 0.5, 
     vertex.label.cex = 0.8, edge.color = "grey")

test_that("rtCorrection", {
    expect_error(rtCorrection(struct_adj[[1]], mat_test, transformations),
        "is not an 'AdjacencyMatrix'")
    tmp <- struct_adj
    assay(tmp, "transformation") <- assay(tmp, "binary")
    expect_error(rtCorrection(tmp, mat_test, transformations), 
        "must be character")
    tmp <- struct_adj
    assay(tmp, "mass") <- assay(tmp, "binary")
    expect_error(rtCorrection(tmp, mat_test, transformations), 
        "must be character")
    expect_error(rtCorrection(struct_adj, NULL, transformations),
        "'x' does not contain the column 'rt'")
    expect_error(rtCorrection(struct_adj, mat_test[, -1], transformations),
        "'x' does not contain the column 'mz'")
    expect_error(rtCorrection(struct_adj, mat_test[, -2], transformations),
        "'x' does not contain the column 'rt'")
    expect_error(rtCorrection(struct_adj, mat_test, NULL),
        "'transformation' does not contain the column 'group'")
    expect_error(rtCorrection(struct_adj, mat_test, transformations[, -1]),
        "'transformation' does not contain the column 'group'")
    expect_error(rtCorrection(struct_adj, mat_test, transformations[, -4]),
        "'transformation' does not contain the column 'rt'")
    expect_error(rtCorrection(struct_adj, mat_test,
        transformation = cbind(transformations[, -4], rt = rep("a", 2))),
        "does contain other")
    expect_true(is.matrix(assay(struct_adj_rt, "binary")))
    expect_true(is.numeric(assay(struct_adj_rt, "binary")))
    expect_true(is.matrix(assay(struct_adj_rt, "group")))
    expect_true(is.character(assay(struct_adj_rt, "group")))
    expect_true(is.matrix(assay(struct_adj_rt, "mass")))
    expect_true(is.character(assay(struct_adj_rt, "mass")))
    expect_true(is.matrix(assay(struct_adj_rt, "formula")))
    expect_true(is.character(assay(struct_adj_rt, "formula")))
    expect_equal(colnames(assay(struct_adj_rt, "binary")), 
        paste(mz, rt, sep = "_"))
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
        table(assay(struct_adj_rt, "group")) == c(41, 6, 2)))
    expect_true(all(
        table(assay(struct_adj_rt, "mass")) == c(41, 2, 6)))
    expect_true(all(
        table(assay(struct_adj_rt, "formula")) == c(41, 6, 2)))
    expect_true(all(
        table(assay(struct_adj_rt_dir, "group")) == c(45, 3, 1)))
    expect_true(all(
        table(assay(struct_adj_rt_dir, "mass")) == c(45, 1, 3)))
    expect_true(all(
        table(assay(struct_adj_rt_dir, "formula")) == c(45, 3, 1)))

    expect_equal(directed(struct_adj_rt), FALSE)
    expect_equal(directed(struct_adj_rt_dir), TRUE)
    expect_equal(type(struct_adj_rt), "structural")
    expect_equal(type(struct_adj_rt_dir), "structural")
    expect_equal(thresholded(struct_adj_rt), TRUE)
    expect_equal(thresholded(struct_adj_rt_dir), TRUE)
    
    ## create a double assignment to certain cells
    transformations <- rbind(transformations, 
            data.frame(group = "pseudo Monosaccharide",
               formula = "C6H10O5", mass = 162.05282, rt = "?"))
    struct_adj_pseudo <- structural(mat_test, transformation = transformations,
        var = c("group", "formula", "mass"), ppm = 5, directed = FALSE)
    expect_equal(assay(struct_adj_pseudo, "group")[3, 1],
        "Monosaccharide (-H2O)/pseudo Monosaccharide")
    expect_equal(assay(struct_adj_pseudo, "mass")[3, 1],
        "162.0528234315/162.05282")
    expect_equal(assay(struct_adj_pseudo, "formula")[3, 1],
        "C6H10O5/C6H10O5")
    struct_adj_pseudo_rt <- rtCorrection(struct_adj_pseudo, x = mat_test, 
        transformation = transformations, var = "group")
})
## END unit test rtCorrection ##
