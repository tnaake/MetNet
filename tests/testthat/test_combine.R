## create toy example data set
data("mat_test", package = "MetNet")
colnames(mat_test) <- paste0("intensity_", 1:20)
mz <- c(100, 150, 200, 200, 262.0528, 348.0532, 448.0532)
mat_test <- cbind(mz = mz, mat_test)

## transformations object for structual calculation
transformations <- rbind(
    c("Malonyl group (-H2O)", "C3H2O3", "86.0003939305"),
    c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"))
transformations <- data.frame(group = as.character(transformations[, 1]),
                                formula = as.character(transformations[, 2]),
                                mass = as.numeric(transformations[, 3]))

## create structural network
struct_adj <- structural(mat_test,
        transformation = transformations, ppm = 5)

## create statistical network
stat_adj <- statistical(mat_test[, -1],
    model = c("clr", "aracne", "pearson", "spearman", "bayes"))
stat_adj_thr <- threshold(stat_adj, type = "top2", args = list(n = 10))

## START unit test combine ##
cons_adj <- combine(struct_adj, stat_adj_thr)

test_that("combine", {
    expect_error(combine(NULL, stat_adj_thr), 
        "'am_structural' is not an 'AdjacencyMatrix' object")
    expect_error(combine(struct_adj, NULL),
        "'am_statistical' is not an 'AdjacencyMatrix' object")
    expect_error(combine(struct_adj, stat_adj),
        "'am_statistical' must contain assay 'consensus'")
    expect_true(validObject(cons_adj))
    expect_equal(assayNames(cons_adj), 
        c("clr_coef", "aracne_coef", "pearson_coef", "pearson_pvalue", 
            "spearman_coef", "spearman_pvalue", "bayes_coef", "consensus",
            "binary", "transformation", "mass_difference", "combine_binary",        
            "combine_transformation", "combine_mass_difference"))
    expect_equal(sum(assay(cons_adj, "combine_binary"), na.rm = TRUE), 4)
    expect_equal(dim(assay(cons_adj, "combine_binary")), c(7, 7))
    expect_equal(rownames(assay(cons_adj, "combine_binary")), paste0("x", 1:7))
    expect_equal(rownames(assay(cons_adj, "combine_transformation")), 
        paste0("x", 1:7))
    expect_equal(rownames(assay(cons_adj, "combine_mass_difference")), 
        paste0("x", 1:7))
    expect_equal(colnames(assay(cons_adj, "combine_binary")), paste0("x", 1:7))
    expect_equal(colnames(assay(cons_adj, "combine_transformation")), 
                 paste0("x", 1:7))
    expect_equal(colnames(assay(cons_adj, "combine_mass_difference")), 
                 paste0("x", 1:7))
    expect_equal(rownames(assay(cons_adj, "combine_binary")),
        colnames(assay(cons_adj, "combine_binary")))
    expect_equal(rownames(assay(cons_adj, "combine_transformation")),
        colnames(assay(cons_adj, "combine_transformation")))
    expect_equal(rownames(assay(cons_adj, "combine_mass_difference")),
        colnames(assay(cons_adj, "combine_mass_difference")))
    expect_true(is.numeric(assay(cons_adj, "combine_binary")))
    expect_true(is.character(assay(cons_adj, "combine_transformation")))
    expect_true(is.character(assay(cons_adj, "combine_mass_difference")))

    ## check rownames/colnames
    struct_adj_rev <- structural(mat_test[7:1, ],
        transformation = transformations, ppm = 5)
    
    expect_error(combine(struct_adj_rev, stat_adj_thr),
        "are not identical to those of the receiving\n  AdjacencyMatrix object")
})
## END unit test combine ##
