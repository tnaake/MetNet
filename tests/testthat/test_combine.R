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
struct_adj <- structural(mat_test, transformation = transformations, 
    var = c("group", "formula", "mass"), ppm = 5)

## create statistical network
stat_adj <- statistical(mat_test[, -1],
    model = c("clr", "aracne", "pearson", "spearman", "ggm"))
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
        c("binary", "group", "formula", "mass", 
            "clr_coef", "aracne_coef", "pearson_coef", "pearson_pvalue", 
          "spearman_coef", "spearman_pvalue","ggm_coef", "ggm_pvalue",
          "consensus","combine_binary", "combine_group", 
          "combine_formula",  "combine_mass"))
    expect_equal(sum(assay(cons_adj, "combine_binary"), na.rm = TRUE), 4)
    expect_equal(as.vector(assay(cons_adj, "combine_binary")[1, ]), 
        c(NA, 0, 0, 0, 1, 0, 0))
    expect_equal(as.vector(assay(cons_adj, "combine_group")[1, ]), 
        c(NA, "", "", "", "Monosaccharide (-H2O)", "", ""))
    expect_equal(as.vector(assay(cons_adj, "combine_formula")[1, ]), 
        c(NA, "", "", "", "C6H10O5", "", ""))
    expect_equal(as.vector(assay(cons_adj, "combine_mass")[1, ]), 
        c(NA, "", "", "", "162.0528234315", "", ""))
    expect_equal(dim(assay(cons_adj, "combine_binary")), c(7, 7))
    expect_equal(rownames(assay(cons_adj, "combine_binary")), paste0("x", 1:7))
    expect_equal(rownames(assay(cons_adj, "combine_group")), paste0("x", 1:7))
    expect_equal(rownames(assay(cons_adj, "combine_formula")), paste0("x", 1:7))
    expect_equal(rownames(assay(cons_adj, "combine_mass")), paste0("x", 1:7))
    expect_equal(colnames(assay(cons_adj, "combine_binary")), paste0("x", 1:7))
    expect_equal(colnames(assay(cons_adj, "combine_group")), paste0("x", 1:7))
    expect_equal(colnames(assay(cons_adj, "combine_formula")), paste0("x", 1:7))
    expect_equal(colnames(assay(cons_adj, "combine_mass")), paste0("x", 1:7))
    expect_equal(rownames(assay(cons_adj, "combine_binary")),
        colnames(assay(cons_adj, "combine_binary")))
    expect_equal(rownames(assay(cons_adj, "combine_group")),
        colnames(assay(cons_adj, "combine_group")))
    expect_equal(rownames(assay(cons_adj, "combine_formula")),
        colnames(assay(cons_adj, "combine_formula")))
    expect_equal(rownames(assay(cons_adj, "combine_mass")),
        colnames(assay(cons_adj, "combine_mass")))
    expect_true(is.numeric(assay(cons_adj, "combine_binary")))
    expect_true(is.character(assay(cons_adj, "combine_group")))
    expect_true(is.character(assay(cons_adj, "combine_formula")))
    expect_true(is.character(assay(cons_adj, "combine_mass")))

    ## check rownames/colnames
    struct_adj_rev <- structural(mat_test[7:1, ],
        transformation = transformations, ppm = 5)
    
    expect_error(combine(struct_adj_rev, stat_adj_thr),
        "names of 'am_structural' do not match names of 'am_statistical'")
    
    ## test with NA values
    mat_test_NA <- mat_test[, -1]
    mat_test_NA[1, 5] <- NA
    stat_adj_NA <- statistical(mat_test_NA, model = c("pearson", "ggm"))
    stat_adj_thr_NA <- threshold(stat_adj_NA, type = "top1", 
        args = list(n = 10, abs = TRUE), na.rm = FALSE)
    cons_adj <- combine(struct_adj, stat_adj_thr_NA)
    expect_true(all(is.na(as.vector(assay(cons_adj, "combine_binary")[, 1]))))
    expect_equal(as.vector(assay(cons_adj, "combine_binary")[, 2]), 
        c(NA, NA, 0, 0, 0, 0, 0))
        
})
## END unit test combine ##
