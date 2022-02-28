## create toy example data set
data("x_test", package = "MetNet")
x_test <- as.matrix(x_test)


## transformations object for structual calculation
transformations <- rbind(
  c("Malonyl group (-H2O)", "C3H2O3", "86.0003939305"),
  c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"))
transformations <- data.frame(group = as.character(transformations[, 1]),
                              formula = as.character(transformations[, 2]),
                              mass = as.numeric(transformations[, 3]))

## create structural network
struct_adj <- structural(x_test, transformation = transformations, 
                         var = c("group", "formula", "mass"), ppm = 5)


## create statistical network
stat_adj <- statistical(x_test[, -c(1:2)],
                        model = c("clr", "aracne", "pearson", "spearman", "ggm"))
stat_adj_thr <- threshold(stat_adj, type = "top2", args = list(n = 10))


## load MS2 
data("ms2_test", package = "MetNet")

## START unit test spectral similarity ##
spect_adj <- addSpectSimil(sps_sub, am_structural = struct_adj, 
                      methods = c("ndotproduct"))



test_that("addSpectSimil", {
  expect_error(addSpectSimil(NULL, struct_adj), "'x' is not a 'Spectra' object")
  expect_error(addSpectSimil(sps_sub, NULL), 
               "'am_structural' is not an 'AdjacencyMatrix' object")
  expect_true(validObject(spect_adj))
  expect_equal(assayNames(spect_adj), c("binary","group",
                                        "formula","mass","ndotproduct"))
  expect_equal(sum(assay(spect_adj, "ndotproduct"), na.rm = TRUE), 
               45.16721, tolerance = 5e-01)
  expect_equal(as.vector(assay(spect_adj, "ndotproduct")[7,1:10]), 
               c(0, 0, 0, 0.989282368, 0, 0, 1, 0.003863543, 0, 0))
  expect_equal(rownames(assay(spect_adj, "ndotproduct")),
               colnames(assay(spect_adj, "ndotproduct")))
  expect_equal(rownames(assay(spect_adj, "ndotproduct")),
               colnames(assay(struct_adj, "binary")))
  expect_true(is.numeric(assay(spect_adj, "ndotproduct")))
  
})


## test threshold() and combine()
args_thr <- list(filter = "ndotproduct > 0.1 | 
                 binary == 1 & is.na(ndotproduct)")

thr_thr <- threshold(spect_adj, type = "threshold", args = args_thr)

cons_adj <- combine(thr_thr, stat_adj_thr)


test_that("threshold and combine similarity", {
  expect_error(threshold(spect_adj, type = "top1", args = list(n = 2)), 
        "'am' 'structural' can only be thresholded by 'type = threshold'")
  expect_true(is(thr_thr, "AdjacencyMatrix"))
  expect_equal(assayNames(thr_thr), 
               c("binary","group","formula","mass","ndotproduct"))
  expect_true(is.numeric(assay(thr_thr, "binary")))
  expect_equal(rownames(assay(thr_thr, "binary")), 
               rownames(assay(spect_adj, "binary")))
  expect_equal(colnames(assay(thr_thr, "binary")), 
               colnames(assay(spect_adj, "binary")))
  expect_equal(sum(assay(thr_thr, "binary"), na.rm = TRUE), 98)
  expect_true(all(assay(thr_thr, "binary") %in% c(0, 1, NaN)))
  
  ## test combine
  expect_true(validObject(cons_adj))
  expect_equal(assayNames(cons_adj), 
               c("binary", "group", "formula", "mass", "ndotproduct",
                 "clr_coef", "aracne_coef", "pearson_coef", "pearson_pvalue", 
                 "spearman_coef", "spearman_pvalue","ggm_coef", "ggm_pvalue",
                 "consensus","combine_binary", "combine_group", 
                 "combine_formula",  "combine_mass", "combine_ndotproduct"))
  expect_equal(sum(assay(cons_adj, "combine_binary"), na.rm = TRUE), 4)
  expect_equal(as.vector(assay(cons_adj, "combine_binary")[1, 1:5]), 
               c(NA, 0, 0, 0, 0))
  expect_equal(as.vector(assay(cons_adj, "combine_ndotproduct")[15, 15:20]), 
               c(NA, "0.130802929034727", "", "", "", ""))
  expect_equal(dim(assay(cons_adj, "combine_binary")), c(36, 36))
  expect_equal(rownames(assay(cons_adj, "combine_ndotproduct")), 
               rownames(assay(spect_adj, "ndotproduct")))
  expect_equal(colnames(assay(cons_adj, "combine_ndotproduct")), 
               colnames(assay(spect_adj, "ndotproduct")))
  
})

## END unit test spectral similarity ##
