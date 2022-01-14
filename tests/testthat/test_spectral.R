## create toy example data set
data("x_test", package = "MetNet")
x_test <- as.matrix(x_test)

## transformations object for structural calculation
transformations <- rbind(
  c("Malonyl group (-H2O)", "C3H2O3", "86.0003939305"),
  c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"))
transformations <- data.frame(group = as.character(transformations[, 1]),
                              formula = as.character(transformations[, 2]),
                              mass = as.numeric(transformations[, 3]))

## create structural network
struct_adj <- structural(x_test, transformation = transformations, 
                         var = c("group", "formula", "mass"), ppm = 5)

## load MS2 
data("ms2_test", package = "MetNet")

## START unit test spectral ##
spect_adj <- spectral(sps_sub, am_structural = struct_adj, 
                      methods = c("ndotproduct", "gnps"))

## create structural with different parameters from compareSpectra()
spect_adj_param <- spectral(sps_sub, am_structural = struct_adj, 
                      methods = c("gnps"), 
                      MAPFUN = joinPeaksGnps, 
                      type = "inner", 
                      wheighted = TRUE,
                      ppm = 10)

test_that("spectral", {
  expect_error(spectral(NULL, struct_adj), "'x' is not a 'Spectra' object")
  expect_error(spectral(sps_sub, NULL), 
               "'am_structural' is not an 'AdjacencyMatrix' object")
  expect_true(validObject(spect_adj))
  expect_equal(assayNames(spect_adj), c("ndotproduct", "gnps"))
  expect_equal(sum(assay(spect_adj, "ndotproduct"), na.rm = TRUE), 
               45.16721, tolerance = 5e-01)
  expect_equal(sum(assay(spect_adj, "gnps"), na.rm = TRUE), 
               55.04219, tolerance = 5e-01)
  expect_equal(as.vector(assay(spect_adj, "ndotproduct")[7,1:10]), 
               c(0, 0, 0, 0.989282368, 0, 0, 1, 0.003863543, 0, 0))
  expect_equal(sum(assay(spect_adj_param, "gnps"), na.rm = TRUE), 
               592.4462, tolerance = 5e-01)
  expect_equal(rownames(assay(spect_adj, "ndotproduct")),
               colnames(assay(spect_adj, "ndotproduct")))
  expect_equal(rownames(assay(spect_adj, "gnps")),
               colnames(assay(spect_adj, "gnps")))
  expect_equal(rownames(assay(spect_adj, "ndotproduct")),
               colnames(assay(struct_adj, "binary")))
  expect_true(is.numeric(assay(spect_adj, "ndotproduct")))
  
  
})
## END unit test spectral ##















