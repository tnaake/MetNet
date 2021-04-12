## create toy example data set
data("mat_test", package = "MetNet")
colnames(mat_test) <- paste0("intensity_", 1:20)
mz <- c(100, 150, 262.0528, 262.0528, 262.0528, 348.0532, 448.0532)
rt <- c(100, 100, 50, 150, 150, 150, 150)
mat_test <- cbind(mz = mz, rt = rt, mat_test)

## transformations object for structural calculation
transformations <- rbind(
  c("Malonyl group (-H2O)", "C3H2O3", 86.0003939305, "+"),
  c("Monosaccharide (-H2O)", "C6H10O5", 162.0528234315, "-"))
transformations_neg <- transformations <- data.frame(
  group = transformations[, 1],
  formula = transformations[, 2],
  mass = as.numeric(transformations[, 3]),
  rt = transformations[, 4])
transformations_neg[, 3] <- -1 * transformations_neg[, 3]


## structural calculation
struct_adj <- structural(mat_test,
                         transformation = transformations, ppm = 5, directed = FALSE)

struct_adj_neg <- structural(mat_test,
                             transformation = transformations_neg, ppm = 5, directed = FALSE)

## statistical calculation
stat_adj <- statistical(mat_test[, -1],
                        model = c("clr", "aracne", "pearson", "spearman", "bayes"))
stat_adj_thr <- threshold(stat_adj, type = "top2", args = list(n = 10))

## combine
cons_adj <- combine(struct_adj, stat_adj_thr)


## START unit test mz_summary##

summary_struct_adj <- mz_summary(struct_adj)
summary_cons_adj <- mz_summary(cons_adj)


test_that("mz_summary", {
  expect_error(mz_summary(transformations),
               "'x' is not an 'AdjacencyMatrix' object")
  expect_error(mz_summary(struct_adj_neg),
               "'x' does not contain any mass-differences")
  expect_error(mz_summary(stat_adj),
               "'x' does not contain a mass-difference adjacency matrix")
  expect_error(mz_summary(struct_adj, filter = TRUE),
               "'filter' needs to be numeric or 'FALSE'")
  expect_equal(colnames(summary_struct_adj), 
               c("transformation",  "mass_difference", "counts"))
  expect_equal(length(summary_struct_adj), 3)
  expect_equal(dim(summary_struct_adj), c(2, 3))
  expect_equal(summary_struct_adj$transformation,
               c("Malonyl group (-H2O)", "Monosaccharide (-H2O)"))
  expect_equal(summary_struct_adj$mass_difference,
               c("162.0528234315", "86.0003939305" ))
  expect_equal(summary_struct_adj$counts,
               c(3,3))
  expect_true(is.data.frame(summary_struct_adj))
  expect_equal(length(summary_cons_adj), 3)
  expect_equal(dim(summary_cons_adj), c(2, 3))
  expect_equal(summary_cons_adj$transformation,
               c("Malonyl group (-H2O)", "Monosaccharide (-H2O)"))
  expect_equal(summary_cons_adj$mass_difference,
               c("162.0528234315", "86.0003939305" ))
  expect_equal(summary_cons_adj$counts,
               c(4,4))
  expect_true(is.data.frame(summary_cons_adj))
})
## END unit test mz_summary ##
