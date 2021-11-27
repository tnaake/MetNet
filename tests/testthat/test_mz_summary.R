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
struct_adj <- structural(mat_test, transformation = transformations, 
    var = c("group", "formula", "mass"), ppm = 5, directed = FALSE)
mz_summary(struct_adj, var = c("group", "formula"))

struct_adj_neg <- structural(mat_test, transformation = transformations_neg, 
    var = c("group", "formula", "mass"), ppm = 5, directed = FALSE)

struct_adj_neg_dir <- structural(mat_test, transformation = transformations_neg, 
    var = c("group", "formula", "mass"), ppm = 5, directed = TRUE)

## statistical calculation
stat_adj <- statistical(mat_test[, -1],
    model = c("clr", "aracne", "pearson", "spearman"))
stat_adj_thr <- threshold(stat_adj, type = "top2", args = list(n = 10))

## combine
cons_adj <- combine(struct_adj, stat_adj_thr)

## START unit test mz_summary##
summary_struct_adj <- mz_summary(struct_adj, var = c("group", "formula"))
summary_struct_adj_neg_dir <- mz_summary(struct_adj_neg_dir,
    var = c("group", "formula"))
summary_cons_adj <- mz_summary(cons_adj, var = c("group", "formula"))


test_that("mz_summary", {
    expect_error(mz_summary(transformations, var = c("group", "formula")),
        "'am' is not an 'AdjacencyMatrix' object")
    expect_error(mz_summary(struct_adj_neg, var = c("group", "formula")),
               "assay 'binary' does not contain any mass differences")
    expect_error(mz_summary(stat_adj, var = c("group", "formula")),
               "'am' is not of type 'structural' or 'combine'")
    expect_error(mz_summary(struct_adj, var = c("group", "formula"), 
        filter = TRUE), "'filter' needs to be numeric")
    expect_error(mz_summary(struct_adj, filter = -1),
        "'filter' needs to be 0 or positive numeric")
    expect_equal(colnames(summary_struct_adj), 
        c("group",  "formula", "count"))
    expect_equal(length(summary_struct_adj), 3)
    expect_equal(dim(summary_struct_adj), c(2, 3))
    expect_equal(summary_struct_adj$group,
        c("Malonyl group (-H2O)", "Monosaccharide (-H2O)"))
    expect_equal(summary_struct_adj$formula, c("C3H2O3", "C6H10O5"))
    expect_equal(summary_struct_adj$count, c(3,3))
    expect_true(is.data.frame(summary_struct_adj))
    expect_equal(length(summary_struct_adj_neg_dir), 3)
    expect_equal(dim(summary_struct_adj_neg_dir), c(2, 3))
    expect_equal(summary_struct_adj_neg_dir$group,
        c("Malonyl group (-H2O)", "Monosaccharide (-H2O)"))
    expect_equal(summary_struct_adj_neg_dir$formula, c("C3H2O3", "C6H10O5"))
    expect_equal(summary_struct_adj_neg_dir$count, c(3,3))
    expect_equal(length(summary_cons_adj), 3)
    expect_equal(dim(summary_cons_adj), c(2, 3))
    expect_equal(summary_cons_adj$group,
        c("Malonyl group (-H2O)", "Monosaccharide (-H2O)"))
    expect_equal(summary_cons_adj$formula, c("C3H2O3", "C6H10O5"))
    expect_equal(summary_cons_adj$count, c(3,3))
    expect_true(is.data.frame(summary_cons_adj))
    expect_error(mz_summary(struct_adj, var = "foo"), "assay 'foo' not in 'am'")
    expect_error(mz_summary(struct_adj, var = NULL), 
        "'var' has to be a character of length > 0")
    expect_error(mz_summary(struct_adj, var = character()), 
        "'var' has to be a character of length > 0")
    expect_error(mz_summary(struct_adj, var = "group", filter = -1), 
        "'filter' needs to be 0 or positive numeric")
    expect_error(mz_summary(struct_adj, var = "group", filter = NULL), 
        "'filter' needs to be numeric")
})
## END unit test mz_summary ##


## START unit test mz_vis ##

## create mass difference count test file missing count
transformation <- c("Malonyl group (-H2O)", "Monosaccharide (-H2O)")
mass_difference <- c("162.0528234315", "86.0003939305")
sum_test <- cbind(transformation, mass_difference)
sum_test_df <- as.data.frame(sum_test)

vis <- mz_vis(summary_struct_adj, var = "group")

test_that("mz_vis", {
  expect_error(mz_vis(sum_test), "'df' is not a data.frame")
  expect_error(mz_vis(sum_test_df), "'df' does not contain the column 'count'")
  expect_true(ggplot2::is.ggplot(vis))
  expect_identical(vis$labels$title, "Numbers of determined mass differences")
  expect_identical(vis$labels$y, "count")
  expect_identical(vis$labels$x, "group")
  expect_equal(vis$data, summary_struct_adj)
})
## END unit test mz_vis ##
