## create toy example data set
data("mat_test", package = "MetNet")
colnames(mat_test) <- paste0("intensity_", 1:20)
mz <- c(100, 150, 200, 200, 262.0528, 348.0532, 448.0532)
mat_test <- cbind(mz = mz, mat_test)

## functional_groups object for structual calculation
functional_groups <- rbind(
    c("Malonyl group (–H2O)", "C3H2O3", "86.0003939305"),
    c("Monosaccharide (–H2O)", "C6H10O5", "162.0528234315"))
functional_groups <- data.frame(group = as.character(functional_groups[, 1]),
                                formula = as.character(functional_groups[, 2]),
                                mass = as.numeric(functional_groups[, 3]))

## create statistical network
stat_net <- create_statistical_network(mat_test[, -1], 
        model = c("clr", "aracne","pearson", "spearman", "bayes"))
## create structural network
struct_net <- create_structural_network(mat_test, 
        functional_groups = functional_groups, ppm = 5)

## START unit test combine_structural_statistical ## 
cons_net <- combine_structural_statistical(struct_net[[1]], stat_net)
test_combine_structural_statistical <- function() {
    checkException(combine_structural_statistical(NULL, stat_net))
    checkException(combine_structural_statistical(struct_net[[1]], NULL))
    checkException(combine_structural_statistical(struct_net[[1]], 
        stat_net, threshold = "a"))
    checkEquals(sum(cons_net), 4)
    checkEquals(dim(cons_net), c(7, 7))
    checkEquals(rownames(cons_net), paste0("x", 1:7))
    checkEquals(colnames(cons_net), paste0("x", 1:7))
    checkEquals(rownames(cons_net), colnames(cons_net))
    checkTrue(is.matrix(cons_net))
    checkTrue(is.numeric(cons_net))
}
## END unit test combine_structural_statistical ## 
