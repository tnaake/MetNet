## create toy example data set
data("mat_test", package = "MetNet")
colnames(mat_test) <- paste0("intensity_", 1:20)
mz <- c(100, 150, 200, 200, 262.0528, 348.0532, 448.0532)
mat_test <- cbind(mz = mz, mat_test)

## transformations object for structual calculation
transformations <- rbind(
    c("Malonyl group (–H2O)", "C3H2O3", "86.0003939305"),
    c("Monosaccharide (–H2O)", "C6H10O5", "162.0528234315"))
transformations <- data.frame(group = as.character(transformations[, 1]),
                                formula = as.character(transformations[, 2]),
                                mass = as.numeric(transformations[, 3]))

## create statistical network
stat_adj <- createStatisticalAdjacency(mat_test[, -1],
        model = c("clr", "aracne", "pearson", "spearman", "bayes"))
## create structural network
struct_adj <- createStructuralAdjacency(mat_test,
        transformation = transformations, ppm = 5)

## START unit test combine ##
cons_adj <- combine(struct_adj, stat_adj)

test_combine <- function() {
    checkException(combine(NULL, stat_adj))
    checkException(combine(struct_adj, NULL))
    checkException(combine(struct_adj,
        stat_adj, threshold = "a"))
    checkEquals(sum(cons_adj[[1]]), 4)
    checkEquals(dim(cons_adj[[1]]), c(7, 7))
    checkEquals(rownames(cons_adj[[1]]), paste0("x", 1:7))
    checkEquals(rownames(cons_adj[[2]]), paste0("x", 1:7))
    checkEquals(colnames(cons_adj[[1]]), paste0("x", 1:7))
    checkEquals(colnames(cons_adj[[2]]), paste0("x", 1:7))
    checkEquals(rownames(cons_adj[[1]]), colnames(cons_adj[[1]]))
    checkEquals(rownames(cons_adj[[2]]), colnames(cons_adj[[2]]))
    checkTrue(is.matrix(cons_adj[[1]]))
    checkTrue(is.matrix(cons_adj[[2]]))
    checkTrue(is.numeric(cons_adj[[1]]))
    checkTrue(is.character(cons_adj[[2]]))
    
    mock <- stat_adj
    colnames(mock)[1] <- "foo"
    checkException(combine(struct_adj, mock))
    
    mock <- stat_adj
    rownames(mock)[1] <- "foo"
    checkException(combine(struct_adj, mock))
    
    ## check for structure
    l <- list(1, cons_adj[[2]])
    checkException(combine(structure=l, statistical=stat_adj))
    l <- list(cons_adj[[1]], 1)
    checkException(combine(structure=l, statistical=stat_adj))
}
## END unit test combine ##
