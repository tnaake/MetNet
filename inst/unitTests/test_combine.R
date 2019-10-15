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

## create structural network
struct_adj <- createStructuralAdjacency(mat_test,
        transformation = transformations, ppm = 5)

## create statistical network
stat_adj_l <- statistical(mat_test, 
                          model = c("clr", "aracne", "pearson", "spearman", "bayes"), 
                          correlation_adjust = "bonferroni")
stat_adj <- threshold(stat_adj_l, type = "top2", args = list(n = 10))

## START unit test combine ##
cons_adj <- combine(struct_adj, stat_adj)

test_combine <- function() {
    checkException(combine(NULL, stat_adj), msg = "not a list of length 2")
    checkException(combine(struct_adj, NULL), msg = "not a numeric matrix")
    checkException(combine(struct_adj,
        stat_adj, threshold = "a"), msg = "is not numeric")
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
    
    ## check rownames/colnames
    mock <- struct_adj
    mock_num <- mock[[1]]
    colnames(mock_num)[1] <- "foo"
    mock[[1]] <- mock_num
    checkException(combine(mock, stat_adj), msg = "are not identical to ")
    
    mock <- struct_adj
    mock_num <- mock[[1]]
    rownames(mock_num)[1] <- "foo"
    mock[[1]] <- mock_num
    checkException(combine(mock, stat_adj), msg = "are not identical to")
    
    mock <- stat_adj
    colnames(mock)[1] <- "foo"
    checkException(combine(struct_adj, mock), msg = "are not identical")
    
    mock <- stat_adj
    rownames(mock)[1] <- "foo"
    checkException(combine(struct_adj, mock), msg = "are not identical")
    
    ## check for structure
    l <- list(1, cons_adj[[2]])
    checkException(combine(structure=l, statistical=stat_adj), 
        msg = "not a numeric matrix")
    l <- list(cons_adj[[1]], 1)
    checkException(combine(structure=l, statistical=stat_adj), 
        msg = "not a character ")
}
## END unit test combine ##
