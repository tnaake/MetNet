## load RUnit
library(RUnit)

## create toy example data set
data("mat_test", package = "MetNet")
colnames(mat_test) <- paste0("intensity_", 1:20)
mz <- c(100, 150, 262.0528, 262.0528, 262.0528, 348.0532, 448.0532)
rt <- c(100, 100, 50, 150, 150, 150, 150)
mat_test <- cbind(mz = mz, rt = rt, mat_test)

## transformations object for structual calculation
transformations <- rbind(
    c("Malonyl group (–H2O)", "C3H2O3", 86.0003939305, "?"),
    c("Monosaccharide (–H2O)", "C6H10O5", 162.0528234315, "-"))
transformations <- data.frame(group = transformations[, 1],
                            formula = transformations[, 2],
                            mass = as.numeric(transformations[, 3]),
                            rt = transformations[, 4])

## START unit test structural ##
struct_adj <- structural(mat_test,
        transformation = transformations, ppm = 5)
test_structural <- function() {
    checkException(structural(mat_test[, -1], transformations))
    checkException(structural(NULL, transformations))
    checkException(structural(mat_test, transformations[,-1]))
    checkException(structural(mat_test, transformations[,-3]))
    checkException(structural(mat_test, matrix()))
    checkException(structural(mat_test, transformations,
                                             ppm = "a"))
    checkEquals(length(struct_adj), 2)
    checkEquals(dim(struct_adj[[1]]), c(7, 7))
    checkEquals(dim(struct_adj[[2]]), c(7, 7))
    checkEquals(rownames(struct_adj[[1]]), colnames(struct_adj[[1]]))
    checkEquals(rownames(struct_adj[[2]]), colnames(struct_adj[[2]]))
    checkEquals(rownames(struct_adj[[1]]), rownames(struct_adj[[2]]))
    checkEquals(rownames(struct_adj[[1]]), paste0("x", 1:7))
    checkEquals(sum(struct_adj[[1]]), 12)
    checkEquals(unique(as.vector(struct_adj[[2]])), 
                c("", "Monosaccharide (–H2O)", "Malonyl group (–H2O)" ))
    checkTrue(is.matrix(struct_adj[[1]]))
    checkTrue(is.matrix(struct_adj[[2]]))
    checkTrue(is.numeric(struct_adj[[1]]))
    checkTrue(is.character(struct_adj[[2]]))
}
## END unit test structural ##

## START unit test rtCorrection ##
struct_adj_rt <- rtCorrection(struct_adj, mat_test, transformations)
test_rtCorrection <- function() {
    checkException(rtCorrection(struct_adj[[1]], mat_test, transformations))
    checkException(
        rtCorrection(list(struct_adj[[1]]), mat_test, transformations))
    checkException(rtCorrection(list(struct_adj[[2]], struct_adj[[2]]),
        mat_test, transformations))
    checkException(rtCorrection(list(struct_adj[[1]], struct_adj[[1]]),
        mat_test, transformations))
    checkException(rtCorrection(struct_adj, NULL, transformations))
    checkException(rtCorrection(struct_adj, mat_test[,-1], transformations))
    checkException(rtCorrection(struct_adj, mat_test[,-2], transformations))
    checkException(rtCorrection(struct_adj, mat_test, NULL))
    checkException(rtCorrection(struct_adj, mat_test, transformations[, -1]))
    checkException(rtCorrection(struct_adj, mat_test, transformations[, -3]))
    checkException(rtCorrection(struct_adj, mat_test, transformations[, -4]))
    checkException(rtCorrection(struct_adj, mat_test,
        cbind(transformations[,-4], rt = rep("a", 4))))
    checkTrue(is.matrix(struct_adj_rt[[1]]))
    checkTrue(is.numeric(struct_adj_rt[[1]]))
    checkTrue(is.matrix(struct_adj_rt[[2]]))
    checkTrue(is.character(struct_adj_rt[[2]]))
    checkEquals(colnames(struct_adj_rt[[1]]), paste0("x", 1:7))
    checkEquals(colnames(struct_adj_rt[[1]]), rownames(struct_adj_rt[[1]]))
    checkEquals(colnames(struct_adj_rt[[1]]), colnames(struct_adj_rt[[2]]))
    checkEquals(colnames(struct_adj_rt[[1]]), rownames(struct_adj_rt[[2]]))
    checkTrue(table(struct_adj_rt)[1] == 41)
    checkTrue(table(struct_adj_rt[[1]])[1] == 41)
    
    ## dims of struct_adj[[1]] and struct[[2]], rownames/colnames
    foo_1 <- struct_adj[[1]]
    foo_2 <- struct_adj[[2]]
    foo <- list(foo_1, foo_2[,1:6])
    checkException(rtCorrection(foo, mat_test, transformations), 
        msg = "dim(structural[[1]] is not equal to dim(structural[[2]])")
    
    rownames(foo_1)[1] <- "foo"
    foo <- list(foo_1, foo_2)
    checkException(rtCorrection(foo, mat_test, transformations), 
        msg = "colnames of structural[[1]] are not identical to rownames of ")
    
    foo_1 <- struct_adj[[1]]
    rownames(foo_2)[1] <- "foo"
    foo <- list(foo_1, foo_2)
    checkException(rtCorrection(foo, mat_test, transformations), 
        msg = "colnames of structural[[2]] are not identical to rownames of ")
    
    foo_2 <- struct_adj[[2]]
    foo <- list(foo_1, foo_2)
    checkException(rtCorrection(foo, mat_test[1:6,], transformations), 
                   msg = "rownames(structural[[1]]) do not fit rownames(x)")
    
    foo_1 <- as.data.frame(struct_adj[[1]])
    foo <- list(foo_1, foo_2)
    checkException(rtCorrection(foo, mat_test, transformations), 
                   msg = "structural[[1]]) is not a numeric matrix")
    
    foo_1 <- struct_adj[[1]]
    foo_2 <- as.data.frame(struct_adj[[2]])
    foo <- list(foo_1, foo_2)
    checkException(rtCorrection(foo, mat_test, transformations), 
                   msg = "structural[[2]]) is not a character matrix")
}
## END unit test rtCorrection ##
