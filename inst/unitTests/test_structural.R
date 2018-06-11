## create toy example data set
set.seed(1)
mat_test <- matrix(rnorm(140, mean = 10, sd = 2), nrow = 7)
mat_test[1:3, ] <- t(apply(mat_test[1:3, ], 1, sort))
mat_test[5:7, ] <- t(apply(mat_test[5:7, ], 1, sort, decreasing = TRUE))
rownames(mat_test) <- paste("x", 1:7, sep = "")
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

## START unit test in_range_which ##
test_in_range_which <- function() {
    checkException(in_range_which(m_1 = NULL, m_2 = 100, 
                                   functional_groups = 99.9))
    checkException(in_range_which(m_1 = 99.8, m_2 = NULL, 
                                   functional_groups = 99.9))
    checkException(in_range_which(m_1 = 99.8, m_2 = 100, 
                                   functional_groups = NULL))
    checkEquals(in_range_which(m_1 = 86.0002, m_2 = 86.0004, 
        functional_groups = functional_groups[, "mass"]), 1)
    checkEquals(in_range_which(m_1 = 86.0004, m_2 = 86.0002, 
        functional_groups = functional_groups[, "mass"]), 1)
    checkEquals(in_range_which(m_1 = 162.0527, m_2 = 162.0529, 
        functional_groups = functional_groups[, "mass"]), 2)
    checkEquals(in_range_which(m_1 = 162.0529, m_2 = 162.0527, 
        functional_groups = functional_groups[, "mass"]), 2)
    checkEquals(in_range_which(m_1 = 86.0002, m_2 = 162.0529, 
        functional_groups = functional_groups[, "mass"]), c(1, 2))
}
## END unit test in_range_which ## 

## START unit test create_structural_network ##
struct_net <- create_structural_network(mat_test, 
        functional_groups = functional_groups, ppm = 5)
test_create_structural_network <- function() {
    checkException(create_structural_network(mat_test[, -1], functional_groups))
    checkException(create_structural_network(NULL, functional_groups))
    checkException(create_structural_network(mat_test, functional_groups[,-1]))
    checkException(create_structural_network(mat_test, functional_groups[,-2]))
    checkException(create_structural_network(mat_test, functional_groups[,-3]))
    checkException(create_structural_network(mat_test, matrix()))
    checkException(create_structural_network(mat_test, functional_groups, ppm = "a"))
    checkEquals(length(struct_net), 2)
    checkEquals(dim(struct_net[[1]]), c(7, 7))
    checkEquals(dim(struct_net[[2]]), c(7, 7))
    checkEquals(rownames(struct_net[[1]]), colnames(struct_net[[1]]))
    checkEquals(rownames(struct_net[[2]]), colnames(struct_net[[2]]))
    checkEquals(rownames(struct_net[[1]]), rownames(struct_net[[2]]))
    checkEquals(rownames(struct_net[[1]]), paste0("x", 1:7))
    checkEquals(sum(struct_net[[1]]), 4)
    checkEquals(unique(as.vector(struct_net[[2]])), 
                c("", "Monosaccharide (–H2O)", "Malonyl group (–H2O)" ))
    checkTrue(is.matrix(struct_net[[1]]))
    checkTrue(is.matrix(struct_net[[2]]))
    checkTrue(is.numeric(struct_net[[1]]))
    checkTrue(is.character(struct_net[[2]]))
}
## END unit test create_structural_network ## 
