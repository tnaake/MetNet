## create toy example data set
data("mat_test", package="MetNet")
colnames(mat_test) <- paste0("intensity_", 1:20)
mz <- c(100, 150, 200, 200, 262.0528, 348.0532, 448.0532)
mat_test <- cbind(mz=mz, mat_test)

## functional_groups object for structual calculation
functional_groups <- rbind(
    c("Malonyl group (–H2O)", "C3H2O3", "86.0003939305"),
    c("Monosaccharide (–H2O)", "C6H10O5", "162.0528234315"))
functional_groups <- data.frame(group=as.character(functional_groups[, 1]),
                            formula=as.character(functional_groups[, 2]),
                            mass=as.numeric(functional_groups[, 3]))

## START unit test inRangeWhich ##
test_inRangeWhich <- function() {
    checkException(inRangeWhich(m_1=NULL, m_2=100, 
                                   functional_groups=99.9))
    checkException(inRangeWhich(m_1=99.8, m_2=NULL, 
                                   functional_groups=99.9))
    checkException(inRangeWhich(m_1=99.8, m_2=100, 
                                   functional_groups=NULL))
    checkEquals(inRangeWhich(m_1=86.0002, m_2=86.0004, 
        functional_groups=functional_groups[, "mass"]), 1)
    checkEquals(inRangeWhich(m_1=86.0004, m_2=86.0002, 
        functional_groups=functional_groups[, "mass"]), 1)
    checkEquals(inRangeWhich(m_1=162.0527, m_2=162.0529, 
        functional_groups=functional_groups[, "mass"]), 2)
    checkEquals(inRangeWhich(m_1=162.0529, m_2=162.0527, 
        functional_groups=functional_groups[, "mass"]), 2)
    checkEquals(inRangeWhich(m_1=86.0002, m_2=162.0529, 
        functional_groups=functional_groups[, "mass"]), c(1, 2))
}
## END unit test inRangeWhich ## 

## START unit test createStructuralAdjacency ##
struct_adj <- createStructuralAdjacency(mat_test, 
        functional_groups=functional_groups, ppm=5)
test_createStructuralAdjacency <- function() {
    checkException(createStructuralAdjacency(mat_test[, -1], functional_groups))
    checkException(createStructuralAdjacency(NULL, functional_groups))
    checkException(createStructuralAdjacency(mat_test, functional_groups[,-1]))
    checkException(createStructuralAdjacency(mat_test, functional_groups[,-3]))
    checkException(createStructuralAdjacency(mat_test, matrix()))
    checkException(createStructuralAdjacency(mat_test, functional_groups, ppm="a"))
    checkEquals(length(struct_adj), 2)
    checkEquals(dim(struct_adj[[1]]), c(7, 7))
    checkEquals(dim(struct_adj[[2]]), c(7, 7))
    checkEquals(rownames(struct_adj[[1]]), colnames(struct_adj[[1]]))
    checkEquals(rownames(struct_adj[[2]]), colnames(struct_adj[[2]]))
    checkEquals(rownames(struct_adj[[1]]), rownames(struct_adj[[2]]))
    checkEquals(rownames(struct_adj[[1]]), paste0("x", 1:7))
    checkEquals(sum(struct_adj[[1]]), 4)
    checkEquals(unique(as.vector(struct_adj[[2]])), 
                c("", "Monosaccharide (–H2O)", "Malonyl group (–H2O)" ))
    checkTrue(is.matrix(struct_adj[[1]]))
    checkTrue(is.matrix(struct_adj[[2]]))
    checkTrue(is.numeric(struct_adj[[1]]))
    checkTrue(is.character(struct_adj[[2]]))
}
## END unit test createStructuralAdjacency ## 
