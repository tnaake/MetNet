#' @name mat_test
#' @title Example data for \code{MetNet}: unit tests
#' @description \code{mat_test} contains 7 toy features that were derived
#' from \code{rnorm}. It will be used as an example data set in unit tests.
#' @docType data
#' @usage mat_test
#' @return \code{matrix}
#' @format \code{matrix}
#' @source 
#' set.seed(1)
#' random_numbers <- rnorm(140, mean = 10, sd = 2)
#' mat_test <- matrix(random_numbers, nrow = 7)
#' mat_test[1:3, ] <- t(apply(mat_test[1:3, ], 1, sort))
#' mat_test[5:7, ] <- t(apply(mat_test[5:7, ], 1, sort, decreasing = TRUE))
#' rownames(mat_test) <- paste("x", 1:7, sep = "")
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL

#' @name mat_test_z
#' @title Example data for \code{MetNet}: unit tests
#' @description \code{mat_test_z} contains 7 toy features that were derived
#' from \code{rnorm}. It will be used as an example data set in unit tests.
#' @docType data
#' @usage mat_test_z
#' @return \code{matrix}
#' @format \code{matrix}
#' @source 
#' set.seed(1)
#' random_numbers <- rnorm(140, mean = 10, sd = 2)
#' mat_test <- matrix(random_numbers, nrow = 7)
#' mat_test[1:3, ] <- t(apply(mat_test[1:3, ], 1, sort))
#' mat_test[5:7, ] <- t(apply(mat_test[5:7, ], 1, sort, decreasing = TRUE))
#' rownames(mat_test) <- paste("x", 1:7, sep = "")
#' mat_test_z <- apply(mat_test, 1, function(x) 
#'                     (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL
