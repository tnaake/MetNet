#' @name spectra_matrix
#' 
#' @title Spectra data to test addSpectralSimilarity
#' 
#' @description 
#' \code{spectra_matrix}  contains one selected putative 
#' annotation of \code{x_test}. Missing annotations are filled with `NA`'s. 
#' It will be used as an example annotation in the vignette to show the 
#' functionality of the package.
#' 
#' @docType data
#' 
#' @return \code{matrix}
#' 
#' @format \code{matrix}
#' 
#' @source
#' library(MsCoreUtils)
#' library(Spectra)
#' 
#' f <- system.file("spectra_matrix/ms2_test.RDS", package = "MetNet")
#' sps_sub <- readRDS(f)
#' 
#' adj_spec <- Spectra::compareSpectra(sps_sub, FUN = ndotproduct)
#' colnames(adj_spec) <- sps_sub$id
#' rownames(adj_spec) <- sps_sub$id
#' 
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
NULL
