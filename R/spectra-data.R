#' @name spectra_matrix
#' 
#' @title Spectra data to test addSpectralSimilarity
#' 
#' @description \code{spectra_matrix} is 
#' contains one selected putative annotation of
#' \code{x_test}. Missing annotations are filled with `NA`'s. It will be used as an example annotation in the vignette to
#' show the functionality of the packages.
#' 
#' @docType data
#' 
#' @return \code{matrix}
#' 
#' @format \code{matrix}
#' 
#' @source
#' 
#' library(MsCoreUtils)
#' library(Spectra)
#' 
#' data("ms2_test", package = "MetNet")
#' 
#' adj_spec <- Spectra::compareSpectra(sps_sub,
#'                                     FUN = ndotproduct)
#' colnames(adj_spec) <- sps_sub$id
#' rownames(adj_spec) <- sps_sub$id
#' 
#'
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
NULL
