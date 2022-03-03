#' @name spectra_data
#' 
#' @title Spectra data to test addSpectralSimilarity
#' 
#' @description \code{spectra_data} contains one selected putative annotation of
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
#' data("x_test", package = "MetNet")
#' 
#' x_annotation <- x_test[,1:2]
#' 
#' x_annotation <- cbind(x_annotation,"database_mz" = NA,
#'                      "database_identifier" = NA,
#'                      "chemical_formula" = NA,
#'                      "smiles" = NA,
#'                      "inchi" = NA,
#'                      "inchikey" =  NA,
#'                      "metabolite_identification" = NA,
#'                      "fragmentations" = NA,
#'                      "modifications" = NA,
#'                      "charge" = NA,
#'                      "database" = NA)
#'
#'x1856 <- cbind(x_annotation["x1856", "mz"],
#'               x_annotation["x1856", "rt"],
#'               "database_mz" = 	308.2,
#'               "database_identifier" = "N-caffeoylspermidine",
#'               "chemical_formula" = "C16H25N3O3",
#'               "smiles" = "C=1(C=C(C(=CC1)O)O)/C=C/C(NCCCNCCCCN)=O",
#'               "inchi" = "InChI=1S/C16H25N3O3/c17-8-1-2-9-18-10-3-11-19-16(22)7-5-13-4-6-14(20)15(21)12-13/h4-7,12,18,20-21H,1-3,8-11,17H2,(H,19,22)/b7-5+",
#'               "inchikey" =  "AZSUJBAOTYNFDE-FNORWQNLSA-N",
#'               "metabolite_identification" = NA,
#'               "fragmentations" = NA,
#'               "modifications" = NA,
#'               "charge" = 1,
#'               "database" = NA)
#'
#'x_annotation[rownames(x_annotation) == "x1856",] <- x1856
#'x_annotation <- x_annotation[,-c(1:2)]
#'
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
NULL
