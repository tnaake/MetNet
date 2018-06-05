#' @name peaklist
#' @title Example data for \code{MetNet}: data input
#' @description The object \code{peaklist} is a \code{data.frame}, where rows 
#' are features and the columns are samples (starting with X001-180). 
#' @docType data
#' @usage peaklist
#' @return \code{data.frame}
#' @format \code{data.frame}
#' @source Internal peaklist from metabolite profiling of Nicotiana species 
#' after W+OS and MeJA treatment. The data was processed by \code{xcms} and 
#' \code{CAMERA} scripts. All unncessary information is removed, keeping only 
#' the columns "mz", "rt" and the respective columns containing the 
#' intensity values.
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL
