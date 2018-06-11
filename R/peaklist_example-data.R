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
#' intensity values. All row entries with retention time < 103 s and > 440 s 
#' were removed. Wntries with m/z values < 250 and > 1200 were removed 
#' as well as entries with m/z values between 510 and 600 to reduce the file 
#' size. 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL
