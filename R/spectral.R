#' @name spectral
#'
#' @aliases spectral
#'
#' @title Create an `AdjacencyMatrix` object containing assays of adjacency 
#' matrices from spectral similarity methods
#'
#' @description
#' The function `spectral` infers adjacency matrix topologies from
#' spectral similarity methods and returns matrices of these networks in an
#' `AdjacencyMatrix` object. 
#' The function includes functionality to calculate adjacency matrices based on 
#' spectral similarity methods included in the `Spectra` package.
#' It uses a `Spectra`-object, storing the ms2 information and using a previously 
#' created `AdjacencyMatrix` object from the type "structural" in order to 
#' perform the mapping on the ms1 features. 
#' The function returns an `AdjacencyMatrix` object of adjacency matrices that 
#' are defined by `methods`.
#'
#' @param x 
#' `Spectra` object that contains a unique "id" (see `spectraVariables()`), 
#' matching to the row-/colnames of the structural `AdjacencyMatrix` and storing
#' important information of MS2 data (i.e. mz and intensity).
#' 
#' @param am_structural 
#' `AdjacencyMatrix` of type "structural" that was created using matching MS1 
#' data of the same data set.
#'
#' @param methods 
#' `character` vector containing the methods that will be used. All methods can 
#' be used that are already implemented in the `RforMassSpectrometry` 
#' infrastructure, e.g. (`"ndotproduct"`(default), `"gnps"`). 
#' `methods` are then forwarded to `FUN` in `Spectra::compareSpectra()`.
#'
#' @param ... 
#' parameters passed to the functions  `combineSpectra` from the `Spectra` 
#' package (e.g. `MAPFUN`, `tolerance`, `ppm`, `type`, ...).
#'
#' @details
#' The function `spectral` includes functionality to calculate adjacency
#' matrices based on different spectral similarity measures provided by the 
#' `RforMassSpectrometry` infrastructure. 
#' `spectral` calls the different functions (`FUN` parameter in 
#' `Spectra::compareSpectra()`) as specified by `methods`. The default is the 
#' normalized dotproduct `"ndotproduct"`. Moreover, different parameters (e.g. 
#' `MAPFUN`, `tolerance`, `ppm`, `type`, ...) of the function 
#' `Spectra::compareSpectra()` are forwarded by `...`. `spectral` will create 
#' adjacency matrices using the specified methods and will return an 
#' `AdjacencyMatrix` containing the weighted adjacency matrices in the `assays` 
#' slot.
#'
#' It is very important that features IDs in the MS1 data (i.e. row/colnames of 
#' `am_structural`) are matching to the IDs of the respective MS2 data (i.e. 
#' `"id"` in the `spectraVariables()`). Also, the `Spectra` object `x` is 
#' required to have unique `"id"`s, meaning on representative spectrum per 
#' feature. If features store multiple spectra, spectra consolidation has to be 
#' performed first (e.g. using the function `Spectra::combineSpectra`).
#' 
#'
#' @return `AdjacencyMatrix` of type "spectral" containing the respective 
#' adjacency matrices in the `assay` slot as specified by `methods`
#'
#' @author Liesa Salzer
#'
#' @examples
#' ## 1. load ms1 data
#' ## 2. load ms2 data as spectra object 
#' ## 3. perform structural()
#' ## 4. perform spectral()
#'
#' @export
#' 
#' @import MsCoreUtils

spectral <- function(x, am_structural, methods = c("ndotproduct"), ...) {
  
  ## sanity checks
  if (!is(x, "Spectra")) 
    stop("'x' is not a 'Spectra' object")
  
  if (!is(am_structural, "AdjacencyMatrix")) 
    stop("'am_structural' is not an 'AdjacencyMatrix' object")
  
  if (!validObject(am_structural)) 
    stop("'am_structural' must be a valid 'AdjacencyMatrix' object")
  
  if(!"id" %in% spectraVariables(x)) 
    stop("Spectra does not contain id column!")
 
  if(!length(unique(x$id)) == length(x)) 
    stop("Spectra shall contain only unique entries!")
    
  l <- list()
  
  for(method in methods) {
    ## create empty slots of adjacency matrix based on feature names of MS1 
    ## data and where we can fill then similarity values of MS2 data
    simil = matrix(NA, nrow = nrow(am_structural), ncol = nrow(am_structural))
    rownames(simil) = rownames(am_structural)
    colnames(simil) = colnames(am_structural)
    
    
    adj_spec <- Spectra::compareSpectra(x,
                                        FUN = get(method),
                                        ...)
    
    colnames(adj_spec) <- x$id
    rownames(adj_spec) <- x$id
    
    simil[rownames(adj_spec), colnames(adj_spec)] <- adj_spec
    
    ## add object to l
    new_index <- length(l) + 1
    l[[new_index]] <- simil
    
    ## assign the name to the newly added entry
    names(l)[new_index] <- method
    
  }
  
  
  rD <- DataFrame(names = rownames(l[[1]]))
  rownames(rD) <- rownames(l[[1]])
  
  
  adj <- AdjacencyMatrix(l, rowData = rD, directed = FALSE, 
                         thresholded = FALSE,
                         type = "spectral")
  
 return(adj)
  
  
}
