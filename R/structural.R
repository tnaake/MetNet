#' @name structural
#'
#' @aliases structural
#'
#' @title Create adjacency matrix based on m/z (molecular weight) difference
#'
#' @description
#' The function `structural` infers an unweighted
#' adjacency matrix using differences in m/z values that are matched against a
#' `data.frame` (`transformation` of calculated theoretical differences of
#' loss/addition of functional groups. `structural` returns
#' an `AdjacencyMatrix` object containing
#' the unweighted `numeric` `matrix` (assay `binary`), together with one or 
#' multiple `character` matrices containing e.g. the type of loss/addition 
#' or the mass differences. The creation of the additional `character` matrices
#' is controlled by the `var` argument that specifies the column in 
#' `transformation` as the data source for the adjacency matrices.
#' 
#' @param
#' x `matrix` or `data.frame`, where columns are the samples and the rows are 
#' features (metabolites), cell entries are intensity values. `x` contains the
#' column `"mz"` that has the m/z information (numerical values) for the
#' calculation of mass differences between features
#'
#' @param transformation 
#' `data.frame`, containing the columns `"group"`,
#' and `"mass"` that will be used for detection of transformation of
#' (functional) groups
#' 
#' @param var
#' `character` corresponding to column names in `transformation`
#'
#' @param ppm 
#' `numeric(1)`, mass accuracy of m/z features in parts per million (ppm)
#' 
#' @param directed 
#' `logical(1)`, if `TRUE`, absolute values of m/z differences will be
#' taken to query against `transformation`  (irrespective the sign of `mass`)
#' and undirected adjacency matrices will be returned as the respective
#' assays. This means, if there is a negative mass in 
#' `transformation[, "mass"]`, this negative mass will not be reported. 
#' If `FALSE`, directed
#' adjacency matrices will be returned with links reported that match the
#' transformations defined in `transformation` (respecting the sign of `mass`).
#' The `directed` slot of the returned `AdjacencyMatrix` object will contain
#' the information on `directed`.
#' 
#' @details
#' `structural` accesses the column `"mz"` of
#' `x` to infer structural topologies based on the functional groups
#' defined by `transformation`. To account for the mass accuracy of
#' the dataset `x`, the user can specify the accuracy of m/z features
#' in parts per million (ppm) by the `ppm` argument. The m/z values in the
#' `"mz"` column of `x`" will be converted to m/z ranges according to
#' the `ppm` argument (default `ppm = 5`).
#' 
#' The returned `AdjacencyMatrix` object contains the assays 
#' `binary` and additional adjacency matrices depending on the `var` 
#' parameter. The assay `binary` stores the `numeric`
#' `matrix` with binary edges inferred from mass differences. The `var` 
#' parameter has to be set according to the column names in `transformation`.
#' E.g. if the `transformation` object contains a `"group"` column that stores
#' the name of the transformation, setting `var = "group"` will create an
#' assay `"group"` that contains the adjacency matrices where the entries
#' correspond to the `"group"` information for the respective feature pairs.
#'
#' The `type` slot is set to `structural`. The `directed` slot is set 
#' accordingly to the `directed` argument of the function `structural`.
#' The `thresholded` slot is set to `FALSE`.
#'
#' @return
#' `AdjacencyMatrix` object. The object will store the adjacency matrix/matrices
#' in the assay slot/slots. The numerical (unweighted) adjacency matrix
#' inferred from mass differences is stored as the assay `"binary"`. Depending
#' on the `var` argument, there are additional adjacency matrices stored
#' in the assay slot.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com} and
#' Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' transformation <- rbind(
#'     c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"),
#'     c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851"),
#'     c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945"))
#' transformation <- data.frame(group = transformation[, 1],
#'                                 formula = transformation[, 2],
#'                                 mass = as.numeric(transformation[, 3]))
#' am_struct <- structural(x_test, transformation, var = c("group", "mass"),
#'     ppm = 10, directed = TRUE)
#'
#' @export
structural <- function(x, transformation, var = character(), 
    ppm = 5, directed = FALSE) {

    ## check for integrity of x
    if (!(is.matrix(x) | is.data.frame(x)))
        stop("'x' has to be a matrix or data.frame")
    if (!"mz" %in% colnames(x)) 
        stop("'x' does not contain the column mz")
    
    ## check for integrity of transformation
    if (!is.data.frame(transformation))
        stop("'transformation' is not a data.frame")
    if (!"mass" %in% colnames(transformation))
        stop("'transformation' does not contain the column mass")
    if (!is.character(var)) 
        stop("'var' is not a character vector")
    var_err <- var[!var %in% colnames(transformation)]
    if (length(var_err) > 0)
        stop(sprintf("'transformation' does not contain the column '%s'", 
                     paste(var_err, collapse = "', '")))

    ## check for integrity of ppm
    if (!is.numeric(ppm) | length(ppm) != 1) 
        stop("'ppm' has to be a numeric of length 1")
    if (ppm <= 0) 
        stop("'ppm' has to be a positive value")
    
    ## check for integrity of directed
    if (!is.logical(directed) | length(directed) != 1)
        stop("'directed' has to be a logical vector of length 1")

    mass <- x[, "mz"]
    mat <- matrix(0, nrow = length(mass), ncol = length(mass))
    rownames(mat) <- colnames(mat) <- mass

    ## create matrix which has rownames per row
    mat <- apply(mat, 1, function(x) as.numeric(mass))
    
    ## receive indices of lower triangle
    lt <- lower.tri(mat)

    ## outline of the function:
    ## example: we have the two features M_1 and M_2, mz(M_1) > mz(M_2), 
    ## we calculate the ppm deviations from M_1 and M_2
    ## M_2+ppm, M_2, M_2-ppm
    ## M_1+ppm, M_1, M_1-ppm
    ## we denote as A = (M_2-ppm) - (M_1+ppm) and B = (M_2+ppm) - (M_1-ppm)
    ## we are interested in the distances A and B between each feature pair in 
    ## the data set and check then in the following if A <= transf <= B, 
    ## where transf is a queried transformation, ## e.g. glucose 
    ## transformation (+162)
    
    ## calculate ppm deviation
    ## mat_1 contains lower values than mat, 
    ## it contains the mz values for M - ppm
    mat_1 <- mat / abs(ppm / 10 ^ 6 + 1) 
    ## mat_2 contains higher values than mat
    ## it contains the mz values for M + ppm
    mat_2 <- mat / abs(ppm / 10 ^ 6 - 1) 
    
    ## calculate difference between rownames and colnames
    ## (hypothetically possible differences between features)
    .mat_1 <- mat
    tmp <- t(mat_1) - mat_2
    .mat_1[upper.tri(tmp, diag = TRUE)] <- tmp[upper.tri(tmp, diag = TRUE)]
    tmp <- -1 * (mat_1 - t(mat_2))
    .mat_1[lt] <- tmp[lt]
   
    
    .mat_2 <- mat
    tmp <- t(mat_2) - mat_1
    .mat_2[upper.tri(tmp, diag = TRUE)] <- tmp[upper.tri(tmp, diag = TRUE)]
    tmp <- -1 * (mat_2 - t(mat_1))
    .mat_2[lt] <- tmp[lt]
    
    ## write the A and B values back to mat_1 and mat_2
    mat_1 <- .mat_1 ## max in lower.tri, min in upper.tri
    mat_2 <- .mat_2 ## min in lower.tri, max in upper.tri

    
    if (!directed) {
 
        ## when the m/z differences between two features is small it is 
        ## possible that one of the differences is negative while the other is 
        ## positive, if we have a transformation with low m/z difference 
        ## (e.g. 0) the link will be likely missed, e.g. when we take the 
        ## absolute values of the differences, we mess up with the negatively 
        ## signed difference, e.g. -1.12 will be converted to 1.12 -> we will 
        ## then check in a interval of 1.12 and the upper bound. 
        ## the following block will take this into account (e.g. by setting 
        ## -1.12 to 0) when there are elements of different sign at the same 
        ## corresponding cell, this will keep the links between the features
        sign_mat <- sign(mat_1) * sign(mat_2)
        mat_1 <- ifelse(sign_mat == -1 & mat_1 < 0, 0, mat_1)
        mat_2 <- ifelse(sign_mat == -1 & mat_2 < 0, 0, mat_2)
        
        ## for the undirected case do not take into account the sign of the 
        ## mass difference
        mat_1_abs <- abs(mat_1)
        mat_2_abs <- abs(mat_2)
        mat_1 <- ifelse(mat_1_abs <= mat_2_abs, mat_2_abs, mat_1_abs) ## max
        mat_2 <- ifelse(mat_1_abs > mat_2_abs, mat_2_abs, mat_1_abs) ## min
    }

    ## create matrices to store result, 
    ## binary contains binary information (0/1) if a connection between two
    ## features is present (this matrix is hard-coded and required),
    ## the other matrices are created depending on the var parameter, e.g.
    ## if there is "group" and "mass" in var, the matrices with names "group"
    ## and "mass" will be added to the list l
    var_binary <- c("binary", var)
    l <- lapply(var_binary, 
        function(x) {
            if (x == "binary") fill <- 0 else fill <- ""
            matrix(fill, nrow = length(mass), ncol = length(mass))
        })
    names(l) <- var_binary

    ## iterate through each column and check if the "mass" is in the interval
    ## defined by the m/z value and ppm
    for (i in seq_along(transformation[, "mass"])) {
        
        transf_i <- transformation[i, ]
        
        ## get intersect from the two (indices where "mass" is in the interval)
        ind_hit <- which(
            (mat_1 >= transf_i[["mass"]] & mat_2 <= transf_i[["mass"]]) |
                (mat_1 <= transf_i[["mass"]] & mat_2 >= transf_i[["mass"]]))
        
        ## write to these indices 1 in the case of the binary matrix
        l[["binary"]][ind_hit] <- 1
        
        ## write to these indices the values stores in transf_i for the 
        ## respective column (paste the value to the entry if there is 
        ## already a value in the cell)
        l_var <- lapply(var, function(var_i) {
            l[[var_i]][ind_hit] <- ifelse(l[[var_i]][ind_hit] != "",
                yes = paste(l[[var_i]][ind_hit], transf_i[[var_i]], sep = "/"),
                no = as.character(transf_i[[var_i]]))
            return(l[[var_i]])
        })
        l[var] <- l_var
    }

    ## add the rownames of x to the list entries (matrices)
    l <- lapply(l, function(l_i) {
        rownames(l_i) <- colnames(l_i) <- rownames(x)
        return(l_i)
    })

    ## create the AdjacencyMatrix object
    rD <- DataFrame(names = rownames(l[["binary"]]), 
                    row.names = rownames(l["binary"]))

    am <- AdjacencyMatrix(l, rowData = rD, type = "structural", 
        directed = directed, thresholded = FALSE)

    return(am)
}

#' @name rtCorrection
#'
#' @aliases rtCorrection
#'
#' @title Correct connections in the structural adjacency matrices by
#' retention time
#'
#' @description
#' The function `rtCorrection` corrects the adjacency matrix
#' infered from structural data based on shifts in the retention time. For
#' known chemical modifications (e.g. addition of glycosyl groups) molecules
#' with the moiety should elute at a different time (in the case of glycosyl
#' groups the metabolite should elute earlier in a reverse-phase
#' liquid chromatography system). If the connection for the metabolite does not
#' fit the expected behaviour, the connection will be removed (otherwise
#' sustained).
#'
#' @param 
#' am `AdjacencyMatrix` object returned by the function `structural`. 
#' The object contains the assays `"binary"` and additional assays with 
#' `character` matrices (only the `"binary"` assay is required). 
#' The assay `"binary"` stores the `numeric` matrix with edges inferred by mass
#' differences.
#'
#' @param 
#' x `matrix`, where columns are the samples and the rows are features
#' (metabolites), cell entries are intensity values, `x` contains the
#' columns "`mz`" and `"rt"` that has the m/z and rt information 
#' (numerical values) for the correction of retention time shifts between 
#' features that have a putative connection assigned based on m/z value 
#' difference. 
#'
#' @param 
#' transformation `data.frame`, containing the columns `var`,
#' and `"rt"` that will be used for correction of transformation of
#' (functional) groups based on retention time shifts derived from `x`
#' 
#' @param 
#' var `character(1)`, the key that is used for matching 
#' between the column `var` in `transformation` and the assay `var` in `am` 
#'
#' @details
#' `rtCorrection` is used to correct the (unweighted) adjacency matrices
#' returned by `structural` when information is available
#' about the retention time and shifts when certain transformation occur
#' (it is meant to filter out connections that were created by
#' m/z differences that have by chance the same m/z difference but
#' different/unexpected retention time behaviour).
#'
#' `rtCorrection` accesses the assay `transformation` of 
#' `am` and matches the elements in the `var` column
#' against the character matrix. In case of matches, `rtCorrection`
#' accesses the `"mz"` and `"rt"` columns of `x` and calculates the retention
#' time difference between the features. `rtCorrection` then checks
#' if the observed retention time difference matches the expected behaviour
#' (indicated by `"+"` for a higher retention time of the feature with
#' the putative group, `"-"` for a lower retention time of the feature
#' with the putative group or `"?"` when there is no information
#' available or features with that group should not be checked). 
#' 
#' In case several transformation were assigned to a feature/feature pair,
#' the connection will be removed if there is an inconsistency with any 
#' of the given transformations.
#'
#' @return
#' `AdjacencyMatrix` containing the slots `binary` and additional `character` 
#' adjacency matrices.
#' The slot `directed` is inherited from `am`.
#' 
#' The assay `binary` stores the
#' `numeric` `matrix` with edges inferred mass differences corrected by
#' retention time shifts.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' rownames(x_test) <- paste(round(x_test[, "mz"], 2),
#'     round(x_test[, "rt"]), sep = "_")
#' transformation <- rbind(
#'     c("Hydroxylation (-H)", "O", 15.9949146221, "+"),
#'     c("Malonyl group (-H2O)", "C3H2O3", 86.0003939305, "+"),
#'     c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315", "-"),
#'     c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851", "-"),
#'     c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945", "-"))
#' transformation <- data.frame(group = transformation[, 1],
#'          formula = transformation[, 2],
#'          mass = as.numeric(transformation[, 3]),
#'          rt = transformation[, 4])
#' am_struct <- structural(x = x_test, transformation = transformation, 
#'     var = c("group", "mass"), ppm = 10, directed = FALSE)
#' am_struct_rt <- rtCorrection(am = am_struct, x = x_test,
#'      transformation = transformation)
#'
#' ## visualize the adjacency matrices
#' library(igraph)
#' g <- graph_from_adjacency_matrix(assay(am_struct, "binary"),
#'     mode = "undirected")
#' g_rt <- graph_from_adjacency_matrix(assay(am_struct_rt, "binary"),
#'     mode = "undirected")
#'
#' plot(g, edge.width = 2, edge.arrow.size = 0.5, vertex.label.cex = 0.7)
#' plot(g_rt, edge.width = 2, edge.arrow.size = 0.5, vertex.label.cex = 0.7)
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @export
rtCorrection <- function(am, x, transformation, var = "group") {
    
    ## check for integrity of am
    if (!is(am, "AdjacencyMatrix")) {
        stop("'am_structural' is not an 'AdjacencyMatrix' object")
    }
    
    if (!validObject(am)) {
        stop("'am' must be a valid 'AdjacencyMatrix' object")
    }
    
    if (am@thresholded) {
        stop("'am' has been already thresholded")
    }
    
    if (am@type != "structural")
        stop("'am' has to be of type 'structural'")
    
    ## check for integrity of x
    if (!"rt" %in% colnames(x)) stop("'x' does not contain the column 'rt'")
    if (!"mz" %in% colnames(x)) stop("'x' does not contain the column 'mz'")
    if (!all(rownames(x) == rownames(am))) 
        stop("'rownames(x)' do not match 'rownames(am)'/'colnames(am)'")
    
    ## check for integrity of transformation and var
    if (!is.character(var) | length(var) != 1)
        stop("'var' has to be a character of length 1")
    if (!var %in% colnames(transformation))
        stop(sprintf("'transformation' does not contain the column '%s'", 
                                                                        var))
    if (!"rt" %in% colnames(transformation))
        stop("'transformation' does not contain the column 'rt'")
    if (!all(transformation[, "rt"] %in% c("+", "-", "?")))
        stop(c("'transformation[, 'rt']' does contain other levels than",
               " '+'', '-'' or '?'"))

    ## allocate binary, transformation, and mass_difference to mat_bin,
    ## mat_type, and mat_mass
    .nms <- SummarizedExperiment::assayNames(am)
    l <- lapply(.nms, function(.nms_i) SummarizedExperiment::assay(am, .nms_i))
    names(l) <- .nms

    n <- nrow(am)
    rn_mat_bin <- names(am)
    mat_rt <- matrix(0, nrow = n, ncol = n)
    colnames(mat_rt) <- rownames(mat_rt) <- x[rn_mat_bin, "rt"]
    
    ## create matrix which has rownames per row
    mat_rt <- apply(mat_rt, 1, function(x) as.numeric(rownames(mat_rt)))
    
    ## calculate difference between rownames and colnames
    ## (difference between features)
    mat_rt <- t(mat_rt) - mat_rt
    colnames(mat_rt) <- rownames(mat_rt) <- rn_mat_bin

    mat_mz <- matrix(0, nrow = n, ncol = n)
    colnames(mat_mz) <- rownames(mat_mz) <- x[rn_mat_bin, "mz"]
    
    ## create matrix which has rownames per row
    mat_mz <- apply(mat_mz, 1, function(i) as.numeric(rownames(mat_mz)))
    
    ## calculate difference between rownames and colnames
    ## (difference between features)
    mat_mz <- mat_mz - t(mat_mz)
    colnames(mat_mz) <- rownames(mat_mz) <- rn_mat_bin

    ## get indices of matching items
    ind <- lapply(seq_len(nrow(transformation)), function(i)
        grep(l[[var]], pattern = transformation[i, var], fixed = TRUE))

    ## iterate through transformation rows
    for (i in seq_len(nrow(transformation))) {

        ind_i <- ind[[i]]
        
        ## check if observed rt shift corresponds to expected one and
        ## remove connection if necessary
        l <- lapply(.nms, function(.nms_i) {
            
            ## obtain the list entry at position .nms_i (e.g. the "binary" slot)
            l_i <- l[[.nms_i]]
            
            ## define the value to which the element is set to
            if (mode(l[[.nms_i]]) == "numeric") {
                vals <- 0
            } else { 
                vals <- ""
            }
            
            ## in case there is a shift to larger retention time for the 
            ## addition, check if the mz/rt pairs match with the expected 
            ## behaviour
            if (transformation[i, "rt"] == "+") {
                l_i[ind_i][mat_mz[ind_i] < 0 & mat_rt[ind_i] < 0] <- vals
                l_i[ind_i][mat_mz[ind_i] > 0 & mat_rt[ind_i] > 0] <- vals
            }
            ## in case there is a shift to lower retention time for the 
            ## addition, check if the mz/rt pairs match with the expected 
            ## behaviour
            if (transformation[i, "rt"] == "-") {
                l_i[ind_i][mat_mz[ind_i] < 0 & mat_rt[ind_i] > 0] <- vals
                l_i[ind_i][mat_mz[ind_i] > 0 & mat_rt[ind_i] < 0] <- vals
            }
            
            return(l_i)
        })
        names(l) <- .nms
    }

    ## create the AdjacencyMatrix object
    rD <- DataFrame(names = rn_mat_bin, row.names = rn_mat_bin)
    
    am <- AdjacencyMatrix(l, rowData = rD, type = "structural", 
        directed = am@directed, thresholded = TRUE)
    
    return(am)
}


#' @name addSpectSimil
#'
#' @aliases addSpectSimil
#'
#' @title Adding a spectral similarity matrix to the "structural" 
#' `AdjacencyMatrix`.
#' 
#' @description
#' The function `addSpectSimil` infers adjacency matrix topologies from
#' spectral similarity methods and adds matrices of these networks into the 
#' "structural" `AdjacencyMatrix` object. 
#' The function includes functionality to calculate adjacency matrices based on 
#' spectral similarity methods included in the `Spectra` package.
#' It uses a `Spectra`-object, storing the ms2 information and using a previously 
#' created `AdjacencyMatrix` object from the type "structural" in order to 
#' perform the mapping on the ms1 features. 
#' The function returns an `AdjacencyMatrix` object of adjacency matrices that 
#' are defined by `methods`.
#'
#' @param spectral 
#' `Spectra` object that contains a unique "id" (see `spectraVariables()`), 
#' matching to the row-/colnames of the structural `AdjacencyMatrix` and storing
#' important information of MS2 data (i.e. mz and intensity).
#' 
#' @param am_structural 
#' `AdjacencyMatrix` of type "structural" that was created using matching MS1 
#' data of the same data set. The respective spectral similarity matrices will 
#' be added into am_structural 
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
#' The function `addSpectSimil` includes functionality to calculate adjacency
#' matrices based on different spectral similarity measures provided by the 
#' `RforMassSpectrometry` infrastructure. 
#' `addSpectSimil` calls the different functions (`FUN` parameter in 
#' `Spectra::compareSpectra()`) as specified by `methods`. The default is the 
#' normalized dotproduct `"ndotproduct"`. Moreover, different parameters (e.g. 
#' `MAPFUN`, `tolerance`, `ppm`, `type`, ...) of the function 
#' `Spectra::compareSpectra()` are forwarded by `...`. 
#' `addSpectSimil` will create adjacency matrices using the specified methods 
#' and will return the "structural" `AdjacencyMatrix` containing the added
#' weighted adjacency matrices in the `assays` slot.
#'
#' It is very important that features IDs in the MS1 data (i.e. row/colnames of 
#' `am_structural`) are matching to the IDs of the respective MS2 data (i.e. 
#' `"id"` in the `spectraVariables()`). Also, the `Spectra` object `spectra` is 
#' required to have unique `"id"`s, meaning on representative spectrum per 
#' feature. If features store multiple spectra, spectra consolidation has to be 
#' performed first (e.g. using the function `Spectra::combineSpectra`).
#' 
#'
#' @return `AdjacencyMatrix` of type "structural" containing the respective 
#' adjacency matrices in the `assay` slot as specified by `methods`
#'
#' @author Liesa Salzer
#'
#' @examples
#' 
#' data("x_test", package = "MetNet")
#' transformation <- rbind(
#'     c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"),
#'     c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851"),
#'     c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945"))
#' transformation <- data.frame(group = transformation[, 1],
#'                                 formula = transformation[, 2],
#'                                 mass = as.numeric(transformation[, 3]))
#' am_struct <- structural(x_test, transformation, var = c("group", "mass"),
#'     ppm = 10, directed = TRUE)
#'
#' data("ms2_test", package = "MetNet")
#' 
#' am_struct-spect <- addSpectSimil(spectra = sps_sub, 
#'                            am_structural = am_struct, 
#'                            methods = "ndotproduct", 
#'                            tolerance = 0.05)
#'
#' @export
#' @importFrom MsCoreUtils ndotproduct
#' @importFrom Spectra spectraVariables

addSpectSimil <- function(spectra, am_structural, 
                          methods = c("ndotproduct"), ...) {
  
  ## sanity checks
  if (!is(spectra, "Spectra")) 
    stop("'x' is not a 'Spectra' object")
  
  if (!is(am_structural, "AdjacencyMatrix")) 
    stop("'am_structural' is not an 'AdjacencyMatrix' object")
  
  if (!validObject(am_structural)) 
    stop("'am_structural' must be a valid 'AdjacencyMatrix' object")
  
  if(!"id" %in% spectraVariables(spectra)) 
    stop("Spectra does not contain id column!")
 
  if(!length(unique(spectra$id)) == length(spectra)) 
    stop("Spectra shall contain only unique entries!")
    
  l <- list()
  
  for(method in methods) {
    ## create empty slots of adjacency matrix based on feature names of MS1 
    ## data and where we can fill then similarity values of MS2 data
    am = matrix(NA, nrow = nrow(am_structural), ncol = nrow(am_structural))
    rownames(am) = rownames(am_structural)
    colnames(am) = colnames(am_structural)
    
    
    ## create spectral similarity matrix
    adj_spec <- Spectra::compareSpectra(spectra,
                                        FUN = get(method),
                                        ...)
    
    colnames(adj_spec) <- spectra$id
    rownames(adj_spec) <- spectra$id
    
    am[rownames(adj_spec), colnames(adj_spec)] <- adj_spec
    
    ## assign the spectral similarity matrix to a new slot in structural
    assay(am_structural, method) <- am
    
  }
  
  
 return(am_structural)
  
}
