#' @name structural
#'
#' @aliases structural
#'
#' @title Create adjacency matrix based on m/z (molecular weight) difference
#'
#' @description
#' The function `structural` infers an unweighted
#' adjacency matrix using differences in m/z values that are matched against a
#' `data.frame` of calculated theoretical differences of
#' loss/addition of functional groups. `structural` returns
#' an `AdjacencyMatrix` object containing
#' the unweighted `numeric` `matrix` (assay `binary`), together with a 
#' `character` `matrix` with the type of loss/addition (assay `transformation`), 
#' and the `character` `matrix` with the mass differences (assay 
#' `mass_difference`). 
#' 
#' @param
#' x `matrix`, where columns are the samples and the rows are features
#' (metabolites), cell entries are intensity values. `x` contains the
#' column `"mz"` that has the m/z information (numerical values) for the
#' calculation of mass differences between features
#'
#' @param
#' transformation `data.frame`, containing the columns `"group"`,
#' and `"mass"` that will be used for detection of transformation of
#' (functional) groups
#'
#' @param
#' ppm `numeric`, mass accuracy of m/z features in parts per million (ppm)
#' 
#' @param
#' directed `logical`, if `TRUE` absolute values of m/z differences will be
#' taken to query against `transformation`  (irrespective the sign of `mass`)
#' and undirected adjacency matrices will be returned as the respective
#' assays. if `FALSE` directed
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
#' `binary`, `transformation`, and `mass_difference`. The `type` slot is
#' set to `structural`. The `directed` slot is set accordingly to the 
#' `directed` argument of the function `structural`.
#' The `thresholded` slot is set to `FALSE`
#'
#' @return
#' `AdjacencyMatrix` object. The object will store the adjacency matrices
#' in the assay slots. The first entry stores the `numeric`
#' `matrix` with binary edges inferred from mass differences. The second entry
#' stores the `character` `matrix` with the type (corresponding to the
#' `"group"` column in `transformation`) is stored. The third entry
#' stores the `character` `matrix` with the `mass_difference` information
#' (corresponding to the `"mass"` column in `transformation`).
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com} and 
#' Liesa Salzer \email{liesa.salzer@@helmholtz-muenchen.de}
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
#' am_struct <- structural(x_test, transformation, ppm = 5, directed = TRUE)
#'
#' @export
structural <- function(x, transformation, ppm = 5, directed = FALSE) {

    if (!is.data.frame(transformation))
        stop("transformation is not a data.frame")
    if (!"group" %in% colnames(transformation))
        stop("transformation does not contain the column group")
    if (!"mass" %in% colnames(transformation))
        stop("transformation does not contain the column mass")
    if (!"mz" %in% colnames(x)) stop("x does not contain the column mz")

    if (!is.numeric(ppm)) stop("ppm is not numeric")

    mass <- x[, "mz"]
    mat <- matrix(0, nrow = length(mass), ncol = length(mass))
    rownames(mat) <- colnames(mat) <- mass

    ## create matrix which has rowmames per row
    mat <- apply(mat, 1, function(x) as.numeric(mass))

    ## calculate ppm deviation
    mat_1 <- mat / abs(ppm / 10 ^ 6 + 1)
    mat_2 <- mat / abs(ppm / 10 ^ 6 - 1)

    ## calculate difference between rownames and colnames
    ## (difference between features)
    mat_1 <- mat - t(mat_1) ## max
    mat_2 <- mat - t(mat_2) ## min

    if (!directed) {
        mat_1_abs <- abs(mat_1)
        mat_2_abs <- abs(mat_2)
        mat_1 <- ifelse(mat_1_abs <= mat_2_abs, mat_2_abs, mat_1_abs) ## max
        mat_2 <- ifelse(mat_1_abs > mat_2_abs, mat_2_abs, mat_1_abs) ## min
    }

    ## create two matrices to store result
    mat_bin <- matrix(0, nrow = length(mass), ncol = length(mass))
    mat_type <- matrix("", nrow = length(mass), ncol = length(mass))
    mat_mass <- matrix("", nrow = length(mass), ncol = length(mass))

    ## iterate through each column and check if the "mass" is in the interval
    ## defined by the m/z value and ppm
    for (i in seq_along(transformation[, "mass"])) {
        
        transf_i <- transformation[i, ]
        ind_mat_1 <- which(mat_1 >= transf_i[["mass"]])
        ind_mat_2 <- which(mat_2 <= transf_i[["mass"]])

        ## get intersect from the two (indices where "mass" is in the interval)
        ind_hit <- intersect(ind_mat_1, ind_mat_2)

        ## write to these indices 1, the "group", and the mass 
        ## (paste the value to group and mass if there is already a value in the
        ## cell)
        mat_bin[ind_hit] <- 1
        mat_type[ind_hit] <- ifelse(mat_type[ind_hit] != "",
            yes = paste(mat_type[ind_hit], transf_i[["group"]], sep = "/"),
            no = as.character(transf_i[["group"]]))
        mat_mass[ind_hit] <- ifelse(mat_mass[ind_hit] != "",
            yes = paste(mat_mass[ind_hit], transf_i[["mass"]], sep = "/"),
            no = as.character(transf_i[["mass"]]))
    }

    rownames(mat_bin) <- colnames(mat_bin) <- rownames(x)
    rownames(mat_type) <- colnames(mat_type) <- rownames(x)
    rownames(mat_mass) <- colnames(mat_mass) <- rownames(x)
    
    ## create the AdjacencyMatrix object
    l <- list(binary = mat_bin, transformation = mat_type, 
        mass_difference = mat_mass)
    rD <- DataFrame(names = rownames(mat_bin), row.names = rownames(mat_bin))
    
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
#' with the moiety should elue at a different time (in the case of glycosyl
#' groups the metabolite should elute earlier in a reverse-phase
#' liquid chromatography system). If the connection for the metabolite does not
#' fit the expected behaviour, the connection will be removed (otherwise
#' sustained).
#'
#' @param
#' am `AdjacencyMatrix` object returned by the function `structural`. 
#' The object contains the assays `"binary"`, `"transformation"`, and 
#' `"mass_difference"`. 
#' The assay `"binary"` stores the `numeric` matrix with edges inferred by mass
#' differences. 
#' The assay `"transformation"` stores the `character` matrix with the type 
#' (corresponding to the `"group"` column in `transformation`).
#' The assay `"mass_difference"` stores the `character` matrix with the type 
#' (corresponding to the `"mass"` column in `transformation`).
#'
#' @param
#' x `matrix`, where columns are the samples and the rows are features
#' (metabolites), cell entries are intensity values, `x` contains the
#' column `"rt"` that has the rt information (numerical values) for the
#' correction of retention time shifts between features that
#' have a putative connection assigned based on m/z value difference
#'
#' @param
#' transformation `data.frame`, containing the columns `"group"`,
#' and `"rt"` that will be used for correction of transformation of
#' (functional) groups based on retention time shifts derived from `x`
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
#' `am` and matches the elements in the `"group"` column
#' against the character matrix. In case of matches, `rtCorrection`
#' accesses the `"rt"` column of `x` and calculates the retention
#' time difference between the features. `rtCorrection` then checks
#' if the observed retention time difference matches the expected behaviour
#' (indicated by `"+"` for a higher retention time of the feature with
#' the putative group, `"-"` for a lower retention time of the feature
#' with the putative group or `"?"` when there is no information
#' available or features with that group should not be checked). In case
#' several transformation were assigned to a feature/feature pair connections
#' will always be removed if there is an inconsistency with any of the given
#' transformation.
#'
#' @return
#' `AdjacencyMatrix` containing the slots `binary`, `transformation`, and 
#' `mass_difference`.
#' The slot `directed` is inherited from `am`
#' 
#' . The first entry stores the
#' `numeric` `matrix` with edges inferred mass differences corrected by
#' retention time shifts. The second entry stores the `character` matrix with
#' the type (corresponding to the `"group`" column
#' in `transformation``) is stored.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' transformation <- rbind(
#'     c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315", "-"),
#'     c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851", "-"),
#'     c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945", "-"))
#' transformation <- data.frame(group = transformation[,1 ],
#'          formula = transformation[, 2],
#'          mass = as.numeric(transformation[, 3]),
#'          rt = transformation[, 4])
#' am_struct <- structural(x = x_test, transformation = transformation, ppm = 5)
#' am_struct_rt <- rtCorrection(am = am_struct, x = x_test, 
#      transformation = transformation)
#'
#' @export
rtCorrection <- function(am, x, transformation) {
    
    if (!is(am, "AdjacencyMatrix")) {
        stop("'am_structural' is not an 'AdjacencyMatrix' object")
    }
    
    if (!validObject(am)) {
        stop("'am' must be a valid 'AdjacencyMatrix' object")
    }
    
    if (thresholded(am)) {
        stop("'am' has been already thresholded")
    }

    ## allocate binary, transformation, and mass_difference to mat_bin,
    ## mat_type, and mat_mass
    mat_bin <- assay(am, "binary")
    mat_type <- assay(am, "transformation")
    mat_mass <- assay(am, "mass_difference")

    if (!"rt" %in% colnames(x)) stop("x does not contain the column rt")

    if (!"group" %in% colnames(transformation))
        stop("transformation does not contain the column group")

    if (!"rt" %in% colnames(transformation))
        stop("transformation does not contain the column rt")

    if (!"mass" %in% colnames(transformation))
        stop("transformation does not contain the column mz")

    if (!all(transformation[, "rt"] %in% c("+", "-", "?")))
        stop(c("transformation[, 'rt'] does contain other levels than",
                " '+'', '-'' or '?'"))
    
    n <- nrow(mat_bin)
    rn_mat_bin <- rownames(mat_bin)
    mat_rt <- matrix(0, nrow = n, ncol = n)
    colnames(mat_rt) <- rownames(mat_rt) <- x[rn_mat_bin, "rt"]
    
    ## create matrix which has rownmames per row
    mat_rt <- apply(mat_rt, 1, function(x) as.numeric(rownames(mat_rt)))
    
    ## calculate difference between rownames and colnames
    ## (difference between features)
    mat_rt <- t(mat_rt) - mat_rt
    colnames(mat_rt) <- rownames(mat_rt) <- rn_mat_bin

    mat_mz <- matrix(0, nrow = n, ncol = n)
    colnames(mat_mz) <- rownames(mat_mz) <- x[rn_mat_bin, "mz"]
    
    ## create matrix which has rowmames per row
    mat_mz <- apply(mat_mz, 1, function(x) as.numeric(rownames(mat_mz)))
    
    ## calculate difference between rownames and colnames
    ## (difference between features)
    mat_mz <- mat_mz - t(mat_mz)
    colnames(mat_mz) <- rownames(mat_mz) <- rn_mat_bin

    ## get indices of matching items
    ind <- lapply(seq_len(nrow(transformation)), function(x)
        grep(mat_type, pattern = transformation[x, 1], fixed = TRUE))

    ## iterate through transformation rows
    for (j in seq_len(nrow(transformation))) {

        ## check if observed rt shift corresponds to expected one and
        ## remove connection if necessary
        if (transformation[j, "rt"] == "+") {
            mat_bin[ind[[j]]][mat_mz[ind[[j]]] < 0 & mat_rt[ind[[j]]] > 0] <- 0
            mat_type[ind[[j]]][mat_mz[ind[[j]]] < 0 & mat_rt[ind[[j]]] > 0] <- ""
            mat_mass[ind[[j]]][mat_mz[ind[[j]]] < 0 & mat_rt[ind[[j]]] > 0] <- ""
            mat_bin[ind[[j]]][mat_mz[ind[[j]]] > 0 & mat_rt[ind[[j]]] < 0] <- 0
            mat_type[ind[[j]]][mat_mz[ind[[j]]] > 0 & mat_rt[ind[[j]]] < 0] <- ""
            mat_mass[ind[[j]]][mat_mz[ind[[j]]] > 0 & mat_rt[ind[[j]]] < 0] <- ""
        }

        if (transformation[j, "rt"] == "-") {
            mat_bin[ind[[j]]][mat_mz[ind[[j]]] < 0 & mat_rt[ind[[j]]] < 0] <- 0
            mat_type[ind[[j]]][mat_mz[ind[[j]]] < 0 & mat_rt[ind[[j]]] < 0] <- ""
            mat_mass[ind[[j]]][mat_mz[ind[[j]]] < 0 & mat_rt[ind[[j]]] < 0] <- ""
            mat_bin[ind[[j]]][mat_mz[ind[[j]]] > 0 & mat_rt[ind[[j]]] > 0] <- 0
            mat_type[ind[[j]]][mat_mz[ind[[j]]] > 0 & mat_rt[ind[[j]]] > 0] <- ""
            mat_mass[ind[[j]]][mat_mz[ind[[j]]] > 0 & mat_rt[ind[[j]]] > 0] <- ""
        }
    }

    ## create the AdjacencyMatrix object
    l <- list(binary = mat_bin, transformation = mat_type, 
              mass_difference = mat_mass)
    rD <- DataFrame(names = rn_mat_bin, row.names = rn_mat_bin)
    
    am <- AdjacencyMatrix(l, rowData = rD, type = "structural", 
        directed = directed(am), thresholded = TRUE)
    
    return(am)
}
