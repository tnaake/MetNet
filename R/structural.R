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
#' the unweighted `numeric` `matrix` together with a `character` `matrix` with
#' the type of loss/addition as a list at the specific positions.
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
#' and an undirected adjacency matrix will be returned, if `FALSE` a directed
#' adjacency matrix will be returned with links reported that match the
#' transformations defined in `transformation` (respecting the sign of `mass`)
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
#' @return
#' `list` containing two matrices. The first entry stores the `numeric`
#' `matrix` with edges inferred from mass differences. The second entry
#' stores the `character` `matrix` with the type (corresponding to the
#' `"group"` column in `transformation`) is stored
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
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
#' struct_adj <- structural(x_test, transformation, ppm = 5, directed = TRUE)
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
    mat <- matrix(0, nrow = length(mass), ncol = length(mass))
    mat_type <- matrix("", nrow = length(mass), ncol = length(mass))

    ## iterate through each column and check if the "mass" is in the interval
    ## defined by the m/z value and ppm
    for (i in seq_along(transformation[, "mass"])) {
        
        transformation_i <- transformation[i, ]
        ind_mat_1 <- which(mat_1 >= transformation_i[["mass"]])
        ind_mat_2 <- which(mat_2 <= transformation_i[["mass"]])

        ## get intersect from the two (indices where "mass" is in the interval)
        ind_hit <- intersect(ind_mat_1, ind_mat_2)

        ## write to these indices 1 and the "group"
        mat[ind_hit] <- 1
        mat_type[ind_hit] <- ifelse(nchar(mat_type[ind_hit]) != 0,
            yes = paste(mat_type[ind_hit], transformation_i[["group"]],
                sep = "/"),
            no = as.character(transformation_i[["group"]]))
    }

    rownames(mat) <- colnames(mat) <- rownames(x)
    rownames(mat_type) <- colnames(mat_type) <- rownames(x)

    return(list(mat, mat_type))

}

#' @name rtCorrection
#'
#' @aliases rtCorrection
#'
#' @title Correct connections in the structural adjacency matrix by
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
#' structural `list` returned by the function `structural`. The first entry
#' stores the `numeric` matrix with edges inferred by mass
#' differences. The second entry stores the `character` matrix with the type
#' (corresponding to the `"group"` column in `transformation`).
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
#' `rtCorrection` is used to correct the unweighted adjacency matrix
#' returned by `structural` when information is available
#' about the retention time and shifts when certain transformation occur
#' (it is meant to filter out connections that were created by
#' m/z differences that have by chance the same m/z difference but
#' different/unexpected retention time behaviour).
#'
#' `rtCorrection` accesses the second list element of
#' `structural` and matches the elements in the `"group"` column
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
#' `list` containing two matrices. The first entry stores the
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
#'                                 formula = transformation[,2 ],
#'                                 mass = as.numeric(transformation[,3 ]),
#'                                 rt = transformation[, 4])
#' struct_adj <- structural(x_test, transformation, ppm = 5)
#' struct_adj_rt <- rtCorrection(struct_adj, x_test, transformation)
#'
#' @export
rtCorrection <- function(structural, x, transformation) {

    ## check arguments
    if (!is.list(structural)) stop("structural is not a list")
    if (!length(structural) == 2) stop("structural is not a list of length 2")
    ## allocate structural[[1]] and structural[[2]] to adj and group
    adj <- structural[[1]]
    group <- structural[[2]]

    if (any(dim(adj) != dim(group)))
        stop("dim(structural[[1]] is not equal to dim(structural[[2]])")

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

    if (!all(colnames(adj) == rownames(adj)))
        stop(c("colnames of structural[[1]] are not identical to rownames of",
                " structural[[1]]"))

    if (!all(colnames(group) == rownames(group)))
        stop(c("colnames of structural[[2]] are not identical to",
                " rownames of structural[[2]]"))

    if (!all(rownames(structural[[1]]) %in% rownames(x)))
        stop("rownames(structural[[1]]) do not fit rownames(x)")

    if (!(is.matrix(adj) && is.numeric(adj)))
        stop("structural[[1]] is not a numeric matrix")

    if (!(is.matrix(group) && is.character(group)))
        stop("structural[[2]] is not a character matrix")

    mat_rt <- matrix(0, nrow = nrow(adj), ncol = ncol(adj))
    colnames(mat_rt) <- rownames(mat_rt) <- x[rownames(adj), "rt"]
    ## create matrix which has rownmames per row
    mat_rt <- apply(mat_rt, 1, function(x) as.numeric(rownames(mat_rt)))
    ## calculate difference between rownames and colnames
    ## (difference between features)
    mat_rt <- t(mat_rt) - mat_rt
    colnames(mat_rt) <- rownames(mat_rt) <- colnames(adj)

    mat_mz <- matrix(0, nrow = nrow(adj), ncol = ncol(adj))
    colnames(mat_mz) <- rownames(mat_mz) <- x[rownames(adj), "mz"]
    
    ## create matrix which has rownmames per row
    mat_mz <- apply(mat_mz, 1, function(x) as.numeric(rownames(mat_mz)))
    
    ## calculate difference between rownames and colnames
    ## (difference between features)
    mat_mz <- mat_mz - t(mat_mz)
    colnames(mat_mz) <- rownames(mat_mz) <- colnames(adj)

    ## get indices of matching items
    ind <- lapply(seq_len(dim(transformation)[1]), function(x)
        grep(group, pattern = transformation[x, 1], fixed = TRUE))

    ## iterate through transformation rows
    for (j in seq_len(nrow(transformation))) {

        ## check if observed rt shift corresponds to expected one and
        ## remove connection if necessary
        if (transformation[j, "rt"] == "+") {
            adj[ind[[j]]][mat_mz[ind[[j]]] < 0 & mat_rt[ind[[j]]] > 0] <- 0
            group[ind[[j]]][mat_mz[ind[[j]]] < 0 & mat_rt[ind[[j]]] > 0] <- ""
            adj[ind[[j]]][mat_mz[ind[[j]]] > 0 & mat_rt[ind[[j]]] < 0] <- 0
            group[ind[[j]]][mat_mz[ind[[j]]] > 0 & mat_rt[ind[[j]]] < 0] <- ""
        }

        if (transformation[j, "rt"] == "-") {
            adj[ind[[j]]][mat_mz[ind[[j]]] < 0 & mat_rt[ind[[j]]] < 0] <- 0
            group[ind[[j]]][mat_mz[ind[[j]]] < 0 & mat_rt[ind[[j]]] < 0] <- ""
            adj[ind[[j]]][mat_mz[ind[[j]]] > 0 & mat_rt[ind[[j]]] > 0] <- 0
            group[ind[[j]]][mat_mz[ind[[j]]] > 0 & mat_rt[ind[[j]]] > 0] <- ""
        }
    }

    return(list(adj, group))
}
