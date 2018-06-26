#' @name inRangeWhich
#' @aliases inRangeWhich
#' @title inRangeWhich
#' @description The function \code{inRangeWhich} checks if a m/z feature
#' (given as a range defined by \code{m_1} and \code{m_2}) is present in a 
#' list of m/z features given by the \code{transformation} argument. 
#' @usage inRangeWhich(m_1, m_2, transformation)
#' @param m_1 numeric, value of the first m/z range
#' @param m_2 numeric, value of the second m/z range
#' @param transformation numeric, m/z values of the functional groups
#' @details The function \code{inRangeWhich} is a helper function for 
#' \code{createStructuralAdjacency}. \code{inRangeWhich} returns the index of
#' the element in \code{transformation} that is in the range defined by 
#' \code{m_1} and \code{m_2}. 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @return numeric
#' @examples 
#' m_1 <- 180.155
#' m_2 <- 180.157
#' transformation <- c(164.16, 180.156)
#' inRangeWhich(m_1, m_2, transformation)
#' @export
inRangeWhich <- function(m_1, m_2, transformation) {
    if (!is.numeric(m_1)) stop("m_1 is not numeric")
    if (!is.numeric(m_2)) stop("m_2 is not numeric")
    if (!is.numeric(transformation)) stop("transformation is not numeric")
    lower_n <- if(m_1 < m_2) m_1 else m_2
    higher_n <- if(m_1 >= m_2) m_1 else m_2
    in_range <- which(lower_n <= transformation & higher_n >= transformation)
    return(in_range)
}

#' @name createStructuralAdjacency
#' @aliases createStructuralAdjacency
#' @title Create adjacency matrix based on m/z (molecular weight) difference
#' @description The function \code{createStructuralAdjacency} infers an 
#' adjacency matrix using differences in m/z values that are matched against a 
#' \code{data.frame} of theoretically calculated differences of 
#' loss/addition of functional groups. \code{createStructuralAdjacency} returns 
#' the unweighted adjacency matrix together with a character matrix with the 
#' type of loss/addition as a list at the specific positions.  
#' @usage createStructuralAdjacency(x, transformation, ppm=5)
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values, \code{x} contains the 
#' column \code{'mz'} that has the m/z information (numerical values) for the 
#' calculation of mass differences between features 
#' @param transformation data.frame, containing the columns \code{"group"}, 
#' and \code{'mass'} that will be used for detection of transformation of 
#' (functional) groups
#' @param ppm numeric, mass accuracy of m/z features in parts per million (ppm)
#' @details \code{createStructuralAdjacency} accesses the column \code{'mz'} of 
#' \code{x} to infer structural topologies based on the functional groups 
#' supplied by \code{transformation}. To account for the mass accuracy of 
#' the dataset \code{x}, the user can specify the accuracy of m/z features 
#' in parts per million (ppm) by the \code{ppm} argument. The m/z values in the 
#' \code{'mz'} column of \code{x} will be converted to m/z ranges according to 
#' the \code{ppm} argument (default \code{ppm=5}). 
#' @return list containing two matrices, in the first list entry the 
#' matrix with edges inferred mass differences is stored, in the second list 
#' entry the matrix with the type (corresponding to the \code{"group"} column
#' in \code{transformation}) is stored
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package="MetNet")
#' transformation <- rbind(
#'     c("Hydroxylation (-H)", "O", "15.9949146221"),
#'     c("Malonyl group (-H2O)", "C3H2O3", "86.0003939305"),
#'     c("C6H10O6", "C6H10O6", "178.0477380536"),
#'     c("D-ribose (-H2O) (ribosylation)", "C5H8O4", "132.0422587452"),
#'     c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851"),
#'     c("Glucuronic acid (-H2O)", "C6H8O6", "176.0320879894"),
#'     c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"),
#'     c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945"))
#' transformation <- data.frame(group=transformation[,1],
#'                                 formula=transformation[,2],
#'                                 mass=as.numeric(transformation[,3]))
#' struct_adj <- createStructuralAdjacency(x_test, transformation, ppm=5)
#' @export
createStructuralAdjacency <- function(x, transformation, ppm=5) {
    
    if (!is.data.frame(transformation)) 
        stop("transformation is not a data.frame")
    if (!"group" %in% colnames(transformation))
        stop("transformation doesn't contain the column group")
    if (!"mass" %in% colnames(transformation))
        stop("transformation doesn't contain the column mass")
    
    if (!is.numeric(ppm)) stop("ppm is not numeric")
    
    mass <- x[, "mz"]
    
    ## calculate according to 
    ## ppm=(mass_measured - mass_theoretical) / mass_theoretical * 10^6
    mass_1 <- mass / abs(ppm / 10 ^ 6  + 1 ) 
    mass_2 <- mass / abs(ppm / 10 ^ 6 - 1)
    mat <- matrix(0, nrow=dim(x)[1], ncol=dim(x)[1])
    mat_type <- matrix("", ncol=nrow(mat), nrow(mat))
    
    mass_t <- transformation[, "mass"]
    ## iterate through columns 
    for (i in seq_len(ncol(mat))) {
        transformation_vec <- lapply(seq_len(length(mass_1)), 
            function(x, z=i) {
                m_1 <- abs(mass_1[z] - mass_2[x])
                m_2 <- abs(mass_2[z] - mass_1[x])
                ## use inRangeWhich to find the indices of the mass 
                ## differences, that are in range with the masses 
                presence_l <- inRangeWhich(m_1, m_2, mass_t)
                presence_l <- transformation[presence_l, "group"]
                presence_l <- as.character(presence_l)
                
                if(length(presence_l) > 0) {
                    ## if there are several possible groups mapped, return the
                    ## names separated by "/"
                    if (length(presence_l) > 1) {
                        return(paste(presence_l, collapse="/"))
                    ## if only one possible groups is mapped, return only this 
                    ## name
                    } else {return(presence_l)}
                } else {return("")}
            }
        )
        transformation_vec <- unlist(transformation_vec)
        mat_type[,i] <- transformation_vec
    }
    ## use mat_type to assign 0 or 1 to the adjacency matrix
    mat <- ifelse(mat_type == "", 0, 1)
    rownames(mat) <- colnames(mat) <- rownames(x)
    rownames(mat_type) <- colnames(mat_type) <- rownames(x)
    return(list(mat, mat_type))
}

#' @name rtCorrection
#' @aliases rtCorrection
#' @title Correct connections in the structural adjacency matrix by 
#' retention time 
#' @description The function \code{rtCorrection} corrects the adjacency matrix
#' infered from structural data based on shifts in the retention time. For 
#' known chemical modifications (e.g. addition of glycosyl groups) molecules 
#' with the moiety should elue at a different time (in the case of glycosyl 
#' groups the metabolite should elute earlier in a reverse-phase 
#' liquid chromatography system). If the connection for the metabolite does not
#' fit the expected behaviour, the connection will be removed (otherwise 
#' sustained). 
#' @usage rtCorrection(struct_adj, x, transformation)
#' @param struct_adj list returned by the function 
#' \code{createStructuralAdjacency}, in the first list entry the 
#' matrix with edges inferred mass differences is stored, in the second list 
#' entry the matrix with the type (corresponding to the \code{'group'} column
#' in \code{transformation}) is stored
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values, \code{x} contains the 
#' column \code{'rt'} that has the rt information (numerical values) for the 
#' correction of retention time shifts between features that 
#' have a putative connection assigned based on m/z value difference
#' @param transformation data.frame, containing the columns \code{"group"}, 
#' and \code{'rt'} that will be used for correction of transformation of 
#' (functional) groups based on retention time shifts derived from 
#' \code{x}
#' @details \code{rtCorrection} is used to correct the adjacency matrix 
#' returned by \code{createStructuralAdjacency} when information is available
#' about the retention time and shifts when certain transformation occur
#' (it is meant to filter out connections that were created by 
#' m/z differences that have by chance the same m/z difference but 
#' different/unexpected retention time behaviour). #' 
#' \code{rtCorrection} accesses the second list element of 
#' \code{struct_adj} and matches the elements in the \code{'group'} column 
#' against the character matrix. In case of matches, \code{rtCorrection} 
#' accesses the \code{'rt'} column of \code{x} and calculates the retention 
#' time difference between the features. \code{rtCorrection} then checks 
#' if the observed retention time difference matches the expected behaviour
#' (indicated by \code{'+'} for a higher retention time of the feature with 
#' the putative group, \code{'-'} for a lower retention time of the feature
#' with the putative group or \code{'?'} when there is no information 
#' available or features with that group should not be checked). In case 
#' several transformation were assigned to a feature/feature pair connections 
#' will always be removed if there is an inconsistency with any of the given 
#' transformation. 
#' @return list containing two matrices, in the first list entry the 
#' matrix with edges inferred mass differences corrected by retention time 
#' shifts  is stored, in the second list entry the matrix with the type 
#' (corresponding to the \code{'group'} column
#' in \code{transformation}) is stored
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("x_test", package="MetNet")
#' transformation <- rbind(
#'     c("Hydroxylation (-H)", "O", "15.9949146221", "-"),
#'     c("Malonyl group (-H2O)", "C3H2O3", "86.0003939305", "?"),
#'     c("C6H10O6", "C6H10O6", "178.0477380536", "-"),
#'     c("D-ribose (-H2O) (ribosylation)", "C5H8O4", "132.0422587452", "-"),
#'     c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851", "-"),
#'     c("Glucuronic acid (-H2O)", "C6H8O6", "176.0320879894", "?"),
#'     c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315", "-"),
#'     c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945", "-"))
#' transformation <- data.frame(group=transformation[,1],
#'                                 formula=transformation[,2],
#'                                 mass=as.numeric(transformation[,3]),
#'                                 rt=transformation[,4])
#' struct_adj <- createStructuralAdjacency(x_test, transformation, ppm=5)
#' struct_adj_rt <- rtCorrection(struct_adj, x_test, transformation)
#' @export
rtCorrection <- function(struct_adj, x, transformation) {
    
    ## check arguments 
    if (!is.list(struct_adj)) stop("struct_adj is not a list")
    if (!length(struct_adj) == 2) stop("struct_adj is not a list of length 2")
    ## allocate struct_adj[[1]] and struct_adj[[2]] to adj and group
    adj <- struct_adj[[1]]
    group <- struct_adj[[2]]
    
    if (any(dim(adj) != dim(group))) 
        stop("dim(struct_adj[[1]] is not equal to dim(struct_adj[[2]]")
    if (!"rt" %in% colnames(x)) stop("x does not contain the column rt")
    if (!"group" %in% colnames(transformation)) 
        stop("transformation does not contain the column group")
    if (!"rt" %in% colnames(transformation)) 
        stop("transformation does not contain the column rt")
    if (!"mass" %in% colnames(transformation)) 
        stop("transformation does not contain the column mz")
    if (!all(levels(transformation[, "rt"]) %in% c("+", "-", "?")))
        stop('in transformation[, "rt"] does contain other levels than 
                "+", "-" or "?"' )
    if (!all(colnames(adj) == rownames(adj)))
        stop("colnames of struct_adj[[1]] is not identical to rownames of 
                struct_adj[[1]]")
    if (!all(colnames(struct_adj[[2]]) == rownames(struct_adj[[2]])))
        stop("colnames of struct_adj[[2]] is not identical to 
                rownames of struct_adj[[2]]")
    if (!all(rownames(x) == rownames(struct_adj[[1]]))) 
        stop("rownames(x) do not fit rownames(struct_adj[[1]])")
    if (!all(rownames(x) == colnames(struct_adj[[1]]))) 
        stop("rownames(x) do not fit colnames(struct_adj[[1]])")
    if (!is.matrix(adj)) stop("struct_adj[[1]] is not a matrix")
    if (!is.matrix(group)) stop("struct_adj[[2]] is not a matrix")
    if (!is.numeric(adj)) stop("struct_adj[[1]] is not a numeric matrix")
    if (!is.character(group)) 
        stop("struct_adj[[2]] is not a character matrix")
    
    ## get indices of matching items
    ind <- lapply(1:dim(transformation)[1], function(x) 
        apply(group, c(1,2), function(y) grep(x=as.vector(y), 
                pattern=as.character(transformation[x,1]), fixed=TRUE)))
    ## get row and col indices where matrix == 1
    ind <- lapply(ind, function(x) which(x == 1, arr.ind=TRUE))
    
    ## iterate through transformation rows
    for (j in 1:nrow(transformation)) { 
        ##  iterate within each entry of transformation rows
        if (!is.null(nrow(ind[[j]]))) {
        for (i in 1:nrow(ind[[j]])) { 
            x_1 <- x[rownames(adj)[ind[[j]][i,1]], ]
            x_2 <- x[rownames(adj)[ind[[j]][i,2]], ]
            
            ## check if the groups match to the predicted retention time 
            ## shift and change the connection accordingly
            if (x_1["mz"] - x_2["mz"] > 0) {
                if (transformation[j, "rt"] == "-") {
                    if (x_2["rt"] - x_1["rt"] < 0) {
                        adj[ind[[j]][i,1], ind[[j]][i, 2]] <- 0
                        group[ind[[j]][i,1], ind[[j]][i, 2]] <- ""
                    }
                }
                if (transformation[j, "rt"] == "+") {
                    if (x_1["rt"] - x_2["rt"] < 0) {
                        adj[ind[[j]][i,1], ind[[j]][i, 2]] <- 0
                        group[ind[[j]][i,1], ind[[j]][i, 2]] <- ""
                    }
                }
            }
            if (x_2["mz"] - x_1["mz"] > 0) {
                if (transformation[j, "rt"] == "-") {
                    if (x_1["rt"] - x_2["rt"] < 0) {
                        adj[ind[[j]][i,1], ind[[j]][i, 2]] <- 0
                        group[ind[[j]][i,1], ind[[j]][i, 2]] <- ""
                    }
                }
                if (transformation[j, "rt"] == "+") {
                    if (x_2["rt"] - x_1["rt"] < 0) {
                        adj[ind[[j]][i,1], ind[[j]][i, 2]] <- 0
                        group[ind[[j]][i,1], ind[[j]][i, 2]] <- ""
                    }
                }
            }
        }
        }
    }
    return(list(adj, group))
}


