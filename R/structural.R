#' @name .in_range_which
#' @title .in_range_which
#' @description The function \code{.in_range_which} checks if a m/z feature
#' (given as a range defined by \code{m_1} and \code{m_2}) is present in a 
#' list of m/z features given by the \code{functional_groups} argument. 
#' @usage .in_range_which(m_1, m_2, functional_groups)
#' @param m_1 numeric, value of the first m/z range
#' @param m_2 numeric, value of the second m/z range
#' @param functional_groups numeric, m/z values of the functional groups
#' @details The function \code{.in_range_which} is a helper function for 
#' \code{create_structural_network}. \code{.in_range_which} returns the index of
#' the element in \code{functional_groups} that is in the range defined by 
#' \code{m_1} and \code{m_2}. 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' m_1 <- 180.155
#' m_2 <- 180.157
#' functional_groups <- c(164.16, 180.156)
#' .in_range_which(m_1, m_2, functional_groups)
.in_range_which <- function(m_1, m_2, functional_groups) {
    if (!is.numeric(m_1)) stop("m_1 is not numeric")
    if (!is.numeric(m_2)) stop("m_2 is not numeric")
    if (!is.numeric(functional_groups)) stop("functional_groups is not numeric")
    lower_n <- if(m_1 < m_2) m_1 else m_2
    higher_n <- if(m_1 >= m_2) m_1 else m_2
    in_range <- which(lower_n <= functional_groups & higher_n >= functional_groups)
    return(in_range)
}

#' @name create_structural_network
#' @title Create structural network
#' @description The function \code{create_structural_network} infers network 
#' topologies using differences in m/z values that are matched against a 
#' data.frame of theoretically calculated differences of loss/addition of 
#' functional groups. \code{create_structural_network} returns the 
#' unweighted network together with a character matrix with the type of 
#' loss/addition as a list.  
#' @usage create_structural_network(x, functional_groups, ppm = 5)
#' @param x matrix, where columns are the samples and the rows are features 
#' (metabolites), cell entries are intensity values, \code{x} contains the 
#' column "mz" that has the m/z information (numerical values) for the 
#' calculation of mass differences between features 
#' @param functional_groups data.frame, containing the columns "group", 
#' "formula" and "mass" that will be used for detection of losses/addition of 
#' (functional) groups
#' @param ppm numeric, mass accuracy of m/z features in parts per million (ppm)
#' @details \code{create_structural_network} accesses the column "mz" of 
#' \code{x} to infer structural topologies based on the functional groups 
#' supplied by \code{functional_groups}. To account for the mass accuracy of 
#' the dataset \code{x}, the user can specify the accuracy of m/z features 
#' in parts per million (ppm) by the \code{ppm} argument. The m/z values in the 
#' "mz" column of \code{x} will be converted to m/z ranges according to the 
#' \code{ppm} argument (default \code{ppm = 5}). 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples create_structural_network(x, functional_groups, ppm = 5)
#' @export
create_structural_network <- function(x, functional_groups, ppm = 5) {
    
    if (!is.data.frame(functional_groups)) stop("functional_groups is not a 
                                                data.frame")
    if (!all(colnames(functional_groups) %in% c("group", "formula", "mass"))) {
        stop("functional_groups does not contain the columns group, formula and
             mass")
    }
    
    mass <- x[, "mz"]
    
    ## calculate according to ppm = (mass_measured - mass_theoretical) / mass_theoretical * 10^6
    mass_1 <- mass / abs(ppm / 10 ^ 6  + 1 ) 
    mass_2 <- mass / abs(ppm / 10 ^ 6 - 1)
    mat <- matrix(0, nrow = dim(peaklist)[1], ncol = dim(peaklist)[1])
    mat_type <- matrix(NA, ncol = ncol(mat), nrow(mat))
    
    mass_fg <- functional_groups[, "mass"]
    ## iterate through columns 
    for (i in 1:dim(mat)[2]) {
        functional_group_vec <- lapply(1:length(mass_1), 
            function(x, z = i) {
                m_1 <- abs(mass_1[z] - mass_2[x])
                m_2 <- abs(mass_2[z] - mass_1[x])
                presence_l <- .in_range_which(m_1, m_2, mass_fg)
                presence_l <- functional_groups[presence_l, "group"]
                presence_l <- as.character(presence_l)
                if(length(presence_l) > 0) {
                    if (length(presence_l) > 1) {
                        return(paste(presence_l, collapse = "/"))
                    } else {return(presence_l)}
                } else {return(NA)}
            }
        )
        functional_group_vec <- unlist(functional_group_vec)
        mat_type[,i] <- functional_group_vec
    }
    mat <- ifelse(is.na(mat_type), 0, 1)
    return(list(mat, mat_type))
}

