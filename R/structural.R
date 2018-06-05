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
    in_range <- which(lower_n <= functional_groups & 
                          higher_n >= functional_groups)
    ## only remove one feature when several are in range
    in_range <- in_range[1]
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
#' @examples 
#' data("x_test", package = "MetNet")
#' x <- as.matrix(x)
#' functional_groups <- rbind(
#'     c("Acetylation (–H)", "C2H3O2", "59.0133043405"),
#'     c("Acetylation (–H2O)", "C2H2O",  "42.0105646863"),
#'     c("Phenylalanine", "C9H9NO",  "147.0684139162"),
#'     c("C2H2", "C2H2", "26.0156500642"),
#'     c("Carboxylation", "CO2", "43.9898292442"),
#'     c("CHO2", "CHO2", "44.9976542763"),
#'     c("Condensation/dehydration", "H2O", "18.0105646863"),
#'     c("Hydrogenation/dehydrogenation", "H2", "2.0156500642"),
#'     c("Hydroxylation (–H)", "O", "15.9949146221"),
#'     c("Secondary amine", "NH", "15.0108990373"),
#'     c("Malonyl group (–H2O)", "C3H2O3", "86.0003939305"),
#'     c("Tertiary amine", "N", "14.0030740052"),
#'     c("C6H10O6", "C6H10O6", "178.0477380536"),
#'     c("D-ribose (–H2O) (ribosylation)", "C5H8O4", "132.0422587452"),
#'     c("Disaccharide (–H2O)", "C12H20O11", "340.1005614851"),
#'     c("Glucose-N-phosphate (–H2O)", "C6H11O8P", "242.0191538399"),
#'     c("Glucuronic acid (–H2O)", "C6H8O6", "176.0320879894"),
#'     c("Monosaccharide (–H2O)", "C6H10O5", "162.0528234315"),
#'     c("Trisaccharide (–H2O)", "C18H30O15", "486.1584702945"))
#' struct_net <- create_structural_network(x, functional_groups, ppm = 5)
#' @export
create_structural_network <- function(x, functional_groups, ppm = 5) {
    
    if (!is.data.frame(functional_groups)) 
        stop("functional_groups is not a data.frame")
    if (!all(c("group", "formula", "mass") %in% colnames(functional_groups))) {
        stop("functional_groups does not contain the columns group, formula and mass")
    }
    
    mass <- x[, "mz"]
    
    ## calculate according to ppm = (mass_measured - mass_theoretical) / mass_theoretical * 10^6
    mass_1 <- mass / abs(ppm / 10 ^ 6  + 1 ) 
    mass_2 <- mass / abs(ppm / 10 ^ 6 - 1)
    mat <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[1])
    mat_type <- matrix(NA, ncol = nrow(mat), nrow(mat))
    
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
    rownames(mat) <- colnames(mat) <- rownames(x)
    rownames(mat_type) <- colnames(mat_type) <- rownames(x)
    return(list(mat, mat_type))
}

