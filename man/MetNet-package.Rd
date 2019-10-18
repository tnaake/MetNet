\name{MetNet-package}
\alias{MetNet-package}
\alias{MetNet}
\docType{package}
\title{
\packageTitle{MetNet}
}
\description{
Inferring metabolic networks from untargeted high-resolution mass 
spectrometry data.
}
\details{
The package infers network topologies from quantitative data 
(intensity values) and structural data (m/z values of mass features). 
\code{MetNet} combines these two data sources to a consensus matrix.
}
\author{
Author: \packageAuthor{MetNet}
Maintainer: \packageMaintainer{MetNet}
}
\references{
Breitling, R. et al. Ab initio prediction of metabolic networks using 
Fourier transform mass spectrometry data. 2006. Metabolomics 2: 155--164. 
10.1007/s11306-006-0029-z
}
\keyword{mass spectrometry, metabolomics }

\examples{
data("x_test", package = "MetNet")
x_test <- as.matrix(x_test)
functional_groups <- rbind(
    c("Hydroxylation (-H)", "O", "15.9949146221"),
    c("Malonyl group (-H2O)", "C3H2O3", "86.0003939305"),
    c("C6H10O6", "C6H10O6", "178.0477380536"),
    c("D-ribose (-H2O) (ribosylation)", "C5H8O4", "132.0422587452"),
    c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851"),
    c("Glucuronic acid (-H2O)", "C6H8O6", "176.0320879894"),
    c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"),
    c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945"))
functional_groups <- data.frame(group = functional_groups[,1],
    formula = functional_groups[,2],
    mass = as.numeric(functional_groups[,3]))
struct_adj <- structural(x_test, functional_groups, ppm = 5)

stat_adj_l <- statistical(x_test, 
    model = c("pearson", "spearman","bayes"))
args_top1 <- list(n = 10)
stat_adj <- threshold(stat_adj_l, type = "top2", args = args_top1)
cons_adj <- combine(struct_adj, stat_adj)
}
