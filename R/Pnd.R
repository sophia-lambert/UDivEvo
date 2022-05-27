#' Probability density of node depths.
#'
#' @description Calculating the probability density of node depths of a phylogenetic tree.
#'
#' @param t Numerical vector. A vector specifying the node depths of a phylogeny in millions of years.
#' @param r Numeric. The net diversification rate : \eqn{\lambda - \mu}. The units are in events/lineages/millions years.
#' @param epsi Numeric. The turnover rate : \eqn{\mu / \lambda}. The units are in events/lineages/millions years.
#' @param y Numeric. The sampling probability also called sampling fraction : \eqn{k / N}.
#'
#' @return The probability density of node depths equivalent to \eqn{f\_y(t)} in Lambert 2018 and to \eqn{\lambda * p1} in Stadler 2010 and Lambert et al. 2022.
#'
#' @author Sophia Lambert
#'
#' @references Lambert, A. (2018). The coalescent of a sample from a binary branching process. Theoretical population biology, 122, 30-35. Stadler, T. (2010). Sampling-through-time in birth–death trees. Journal of theoretical biology, 267(3), 396-404.
#'
#' @export

#@seealso \code{\link{likelihood_bdRho}}, \code{\link{fitMCMC_bdRho}}

Pnd <- function(t, r, epsi, y){
  y*r^3*exp(-r*t)/
    (1-epsi)/
    (y*r/(1-epsi)+
       (r-y*r/(1-epsi))*exp(-r*t))^2
}
