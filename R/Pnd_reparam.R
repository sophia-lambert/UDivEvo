#' Probability density of node depths reparametrised.
#'
#' @description Calculating the probability density of node depths of a phylogenetic tree parametrised on \eqn{y*\lambda} and \eqn{r}.
#'
#' @param t Numerical vector. A vector specifying the node depths of a phylogeny in millions of years.
#' @param yl Numeric. The sampling probability times the speciation rate : \eqn{y*\lambda}.
#' @param r Numeric. The net diversification rate : \eqn{\lambda - \mu}. The units are in events/lineages/millions years.
#'
#' @return The probability density of node depths equivalent to \eqn{f\_y(t)} in Lambert 2018 and to \eqn{\lambda*p1} in Stadler 2010 and Lambert et al. 2022 but parametrised on \eqn{y*\lambda} and \eqn{r} see also Stadler 2009.
#'
#' @author Sophia Lambert
#'
#' @references Lambert, A. (2018). The coalescent of a sample from a binary branching process. Theoretical population biology, 122, 30-35.
#' Stadler, T. (2010). Sampling-through-time in birth–death trees. Journal of theoretical biology, 267(3), 396-404.
#' Stadler, T. (2009). On incomplete sampling under birth–death models and connections to the sampling-based coalescent. Journal of theoretical biology, 261(1), 58-66.
#'
#' @export

#@seealso \code{\link{likelihood_bdRho}}, \code{\link{fitMCMC_bdRho}}

Pnd_reparam <- function(t, yl, r){
  yl*r^2*exp(-r*t)/
    (yl+
       (r-yl)*exp(-r*t))^2
}
