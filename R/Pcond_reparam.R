#' Conditioning the birth-death process on the age and the survival of the process reparametrised.
#'
#' @description Calculating the probability of conditioning the birth-death process on the age of the process and its survival parametrised on \eqn{y*\lambda} and \eqn{r}.
#'
#' @param t Numeric. The time of the start of the birth-death process in millions of years.
#' @param yl Numeric. The sampling probability times the speciation rate: \eqn{y*\lambda}.
#' @param r Numeric. The net diversification rate: \eqn{\lambda - \mu}. The units are in events/lineages/millions years.
#'
#' @return The probability of the birth-death process starting with one lineage at time t and surviving at present parametrised on \eqn{y*\lambda} and \eqn{r}.
#'
#' @author Sophia Lambert
#'
#' @export
#'
#' @seealso \code{\link{likelihood_bdRho}}, \code{\link{fitMCMC_bdRho}}

Pcond_reparam <- function(t, yl, r){
  1/
    (1+
       yl/r*(exp(r*t)-1))
}
