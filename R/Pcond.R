#' Calculating the conditioning of the birth-death process on the age of the process and its survival.
#'
#' @param t Numeric. The time of the start of the birth-death process in millions of years.
#' @param r Numeric. The net diversification rate : \eqn{\lambda - \mu}. The units are in events/lineages/millions years.
#' @param epsi Numeric. The turnover rate : \eqn{\mu / \lambda}. The units are in events/lineages/millions years.
#' @param y Numeric. The sampling probability also called sampling fraction : \eqn{k / N}.
#'
#' @return The probability of the birth-death process starting with one lineage at time t and surviving at present.
#'
#' @author Sophia Lambert
#'

#@seealso \code{\link{likelihood_bdRho}}, \code{\link{fitMCMC_bdRho}}


Pcond <- function(t, r, epsi, y){
  1/
    (1+
       y/(1-epsi)*(exp(r*t)-1))
}
