#' Conditioning the birth-death process on a certain time and the survival of the process.
#'
#' @description Calculating the probability of conditioning the birth-death process on the time \eqn{t} and its survival to present day.
#'
#' @param t Numeric. The specific time in the past in millions of years. Typically, the age of the clade is used to condition on the start of the process.
#' @param r Numeric. The net diversification rate: \eqn{\lambda - \mu}. The units are in events/lineages/millions years.
#' @param epsi Numeric. The turnover rate: \eqn{\mu / \lambda}. The units are in events/lineages/millions years.
#' @param y Numeric. The sampling probability also called sampling fraction: \eqn{k / N}.
#'
#' @return The probability of the birth-death process starting with one lineage at time \eqn{t} and surviving at present.
#'
#' @author Sophia Lambert
#'
#' @export
#'
#' @seealso \code{\link{likelihood_bdRho}}, \code{\link{fitMCMC_bdRho}}

Pcond <- function(t, r, epsi, y){
  1/
    (1+
       y/(1-epsi)*(exp(r*t)-1))
}
