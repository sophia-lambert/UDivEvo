#' Probabillity that all speciation and extinction events happened after a certain time.
#'
#' @description Calculating the probabillity that all speciation and extinction events happened after a certain time \eqn{t}.
#'
#' @param t Numeric The specific time in the past after which all the events happened in millions of years. Typically, the age of the clade is used to condition on the start of the process.
#' @param r Numeric. The net diversification rate: \eqn{\lambda - \mu}. The units are in events/lineages/millions years.
#' @param epsi Numeric. The turnover rate: \eqn{\mu / \lambda}. The units are in events/lineages/millions years.
#'
#' @return The probability that all birth and death events of the birth-death process happened after a time \eqn{t}.
#'
#' @author Sophia Lambert
#'
#' @export
#'
#'@seealso \code{\link{likelihood_bdK}}, \code{\link{fitMCMC_bdK}}

Pa <- function(t, r, epsi){
  1-1/(1+(exp(r*t)-1)/(1-epsi))
}
