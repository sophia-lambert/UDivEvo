#' Integrate function
#'
#' @description Integrate function either from \code{stats::integrate} while overpassing errors or manual integration.
#'
#' @param dt Numeric. The length of the intervals on which the function is assumed to be constant. If \code{dt = 0} the function is integrated using  \code{stats::integrate} while overpassing errors. If a positive \code{dt} is given as argument, the integral is computed using a piece-wise constant approximation.
#' @param f Function. An \code{R} function of one variable on which to integrate.
#' @param lower Numeric. Lower limit of integration.
#' @param upper Numeric. Upper limit of integration.
#' @param rel.tol Numeric. Relative accuracy requested.
#'
#' @details If \code{dt = 0} any other arguments of \code{stats::integrate} can be specified. This function is a closure, it takes \code{dt} as argument and creates another function with arguments \code{f}, \code{lower}, \code{upper} and \code{rel.tol}.
#'
#' @return The final estimate of the integral.
#'
#' @author Sophia Lambert
#'
#' @export

#@seealso \code{\link{likelihood_bdRho}}, \code{\link{fitMCMC_bdRho}}

int_Rho <- function(dt){
  if(dt==0){
    function(...)
    {
      op <- getOption("show.error.messages")
      options(show.error.messages=FALSE)
      res <- try(stats::integrate(...));
      options(show.error.messages=op)
      if (res[1] == "Error in stats::integrate(...) : non-finite function value\n"| # always Inf
          res[1] =="Error in stats::integrate(...) : roundoff error was detected\n"| # always Inf
          res[1] == Inf|
          res[1] == "Error in stats::integrate(...) : the integral is probably divergent\n"| # maybe we should just adjust the rel.tol because it could just be very low values close to 0 (and thus actually should return 0 and not Inf so that the log would be -Inf and not Inf, we have tested that and return 0 and had actually more non finite values at the end) but also Inf empirically
          res[1] == "Error in stats::integrate(...) : maximum number of subdivisions reached\n"| # we could also just adjust the rel.tol or increase the number of subdivisions
          res[1] == "Error in stats::integrate(...) : extremely bad integrand behaviour\n") # this has been checked and it works when gives infinite and not 0 # this should not give an Infinite value so probably just 0 problem due to bad rel.tol (increase it)
        # before we had class(res) == "try-error"
      {
        # cat("We are overtaking the error message of integrate that gives non-finite function value and give an Inf value to the likelihood instead since it is a problem link to Inf")
        return(Inf)
      }
      else
      {
        return(res$value)
      }
    }
  }
  else if(dt>0)
  {
    function(f, lower, upper, rel.tol = NULL){
      int <- 1/dt
      X <- seq(lower, upper, length.out = int + 1)
      res_intr <- sum(sapply(X, f)) / int # could also be cumsum()[int+1]
      return(res_intr)
    }
  }
}
