#' Bayesian fit of the birth-death rho model on a phylogeny or a set of phylogenies
#'
#' @description Fits the birth-death rho model (a constant-time homogeneous birth-death model under Bernoulli sampling) to a rooted ultrametric phylogeny using Bayesiann inference. The birth-death process is conditioned on the starting time of the process \code{tot_time} and the survival of the process at present time. The inference can be done specifying the sampling probability or integrating over it according to a specified sampling probability distribution (either uniform: \code{unif = TRUE} or beta distribution: \code{beta = TRUE}). This function can fit the birth-death rho model on a stem or crown phylogeny or a set of phylogenies assuming common or specific diversification rates. It is by default parametrised on the net diversification rate and the turnover rate but can be reparametrised on the product \eqn{y*\lambda} and the net diversification rate \eqn{r}. This function is specifically adapted for diversification analysis on phylogenies on which the sampling probability is unknown.
#'
#' @param phylo Object of class \code{phylo} or \code{multiPhylo}. A rooted ultrametric phylogeny of class \code{phylo} or a set of rooted ultrametric phylogenies of class \code{multiPhylo}.
#' @param tot_time Numeric vector. The stem or crown age (also called MRCA) of the phylogenie(s) depending on the conditioning of the process specified (see \code{cond} argument) and the phylogenie(s) used accordingly. The length of the numeric vector equals the number of phylogenie(s) used. It is of length 1 if a unique phylogeny is used. The stem age of the phylogeny can be computed using max(TreeSim::getx(phylo))+phylo$root.edge (note that the phylo$root.edge needs to be known) and the crown age of the phylogeny can be computed using max(TreeSim::getx(phylo)). If multiple phylogenies are used the following can be used for calculating the stem age: \code{sapply(seq_along(multiPhylo), function(i) max(TreeSim::getx(multiPhylo[[i]]))+multiPhylo[[i]]$root.edge)} and the crown age: \code{sapply(seq_along(multiPhylo), function(i) max(TreeSim::getx(multiPhylo[[i]])))}. Note that if multiple phylogenies are used, the phylogenies do not need to have the same root age but they need to be conditioned the same way (for all phylogenies either their stem or crown age).
#' @param y Numeric vector. The sampling probabilitie(s) typically calculated as \eqn{k/N} where \code{k} is the number of extant sampled tips and \eqn{N} is the global diversity of the clade. The length of the numeric vector equals the number of phylogenie(s) used. If \code{NULL} (default option) and \code{reparam = FALSE}, the sampling probabilitie(s) are integrated according to the specified sampling probability distribution (corresponding to the model \eqn{birth-death∫\rho} in Lambert et al. 2022).
#' @param reparam Logical. If \code{FALSE} (the default option), the log likelihood is calculated using the parameters \eqn{r} and \eqn{\epsilon}. If \code{TRUE} and \code{yj = NULL}, the log likelihood is parametrised using the product \eqn{y*\lambda} and the net diversification rate \eqn{r}.
#' @param common Logical. This argument is only used when a set of phylogenies are provided in the \code{phylo} argument. If \code{TRUE} (default option), common diversification rates are inferred for the set of phylogenies used. If \code{FALSE}, each phylogeny will have its own specific diversification rates inferred.
#' @param beta Logical. This argument is only used if \code{y = NULL} and \code{reparam = FALSE}. If \code{TRUE} a beta distribution is assumed on the sampling probabilitie(s). Note that the parameters of the beta distribution can be fixed or inferred.
#' @param unif Logical. This argument is only used if \code{y = NULL} and \code{reparam = FALSE}. If \code{TRUE} (default option) a uniform distribution is assumed on the sampling probabilitie(s). Note that the parameters of the uniform distribution can be fixed or inferred.
#' @param a Numeric. This argument is only used if \code{y = NULL}, \code{reparam = FALSE} and \code{afix = TRUE}. It corresponds to the value of \eqn{\alpha} (\eqn{\alpha>0}) or \eqn{a} (the lower bound \eqn{0\lea<1}) respectively for the beta or the uniform distribution on the sampling probabilitie(s).
#' @param b Numeric. This argument is only used if \code{y = NULL}, \code{reparam = FALSE} and \code{bfix = TRUE}. It corresponds to the value of \eqn{\beta} (\eqn{\beta>0}) or \eqn{b} (the higher bound (\eqn{0<b\le1} and \eqn{b>a})) respectively for the beta and the uniform distribution on the sampling probabilitie(s).
#' @param afix Logical. This argument is only used if \code{y = NULL} and \code{reparam = FALSE}. If \code{TRUE} (the default option), the hyperparameter \eqn{a} of the model is fixed. If \code{FALSE}, the hyperparameter \eqn{a} of the model is inferred.
#' @param bfix Logical. This argument is only used if \code{y = NULL} and \code{reparam = FALSE}. If \code{TRUE} (the default option), the hyperparameter \eqn{b} of the model is fixed. If \code{FALSE}, the hyperparameter \eqn{b} of the model is inferred.
#' @param cond Character. Specifying the conditioning of the birth-death process. Two conditioning are available, either \code{cond = "crown"} (the default option) if the phylogeny used is a crown phylogeny or \code{cond = "stem"} if the phylogeny used is a stem phylogeny. Note that if a set of phylogenies are used, they will be conditioned the same way according to this argument.
#' @param YULE Logical. If \code{TRUE}, the extinction rate \eqn{\mu} thus the turnover rate \eqn{\epsilon} are fixed to \eqn{0} and the net diversification rate \eqn{r} equals the speciation rate \eqn{\lambda}. If \code{FALSE} (the default option), the turnover rate \eqn{\epsilon} is not fixed to \eqn{0} and is thus inferred. This option is not available if the model is reparametrised (\code{reparam = TRUE}).
#' @param dt Numeric. This argument is only used if \code{y = NULL} and \code{reparam = FALSE}. If \code{dt = 0}, the integral on the sampling probabilitie(s) is computed using the R \code{stats::integrate} function. If \eqn{dt\ge0}, the integral of the sampling probabilitie(s) is performed manually using a piece-wise constant approximation. \code{dt} represents the length of the interval on which the function integrated is assumed to be constant. For manual integral, advised value of \code{dt} are \eqn{1e-3} to \eqn{1e-5}.
#' @param rel.tol Numeric. This argument is only used if \code{y = NULL}, \code{reparam = FALSE} and \code{dt = 0}. This represents the relative accuracy requested when the integral is performed using the \code{stats::integrate} function. Typically \code{.Machine$double.eps^0.25} is used but a value of \eqn{1e-10} (the default value) has been tested and performs well.
#' @param tuned_dichotomy Logical. This argument is only used if \code{y = NULL} and \code{reparam = FALSE}. If \code{TRUE}, when the log likelihood of the model is equal to non finite value due to approximations, a dichotomy search is performed to find a tuning parameter that will be used for getting a finite value of the log likelihood. If \code{TRUE}, the log likelihood will take longer to calculate. Else if \code{FALSE}, no dichotomy search is performed; if the log likelihood is equal to non finite value due to approximations, the log likelihood will take this non finite value for the corresponding parameters.
#' @param brk Numeric. This argument is only used if \code{y = NULL}, \code{reparam = FALSE} and \code{tuned_dichotomy = TRUE}. The number of steps used in the dichotomy search. Typically the value \eqn{200} is sufficient to avoid non finite values. In some case if the log likelihood is still equal to non finite value, the \code{brk} value \eqn{2000} will be required for more tuning but it will rarely take a larger value.
#' @param savedBayesianSetup BayesianOutput. A BayesianOutput created by \code{fitMCMC_bdRho}. If \code{NULL} (the default option), no previous MCMC run is continued and the Bayesian inference start from scratch. If a \code{BayesianOutput} is provided the Bayesian inference continue the previous MCMC run.
#' @param mcmcSettings List. A list of settings for the Bayesian inference using the sampler \code{DEzs} of the package \code{BayesianTools}. Typically, the number of iterations and the starting values will be specified as the following example: \code{mcmcSettings = list(iterations = 3*nbIter, startValue = startValueMatrix)} where \eqn{3} is the number of chains, \code{nbIter} is the number of iterations and, \code{startValueMatrix} is a matrix containing parameters starting values for the MCMC chains. In this example this matrix takes \eqn{3} rows (one for each chain) and the number of columns equals to the number of parameters to infer. Check \code{\link[BayesianTools:runMCMC]{BayesianTools::runMCMC()}} for more details on the \code{settings} options for the sampler \code{DEzs}.
#' @param prior Prior or function. Either a prior class (see \code{\link[BayesianTools:createPrior]{BayesianTools::createPrior()}}) or a log prior density function.
#' @param parallel Numeric or logical. If \code{FALSE} (the default option), the calculation of the likelihood is not parallelised. If \code{>1}, the calculation of the likelihood is parallelised. Note that parallelising the computation is not always faster. This should be checked and depends on the number of cores used for the parallelisation.
#' @param save_inter Numeric vector. A vector specifying the timings at which the MCMC chains should be saved for checkpointing. This can be computed using the following example : \code{c(seq(from = proc.time()[3], to = proc.time()[3]+maxTime, by = freqTime),stopTime)}. It is particularly useful when launching the inference on a cluster where some time restrictions exist.
#' @param index_saving Factor. A factor specifying the name of the MCMC chains saved during the checkpointing. The MCMC chains will be saved as a RDS file in your working directory and will have the following syntax \code{chainMWindex_saving.RDS}.
#'
#' @details This function will fit different birth-death models depending on the arguments chosen:
#' \itemize{
#'   \item If a unique phylogeny is used and the corresponding sampling probability is given in \code{y}, then the classical birth-death-sampling model is used and the function will infer the net diversification rate and the turnover rate.
#'   \item If a unique phylogeny is used, \code{y = NULL} and the model is set for being reparametrised \code{reparam = TRUE}, then the reparametrised birth-death-sampling model is used and the function will infer the net diversification rate and the product of the sampling probability and speciation rate \eqn{y*\lambda}.
#'   \item If a unique phylogeny is used, \code{y = NULL} and \code{reparam = FALSE}, then the \eqn{birth-death∫\rho} model is used and the function will infer the net diversification rate, the turnover rate, and hyperparameters of the sampling probability distribution \code{a} and \code{b} depending if they are set to be fixed or not (see \code{afix} and \code{bfix} arguments). Make sure the desired sampling probability distribution is chosen (see \code{beta} and \code{unif} arguments). See \code{\link{phi}} for more details).
#'   \item If a unique phylogeny is used and \code{YULE = TRUE}, then the corresponding model will be used with one parameter less since the turnover rate will be fixed to \eqn{0} and thus will not be inferred. Note that this option is not available if the model is reparametrised.
#'   \item If a set of phylogenies are used for fitting the model and \code{common = TRUE}, then the \eqn{birth-death∫\rho_mult} model is used and the function will infer the same number of parameters as \eqn{birth-death∫\rho} depending on whether the hyperparameters are set to be fixed or not (see \code{afix} and \code{bfix} arguments).
#'   \item If a set of phylogenies are used for fitting the model and \code{common = FALSE}, then the \eqn{birth-death∫\rho_mult_x} model is used and the function will infer the specific net diversification rate per phylogenies, the turnover rate per phylogenies and the hyperparameters depending on whether they are set to be fixed or not (see \code{afix} and \code{bfix} arguments).
#' }
#' @details Note that depending on the model chosen, the number of parameters inferred can vary thus the \code{prior} and the \code{mcmcSettings} should be adapted to the number of parameters and the parameters should be ordered as described above. This function is specifically intended to be used on phylogenies with unknown or highly uncertain global diversity estimates (the sampling probability is not known with accuracy). Note that the sampling probability is never estimated and that this function is not able to evaluate negative rates.
#'
#' @return Returns an object of class \code{MCMC_bd}. This \code{MCMC_bd} object is a list containing the name of the birth-death model performed and an object of class \code{"mcmcSampler"    "bayesianOutput"} (see the output of \code{\link[BayesianTools:runMCMC]{BayesianTools::runMCMC()}}). This second object contains the MCMC chains and the information about the MCMC run. For analysis of the chains, it can be converted to a \code{coda} object (\code{\link[BayesianTools:getSample]{BayesianTools::getSample()}}) or used in line with the appropriate functions e.g. \code{\link[BayesianTools:MAP]{BayesianTools::MAP()}}.
#'
#' @author Sophia Lambert
#'
#' @example /inst/examples/ExfitMCMC_bdRho.R
#'
#' @export
#'
#' @seealso \code{\link{likelihood_bdRho}} and \code{\link{fitMCMC_bdK}}


fitMCMC_bdRho <- function(phylo, tot_time, y = NULL,
                          reparam = FALSE, common = TRUE,
                          beta = FALSE, unif = TRUE,
                          a = 0, b = 1,
                          afix = TRUE, bfix =TRUE,
                          cond = "crown",
                          YULE = FALSE, dt = 0,
                          rel.tol = 1e-10,
                          tuned_dichotomy = TRUE,
                          brk = 2000,
                          savedBayesianSetup = NULL,
                          mcmcSettings = NULL,
                          prior = NULL,
                          parallel = FALSE,
                          save_inter = NULL,
                          index_saving = NULL){

  # TO DO : rel.tol, if the inference stop by itself it means that the rel.tol is probably too low maybe add an option if the integrate is NULL and does not exist
  # same thing for save_inter with the checkpointing

  if (!inherits(phylo, c("phylo", "multiPhylo")))
    stop("object \"phylo\" is not of class \"phylo\" or \"multiPhylo\"") # check the stops

  # add if common only one starting value
  # if beta = T cannot start with a = 0 or b = 0

  if(class(phylo)=="phylo"){
    phylo <- list(phylo)
  }
  seqphy <- seq_along(phylo)
  len <- length(seqphy)
  ## put a warning message if seqphy is not the same length as tot_time, y, r and epsi, yl
  nbtips <- sapply(phylo, FUN = ape::Ntip)
  from_past <- lapply(seqphy, function (i)
    cbind(phylo[[i]]$edge, picante::node.age(phylo[[i]])$ages))
  ages <- lapply(seqphy, function(i)
    rbind(from_past[[i]][, 2:3], c(nbtips[i] + 1, 0)))
  ages <- lapply(seqphy, function(i)
    ages[[i]][order(ages[[i]][, 1]), ])
  age <- sapply(seqphy, function(i)
    max(ages[[i]][, 2]))
  if (cond=="crown")
  {node_start <- 2
  root <- 1}
  else if (cond=="stem")
  {node_start <- 1
  root <- 0}
  j <- lapply(seqphy, function(i)
    node_start:(as.numeric(nbtips[i]) - 1))
  node <- lapply(seqphy, function(i)
    (as.numeric(nbtips[i]) + j[[i]]))
  tj <- lapply(seqphy, function(i)
    age[i] - ages[[i]][node[[i]], 2])
  chainMWvector <- c()

  likeli <- likelihood_bdRho(tottime = tot_time, nbtips = nbtips,
                             tj = tj, yj = y, reparam = reparam,
                             beta = beta, unif = unif, root = root,
                             dt = dt, rel.tol = rel.tol,
                             tuned_dichotomy = tuned_dichotomy, brk = brk)


  # when we know the sampling fraction

  if(class(y)=="numeric"){
    message("known sampling fraction(s) used")
    if(reparam==TRUE | beta==TRUE | unif==TRUE | class(a)=="numeric" | class(b)=="numeric")
      message("no arguments related to unknown sampling fraction are used here")
    if (YULE==FALSE)
    {
      p <- length(prior$lower)

      # and we want to force multi phylo to have the same values

      if(common == TRUE){

        if(p!=2)
          stop("if common = TRUE only two parameters *r* and *epsi* should be inferred and thus have its associated prior")

        optimLH <- function(init)
        {
          r <- rep(init[1:(p/2)], len)
          epsi <- rep(init[(1+(p/2)):p], len)
          LH <- likeli(div = r, turn = epsi)
          return(LH)
        }

        model.bd <- "bd"

        # or not the same

      }else if (common == FALSE){

        if(p!=c(len*2))
          stop("if common = FALSE each phylogeny should have two parameters *r* and *epsi* thus have its associated prior")

        optimLH <- function(init)
        {
          r <- init[1:(p/2)]
          epsi <- init[(1+(p/2)):p]
          LH <- likeli(div = r, turn = epsi)
          return(LH)
        }

        model.bd <- "bd.mult"

      }
      names = c(paste("div rate", 1:(p/2)), paste("epsi rate", 1:(p/2)))
    }

    # if you know the sampling fraction and you fix the turnover rate to be 0 (Yule model)

    else # YULE = TRUE means fix.mu = TRUE and mu = 0 the r parameter becomes the lambda
    {
      message('extinction rate(s) and turnover rate(s) "epsi" are fixed to 0')
      p <- length(prior$lower)
      # and we want to force multi phylo to have the same values

      if(common == TRUE){

        if(p!=1)
          stop("if YULE = TRUE and common = TRUE only one parameter *r* should be inferred and thus have its associated prior")

        optimLH <- function(init)
        {
          LH <- likeli(div = rep(init, len), turn = rep(0, len))
          return(LH)
        }

        model.bd <- "Yule"

        # or not the same

      }else if (common == FALSE){

        if(p!=c(len))
          stop("if YULE = TRUE and common = FALSE each phylogeny should have one parameter *r* thus have its associated prior")

        optimLH <- function(init)
        {
          LH <- likeli(div = init, turn = rep(0, p))
          return(LH)
        }

        model.bd <- "Yule.mult"

      }
      names = c(paste("div rate", 1:p))
    }
  }

  # if you want to use the model with reparametrisation

  else if(reparam == TRUE){
    if (YULE==FALSE)
    {
      p <- length(prior$lower)

      # and we want to force multi phylo to have the same values

      if(common == TRUE){

        if(p!=2)
          stop("if common = TRUE only two parameters *r* and \"ylamb\" should be inferred and thus have its associated prior")

        optimLH <- function(init)
        {
          r <- rep(init[1:(p/2)], len)
          yl <- rep(init[(1+(p/2)):p], len)
          LH <- likeli(div = r, ylamb = yl)
          return(LH)
        }

        model.bd <- "bd.reparam"

        # or not the same

      }else if (common == FALSE){

        if(p!=c(len*2))
          stop("if common = FALSE each phylogeny should have two parameters *r* and \"ylamb\" thus have its associated prior")

        optimLH <- function(init)
        {
          r <- init[1:(p/2)]
          yl <- init[(1+(p/2)):p]
          LH <- likeli(div = r, ylamb = yl)
          return(LH)
        }

        model.bd <- "bd.mult.reparam"

      }
    }
    else # YULE = TRUE means fix.mu = TRUE and mu = 0 the r parameter becomes the lambda
      stop("YULE model for reparam = TRUE is not available") # check the stop

    names = c(paste("div rate", 1:(p/2)),paste("ylamb", 1:(p/2)))
  }


  # if we want to use the model with rho integrated
  else{
    tun.init <- rep(log(1), len)
    iter <- 0
    # burnin <- 2000

    if (YULE==FALSE)
    {

      # and we are fixing a and b

      if(afix == TRUE & bfix == TRUE){
        p <- length(prior$lower)

        # if(p<1)
        #   stop(*r* and *epsi* priors value should be defined for this model")

        # and we want to force multi phylo to have the same values

        if(common == TRUE){

          if(p!=2)
            stop("if common = TRUE only two parameters *r* and *epsi* should be inferred and thus have its associated prior")

          optimLH <- function(init)
          {
            iter <<- iter+1
            r <- rep(init[1:(p/2)], len)
            epsi <- rep(init[(1+(p/2)):p], len)
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = epsi, a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            # if(iter>burnin) tun.init <<- bfLH$tun # memoryless tuning since the dichotomy search is faster
            return(LH)
          }

          model.bd <- "bd.Rho_fix_a_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len*2))
            stop("if common = FALSE each phylogeny should have two parameters *r* and *epsi* thus have its associated prior")

          optimLH <- function(init)
          {
            iter <<- iter+1
            r <- init[1:(p/2)]
            epsi <- init[(1+(p/2)):p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = epsi, a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            # if(iter>burnin) tun.init <<- bfLH$tun # modified
            return(LH)
          }

          model.bd <- "bd.mult.Rho_fix_a_b"

        }
        names = c(paste("div rate", 1:(p/2)), paste("epsi rate", 1:(p/2)))
      }

      # and we are inferring a

      else if(afix == FALSE & bfix == TRUE){
        p <- length(prior$lower)

        # and we want to force multi phylo to have the same values

        if(common == TRUE){

          if(p!=3)
            stop("if common = TRUE and afix = FALSE three parameters *r*, *epsi* and *a* should be inferred and thus have a prior")

          optimLH <- function(init)
          {
            iter <<- iter+1
            r <- rep(init[1:((p-1)/2)], len)
            epsi <- rep(init[(1+((p-1)/2)):(p-1)], len)
            a <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = epsi, a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }

          model.bd <- "bd.Rho_fix_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len*2+1))
            stop("if common = FALSE and afix = FALSE each phylogeny should have two parameters *r*, *epsi* and one common *a* and thus have its associated prior")

          optimLH <- function(init)
          {
            r <- init[1:((p-1)/2)]
            epsi <- init[(1+((p-1)/2)):(p-1)]
            a <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = epsi, a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }

          model.bd <- "bd.mult.Rho_fix_b"

        }
        names = c(paste("div rate", 1:((p-1)/2)), paste("epsi rate", 1:((p-1)/2)), paste("a", 1))
      }

      # and if we are infering b

      else if(afix == TRUE & bfix == FALSE){
        p <- length(prior$lower)
        # and we want to force multi phylo to have the same values

        if(common == TRUE){

          if(p!=3)
            stop("if common = TRUE and bfix = FALSE three parameters *r*, *epsi* and *b* should be inferred and thus have a prior")

          optimLH <- function(init)
          {
            r <- rep(init[1:((p-1)/2)], len)
            epsi <- rep(init[(1+((p-1)/2)):(p-1)], len)
            b <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = epsi, a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }

          model.bd <- "bd.Rho_fix_a"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len*2+1))
            stop("if common = FALSE and bfix = FALSE each phylogeny should have two parameters *r*, *epsi* and one common *b* and thus have its associated prior")

          optimLH <- function(init)
          {
            r <- init[1:((p-1)/2)]
            epsi <- init[(1+((p-1)/2)):(p-1)]
            b <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = epsi, a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }

          model.bd <- "bd.mult.Rho_fix_a"

        }
        names = c(paste("div rate", 1:((p-1)/2)), paste("epsi rate", 1:((p-1)/2)), paste("b", 1))
      }

      # if we are infering both a and b

      else if(afix == FALSE & bfix == FALSE){
        p <- length(prior$lower)

        # and we want to force multi phylo to have the same values

        if(common == TRUE){

          if(p!=4)
            stop("if common = TRUE, afix = FALSE and bfix = FALSE four parameters *r*, *epsi*, *a* and *b* should be inferred and thus have a prior")

          optimLH <- function(init)
          {
            r <- rep(init[1:((p-2)/2)], len)
            epsi <- rep(init[(1+((p-2)/2)):(p-2)], len)
            a <- init[p-1]
            b <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = epsi, a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }

          model.bd <- "bd.Rho_inf_a_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len*2+2))
            stop("if common = FALSE, afix = FALSE and bfix = FALSE each phylogeny should have two parameters *r*, *epsi* and one common *a* and *b* and thus have its associated prior")

          optimLH <- function(init)
          {
            r <- init[1:((p-2)/2)]
            epsi <- init[(1+((p-2)/2)):(p-2)]
            a <- init[p-1]
            b <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = epsi, a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }

          model.bd <- "bd.mult.Rho_inf_a_b"

        }
        names = c(paste("div rate", 1:((p-2)/2)), paste("epsi rate", 1:((p-2)/2)), paste("a", 1), paste("b", 1))
      }
    }

    # if we do not know the sampling fraction and we force the turnover rate to be 0

    else # YULE = TRUE means fix.mu = TRUE and mu = 0 the d parameter becomes the lambda
    {
      message('extinction rate(s) and turnover rate(s) "epsi" are fixed to 0')

      # if we fix a and b

      if(afix == TRUE & bfix == TRUE){
        p <- length(prior$lower)

        # and we want to force multi phylo to have the same values

        if(common == TRUE){

          if(p!=1)
            stop("if YULE = TRUE and common = TRUE only one parameter *r* should be inferred and thus have its associated prior")

          optimLH <- function(init)
          {
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = rep(init, len), turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.Rho_fix_a_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len))
            stop("if YULE = TRUE and common = FALSE each phylogeny should have one parameter *r* and thus have its associated prior")

          optimLH <- function(init)
          {
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = init, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }

          model.bd <- "Yule.mult.Rho_fix_a_b"

          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")
        }
        names = c(paste("div rate", 1:p))
      }

      # if a is inferred

      else if(afix == FALSE & bfix == TRUE){

        if(p!=2)
          stop("if YULE = TRUE, common = TRUE and afix = FALSE two parameters *r* and *a* should be inferred and thus have a prior")

        p <- length(prior$lower)

        # and we want to force multi phylo to have the same values

        if(common == TRUE){
          optimLH <- function(init)
          {
            r <- rep(init[1:(p-1)], len)
            a <- rep(init[p], len)
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.Rho_fix_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len+1))
            stop("if YULE = TRUE, common = FALSE and afix = FALSE each phylogeny should have one parameter *r* and one common *a* and thus have its associated prior")

          optimLH <- function(init)
          {
            r <- init[1:(p-1)]
            a <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.mult.Rho_fix_b"

        }
        names = c(paste("div rate", 1:(p-1)), paste("a", 1))
      }

      # if b is inferred

      else if(afix == TRUE & bfix == FALSE){
        p <- length(prior$lower)

        # and we want to force multi phylo to have the same values

        if(common == TRUE){

          if(p!=2)
            stop("if YULE = TRUE, common = TRUE and bfix = FALSE two parameters *r* and *b* should be inferred and thus have a prior")

          optimLH <- function(init)
          {
            r <- rep(init[1:(p-1)], len)
            b <- rep(init[p], len)
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.Rho_fix_a"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len+1))
            stop("if YULE = TRUE, common = FALSE and bfix = FALSE each phylogeny should have one parameter *r* and one common *b* and thus have its associated prior")

          optimLH <- function(init)
          {
            r <- init[1:(p-1)]
            b <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.mult.Rho_fix_a"

        }
        names = c(paste("div rate", 1:(p-1)), paste("b", 1))
      }

      # if we infer both a and b

      else if(afix == FALSE & bfix == FALSE){
        p <- length(prior$lower)

        # and we want to force multi phylo to have the same values

        if(common == TRUE){

          if(p!=3)
            stop("if YULE = TRUE, common = TRUE, afix = FALSE and bfix = FALSE three parameters *r*, *a* and *b* should be inferred and thus have a prior")

          optimLH <- function(init)
          {
            r <- rep(init[1:(p-2)], len)
            a <- rep(init[p-1], len)
            b <- rep(init[p], len)
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.Rho_inf_a_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len+2))
            stop("if YULE = TRUE, common = FALSE, afix = FALSE and bfix = FALSE each phylogeny should have one parameter *r* and one common *a* and *b* and thus have its associated prior")

          optimLH <- function(init)
          {
            r <- init[1:(p-2)]
            a <- init[p-1]
            b <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.mult.Rho_inf_a_b"

        }
        names = c(paste("div rate", 1:(p-2)), paste("a", 1), paste("b", 1))
      }
    }
  }

  LL <- function(par) {optimLH(par)}
  if(is.null(prior)) prior = BayesianTools::createUniformPrior(lower = rep(0,p), upper = rep(10,p))
  if(is.null(savedBayesianSetup)) savedBayesianSetup = BayesianTools::createBayesianSetup(likelihood = LL, prior = prior, names = names, parallel = parallel)
  if(("bayesianOutput" %in% class(savedBayesianSetup)) & (length(class(savedBayesianSetup))==1)) {
    chain <- coda::as.mcmc.list(lapply(1:length(savedBayesianSetup$saveChain),function(i) coda::as.mcmc(savedBayesianSetup$saveChain[[i]])))
    savedBayesianSetup = list(
      setup = savedBayesianSetup$setup,
      settings = savedBayesianSetup$settings,
      chain = chain,
      codaChain = coda::mcmc(chain),
      X = savedBayesianSetup$X,
      Z = savedBayesianSetup$Z,
      sampler = savedBayesianSetup$sampler)
    class(savedBayesianSetup) <- c("mcmcSampler", "bayesianOutput")
  }
  MCMCres <- runMCMC(bayesianSetup = savedBayesianSetup, sampler = "DEzs", settings = mcmcSettings, save_inter = save_inter, index_saving = index_saving)
  res <- list(model = model.bd, mcmc = MCMCres)
  class(res) <- "MCMC_bd"
  return(res)
}
