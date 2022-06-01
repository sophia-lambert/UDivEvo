#' Bayesian fit of the birth-death k model on a phylogeny
#'
#' @description Fits the birth-death k model (a constant-time homogeneous birth-death model under k-sampling) to a rooted ultrametric phylogeny using Bayesiann inference. The birth-death process is conditioned on the starting time of the process \code{tot_time} and the survival of the process at present time as well as having \eqn{k} extant sampled tips at present. This function can fit the birth-death k model on a stem or crown phylogeny. This function is specifically adapted for diversification analysis on phylogenies on which the sampling probability is unknown.
#'
#' @param phylo Object of class \code{phylo}. A rooted ultrametric phylogeny of class \code{phylo}.
#' @param tot_time Numeric. The stem or crown age (also called MRCA) of the phylogeny depending on the conditioning of the process specified (see \code{cond} argument) and the phylogeny used accordingly. The stem age of the phylogeny can be computed using max(TreeSim::getx(phylo))+phylo$root.edge (note that the phylo$root.edge needs to be known) and the crown age of the phylogeny can be computed using max(TreeSim::getx(phylo)).
#' @param cond Character. Specifying the conditioning of the birth-death process. Two conditioning are available, either \code{cond = "crown"} (the default option) if the phylogeny used is a crown phylogeny or \code{cond = "stem"} if the phylogeny used is a stem phylogeny.
#' @param YULE Logical. If \code{TRUE}, the extinction rate \eqn{\mu} thus the turnover rate \eqn{\epsilon} are fixed to \eqn{0} and the net diversification rate \eqn{r} equals the speciation rate \eqn{\lambda}. If \code{FALSE} (the default option), the turnover rate \eqn{\epsilon} is not fixed to \eqn{0} and is thus inferred.
#' @param dt Numeric. If \code{dt = 0}, the integral on the sampling probability is computed using the R \code{stats::integrate} function. If \eqn{dt\ge0}, the integral of the sampling probability is performed manually using a piece-wise constant approximation. \code{dt} represents the length of the interval on which the function integrated is assumed to be constant. For manual integral, advised value of \code{dt} are \eqn{1e-3} to \eqn{1e-5}.
#' @param rel.tol Numeric. This argument is only used if \code{dt = 0}. This represents the relative accuracy requested when the integral is performed using the \code{stats::integrate} function. Typically \code{.Machine$double.eps^0.25} is used but a value of \eqn{1e-10} (the default value) has been tested and performs well.
#' @param tuned_dichotomy Logical. If \code{TRUE}, when the log likelihood of the model is equal to non finite value due to approximations, a dichotomy search is performed to find a tuning parameter that will be used for getting a finite value of the log likelihood. If \code{TRUE}, the log likelihood will take longer to calculate. Else if \code{FALSE}, no dichotomy search is performed; if the log likelihood is equal to non finite value due to approximations, the log likelihood will take this non finite value for the corresponding parameters.
#' @param brk Numeric. This argument is only used if \code{tuned_dichotomy = TRUE}. The number of steps used in the dichotomy search. Typically the value \eqn{200} is sufficient to avoid non finite values. In some case if the log likelihood is still equal to non finite value, the \code{brk} value \eqn{2000} will be required for more tuning but it will rarely take a larger value.
#' @param savedBayesianSetup BayesianOutput. A BayesianOutput created by \code{fitMCMC_bdK}. If \code{NULL} (the default option), no previous MCMC run is continued and the Bayesian inference start from scratch. If a \code{BayesianOutput} is provided the Bayesian inference continue the previous MCMC run.
#' @param mcmcSettings List. A list of settings for the Bayesian inference using the sampler \code{DEzs} of the package \code{BayesianTools}. Typically, the number of iterations and the starting values will be specified as the following example: \code{mcmcSettings = list(iterations = 3*nbIter, startValue = startValueMatrix)} where \eqn{3} is the number of chains, \code{nbIter} is the number of iterations and, \code{startValueMatrix} is a matrix containing parameters starting values for the MCMC chains. In this example this matrix takes \eqn{3} rows (one for each chain) and the number of columns equals to the number of parameters to infer (here \eqn{2} parameters). Check \code{\link[BayesianTools:runMCMC]{BayesianTools::runMCMC()}} for more details on the \code{settings} options for the sampler \code{DEzs}.
#' @param prior Prior or function. Either a prior class (see \code{\link[BayesianTools:createPrior]{BayesianTools::createPrior()}}) or a log prior density function.
#' @param parallel Numeric or logical. If \code{FALSE} (the default option), the calculation of the likelihood is not parallelised. If \code{>1}, the calculation of the likelihood is parallelised. Note that parallelising the computation is not always faster. This should be checked and depends on the number of cores used for the parallelisation.
#' @param save_inter Numeric vector. A vector specifying the timings at which the MCMC chains should be saved for checkpointing. This can be computed using the following example : \code{c(seq(from = proc.time()[3], to = proc.time()[3]+maxTime, by = freqTime),stopTime)}. It is particularly useful when launching the inference on a cluster where some time restrictions exist.
#' @param index_saving Factor. A factor specifying the name of the MCMC chains saved during the checkpointing. The MCMC chains will be saved as a RDS file in your working directory and will have the following syntax \code{chainMWindex_saving.RDS}.
#'
#' @details This function will fit the birth-death k-sampling model and the function will infer the net diversification rate \eqn{r} and the turnover rate \eqn{\epsilon}. Note that the \code{prior} and the \code{mcmcSettings} should be adapted to the number of parameters (here \eqn{2}, \eqn{r} and \eqn{\epsilon}). This function is specifically intended to be used on phylogenies with unknown or highly uncertain global diversity estimates (the sampling probability is not known with accuracy). Note that the sampling probability is never estimated and that this function is not able to evaluate negative rates.
#'
#' @return Returns an object of class \code{MCMC_bdK}. This \code{MCMC_bdK} object is a list containing the name of the birth-death model performed and an object of class \code{"mcmcSampler"    "bayesianOutput"} (see the output of \code{\link[BayesianTools:runMCMC]{BayesianTools::runMCMC()}}). This second object contains the MCMC chains and the information about the MCMC run. For analysis of the chains, it can be converted to a \code{coda} object (\code{\link[BayesianTools:getSample]{BayesianTools::getSample()}}) or used in line with the appropriate functions e.g. \code{\link[BayesianTools:MAP]{BayesianTools::MAP()}}.
#'
#' @author Sophia Lambert
#'
#' @export
#'
#' @seealso \code{\link{likelihood_bdK}} and \code{\link{fitMCMC_bdRho}}


fitMCMC_bdK <- function(phylo, tot_time,
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

  if (!inherits(phylo, "phylo"))
    stop("object \"phylo\" is not of class \"phylo\"") # check stops

  nbtip <- ape::Ntip(phylo)
  from_past <- cbind(phylo$edge, picante::node.age(phylo)$ages)
  ages <- rbind(from_past[, 2:3], c(nbtip + 1, 0))
  age <- max(ages[, 2])
  tj <- list()
  if (cond=="crown")
  {
    root <- 2
    ntot <- nbtip-1+nbtip
    select <- which(from_past[,1]==nbtip+1)[2]
    nbtips <- c()
    nbtips[1] <- ((length(ages[1:select-1,])/2)+1)/2
    nbtips[2] <- ((length(ages[select:(ntot-1),])/2)+1)/2
    out <-  which(ages[,1]%in%c(1:nbtip))
    ages_nodes <- ages[-out,]
    node <- list()
    node[[1]] <- 1:(nbtips[1]-1)
    node[[2]] <- nbtips[1]:(nbtips[1]+nbtips[2]-2)
    node[sapply(node, is.unsorted)] <- 0 # essential : if there is only one tip in a subphylo there is no tj for it
    tj[[1]] <- age - ages_nodes[node[[1]], 2]
    tj[[2]] <- age - ages_nodes[node[[2]], 2]
  }

  else if (cond=="stem")
  {
    root <- 1
    nbtips <- nbtip
    j <- 1:(as.numeric(nbtips) - 1)
    node <- as.numeric(nbtips) + j
    ages <- ages[order(ages[, 1]), ]
    tj[[1]] <- age - ages[node, 2]
  }

  likeli <- likelihood_bdK(tottime = tot_time, nbtips, tj,
                           root, dt = dt, rel.tol = rel.tol,
                           tuned_dichotomy = tuned_dichotomy,
                           brk = brk)

  # if we want to use the model with K sampling
  tun.init <- rep(log(1), root)
  if (YULE==FALSE)
  {
    p <- length(prior$lower)

    optimLH <- function(init)
    {
      r <- init[1:length(p/2)]
      epsi <- init[(1 + length(p/2)):length(init)]
      tuningLH <- likeli(tun.init, seqphy.init = 1:root)
      bfLH <- tuningLH(div = r, turn = epsi)
      LH <- bfLH$logLik
      # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
      return(LH)
    }
    names = c(paste("div rate", 1:length(p/2)), paste("epsi rate", 1:length(p/2)))
  }

  # if we want to use the model with K sampling and we force the turnover rate to be 0

  else # YULE = TRUE means fix.mu = TRUE and mu = 0 the d parameter becomes the lambda
  {
    message('extinction rate(s) and turnover rate(s) "epsi" are fixed to 0')

    p <- length(prior$lower)

    optimLH <- function(init)
    {
      tuningLH <- likeli(tun.init, seqphy.init = 1:root)
      bfLH <- tuningLH(div = init, turn = 0)
      LH <- bfLH$logLik
      # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
      return(LH)
    }
    names = c(paste("div rate", 1:length(p)))
  }
  LL <- function(par) {optimLH(par)}
  if(is.null(prior)) prior = BayesianTools::createUniformPrior(lower = rep(0,p), upper = rep(10,p))
  if(is.null(savedBayesianSetup)) savedBayesianSetup = BayesianTools::createBayesianSetup(likelihood = LL, prior = prior, names = names, parallel = parallel) # modified
  if(("bayesianOutput" %in% class(savedBayesianSetup)) & (length(class(savedBayesianSetup))==1)) {
    chain <- coda::as.mcmc.list(lapply(1:length(savedBayesianSetup$saveChain),function(i) coda::as.mcmc(savedBayesianSetup$saveChain[[i]]))) # change Npop
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
  res <- list(model = "bd.K", mcmc = MCMCres)
  class(res) <- "MCMC_bdK"
  return(res)
}

