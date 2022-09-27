#' Maximum likelihood fit of the birth-death rho model on a phylogeny or a set of phylogenies
#'
#' @description Fits the birth-death rho model (a constant-time homogeneous birth-death model under Bernoulli sampling) to a rooted ultrametric phylogeny using Maximum likelihood estimation. Precisely, it uses the Nelder-Mead algorithm to optimise the likelihood. The birth-death process is conditioned on the starting time of the process \code{tot_time} and the survival of the process at present time. The inference can be done specifying the sampling probability or integrating over it according to a specified sampling probability distribution (either uniform: \code{unif = TRUE} or beta distribution: \code{beta = TRUE}). This function can fit the birth-death rho model on a stem or crown phylogeny or a set of phylogenies assuming common or specific diversification rates. It is by default parametrised on the net diversification rate and the turnover rate but can be reparametrised on the product \eqn{y*\lambda} and the net diversification rate \eqn{r}. This function is specifically adapted for diversification analysis on phylogenies on which the sampling probability is unknown.
#'
#' @param phylo Object of class \code{phylo} or \code{multiPhylo}. A rooted ultrametric phylogeny of class \code{phylo} or a set of rooted ultrametric phylogenies of class \code{multiPhylo}. The rooted ultrametric phylogenie(s) can have polytomie(s) (i.e. non binary tree).
#' @param tot_time Numeric vector. The stem or crown age (also called MRCA) of the phylogenie(s) depending on the conditioning of the process specified (see \code{cond} argument) and the phylogenie(s) used accordingly. The length of the numeric vector equals the number of phylogenie(s) used. It is of length 1 if a unique phylogeny is used. The stem age of the phylogeny can be computed using max(TreeSim::getx(phylo))+phylo$root.edge (note that the phylo$root.edge needs to be known) and the crown age of the phylogeny can be computed using max(TreeSim::getx(phylo)). If multiple phylogenies are used the following can be used for calculating the stem age: \code{sapply(seq_along(multiPhylo), function(i) max(TreeSim::getx(multiPhylo[[i]]))+multiPhylo[[i]]$root.edge)} and the crown age: \code{sapply(seq_along(multiPhylo), function(i) max(TreeSim::getx(multiPhylo[[i]])))}. Note that if multiple phylogenies are used, the phylogenies do not need to have the same root age but they need to be conditioned the same way (for all phylogenies either their stem or crown age).
#' @param r Numeric vector. The net diversification rate \eqn{r} starting value(s).
#' @param epsi Numeric vector. The turnover rate \eqn{\epsilon} starting value(s).
#' @param y Numeric vector. The sampling probabilitie(s) typically calculated as \eqn{k/N} where \code{k} is the number of extant sampled tips and \eqn{N} is the global diversity of the clade. The length of the numeric vector equals the number of phylogenie(s) used. If \code{NULL} (default option) and \code{reparam = FALSE}, the sampling probabilitie(s) are integrated according to the specified sampling probability distribution (corresponding to the model \eqn{birth-death∫\rho} in Lambert et al. 2022).
#' @param ylamb Numeric vector. The weighted speciation rate \eqn{\y*\lambda} starting value(s).
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
#' @details Note that depending on the model chosen, the number of parameters inferred can vary thus the starting values should be adapted to the number of parameters and the parameters should be ordered as described above. This function is specifically intended to be used on phylogenies with unknown or highly uncertain global diversity estimates (the sampling probability is not known with accuracy). Note that the sampling probability is never estimated and that this function is not able to evaluate negative rates.
#'
#' @return Returns an object of class \code{MLE_bd}. This \code{MLE_bd} object is a list containing the name of the birth-death model performed, the log likelihood of the data knowing the parameters, the akaike information criterion corrected and the inferred parameters.
#'
#' @author Sophia Lambert
#'
#' @example /inst/examples/ExfitMLE_bdRho.R
#'
#' @export
#'
#' @seealso \code{\link{likelihood_bdRho}} and \code{\link{fitMCMC_bdRho}}

fitMLE_bdRho <- function(phylo, tot_time, r = NULL, epsi = NULL,
                         y = NULL, ylamb = NULL,
                         reparam = FALSE, common = TRUE,
                         beta = FALSE, unif = TRUE,
                         a = 0, b = 1,
                         afix = TRUE, bfix =TRUE,
                         cond = "crown",
                         YULE = FALSE, dt = 0,
                         rel.tol = 1e-10,
                         tuned_dichotomy = TRUE,
                         brk = 2000){

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
  counts <- lapply(seqphy, function (i)
    table(phylo[[i]]$edge[,1])) # handling polytomies
  polytom_nodes <- lapply(seqphy, function (i)
    as.numeric(names(counts[[i]])[counts[[i]] > 2]))
  polytomTimes <- lapply(seqphy, function (i)
    counts[[i]][counts[[i]] > 2]-1)
  ages <- lapply(seqphy, function (i)
    cbind(ages[[i]], 1))
  ages <- lapply(seqphy, function(i){
    ages[[i]][polytom_nodes[[i]],3] <- polytomTimes[[i]]; ages[[i]]})
  ages_polytom <- lapply(seqphy, function (i)
    as.data.frame(lapply(as.data.frame(ages[[i]]), rep, ages[[i]][,3])))
  age <- sapply(seqphy, function(i)
    max(ages_polytom[[i]][, 2]))
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
    age[i] - ages_polytom[[i]][node[[i]], 2])
  chainMWvector <- c()

  likeli <- likelihood_bdRho(tottime = tot_time, nbtips = nbtips,
                             tj = tj, yj = y, reparam = reparam,
                             beta = beta, unif = unif, root = root,
                             dt = dt, rel.tol = rel.tol,
                             tuned_dichotomy = tuned_dichotomy, brk = brk)


  # when we know the sampling fraction

  if(class(y)=="numeric"){
    message("known sampling fraction(s) used")
    if(reparam==TRUE | class(ylamb)=="numeric" | beta==TRUE | unif==TRUE | class(a)=="numeric" | class(b)=="numeric")
      message("no arguments related to unknown sampling fraction are used here")
    if (YULE==FALSE)
    {
      if(is.null(r)|is.null(epsi))
        stop("\"r\" and \"epsi\" starting value should be specified") # checking if the stop actually stops

      startV <- c(r, epsi)
      p <- length(startV)

      # and we want to force multi phylo to have the same values

      if(common == TRUE){

        if(p!=2)
          stop("if common = TRUE only two parameters *r* and *epsi* should be inferred and thus have its associated starting value")

        optimLH <- function(init)
        {
          r <- rep(init[1:(p/2)], len)
          epsi <- rep(init[(1+(p/2)):p], len)
          LH <- likeli(div = r, turn = epsi)
          return(-LH)
        }

        model.bd <- "bd"

        # or not the same

      }else if (common == FALSE){

        if(p!=c(len*2))
          stop("if common = FALSE each phylogeny should have two parameters *r* and *epsi* thus have its associated starting value")

        optimLH <- function(init)
        {
          r <- init[1:(p/2)]
          epsi <- init[(1+(p/2)):p]
          LH <- likeli(div = r, turn = epsi)
          return(-LH)
        }

        model.bd <- "bd.mult"

      }
      names = c(paste("div rate", 1:(p/2)), paste("epsi rate", 1:(p/2)))
    }

    # if you know the sampling fraction and you fix the turnover rate to be 0 (Yule model)

    else # YULE = TRUE means fix.mu = TRUE and mu = 0 the r parameter becomes the lambda
    {
      message('extinction rate(s) and turnover rate(s) "epsi" are fixed to 0')
      if(is.null(r))
        stop("\"r\" starting value (becoming \"lambda\" when YULE = TRUE) should be specified") # checking if the stop actually stops

      startV <- c(r)
      p <- length(startV)
      # and we want to force multi phylo to have the same values

      if(common == TRUE){

        if(p!=1)
          stop("if YULE = TRUE and common = TRUE only one parameter *r* should be inferred and thus have its associated starting value")

        optimLH <- function(init)
        {
          LH <- likeli(div = rep(init, len), turn = rep(0, len))
          return(-LH)
        }

        model.bd <- "Yule"

        # or not the same

      }else if (common == FALSE){

        if(p!=c(len))
          stop("if YULE = TRUE and common = FALSE each phylogeny should have one parameter *r* thus have its associated starting value")

        optimLH <- function(init)
        {
          LH <- likeli(div = init, turn = rep(0, p))
          return(-LH)
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
      if(is.null(r)|is.null(ylamb))
        stop("\"r\" and \"ylamb\" starting value should be specified") # checking if the stop actually stops

      startV <- c(r, ylamb)
      p <- length(startV)

      # and we want to force multi phylo to have the same values

      if(common == TRUE){

        if(p!=2)
          stop("if common = TRUE only two parameters *r* and \"ylamb\" should be inferred and thus have its associated starting value")

        optimLH <- function(init)
        {
          r <- rep(init[1:(p/2)], len)
          yl <- rep(init[(1+(p/2)):p], len)
          LH <- likeli(div = r, ylamb = yl)
          return(-LH)
        }

        model.bd <- "bd.reparam"

        # or not the same

      }else if (common == FALSE){

        if(p!=c(len*2))
          stop("if common = FALSE each phylogeny should have two parameters *r* and \"ylamb\" thus have its associated starting value")

        optimLH <- function(init)
        {
          r <- init[1:(p/2)]
          yl <- init[(1+(p/2)):p]
          LH <- likeli(div = r, ylamb = yl)
          return(-LH)
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

        if(is.null(r)|is.null(epsi)|is.null(a)|is.null(b))
          stop("\"r\", \"epsi\" starting value as well as \"a\" and \"b\" fixed value should be specified") # checking if the stop actually stops

        startV <- c(r, epsi)
        p <- length(startV)

        # if(p<1)
        #   stop(*r* and *epsi* priors value should be defined for this model")

        # and we want to force multi phylo to have the same values

        if(common == TRUE){

          if(p!=2)
            stop("if common = TRUE only two parameters *r* and *epsi* should be inferred and thus have its associated starting value")

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
            return(-LH)
          }

          model.bd <- "bd.Rho_fix_a_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len*2))
            stop("if common = FALSE each phylogeny should have two parameters *r* and *epsi* thus have its associated starting value")

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
            return(-LH)
          }

          model.bd <- "bd.mult.Rho_fix_a_b"

        }
        names = c(paste("div rate", 1:(p/2)), paste("epsi rate", 1:(p/2)))
      }

      # and we are inferring a

      else if(afix == FALSE & bfix == TRUE){

        if(is.null(r)|is.null(epsi)|is.null(a)|is.null(b))
          stop("\"r\", \"epsi\", \"a\" starting value as well as \"b\" fixed value should be specified") # checking if the stop actually stops

        startV <- c(r, epsi, a)
        p <- length(startV)

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
            return(-LH)
          }

          model.bd <- "bd.Rho_fix_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len*2+1))
            stop("if common = FALSE and afix = FALSE each phylogeny should have two parameters *r*, *epsi* and one common *a* and thus have its associated starting value")

          optimLH <- function(init)
          {
            r <- init[1:((p-1)/2)]
            epsi <- init[(1+((p-1)/2)):(p-1)]
            a <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = epsi, a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(-LH)
          }

          model.bd <- "bd.mult.Rho_fix_b"

        }
        names = c(paste("div rate", 1:((p-1)/2)), paste("epsi rate", 1:((p-1)/2)), paste("a", 1))
      }

      # and if we are infering b

      else if(afix == TRUE & bfix == FALSE){

        if(is.null(r)|is.null(epsi)|is.null(a)|is.null(b))
          stop("\"r\", \"epsi\", \"b\" starting value as well as \"a\" fixed value should be specified") # checking if the stop actually stops

        startV <- c(r, epsi, b)
        p <- length(startV)
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
            return(-LH)
          }

          model.bd <- "bd.Rho_fix_a"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len*2+1))
            stop("if common = FALSE and bfix = FALSE each phylogeny should have two parameters *r*, *epsi* and one common *b* and thus have its associated starting value")

          optimLH <- function(init)
          {
            r <- init[1:((p-1)/2)]
            epsi <- init[(1+((p-1)/2)):(p-1)]
            b <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = epsi, a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(-LH)
          }

          model.bd <- "bd.mult.Rho_fix_a"

        }
        names = c(paste("div rate", 1:((p-1)/2)), paste("epsi rate", 1:((p-1)/2)), paste("b", 1))
      }

      # if we are infering both a and b

      else if(afix == FALSE & bfix == FALSE){

        if(is.null(r)|is.null(epsi)|is.null(a)|is.null(b))
          stop("\"r\", \"epsi\", \"a\" and \"b\" starting value should be specified") # checking if the stop actually stops

        startV <- c(r, epsi, a, b)
        p <- length(startV)

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
            return(-LH)
          }

          model.bd <- "bd.Rho_inf_a_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len*2+2))
            stop("if common = FALSE, afix = FALSE and bfix = FALSE each phylogeny should have two parameters *r*, *epsi* and one common *a* and *b* and thus have its associated starting value")

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
            return(-LH)
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

        if(is.null(r)|is.null(a)|is.null(b))
          stop("\"r\" (becoming \"lambda\" when YULE = TRUE) starting value as well as \"a\" and \"b\" fixed value should be specified") # checking if the stop actually stops

        startV <- c(r)
        p <- length(startV)

        # and we want to force multi phylo to have the same values

        if(common == TRUE){

          if(p!=1)
            stop("if YULE = TRUE and common = TRUE only one parameter *r* should be inferred and thus have its associated starting value")

          optimLH <- function(init)
          {
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = rep(init, len), turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(-LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.Rho_fix_a_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len))
            stop("if YULE = TRUE and common = FALSE each phylogeny should have one parameter *r* and thus have its associated starting value")

          optimLH <- function(init)
          {
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = init, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(-LH)
          }

          model.bd <- "Yule.mult.Rho_fix_a_b"

          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")
        }
        names = c(paste("div rate", 1:p))
      }

      # if a is inferred

      else if(afix == FALSE & bfix == TRUE){

        if(is.null(r)|is.null(a)|is.null(b))
          stop("\"r\" (becoming \"lambda\" when YULE = TRUE) and \"a\" starting value as well as \"b\" fixed value should be specified") # checking if the stop actually stops

        startV <- c(r, a)
        p <- length(startV)

        if(p!=2)
          stop("if YULE = TRUE, common = TRUE and afix = FALSE two parameters *r* and *a* should be inferred and thus have a prior")

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
            return(-LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.Rho_fix_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len+1))
            stop("if YULE = TRUE, common = FALSE and afix = FALSE each phylogeny should have one parameter *r* and one common *a* and thus have its associated starting value")

          optimLH <- function(init)
          {
            r <- init[1:(p-1)]
            a <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(-LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.mult.Rho_fix_b"

        }
        names = c(paste("div rate", 1:(p-1)), paste("a", 1))
      }

      # if b is inferred

      else if(afix == TRUE & bfix == FALSE){

        if(is.null(r)|is.null(a)|is.null(b))
          stop("\"r\" (becoming \"lambda\" when YULE = TRUE) and \"b\" starting value as well as \"a\" fixed value should be specified") # checking if the stop actually stops

        startV <- c(r, b)
        p <- length(startV)

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
            return(-LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.Rho_fix_a"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len+1))
            stop("if YULE = TRUE, common = FALSE and bfix = FALSE each phylogeny should have one parameter *r* and one common *b* and thus have its associated starting value")

          optimLH <- function(init)
          {
            r <- init[1:(p-1)]
            b <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(-LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.mult.Rho_fix_a"

        }
        names = c(paste("div rate", 1:(p-1)), paste("b", 1))
      }

      # if we infer both a and b

      else if(afix == FALSE & bfix == FALSE){

        if(is.null(r)|is.null(a)|is.null(b))
          stop("\"r\" (becoming \"lambda\" when YULE = TRUE), \"a\" and \"b\" starting value should be specified") # checking if the stop actually stops

        startV <- c(r, a, b)
        p <- length(startV)

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
            return(-LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.Rho_inf_a_b"

          # or not the same

        }else if (common == FALSE){

          if(p!=c(len+2))
            stop("if YULE = TRUE, common = FALSE, afix = FALSE and bfix = FALSE each phylogeny should have one parameter *r* and one common *a* and *b* and thus have its associated starting value")

          optimLH <- function(init)
          {
            r <- init[1:(p-2)]
            a <- init[p-1]
            b <- init[p]
            tuningLH <- likeli(tun.init, seqphy.init = seqphy)
            bfLH <- tuningLH(div = r, turn = rep(0, len), a = a, b = b)
            LH <- bfLH$logLik
            # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
            return(-LH)
          }
          # if(is.infinite(temp$value)) warning("Yule model is highly improbable")

          model.bd <- "Yule.mult.Rho_inf_a_b"

        }
        names = c(paste("div rate", 1:(p-2)), paste("a", 1), paste("b", 1))
      }
    }
  }

  temp <- stats::optim(startV, optimLH, method = "Nelder-Mead")
  df_params <- data.frame(temp$par, row.names = names)
  res <- list(model = model.bd, LH = -temp$value,
              aicc=2 * temp$value + 2 * p + (2 * p * (p + 1)) / (sum(nbtips) - p - 1),
              parameters = df_params)

  class(res) <- "MLE_bd"
  return(res)
}
