# A modified version of the BayesianTools::DEzs

DEzs <- function(bayesianSetup, save_inter = NULL, index_saving = NULL,
                 settings = list(iterations=10000,
                                 Z = NULL,
                                 startValue = NULL,
                                 pSnooker = 0.1,
                                 burnin = 0,
                                 thin = 1,
                                 f = 2.38,
                                 eps = 0,
                                 parallel = NULL,
                                 pGamma1 = 0.1,
                                 eps.mult =0.2,
                                 eps.add = 0,
                                 consoleUpdates = 100,
                                 zUpdateFrequency = 1,
                                 currentChain = 1,
                                 blockUpdate = list("none", k = NULL, h = NULL, pSel = NULL, pGroup = NULL,
                                                    groupStart = 1000, groupIntervall = 1000)
                                 ,message = TRUE))
{


  #  X = startValue


  # print("I got inside the function")
  nbSave <- 1

  if("bayesianOutput" %in% class(bayesianSetup)){
    restart <- TRUE
  } else restart <- FALSE


  if(restart){
    if(is.null(settings)) settings <- bayesianSetup$settings
    else  settings <- BayesianTools::applySettingsDefault(settings = settings, sampler = "DEzs")
  }else{
    # If nothing provided use default settings
    settings <- BayesianTools::applySettingsDefault(settings = settings, sampler = "DEzs")
  }

  if(!restart){
    setup <- bayesianSetup
  } else setup <- bayesianSetup$setup

  setup <- BayesianTools::checkBayesianSetup(setup, parallel = settings$parallel) # calling parallel will check if requested parallelization in settings is provided by the BayesianSetup
  if(is.null(settings$parallel)) settings$parallel = setup$parallel # checking back - if no parallelization is provided, we use the parallelization in the BayesianSetup. We could also set parallel = F, but I feel it makes more sense to use the Bayesiansetup as default

  if(!restart){
    if(is.null(settings$startValue)){
      parLen = length(bayesianSetup$prior$sampler(1))
      X = bayesianSetup$prior$sampler(3)
    }
    if(is.function(settings$startValue)){
      X = settings$startValue()
    }
    if(class(settings$startValue)[1] == "numeric"){
      X = bayesianSetup$prior$sampler(settings$startValue)
    }

    if(is.matrix(settings$startValue)) X <- settings$startValue

    if(is.null(settings$Z)){
      parLen = length(bayesianSetup$prior$sampler(1))
      Z = bayesianSetup$prior$sampler(parLen * 10)
    }
    if(is.function(settings$Z)){
      Z = settings$Z()
    }

    if(class(settings$Z)[1] == "numeric"){
      Z = bayesianSetup$prior$sampler(settings$Z)
    }
    if(is.matrix(settings$Z)) Z <- settings$Z

  }else{
    X <- bayesianSetup$X
    Z <- bayesianSetup$Z
    if(is.vector(Z)) Z = as.matrix(Z)
  }


  if (! is.matrix(X)) stop("wrong starting values")
  if (! is.matrix(Z)) stop("wrong Z values")


  FUN = setup$posterior$density

  if(is.null(settings$parallel)) parallel = setup$parallel else parallel <- settings$parallel
  if(parallel == T & setup$parallel == F) stop("parallel = T requested in DEzs but BayesianSetup does not support parallelization. See help of BayesianSetup on how to enable parallelization")

  ## Initialize blockUpdate parameters and settings
  blockdefault <- list("none", k = NULL, h = NULL, pSel = NULL, pGroup = NULL,
                       groupStart = 1000, groupIntervall = 1000)

  if(!is.null(settings$blockUpdate)){
    blockUpdate <- utils::modifyList(blockdefault, settings$blockUpdate)
    blockUpdate[[1]] <- settings$blockUpdate[[1]] # to catch first argument
    if(blockUpdate[[1]] == "none"){
      blockUpdateType <- "none"
      blocks = FALSE
      BlockStart = FALSE
    }else{
      groupStart <- blockUpdate$groupStart
      groupIntervall <- blockUpdate$groupIntervall
      blockUpdateType = blockUpdate[[1]]
      blocks = TRUE
      ## Initialize BlockStart
      BlockStart = FALSE
      Bcount = 0
    }
  }else{
    blockUpdateType <- "none"
    blocks = FALSE
    BlockStart = FALSE
  }


  # Initialize parameter values. Because they are called in
  # the loop this saves time in comparison to referencing them
  # every iteration using settings$...
  iterations <- settings$iterations
  consoleUpdates <- settings$currentChain
  currentChain <- settings$currentChain
  pSnooker <- settings$pSnooker
  zUpdateFrequency <- settings$zUpdateFrequency
  pGamma1 <- settings$pGamma1
  eps.mult <- settings$eps.mult
  eps.add <- settings$eps.add

  # Initialization of previous chain length (= 0 if restart = F)
  lChainOld <- 0

  Npar <- ncol(X)
  Npar12 <- (Npar - 1)/2 # factor for Metropolis ratio DE Snooker update

  # M0 is initial population size of Z is the size of Z, it's the same number, only kept 2 to stay consistent with the ter Brakk & Vrugt 2008
  M = M0 = nrow(Z)
  Npop <- nrow(X)

  F2 = settings$f/sqrt(2*Npar)
  F1 = 1.0
  rr = NULL
  r_extra = 0

  #if(burnin != 0) stop("burnin option is currently not implemented")

  burnin <- settings$burnin/Npop
  n.iter <- ceiling(settings$iterations/Npop)
  if (n.iter < 2) stop ("The total number of iterations must be greater than 3")

  lChain <- ceiling((n.iter - burnin)/settings$thin)+1

  pChain <- array(NA, dim=c(lChain, Npar+3, Npop))

  colnames(pChain) <- c(setup$names, "LP", "LL", "LPr")


  # Print adjusted iterations
  #  cat("Iterations adjusted to", n.iter*Npop,"to fit settings", "\n")


  # assign memory for Z
  Zold <- Z
  Z <- matrix(NA, nrow= M0 + floor((n.iter-1) /zUpdateFrequency) * Npop, ncol=Npar)

  Z[1:M,] <- Zold


  counter <- 1
  counterZ <- 0

  # accept.prob <- 0
  logfitness_X <- FUN(X, returnAll = T)


  # Write first values in chain
  pChain[1,,] <- t(cbind(X,logfitness_X))



  for (iter in 2:n.iter) {
    f <- ifelse(iter%%10 == 0, 0.98, F1)
    #accept <- 0


    if(blocks){
      ### Update the groups.
      if(iter == groupStart+ Bcount*groupIntervall){
        blockSettings <- updateGroups(chain = pChain[1:counter,, ], blockUpdate)
        BlockStart <- TRUE
        Bcount <- Bcount + 1
      }
    }


    if(parallel == TRUE | parallel == "external"){
      x_prop <- matrix(NA, nrow= Npop, ncol=Npar)
      r_extra <- numeric(Npop)


      for(i in 1:Npop){
        # select to random different individuals (and different from i) in rr, a 2-vector
        rr <- sample.int(M, 3, replace = FALSE)
        if(stats::runif(1) < pSnooker) {
          z <- Z[rr[3],]
          x_z <- X[i,] - z
          D2 <- max(sum(x_z*x_z), 1.0e-300)
          projdiff <- sum((Z[rr[1],] -Z[rr[2],]) * x_z)/D2 # inner_product of difference with x_z / squared norm x_z
          gamma_snooker <- stats::runif(1, min=1.2,max=2.2)

          x_prop[i,] <-  X[i,] + gamma_snooker * projdiff * x_z
          x_z <- x_prop[i,] - z
          D2prop <- max(sum(x_z*x_z), 1.0e-300)
          r_extra[i] <- Npar12 * (log(D2prop) - log(D2))

        } else {
          if ( stats::runif(1)< pGamma1 ) { gamma_par = F1 # to be able to jump between modes
          } else {
            gamma_par = F2 * stats::runif(Npar, min=1-eps.mult, max=1+eps.mult)    # multiplicative error to be applied to the difference
            # gamma_par = F2
          }
          rr = sample.int(M, 2, replace = FALSE)
          if (eps.add ==0) {  # avoid generating normal random variates if possible
            x_prop[i,] = X[i,] + gamma_par * (Z[rr[1],]-Z[rr[2],])
          } else {
            x_prop[i,] = X[i,] + gamma_par * (Z[rr[1],]-Z[rr[2],])  +  eps.add*stats::rnorm(Npar,0,1)
          }
          r_extra = rep(0, Npop)
        }
      }
      # end proposal creation

      if(BlockStart){
        # Get the current group and update the proposal accordingly
        Member <- getBlock(blockSettings)
        x_prop[,-Member] <- X[,-Member]
        ####
      }


      # run proposals
      logfitness_x_prop <- FUN(x_prop, returnAll = T)

      # evaluate acceptance
      for(i in 1:Npop){
        if(!is.na(logfitness_x_prop[i,1] - logfitness_X[i,1])){
          if ((logfitness_x_prop[i,1] - logfitness_X[i,1] + r_extra[i]) > log(stats::runif(1))){
            # accept <- accept + 1
            X[i,] <- x_prop[i,]
            logfitness_X[i,] <- logfitness_x_prop[i,]
          }
        }
      }

    } else{
      # if not parallel

      for (i in 1:Npop){
        # select to random different individuals (and different from i) in rr, a 2-vector
        rr <- sample.int(M, 3, replace = FALSE)
        if(stats::runif(1) < pSnooker) {
          z <- Z[rr[3],]
          x_z <- X[i,] - z
          D2 <- max(sum(x_z*x_z), 1.0e-300)
          projdiff <- sum((Z[rr[1],] -Z[rr[2],]) * x_z)/D2 # inner_product of difference with x_z / squared norm x_z
          gamma_snooker <- stats::runif(1, min=1.2,max=2.2)
          x_prop <- X[i,] + gamma_snooker * projdiff * x_z
          x_z <- x_prop - z
          D2prop <- max(sum(x_z*x_z), 1.0e-300)
          r_extra <- Npar12 * (log(D2prop) - log(D2))
        } else {

          if ( stats::runif(1)< pGamma1 ) { gamma_par = F1 # to be able to jump between modes
          } else {
            gamma_par = F2 * stats::runif(Npar, min=1-eps.mult, max=1+eps.mult)    # multiplicative error to be applied to the difference
            # gamma_par = F2
          }
          rr = sample.int(M, 2, replace = FALSE)
          if (eps.add ==0) {  # avoid generating normal random variates if possible
            x_prop = X[i,] + gamma_par * (Z[rr[1],]-Z[rr[2],]) } else {
              x_prop = X[i,] + gamma_par * (Z[rr[1],]-Z[rr[2],])  +  eps.add*stats::rnorm(Npar,0,1)
            }
          r_extra = 0

        }
        if(BlockStart){
          # Get the current group and update the proposal accordingly
          Member <- getBlock(blockSettings)
          x_prop[-Member] <- X[i,-Member]
          ####
        }


        # evaluate proposal - can this be mixed with the parallel above?
        logfitness_x_prop <- FUN(x_prop, returnAll = T)

        # evaluate acceptance
        if(!is.na(logfitness_x_prop[1] - logfitness_X[i,1])){
          if ((logfitness_x_prop[1] - logfitness_X[i,1] + r_extra) > log(stats::runif(1))){
            # accept <- accept + 1
            X[i,] <- x_prop
            logfitness_X[i,] <- logfitness_x_prop
          }
        }
      } # for Npop


    }

    if ((iter > burnin) && (iter %% settings$thin == 0) ) { # retain sample
      counter <- counter+1
      pChain[counter,,] <- t(cbind(X,logfitness_X))

    }

    if (iter%%zUpdateFrequency == 0) { # update history

      Z[( M0 + (counterZ*Npop) + 1 ):( M0 + (counterZ+1)*Npop),] <- X
      counterZ <- counterZ +1
      M <- M + Npop
    }

    # Console update

    if(settings$message){
      if( (iter %% settings$consoleUpdates == 0) | (iter == n.iter)) cat("\r","Running DEzs-MCMC, chain ", currentChain,
                                                                         "iteration" ,iter*Npop,"of",n.iter*Npop,". Current logp ",
                                                                         logfitness_X[,1],". Please wait!","\r")
      utils::flush.console()
    }

    if(proc.time()[3]>save_inter[nbSave]){ # before we were saving at every iterations any(save_inter %in% iter)

      nbSave <- nbSave +1
      cat("\r", proc.time()[3])
      # cat("\r", iter)

      saveChain <- pChain

      saveChain <- lapply(seq(dim(saveChain)[3]), function(x) saveChain[ , , x])

      if(restart) saveChain <- lapply(seq_along(bayesianSetup$chain), function(x){rbind(as.array(bayesianSetup$chain[[x]]),pChain[,,x])})

      saveChain <- lapply(seq_along(saveChain), function(x) saveChain[[x]][rowSums(is.na(saveChain[[x]])) != ncol(saveChain[[x]]), ])
      saveZ <- Z[rowSums(is.na(Z)) != ncol(Z),]

      mcmcSampler = list(
        setup = setup,
        settings = settings,
        saveChain = saveChain,
        X = as.matrix(X[,1:Npar]),
        Z = saveZ,
        sampler = "DEzs"
      )

      class(mcmcSampler) <- "bayesianOutput"

      saveRDS(mcmcSampler, file = paste0("chainMW", index_saving, ".RDS"), compress = FALSE)
    }

  } # n.iter


  pChain <- pChain[1:counter,,]

  if(restart){ # Combine chains
    newchains <- array(NA, dim = c((counter+nrow(bayesianSetup$chain[[1]])), (Npar+3), Npop))

    for(i in 1:Npop){
      for(k in 1:(Npar+3)){
        newchains[,k,i] <- c(bayesianSetup$chain[[i]][,k],pChain[,k,i])
      }
    }
    pChain <- newchains
  }

  pChain<- coda::as.mcmc.list(lapply(1:Npop,function(i) coda::as.mcmc(pChain[,1:(Npar+3),i])))



  list(Draws = pChain,  X = as.matrix(X[,1:Npar]), Z = Z)

}




# getBlock function not exported in BayesianTools

getBlock <- function(blockSettings){
  groups <- blockSettings$cT
  pGroup <- blockSettings$pGroup
  pSel <- blockSettings$pSel


  nGroups = max(groups)
  if(nGroups == 1) return(1:length(groups))
  if (is.null(pGroup)) pGroup = rep(1,nGroups)
  if(length(pSel) > nGroups) pSel <- pSel[1:nGroups]
  pSel = c(pSel, rep(0,nGroups - length(pSel)))
  groupsToSample = sample.int(nGroups, 1, prob = pSel)

  selectedGroups = sample.int(nGroups,groupsToSample, prob = pGroup[1:nGroups])
  GroupMember  <- which(is.element(groups,selectedGroups))
  return(GroupMember)

}




# updateGroups function not exported in BayesianTools

updateGroups <- function(chain,blockSettings){

  settings <- getBlockSettings(blockSettings)
  blockUpdateType  <- settings$blockUpdateType

  switch(blockUpdateType,
         "correlation" = {
           ## (Pair wise) Correlation in the parameters
           cormat <- abs(stats::cor(chain[,1:(ncol(chain)-3),sample(1:dim(chain)[3],1)]))
           diag(cormat) <- 0
           # Correct for NA and Inf values as this could cause error in as.dist()
           cormat[c(which(is.na(cormat)),which(cormat == Inf),which(cormat == -Inf)) ] <- 0
           tree <- stats::hclust(stats::as.dist(1-cormat))  # get tree based on distance(dissimilarity = 1-cor).
           cT <- stats::cutree(tree, k = settings$k, h = settings$h) # get groups. With h we can manipulate the strength of the interaction.
         },
         "user" = {
           cT <-  settings$groups
         },
         "random" = {
           pool <- c(1:settings$k, sample(1:settings$k, (ncol(chain)-3-settings$k)))
           cT <- sample(pool)
         }
  )

  pSel <- settings$pSel
  if(is.null(pSel) && is.null(settings$pGroup)) pSel = rep(1,ncol(chain)-3)
  return(list(cT = cT, pGroup = settings$pGroup, pSel = pSel))
}




# getBlockSettings function not exported in BayesianTools

getBlockSettings <- function(blockUpdate){

  h <- k <- pSel <- pGroup <- groups <- NULL
  blockUpdateType <- blockUpdate[[1]]

  switch(blockUpdateType,
         "correlation" = {
           h <- blockUpdate$h
           k <- blockUpdate$k
           pSel <- blockUpdate$pSel
           pGroup <- blockUpdate$pGroup
         },
         "random"={
           k <- blockUpdate$k
         },
         "user"= {
           groups <- blockUpdate$groups
           pSel <- blockUpdate$pSel
           pGroup <- blockUpdate$pGroup
         })

  return(list(blockUpdateType = blockUpdateType, h = h, k = k, pSel = pSel,
              pGroup = pGroup, groups = groups))
}
