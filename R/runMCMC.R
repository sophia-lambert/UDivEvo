# A modified version of the BayesianTools::runMCMC

runMCMC <- function(bayesianSetup , sampler = "DEzs", settings = NULL,  save_inter = NULL, index_saving = NULL){

  options(warn = 0)

  ptm <- proc.time()

  ####### RESTART ##########

  if("bayesianOutput" %in% class(bayesianSetup)){

    # TODO - the next statements should have assertions in case someone overwrites the existing setting or similar

    previousMcmcSampler <- bayesianSetup


    # Catch the settings in case of nrChains > 1
    if(!("mcmcSamplerList" %in% class(previousMcmcSampler) |  "smcSamplerList" %in% class(previousMcmcSampler) )){
      if(is.null(settings)) settings <- previousMcmcSampler$settings
      setup <- previousMcmcSampler$setup
      sampler <- previousMcmcSampler$settings$sampler
      previousSettings <- previousMcmcSampler$settings
    } else{
      if(is.null(settings)) settings <- previousMcmcSampler[[1]]$settings
      settings$nrChains <- length(previousMcmcSampler)
      setup <- previousMcmcSampler[[1]]$setup
      sampler <- previousMcmcSampler[[1]]$settings$sampler
      previousSettings <- previousMcmcSampler[[1]]$settings
    }

    # Set settings$sampler (only needed if new settings are supplied)
    settings$sampler <- sampler

    # overwrite new settings
    for(name in names(settings)) previousSettings[[name]] <- settings[[name]]

    settings <- previousSettings

    # Check if previous settings will be new default

    previousMcmcSampler$settings <- BayesianTools::applySettingsDefault(settings = settings, sampler = settings$sampler, check = TRUE)

    restart <- TRUE


    ## NOT RESTART STARTS HERE ###################

  }else if(class(bayesianSetup) == "BayesianSetup"){
    restart <- FALSE

    if(is.null(settings$parallel)) settings$parallel <- bayesianSetup$parallel
    if(is.numeric(settings$parallel)) settings$parallel <- TRUE

    setup <- BayesianTools::checkBayesianSetup(bayesianSetup, parallel = settings$parallel)
    settings <- BayesianTools::applySettingsDefault(settings = settings, sampler = sampler, check = TRUE)
  } else stop("runMCMC requires a class of type BayesianOutput or BayesianSetup")

  ###### END RESTART ##############


  # TODO - the following statement should be removed once all further functions access settings$sampler instead of sampler
  # At the moment only the same sampler can be used to restart sampling.
  sampler = settings$sampler

  #### Assertions
  if(!restart && setup$numPars == 1) if(!BayesianTools::getPossibleSamplerTypes()$univariate[which(BayesianTools::getPossibleSamplerTypes()$BTname == settings$sampler)]) stop("This sampler can not be applied to a univariate distribution")

  if(restart == T) if(!BayesianTools::getPossibleSamplerTypes()$restartable[which(BayesianTools::getPossibleSamplerTypes()$BTname == settings$sampler)]) stop("This sampler can not be restarted")

  ########### Recursive call in case multiple chains are to be run
  if(settings$nrChains >1){

    # Initialize output list
    out<- list()

    # Run several samplers
    for(i in 1:settings$nrChains){

      settingsTemp <- settings
      settingsTemp$nrChains <- 1 # avoid infinite loop
      settingsTemp$currentChain <- i

      if(restart){
        out[[i]] <- runMCMC(bayesianSetup = previousMcmcSampler[[i]], settings = settingsTemp) # is this a problem runMCMC ?
      }else{
        if(is.list(settings$startValue)) settingsTemp$startValue = settings$startValue[[i]]
        out[[i]] <- runMCMC(bayesianSetup = setup, sampler = settings$sampler, settings = settingsTemp)
      }
    }
    if(settings$sampler == "SMC") class(out) = c("smcSamplerList", "bayesianOutput")
    else class(out) = c("mcmcSamplerList", "bayesianOutput")
    return(out)

    ######### END RECURSIVE CALL
    # MAIN RUN FUNCTION HERE
  }else{

    # check start values
    setup$prior$checkStart(settings$startValue)


    if (sampler == "Metropolis" || sampler == "AM" || sampler == "DR" || sampler == "DRAM"){
      if(restart == FALSE){
        mcmcSampler <- BayesianTools::Metropolis(bayesianSetup = setup, settings = settings)
        mcmcSampler <- BayesianTools::sampleMetropolis(mcmcSampler = mcmcSampler, iterations = settings$iterations)
      } else {
        mcmcSampler <- BayesianTools::sampleMetropolis(mcmcSampler = previousMcmcSampler, iterations = settings$iterations)
      }
    }



    ############## Differential Evolution #####################
    if (sampler == "DE"){

      if(restart == F) out <- BayesianTools::DE(bayesianSetup = setup, settings = settings)
      else out <- BayesianTools::DE(bayesianSetup = previousMcmcSampler, settings = settings)

      #out <- DE(bayesianSetup = bayesianSetup, settings = list(startValue = NULL, iterations = settings$iterations, burnin = settings$burnin, eps = settings$eps, parallel = settings$parallel, consoleUpdates = settings$consoleUpdates,
      #            blockUpdate = settings$blockUpdate, currentChain = settings$currentChain))

      mcmcSampler = list(
        setup = setup,
        settings = settings,
        chain = out$Draws,
        codaChain = coda::mcmc(out$Draws),
        X = out$X,
        sampler = "DE"
      )
    }

    ############## Differential Evolution with snooker update
    if (sampler == "DEzs"){
      # check z matrix
      if(!is.null(settings$Z)) setup$prior$checkStart(settings$Z,z = TRUE)

      if(restart == F) out <- DEzs(bayesianSetup = setup, settings = settings, save_inter = save_inter, index_saving = index_saving)
      else out <- DEzs(bayesianSetup = previousMcmcSampler, settings = settings, save_inter = save_inter, index_saving = index_saving)

      mcmcSampler = list(
        setup = setup,
        settings = settings,
        chain = out$Draws,
        codaChain = coda::mcmc(out$Draws),
        X = out$X,
        Z = out$Z,
        sampler = "DEzs"
      )
    }

    ############## DREAM
    if (sampler == "DREAM"){

      if(restart == F) out <- BayesianTools::DREAM(bayesianSetup = setup, settings = settings)
      else out <- BayesianTools::DREAM(bayesianSetup = previousMcmcSampler, settings = settings)

      mcmcSampler = list(
        setup = setup,
        settings = settings,
        chain = out$chains,
        pCR = out$pCR,
        sampler = "DREAM",
        lCR = out$lCR,
        X = out$X,
        delta = out$delta
      )
    }

    ############## DREAMzs
    if (sampler == "DREAMzs"){
      # check z matrix
      if(!is.null(settings$Z)) setup$prior$checkStart(settings$Z,z = TRUE)

      if(restart == F) out <- BayesianTools::DREAMzs(bayesianSetup = setup, settings = settings)
      else out <- BayesianTools::DREAMzs(bayesianSetup = previousMcmcSampler, settings = settings)

      mcmcSampler = list(
        setup = setup,
        settings = settings,
        chain = out$chains,
        pCR = out$pCR,
        sampler = "DREAMzs",
        JumpRates = out$JumpRates,
        X = out$X,
        Z = out$Z
      )

    }

    if(sampler == "Twalk"){
      warning("At the moment using T-walk is discouraged: numeric instability")
      if(!restart){
        if(is.null(settings$startValue)){
          settings$startValue = bayesianSetup$prior$sampler(2)
        }
        mcmcSampler <- BayesianTools::Twalk(bayesianSetup = setup, settings = settings)
      }else{
        mcmcSampler <- BayesianTools::Twalk(bayesianSetup = previousMcmcSampler, settings = settings)
      }
      mcmcSampler$setup <- setup
      mcmcSampler$sampler <- "Twalk"
    }


    if ((sampler != "SMC")){
      class(mcmcSampler) <- c("mcmcSampler", "bayesianOutput")
    }

    ############# SMC #####################

    if (sampler == "SMC"){

      mcmcSampler <- BayesianTools::smcSampler(bayesianSetup = bayesianSetup, initialParticles = settings$initialParticles, iterations = settings$iterations, resampling = settings$resampling, resamplingSteps = settings$resamplingSteps, proposal = settings$proposal, adaptive = settings$adaptive, proposalScale = settings$proposalScale )
      mcmcSampler$settings = settings
    }

    mcmcSampler$settings$runtime = mcmcSampler$settings$runtime + proc.time() - ptm
    if(is.null(settings$message) || settings$message == TRUE){
      message("runMCMC terminated after ", mcmcSampler$settings$runtime[3], "seconds")
    }
    return(mcmcSampler)
  }
}
