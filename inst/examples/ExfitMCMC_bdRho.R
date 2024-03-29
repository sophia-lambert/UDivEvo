# Creating a phylogeny with 0.05 net diversification rate and 0.5 turnover rate.

set.seed(1234)
tree1 <- TESS::tess.sim.age(1, 100, lambda = 0.1, mu = 0.05, MRCA = TRUE, samplingProbability = 0.5)[[1]]
plot(tree1, root.edge = TRUE)

# Creating variables to give to arguments

tot_time <- max(TreeSim::getx(tree1))
Ntips <- ape::Ntip(tree1)
lamb_moments <- log(Ntips)/tot_time
ysim <- 0.5

# Creating setting for MCMC

densityTest4 = function(x) {
  sum(dunif(x[1], min = 0, max = 1, log =TRUE)) + sum(dunif(x[2], min = 0, max = 1, log =TRUE))
}
samplerTest4 = function(n=1){
  s1 = runif(n, min = 0, max = 1)
  s2 = runif(n, min = 0, max = 1)
  return(cbind(s1,s2))
}
priorTest4 <- BayesianTools::createPrior(density = densityTest4, sampler = samplerTest4,
                                         lower = c(0,0), upper = c(1,1), best = NULL)
StartValueDTest4 = c(lamb_moments, runif(2, min = 0, max = 0.1))
StartValueEpsiTest4 = runif(3, min = 0, max = 1)
startValueTest4 = matrix(data = c(StartValueDTest4, StartValueEpsiTest4), nrow = 3, ncol = 2)

# Parameters for the checkpointing

nbIter <- 20000
maxTime <- 60*60*19.3 # 20 hours max (tiny less because of some processing issues)
stopTime <- 60*60*20
freqTime <- 60*60*3.21 # save every 3 hours
previousMCMC = NULL

# Fitting the birth-death∫rho model

res_fitMCMC_M1 <- fitMCMC_bdRho(phylo = tree1,
                                tot_time = tot_time, y = NULL,
                                reparam = FALSE, common = FALSE,
                                beta = FALSE, unif = TRUE,
                                a = 0, b = 1,
                                afix = TRUE, bfix =TRUE,
                                cond = "crown", YULE = FALSE,
                                dt = 0, rel.tol = 1e-10,
                                tuned_dichotomy = TRUE,
                                brk = 2000,
                                savedBayesianSetup = previousMCMC,
                                mcmcSettings = list(iterations = 3*nbIter,
                                                    startValue = startValueTest4),
                                prior = priorTest4,
                                parallel = FALSE, save_inter =
                                  c(seq(from = proc.time()[3], to = maxTime, by = freqTime),stopTime),
                                index_saving = as.factor("M1_tree1"))

plot(res_fitMCMC_M1$mcmc)

# Fitting the classical birth-death-sampling model

res_fitMCMC_M5 <- fitMCMC_bdRho(phylo = tree1,
                                tot_time = tot_time, y = ysim,
                                reparam = FALSE, common = FALSE,
                                beta = FALSE, unif = FALSE,
                                a = NULL, b = NULL,
                                afix = NULL, bfix =NULL,
                                cond = "crown", YULE = FALSE,
                                dt = 0, rel.tol = 1e-10,
                                tuned_dichotomy = TRUE,
                                brk = 2000,
                                savedBayesianSetup = previousMCMC,
                                mcmcSettings = list(iterations = 3*nbIter,
                                                    startValue = startValueTest4),
                                prior = priorTest4,
                                parallel = FALSE, save_inter =
                                  c(seq(from = proc.time()[3], to = maxTime, by = freqTime),stopTime),
                                index_saving = as.factor("M5_tree1"))
plot(res_fitMCMC_M5$mcmc)
