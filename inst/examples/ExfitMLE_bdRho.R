# Creating a phylogeny with 0.05 net diversification rate and 0.5 turnover rate.

set.seed(1234)
tree1 <- TESS::tess.sim.age(1, 100, lambda = 0.1, mu = 0.05, MRCA = TRUE, samplingProbability = 0.5)[[1]]
plot(tree1, root.edge = TRUE)

# Creating variables to give to arguments

tot_time <- max(TreeSim::getx(tree1))
Ntips <- ape::Ntip(tree1)
lamb_moments <- log(Ntips)/tot_time
ysim <- 0.5

# Fitting the birth-death∫rho model

res_fitMLE_M1 <- fitMLE_bdRho(phylo = tree1,
                              tot_time = tot_time, r = lamb_moments,
                              epsi = 0, y = NULL, ylamb = NULL,
                              reparam = FALSE, common = FALSE,
                              beta = FALSE, unif = TRUE,
                              a = 0, b = 1,
                              afix = TRUE, bfix =TRUE,
                              cond = "crown", YULE = FALSE,
                              dt = 0, rel.tol = 1e-10,
                              tuned_dichotomy = TRUE,
                              brk = 2000)

res_fitMLE_M1$parameters

# Fitting the classical birth-death-sampling model

res_fitMLE_M5 <- fitMLE_bdRho(phylo = tree1,
                              tot_time = tot_time, r = lamb_moments,
                              epsi = 0, y = ysim, ylamb = NULL,
                              reparam = FALSE, common = FALSE,
                              beta = FALSE, unif = FALSE,
                              a = NULL, b = NULL,
                              afix = NULL, bfix =NULL,
                              cond = "crown", YULE = FALSE,
                              dt = 0, rel.tol = 1e-10,
                              tuned_dichotomy = TRUE,
                              brk = 2000)

res_fitMLE_M5$parameters
