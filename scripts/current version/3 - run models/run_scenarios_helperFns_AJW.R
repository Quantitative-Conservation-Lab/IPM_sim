# this file contains a function for running each model type


#### IPM ####

runIPMmod <- function(nb, ni, nt, nc,
                      detect,
                      comb) {
  #### DATA ####
  dat1 <- list(y = surv_cnts,
               marr.a = marr.a,
               marr.j = marr.j,
               R.j = rowSums(marr.j), 
               R.a = rowSums(marr.a),
               OBS_nestlings = obs_nestlings,
               R_obs = obs_nests)
  
  #### CONSTANTS ####
  const1 <- list(nyears = nyears,
                 maxcount = maxcount,
                 stable = stable,
                 n.sam = n.sam)
  
  #### INITIAL VALUES ####
  inits1 <- list(
    mean.phi = c(comb$phi1, comb$phiad),
    mean.p = det.MR,
    p.surv = det.abund,
    fec = comb$fec,
    n1.start = pop1$N[1,1]+pop1$N[1,1]*0.25,
    nad.start = pop1$N[2,1]+pop1$N[2,1]*0.25,
    N1 = as.numeric(round(pop1$N[1,]+pop1$N[1,]*0.25)),
    Nad = as.numeric(round(pop1$N[2,]+pop1$N[2,]*0.25))
  )
  
  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda","Ntot")
  
  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = IPMmod, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
  
  #### RUN MCMC ####
  outIPM <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1, thin=nt,
                    setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  
  return(outIPM)
}

#### NO NESTS ####

runnonests <- function(nb, ni, nt, nc,
                       comb, detect) {
  
  dat1 <- list(y = surv_cnts,
               marr.a = marr.a,
               marr.j = marr.j,
               R.j = rowSums(marr.j), 
               R.a = rowSums(marr.a))
  
  
  #### CONSTANTS ####
  const1 <- list(nyears = nyears,
                 maxcount = maxcount,
                 stable = stable,
                 n.sam = n.sam)
  
  #### INITIAL VALUES ####
  inits1 <- list(
    mean.phi = c(comb$phi1, comb$phiad),
    mean.p = det.MR,
    p.surv = det.abund,
    fec = comb$fec,
    n1.start = pop1$N[1,1]+pop1$N[1,1]*0.25,
    nad.start = pop1$N[2,1]+pop1$N[2,1]*0.25,
    N1 = as.numeric(round(pop1$N[1,]+pop1$N[1,]*0.25)),
    Nad = as.numeric(round(pop1$N[2,]+pop1$N[2,]*0.25))
  )
  
  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda","Ntot")
  
  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = nonests, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
  
  #### RUN MCMC ####
  outnonests <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,thin=nt,
                        setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  
  return(outnonests)
  
}

#### NO MR ####

runnomr <- function(nb, ni, nt, nc,
                    comb, detect) {
  
  #### DATA ####
  dat1 <- list(y = surv_cnts,
               OBS_nestlings = obs_nestlings,
               R_obs = obs_nests)
  
  
  #### CONSTANTS ####
  const1 <- list(nyears = nyears,
                 maxcount = maxcount,
                 stable = stable,
                 n.sam = n.sam)
  
  #### INITIAL VALUES ####
  inits1 <- list(
    mean.phi = c(comb$phi1, comb$phiad),
    p.surv = det.abund,
    fec = comb$fec,
    n1.start = pop1$N[1,1]+pop1$N[1,1]*0.25,
    nad.start = pop1$N[2,1]+pop1$N[2,1]*0.25,
    N1 = as.numeric(round(pop1$N[1,]+pop1$N[1,]*0.25)),
    Nad = as.numeric(round(pop1$N[2,]+pop1$N[2,]*0.25))
  )
  
  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi", "fec", "lambda","Ntot")
  
  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = nomr, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
  
  #### RUN MCMC ####
  outnomr <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,thin=nt,
                     setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  
  return(outnomr)
  
}

#### ABUND ONLY ####

runabundonly <- function(nb, ni, nt, nc,
                         comb, detect) {
  
  #### DATA ####
  dat1 <- list(y = surv_cnts)
  
  
  #### CONSTANTS ####
  const1 <- list(nyears = nyears,
                 maxcount = maxcount,
                 stable = stable,
                 n.sam = n.sam)
  
  #### INITIAL VALUES ####
  inits1 <- list(
    mean.phi = c(comb$phi1, comb$phiad),
    mean.p = det.MR,
    p.surv = det.abund,
    fec = comb$fec,
    n1.start = pop1$N[1,1]+pop1$N[1,1]*0.25,
    nad.start = pop1$N[2,1]+pop1$N[2,1]*0.25,
    N1 = as.numeric(round(pop1$N[1,]+pop1$N[1,]*0.25)),
    Nad = as.numeric(round(pop1$N[2,]+pop1$N[2,]*0.25))
  )
  
  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", 
               "mean.phi", "fec",
               "lambda","Ntot")
  
  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = abundonly, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
  
  #### RUN MCMC ####
  outabund <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,thin=nt,
                      setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  
  return(outabund)
  
}
