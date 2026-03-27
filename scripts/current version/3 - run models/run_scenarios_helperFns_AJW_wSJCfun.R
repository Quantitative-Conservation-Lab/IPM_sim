# this file contains a function for running each model type


#### IPM ####

runIPMmod <- function(nb, ni, nt, nc,
                      detect,
                      comb) {
  #### DATA ####
  dat1 <- list(marr.j = marray[,,1], marr.a = marray[,,2], 
               fledge = obs_nestlings, broods = obs_broods, count = count)

  #### CONSTANTS ####
  pNinit <- dUnif(10,800)
  lpNinit <- length(pNinit)
  const1 <- list(n.occasions = nyears, rel.j = rowSums(marray[,,1]), 
                 rel.a = rowSums(marray[,,2]), nsurvs = dim(count)[2], 
                 pNinit = pNinit, stable = stable, pNinit = pNinit,
                 lpNinit = lpNinit)

  #### INITIAL VALUES ####
  N.st <- matrix(NA,nrow = 2,ncol = nyears)
  N.st[1,] <- apply(count,1,max)+30
  N.st[2,] <- apply(count,1,max)+30

  #Initial values 
  inits1 <- list(mean.phi = phi, #mean.phi=runif(2),
                 fec = fec,
                 mean.p=runif(1),
                 p.count = det.abund,
                 N = N.st)

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
  outIPM <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1, thin = nt,
                    setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  return(outIPM)
}

#### NO NESTS ####

runnonests <- function(nb, ni, nt, nc,
                      comb, detect) {

  #### DATA ####
  dat1 <- list(marr.j = marray[,,1], marr.a = marray[,,2], 
               fledge = obs_nestlings, broods = obs_broods, count = count)
  
  #### CONSTANTS ####
  pNinit <- dUnif(10,800)
  lpNinit <- length(pNinit)
  const1 <- list(n.occasions = nyears, rel.j = rowSums(marray[,,1]), 
                 rel.a = rowSums(marray[,,2]), nsurvs = dim(count)[2], 
                 pNinit = pNinit, stable = stable, pNinit = pNinit,
                 lpNinit = lpNinit)
  
  #### INITIAL VALUES ####
  N.st <- matrix(NA,nrow = 2,ncol = nyears)
  N.st[1,] <- apply(count,1,max)+30
  N.st[2,] <- apply(count,1,max)+30

  inits1 <- list(mean.phi = phi, #mean.phi = runif(2),
                 fec = fec,
                 mean.p = runif(1),
                 p.count = det.abund,
                 N = N.st)
  
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
  outnonests <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1, thin = nt,
                        setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  return(outnonests)

}

#### NO MR ####

runnomr <- function(nb, ni, nt, nc,
                      comb, detect) {

  #### DATA ####
  dat1 <- list(marr.j = marray[,,1], marr.a = marray[,,2], 
               fledge = obs_nestlings, broods = obs_broods, count = count)
  
  #### CONSTANTS ####
  pNinit <- dUnif(10,800)
  lpNinit <- length(pNinit)
  const1 <- list(n.occasions = nyears, rel.j = rowSums(marray[,,1]), 
                 rel.a = rowSums(marray[,,2]), nsurvs = dim(count)[2], 
                 pNinit = pNinit, stable = stable, pNinit = pNinit,
                 lpNinit = lpNinit)
  
  #### INITIAL VALUES ####
  N.st <- matrix(NA,nrow = 2,ncol = nyears)
  N.st[1,] <- apply(count,1,max)+30
  N.st[2,] <- apply(count,1,max)+30

  inits1 <- list(mean.phi = phi, #mean.phi = runif(2),
                 fec = fec,
                 mean.p = runif(1),
                 p.count = det.abund,
                 N = N.st)
  
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
  outnomr <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1, thin=nt,
                     setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  return(outnomr)

}

#### ABUND ONLY ####

runabundonly <- function(nb, ni, nt, nc,
                      comb, detect) {

  #### DATA ####
  dat1 <- list(marr.j = marray[,,1], marr.a = marray[,,2], 
               fledge = obs_nestlings, broods = obs_broods, count = count)
  
  #### CONSTANTS ####
  pNinit <- dUnif(10,800)
  lpNinit <- length(pNinit)
  const1 <- list(n.occasions = nyears, rel.j = rowSums(marray[,,1]), 
                 rel.a = rowSums(marray[,,2]), nsurvs = dim(count)[2], 
                 pNinit = pNinit, stable = stable, pNinit = pNinit,
                 lpNinit = lpNinit)
  
  #### INITIAL VALUES ####
  N.st <- matrix(NA,nrow = 2,ncol = nyears)
  N.st[1,] <- apply(count,1,max)+30
  N.st[2,] <- apply(count,1,max)+30

  inits1 <- list(mean.phi = phi, #mean.phi = runif(2),
                 fec = fec,
                 mean.p = runif(1),
                 p.count = det.abund,
                 N = N.st)
  
  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi", "mean.p", "fec", "lambda","Ntot")
  
  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = IPMmod, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)

  #### RUN MCMC ####
  outabund <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1, thin = nt,
                      setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  return(outabund)

}
