# this file contains a function for running each model type

# TODO - any updates needed here? HAS tested already
# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,(n.occasions+1)] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

#### IPM ####

runIPMmod <- function(nb, ni, nt, nc,
                      #popDat, #popTraj,
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
                 n.sam = n.sam)

  #### INITIAL VALUES ####
  inits1 <- list(
    mean.phi = c(comb$phi1, comb$phiad),
    mean.p = det.MR,
    p.surv = det.abund,
    fec = comb$fec,
    n1.start = pop1$N[1,1],
    nad.start = pop1$N[2,1]
  )

  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda","Ntot")#,"N1","Nad","f","rho")#0.3764911

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
                      #popDat, #popTraj,
                      comb, detect) {

  dat1 <- list(y = surv_cnts,
               marr.a = marr.a,
               marr.j = marr.j,
               R.j = rowSums(marr.j), 
               R.a = rowSums(marr.a))


  #### CONSTANTS ####

  const1 <- list(nyears = nyears,
                 n.sam = n.sam)

  #### INITIAL VALUES ####
  #z.state <- state.data(popDat$ch)

  inits1 <- list(
    mean.phi = c(comb$phi1, comb$phiad),
    mean.p = det.MR,
    p.surv = det.abund,
    fec = comb$fec,
    n1.start = pop1$N[1,1],
    nad.start = pop1$N[2,1]
  )

  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda","Ntot")#,"N1","Nad","f","rho")#0.3764911

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
                      #popDat, #popTraj,
                      comb, detect) {

  #### DATA ####
  dat1 <- list(y = surv_cnts,
               OBS_nestlings = obs_nestlings,
               R_obs = obs_nests)


  #### CONSTANTS ####

  const1 <- list(nyears = nyears,
                 n.sam = n.sam)

  #### INITIAL VALUES ####
  #z.state <- state.data(popDat$ch)

  inits1 <- list(
    mean.phi = c(comb$phi1, comb$phiad),
    mean.p = det.MR,
    p.surv = det.abund,
    fec = comb$fec,
    n1.start = pop1$N[1,1],
    nad.start = pop1$N[2,1]
  )

  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi", "fec", "lambda","Ntot")#,"N1","Nad","f","rho")#0.3764911

  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = nomr, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)

  #### RUN MCMC ####
  #sink("sad_output.txt")
  outnomr <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,thin=nt,
                     setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  return(outnomr)

}

#### ABUND ONLY ####

runabundonly <- function(nb, ni, nt, nc,
                      #popDat, #popTraj,
                      comb, detect) {

  #### DATA ####
  dat1 <- list(y = surv_cnts)


  #### CONSTANTS ####

  const1 <- list(nyears = nyears,
                 n.sam = n.sam)

  #### INITIAL VALUES ####
  #z.state <- state.data(popDat$ch)

  inits1 <- list(
    mean.phi = c(comb$phi1, comb$phiad),
    mean.p = det.MR,
    p.surv = det.abund,
    fec = comb$fec,
    n1.start = pop1$N[1,1],
    nad.start = pop1$N[2,1]
  )

  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi", "fec", "lambda","Ntot")#,"N1","Nad","f","rho")#0.3764911

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
