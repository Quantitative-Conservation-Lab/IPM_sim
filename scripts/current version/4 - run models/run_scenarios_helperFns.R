# this should have a function for
# running each model type

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
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

#### IPM ####


#### MCMC SETTINGS ####
nb <- 100000#0 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

# TODO
#change sutff

runIPMmod <- function(nb, ni, nt, nc,
                      popDat, popTraj,
                      comb, detect) {
  #### DATA ####
  dat1 <- list(y = popDat$SUR,
               marr = marray(popDat$ch),
               R = rowSums(marray(popDat$ch)),
               OBS_nestlings = popDat$OBS_nestlings,
               R_obs = popDat$R_obs
  )


  #### CONSTANTS ####

  const1 <- list(nyears = ncol(popDat$ch),
                 n.sam = nrow(popDat$SUR)#,
                 #first = popDat$firstobs
                 )

  #### INITIAL VALUES ####
  #z.state <- state.data(popDat$ch)

  inits1 <- list(
    mean.phi = c(comb$phi1, comb$phiad),#c(detect.h, detect.h),
    mean.p = detect[1],
    p.surv = detect[2],
    fec = comb$fec,#detect.h,
    #mean.phi = runif(2,0,1),#c(detect.h, detect.h),
    #mean.p = runif(1,0,1),#detect.h,
    #fec = runif(1,0,5),#detect.h,
    #z=z.state,
    n1.start=popTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
    nad.start=popTraj$Nouts[2,1]
  )

  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda") #,"Ntot","N1","Nad","f","rho")#0.3764911
  #dont need all of these later, but wanted to look at them now

  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = IPMmod, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  # lots of initial model checking you can do by exploring conf1
  # if you wanted to change samplers this is where you would do that
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)

  #### RUN MCMC ####
  outIPM <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
                    setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  return(outIPM)
}

#### NO NESTS ####

runnonests <- function(nb, ni, nt, nc,
                      popDat, popTraj,
                      comb, detect) {

  dat1 <- list(y = popDat$SUR,
               marr = marray(popDat$ch),
               R = rowSums(marray(popDat$ch))
  )


  #### CONSTANTS ####

  const1 <- list(nyears = ncol(popDat$ch),
                 n.sam = nrow(popDat$SUR)#,
                 #n.ind = nrow(popDat$ch),
                 #first = popDat$firstobs
                 )

  #### INITIAL VALUES ####
  #z.state <- state.data(popDat$ch)

  inits1 <- list(mean.phi = c(comb$phi1, comb$phiad),#c(detect.h, detect.h),
                 mean.p = detect[1],
                 p.surv = detect[2],
                 fec = comb$fec,#detect.h,
                 #mean.phi = runif(2,0,1),#c(detect.h, detect.h),
                 #mean.p = runif(1,0,1),#detect.h,
                 #fec = runif(1,0,5),#detect.h,
                 #z=z.state,
                 n1.start=popTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
                 nad.start=popTraj$Nouts[2,1]
  )

  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda") #,"Ntot","N1","Nad","f","rho")#0.3764911
  #dont need all of these later, but wanted to look at them now

  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = nonests, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  # lots of initial model checking you can do by exploring conf1
  # if you wanted to change samplers this is where you would do that
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)

  #### RUN MCMC ####
  #sink("sad_output.txt")
  #changed to checking to just see a matrix, since it is working!
  outnonests <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
                        setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  return(outnonests)

}

#### NO MR ####

runnomr <- function(nb, ni, nt, nc,
                      popDat, popTraj,
                      comb, detect) {

  #### DATA ####
  dat1 <- list(y = popDat$SUR,
               OBS_nestlings = popDat$OBS_nestlings,
               R_obs = popDat$R_obs
  )


  #### CONSTANTS ####

  const1 <- list(nyears = ncol(popDat$ch),
                 n.sam = nrow(popDat$SUR))

  #### INITIAL VALUES ####
  #z.state <- state.data(popDat$ch)

  inits1 <- list(mean.phi = c(comb$phi1, comb$phiad),#c(detect.h, detect.h),
                 #mean.p = detect.h,
                 p.surv = detect[2],
                 fec = comb$fec,#detect.h,
                 #mean.phi = runif(2,0,1),#c(detect.h, detect.h),
                 #mean.p = runif(1,0,1),#detect.h,
                 #fec = runif(1,0,5),#detect.h,
                 #z=z.state,
                 n1.start=popTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
                 nad.start=popTraj$Nouts[2,1]
  )

  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi", "fec", "lambda") #,"Ntot","N1","Nad","f","rho")#0.3764911
  #dont need all of these later, but wanted to look at them now

  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = nomr, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  # lots of initial model checking you can do by exploring conf1
  # if you wanted to change samplers this is where you would do that
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)

  #### RUN MCMC ####
  #sink("sad_output.txt")
  outnomr <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
                     setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  return(outnomr)

}

#### ABUND ONLY ####

runabundonly <- function(nb, ni, nt, nc,
                      popDat, popTraj,
                      comb, detect) {

  #### DATA ####
  dat1 <- list(y = popDat$SUR
  )


  #### CONSTANTS ####

  const1 <- list(nyears = ncol(popDat$ch),
                 n.sam = nrow(popDat$SUR))

  #### INITIAL VALUES ####
  #z.state <- state.data(popDat$ch)

  inits1 <- list(mean.phi = c(comb$phi1, comb$phiad),#c(detect.h, detect.h),
                 #mean.p = detect.h,
                 p.surv = detect[2],
                 fec = comb$fec,#detect.h,
                 #mean.phi = runif(2,0,1),#c(detect.h, detect.h),
                 #mean.p = runif(1,0,1),#detect.h,
                 #fec = runif(1,0,5),#detect.h,
                 #z=z.state,
                 n1.start=popTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
                 nad.start=popTraj$Nouts[2,1]
  )

  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi", "fec", "lambda") #,"Ntot","N1","Nad","f","rho")#0.3764911
  #dont need all of these later, but wanted to look at them now

  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = abundonly, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  # lots of initial model checking you can do by exploring conf1
  # if you wanted to change samplers this is where you would do that
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)

  #### RUN MCMC ####
  outabund <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
                      setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  return(outabund)

}
