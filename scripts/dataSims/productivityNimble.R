# NEST SURV NIMBLE TEST

# we need to make sure we can recover simulated parameters 
# with a simple nest survival model

# load libraries #####
library(nimble)
library(here)

# source files #####
source(here("scripts", "dataSims", "productivityDataSim.R"))

# Nimble model #####
code1 <- nimbleCode({
  
  # Priors and constraints #####
  lambdaf ~ dunif(0, 10)
  
  phi ~ dunif(0, 1)
  
  # END priors and constraints
  
  # Likelihood #####
  
  for (n in 1:n.nests) {
    for(d in (first[n] + 1):last[n]) {
      H[n, d] ~ dbern(phi * H[n, d - 1])
    } # d
  } # n
  
  for (c in 1:n.succ.nests) {
    Fledged[c] ~ dpois(lambdaf) 
  } # c
  
  # END likelihood
  
  # Derived quantities #####

  fec <- 1/2 * lambdaf * phi^max.nest.age
  
  # END derived quantities
})

# DATA ####
dat1 <- list(
  H = observed.nest.status, 
  Fledged = clutch.sizes
) # Stage history for each nest

# CONSTANTS ####
const1 <- list(
  n.nests = N.nests.found, 
  n.succ.nests = N.nests.successful, 
  first = first, 
  last = last, 
  max.nest.age = max.nest.age
)

# INITIAL VALUES ####

# function for H inits - set every NA to 1
Hinits <- observed.nest.status
Hinits[!is.na(observed.nest.status)] <- NA
Hinits[is.na(observed.nest.status)] <- 1


inits1 <- list(
  phi = runif(1, 0, 1), 
  lambdaf = runif(1, 0, 10),
  H = Hinits
)

# PARAMETERS ####
params1 <- c(
  "phi", 
  "lambdaf"
)

# MCMC SETTINGS ####
nb <- 5000 #burn-in
ni <- 10000 + nb #total iterations
nt <- 1  #thin
nc <- 3  #chains

# COMPILE CONFIGURE AND BUILD ####
Rmodel1 <- nimbleModel(code = code1, constants = const1, data = dat1, 
                       check = FALSE, calculate = FALSE, inits = inits1)
conf1 <- configureMCMC(Rmodel1, monitors = params1, thin = nt) # takes foreverrrrr 
#conf1$printSamplers(1:500) # check sampler defaults
Rmcmc1 <- buildMCMC(conf1)  
Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE) 
Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)

# RUN MCMC ####
t.start <- Sys.time()
out1 <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
                setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
t.end <- Sys.time()
(runTime <- t.end - t.start)
beep(sound = 1)

summary(out1)
library(coda)
plot(out1)
gelman.diag(out1)

