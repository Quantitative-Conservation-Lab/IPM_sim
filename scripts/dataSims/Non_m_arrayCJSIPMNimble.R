# TODO
# abby put nest stuff in here

#IPM model for sim data
library(here)
#source(here("scripts", "dataSims","BPA_IPMsimData.R")) #full on BPA WinBUGS code
source(here("scripts", "dataSims","IPMsimData_nonest.R")) #changed the MR data

# 11.3. Example of a simple IPM (counts, capture-recapture, reproduction)
# 11.3.1. Load data
# Population counts (from years 1 to n.years)
y <- SUR#df$SUR

# Capture-recapture data (in m-array format, from years 1 to n.years)
m <- ch#df$ch


#Productivity/nest success not in the data sim right now
#these will need to be changed once we have the nest success data
# Productivity data (from years 1 to n.years-1)
J<-nestlings#df$nestlings
R<-R#df$R


#Use the IPM from WinBugs, but put in nimble
# 11.3.2. Analysis of the model
# Specify model in BUGS language
library("nimble")
IPMmod<-nimbleCode({
  #-------------------------------------------------
  #  Integrated population model
  #  - Age structured model with 2 age classes: 
  #		1-year old and adult (at least 2 years old)
  #  - Age at first breeding = 1 year
  #  - Prebreeding census, female-based
  #  - All vital rates assumed to be constant
  #-------------------------------------------------
  
  #-------------------------------------------------
  # 1. Define the priors for the parameters
  #-------------------------------------------------
  # Observation error
  tauy <- pow(sigma.y, -2)
  sigma.y ~ dunif(0, 50)
  sigma2.y <- pow(sigma.y, 2)
  
  # Initial population sizes
  # Note that this part is different than in BUGS:
  # 1. JAGS seems to be very sensitive to the choice of the prior distribution for the initial population sizes and of the choice of the observation model (priors)
  # 2. Since the initial population sizes are used in binomial distributions, the numbers must be integers, otherwise JAGS does not run.
  # The following specification seems to work very well and produces similar results as BUGS:
  
  n1 ~ T(dnorm(25, tauy),0,)     # 1-year
  nad ~ T(dnorm(25, tauy),0,)    # Adults
  N1[1] <- round(n1)
  Nad[1] <- round(nad)
  
  # Survival and recapture probabilities, as well as productivity
  for (t in 1:(nyears-1)){
    f[t] <- mean.fec
  }
  

  mean.fec ~ dunif(0, 20)
  
  #-------------------------------------------------
  # 2. Derived parameters
  #-------------------------------------------------
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / Ntot[t]
  }
  
  #-------------------------------------------------
  # 3. The likelihoods of the single data sets
  #-------------------------------------------------
  # 3.1. Likelihood for population population count data (state-space model)
  # 3.1.1 System process
  for (t in 2:nyears){
    mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
    N1[t] ~ dpois(mean1[t])
    Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
  }
  for (t in 1:nyears){
    Ntot[t] <- Nad[t] + N1[t]
  }
  
  # 3.1.2 Observation process
  for (t in 1:nyears){
    y[t] ~ dnorm(Ntot[t], tauy)
  }
  
  # 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
  for(i in 1:n.ind){
    for(t in (first[i]+1):(nyears)){
      z[i,t]~dbern(z[i,t-1]*mean.phi)
      ch.y[i,t]~dbern(z[i,t]*mean.p)
    }
  }
  #priors for CJS
  # for(i in 1:n.ind){
  #   for(t in 1:(nyears-1)){
  #     phi[i,t]<-mean.phi[i,t]
  #     p[i,t]<-mean.p
  #   }
  # }
  mean.phi~dunif(0,1)
 #mean.phi[1]~dunif(0,1) #surv 1 year olds
# mean.phi[2]~dunif(0,1) #surv adults
 mean.p~dunif(0,1) #resight prob
 #mean.m~dunif(0,1) #mark probability??????
 
  
  # 3.3. Likelihood for productivity data: Poisson regression
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*f[t]
  }
})


# data

datipm <- list(ch.y = m, y = y, J = J, R = R)
constants<-list(nyears = length(y), 
                n.ind=nrow(m), first=first)
                
#Provide known state as data  
state.data <- function(EH){
  state <- EH
  for (i in 1:dim(EH)[1]){
    first <- min(which(EH[i,]==1))
    last <- max(which(EH[i,]==1))
    state[i,first:last] <- 1
    #   state[i,first] <- NA
  }
  state[state==0] <- NA
  return(state)
}

z.state <- state.data(m)
# Initial values
inits <- list(mean.phi=runif(1,0,1),
              mean.p = runif(1, 0, 1), 
              mean.fec = runif(1, 0, 10), 
              sigma.y = runif(1, 0, 1), 
              n1 = rpois(1, 100), 
              nad = rpois(1, 100),
              z=z.state)
parameters <- c("mean.phi", 
                "mean.p", "mean.fec", "N1", "Nad", 
                "Ntot", "sigma2.y", "lambda")

samples <- nimbleMCMC(
  code = IPMmod,
  constants = constants,#list(nyears = dim(m)[2], r = rowSums(m)),
  data =  datipm,#list(m = m, y = y, J = J, R = R), 
  summary=TRUE,
  inits = inits,
  monitors = parameters,
  niter = 1000, nburnin = 500, thin = 1)
