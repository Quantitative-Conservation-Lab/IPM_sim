# 11.3. Example of a simple IPM (counts, capture-recapture, reproduction)
# 11.3.1. Load data
# Population counts (from years 1 to n.years)

###### FULL IPM ######

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
  
  #population projection part
  # Initial population sizes
  n1.start ~ dunif(0, 1000)
  nad.start ~ dunif(0, 1000)
  N1[1] <- n1.start
  Nad[1] <- nad.start
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / Ntot[t]
  }
  #prior for survey detectin probability
  p.surv~dunif(0,1)
  #  productivity
  for (t in 1:(nyears-1)){
    f[t] <- fec
  }
  #mean.fec ~ dunif(0, 20)

  # 3.1. Likelihood for population population count data (state-space model)
  # 3.1.1 System process
  for (t in 2:nyears){
    #I think this was causing the problem, it was a little mis-specified
    #Double check that this is right, but it doesnt cause an error in
    #Slice samplers this way
    #mean1[t]<-((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1]
    #for reference, this is how it was:
    #mean1[t] <- f[t-1] * mean.phi[1] * Ntot[t-1]
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    #N1[t] ~ dpois(mean1[t])
    #Also changed that the total is the sum of N1 and Nad, it was just Ntot before
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }
  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }

  # 3.1.2 Observation process
  for(n in 1:n.sam){
    for (t in 1:nyears){
      y[n,t] ~ dbin(p.surv,Ntot[t])
    }
  }

  # 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
  for(i in 1:n.ind){
    z[i,first[i]]<-1
    for(t in (first[i]+1):(nyears)){
      z[i,t]~dbern(z[i,t-1]*mean.phi[age[i,t]]) #is this right? Because they survived as YOY
      ch.y[i,t]~dbern(z[i,t]*mean.p)
    }
  }
  #priors for resight and adult or yoy survival
  mean.phi[1]~dunif(0,1) #surv 1 year olds
  mean.phi[2]~dunif(0,1) #surv adults
  mean.p~dunif(0,1) #resight prob
  
  
  ##################################
  #################################
  #nest success model goes here
  
  # Priors and constraints #####
  lambdaf ~ dunif(0, 5)
  
  #changed this prior, since if it is too low the population goes negative
  #and I gather from looking that this may be reasonable?
  phi.nest ~ dunif(0.925, 1)
   
  for (n in 1:n.nests) {
    for(d in (first.nest[n] + 1):last.nest[n]) {
      H[n, d] ~ dbern(phi.nest * H[n, d - 1])
    } # d
  } # n
  
  for (c in 1:n.succ.nests) {
    # AEB note - should be using zero truncated poisson
    #Fledged[c] ~ T(dpois(lambdaf), 1, Inf)
    Fledged[c] ~ dpois(lambdaf)
  } # c

  # Derived quantities #####
  fec <- 1/2 * lambdaf * (phi.nest^max.nest.age)
  
  # END derived quantities
  
})

###### NO NEST ######

# create model versions that omit some of the data
noNests<-nimbleCode({
  #-------------------------------------------------
  #  Integrated population model
  #  - Age structured model with 2 age classes: 
  #		1-year old and adult (at least 2 years old)
  #  - Age at first breeding = 1 year
  #  - Prebreeding census, female-based
  #  - All vital rates assumed to be constant
  #-------------------------------------------------
  
  #population projection part
  # Initial population sizes
  n1.start ~ dunif(0, 1000)
  nad.start ~ dunif(0, 1000)
  N1[1] <- n1.start
  Nad[1] <- nad.start
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / Ntot[t]
  }
  #prior for survey detectin probability
  p.surv~dunif(0,1)
  #  productivity
  for (t in 1:(nyears-1)){
    f[t] <- fec
  }
  #mean.fec ~ dunif(0, 20)
  
  # 3.1. Likelihood for population population count data (state-space model)
  # 3.1.1 System process
  for (t in 2:nyears){
    #I think this was causing the problem, it was a little mis-specified
    #Double check that this is right, but it doesnt cause an error in
    #Slice samplers this way
    #mean1[t]<-((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1]
    #for reference, this is how it was:
    #mean1[t] <- f[t-1] * mean.phi[1] * Ntot[t-1]
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    #N1[t] ~ dpois(mean1[t])
    #Also changed that the total is the sum of N1 and Nad, it was just Ntot before
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }
  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }
  
  # 3.1.2 Observation process
  for(n in 1:n.sam){
    for (t in 1:nyears){
      y[n,t] ~ dbin(p.surv,Ntot[t])
    }
  }
  
  # 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
  for(i in 1:n.ind){
    z[i,first[i]]<-1
    for(t in (first[i]+1):(nyears)){
      z[i,t]~dbern(z[i,t-1]*mean.phi[age[i,t]]) #is this right? Because they survived as YOY
      ch.y[i,t]~dbern(z[i,t]*mean.p)
    }
  }
  #priors for resight and adult or yoy survival
  mean.phi[1]~dunif(0,1) #surv 1 year olds
  mean.phi[2]~dunif(0,1) #surv adults
  mean.p~dunif(0,1) #resight prob
  
  
  ##################################
  #################################
  #nest success model goes here
  
  # AEB note - this might not be a good way of doing it
  
  # Priors and constraints #####
  lambdaf ~ dunif(0, 5)
  
  #changed this prior, since if it is too low the population goes negative
  #and I gather from looking that this may be reasonable?
  phi.nest ~ dunif(0.925, 1)
  
  fec <- 1/2 * lambdaf * (phi.nest^max.nest.age)
  
})

###### NO MR ######

noMR<-nimbleCode({
  #-------------------------------------------------
  #  Integrated population model
  #  - Age structured model with 2 age classes: 
  #		1-year old and adult (at least 2 years old)
  #  - Age at first breeding = 1 year
  #  - Prebreeding census, female-based
  #  - All vital rates assumed to be constant
  #-------------------------------------------------
  
  #population projection part
  # Initial population sizes
  n1.start ~ dunif(0, 1000)
  nad.start ~ dunif(0, 1000)
  N1[1] <- n1.start
  Nad[1] <- nad.start
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / Ntot[t]
  }
  #prior for survey detection probability
  p.surv~dunif(0,1)
  #  productivity
  for (t in 1:(nyears-1)){
    f[t] <- fec
  }
  #mean.fec ~ dunif(0, 20)
  
  # 3.1. Likelihood for population population count data (state-space model)
  # 3.1.1 System process
  for (t in 2:nyears){
    #I think this was causing the problem, it was a little mis-specified
    #Double check that this is right, but it doesnt cause an error in
    #Slice samplers this way
    #mean1[t]<-((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1]
    #for reference, this is how it was:
    #mean1[t] <- f[t-1] * mean.phi[1] * Ntot[t-1]
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    #N1[t] ~ dpois(mean1[t])
    #Also changed that the total is the sum of N1 and Nad, it was just Ntot before
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }
  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }
  
  # 3.1.2 Observation process
  for(n in 1:n.sam){
    for (t in 1:nyears){
      y[n,t] ~ dbin(p.surv,Ntot[t])
    }
  }
  
  # MR stuff
  #priors for resight and adult or yoy survival
  mean.phi[1]~dunif(0,1) #surv 1 year olds
  mean.phi[2]~dunif(0,1) #surv adults
  
  ##################################
  #################################
  #nest success model goes here
  
  # Priors and constraints #####
  lambdaf ~ dunif(0, 5)
  
  #changed this prior, since if it is too low the population goes negative
  #and I gather from looking that this may be reasonable?
  phi.nest ~ dunif(0.925, 1)
  
  for (n in 1:n.nests) {
    for(d in (first.nest[n] + 1):last.nest[n]) {
      H[n, d] ~ dbern(phi.nest * H[n, d - 1])
    } # d
  } # n
  
  for (c in 1:n.succ.nests) {
    Fledged[c] ~ dpois(lambdaf) 
  } # c
  
  # Derived quantities #####
  fec <- 1/2 * lambdaf * (phi.nest^max.nest.age)
  
  # END derived quantities
  
})
