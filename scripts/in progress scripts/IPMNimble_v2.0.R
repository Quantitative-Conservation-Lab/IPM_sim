# TODO
# load data
# toggle between right model

library("nimble")

# initial value functions
source(here("scripts", "in progress scripts", "IPMinitvalues.R"))

###### FULL IPM ######
IPMmod<-nimbleCode({

  # COUNTS #####
  
  # System process
  
  # Initial population sizes
  n1.start ~ dunif(0, 500)
  nad.start ~ dunif(0, 500)
  N1[1] <- n1.start
  Nad[1] <- nad.start
  
  for (t in 2:nyears){
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }
  
  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }

  # Observation process
  
  #prior for survey detection probability
  p.surv~dunif(0,1)
  
  for(n in 1:n.sam){
    for (t in 1:nyears){
      y[n,t] ~ dbin(p.surv,Ntot[t])
    }
  }

  # CAPTURE RECAPTURE #####
  
  # System process
  for(i in 1:n.ind){
    z[i,first[i]]<-1
    for(t in (first[i]+1):(nyears)){
      # AEB note - changed this from age because the thing wasnt working
      # use this if we mark chicks at all
      #z[i,t]~dbern(z[i,t-1]*mean.phi[age[i, t-1] + 1]) # was mean.phi[age[i, t-1]]
      # use this if we mark adults or YOY only
      z[i,t]~dbern(z[i,t-1]*mean.phi[2])
      ch.y[i,t]~dbern(z[i,t]*mean.p)
    }
  }
  
  #priors for resight and adult or yoy survival
  mean.phi[1]~dunif(0,1) #surv 1 year olds
  mean.phi[2]~dunif(0,1) #surv adults
  mean.p~dunif(0,1) #resight prob
  
  # PRODUCTIVITY #####

  for (t in 1:(nyears-1)){
     OBS_nestlings[t] ~ dpois(rho[t])
     rho[t] <- R_obs[t]*f[t]
  }
  
  # priors
  for (t in 1:(nyears-1)){
    f[t] <- fec
  }
  
  fec ~ dunif(0, 5)
  
  # DERIVED QUANTITIES #####
  
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- (Ntot[t+1] + 1e-8) / (Ntot[t] + 1e-8) # adding tiny number to avoid Nan
  }
  # END derived quantities
  
})

##### NO NESTS ##### 

nonests<-nimbleCode({
  
  # COUNTS #####
  
  # System process
  
  # Initial population sizes
  n1.start ~ dunif(0, 1000)
  nad.start ~ dunif(0, 1000)
  N1[1] <- n1.start
  Nad[1] <- nad.start
  
  for (t in 2:nyears){
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }
  
  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }
  
  # Observation process
  
  #prior for survey detection probability
  p.surv~dunif(0,1)
  
  for(n in 1:n.sam){
    for (t in 1:nyears){
      y[n,t] ~ dbin(p.surv,Ntot[t])
    }
  }
  
  # CAPTURE RECAPTURE #####
  
  # System process
  for(i in 1:n.ind){
    z[i,first[i]]<-1
    for(t in (first[i]+1):(nyears)){
      z[i,t]~dbern(z[i,t-1]*mean.phi[age[i,t]]) 
      ch.y[i,t]~dbern(z[i,t]*mean.p)
    }
  }
  
  #priors for resight and adult or yoy survival
  mean.phi[1]~dunif(0,1) #surv 1 year olds
  mean.phi[2]~dunif(0,1) #surv adults
  mean.p~dunif(0,1) #resight prob
  
  # PRODUCTIVITY #####
  
  # priors
  for (t in 1:(nyears-1)){
    f[t] <- fec
  }
  fec ~ dunif(0, 5)
  
  # DERIVED QUANTITIES #####
  
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / Ntot[t]
  }
  # END derived quantities
  
})

##### NO MR ##### 

nomr<-nimbleCode({
  
  # COUNTS #####
  
  # System process
  
  # Initial population sizes
  n1.start ~ dunif(0, 1000)
  nad.start ~ dunif(0, 1000)
  N1[1] <- n1.start
  Nad[1] <- nad.start
  
  for (t in 2:nyears){
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }
  
  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }
  
  # Observation process
  
  #prior for survey detection probability
  p.surv~dunif(0,1)
  
  for(n in 1:n.sam){
    for (t in 1:nyears){
      y[n,t] ~ dbin(p.surv,Ntot[t])
    }
  }
  
  # CAPTURE RECAPTURE #####

  #priors for resight and adult or yoy survival
  mean.phi[1]~dunif(0,1) #surv 1 year olds
  mean.phi[2]~dunif(0,1) #surv adults
  
  # PRODUCTIVITY #####
  
  for (t in 1:(nyears-1)){
    OBS_nestlings[t] ~ dpois(rho[t])
    rho[t] <- R_obs[t]*f[t]
  }
  
  # priors
  for (t in 1:(nyears-1)){
    f[t] <- fec
  }
  fec ~ dunif(0, 5)
  
  # DERIVED QUANTITIES #####
  
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / Ntot[t]
  }
  # END derived quantities
  
})

##### ABUND ONLY ##### 

abundonly<-nimbleCode({
  
  # COUNTS #####
  
  # System process
  
  # Initial population sizes
  n1.start ~ dunif(0, 1000)
  nad.start ~ dunif(0, 1000)
  N1[1] <- n1.start
  Nad[1] <- nad.start
  
  for (t in 2:nyears){
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }
  
  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }
  
  # Observation process
  
  #prior for survey detection probability
  p.surv~dunif(0,1)
  
  for(n in 1:n.sam){
    for (t in 1:nyears){
      y[n,t] ~ dbin(p.surv,Ntot[t])
    }
  }
  
  # CAPTURE RECAPTURE #####
  
  #priors for resight and adult or yoy survival
  mean.phi[1]~dunif(0,1) #surv 1 year olds
  mean.phi[2]~dunif(0,1) #surv adults
  
  # PRODUCTIVITY #####
  
  # priors
  for (t in 1:(nyears-1)){
    f[t] <- fec
  }
  fec ~ dunif(0, 5)
  
  # DERIVED QUANTITIES #####
  
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / Ntot[t]
  }
  # END derived quantities
  
})

