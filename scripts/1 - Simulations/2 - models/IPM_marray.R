# this file contains all of the model code for each combination of available datasets

# load libraries
library("here")
library("nimble")

# source initial value functions
source(here("scripts", "current version", "2 - models", "IPMinitvalues.R"))

###### FULL IPM ######
IPMmod<-nimbleCode({

  # COUNTS #####

  # System process
  n1.start ~ dunif(0, 500)
  nad.start ~ dunif(0, 500)
  N1[1] <- round(n1.start)
  Nad[1] <- round(nad.start)

  for (t in 2:nyears){
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }

  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }

  # Observation process
  p.surv~dunif(0,1)
  for(n in 1:n.sam){
    for (t in 1:nyears){
      #y[n,t]~dnorm(Ntot[t],sd=p.surv)
      y[n,t] ~ dbin(p.surv,Ntot[t])
    }
  }

  # CAPTURE RECAPTURE #####
  
  #m-array, multinomial likelihood
  for(t in 1:(nyears-1)){
    marr[t,1:nyears]~dmulti(pr[t,1:nyears],R[t])
  }
  for(t in 1:(nyears-1)){
    q[t]<-1-p[t]
    pr[t,t]<-phi.ad[t]*p[t]
    for(j in (t+1):(nyears-1)){
      pr[t,j]<-prod(phi.ad[t:j])*prod(q[t:(j-1)])*p[j]
    }
    for(j in 1:(t-1)){
      pr[t,j]<-0
    }
  }
  for(t in 1:(nyears-1)){
    pr[t,nyears]<-1-sum(pr[t,1:(nyears-1)])
  }
  for(t in 1:(nyears)){
    phi.ad[t]<-mean.phi[2]
    p[t]<-mean.p
  }

  #priors for resight and adult or 1yo survival
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
  n1.start ~ dunif(0, 500)
  nad.start ~ dunif(0, 500)
  N1[1] <- round(n1.start)
  Nad[1] <- round(nad.start)

  for (t in 2:nyears){
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }

  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }

  # Observation process

  #prior for survey detection probability
  #p.surv ~ T(dnorm(0.5, sd = 0.3), 0, Inf)
  p.surv~dunif(0,1)
  for(n in 1:n.sam){
    for (t in 1:nyears){
      y[n,t] ~ dbin(p.surv,Ntot[t])
    }
  }

  # CAPTURE RECAPTURE #####
  #m-array, multinomial likelihood
  for(t in 1:(nyears-1)){
    marr[t,1:nyears]~dmulti(pr[t,1:nyears],R[t])
  }
  for(t in 1:(nyears-1)){
    q[t]<-1-p[t]
    pr[t,t]<-phi.ad[t]*p[t]
    for(j in (t+1):(nyears-1)){
      pr[t,j]<-prod(phi.ad[t:j])*prod(q[t:(j-1)])*p[j]
    }
    for(j in 1:(t-1)){
      pr[t,j]<-0
    }
  }
  for(t in 1:(nyears-1)){
    pr[t,nyears]<-1-sum(pr[t,1:(nyears-1)])
  }
  for(t in 1:(nyears)){
    phi.ad[t]<-mean.phi[2]
    p[t]<-mean.p
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
    lambda[t] <- (Ntot[t+1] + 1e-8) / (Ntot[t] + 1e-8)
  }
  # END derived quantities

})

##### NO MR #####

nomr<-nimbleCode({

  # COUNTS #####

  # System process

  # Initial population sizes
  n1.start ~ dunif(0, 500)
  nad.start ~ dunif(0, 500)
  N1[1] <- round(n1.start)
  Nad[1] <- round(nad.start)

  for (t in 2:nyears){
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }

  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }

  # Observation process

  #prior for survey detection probability
  #p.surv ~ T(dnorm(0.5, sd = 0.3), 0, Inf)
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
    lambda[t] <- (Ntot[t+1] + 1e-8) / (Ntot[t] + 1e-8)
  }
  # END derived quantities

})

##### ABUND ONLY #####

abundonly<-nimbleCode({

  # COUNTS #####

  # System process

  # Initial population sizes
  n1.start ~ dunif(0, 500)
  nad.start ~ dunif(0, 500)
  N1[1] <- round(n1.start)
  Nad[1] <- round(nad.start)

  for (t in 2:nyears){
    N1[t] ~ dpois(((f[t-1]*N1[t-1])+(f[t-1]*Nad[t-1]))*mean.phi[1])
    Nad[t] ~ dbin(mean.phi[2],(Ntot[t-1]))
  }

  for (t in 1:nyears){
    Ntot[t] <- round(Nad[t] + N1[t])
  }

  # Observation process

  #prior for survey detection probability
  #p.surv ~ T(dnorm(0.5, sd = 0.3), 0, Inf)
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
    lambda[t] <- (Ntot[t+1] + 1e-8) / (Ntot[t] + 1e-8)
  }
  # END derived quantities

})

