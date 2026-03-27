# this file contains all of the model code for each combination of available datasets

# load libraries
library(here)
library(nimble)

###### FULL IPM ######
IPMmod<-nimbleCode({

  # Model for the initial population size
  
  N.first ~ dcat(pNinit[1:lpNinit])
  
  N[1,1] ~ dbin(stable[1],N.first)
  N[2,1] <- N.first - N[1,1]
  
  # Process model over time
  for (t in 1:(n.occasions-1)){  
    N[1,t+1] ~ dpois(sj[t] * f[t] * (N[1,t] + N[2,t])) 
    N[2,t+1] ~ dbin(sa[t], (N[1,t] + N[2,t]))
  }
  
  #COUNT MODEL 
  for(t in 1:n.occasions){  
    Ntot[t] <- N[1,t] + N[2,t]
    for(i in 1:nsurvs){
      count[t,i] ~ dbin(p.count,Ntot[t])
    }
  }
  
  # p.count ~ dbeta(1,1)
  p.count ~ dunif(0,1)
  p.surv <- p.count #just rename to align with processing
  
  #REPRO MODEL  
  for(t in 1:n.occasions){  
    fledge[t] ~ dpois(rho[t])
    rho[t] <- broods[t]*f[t]
    f[t] <- fec     
  } 
  
  fec ~ dunif(0,5)
  
  #CMR MODEL 
  # priors for mark-recapture model 
  mean.phi[1] ~ dbeta(1,1) 
  mean.phi[2] ~ dbeta(1,1) 
  mean.p ~ dbeta(1,1)
  
  # constraints for mark-recapture model 
  for (t in 1:(n.occasions-1)){ 
    sj[t] <- mean.phi[1]
    sa[t] <- mean.phi[2]
    p[t] <- mean.p
  }
  
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,1:n.occasions], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,1:n.occasions], rel.a[t])
  }
  # Define the cell probabilities of the m-array
  for (t in 1:(n.occasions-1)){
    # Define probability of non-recapture (for compactness)
    q[t] <- 1-p[t]                
    
    # Main diagonal - captured at occasion after release
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    
    # Above main diagonal - captured again with some internal 0s
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    
    # Below main diagonal - not possible 
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  } #t  
  
  # DERIVED QUANTITIES #####

  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- (Ntot[t+1] + 1e-8) / (Ntot[t] + 1e-8) # adding tiny number to avoid Nan
  }
  # END derived quantities

})



##### NO NESTS #####
nonests <- nimbleCode({

  # Model for the initial population size
  
  N.first ~ dcat(pNinit[1:lpNinit])
  
  N[1,1] ~ dbin(stable[1],N.first)
  N[2,1] <- N.first - N[1,1]
  
  # Process model over time
  for (t in 1:(n.occasions-1)){  
    N[1,t+1] ~ dpois(sj[t] * f[t] * (N[1,t] + N[2,t])) 
    N[2,t+1] ~ dbin(sa[t], (N[1,t] + N[2,t]))
  }
  
  #COUNT MODEL 
  for(t in 1:n.occasions){  
    Ntot[t] <- N[1,t] + N[2,t]
    for(i in 1:nsurvs){
      count[t,i] ~ dbin(p.count,Ntot[t])
    }
  }
  
  p.count ~ dbeta(1,1)
  p.surv <- p.count #just rename to align with processing
  
  #REPRO MODEL  
  for(t in 1:n.occasions){  
    f[t] <- fec     
  } 
  
  fec ~ dunif(0,5)
  
  #CMR MODEL 
  # priors for mark-recapture model 
  mean.phi[1] ~ dbeta(1,1) 
  mean.phi[2] ~ dbeta(1,1) 
  mean.p ~ dbeta(1,1)
  
  # constraints for mark-recapture model 
  for (t in 1:(n.occasions-1)){ 
    sj[t] <- mean.phi[1]
    sa[t] <- mean.phi[2]
    p[t] <- mean.p
  }
  
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,1:n.occasions], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,1:n.occasions], rel.a[t])
  }
  # Define the cell probabilities of the m-array
  for (t in 1:(n.occasions-1)){
    # Define probability of non-recapture (for compactness)
    q[t] <- 1-p[t]                
    
    # Main diagonal - captured at occasion after release
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    
    # Above main diagonal - captured again with some internal 0s
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    
    # Below main diagonal - not possible 
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  } #t  
  
  # DERIVED QUANTITIES #####
  
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- (Ntot[t+1] + 1e-8) / (Ntot[t] + 1e-8) # adding tiny number to avoid Nan
  }
  # END derived quantities

})

##### NO MR #####

nomr<-nimbleCode({

  # Model for the initial population size
  
  N.first ~ dcat(pNinit[1:lpNinit])
  
  N[1,1] ~ dbin(stable[1],N.first)
  N[2,1] <- N.first - N[1,1]
  
  # Process model over time
  for (t in 1:(n.occasions-1)){  
    N[1,t+1] ~ dpois(sj[t] * f[t] * (N[1,t] + N[2,t])) 
    N[2,t+1] ~ dbin(sa[t], (N[1,t] + N[2,t]))
  }
  
  #COUNT MODEL 
  for(t in 1:n.occasions){  
    Ntot[t] <- N[1,t] + N[2,t]
    for(i in 1:nsurvs){
      count[t,i] ~ dbin(p.count,Ntot[t])
    }
  }
  
  p.count ~ dbeta(1,1)
  p.surv <- p.count #just rename to align with processing
  
  #REPRO MODEL  
  for(t in 1:n.occasions){  
    fledge[t] ~ dpois(rho[t])
    rho[t] <- broods[t]*f[t]
    f[t] <- fec     
  } 
  
  fec ~ dunif(0,5)
  
  #CMR MODEL 
  # priors for mark-recapture model 
  mean.phi[1] ~ dbeta(1,1) 
  mean.phi[2] ~ dbeta(1,1) 
  
  # DERIVED QUANTITIES #####
  
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- (Ntot[t+1] + 1e-8) / (Ntot[t] + 1e-8) # adding tiny number to avoid Nan
  }
  # END derived quantities

})

##### ABUND ONLY #####

abundonly<-nimbleCode({

  # Model for the initial population size
  
  N.first ~ dcat(pNinit[1:lpNinit])
  
  N[1,1] ~ dbin(stable[1],N.first)
  N[2,1] <- N.first - N[1,1]
  
  # Process model over time
  for (t in 1:(n.occasions-1)){  
    N[1,t+1] ~ dpois(sj[t] * f[t] * (N[1,t] + N[2,t])) 
    N[2,t+1] ~ dbin(sa[t], (N[1,t] + N[2,t]))
  }
  
  #COUNT MODEL 
  for(t in 1:n.occasions){  
    Ntot[t] <- N[1,t] + N[2,t]
    for(i in 1:nsurvs){
      count[t,i] ~ dbin(p.count,Ntot[t])
    }
  }
  
  p.count ~ dbeta(1,1)
  p.surv <- p.count #just rename to align with processing
  
  #REPRO MODEL  
  for(t in 1:n.occasions){  
    f[t] <- fec     
  } 
  
  fec ~ dunif(0,5)
  
  #CMR MODEL 
  # priors for mark-recapture model 
  mean.phi[1] ~ dbeta(1,1) 
  mean.phi[2] ~ dbeta(1,1) 
  
  # DERIVED QUANTITIES #####
  
  # Population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- (Ntot[t+1] + 1e-8) / (Ntot[t] + 1e-8) # adding tiny number to avoid Nan
  }
  # END derived quantities

})

