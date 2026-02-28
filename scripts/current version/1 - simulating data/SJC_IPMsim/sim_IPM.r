####Version 2/7/2018

library(here)
library(nimble)
library(MCMCvis)
library(coda)

rm(list=ls())

source(here("scripts/current version/1 - simulating data/SJC_IPMsim","create.pop.r"))
source(here("scripts/current version/1 - simulating data/SJC_IPMsim","create.eh.r"))
source(here("scripts/current version/1 - simulating data/SJC_IPMsim","create.survey.r"))
source(here("scripts/current version/1 - simulating data/SJC_IPMsim","create.repro.r"))

######################################################################
#                                                                    #  
#                        Simulation Parameters                       #
#                                                                    #    
######################################################################
sims <- 2

# Number of years
yrs <- 10 

#Proportion of breeding population in sample
pbrood <- 0.05

#Probability of detection for abundance surveys
pobs <- 0.2

######################################################################
#                                                                    #  
#                        Create Population                           #
#                                                                    #    
######################################################################

# Record time 
start.time <- Sys.time()

all.results <- array(NA,dim = c(sims,25,5))
dimnames(all.results)[[2]] <- c("N[1, 1]","N[2, 1]","N[1, 2]","N[2, 2]","N[1, 3]","N[2, 3]",
                                "N[1, 4]","N[2, 4]","N[1, 5]","N[2, 5]","N[1, 6]","N[2, 6]",    
                                "N[1, 7]","N[2, 7]","N[1, 8]","N[2, 8]","N[1, 9]","N[2, 9]",  
                                "N[1, 10]","N[2, 10]","mean.f","mean.p","mean.phi[1]","mean.phi[2]",
                                "p.count")
dimnames(all.results)[[3]] <- c("Mean","Median","St.Dev.","95%CI_low","95%CI_upp")

# Start simulations 
for(s in 1:sims){

  # Age specific survival probabilities (juv, adult)
  sj <- 0.33
  sa <- 0.52

  # Fecundity parameters 
  f <- 1 
  
  # Survival and fecundity   
  surv.mat <- matrix(NA,nrow=2,ncol=yrs)
  surv.mat[1,1:yrs] <- rep(sj,yrs)
  surv.mat[2,1:yrs] <- rep(sa,yrs)
  surv.mat <- surv.mat[,-yrs]
  
  fec.mat <- matrix(f,nrow=nrow(surv.mat),ncol=yrs)
  
  # Initial population size per age class 
  Ni <- c(10000, 10000)
      
  # Create the true population for this simulation and year  
  ind <- create.population(phi = surv.mat, f = fec.mat, Im = rep(0, yrs), Ni = Ni)
      
  ######################################################################
  #                                                                    #  
  #                        Create Survey Data                          #
  #                                                                    #    
  ######################################################################
      
  # Detection probability for the population survey
  pobs <- pobs
  nsurvs <- 2 #number of surveys per year 
  
  # Create the population survey data
  count <- create.survey.bin(ind$Nu["Total",], rep(pobs, yrs), nsurvs = nsurvs)
  
  ######################################################################
  #                                                                    #  
  #                        Create Capture Data                         #
  #                                                                    #    
  ######################################################################
      
  # Capture and recapture probabilities
  cjuv <- 0.3           # initial capture probability of juveniles
  cad <- 0.3            # initial capture probability of adults
  prec <- 0.6           # recapture probability 
      
  # Create the capture histories and the corresponding m-arrays
    
  ch <- create.capturehistory(ind$IND, c = matrix(c(rep(cjuv, yrs), rep(cad, yrs)), nrow = 2, byrow = TRUE), p = matrix(c(rep(prec, yrs-1), rep(prec, yrs-1)), nrow = 2, byrow = TRUE))
  
  ch.a <- ch$ch[,1:yrs]
  incl <- which(rowSums(ch.a)>=1)
  ch.a <- ch.a[incl,]
  ch.b <- ch$age[incl]
      
  marray <- marray.age(ch.a, ch.b)
      
  ######################################################################
  #                                                                    #  
  #                        Create Repro Data                           #
  #                                                                    #    
  ######################################################################
      
  # Probability to find a brood whose reproductive ouput is recorded
  pbrood <- pbrood
      
  # Create productivity data
  P <- create.reproduction(ind$IND, rep(pbrood, yrs))
      
  ######################################################################
  #                                                                    #  
  #                       Specify IPM Code                             #
  #                                                                    #    
  ######################################################################
      
  # Write BUGS model file
  IPM.marr <- nimbleCode( { 
      
    # Priors and constraints
      
    # priors for mark-recapture model 
    mean.phi[1] ~ dbeta(1,1) 
    mean.phi[2] ~ dbeta(1,1) 
    mean.p ~ dbeta(1,1)
    
    #prior for repro model       
    mean.f ~ dunif(0,5)
    
    # prior for count model  
    p.count ~ dbeta(1,1) 
    
    # constraints for mark-recapture model 
    for (t in 1:(n.occasions-1)){ 
      sj[t] <- mean.phi[1]
      sa[t] <- mean.phi[2]
      p[t] <- mean.p
    }
    # constraints for repro model   
    for (t in 1:(n.occasions-1)){
      f[t] <- mean.f
    }
      
    # Likelihood 
      
    # State-space model for count data
      
    # Model for the initial population size
    N[1,1] ~ dcat(pNinit[])
    N[2,1] ~ dcat(pNinit[])
      
    # Process model over time
    for (t in 1:(n.occasions-1)){  
      N[1,t+1] ~ dpois(sj[t] * f[t] * (N[1,t] + N[2,t])) 
      N[2,t+1] ~ dbin(sa[t], (N[1,t] + N[2,t]))
    }
      
    # Observation model for count data 
    for(t in 1:n.occasions){  
      Ntot[t] <- N[1,t] + N[2,t]
      for(i in 1:nsurvs){
        count[t,i] ~ dbin(p.count,Ntot[t])
      }
    }
    
    for(t in 1:(n.occasions-1)){  
      fledge[t] ~ dpois(rho[t])
      rho[t] <- broods[t]*f[t]
    }  
    
    # Capture-recapture model (multinomial likelihood)
    # Define the multinomial likelihood
    
    
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
    
  })
  
  disc.unif <- function(A, B){
    pprob <- c(rep(0, A-1), rep(1/(B-A+1), (B-A+1)))
    return(pprob)
  }
      
  constants <- list(n.occasions = yrs, nsurvs=nsurvs, pNinit = disc.unif(1, 10000)) 
  
  # Bundle data                         
  data <- list(marr.j = marray[,,1], marr.a = marray[,,2], rel.j = rowSums(marray[,,1]), rel.a = rowSums(marray[,,2]), count = count, fledge = P$rep.agg[-c(yrs),1], broods = P$rep.agg[-c(yrs),2])
      
  # Initial values
  count.st <- count
  N.st <- matrix(NA,nrow=2,ncol=yrs)
  N.st[1,] <- apply(count.st,1,max)+1
  N.st[2,] <- apply(count.st,1,max)+1
  inits <- function(){list(N = N.st)}
      
  params <- c("mean.phi","mean.f","mean.p","p.count","N") 
  inits =  function() {list(mean.phi=runif(2),mean.p=runif(1))} 
  
  ni <- 50
  nb <- 20
  nc <- 3
  nt <- 1
  
  ## run model 
  samples.IPM.marr <- nimbleMCMC(
    code = IPM.marr,  
    data=data,
    constants = constants, 
    inits = inits,
    monitors = params,
    niter = ni,
    nburnin = nb,
    nchains = nc,
    thin = nt,
    summary = TRUE,
    samplesAsCodaMCMC = TRUE)
  
  ########################################################
  #nimble output 
  ########################################################
  
  #check diagnostics
  #MCMCtrace(object = samples.IPM.marr,
  #         pdf = FALSE,
  #         ind = TRUE)
  
  #summary 
  samples.IPM.marr$summary$all.chains 
  
  #Gelman-Rubin diagnostic
  #gelman.diag(samples.IPM.marr)
  
  all.results[s,,] <- samples.IPM.marr$summary$all.chains

  print(s)
}

