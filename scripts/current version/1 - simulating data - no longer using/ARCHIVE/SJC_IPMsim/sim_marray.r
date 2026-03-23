####Version 2/7/2018

library(here)
library(nimble)
library(MCMCvis)
library(coda)

rm(list=ls())

source(here("scripts/current version/1 - simulating data/SJC_IPMsim","create.pop.r"))
source(here("scripts/current version/1 - simulating data/SJC_IPMsim","create.eh.r"))

######################################################################
#                                                                    #  
#                        Simulation Parameters                       #
#                                                                    #    
######################################################################
sims <- 25
sims <- 1

# Number of years
yrs <- 10 

######################################################################
#                                                                    #  
#                        Create Population                           #
#                                                                    #    
######################################################################

# Record time 
start.time <- Sys.time()

all.results <- array(NA,dim = c(3,5,sims))
dimnames(all.results)[[1]] <- c("mean.p","mean.phi[1]","mean.phi[2]")
dimnames(all.results)[[2]] <- c("Mean","Median","St.Dev.","95%CI_low","95%CI_upp")

# Start simulations 
for(s in 1:sims){

  # Age specific survival probabilities (juv, adult)
  sj <- 0.35
  sa <- 0.65
  # Fecundity parameters 
  f <- 1 
  
  # Survival and fecundity   
  surv.mat <- matrix(NA,nrow=2,ncol=yrs)
  surv.mat[1,1:yrs] <- rep(sj,yrs)
  surv.mat[2,1:yrs] <- rep(sa,yrs)
  surv.mat <- surv.mat[,-yrs]
  
  fec.mat <- matrix(f,nrow=nrow(surv.mat),ncol=yrs)
  
  # Initial population size per age class 
  Ni <- c(3500, 6500)
      
  # Create the true population for this simulation and year  
  ind <- create.population(phi = surv.mat, f = fec.mat, Im = rep(0, yrs), Ni = Ni)
      
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
  
  EH <- ch$ch[,1:yrs]
  incl <- which(rowSums(EH)>=1)
  EH <- EH[incl,]
  age <- ch$age[incl]
      
  marray <- marray.age(EH,age)
      
  # Write model file
  marr.age <- nimbleCode( { 
      
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
    
  })
  
  # Bundle data                         
  constants <- list(n.occasions = yrs, rel.j = rowSums(marray[,,1]), rel.a = rowSums(marray[,,2])) 
  data <- list(marr.j = marray[,,1], marr.a = marray[,,2])
  
  #Parameters
  params <- c("mean.phi","mean.p") 
  
  #Initial values 
  inits =  function() {list(mean.phi=runif(2),mean.p=runif(1))} 
  
  ni <- 15000
  nb <- 2500
  nc <- 4
  nt <- 1
  
  ## run model 
  samples.marr.age <- nimbleMCMC(
    code = marr.age,  
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
  MCMCtrace(object = samples.marr.age$samples,
            pdf = FALSE,
            ind = TRUE)
  
  #summary 
  samples.marr.age$summary$all.chains 
  
  #Gelman-Rubin diagnostic
  gelman.diag(samples.marr.age$samples)
  
  all.results[,,s] <- samples.marr.age$summary$all.chains

  print(s)
}

apply(all.results,c(1,2),mean)
