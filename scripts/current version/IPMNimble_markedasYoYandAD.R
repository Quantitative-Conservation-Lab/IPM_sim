# TODO
# abby put nest stuff in here

#IPM model for sim data
library(here)
#source(here("scripts", "dataSims","BPA_IPMsimData.R")) #full on BPA WinBUGS code
source(here("scripts", "current version","IPMsimData_markasYoYandAd.R")) #changed the MR data
source(here("scripts", "current version","productivityDataSim.R"))

# 11.3. Example of a simple IPM (counts, capture-recapture, reproduction)
# 11.3.1. Load data
# Population counts (from years 1 to n.years)

# AEB NOTE - if 1000 individuals, why are the counts so low?
y <- df$SUR
n.sam <- df$n.sam

# Capture-recapture data (in m-array format, from years 1 to n.years)
m <- df$ch
first <- df$first
age_ch <- df$age_ch

#Productivity/nest success not in the data sim right now
#these will need to be changed once we have the nest success data
# Productivity data (from years 1 to n.years-1)
#J<-nestlings#df$nestlings
#R<-R#df$R
 
H <- prod$observed.nest.status
Fledged <- prod$clutch.sizes
first.nest <- prod$first.nest
last.nest <- prod$last.nest
max.nest.age <- prod$max.nest.age
n.nests <- prod$N.nests.found
n.succ.nests <- prod$N.nests.successful

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
  N1[1]<-n1.start
  Nad[1]<-nad.start
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
    Fledged[c] ~ dpois(lambdaf) 
  } # c

  # Derived quantities #####
  fec <- 1/2 * lambdaf * (phi.nest^max.nest.age)
  
  # END derived quantities
  
})




#Provide known state as data  
state.data <- function(EH){
  state <- EH
  for (i in 1:dim(EH)[1]){
    first <- min(which(EH[i,]==1))
    last <- max(which(EH[i,]==1))
    state[i,first:last] <- 1
  }
  state[state==0] <- NA
  return(state)
}
#age for age specific survival probabilities
ageunknown<-function(age_ch){
  allages<-age_ch
  f1<-l1<-numeric(dim(allages)[1])
  t<-dim(age_ch)[2]
  for(i in 1:dim(allages)[1]){
    f1[i]<-min(which(allages[i,]>=2))
    l1[i]<-max(which(allages[i,]>0))
    if((f1[i]+1)<=t){
      allages[i,(f1[i]+1)]<-3
    } else{}
    if((f1[i]+2)<=t){
      allages[i,(f1[i]+2):t]<-3
    } else{}
  }
  allages<-allages-1
  allages[which(allages==-1)]<-0
  #allages[which(allages==1)]<-0
  # allages[which(allages==2)]<-1
  # allages[which(allages==3)]<-2
  return(allages)
}

age<-ageunknown(age_ch)

z.state <- state.data(m)

# data
#HAS- can probably clean all of this up when the data sim is ready to go, 
#so it automatically inputs the right data and constants 
datipm <- list(ch.y = m, y = y, 
               H = H, 
               Fledged = Fledged)
constants<-list(nyears = ncol(m), 
                n.ind=nrow(m), first=first, age=age, n.sam=n.sam,
                n.nests = n.nests, 
                n.succ.nests = n.succ.nests, 
                first.nest = first.nest, 
                last.nest = last.nest, 
                max.nest.age = max.nest.age)
# Initial values

# function for H inits - set every NA to 1
# ugly but whatever
Hinits <- H
Hinits[!is.na(H)] <- NA
Hinits[is.na(H)] <- 1

inits <- list(mean.phi=c(0.4, 0.77),
              mean.p = 0.5, 
              #mean.fec = runif(1, 0, 10), 
              p.surv=0.5,
              z=z.state,
              n1.start=1000,#sample(1:30,1),#super sensitive to these values, tried rpois(1,30) and it dodnt work
              nad.start=1000,#sample(1:30,1),
              phi.nest = 0.975, 
              lambdaf = 2.5 ,
              H = Hinits
              )
parameters <- c("mean.phi", 
                "mean.p", "fec", 
                "p.surv", 
                "phi.nest", "lambdaf", 
                "N1", "Nad")
mod<-nimbleModel(IPMmod, constants=constants, data=datipm, inits=inits)
conf<-configureMCMC(mod)
conf$addMonitors(parameters)
Rmcmc<-buildMCMC(conf)
Cmodel<-compileNimble(mod)#, maxContractions = 1000)
#for re-running without compiling all the above
Cmodel$setInits(inits)
#Cmodel$setData(newdata)
Cmcmc<-compileNimble(Rmcmc, project=Cmodel)
Cmcmc$run(thin=10, reset=T, niter=10000, nburnin=5000)
out<-as.data.frame(as.matrix(Cmcmc$mvSamples))

# TODO
# slice sampler reached maximum number of contractions for Ns :(

