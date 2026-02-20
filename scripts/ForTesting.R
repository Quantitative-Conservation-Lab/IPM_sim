library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)

scenarios <- readRDS(here("data", "demographic_scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename(
    "phi1" = "S.J", 
    "phiad" = "S.A",
    "fec" = "f"
  )
# low.lam.combos <- readRDS(here("data","low.lam.params.RDS"))
# med.lam.combos <- readRDS(here("data","med.lam.params.RDS"))
# high.lam.combos <- readRDS(here("data","high.lam.params.RDS"))


#HAS: just picked on, could do whichever

# low.lam.params <- scenarios %>% 
#   filter(trend == "decline")
med.lam.params <- scenarios %>% 
  filter(trend == "stable")
# high.lam.params <- scenarios %>% 
#   filter(trend == "increase")


source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_3.0function.R"))
source(here("scripts", "current version",
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "3 - run models", "run_scenarios_helperFns.R"))


# detection levels
detect.l <- 0.3
detect.m <- 0.5
detect.h <- 0.8

detect <- c(detect.l, detect.m, detect.h)

data_scenarios <- readRDS(here("data", "data_scenarios.RDS"))

# MCMC settings #######

nb <- 100000#0 #burn-in # TODO play with this
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

#want to check data_scenarios 1 (lll), 17 (mmm), 33 (hhh)
# see if there is anything there

#just look at moderate LH (similar trend among fast and slow)
#i is scenarios picked
#j is sims per

i<-1
j<-4

#HAS: here I'm not using the already created true pops
#but can turn that back on
#medpopTraj <- readRDS(here("data", "medTrajectories", paste("medpopTraj", "-", i, "-", j, ".RDS", sep = "")))

#HAS: here I am setting phi1 and phiad equal to see if that 
# fixes the issues - if it does then may need to redo some of the 
#sim code

temptrue<-simPopTrajectory(n.years=15,
                           age.init=c(150,150),
                           phi.1=0.35,phi.ad=0.35,
                           f=2)
#check whether phi1=phiad removes the bias

#below, just picking detection levels for the non-mr data
d<-17
# translate detection levels into numbes
det.levels <- data_scenarios[d, 1:3]
det.numeric <- det.levels[1:3]
det.numeric[which(det.numeric == "L")] <- detect.l
det.numeric[which(det.numeric == "M")] <- detect.m
det.numeric[which(det.numeric== "H")] <- detect.h
det.numeric[which(det.numeric== "NA")] <- NA

#only look at full IPM

#Simulating data collection from 'truth'
medpopDat <- simData (indfates = temptrue$indfates,# medpopTraj$indfates,
                      n.years = 15,
                      #n.data.types = c(0.25,0.25,0.25),
                      #remove n.data.types (not needed in newest)
                      #ADonly = F, 
                      #HAS: removed ADonly
                      p.1 =1,# as.numeric(det.numeric[2]), #
                      p.ad = 1,#as.numeric(det.numeric[2]), #
                      p.count = as.numeric(det.numeric[1]), #
                      p.prod = as.numeric(det.numeric[3]), #
                      BinMod = T,
                      n.sam = 3,
                      sig = 0,
                      productivity = T)

#HAS: lazy coding here, just changing the 'truth'
comb <- med.lam.params[i,]
comb$phi1<-0.35
comb$phiad<-0.35
comb$fec<-2

runIPMmod <- function(nb, ni, nt, nc,
                      popDat, popTraj,
                      comb, detect) {
  #### DATA ####
  dat1 <- list(y = popDat$SUR,
               marr.a = marray(popDat$ch.a), 
               marr.j=marray(popDat$ch.j), # TODO note hardcode
               R.j=rowSums(marray(popDat$ch.j)), # TODO could try making same as adults for testig
               R.a = rowSums(marray(popDat$ch.a)),
               OBS_nestlings = popDat$OBS_nestlings,
               R_obs = popDat$R_obs
  )
  
  
  #### CONSTANTS ####
  
  const1 <- list(nyears = 15,
                 n.sam = nrow(popDat$SUR))
  
  #### INITIAL VALUES ####
  #z.state <- state.data(popDat$ch)
  
  inits1 <- list(
    mean.phi = c(comb$phi1, comb$phiad),#c(detect.h, detect.h),
    mean.p = detect[2],
    p.surv = detect[1],
    fec = comb$fec,#detect.h,
    #mean.phi = runif(2,0,1),#c(detect.h, detect.h),
    #mean.p = runif(1,0,1),#detect.h,
    #fec = runif(1,0,5),#detect.h,
    #z=z.state,
    n1.start=popTraj$Nouts[1,1], 
    nad.start=popTraj$Nouts[2,1]
  )
  
  #### PARAMETERS TO MONITOR ####
  params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda","Ntot","N1","Nad")#0.3764911
  
  #### COMPILE CONFIGURE AND BUILD ####
  Rmodel1 <- nimbleModel(code = IPMmod, constants = const1, data = dat1,
                         check = FALSE, calculate = FALSE, inits = inits1)
  conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
  #control = list(maxContractions = 1000))
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
  
  #### RUN MCMC ####
  outIPM <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1, thin=nt,
                    setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  
  return(outIPM)
}

popDat <- medpopDat
#HAS: for popTraj, can change back to those already created sims
#here, just using that phi1=phiad trajectory from above
popTraj <- temptrue#medpopTraj


# run model and save results ####
medout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, 
                    popDat, popTraj, comb, detect = as.numeric(det.numeric))


medout1<-data.frame(do.call(rbind, medout))


hist(medout1$N1.5.)
#abline(v=medpopTraj$Nouts[1,5], col="red")
abline(v=temptrue$Nouts[1,5])

hist(medout1$mean.phi.1.)
abline(v=comb$phi1, col="red")

hist(medout1$mean.phi.2.)
abline(v=comb$phiad, col="red")

hist(medout1$Nad.10.)
#abline(v=medpopTraj$Nouts[2,10], col="red")
abline(v=temptrue$Nouts[2,10], col="red")

