# load libraries
library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)
library(popbio)

source(here("scripts/current version/1 - simulating data/SJC_IPMsim","create.pop.r"))
source(here("scripts/current version/1 - simulating data/SJC_IPMsim","create.eh.r"))
source(here("scripts/current version/1 - simulating data/SJC_IPMsim","create.repro.r"))
source(here("scripts/current version/1 - simulating data/SJC_IPMsim","create.survey.r"))

# load data and dem/survey scenarios
surv_scenarios <- readRDS(here("data", "data_scenarios.RDS"))

surv_scenarios_num <- surv_scenarios %>%
   transform(det.MR = ifelse(det.MR == 'L', 0.3,
                          ifelse(det.MR == 'M', 0.5, 
                                 ifelse(det.MR == 'H', 0.8, NA)))) %>%
  transform(det.prod = ifelse(det.prod == 'L', 0.3,
                              ifelse(det.prod == 'M', 0.5, 
                                     ifelse(det.prod == 'H', 0.8, NA)))) %>%
  transform(det.abund = ifelse(det.abund == 'L', 0.3,
                                ifelse(det.abund == 'M', 0.5, 0.8)))


dem_scenarios <- readRDS(here("data", "demographic_scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename("phi1" = "S.J", 
         "phiad" = "S.A",
         "fec" = "f")
 
Ni <- c(150, 150)
nyears <- 15

# source functions
source(here("scripts", "current version",
            "2 - models", "IPM_marray_vSJC.R"))
source(here("scripts", "current version",
            "3 - run models", "run_scenarios_helperFns_AJW_wSJCfun.R"))


# MCMC settings #######
nb <- 200000 #burn-in 
ni <- 350000 #total iterations
nt <- 10  #thin
nc <- 3  #chains
#testing
nb <- 40000 #burn-in
ni <- 80000 #total iterations
# nc <- 2

sims.per <- 20

cores = detectCores()
cl <- makeCluster(nrow(dem_scenarios), setup_strategy = "sequential") #not to overload your computer
registerDoParallel(cl)


#simulation replicates
foreach(i = 1:sims.per) %dopar% { # loop over replicate sims  #####
    
  library(here)
  library(nimble)
  # library(IPMbook)
  
  #make true population trajectory data across demographic scenarios
  for (d in 1:nrow(dem_scenarios)) { #
    
    # Age specific survival probabilities (juv, adult)
    #double names are in here to most quickly use Sarah's functions and amanda's
    #initial values/model set-up and processing, can consolidate latter
    phi <- as.numeric(c(dem_scenarios[d,'phi1'], dem_scenarios[d,'phiad']))
    sj <- phi[1]
    sa <- phi[2]
    f <- as.numeric(dem_scenarios[d,'fec'])*2 #check; function assumes both sexes
    fec <- f
    
    pop.mat <- matrix(NA,nrow=2,ncol=2)
    pop.mat[1,] <- c(sj*f/2,sj*f/2)
    pop.mat[2,] <- c(sa,sa)
    
    stable <- eigen.analysis(pop.mat)$stable
    
    # Survival and fecundity   
    surv.mat <- matrix(NA,nrow=2,ncol=nyears)
    surv.mat[1,1:nyears] <- rep(sj,nyears)
    surv.mat[2,1:nyears] <- rep(sa,nyears)
    surv.mat <- surv.mat[,-nyears]
    
    fec.mat <- matrix(f/2,nrow=nrow(surv.mat),ncol=nyears)
    
    # Create the true population for this simulation and year  
    ind <- create.population(phi = surv.mat, f = fec.mat, Im = rep(0, nyears), Ni = Ni)
    
    #generate data based on different survey scenarios
    for (s in 1:nrow(surv_scenarios)) {
      
      det.abund <- surv_scenarios_num[s,'det.abund']
      det.prod <- surv_scenarios_num[s,'det.prod']
      det.MR <- surv_scenarios_num[s,'det.MR']
      
      comb <- dem_scenarios[d,]
    
    # generate capture and recapture data
    cjuv <- det.MR           # initial capture probability of juveniles
    cad <- det.MR            # initial capture probability of adults
    prec <- det.MR           # recapture probability 
    
    # Create the capture histories and the corresponding m-arrays
    if (!is.na(det.MR)) {
      
      ch <- create.capturehistory(ind$IND, c = matrix(c(rep(cjuv, nyears), rep(cad, nyears)), nrow = 2, byrow = TRUE), 
                                  p = matrix(c(rep(prec, nyears-1), rep(prec, nyears-1)), nrow = 2, byrow = TRUE))
      
      EH <- ch$ch[,1:nyears]
      incl <- which(rowSums(EH)>=1)
      EH <- EH[incl,]
      age <- ch$age[incl]
      
      marray <- marray.age(EH,age)
      marr.j <- as.matrix(marray[,,1])
      marr.a <- as.matrix(marray[,,2])
    }
    
    #create reproductive data
    if (!is.na(det.prod)) { 

    # Proportion of breeding population in sample
    pbrood <- det.prod
    
    # Create productivity data
    nest_dat <- create.reproduction(ind$IND, rep(pbrood, nyears))
    
    obs_nestlings <- nest_dat$rep.agg[,1]
    obs_broods <- nest_dat$rep.agg[,2]
    }
    

    #create survey data
    # Proportion of breeding population in sample
    pcount <- det.abund
    psurv <- rep(pcount,nyears)
    n.sam <- 3
    secondaries <- n.sam
    
    # Create productivity data
    count <- create.survey.bin(ind$Nu['Total',], psur = psurv, sec = secondaries)
    
      ##run models
      #abundance data only
      if (is.na(det.prod) & is.na(det.MR)) { 
        
        out_abundOnly <- runabundonly(nb = nb*2, ni = ni*2, nt = nt, nc = nc, 
                               comb, detect = det.abund)
        saveRDS(out_abundOnly, here("results", 'vSJC', paste("out_abundOnly-",d,"-",s,"-",i,".RDS", sep = "")))
        rm(out_abundOnly) 
        
      } #abund only
        
      #missing productivity data
        else if (is.na(det.prod)) { 
          out_noProd <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, 
                                   comb, detect = det.abund)
          saveRDS(out_noProd, here("results", 'vSJC', paste("out_noProd-",d,"-",s,"-",i,".RDS", sep = "")))
          rm(out_noProd)
          
        } #missing prod
      
      #missing MR data
        else if (is.na(det.MR)) { 
          out_noMR <- runnomr(nb = nb*2, ni = ni*2, nt = nt, nc = nc, 
                              comb, detect = det.abund)
          saveRDS(out_noMR, here("results", 'vSJC', paste("out_noMR-",d,"-",s,"-",i,".RDS", sep = "")))
          rm(out_noMR)
          
        } #missing MR
        
      #full IPM
          else { 
            out_IPM <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, 
                                 comb, detect = det.abund)
            saveRDS(out_IPM, here("results", 'vSJC', paste("out_IPM-",d,"-",s,"-",i,".RDS", sep = "")))
            rm(out_IPM)
          }
    } #s
  } #d
} #i
    
  
stopCluster(cl)
      
# MCMCsummary(out_IPM, Rhat = 1.1)
# gelman.diag(out_IPM)
