# load libraries
library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)

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
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "3 - run models", "run_scenarios_helperFns_AJW.R"))


# MCMC settings #######
nb <- 200000 #burn-in 
ni <- 350000 #total iterations
nt <- 10  #thin
nc <- 3  #chains
# nb <- 2000 #burn-in 
# ni <- 8000 #total iterations
# nc <- 2

sims.per <- 40
n.sam <- 3 #number of annual surveys for n-mixture abundance model

cores = detectCores()
cl <- makeCluster(nrow(dem_scenarios), setup_strategy = "sequential") #not to overload your computer
registerDoParallel(cl)


#simulation replicates
foreach(i = 21:sims.per) %dopar% { # loop over replicate sims  #####
    
  library(here)
  library(nimble)
  library(IPMbook)
  
  #make true population trajectory data across demographic scenarios
  for (d in 1:nrow(dem_scenarios)) { #
    
    phi <- as.numeric(c(dem_scenarios[d,'phi1'], dem_scenarios[d,'phiad']))
    fec <- as.numeric(dem_scenarios[d,'fec'])*2 #function assumes both sexes
    
    pop1 <- simPop(Ni = Ni, phi = phi, f = fec, nYears = nyears)
    pop2 <- simPop(Ni = Ni, phi = phi, f = fec, nYears = nyears)
    pop3 <- simPop(Ni = Ni, phi = phi, f = fec, nYears = nyears)
    
    #survey the real populations across survey scenarios and run models
    for (s in 1:nrow(surv_scenarios)) {
      
      det.abund <- surv_scenarios_num[s,'det.abund']
      det.prod <- surv_scenarios_num[s,'det.prod']
      det.MR <- surv_scenarios_num[s,'det.MR']
      
      comb <- dem_scenarios[d,]

      # population survey data 
      Nad_count <- simCountBin(N=pop1$N[2,], pDetect = det.abund)
      N1_count <- simCountBin(N=pop1$N[1,], pDetect = det.abund)
      
      tot_count1 <- simCountBin(N = pop1$totAdults, pDetect = det.abund)
      tot_count2 <- simCountBin(N = pop1$totAdults, pDetect = det.abund)
      tot_count3 <- simCountBin(N = pop1$totAdults, pDetect = det.abund)
      
      surv_cnts <- rbind(tot_count1$count,
                             tot_count2$count,
                             tot_count3$count)
      
      # capture histories
      #check; not sure what to put for cap
      if (!is.na(det.MR)) {
      ch <- simCapHist(state=pop2$state, cap=1, recap=det.MR, maxAge=2, verbose = F)
      
      # m-arrays
      marr <- marrayAge(ch$ch, ch$age)
      marr.j <- as.matrix(marr[,,1])
      marr.a <- as.matrix(marr[,,2])
      }
      
      # productivity data; check females.only
      if (!is.na(det.prod)) { 
      nest_dat <- simProd(reprod = pop3$reprod, pInclude = det.prod, females.only = TRUE,
                     verbose = F)
      obs_nestlings <- as.numeric(nest_dat$prod.agg[,'Juveniles'])
      obs_nests <- as.numeric(nest_dat$prod.agg[,'Surveyed broods'])
      }
      
      ##run models
      #abundance data only
      if (is.na(det.prod) & is.na(det.MR)) { 
        
        out_abundOnly <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, 
                               comb, detect = det.abund)
        saveRDS(out_abundOnly, here("results", paste("out_abundOnly-",d,"-",s,"-",i,".RDS", sep = "")))
        rm(out_abundOnly) 
        
      } #abund only
        
      #missing productivity data
        else if (is.na(det.prod)) { 
          out_noProd <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, 
                                   comb, detect = det.abund)
          saveRDS(out_noProd, here("results", paste("out_noProd-",d,"-",s,"-",i,".RDS", sep = "")))
          rm(out_noProd)
          
        } #missing prod
      
      #missing MR data
        else if (is.na(det.MR)) { 
          out_noMR <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, 
                              comb, detect = det.abund)
          saveRDS(out_noMR, here("results", paste("out_noMR-",d,"-",s,"-",i,".RDS", sep = "")))
          rm(out_noMR)
          
        } #missing MR
        
      #full IPM
          else { 
            out_IPM <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, 
                                 comb, detect = det.abund)
            saveRDS(out_IPM, here("results", paste("out_IPM-",d,"-",s,"-",i,".RDS", sep = "")))
            rm(out_IPM)
          }
    } #s
  } #d
} #i
    
  
stopCluster(cl)
      
# MCMCsummary(out_IPM, Rhat = 1.1)
# gelman.diag(out_IPM)
