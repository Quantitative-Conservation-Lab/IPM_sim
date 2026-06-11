# load libraries
library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)

# load data and dem/survey scenarios
surv_scenarios <- readRDS(here("data", "data_scenarios.RDS")) %>%
  dplyr::filter(det.abund != 'L')
  
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
 
# Ni <- c(150, 150)
# Ni <- c(300, 300)

nyears <- 20
Nst.tot <- 300

# source functions
source(here("scripts", "current version",
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "3 - run models", "run_scenarios_helperFns_AJW.R"))

# MCMC settings #######
nb <- 100000 #burn-in
ni <- 250000 #total iterations
nc <- 4
# nb <- 125000
# ni <- 250000
# nc <- 3  #chains
nt <- 10  #thin

sims.per <- 120 #10#0 

cores = detectCores()
cl <- makeCluster(nrow(dem_scenarios), setup_strategy = "sequential") #not to overload your computer
registerDoParallel(cl)


#simulation replicates
foreach(i = 101:120) %dopar% { # loop over replicate sims  #####
  library(here)
  library(nimble)
  library(IPMbook)
  library(popbio)
  
  #make true population trajectory data across demographic scenarios
  for (d in 1:nrow(dem_scenarios)) { #
    
    phi <- as.numeric(c(dem_scenarios[d,'phi1'], dem_scenarios[d,'phiad']))
    fec <- as.numeric(dem_scenarios[d,'fec'])*2 #check; function assumes both sexes
    
    sj <- phi[1]
    sa <- phi[2]
    
    pop.mat <- matrix(NA,nrow=2,ncol=2)
    pop.mat[1,] <- c(sj*fec/2,sj*fec/2)
    pop.mat[2,] <- c(sa,sa)
    
    stable <- eigen.analysis(pop.mat)$stable
    
    Ni <- round(c(Nst.tot*stable))
    
    pop1 <- simPop(Ni = Ni, phi = phi, f = fec, nYears = nyears)
    pop2 <- simPop(Ni = Ni, phi = phi, f = fec, nYears = nyears)
    pop3 <- simPop(Ni = Ni, phi = phi, f = fec, nYears = nyears)
    
    #survey the real populations across survey scenarios and run models
    for (s in 1:nrow(surv_scenarios)) {
      
      det.abund <- surv_scenarios_num[s,'det.abund']
      det.abund <- surv_scenarios_num[s,'det.abund']
      det.prod <- surv_scenarios_num[s,'det.prod']
      det.MR <- surv_scenarios_num[s,'det.MR']
      
      comb <- dem_scenarios[d,]

      # population survey data 
      tot_count1 <- simCountBin(N = pop1$totAdults, pDetect = det.abund)
      tot_count2 <- simCountBin(N = pop1$totAdults, pDetect = det.abund)
      tot_count3 <- simCountBin(N = pop1$totAdults, pDetect = det.abund)
      tot_count4 <- simCountBin(N = pop1$totAdults, pDetect = det.abund)
      tot_count5 <- simCountBin(N = pop1$totAdults, pDetect = det.abund)
      
      surv_cnts <- rbind(tot_count1$count,
                        tot_count2$count,
                        tot_count4$count,
                        tot_count5$count,
                        tot_count3$count)
      
      n.sam <- dim(surv_cnts)[1]
      
      maxcount <- max(surv_cnts[,1])
      
      # capture histories
      if (!is.na(det.MR)) {
      ch <- simCapHist(state=pop2$state, cap=det.MR, recap=det.MR, maxAge=2, verbose = F)
      
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
        saveRDS(out_abundOnly, here("results", 'nmix', 'final', paste("out_abundOnly-",d,"-",s,"-",i,".RDS", sep = "")))
        rm(out_abundOnly)

      } #abund only
      #   
      # #missing productivity data
        else if (is.na(det.prod)) {
          out_noProd <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc,
                                   comb, detect = det.abund)
          saveRDS(out_noProd, here("results", 'nmix', 'final', paste("out_noProd-",d,"-",s,"-",i,".RDS", sep = "")))
          rm(out_noProd)

        } #missing prod
      # 
      # #missing MR data
        else if (is.na(det.MR)) {
          out_noMR <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc,
                              comb, detect = det.abund)
          saveRDS(out_noMR, here("results", 'nmix', 'final', 
          paste("out_noMR-",d,"-",s,"-",i,".RDS", sep = "")))
          rm(out_noMR)

        } #missing MR
      #   
      # #full IPM
          else {
            out_IPM <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, 
                                 comb, detect = det.abund)
            saveRDS(out_IPM, here("results", 'nmix', 'final', paste("out_IPM-",d,"-",s,"-",i,".RDS", sep = "")))
            rm(out_IPM)
          }
    } #s
  } #d
} #i
    
  
stopCluster(cl)
      
# MCMCsummary(out_IPM, Rhat = 1.1)
# gelman.diag(out_IPM)
