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
 
low.lam.params <- dem_scenarios %>% 
  filter(trend == "decline")
med.lam.params <- dem_scenarios %>% 
  filter(trend == "stable")
high.lam.params <- dem_scenarios %>% 
  filter(trend == "increase")

Ni <- c(150, 150)
nyears <- 15

# source functions
# source(here("scripts", "current version",
#             "1 - simulating data", "IPM_sim_3.0function.R"))

source(here("scripts", "current version",
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "3 - run models", "run_scenarios_helperFns.R"))


# MCMC settings #######
nb <- 120000#0 #burn-in 
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

sims.per <- 100
sims.per <- 2

cores=detectCores()
cl <- makeCluster(nrow(low.lam.params), setup_strategy = "sequential") #not to overload your computer
registerDoParallel(cl)

# from scenario -> trajectory -> generate data -> fit all models
# i is the unique trajectory (within trend)
# d is scenario number

#simulation replicates
foreach(i = 1:sims.per) %dopar% { # loop over replicate sims  #####
    
  library(here)
  library(nimble)
  
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
      
      # population survey data 
      Ncount <- simCountBin(N=pop1$totB, pDetect = det.abund)
      
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
      pro <- simProd(reprod=pop3$reprod, pInclude = det.prod, females.only = TRUE,
                     verbose = F)
      }
      
      ##run models
      
      comb <- dem_scenarios[d,]
      
      #abundance data only
      if (is.na(det.prod) & is.na(det.MR)) { 
        
        out_abundOnly <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, 
                               popDat, popTraj, comb, detect = det.abund)
        saveRDS(out_abundOnly, here("results",paste("out_abundOnly-",d,"-",s,"-",i,".RDS", sep = "")))
        rm(out_abundOnly) 
        
      #missing productivity data
        else if (is.na(det.prod)) { 
          out_noProd <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, 
                                   popDat, popTraj, comb, detect = det.abund)
          saveRDS(out_noProd, here("results", paste("out_noProd-",d,"-",s,"-",i,".RDS", sep = "")))
          rm(out_noProd)
      
      #missing MR data
        else if (is.na(det.MR)) { 
          out_noMR <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, 
                              popDat, popTraj, comb, detect = det.abund)
          saveRDS(out_noMR, here("results", paste("out_noMR-",d,"-",s,"-",i,".RDS", sep = "")))
          rm(out_noMR)
        
      #full IPM
          else { 
            out_IPM <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, 
                                 popDat, popTraj, comb, detect = det.abund)
            saveRDS(out_IPM, here("results", paste("out_IPM-",d,"-",s,"-",i,".RDS", sep = "")))
            rm(out_IPM)
             }
        }
        }
      }
    } #s
  } #d
} #i
    
  
stopCluster(cl)
      
      
      
#       
#     
#     # load relevant population trajectories ####
#     lowpopTraj <- readRDS(here("data", "lowTrajectories", paste("lowpopTraj", "-", i, "-", j, ".RDS", sep = "")))
#     medpopTraj <- readRDS(here("data", "medTrajectories", paste("medpopTraj", "-", i, "-", j, ".RDS", sep = "")))
#     highpopTraj <- readRDS(here("data", "highTrajectories", paste("highpopTraj", "-", i, "-", j, ".RDS", sep = "")))
#     
#     for (d in 1:dim(surv_scenarios)[1]) { # loop over model scenario  #####
#       # print(paste0(d, " out of ", 10))
#       
#       # translate detection levels into numbes
#       det.levels <- data_scenarios[d, 1:3]
#       det.numeric <- det.levels[1:3]
#       det.numeric[which(det.numeric == "L")] <- detect.l
#       det.numeric[which(det.numeric == "M")] <- detect.m
#       det.numeric[which(det.numeric== "H")] <- detect.h
#       det.numeric[which(det.numeric== "NA")] <- NA
#       if (is.na(det.numeric[2]) & is.na(det.numeric[3])) { # ABUNDANCE ONLY #####
#         #if (det.levels[4] == "L") { # simulate low trajectory data #####
#           lowpopDat <- simData (indfates = lowpopTraj$indfates,
#                                 n.years = 15,
#                                 #n.data.types = c(0.25,0.25,0.25),
#                                 #remove n.data.types (not needed in newest)
#                                 #ADonly = F, 
#                                 #HAS: removed ADonly
#                                 p.1 = as.numeric(det.numeric[2]), #
#                                 p.ad = as.numeric(det.numeric[2]), #
#                                 p.count = as.numeric(det.numeric[1]), #
#                                 p.prod = as.numeric(det.numeric[3]), #
#                                 BinMod = T, # TODO could tidy these up as won't change
#                                 n.sam = 3,
#                                 sig = 0,
#                                 productivity = T)
#           popDat <- lowpopDat
#           popTraj <- lowpopTraj
#           comb <- low.lam.params[i,]
#           # run model and save results ####
#           lowout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(lowout, here("results",paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(lowout)
#         #} else if (det.levels[4] == "M") { # simulate medium trajectory data #####
#           medpopDat <- simData (indfates = medpopTraj$indfates,
#                                 n.years = 15,
#                                 #n.data.types = c(0.25,0.25,0.25),
#                                 #remove n.data.types (not needed in newest)
#                                 #ADonly = F, 
#                                 #HAS: removed ADonly
#                                 p.1 = as.numeric(det.numeric[2]), #
#                                 p.ad = as.numeric(det.numeric[2]), #
#                                 p.count = as.numeric(det.numeric[1]), #
#                                 p.prod = as.numeric(det.numeric[3]), #
#                                 BinMod = T,
#                                 n.sam = 3,
#                                 sig = 0,
#                                 productivity = T)
#           popDat <- medpopDat
#           popTraj <- medpopTraj
#           comb <- med.lam.params[i,]
#           # run model and save results ####
#           medout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(medout)
#         #} else if (det.levels[4] == "H") { # simulate high trajectory data #####
#           highpopDat <- simData (indfates = highpopTraj$indfates,
#                                  n.years = 15,
#                                  #n.data.types = c(0.25,0.25,0.25),
#                                  #remove n.data.types (not needed in newest)
#                                  #ADonly = F, 
#                                  #HAS: removed ADonly
#                                  p.1 = as.numeric(det.numeric[2]), #
#                                  p.ad = as.numeric(det.numeric[2]), #
#                                  p.count = as.numeric(det.numeric[1]), #
#                                  p.prod = as.numeric(det.numeric[3]), #
#                                  BinMod = T,
#                                  n.sam = 3,
#                                  sig = 0,
#                                  productivity = T)
#           popDat <- highpopDat
#           popTraj <- highpopTraj
#           comb <- high.lam.params[i,]
#           # run model and save results ####
#           highout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(highout)
#         #}
#       } else if (is.na(det.numeric[2])) { # NO MARK RECAPTURE ######
#         #if (det.levels[4] == "L") { # simulate low trajectory data #####
#           lowpopDat <- simData (indfates = lowpopTraj$indfates,
#                                 n.years = 15,
#                                 #n.data.types = c(0.25,0.25,0.25),
#                                 #remove n.data.types (not needed in newest)
#                                 #ADonly = F, 
#                                 #HAS: removed ADonly
#                                 p.1 = as.numeric(det.numeric[2]), #
#                                 p.ad = as.numeric(det.numeric[2]), #
#                                 p.count = as.numeric(det.numeric[1]), #
#                                 p.prod = as.numeric(det.numeric[3]), #
#                                 BinMod = T,
#                                 n.sam = 3,
#                                 sig = 0,
#                                 productivity = T)
#           popDat <- lowpopDat
#           popTraj <- lowpopTraj
#           comb <- low.lam.params[i,]
#           # run model and save results ####
#           lowout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(lowout, here("results", paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(lowout)
#         #} else if (det.levels[4] == "M") { # simulate medium trajectory data #####
#           medpopDat <- simData (indfates = medpopTraj$indfates,
#                                 n.years = 15,
#                                 #n.data.types = c(0.25,0.25,0.25),
#                                 #remove n.data.types (not needed in newest)
#                                 #ADonly = F, 
#                                 #HAS: removed ADonly
#                                 p.1 = as.numeric(det.numeric[2]), #
#                                 p.ad = as.numeric(det.numeric[2]), #
#                                 p.count = as.numeric(det.numeric[1]), #
#                                 p.prod = as.numeric(det.numeric[3]), #
#                                 BinMod = T,
#                                 n.sam = 3,
#                                 sig = 0,
#                                 productivity = T)
#           popDat <- medpopDat
#           popTraj <- medpopTraj
#           comb <- med.lam.params[i,]
#           # run model and save results ####
#           medout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(medout)
#         #} else if (det.levels[4] == "H") { # simulate high trajectory data #####
#           highpopDat <- simData (indfates = highpopTraj$indfates,
#                                  n.years = 15,
#                                  #n.data.types = c(0.25,0.25,0.25),
#                                  #remove n.data.types (not needed in newest)
#                                  #ADonly = F, 
#                                  #HAS: removed ADonly
#                                  p.1 = as.numeric(det.numeric[2]), #
#                                  p.ad = as.numeric(det.numeric[2]), #
#                                  p.count = as.numeric(det.numeric[1]), #
#                                  p.prod = as.numeric(det.numeric[3]), #
#                                  BinMod = T,
#                                  n.sam = 3,
#                                  sig = 0,
#                                  productivity = T)
#           popDat <- highpopDat
#           popTraj <- highpopTraj
#           comb <- high.lam.params[i,]
#           # run model and save results ####
#           highout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(highout)
#         #}
#       } else if (is.na(det.numeric[3])) { # NO NEST SURVIVAL ######
#         #if (det.levels[4] == "L") { # simulate low trajectory data #####
#           lowpopDat <- simData (indfates = lowpopTraj$indfates,
#                                 n.years = 15,
#                                 #n.data.types = c(0.25,0.25,0.25),
#                                 #remove n.data.types (not needed in newest)
#                                 #ADonly = F, 
#                                 #HAS: removed ADonly
#                                 p.1 = as.numeric(det.numeric[2]), #
#                                 p.ad = as.numeric(det.numeric[2]), #
#                                 p.count = as.numeric(det.numeric[1]), #
#                                 p.prod = as.numeric(det.numeric[3]), #
#                                 BinMod = T,
#                                 n.sam = 3,
#                                 sig = 0,
#                                 productivity = T)
#           popDat <- lowpopDat
#           popTraj <- lowpopTraj
#           
#           # TODO change here
#           # TODO will need to change the loop structure to deal w this
#           #comb <- low.lam.params[i,]
#           comb <- low.lam.params[i,]
#           
#           # run model and save results ####
#           lowout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(lowout, here("results", paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(lowout)
#         #} else if (det.levels[4] == "M") { # simulate medium trajectory data #####
#           medpopDat <- simData (indfates = medpopTraj$indfates,
#                                 n.years = 15,
#                                 #n.data.types = c(0.25,0.25,0.25),
#                                 #remove n.data.types (not needed in newest)
#                                 #ADonly = F, 
#                                 #HAS: removed ADonly
#                                 p.1 = as.numeric(det.numeric[2]), #
#                                 p.ad = as.numeric(det.numeric[2]), #
#                                 p.count = as.numeric(det.numeric[1]), #
#                                 p.prod = as.numeric(det.numeric[3]), #
#                                 BinMod = T,
#                                 n.sam = 3,
#                                 sig = 0,
#                                 productivity = T)
#           popDat <- medpopDat
#           popTraj <- medpopTraj
#           comb <- med.lam.params[i,]
#           # run model and save results ####
#           medout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(medout)
#         #} else if (det.levels[4] == "H") { # simulate high trajectory data #####
#           highpopDat <- simData (indfates = highpopTraj$indfates,
#                                  n.years = 15,
#                                  #n.data.types = c(0.25,0.25,0.25),
#                                  #remove n.data.types (not needed in newest)
#                                  #ADonly = F, 
#                                  #HAS: removed ADonly
#                                  p.1 = as.numeric(det.numeric[2]), #
#                                  p.ad = as.numeric(det.numeric[2]), #
#                                  p.count = as.numeric(det.numeric[1]), #
#                                  p.prod = as.numeric(det.numeric[3]), #
#                                  BinMod = T,
#                                  n.sam = 3,
#                                  sig = 0,
#                                  productivity = T)
#           popDat <- highpopDat
#           popTraj <- highpopTraj
#           comb <- high.lam.params[i,]
#           # run model and save results ####
#           highout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(highout)
#        # }
#       } else { # FULL IPM ########
#        # if (det.levels[4] == "L") { # simulate low trajectory data #####
#           lowpopDat <- simData (indfates = lowpopTraj$indfates,
#                                 n.years = 15,
#                                 #n.data.types = c(0.25,0.25,0.25),
#                                 #remove n.data.types (not needed in newest)
#                                 #ADonly = F,
#                                 #HAS: removed ADonly
#                                 p.1 = as.numeric(det.numeric[2]), #
#                                 p.ad = as.numeric(det.numeric[2]), #
#                                 p.count = as.numeric(det.numeric[1]), #
#                                 p.prod = as.numeric(det.numeric[3]), #
#                                 BinMod = T,
#                                 n.sam = 3,
#                                 sig = 0,
#                                 productivity = T)
#           popDat <- lowpopDat
#           popTraj <- lowpopTraj
# 
#           comb <- low.lam.params[i,] # TODO
# 
#           # run model and save results ####
#           lowout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(lowout, here("results", paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(lowout)
#        # } else if (det.levels[4] == "M") { # simulate medium trajectory data #####
#           medpopDat <- simData (indfates = medpopTraj$indfates,
#                                 n.years = 15,
#                                 #n.data.types = c(0.25,0.25,0.25),
#                                 #remove n.data.types (not needed in newest)
#                                 #ADonly = F, 
#                                 #HAS: removed ADonly
#                                 p.1 = as.numeric(det.numeric[2]), #
#                                 p.ad = as.numeric(det.numeric[2]), #
#                                 p.count = as.numeric(det.numeric[1]), #
#                                 p.prod = as.numeric(det.numeric[3]), #
#                                 BinMod = T,
#                                 n.sam = 3,
#                                 sig = 0,
#                                 productivity = T)
#           popDat <- medpopDat
#           popTraj <- medpopTraj
#           
#           comb <- med.lam.params[i,]
#           
#           # run model and save results ####
#           medout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#           saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
#           rm(medout)
#        # } else if (det.levels[4] == "H") { # simulate high trajectory data #####
#        highpopDat <- simData (indfates = highpopTraj$indfates,
#                               n.years = 15,
#                               #n.data.types = c(0.25,0.25,0.25),
#                               #remove n.data.types (not needed in newest)
#                               #ADonly = F,
#                               #HAS: removed ADonly
#                               p.1 = as.numeric(det.numeric[2]), #
#                               p.ad = as.numeric(det.numeric[2]), #
#                               p.count = as.numeric(det.numeric[1]), #
#                               p.prod = as.numeric(det.numeric[3]), #
#                               BinMod = T,
#                               n.sam = 3,
#                               sig = 0,
#                               productivity = T)
#        popDat <- highpopDat
#        popTraj <- highpopTraj
# 
#        comb <- high.lam.params[i,]
# 
#        # run model and save results ####
#        highout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
#        saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
#        rm(highout)
#        # }
#       } # else
#     } # scenarios row (d)
#   } # sims per (j)
# } # foreach - scenarios picked (i)

