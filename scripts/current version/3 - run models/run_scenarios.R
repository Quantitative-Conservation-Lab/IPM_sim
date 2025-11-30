# load libraries
library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)

# load data
# TODO - adjust as in previous file
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
# TODO fix hard coding
low.lam.params <- scenarios %>% 
  filter(trend == "decline")
med.lam.params <- scenarios %>% 
  filter(trend == "stable")
high.lam.params <- scenarios %>% 
  filter(trend == "increase")

# source functions
# TODO - don't source from here? or do? what is easier to avoid hardcoding
# source(here("scripts", "current version",
#             "0 - preparing scenarios", "compute_time_calc.R"))

# source(here("scripts", "current version",
#             "1 - simulating data", "IPM_sim_2.0function.R"))
#changed to latest simulate script
source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_3.0function.R"))
source(here("scripts", "current version",
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "3 - run models", "run_scenarios_helperFns.R"))

# TODO cut this bit
# determine priority score for scenarios
# scenarios %>% mutate(priority = NA_integer_)
# for (i in 1:nrow(scenarios)) {
#   tmp <- scenarios[i, 1:3]
#   tmp <- tmp[!is.na(tmp)]
#   scenarios[i, "priority"] <- length(unique(tmp))
# }
# scenarios <- scenarios %>% arrange(priority) %>% # save in prioritized order
#   transform(simscenarios = 1:144)

# write.csv(scenarios, here::here('data', 'scenario_ID.csv'), row.names = F)

# which.prio.1 <- which(scenarios$priority == 1)
# which.prio.2 <- which(scenarios$priority == 2)
# which.prio.3 <- which(scenarios$priority == 3)

# TODO note hardcoding here. this seems like an ok place to define it though
# detection levels
detect.l <- 0.3
detect.m <- 0.5
detect.h <- 0.8

detect <- c(detect.l, detect.m, detect.h)

data_scenarios <- readRDS(here("data", "data_scenarios.RDS"))

# MCMC settings #######
# TODO this is defined differently elsewhere
# does it need to be longer for convergence, especially given stricter cutoff
nb <- 100000#0 #burn-in # TODO play with this
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

cores=detectCores()
cl <- makeCluster(scenarios.picked, setup_strategy = "sequential") #not to overload your computer
registerDoParallel(cl)

# TODO rename these to more sensible things (pop scenario, data scenario, sims per?)
# this one needs some rewriting
# TODO - is it more streamlined / is it painful recoding to have one pipeline 
# from scenario -> trajectory -> generate data -> fit all models
# i is the unique trajectory (within trend)
# j is replicate
# d is scenario number

# TODO tidy this up and create a globals doc
scenarios.picked <- nrow(low.lam.params) # TODO note change here
sims.per <- 50 # TODO - is this what we decided? didn't we discuss either less

# TODO this for loop requires some recoding given new structure of 
# scenarios and data scenarios
#foreach(i = 1:scenarios.picked) %dopar% { # loop over population trajectory  #####
foreach(j = 1:sims.per) %dopar% { # loop over population trajectory  #####
    
  library(here)
  library(nimble)
  #for (j in 1:sims.per) { # loop over replicate  #####
  for (i in 1:scenarios.picked) { # loop over replicate  #####
    # load relevant population trajectories ####
    lowpopTraj <- readRDS(here("data", "lowTrajectories", paste("lowpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    medpopTraj <- readRDS(here("data", "medTrajectories", paste("medpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    highpopTraj <- readRDS(here("data", "highTrajectories", paste("highpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    for (d in 1:dim(data_scenarios)[1]) { # loop over model scenario  #####
    
    # TESTING
    # for (d in 1:10) { # loop over model scenario  #####
    #   print(paste0(d, " out of ", 10))
      
      # translate detection levels into numbes
      det.levels <- data_scenarios[d, 1:3]
      det.numeric <- det.levels[1:3]
      det.numeric[which(det.numeric == "L")] <- detect.l
      det.numeric[which(det.numeric == "M")] <- detect.m
      det.numeric[which(det.numeric== "H")] <- detect.h
      det.numeric[which(det.numeric== "NA")] <- NA
      if (is.na(det.numeric[2]) & is.na(det.numeric[3])) { # ABUNDANCE ONLY #####
        #if (det.levels[4] == "L") { # simulate low trajectory data #####
          lowpopDat <- simData (indfates = lowpopTraj$indfates,
                                n.years = 15,
                                #n.data.types = c(0.25,0.25,0.25),
                                #remove n.data.types (not needed in newest)
                                #ADonly = F, 
                                #HAS: removed ADonly
                                p.1 = as.numeric(det.numeric[2]), #
                                p.ad = as.numeric(det.numeric[2]), #
                                p.count = as.numeric(det.numeric[1]), #
                                p.prod = as.numeric(det.numeric[3]), #
                                BinMod = T, # TODO could tidy these up as won't change
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- lowpopDat
          popTraj <- lowpopTraj
          comb <- low.lam.params[i,]
          # run model and save results ####
          lowout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(lowout, here("results",paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(lowout)
        #} else if (det.levels[4] == "M") { # simulate medium trajectory data #####
          medpopDat <- simData (indfates = medpopTraj$indfates,
                                n.years = 15,
                                #n.data.types = c(0.25,0.25,0.25),
                                #remove n.data.types (not needed in newest)
                                #ADonly = F, 
                                #HAS: removed ADonly
                                p.1 = as.numeric(det.numeric[2]), #
                                p.ad = as.numeric(det.numeric[2]), #
                                p.count = as.numeric(det.numeric[1]), #
                                p.prod = as.numeric(det.numeric[3]), #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- medpopDat
          popTraj <- medpopTraj
          comb <- med.lam.params[i,]
          # run model and save results ####
          medout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(medout)
        #} else if (det.levels[4] == "H") { # simulate high trajectory data #####
          highpopDat <- simData (indfates = highpopTraj$indfates,
                                 n.years = 15,
                                 #n.data.types = c(0.25,0.25,0.25),
                                 #remove n.data.types (not needed in newest)
                                 #ADonly = F, 
                                 #HAS: removed ADonly
                                 p.1 = as.numeric(det.numeric[2]), #
                                 p.ad = as.numeric(det.numeric[2]), #
                                 p.count = as.numeric(det.numeric[1]), #
                                 p.prod = as.numeric(det.numeric[3]), #
                                 BinMod = T,
                                 n.sam = 3,
                                 sig = 0,
                                 productivity = T)
          popDat <- highpopDat
          popTraj <- highpopTraj
          comb <- high.lam.params[i,]
          # run model and save results ####
          highout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(highout)
        #}
      } else if (is.na(det.numeric[2])) { # NO MARK RECAPTURE ######
        #if (det.levels[4] == "L") { # simulate low trajectory data #####
          lowpopDat <- simData (indfates = lowpopTraj$indfates,
                                n.years = 15,
                                #n.data.types = c(0.25,0.25,0.25),
                                #remove n.data.types (not needed in newest)
                                #ADonly = F, 
                                #HAS: removed ADonly
                                p.1 = as.numeric(det.numeric[2]), #
                                p.ad = as.numeric(det.numeric[2]), #
                                p.count = as.numeric(det.numeric[1]), #
                                p.prod = as.numeric(det.numeric[3]), #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- lowpopDat
          popTraj <- lowpopTraj
          comb <- low.lam.params[i,]
          # run model and save results ####
          lowout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(lowout, here("results", paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(lowout)
        #} else if (det.levels[4] == "M") { # simulate medium trajectory data #####
          medpopDat <- simData (indfates = medpopTraj$indfates,
                                n.years = 15,
                                #n.data.types = c(0.25,0.25,0.25),
                                #remove n.data.types (not needed in newest)
                                #ADonly = F, 
                                #HAS: removed ADonly
                                p.1 = as.numeric(det.numeric[2]), #
                                p.ad = as.numeric(det.numeric[2]), #
                                p.count = as.numeric(det.numeric[1]), #
                                p.prod = as.numeric(det.numeric[3]), #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- medpopDat
          popTraj <- medpopTraj
          comb <- med.lam.params[i,]
          # run model and save results ####
          medout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(medout)
        #} else if (det.levels[4] == "H") { # simulate high trajectory data #####
          highpopDat <- simData (indfates = highpopTraj$indfates,
                                 n.years = 15,
                                 #n.data.types = c(0.25,0.25,0.25),
                                 #remove n.data.types (not needed in newest)
                                 #ADonly = F, 
                                 #HAS: removed ADonly
                                 p.1 = as.numeric(det.numeric[2]), #
                                 p.ad = as.numeric(det.numeric[2]), #
                                 p.count = as.numeric(det.numeric[1]), #
                                 p.prod = as.numeric(det.numeric[3]), #
                                 BinMod = T,
                                 n.sam = 3,
                                 sig = 0,
                                 productivity = T)
          popDat <- highpopDat
          popTraj <- highpopTraj
          comb <- high.lam.params[i,]
          # run model and save results ####
          highout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(highout)
        #}
      } else if (is.na(det.numeric[3])) { # NO NEST SURVIVAL ######
        #if (det.levels[4] == "L") { # simulate low trajectory data #####
          lowpopDat <- simData (indfates = lowpopTraj$indfates,
                                n.years = 15,
                                #n.data.types = c(0.25,0.25,0.25),
                                #remove n.data.types (not needed in newest)
                                #ADonly = F, 
                                #HAS: removed ADonly
                                p.1 = as.numeric(det.numeric[2]), #
                                p.ad = as.numeric(det.numeric[2]), #
                                p.count = as.numeric(det.numeric[1]), #
                                p.prod = as.numeric(det.numeric[3]), #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- lowpopDat
          popTraj <- lowpopTraj
          
          # TODO change here
          # TODO will need to change the loop structure to deal w this
          #comb <- low.lam.params[i,]
          comb <- low.lam.params[i,]
          
          # run model and save results ####
          lowout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(lowout, here("results", paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(lowout)
        #} else if (det.levels[4] == "M") { # simulate medium trajectory data #####
          medpopDat <- simData (indfates = medpopTraj$indfates,
                                n.years = 15,
                                #n.data.types = c(0.25,0.25,0.25),
                                #remove n.data.types (not needed in newest)
                                #ADonly = F, 
                                #HAS: removed ADonly
                                p.1 = as.numeric(det.numeric[2]), #
                                p.ad = as.numeric(det.numeric[2]), #
                                p.count = as.numeric(det.numeric[1]), #
                                p.prod = as.numeric(det.numeric[3]), #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- medpopDat
          popTraj <- medpopTraj
          comb <- med.lam.params[i,]
          # run model and save results ####
          medout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(medout)
        #} else if (det.levels[4] == "H") { # simulate high trajectory data #####
          highpopDat <- simData (indfates = highpopTraj$indfates,
                                 n.years = 15,
                                 #n.data.types = c(0.25,0.25,0.25),
                                 #remove n.data.types (not needed in newest)
                                 #ADonly = F, 
                                 #HAS: removed ADonly
                                 p.1 = as.numeric(det.numeric[2]), #
                                 p.ad = as.numeric(det.numeric[2]), #
                                 p.count = as.numeric(det.numeric[1]), #
                                 p.prod = as.numeric(det.numeric[3]), #
                                 BinMod = T,
                                 n.sam = 3,
                                 sig = 0,
                                 productivity = T)
          popDat <- highpopDat
          popTraj <- highpopTraj
          comb <- high.lam.params[i,]
          # run model and save results ####
          highout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(highout)
       # }
      } else { # FULL IPM ########
       # if (det.levels[4] == "L") { # simulate low trajectory data #####
          lowpopDat <- simData (indfates = lowpopTraj$indfates,
                                n.years = 15,
                                #n.data.types = c(0.25,0.25,0.25),
                                #remove n.data.types (not needed in newest)
                                #ADonly = F, 
                                #HAS: removed ADonly
                                p.1 = as.numeric(det.numeric[2]), #
                                p.ad = as.numeric(det.numeric[2]), #
                                p.count = as.numeric(det.numeric[1]), #
                                p.prod = as.numeric(det.numeric[3]), #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- lowpopDat
          popTraj <- lowpopTraj
          
          comb <- low.lam.params[i,] # TODO
          
          # run model and save results ####
          lowout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(lowout, here("results", paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(lowout)
       # } else if (det.levels[4] == "M") { # simulate medium trajectory data #####
          medpopDat <- simData (indfates = medpopTraj$indfates,
                                n.years = 15,
                                #n.data.types = c(0.25,0.25,0.25),
                                #remove n.data.types (not needed in newest)
                                #ADonly = F, 
                                #HAS: removed ADonly
                                p.1 = as.numeric(det.numeric[2]), #
                                p.ad = as.numeric(det.numeric[2]), #
                                p.count = as.numeric(det.numeric[1]), #
                                p.prod = as.numeric(det.numeric[3]), #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- medpopDat
          popTraj <- medpopTraj
          
          comb <- med.lam.params[i,]
          
          # run model and save results ####
          medout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(medout)
       # } else if (det.levels[4] == "H") { # simulate high trajectory data #####
          highpopDat <- simData (indfates = highpopTraj$indfates,
                                 n.years = 15,
                                 #n.data.types = c(0.25,0.25,0.25),
                                 #remove n.data.types (not needed in newest)
                                 #ADonly = F, 
                                 #HAS: removed ADonly
                                 p.1 = as.numeric(det.numeric[2]), #
                                 p.ad = as.numeric(det.numeric[2]), #
                                 p.count = as.numeric(det.numeric[1]), #
                                 p.prod = as.numeric(det.numeric[3]), #
                                 BinMod = T,
                                 n.sam = 3,
                                 sig = 0,
                                 productivity = T)
          popDat <- highpopDat
          popTraj <- highpopTraj
          
          comb <- high.lam.params[i,]
          
          # run model and save results ####
          highout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
          saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(highout)
       # }
      } # else
    } # scenarios row (d)
  } # sims per (j)
} # foreach - scenarios picked (i)
stopCluster(cl)
