# aeb
# april 2, 2020

# OUTLINE ####

# NOTES ON WHAT HAS RUN AND WHERE IT IS LOCATED
#### ii ####
# first half
# second half
#############

library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
low.lam.combos <- readRDS(here("data","low.lam.params.RDS"))
med.lam.combos <- readRDS(here("data","med.lam.params.RDS"))
high.lam.combos <- readRDS(here("data","high.lam.params.RDS"))

# functions
source(here("scripts", "current version",
            "0 - preparing scenarios", "compute_time_calc.R"))
source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_2.0function.R"))
source(here("scripts", "current version",
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "4 - run models", "run_scenarios_helperFns.R"))

# determine priority score for scenarios
scenarios %>% mutate(priority = NA_integer_)
for (i in 1:nrow(scenarios)) {
  tmp <- scenarios[i, 1:3]
  tmp <- tmp[!is.na(tmp)]
  scenarios[i, "priority"] <- length(unique(tmp))
}
scenarios <- scenarios %>% arrange(priority) # save in prioritized order
which.prio.1 <- which(scenarios$priority == 1)
which.prio.2 <- which(scenarios$priority == 2)
which.prio.3 <- which(scenarios$priority == 3)

# simulate data
detect.l <- 0.3
detect.m <- 0.5
detect.h <- 0.8

detect <- c(detect.l, detect.m, detect.h)

nb <- 200000#0 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

cores=detectCores()
cl <- makeCluster(scenarios.picked, setup_strategy = "sequential") #not to overload your computer
registerDoParallel(cl)

# i is the unique trajectory (within trend)
# j is replicate
# d is scenario number


toRerun <- readRDS(here("toRerun.RDS"))

foreach(ii = 1:((dim(rerun)[1])/2)) %dopar% { #scenarios picked
  library(here)
  library(nimble)
  
  i <- toRerun[ii, 1]
  j <- toRerun[ii, 2]
  d <- toRerun[ii, 3]
  
  lowpopTraj <- readRDS(here("data", "lowTrajectories", paste("lowpopTraj", "-", i, "-", j, ".RDS", sep = "")))
  medpopTraj <- readRDS(here("data", "medTrajectories", paste("medpopTraj", "-", i, "-", j, ".RDS", sep = "")))
  highpopTraj <- readRDS(here("data", "highTrajectories", paste("highpopTraj", "-", i, "-", j, ".RDS", sep = "")))
  
  det.levels <- scenarios[d, 1:4]
  det.numeric <- det.levels[1:3]
  det.numeric[which(det.numeric == "L")] <- detect.l
  det.numeric[which(det.numeric == "M")] <- detect.m
  det.numeric[which(det.numeric== "H")] <- detect.h
  det.numeric[which(det.numeric== "NA")] <- NA
  if (is.na(det.numeric[2]) & is.na(det.numeric[3])) { # ABUNDANCE ONLY
    if (det.levels[4] == "L") {
      lowpopDat <- simData (indfates = lowpopTraj$indfates,
                            n.years = 15,
                            n.data.types = c(0.25,0.25,0.25),
                            ADonly = T,
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
      comb <- low.lam.combos[i,]
      lowout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("lowout-",i,"-",j,"-",d, sep = ""), lowout)
      saveRDS(lowout, here("results",paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(lowout)
    } else if (det.levels[4] == "M") {
      medpopDat <- simData (indfates = medpopTraj$indfates,
                            n.years = 15,
                            n.data.types = c(0.25,0.25,0.25),
                            ADonly = T,
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
      comb <- med.lam.combos[i,]
      medout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("medout-",i,"-",j,"-",d, sep = ""), medout)
      saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(medout)
    } else if (det.levels[4] == "H") {
      highpopDat <- simData (indfates = highpopTraj$indfates,
                             n.years = 15,
                             n.data.types = c(0.25,0.25,0.25),
                             ADonly = T,
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
      comb <- high.lam.combos[i,]
      highout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("highout-",i,"-",j,"-",d, sep = ""), highout)
      saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(highout)
    }
  } else if (is.na(det.numeric[2])) { # NO MARK RECAPTURE
    if (det.levels[4] == "L") {
      lowpopDat <- simData (indfates = lowpopTraj$indfates,
                            n.years = 15,
                            n.data.types = c(0.25,0.25,0.25),
                            ADonly = T,
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
      comb <- low.lam.combos[i,]
      lowout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("lowout-",i,"-",j,"-",d, sep = ""), lowout)
      saveRDS(lowout, here("results", paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(lowout)
    } else if (det.levels[4] == "M") {
      medpopDat <- simData (indfates = medpopTraj$indfates,
                            n.years = 15,
                            n.data.types = c(0.25,0.25,0.25),
                            ADonly = T,
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
      comb <- med.lam.combos[i,]
      medout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("medout-",i,"-",j,"-",d, sep = ""), medout)
      saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(medout)
    } else if (det.levels[4] == "H") {
      highpopDat <- simData (indfates = highpopTraj$indfates,
                             n.years = 15,
                             n.data.types = c(0.25,0.25,0.25),
                             ADonly = T,
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
      comb <- high.lam.combos[i,]
      highout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("highout-",i,"-",j,"-",d, sep = ""), highout)
      saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(highout)
    }
  } else if (is.na(det.numeric[3])) { # NO NEST SURVIVAL
    if (det.levels[4] == "L") {
      lowpopDat <- simData (indfates = lowpopTraj$indfates,
                            n.years = 15,
                            n.data.types = c(0.25,0.25,0.25),
                            ADonly = T,
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
      comb <- low.lam.combos[i,]
      lowout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("lowout-",i,"-",j,"-",d, sep = ""), lowout)
      saveRDS(lowout, here("results", paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(lowout)
    } else if (det.levels[4] == "M") {
      medpopDat <- simData (indfates = medpopTraj$indfates,
                            n.years = 15,
                            n.data.types = c(0.25,0.25,0.25),
                            ADonly = T,
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
      comb <- med.lam.combos[i,]
      medout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("medout-",i,"-",j,"-",d, sep = ""), medout)
      saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(medout)
    } else if (det.levels[4] == "H") {
      highpopDat <- simData (indfates = highpopTraj$indfates,
                             n.years = 15,
                             n.data.types = c(0.25,0.25,0.25),
                             ADonly = T,
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
      comb <- high.lam.combos[i,]
      highout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("highout-",i,"-",j,"-",d, sep = ""), highout)
      saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(highout)
    }
  } else { # FULL IPM
    if (det.levels[4] == "L") {
      lowpopDat <- simData (indfates = lowpopTraj$indfates,
                            n.years = 15,
                            n.data.types = c(0.25,0.25,0.25),
                            ADonly = T,
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
      comb <- low.lam.combos[i,]
      lowout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("lowout-",i,"-",j,"-",d, sep = ""), lowout)
      saveRDS(lowout, here("results", paste("lowout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(lowout)
    } else if (det.levels[4] == "M") {
      medpopDat <- simData (indfates = medpopTraj$indfates,
                            n.years = 15,
                            n.data.types = c(0.25,0.25,0.25),
                            ADonly = T,
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
      comb <- med.lam.combos[i,]
      medout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("medout-",i,"-",j,"-",d, sep = ""), medout)
      saveRDS(medout, here("results", paste("medout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(medout)
    } else if (det.levels[4] == "H") {
      highpopDat <- simData (indfates = highpopTraj$indfates,
                             n.years = 15,
                             n.data.types = c(0.25,0.25,0.25),
                             ADonly = T,
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
      comb <- high.lam.combos[i,]
      highout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = as.numeric(det.numeric))
      #assign(paste("highout-",i,"-",j,"-",d, sep = ""), highout)
      saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
      rm(highout)
    }
  } # else
  # assign(paste("highout-",i,"-",j,"-",d, sep = ""), highout)
  # saveRDS(highout, paste("highout-",i,"-",j,"-",d,".RDS", sep = ""))
  # rm(highout)
  # assign(paste("medout-",i,"-",j,"-",d, sep = ""), medout)
  # saveRDS(medout, paste("medout-",i,"-",j,"-",d,".RDS", sep = ""))
  # rm(medout)
  # assign(paste("lowout-",i,"-",j,"-",d, sep = ""), lowout)
  # saveRDS(lowout, paste("lowout-",i,"-",j,"-",d,".RDS", sep = ""))
  # rm(lowout)
  
  # saveRDS(highout, paste("highout-",i,"-",j,"-",d,".RDS", sep = ""))
  # saveRDS(medout, paste("medout-",i,"-",j,"-",d,".RDS", sep = ""))
  # saveRDS(lowout, paste("lowout-",i,"-",j,"-",d,".RDS", sep = ""))
  
  #rm(list = c("highout", "medout", "lowout"))
} # foreach - scenarios picked (i)
stopCluster(cl)