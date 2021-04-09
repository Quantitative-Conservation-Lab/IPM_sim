# aeb
# april 2, 2020

# OUTLINE ####

# can we write this more efficiently??

# if m-array speeds things significantly
  # add more iterations to all
  # create a flag for things that didn't converge -
  # could just rerun these, or restart chains

# HAS idea
# add more iterations to models that we think should converge slowest
  # e.g. fewer datasets and low detection

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

# functions
source(here("scripts", "current version",
            "0 - preparing scenarios", "compute_time_calc.R"))
source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_2.0function.R"))
source(here("scripts", "current version",
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "4 - run models", "run_scenarios_helperFns.R"))

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
cl <- makeCluster(10, setup_strategy = "sequential") #not to overload your computer
registerDoParallel(cl)

foreach(i = 1:length(which.prio.1)) %dopar% { #scenarios picked
  library(here)
  library(nimble)
  for (j in 1:sims.per) {
    lowpopTraj <- readRDS(here("data", "lowTrajectories", paste("lowpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    medpopTraj <- readRDS(here("data", "medTrajectories", paste("medpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    highpopTraj <- readRDS(here("data", "highTrajectories", paste("highpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    for (d in 1:1){#nrow(scenarios)) {
      det.levels <- scenarios[d, 1:4]
      det.numeric <- det.levels[1:3]
      det.numeric[which(det.numeric == "L")] <- detect.l
      det.numeric[which(det.numeric == "M")] <- detect.m
      det.numeric[which(det.numeric== "H")] <- detect.h
      if (is.na(det.levels[2]) & is.na(det.levels[3])) { # ABUNDANCE ONLY
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
          lowout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
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
          medout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
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
          highout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
          #assign(paste("highout-",i,"-",j,"-",d, sep = ""), highout)
          saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(highout)
        }
      } else if (is.na(det.levels[2])) { # NO MARK RECAPTURE
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
          lowout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
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
          medout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
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
          highout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
          #assign(paste("highout-",i,"-",j,"-",d, sep = ""), highout)
          saveRDS(highout, here("results", paste("highout-",i,"-",j,"-",d,".RDS", sep = "")))
          rm(highout)
        }
      } else if (is.na(det.levels[3])) { # NO NEST SURVIVAL
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
          lowout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
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
          medout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
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
          highout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
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
          lowout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
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
          medout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
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
          highout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
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
    } # scenarios row (d)
  } # sims per (j)
} # foreach - scenarios picked (i)
stopCluster(cl)


outarray<-array(dim=c(10,19,6))
gl<-numeric(10)
#colnames(outarray)<-c("mean","2.5","25","50","75","97.5")
for(i in 1:10){
  outarray[i,,1]<-summary(readRDS(here(paste("lowout-1-",i,"-1.RDS", sep=""))))[[1]][,1]
  outarray[i,,2:6]<-summary(readRDS(here(paste("lowout-1-",i,"-1.RDS", sep=""))))[[2]]
  gl[i]<-gelman.diag(readRDS(here(paste("lowout-1-",i,"-1.RDS", sep=""))))[[2]]
}
