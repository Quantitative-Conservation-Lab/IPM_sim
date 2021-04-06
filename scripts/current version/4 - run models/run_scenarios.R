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
low.lam.combos <- readRDS(here("data", "low.lam.params.RDS"))
med.lam.combos <- readRDS(here("data", "med.lam.params.RDS"))
high.lam.combos <- readRDS(here("data", "high.lam.params.RDS"))

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

nb <- 100000#0 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

cores=detectCores()
cl <- makeCluster(cores-2, setup_strategy = "sequential") #not to overload your computer
registerDoParallel(cl)

foreach(i = 1:scenarios.picked) %dopar% { #scenarios picked
  library(here)
  library(nimble)
  for (j in 1:sims.per) {
    lowpopTraj <- readRDS(here("data", "lowTrajectories", paste("lowpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    medpopTraj <- readRDS(here("data", "medTrajectories", paste("medpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    highpopTraj <- readRDS(here("data", "highTrajectories", paste("highpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    for (d in 1:nrow(scenarios)) {
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
                                p.1 = det.numeric[2], #
                                p.ad = det.numeric[2], #
                                p.count = det.numeric[1], #
                                p.prod = det.numeric[3], #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- lowpopDat
          popTraj <- lowpopTraj
          comb <- low.comb
          lowout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        } else if (det.levels[4] == "M") {
          medpopDat <- simData (indfates = medpopTraj$indfates,
                                n.years = 15,
                                n.data.types = c(0.25,0.25,0.25),
                                ADonly = T,
                                p.1 = det.numeric[2], #
                                p.ad = det.numeric[2], #
                                p.count = det.numeric[1], #
                                p.prod = det.numeric[3], #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- medpopDat
          popTraj <- medpopTraj
          comb <- med.comb
          medout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        } else if (det.levels[4] == "H") {
          highpopDat <- simData (indfates = highpopTraj$indfates,
                                 n.years = 15,
                                 n.data.types = c(0.25,0.25,0.25),
                                 ADonly = T,
                                 p.1 = det.numeric[2], #
                                 p.ad = det.numeric[2], #
                                 p.count = det.numeric[1], #
                                 p.prod = det.numeric[3], #
                                 BinMod = T,
                                 n.sam = 3,
                                 sig = 0,
                                 productivity = T)
          popDat <- highpopDat
          popTraj <- highpopTraj
          comb <- high.comb
          highout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        }
      } else if (is.na(det.levels[2])) { # NO MARK RECAPTURE
        if (det.levels[4] == "L") {
          lowpopDat <- simData (indfates = lowpopTraj$indfates,
                                n.years = 15,
                                n.data.types = c(0.25,0.25,0.25),
                                ADonly = T,
                                p.1 = det.numeric[2], #
                                p.ad = det.numeric[2], #
                                p.count = det.numeric[1], #
                                p.prod = det.numeric[3], #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- lowpopDat
          popTraj <- lowpopTraj
          comb <- low.comb
          lowout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        } else if (det.levels[4] == "M") {
          medpopDat <- simData (indfates = medpopTraj$indfates,
                                n.years = 15,
                                n.data.types = c(0.25,0.25,0.25),
                                ADonly = T,
                                p.1 = det.numeric[2], #
                                p.ad = det.numeric[2], #
                                p.count = det.numeric[1], #
                                p.prod = det.numeric[3], #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- medpopDat
          popTraj <- medpopTraj
          comb <- med.comb
          medout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        } else if (det.levels[4] == "H") {
          highpopDat <- simData (indfates = highpopTraj$indfates,
                                 n.years = 15,
                                 n.data.types = c(0.25,0.25,0.25),
                                 ADonly = T,
                                 p.1 = det.numeric[2], #
                                 p.ad = det.numeric[2], #
                                 p.count = det.numeric[1], #
                                 p.prod = det.numeric[3], #
                                 BinMod = T,
                                 n.sam = 3,
                                 sig = 0,
                                 productivity = T)
          popDat <- highpopDat
          popTraj <- highpopTraj
          comb <- high.comb
          highout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        }
      } else if (is.na(det.levels[3])) { # NO NEST SURVIVAL
        if (det.levels[4] == "L") {
          lowpopDat <- simData (indfates = lowpopTraj$indfates,
                                n.years = 15,
                                n.data.types = c(0.25,0.25,0.25),
                                ADonly = T,
                                p.1 = det.numeric[2], #
                                p.ad = det.numeric[2], #
                                p.count = det.numeric[1], #
                                p.prod = det.numeric[3], #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- lowpopDat
          popTraj <- lowpopTraj
          comb <- low.comb
          lowout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        } else if (det.levels[4] == "M") {
          medpopDat <- simData (indfates = medpopTraj$indfates,
                                n.years = 15,
                                n.data.types = c(0.25,0.25,0.25),
                                ADonly = T,
                                p.1 = det.numeric[2], #
                                p.ad = det.numeric[2], #
                                p.count = det.numeric[1], #
                                p.prod = det.numeric[3], #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- medpopDat
          popTraj <- medpopTraj
          comb <- med.comb
          medout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        } else if (det.levels[4] == "H") {
          highpopDat <- simData (indfates = highpopTraj$indfates,
                                 n.years = 15,
                                 n.data.types = c(0.25,0.25,0.25),
                                 ADonly = T,
                                 p.1 = det.numeric[2], #
                                 p.ad = det.numeric[2], #
                                 p.count = det.numeric[1], #
                                 p.prod = det.numeric[3], #
                                 BinMod = T,
                                 n.sam = 3,
                                 sig = 0,
                                 productivity = T)
          popDat <- highpopDat
          popTraj <- highpopTraj
          comb <- high.comb
          highout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        }
      } else { # FULL IPM
        if (det.levels[4] == "L") {
          lowpopDat <- simData (indfates = lowpopTraj$indfates,
                                n.years = 15,
                                n.data.types = c(0.25,0.25,0.25),
                                ADonly = T,
                                p.1 = det.numeric[2], #
                                p.ad = det.numeric[2], #
                                p.count = det.numeric[1], #
                                p.prod = det.numeric[3], #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- lowpopDat
          popTraj <- lowpopTraj
          comb <- low.comb
          lowout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        } else if (det.levels[4] == "M") {
          medpopDat <- simData (indfates = medpopTraj$indfates,
                                n.years = 15,
                                n.data.types = c(0.25,0.25,0.25),
                                ADonly = T,
                                p.1 = det.numeric[2], #
                                p.ad = det.numeric[2], #
                                p.count = det.numeric[1], #
                                p.prod = det.numeric[3], #
                                BinMod = T,
                                n.sam = 3,
                                sig = 0,
                                productivity = T)
          popDat <- medpopDat
          popTraj <- medpopTraj
          comb <- med.comb
          medout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        } else if (det.levels[4] == "H") {
          highpopDat <- simData (indfates = highpopTraj$indfates,
                                 n.years = 15,
                                 n.data.types = c(0.25,0.25,0.25),
                                 ADonly = T,
                                 p.1 = det.numeric[2], #
                                 p.ad = det.numeric[2], #
                                 p.count = det.numeric[1], #
                                 p.prod = det.numeric[3], #
                                 BinMod = T,
                                 n.sam = 3,
                                 sig = 0,
                                 productivity = T)
          popDat <- highpopDat
          popTraj <- highpopTraj
          comb <- high.comb
          highout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
        }
      } # else
      assign(paste("highout-",i,"-",j,"-",d, sep = ""), highout)
      assign(paste("medout-",i,"-",j,"-",d, sep = ""), medout)
      assign(paste("lowout-",i,"-",j,"-",d, sep = ""), lowout)

      saveRDS(highout, paste("highout-",i,"-",j,"-",d,".RDS", sep = ""))
      saveRDS(medout, paste("medout-",i,"-",j,"-",d,".RDS", sep = ""))
      saveRDS(lowout, paste("lowout-",i,"-",j,"-",d,".RDS", sep = ""))

      rm(list = c("highout", "medout", "lowout"))
    } # scenarios row (d)
  } # sims per (j)
} # foreach - scenarios picked (i)

