# aeb
# march 2, 2020
# example post processing
library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)
library(beepr)

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
low.lam.combos <- readRDS(here("data", "low.lam.combos.RDS"))
med.lam.combos <- readRDS(here("data", "med.lam.combos.RDS"))
high.lam.combos <- readRDS(here("data", "high.lam.combos.RDS"))

# functions
source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_2.0function.R"))
source(here("scripts", "current version",
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "4 - run models", "run_scenarios_helperFns.R"))

low.comb <- low.lam.combos[sample(1:5000, 1), 1:3]
med.comb <- med.lam.combos[sample(1:5000, 1), 1:3]
high.comb <- high.lam.combos[sample(1:5000, 1), 1:3]

lowpopTraj <- simPopTrajectory(n.years=15,
                               n.data.types=c(0.25,0.25,0.25),
                               age.init=c(150,150),
                               phi.1=as.numeric(low.comb[2]),
                               phi.ad=as.numeric(low.comb[3]),
                               f=as.numeric(low.comb[1]))

medpopTraj <- simPopTrajectory(n.years=15,
                               n.data.types=c(0.25,0.25,0.25),
                               age.init=c(150,150),
                               phi.1=as.numeric(med.comb[2]),
                               phi.ad=as.numeric(med.comb[3]),
                               f=as.numeric(med.comb[1]))

highpopTraj <- simPopTrajectory(n.years=15,
                                n.data.types=c(0.25,0.25,0.25),
                                age.init=c(150,150),
                                phi.1=as.numeric(high.comb[2]),
                                phi.ad=as.numeric(high.comb[3]),
                                f=as.numeric(high.comb[1]))

# simulate data
detect.l <- 0.3
detect.m <- 0.5
detect.h <- 0.8

detect <- c(detect.l, detect.m, detect.h)

nb <- 100000#0 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

beep(sound = 8)
cores=detectCores()
cl <- makeCluster(3, setup_strategy = "sequential") #not to overload your computer
registerDoParallel(cl)
foreach(d = 1:3) %dopar% { # detection level
#for (d in 1:3) {
  library(nimble)
  print(paste("detection level ", detect[d]))
  lowpopDat <- simData (indfates = lowpopTraj$indfates,
                        n.years = 15,
                        n.data.types = c(0.25,0.25,0.25),
                        ADonly = T,
                        p.1 = detect[d],
                        p.ad = detect[d],
                        p.count = detect[d],
                        p.prod = detect[d],
                        BinMod = T,
                        n.sam = 3,
                        sig = 0,
                        productivity = T)
  medpopDat <- simData (indfates = medpopTraj$indfates,
                        n.years = 15,
                        n.data.types = c(0.25,0.25,0.25),
                        ADonly = T,
                        p.1 = detect[d],
                        p.ad = detect[d],
                        p.count = detect[d],
                        p.prod = detect[d],
                        BinMod = T,
                        n.sam = 3,
                        sig = 0,
                        productivity = T)
  highpopDat <- simData (indfates = highpopTraj$indfates,
                         n.years = 15,
                         n.data.types = c(0.25,0.25,0.25),
                         ADonly = T,
                         p.1 = detect[d],
                         p.ad = detect[d],
                         p.count = detect[d],
                         p.prod = detect[d],
                         BinMod = T,
                         n.sam = 3,
                         sig = 0,
                         productivity = T)
  for (m in 1:4) {
    print(paste("model ", m))
    if(m == 1) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      highout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      medout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      lowout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    } else if (m == 2) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      highout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      medout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      lowout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    } else if (m == 3) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      highout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      medout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      lowout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    } else if (m == 4) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      highout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      medout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      lowout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    }
    assign(paste("highout-",d,"-",m, sep = ""), highout)
    assign(paste("medout-",d,"-",m, sep = ""), medout)
    assign(paste("lowout-",d,"-",m, sep = ""), lowout)

    saveRDS(highout, paste("highout-",d,"-",m, ".RDS", sep = ""))
    saveRDS(medout, paste("medout-",d,"-",m,  ".RDS", sep = ""))
    saveRDS(lowout, paste("lowout-",d,"-",m,  ".RDS", sep = ""))

    rm(list = c("highout", "medout", "lowout"))
  }
}
stopCluster(cl)
beep(sound = 8)

# NOW PLAY WITH THE OUTPUT

high.comb
summary(`highout-1-1`)
summary(`highout-1-2`)
summary(`highout-1-3`)
summary(`highout-1-4`)
summary(`highout-1-1`)
summary(`highout-2-1`)
summary(`highout-2-2`)
summary(`highout-2-3`)
summary(`highout-2-4`)
summary(`highout-3-1`)
summary(`highout-3-2`)
summary(`highout-3-3`)
summary(`highout-3-4`)
