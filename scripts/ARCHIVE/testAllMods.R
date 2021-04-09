library(here)
library(nimble)

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

# TODO #####
# Fix this section
assign("test8", outIPM)

saveRDS(test7, "test7.RDS")
saveRDS(test8, "test8.RDS")
#######

for (d in 1:3) {
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
  medpopDat <- simData (indfates = lowpopTraj$indfates,
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
  highpopDat <- simData (indfates = lowpopTraj$indfates,
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
    if(m == 1) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    } else if (m == 2) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    } else if (m == 3) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    } else if (m == 4) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    }
  }
}
