# TODO
# parallelize data

library(here)

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
low.lam.combos <- readRDS(here("data", "low.lam.combos.RDS"))
med.lam.combos <- readRDS(here("data", "med.lam.combos.RDS"))
high.lam.combos <- readRDS(here("data", "high.lam.combos.RDS"))

# functions
source(here("scripts", "in progress scripts", "simHelperFns.R"))
source(here("scripts", "in progress scripts", "IPM_sim_2.0function.R"))
source(here("scripts", "in progress scripts", "compute_time_calc.R"))

# DRAW LAMBDA COMBINATIONS
# PLOT THEM
# SAVE THEM
#low.lam.params <- sample_n(low.lam.combos, scenarios.picked)
low.lam.params <- readRDS("low.lam.params.RDS")

par(mfrow = c(3, 1))
plot(low.lam.params$fec, low.lam.params$phi1)
plot(low.lam.params$fec, low.lam.params$phiad)
plot(low.lam.params$phi1, low.lam.params$phiad)

#saveRDS(low.lam.params, "low.lam.params.RDS")

#med.lam.params <- sample_n(med.lam.combos, scenarios.picked)
med.lam.params <- readRDS("med.lam.params.RDS")

par(mfrow = c(3, 1))
plot(med.lam.params$fec, med.lam.params$phi1)
plot(med.lam.params$fec, med.lam.params$phiad)
plot(med.lam.params$phi1, med.lam.params$phiad)

#saveRDS(med.lam.params, "med.lam.params.RDS")

#high.lam.params <- sample_n(high.lam.combos, scenarios.picked)
high.lam.params <- readRDS("high.lam.params.RDS")

par(mfrow = c(3, 1))
plot(high.lam.params$fec, high.lam.params$phi1)
plot(high.lam.params$fec, high.lam.params$phiad)
plot(high.lam.params$phi1, high.lam.params$phiad)

#saveRDS(med.lam.params, "med.lam.params.RDS")

# SET DETECTION PARAMETERS

sigma.detect <- 0.1
detect.l <- 0.3
detect.m <- 0.5
detect.h <- 0.8

# TODO
# split simulation function into two parts
# done, but now need to test

# simulate populations

for (i in 1:scenarios.picked) {
  
  for (j in 1:sims.per) {
    lowpopTraj <- simPopTrajectory(n.years=15, 
                                   n.data.types=c(0.25,0.25,0.25),
                                   age.init=c(150,150), 
                                   phi.1=low.lam.params$phi1[i], 
                                   phi.ad=low.lam.params$phiad[i], 
                                   f=low.lam.params$fec[i])
    
    medpopTraj <- simPopTrajectory(n.years=15, 
                                   n.data.types=c(0.25,0.25,0.25),
                                   age.init=c(150,150), 
                                   phi.1=med.lam.params$phi1[i], 
                                   phi.ad=med.lam.params$phiad[i], 
                                   f=med.lam.params$fec[i])
    
    highpopTraj <- simPopTrajectory(n.years=15, 
                                    n.data.types=c(0.25,0.25,0.25),
                                    age.init=c(150,150), 
                                    phi.1=high.lam.params$phi1[i], 
                                    phi.ad=high.lam.params$phiad[i], 
                                    f=high.lam.params$fec[i])
    
    assign(paste("lowpopTraj", "-", i, "-", j, sep = ""), lowpopTraj)
    saveRDS(lowpopTraj, here("data", "lowTrajectories", paste("lowpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    rm(lowpopTraj)
    
    assign(paste("medpopTraj", "-", i, "-", j, sep = ""), medpopTraj)
    saveRDS(medpopTraj, here("data", "medTrajectories", paste("medpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    rm(medpopTraj)
    
    assign(paste("highpopTraj", "-", i, "-", j, sep = ""), highpopTraj)
    saveRDS(highpopTraj, here("data", "highTrajectories", paste("highpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    rm(highpopTraj)  
  }
  
}

# simulate observations 
# simulate data on the fly

