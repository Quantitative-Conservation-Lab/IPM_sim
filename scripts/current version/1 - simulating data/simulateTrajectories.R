# TODO
# parallelize data

library(here)

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
low.lam.combos <- readRDS(here("data", "low.lam.combos.RDS"))
med.lam.combos <- readRDS(here("data", "med.lam.combos.RDS"))
high.lam.combos <- readRDS(here("data", "high.lam.combos.RDS"))

# functions
source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_2.0function.R"))
source(here("scripts", "current version",
            "O - preparing scenarios", "compute_time_calc.R"))

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

