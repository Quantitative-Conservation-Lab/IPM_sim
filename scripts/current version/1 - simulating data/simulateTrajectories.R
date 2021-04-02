library(here)
library(foreach)
library(doParallel)

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
low.lam.combos <- readRDS(here("data", "low.lam.combos.RDS"))
med.lam.combos <- readRDS(here("data", "med.lam.combos.RDS"))
high.lam.combos <- readRDS(here("data", "high.lam.combos.RDS"))

# functions
source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_2.0function.R"))
source(here("scripts", "current version",
            "0 - preparing scenarios", "compute_time_calc.R"))

set.seed(1234)
low.rows <- sample(nrow(low.lam.combos), scenarios.picked, replace = FALSE)
med.rows <- sample(nrow(med.lam.combos), scenarios.picked, replace = FALSE)
high.rows <- sample(nrow(high.lam.combos), scenarios.picked, replace = FALSE)

low.lam.params <- low.lam.combos[low.rows, ]
med.lam.params <- med.lam.combos[med.rows, ]
high.lam.params <- high.lam.combos[high.rows, ]

# simulate populations
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2, setup_strategy = "sequential") #not to overload your computer
registerDoParallel(cl)
foreach(i = 1:scenarios.picked) %dopar% { #scenarios picked
  library(here)

#for (i in 1:scenarios.picked) {

  for (j in 1:sims.per) { # sims per
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
  } # sims per

#} # scenarios picked

} # foreach

stopCluster(cl)

