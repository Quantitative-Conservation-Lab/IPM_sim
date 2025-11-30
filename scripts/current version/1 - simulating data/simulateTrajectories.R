library(here)
library(foreach)
library(doParallel)
library(tidyverse)

# load data
# TODO - this should maybe go to the end of the previous file
scenarios <- readRDS(here("data", "scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename(
    "phi1" = "S.J", 
    "phiad" = "S.A",
    "fec" = "f"
  )
# low.lam.combos <- readRDS(here("data", "low.lam.combos.RDS"))
# med.lam.combos <- readRDS(here("data", "med.lam.combos.RDS"))
# high.lam.combos <- readRDS(here("data", "high.lam.combos.RDS"))

# load functions and namespace
source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_3.0function.R"))
# source(here("scripts", "current version",
#             "0 - preparing scenarios", "compute_time_calc.R"))

# TODO - don't need scenarios picked because we are going to run all
# correct? yes
# pick param combinations to generate population trajectories from
# set.seed(1234)
# low.rows <- sample(nrow(low.lam.combos), scenarios.picked, replace = FALSE)
# med.rows <- sample(nrow(med.lam.combos), scenarios.picked, replace = FALSE)
# high.rows <- sample(nrow(high.lam.combos), scenarios.picked, replace = FALSE)
# low.lam.params <- low.lam.combos[low.rows, ]
# med.lam.params <- med.lam.combos[med.rows, ]
# high.lam.params <- high.lam.combos[high.rows, ]

# TODO
  # rename low lam params etc. get called below
low.lam.params <- scenarios %>% 
  filter(trend == "decline")
med.lam.params <- scenarios %>% 
  filter(trend == "stable")
high.lam.params <- scenarios %>% 
  filter(trend == "increase")

# TODO discussed creating a global constants file, decide if useful
scenarios.picked <- nrow(scenarios)
sims.per <- 50 # TODO - is this what we decided? didn't we discuss either less
# demog stochasticity or less sampling stochasticity?
# if we dropped demog stochasticity we'd have less to save which would save memory among other things

# TODO - is there an advantage to using furrr? or leave as is
# simulate populations
# and save the trajectories in data files
cores=detectCores()
cl <- makeCluster(cores[1]-2, setup_strategy = "sequential") 

# TODO remove
# TESTING 
sims.per <- 3
cl <- makeCluster(3, setup_strategy = "sequential")
# TESTING

registerDoParallel(cl)
foreach(i = 1:scenarios.picked) %dopar% { #scenarios picked
  library(here)
  for (j in 1:sims.per) { # sims per
    
  lowpopTraj <- simPopTrajectory(n.years=15,
                                   #n.data.types=c(0.25,0.25,0.25),
                                   age.init=c(150,150),
                                   phi.1=low.lam.params$phi1[i],
                                   phi.ad=low.lam.params$phiad[i],
                                   f=low.lam.params$fec[i])

    medpopTraj <- simPopTrajectory(n.years=15,
                                   #n.data.types=c(0.25,0.25,0.25),
                                   age.init=c(150,150),
                                   phi.1=med.lam.params$phi1[i],
                                   phi.ad=med.lam.params$phiad[i],
                                   f=med.lam.params$fec[i])

    highpopTraj <- simPopTrajectory(n.years=15,
                                    #n.data.types=c(0.25,0.25,0.25),
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

} # foreach

stopCluster(cl)

