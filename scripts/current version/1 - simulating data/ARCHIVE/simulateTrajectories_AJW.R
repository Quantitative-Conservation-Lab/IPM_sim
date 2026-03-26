library(here)
library(foreach)
library(doParallel)
library(tidyverse)

# load demographic scenarios
scenarios <- readRDS(here("data", "demographic_scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename(
    "phi1" = "S.J", 
    "phiad" = "S.A",
    "fec" = "f"
  )

# load functions and namespace
# source(here("scripts", "current version",
#             "1 - simulating data", "IPM_sim_3.0function.R"))

low.lam.params <- scenarios %>% 
  filter(trend == "decline")
med.lam.params <- scenarios %>% 
  filter(trend == "stable")
high.lam.params <- scenarios %>% 
  filter(trend == "increase")

scenarios.picked <- nrow(low.lam.params) # TODO note change here
sims.per <- 100

cores=detectCores()
cl <- makeCluster(sims.per/2, setup_strategy = "sequential") 

t.start <- Sys.time()
registerDoParallel(cl)

foreach(j = 1:sims.per) %dopar% { 
  library(here)
  for (i in 1:scenarios.picked) { 
    
  lowpopTraj <- simPopTrajectory(n.years=15,
                                   age.init=c(150,150),
                                   phi.1=low.lam.params$phi1[i],
                                   phi.ad=low.lam.params$phiad[i],
                                   f=low.lam.params$fec[i])
  
  assign(paste("lowpopTraj", "-", i, "-", j, sep = ""), lowpopTraj)
  saveRDS(lowpopTraj, here("data", "lowTrajectories", paste("lowpopTraj", "-", i, "-", j, ".RDS", sep = "")))
  rm(lowpopTraj)

    medpopTraj <- simPopTrajectory(n.years=15,
                                   age.init=c(150,150),
                                   phi.1=med.lam.params$phi1[i],
                                   phi.ad=med.lam.params$phiad[i],
                                   f=med.lam.params$fec[i])
    
    assign(paste("medpopTraj", "-", i, "-", j, sep = ""), medpopTraj)
    saveRDS(medpopTraj, here("data", "medTrajectories", paste("medpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    rm(medpopTraj)

    highpopTraj <- simPopTrajectory(n.years=15,
                                    age.init=c(150,150),
                                    phi.1=high.lam.params$phi1[i],
                                    phi.ad=high.lam.params$phiad[i],
                                    f=high.lam.params$fec[i])

    assign(paste("highpopTraj", "-", i, "-", j, sep = ""), highpopTraj)
    saveRDS(highpopTraj, here("data", "highTrajectories", paste("highpopTraj", "-", i, "-", j, ".RDS", sep = "")))
    rm(highpopTraj)
  } # sims per

} # foreach

stopCluster(cl)
t.end <- Sys.time()
t.end - t.start

