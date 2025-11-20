library(tidyverse)
library(stringr)

getNviable <- function() {
  # 1 - figure out how many viable combos there are - then add them to the scenarios table
  # 2 - what ARE the viable combinations

  # parameter combination
  coms <- expand.grid(fec = fec, phi1 = surv, phiad = surv)
  #output vector
  lams<-numeric(length(coms[,1]))
  #leslie matrix for each
  mats<-array(dim=c(length(coms[,1]), 2,2))

  for(i in 1:length(coms[,1])){ #loop over combinations
    #leslie matrix for each combination
    mats[i,1,]<-c(coms[i,"phi1"]*coms[i,"fec"],coms[i,"phi1"]*coms[i,"fec"])
    mats[i,2,]<-c(coms[i,"phiad"], coms[i,"phiad"])
    #Get the eigenvalue for each leslie matrix, which is the lambda
    lams[i]<-eigen(mats[i,,])$values[1]
  }

  return(cbind(coms, lams))
}

# construct matrix of scenarios #####
det.abund <- factor(x = c("L", "M", "H"))
det.MR <- factor(x = c("L", "M", "H", "NA"))
det.prod <- factor(x = c("L", "M", "H", "NA"))
lambda <- factor(x = c("L", "M", "H"))
scenarios <- expand.grid(det.abund = det.abund, det.MR = det.MR, det.prod = det.prod,
                         lambda = lambda)

# run function to generate possible parameter combinations #####
fec <- seq(1, 6, 0.05)/2 # realistic demographic parameters
surv <- seq(0.1, 0.98, 0.01) # realistic demographic parameters
lams <- getNviable() # NOTE this is kind of slow

# filter into decreasing, stable, and increasing trajectories ####
lamlowcombs <- lams %>% filter(lams >= 0.94 & lams <= 0.96)
lammedcombs <- lams %>% filter(lams >= 0.99 & lams <= 1.01)
lamhighcombs <- lams %>% filter(lams >= 1.04 & lams <= 1.06)

# how many viable parameter combinations are there to choose from?
n.viable.low <- dim(lamlowcombs)[1]
n.viable.med <- dim(lammedcombs)[1]
n.viable.high <- dim(lamhighcombs)[1]

scenarios <- scenarios %>%
  mutate(n.viable.combinations = case_when(
    lambda == "L" ~ n.viable.low,
    lambda == "M" ~ n.viable.med,
    lambda == "H" ~ n.viable.high
    )
  )

# save objects 
saveRDS(scenarios, here("data", "scenarios.RDS"))
saveRDS(lamlowcombs, here("data", "low.lam.combos.RDS"))
saveRDS(lammedcombs, here("data", "med.lam.combos.RDS"))
saveRDS(lamhighcombs, here("data", "high.lam.combos.RDS"))


