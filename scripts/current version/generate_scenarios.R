# GOAL
# use april and may to run whatever hannah presents at euring

# PLAN
# for each scenario
# each simulation should have data simulated from a randomly selected viable combination
# should have detection parameters drawn from a distribution

# TODO
# write code in parallel so we don't have to stop and start 100 instances at a time
# potentially a lot of big files - do we need all? 
# could we do preprocessing and save the summary
# create flags for sims that need attention
# test convergence on some scenarios to figure out what burnin iterations and thinning should be

# standardize language
# scenario - parameter combination that we hope to recover
# observation structure/version - which datasets are included and detection levels
# simulation - unique model run

# generate 25 scenarios per lambda 
# for every scenario generate 50 populations
# then apply each detection level to generate 3x3x3= 27 datasets representing each
# level of detection for abundance, productivity, and mark-resight
# then - do like leave one out approach fitting 4 different models to the same dataset
# 27 * 4 = 108 versions per level of lambda 
# lucidchart figure for this

# assuming this
# generate 1250 populations
# 1250 * 27 datasets
# fit 4 models to each of those (135000)

library(tidyverse)
library(stringr)

# FUNCTION FROM HANNAH CODE
getNviable <- function() {
  # TODO
  # 1 - figure out how many viable combos there are - then add them to the scenarios table
  # 2 - what ARE the viable combinations
  
  # NOW THE EXCITING PART - how many valid combos
  #creat a matrix for all possible combinations of these 3
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

# TODO
# for these ones, draw detection values from
# a standard normal
det.abund <- factor(x = c("L", "M", "H"))
det.MR <- factor(x = c("L", "M", "H", "NA"))
det.prod <- factor(x = c("L", "M", "H", "NA"))

# can play around with this
# l - 0.98 pm 0.01
# m - 1 pm 0.01
# h - 1.02 pm 0.01
lambda <- factor(x = c("L", "M", "H"))

# TODO
# run best case scenarios with these trends to make sure model works well
# if can't capture well - increase those lambdas

# Fec levels
# can play around with this

# TODO
# up range to 6 for example
# what is rel between fecundity and adult survival and lambda - are we capturing that 2d space well
fec <- seq(1, 6, 0.05)/2

# surviv levels
# can play around with this

# TODO
# up this range too - 0.975
surv <- seq(0.1, 0.98, 0.01)

# TODO
# could bin these into HML and sample the lambda scenarios makig sure we represent the 2d space well

scenarios <- expand.grid(det.abund = det.abund, det.MR = det.MR, det.prod = det.prod, 
                         lambda = lambda)

priority_score <- apply(scenarios[, 1:3], 1, function(x) sum(str_detect(x, "[HL]")))

lams <- getNviable() # this is kind of slow

lamlowcombs <- lams %>% filter(lams >= 0.96 & lams <= 0.98)
lammedcombs <- lams %>% filter(lams >= 0.99 & lams <= 1.01)
lamhighcombs <- lams %>% filter(lams >= 1.02 & lams <= 1.04)

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

scenarios <- cbind(scenarios, 
                   scenarios.picked,
                   priority_score) %>% 
  arrange(desc(priority_score)) %>% 
  mutate(batch = ceiling(row_number() / 24))

# STUFF TO SAVE
saveRDS(scenarios, "scenarios.RDS")
saveRDS(lamlowcombs, "low.lam.combos.RDS")
saveRDS(lammedcombs, "med.lam.combos.RDS")
saveRDS(lamhighcombs, "high.lam.combos.RDS")


