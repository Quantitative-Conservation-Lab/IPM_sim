
# standardize language
# scenario - parameter combination that we hope to recover
# observation structure/version - which datasets are included and detection levels
# simulation - unique model run

# AEB note - this is no longer correct but not fixing right now
# generate 25 scenarios per lambda
# for every scenario generate 50 populations
# then apply each detection level to generate 3x3x3= 27 datasets representing each
# level of detection for abundance, productivity, and mark-resight
# then - do like leave one out approach fitting 4 different models to the same dataset
# 27 * 4 = 108 versions per level of lambda
# see lucidchart figure for this
# assuming this
# generate 1250 populations
# 1250 * 27 datasets
# fit 4 models to each of those (135000)

library(tidyverse)
library(stringr)

# FUNCTION FROM HANNAH CODE
getNviable <- function() {
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


det.abund <- factor(x = c("L", "M", "H"))
det.MR <- factor(x = c("L", "M", "H", "NA"))
det.prod <- factor(x = c("L", "M", "H", "NA"))

lambda <- factor(x = c("L", "M", "H"))

# realistic demographic parameters
fec <- seq(1, 6, 0.05)/2
surv <- seq(0.1, 0.98, 0.01)

scenarios <- expand.grid(det.abund = det.abund, det.MR = det.MR, det.prod = det.prod,
                         lambda = lambda)

priority_score <- apply(scenarios[, 1:3], 1, function(x) sum(str_detect(x, "[HL]")))

lams <- getNviable() # this is kind of slow

lamlowcombs <- lams %>% filter(lams >= 0.94 & lams <= 0.96)
lammedcombs <- lams %>% filter(lams >= 0.99 & lams <= 1.01)
lamhighcombs <- lams %>% filter(lams >= 1.04 & lams <= 1.06)

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
saveRDS(scenarios, here("data", "scenarios.RDS"))
saveRDS(lamlowcombs, here("data", "low.lam.combos.RDS"))
saveRDS(lammedcombs, here("data", "med.lam.combos.RDS"))
saveRDS(lamhighcombs, here("data", "high.lam.combos.RDS"))


