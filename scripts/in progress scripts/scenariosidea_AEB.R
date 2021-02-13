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

# simulations per scenario
sims.per <- 100 # 50 for euring, 100 for manuscript

incl.abund <- factor(x = c("Y"))
incl.MR <- factor(x = c("Y", "N"))
incl.prod <- factor(x = c("Y", "N"))

# TODO
# for these ones, draw detection values from
# a standard normal
det.abund <- factor(x = c("L", "M", "H"))
det.MR <- factor(x = c("L", "M", "H"))
det.prod <- factor(x = c("L", "M", "H"))

# can play around with this
# l - 0.98 pm 0.01
# m - 1 pm 0.01
# h - 1.02 pm 0.01
lambda <- factor(x = c("L", "M", "H"))

# Fec levels
# can play around with this
fec <- seq(1, 3, 0.05)/2

# surviv levels
# can play around with this
surv <- seq(0.1, 0.9, 0.01)

scenarios <- expand.grid(sims.per = sims.per, 
                         incl.abund = incl.abund, incl.MR = incl.MR, incl.prod = incl.prod, 
                         det.abund = det.abund, det.MR = det.MR, det.prod = det.prod, 
                         lambda = lambda)

priority_score <- apply(scenarios[, 5:8], 1, function(x) sum(str_detect(x, "[HL]")))

lams <- getNviable() # this is kind of slow

lamlowcombs <- lams %>% filter(lams >= 0.96 & lams <= 0.98)
lammedcombs <- lams %>% filter(lams >= 0.99 & lams <= 1.01)
lamhighcombs <- lams %>% filter(lams >= 1.02 & lams <= 1.04)

n.viable.low <- dim(lamlowcombs)[1]
n.viable.med <- dim(lammedcombs)[1]
n.viable.high <- dim(lamhighcombs)[1]

scenarios.picked <- 10 # there are a bajillion of each scenario viable, so sample some

scenarios <- scenarios %>% 
  mutate(n.viable.combinations = case_when(
    lambda == "L" ~ n.viable.low,
    lambda == "M" ~ n.viable.med,
    lambda == "H" ~ n.viable.high
    )
  )

scenarios <- cbind(scenarios, 
                   priority_score, 
                   scenarios.picked)

# let's say we can run each incl.X * det.X combination on a separate r instance - 108 instances
# would require multiple machines - perhaps loon + delphine?
n.instances <- dim(expand.grid(incl.abund = incl.abund, incl.MR = incl.MR, incl.prod = incl.prod, 
                               det.abund = det.abund, det.MR = det.MR, det.prod = det.prod))[1]

total.sims <- sims.per * dim(scenarios)[1] * scenarios.picked
total.sims

# let's say length per dat sim + model run is 5 mins
time.per <- 10 # on average

# what is reasonable here...
# how should we cap this for now
total.time <- time.per * total.sims / n.instances/ 60 / 24  # in days, 
# NOTE this assumes total use of the computers
total.time

# STUFF TO SAVE
saveRDS(scenarios, "scenarios.RDS")
saveRDS(lamlowcombs, "low.lam.combos.RDS")
saveRDS(lammedcombs, "med.lam.combos.RDS")
saveRDS(lamhighcombs, "high.lam.combos.RDS")







