library(here)

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
low.lam.combos <- readRDS(here("data", "low.lam.combos.RDS"))
med.lam.combos <- readRDS(here("data", "med.lam.combos.RDS"))
high.lam.combos <- readRDS(here("data", "high.lam.combos.RDS"))

# functions
source(here("scripts", "in progress scripts", "simHelperFns.R"))
source(here("scripts", "in progress scripts", "IPM_sim_2.0function.R"))

# df<-IPMSimFunction(n.years=10, n.data.types=c(0.25,0.25,0.25), 
#                    age.init=c(150,150), phi.1=0.3, phi.ad=0.3, f=0.5, max.nest.age=NA,
#                    mean.clutch.size=NA, phi.nest=NA, ADonly=T,p.1=NA,p.ad=0.8,
#                    BinMod=T,n.sam=3,p.count=0.55,sig=NA,
#                    productivity=T,p.prod=0.65)

# PLAN
# for each scenario
# each simulation should have data simulated from a randomly selected viable combination
# should have detection parameters drawn from a distribution