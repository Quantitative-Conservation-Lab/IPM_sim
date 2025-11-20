library(here)


# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
data_scenarios <- readRDS(here("data", "data_scenarios.RDS"))

# simulations per scenario
sims.per <- 50
scenarios.picked <- nrow(scenarios) # there are a bajillion of each scenario viable, so sample some
data.scenarios.picked <- nrow(data_scenarios)
cores.per.computer <- 24 # cores at a time
computers <- 1 
n.instances <- cores.per.computer * computers

total.sims <- sims.per * scenarios.picked * data.scenarios.picked
total.sims

# let's say length per dat sim + model run is x minutes
time.per <- 10 # on average

total.time <- time.per * total.sims / n.instances/ 60 / 24  # in days,
total.time
# 1 week on loon, just set and forget? over thanksgiving break
# set this as a goal




