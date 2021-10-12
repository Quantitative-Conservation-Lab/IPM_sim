library(here)

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))

# simulations per scenario
sims.per <- 25 
scenarios.picked <- 25 # there are a bajillion of each scenario viable, so sample some
cores.per.computer <- 24 # cores at a time
computers <- 3 
n.instances <- cores.per.computer * computers

total.sims <- sims.per * scenarios.picked * dim(scenarios)[1]
total.sims

# let's say length per dat sim + model run is x minutes
time.per <- 10 # on average

total.time <- time.per * total.sims / n.instances/ 60 / 24  # in days,
# NOTE this assumes total use of the computers
total.time





