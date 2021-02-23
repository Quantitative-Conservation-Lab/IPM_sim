library(here)

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))

# simulations per scenario
sims.per <- 50 # 50 for euring, 100 for manuscript
scenarios.picked <- 25 # there are a bajillion of each scenario viable, so sample some
cores.per.computer <- 24 # cores at a time
computers <- 3 # loon, delphine, ursus
n.instances <- cores.per.computer * computers

total.sims <- sims.per * scenarios.picked * dim(scenarios)[1]
total.sims

# let's say length per dat sim + model run is 5 mins
time.per <- 10 # on average

# what is reasonable here...
# how should we cap this for now
total.time <- time.per * total.sims / n.instances/ 60 / 24  # in days, 
# NOTE this assumes total use of the computers
total.time





