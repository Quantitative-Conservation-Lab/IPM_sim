# AEB 
# Nov 6

# jotting some stuff down for data sim

# TODO 
# interp 1s inbetween consecutive visits where nest was active
# make sure that a basic nest survival model can resolve these parameters
# come up with inits

# load libraries
library(nimble)

# parameters ######

# set season length and stuff
n.initiation.dates <- 120
max.nest.age <- 30
first.initiation.date <- 1
last.fledge.date <- n.initiation.dates + max.nest.age
# AEB note - 2 extra days for monitoring
season.length <- last.fledge.date - first.initiation.date + 1 + 2

# nests for population
N.nests.total <- 100
prop.nests.found <- 0.8

# mean clutch size
mean.clutch.size <- 5
# AEB note - additional assumption that all eggs hatch and all nestlings fledge

# daily nest survival
phi <- 0.975

# true fecundity (female chicks produced per nest)
# assuming 50/50 sex ratio at hatching
true.fec <- 1/2 * mean.clutch.size * phi^max.nest.age

# observation params
visit.interval <- 3

# start data sim ####

# all nests - this is latent
total.nests.age <- matrix(NA, nrow = N.nests.total, ncol = season.length)
total.nests.status <- matrix(NA, nrow = N.nests.total, ncol = season.length)

# initiation date
# uniform across season length
# sort by initiation date
init.dates <- sort(rcat(N.nests.total, rep(1/n.initiation.dates, length.out = n.initiation.dates)))

for (i in 1:N.nests.total) {
  # true age of nest 
  total.nests.age[i, init.dates[i]:(init.dates[i] + max.nest.age - 1)] <- 1:max.nest.age
  
  # true status of nest
  total.nests.status[i, init.dates[i]] <- 1 
  for (a in 1:(max.nest.age-1)) {
    total.nests.status[i, init.dates[i] + a] <- rbinom(1, 1, total.nests.status[i, init.dates[i] + a - 1] * phi)
  }
}

# which of these nests were found - this is also latent
which.found <- sort(sample(1:N.nests.total, prop.nests.found * N.nests.total))
init.dates.found <- init.dates[which.found]
N.nests.found <- length(which.found)

found.nests.age <- total.nests.age[which.found, ]
found.nests.status <- total.nests.status[which.found, ]

# detection process

# age at first detection
# nests are checked every 3 days after that
# AEB note - currently equally likely to be observed at any age
age.when.found <- rcat(N.nests.found, rep(1/max.nest.age, length.out = max.nest.age))
ages.visited <- list(NULL)
for (i in 1:N.nests.found) {
  ages.visited[[i]] <- unique(c(seq(age.when.found[i], max.nest.age, by = visit.interval), max.nest.age))
}

# TODO
# if nest failed before found
# it's not actually a nest for our purp

# create observation matrix
observed.nest.age <- matrix(NA, nrow = N.nests.found, ncol = season.length)
observed.nest.status <- matrix(NA, nrow = N.nests.found, ncol = season.length)

for (i in 1:N.nests.found) {
  # age of nest
  observed.nest.age[i, init.dates.found[i]+ages.visited[[i]]-1] <- ages.visited[[i]]
    
  # status of nest
  observed.nest.status[i, init.dates.found[i]+ages.visited[[i]]-1] <- found.nests.status[i, init.dates.found[i]+ages.visited[[i]]-1]

  # AEB note - we assume we can age the nest
  # and that our observations are perfect
  # so if it fledged we know the date it fledged or *should* have fledged
}

# HAVE TO REMOVE NESTS THAT FAILED BEFORE FOUND

is.nest <- apply(observed.nest.status, 1, function(x) x[min(which(!is.na(x)))] == 1)
  
observed.nest.age <- observed.nest.age[is.nest, ]
observed.nest.status <- observed.nest.status[is.nest, ]

N.nests.found <- dim(observed.nest.age)[1]

# how many were successful
# is last status a 1 or a 0
last.status <- apply(observed.nest.status, 1, function(x) x[max(which(!is.na(x)))])

which.successful <- which(last.status == 1)

N.nests.successful <- length(which.successful)

first <- apply(observed.nest.status, 1, function(x) min(which(!is.na(x))))

# if nest fledged, last is fledge day
# if nest failed, last is first zero
last <- rep(NA, length.out = N.nests.found)
last[which.successful] <- apply(observed.nest.status[which.successful, ], 1, function(x) max(which(!is.na(x))))
last[-which.successful] <- apply(observed.nest.status[-which.successful, ], 1, function(x) min(which(!is.na(x) & x == 0)))

# clutch size
clutch.sizes <- rpois(N.nests.successful, mean.clutch.size)


