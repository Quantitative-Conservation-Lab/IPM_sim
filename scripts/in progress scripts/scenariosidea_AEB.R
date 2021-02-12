# GOAL
# use april and may to run whatever hannah presents at euring

# simulations per scenario
sims.per <- 50 # 50 for euring, 100 for manuscript

incl.abund <- factor(x = c("Y"))
incl.MR <- factor(x = c("Y", "N"))
incl.prod <- factor(x = c("Y", "N"))

det.abund <- factor(x = c("L", "M", "H"))
det.MR <- factor(x = c("L", "M", "H"))
det.prod <- factor(x = c("L", "M", "H"))

# can play around with this
# l - 0.98 pm 0.01
# m - 1 pm 0.01
# h - 1.02 pm 0.01
lambda <- factor(x = c("L", "M", "H"))

# pick h/m/l levels for each param
# figure out how many combinations work for h/m/l lambda

n.viable.combinations <- 10 #e.g. 5 on average

# Fec levels
# l - 
# m -
# h - 

# surviv levels
# l - 
# m - 
# h - 

#demog.abund <- factor(x = c("L", "M", "H"))
#demog.MR <- factor(x = c("L", "M", "H"))
#demog.prod <- factor(x = c("L", "M", "H"))

scenarios <- cbind(expand.grid(sims.per = sims.per, 
                         incl.abund = incl.abund, incl.MR = incl.MR, incl.prod = incl.prod, 
                         det.abund = det.abund, det.MR = det.MR, det.prod = det.prod, 
                         lambda = lambda
                         ), 
                   n.viable.combinations)

# which things to run first, for euring
# things with more highs and lows vs mediums
priority.score <- 

total.sims <- sims.per * n.viable.combinations * dim(scenarios)[1]

# let's say length per dat sim + model run is 5 mins
time.per <- 5 # on average

# let's say we can run each incl.X * det.X combination on a separate r instance - 108 instances
# would require multiple machines - perhaps loon + delphine?
n.instances <- dim(expand.grid(incl.abund = incl.abund, incl.MR = incl.MR, incl.prod = incl.prod, 
                               det.abund = det.abund, det.MR = det.MR, det.prod = det.prod))[1]

# what is reasonable here...
# how should we cap this for now
total.time <- time.per * total.sims / n.instances/ 60 / 24  # in days, 
# NOTE this assumes total use of the computers

# need to know
# ram.per.inst
# cpu.per.inst

# TODO
# write code in parallel so we don't have to stop and start 100 instances at a time
# potentially a lot of big files - do we need all? 
  # could we do preprocessing and save the summary
  # create flags for sims that need attention
# test convergence on some scenarios to figure out what burnin iterations and thinning should be

# vs
# loon specs
# ursus specs

