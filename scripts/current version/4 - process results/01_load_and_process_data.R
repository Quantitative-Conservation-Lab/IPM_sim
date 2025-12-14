library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)
library(coda)
library(nlist)
library(beepr)

# TODO
# needs reassessed
# are there unpushed changes here? don't see the change beth made from 1.2 -> 1.1
# count converged ones

# scenarios <- readRDS(here("data", "scenarios.RDS"))
# low.lam.combos <- readRDS(here("data","low.lam.params.RDS"))
# med.lam.combos <- readRDS(here("data","med.lam.params.RDS"))
# high.lam.combos <- readRDS(here("data","high.lam.params.RDS"))

# load data
# TODO - adjust as in previous file
scenarios <- readRDS(here("data", "demographic_scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename(
    "phi1" = "S.J", 
    "phiad" = "S.A",
    "fec" = "f"
  )
# low.lam.combos <- readRDS(here("data","low.lam.params.RDS"))
# med.lam.combos <- readRDS(here("data","med.lam.params.RDS"))
# high.lam.combos <- readRDS(here("data","high.lam.params.RDS"))
# TODO fix hard coding
low.lam.params <- scenarios %>% 
  filter(trend == "decline")
med.lam.params <- scenarios %>% 
  filter(trend == "stable")
high.lam.params <- scenarios %>% 
  filter(trend == "increase")

# # functions
# source(here("scripts", "current version",
#             "0 - preparing scenarios", "compute_time_calc.R"))
# 
# # determine priority score for scenarios
# scenarios %>% mutate(priority = NA_integer_)
# for (i in 1:nrow(scenarios)) {
#   tmp <- scenarios[i, 1:3]
#   tmp <- tmp[!is.na(tmp)]
#   scenarios[i, "priority"] <- length(unique(tmp))
# }
# scenarios <- scenarios %>% arrange(priority) # save in prioritized order

# set working directory
setwd("C:/Users/AbbyBratt/Desktop/IPM SIM results/")

# scenarios.picked <- 25
# sims.per <- 25
# scenario <- dim(scenarios)[1]

scenarios.picked <- nrow(low.lam.params) # TODO note change here
sims.per <- 50 # TODO - is this what we decided? didn't we discuss either less
data_scenarios <- readRDS(here("data", "data_scenarios.RDS"))

# SLOW TO RUN 
for (i in 1:scenarios.picked) { #scenarios picked
  for (j in 1:sims.per) { # sims per
    for (k in 1:nrow(data_scenarios)) { # simulation scenario
      #if (i == 25 & j == 25 & k == 63) next
      print(paste(i, j, k), sep = " ")
      #if (scenarios[k, "lambda"] == "L") {
        out_L <- readRDS(paste("lowout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
        tmp <- max(gelman.diag(out_L, multivariate = FALSE)[[1]][, 1])
        # TODO - consider whether *everything* needs to have converged
        # TODO - for things that didn't converge, should we report on what didn't converge?
        # TODO - do we need to run longer? in that case need to fix that DLL error
        if (!is.na(tmp) & tmp <= 1.1) {
          out_L <- out_L %>% 
            collapse_chains() %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            filter(row_number() %% 60 == 1) %>% # thin chains
            mutate(scenario= i) %>% 
            mutate(sims = j) %>% 
            mutate(simscenarios = k)
          assign(paste("lowout", "-", i, "-", j, "-", k, sep = ""), out_L)
        }
        out_M <- readRDS(paste("medout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
        tmp <- max(gelman.diag(out_M, multivariate = FALSE)[[1]][, 1])
        if (!is.na(tmp) & tmp <= 1.1) {
          out_M <- out_M %>% 
            collapse_chains() %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            filter(row_number() %% 60 == 1) %>% # thin chains
            mutate(scenario= i) %>% 
            mutate(sims = j) %>% 
            mutate(simscenarios = k)
          assign(paste("medout", "-", i, "-", j, "-", k, sep = ""), out_M)
        }
        out_H <- readRDS(paste("highout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
        tmp <- max(gelman.diag(out_H, multivariate = FALSE)[[1]][, 1])
        if (!is.na(tmp) & tmp <= 1.1) {
          out_H <- out_H %>% 
            collapse_chains() %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            filter(row_number() %% 60 == 1) %>% # thin chains
            mutate(scenario= i) %>% 
            mutate(sims = j) %>% 
            mutate(simscenarios = k)
          assign(paste("highout", "-", i, "-", j, "-", k, sep = ""), out_H)
        }
      #} else if (scenarios[k, "lambda"] == "M") {
      #   out <- readRDS(paste("medout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
      #   tmp <- max(gelman.diag(out, multivariate = FALSE)[[1]][, 1])
      #   if (!is.na(tmp) & tmp <= 1.2) {
      #     out <- out %>% 
      #       collapse_chains() %>% 
      #       as.matrix() %>% 
      #       as.data.frame() %>% 
      #       filter(row_number() %% 60 == 1) %>% # thin chains
      #       mutate(scenario= i) %>% 
      #       mutate(sims = j) %>% 
      #       mutate(simscenarios = k)
      #     assign(paste("medout", "-", i, "-", j, "-", k, sep = ""), out)
      #   }
      # } else if (scenarios[k, "lambda"] == "H") {
      #   out <- readRDS(paste("highout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
      #   tmp <- max(gelman.diag(out, multivariate = FALSE)[[1]][, 1])
      #   if (!is.na(tmp) & tmp <= 1.2) {
      #     out <- out %>% 
      #       collapse_chains() %>% 
      #       as.matrix() %>% 
      #       as.data.frame() %>% 
      #       filter(row_number() %% 60 == 1) %>% # thin chains
      #       mutate(scenario= i) %>% 
      #       mutate(sims = j) %>% 
      #       mutate(simscenarios = k)
      #     assign(paste("highout", "-", i, "-", j, "-", k, sep = ""), out)
      #   }
      #}
    }
  }
}

rm(list=grep("highout|medout|lowout",ls(),value=TRUE,invert=TRUE))
#save.image("processedIPMOutput.RData") # doesn't work - stack overflow error

row.low <- do.call(bind_rows, lapply( ls(patt="lowout"), get) )
row.med <- do.call(bind_rows, lapply( ls(patt="medout"), get) )
row.high <- do.call(bind_rows, lapply( ls(patt="highout"), get) )

write.csv(row.high, file = "highout.csv")
write.csv(row.med, file = "medout.csv")
write.csv(row.low, file = "lowout.csv")

#rm(list = ls())

# TODO
  # count number that converged 
nlow <- length(ls(patt="lowout"))
nmed <- length(ls(patt="medout"))
nhigh <- length(ls(patt="highout"))

# TODO 
  # should i run at 1.2 just for a start
  # take converged ones
  # make smaller
  # bundle
  # upload?? 

# going to need some help diagnosing what is wrong

# MISSING ONES
# med 2 44 5 

# should 
