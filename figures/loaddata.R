# LOAD data

library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)
library(coda)
library(nlist)
library(beepr)

scenarios <- readRDS(here("data", "scenarios.RDS"))
low.lam.combos <- readRDS(here("data","low.lam.params.RDS"))
med.lam.combos <- readRDS(here("data","med.lam.params.RDS"))
high.lam.combos <- readRDS(here("data","high.lam.params.RDS"))

# functions
source(here("scripts", "current version",
            "0 - preparing scenarios", "compute_time_calc.R"))

# determine priority score for scenarios
scenarios %>% mutate(priority = NA_integer_)
for (i in 1:nrow(scenarios)) {
  tmp <- scenarios[i, 1:3]
  tmp <- tmp[!is.na(tmp)]
  scenarios[i, "priority"] <- length(unique(tmp))
}
scenarios <- scenarios %>% arrange(priority) # save in prioritized order

# set working directory
setwd("D:/all_results")

scenarios.picked <- 25
sims.per <- 25
scenario <- dim(scenarios)[1]

# SLOW TO RUN BUT DOES RUN
for (i in 1:scenarios.picked) { #scenarios picked
  for (j in 1:sims.per) { # sims per
    for (k in 1:scenario) { # simulation scenario
      print(paste(i, j, k), sep = " ")
      if (scenarios[k, "lambda"] == "L") {
        out <- readRDS(paste("lowout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
        tmp <- max(gelman.diag(out, multivariate = FALSE)[[1]][, 1])
        if (!is.na(tmp) & tmp <= 1.2) {
          out <- out %>% 
            collapse_chains() %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            filter(row_number() %% 60 == 1) %>% 
            mutate(scenario= i) %>% 
            mutate(sims = j) %>% 
            mutate(simscenarios = k)
          assign(paste("lowout", "-", i, "-", j, "-", k, sep = ""), out)
        }
      } else if (scenarios[k, "lambda"] == "M") {
        out <- readRDS(paste("medout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
        tmp <- max(gelman.diag(out, multivariate = FALSE)[[1]][, 1])
        if (!is.na(tmp) & tmp <= 1.2) {
          out <- out %>% 
            collapse_chains() %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            filter(row_number() %% 60 == 1) %>% 
            mutate(scenario= i) %>% 
            mutate(sims = j) %>% 
            mutate(simscenarios = k)
          assign(paste("medout", "-", i, "-", j, "-", k, sep = ""), out)
        }
      } else if (scenarios[k, "lambda"] == "H") {
        out <- readRDS(paste("highout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
        tmp <- max(gelman.diag(out, multivariate = FALSE)[[1]][, 1])
        if (!is.na(tmp) & tmp <= 1.2) {
          out <- out %>% 
            collapse_chains() %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            filter(row_number() %% 60 == 1) %>% 
            mutate(scenario= i) %>% 
            mutate(sims = j) %>% 
            mutate(simscenarios = k)
          assign(paste("highout", "-", i, "-", j, "-", k, sep = ""), out)
        }
      }
    }
  }
}

rm(list=grep("highout|medout|lowout",ls(),value=TRUE,invert=TRUE))
save.image("processedIPMOutput.RData")
rm(list = ls())
