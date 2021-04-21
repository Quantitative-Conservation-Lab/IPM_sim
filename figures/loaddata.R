# LOAD data

library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)
library(coda)
library(nlist)
library(beepr)

setwd("/Users/aebratt/Desktop/IPMEURING")

# SLOW TO RUN BUT DOES RUN
for (i in 1:10) { #scenarios picked
  for (j in 1:25) { # sims per
    for (k in 1:9) { # simulation scenario
      print(paste(i, j, k), sep = " ")
      if (k %in% 1:3) {
        out <- readRDS(paste("lowout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
        tmp <- max(gelman.diag(out, multivariate = FALSE)[[1]][, 1])
        if (!is.na(tmp) & tmp <= 1.2) {
          out <- out %>% 
            collapse_chains() %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            filter(row_number() %% 60 == 1)
          assign(paste("lowout", "-", i, "-", j, "-", k, sep = ""), out)
        }
      } else if (k %in% 4:6) {
        out <- readRDS(paste("medout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
        tmp <- max(gelman.diag(out, multivariate = FALSE)[[1]][, 1])
        if (!is.na(tmp) & tmp <= 1.2) {
          out <- out %>% 
            collapse_chains() %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            filter(row_number() %% 60 == 1)
          assign(paste("medout", "-", i, "-", j, "-", k, sep = ""), out)
        }
      } else if (k %in% 7:9) {
        out <- readRDS(paste("highout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
        tmp <- max(gelman.diag(out, multivariate = FALSE)[[1]][, 1])
        if (!is.na(tmp) & tmp <= 1.2) {
          out <- out %>% 
            collapse_chains() %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            filter(row_number() %% 60 == 1)
          assign(paste("highout", "-", i, "-", j, "-", k, sep = ""), out)
        }
      }
    }
  }
}

# presumably all of these are good to go, so let's collapse the chains 
# and thin the

rm(list=grep("highout|medout|lowout",ls(),value=TRUE,invert=TRUE))
save.image("processedIPMOutput.RData")

# repeat this for case study
# SAVE THIS AS RDATA
setwd("/Users/aebratt/Desktop/IPMEURING")

# LOAD data

rm(list = ls())
setwd("/Users/aebratt/Desktop/OOPSIPM")

# SLOW TO RUN BUT DOES RUN
for (i in 1:10) { #scenarios picked
  for (j in 1:25) { # sims per
    for (k in 1:9) { # simulation scenario
      print(paste(i, j, k), sep = " ")
        out <- readRDS(paste("lowout", "-", i, "-", j, "-", k, ".RDS", sep = ""))
        tmp <- max(gelman.diag(out, multivariate = FALSE)[[1]][, 1])
        if (!is.na(tmp) & tmp <= 1.2) {
          out <- out %>% 
            collapse_chains() %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            filter(row_number() %% 60 == 1)
          assign(paste("lowout", "-", i, "-", j, "-", k, sep = ""), out)
        }
    }
  }
}

# presumably all of these are good to go, so let's collapse the chains 
# and thin the

rm(list=grep("highout|medout|lowout",ls(),value=TRUE,invert=TRUE))
save.image("processedIPMcasestudy.RData")
