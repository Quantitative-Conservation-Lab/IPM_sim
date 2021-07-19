# figure out which scenarios did not run

# output files

library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)
library(coda)
library(nlist)
library(beepr)

setwd("D:/all_results")

files <- list.files() 

