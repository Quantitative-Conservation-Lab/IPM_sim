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
library(purrr)

setwd("D:/all_results")

files <- list.files()

files <- readRDS("~/Desktop/IPM_sim/scripts/current version/files.RDS")
files <- files[str_detect(files, "RDS")] %>% 
  str_split("-", 2) %>% 
  sapply("[", 2) %>% 
  str_remove(".RDS") %>% 
  sort()
  

scenarios <- readRDS(here("data", "scenarios.RDS"))
# functions

# determine priority score for scenarios
scenarios %>% mutate(priority = NA_integer_)
expected.number <- 25 * 25 * dim(scenarios)[1]

scenarios.picked <- 1:25
sims.per <- 1:25
scenario <- 1:dim(scenarios)[1]

rows <- expand_grid(scenarios.picked, sims.per, scenario) %>% 
  unite(col = "string", sep = "-") %>% 
  as_vector() %>% 
  unname()

rerun <- rows[which(!(rows %in% files))] %>% 
  as.data.frame %>% 
  separate(".", into = c("c1", "c2", "c3")) %>% 
  sapply(as.numeric ) %>% 
  as.data.frame() %>% 
  arrange(c3, c2, c1)

# aaaaa <- rerun %>% 
#   filter((!c3 %in% c(1:9, 100:108, 118:126, 37:45, 46:54, 55:63)))

bbb <- rerun %>% 
  filter(c3 %in% 1:9)

# 1:9 - these are on loon - need to transfer over

# ONLY CERTAIN PARAM COMBINATIONS DID NOT CONVERGE - discuss this with the group
# 3-9, 25
# 37-45 - ALL
# 46:54 - ALL
# 55:63 - ALL

# amanda noticed these got fucked up
# 100:108 - 14:25
# 118:126 - 1-13

