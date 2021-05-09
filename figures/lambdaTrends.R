#### Processing #######

library(tidyverse)
library(ggplot2)
library(gtable)
library(RColorBrewer)
library(wesanderson)
library(coda)

pal <- rev(wes_palette("Zissou1", 3, type = "continuous"))

# row.low <- do.call(rbind, lapply( ls(patt="lowout"), get) ) %>%
#   select(contains("lambda") | contains("sims") | contains("scenario"))
# row.med <- do.call(rbind, lapply( ls(patt="medout"), get) ) %>%
#   select(contains("lambda") | contains("sims") | contains("scenario"))
# row.high <- do.call(rbind, lapply( ls(patt="highout"), get) ) %>%
#   select(contains("lambda") | contains("sims") | contains("scenario"))
# 
# rm(list=grep("highout|medout|lowout",ls(),value=TRUE,invert=FALSE))

# TODO
# generate geometric means

# row 1
# declining population

# for (i in 1:14) {
#   row.low <- row.low %>%
#     mutate("geomean.{i}" :=  NA_real_)
#   row.med <- row.med %>%
#     mutate("geomean.{i}" :=  NA_real_)
#   row.high <- row.high %>%
#     mutate("geomean.{i}" :=  NA_real_)
# }
# 
# # super slow
# for(i in 1:dim(row.low)[1]) {
#   print(paste("row", i))
#   for(j in 1:((ncol(row.low) - 3 )/2)) {
#     row.low[i, ((ncol(row.low)+3)/2 + j)] <- exp(mean(unlist(log(row.low[i, 1:j]))))
#   }
# }
# 
# # super slow
# for(i in 1:dim(row.med)[1]) {
#   print(paste("row", i))
#   for(j in 1:((ncol(row.med) - 3 )/2)) {
#     row.med[i, ((ncol(row.med)+3)/2 + j)] <- exp(mean(unlist(log(row.med[i, 1:j]))))
#   }
# }
# 
# # super slow
# for(i in 1:dim(row.high)[1]) {
#   print(paste("row", i))
#   for(j in 1:((ncol(row.high) - 3 )/2)) {
#     row.high[i, ((ncol(row.high) + 3)/2 + j)] <- exp(mean(unlist(log(row.high[i, 1:j]))))
#   }
# }
# 
# # save objects
# 
# #saveRDS(row.low, "row.low.RDS")

## load RDS
row.low <- readRDS(file = here::here('data', 'EuringPosterReg', 'row.low.rds'))
row.med <- readRDS(file = here::here('data', 'EuringPosterReg', 'row.med.rds'))
row.high <- readRDS(file = here::here('data', 'EuringPosterReg', 'row.high.rds'))

# Reformat for plotting
toplot1 <- row.low %>%
  select(contains("geomean"), scenario, sims, simscenarios) %>%
  #group_by(model, detection)
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(scenario, sims, simscenarios, Year) %>% # checked through here - TODO
  # summarise(low = quantile(value, 0.025),
  #           med = quantile(value, 0.5),
  #           high = quantile(value, 0.975)) %>%
  ungroup() %>%
  mutate(scenario = as.factor(scenario),
         sims = as.factor(sims), 
         simscenarios = as.factor(simscenarios)) 

toplot2 <- row.med %>%
  select(contains("geomean"), scenario, sims, simscenarios) %>%
  #group_by(model, detection)
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(scenario, sims, simscenarios, Year) %>% # checked through here - TODO
  # summarise(low = quantile(value, 0.025),
  #           med = quantile(value, 0.5),
  #           high = quantile(value, 0.975)) %>%
  ungroup() %>%
  mutate(scenario = as.factor(scenario),
         sims = as.factor(sims), 
         simscenarios = as.factor(simscenarios))

toplot3 <- row.high %>%
  select(contains("geomean"), scenario, sims, simscenarios) %>%
  #group_by(model, detection)
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(scenario, sims, simscenarios, Year) %>% # checked through here - TODO
  # summarise(low = quantile(value, 0.025),
  #           med = quantile(value, 0.5),
  #           high = quantile(value, 0.975)) %>%
  ungroup() %>%
  mutate(scenario = as.factor(scenario),
         sims = as.factor(sims), 
         simscenarios = as.factor(simscenarios)) 

# toplot2 <- row.med %>%
#   select(contains("geomean"), model, detection) %>%
#   #group_by(model, detection)
#   pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
#   filter(!is.na(value)) %>%
#   mutate(Year = str_remove(Year, "geomean\\.")) %>%
#   mutate(Year = as.numeric(Year)) %>%
#   group_by(model, detection, Year) %>%
#   summarise(low = quantile(value, 0.025),
#             med = quantile(value, 0.5),
#             high = quantile(value, 0.975)) %>%
#   ungroup() %>%
#   mutate(model = as.factor(model),
#          detection = as.factor(detection)) %>%
#   mutate(high = if_else(high > 1.3, 1.299, high))
# 
# toplot3 <- row.high %>%
#   select(contains("geomean"), model, detection) %>%
#   #group_by(model, detection)
#   pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
#   filter(!is.na(value)) %>%
#   mutate(Year = str_remove(Year, "geomean\\.")) %>%
#   mutate(Year = as.numeric(Year)) %>%
#   group_by(model, detection, Year) %>%
#   summarise(low = quantile(value, 0.025),
#             med = quantile(value, 0.5),
#             high = quantile(value, 0.975)) %>%
#   ungroup() %>%
#   mutate(model = as.factor(model),
#          detection = as.factor(detection))
