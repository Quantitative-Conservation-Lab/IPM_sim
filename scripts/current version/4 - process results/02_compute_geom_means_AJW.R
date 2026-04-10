# load libraries ####
library(tidyverse)
library(ggplot2)
library(gtable)
library(coda)
library(readr)
library(here)
library(purrr)

# function to compute quantiles of every variable in dataframe ####
p <- c(0.025, 0.5, 0.975)
p_names <- map_chr(p, ~paste0(.x*100, "%"))
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

dem_scenarios <- readRDS(here("data", "demographic_scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename("phi1" = "S.J", 
         "phiad" = "S.A",
         "fec" = "f")  %>%
  mutate(dem_scenario = row_number())

#load all processed mcmc results
# results_lam <- readRDS(file = here('results', 'processed', 'results_all_vSJC.RDS')) %>%
  results_lam <- readRDS(file = here('results', 'processed', 'results_all_ind300_nsam5.RDS')) %>%
  mutate(iter = row_number()) %>%
  select(contains("lambda") | contains("sim_rep") | contains("scenario") | contains('type')) 

lambda_dat <- results_lam %>%
  group_by(sim_rep, model_type, surv_scenario, dem_scenario) %>%
  summarise(across(everything(), p_funs)) %>% # deprecated dplyr code above
  pivot_longer( # begin reshaping
    cols = contains("%"), 
    names_to = "quantile",
    values_to = "lambda",
  ) %>% 
  mutate(Year = str_extract(quantile, "\\[\\d+\\]")) %>% 
  mutate(Year = str_extract(quantile, "\\d+")) %>% 
  mutate(Quantile = str_extract(quantile, "\\d+\\.?\\d?%")) %>% 
  select(-quantile) %>% 
  # select(c(1:3, 6, 5, 4)) %>% #AJW not sure what this did
  arrange(sim_rep, surv_scenario, dem_scenario, model_type, Quantile, Year) %>% 
  pivot_wider(names_from = Year, values_from = lambda, names_prefix = "Year_") %>% 
  select(c(sim_rep, surv_scenario, dem_scenario, model_type, Quantile, "Year_1", "Year_2", "Year_3", "Year_4", "Year_5", "Year_6", 
                "Year_7", "Year_8", "Year_9", "Year_10", "Year_11", "Year_12", 
                "Year_13", "Year_14")) # end reshaping


# create new variables for geometric mean by year
for (i in 1:14) { # number of years
  lambda_dat <- lambda_dat %>%
    mutate("geomean.{i}" :=  NA_real_)
}

# compute geometric means #######
num_id_cols <- 5
for(i in 1:dim(lambda_dat)[1]) {
  print(paste("row", i))
  for(j in 1:((ncol(lambda_dat) - num_id_cols)/2)) {
    lambda_dat[i,((ncol(lambda_dat)+num_id_cols)/2+j)] <- exp(mean(unlist(log(lambda_dat[i,num_id_cols+1:j]))))
  }
}

# save objects  #######
# write_csv(lambda_dat, here('results', 'processed', "lambda_geo_vSJC.csv"))
# write_csv(lambda_dat, here('results', 'processed', "lambda_geo_ind300_nsam5.csv"))


