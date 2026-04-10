## IPM Sim Figures
## Adapted from PaperFigures.RMD
## A DuVall
## 15 Oct 2021
## Updated by A Bratt
## 8 Jan 2026
## Updated by A Warlick
## 8 Jan 2026, 3/24/2026

##################### Libraries ###########################

library(tidyverse)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(coda)
library(knitr)
library(reshape2)
library(here)
library(RColorBrewer)
library(colorspace)
library(ggh4x)
library(Hmisc)
library(patchwork)
library(tidybayes)

rainbow2 <- c("violetred4", "dodgerblue3", 'deepskyblue1', "#4aaaa5", "#a3d39c", "#f6b61c", "chocolate2", "red3")

#demographic scenarios: lambda, life history, vital rate combos
dem_scenarios <- readRDS(here("data", "demographic_scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename("phi1" = "S.J", 
    "phiad" = "S.A",
    "fec" = "f")  %>%
  mutate(dem_scenario = row_number()) 

det.abund <- factor(x = c("L", "M", "H"))
det.MR <- factor(x = c("L", "M", "H", "NA"))
det.prod <- factor(x = c("L", "M", "H", "NA"))
lambda <- factor(x = c("L", "M", "H"))

surv_scenarios <- readRDS(here("data", "data_scenarios.RDS")) %>%
  mutate(surv_scenario = row_number()) %>%
  transform(p.surv = ifelse(det.abund == 'L', 0.3, 
                                ifelse(det.abund == 'M', 0.5, 
                                       ifelse(det.abund == 'H', 0.8, NA)))) %>%
  transform(mean.p = ifelse(det.MR == 'L', 0.3, 
                                ifelse(det.MR == 'M', 0.5, 
                                       ifelse(det.MR == 'H', 0.8, NA)))) %>%
  transform(det.prod = ifelse(det.prod == 'L', 0.3, 
                            ifelse(det.prod == 'M', 0.5, 
                                   ifelse(det.prod == 'H', 0.8, NA)))) %>%
  dplyr::select(-c(det.abund, det.MR))

#character version to merge with stats below
surv_scenarios_char <- readRDS(here("data", "data_scenarios.RDS")) %>%
  mutate(surv_scenario = row_number())

##################### Data Prep ###########################

# results <- readRDS(file = here('results', 'processed', 'results_all_ind400.RDS')) %>%
  # results <- readRDS(file = here('results', 'processed', 'results_all_normObs.RDS')) %>%
# results <- readRDS(file = here('results', 'processed', 'results_all_ind300_nsam5.RDS')) %>%
results <- readRDS(file = here('results', 'processed', 'results_all_vSJC.RDS')) %>%
  dplyr::rename(phi1 = `mean.phi[1]`, phiad = `mean.phi[2]`) %>%
  inner_join(dem_scenarios, by = 'dem_scenario', suffix = c('.obs', '.true')) %>%
  inner_join(surv_scenarios, by = 'surv_scenario', suffix = c('.obs', '.true')) 

# summarize medians, sd, cvs, relative bias, and error at the model (sim) level
all.stats.sims <- results %>%
  dplyr::select('phi1.obs', 'phi1.true', 'phiad.obs', 'phiad.true', 'fec.obs', 'fec.true',
                "p.surv.obs", "p.surv.true", "mean.p.obs", "mean.p.true",
                'dem_scenario', 'surv_scenario', 'model_type', 'trend', 'life_hist', 'sim_rep') %>%
  group_by(trend, model_type, life_hist, surv_scenario, dem_scenario, sim_rep) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(-c(trend, model_type, life_hist, surv_scenario, dem_scenario, sim_rep, iter)) %>% 
  transform(name = ifelse(name == 'p.surv.obs', 'psurv.obs',
                          ifelse(name == 'p.surv.true', 'psurv.true',
                                 ifelse(name == 'mean.p.obs', 'meanp.obs', 
                                        ifelse(name == 'mean.p.true', 'meanp.true', name))))) %>%
  separate_wider_delim(name, delim = ".", names = c("param", "type")) %>% 
  
  reshape2::dcast(trend + model_type + life_hist + surv_scenario + dem_scenario + param + 
                   sim_rep + iter ~ type, value.var = 'value') %>%
  #take mean over iters but keep sim_rep
  mutate(error = (obs - true)) %>% 
  group_by(trend, model_type, life_hist, surv_scenario, dem_scenario, param, sim_rep) %>% 
  dplyr::summarize(
    median = median(obs), 
    sd = sd(obs), 
    cv = sd(obs)/mean(obs),
    rb = (mean(obs) - mean(true)) / mean(true), # NOTE using mean(true) tho it doesnt vary within group
    rmse = sqrt(mean(error^2))
  )

# average performance stats across sim replicates
all.stats <- all.stats.sims %>% 
  select(trend, model_type, life_hist, surv_scenario, dem_scenario, param, sim_rep, cv, rb, rmse) %>% 
  group_by(trend, model_type, life_hist, surv_scenario, dem_scenario, param) %>% 
  dplyr::summarize(
    across(
      c(cv, rb, rmse),
      list(
        mean = ~ mean(., na.rm = T),
        lower = ~ quantile(., 0.025, na.rm = TRUE),
        upper = ~ quantile(., 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

##################### Bias ###########################

# rel.bias
rel.bias <- all.stats %>%
  inner_join(surv_scenarios_char, by = "surv_scenario") %>%
  dplyr::select('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type',
                "param", "rb_mean", 
                #'rb_lower', 'rb_upper',
                'det.MR', 'det.abund', 'det.prod') %>%
  reshape2::melt(id.vars = c('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type', 
                             "param", 'det.MR', 'det.abund', 'det.prod')) %>%
  # mean across life history strategy (dem_scenarios);  
  group_by(trend, model_type, det.MR, det.abund, det.prod, param) %>%
  ##take mean of upper and lower if add above
  dplyr::summarize(bias = mean(value), .groups = 'keep') %>% 
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec', 'psurv', 'meanp'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity', 'Count survey detection', 'MR detection'))) %>%
  transform(trend = factor(trend, levels = c('decline', 'stable', 'increase'),
                           labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = factor(model_type, levels = c('out_IPM', 'out_noProd', 'out_noMR', 'out_abundOnly'),
                             labels = c('Full IPM', 'Abundance & Survival', 
                                        'Abundance & Productivity', 'Abundance Only'))) 

# average over two layers of detection (det.MR and det.prod)
obs.pars <- c('MR detection', 'Count survey detection')
rel.bias.few <- rel.bias %>%
  group_by(variable, det.abund, trend, dataset) %>%
  dplyr::summarize(bias = mean(bias), .groups = 'keep') %>% # taking average
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) 

rel.bias.sc <- all.stats %>%
  inner_join(surv_scenarios_char, by = "surv_scenario") %>%
  dplyr::select('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type',
                "param", "rb_mean", 
                #'rb_lower', 'rb_upper',
                'det.MR', 'det.abund', 'det.prod') %>%
  reshape2::melt(id.vars = c('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type', 
                             "param", 'det.MR', 'det.abund', 'det.prod')) %>%
  #keep dem_scenario and life history, I think....? maybe just life_hist?
  group_by(trend, model_type, dem_scenario, life_hist, det.MR, det.abund, det.prod, param) %>%
  ##take mean of upper and lower if add above
  dplyr::summarize(bias = mean(value), .groups = 'keep') %>% 
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec', 'psurv', 'meanp'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity', 'Count survey detection', 'MR detection'))) %>%
  transform(trend = factor(trend, levels = c('decline', 'stable', 'increase'),
                           labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = factor(model_type, levels = c('out_IPM', 'out_noProd', 'out_noMR', 'out_abundOnly'),
                             labels = c('Full IPM', 'Abundance & Survival', 
                                        'Abundance & Productivity', 'Abundance Only')))

rel.bias.dem <- rel.bias.sc %>%
  group_by(trend, dataset, life_hist, det.MR, det.abund, det.prod, variable) %>%
  dplyr::summarize(bias = mean(bias), .groups = 'keep')

# facet by both fec and juv true vals; average over all detection levels
bias.plot.vals <- rel.bias.dem %>%
  group_by(variable, life_hist, trend, dataset) %>% 
  dplyr::summarize(value = mean(bias), .groups = 'keep') 

##################### RMSE ###########################

rmse.vals <- all.stats %>%
  inner_join(surv_scenarios_char, by = "surv_scenario") %>%
  dplyr::select('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type',
                "param", "rmse_mean", 
                #'rmse_lower', 'rmse_upper',
                'det.MR', 'det.abund', 'det.prod') %>%
  reshape2::melt(id.vars = c('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type', 
                             "param", 'det.MR', 'det.abund', 'det.prod')) %>%
  # mean across life history strategy (dem_scenarios);  
  group_by(trend, model_type, det.MR, det.abund, det.prod, param) %>%
  ##take mean of upper and lower if add above
  dplyr::summarize(rmse = mean(value), .groups = 'keep') %>% 
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec', 'psurv', 'meanp'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity', 'Count survey detection', 'MR detection'))) %>%
  transform(trend = factor(trend, levels = c('decline', 'stable', 'increase'),
                           labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = factor(model_type, levels = c('out_IPM', 'out_noProd', 'out_noMR', 'out_abundOnly'),
                             labels = c('Full IPM', 'Abundance & Survival', 
                                        'Abundance & Productivity', 'Abundance Only'))) 

# average over two layers of detection (det.MR and det.prod)
obs.pars <- c('MR detection', 'Count survey detection')
rmse.few <- rmse.vals %>%
  group_by(variable, det.abund, trend, dataset) %>%
  dplyr::summarize(rmse = mean(rmse), .groups = 'keep') %>% # taking average
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) 

rmse.sc <- all.stats %>%
  inner_join(surv_scenarios_char, by = "surv_scenario") %>%
  dplyr::select('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type',
                "param", "rb_mean", 
                #'rb_lower', 'rb_upper',
                'det.MR', 'det.abund', 'det.prod') %>%
  reshape2::melt(id.vars = c('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type', 
                             "param", 'det.MR', 'det.abund', 'det.prod')) %>%
  #keep dem_scenario and life history, I think....? maybe just life_hist?
  group_by(trend, model_type, dem_scenario, life_hist, det.MR, det.abund, det.prod, param) %>%
  ##take mean of upper and lower if add above
  dplyr::summarize(rmse = mean(value), .groups = 'keep') %>% 
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec', 'psurv', 'meanp'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity', 'Count survey detection', 'MR detection'))) %>%
  transform(trend = factor(trend, levels = c('decline', 'stable', 'increase'),
                           labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = factor(model_type, levels = c('out_IPM', 'out_noProd', 'out_noMR', 'out_abundOnly'),
                             labels = c('Full IPM', 'Abundance & Survival', 
                                        'Abundance & Productivity', 'Abundance Only')))

rmse.dem <- rmse.sc %>%
  group_by(trend, dataset, life_hist, det.MR, det.abund, det.prod, variable) %>%
  dplyr::summarize(rmse = mean(rmse), .groups = 'keep')

# facet by both fec and juv true vals; average over all detection levels
rmse.plot.vals <- rmse.dem %>%
  group_by(variable, life_hist, trend, dataset) %>% 
  dplyr::summarize(value = mean(rmse), .groups = 'keep') 
# 
# ##################### CV ###########################
# 
cv.vals <- all.stats %>%
  inner_join(surv_scenarios_char, by = "surv_scenario") %>%
  dplyr::select('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type',
                "param", "cv_mean", 
                #'cv_lower', 'cv_upper',
                'det.MR', 'det.abund', 'det.prod') %>%
  reshape2::melt(id.vars = c('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type', 
                             "param", 'det.MR', 'det.abund', 'det.prod')) %>%
  # mean across life history strategy (dem_scenarios);  
  group_by(trend, model_type, det.MR, det.abund, det.prod, param) %>%
  ##take mean of upper and lower if add above
  dplyr::summarize(cv = mean(value), .groups = 'keep') %>% 
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec', 'psurv', 'meanp'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity', 'Count survey detection', 'MR detection'))) %>%
  transform(trend = factor(trend, levels = c('decline', 'stable', 'increase'),
                           labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = factor(model_type, levels = c('out_IPM', 'out_noProd', 'out_noMR', 'out_abundOnly'),
                             labels = c('Full IPM', 'Abundance & Survival', 
                                        'Abundance & Productivity', 'Abundance Only'))) 

# average over two layers of detection (det.MR and det.prod)
obs.pars <- c('MR detection', 'Count survey detection')
cv.few <- cv.vals %>%
  group_by(variable, det.abund, trend, dataset) %>%
  dplyr::summarize(cv = mean(cv), .groups = 'keep') %>% # taking average
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) 

cv.sc <- all.stats %>%
  inner_join(surv_scenarios_char, by = "surv_scenario") %>%
  dplyr::select('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type',
                "param", "rb_mean", 
                #'rb_lower', 'rb_upper',
                'det.MR', 'det.abund', 'det.prod') %>%
  reshape2::melt(id.vars = c('trend', 'surv_scenario', 'dem_scenario', 'life_hist', 'model_type', 
                             "param", 'det.MR', 'det.abund', 'det.prod')) %>%
  #keep dem_scenario and life history, I think....? maybe just life_hist?
  group_by(trend, model_type, dem_scenario, life_hist, det.MR, det.abund, det.prod, param) %>%
  ##take mean of upper and lower if add above
  dplyr::summarize(cv = mean(value), .groups = 'keep') %>% 
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec', 'psurv', 'meanp'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity', 'Count survey detection', 'MR detection'))) %>%
  transform(trend = factor(trend, levels = c('decline', 'stable', 'increase'),
                           labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = factor(model_type, levels = c('out_IPM', 'out_noProd', 'out_noMR', 'out_abundOnly'),
                             labels = c('Full IPM', 'Abundance & Survival', 
                                        'Abundance & Productivity', 'Abundance Only')))

cv.dem <- cv.sc %>%
  group_by(trend, dataset, life_hist, det.MR, det.abund, det.prod, variable) %>%
  dplyr::summarize(cv = mean(cv), .groups = 'keep')

# facet by both fec and juv true vals; average over all detection levels
cv.plot.vals <- cv.dem %>%
  group_by(variable, life_hist, trend, dataset) %>% 
  dplyr::summarize(value = mean(cv), .groups = 'keep')

##################### LAMBDA #########################

# lambda_dat <- read.csv(file = here('results', 'processed', 'lambda_geo_vSJC.csv'), header = T, 
#                        stringsAsFactors = F)
lambda_dat <- read.csv(file = here('results', 'processed', 'lambda_geo_ind300_nsam5.csv'), header = T, 
                       stringsAsFactors = F)

#reformat for plotting
lam_dat <- lambda_dat %>%
  select(contains("geomean"), sim_rep, surv_scenario, dem_scenario, model_type, Quantile) %>%
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(sim_rep, surv_scenario, dem_scenario, model_type, Year, Quantile) %>% 
  ungroup() %>%
  inner_join(surv_scenarios_char, by = "surv_scenario") %>%
  inner_join(dem_scenarios, by = 'dem_scenario') %>% 
  rename(phi1.true = phi1, 
                      phiad.true = phiad, 
                      fec.true = fec) %>% 
  mutate(
    det.MR = na_if(as.character(det.MR), "NA"), 
    det.prod = na_if(as.character(det.prod), "NA"),
    det.abund = na_if(as.character(det.abund), "NA")) %>% 
  transform(dataset = ifelse(is.na(det.MR)&!is.na(det.prod), 'Abundance & Productivity', 
                             ifelse(!is.na(det.MR)&is.na(det.prod), 'Abundance & Survival',
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM')))) %>%
  group_by(Quantile, Year, det.abund, det.MR, det.prod, trend, lambda, life_hist, dataset) %>% 
  dplyr::summarize(value = mean(value), .groups = "drop") %>% 
  ungroup() %>% 
  mutate(Quantile = str_remove(Quantile, "\\%")) %>% 
  mutate(Quantile = paste("X", Quantile, sep = "")) %>% 
  reshape2::dcast(dataset + Year + det.MR + det.prod + det.abund + lambda + trend + life_hist ~ Quantile, value.var = "value") %>% 
  mutate(Year = Year + 1) %>% 
  filter(Year %in% c(15)) %>% 
  mutate(Year = factor(Year)) %>% 
  mutate(det.abund = factor(det.abund, levels = c("L", "M", "H"), labels = c("Low", "Medium", "High"))) %>% 
  mutate(det.prod = factor(det.prod, levels = c("L", "M", "H"), labels = c("Low", "Medium", "High"))) %>% 
  mutate(det.MR = factor(det.MR, levels = c("L", "M", "H"), labels = c("Low", "Medium", "High"))) %>% 
  transform(trend = factor(trend, levels = c("decline", "stable", "increase"),
                                     labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'))) 

# average over two layers of detection (det.MR and det.prod) and life history type
plot_dat_lam_few <- lam_dat %>%
  group_by(dataset, Year, det.abund, trend, lambda) %>%
  dplyr::summarize(
    `X2.5` = mean(`X2.5`), 
    `X50` = mean(`X50`),
    `X97.5` = mean(`X97.5`),
    .groups = 'drop') #%>%
  # transform(trend = factor(trend, levels = c("Decreasing", "Stable", "Increasing"))) %>%
  # transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
  #                            labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))

# average over all layers of detection and NOT life history type
plot_dat_lam_dem <- lam_dat %>%
  group_by(dataset, Year, life_hist, trend, lambda) %>%
  dplyr::summarize(
    `X2.5` = mean(`X2.5`), 
    `X50` = mean(`X50`),
    `X97.5` = mean(`X97.5`),
    .groups = 'drop') %>%
  # transform(lambda.scenario = factor(lambda.scenario,
  #                                    levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(life_hist = factor(life_hist, levels = c("slow", "mod", "fast"), 
                               labels = c("Slow", "Moderate", "Fast"))) #%>%
  # transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
  #                            labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))

################### FIGURES ##########################

#### Figure 3: RMSE and bias ecological paramters x count survey detection ####
dataset.labs <- c("Full IPM", "Abundance & Prod.", "Abundance & Surv.", "Abundance Only")
names(dataset.labs) <- c("Full IPM", "Abundance & Productivity", "Abundance & Survival", "Abundance Only")
lambda.labs <- c("Decrease", "Stable", "Increase")
names(lambda.labs) <- c("Decreasing", "Stable", "Increasing")

##### Bias dot plot ####
a1 <- ggplot(rel.bias.few %>% filter(variable %nin% obs.pars), 
             aes(x = det.abund, y = bias, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + 
  geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  #ylim(c(-1.75, 1.75)) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  xlab('Count survey detection') + 
  ylab('Relative bias') +
  facet_grid(dataset~trend, scales = 'free', 
             labeller = labeller(dataset = dataset.labs, trend = lambda.labs)) +
  # scale_y_continuous(limits = c(-1.2, 1.2), breaks = c(-1,0,1)) +
  theme_bw() +
  theme(legend.position = 'top',
        #plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        #strip.text = element_text(color = "black"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 10, vjust = 0.75),
        axis.title = element_text(size = 10, vjust = 0.75),
        strip.text = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        #panel.spacing.x = unit(0.75, "line")) +
        #scale_color_manual(values = rainbow2[-c(1,4)], name = '') +
        #scale_shape_manual(values = c(15, 16, 17), name = '')
        panel.spacing.x = unit(0.75, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = rainbow2[-c(1,4)], name = '', 
                     labels = c(expression(phi["2"]), 
                                expression(phi["1"]), expression(f))) +
  scale_shape_manual(values = c(15, 16, 17), name = '',
                     labels = c(expression(phi["2"]),
                                expression(phi["1"]), expression(f)))
a1

### RMSE dot plot 
a2 <- ggplot(rmse.few %>% filter(variable %nin% obs.pars), 
             aes(x = det.abund, y = rmse, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + geom_line() +
  # geom_tile(color = 'grey50') +
  xlab('Count survey detection') + ylab('RMSE') +
  facet_grid(dataset ~ trend, drop = T, scales = 'free_x', 
             labeller = labeller(dataset = dataset.labs, 
                                 trend = lambda.labs)) +
  scale_fill_gradient2(name = "RMSE",
                       #mid = "white", high = rainbow2[2], midpoint = 0) +
                       #low = "white", mid = rainbow2[3], high = rainbow2[2]) + #,
                       low = "white", high = rainbow2[2]) + #,
  #midpoint = 0.5) + # TODO - note change here 
  theme_bw() +
  theme(legend.position = 'top',
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 10, vjust = 0.75),
        axis.title = element_text(size = 10, vjust = 0.75),
        strip.text = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        #panel.spacing.x = unit(0.75, "line")) +
        #scale_color_manual(values = rainbow2[-c(1,4)], name = '') +
        #scale_shape_manual(values = c(15, 16, 17), name = '')
        panel.spacing.x = unit(0.75, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  scale_color_manual(values = rainbow2[-c(1,4)], name = '', 
                     labels = c(expression(phi["2"]), 
                                expression(phi["1"]), expression(f))) +
  scale_shape_manual(values = c(15, 16, 17), name = '',
                     labels = c(expression(phi["2"]),
                                expression(phi["1"]), expression(f)))
a2

### CV dot plot 
a3 <- ggplot(cv.few %>% filter(variable %nin% obs.pars), 
             aes(x = det.abund, y = cv, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + geom_line() +
  # geom_tile(color = 'grey50') +
  xlab('Count survey detection') + ylab('Coefficient of variation (CV)') +
  #facet_grid(dataset ~ lambda.scenario, drop = T, scales = 'free_x', labeller = label_wrap_gen()) +
  facet_grid(dataset ~ trend, drop = T, scales = 'free_x', 
             labeller = labeller(dataset = dataset.labs, 
                                 trend = lambda.labs)) +
  scale_fill_gradient2(name = "CV",
                       #mid = "white", high = rainbow2[2], midpoint = 0) +
                       #low = "white", mid = rainbow2[3], high = rainbow2[2]) + #,
                       low = "white", high = rainbow2[2]) + #,
  theme_bw() +
  theme(legend.position = 'top',
        #plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        #strip.text = element_text(color = "black"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 10, vjust = 0.75),
        axis.title = element_text(size = 10, vjust = 0.75),
        strip.text = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        #panel.spacing.x = unit(0.75, "line")) +
        #scale_color_manual(values = rainbow2[-c(1,4)], name = '') +
        #scale_shape_manual(values = c(15, 16, 17), name = '')
        panel.spacing.x = unit(0.75, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  scale_color_manual(values = rainbow2[-c(1,4)], name = '', 
                     labels = c(expression(phi["2"]), 
                                expression(phi["1"]), expression(f))) +
  scale_shape_manual(values = c(15, 16, 17), name = '',
                     labels = c(expression(phi["2"]),
                                expression(phi["1"]), expression(f)))
a3

##### Combine ####
# plot_grid(a1, a2, a3, nrow = 3, labels = "AUTO",
#           align = "hv", label_size = 12)
# ggsave(width = 6.5, height = 18, here("figures", 'final', "fig3.png"))

plot_grid(a1, a2, a3, nrow = 2, labels = "AUTO",
          align = "hv", label_size = 12)
# ggsave(width = 8, height = 12, here("figures", 'final', "fig3.png"))


#### Figure 4: RMSE and bias ecological parameters x true fecundity ####

phi1_cat_lab <- c("True first-year\nφ: Low", 
                  "True first-year\nφ: Med",
                  "True first-year\nφ: High")

names(phi1_cat_lab) <- c("True first-year survival: Low", 
                         "True first-year survival: Medium", 
                         "True first-year survival: High")

##### Bias dot plot ####
b1 <- ggplot(bias.plot.vals %>% filter(variable %nin% obs.pars), 
             aes(x = life_hist, y = value, col = factor(variable), group = factor(variable),
                            shape = factor(variable))) +
  geom_point() + geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  #scale_x_discrete(labels = c("L", "M", "H")) +
  xlab('Life history type') + ylab('Relative bias') +
  facet_grid(dataset~trend, scales = 'free') + #, 
  #labeller = labeller(dataset = dataset.labs, phi1_cat = phi1_cat_lab)) +
  #labeller = label_wrap_gen()) +
  #ylim(c(-1.75, 1.75)) +
  # scale_y_continuous(limits = c(-1.75, 1.75), breaks = c(-1.5,0,1.5)) +
  scale_color_manual(values = rainbow2[-c(1,4)], name = '',
                     labels = c(expression(phi["2"]), expression(phi["1"]), expression(f))) +
  scale_shape_manual(values = c(15, 16, 17), name = '',
                     labels = c(expression(phi["2"]), expression(phi["1"]), expression(f))) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 10, vjust = 0.75),
        axis.title = element_text(size = 10, vjust = 0.75),
        #plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #+
#scale_color_manual(values = rainbow2[-c(1,4)], name = '') +
#scale_shape_manual(values = c(15, 16, 17), name = '')
b1

#dot plot
b2 <- ggplot(rmse.plot.vals %>% filter(variable %nin% obs.pars), 
             aes(x = life_hist, y = value, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + geom_line() +
  # geom_tile(color = 'grey50') +
  xlab('Life history type') + ylab('RMSE') +
  #facet_grid(dataset ~ trend, drop = T, scales = 'free_x', labeller = label_wrap_gen()) +
  facet_grid(dataset ~ trend, drop = T, scales = 'free', 
             labeller = labeller(dataset = dataset.labs, 
                                 trend = lambda.labs)) +
  scale_fill_gradient2(name = "RMSE",
                       #mid = "white", high = rainbow2[2], midpoint = 0) +
                       #low = "white", mid = rainbow2[3], high = rainbow2[2]) + #,
                       low = "white", high = rainbow2[2]) + #,
  #midpoint = 0.5) + # TODO - note change here 
  scale_color_manual(values = rainbow2[-c(1,4)], name = '',
                     labels = c(expression(phi["2"]), expression(phi["1"]), expression(f))) +
  scale_shape_manual(values = c(15, 16, 17), name = '',
                     labels = c(expression(phi["2"]), expression(phi["1"]), expression(f))) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 10, vjust = 0.75),
        axis.title = element_text(size = 10, vjust = 0.75),
        #plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
b2

## dot plot
b3 <- ggplot(cv.plot.vals %>% filter(variable %nin% obs.pars), 
             aes(x = life_hist, y = value, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + geom_line() +
  # geom_tile(color = 'grey50') +
  xlab('Life history type') + ylab('Coefficient of variation (CV)') +
  #facet_grid(dataset ~ trend, drop = T, scales = 'free_x', labeller = label_wrap_gen()) +
  facet_grid(dataset ~ trend, drop = T, scales = 'free', 
             labeller = labeller(dataset = dataset.labs, 
                                 trend = lambda.labs)) +
  scale_fill_gradient2(name = "CV",
                       #mid = "white", high = rainbow2[2], midpoint = 0) +
                       #low = "white", mid = rainbow2[3], high = rainbow2[2]) + #,
                       low = "white", high = rainbow2[2]) + #,
  #midpoint = 0.5) + # TODO - note change here 
  scale_color_manual(values = rainbow2[-c(1,4)], name = '',
                     labels = c(expression(phi["2"]), expression(phi["1"]), expression(f))) +
  scale_shape_manual(values = c(15, 16, 17), name = '',
                     labels = c(expression(phi["2"]), expression(phi["1"]), expression(f))) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 10, vjust = 0.75),
        axis.title = element_text(size = 10, vjust = 0.75),
        #plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
b3

##### Combine ####
# plot_grid(b1, b2, b3, nrow = 3, labels = "AUTO", align = "hv", label_size = 12)
# ggsave(width = 6.5, height = 18, here("figures", "fig4_NEW.png"))

plot_grid(b1, b2, b3, nrow = 2, labels = "AUTO", align = "hv", label_size = 12)
# ggsave(width = 11, height = 15, here("figures", 'final', "fig4.png"))


#### Figure 5: Lambda trends ####

##### WRT count survey detection ####
c1 <- ggplot(plot_dat_lam_few) +
  geom_point(aes(x = Year, y = X50, col = det.abund, group = det.abund, shape = det.abund), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(x = Year, ymin = X2.5, ymax = X97.5, col = det.abund, group = det.abund,
                     shape = det.abund), position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = as.numeric(as.character(lambda))), linetype = 'dotted') +
  geom_hline(aes(yintercept = 1.0), linetype = 'solid') +
  xlab('Final year (t=15)') + #renamed the axis - TODO is this correct
  # oh right, yes, because we should have stopped the model after each year 
  # (to 'hide' the full time series from being used to estimate trend)
  ylab(expression(lambda)) + 
  facet_grid(dataset ~ trend, scales = 'free', labeller = label_wrap_gen()) +
  theme_bw() +
  theme(legend.position = 'top',
        plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        #HAS: removed axis labels
        axis.text.x = element_blank(),#element_text(angle = 0, vjust = 1.5),
        axis.ticks.x=element_blank(),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line")) +
  scale_color_manual(values = rainbow2[-c(1,4)], name = 'Count survey detection level') +
  scale_shape_manual(values = c(15, 16, 17), name = 'Count survey detection level')
c1

##### WRT life history type ####
c2 <- ggplot(plot_dat_lam_dem) +
  geom_point(aes(x = Year, y = X50, col = life_hist, group = life_hist, shape = life_hist), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(x = Year, ymin = X2.5, ymax = X97.5, col = life_hist, group = life_hist,
                     shape = life_hist), position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = as.numeric(as.character(lambda))), linetype = 'dotted') +
  geom_hline(aes(yintercept = 1.0), linetype = 'solid') +
  xlab('Final year (t=15)') + #renamed the axis - TODO is this correct
  # oh right, yes, because we should have stopped the model after each year 
  # (to 'hide' the full time series from being used to estimate trend)
  ylab(expression(lambda)) + 
  facet_grid(dataset ~ trend, scales = 'free', labeller = label_wrap_gen()) +
  theme_bw() +
  theme(legend.position = 'top',
        plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        #HAS: removed axis labels
        axis.text.x = element_blank(),#element_text(angle = 0, vjust = 1.5),
        axis.ticks.x=element_blank(),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line")) +
  scale_color_manual(values = rainbow2[-c(1,4)], name = 'Life History Type') +
  scale_shape_manual(values = c(15, 16, 17), name = 'Life History Type')
c2

##### Combine ####
plot_grid(c1, c2, ncol = 2, labels = "AUTO", align = "hv", label_size = 12)
# ggsave(width = 10, height = 7.5, here("figures", 'final', "fig5.png"))

################### Appendix #########################

#### Figure 6: RMSE, CV, bias of estimated observation parameters x count survey detection ####

##### Bias dot plot ####
d1 <- ggplot(rel.bias.few  %>% filter(variable %in% obs.pars), 
             aes(x = det.abund, y = bias, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  xlab('Count survey detection') + ylab('Relative bias') +
  facet_grid(dataset~trend, scales = 'free_x', labeller = label_wrap_gen()) +
  ylim(c(-0.3, 0.3)) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  theme_bw() +
  theme(legend.position = 'top',
        plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line")) +
  scale_color_manual(values = rainbow2[c(2,3)], name = '') +
  scale_shape_manual(values = c(15,16), name = '')
# d1

##### RMSE dot ####
d2 <- ggplot(rmse.few  %>% filter(variable %in% obs.pars), 
             aes(x = det.abund, y = rmse, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  xlab('Count survey detection') + ylab('RMSE') +
  facet_grid(dataset~trend, scales = 'free_x', labeller = label_wrap_gen()) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  theme_bw() +
  theme(legend.position = 'top',
        plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line")) +
  scale_color_manual(values = rainbow2[c(2,3)], name = '') +
  scale_shape_manual(values = c(15,16), name = '')
# d2

##### CV dot ####
d3 <- ggplot(cv.few  %>% filter(variable %in% obs.pars), 
             aes(x = det.abund, y = cv, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  xlab('Count detection survey') + ylab('Coefficient of variation (CV)') +
  facet_grid(dataset~trend, scales = 'free_x', labeller = label_wrap_gen()) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  theme_bw() +
  theme(legend.position = 'top',
        plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line")) +
  scale_color_manual(values = rainbow2[c(2,3)], name = '') +
  scale_shape_manual(values = c(15,16), name = '')
# d3

##### Combine ####
plot_grid(d1, d2, d3, nrow = 2, labels = "AUTO", align = "hv", label_size = 12)
# ggsave(width = 12, height = 15, here("figures", 'final', "fig6.png"))

#### Amanda is trying things ####
## curious about emphasizing model structure/datasets included

dataset.labs <- c("Full IPM", "Abundance & Prod.", "Abundance & Surv.", "Abundance Only")
names(dataset.labs) <- c("Full IPM", "Abundance & Productivity", "Abundance & Survival", "Abundance Only")
lambda.labs <- c("Decrease", "Stable", "Increase")
names(lambda.labs) <- c("Decreasing", "Stable", "Increasing")

### what about the 'power' of MR.det?
### and now trying to understand why bias in phi1 *increases* with increasing det.MR? 
test.bias <- rel.bias.sc  %>% 
  group_by(trend, life_hist, variable, det.MR, det.abund, det.prod, dataset) %>%
  dplyr::summarize(value = mean(bias)) %>%
  filter(trend == 'Stable') #%>%
  # filter(variable %nin% obs.pars) #%>%
  # filter(dataset == 'Full IPM')
  #phi1 has the most/only bias in full IPM - exploring
  # filter(dataset == 'Full IPM' & det.prod == 'High') #%>%
  #check this
  # transform(scenario = factor(scenario, levels = c(1,3,2),
  #                             labels = c('fast', 'mod', 'slow')))

#doesn't seem to vary over det.prod, so didn't visualize             
ggplot(test.bias %>% filter(dataset == 'Full IPM' & !is.na(det.prod) & !is.na(det.MR)), 
       aes(x = det.MR, y = value, col = variable, 
           group = variable)) +
  geom_point() + 
  geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  #ylim(c(-1.75, 1.75)) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  xlab('MR detection') + 
  ylab('Relative bias') +
  facet_nested(det.prod~det.abund + life_hist, scales = 'free_x')
  # facet_nested(scenario~det.abund, scales = 'free_x') +
  theme_bw()
  
  ggplot(test.bias %>% filter(dataset == 'Full IPM'  & !is.na(det.prod) & !is.na(det.MR)
                              & variable %nin% obs.pars), 
         aes(x = det.MR, y = value, col = variable, 
             group = variable)) +
    geom_point() + 
    geom_line() +
    geom_hline(aes(yintercept = 0), linetype = 'dotted') +
    #ylim(c(-1.75, 1.75)) +
    scale_x_discrete(labels = c("L", "M", "H")) +
    xlab('MR detection') + 
    ylab('Relative bias') +
    facet_nested(det.prod~det.abund + life_hist, scales = 'free_x')
  # facet_nested(scenario~det.abund, scales = 'free_x') +
  theme_bw()
  
  #no nests
  ggplot(test.bias %>% filter(dataset == 'Abundance & Survival'), 
         aes(x = det.MR, y = value, col = variable, 
             group = variable)) +
    geom_point() + 
    geom_line() +
    geom_hline(aes(yintercept = 0), linetype = 'dotted') +
    #ylim(c(-1.75, 1.75)) +
    scale_x_discrete(labels = c("L", "M", "H")) +
    xlab('MR detection') + 
    ylab('Relative bias') +
    facet_nested(det.abund~life_hist, scales = 'free_x')
  # facet_nested(scenario~det.abund, scales = 'free_x') +
  theme_bw()
  
  #no MR
  ggplot(test.bias %>% filter(dataset == 'Abundance & Productivity'), 
         aes(x = det.prod, y = value, col = variable, 
             group = variable)) +
    geom_point() + 
    geom_line() +
    geom_hline(aes(yintercept = 0), linetype = 'dotted') +
    #ylim(c(-1.75, 1.75)) +
    scale_x_discrete(labels = c("L", "M", "H")) +
    xlab('Nest detection') + 
    ylab('Relative bias') +
    facet_nested(det.abund~life_hist, scales = 'free_x')
  # facet_nested(scenario~det.abund, scales = 'free_x') +
  theme_bw()
  
  #abund only
  ggplot(test.bias %>% filter(dataset == 'Abundance Only'), 
         aes(x = det.abund, y = value, col = variable, 
             group = variable)) +
    geom_point() + 
    geom_line() +
    geom_hline(aes(yintercept = 0), linetype = 'dotted') +
    #ylim(c(-1.75, 1.75)) +
    scale_x_discrete(labels = c("L", "M", "H")) +
    xlab('Count detection') + 
    ylab('Relative bias') +
    facet_nested(. ~ life_hist, scales = 'free_x')
  # facet_nested(scenario~det.abund, scales = 'free_x') +
  theme_bw()

###what about not in the full IPM? -- didn't re-tool this code to match new stuff
# test.bias2 <- rel.bias.sc  %>% 
#   group_by(lambda.scenario, scenario, variable, det.MR, det.abund, det.prod, dataset) %>%
#   dplyr::summarize(value = mean(value)) %>%
#   #phi1 has the most/only bias in full IPM - exploring
#   filter(scenario == 2) 
# 
# ggplot(test.bias2 %>% filter(lambda.scenario == 'Stable') %>%
#          filter(det.prod == 'Medium'), 
#        aes(x = det.MR, y = value, col = factor(variable), 
#            group = factor(variable),
#            shape = factor(variable))) +
#   geom_point() + 
#   geom_line() +
#   geom_hline(aes(yintercept = 0), linetype = 'dotted') +
#   #ylim(c(-1.75, 1.75)) +
#   scale_x_discrete(labels = c("L", "M", "H")) +
#   xlab('MR detection') + 
#   ylab('Relative bias') +
#   # facet_nested(det.prod~det.abund + scenario, scales = 'free_x')
#   facet_nested(dataset~det.abund, scales = 'free') +
#   theme_bw()
# 
# test.cv2 <- cv.vals.sc  %>% 
#   group_by(lambda.scenario, scenario, variable, det.MR, det.abund, det.prod, dataset) %>%
#   dplyr::summarize(value = mean(value)) %>%
#   #phi1 has the most/only bias in full IPM - exploring
#   filter(scenario == 2) 
# 
# ggplot(test.cv2 %>% filter(lambda.scenario == 'Stable') %>%
#          filter(det.prod == 'Medium'), 
#        aes(x = det.MR, y = value, col = factor(variable), 
#            group = factor(variable),
#            shape = factor(variable))) +
#   geom_point() + 
#   geom_line() +
#   geom_hline(aes(yintercept = 0), linetype = 'dotted') +
#   #ylim(c(-1.75, 1.75)) +
#   scale_x_discrete(labels = c("L", "M", "H")) +
#   xlab('MR detection') + 
#   ylab('Relative bias') +
#   # facet_nested(det.prod~det.abund + scenario, scales = 'free_x')
#   facet_nested(dataset~det.abund, scales = 'free') +
#   theme_bw()

# test.rmse <- rmse.vals.sc  %>% 
#   group_by(lambda.scenario, scenario, variable, det.MR, det.abund, det.prod, dataset) %>%
#   dplyr::summarize(value = mean(value)) %>%
#   filter(dataset == 'Full IPM' & scenario == 2)
# 
# ggplot(test.rmse, 
#        aes(x = det.MR, y = value, col = factor(variable), group = factor(variable),
#            shape = factor(variable))) +
#   geom_point() + 
#   geom_line() +
#   geom_hline(aes(yintercept = 0), linetype = 'dotted') +
#   #facet_grid(dataset~lambda.scenario, scales = 'free_x', labeller = label_wrap_gen()) +
#   #ylim(c(-1.75, 1.75)) +
#   scale_x_discrete(labels = c("L", "M", "H")) +
#   xlab('MR detection') + 
#   ylab('RMSE') +
#   facet_nested(lambda.scenario + det.prod~det.abund, scales = 'free_x')
# 
# test.cv <- cv.vals.sc  %>% 
#   group_by(lambda.scenario, scenario, variable, det.MR, det.abund, det.prod, dataset) %>%
#   dplyr::summarize(value = mean(value)) %>%
#   filter(dataset == 'Full IPM' & scenario == 2)
# 
# ggplot(test.cv, 
#        aes(x = det.MR, y = value, col = factor(variable), group = factor(variable),
#            shape = factor(variable))) +
#   geom_point() + 
#   geom_line() +
#   geom_hline(aes(yintercept = 0), linetype = 'dotted') +
#   #facet_grid(dataset~lambda.scenario, scales = 'free_x', labeller = label_wrap_gen()) +
#   #ylim(c(-1.75, 1.75)) +
#   scale_x_discrete(labels = c("L", "M", "H")) +
#   xlab('MR detection') + 
#   ylab('CV') +
#   facet_nested(lambda.scenario + det.prod~det.abund, scales = 'free_x')

#
