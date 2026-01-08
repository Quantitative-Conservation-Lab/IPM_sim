## IPM Sim Figures
## Adapted from PaperFigures.RMD
## A DuVall
## 15 Oct 2021
## Updated by A Bratt
## 8 Jan 2026

##################### Libraries ###########################

library(tidyverse)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(coda)
library(captioner)
library(knitr)
library(reshape2)
library(here)
library(RColorBrewer)
library(colorspace)
library(ggh4x)
library(Hmisc)
library(patchwork)
library(tidybayes)

##################### Constants ###########################

# TODO - clean up globals

rainbow2 <- c("violetred4", "dodgerblue3", 'deepskyblue1', "#4aaaa5", "#a3d39c", "#f6b61c", "chocolate2", "red3")

scenarios <- readRDS(here("data", "demographic_scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename(
    "phi1" = "S.J", 
    "phiad" = "S.A",
    "fec" = "f"
  )

low.lam.params <- scenarios %>% 
  filter(trend == "decline")
med.lam.params <- scenarios %>% 
  filter(trend == "stable")
high.lam.params <- scenarios %>% 
  filter(trend == "increase")

det.abund <- factor(x = c("L", "M", "H"))
det.MR <- factor(x = c("L", "M", "H", "NA"))
det.prod <- factor(x = c("L", "M", "H", "NA"))
lambda <- factor(x = c("L", "M", "H"))
data_scenarios <- readRDS(here("data", "data_scenarios.RDS"))

##################### Data Prep ###########################

# TODO change file path for user
# NOTE - this is slow, files are large
row.low <- read.csv(file = "C:/Users/AbbyBratt/Desktop/IPM SIM results/lowout.csv")
row.med <- read.csv(file = "C:/Users/AbbyBratt/Desktop/IPM SIM results/medout.csv")
row.high <- read.csv(file = "C:/Users/AbbyBratt/Desktop/IPM SIM results/highout.csv")

low.dat <- row.low %>%
  dplyr::rename(phi1 = `mean.phi.1.`, phiad = `mean.phi.2.`) %>% 
  inner_join(low.lam.params %>% mutate(scenario = row_number()), 
             by = 'scenario', suffix = c('.obs', '.true')) %>%
  transform(lambda.scenario = 'Decreasing')

med.dat <- row.med %>%
  dplyr::rename(phi1 = `mean.phi.1.`, phiad = `mean.phi.2.`) %>% 
  inner_join(med.lam.params %>% mutate(scenario = row_number()), 
             by = 'scenario', suffix = c('.obs', '.true')) %>%
  transform(lambda.scenario = 'Stable')

high.dat <- row.high %>%
  dplyr::rename(phi1 = `mean.phi.1.`, phiad = `mean.phi.2.`) %>% 
  inner_join(high.lam.params %>% mutate(scenario = row_number()),
             by = 'scenario', suffix = c('.obs', '.true')) %>%
  transform(lambda.scenario = 'Increasing')

all.dat <- bind_rows(low.dat, med.dat, high.dat)

# summarize medians, sd, cvs, relative bias, and error at the model (sim) level
all.stats.sims <- all.dat %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  transform(psurv.true = ifelse(det.abund == 'L', 0.3, 
                                ifelse(det.abund == 'M', 0.5, 
                                       ifelse(det.abund == 'H', 0.8, NA)))) %>%
  transform(meanp.true = ifelse(det.MR == 'L', 0.3, 
                                ifelse(det.MR == 'M', 0.5, 
                                       ifelse(det.MR == 'H', 0.8, NA)))) %>%
  rename("psurv.obs" = "p.surv", "meanp.obs" = "mean.p") %>% 
  dplyr::select('phi1.obs', 'phi1.true', 'phiad.obs', 'phiad.true', 'fec.obs', 'fec.true',
                "psurv.obs", "psurv.true", "meanp.obs", "meanp.true",
                'sims', 'scenario', 'simscenarios', 'lambda.scenario') %>%
  group_by(lambda.scenario, scenario, simscenarios, sims) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(-c(lambda.scenario, scenario, simscenarios, sims, iter)) %>% 
  separate_wider_delim(name, delim = ".", names = c("param", "type")) %>% 
  pivot_wider(id_cols = c(lambda.scenario, scenario, simscenarios, sims, param, iter), 
              names_from = type, values_from = value) %>% 
  mutate(error = (obs - true)) %>% 
  group_by(lambda.scenario, scenario, simscenarios, sims, param) %>% 
  dplyr::summarize(
    median = median(obs), 
    sd = sd(obs), 
    cv = sd(obs)/mean(obs),
    rb = (mean(obs) - mean(true)) / mean(true), # NOTE using mean(true) tho it doesnt vary within group
    rmse = sqrt(mean(error^2))
  )

# average performance stats across sim replicates
all.stats <- all.stats.sims %>% 
  select(lambda.scenario, scenario, simscenarios, sims, param, cv, rb, rmse) %>% 
  group_by(lambda.scenario, scenario, simscenarios, param) %>% 
  dplyr::summarize(
    across(
      c(cv, rb, rmse),
      list(
        mean = ~ mean(.),
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
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()),
             by = "simscenarios") %>%
  dplyr::select('lambda.scenario', 'scenario', 'simscenarios',
                "param", "rb_mean",
                'det.MR', 'det.abund', 'det.prod') %>%
  reshape2::melt(id.vars = c('lambda.scenario', 'scenario', "simscenarios", "param", 'det.MR', 'det.abund', 'det.prod')) %>%
  group_by(lambda.scenario, det.MR, det.abund, det.prod, param) %>%
  dplyr::summarize(bias = mean(value), .groups = 'keep') %>% # taking average across scenario + simscenarios
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec', 'psurv', 'meanp'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity', 'Count survey detection', 'MR detection'))) %>%
  transform(lambda.scenario = factor(lambda.scenario,
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = ifelse(is.na(det.MR)&!is.na(det.prod), 'Abundance & Productivity', 
                             ifelse(!is.na(det.MR)&is.na(det.prod), 'Abundance & Survival',
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM'))))

# average over two layers of detection (det.MR and det.prod)
obs.pars <- c('MR detection', 'Count survey detection')
rel.bias.few <- rel.bias %>%
  group_by(variable, det.abund, lambda.scenario, dataset) %>%
  dplyr::summarize(bias = mean(bias), .groups = 'keep') %>% # taking average
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))

# just bias; same as above, keep 'scenario'
rel.bias.sc <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  filter(param %nin% c("meanp", "psurv")) %>% 
  dplyr::select('lambda.scenario', 'scenario', "simscenarios", "param", "rb_mean",
                'det.MR', 'det.abund', 'det.prod') %>%
  transform(lambda.scenario = factor(lambda.scenario,
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = ifelse(is.na(det.MR)&!is.na(det.prod), 'Abundance & Productivity', 
                             ifelse(!is.na(det.MR)&is.na(det.prod), 'Abundance & Survival',
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM')))) %>% 
  rename("variable" = 'param', "value" = "rb_mean")

# merge back onto bias df so can have something to name/view the 'scenarios' 
rel.bias.dem <- rel.bias.sc %>%
  inner_join(scenarios %>% 
               arrange(trend) %>% 
               mutate(scenario = c(1:3, 1:3, 1:3), #messy but works for now
                      lambda.scenario = case_when(
                        trend == "increase" ~ "Increasing",
                        trend == "decline" ~ "Decreasing",
                        TRUE ~ "Stable"
                      )) %>% 
               rename(phi1.true = phi1, 
                      phiad.true = phiad, 
                      fec.true = fec), 
             by = c("lambda.scenario", "scenario")) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity'))) %>% 
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(life_hist = factor(life_hist, levels = c("slow", "mod", "fast"), 
                               labels = c("Slow", "Moderate", "Fast"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))

# facet by both fec and juv true vals
plot.vals <- rel.bias.dem %>%
  group_by(variable, life_hist, lambda.scenario, dataset) %>% 
  dplyr::summarize(value = mean(value), .groups = 'keep') 

##################### RMSE ###########################

rmse.vals <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  dplyr::select(lambda.scenario, scenario, simscenarios, 
                param, rmse_mean,
                det.MR, det.abund, det.prod) %>%
  reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'simscenarios', 'param', 'det.MR', 'det.abund', 'det.prod')) %>%
  group_by(lambda.scenario, det.MR, det.abund, det.prod, param) %>%
  dplyr::summarize(mean.rmse = mean(value), .groups = 'keep') %>% 
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec', 
                                                   'psurv', 'meanp'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity',
                                         'Count survey detection', 'MR detection'))) %>%
  transform(lambda.scenario = factor(lambda.scenario,
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'),labels = c('Low', 'Medium', 'High'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = ifelse(is.na(det.MR)&!is.na(det.prod), 'Abundance & Productivity', 
                             ifelse(!is.na(det.MR)&is.na(det.prod), 'Abundance & Survival',
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM'))))

# average over two layers of detection (det.MR and det.prod)
rmse.few <- rmse.vals %>%
  group_by(variable, det.abund, lambda.scenario, dataset) %>%
  dplyr::summarize(rmse = mean(mean.rmse), .groups = 'keep') %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))

# data-generating values
rmse.vals.sc <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  filter(param %nin% c("meanp", "psurv")) %>% 
  dplyr::select(lambda.scenario, scenario, simscenarios, 
                param, rmse_mean, 
                det.MR, det.abund, det.prod) %>%
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity'))) %>%
  transform(lambda.scenario = factor(lambda.scenario,
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = ifelse(is.na(det.MR)&!is.na(det.prod), 'Abundance & Productivity', 
                             ifelse(!is.na(det.MR)&is.na(det.prod), 'Abundance & Survival',
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM')))) %>% 
  rename("value" = "rmse_mean")

rmse.dem <- rmse.vals.sc %>%
  inner_join(scenarios %>% 
               arrange(trend) %>% 
               mutate(scenario = c(1:3, 1:3, 1:3), 
                      lambda.scenario = case_when(
                        trend == "increase" ~ "Increasing",
                        trend == "decline" ~ "Decreasing",
                        TRUE ~ "Stable"
                      )) %>% 
               rename(phi1.true = phi1, 
                      phiad.true = phiad, 
                      fec.true = fec), 
             by = c("lambda.scenario", "scenario")) %>% 
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(life_hist = factor(life_hist, levels = c("slow", "mod", "fast"), 
                               labels = c("Slow", "Moderate", "Fast"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))


# facet by both fec and juv true vals
plot.vals.rmse <- rmse.dem %>%
  group_by(variable, life_hist, lambda.scenario, dataset) %>% 
  dplyr::summarize(value = mean(value), .groups = 'keep')

##################### CV ###########################

cv.vals <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  dplyr::select(lambda.scenario, scenario, simscenarios, 
                param, cv_mean,
                det.MR, det.abund, det.prod) %>%
  reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'simscenarios', 'param', 'det.MR', 'det.abund', 'det.prod')) %>%
  group_by(lambda.scenario, det.MR, det.abund, det.prod, param) %>%
  dplyr::summarize(mean.cv = mean(value), .groups = 'keep') %>% 
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec', 
                                                   'psurv', 'meanp'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity',
                                         'Count survey detection', 'MR detection'))) %>%
  transform(lambda.scenario = factor(lambda.scenario,
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'),labels = c('Low', 'Medium', 'High'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = ifelse(is.na(det.MR)&!is.na(det.prod), 'Abundance & Productivity', 
                             ifelse(!is.na(det.MR)&is.na(det.prod), 'Abundance & Survival',
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM'))))

# average over two layers of detection (det.MR and det.prod)
cv.few <- cv.vals %>%
  group_by(variable, det.abund, lambda.scenario, dataset) %>%
  dplyr::summarize(cv = mean(mean.cv), .groups = 'keep') %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))

# data-generating values
cv.vals.sc <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  filter(param %nin% c("meanp", "psurv")) %>% 
  dplyr::select(lambda.scenario, scenario, simscenarios, 
                param, cv_mean, 
                det.MR, det.abund, det.prod) %>%
  rename(variable = param) %>% 
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity'))) %>%
  transform(lambda.scenario = factor(lambda.scenario,
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(missing.MR = ifelse(is.na(det.MR), 1, 0),
            missing.prod = ifelse(is.na(det.prod), 1, 0)) %>%
  transform(num.miss = missing.MR + missing.prod) %>%
  transform(dataset = ifelse(is.na(det.MR)&!is.na(det.prod), 'Abundance & Productivity', 
                             ifelse(!is.na(det.MR)&is.na(det.prod), 'Abundance & Survival',
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM')))) %>% 
  rename("value" = "cv_mean")

cv.dem <- cv.vals.sc %>%
  inner_join(scenarios %>% 
               arrange(trend) %>% 
               mutate(scenario = c(1:3, 1:3, 1:3), #janky but try for now 
                      lambda.scenario = case_when(
                        trend == "increase" ~ "Increasing",
                        trend == "decline" ~ "Decreasing",
                        TRUE ~ "Stable"
                      )) %>% 
               rename(phi1.true = phi1, 
                      phiad.true = phiad, 
                      fec.true = fec), 
             by = c("lambda.scenario", "scenario")) %>% 
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(life_hist = factor(life_hist, levels = c("slow", "mod", "fast"), 
                               labels = c("Slow", "Moderate", "Fast"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))

# facet by both fec and juv true vals
plot.vals.cv <- cv.dem %>%
  group_by(variable, life_hist, lambda.scenario, dataset) %>% 
  dplyr::summarize(value = mean(value), .groups = 'keep')

##################### LAMBDA #########################

# TODO change this for each user
row.low <- read_csv(file = "C:/Users/AbbyBratt/Desktop/IPM SIM results/row_low.csv") %>% 
  transform(lambda.scenario = 'Decreasing')
row.med <- read_csv(file = "C:/Users/AbbyBratt/Desktop/IPM SIM results/row_med.csv") %>%
  transform(lambda.scenario = 'Stable')
row.high <- read_csv(file = "C:/Users/AbbyBratt/Desktop/IPM SIM results/row_high.csv") %>% 
  transform(lambda.scenario = 'Increasing')

# Reformat for plotting
toplot1 <- row.low %>%
  select(contains("geomean"), scenario, sims, simscenarios, Quantile, lambda.scenario) %>%
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(scenario, sims, simscenarios, Year, Quantile) %>% 
  ungroup() %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()),
             by = "simscenarios") %>%
  inner_join(scenarios %>% 
               arrange(trend) %>% 
               mutate(scenario = c(1:3, 1:3, 1:3), 
                      lambda.scenario = case_when(
                        trend == "increase" ~ "Increasing",
                        trend == "decline" ~ "Decreasing",
                        TRUE ~ "Stable"
                      )) %>% 
               rename(phi1.true = phi1, 
                      phiad.true = phiad, 
                      fec.true = fec), 
             by = c("lambda.scenario", "scenario")) %>% 
  mutate(scenario = as.factor(scenario),
         sims = as.factor(sims), 
         simscenarios = as.factor(simscenarios))

toplot2 <- row.med %>%
  select(contains("geomean"), scenario, sims, simscenarios, Quantile, lambda.scenario) %>%
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(scenario, sims, simscenarios, Year, Quantile) %>% 
  ungroup() %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()),
             by = "simscenarios") %>%
  inner_join(scenarios %>% 
               arrange(trend) %>% 
               mutate(scenario = c(1:3, 1:3, 1:3),
                      lambda.scenario = case_when(
                        trend == "increase" ~ "Increasing",
                        trend == "decline" ~ "Decreasing",
                        TRUE ~ "Stable"
                      )) %>% 
               rename(phi1.true = phi1, 
                      phiad.true = phiad, 
                      fec.true = fec), 
             by = c("lambda.scenario", "scenario")) %>% 
  mutate(scenario = as.factor(scenario),
         sims = as.factor(sims), 
         simscenarios = as.factor(simscenarios))

toplot3 <- row.high %>%
  select(contains("geomean"), scenario, sims, simscenarios, Quantile, lambda.scenario) %>%
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(scenario, sims, simscenarios, Year, Quantile) %>% 
  ungroup() %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()),
             by = "simscenarios") %>%
  inner_join(scenarios %>% 
               arrange(trend) %>% 
               mutate(scenario = c(1:3, 1:3, 1:3), 
                      lambda.scenario = case_when(
                        trend == "increase" ~ "Increasing",
                        trend == "decline" ~ "Decreasing",
                        TRUE ~ "Stable"
                      )) %>% 
               rename(phi1.true = phi1, 
                      phiad.true = phiad, 
                      fec.true = fec), 
             by = c("lambda.scenario", "scenario")) %>% 
  mutate(scenario = as.factor(scenario),
         sims = as.factor(sims), 
         simscenarios = as.factor(simscenarios))

toplot <- bind_rows(toplot1, toplot2, toplot3) %>% 
  mutate(
    det.MR = na_if(as.character(det.MR), "NA"), 
    det.prod = na_if(as.character(det.prod), "NA"),
    det.abund = na_if(as.character(det.abund), "NA"),
    lambda = factor(lambda)
  ) %>% 
  transform(dataset = ifelse(is.na(det.MR)&!is.na(det.prod), 'Abundance & Productivity', 
                             ifelse(!is.na(det.MR)&is.na(det.prod), 'Abundance & Survival',
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM')))) %>%
  group_by(Quantile, Year, det.abund, lambda, lambda.scenario, life_hist, dataset, det.MR, det.prod) %>% 
  dplyr::summarize(value = mean(value), .groups = "drop") %>% 
  ungroup() %>% 
  mutate(Quantile = str_remove(Quantile, "\\%")) %>% 
  mutate(Quantile = paste("X", Quantile, sep = "")) %>% 
  reshape2::dcast(dataset + Year +  det.MR + det.prod + det.abund + lambda + lambda.scenario + life_hist ~ Quantile, value.var = "value") %>% 
  mutate(Year = Year + 1) %>% 
  filter(Year %in% c(15)) %>% 
  mutate(Year = factor(Year)) %>% 
  mutate(det.abund = factor(det.abund, levels = c("L", "M", "H"), labels = c("Low", "Medium", "High"))) %>% 
  mutate(det.prod = factor(det.prod, levels = c("L", "M", "H"), labels = c("Low", "Medium", "High"))) %>% 
  mutate(det.MR = factor(det.MR, levels = c("L", "M", "H"), labels = c("Low", "Medium", "High"))) %>% 
  transform(lambda.scenario = factor(lambda,
                                     labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'))) 

# average over two layers of detection (det.MR and det.prod) and life history type
toplot.few <- toplot %>%
  group_by(dataset, Year, det.abund, lambda, lambda.scenario) %>%
  dplyr::summarize(
    `X2.5` = mean(`X2.5`), 
    `X50` = mean(`X50`),
    `X97.5` = mean(`X97.5`),
    .groups = 'drop') %>%
  transform(lambda.scenario = factor(lambda.scenario,
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))

# average over all layers of detection and NOT life history type
toplot.lh <- toplot %>%
  group_by(dataset, Year, life_hist, lambda, lambda.scenario) %>%
  dplyr::summarize(
    `X2.5` = mean(`X2.5`), 
    `X50` = mean(`X50`),
    `X97.5` = mean(`X97.5`),
    .groups = 'drop') %>%
  transform(lambda.scenario = factor(lambda.scenario,
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(life_hist = factor(life_hist, levels = c("slow", "mod", "fast"), 
                               labels = c("Slow", "Moderate", "Fast"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))

################### FIGURES ##########################

# TODO might need some adjustments

#### Figure 3: RMSE and bias ecological paramters x count survey detection ####
dataset.labs <- c("Full IPM", "Abundance & Prod.", "Abundance & Surv.", "Abundance Only")
names(dataset.labs) <- c("Full IPM", "Abundance & Productivity", "Abundance & Survival", "Abundance Only")
lambda.labs <- c("Decrease", "Stable", "Increase")
names(lambda.labs) <- c("Decreasing", "Stable", "Increasing")

##### Bias dot plot ####
a1 <- ggplot(rel.bias.few  %>% filter(variable %nin% obs.pars), 
             aes(x = det.abund, y = bias, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + 
  geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  #xlab('Count survey detection') + ylab('Relative bias') +
  #facet_grid(dataset~lambda.scenario, scales = 'free_x', labeller = label_wrap_gen()) +
  #ylim(c(-1.75, 1.75)) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  xlab('Count survey detection') + 
  ylab('Relative bias') +
  facet_grid(dataset~lambda.scenario, scales = 'free_x', 
             labeller = labeller(dataset = dataset.labs, lambda.scenario = lambda.labs)) +
  scale_y_continuous(limits = c(-1.2, 1.2), breaks = c(-1,0,1)) +
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


##### RMSE heat map ####
a2 <- ggplot(rmse.few %>% filter(variable %nin% obs.pars), aes(x = factor(det.abund), y = variable, fill = rmse)) +
  geom_tile(color = 'grey50') +
  xlab('Count survey detection') + ylab('') +
  #facet_grid(dataset ~ lambda.scenario, drop = T, scales = 'free_x', labeller = label_wrap_gen()) +
  facet_grid(dataset ~ lambda.scenario, drop = T, scales = 'free_x', 
             labeller = labeller(dataset = dataset.labs, 
                                 lambda.scenario = lambda.labs)) +
  scale_fill_gradient2(name = "RMSE",
                       #mid = "white", high = rainbow2[2], midpoint = 0) +
                       #low = "white", mid = rainbow2[3], high = rainbow2[2]) + #,
                       low = "white", high = rainbow2[2]) + #,
  #midpoint = 0.5) + # TODO - note change here 
  theme_light() +
  scale_y_discrete(labels = c(expression(phi["2"]), expression(phi["1"]), expression(f))) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  theme(legend.position = 'top',
        #legend.title = element_text(size = 11, vjust = 0.75),
        legend.title = element_text(size = 10, vjust = 0.75),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 10, vjust = 0.75),
        axis.title = element_text(size = 10, vjust = 0.75),
        strip.text = element_text(color = "black", size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line"))
a2

##### CV heat map ####
a3 <- ggplot(cv.few %>% filter(variable %nin% obs.pars), aes(x = factor(det.abund), y = variable, fill = cv)) +
  geom_tile(color = 'grey50') +
  xlab('Count survey detection') + ylab('') +
  #facet_grid(dataset ~ lambda.scenario, drop = T, scales = 'free_x', labeller = label_wrap_gen()) +
  facet_grid(dataset ~ lambda.scenario, drop = T, scales = 'free_x', 
             labeller = labeller(dataset = dataset.labs, 
                                 lambda.scenario = lambda.labs)) +
  scale_fill_gradient2(name = "CV",
                       #mid = "white", high = rainbow2[2], midpoint = 0) +
                       low = "white", high = rainbow2[2]) + #,
  #midpoint = 0.5) + # TODO - note change here 
  theme_light() +
  scale_y_discrete(labels = c(expression(phi["2"]), expression(phi["1"]), expression(f))) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  theme(legend.position = 'top',
        #legend.title = element_text(size = 11, vjust = 0.75),
        legend.title = element_text(size = 10, vjust = 0.75),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 10, vjust = 0.75),
        axis.title = element_text(size = 10, vjust = 0.75),
        strip.text = element_text(color = "black", size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line"))
a3

##### Combine ####
plot_grid(a1, a2, a3, nrow = 3, labels = "AUTO",
          align = "hv", label_size = 12)
ggsave(width = 6.5, height = 18, here("figures", "fig3_NEW.png"))


#### Figure 4: RMSE and bias ecological parameters x true fecundity ####

phi1_cat_lab <- c("True first-year\nφ: Low", 
                  "True first-year\nφ: Med",
                  "True first-year\nφ: High")

names(phi1_cat_lab) <- c("True first-year survival: Low", 
                         "True first-year survival: Medium", 
                         "True first-year survival: High")

##### Bias dot plot ####
b1 <- ggplot(plot.vals, aes(x = life_hist, y = value, col = factor(variable), group = factor(variable),
                            shape = factor(variable))) +
  geom_point() + geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  #scale_x_discrete(labels = c("L", "M", "H")) +
  xlab('Life History Type') + ylab('Relative bias') +
  facet_grid(dataset~lambda.scenario, scales = 'free_x') + #, 
  #labeller = labeller(dataset = dataset.labs, phi1_cat = phi1_cat_lab)) +
  #labeller = label_wrap_gen()) +
  #ylim(c(-1.75, 1.75)) +
  scale_y_continuous(limits = c(-1.75, 1.75), breaks = c(-1.5,0,1.5)) +
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
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #+
#scale_color_manual(values = rainbow2[-c(1,4)], name = '') +
#scale_shape_manual(values = c(15, 16, 17), name = '')
b1

##### RMSE heat map ####
b2 <- ggplot(plot.vals.rmse, aes(x = life_hist, y = factor(variable), fill = value)) +
  geom_tile(color = 'grey50') +
  xlab('Life History Type') + ylab('') +
  facet_grid(dataset~lambda.scenario, drop = T, scales = 'free_x'#, 
             #labeller = labeller(dataset = dataset.labs, phi1_cat = phi1_cat_lab)
  ) + #label_wrap_gen()) +
  #scale_fill_gradient2(name = "RMSE", mid = "white", high = rainbow2[2], midpoint = 0) +
  scale_fill_gradient2(name = "RMSE", low = "white", high = rainbow2[2],
                       n.breaks = 4) +
  #midpoint = 1, limits = c(0, 2), breaks = c(0, 0.5, 1, 1.5, 2)) +
  theme_light() +
  scale_y_discrete(labels = c(expression(φ["2"]), expression(φ["1"]), expression(f))) +
  #scale_x_discrete(labels = c("L", "M", "H")) +
  theme(legend.position = 'top',
        legend.title = element_text(size = 10, vjust = 0.75),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10, vjust = 0.75),
        axis.title = element_text(size = 10, vjust = 0.75),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line"))
b2

##### CV heat map ####
b3 <- ggplot(plot.vals.cv, aes(x = life_hist, y = factor(variable), fill = value)) +
  geom_tile(color = 'grey50') +
  xlab('Life History Type') + ylab('') +
  facet_grid(dataset~lambda.scenario, drop = T, scales = 'free_x'#, 
             #labeller = labeller(dataset = dataset.labs, phi1_cat = phi1_cat_lab)
  ) + #label_wrap_gen()) +
  #scale_fill_gradient2(name = "RMSE", mid = "white", high = rainbow2[2], midpoint = 0) +
  scale_fill_gradient2(name = "CV", low = "white", high = rainbow2[2]) + #,
  # midpoint = 0.5) +
  theme_light() +
  scale_y_discrete(labels = c(expression(φ["2"]), expression(φ["1"]), expression(f))) +
  #scale_x_discrete(labels = c("L", "M", "H")) +
  theme(legend.position = 'top',
        legend.title = element_text(size = 10, vjust = 0.75),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10, vjust = 0.75),
        axis.title = element_text(size = 10, vjust = 0.75),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line"))
b3

##### Combine ####
plot_grid(b1, b2, b3, nrow = 3, labels = "AUTO", align = "hv", label_size = 12)
ggsave(width = 6.5, height = 18, here("figures", "fig4_NEW.png"))


#### Figure 5: Lambda trends ####

##### WRT count survey detection ####
c1 <- ggplot(toplot.few) +
  geom_point(aes(x = Year, y = X50, col = det.abund, group = det.abund, shape = det.abund), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(x = Year, ymin = X2.5, ymax = X97.5, col = det.abund, group = det.abund,
                     shape = det.abund), position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = as.numeric(as.character(lambda))), linetype = 'dotted') +
  geom_hline(aes(yintercept = 1.0), linetype = 'solid') +
  xlab('Final year (t=15)') + #renamed the axis - TODO is this correct
  # oh right, yes, because we should have stopped the model after each year 
  # (to 'hide' the full time series from being used to estimate trend)
  ylab(expression(lambda)) + 
  facet_grid(dataset ~ lambda.scenario, scales = 'free', labeller = label_wrap_gen()) +
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
c2 <- ggplot(toplot.lh) +
  geom_point(aes(x = Year, y = X50, col = life_hist, group = life_hist, shape = life_hist), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(x = Year, ymin = X2.5, ymax = X97.5, col = life_hist, group = life_hist,
                     shape = life_hist), position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = as.numeric(as.character(lambda))), linetype = 'dotted') +
  geom_hline(aes(yintercept = 1.0), linetype = 'solid') +
  xlab('Final year (t=15)') + #renamed the axis - TODO is this correct
  # oh right, yes, because we should have stopped the model after each year 
  # (to 'hide' the full time series from being used to estimate trend)
  ylab(expression(lambda)) + 
  facet_grid(dataset ~ lambda.scenario, scales = 'free', labeller = label_wrap_gen()) +
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
ggsave(width = 15, height = 6, here("figures", "fig5_NEW.png"))

################### Appendix #########################

#### Figure 6: RMSE and bias observation parameters x count survey detection ####

##### Bias dot plot ####
d1 <- ggplot(rel.bias.few  %>% filter(variable %in% obs.pars), 
             aes(x = det.abund, y = bias, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  xlab('Count survey detection') + ylab('Relative bias') +
  facet_grid(dataset~lambda.scenario, scales = 'free_x', labeller = label_wrap_gen()) +
  ylim(c(-0.3, 0.3)) +
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
d1


##### RMSE heat map ####
d2 <- ggplot(rmse.few %>% filter(variable %in% obs.pars), aes(x = factor(det.abund), y = variable, fill = rmse)) +
  geom_tile(color = 'grey50') +
  xlab('Count survey detection') + ylab('') +
  facet_grid(dataset ~ lambda.scenario, drop = T, scales = 'free_x', labeller = label_wrap_gen()) +
  scale_fill_gradient2(name = "RMSE",
                       mid = "white", high = rainbow2[2], n.breaks = 3) +
  # breaks = c(-0.2, 0, 0.2), n.breaks = 3, labels = c("-0.2", "0", "0.2")) +
  theme_light() +
  theme(legend.position = 'top',
        legend.title = element_text(size = 11, vjust = 0.75),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line"))
d2

##### CV heat map ####
d3 <- ggplot(cv.few %>% filter(variable %in% obs.pars), aes(x = factor(det.abund), y = variable, fill = cv)) +
  geom_tile(color = 'grey50') +
  xlab('Count survey detection') + ylab('') +
  facet_grid(dataset ~ lambda.scenario, drop = T, scales = 'free_x', labeller = label_wrap_gen()) +
  scale_fill_gradient2(name = "CV",
                       mid = "white", high = rainbow2[2], n.breaks = 4) +
  # breaks = c(-0.2, 0, 0.2), n.breaks = 3, labels = c("-0.2", "0", "0.2")) +
  theme_light() +
  theme(legend.position = 'top',
        legend.title = element_text(size = 11, vjust = 0.75),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line"))
d3

##### Combine ####
plot_grid(d1, d2, d3, nrow = 3, labels = "AUTO", align = "hv", label_size = 12)
ggsave(width = 6.5, height = 18, here("figures", "fig6_NEW.png"))
