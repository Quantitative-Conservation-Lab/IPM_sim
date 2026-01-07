## IPM Sim Figures
## Adapted from PaperFigures.RMD
## A DuVall
## 15 Oct 2021

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

# TODO ab to tidy comments

rainbow2 <- c("violetred4", "dodgerblue3", 'deepskyblue1', "#4aaaa5", "#a3d39c", "#f6b61c", "chocolate2", "red3")

#scenarios <- read.csv(file = here::here('data', 'scenario_ID.csv'), header = T, stringsAsFactors = F)
# TODO - adjust as in previous file
scenarios <- readRDS(here("data", "demographic_scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename(
    "phi1" = "S.J", 
    "phiad" = "S.A",
    "fec" = "f"
  )

# low.params <- readRDS(file = here::here('data', 'low.lam.params.rds')) %>%
#   transform(scenario = 1:25)
# 
# med.params <- readRDS(file = here::here('data', 'med.lam.params.rds')) %>%
#   transform(scenario = 1:25)
# 
# high.params <- readRDS(file = here::here('data', 'high.lam.params.rds')) %>%
#   transform(scenario = 1:25)

low.lam.params <- scenarios %>% 
  filter(trend == "decline")
med.lam.params <- scenarios %>% 
  filter(trend == "stable")
high.lam.params <- scenarios %>% 
  filter(trend == "increase")

# TODO change file path for user
# NOTE - this is slow, files are large
row.low <- read.csv(file = "C:/Users/AbbyBratt/Desktop/IPM SIM results/lowout.csv")
row.med <- read.csv(file = "C:/Users/AbbyBratt/Desktop/IPM SIM results/medout.csv")
row.high <- read.csv(file = "C:/Users/AbbyBratt/Desktop/IPM SIM results/highout.csv")

# TODO where does this get made
#combine to get true values integrated in results
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

# all.meds <- all.dat %>%
#   dplyr::select('phi1.obs', 'phi1.true', 'phiad.obs', 'phiad.true', 'fec.obs', 'fec.true',
#                 "p.surv", "mean.p",
#                 'sims', 'scenario', 'simscenarios', 'lambda.scenario') %>%
#   reshape2::melt(id.vars = c('lambda.scenario', 'sims', 'scenario', 'simscenarios')) %>%
#   group_by(lambda.scenario, scenario, simscenarios, variable) %>%
#   dplyr::summarize(median = median(value), .groups = 'keep') %>%
#   reshape2::dcast(lambda.scenario + scenario + simscenarios ~ variable, value.var = 'median') 
# NOTE stopping here for now - returning to manuscript figs

# first chunk in figures.Rmd - move over here

# all.meds <- read.csv(file = here::here('figures', 'Processed csvs', 'all.meds.csv'), 
#                      header = T, stringsAsFactors = F)

# summarize
  # medians
  # sds
  # cvs
  # error - do this separately
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
# ok - so this has the stats per model run (i.e., replicate/sim)

# now we want to get these things averaged across replicates 
# (so, average cv, average rb, average rmse) 
# also keeping intervals for these 
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

# Now we want to grab the average RB, RMSE, and CV across replicates

######################################################
##################### Bias ###########################
######################################################

# TODO - somewhere here we need to add precisions

# TODO remove hardcoding on det.abund etc and make global constants file
det.abund <- factor(x = c("L", "M", "H"))
det.MR <- factor(x = c("L", "M", "H", "NA"))
det.prod <- factor(x = c("L", "M", "H", "NA"))
lambda <- factor(x = c("L", "M", "H"))
data_scenarios <- readRDS(here("data", "data_scenarios.RDS"))

# rel.bias
rel.bias <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()),
             by = "simscenarios") %>%
  # transform(p.surv.true = ifelse(det.abund == 'L', 0.3, 
  #                                ifelse(det.abund == 'M', 0.5, 
  #                                       ifelse(det.abund == 'H', 0.8, NA)))) %>%
  # transform(mean.p.true = ifelse(det.MR == 'L', 0.3, 
  #                                ifelse(det.MR == 'M', 0.5, 
  #                                       ifelse(det.MR == 'H', 0.8, NA)))) %>%
  #calculate relative bias, mean across scenarios
  # transform(phi1.bias = (phi1.obs-phi1.true)/phi1.true,
  #           phiad.bias = (phiad.obs-phiad.true)/phiad.true,
  #           fec.bias = (fec.obs-fec.true)/fec.true,
  #           p.surv.bias = (p.surv-p.surv.true)/p.surv.true,
  #           mean.p.bias = (mean.p-mean.p.true)/mean.p.true) %>%
  dplyr::select('lambda.scenario', 'scenario', 'simscenarios',
                "param", "rb_mean",
                'det.MR', 'det.abund', 'det.prod') %>%
  reshape2::melt(id.vars = c('lambda.scenario', 'scenario', "simscenarios", "param", 'det.MR', 'det.abund', 'det.prod')) %>%
  group_by(lambda.scenario, det.MR, det.abund, det.prod, param) %>%
  dplyr::summarize(bias = mean(value), .groups = 'keep') %>%
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

## average over two layers of detection (det.MR and det.prod)
obs.pars <- c('MR detection', 'Count survey detection')
rel.bias.few <- rel.bias %>%
  group_by(variable, det.abund, lambda.scenario, dataset) %>%
  dplyr::summarize(bias = mean(bias), .groups = 'keep') %>%
  transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))
  # transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
  #                            labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))

# TODO - go back to figures.RMD to make this object
# all.meds.sc <- read.csv(file = here::here('figures', 'Processed csvs', 'all.meds.sc.csv'), 
#                         header = T, stringsAsFactors = F)
#get posterior medians; same as above except keep 'scenario' in group_by
# all.meds.sc <- all.dat %>%
#   dplyr::select('phi1.obs', 'phi1.true', 'phiad.obs', 'phiad.true', 'fec.obs', 'fec.true',
#                 'sims', 'scenario', 'simscenarios', 'lambda.scenario') %>%
#   reshape2::melt(id.vars = c('sims', 'scenario', 'simscenarios', 'lambda.scenario')) %>%
#   group_by(lambda.scenario, scenario, simscenarios, variable) %>%
#   dplyr::summarize(median = median(value)) %>%
#   reshape2::dcast(lambda.scenario + scenario + simscenarios ~ variable, value.var = 'median') %>%
#   #calculate relative bias
#   transform(phi1.bias = (phi1.obs-phi1.true)/phi1.true,
#             phiad.bias = (phiad.obs-phiad.true)/phiad.true,
#             fec.bias = (fec.obs-fec.true)/fec.true)
# FROM HERE - back to manuscript figs

## just bias; same as above, keep 'scenario'
#rel.bias.sc <- all_meds_sc %>%
rel.bias.sc <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  filter(param %nin% c("meanp", "psurv")) %>% 
  # transform(phi1.bias = (phi1.obs-phi1.true)/phi1.true,
  #           phiad.bias = (phiad.obs-phiad.true)/phiad.true,
  #           fec.bias = (fec.obs-fec.true)/fec.true) %>%
  dplyr::select('lambda.scenario', 'scenario', "simscenarios", "param", "rb_mean",
                'det.MR', 'det.abund', 'det.prod') %>%
  #reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'det.MR', 'det.abund', 'det.prod')) %>%
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

# TODO - AEB doesnt think we need this section here because we have fewer scenario ####
#load saved true values
# true_vals <- read.csv(file = here::here('data', 'true.vals.csv'), header = T, stringsAsFactors = F)
# 
# #for categories: 
# fec_lims <- quantile(true_vals$fec.true, probs = c(0.33, 0.7), names = F)
# phiad_lims <- quantile(true_vals$phiad.true, probs = c(0.33, 0.7), names = F)
# phi1_lims <- quantile(true_vals$phi1.true, probs = c(0.33, 0.7), names = F)
#### ----

#merge back onto bias df so can have something to name/view the 'scenarios' 
rel.bias.dem <- rel.bias.sc %>%
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
  #merge(true_vals, by = c('lambda.scenario', 'scenario')) %>%
  transform(variable = factor(variable, levels = c('phiad', 'phi1', 'fec'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity'))) %>% 
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  
  # TODO - these should now be fast slow mod life histories
  # transform(fec_cat = ifelse(fec.true < fec_lims[1], 'L', ifelse(fec.true > fec_lims[2], 'H', 'M')),
  #           phiad_cat = ifelse(phiad.true < phiad_lims[1], 'L', ifelse(phiad.true > phiad_lims[2], 'H', 'M')),
  #           phi1_cat = ifelse(phi1.true < phi1_lims[1], 'L', ifelse(phi1.true > phi1_lims[2], 'H', 'M'))) %>%
  # transform(fec_cat = factor(fec_cat, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High')),
  #           phiad_cat = factor(phiad_cat, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High')),
  #           phi1_cat = factor(phi1_cat, levels = c('L', 'M', 'H'),
  #                             labels = c('True first-year survival: Low', 
  #                                        'True first-year survival: Medium', 
  #                                        'True first-year survival: High'))) %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(life_hist = factor(life_hist, levels = c("slow", "mod", "fast"), 
                                     labels = c("Slow", "Moderate", "Fast"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))
# transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
#                            labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))

#facet by both fec and juv true vals
plot.vals <- rel.bias.dem %>%
  group_by(variable, life_hist, lambda.scenario, dataset) %>% # TODO changed here
  dplyr::summarize(value = mean(value), .groups = 'keep') 



######################################################
##################### RMSE ###########################
######################################################

# TODO - repeat changes above here. LATER

rmse.vals <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  # transform(p.surv.true = ifelse(det.abund == 'L', 0.3, 
  #                                ifelse(det.abund == 'M', 0.5, 
  #                                       ifelse(det.abund == 'H', 0.8, NA)))) %>%
  # transform(mean.p.true = ifelse(det.MR == 'L', 0.3, 
  #                                ifelse(det.MR == 'M', 0.5, 
  #                                       ifelse(det.MR == 'H', 0.8, NA)))) %>%
  # transform(fec.rmse = (fec.obs-fec.true)^2,
  #           phiad.rmse = (phiad.obs-phiad.true)^2,
  #           phi1.rmse = (phi1.obs-phi1.true)^2,
  #           mean.p.rmse = (mean.p-mean.p.true)^2,
  #           p.surv.rmse = (p.surv-p.surv.true)^2) %>%
  dplyr::select(lambda.scenario, scenario, simscenarios, 
                param, rmse_mean,
                det.MR, det.abund, det.prod) %>%
  reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'simscenarios', 'param', 'det.MR', 'det.abund', 'det.prod')) %>%
  group_by(lambda.scenario, det.MR, det.abund, det.prod, param) %>%
  dplyr::summarize(mean.rmse = mean(value), .groups = 'keep') %>% # TODO - added sqrt, check
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

###average over two layers of detection (det.MR and det.prod)
rmse.few <- rmse.vals %>%
  group_by(variable, det.abund, lambda.scenario, dataset) %>%
  dplyr::summarize(rmse = mean(mean.rmse), .groups = 'keep') %>%
  #transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))
# transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
#                            labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))

########data-generating values
rmse.vals.sc <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  filter(param %nin% c("meanp", "psurv")) %>% 
  # transform(fec.rmse = (fec.obs-fec.true)^2,
  #           phiad.rmse = (phiad.obs-phiad.true)^2,
  #           phi1.rmse = (phi1.obs-phi1.true)^2) %>%
  dplyr::select(lambda.scenario, scenario, simscenarios, 
                param, rmse_mean, 
                det.MR, det.abund, det.prod) %>%
  #reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'det.MR', 'det.abund', 'det.prod')) %>%
  rename(variable = param) %>% 
  #group_by(lambda.scenario, scenario, simscenarios, det.MR, det.abund, det.prod, variable) %>%
  #dplyr::summarize(mean.rmse = mean(value), .groups = 'keep') %>% # TODO changed here, check
  #dplyr::summarize(cv = sd(value)/value, .groups = 'keep') %>% # CV version
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
  # merge(true_vals, by = c('lambda.scenario', 'scenario')) %>%
  # transform(variable = factor(variable, levels = c('phiad.rmse', 'phi1.rmse', 'fec.rmse'),
  #                             labels = c('Adult survival', 'Juv survival', 'Fecundity'))) %>% 
  # transform(lambda.scenario = factor(lambda.scenario, 
  #                                    levels = c("Decreasing", "Stable", "Increasing"))) %>%
  # transform(fec_cat = ifelse(fec.true < fec_lims[1], 'L', ifelse(fec.true > fec_lims[2], 'H', 'M')),
  #           phiad_cat = ifelse(phiad.true < phiad_lims[1], 'L', ifelse(phiad.true > phiad_lims[2], 'H', 'M')),
  #           phi1_cat = ifelse(phi1.true < phi1_lims[1], 'L', ifelse(phi1.true > phi1_lims[2], 'H', 'M'))) %>%
  # transform(fec_cat = factor(fec_cat, levels = c('L', 'M', 'H'),  labels = c('Low', 'Medium', 'High')),
  #           phiad_cat = factor(phiad_cat, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High')),
  #           phi1_cat = factor(phi1_cat, levels = c('L', 'M', 'H'),
  #                             labels = c('True first-year survival: Low', 
  #                                        'True first-year survival: Medium', 
  #                                        'True first-year survival: High'))) %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(life_hist = factor(life_hist, levels = c("slow", "mod", "fast"), 
                               labels = c("Slow", "Moderate", "Fast"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))
# transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
#                            labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))


#facet by both fec and juv true vals
plot.vals.rmse <- rmse.dem %>%
  group_by(variable, life_hist, lambda.scenario, dataset) %>% # TODO changed here
  dplyr::summarize(value = mean(value), .groups = 'keep')

######################################################
##################### CV ###########################
######################################################

cv.vals <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  # transform(p.surv.true = ifelse(det.abund == 'L', 0.3, 
  #                                ifelse(det.abund == 'M', 0.5, 
  #                                       ifelse(det.abund == 'H', 0.8, NA)))) %>%
  # transform(mean.p.true = ifelse(det.MR == 'L', 0.3, 
  #                                ifelse(det.MR == 'M', 0.5, 
  #                                       ifelse(det.MR == 'H', 0.8, NA)))) %>%
  # transform(fec.cv = (fec.obs-fec.true)^2,
  #           phiad.cv = (phiad.obs-phiad.true)^2,
  #           phi1.cv = (phi1.obs-phi1.true)^2,
  #           mean.p.cv = (mean.p-mean.p.true)^2,
  #           p.surv.cv = (p.surv-p.surv.true)^2) %>%
  dplyr::select(lambda.scenario, scenario, simscenarios, 
                param, cv_mean,
                det.MR, det.abund, det.prod) %>%
  reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'simscenarios', 'param', 'det.MR', 'det.abund', 'det.prod')) %>%
  group_by(lambda.scenario, det.MR, det.abund, det.prod, param) %>%
  dplyr::summarize(mean.cv = mean(value), .groups = 'keep') %>% # TODO - added sqrt, check
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

###average over two layers of detection (det.MR and det.prod)
cv.few <- cv.vals %>%
  group_by(variable, det.abund, lambda.scenario, dataset) %>%
  dplyr::summarize(cv = mean(mean.cv), .groups = 'keep') %>%
  #transform(det.abund = factor(det.abund, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High'))) %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))
# transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
#                            labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))

########data-generating values
cv.vals.sc <- all.stats %>%
  inner_join(data_scenarios %>% mutate(simscenarios = row_number()), 
             by = "simscenarios") %>% 
  filter(param %nin% c("meanp", "psurv")) %>% 
  # transform(fec.cv = (fec.obs-fec.true)^2,
  #           phiad.cv = (phiad.obs-phiad.true)^2,
  #           phi1.cv = (phi1.obs-phi1.true)^2) %>%
  dplyr::select(lambda.scenario, scenario, simscenarios, 
                param, cv_mean, 
                det.MR, det.abund, det.prod) %>%
  #reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'det.MR', 'det.abund', 'det.prod')) %>%
  rename(variable = param) %>% 
  #group_by(lambda.scenario, scenario, simscenarios, det.MR, det.abund, det.prod, variable) %>%
  #dplyr::summarize(mean.cv = mean(value), .groups = 'keep') %>% # TODO changed here, check
  #dplyr::summarize(cv = sd(value)/value, .groups = 'keep') %>% # CV version
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
  # merge(true_vals, by = c('lambda.scenario', 'scenario')) %>%
  # transform(variable = factor(variable, levels = c('phiad.cv', 'phi1.cv', 'fec.cv'),
  #                             labels = c('Adult survival', 'Juv survival', 'Fecundity'))) %>% 
  # transform(lambda.scenario = factor(lambda.scenario, 
  #                                    levels = c("Decreasing", "Stable", "Increasing"))) %>%
  # transform(fec_cat = ifelse(fec.true < fec_lims[1], 'L', ifelse(fec.true > fec_lims[2], 'H', 'M')),
  #           phiad_cat = ifelse(phiad.true < phiad_lims[1], 'L', ifelse(phiad.true > phiad_lims[2], 'H', 'M')),
  #           phi1_cat = ifelse(phi1.true < phi1_lims[1], 'L', ifelse(phi1.true > phi1_lims[2], 'H', 'M'))) %>%
  # transform(fec_cat = factor(fec_cat, levels = c('L', 'M', 'H'),  labels = c('Low', 'Medium', 'High')),
  #           phiad_cat = factor(phiad_cat, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High')),
  #           phi1_cat = factor(phi1_cat, levels = c('L', 'M', 'H'),
  #                             labels = c('True first-year survival: Low', 
  #                                        'True first-year survival: Medium', 
  #                                        'True first-year survival: High'))) %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(life_hist = factor(life_hist, levels = c("slow", "mod", "fast"), 
                               labels = c("Slow", "Moderate", "Fast"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only')))
# transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
#                            labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))


#facet by both fec and juv true vals
plot.vals.cv <- cv.dem %>%
  group_by(variable, life_hist, lambda.scenario, dataset) %>% # TODO changed here
  dplyr::summarize(value = mean(value), .groups = 'keep')

######################################################
##################### Lambda #########################
######################################################

# TODO - repeat changes above down here
# TODO - actually this one needs more attention, and also to have computed geom means first

row.low <- read_csv(here("results", "row_low.csv"))
row.med <- read_csv(here("results", "row_med.csv"))
row.high <- read_csv(here("results", "row_high.csv"))

scenario_ID <- read_csv(here("data", "scenario_ID.csv")) %>% 
  select(-c(n.viable.combinations, priority))

# Reformat for plotting
toplot1 <- row.low %>%
  select(contains("geomean"), scenario, sims, simscenarios, Quantile) %>%
  #group_by(model, detection)
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(scenario, sims, simscenarios, Year, Quantile) %>% # checked through here - TODO
  # summarise(low = quantile(value, 0.025),
  #           med = quantile(value, 0.5),
  #           high = quantile(value, 0.975)) %>%
  ungroup() %>%
  left_join(scenario_ID, by = "simscenarios") %>% 
  mutate(scenario = as.factor(scenario),
         sims = as.factor(sims), 
         simscenarios = as.factor(simscenarios))

# use reshape cast to spread the quantiles out into lower middle upper for line plot

toplot2 <- row.med %>%
  select(contains("geomean"), scenario, sims, simscenarios, Quantile) %>%
  #group_by(model, detection)
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(scenario, sims, simscenarios, Year, Quantile) %>% # checked through here - TODO
  # summarise(low = quantile(value, 0.025),
  #           med = quantile(value, 0.5),
  #           high = quantile(value, 0.975)) %>%
  ungroup() %>%
  left_join(scenario_ID, by = "simscenarios") %>% 
  mutate(scenario = as.factor(scenario),
         sims = as.factor(sims), 
         simscenarios = as.factor(simscenarios))

toplot3 <- row.high %>%
  select(contains("geomean"), scenario, sims, simscenarios, Quantile) %>%
  #group_by(model, detection)
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(scenario, sims, simscenarios, Year, Quantile) %>% # checked through here - TODO
  # summarise(low = quantile(value, 0.025),
  #           med = quantile(value, 0.5),
  #           high = quantile(value, 0.975)) %>%
  ungroup() %>%
  left_join(scenario_ID, by = "simscenarios") %>% 
  mutate(scenario = as.factor(scenario),
         sims = as.factor(sims), 
         simscenarios = as.factor(simscenarios))

toplot <- bind_rows(toplot1, toplot2, toplot3) %>% 
  select(-simscenarios) %>% 
  transform(dataset = ifelse(is.na(det.MR)&!is.na(det.prod), 'Abundance & Productivity', 
                             ifelse(!is.na(det.MR)&is.na(det.prod), 'Abundance & Survival',
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM')))) %>%
  group_by(Quantile, Year, det.abund, lambda, dataset, det.MR, det.prod) %>% 
  # took mean over demographic scenario (n = 25) and replicate (n = 50)
  # and mark recapture detection and fecundity
  # AEB - is it ok to take mean of quantiles? review here ######
dplyr::summarize(value = mean(value), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(Quantile = str_remove(Quantile, "\\%")) %>% 
  mutate(Quantile = paste("X", Quantile, sep = "")) %>% 
  reshape2::dcast(dataset + Year +  det.MR + det.prod + det.abund   + lambda ~ Quantile, value.var = "value") %>% 
  mutate(Year = Year + 1) %>% 
  filter(Year %in% c(15)) %>% 
  mutate(Year = factor(Year)) %>% 
  mutate(det.abund = factor(det.abund, levels = c("L", "M", "H"), labels = c("Low", "Medium", "High"))) %>% 
  mutate(det.prod = factor(det.prod, levels = c("L", "M", "H"), labels = c("Low", "Medium", "High"))) %>% 
  mutate(det.MR = factor(det.MR, levels = c("L", "M", "H"), labels = c("Low", "Medium", "High"))) %>% 
  transform(lambda = factor(lambda, levels = c("L", "M", "H"), 
                            labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Survival', 'Abundance & Productivity', 'Abundance Only'))) %>% 
# transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
#                            labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))
  mutate(intercept = case_when(
    lambda == "Decreasing" ~ 0.95, 
    lambda == "Stable" ~ 1,
    lambda == "Increasing" ~ 1.05))

######################################################
################### Figures ##########################
######################################################

# TODO might need some adjustments


#### Figure 3: RMSE and bias ecological paramters x count survey detection ####

dataset.labs <- c("Full IPM", "Abundance & Prod.", "Abundance & Surv.", "Abundance Only")
names(dataset.labs) <- c("Full IPM", "Abundance & Productivity", "Abundance & Survival", "Abundance Only")
lambda.labs <- c("Decrease", "Stable", "Increase")
names(lambda.labs) <- c("Decreasing", "Stable", "Increasing")

## bias dot plot
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


## RMSE heat map
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

## CV heat map
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

## combine
plot_grid(a1, a2, a3, nrow = 3, labels = "AUTO",
          align = "hv", label_size = 12)
ggsave(width = 6.5, height = 18, here("figures", "fig3_NEW.png"))


#### Figure 4: RMSE and bias ecological parameters x true fecundity ####

# TODO - revise from here

phi1_cat_lab <- c("True first-year\nφ: Low", 
                  "True first-year\nφ: Med",
                  "True first-year\nφ: High")

names(phi1_cat_lab) <- c("True first-year survival: Low", 
                         "True first-year survival: Medium", 
                         "True first-year survival: High")

## bias dot plot
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

## RMSE heat map
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

## CV heat map
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

## combine
plot_grid(b1, b2, b3, nrow = 3, labels = "AUTO", align = "hv", label_size = 12)
ggsave(width = 6.5, height = 18, here("figures", "fig4_NEW.png"))


#### Figure 5: Lambda trends ####

# TODO edit

c1 <- ggplot(toplot) +
  geom_point(aes(x = Year, y = X50, col = det.abund, group = det.abund, shape = det.abund), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(x = Year, ymin = X2.5, ymax = X97.5, col = det.abund, group = det.abund,
                     shape = det.abund), position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = intercept), linetype = 'dotted') +
  geom_hline(aes(yintercept = 1.0), linetype = 'solid') +
  xlab('Final year (t=15)') + #renamed the axis
  ylab(expression(lambda)) + 
  facet_grid(dataset ~ lambda, scales = 'free', labeller = label_wrap_gen()) +
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

ggsave(width = 6, height = 8, here("figures", "fig5.png"))

######################################################
################### Appendix #########################
######################################################

#### Figure 6: RMSE and bias observation parameters x count survey detection ####

## bias dot plot
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


## RMSE heat map
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

## CV heat map
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

## combine
plot_grid(d1, d2, d3, nrow = 3, labels = "AUTO", align = "hv", label_size = 12)
ggsave(width = 6.5, height = 18, here("figures", "fig6_NEW.png"))

# TODO - the widths are not working on the figures
# need some aesthetic adjustments