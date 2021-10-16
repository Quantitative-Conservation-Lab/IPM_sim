## IPM Sim Figures
## Adapted from PaperFigures.RMD
## A DuVall
## 15 Oct 2021

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

rainbow2 <- c("violetred4", "dodgerblue3", 'deepskyblue1', "#4aaaa5", "#a3d39c", "#f6b61c", "chocolate2", "red3")

scenarios <- read.csv(file = here::here('data', 'scenario_ID.csv'), header = T, stringsAsFactors = F)

low.params <- readRDS(file = here::here('data', 'low.lam.params.rds')) %>%
  transform(scenario = 1:25)

med.params <- readRDS(file = here::here('data', 'med.lam.params.rds')) %>%
  transform(scenario = 1:25)

high.params <- readRDS(file = here::here('data', 'high.lam.params.rds')) %>%
  transform(scenario = 1:25)

all.meds <- read.csv(file = here::here('figures', 'Processed csvs', 'all.meds.csv'), 
                     header = T, stringsAsFactors = F)

######################################################
##################### Bias ###########################
######################################################

rel.bias <- all.meds %>%
  transform(p.surv.true = ifelse(det.abund == 'L', 0.3, 
                                 ifelse(det.abund == 'M', 0.5, 
                                        ifelse(det.abund == 'H', 0.8, NA)))) %>%
  transform(mean.p.true = ifelse(det.MR == 'L', 0.3, 
                                 ifelse(det.MR == 'M', 0.5, 
                                        ifelse(det.MR == 'H', 0.8, NA)))) %>%
  #calculate relative bias, mean across scenarios
  transform(phi1.bias = (phi1.obs-phi1.true)/phi1.true,
            phiad.bias = (phiad.obs-phiad.true)/phiad.true,
            fec.bias = (fec.obs-fec.true)/fec.true,
            p.surv.bias = (p.surv-p.surv.true)/p.surv.true,
            mean.p.bias = (mean.p-mean.p.true)/mean.p.true) %>%
  dplyr::select('lambda.scenario', 'scenario',
                'phi1.bias', 'phiad.bias', 'fec.bias', 'det.MR', 'det.abund', 'det.prod', 'mean.p.bias', 'p.surv.bias') %>%
  reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'det.MR', 'det.abund', 'det.prod')) %>%
  group_by(lambda.scenario, det.MR, det.abund, det.prod, variable) %>%
  dplyr::summarize(bias = mean(value), .groups = 'keep') %>%
  transform(variable = factor(variable, levels = c('phiad.bias', 'phi1.bias', 'fec.bias', 'p.surv.bias', 'mean.p.bias'),
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
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))

all.meds.sc <- read.csv(file = here::here('figures', 'Processed csvs', 'all.meds.sc.csv'), 
                        header = T, stringsAsFactors = F)

## just bias; same as above, keep 'scenario'
rel.bias.sc <- all.meds.sc %>%
  transform(phi1.bias = (phi1.obs-phi1.true)/phi1.true,
            phiad.bias = (phiad.obs-phiad.true)/phiad.true,
            fec.bias = (fec.obs-fec.true)/fec.true) %>%
  dplyr::select('lambda.scenario', 'scenario', 
                'phi1.bias', 'phiad.bias', 'fec.bias', 'det.MR', 'det.abund', 'det.prod') %>%
  reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'det.MR', 'det.abund', 'det.prod')) %>%
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
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM'))))

#load saved true values
true_vals <- read.csv(file = here::here('data', 'true.vals.csv'), header = T, stringsAsFactors = F)

#for categories: 
fec_lims <- quantile(true_vals$fec.true, probs = c(0.33, 0.7), names = F)
phiad_lims <- quantile(true_vals$phiad.true, probs = c(0.33, 0.7), names = F)
phi1_lims <- quantile(true_vals$phi1.true, probs = c(0.33, 0.7), names = F)

#merge back onto bias df so can have something to name/view the 'scenarios' 
rel.bias.dem <- rel.bias.sc %>%
  merge(true_vals, by = c('lambda.scenario', 'scenario')) %>%
  transform(variable = factor(variable, levels = c('phiad.bias', 'phi1.bias', 'fec.bias'),
                              labels = c('Adult survival', 'First-year survival', 'Fecundity'))) %>% 
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(fec_cat = ifelse(fec.true < fec_lims[1], 'L', ifelse(fec.true > fec_lims[2], 'H', 'M')),
            phiad_cat = ifelse(phiad.true < phiad_lims[1], 'L', ifelse(phiad.true > phiad_lims[2], 'H', 'M')),
            phi1_cat = ifelse(phi1.true < phi1_lims[1], 'L', ifelse(phi1.true > phi1_lims[2], 'H', 'M'))) %>%
  transform(fec_cat = factor(fec_cat, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High')),
            phiad_cat = factor(phiad_cat, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High')),
            phi1_cat = factor(phi1_cat, levels = c('L', 'M', 'H'),
                              labels = c('True first-year survival: Low', 
                                         'True first-year survival: Medium', 
                                         'True first-year survival: High'))) %>%
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))

#facet by both fec and juv true vals
plot.vals <- rel.bias.dem %>%
  group_by(variable, fec_cat, phi1_cat, dataset) %>%
  dplyr::summarize(value = mean(value), .groups = 'keep') 



######################################################
##################### RMSE ###########################
######################################################
rmse.vals <- all.meds %>%
  transform(p.surv.true = ifelse(det.abund == 'L', 0.3, 
                                 ifelse(det.abund == 'M', 0.5, 
                                        ifelse(det.abund == 'H', 0.8, NA)))) %>%
  transform(mean.p.true = ifelse(det.MR == 'L', 0.3, 
                                 ifelse(det.MR == 'M', 0.5, 
                                        ifelse(det.MR == 'H', 0.8, NA)))) %>%
  transform(fec.rmse = (fec.obs-fec.true)^2,
            phiad.rmse = (phiad.obs-phiad.true)^2,
            phi1.rmse = (phi1.obs-phi1.true)^2,
            mean.p.rmse = (mean.p-mean.p.true)^2,
            p.surv.rmse = (p.surv-p.surv.true)^2) %>%
  dplyr::select(lambda.scenario, scenario, fec.rmse, phiad.rmse, phi1.rmse, p.surv.rmse, mean.p.rmse,
                det.MR, det.abund, det.prod) %>%
  reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'det.MR', 'det.abund', 'det.prod')) %>%
  group_by(lambda.scenario, det.MR, det.abund, det.prod, variable) %>%
  dplyr::summarize(mean.rmse = mean(value), .groups = 'keep') %>%
  transform(variable = factor(variable, levels = c('phiad.rmse', 'phi1.rmse', 'fec.rmse', 
                                                   'p.surv.rmse', 'mean.p.rmse'),
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
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))

########data-generating values
rmse.vals.sc <- all.meds %>%
  transform(fec.rmse = (fec.obs-fec.true)^2,
            phiad.rmse = (phiad.obs-phiad.true)^2,
            phi1.rmse = (phi1.obs-phi1.true)^2) %>%
  dplyr::select(lambda.scenario, scenario, fec.rmse, phiad.rmse, phi1.rmse, 
                det.MR, det.abund, det.prod) %>%
  reshape2::melt(id.vars = c('lambda.scenario', 'scenario', 'det.MR', 'det.abund', 'det.prod')) %>%
  group_by(lambda.scenario, scenario, det.MR, det.abund, det.prod, variable) %>%
  dplyr::summarize(mean.rmse = mean(value), .groups = 'keep') %>%
  transform(variable = factor(variable, levels = c('phiad.rmse', 'phi1.rmse', 'fec.rmse'),
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
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Abundance Only', 'Full IPM'))))

rmse.dem <- rmse.vals.sc %>%
  merge(true_vals, by = c('lambda.scenario', 'scenario')) %>%
  # transform(variable = factor(variable, levels = c('phiad.rmse', 'phi1.rmse', 'fec.rmse'),
  #                             labels = c('Adult survival', 'Juv survival', 'Fecundity'))) %>% 
  transform(lambda.scenario = factor(lambda.scenario, 
                                     levels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(fec_cat = ifelse(fec.true < fec_lims[1], 'L', ifelse(fec.true > fec_lims[2], 'H', 'M')),
            phiad_cat = ifelse(phiad.true < phiad_lims[1], 'L', ifelse(phiad.true > phiad_lims[2], 'H', 'M')),
            phi1_cat = ifelse(phi1.true < phi1_lims[1], 'L', ifelse(phi1.true > phi1_lims[2], 'H', 'M'))) %>%
  transform(fec_cat = factor(fec_cat, levels = c('L', 'M', 'H'),  labels = c('Low', 'Medium', 'High')),
            phiad_cat = factor(phiad_cat, levels = c('L', 'M', 'H'), labels = c('Low', 'Medium', 'High')),
            phi1_cat = factor(phi1_cat, levels = c('L', 'M', 'H'),
                              labels = c('True first-year survival: Low', 
                                         'True first-year survival: Medium', 
                                         'True first-year survival: High'))) %>%
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only')))


#facet by both fec and juv true vals
plot.vals.rmse <- rmse.dem %>%
  group_by(variable, fec_cat, phi1_cat, dataset) %>%
  dplyr::summarize(value = mean(mean.rmse), .groups = 'keep')

######################################################
##################### Lambda #########################
######################################################

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
  transform(dataset = factor(dataset, levels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'),
                             labels = c('Full IPM', 'Abundance & Productivity', 'Abundance & Survival', 'Abundance Only'))) %>% 
  mutate(intercept = case_when(
    lambda == "Decreasing" ~ 0.95, 
    lambda == "Stable" ~ 1,
    lambda == "Increasing" ~ 1.05))

######################################################
################### Figures ##########################
######################################################


#### Figure 3: RMSE and bias ecological paramters x count survey detection ####

## bias dot plot
a1 <- ggplot(rel.bias.few  %>% filter(variable %nin% obs.pars), 
             aes(x = det.abund, y = bias, col = factor(variable), group = factor(variable),
                 shape = factor(variable))) +
  geom_point() + geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  xlab('Count survey detection') + ylab('Relative bias') +
  facet_grid(dataset~lambda.scenario, scales = 'free_x', labeller = label_wrap_gen()) +
  ylim(c(-1.75, 1.75)) +
  theme_bw() +
  theme(legend.position = 'top',
        plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line")) +
  scale_color_manual(values = rainbow2[-c(1,4)], name = '') +
  scale_shape_manual(values = c(15, 16, 17), name = '')
a1

## RMSE heat map
a2 <- ggplot(rmse.few %>% filter(variable %nin% obs.pars), aes(x = factor(det.abund), y = variable, fill = rmse)) +
  geom_tile(color = 'grey50') +
  xlab('Count survey detection') + ylab('') +
  facet_grid(dataset ~ lambda.scenario, drop = T, scales = 'free_x', labeller = label_wrap_gen()) +
  scale_fill_gradient2(name = "RMSE",
                       mid = "white", high = rainbow2[2], midpoint = 0) +
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
a2

## combine
plot_grid(a1, a2, ncol = 2, labels = "AUTO", align = "h", label_size = 16)
ggsave(width = 15, height = 8, here("figures", "fig3.png"))

#### Figure 4: RMSE and bias ecological parameters x true fecundity ####

## bias dot plot
b1 <- ggplot(plot.vals, aes(x = fec_cat, y = value, col = factor(variable), group = factor(variable),
                            shape = factor(variable))) +
  geom_point() + geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  xlab('True fecundity') + ylab('Relative bias') +
  facet_grid(dataset~phi1_cat, scales = 'free_x', labeller = label_wrap_gen()) +
  ylim(c(-1.75, 1.75)) +
  theme_bw() +
  theme(legend.position = 'top',
        plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
        panel.border = element_rect(color = "black", fill = NA),  
        panel.spacing.x = unit(0.75, "line")) +
  scale_color_manual(values = rainbow2[-c(1,4)], name = '') +
  scale_shape_manual(values = c(15, 16, 17), name = '')
b1

## RMSE heat map
b2 <- ggplot(plot.vals.rmse, aes(x = factor(fec_cat), y = factor(variable), fill = value)) +
  geom_tile(color = 'grey50') +
  xlab('True fecundity') + ylab('') +
  facet_grid(dataset~phi1_cat, drop = T, scales = 'free_x', labeller = label_wrap_gen()) +
  scale_fill_gradient2(name = "RMSE", mid = "white", high = rainbow2[2], midpoint = 0) +
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
b2

## combine
plot_grid(b1, b2, ncol = 2, labels = "AUTO", align = "h", label_size = 16)
ggsave(width = 15, height = 8, here("figures", "fig4.png"))

#### Figure 5: Lambda trends ####

c1 <- ggplot(toplot) +
  geom_point(aes(x = Year, y = X50, col = det.abund, group = det.abund, shape = det.abund), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(x = Year, ymin = X2.5, ymax = X97.5, col = det.abund, group = det.abund,
                     shape = det.abund), position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = intercept), linetype = 'dotted') +
  geom_hline(aes(yintercept = 1.0), linetype = 'solid') +
  xlab('Years') + 
  ylab(expression(lambda)) + 
  facet_grid(dataset ~ lambda, scales = 'free', labeller = label_wrap_gen()) +
  theme_bw() +
  theme(legend.position = 'top',
        plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 1.5),
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
  ylim(c(-1.75, 1.75)) +
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
                       mid = "white", high = rainbow2[2], midpoint = 0) +
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

## combine
plot_grid(d1, d2, ncol = 2, labels = "AUTO", align = "h", label_size = 16)
ggsave(width = 15, height = 8, here("figures", "fig6.png"))
