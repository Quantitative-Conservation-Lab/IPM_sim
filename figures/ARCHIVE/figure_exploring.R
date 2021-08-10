library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(coda)
library(captioner)
library(knitr)
library(reshape2)
library(here)

rm(list=ls())

#true vals
truth <- readRDS(here::here('data', 'ARCHIVE', 'scenarios.Rdata')) 
colnames(truth) <- gsub(' ', '', colnames(truth))
true_vals <- truth %>%
  dplyr::select(ScenarioNumber, MRdetection, Abunddetection, AdultSurv, JuvSurv, Dailynestsurvival, Fec) %>%
  dplyr::rename(scenario = ScenarioNumber) %>%
  reshape2::melt(id.vars = c('scenario'), value.name = 'Truth') %>%
  transform(variable = factor(variable,
                              levels = c('JuvSurv', 'AdultSurv', 'Fec', 'Dailynestsurvival', 'MRdetection', 'Abunddetection'),
                              labels = c('Juvenile survival', 'Adult survival', 'Fecundity', 'Daily nest survival', 
                                         'MR detection', 'Survey detection'))) %>%
  transform(Data = ifelse(scenario < 6, 'All', 
                          ifelse(scenario < 11 & scenario > 5, 'Counts + Productivity', 'Counts + MR'))) %>%
  transform(Truth = ifelse(Data == 'Counts + Productivity' & variable == 'MR detection', NA, Truth))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~ #

### assessing convergence

# test <-readRDS(file = here::here('data', paste0('out_3-1.Rdata')))[3500:4000,]


sc <- c(1:5, 11:15)
# sc <- c(1:5)
nsim <- 25
pars <- c('fec', 'lambdaf', 'mean.p', 'mean.phi[1]', 'mean.phi[2]', 'p.surv', 'phi.nest')

out_meds1 <- data.frame() #storage
for (s in sc) {
  for (i in 1:nsim) { 
    chains <- as.mcmc.list(list(
      as.mcmc(readRDS(file = here::here('data', 'ARCHIVE', paste0('out', s, '_', i, '-', 1, '.Rdata')))[3500:4000,]), 
      as.mcmc(readRDS(file = here::here('data', 'ARCHIVE', paste0('out', s, '_', i, '-', 2, '.Rdata')))[3500:4000,]), 
      as.mcmc(readRDS(file = here::here('data', 'ARCHIVE', paste0('out', s, '_', i, '-', 3, '.Rdata')))[3500:4000,])))
    
    max_rhat <- max(gelman.diag(chains[,pars[-6]], multivariate=F, autoburnin=F)$psrf[,1], na.rm = T)
    
    ## posterior medians
    out_temp <- data.frame(t(c(summary(chains)$q[pars,"50%"], max_rhat, i, s)))
    colnames(out_temp) <- c(names(summary(chains)$q[pars,"50%"]), 'max_rhat', 'sim', 'scenario')
    
    out_meds1 <- rbind(out_meds1, out_temp)
    
  } #sims
} #scenarios

#scenarios without p.surv
sc <- c(6:10)
# sc <- c(6:9)

pars <- c('fec', 'lambdaf', 'mean.phi[1]', 'mean.phi[2]', 'p.surv', 'phi.nest')

out_meds2 <- data.frame() #storage
for (s in sc) {
  for (i in 1:nsim) { 
    chains <- as.mcmc.list(list(
      as.mcmc(readRDS(file = here::here('data', 'ARCHIVE', paste0('out', s, '_', i, '-', 1, '.Rdata')))[3500:4000,]), 
      as.mcmc(readRDS(file = here::here('data', 'ARCHIVE', paste0('out', s, '_', i, '-', 2, '.Rdata')))[3500:4000,]), 
      as.mcmc(readRDS(file = here::here('data', 'ARCHIVE', paste0('out', s, '_', i, '-', 3, '.Rdata')))[3500:4000,])))
    
    max_rhat <- max(gelman.diag(chains[,pars], multivariate=F, autoburnin=F)$psrf[,1], na.rm = T)
    
    ## posterior medians
    out_temp <- data.frame(t(c(summary(chains)$q[pars,"50%"], max_rhat, i, s)))
    colnames(out_temp) <- c(names(summary(chains)$q[pars,"50%"]), 'max_rhat', 'sim', 'scenario')
    
    out_meds2 <- rbind(out_meds2, out_temp)
    
  } #sims
} #scenarios

out_meds <- bind_rows(out_meds1, out_meds2)

post_meds <- out_meds %>% 
  filter(max_rhat < 1.8) %>%
  dplyr::select(-c(max_rhat, lambdaf)) %>%
  reshape2::melt(id.vars = c('scenario', 'sim')) %>%
  filter(variable %in% c('mean.phi[1]', 'mean.phi[2]', 'fec', 'phi.nest', 'p.surv', 'mean.p')) %>%
  transform(variable = factor(variable,
                              levels = c('mean.phi[1]', 'mean.phi[2]', 'fec', 'phi.nest', 'mean.p', 'p.surv'),
                              labels = c('Juvenile survival', 'Adult survival', 'Fecundity', 'Daily nest survival', 
                                         'MR detection', 'Survey detection'))) %>%
  transform(Data = ifelse(scenario < 6, 'All', 
                          ifelse(scenario < 11 & scenario > 5, 'Counts + Productivity', 'Counts + MR')))

meds_plot <- ggplot(post_meds, aes(y = value, x = scenario, group = scenario)) +
  geom_boxplot(aes(color = Data)) +
  geom_point(size = 0.7, col = 'darkgrey') + 
  geom_point(data = true_vals, aes(y = Truth, x = scenario, group = scenario), col = 'darkred') +
  xlab('Scenario') + ylab('Posterior medians') +
  facet_wrap(~variable, scales = 'free_y', nrow = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'top') +
  scale_x_continuous(breaks = 1:15)
 # scale_color_manual(values = cols[c(2,3,5)])

## try ridgeline
library(ggridges)
library(ggplot2)
library(viridis)

ggplot(post_meds) +
  geom_density_ridges(mapping = aes(x = value, y = as.factor(scenario), fill = as.factor(scenario))) +
  theme_minimal() +
  facet_wrap(~variable, scales = 'free_y', nrow = 3) +
  xlab("Scenario") + ylab("") +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  geom_point(data = true_vals, aes(y = as.factor(scenario), x = Truth, group = scenario), 
             col = 'black')

## also figure on over, under-estimaitng posteriors

#relative bias
bias_vals <- post_meds %>% 
  merge(true_vals %>% dplyr::select(scenario, variable, Truth), by = c('scenario', 'variable')) %>%
  transform(bias = value-Truth, rel_bias = (value-Truth)/Truth) 

bias_plot <- ggplot(bias_vals, aes(y = rel_bias, x = scenario, group = scenario)) +
  geom_boxplot(aes(color = Data)) +
  geom_point(size = 0.7, col = 'darkgrey') + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed', col = 'darkred') +
  xlab('Scenario') + ylab('Relative bias') +
  facet_wrap(~variable, scales = 'free_y', nrow = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'top') +
  scale_x_continuous(breaks = 1:15) +
  scale_color_manual(values = cols[c(2,3,5)])



#rmse
rmse_vals <- post_meds %>% 
  merge(true_vals %>% dplyr::select(scenario, variable, Truth), by = c('scenario', 'variable')) %>%
  transform(rmse.temp = (value-Truth)^2) %>%
  group_by(Data, scenario, variable) %>%
  dplyr::summarize(rmse = sqrt(mean(rmse.temp)))

# rmse_plot <- ggplot(rmse_vals, aes(y = rmse, x = scenario, group = scenario)) +
#   geom_boxplot(aes(color = Data)) +
#   geom_point(size = 0.7, col = 'darkgrey') +
#   # geom_hline(aes(yintercept = 0), linetype = 'dashed', col = 'darkred') +
#   xlab('Scenario') + ylab('RMSE') +
#   facet_wrap(~variable, scales = 'free_y', nrow = 3) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         legend.position = 'top') +
#   scale_x_continuous(breaks = 1:15) +
#   scale_color_manual(values = cols[c(2,3,5)])
# 
# figs("rmse_plot", caption = "Figure of RMSE per scenario.")

posterior_summary <- post_meds %>%
  group_by(scenario, variable) %>%
  dplyr::summarize(post_med = round(median(value), 2)) %>%
  dcast(scenario ~ variable, value.var = 'post_med')

bias_summary <- bias_vals %>%
  group_by(scenario, variable) %>%
  dplyr::summarize(mean_bias = round(mean(rel_bias), 2)) %>%
  dcast(scenario ~ variable, value.var = 'mean_bias')
write.csv(bias_summary, file = here::here('results', 'bias_summary.csv'))

rmse_summary <- rmse_vals %>%
  group_by(scenario, variable) %>%
  dplyr::summarize(mean_rmse = round(mean(rmse), 2)) %>%
  dcast(scenario ~ variable, value.var = 'mean_rmse')
write.csv(rmse_summary, file = here::here('results', 'rmse_summary.csv'))

summary_tab <- bias_summary %>%
  merge(rmse_summary, by = 'scenario', suffixes = c('.bias', '.rmse')) %>%
  merge(true_vals %>% dplyr::select(scenario, Data) %>% distinct(), by = 'scenario') %>%
  dplyr::select(Data, scenario, JuvSurv.bias, JuvSurv.rmse, AdultSurv.bias, AdultSurv.rmse, 
                Fec.bias, Fec.rmse, Dailynestsurvival.bias, Dailynestsurvival.rmse, 
                Abunddetection.bias, Abunddetection.rmse, MRdetection.bias, MRdetection.rmse)

