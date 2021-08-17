# AEB NOTE 
# this needs to be modified

# 9 panel plot
# rows - low, med, high lambda
# columns - low, med, high detection

#or 3? column is model type?

library(tidybayes)
library(ggplot2)
library(wesanderson)
library(here)
library(tidyverse)
library(reshape2)

pal <- rev(wes_palette("Zissou1", 3, type = "continuous"))

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
  transform(dataset = ifelse(is.na(det.MR)&!is.na(det.prod), 'Counts+Prod', 
                             ifelse(!is.na(det.MR)&is.na(det.prod), 'Counts+MR',
                                    ifelse(is.na(det.MR)&is.na(det.prod), 'Counts', 'Full')))) %>% 
  group_by(Quantile, Year, det.abund, lambda, dataset, det.MR, det.prod) %>% 
  # took mean over demographic scenario (n = 25) and replicate (n = 50)
  # and mark recapture detection and fecundity
  # AEB - is it ok to take mean of quantiles? review here ######
  summarize(value = mean(value), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(Quantile = str_remove(Quantile, "\\%")) %>% 
  mutate(Quantile = paste("X", Quantile, sep = "")) %>% 
  reshape2::dcast(dataset + Year +  det.MR + det.prod + det.abund   + lambda ~ Quantile, value.var = "value") %>% 
  filter(Year %in% c(1:5, 10)) %>% 
  mutate(Year = factor(Year)) %>% 
  mutate(det.abund = factor(det.abund, levels = c("L", "M", "H"))) %>% 
  mutate(det.prod = factor(det.prod, levels = c("L", "M", "H"))) %>% 
  mutate(det.MR = factor(det.MR, levels = c("L", "M", "H"))) %>% 
  transform(lambda = factor(lambda, levels = c("L", "M", "H"), 
                            labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(dataset = factor(dataset, levels = c('Full', 'Counts+Prod', 'Counts+MR', 'Counts'),
                             labels = c('Full', 'Counts+Prod', 'Counts+MR', 'Counts')))


pdf(here("figures", "lambdaTrends.pdf"), width = 12, height = 8)
rainbow2 <- c("violetred4", "dodgerblue3", 'deepskyblue1', "#4aaaa5", "#a3d39c", "#f6b61c", "chocolate2", "red3")
ggplot(toplot) +
  geom_point(aes(x = Year, y = X50, col = det.abund, group = det.abund), position = position_dodge(width = 0.5)) +
  #geom_linerange(aes(ymin = X2.5, ymax = X97.5, x = Year), position = position_dodge(width = 0.5)) +
  #geom_hline(aes(yintercept = 1.0), linetype = 'dotted') +
  xlab('Year') + ylab('Lambda') + 
  ylim(0.87, 1.14) +
  facet_grid(dataset ~ lambda, scales = 'free') +
  theme_bw() +
  theme(legend.position = 'top',
        plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1)) +
  scale_color_manual(values = rainbow2[-c(1,4)], name = 'Abundance detection level') + 
  geom_rect(aes(ymin=-Inf,
                ymax=Inf,
                xmin=which(levels(Year) == "10")-0.65,
                xmax=which(levels(Year) == "10")+0.65),
            fill="grey85", alpha=0.25, col = NA) +
  geom_hline(aes(yintercept = 1.0), linetype = 'dotted') +
  geom_linerange(aes(ymin = X2.5, ymax = X97.5, x = Year, col = det.abund, group = det.abund), position = position_dodge(width = 0.5)) +
  geom_point(aes(x = Year, y = X50, col = det.abund, group = det.abund), position = position_dodge(width = 0.5))
dev.off()
  
######

# ggplot(toplot, aes(x = Year, y = X50, col = factor(dataset), group = factor(dataset))) +
#   geom_point(position = position_dodge(width = 0.75)) + geom_linerange(aes(ymin = X2.5, ymax = X97.5, x = Year), position = position_dodge(width = 0.75)) +
#   geom_hline(aes(yintercept = 1.0), linetype = 'dotted') +
#   xlab('Year') + ylab('Lambda') +
#   facet_grid(det.abund ~ lambda, scales = 'free', labeller = "label_both") +
#   theme_bw() +
#   theme(legend.position = 'top',
#         plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1)) #+
#   #scale_color_manual(values = rainbow2[-c(1,4)], name = '')
# 
# ggplot(toplot, aes(x = Year, y = X50, col = factor(dataset), group = factor(dataset))) +
#   geom_point(position = position_dodge(width = 0.75)) + geom_linerange(aes(ymin = X2.5, ymax = X97.5, x = Year), position = position_dodge(width = 0.75)) +
#   geom_hline(aes(yintercept = 1.0), linetype = 'dotted') +
#   xlab('Year') + ylab('Lambda') +
#   facet_grid(det.MR ~ lambda, scales = 'free', labeller = "label_both") +
#   theme_bw() +
#   theme(legend.position = 'top',
#         plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1))
# 
# ggplot(toplot, aes(x = Year, y = X50, col = factor(dataset), group = factor(dataset))) +
#   geom_point(position = position_dodge(width = 0.75)) + geom_linerange(aes(ymin = X2.5, ymax = X97.5, x = Year), position = position_dodge(width = 0.75)) +
#   geom_hline(aes(yintercept = 1.0), linetype = 'dotted') +
#   xlab('Year') + ylab('Lambda') +
#   facet_grid(det.prod ~ lambda, scales = 'free', labeller = "label_both") +
#   theme_bw() +
#   theme(legend.position = 'top',
#         plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1))
# 
# ggplot(toplot %>%  filter(lambda == "L"), aes(x = Year, y = X50, col = factor(det.abund), group = factor(det.abund))) +
#   geom_point(position = position_dodge(width = 0.75)) + geom_linerange(aes(ymin = X2.5, ymax = X97.5, x = Year), position = position_dodge(width = 0.75)) +
#   geom_hline(aes(yintercept = 1.0), linetype = 'dotted') +
#   xlab('Year') + ylab('Lambda') +
#   facet_grid(det.prod ~ det.MR, scales = 'free', labeller = "label_both") +
#   theme_bw() +
#   theme(legend.position = 'top',
#         plot.subtitle = element_text(size = 10, hjust = 0.5, vjust = 1))
# 
# 
# library(patchwork)
# library(cowplot)
