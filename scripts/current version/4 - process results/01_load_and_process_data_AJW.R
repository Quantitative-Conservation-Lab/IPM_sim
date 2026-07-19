library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)
library(coda)
library(nlist)
library(ggplot2)
# library(beepr)

rainbow2 <- c("violetred4", "dodgerblue3", 'deepskyblue1', "#4aaaa5", "#a3d39c", "#f6b61c", "chocolate2", "red3")


# load scenarios
dem_scenarios <- readRDS(here("data", "demographic_scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename("phi1" = "S.J", 
    "phiad" = "S.A",
    "fec" = "f")

surv_scenarios <- readRDS(here("data", "data_scenarios.RDS")) %>%
  dplyr::filter(det.abund != 'L')

sims.per <- 3000 

### new option
library(dplyr)
library(coda)

# storage
results_list <- list()
missing_files <- data.frame(type=character(), d=integer(), s=integer(), i=integer())
not_converged <- data.frame(type=character(), d=integer(), s=integer(), 
                                 i=integer(), gelman=numeric(), max_param=character())

# processing function (returns the data frame or NULL)
process_model <- function(prefix, d, s, i) {
  
  file_name <- file.path("D:/IPM sim results", 
                         paste0(prefix, "-", d, "-", s, "-", i, ".RDS"))
  
  # 1. Check if file physically exists
  if (!file.exists(file_name)) {
    missing_files <<- bind_rows(missing_files, data.frame(type=prefix, d=d, s=s, i=i))
    return(NULL)
  }
  
  # 2. Try to read the file safely (handles corrupted/incomplete files)
  out_temp <- tryCatch({
    readRDS(file_name)
  }, error = function(e) {
    message(paste("Skipping corrupted file at i =", i, "s =", s, "d =", d))
    missing_files <<- bind_rows(missing_files, data.frame(type=prefix, d=d, s=s, i=i))
    return(NULL)
  })
  
  # If tryCatch caught an error and returned NULL, exit the function early
  if (is.null(out_temp)) return(NULL)
  
  # 3. Convergence checks
  diag_result <- gelman.diag(out_temp, multivariate = FALSE)[[1]]
  # Find the value and the name of the worst offender
  max_idx   <- which.max(diag_result[, 1])
  max_val   <- diag_result[max_idx, 1]
  max_param <- rownames(diag_result)[max_idx]
  
  # compile convergence info
  if (is.na(max_val) || max_val > 1.1) {
    not_converged <<- bind_rows(not_converged, 
                                data.frame(
                                  type=prefix, d=d, s=s, i=i, 
                                  gelman=max_val, 
                                  max_param=max_param
                                ))
    return(NULL)
  }
  
  # 4. Process and add indices (only runs if file is good AND converged)
  out <- out_temp %>% 
    collapse_chains() %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    filter(row_number() %% 60 == 1) %>% 
    mutate(model_type = prefix,
           dem_scenario = d, 
           surv_scenario = s, 
           sim_rep = i)
  
  return(out)
}

# process
model_types <- c("out_IPM", "out_noMR", "out_noProd", "out_abundOnly")
# model_types <- c("out_IPM")

for (i in 1:sims.per) { # sims per
  for (s in 1:nrow(surv_scenarios)) { #scenarios picked
    for (d in 1:nrow(dem_scenarios)) { # simulation scenario
      for (type in model_types) {
        
        result <- process_model(type, d, s, i)
        
        if (!is.null(result)) {
          # Unique key for the list entry
          item_name <- paste(type, d, s, i, sep = "_")
          results_list[[item_name]] <- result
        }
      }
    }
  }
}

# combine and save
results_all <- bind_rows(results_list)

#just take 100 of each type x d x s since some combo's have accumulated 100+
results_all_100 <- results_all %>%
  group_by(model_type, dem_scenario, surv_scenario) %>%
  slice_sample(n = 100, replace = F) 

saveRDS(results_all_100, file = here('results', 'processed', "results_all_final_ursus.RDS"))

# convergence_summary <- results_all %>%
#   distinct(model_type, dem_scenario, surv_scenario, sim_rep) %>%
#   group_by(model_type) %>%
#   dplyr::summarize(
#     successful_sims = n(), 
#     # Optional: Calculate percentage based on your 'sims.per' variable
#     # percent_converged = (n() / sims.per) * 100,
#     .groups = "drop")
# 
# library(dplyr)

# count successful runs per scenario 
success_counts <- results_all %>%
  distinct(model_type, dem_scenario, surv_scenario, sim_rep) %>%
  group_by(model_type, dem_scenario, surv_scenario) %>%
  summarize(converged_reps = n(), .groups = "drop") %>%
  # Rename columns to match the trackers
  rename(type = model_type, d = dem_scenario, s = surv_scenario)

# count runs that failed convergence per scenario 
fail_counts <- not_converged %>%
  group_by(type, d, s) %>%
  summarize(failed_reps = n(), .groups = "drop")

# combine to get the real denominator (total files found)
convergence_summary <- full_join(success_counts, fail_counts, by = c("type", "d", "s")) %>%
  mutate(
    # Replace NAs with 0 if a scenario had 100% success or 100% failure
    converged_reps = coalesce(converged_reps, 0L),
    failed_reps    = coalesce(failed_reps, 0L),
    
    # Total files actually found for this scenario
    total_found    = converged_reps + failed_reps,
    
    # Calculate percentage based ONLY on existing files
    percent_converged = (converged_reps / total_found) * 100
  ) %>%
  arrange(type, d, s)

needs_more <- convergence_summary %>%
  filter(converged_reps < 90)

#total "successful sims" expected
#100 replicate simulations
#9 dem scenarios
#32 surv scenarios

#IPM: 18 surv scenarios * 9 dem scenarios = 162, need 16,200 success 
#no_MR and no_nest: 6 surv scenarios * 9 dem scenarios = 54, need 5400 success
#abund_only: 2 surv scenarios * 9 dem scenarios = 18, need 1800 success
# 162+108+18 = 288 total * 100 = 28800 files total expected


ggplot(convergence_summary, aes(x = factor(s), y = percent_converged, 
                                col = factor(d), group = s)) +
  geom_boxplot(position = position_dodge(0.5)) +
  facet_wrap(~type, scales = 'free') +
  theme_bw()

ggplot(convergence_summary, aes(x = factor(s), y = percent_converged, 
                                col = factor(d), group = s)) +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~type, scales = 'free') +
  theme_bw()

ggplot(convergence_summary, aes(x = factor(s), y = converged_reps,
                                col = factor(d), group = s)) +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~type, scales = 'free') +
  theme_bw()

#needs more
ggplot(needs_more %>% filter(type == 'out_noMR'), aes(x = factor(s), y = converged_reps,
                                col = factor(d), group = s)) +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~type, scales = 'free') +
  theme_bw()

ggplot(convergence_summary, aes(x = factor(s), y = total_found,
                                col = factor(d), group = s)) +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~type, scales = 'free') +
  theme_bw()

#add back survey and demographic features
conv_scenarios <- convergence_summary %>%
  merge(surv_scenarios %>%
          mutate(s = row_number()), by = 's') %>%
  merge(dem_scenarios %>%
          mutate(d = row_number()), by = 'd') %>%
  dplyr::select(type, percent_converged, det.abund, det.MR, det.prod, life_hist, trend) %>%
  transform(trend = factor(trend, levels = c('decline', 'stable', 'increase'),
                           labels = c("Decreasing", "Stable", "Increasing"))) %>%
  transform(det.MR = factor(det.MR, levels = c('L', 'M', 'H'))) %>%
  transform(det.abund = factor(det.abund, levels = c('M', 'H'),
                               labels = c('L', 'H'))) %>%
  transform(det.prod = factor(det.prod, levels = c('L', 'M', 'H'))) %>%
  transform(type = factor(type, levels = c('out_IPM', 'out_noProd', 'out_noMR', 'out_abundOnly'),
                             labels = c('Full IPM', 'Abundance & Survival', 
                                        'Abundance & Productivity', 'Abundance Only'))) %>%
  transform(life_hist = factor(life_hist, levels = c('slow', 'mod', 'fast'),
                               labels = c('Slow', 'Moderate', 'Fast')))


#by count survey detection and lambda
conv_plot1 <- ggplot(conv_scenarios, aes(x = det.abund, y = percent_converged, 
                                col = factor(det.abund), group = det.abund)) +
  geom_boxplot(position = position_dodge(0.5)) +
  xlab('Count survey detection') + ylab('Percent converged') +
  facet_grid(type~trend) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_manual(values = rainbow2[c(2,5)])

#by life history and lambda
conv_plot2 <- ggplot(conv_scenarios, aes(x = life_hist, y = percent_converged, 
                           col = factor(life_hist), group = life_hist)) +
  geom_boxplot(position = position_dodge(0.5)) +
  xlab('Life history') + ylab('Percent converged') +
  facet_grid(type~trend) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_manual(values = rainbow2[c(2,5,6)])

plot_grid(conv_plot1, conv_plot2, nrow = 1, labels = "AUTO", align = "hv", label_size = 12)
# ggsave(width = 12, height = 10, here("figures", 'final', "fig7.png"))

##which params are failing most? 
conv_dat <- not_converged %>%
  merge(convergence_summary %>% dplyr::select(d,s,type,total_found, failed_reps), 
        by = c('d', 's', 'type')) %>%
  group_by(d, s, type, max_param, total_found) %>%
  dplyr::summarize(n_fail = n_distinct(gelman)) %>%
  filter(n_fail > 10) %>%
  transform(prop_fail = n_fail/total_found) %>%
  transform(max_param = ifelse(max_param == 'mean.phi[2]', 'Adult survival',
                               ifelse(max_param == 'mean.phi[1]', 'Juvenile survival',
                                      ifelse(max_param == 'fec', 'Fecundity',
                                             ifelse(max_param == 'p.surv', 'Count surv detection', max_param))))) %>%
  transform(max_param = factor(max_param, levels = c('Adult survival', 'Juvenile survival',
                                                     'Fecundity', 'Count surv detection'))) %>%
  transform(type = factor(type, levels = c('out_IPM', 'out_noProd', 'out_noMR', 'out_abundOnly'),
                          labels = c('Full IPM', 'Abundance & Survival', 
                                     'Abundance & Productivity', 'Abundance Only')))
  

ggplot(conv_dat, aes(x = max_param, y = prop_fail)) +
  geom_boxplot() +
  xlab('') + ylab('Proportion not converged') +
  facet_wrap(~type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# parameters that are failing to converge
ggplot(not_converged, aes(x = reorder(max_param, max_param, function(x) -length(x)))) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  facet_wrap(.~type) +
  labs(
    title = "Which parameters fail to converge most often?",
    subtitle = "Counts of highest R-hat (> 1.1) per model run",
    x = "Parameter Name",
    y = "Number of Failures"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#rhat values
ggplot(not_converged, aes(x = max_param, y = gelman, color = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
  facet_wrap(.~type) +
  theme_minimal() +
  labs(
    x = "Parameter",
    y = "Gelman-Rubin Diagnostic (R-hat)",
    color = "Model Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggplot(not_converged %>% filter(gelman < 20), aes(x = max_param, y = gelman, color = type)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   geom_jitter(width = 0.2, alpha = 0.7) +
#   geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
#   facet_wrap(.~type) +
#   theme_minimal() +
#   labs(
#     x = "Parameter",
#     y = "Gelman-Rubin Diagnostic (R-hat)",
#     color = "Model Type"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(not_converged, aes(x = max_param, y = gelman, color = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
  facet_wrap(.~d) +
  theme_minimal() +
  labs(
    x = "Parameter",
    y = "Gelman-Rubin Diagnostic (R-hat)",
    color = "Model Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(not_converged, aes(x = max_param, y = gelman, color = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
  facet_wrap(.~s) +
  theme_minimal() +
  labs(
    x = "Parameter",
    y = "Gelman-Rubin Diagnostic (R-hat)",
    color = "Model Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#looking at individual runs
#struggle params for IPM are psurv and Ntot's
testIPM.b <- readRDS(file = here('results', 'normObs',
                                 # 'ind400',
                                 'out_IPM-1-12-7.RDS'))
testIPM.w <- readRDS(file = here('results', 'nmix', 'ind300_nsam5', 'out_IPM-1-13-2.RDS')) #best
testIPM.w <- readRDS(file = here('results', 'vSJC', 'out_IPM-1-13-3.RDS')) #best


plot(testIPM.b[,'mean.phi[1]'])
plot(testIPM.b[,'p.surv'])
plot(testIPM.b[,'Ntot[6]'])
plot(testIPM.w[,'mean.phi[1]'])
plot(testIPM.w[,'p.surv'])
plot(testIPM.w[,'mean.p'])
plot(testIPM.w[,'Ntot[10]'])

#struggle params - Ntots and psurv (zero vital rates)
testNoprod.b <- readRDS(file = here('results', 'ind400', 'out_noProd-4-37-8.RDS')) #3-44-2 best
testNoprod.w <- readRDS(file = here('results', 'ind400', 'out_noProd-1-40-7.RDS'))

plot(testNoprod.b[,'mean.phi[1]'])
plot(testNoprod.b[,'fec'])
plot(testNoprod.b[,'p.surv'])
plot(testNoprod.b[,'Ntot[10]'])
plot(testNoprod.w[,'mean.phi[1]'])
plot(testNoprod.w[,'p.surv'])
plot(testNoprod.w[,'Ntot[10]'])

#struggle params - phi's only, very little issues with Ntot and psurv
#also look at 3-23-9; 1-24-4 (worst); and 3-36-7 (best, 1.1005)
testNoMR.w <- readRDS(file = here('results', 'normObs', 'out_noMR-1-23-8.RDS'))
testNoMR.w <- readRDS(file = here('results', 'ind400', 'out_noMR-3-11-4.RDS'))
testNoMR.b <- readRDS(file = here('results', 'ind400', 'out_noMR-5-10-5.RDS'))
testNoMR.b <- readRDS(file = here('results', 'ind400', 'out_noMR-2-11-7.RDS'))

plot(testNoMR.w[,'mean.phi[1]'])
plot(testNoMR.w[,'mean.phi[2]'])
plot(testNoMR.w[,'fec'])
plot(testNoMR.b[,'mean.phi[1]'])
plot(testNoMR.b[,'mean.phi[2]'])
plot(testNoMR.b[,'p.surv'])

#struggle params - all vitals, nothing else; 1-47-1 best
testabund.w  <- readRDS(file = here('results', 'ind400', 'out_abundOnly-4-46-1.RDS'))
testabund.b  <- readRDS(file = here('results', 'ind400', 'out_abundOnly-2-48-5.RDS'))

plot(testabund.b[,'fec'])
plot(testabund.b[,'mean.phi[2]'])
plot(testabund.b[,'p.surv'])
plot(testabund.b[,'Ntot[3]'])
plot(testabund.w[,'mean.phi[1]'])
plot(testabund.w[,'mean.phi[2]'])
plot(testabund.w[,'p.surv'])

#### archive - older version ####
# # takes times
# for (i in 1:sims.per) { # sims per
#   for (s in 1:nrow(surv_scenarios)) { #scenarios picked
#     for (d in 1:nrow(dem_scenarios)) { # simulation scenario
#       
#       out_IPM <- readRDS(paste("out_IPM", "-", d, "-", s, "-", i, ".RDS", sep = ""))
#       tmp <- max(gelman.diag(out_IPM, multivariate = FALSE)[[1]][, 1])
#       
#       if (!is.na(tmp) & tmp <= 1.1) {
#         out_IPM <- out_IPM %>% 
#           collapse_chains() %>% 
#           as.matrix() %>% 
#           as.data.frame() %>% 
#           filter(row_number() %% 60 == 1) %>% # thin chains
#           mutate(dem_scenario = d) %>% 
#           mutate(sim_rep = i) %>% 
#           mutate(surv_scenario = s)
#         assign(paste("out_IPM", "-", d, "-", s, "-", i, sep = ""), out_IPM)
#       }
#       
#       out_noMR <- readRDS(paste("out_noMR", "-", d, "-", s, "-", i, ".RDS", sep = ""))
#       tmp <- max(gelman.diag(out_noMR, multivariate = FALSE)[[1]][, 1])
#       
#       if (!is.na(tmp) & tmp <= 1.1) {
#         out_noMR <- out_noMR %>% 
#           collapse_chains() %>% 
#           as.matrix() %>% 
#           as.data.frame() %>% 
#           filter(row_number() %% 60 == 1) %>% # thin chains
#           mutate(dem_scenario = d) %>% 
#           mutate(sim_rep = i) %>% 
#           mutate(surv_scenario = s)
#         assign(paste("out_noMR", "-", d, "-", s, "-", i, sep = ""), out_noMR)
#       }
#       
#       out_noProd <- readRDS(paste("out_noProd", "-", d, "-", s, "-", i, ".RDS", sep = ""))
#       tmp <- max(gelman.diag(out_noProd, multivariate = FALSE)[[1]][, 1])
#       
#       if (!is.na(tmp) & tmp <= 1.1) {
#         out_noProd <- out_noProd %>% 
#           collapse_chains() %>% 
#           as.matrix() %>% 
#           as.data.frame() %>% 
#           filter(row_number() %% 60 == 1) %>% # thin chains
#           mutate(dem_scenario = d) %>% 
#           mutate(sim_rep = i) %>% 
#           mutate(surv_scenario = s)
#         assign(paste("out_noProd", "-", d, "-", s, "-", i, sep = ""), out_noProd)
#       }
#       
#       out_abundOnly <- readRDS(paste("out_abundOnly", "-", d, "-", s, "-", i, ".RDS", sep = ""))
#       tmp <- max(gelman.diag(out_abundOnly, multivariate = FALSE)[[1]][, 1])
#       
#       if (!is.na(tmp) & tmp <= 1.1) {
#         out_abundOnly <- out_abundOnly %>% 
#           collapse_chains() %>% 
#           as.matrix() %>% 
#           as.data.frame() %>% 
#           filter(row_number() %% 60 == 1) %>% # thin chains
#           mutate(dem_scenario = d) %>% 
#           mutate(sim_rep = i) %>% 
#           mutate(surv_scenario = s)
#         assign(paste("out_abundOnly", "-", d, "-", s, "-", i, sep = ""), out_abundOnly)
#       }
#       
#     } #d
#   } #s
# } #i
# 
# 
# fullIPM <- do.call(bind_rows, lapply( ls(patt="out_IPM"), get))
# noMR <- do.call(bind_rows, lapply( ls(patt="out_noMR"), get))
# noProd <- do.call(bind_rows, lapply( ls(patt="out_noProd"), get))
# abundOnly <- do.call(bind_rows, lapply( ls(patt="out_abundOnly"), get))
# 
# rm(list=grep("out_IPM|out_noMR|out_abundOnly|out_noProd",ls(),value=TRUE,invert=TRUE))
# 
# #rm(list = ls())


