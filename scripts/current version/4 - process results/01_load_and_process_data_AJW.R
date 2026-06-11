library(tidyverse)
library(here)
library(nimble)
library(foreach)
library(doParallel)
library(coda)
library(nlist)
# library(beepr)

# load scenarios
dem_scenarios <- readRDS(here("data", "demographic_scenarios.RDS")) %>% 
  separate_wider_delim(cols = scenario, delim = ",", 
                       names = c("life_hist", "trend")) %>% 
  rename("phi1" = "S.J", 
    "phiad" = "S.A",
    "fec" = "f")

surv_scenarios <- readRDS(here("data", "data_scenarios.RDS"))
sims.per <- 100


# storage
results_list <- list()
missing_files <- data.frame(type=character(), d=integer(), s=integer(), i=integer())
not_converged <- data.frame(type=character(), d=integer(), s=integer(), 
                                 i=integer(), gelman=numeric(), max_param=character())

# processing function (returns the data frame or NULL)
process_model <- function(prefix, d, s, i) {
  # file_name <- here('results', 'ind300', paste0(prefix, "-", d, "-", s, "-", i, ".RDS"))
  # file_name <- here('results', 'ind400', paste0(prefix, "-", d, "-", s, "-", i, ".RDS"))
  # file_name <- here('results', 'nmix', 'ind300_nsam5',
  file_name <- here('results', 'nmix', 'final', #'ind300_nsam5',
                    paste0(prefix, "-", d, "-", s, "-", i, ".RDS"))
  
  if (!file.exists(file_name)) {
    # <<- updates tracking items outside the function
    missing_files <<- bind_rows(missing_files, data.frame(type=prefix, d=d, s=s, i=i))
    return(NULL)
  }
  
  out_temp <- readRDS(file_name)
  diag_result <- gelman.diag(out_temp, multivariate = FALSE)[[1]]
  # Find the value and the name of the worst offender
  max_idx   <- which.max(diag_result[, 1])
  max_val   <- diag_result[max_idx, 1]
  max_param <- rownames(diag_result)[max_idx]
  
  # compile convergence info
  if (is.na(max_val) || max_val > 1.12) {
    not_converged <<- bind_rows(not_converged, 
                                     data.frame(
                                       type=prefix, d=d, s=s, i=i, 
                                       gelman=max_val, 
                                       max_param=max_param
                                     )) %>% distinct()
    return(NULL)
  }
  
  # process
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

saveRDS(results_all, file = here('results', 'processed', "results_all_final_delphine_batch1.RDS"))
# saveRDS(results_all, file = here('results', 'processed', "results_all_vSJC.RDS"))

#this still isn't quite right.... something like this but accounting for everything in missing_files
var_reps <- n_distinct(results_all$sim_rep)*n_distinct(results_all$dem_scenario)*n_distinct(results_all$surv_scenario)

convergence_summary <- results_all %>%
  distinct(model_type, dem_scenario, surv_scenario, sim_rep) %>%
  group_by(model_type) %>%
  dplyr::summarize(
    successful_sims = n(), 
    # calculate percentage based on 'sims.per' variable
    percent_converged = (n()/var_reps)*100,
    .groups = "drop")

# not_conv_summary <- not_converged %>%
#   group_by(type) %>%
#   dplyr::summarize(n = nrow(gelman))

# abund_conv <- results_all %>%
#   filter(model_type == 'out_abundOnly') %>%
#   distinct(surv_scenario, dem_scenario, sim_rep)

library(ggplot2)

#total scenarios finished or proportion per model type
conv <- sum(convergence_summary$successful_sims)

prop_not_conv <- dim(not_converged)[1]/conv

# parameters that are failing to converge
ggplot(not_converged, aes(x = reorder(max_param, max_param, function(x) -length(x)))) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  facet_grid(.~type) +
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
testIPM.b <- readRDS(file = here('results', 'nmix', 'final',
                                 # 'ind400',
                                 'out_IPM-6-25-7.RDS'))
testIPM.w <- readRDS(file = here('results', 'nmix', 'final', 'out_IPM-4-17-4.RDS')) #best
testIPM.w <- readRDS(file = here('results', 'vSJC', 'out_IPM-1-13-6.RDS')) #best


plot(testIPM.b[,'mean.phi[1]'])
plot(testIPM.b[,'p.surv'])
plot(testIPM.b[,'Ntot[4]'])
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
testNoMR.w <- readRDS(file = here('results', 'nmix', 'final', 'out_noMR-6-15-5.RDS'))


plot(testNoMR.w[,'mean.phi[1]'])
plot(testNoMR.w[,'mean.phi[2]'])
plot(testNoMR.w[,'fec'])
plot(testNoMR.w[,'p.surv'])
plot(testNoMR.w[,'Ntot[10]'])

#struggle params - all vitals, nothing else; 1-47-1 best
testabund.w  <- readRDS(file = here('results', 'nmix', 'final', 'out_abundOnly-6-46-1.RDS'))
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


