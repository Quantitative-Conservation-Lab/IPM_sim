# load libraries ####

library(here)
library(dplyr)
library(readxl)
library(lubridate)
library(hms)
library(stringr)
library(tidyr)

# source files ######

# source data simulation function
# source model function

# load scenarios ####

# site
paramlevels <- read_excel(here("scripts", "dataSims", "IPM sim spreadsheet.xlsx"), range = "E1:K4")
paramlevels[, 1] <- c("L", "M", "H")
colnames(paramlevels)[1] <- "Level"
scenarios <- read_excel(here("scripts", "dataSims", "IPM sim spreadsheet.xlsx"), range = "A5:K22")
scenarios <- scenarios[-c(6, 12), ]

n.scenarios <- max(scenarios$`Scenario Number`)
sims.per <- max(scenarios$`Sims per`)

scenarios <- scenarios %>% 
  mutate(`MR detection` = case_when(`MR detection` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "MR detection"]),
                                    `MR detection` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "MR detection"]),
                                    `MR detection` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "MR detection"]))
         ) %>% 
  mutate(`Abund detection` = case_when(`Abund detection` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "Abund detection"]),
                                    `Abund detection` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "Abund detection"]),
                                    `Abund detection` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "Abund detection"]))
  ) %>% 
  mutate(`Adult Surv` = case_when(`Adult Surv` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "Adult Surv"]),
                                    `Adult Surv` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "Adult Surv"]),
                                    `Adult Surv` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "Adult Surv"]))
  ) %>% 
  mutate(`Juv Surv` = case_when(`Juv Surv` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "Juv Surv"]),
                                    `Juv Surv` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "Juv Surv"]),
                                    `Juv Surv` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "Juv Surv"]))
  ) %>% 
  mutate(`Mean Clutch Size` = case_when(`Mean Clutch Size` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "Mean Clutch Size"]),
                                    `Mean Clutch Size` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "Mean Clutch Size"]),
                                    `Mean Clutch Size` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "Mean Clutch Size"]))
  ) %>% 
  mutate(`Daily nest survival` = case_when(`Daily nest survival` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "Daily nest survival"]),
                                    `Daily nest survival` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "Daily nest survival"]),
                                    `Daily nest survival` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "Daily nest survival"]))
  ) %>% 
  mutate(`MR Included` = ifelse(`MR Included` == "Y",1,0)) %>% 
  mutate(`Abund Included` = ifelse(`Abund Included` == "Y",1,0)) %>% 
  mutate(`Nests Included` = ifelse(`Nests Included` == "Y",1,0)) %>%
  select(-`Sims per`)

# run simulations ######

for (s in 1:n.scenarios) {
  for (i in 1:sims.per) {
    if (scenarios[s, `MR Included`] == 1 & scenarios[s, `Nests Includedd`] == 1) {
      # run full model
    } else if (scenarios[s, `MR Included`] == 1 & scenarios[s, `Nests Includedd`] == 0) {
      # run no nest model
    } else if (scenarios[s, `MR Included`] == 0 & scenarios[s, `Nests Includedd`] == 1) {
      # run no resight model
    }
    
    # save results to file and to environment
    outcopy <- out
    assign(paste("out", s, "_",  i, sep = ""), outcopy)
    saveRDS(out, paste("out", s, "_",  i, ".Rdata", sep = ""))
    rm(out, outcopy)
  }
}

# visualize results #####


# create tables #####


