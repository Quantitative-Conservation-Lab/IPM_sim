# load libraries ####

library(here)
library(dplyr)
#library(readxl)
#library(lubridate)
#library(hms)
#library(stringr)
library(tidyr)
library(nimble)

# source files ######

# source data simulation functions
source(here("scripts", "current version","IPMsimData_markasYoYandAd.R")) #changed the MR data
source(here("scripts", "current version","productivityDataSim.R"))
# source model functions
source(here("scripts", "current version","IPMNimble_markedasYoYandAD.R")) #changed the MR data
# source initial value functions
source(here("scripts", "current version","IPMinitvalues.R")) #changed the MR data

# load scenarios ####

# # site
# paramlevels <- read_excel(here("data", "IPM sim spreadsheet.xlsx"), range = "E1:K4")
# paramlevels[, 1] <- c("L", "M", "H")
# colnames(paramlevels)[1] <- "Level"
# scenarios <- read_excel(here("data", "IPM sim spreadsheet.xlsx"), range = "A5:K22")
# scenarios <- scenarios[-c(6, 12), ]
# 
# 
# scenarios <- scenarios %>% 
#   mutate(`MR detection` = case_when(`MR detection` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "MR detection"]),
#                                     `MR detection` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "MR detection"]),
#                                     `MR detection` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "MR detection"]))
#          ) %>% 
#   mutate(`Abund detection` = case_when(`Abund detection` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "Abund detection"]),
#                                     `Abund detection` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "Abund detection"]),
#                                     `Abund detection` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "Abund detection"]))
#   ) %>% 
#   mutate(`Adult Surv` = case_when(`Adult Surv` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "Adult Surv"]),
#                                     `Adult Surv` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "Adult Surv"]),
#                                     `Adult Surv` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "Adult Surv"]))
#   ) %>% 
#   mutate(`Juv Surv` = case_when(`Juv Surv` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "Juv Surv"]),
#                                     `Juv Surv` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "Juv Surv"]),
#                                     `Juv Surv` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "Juv Surv"]))
#   ) %>% 
#   mutate(`Mean Clutch Size` = case_when(`Mean Clutch Size` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "Mean Clutch Size"]),
#                                     `Mean Clutch Size` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "Mean Clutch Size"]),
#                                     `Mean Clutch Size` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "Mean Clutch Size"]))
#   ) %>% 
#   mutate(`Daily nest survival` = case_when(`Daily nest survival` == "L" ~ as.numeric(paramlevels[paramlevels$Level == "L", "Daily nest survival"]),
#                                     `Daily nest survival` == "M" ~ as.numeric(paramlevels[paramlevels$Level == "M", "Daily nest survival"]),
#                                     `Daily nest survival` == "H" ~ as.numeric(paramlevels[paramlevels$Level == "H", "Daily nest survival"]))
#   ) %>% 
#   mutate("MR Included" = ifelse("MR Included" == "Y",1,0)) %>% 
#   mutate(`Abund Included` = ifelse(`Abund Included` == "Y",1,0)) %>% 
#   mutate("Nests Included" = ifelse("Nests Included" == "Y",1,0)) %>%
#   select(-`Sims per`) %>%
#   mutate(`Fec` = 1/2 * `Mean Clutch Size` * `Daily nest survival`^30)
# 
# lams <- apply(scenarios, 1, function(x) {eigen(matrix(data = c(x[8] * x[11], x[8] * x[11], x[7], x[7]), 
#                                               byrow = TRUE, nrow = 2, ncol = 2))$values[1]})
# 
# scenarios <- cbind(scenarios, lams)
# #View(scenarios)
# saveRDS(scenarios, here("data", "scenarios.Rdata"))
scenarios <- readRDS(here("data", "scenarios.Rdata"))

n.scenarios <- max(scenarios$`Scenario Number`)
sims.per <- 25

# AEB note
# parameter values picked such that
# parameter with most uncertainty - juv survival (one we observe most indirectly)
  # how good at we are recovering that 

# run simulations ######

# THINGS THAT NEVER CHANGE
phi.nest <- 0.975
mean.clutch.size <- 2.5
n.sam <- 3
max.nest.age <- 30
n.years=10
n.data=c(200,1000)
init.age = c(1000,1000)
phi.ad = 0.77
p.1 = 0.98

for (s in 15:15) {
  # THINGS THAT DO CHANGE
  phi.1 = scenarios[s, "Juv Surv"]
  p.ad = scenarios[s, "MR detection"]
  p.sur = scenarios[s, "Abund detection"]
  for (i in 1:sims.per) {
    print(paste("scenario: ", s, "; simulation number: ", i, sep = ""))
    if (scenarios[s, "MR Included"] == 1 & scenarios[s, "Nests Included"] == 1) {
      # simulate datasets
      df<-simIPMdata(n.years, n.data, init.age, phi.1, phi.ad, p.1, p.ad, p.sur,
                     max.nest.age, mean.clutch.size, phi.nest, n.sam) 
      prod <- getNestDat()
      for (c in 1:3) {
        # abund data
        y <- df$SUR
        n.sam <- df$n.sam
        # Capture-recapture data (in m-array format, from years 1 to n.years)
        m <- df$ch
        first <- df$first
        age_ch <- df$age_ch
        # Nest data
        H <- prod$observed.nest.status
        Fledged <- prod$clutch.sizes
        first.nest <- prod$first.nest
        last.nest <- prod$last.nest
        max.nest.age <- prod$max.nest.age
        n.nests <- prod$N.nests.found
        n.succ.nests <- prod$N.nests.successful
        # initial values
        age<-ageunknown(age_ch)
        z.state <- state.data(m)
        Hinits <- getHinits(H)
        
        datipm <- list(ch.y = m, y = y, 
                       H = H, 
                       Fledged = Fledged)
        constants<-list(nyears = ncol(m), 
                        n.ind=nrow(m), first=first, age=age, n.sam=n.sam,
                        n.nests = n.nests, 
                        n.succ.nests = n.succ.nests, 
                        first.nest = first.nest, 
                        last.nest = last.nest, 
                        max.nest.age = max.nest.age)
        inits <- list(mean.phi=c(0.4, 0.77),
                      mean.p = 0.5, 
                      #mean.fec = runif(1, 0, 10), 
                      p.surv=0.9,
                      z=z.state,
                      n1.start=round(mean(y[,1]) * 1.5),#sample(1:30,1),#super sensitive to these values, tried rpois(1,30) and it dodnt work
                      nad.start=round(mean(y[,1]) * 1.5),#sample(1:30,1),
                      phi.nest = 0.975, 
                      lambdaf = 2.5 ,
                      H = Hinits
        )
        # run full model
        # IPM
        # AEB note - weird stuff happening with p.surv
        parameters <- c("mean.phi", "mean.p", "fec", "p.surv", "phi.nest", "lambdaf", "lambda")
        modIPM<-nimbleModel(IPMmod, constants=constants, data=datipm, inits=inits)
        confIPM<-configureMCMC(modIPM) 
        confIPM$addMonitors(parameters)
        RmcmcIPM<-buildMCMC(confIPM) 
        CmodelIPM<-compileNimble(modIPM)
        CmcmcIPM<-compileNimble(RmcmcIPM, project=CmodelIPM)
        CmcmcIPM$run(thin=10, reset=T, niter=45000, nburnin=5000)
        out<-as.data.frame(as.matrix(CmcmcIPM$mvSamples))
        
        # save results to file and to environment
        outcopy <- out
        assign(paste("out", s, "_",  i, "-", c, sep = ""), outcopy)
        saveRDS(out, here("data", paste("out", s, "_",  i, "-", c, ".Rdata", sep = "")))
        rm(out, outcopy)
      }
    } else if (scenarios[s, "MR Included"] == 1 & scenarios[s, "Nests Included"] == 0) {
      # simulate datasets
      df<-simIPMdata(n.years, n.data, init.age, phi.1, phi.ad, p.1, p.ad, p.sur,
                     max.nest.age, mean.clutch.size, phi.nest, n.sam) 
      prod <- getNestDat()
      for (c in 1:3) {
        # abund data
        y <- df$SUR
        n.sam <- df$n.sam
        # Capture-recapture data (in m-array format, from years 1 to n.years)
        m <- df$ch
        first <- df$first
        age_ch <- df$age_ch
        # Nest data
        H <- prod$observed.nest.status
        Fledged <- prod$clutch.sizes
        first.nest <- prod$first.nest
        last.nest <- prod$last.nest
        max.nest.age <- prod$max.nest.age
        n.nests <- prod$N.nests.found
        n.succ.nests <- prod$N.nests.successful
        # initial values
        age<-ageunknown(age_ch)
        z.state <- state.data(m)
        Hinits <- getHinits(H)
        
        datipm <- list(ch.y = m, y = y, 
                       H = H, 
                       Fledged = Fledged)
        constants<-list(nyears = ncol(m), 
                        n.ind=nrow(m), first=first, age=age, n.sam=n.sam,
                        n.nests = n.nests, 
                        n.succ.nests = n.succ.nests, 
                        first.nest = first.nest, 
                        last.nest = last.nest, 
                        max.nest.age = max.nest.age)
        inits <- list(mean.phi=c(0.4, 0.77),
                      mean.p = 0.5, 
                      #mean.fec = runif(1, 0, 10), 
                      p.surv=0.9,
                      z=z.state,
                      n1.start=round(mean(y[,1]) * 1.5),#sample(1:30,1),#super sensitive to these values, tried rpois(1,30) and it dodnt work
                      nad.start=round(mean(y[,1]) * 1.5),#sample(1:30,1),
                      phi.nest = 0.975, 
                      lambdaf = 2.5 ,
                      H = Hinits
        )
        # run no nest model
        # NO NEST
        parameters <- c("mean.phi", "mean.p", "fec", "p.surv", "phi.nest", "lambdaf", "lambda")
        modnonest<-nimbleModel(noNests, constants=constants, data=datipm, inits=inits)
        confnonest<-configureMCMC(modnonest) 
        confnonest$addMonitors(parameters)
        Rmcmcnonest<-buildMCMC(confnonest) 
        Cmodelnonest<-compileNimble(modnonest)
        Cmcmcnonest<-compileNimble(Rmcmcnonest, project=Cmodelnonest)
        Cmcmcnonest$run(thin=10, reset=T, niter=10000, nburnin=5000)
        out<-as.data.frame(as.matrix(Cmcmcnonest$mvSamples))
        
        # save results to file and to environment
        outcopy <- out
        assign(paste("out", s, "_",  i, "-", c, sep = ""), outcopy)
        saveRDS(out, here("data", paste("out", s, "_",  i, "-", c, ".Rdata", sep = "")))
        rm(out, outcopy)
      }
    } else if (scenarios[s, "MR Included"] == 0 & scenarios[s, "Nests Included"] == 1) {
      # simulate datasets
      df<-simIPMdata(n.years, n.data, init.age, phi.1, phi.ad, p.1, p.ad, p.sur,
                     max.nest.age, mean.clutch.size, phi.nest, n.sam) 
      prod <- getNestDat()
      for (c in 1:3) {
        # abund data
        y <- df$SUR
        n.sam <- df$n.sam
        # Capture-recapture data (in m-array format, from years 1 to n.years)
        m <- df$ch
        first <- df$first
        age_ch <- df$age_ch
        # Nest data
        H <- prod$observed.nest.status
        Fledged <- prod$clutch.sizes
        first.nest <- prod$first.nest
        last.nest <- prod$last.nest
        max.nest.age <- prod$max.nest.age
        n.nests <- prod$N.nests.found
        n.succ.nests <- prod$N.nests.successful
        # initial values
        age<-ageunknown(age_ch)
        z.state <- state.data(m)
        Hinits <- getHinits(H)
        
        datipm <- list(ch.y = m, y = y, 
                       H = H, 
                       Fledged = Fledged)
        constants<-list(nyears = ncol(m), 
                        n.ind=nrow(m), first=first, age=age, n.sam=n.sam,
                        n.nests = n.nests, 
                        n.succ.nests = n.succ.nests, 
                        first.nest = first.nest, 
                        last.nest = last.nest, 
                        max.nest.age = max.nest.age)
        inits <- list(mean.phi=c(0.4, 0.77),
                      mean.p = 0.5, 
                      #mean.fec = runif(1, 0, 10), 
                      p.surv=0.9,
                      z=z.state,
                      n1.start=round(mean(y[,1]) * 1.5),#sample(1:30,1),#super sensitive to these values, tried rpois(1,30) and it dodnt work
                      nad.start=round(mean(y[,1]) * 1.5),#sample(1:30,1),
                      phi.nest = 0.975, 
                      lambdaf = 2.5 ,
                      H = Hinits
        )
        # NO MR
        parameters <- c("mean.phi", "fec", "p.surv", "phi.nest", "lambdaf", "lambda")
        modnoMR<-nimbleModel(noMR, constants=constants, data=datipm, inits=inits)
        confnoMR<-configureMCMC(modnoMR) 
        confnoMR$addMonitors(parameters)
        RmcmcnoMR<-buildMCMC(confnoMR) 
        CmodelnoMR<-compileNimble(modnoMR)
        CmcmcnoMR<-compileNimble(RmcmcnoMR, project=CmodelnoMR)
        CmcmcnoMR$run(thin=10, reset=T, niter=10000, nburnin=5000)
        out<-as.data.frame(as.matrix(CmcmcnoMR$mvSamples))
        
        # save results to file and to environment
        outcopy <- out
        assign(paste("out", s, "_",  i, "-", c, sep = ""), outcopy)
        saveRDS(out, here("data", paste("out", s, "_",  i, "-", c, ".Rdata", sep = "")))
        rm(out, outcopy)
      }  
    }
  }
}

# s <- 1 # for full test
# s <- 11 # for no nests test
# s <- 6 # for no nests test

# TODO
# try one of each
# put up on loon...
# could start many instances of r
# 15 instances of R - 1 for each simulation scenario

# visualize results #####
# TODO
# fill this in

# create tables #####
# TODO
# fill this in

