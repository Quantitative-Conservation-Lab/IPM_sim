# load libraries ####

library(here)
library(dplyr)
#library(readxl)
#library(lubridate)
#library(hms)
#library(stringr)
library(tidyr)
library(nimble)

# source functions

# load simulated data ####

# TODO
# set this code up to be parallelized

# 324 different scenarios
# how to split them up

# Most efficient way to run things
# loon has 20 cores

# using 12 at a time
  # which datasets present x lambda level

# https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html
# ask hannah for example 
# also ask hannah for example of running models until all converge

# then have to run all detection levels - 27



for (s in 2:2) {
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
        Cmcmcnonest$run(thin=10, reset=T, niter=45000, nburnin=5000)
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
        CmcmcnoMR$run(thin=10, reset=T, niter=45000, nburnin=5000)
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

