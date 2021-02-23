#example for how to check for convergence and update model runs
#using an occupancy model and a bit of clunky coding!

library(nimble)
################simulate data
#simulate some occupancy data for a super simple occupany model
#50 sites, 3 surveys
true_z<-rbinom(50,1,0.8) #true state of each site
obs_y<-matrix(nrow=50,ncol=3)
for(i in 1:50){
  #observed state of each site at each survey
  #0.6 is the probability of detecting whatever species during a survey
  obs_y[i,]<-rbinom(3,1, true_z[i]*0.6)
}

#make the observation data
dat_sim<-list(y=obs_y)

############nimble stuff
#nimble model
mod<-nimbleCode({
  p~dunif(0,1)
  phi~dunif(0,1)
  for(i in 1:N){
    z[i]~dbern(phi)
  }
  for(i in 1:N){
    for(j in 1:J){
      zp[i,j]<-z[i]*p
      y[i,j]~dbern(zp[i,j])
    }
  }
})

#function to grab initial values
initsfn<-function(){
  p<-runif(1,0,1)
  phi<-runif(1,0,1)
  return(list(p=p,phi=phi))
}

##################################
#Here is a clunky version
#these commanands are for all chains
constants<-list(N=50, J=3)
Rmodel<-nimbleModel(code=mod, constants=constants)
conf<-configureMCMC(Rmodel)
conf$addMonitors(c("phi","p"))

#then need to set it up to run different chains with different initial values
#so one for each chain
Rmcmc1<-buildMCMC(conf)
Rmcmc2<-buildMCMC(conf)

Cmodel1<-compileNimble(Rmodel)
Cmodel2<-compileNimble(Rmodel)

Cmodel1$setData(dat_sim)
Cmodel2$setData(dat_sim)

Cmodel1$setInits(inits=initsfn())
Cmodel2$setInits(inits=initsfn())

Cmcmc1<-compileNimble(Rmcmc1, project=Cmodel1)
#run the chain
Cmcmc1$run(thin=10, reset=T, niter=1000,nburnin=500)
#save output for this chain
x11<-(as.data.frame(as.matrix(Cmcmc1$mvSamples)))

#same as above
Cmcmc2<-compileNimble(Rmcmc2, project=Cmodel2)
Cmcmc2$run(thin=10, reset=T, niter=1000,nburnin=500)

x22<-(as.data.frame(as.matrix(Cmcmc2$mvSamples)))

library(coda)
#then check the rhat value for what was run
xone<-as.mcmc(x11)
xtwo<-as.mcmc(x22)

out.mcmc <- as.mcmc.list(list(xone,xtwo))
v<-gelman.diag(out.mcmc, multivariate = T, transform = F)
#check it out
v

#so then if they havent converged pick up running the chain from where they left off
#to do this, you do 'reset=F' and it will pick up where it left off
Cmcmc1$run(thin=10, reset=F, niter=7000)
#then save the updated chain
x11<-(as.data.frame(as.matrix(Cmcmc1$mvSamples)))

#same as above
Cmcmc2$run(thin=10, reset=F, niter=7000)
x22<-(as.data.frame(as.matrix(Cmcmc2$mvSamples)))

#then recheck
xone<-as.mcmc(x11)
xtwo<-as.mcmc(x22)

out.mcmc <- as.mcmc.list(list(xone,xtwo))
v<-gelman.diag(out.mcmc)
v
#if it hasnt converged then do it again with longer run

##########################
#so you could set this up this way too, which would just automate checking
#it is still a little clunky and you have to run each for each chain this way
constants<-list(N=50, J=3)
Rmodel<-nimbleModel(code=mod, constants=constants)
conf<-configureMCMC(Rmodel)
conf$addMonitors(c("phi","p"))

#then need to set it up to run different initial values
#so one for each chain
Rmcmc1<-buildMCMC(conf)
Rmcmc2<-buildMCMC(conf)

Cmodel1<-compileNimble(Rmodel)
Cmodel2<-compileNimble(Rmodel)

Cmodel1$setData(dat_sim)
Cmodel2$setData(dat_sim)

Cmodel1$setInits(inits=initsfn())
Cmodel2$setInits(inits=initsfn())

Cmcmc1<-compileNimble(Rmcmc1, project=Cmodel1)
Cmcmc2<-compileNimble(Rmcmc2, project=Cmodel2)
library(coda)

#run the chains
Cmcmc1$run(thin=10, reset=T, niter=1000,nburnin=500)
Cmcmc2$run(thin=10, reset=T, niter=1000,nburnin=500)
#save for this chain
x11<-(as.data.frame(as.matrix(Cmcmc1$mvSamples)))
x22<-(as.data.frame(as.matrix(Cmcmc2$mvSamples)))

xone<-as.mcmc(x11)
xtwo<-as.mcmc(x22)

out.mcmc <- as.mcmc.list(list(xone,xtwo))
v<-gelman.diag(out.mcmc, multivariate = T, transform = F)
v
#hasnt converged, so update until it does
while(v$mpsrf>1.1){
  Cmcmc1$run(thin=10, reset=F, niter=1000)
  Cmcmc2$run(thin=10, reset=F, niter=1000)
  x11<-(as.data.frame(as.matrix(Cmcmc1$mvSamples)))
  x22<-(as.data.frame(as.matrix(Cmcmc2$mvSamples)))
  
  xone<-as.mcmc(x11)
  xtwo<-as.mcmc(x22)
  
  out.mcmc <- as.mcmc.list(list(xone,xtwo))
  v<-gelman.diag(out.mcmc, multivariate = T, transform = F)
}
#again, there are probably 'slicker' ways to do this and wouldn't be ideal if 
#model runs are substantial in time or memory.
#it would be cool to figure out how to do it in parallel.
#What I do is just do one instance per chain and then check in a separate
#r instance while keeping the other instances open
#then if they haven't converged, I do 'Cmcmc$run(thin=10, reset=F, niter=1000)'
