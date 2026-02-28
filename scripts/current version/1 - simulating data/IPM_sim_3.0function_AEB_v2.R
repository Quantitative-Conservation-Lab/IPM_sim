# function to simulate the trajectory of a population
# adapted from Kery and Schaub, 2011

# Function input definitions
# n.years = number of years to simulate,
# n.data.types = vector with proportion of pop observed for each dataset,
# age.init = vector of starting age structure/number in each age class, i.e.,
#           = c(1st years, adults), 
# phi.1 = juv survival,  
# phi.ad = adult survival,
# f = fecundity

simPopTrajectory <- function(n.years, age.init,
                             phi.1, phi.ad, f){
  
  ####### HELPER FUNCTIONS  ####### 
  # states: [1yrolds, Adult, chicks, Dead]
  # indfates is an array [states,year,individual]
  
  # fate of juveniles; question/check for Hannah - line 25 is 'time' and all the others are 'time+1'
  oneyr_fatefn<-function(ind,time){
    indfates[3,time,ind]<-rpois(1,fec) # number of chicks produced by ind. in current year
    zsurv1<-rbinom(1,1,phi.ad)
    indfates[2,time+1,ind]<-ifelse(zsurv1==1,1,NA) #do they survive?
    indfates[4,time+1,ind]<-ifelse(zsurv1==0, 1,NA) #or die?
    indfates[1,time+1,ind]<-NA #cant stay 1yearold
    return(indfates[,time:(time+1),ind])
  }
  
  # fate of adults
  adfatefn<-function(ind,time){
    indfates[3,time,ind]<-rpois(1,fec) # number of chicks produced by ind. in current year
    zsurv1<-rbinom(1,1,phi.ad)
    indfates[2,time+1,ind]<-ifelse(zsurv1==1,1,NA) #do they survive?
    indfates[4,time+1,ind]<-ifelse(zsurv1==0, 1,NA) #or die?
    indfates[1,time+1,ind]<-NA #can't age backwards
    return(indfates[,time:(time+1),ind])
  }
  
  # if you are dead you stay dead
  deadfn<-function(ind,time){
    indfates[4,time+1,ind]<-1
    indfates[1:3,time+1,ind]<-NA
    return(indfates[,(time+1),ind])
  }
  
  # did the chicks produced above survive to the next time step?
  chickfatefn<-function(ind, time){
    zsurv1<-rbinom(1,1,phi.1)
    indfates[1,time+1,ind]<-ifelse(zsurv1==1,1,NA) #do they survive to become 1yrolds
    indfates[4,time+1,ind]<-ifelse(zsurv1==0,1,NA) #or die?
    indfates[2,time+1,ind]<-NA #not adult yet
    indfates[3,time+1,ind]<-NA
    return(indfates[,(time+1),ind])
  }
  
  ####### END HELPER FUNCTIONS  ####### 
  
  nplus1.years<-n.years+1
  fec <- f 
  
  #Leslie matrix for population size
  N<-matrix(nrow=length(age.init), ncol=nplus1.years)
  N[,1]<-age.init
  lesmat<-matrix(nrow=2, ncol=2)
  lesmat[1,]<-c(phi.1*fec, phi.1*fec)
  lesmat[2,]<-c(phi.ad, phi.ad)
  # populate Leslie matrix
  for(t in 1:n.years){
    N[,(t+1)]<-lesmat%*%N[,t]
  }
  
  # give a warning if param values are outside of desired pop. growth rates
  el<-eigen(lesmat)$values[1]
  
  no.animals<-sum(N) #number of animals ever in the system at anytime
  no.ani.max<-round(no.animals*5) #include more for simulation of offspring
  ####stable age distribution is:
  #only want to look at the end distribution
  sad<-N[,nplus1.years]/sum(N[,nplus1.years]) # the proportion of the population that is in each age class
  #lambda from N, matches eigenvalue
  #note that this is without demographic stochasticity and
  #when fates are simulated the actual lambda will change slightly
  
  ######
  #simulate what happens to each individual:
  ######
  
  #dont run if there arent any individuals given
  if(N[1]==0 || N[2]==0){print("no individuals!!!!!!!!!")}
  
  #make an array for each individual, for each year, what state it is in
  # states: [1yrolds, Adult, chicks, Dead]
  # Oneyearolds are the 1st years, chicks are how many chicks per individual,
  #dead is when they die
  
  #simulate their fates
  #set up array for individuals
  #first dimension is age groups; juvenile, adult, chicks, dead
  indfates<-array(dim=c(4, nplus1.years, no.ani.max))
  #add in the stable age distribution to start from, given our total starting population
  age1<-round(sum(age.init)*sad[1])
  age2<-round(sum(age.init)*sad[2])
  indfates[1,1,1:age1]<-1 #number of one year olds in year 1
  indfates[2,1,(age1+1):sum(age1+age2)]<-1 #number of adults in year 1
  
  inpop<-numeric(n.years) #track the population over time
  inpop[1]<-sum(age1+age2) #first year total, stable age distribution way to do it
  
  # simulate ind. trajectories with dem stochasticity
  chickst<-numeric(n.years) # tracking number of chicks produced each year
  tempstep<-matrix(nrow=no.ani.max, ncol=(n.years+1))
  for(t in 1:n.years){
    for(i in 1:inpop[t]){
      tempstep[i,t]<-which(indfates[c(1,2,4),t,i]==1)
      if(tempstep[i,t]==1){ #if first years
        indfates[,(t:(t+1)),i]<-oneyr_fatefn(ind = i, time = t)#[,(t:(t+1))]
      }else if(tempstep[i,t]==2){ #if adults
        indfates[,(t:(t+1)),i]<-adfatefn(ind = i, time = t)#[,(t:(t+1))]
      } else if(tempstep[i,t]==3){ #if dead
        indfates[,(t+1),i]<-deadfn(ind = i, time = t)#[,(t:(t+1))]
      }
    } #i loop over population size at time t
    
    #see how many chicks were produced and add them to the population to track
    chickst[t]<-sum(indfates[3,t,], na.rm=T)
    inpop[t+1]<-inpop[t]+chickst[t]
    if(inpop[t+1]>inpop[t]){ # if pop is growing
      indfates[3,t,(inpop[t]+1):inpop[t+1]]<-1 # add chicks to pop
      for(i in (inpop[t]+1):inpop[t+1]){
        #if chicks were produced, do they survive
        indfates[,(t+1),i]<-chickfatefn(ind = i, time = t)#[,(t:(t+1))]
      }
    }else{}
  }
  
  #check observed lambda
  ad<-j<-numeric(n.years)
  ch<-dead<-numeric(n.years)
  for(i in 1:n.years){
    ad[i]<-sum(indfates[2,i,], na.rm=T)
    ch[i]<-sum(indfates[3,i,], na.rm=T)
    j[i]<-sum(indfates[1,i,], na.rm=T)
    dead[i]<-sum(indfates[4,i,], na.rm=T)
  }
  
  adj<-numeric(length(ad))
  for(i in 2:length(adj)){
    adj[i]<-sum(ad[i]+j[i])/sum(ad[i-1]+j[i-1])
  }
  #mean(adj[2:10])
  #close *enough* to actual lambda, again demographic stochasticity is the cause!
  
  #we dont want to output the leslie matrix N, we want to output the N
  #from the indfates
  Nouts<-matrix(nrow=2,ncol=n.years)
  for(t in 1:n.years){
    Nouts[1,t]<-sum(indfates[1,t,], na.rm=T)
    Nouts[2,t]<-sum(indfates[2,t,],na.rm=T)
  }
  #so we would compare the model estimated Nouts to see if it is tracking it
  #we probably dont care about this for actual analysis, but for checking it is good to have on hand
  return(list(indfates = indfates, Nouts = Nouts))
  
}


# Definitions ####
# indfates = matrix of individual fates, generated using simPopTrajectory above
# n.years = n years to simulate,
# n.data.types = proportion of pop observed for each dataset,
# ADonly = binary variable indicating whether only adults are marked or not
# p.1 = detection of 1 year old birds
# p.ad = detection of adult birds
# BinMod = binary variable indicating whether to use a binomial observation model on repeated counts
# n.sam = number of repeated counts to take
# p.count = probability of detecting individuals on repeated counts
# sig = used if BinMod = F, sd parameter for normal observation model on repeated counts
# productivity = binary variable indicating whether or not to use a nest survival model. deprecated, always T. 
# p.prod = probability of detection for productivity model
simData <- function(indfates, n.years,
                    p.1, p.ad,
                    BinMod, n.sam, p.count, sig, productivity, p.prod){
  
  ################### For output of data  ###################
  #randomly sampling from whole indfates based on detection probabilities or 
  #other sampling probabilities (no longer dividing up the population into 
  #distinct sets for each data type)
  
  IND_MR <- indfates
  IND_Count <- indfates
  IND_Nest <- indfates
  
  ###################Create Mark-Resight data  ###################
  
  #Marking adults and 1 year olds, which will have the same
  #for adults and 1yr olds only
  
  # TODO ugh something is wrong here
  
  #for adults, 1 year olds and chicks
  if(!is.na(p.1) & !is.na(p.ad)) {
    
    # ok yeah, this is just grabbing new and established breeders
    ind_mr.a<-IND_MR[1:3,1:n.years,] 
    ind_mr.j<-IND_MR[1:3,1:n.years,] 
    
    rm.a<-numeric(dim(ind_mr.a)[3]) #vector to remove individuals
    rm.j<-numeric(dim(ind_mr.j)[3]) #vector to remove individuals
    for(i in 1:dim(ind_mr.a)[3]){
      
      # at least one occasion where there is an entry for chick and no entry for 1:2
      if(sum(!is.na(ind_mr.j[3,,i]) & ind_mr.j[3,,i]==1 & apply(ind_mr.j[1:2,,i], 2, function(x)all(is.na(x))))){
        if (rbinom(1,1,p.prod) == 1) {
          rm.j[i]<-0 # don't remove, up for marking as a chick
        } else {
          rm.j[i] <- 1
        }
      }else{
        rm.j[i]<-1 # otherwise remove
      }
      # entry either as new or established breeder, 
      # and is not in the juvenile array
      if(length(which(!is.na(ind_mr.a[1:2,,i])))>0 & rm.j[i]==1){
        rm.a[i]<-0
      }else{
        rm.a[i]<-1 # otherwise remove
      }
    }
    
    # sample subset
    if(sum(rm.a>0)){
      ind_mr.a<-ind_mr.a[,,-which(rm.a==1)]
      # use some prob of capture same as subsequent resight...
      samp.a<-sample(1:(dim(ind_mr.a)[3]), round(p.ad*(dim(ind_mr.a)[3]))) %>% sort()
      ind_mr.a<-ind_mr.a[,,samp.a]
    }else{}
    if(sum(rm.j>0)){
      ind_mr.j<-ind_mr.j[,,-which(rm.j==1)]
      # AEB note - let's band from nests...
      # samp.j<-sample(1:(dim(ind_mr.j)[3]), round(p.prod*(dim(ind_mr.j)[3]))) %>% sort()
      # ind_mr.j<-ind_mr.j[,,samp.j]
    }else{}
    
    age.a<-first.a<-last.a<-numeric() # age, first and last encounters
    age.j<-first.j<-last.j<-numeric() # age, first and last encounters
    
    mr_t.a<-mr_t.j<-dim(ind_mr.a)[2]
    mr_ind.a<-dim(ind_mr.a)[3]
    mr_ind.j<-dim(ind_mr.j)[3]
    for(i in 1:mr_ind.a){
      g <- which(!is.na(ind_mr.a[1:2,,i]), arr.ind = TRUE) # ind that were seen
      age.a[i] <- g[1,1] # age at marking
      first.a[i] <- g[1,2] # first time seen
      h <- which(ind_mr.a[1:2,,i]==1, arr.ind = TRUE) # last time seen
      last.a[i] <- max(h[,2])
    }
    ch.true.a<-matrix(0,ncol=mr_t.a, nrow=mr_ind.a)
    for(i in 1:mr_ind.a){
      ch.true.a[i,first.a[i]:last.a[i]]<-1
    }
    for(i in 1:mr_ind.j){
      g <- which(!is.na(ind_mr.j[3,,i]))
      age.j[i] <- 3 #at marking (1=1year olds, 2=adults, 3= chicks)
      first.j[i] <- g[1] 
      h <- which(ind_mr.j[3,1:n.years,i]==1) 
      last.j[i] <- max(h)
    }
    ch.true.j<-matrix(0,ncol=mr_t.j, nrow=mr_ind.j)
    for(i in 1:mr_ind.j){
      ch.true.j[i,first.j[i]:last.j[i]]<-1
    }
    
    #detection of true marked individuals:
    in.mark.a<-ch.a<-matrix(0,nrow=mr_ind.a, ncol=mr_t.a)
    for(i in 1:mr_ind.a){
      # AEB change - conditioning on first capture, no p.ad needed here
      in.mark.a[i,first.a[i]]<-rbinom(1,1,ch.true.a[i,first.a[i]])
      if(first.a[i]==mr_t.a | first.a[i] == last.a[i]) {
        ch.a[i,]<-in.mark.a[i,]
      } else {
        for(t in (first.a[i]+1):last.a[i]){
          in.mark.a[i,t]<-rbinom(1,1,p.ad*ch.true.a[i,t])
        }
        ch.a[i,]<-in.mark.a[i,]
      }
    }
    in.mark.j<-ch.j<-matrix(0,nrow=mr_ind.j, ncol=mr_t.j)
    for(i in 1:mr_ind.j){
      in.mark.j[i,first.j[i]]<-rbinom(1,1,ch.true.j[i,first.j[i]])
      if(first.j[i]==mr_t.j | first.j[i] == last.j[i]) {
        ch.j[i,]<-in.mark.j[i,]
      } else {
        for(t in (first.j[i]+1):last.j[i]){
          in.mark.j[i,t]<-rbinom(1,1,p.ad*ch.true.j[i,t])
        }
        ch.j[i,]<-in.mark.j[i,]
      }
    }
    
    #code to track the ages
    add_age_chtrue.a<-age_ch.a<-matrix(0,nrow=mr_ind.a, ncol=mr_t.a)
    for(i in 1:mr_ind.a){
      for(x in (first.a[i]):last.a[i]){
        add_age_chtrue.a[i,x]<-as.numeric(which(!is.na(ind_mr.a[1:2,x,i])))
        age_ch.a[i,x]<-add_age_chtrue.a[i,x]*ch.a[i,x]
      }
    }
    # AEB - cutting section removing those we never saw - should see everyone at first cap
    firstobs.a<-lastobs.a<-numeric(length(ch.a[,1]))
    for(i in 1:length(ch.a[,1])){
      firstobs.a[i]<-min(which(ch.a[i,]==1))
      lastobs.a[i]<-max(which(ch.a[i,]==1))
    }
    add_age_chtrue.j<-age_ch.j<-matrix(0,nrow=mr_ind.j, ncol=mr_t.j)
    for(i in 1:mr_ind.j){
      add_age_chtrue.j[i,first.j[i]] <- 3
      if(first.j[i]==mr_t.j) next
      add_age_chtrue.j[i,first.j[i]+1] <- 1
      if((first.j[i]+1)==mr_t.j) next
      add_age_chtrue.j[i,(first.j[i]+2):mr_t.j] <- 2
      for(x in (first.j[i]):last.j[i]){
        #add_age_chtrue[i,x]<-as.numeric(which(!is.na(ind_mr[,x,i])))[1]
        age_ch.j[i,x]<-add_age_chtrue.j[i,x]*ch.j[i,x]
      }
    }
    firstobs.j<-lastobs.j<-numeric(length(ch.j[,1]))
    for(i in 1:length(ch.j[,1])){
      firstobs.j[i]<-min(which(ch.j[i,]==1))
      lastobs.j[i]<-max(which(ch.j[i,]==1))
    }
    
  } else {
    ch.a=NULL
    ch.j=NULL
    age_ch.a=NULL
    age_ch.j=NULL
    firstobs.a=NULL
    firstobs.j=NULL
    lastobs.a=NULL
    lastobs.j=NULL
  }
  
  
  ################### Create population survey data ###################
  #n.sam is the number of times that the population was sampled in a year
  #Need to choose in the function argument which model to use, BinMod=T if binom, F if normal
  TRUE_Count<-matrix(nrow=2, ncol=n.years) #first row, number of one year-olds; second row is adults
  SUR<-matrix(nrow=n.sam, ncol=n.years)#matrix for surveys by survey and year
  for(u in 1:n.years){
    #Some discussion about should we include 1yearolds
    #if so, remove the next line
    TRUE_Count[1,u]<-sum(IND_Count[1,u,], na.rm = T)
    TRUE_Count[2,u]<-sum(IND_Count[2,u,], na.rm = T)
    if(BinMod==T){
      SUR[,u]<-rbinom(n.sam, sum(TRUE_Count[,u]), p.count)
    }else{
      SUR[,u]<-rnorm(n.sam, TRUE_Count[,u], sig)
    }
    
  }
  
  ################### Create reproductive success data ###################
  TRUE_nestlings<-numeric(n.years)
  R_true<-R_obs <- numeric(n.years) # number of pairs whose productivity was observed
  OBS_nestlings <- numeric(n.years)  # total number of nestlings recorded in a year
  tomark<-matrix(0,nrow=n.years, ncol=dim(indfates)[3])
  for(t in 1:n.years){
    #true number of nestlings
    TRUE_nestlings[t]<-sum(IND_Nest[3,t,], na.rm=T)
    #true number of reproducing females
    R_true[t]<-length(which(IND_Nest[3,t,]>0))
    #track which are observed, to mark
    
    if(!is.na(p.prod)) {
      #observation, who dont we see and who do we see
      #Assumes we can count the number of chicks perfectly if we see the nest
      #HAS: corrected this
      #before it may have been double counting
      for(i in 1:length(IND_Nest[1,1,])){
        if(!is.na(IND_Nest[3,t,i]) & length(which(!is.na(IND_Nest[1:2,t,i])))>=1){
          obstemp<-rbinom(1,1,p.prod)
          if(obstemp==1){
            R_obs[t]<-R_obs[t]+1
            OBS_nestlings[t]<-IND_Nest[3,t,i]+OBS_nestlings[t]
            
          }
        }
      }
    }
  }
  
  
  
  return(list(ch.a=ch.a, ch.j=ch.j,
              age_ch.a=age_ch.a,
              age_ch.j=age_ch.j,
              firstobs.j=firstobs.j, lastobs.j=lastobs.j,
              firstobs.a=firstobs.a, lastobs.a=lastobs.a,
              SUR=SUR,
              n.sam = n.sam, R_obs=R_obs, OBS_nestlings=OBS_nestlings))
  
}

temptrue.stable<-simPopTrajectory(n.years=15,
                           age.init=c(200,200),
                           phi.1=0.3,phi.ad=0.4,
                           f=2)
temp.data<-simData(indfates = temptrue.stable$indfates, 
                   n.years=15,
                   p.1=1,
                   p.ad=1,
                   p.prod=1,
                   n.sam=3,
                   p.count=0.5,
                   BinMod=T)
source(here("scripts", "current version",
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "3 - run models", "run_scenarios_helperFns.R"))

nb <- 10000#0 #burn-in # TODO play with this
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

comb <- c(0.3, 0.4, 2) %>% t() %>% as.data.frame()
colnames(comb) <- c("phi1", "phiad", "fec") 
detect <- c(0.5, 1) #p.surv, mean.p
out <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, 
                    popDat = temp.data, 
                    popTraj = temptrue.stable, 
                    comb, 
                    detect)

library(postpack)
summ <- post_summ(out, get_params(out), Rhat = TRUE, neff = TRUE) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
# TODO seems to be confounding between phi1 and phi2?????

# trying again without perfect detection
temp.data<-simData(indfates = temptrue.stable$indfates, 
                   n.years=15,
                   p.1=0.8,
                   p.ad=0.8,
                   p.prod=0.5,
                   n.sam=3,
                   p.count=0.5,
                   BinMod=T)

nb <- 10000#0 #burn-in # TODO play with this
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

comb <- c(0.3, 0.4, 2) %>% t() %>% as.data.frame()
colnames(comb) <- c("phi1", "phiad", "fec") 
detect <- c(0.5, 0.8) #p.surv, mean.p
out2 <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, 
                 popDat = temp.data, 
                 popTraj = temptrue.stable, 
                 comb, 
                 detect)

library(postpack)
summ2 <- post_summ(out2, get_params(out2), Rhat = TRUE, neff = TRUE) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
# huh, now this looks better???

# TODO maybe just need to try with larger sample sizes??
# nope, that didn't work
# TODO - also try running just the marray model

marray(temp.data$ch.j) # seems to be increasing suspiciously...
marray(temp.data$ch.a)

marrayonly<-nimbleCode({
  
  # CAPTURE RECAPTURE #####
  
  #m-array, multinomial likelihood
  for(t in 1:(nyears-1)){
    marr.j[t,1:nyears]~dmulti(pr.j[t,1:nyears],R.j[t])
    #marr.a[t,1:nyears]~dmulti(pr.a[t,1:nyears],R.a[t])
  }
  
  # TODO
  # monitor pr.j and pr.a and see what the structure is
  
  # diagonal
  for(t in 1:(nyears-1)){ # 1:14, 1:14
    q[t]<-1-p[t]
    pr.j[t,t]<-phi.j[t]*p[t]
    pr.a[t,t]<-phi.a[t]*p[t]
  }
  #upper triangle
  for(t in 1:(nyears-2)){ #1:13, 2:14
    for(j in (t+1):(nyears-1)){
      pr.j[t,j]<-phi.j[t]*prod(phi.a[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j]<-prod(phi.a[t:j])*prod(q[t:(j-1)])*p[j]
    }
  }
  #lower triangle
  for (t in 2:(nyears-1)){ #2:14, 1:13
    for(j in 1:(t-1)){
      pr.j[t,j]<-0
      pr.a[t,j]<-0
    }
  }
  for(t in 1:(nyears-1)){
    pr.a[t,nyears]<-1-sum(pr.a[t,1:(nyears-1)])
    pr.j[t,nyears]<-1-sum(pr.j[t,1:(nyears-1)])
  }
  for(t in 1:(nyears)){
    phi.j[t]<-mean.phi[1]
    phi.a[t]<-mean.phi[2]
    p[t]<-mean.p
  }
  
  #priors for resight and adult or 1yo survival
  mean.phi[1]~dunif(0,1) #surv 1 year olds
  mean.phi[2]~dunif(0,1) #surv adults
  mean.p~dunif(0,1) #resight prob
  
})

#### DATA ####
dat1 <- list(#marr.a = marray(temp.data$ch.a), 
             marr.j=marray(temp.data$ch.j),
             R.j=rowSums(marray(temp.data$ch.j))#, 
             #R.a = rowSums(marray(temp.data$ch.a))
)


#### CONSTANTS ####

const1 <- list(nyears = 15)

#### INITIAL VALUES ####
#z.state <- state.data(popDat$ch)

inits1 <- list(
  mean.phi = c(comb$phi1, comb$phiad),#c(detect.h, detect.h),
  mean.p = detect[2]
)

#### PARAMETERS TO MONITOR ####
params1 <- c("mean.phi","mean.p")#,"N1","Nad","f","rho")#0.3764911

#### COMPILE CONFIGURE AND BUILD ####
Rmodel1 <- nimbleModel(code = marrayonly, constants = const1, data = dat1,
                       check = FALSE, calculate = FALSE, inits = inits1)
conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
#control = list(maxContractions = 1000))
Rmcmc1 <- buildMCMC(conf1)
Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)

#### RUN MCMC ####
out_marray <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1, thin=nt,
                  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

post_summ(out_marray, get_params(out_marray), Rhat = TRUE, neff = TRUE) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column()

# so these estimates are not right, but they are not right in a different way than the fullIPM...
# the adult survival looks ok though...
# TODO let's test with just the adult m array
# ok this one looks pretty good for both mean p (0.8) and adult survival (0.4)
# TODO then just the juvenile m array
# ok we're underestimating p for one thing
# phi1 looks ok actually, if uncertain. oh wait - not converging actually...nvm
