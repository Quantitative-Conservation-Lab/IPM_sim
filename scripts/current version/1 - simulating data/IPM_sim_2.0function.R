# n.years = n years to simulate, n.data.types = proportion observed for each dataset, 
# age.init = starting age structure, phi.1 = juv survival,  phi.ad = adult survival,
# f = fecundity
simPopTrajectory <- function(n.years, n.data.types, age.init,
                             phi.1, phi.ad, f){
  
  #######
  # HELPER FUNCTIONS -- transition matrices
  # states: [1yrolds, Adult, chicks, Dead]
  # fate of juveniles
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
  
  deadfn<-function(ind,time){
    indfates[4,time+1,ind]<-1
    indfates[1:3,time+1,ind]<-NA
    return(indfates[,(time+1),ind])
  }
  
  # did the chicks produced above survive to the next time step?
  chickfatefn<-function(ind, time){
    zsurv1<-rbinom(1,1,phi.1)
    indfates[1,time+1,ind]<-ifelse(zsurv1==1,1,NA) #do they survive to become 1yrolds
    indfates[4,time+1,ind]<-ifelse(zsurv1==0, 1,NA) #or die?
    indfates[2,time+1,ind]<-NA #not adult yet
    indfates[3,time+1,ind]<-NA
    return(indfates[,(time+1),ind])
  }
  ######
  
  #nmin1.years<-n.years-1 # review again - is this still being used?
  nplus1.years<-n.years+1
  
  #derive true fecundity from the parameters, for productivity use f
  #if(productivity==T){
  fec <- f # AEB - do not divide by 2 
  #}else{
  #  fec<-1/2*mean.clutch.size*phi.nest^max.nest.age
  #}
  
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
  if(el<0.90 || el>1.10){
    print("Lambda is either too low or too high, revisit parameters")
    print(el)
  }
  #eigen(lesmat) #for lambda, if we need to check
  no.animals<-sum(N) #number of animals ever in the system at anytime
  no.ani.max<-round(no.animals*2.5) #include more for simulation of offspring
  ####stable age distribution is: 
  #only want to look at the end distribution
  sad<-N[,nplus1.years]/sum(N[,nplus1.years]) # the proportion of the population that is in each age class
  #lambda from N, matches eigenvalue
  #nlam<-N[,nplus1.years]/N[,n.years]
  #note that this is without demographic stochasticity and 
  #when we simulate the fates the actual lambda will change based on that
  
  
  #simulate what happens to each individual:
  if(N[1]==0 || N[2]==0){print("noindividuals!!!!!!!!!")}
  #dont run if there arent any individuals given 
  
  #make an array for each individual, for each year, what state it is in
  # states: [1yrolds, Adult, chicks, Dead]
  # Oneyearolds are the 1st years, chicks are how many chicks per individual,
  #dead is when they die
  
  #simulate their fates 
  
  #set up array for individuals
  indfates<-array(dim=c(4, nplus1.years, no.ani.max))
  #add in the stable age distribution to start from, given our total starting population
  age1<-round(sum(age.init)*sad[1])
  age2<-round(sum(age.init)*sad[2])
  indfates[1,1,1:age1]<-1 #one year olds in year 1
  indfates[2,1,(age1+1):sum(age1+age2)]<-1 #adults in year 1
  
  #this commented out bit is for if we dont want to start at stable age distribution
  #indfates[1,1,1:age.init[1]]<-1
  #indfates[2,1,(age.init[1]+1):sum(age.init)]<-1
  
  inpop<-numeric(n.years) #track the population over time
  #inpop[1]<-sum(age.init) alternative way to do it
  inpop[1]<-sum(age1+age2) #stable age distribution way to do it
  
  # simulate ind. trajectories with dem stochasticity  
  chickst<-numeric(n.years) # tracking number of chicks produced each year
  tempstep<-matrix(nrow=no.ani.max, ncol=(n.years+1))
  for(t in 1:n.years){
    for(i in 1:inpop[t]){
      tempstep[i,t]<-which(indfates[c(1,2,4),t,i]==1)
      if(tempstep[i,t]==1){
        indfates[,(t:(t+1)),i]<-oneyr_fatefn(ind = i, time = t)#[,(t:(t+1))]
      }else if(tempstep[i,t]==2){
        indfates[,(t:(t+1)),i]<-adfatefn(ind = i, time = t)#[,(t:(t+1))]
      } else if(tempstep[i,t]==3){
        indfates[,((t+1)),i]<-deadfn(ind = i, time = t)#[,(t:(t+1))]
      }
    } #i loop over population size at time t
    #see how many chicks were produced and add them to the population to track
    chickst[t]<-sum(indfates[3,t,], na.rm=T)
    inpop[t+1]<-inpop[t]+chickst[t]
    if(inpop[t+1]>inpop[t]){
      indfates[3,t,(inpop[t]+1):inpop[t+1]]<-1
      for(i in (inpop[t]+1):inpop[t+1]){
        indfates[,((t+1)),i]<-chickfatefn(ind = i, time = t)#[,(t:(t+1))]
      }
    }else{}
  } 
  
  #check observed lambda
  ad<-j<-numeric(n.years)
  ch<-dead<-numeric(n.years)
  for(i in 1:n.years){
    ad[i]<-sum(indfates[2,i,], na.rm=T)
    #ch[i]<-sum(indfates[3,i,], na.rm=T)
    j[i]<-sum(indfates[1,i,], na.rm=T)
    #dead[i]<-sum(indfates[4,i,], na.rm=T)
  }
  
  adj<-numeric(length(ad))
  for(i in 2:length(adj)){
    adj[i]<-sum(ad[i]+j[i])/sum(ad[i-1]+j[i-1])
  }
  mean(adj[2:10])
  #close *enough* to actual lambda, again demographic stochasticity is the cause!
  
  #HAS:
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
  
simData <- function(indfates, n.years, n.data.types, 
                    ADonly, p.1, p.ad,
                    BinMod, n.sam, p.count, sig, productivity, p.prod){
  
  ###################For output of data ##############
  # Three independent samples from the individuals we just simulated
  Ntot <- dim(indfates)[3]
  xt1 <- matrix(data = seq(1:Ntot), ncol = 1) # individual ID
  resamp1 <- resamp2 <- resamp3 <- numeric()
  nds<-numeric(length(n.data.types))
  if(sum(n.data.types)>1){ # is n.data.types numbers
    nds<-n.data.types # number of individuals we want in each independent dataset
  }else{ # or proportions
    nds<-Ntot*n.data.types
  }
  # Sample 1
  resamp1 <- sample(xt1, nds[1], replace = F)
  # Sample 2
  #resamp2 <- sample(xt1[-resamp1,], nds[2], replace = F)
  resamp2 <- xt1[,1] # observe all the individuals in the population
  
  resamp3<-sample(xt1[c(-resamp1,-resamp2),], nds[3], replace = F)
  
  # Individuals selected for capture-recapture, survey and reproductive success data
  IND_MR <- indfates[,,resamp1]
  IND_Count <- indfates[,,resamp2]
  IND_Nest <- indfates[,,resamp3]
  
  ######################################################
  # Create Mark-Resight data
  ######################################################
  #HAS: Do we only want to mark adults? and never see 1year olds
  #Or do we want to be able to mark anyone we see? then we can get at 1yearoldSurv
  #We did just adults, but I added in an option for 
  #also marking any 1 year olds we see
  
  
  #Marking adults and 1 year olds, which will have the same
  #for adults and 1yr olds only
  #AD_only=T
  
  #for adults, 1 year olds and chicks
  #AD_only=F
  
  if(ADonly==T){
    #using the individual population array for MR data
    mr_classes<-dim(IND_MR)[1] - 3 #since we wont see Dead and the reproduction doesnt matter
    #so we have 2 classes, we care about in marking: 1year olds, and adults
    ind_mr<-IND_MR[c(1,2),1:n.years,]
    #ind_mr[1,,]<-NA #if not banding 1yearold/s, then remove
    rm<-numeric(dim(ind_mr)[3]) # remove dead individuals
    for(i in 1:dim(ind_mr)[3]){
      if(length(which(!is.na(ind_mr[1:2,,i])))==0){
        rm[i]<-1
      }else{
        rm[i]<-0
      }
    }
    if(sum(rm>0)){
      ind_mr<-ind_mr[,,-which(rm==1)]
    }else{}
    
    
    age<-first<-last<-numeric() # age, first and last encounters
    mr_t<-dim(ind_mr)[2]
    mr_ind<-dim(ind_mr)[3] 
    for(i in 1:mr_ind){
      g <- which(!is.na(ind_mr[1:2,,i]), arr.ind = TRUE) # ind that were seen
      age[i] <- g[1,1] # age at marking
      first[i] <- g[1,2] # first time seen
      h <- which(ind_mr[1:2,,i]==1, arr.ind = TRUE) # last time seen
      last[i] <- max(h[,2])
    }  
    #remove those that we never banded as a chick
    ch.true<-matrix(0,ncol=mr_t, nrow=mr_ind)
    for(i in 1:mr_ind){
      ch.true[i,first[i]:last[i]]<-1
    }  
    
    #detection of true marked individuals:
    #since inital marking is constant prob
    #and resight is constant prob
    #phi.1 #probability you were marked as 1 year old
    #p.ad #probability you were marked as adult
    in.mark<-ch<-matrix(0,nrow=mr_ind, ncol=mr_t)
    for(i in 1:mr_ind){
      # if(age[i]==1){
      #   in.mark[i,first[i]]<-rbinom(1,1,p.1*ch.true[i,first[i]])
      # } else{
      in.mark[i,first[i]]<-rbinom(1,1,p.ad*ch.true[i,first[i]]) 
      #}
      if(first[i]==mr_t) next
      for(t in (first[i]+1):last[i]){
        #if(age[i]==1){
        #  in.mark[i,t]<-rbinom(1,1,p.1*ch.true[i,t])
        #}
        in.mark[i,t]<-rbinom(1,1,p.ad*ch.true[i,t])
      }
      ch[i,]<-in.mark[i,]
    }
    
    #code to track the ages, 1 is for chicks, 2 1years, 3 adults (code not currently in use)
    add_age_chtrue<-age_ch<-matrix(0,nrow=mr_ind, ncol=mr_t)
    for(i in 1:mr_ind){
      for(x in (first[i]):last[i]){
        add_age_chtrue[i,x]<-as.numeric(which(!is.na(ind_mr[,x,i])))
        age_ch[i,x]<-add_age_chtrue[i,x]*ch[i,x]
      }
    }
    rm_2<-numeric(mr_ind) # remove those we never saw
    for(i in 1:mr_ind){
      if(sum(ch[i,])==0){
        rm_2[i]<-1
      }else {
        rm_2[i]<-0
      }
      #rm_2[i]<-sum(ch[i,])
    }
    if(sum(rm_2)>0){
      ch<-ch[-(which(rm_2==1)),]
      age_ch<-age_ch[-(which(rm_2==1)),]
      add_age_chtrue[-(which(rm_2==1)),]
    } else {} 
    firstobs<-lastobs<-numeric(length(ch[,1]))
    for(i in 1:length(ch[,1])){
      firstobs[i]<-min(which(ch[i,]==1))
      lastobs[i]<-max(which(ch[i,]==1))
    }
    
  }else{ 
    #if we mark chicks, 1 years, adults
    #!!
    #using the individual population array for MR data
    #mr_classes<-dim(IND_MR)[1] - 3 #since we wont see Dead and the reproduction doesnt matter
    #so we have 2 classes, we care about in marking: 1year olds, and adults
    ind_mr<-IND_MR[c(1,2,3),1:n.years,]
    #ind_mr[1,,]<-NA #if not banding 1yearold/s, then remove
    rm<-numeric(dim(ind_mr)[3])
    for(i in 1:dim(ind_mr)[3]){
      if(length(which(!is.na(ind_mr[1:3,,i])))==0){
        rm[i]<-1
      }else{
        rm[i]<-0
      }
    }
    if(sum(rm>0)){
      ind_mr<-ind_mr[,,-which(rm==1)]
    }else{}
    
    
    age<-first<-last<-numeric()
    mr_t<-dim(ind_mr)[2]
    mr_ind<-dim(ind_mr)[3] 
    for(i in 1:mr_ind){
      g <- which(!is.na(ind_mr[1:3,,i]), arr.ind = TRUE)
      age[i] <- g[1,1] #at marking
      first[i] <- g[1,2]
      h <- which(ind_mr[1:3,,i]==1, arr.ind = TRUE)
      last[i] <- max(h[,2])
    }  
    ch.true<-matrix(0,ncol=mr_t, nrow=mr_ind)
    for(i in 1:mr_ind){
      ch.true[i,first[i]:last[i]]<-1
    }  
    
    #since inital marking is constant prob
    #and resight is constant prob
    #phi.1 #probability you were marked as 1 year old
    #p.ad #probability you were marked as adult
    in.mark<-ch<-matrix(0,nrow=mr_ind, ncol=mr_t)
    for(i in 1:mr_ind){
      if(age[i]==1){
        in.mark[i,first[i]]<-rbinom(1,1,p.1*ch.true[i,first[i]])
      } else{
        in.mark[i,first[i]]<-rbinom(1,1,p.ad*ch.true[i,first[i]])
      }
      if(first[i]==mr_t) next
      for(t in (first[i]+1):last[i]){
        if(age[i]==1){
          in.mark[i,t]<-rbinom(1,1,p.1*ch.true[i,t])
        }
        in.mark[i,t]<-rbinom(1,1,p.ad*ch.true[i,t])
      }
      ch[i,]<-in.mark[i,]
    }
    
    #code to track the ages, 1 is for chicks, 2 1years, 3 adults
    add_age_chtrue<-age_ch<-matrix(0,nrow=mr_ind, ncol=mr_t)
    for(i in 1:mr_ind){
      for(x in (first[i]):last[i]){
        add_age_chtrue[i,x]<-as.numeric(which(!is.na(ind_mr[,x,i])))[1]
        age_ch[i,x]<-add_age_chtrue[i,x]*ch[i,x]
      }
    }#here 3 is chick, 1 is 1yearold, 2 is adult
    rm_2<-numeric(mr_ind)
    for(i in 1:mr_ind){
      if(sum(ch[i,])==0){
        rm_2[i]<-1
      }else {
        rm_2[i]<-0
      }
      #rm_2[i]<-sum(ch[i,])
    }
    if(sum(rm_2)>0){
      ch<-ch[-(which(rm_2==1)),]
      age_ch<-age_ch[-(which(rm_2==1)),]
      add_age_chtrue[-(which(rm_2==1)),]
    } else {}
    firstobs<-lastobs<-numeric(length(ch[,1]))
    for(i in 1:length(ch[,1])){
      firstobs[i]<-min(which(ch[i,]==1))
      lastobs[i]<-max(which(ch[i,]==1))
    }
  }
  
  
  ######################################################
  # Create population survey data
  ######################################################
  #n.sam is the number of times that the population was sampled in a year
  #Need to choose in the function argument which model to use, BinMod=T if binom, F if normal
  #HAS - should output the true count data for comparison
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
  
  ######################################################
  # Create reproductive success data
  ######################################################
  #if productivity is T, then this nest success model
  #Abby can add in if we want to include it - AEB: let's omit for now
  #if F, then it is the below model
  #IND_Nest <- indfates[,,resamp3]
  TRUE_nestlings<-numeric(n.years)
  R_true<-R_obs <- numeric(n.years) # number of pairs whose productivity was observed 
  OBS_nestlings <- numeric(n.years)  # total number of nestlings recorded in a year
  for(t in 1:n.years){
    #true number of nestlings
    TRUE_nestlings[t]<-sum(IND_Nest[3,t,], na.rm=T)
    #true number of reproducing females
    R_true[t]<-length(which(IND_Nest[3,t,]>0)) 
    
    #observation, who dont we see and who do we see
    #Assumes we can count the number of chicks perfectly if we see the nest
    for(i in 1:length(IND_Nest[1,1,])){
      if(!is.na(IND_Nest[3,t,i])){
        obstemp<-rbinom(1,1,p.prod)
        if(obstemp==1){
          R_obs[t]<-R_obs[t]+1
          OBS_nestlings[t]<-IND_Nest[3,t,i]+OBS_nestlings[t]
        }
      }
    }
    
  }
  
  
  return(list(ch=ch, SUR=SUR, age_ch=age_ch, firstobs=firstobs, lastobs=lastobs,
              n.sam = n.sam, R_obs=R_obs, OBS_nestlings=OBS_nestlings))
  
}


# ARCHIVE #######

#IPM simulation code

#INPUT FOR DATA SIMULATION FUNCTION

# TODO
# could do a T/F for stable age dist or not

# #1.) YEARS: number of years to simulate
# n.years<-15 # first 5 are burn-in (population deviates from stable age dist), last 10 are fit
# 
# #2.) OBSERVED DATA OUTPUT: a vector specified as c(Mark-resight, Count, Productivity)
# # Can be either the proportion of the simulated individuals that are observed for 
# # each data type.
# n.data.types<-c(0.25,0.3,0.25)
# #OR the number of individuals in the population that we observe
# #Note that this is the maximum and especially for MR data, if the individual died, then
# # there will be no observation of them, so observing 100 individuals is not guaranteed
# #same goes for count and nest data
# n.data.types<-c(100,100,100)
# 
# #3.) TRUE STARTING POPULATION/INITIAL AGE DISTRIBUTION
# #Can specify it as the number of 1year olds and adults:
# age.init<-c(500,500)
# 
# #The stable age distribution from the leslie matrix will be used to initialize 
# #the initial starting population, so it will always start with 1000 individuals, but
# # the distribution will depend on the stable age distribution
# 
# 
# #4.) Constant Parameter Values to simulate population from
# 
#   #a.) 1year old survival
# phi.1<-0.5
#   #b.) adult survival
# phi.ad<-0.7
#   #c.) fecundity
#   #where fecundity is the same for 1year olds and adults
#   #can be specified as a single rate with 'productivity=T' or as in the nest success model
# productivity=T
# f<-1.4 
# 
# #5.) PARAMETERS FOR DETECTION AND MODEL SELECTION
# 
# # TODO
# # AEB - this section will change
# 
# #a.) Mark Resight model
# #If we only mark adults and resight adults, then 
# ADonly = T
# p.ad <- 0.7
# #If we mark 1yearolds and adults
# ADonly=F
# p.1<-0.5
# 
# #b.) Count model
# #define the number of times we go and count
# n.sam<-3
# #if we want to model the observation process for the count data as binomial, 
# # where counts~Binomial(truecount,p.count) then:
# BinMod=T
# p.count<-0.55
# #If we want to model it as normal, where counts~normal(truecounts, sig)
# BinMod=F
# sig<-2
# 
# #c) Productivity/nest success
# #we defined which we want to use above, with productivity=T/F
# #If we want a Poisson process to model the productivity, then add in observation
# p.prod<-0.4

IPMSimFunction<-function(n.years, n.data.types, age.init, phi.1, phi.ad, f, max.nest.age,
                         mean.clutch.size, phi.nest, ADonly, p.1, p.ad,
                         BinMod, n.sam, p.count, sig, productivity, p.prod){
  
  nmin1.years<-n.years-1
  nplus1.years<-n.years+1
  
  #derive true fecundity from the parameters, for productivity use f
  if(productivity==T){
    fec<-f/2 
  }else{
    fec<-1/2*mean.clutch.size*phi.nest^max.nest.age
  }

#Leslie matrix for population size
N<-matrix(nrow=length(age.init), ncol=nplus1.years)
N[,1]<-age.init
lesmat<-matrix(nrow=2, ncol=2)
lesmat[1,]<-c(phi.1*fec, phi.1*fec)
lesmat[2,]<-c(phi.ad, phi.ad)
for(t in 1:n.years){
  N[,(t+1)]<-lesmat%*%N[,t]
}
el<-eigen(lesmat)$values[1]
if(el<0.95 || el>1.05){
  print("Lambda is either too low or too high, revisit parameters")
  print(el)
}
#eigen(lesmat) #for lambda, if we need to check
no.animals<-sum(N) #number of animals ever in the system at anytime
no.ani.max<-round(no.animals*2.5) #include more for simulation of offspring
####stable age distribution is:
sad<-N[,nplus1.years]/sum(N[,nplus1.years]) # the proportion of the population that is in each age class
#lambda from N, matches eigenvalue
#nlam<-N[,nplus1.years]/N[,n.years]
#note that this is without demographic stochasticity and 
#when we simulate the fates the actual lambda will change based on that


#simulate what happens to each individual:
if(N[1]==0 || N[2]==0){print("noindividuals!!!!!!!!!")}
#dont run if there arent any individuals given 

#make an array for each individual, for each year, what state it is in
# states: [1yrolds, Adult, chicks, Dead]
# Oneyearolds are the 1st years, chicks are how many chicks per individual,
#dead is when they die

#simulate their fates 

#set up array for individuals
indfates<-array(dim=c(4, nplus1.years, no.ani.max))
#add in the stable age distribution to start from, given our total starting population
age1<-round(sum(age.init)*sad[1])
age2<-round(sum(age.init)*sad[2])
indfates[1,1,1:age1]<-1 #one year olds in year 1
indfates[2,1,(age1+1):sum(age1+age2)]<-1 #adults in year 1

#this commented out bit is for if we dont want to start at stable age distribution
#indfates[1,1,1:age.init[1]]<-1
#indfates[2,1,(age.init[1]+1):sum(age.init)]<-1

inpop<-numeric(n.years) #track the population over time
#inpop[1]<-sum(age.init) alternative way to do it
inpop[1]<-sum(age1+age2) #stable age distribution way to do it

chickst<-numeric(n.years)
tempstep<-matrix(nrow=no.ani.max, ncol=(n.years+1))
for(t in 1:n.years){
  for(i in 1:inpop[t]){
    tempstep[i,t]<-which(indfates[c(1,2,4),t,i]==1)
    if(tempstep[i,t]==1){
      indfates[,(t:(t+1)),i]<-oneyr_fatefn(ind=i, time=t)#[,(t:(t+1))]
    }else if(tempstep[i,t]==2){
      indfates[,(t:(t+1)),i]<-adfatefn(ind=i, time=t)#[,(t:(t+1))]
    } else if(tempstep[i,t]==3){
      indfates[,((t+1)),i]<-deadfn(ind=i, time=t)#[,(t:(t+1))]
    }
  } #i loop over population size at time t
  #see how many chicks were produced and add them to the population to track
  chickst[t]<-sum(indfates[3,t,], na.rm=T)
  inpop[t+1]<-inpop[t]+chickst[t]
  if(inpop[t+1]>inpop[t]){
  indfates[3,t,(inpop[t]+1):inpop[t+1]]<-1
  for(i in (inpop[t]+1):inpop[t+1]){
    indfates[,((t+1)),i]<-chickfatefn(ind=i, time=t)#[,(t:(t+1))]
  }
  }else{}
} 

#check observed lambda
ad<-j<-numeric(n.years)
ch<-dead<-numeric(n.years)
for(i in 1:n.years){
  ad[i]<-sum(indfates[2,i,], na.rm=T)
  #ch[i]<-sum(indfates[3,i,], na.rm=T)
  j[i]<-sum(indfates[1,i,], na.rm=T)
  #dead[i]<-sum(indfates[4,i,], na.rm=T)
}

adj<-numeric(length(ad))
for(i in 2:length(adj)){
  adj[i]<-sum(ad[i]+j[i])/sum(ad[i-1]+j[i-1])
}
mean(adj[2:10])
#close *enough* to actual lambda, again demographic stochasticity is the cause!


###################For output of data ##############
# Three independent samples from the individuals we just simulated
Ntot <- dim(indfates)[3]
xt1 <- matrix(data = seq(1:Ntot), ncol = 1)
resamp1 <- resamp2 <- resamp3 <- numeric()
nds<-numeric(length(n.data.types))
if(sum(n.data.types)>1){ # is n.data.types numbers
  nds<-n.data.types
}else{ # or proportions
  nds<-Ntot*n.data.types
}
# Sample 1
resamp1 <- sample(xt1, nds[1], replace = F)
# Sample 2
resamp2 <- sample(xt1[-resamp1,], nds[2], replace = F)

resamp3<-sample(xt1[c(-resamp1,-resamp2),], nds[3], replace = F)

# Individuals selected for capture-recapture, survey and reproductive success data
IND_MR <- indfates[,,resamp1]
IND_Count <- indfates[,,resamp2]
IND_Nest <- indfates[,,resamp3]

######################################################
# Create Mark-Resight data
######################################################
#HAS: Do we only want to mark adults? and never see 1year olds
#Or do we want to be able to mark anyone we see? then we can get at 1yearoldSurv
#We did just adults, but I added in an option for 
#also marking any 1 year oldss we see


#Marking adults and 1 year olds, which will have the same
#for adults and 1yr olds only
#AD_only=T

#for adults, 1 year olds and chicks
#AD_only=F

if(ADonly==T){
  #using the individual population array for MR data
  mr_classes<-dim(IND_MR)[1] - 3 #since we wont see Dead and the reproduction doesnt matter
  #so we have 2 classes, we care about in marking: 1year olds, and adults
  ind_mr<-IND_MR[c(1,2),1:n.years,]
  #ind_mr[1,,]<-NA #if not banding 1yearold/s, then remove
  rm<-numeric(dim(ind_mr)[3])
  for(i in 1:dim(ind_mr)[3]){
    if(length(which(!is.na(ind_mr[1:2,,i])))==0){
      rm[i]<-1
    }else{
      rm[i]<-0
    }
  }
  if(sum(rm>0)){
    ind_mr<-ind_mr[,,-which(rm==1)]
  }else{}
  
  
  age<-first<-last<-numeric()
  mr_t<-dim(ind_mr)[2]
  mr_ind<-dim(ind_mr)[3] 
  for(i in 1:mr_ind){
    g <- which(!is.na(ind_mr[1:2,,i]), arr.ind = TRUE)
    age[i] <- g[1,1] #at marking
    first[i] <- g[1,2]
    h <- which(ind_mr[1:2,,i]==1, arr.ind = TRUE)
    last[i] <- max(h[,2])
  }  
  #remove those that we never banded as a chick
  ch.true<-matrix(0,ncol=mr_t, nrow=mr_ind)
  for(i in 1:mr_ind){
    ch.true[i,first[i]:last[i]]<-1
  }  
  
  #since inital marking is constant prob
  #and resight is constant prob
  #phi.1 #probability you were marked as 1 year old
  #p.ad #probability you were marked as adult
  in.mark<-ch<-matrix(0,nrow=mr_ind, ncol=mr_t)
  for(i in 1:mr_ind){
    # if(age[i]==1){
    #   in.mark[i,first[i]]<-rbinom(1,1,p.1*ch.true[i,first[i]])
    # } else{
      in.mark[i,first[i]]<-rbinom(1,1,p.ad*ch.true[i,first[i]])
    #}
    if(first[i]==mr_t) next
    for(t in (first[i]+1):last[i]){
      #if(age[i]==1){
      #  in.mark[i,t]<-rbinom(1,1,p.1*ch.true[i,t])
      #}
      in.mark[i,t]<-rbinom(1,1,p.ad*ch.true[i,t])
    }
    ch[i,]<-in.mark[i,]
  }
  
  #code to track the ages, 1 is for chicks, 2 1years, 3 adults
  add_age_chtrue<-age_ch<-matrix(0,nrow=mr_ind, ncol=mr_t)
  for(i in 1:mr_ind){
    for(x in (first[i]):last[i]){
      add_age_chtrue[i,x]<-as.numeric(which(!is.na(ind_mr[,x,i])))
      age_ch[i,x]<-add_age_chtrue[i,x]*ch[i,x]
    }
  }
  rm_2<-numeric(mr_ind)
  for(i in 1:mr_ind){
    if(sum(ch[i,])==0){
      rm_2[i]<-1
    }else {
      rm_2[i]<-0
    }
    #rm_2[i]<-sum(ch[i,])
  }
  if(sum(rm_2)>0){
    ch<-ch[-(which(rm_2==1)),]
    age_ch<-age_ch[-(which(rm_2==1)),]
    add_age_chtrue[-(which(rm_2==1)),]
  } else {}
  firstobs<-lastobs<-numeric(length(ch[,1]))
  for(i in 1:length(ch[,1])){
    firstobs[i]<-min(which(ch[i,]==1))
    lastobs[i]<-max(which(ch[i,]==1))
  }
  
}else{
#if we mark chicks, 1 years, adults
  #!!
#using the individual population array for MR data
#mr_classes<-dim(IND_MR)[1] - 3 #since we wont see Dead and the reproduction doesnt matter
#so we have 2 classes, we care about in marking: 1year olds, and adults
ind_mr<-IND_MR[c(1,2,3),1:n.years,]
#ind_mr[1,,]<-NA #if not banding 1yearold/s, then remove
rm<-numeric(dim(ind_mr)[3])
for(i in 1:dim(ind_mr)[3]){
  if(length(which(!is.na(ind_mr[1:3,,i])))==0){
    rm[i]<-1
  }else{
    rm[i]<-0
  }
}
if(sum(rm>0)){
  ind_mr<-ind_mr[,,-which(rm==1)]
}else{}


age<-first<-last<-numeric()
mr_t<-dim(ind_mr)[2]
mr_ind<-dim(ind_mr)[3] 
for(i in 1:mr_ind){
  g <- which(!is.na(ind_mr[1:3,,i]), arr.ind = TRUE)
  age[i] <- g[1,1] #at marking
  first[i] <- g[1,2]
  h <- which(ind_mr[1:3,,i]==1, arr.ind = TRUE)
  last[i] <- max(h[,2])
}  
ch.true<-matrix(0,ncol=mr_t, nrow=mr_ind)
for(i in 1:mr_ind){
  ch.true[i,first[i]:last[i]]<-1
}  

#since inital marking is constant prob
#and resight is constant prob
#phi.1 #probability you were marked as 1 year old
#p.ad #probability you were marked as adult
in.mark<-ch<-matrix(0,nrow=mr_ind, ncol=mr_t)
for(i in 1:mr_ind){
  if(age[i]==1){
    in.mark[i,first[i]]<-rbinom(1,1,p.1*ch.true[i,first[i]])
  } else{
    in.mark[i,first[i]]<-rbinom(1,1,p.ad*ch.true[i,first[i]])
  }
  if(first[i]==mr_t) next
  for(t in (first[i]+1):last[i]){
    if(age[i]==1){
      in.mark[i,t]<-rbinom(1,1,p.1*ch.true[i,t])
    }
    in.mark[i,t]<-rbinom(1,1,p.ad*ch.true[i,t])
  }
  ch[i,]<-in.mark[i,]
}

#code to track the ages, 1 is for chicks, 2 1years, 3 adults
add_age_chtrue<-age_ch<-matrix(0,nrow=mr_ind, ncol=mr_t)
for(i in 1:mr_ind){
  for(x in (first[i]):last[i]){
    add_age_chtrue[i,x]<-as.numeric(which(!is.na(ind_mr[,x,i])))[1]
    age_ch[i,x]<-add_age_chtrue[i,x]*ch[i,x]
  }
}#here 3 is chick, 1 is 1yearold, 2 is adult
rm_2<-numeric(mr_ind)
for(i in 1:mr_ind){
  if(sum(ch[i,])==0){
    rm_2[i]<-1
  }else {
    rm_2[i]<-0
  }
  #rm_2[i]<-sum(ch[i,])
}
if(sum(rm_2)>0){
  ch<-ch[-(which(rm_2==1)),]
  age_ch<-age_ch[-(which(rm_2==1)),]
  add_age_chtrue[-(which(rm_2==1)),]
} else {}
firstobs<-lastobs<-numeric(length(ch[,1]))
for(i in 1:length(ch[,1])){
  firstobs[i]<-min(which(ch[i,]==1))
  lastobs[i]<-max(which(ch[i,]==1))
}
}


######################################################
# Create population survey data
######################################################
#n.sam is the number of times that the population was sampled in a year
#Need to choose in the function argument which model to use, BinMod=T if binom, F if normal
#HAS - should output the true count data for comparison
TRUE_Count<-matrix(nrow=2, ncol=n.years) #first row, number of YOY; second row is adults
SUR<-matrix(nrow=n.sam, ncol=n.years)#matrix for surveys by survey and year
for(u in 1:n.years){
  #Some discussion about should we include ?1yearolds
  #if so, remove the next line
  TRUE_Count[1,u]<-sum(IND_Count[1,u,], na.rm = T)
  TRUE_Count[2,u]<-sum(IND_Count[2,u,], na.rm = T)
  if(BinMod==T){
    SUR[,u]<-rbinom(n.sam, sum(TRUE_Count[,u]), p.count)
  }else{
    SUR[,u]~rnorm(n.sam, TRUE_count[,u], sig)
  }
}

######################################################
# Create reproductive success data
######################################################
#if producivity is T, then this nest success model
#Abby can add in if we want to include it - AEB: let's omit for now
#if F, then it is the below model
#IND_Nest <- indfates[,,resamp3]
TRUE_nestlings<-numeric(n.years)
R_true<-R_obs <- numeric(n.years) # number of pairs whose productivity was observed 
OBS_nestlings <- numeric(n.years)  # total number of nestlings recoded in a year
for(t in 1:n.years){
  #true number of nestlings
  TRUE_nestlings[t]<-sum(IND_Nest[3,t,], na.rm=T)
  #true number of reproducing females
  R_true[t]<-length(which(IND_Nest[3,t,]>0)) 
  
  #observation, who dont we see and who do we see
  #Assumes we can count the number of chicks perfectly if we see the nest
  for(i in 1:length(IND_Nest[1,1,])){
    if(!is.na(IND_Nest[3,t,i])){
      obstemp<-rbinom(1,1,p.prod)
      if(obstemp==1){
        R_obs[t]<-R_obs[t]+1
        OBS_nestlings[t]<-IND_Nest[3,t,i]+OBS_nestlings[t]
      }
    }
  }
  
}
  

  return(list(ch=ch, SUR=SUR, age_ch=age_ch, firstobs=firstobs, lastobs=lastobs,
              n.sam, fec=fec, N=N, R_obs=R_obs, OBS_nestlings=OBS_nestlings))

}

# df<-IPMSimFunction(n.years=10, n.data.types=c(0.25,0.25,0.25),
#                    age.init=c(150,150), phi.1=0.3, phi.ad=0.3, f=0.5, max.nest.age=NA,
#                    mean.clutch.size=NA, phi.nest=NA, ADonly=T,p.1=NA,p.ad=0.8,
#                    BinMod=T,n.sam=3,p.count=0.55,sig=NA,
#                    productivity=T,p.prod=0.65)
# n.years=10
# n.data.types=c(0.25,0.25,0.25) 
# age.init=c(50,50)
# phi.1=0.2
# phi.ad=0.2
# f=0.5
# max.nest.age=NA
# mean.clutch.size=NA
# phi.nest=NA
# ADonly=F
# p.1=0.3
# p.ad=0.2
# BinMod=F
# n.sam=3
# p.count=0.2
# sig=0.3
# productivity=T
# p.prod=0.65
