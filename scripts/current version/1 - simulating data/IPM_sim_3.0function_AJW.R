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

# Pick values for the function arguments
nyears <- 20                                 # Number of years
phi <- c(0.3, 0.55)                     # Age specific survival probabilities (juv, adult)
f <- c(2.5, 2.5)    # Age-specific productivity (1y, older)
Ni <- c(50, 50)                         # Initial pop. size for each age class (1y, older)

# Apply the function and produce data overview
# set.seed(111167)                        # To initialize the RNGs at the same place
pop <- simPop(Ni=Ni, phi=phi, f=f, nYears=nyears)
# str(pop)

pop1 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
pop2 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
pop3 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)

pop1$totBreeders

# Pick a value of the observation error (SD) for the population survey
# sigma <- 10
pDetect <- 0.8

# Create the population survey data and produce data overview
# count <- simCountNorm(N=pop1$totB, sigma=sigma)
Ncount <- simCountBin(N=pop1$totB, pDetect)
str(count)

# Pick values for capture and recapture probabilities
cap <- 0.4                              # Initial capture probability (same for juv. and adults)
recap <- 0.6                            # Recapture probability

# Create the capture histories and produce data overview
ch <- simCapHist(state=pop2$state, cap=cap, recap=recap, maxAge=2, verbose = F)

# Create m-arrays
marr <- marrayAge(ch$ch, ch$age)

# Pick probability to find a brood
pprod <- 0.3

# Create productivity data; check females.only
pro <- simProd(reprod=pop3$reprod, pInclude=pprod, females.only = T)
str(pro)















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
  # if(el<0.90 || el>1.10){
  #   print("Lambda is either too low or too high, revisit parameters")
  #   print(el)
  # }
  #eigen(lesmat) #for lambda, if we need to check
  
  no.animals<-sum(N) #number of animals ever in the system at anytime
  no.ani.max<-round(no.animals*5) #include more for simulation of offspring
  ####stable age distribution is:
  #only want to look at the end distribution
  sad<-N[,nplus1.years]/sum(N[,nplus1.years]) # the proportion of the population that is in each age class
  #lambda from N, matches eigenvalue
  #nlam<-N[,nplus1.years]/N[,n.years]
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
  #first dimension is age groups; j, adult, chicks, dead
  indfates<-array(dim=c(4, nplus1.years, no.ani.max))
  #add in the stable age distribution to start from, given our total starting population
  age1<-round(sum(age.init)*sad[1])
  age2<-round(sum(age.init)*sad[2])
  indfates[1,1,1:age1]<-1 #number of one year olds in year 1
  indfates[2,1,(age1+1):sum(age1+age2)]<-1 #number of adults in year 1

  #this commented out bit is for if we dont want to start at stable age distribution
  #indfates[1,1,1:age.init[1]]<-1
  #indfates[2,1,(age.init[1]+1):sum(age.init)]<-1

  inpop<-numeric(n.years) #track the population over time
  #inpop[1]<-sum(age.init) alternative way to do it
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
    #question for Hannah.... this counts 'if you've ever been alive'? 
    #if 50 in first year, add chicks, then.... does this get rid of chicks that die? 
    chickst[t]<-sum(indfates[3,t,], na.rm=T)
    inpop[t+1]<-inpop[t]+chickst[t]
    #question for hannah; what's this 'if' statement accomplishing?
    if(inpop[t+1]>inpop[t]){
      indfates[3,t,(inpop[t]+1):inpop[t+1]]<-1
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
    #ch[i]<-sum(indfates[3,i,], na.rm=T)
    j[i]<-sum(indfates[1,i,], na.rm=T)
    #dead[i]<-sum(indfates[4,i,], na.rm=T)
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

  #for adults, 1 year olds and chicks
  if(!is.na(p.1) & !is.na(p.ad)) {
    
    # #if(ADonly==T){
      #using the population array for MR data
      mr_classes.a<-dim(IND_MR)[1] - 3 #since we wont see Dead and the reproduction doesnt matter
      #so we have 2 classes, we care about in marking: 1year olds, and adults
      
      ind_mr.a<-IND_MR[c(1,2),1:n.years,] 
      
      #ind_mr[1,,]<-NA #if not banding 1yearold/s, then remove
      rm<-numeric(dim(ind_mr.a)[3]) #vector to remove dead individuals
      for(i in 1:dim(ind_mr.a)[3]){
        if(length(which(!is.na(ind_mr.a[1:2,,i])))==0){
          rm[i]<-1
        }else{
          rm[i]<-0
        }
      }
      if(sum(rm>0)){
        ind_mr.a<-ind_mr.a[,,-which(rm==1)]
      }else{}
      
      
      age.a<-first.a<-last.a<-numeric() # age, first and last encounters
      mr_t.a<-dim(ind_mr.a)[2]
      mr_ind.a<-dim(ind_mr.a)[3]
      for(i in 1:mr_ind.a){
        g <- which(!is.na(ind_mr.a[1:2,,i]), arr.ind = TRUE) # ind that were seen
        age.a[i] <- g[1,1] # age at marking
        first.a[i] <- g[1,2] # first time seen
        h <- which(ind_mr.a[1:2,,i]==1, arr.ind = TRUE) # last time seen
        last.a[i] <- max(h[,2])
      }
      #remove those that we never banded as a chick
      ch.true.a<-matrix(0,ncol=mr_t.a, nrow=mr_ind.a)
      for(i in 1:mr_ind.a){
        ch.true.a[i,first.a[i]:last.a[i]]<-1
      }
      
      #detection of true marked individuals:
      #since inital marking is constant prob
      #and resight is constant prob
      #p.1 #probability you were marked as 1 year old
      #p.ad #probability you were marked as adult
      in.mark.a<-ch.a<-matrix(0,nrow=mr_ind.a, ncol=mr_t.a)
      for(i in 1:mr_ind.a){
        # if(age[i]==1){
        #   in.mark[i,first[i]]<-rbinom(1,1,p.1*ch.true[i,first[i]])
        # } else{
        in.mark.a[i,first.a[i]]<-rbinom(1,1,p.ad*ch.true.a[i,first.a[i]])
        #}
        if(first.a[i]==mr_t.a) next
        for(t in (first.a[i]+1):last.a[i]){
          #if(age[i]==1){
          #  in.mark[i,t]<-rbinom(1,1,p.1*ch.true[i,t])
          #}
          in.mark.a[i,t]<-rbinom(1,1,p.ad*ch.true.a[i,t])
        }
        ch.a[i,]<-in.mark.a[i,]
      }
      
      #code to track the ages, 1 is for chicks, 2 1years, 3 adults (code not currently in use)
      add_age_chtrue.a<-age_ch.a<-matrix(0,nrow=mr_ind.a, ncol=mr_t.a)
      for(i in 1:mr_ind.a){
        for(x in (first.a[i]):last.a[i]){
          add_age_chtrue.a[i,x]<-as.numeric(which(!is.na(ind_mr.a[,x,i])))
          age_ch.a[i,x]<-add_age_chtrue.a[i,x]*ch.a[i,x]
        }
      }
      rm_2<-numeric(mr_ind.a) # remove those we never saw
      for(i in 1:mr_ind.a){
        if(sum(ch.a[i,])==0){
          rm_2[i]<-1
        }else {
          rm_2[i]<-0
        }
        #rm_2[i]<-sum(ch[i,])
      }
      if(sum(rm_2)>0){
        ch.a<-ch.a[-(which(rm_2==1)),]
        age_ch.a<-age_ch.a[-(which(rm_2==1)),]
        add_age_chtrue.a[-(which(rm_2==1)),]
      } else {}
      firstobs.a<-lastobs.a<-numeric(length(ch.a[,1]))
      for(i in 1:length(ch.a[,1])){
        firstobs.a[i]<-min(which(ch.a[i,]==1))
        lastobs.a[i]<-max(which(ch.a[i,]==1))
      }
      
      #for marking chicks
      tomark.j<-matrix(0,nrow=n.years, ncol=dim(indfates)[3])
      for(t in 1:n.years){
        #track which are observed, to mark
          for(i in 1:length(IND_MR[1,1,])){
            if(!is.na(IND_Nest[3,t,i]) & length(which(is.na(IND_Nest[1:2,t,i])))==2){
              tomark.j[t,i]<-1
            }  
          }
        }
      #mark chicks that were observed in nests
      #only want those that are new additions
      #those what were new to the pop:
      indsmarked.j<-which((tomark.j==1), arr.ind=T)[,2]
      
      age.j<-first.j<-last.j<-numeric()
      mr_t.j<-dim(indfates)[2] - 1 # Note change to prevent dimension mismatch
      mr_ind.j<-length(indsmarked.j)
      for(i in 1:mr_ind.j){
        g <- which(!is.na(indfates[1:3,,indsmarked.j[i]]), arr.ind = TRUE)
        age.j[i] <- g[1,1] #at marking (1=1year olds, 2=adults, 3= chicks)
        first.j[i] <- g[1,2]
        h <- which(indfates[1:3,1:n.years,indsmarked.j[i]]==1, arr.ind = TRUE) 
        last.j[i] <- max(h[,2])
      }
      ch.true.j<-matrix(0,ncol=mr_t.j, nrow=mr_ind.j)
      for(i in 1:mr_ind.j){
        ch.true.j[i,first.j[i]:last.j[i]]<-1
      }
      
      #code to track the ages, 3 is for chicks, 1 1years, 2 adults
      add_age_chtrue.j<-matrix(0,nrow=mr_ind.j, ncol=mr_t.j)
      for(i in 1:mr_ind.j){
        for(x in (first.j[i]):last.j[i]){
          add_age_chtrue.j[i,x]<-as.numeric(which(!is.na(indfates[1:3,x,indsmarked.j[i]])))[1]
        }
      }#where 3 is chick, 1 is 1yearold, 2 is adult in add_age_chtrue
      #since inital marking is constant prob
      #and resight is constant prob
      #p.ad #probability you were marked as adult
      #p.1 #probability of....? 
      in.mark.j<-ch.j<-matrix(0,nrow=mr_ind.j, ncol=mr_t.j)
      for(i in 1:mr_ind.j){
        if(age.j[i]==3){ #hatchlings
          in.mark.j[i,first.j[i]]<-rbinom(1,1,p.1*ch.true.j[i,first.j[i]])
        } #else{ #breeders
          # in.mark.j[i,first.j[i]]<-rbinom(1,1,p.ad*ch.true.j[i,first.j[i]])
        #}
        if(first.j[i]==mr_t.j) next
        for(t in (first.j[i]+1):last.j[i]){
          # if(add_age_chtrue.j[i,t]==3){ 
          #   in.mark.j[i,t]<-rbinom(1,1,p.ad*ch.true.j[i,t])
          #   # in.mark.j[i,t]<-rbinom(1,1,p.1*ch.true.j[i,t])
          # }
          in.mark.j[i,t]<-rbinom(1,1,p.ad*ch.true.j[i,t])
        }
        ch.j[i,]<-in.mark.j[i,]
      }
      
      #code to track the ages, 1 is for chicks, 2 1years, 3 adults
      age_ch.j<-matrix(0,nrow=mr_ind.j, ncol=mr_t.j)
      for(i in 1:mr_ind.j){
        for(x in (first.j[i]):last.j[i]){
          #add_age_chtrue[i,x]<-as.numeric(which(!is.na(ind_mr[,x,i])))[1]
          age_ch.j[i,x]<-add_age_chtrue.j[i,x]*ch.j[i,x]
        }
      }#where 3 is chick, 1 is 1yearold, 2 is adult in add_age_chtrue
      rm_2<-numeric(mr_ind.j)
      for(i in 1:mr_ind.j){
        if(sum(ch.j[i,])==0){
          rm_2[i]<-1
        }else {
          rm_2[i]<-0
        }
        #rm_2[i]<-sum(ch[i,])
      }
      if(sum(rm_2)>0){
        ch.j<-ch.j[-(which(rm_2==1)),]
        age_ch.j<-age_ch.j[-(which(rm_2==1)),]
        add_age_chtrue.j[-(which(rm_2==1)),]
      } else {}
      
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

temptrue<-simPopTrajectory(n.years=10,
                           age.init=c(20,20),
                           phi.1=0.5,phi.ad=0.6,
                           f=1)
temp.high<-simData(indfates = temptrue$indfates, 
              n.years=10,
              p.1=0.85,
              p.ad=0.85,
              p.prod=0.6,
              n.sam=0.6,
              p.count=0.6,
              BinMod=T)

#might need to divorce nest from mr, since can't fully track
# then just put above and clean up