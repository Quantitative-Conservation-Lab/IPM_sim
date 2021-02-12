#rewriting IPM simulation code

#parameters to put in
#initial set up
#Leslie matrix
#array for age classes, then specific age class demostoch
#divide them into data
#add in observation process for each data type

#INPUT FOR FUNCTION
#1.) number of years to simulate
#2.) vector for the number individuals to use to simulate data w or w/o nest model
#3.) vector for initial age distribution of the data, 
#   true number of individuals in each age class in the initial year
#4.) constant parameter values to use
  #a.) 1year old survival
  #b.) adult survival
  #c.) 1 year old survival detection (to include?)
  #d.) adult survival detection
  #e.) count survey detection
  #f.) maximum nest age
  #g.) mean clutch size
  #h.) nest survival
  #i.) number of nest surveys
#5.) additional arguments to maybe include: 
      #a.) T/F for which nest survival/productivity data to simulate
      #b.) T/F for which count survey model to use
      #c.) T/F for whether to include 1 year old survival as an observed parameter

n.years<-10

#for data, I think it might work better to specify the proportion of 
# the population that we want for each; MR, Count, Productivity, i.e.:
n.data.types<-c(1/3,1/3,1/3)
#which will be multiplied by the total number of individuals

#can also specify it this way, function will take either one
n.data.types<-c(100,100,100) #MR, Count, Productivity, number if 
#individuals to grab from the total population for data points
#does not specify if they survived very long, so can limit the outpute

age.init<-c(150,150) #1 year olds, adults
#note that this is just to inform the leslie matrix
#the age distribution we will start with will be the stable age distribution
#could probably remove, since it really doesnt matter
#what we want will be the total number of individuals to start from
#then can pull out the stable age distribution from that, so the input would be
age.init<-300 #then split equally 1 year olds and adults

#MR model and population projection
phi.1<-0.5 
phi.ad<-0.7
p.1<-0.7 #for if 1 year olds are detected in mr
p.ad<-0.9
#count model
p.count<-0.55 #for binomial: y_count~bin(truecount, p.count)
#alternative: for normal count model
sig<-0.5 #observation error for normal count model
# where for normal: y_count~norm(truecount, sig)
#if you want binomial
BinMod=T #tell it which you want, if false it is the normal
#if you want normal:
BinMod=F

#productivity
#if you want nest survival:
productivity=T
#if simple productivity model
productivity=F

max.nest.age<-30
mean.clutch.size <- 2.5
phi.nest <- 0.975
n.sam <- 3
#alternative for if different model for nest success
#fecundity for 1yearolds and ads
f<-1.2
p.prod<-0.4 #probability of detection 

#In function add in T/F for which sims to do
#DO FUNCTION STUFF HERE simIPMdat<-function(){}
productivity=F

nmin1.years<-n.years-1
nplus1.years<-n.years+1

#derive true fecundity from the parameters for productivity or use other f
if(productivity==T){
  fec<-1/2*mean.clutch.size*phi.nest^max.nest.age
}else{
  fec<-f/2 
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
eigen(lesmat) #for lambda
no.animals<-sum(N) #number of animals ever in the system at anytime
no.ani.max<-round(no.animals*2.5) #include more for simulation of offspring
####stable age distribution is:
sad<-N[,nplus1.years]/sum(N[,nplus1.years]) # the proportion of the popualtion that is in each age class
#lambda from N, matches eigenvalue
nlam<-N[,nplus1.years]/N[,n.years]
#note that this is without demographic stochasticity and 
#when we simulate the fates the actual lambda will change based on that

#simulate what happens to each individual:
if(N[1]==0 || N[2]==0){print("noindividuals!!!!!!!!!")}
#dont run if there arent any individuals given 

#make an array for each individual, for each year, what state it is in
# states: [1yrolds, Adult, chicks, Dead]
# Oneyearolds are the 1st years, chicks are how many chicks per individual,
#dead is when they die, duh

#simulate their fates 
#1yearolds - 
#we want to track whether they survive and become adults in each year [2,t+1,]
#AND if they reproduce, how many chicks do they reproduce [3,t,] and add individuals
#that were reproduced into the chicks [4,t,]
# track if they did not survive 


#do individual fates for each per time step, track the number of chicks
indfates<-array(dim=c(4, nplus1.years, no.ani.max))
#add in the stable age distribution to start from
age1<-round(sum(age.init)*sad[1])
age2<-round(sum(age.init)*sad[2])
indfates[1,1,1:age1]<-1
indfates[2,1,(age1+1):sum(age1+age2)]<-1

#this commented out bit is for if we dont want to start at stable age distribution
#indfates[1,1,1:age.init[1]]<-1
#indfates[2,1,(age.init[1]+1):sum(age.init)]<-1
#source in these functions, they are annoying here
oneyr_fatefn<-function(ind,time){
  #if(which(indfates[,time,ind]==1)){
    indfates[3,time,ind]<-rpois(1,fec)
    #indfates[3,time,ind]<-ifelse(indfates[4,1,ind]>0,1,0) #indicates that they reproduced
    zsurv1<-rbinom(1,1,phi.ad)
    indfates[2,time+1,ind]<-ifelse(zsurv1==1,1,NA) #do they survive?
    indfates[4,time+1,ind]<-ifelse(zsurv1==0, 1,NA) #or die?
  return(indfates[,,ind])
}

adfatefn<-function(ind,time){
    indfates[3,time,ind]<-rpois(1,fec)
    #indfates[3,time,ind]<-ifelse(indfates[4,1,ind]>0,1,0) #indicates that they reproduced
    zsurv1<-rbinom(1,1,phi.ad)
    indfates[2,time+1,ind]<-ifelse(zsurv1==1,1,NA) #do they survive?
    indfates[4,time+1,ind]<-ifelse(zsurv1==0, 1,NA) #or die?
    return(indfates[,,ind])
}

deadfn<-function(ind,time){
  indfates[4,time+1,ind]<-1
  return(indfates[,,ind])
}

chickfatefn<-function(ind, time){
  zsurv1<-rbinom(1,1,phi.1)
  indfates[1,time+1,ind]<-ifelse(zsurv1==1,1,NA) #do they survive to become 1yrolds
  indfates[4,time+1,ind]<-ifelse(zsurv1==0, 1,NA) #or die?
  return(indfates[,,ind])
}


inpop<-numeric(n.years)
#inpop[1]<-sum(age.init)
inpop[1]<-sum(age1+age2)

chickst<-numeric(n.years)
tempstep<-matrix(nrow=no.ani.max, ncol=(n.years+1))
for(t in 1:n.years){
  for(i in 1:inpop[t]){
    tempstep[i,t]<-which(indfates[c(1,2,4),t,i]==1)
    if(tempstep[i,t]==1){
      indfates[,(t:(t+1)),i]<-oneyr_fatefn(ind=i, time=t)[,(t:(t+1))]
    }else if(tempstep[i,t]==2){
        indfates[,(t:(t+1)),i]<-adfatefn(ind=i, time=t)[,(t:(t+1))]
      # }else if (tempstep==3){
      # indfates[,(t:(t+1)),i]<-chickfatefn(ind=i, time=t)[,(t:(t+1))]
    } else if(tempstep[i,t]==3){
      indfates[,(t:(t+1)),i]<-deadfn(ind=i, time=t)[,(t:(t+1))]
    # }else{
     }
  } #i loop over population size at time t
  #see how many chicks were produced and add them to the population
  chickst[t]<-sum(indfates[3,t,], na.rm=T)
  inpop[t+1]<-inpop[t]+chickst[t]
  indfates[3,t,(inpop[t]+1):inpop[t+1]]<-1
  for(i in (inpop[t]+1):inpop[t+1]){
    indfates[,(t:(t+1)),i]<-chickfatefn(ind=i, time=t)[,(t:(t+1))]
  }
} 

ad<-ch<-j<-dead<-numeric(n.years)
for(i in 1:n.years){
  ad[i]<-sum(indfates[2,i,], na.rm=T)
  ch[i]<-sum(indfates[3,i,], na.rm=T)
  j[i]<-sum(indfates[1,i,], na.rm=T)
  dead[i]<-sum(indfates[4,i,], na.rm=T)
  
}
plot(dead, pch=16)
points(ad, col="red", pch=16)
points(ch,col="blue", pch=16)
points(j,col="purple", pch=16)

adj<-numeric(length(ad))
for(i in 2:length(adj)){
  adj[i]<-sum(ad[i]+j[i])/sum(ad[i-1]+j[i-1])
}
mean(adj[2:10])
#close *enough* to actual lambda, again demographic stochasticity is the cause!


# Three independent samples
Ntot <- dim(indfates)[3]
xt1 <- matrix(data = seq(1:Ntot), ncol = 1)
resamp1 <- resamp2 <- resamp3 <- numeric()
nds<-numeric(length(n.data.types))
if(sum(n.data.types)>1){
  nds<-n.data.types
}else{
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

#for adults only
#AD_only=T
#for adults and 1yearolds 
#AD_only=F


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
#Abby can add in if we want to include it
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
#for the model it would be:
#where the data to give it is the observed OBS_nestlings and R_obs
# for (t in 1:(nyears-1)){
#   OBS_nestlings[t] ~ dpois(rho[t])
#   rho[t] <- R_obs[t]*fecundity[t]
# }







