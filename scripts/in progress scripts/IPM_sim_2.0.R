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
  #a.) juvenile survival
  #b.) adult survival
  #c.) juvenile survival detection (to include?)
  #d.) adult survival detection
  #e.) count survey detection
  #f.) maximum nest age
  #g.) mean clutch size
  #h.) nest survival
  #i.) number of nest surveys
#5.) additional arguments to maybe include: 
      #a.) T/F for which nest survival/productivity data to simulate
      #b.) T/F for which count survey model to use
      #c.) T/F for whether to include juvenile survival as an observed parameter

n.years<-10
n.data.types<-c(10,10,10) #MR, Count, Productivity
age.init<-c(10,10) #juveniles, adults
#MR model and population projection
phi.1<-0.7
phi.ad<-0.85
p.1<-0.76 #for if juveniles are detected in mr
p.ad<-0.8
#count model
p.count<-0.55 #for binomial: y_count~bin(truecount, p.count)
#alternative: for normal count model
sig<-0.5 #observation error for normal count model
# where for normal: y_count~norm(truecount, sig)
#productivity
max.nest.age<-30
mean.clutch.size <- 2.5
phi.nest <- 0.975
n.sam <- 3
#alternative for if different model for nest success
f<-1 #just adults reproduce?
p.prod<-0.4 #probability of detection the reproduction

#In function add in T/F for which sims to do
#DO FUNCTION STUFF HERE simIPMdat<-function(){}
productivity=F

nmin1.years<-n.years-1
nplus1.years<-n.years+1

#derive true fecundity from the parameters for productivity or use other f
if(productivity==T){
  fec<-1/2*mean.clutch.size*phi.nest^max.nest.age
}else{
  fec<-f #?
}

#Leslie matrix for population size
N<-matrix(nrow=length(age.init), ncol=nplus1.years)
N[,1]<-age.init
lesmat<-matrix(nrow=2, ncol=2)
lesmat[1,]<-c(phi.1*f, phi.1*f)
lesmat[2,]<-c(phi.ad, phi.ad)
for(t in 1:n.years){
  N[,(t+1)]<-lesmat%*%N[,t]
}
eigen(lesmat) #for lambda
no.animals<-sum(N) #number of animals ever in the system at anytime
no.ani.max<-round(no.animals*2.5) #include more for simulation of offspring

#simulate what happens to each individual

#make an array for each individual, for each year, what state it is in
# states: [Juv1, Adult, Rep, Newborns, Dead]
# Juv1 are the 1st years, Rep is how many chicks were produced
#Newborns are the chicks that were produced and which go into the population
#dead is when they die, duh
indfates<-array(dim=c(5, nplus1.years, no.ani.max))

#simulate their fates 

#Juveniles/1yearolds - 
#we want to track whether they survive and become adults in each year [2,t+1,]
#AND if they reproduce, how many chicks do they reproduce [3,t,] and add individuals
#that were reproduced into the chicks [4,t,]
# track if they did not survive 
if(N[1]==0 || N[2]==0){print("noindividuals!!!!!!!!!")}
for(i in 1:age.init[1]){
  #what happens in the first year
  indfates[1,1,i]<-1 #initial time step, they are there
  #do they reproduce chicks in this initial year
  indfates[4,1,i]<-rpois(1,fec)
  indfates[3,1,i]<-ifelse(indfates[4,1,i]>0,1,0) #indicates that they reproduced
  zsurv1<-rbinom(1,1,phi.ad)
  indfates[2,2,i]<-ifelse(zsurv1==1,1,NA) #do they survive?
  indfates[5,2,i]<-ifelse(zsurv1==0, 1,NA) #or die?
  for(t in 2:n.years){
    if(is.na(indfates[5,t,i]) && indfates[2,t,i]==1){ #not dead
      indfates[4,t,i]<-rpois(1,fec) #no. offspring
      indfates[3,t,i]<-ifelse(indfates[4,t,i]>0,1,0)
      zsurv1<-rbinom(1,1,phi.ad)
      indfates[2,t+1,i]<-ifelse(zsurv1==1,1,NA)
      indfates[5,t+1,i]<-ifelse(zsurv1==0, 1,NA)
    }else{
      indfates[5,t+1,i]<-1 #if they are dead, they are still dead
    }
  }
}

#now for adults
#we filled in the individuals that were juveniles, so now the 
#next group of individuals are adults 
for(i in ((age.init[1]+1):(age.init[1]+age.init[2]))){
  indfates[2,1,i]<-1
  #do they reproduce chicks in this initial year
  indfates[4,1,i]<-rpois(1,fec)
  indfates[3,1,i]<-ifelse(indfates[4,1,i]>0,1,0) #indicates that they reproduced
  zsurv1<-rbinom(1,1,phi.ad)
  indfates[2,2,i]<-ifelse(zsurv1==1,1,NA) #do they survive?
  indfates[5,2,i]<-ifelse(zsurv1==0, 1,NA) #or die?
  for(t in 2:n.years){
    if(is.na(indfates[5,t,i]) && indfates[2,t,i]==1){ #not dead
      indfates[4,t,i]<-rpois(1,fec) #no. offspring
      indfates[3,t,i]<-ifelse(indfates[4,t,i]>0,1,0)
      zsurv1<-rbinom(1,1,phi.ad)
      indfates[2,t+1,i]<-ifelse(zsurv1==1,1,NA)
      indfates[5,t+1,i]<-ifelse(zsurv1==0, 1,NA)
    }else{
      indfates[5,t+1,i]<-1 #if they are dead, they are still dead
    }
  }
}

#simulate the fates of the newborn chicks
#add in the chicks that survive to the population
chicks<-numeric(n.years)
for(t in 1:n.years){
  chicks[t]<-sum(indfates[4,t,], na.rm=T)
  
}
#get the total number to use in the array
sum.chicks.tot<-sum(chicks)
for(t in 1:n.years){
for(i in 1:chicks[t]){
  if(chicks[t]==0){
    indfates[4,t,(sum(age.init)+1):chicks[t]]<-NA
  }else{
    indfates[4,1,(sum(age.init))+i]<-1
    zsurv1<-rbinom(1,1,phi.1)
    indfates[1,2,(sum(age.init))+i]<-ifelse(zsurv1==1,1,NA)
    indfates[5,2,(sum(age.init))+i]<-ifelse(zsurv1==0,1,NA)
    for(n in 2:n.years){
      if(is.na(indfates[5,n,(sum(age.init))+i]) && (indfates[1,n,(sum(age.init))+i]==1 || indfates[2,n,(sum(age.init))+i]==1)){ #not dead
        indfates[4,n,(sum(age.init))+i]<-rpois(1,fec) #no. offspring
        indfates[3,n,(sum(age.init))+i]<-ifelse(indfates[4,n,(sum(age.init))+i]>0,1,0)
        zsurv1<-rbinom(1,1,phi.ad)
        indfates[2,n+1,(sum(age.init))+i]<-ifelse(zsurv1==1,1,NA)
        indfates[5,n+1,(sum(age.init))+i]<-ifelse(zsurv1==0,1,NA)
      }else{
        indfates[5,n+1,(sum(age.init))+i]<-1 #if they are dead, they are still dead
      }
    }
  }
}
}
#do individual fates for each per time step, track the number of chicks



