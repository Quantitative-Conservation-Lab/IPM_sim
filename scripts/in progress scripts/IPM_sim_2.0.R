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

n.years<-20
n.data.types<-c(10,10,10) #MR, Count, Productivity
age.init<-c(50,50) #juveniles, adults
#MR model and population projection
phi.1<-0.65
phi.ad<-0.75
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
  fec<-f/2 #?
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
# states: [Juv1, Adult, chicks, Dead]
# Juv1 are the 1st years, chicks are how many chicks per individual,
#dead is when they die, duh
#indfates<-array(dim=c(4, nplus1.years, no.ani.max))

#simulate their fates 

#Juveniles/1yearolds - 
#we want to track whether they survive and become adults in each year [2,t+1,]
#AND if they reproduce, how many chicks do they reproduce [3,t,] and add individuals
#that were reproduced into the chicks [4,t,]
# track if they did not survive 
if(N[1]==0 || N[2]==0){print("noindividuals!!!!!!!!!")}
# for(i in 1:age.init[1]){
#   #what happens in the first year
#   indfates[1,1,i]<-1 #initial time step, they are there
#   #do they reproduce chicks in this initial year
#   indfates[4,1,i]<-rpois(1,fec)
#   indfates[3,1,i]<-ifelse(indfates[4,1,i]>0,1,0) #indicates that they reproduced
#   zsurv1<-rbinom(1,1,phi.ad)
#   indfates[2,2,i]<-ifelse(zsurv1==1,1,NA) #do they survive?
#   indfates[5,2,i]<-ifelse(zsurv1==0, 1,NA) #or die?
#   for(t in 2:n.years){
#     if(is.na(indfates[5,t,i]) && indfates[2,t,i]==1){ #not dead
#       indfates[4,t,i]<-rpois(1,fec) #no. offspring
#       indfates[3,t,i]<-ifelse(indfates[4,t,i]>0,1,0)
#       zsurv1<-rbinom(1,1,phi.ad)
#       indfates[2,t+1,i]<-ifelse(zsurv1==1,1,NA)
#       indfates[5,t+1,i]<-ifelse(zsurv1==0, 1,NA)
#     }else{
#       indfates[5,t+1,i]<-1 #if they are dead, they are still dead
#     }
#   }
# }
# 
# #now for adults
# #we filled in the individuals that were juveniles, so now the 
# #next group of individuals are adults 
# for(i in ((age.init[1]+1):(age.init[1]+age.init[2]))){
#   indfates[2,1,i]<-1
#   #do they reproduce chicks in this initial year
#   indfates[4,1,i]<-rpois(1,fec)
#   indfates[3,1,i]<-ifelse(indfates[4,1,i]>0,1,0) #indicates that they reproduced
#   zsurv1<-rbinom(1,1,phi.ad)
#   indfates[2,2,i]<-ifelse(zsurv1==1,1,NA) #do they survive?
#   indfates[5,2,i]<-ifelse(zsurv1==0, 1,NA) #or die?
#   for(t in 2:n.years){
#     if(is.na(indfates[5,t,i]) && indfates[2,t,i]==1){ #not dead
#       indfates[4,t,i]<-rpois(1,fec) #no. offspring
#       indfates[3,t,i]<-ifelse(indfates[4,t,i]>0,1,0)
#       zsurv1<-rbinom(1,1,phi.ad)
#       indfates[2,t+1,i]<-ifelse(zsurv1==1,1,NA)
#       indfates[5,t+1,i]<-ifelse(zsurv1==0, 1,NA)
#     }else{
#       indfates[5,t+1,i]<-1 #if they are dead, they are still dead
#     }
#   }
# }
# 
# #simulate the fates of the newborn chicks
# #add in the chicks that survive to the population
# chicks<-numeric(n.years)
# for(t in 1:n.years){
#   chicks[t]<-sum(indfates[4,t,], na.rm=T)
#   
# }
# #get the total number to use in the array
# sum.chicks.tot<-sum(chicks)
# for(t in 1:n.years){
# for(i in 1:chicks[t]){
#   if(chicks[t]==0){
#     indfates[4,t,(sum(age.init)+1):chicks[t]]<-NA
#   }else{
#     indfates[4,1,(sum(age.init))+i]<-1
#     zsurv1<-rbinom(1,1,phi.1)
#     indfates[1,2,(sum(age.init))+i]<-ifelse(zsurv1==1,1,NA)
#     indfates[5,2,(sum(age.init))+i]<-ifelse(zsurv1==0,1,NA)
#     for(n in 3:n.years){
#       if(is.na(indfates[5,n-1,(sum(age.init))+i])){ #not dead
#         indfates[4,n-1,(sum(age.init))+i]<-rpois(1,fec) #no. offspring
#         indfates[3,n-1,(sum(age.init))+i]<-ifelse(indfates[4,n-1,(sum(age.init))+i]>0,1,0)
#         zsurv1<-rbinom(1,1,phi.ad)
#         indfates[2,n,(sum(age.init))+i]<-ifelse(zsurv1==1,1,NA)
#         indfates[5,n,(sum(age.init))+i]<-ifelse(zsurv1==0,1,NA)
#       }else{
#         indfates[5,n,(sum(age.init))+i]<-1 #if they are dead, they are still dead
#       }
#     }
#   }
# }
#}
#do individual fates for each per time step, track the number of chicks
indfates<-array(dim=c(4, nplus1.years, no.ani.max))

indfates[1,1,1:age.init[1]]<-1
indfates[2,1,(age.init[1]+1):sum(age.init)]<-1
juvfatefn<-function(ind,time){
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
  indfates[1,time+1,ind]<-ifelse(zsurv1==1,1,NA) #do they survive to become juvs
  indfates[4,time+1,ind]<-ifelse(zsurv1==0, 1,NA) #or die?
  return(indfates[,,ind])
}

#indfates[,1:2,1]<-juvfatefn(ind=1,time=1)[,1:2]
##indfates[,2:3,1]<-adfatefn(ind=1, time=2)[,2:3]
#indfates[,1:2,21]<-chickfatefn(ind=21, time=1)[,1:2]
#indfates[,2:3,21]<-juvfatefn(ind=21,time=2)[,2:3]
#indfates[,3:4,21]<-adfatefn(ind=21, time=3)[,3:4]

inpop<-numeric(n.years)
inpop[1]<-sum(age.init)
chickst<-numeric(n.years)
tempstep<-matrix(nrow=no.ani.max, ncol=(n.years+1))
for(t in 1:n.years){
  for(i in 1:inpop[t]){
    tempstep[i,t]<-which(indfates[c(1,2,4),t,i]==1)
    if(tempstep[i,t]==1){
      indfates[,(t:(t+1)),i]<-juvfatefn(ind=i, time=t)[,(t:(t+1))]
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
