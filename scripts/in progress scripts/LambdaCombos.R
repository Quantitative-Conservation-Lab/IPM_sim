
#Here is a little code to pick scenarios of parameter values
# that align with lambda

getNviable <- function() {
  # TODO
  # 1 - figure out how many viable combos there are - then add them to the scenarios table
  # 2 - what ARE the viable combinations
  
  # NOW THE EXCITING PART - how many valid combos
  #creat a matrix for all possible combinations of these 3
  coms <- expand.grid(fec = fec, phi1 = surv, phiad = surv)
  #output vector
  lams<-numeric(length(coms[,1]))
  #leslie matrix for each
  mats<-array(dim=c(length(coms[,1]), 2,2))
  
  for(i in 1:length(coms[,1])){ #loop over combinations
    #leslie matrix for each combination
    mats[i,1,]<-c(coms[i,"phi1"]*coms[i,"fec"],coms[i,"phi1"]*coms[i,"fec"])
    mats[i,2,]<-c(coms[i,"phiad"], coms[i,"phiad"])
    #Get the eigenvalue for each leslie matrix, which is the lambda
    lams[i]<-eigen(mats[i,,])$values[1]
  }
  
  return(cbind(coms, lams))
}  

#saveRDS(scenarios, scenarios.RDS)

      # #look at how many are within these arbitrary bounds
      # length(which(lams>=0.95 & lams<=1.05))
      # #put those that we care about in a matrix and pull out the values for the parameters
      # intv<-which(lams>=0.95 & lams<=1.05)
      # scenmat<-matrix(nrow=length(intv), ncol=3)
      # for(i in 1:length(scenmat[,1])){
      #   scenmat[i,1]<-phi1[coms[intv[i],1]]
      #   scenmat[i,2]<-fec[coms[intv[i],2]]
      #   scenmat[i,3]<-phiad[coms[intv[i],3]]

#HANNAH STUFF ########

#make a sequence of adult and juvenile survival
phiad<-phi1<-seq(0.05,0.95,by=0.05)
#make the same length for fecundity
fec<-seq(0.25,4, length.out=length(phiad))/2
#creat a matrix for all possible combinations of these 3
coms<-permutations(length(phiad),3,replace=T)
#output vector
lams<-numeric(length(coms[,1]))
#leslie matrix for each
mats<-array(dim=c(length(coms[,1]), 2,2))

for(i in 1:length(coms[,1])){ #loop over combinations
  #leslie matrix for each combination
  mats[i,1,]<-c(phi1[coms[i,1]]*fec[coms[i,2]],phi1[coms[i,1]]*fec[coms[i,2]])
  mats[i,2,]<-c(phiad[coms[i,3]], phiad[coms[i,3]])
  #Get the eigenvalue for each leslie matrix, which is the lambda
  lams[i]<-eigen(mats[i,,])$values[1]
}

#look at how many are within these arbitrary bounds
length(which(lams>=0.95 & lams<=1.05))
#put those that we care about in a matrix and pull out the values for the parameters
intv<-which(lams>=0.95 & lams<=1.05)
scenmat<-matrix(nrow=length(intv), ncol=3)
for(i in 1:length(scenmat[,1])){
  scenmat[i,1]<-phi1[coms[intv[i],1]]
  scenmat[i,2]<-fec[coms[intv[i],2]]
  scenmat[i,3]<-phiad[coms[intv[i],3]]
}

#look at the coverage of parameters and lambdas
#just comparing across for curiousity
hist(scenmat[,1], breaks=seq(0,1,by=0.1))
plot(sort(scenmat[,1]), sort(lams[intv]))
plot(scenmat[,1], scenmat[,2])
plot(scenmat[,3])
plot(scenmat[,1]*scenmat[,2], scenmat[,3])
plot(scenmat[,1],scenmat[,3])
plot(scenmat[,1]*scenmat[,2], lams[intv])
plot(scenmat[,1], lams[intv])
plot(scenmat[,2], lams[intv])
plot(scenmat[,3], lams[intv])


####Function to do this
#inputs for phi1, phiad, fec, are sequence or vector of values to use
#make sure they are all the same length, can change this later if we need to

#lam_range is a vector for low and high, where we will choose lambdas>=low and <=high
#can change the lambda to do discrete intervals or just use function multiple times
lambdacombosFN<-function(phi1, phiad, fec, lam_range){
  if(length(phi1)==length(phiad) && length(phi1)==length(fec) && length(phiad)==length(fec)){
  #creat a matrix for all possible combinations of these 3
  coms<-permutations(length(phiad),3,replace=T)
  #output vector
  lams<-numeric(length(coms[,1]))
  #leslie matrix for each
  mats<-array(dim=c(length(coms[,1]), 2,2))
  
  for(i in 1:length(coms[,1])){ #loop over combinations
    #leslie matrix for each combination
    mats[i,1,]<-c(phi1[coms[i,1]]*fec[coms[i,2]],phi1[coms[i,1]]*fec[coms[i,2]])
    mats[i,2,]<-c(phiad[coms[i,3]], phiad[coms[i,3]])
    #Get the eigenvalue for each leslie matrix, which is the lambda
    lams[i]<-eigen(mats[i,,])$values[1]
  }
  
  #look at how many are within these arbitrary bounds
  x<-length(which(lams>=lam_range[1] & lams<=lam_range[2]))
  #put those that we care about in a matrix and pull out the values for the parameters
  intv<-which(lams>=lam_range[1] & lams<=lam_range[2])
  scenmat<-matrix(nrow=length(intv), ncol=3)
  for(i in 1:length(scenmat[,1])){
    scenmat[i,1]<-phi1[coms[intv[i],1]]
    scenmat[i,2]<-fec[coms[intv[i],2]]
    scenmat[i,3]<-phiad[coms[intv[i],3]]
  }
  return(list(scenarios=scenmat, lambdas=lams[intv], numscen=x))
  }else{
    print("Error in input: Parameter input lengths must be the same")
    #can change this to not be true later if we need to, but seems like we would want
    #to look at the same number of parameters maybe
  }
}

lambdacombosFN(phi1=phi1, phiad=phiad, fec=fec, lam_range = c(0.95,1.05))
