
#Here is a little code to pick scenarios of parameter values
# that align with lambda

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

