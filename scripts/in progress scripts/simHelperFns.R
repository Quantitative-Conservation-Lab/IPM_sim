oneyr_fatefn<-function(ind,time){
  indfates[3,time,ind]<-rpois(1,fec)
  zsurv1<-rbinom(1,1,phi.ad)
  indfates[2,time+1,ind]<-ifelse(zsurv1==1,1,NA) #do they survive?
  indfates[4,time+1,ind]<-ifelse(zsurv1==0, 1,NA) #or die?
  indfates[1,time+1,ind]<-NA #cant stay 1yearold
  return(indfates[,time:(time+1),ind])
}

adfatefn<-function(ind,time){
  indfates[3,time,ind]<-rpois(1,fec)
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

chickfatefn<-function(ind, time){
  zsurv1<-rbinom(1,1,phi.1)
  indfates[1,time+1,ind]<-ifelse(zsurv1==1,1,NA) #do they survive to become 1yrolds
  indfates[4,time+1,ind]<-ifelse(zsurv1==0, 1,NA) #or die?
  indfates[2,time+1,ind]<-NA #not adult yet
  indfates[3,time+1,ind]<-NA
  return(indfates[,(time+1),ind])
}
