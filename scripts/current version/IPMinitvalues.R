#Provide known state as data  
state.data <- function(EH){
  state <- EH
  for (i in 1:dim(EH)[1]){
    first <- min(which(EH[i,]==1))
    last <- max(which(EH[i,]==1))
    state[i,first:last] <- 1
  }
  state[state==0] <- NA
  return(state)
}

#age for age specific survival probabilities
ageunknown<-function(age_ch){
  allages<-age_ch
  f1<-l1<-numeric(dim(allages)[1])
  t<-dim(age_ch)[2]
  for(i in 1:dim(allages)[1]){
    f1[i]<-min(which(allages[i,]>=2))
    l1[i]<-max(which(allages[i,]>0))
    if((f1[i]+1)<=t){
      allages[i,(f1[i]+1)]<-3
    } else{}
    if((f1[i]+2)<=t){
      allages[i,(f1[i]+2):t]<-3
    } else{}
  }
  allages<-allages-1
  allages[which(allages==-1)]<-0
  #allages[which(allages==1)]<-0
  # allages[which(allages==2)]<-1
  # allages[which(allages==3)]<-2
  return(allages)
}

getHinits <- function(H) {
  Hinits <- H
  Hinits[!is.na(H)] <- NA
  Hinits[is.na(H)] <- 1
  return(Hinits)
}
