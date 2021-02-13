#HAS - need to comment throughout, but got the MR data (I think) how it should be


#altered IPM sim with placeholder for nest success
#population projections and count code is from BPA WinBUGS appendix 2
#code adapted from WinBugsBPA Appendix 2, code to simulate IPM
# link to code is here: https://www.vogelwarte.ch/assets/files/publications/BPA/Web%20appendix%202.txt

library(nimble)

simIPMdata<-function(n.years, n.data, init.age, phi.1, phi.ad, p.1, p.ad, p.sur,
 max.nest.age, mean.clutch.size, phi.nest, n.sam){
  
  ti<-n.years
  ni<-n.years-1
  nd<-n.data
  Ni<-init.age
  phi1<-phi.1
  phi2<-phi.ad
  pjuv<-p.1
  psur<-p.sur
  padults<-p.ad
  max.nest.age <- max.nest.age
  mean.clutch.size <- mean.clutch.size
  phi.nest <- phi.nest
  n.sam <- n.sam
  
  # DERIVE TRUE FECUNDITY
  true.fec <- 1/2 * mean.clutch.size * phi.nest^max.nest.age
  
  # Combine the values to vectors
  PHI1 <- rep(phi1, ti)
  PHI2 <- rep(phi2, ti)
  PJUV <- rep(pjuv, ti)
  PADULTS <- rep(padults, ti)
  PSUR <- rep(psur, ti)
  FL1 <- rep(true.fec, ti+1) # THIS MATCHES WITH PROD DATA SIM
  FL2 <- rep(true.fec, ti+1) # THIS MATCHES WITH PROD DATA SIM
  
  ###########################
  # Population simulation
  ##########################
  
  # 1. Matrix model to estimate roughly the population size with a Leslie matrix
  
  N <- matrix(data = NA, nrow = length(Ni), ncol = ti+1)
  N[,1] <- Ni
  for (i in 1:ti){
    les <- matrix(data = c(PHI1[i]*FL1[i], PHI1[i]*FL2[i], PHI2[i], PHI2[i]), byrow=T, nrow = length(Ni), ncol = length(Ni))
    N[,(i+1)]<-les%*%N[,i]
    les
    eigen(les)
  }
  no.anim <- sum(N)
  no.ani <- round(no.anim*2.5)
  
  # 2. Define array for each individual
  ind <- array(data = NA, dim = c(5, ti+1, no.ani))   # infn about [1ye, adu, dead, rep, juv]
  
  # 3. Simulate the fates of initial individuals
  
  # 1-Year old individuals
  for (i in 1:Ni[1]){
    if(Ni[1]==0) break
    ind[1,1,i] <- 1
    ind[4,1,i] <- rpois(1, FL1[1])
    j <- 1
    z1 <- rbinom(1, 1, PHI2[j])           # survived or not
    if(z1==0){ind[3,j+1,i] <- 1
    next}                               # if dead
    ind[2,j+1,i] <- 1
    ind[4,j+1,i] <- rpois(1, FL2[j+1])  # no.of fledglings produced
    for (j in 2:ti){
      z2 <- rbinom(1, 1, PHI2[j])
      if(z2==0){
        ind[3,j+1,i] <- 1
        break}
      ind[2,j+1,i] <- 1
      ind[4,j+1,i] <- rpois(1, FL2[j+1])
    }
  }
  
  # Adults
  for (i in 1:Ni[2]){
    if(Ni[2]==0) break
    m <- Ni[1]+i
    ind[2,1,m] <- 1
    ind[4,1,m] <- rpois(1, FL2[1])
    for (j in 1:ti){
      z<-rbinom(1, 1, PHI2[j])
      if(z==0){
        ind[3,j+1,m] <- 1
        break}
      ind[2,j+1,m] <- 1
      ind[4,j+1,m] <- rpois(1, FL2[j+1])
    }
  }
  
  
  # 4. Simulate the fates of newborn
  juv <- rep(0,ti)
  for (g in 1:ti){
    m <- sum(Ni)+sum(juv)
    juv[g] <- sum(ind[4,g,], na.rm = T)
    if(juv[g]==0) next
    for (i in 1:juv[g]){
      ind[5,g,m+i] <- 1
      # first time step
      z1 <- rbinom(1, 1, PHI1[g])
      if(z1==0){ind[3,g+1,m+i] <- 1
      next}
      ind[1,g+1,m+i] <- 1
      ind[4,g+1,m+i] <- rpois(1, FL1[g+1])
      if(g>=ti) next                       # if cycle too long
      # second time step
      z2 <- rbinom(1, 1, PHI2[g+1])
      if(z2==0){ind[3,g+2,m+i] <- 1
      next}
      ind[2,g+2,m+i] <- 1
      ind[4,g+2,m+i] <- rpois(1, FL2[g+2])
      # third time step and above
      if(g>=(ti-1)) next
      for (k in (g+2):ti){
        z3 <- rbinom(1, 1, PHI2[k])
        if(z3==0){ind[3,k+1,m+i] <- 1
        break}
        ind[2,k+1,m+i] <- 1
        ind[4,k+1,m+i] <- rpois(1, FL2[k+1])
      }
    }
  }
  
  # 5. Total number of animals
  Ntotal <- sum(Ni)+sum(ind[4,1:ti,], na.rm = T)
  
  # Reshape the array, remove empty cells
  IND <- ind[,,1:Ntotal]
  rownames(IND) <- c("1-Year", "Adu", "Dead", "Rep", "Chicks")
  
  ###################################################
  #  Create completely independent samples
  ###################################################
  
  # Three independent samples
  Ntot <- dim(IND)[3]
  xt1 <- matrix(data = seq(1:Ntot), ncol = 1)
  resamp1 <- resamp2 <- resamp3 <- numeric()
  
  # Sample 1
  resamp1 <- sample(xt1, nd[1], replace = F)
  # Sample 2
  resamp2 <- sample(xt1[-resamp1,], nd[2], replace = F)
  
  # Individuals selected for capture-recapture, survey and reproductive success data
  IND_MR <- IND[,,resamp1]
  IND_Count <- IND[,,resamp2]
  IND_Nest <- IND[,,resamp3]
  #only use the yoy and adults that are alive
  
  #######Mark resight data
  #pretend marking and resighting:
  
  #using the individual population array for MR data
  mr_classes<-dim(IND_MR)[1] - 3 #since we wont see Dead and the reproduction doesnt matter
  #so we have 3 classes, chicks, 1year olds, and adults
  ind_mr<-IND_MR[c(5,1,2),1:ti,]
  ind_mr[1,,]<-NA #not banding chicks, so remove them
  rm<-numeric(dim(ind_mr)[3])
   for(i in 1:dim(ind_mr)[3]){
     if(length(which(!is.na(ind_mr[2:3,,i])))==0){
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
    g <- which(!is.na(ind_mr[2:(mr_classes+1),,i]), arr.ind = TRUE)
    age[i] <- g[1,1]+1
    first[i] <- g[1,2]
    h <- which(ind_mr[2:(mr_classes+1),,i]==1, arr.ind = TRUE)
    last[i] <- max(h[,2])
  }  
  #remove those that we never banded as a chick
  ch.true<-matrix(0,ncol=mr_t, nrow=mr_ind)
  for(i in 1:mr_ind){
    ch.true[i,first[i]:last[i]]<-1
  }  
  
  #since inital marking is constant prob
  #and resight is constant prob
  p.juv<-PJUV[1] #probability you were marked as YoY
  p.ad<-PADULTS[1] #probability you were marked as adult
  in.mark<-ch<-matrix(0,nrow=mr_ind, ncol=mr_t)
  for(i in 1:mr_ind){
    if(age[i]==2){
    in.mark[i,first[i]]<-rbinom(1,1,p.juv*ch.true[i,first[i]])
    } else{
      in.mark[i,first[i]]<-rbinom(1,1,p.ad*ch.true[i,first[i]])
    }
    if(first[i]==mr_t) next
    for(t in (first[i]+1):last[i]){
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
  first<-last<-numeric(length(ch[,1]))
  for(i in 1:length(ch[,1])){
    first[i]<-min(which(ch[i,]==1))
    last[i]<-max(which(ch[i,]==1))
  }
  #############################
  ######################################################
  # Create population survey data
  ######################################################
  #n.sam is the number of times that the population was sampled in a year
  # I rewrote this, it just needs to be rbinom(n.sam, TrueCount,p_sur)
  #HAS - should output the true count data for comparison
  TRUE_Count<-matrix(nrow=2, ncol=ti) #first row, number of YOY; second row is adults
  SUR<-matrix(nrow=n.sam, ncol=ti)#matrix for surveys by survey and year
  #SUR <- rep(0,ti)
  for(u in 1:ti){
    TRUE_Count[1,u]<-sum(IND_Count[1,u,], na.rm = T)
    TRUE_Count[2,u]<-sum(IND_Count[2,u,], na.rm = T)
    SUR[,u]<-rbinom(n.sam, sum(TRUE_Count[,u]), PSUR[u])
  }
  
  ###########################
  ##########################
  #ABBY NEST MODEL could go here
  
  # AEB NOTE - this was super buggy so the nest survival stuff currently comes from
  # productivityDataSim.R and we just assume independence
  
  #Get all the data together, just MR and Survey Counts until nest in here
  return(list(ch=ch, SUR=SUR, age_ch=age_ch, first=first, last=last,
              n.sam = n.sam, true.fec=true.fec, 
              N = N))
}
