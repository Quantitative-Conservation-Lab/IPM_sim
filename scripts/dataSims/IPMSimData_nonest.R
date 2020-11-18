#HAS - need to comment throughout, but got the MR data (I think) how it should be


#altered IPM sim with placeholder for nest success
#population projections and count code is from BPA WinBUGS appendix 2
#code adapted from WinBugsBPA Appendix 2, code to simulate IPM
# link to code is here: https://www.vogelwarte.ch/assets/files/publications/BPA/Web%20appendix%202.txt

n.years=10; n.data=c(50,50,50);init.age = c(100,100); 
phi.1=0.7; phi.ad=0.76;p.1=0.98; p.ad=0.65;
f.1=0.8;f.ad=0.8; p.sur=0.8; p.prod=0.77

simIPMdata<-function(n.years, n.data, init.age, phi.1, phi.ad, p.1,p.ad,
 f.1, f.ad, p.sur, p.prod){
  ti<-n.years
  ni<-n.years-1
  nd<-n.data
  Ni<-init.age
  phi1<-phi.1
  phi2<-phi.ad
  fl1<-f.1
  fl2<-f.ad
  pjuv<-p.1
  padults<-p.ad
  psur<-p.sur
  pprod<-p.prod
  
  # Combine the values to vectors
  PHI1 <- rep(phi1, ti)
  PHI2 <- rep(phi2, ti)
  PJUV <- rep(pjuv, ti)
  PADULTS <- rep(padults, ti)
  PSUR <- rep(psur, ti)
  PPROD <- rep(pprod, ti)
  FL1 <- rep(fl1, ti+1)
  FL2 <- rep(fl2, ti+1)
  
  ###########################
  # Population simulation
  ##########################
  
  # 1. Matrix model to estimate roughly the population size with a Leslie matrix
  
  N <- matrix(data = NA, nrow = length(Ni), ncol = ti+1)
  N[,1] <- Ni
  for (i in 1:ti){
    les <- matrix(data = c(PHI1[i]*FL1[i], PHI1[i]*FL2[i], PHI2[i], PHI2[i]), byrow=T, nrow = length(Ni), ncol = length(Ni))
    N[,(i+1)]<-les%*%N[,i]
  }
  no.anim <- sum(N)
  no.ani <- round(no.anim*5)
  
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
  # Sample 3
  resamp3 <- sample(xt1[c(-resamp1, -resamp2),], nd[3], replace = F)
  
  # Individuals selected for capture-recapture, survey and reproductive success data
  IND_MR <- IND[,,resamp1]
  IND_Count <- IND[,,resamp2]
  IND_Nest <- IND[,,resamp3]
  
  #######Mark resight data
  #pretend marking and resighting:
  #Chicks are banded in year 1
  # years>1, chicks are banded and older are resighted
  
  #using the individual population array for MR data
  mr_classes<-dim(IND_MR)[1] - 3 #since we wont see Dead and the reproduction doesnt matter
  #so we have 3 classes, chicks, 1year olds, and adults
  ind_mr<-IND_MR[c(5,1,2),,]
  rm<-numeric(dim(ind_mr)[3])
  for(i in 1:dim(ind_mr)[3]){
    if(!is.na(ind_mr[2,1,i])){
      rm[i]<-1
      ind_mr[2,1,i]<-NA
    }else{}
    if(!is.na(ind_mr[3,1,i])){
      ind_mr[3,1,i]<-NA
      rm[i]<-1
    }else{}
  }
  if(sum(rm>0)){
    ind_mr<-ind_mr[,,-which(rm==1)]
  }else{}
  age<-first<-last<-numeric()
  mr_t<-dim(ind_mr)[2]
  mr_ind<-dim(ind_mr)[3] 
  for(i in 1:mr_ind){
    g <- which(!is.na(ind_mr[1:(mr_classes+1),,i]), arr.ind = TRUE)
    age[i] <- g[1,1]
    first[i] <- g[1,2]
    h <- which(ind_mr[1:(mr_classes+1),,i]==1, arr.ind = TRUE)
    last[i] <- max(h[,2])
  }  
  #remove those that we never banded as a chick
  ch.true<-matrix(0,ncol=mr_t, nrow=mr_ind)
  for(i in 1:mr_ind){
    ch.true[i,first[i]:last[i]]<-1
  }  
  
  #since inital marking is constant prob
  #and resight is constant prob
  p.juv<-PJUV[1]
  p.ad<-PADULTS[1]
  in.mark<-ch<-matrix(0,nrow=mr_ind, ncol=mr_t)
  for(i in 1:mr_ind){
    in.mark[i,first[i]]<-rbinom(1,1,p.juv*ch.true[i,first[i]])
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
  } else {}
  #############################
  ######################################################
  # Create population survey data
  ######################################################
  SUR <- rep(0,ti)
  for (i in 1:nd[2]){
    for (u in 1:ti){
      if(!is.na(IND_Count[1,u,i])){
        y3 <- rbinom(1, 1, PSUR[u])
        if(y3==1){
          SUR[u] <- SUR[u]+1
        }
      }
      if(!is.na(IND_Count[2,u,i])){
        y3 <- rbinom(1, 1, PSUR[u])
        if(y3==1){
          SUR[u] <- SUR[u]+1
        }
      }
    }
  }
  
  
  ###########################
  ##########################
  #ABBY NEST MODEL could go here
  
  #The individuals that are separated from the other data can be found in:
  #IND_Nest
  
  #I am using this just as a placeholder to make sure everything works
  R <- rep(0, ti) # number of pairs whose productivity was observed 
  nestlings <- rep(0, ti)  # total number of nestlings recoded in a year
  for (i in 1:nd[3]){
    for (v in 1:ti){
      if(!is.na(IND_Nest[4,v,i])){
        y <- rbinom(1, 1, PPROD[v])
        if(y==1){
          R[v] <- R[v]+1
          nestlings[v] <- nestlings[v]+IND_Nest[4,v,i]
        }
      }
    }
  }
  
  
  #Get all the data together, just MR and Survey Counts until nest in here
  return(list(ch=ch, SUR=SUR, age_ch=age_ch, first=first, last=last,
              R=R, nestlings=nestlings*2))
}
df<-simIPMdata(n.years=10, n.data=c(50,50,50), init.age = c(100,100), 
               phi.1=0.7, phi.ad=0.76,p.1=0.98, p.ad=0.65,
               f.1=0.8,f.ad=0.8, p.sur=0.8, p.prod=0.77)  
