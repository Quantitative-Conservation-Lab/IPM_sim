#code adapted from WinBugsBPA Appendix 2, code to simulate IPM
# link to code is here: https://www.vogelwarte.ch/assets/files/publications/BPA/Web%20appendix%202.txt


#function for simulation of data - just using the WinBugs Code 
# from the web appendix
#their IPM is basicially the same basic structure
#HAS: we will have to change some aspects, but would be good to go through 
#as a group to do so?
#for example, does this work for our Mark-Resight model, do we not want to do m-array

#Abby to add in how she simulates Nest Success data
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
  rownames(IND) <- c("1-Year", "Adu", "Dead", "Rep", "Juv")
  
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
  INDC <- IND[,,resamp1]
  INDB <- IND[,,resamp2]
  INDR <- IND[,,resamp3]
  
  
  #################################
  # Create Capture-Recapture Data
  #################################
  
  # 1. Construct capture histories
  #year 1 - go out and band all chicks 
  #year 2+ go out and band chick, resight others that are banded
  
  CR <- matrix(data = rep(0, nd[1]*(ti+3)), ncol = ti+3, nrow = nd[1])
  for ( i in 1:nd[1]){
    for (r in 1:ti){
      # Juveniles
      if(!is.na(INDC[5,r,i])){
        y1 <- rbinom(1, 1, PJUV[r])
        if(y1==1){
          CR[i,r] <- 1
          CR[i,ti+1] <- 1
        }
      }
      # 1-year old
      if(!is.na(INDC[1,r,i])){
        y2 <- rbinom(1, 1, PADULTS[r])
        if(y2==1){CR[i,r] <- 2
        if(CR[i,ti+1]==0) CR[i,ti+2] <- 1
        }
      }
      # Adults
      if(!is.na(INDC[2,r,i])){
        y3 <- rbinom(1, 1, PADULTS[r])
        if(y3==1){CR[i,r] <- 3
        if(CR[i,ti+1]==0&CR[i,ti+2]==0) CR[i,ti+3] <- 1
        }
      }
    }
  }
  colnames(CR)<-c(seq(1,ti), "Juveniles", "1-year", "Adults")
  
  
  # 2. Create m-array for the capture-histories
  ####NOT M-ARRAY, change****
  ####
  # 2.1. Two functions for data manipulation:
  # 2.1.1. Function to remove individuals without histories
  clean <- function(crma){
    su <- rep(NA, nrow(crma))
    for (i in 1:nrow(crma)){
      su[i] <- sum(crma[i,])}
    final <- crma[which(su!=0),]
  }
  
  # 2.1.2. Function to remove all first captures
  rmfirst <- function(crma){
    newma <- matrix(data = 0, ncol = ncol(crma), nrow = nrow(crma))
    for (i in 1:nrow(crma)){
      if (sum(crma[i,])==1) next  # if the individual is captured only once, move to next individual
      k <- min(which(crma[i,]==1))
      newma[i,] <- crma[i,]
      newma[i,k] <- 0
    }
    final <- clean(newma)
    if (is.null(ncol(final))==TRUE) {final <- matrix(final, ncol = length(final), nrow = 1)}
    final <- final
  }
  
  # 2.2. Recode capture-recapture matrix
  cr2 <- replace(CR,which(CR[,]==2),1)
  cr3 <- replace(cr2,which(cr2[,]==3),1)
  cr <- clean(cr3)
  m.juv <- matrix(data = 0, ncol = ti, nrow = ti)  # juv m-array
  m.adu <- matrix(data = 0, ncol = ti, nrow = ti)  # 1y + adults m-array
  cr.juv <- cr[which(cr[,(ti+1)]==1),1:ti]         # indivs captured as juv
  cr.adu <- cr[c(which(cr[,(ti+2)]==1), which(cr[,(ti+3)]==1)), 1:ti]  # indivs captured as 1y + adults
  
  # 2.3. Calculate m-array for juveniles
  for (tt in 1:(ti-1)){
    cr.h <- matrix(data = 0, ncol = ti, nrow = nrow(cr.juv))
    if(sum(cr.juv[,tt])==0)next
    for (i in 1:nrow(cr.juv)){
      if (cr.juv[i,tt]==1){
        cr.h[i,] <- cr.juv[i,]
        cr.juv[i,] <- 0
      }
    }
    
    m.juv[tt,1] <- sum(cr.h[,tt])
    
    # Remove all first captures in cr.h
    if(sum(cr.h)==0) next
    # cr.h1<-rmfirst(clean(cr.h))
    cr.h1 <- clean(cr.h)
    cr.h1 <- matrix(cr.h1, ncol = (ti))
    cr.h1 <- rmfirst(cr.h1)
    if(sum(cr.h1)==0) next
    # Determine when the birds were first recaptured
    po<-rep(0, nrow(cr.h1))
    for (j in 1:nrow(cr.h1)){
      po[j] <- min(which(cr.h1[j,]==1))
    }
    
    k <- as.data.frame(table(po))    #  a table of time first recaptured (po) by #indivs (freq)
    if(length(k$Freq)>=1){
      for (i in 1:length(k$Freq)){
        m.juv[tt,as.numeric(as.vector(k[i,1]))] <- k[i,2]   #give the m-array of juveniles
      }
    }
    cr.adu <- rbind(cr.adu, cr.h1[which(po>=(tt+1)),])  #  merge adults cap hist and juvs recaptured as adults(1y & adults) following first release
  }
  
  # 2.4. m-array for adults
  for (tt in 1:(ti-1)){
    cr.h <- matrix(data = 0, ncol = ti, nrow = nrow(cr.adu))
    if(sum(cr.adu[,tt])==0) next
    for (i in 1:nrow(cr.adu)){
      if (cr.adu[i,tt]==1){
        cr.h[i,] <- cr.adu[i,]
        cr.adu[i,] <- 0
      }
    }
    m.adu[tt,1] <- sum(cr.h[,tt])
    # Remove all first captures in cr.h
    if(sum(cr.h)==0) next
    # cr.h1 <- rmfirst(clean(cr.h))
    cr.h1 <- clean(cr.h)
    cr.h1 <- matrix(cr.h1, ncol = (ti))
    cr.h1 <- rmfirst(cr.h1)
    if(sum(cr.h1)==0) next
    # Determine when the birds were first recaptured
    po <- rep(0,nrow(cr.h1))
    for (j in 1:nrow(cr.h1)){
      po[j] <- min(which(cr.h1[j,]==1))
    }
    k <- as.data.frame(table(po))
    if(length(k$Freq)>=1){
      for (i in 1:length(k$Freq)){
        m.adu[tt,as.numeric(as.vector(k[i,1]))] <- k[i,2]
      }
    }
    cr.adu <- rbind(cr.adu, cr.h1)
  }
  
  # 2.5. Rearrange the m-array
  m <- matrix(NA, ncol = ti, nrow = 2*(ti-1))
  # Juv
  for (i in 1:(ti-1)){
    for (j in 1:(ti-1)){
      m[i,j] <- m.juv[i,j+1]
      m[i,ti] <- m.juv[i,1]-sum(m.juv[i,2:ti])
    }
  }
  
  # Adults
  for (i in 1:(ti-1)){
    for (j in 1:(ti-1)){
      m[i+ti-1,j] <- m.adu[i,j+1]
      m[i+ti-1,ti] <- m.adu[i,1]-sum(m.adu[i,2:ti])
    }
  }
  
  ######################################################
  # Create population survey data
  ######################################################
  SUR <- rep(0,ti)
  for (i in 1:nd[2]){
    for (u in 1:ti){
      if(!is.na(INDB[1,u,i])){
        y3 <- rbinom(1, 1, PSUR[u])
        if(y3==1){
          SUR[u] <- SUR[u]+1
        }
      }
      if(!is.na(INDB[2,u,i])){
        y3 <- rbinom(1, 1, PSUR[u])
        if(y3==1){
          SUR[u] <- SUR[u]+1
        }
      }
    }
  }
  
  
  #############################################
  # Create reproductive success data
  #############################################
  
  R <- rep(0, ti) # number of pairs whose productivity was observed 
  nestlings <- rep(0, ti)  # total number of nestlings recoded in a year
  for (i in 1:nd[3]){
    for (v in 1:ti){
      if(!is.na(INDR[4,v,i])){
        y <- rbinom(1, 1, PPROD[v])
        if(y==1){
          R[v] <- R[v]+1
          nestlings[v] <- nestlings[v]+INDR[4,v,i]
        }
      }
    }
  }
  
  #######################
  # Bundle data
  #######################
  
  return(list(m = m, SUR = SUR, nestlings = nestlings*2, R = R))
  
}

df<-simIPMdata(n.years=10, n.data=c(50,50,50), init.age = c(100,100), 
               phi.1=0.7, phi.ad=0.76,p.1=0.98, p.ad=0.65,
               f.1=0.8,f.ad=0.8, p.sur=0.8, p.prod=0.77)
df
