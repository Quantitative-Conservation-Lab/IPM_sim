#HAS - need to comment throughout, but got the MR data (I think) how it should be


#altered IPM sim with placeholder for nest success
#population projections and count code is from BPA WinBUGS appendix 2
#code adapted from WinBugsBPA Appendix 2, code to simulate IPM
# link to code is here: https://www.vogelwarte.ch/assets/files/publications/BPA/Web%20appendix%202.txt

library(nimble)

n.years=10; n.data=c(50,50,100);init.age = c(100,100); 
phi.1=0.4; phi.ad=0.76;p.1=0.98; p.ad=0.65;
p.prod=0.77; p.sur=0.8

n.initiation.dates <- 31

max.nest.age <- 30
first.initiation.date <- 1
last.fledge.date <- n.initiation.dates + max.nest.age
season.length <- last.fledge.date - first.initiation.date + 1 + 2

# mean clutch size
mean.clutch.size <- 3

# daily nest survival
phi.nest <- 0.975

prop.nests.found <- 0.8

n.sam<-3

visit.interval <- 3

simIPMdata<-function(n.years, n.data, init.age, phi.1, phi.ad, p.1,p.ad,p.sur,
 n.initiation.dates, max.nest.age, first.initiation.date, last.fledge.date, 
 season.length, mean.clutch.size, phi.nest, prop.nests.found, visit.interval, n.sam){
  
  ti<-n.years
  ni<-n.years-1
  nd<-n.data
  Ni<-init.age
  phi1<-phi.1
  phi2<-phi.ad
  pjuv<-p.1
  psur<-p.sur
  padults<-p.ad
  n.initiation.dates <- n.initiation.dates
  max.nest.age <- max.nest.age
  first.initiation.date <- first.initiation.date
  last.fledge.date <- n.initiation.dates + max.nest.age
  season.length <- last.fledge.date - first.initiation.date + 1 + 2
  mean.clutch.size <- mean.clutch.size
  phi.nest <- phi.nest
  prop.nests.found <- prop.nests.found
  visit.interval <- visit.interval
  n.sam <- n.sam
  
  # DERIVE TRUE FECUNDITY
  true.fec <- 1/2 * mean.clutch.size * phi.nest^max.nest.age
  
  # Combine the values to vectors
  PHI1 <- rep(phi1, ti)
  PHI2 <- rep(phi2, ti)
  PJUV <- rep(pjuv, ti)
  PADULTS <- rep(padults, ti)
  PSUR <- rep(psur, ti)
  FL1 <- rep(true.fec, ti+1)
  FL2 <- rep(true.fec, ti+1)
  
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
  
  # LINK YOY AND ADULTS TO DAILY NEST SUCCESS
  nests <- array(data = NA, dim = c(2, ti+1, Ntotal, season.length))
  for (age in 1:2) { # age classes
    for (year in 1:(ti+1)) { # year
      for (i in 1:sum(IND[age, year, ], na.rm = TRUE)) { # every alive individual
        
        # assume everybody attempts a nest
        N.nests.total <- sum(IND[age, year, ], na.rm = TRUE)
        
        inds.nesting <- which(!is.na(IND[age, year, ]))
        
        total.nests.age <- matrix(NA, nrow = Ntotal, ncol = season.length)
        total.nests.status <- matrix(NA, nrow = Ntotal, ncol = season.length)
        
        # initiation date
        # uniform across season length
        # sort by initiation date
        init.dates <- sort(rcat(N.nests.total, rep(1/n.initiation.dates, length.out = n.initiation.dates)))
        
        for (i in 1:N.nests.total) {
          # true age of nest 
          total.nests.age[inds.nesting[i], init.dates[i]:(init.dates[i] + max.nest.age - 1)] <- 1:max.nest.age
          
          # true status of nest
          total.nests.status[inds.nesting[i], init.dates[i]] <- 1 
          for (a in 1:(max.nest.age-1)) {
            total.nests.status[inds.nesting[i], init.dates[i] + a] <- rbinom(1, 1, total.nests.status[inds.nesting[i], init.dates[i] + a - 1] * phi.nest)
          }
        }
      }
      
      nests[age, year, inds.nesting, ] <- total.nests.status[inds.nesting, ]
    }
  }
  
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
  ind_mr<-IND_MR[c(5,1,2),1:ti,]
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
  # for (i in 1:nd[2]){
  #   for (u in 1:ti){
  #     TRUE_Count[1,u]<-sum(IND_Count[1,u,], na.rm = T)
  #     TRUE_Count[2,u]<-sum(IND_Count[2,u,], na.rm = T)
  #     if(!is.na(IND_Count[1,u,i])){
  #       y3 <- rbinom(1, 1, PSUR[u])
  #       if(y3==1){
  #         SUR[u] <- SUR[u]+1
  #       }
  #     }
  #     if(!is.na(IND_Count[2,u,i])){
  #       y3 <- rbinom(1, 1, PSUR[u])
  #       if(y3==1){
  #         SUR[u] <- SUR[u]+1
  #       }
  #     }
  #   }
  # }
  
  
  ###########################
  ##########################
  #ABBY NEST MODEL could go here
  
  inds.nests.available <- resamp3
  
  n.found <- length(inds.nests.available) #* prop.nests.found
  
  H <- array(data = NA, dim = c(ti+1, n.found, season.length))
  Hlatent <- array(data = NA, dim = c(ti+1, n.found, season.length))
  which.found <- matrix(data = NA, nrow = ti+1, ncol = n.found)  
  
  Fledged <- matrix(data = NA, nrow = ti+1, ncol = n.found)  
  
  is.nest <- matrix(data = NA, nrow = ti+1, ncol = n.found)
  last.status <- matrix(data = NA, nrow = ti+1, ncol = n.found)
  first.nest <- matrix(data = NA, nrow = ti+1, ncol = n.found)
  last.nest <- matrix(data = NA, nrow = ti+1, ncol = n.found)
  for (year in 1:(ti+1)) { # year
    # which of these nests were found - this is also latent
    which.found[year, ] <- sort(sample(resamp3, n.found))
    init.dates <- numeric(n.found)
    age.when.found <- numeric(n.found)
    ages.visited <- list(NULL)
    for (i in 1:n.found) {
      if (!is.na(IND[1, year, which.found[year, ][i]])) {
        if (IND[1, year, which.found[year, ][i]] > 0) {
          Hlatent[year, i, ] <- nests[1, year, which.found[year, ][i], ]
          init.dates[i] <- min(which(!is.na(Hlatent[year, i, ])))
          age.when.found[i] <- rcat(1, rep(1/max.nest.age, length.out = max.nest.age))
          ages.visited[[i]] <- unique(c(seq(age.when.found[i], max.nest.age, by = visit.interval), max.nest.age))
          
          H[year, i, init.dates[i]+ages.visited[[i]]-1] <- Hlatent[year, i, init.dates[i]+ages.visited[[i]]-1]
          
          is.nest[year, i] <- H[year, i, min(which(!is.na(H[year, i, ])))] == 1
          last.status[year, i] <- H[year, i, max(which(!is.na(H[year, i, ])))]
          
          if (last.status[year, i]) {
            Fledged[year, i] <- IND[4, year, which.found[year, ][i]]
          }
          
          first.nest[year, i] <- min(which(!is.na(H[year, i, ])))
          if (is.nest[year, i] & last.status[year, i]) { # found active and succeeded
            last.nest[year, i] <- max(which(!is.na(H[year, i, ])))
          } else if (is.nest[year, i] & !last.status[year, i]) { # found active but then failed
            last.nest[year, i] <- min(which(!is.na(H[year, i, ]) & H[year, i, ] == 0))
          } else { # found inactive
            last.nest[year, i] <- first.nest[year, i]
          }
            
        } else if (IND[2, year, which.found[year, ][i]] > 0) {
          Hlatent[year, i, ] <- nests[2, year, which.found[year, ][i], ]
          init.dates[i] <- min(which(!is.na(Hlatent[year, i, ])))
          age.when.found[i] <- rcat(1, rep(1/max.nest.age, length.out = max.nest.age))
          ages.visited[[i]] <- unique(c(seq(age.when.found[i], max.nest.age, by = visit.interval), max.nest.age))
          
          H[year, i, init.dates[i]+ages.visited[[i]]-1] <- Hlatent[year, i, init.dates[i]+ages.visited[[i]]-1]
          
          is.nest[year, i] <- H[year, i, min(which(!is.na(H[year, i, ])))] == 1
          last.status[year, i] <- H[year, i, max(which(!is.na(H[year, i, ])))]
          
          if (last.status[year, i]) {
            Fledged[year, i] <- IND[4, year, which.found[year, ][i]]
          }
          
          first.nest[year, i] <- min(which(!is.na(H[year, i, ])))
          if (is.nest[year, i] & last.status[year, i]) { # found active and succeeded
            last.nest[year, i] <- max(which(!is.na(H[year, i, ])))
          } else if (is.nest[year, i] & !last.status[year, i]) { # found active but then failed
            last.nest[year, i] <- min(which(!is.na(H[year, i, ]) & H[year, i, ] == 0))
          } else { # found inactive
            last.nest[year, i] <- first.nest[year, i]
          }
          
        } # else if
      } # if not NA
    } # found
  } # years
  
  N.nests.found <- rowSums(is.nest, na.rm = TRUE)
  N.nests.successful <- rowSums(last.status, na.rm = TRUE)
  
  first.nesttidy <- matrix(data = NA, nrow = ti+1, max(N.nests.found)) 
  
  last.nesttidy <- matrix(data = NA, nrow = ti+1, max(N.nests.found)) 
  
  Htidy <- array(data = NA, dim = c(ti+1, max(N.nests.found), season.length))
  Fledgedtidy <- matrix(data = NA, nrow = ti+1, max(N.nests.successful)) 
  for (year in 1:(ti+1)) { # year
    
    is.nest.this.year <- which(is.nest[year, ])
    
    if (N.nests.found[year] > 0) {
      for (n in 1:N.nests.found[year]) {
        Htidy[year, n, ] <- H[year, is.nest.this.year[n], ]
        first.nesttidy[year, n] <- first.nest[year, is.nest.this.year[n]]
        last.nesttidy[year, n] <- last.nest[year, is.nest.this.year[n]]
      }  
    }
    
    if (!all(is.na(Fledged[year, ]))) {
      Fledgedtidy[year, 1:N.nests.successful[year]] <- na.omit(Fledged[year, ])  
    }
  }
  
  #Get all the data together, just MR and Survey Counts until nest in here
  return(list(ch=ch, SUR=SUR, age_ch=age_ch, first=first, last=last,
              H=Htidy, Fledged=Fledgedtidy, 
              first.nest = first.nesttidy, last.nest = last.nesttidy, 
              max.nest.age = max.nest.age, 
              n.nests = N.nests.found, n.succ.nests = N.nests.successful, 
              n.sam = n.sam))
}
df<-simIPMdata(n.years, n.data, init.age, phi.1, phi.ad, p.1,p.ad,p.sur,
               n.initiation.dates, max.nest.age, first.initiation.date, last.fledge.date, 
               season.length, mean.clutch.size, phi.nest, prop.nests.found, visit.interval, n.sam)  
