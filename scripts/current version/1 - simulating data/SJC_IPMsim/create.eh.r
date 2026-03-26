#####################################################################################################
#
# Function to generate capture histories from the population (ind) 
#
# Input variables
#    ind: array with the population
#    c: matrix with age- and time-specific capture probabilities (probability of first capture)
#    p: matrix with age- and time-specific REcapture probabilities
#    maxAge: maximal number of age classes that can be identified when the individuals are captured for the first time
#
# Output
#    ch: matrix with the capture histories
#    age: vector with the age class at first capture for each individual
#
# Written: 11.3.2016, M.Schaub
#
# Last up-date: 19.5.2017
#
#####################################################################################################


#####################################################################################################
#
# Function to generate capture histories from the population (ind) 
#
# Input variables
#    ind: array with the population
#    c: matrix with age- and time-specific capture probabilities (probability of first capture)
#    p: matrix with age- and time-specific REcapture probabilities
#    maxAge: maximal number of age classes that can be identified when the individuals are captured for the first time
#
# Output
#    ch: matrix with the capture histories
#    age: vector with the age class at first capture for each individual
#
# Written: 11.3.2016, M.Schaub
#
# Last up-date: 19.5.2017
#
#####################################################################################################


create.capturehistory <- function(ind, c, p, maxAge = 2, seed = NA){
  
  if (!is.na(seed)) {set.seed(seed)}
  T <- dim(ind)[2]
  nind <- dim(ind)[3]
  nstage <- dim(ind)[1]
  aclasses <- nstage-3
  age <- first <- last <- numeric()
  
  for (i in 1:nind){
    g <- which(!is.na(ind[1:(aclasses+1),,i]), arr.ind = TRUE)
    age[i] <- g[1,1]
    first[i] <- g[1,2]
    h <- which(ind[1:(aclasses+1),,i]==1, arr.ind = TRUE)
    last[i] <- max(h[,2])
  } # i
  
  ch.true <- ch <- in.cap <- matrix(0, ncol = T, nrow = nind)
  for (i in 1:nind){
    ch.true[i,first[i]:last[i]] <- 1
  } # i
  # Recode age
  age[age > maxAge] <- maxAge
  
  # Sampling
  # Expand c and p (to higher age classes)
  C <- matrix(0, ncol = T, nrow = max(c(maxAge, nrow(c))) + T)
  C[1:nrow(c),] <- c
  u <- max(c(maxAge, nrow(c))) + T - nrow(c)
  if (u > 0){
    for (j in 1:u){
      C[nrow(c)+j,] <- c[nrow(c),]
    } # j
  } # if
  
  P <- matrix(0, ncol = T-1, nrow = max(c(maxAge, nrow(p))) + T)
  P[1:nrow(p),] <- p
  u <- max(c(maxAge, nrow(p))) + T - nrow(p)
  if (u > 0){
    for (j in 1:u){
      P[nrow(p)+j,] <- p[nrow(p),]
    } # j
  } # if
  
  for (i in 1:nind){
    # First capture
    # When captured for the first time?
    for (t in first[i]:last[i]){
      in.cap[i,t] <- rbinom(1, 1, C[age[i]+t-first[i],t])
    }
    ch[i,] <- in.cap[i,]
    # up-date the age at first capture
    if (sum(in.cap[i,]) >= 1){
      k <- min(which(in.cap[i,] == 1))
      age[i] <- age[i] + k - first[i]
      if (age[i] > maxAge) {age[i] <- maxAge}
    }  # if
    if (sum(in.cap[i,]) > 1){
      h <- which(in.cap[i,] == 1)
      ch[i,h[2:length(h)]] <- 0
    }
    if (first[i]==last[i]) next
    
    if (sum(in.cap[i,])!=0){
      # Recapture (conditional on first capture)
      cap.occ <- min(which(in.cap[i,] == 1))
      if (cap.occ==last[i]) next
      for (t in (cap.occ+1):last[i]){
        ch[i,t] <- rbinom(1, 1, P[age[i]+t-first[i],t-1])
      } # t
    } # if
  } # i
  
  # Compute the size of the population to check whether the simulated population went extinct
  Ptot <- numeric()
  for (t in 1:T){
    Ptot[t] <- sum(ind[2:(aclasses+1),t,], na.rm = TRUE)
  }
  
  # Remove individuals that have never been captured/marked and occasions after population extinction
  incl <- which(rowSums(ch)>=1)
  ch <- ch[incl,]
  age <- age[incl]
  incl <- which(Ptot>0)
  ch <- ch[,incl]
  if (length(incl)<T) print (paste("Population extinct after ",length(incl)," years"))
  
  return(list(ch = ch, age = age))  
}




#####################################################################################################
#
# Function to create age-dependent m-arrays 
#
# Input variables
#    ch: matrix with capture histories. Note, this is a single file including all age classes
#    age: vector with the age for each individual at first capture
#    mAge: maximal number of age classes for which m-arrays are constructed. Input is optional and only required if the age matrix has fewer age classes as we want to separate (e.g. CH contains only individuals marked as juveniles, and we want 2 age classes)
#
# Output
#    marr: 3-d array with the m-array. The third dimension is the age class. The last column of each m-array is the number of released individuals that were never recaptured. Thus, the total number of released individuals per occasion is the row sum of each m-array.
#
# Written: 14.3.2016, M.Schaub
#
# Last up-date:
#
#####################################################################################################

marray.age <- function(ch, age, mAge = 1){
  
  # 1. Helper functions
  # 1.1. Function to create a m-array based on capture-histories (ch)
  marray <- function(ch){
    nind <- nrow(ch)
    n.occasions <- ncol(ch)
    m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
    # Calculate the number of released individuals at each time period
    m.array[,1] <- colSums(ch)
    for (i in 1:nind){
      pos <- which(ch[i,]==1)
      g <- length(pos)
      if (g==1) next
      for (z in 1:(g-1)){
        m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
      } # z
    } # i
    # Calculate the number of individuals never recaptured
    for (t in 1:n.occasions){
      m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
    } # t
    out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
    return(out)
  }
  
  # 1.2. Function to remove histories without any capture from a capture-recapture matrix
  clean.ch <- function(ch){
    incl <- which(rowSums(ch)>=1)
    ch <- ch[incl,]
    return(ch)
  }
  
  # 1.3. Function to remove the first capture in a capture-recapture matrix
  rm.first <- function(ch) {
    get.first <- function(x) min(which(x==1))
    first <- apply(ch, 1, get.first)
    for (i in 1:nrow(ch)){
      ch[i,first[i]] <- 0
    }
    return(ch)
  }
  
  # 1.4. Function to calculate the occasion of first capture
  get.first <- function(x) min(which(x==1))
  
  
  # 2. Calculations   
  if (is.matrix(ch)==FALSE) ch <- matrix(ch, nrow = 1)   
  maxAge <- max(c(max(age), mAge))
  nind <- nrow(ch)
  n.occasions <- ncol(ch)
  
  first <- apply(ch, 1, get.first)
  age.matrix <- matrix(0, ncol = n.occasions, nrow = nind)
  for (i in 1:nind){
    age.matrix[i,first[i]:n.occasions] <- 1:(n.occasions-first[i]+1)+(age[i]-1)
  }
  age.matrix[age.matrix > maxAge] <- maxAge
  
  # Recode capture history
  ch.rec <- ch
  for (i in 1:nind){
    h <- which(ch.rec[i,]==1)
    for (j in 1:length(h)){
      ch.rec[i,h[j]] <- j
    } # j
  } # i
  ch.rec[ch.rec > maxAge] <- maxAge
  
  ch.split <- array(0, dim = c(nrow(ch), ncol(ch), maxAge))
  for (a in 1:maxAge){
    for (i in 1:nind){
      j <- which(ch.rec[i,]==a | ch.rec[i,]==(a+1))
      if (length(j)==0) next
      ch.split[i,j[1:2],age.matrix[i,j[1]]] <- 1
      if (length(j)>1){
        ch.split[i,j[2:length(j)],age.matrix[i,j[2]]] <- 1
      }
    } # i
  } # a
  
  marr <- array(0, dim = c(n.occasions-1, n.occasions, maxAge))
  for (a in 1:(maxAge-1)){
    for (i in 1:nind){
      u <- which(ch.split[i,,a]==1)
      if (length(u)==0) next
      if (u[1]==n.occasions) next
      if (length(u)==1) marr[u,n.occasions,a] <- marr[u,n.occasions,a] + 1
      if (length(u)==2) marr[u[1],u[2]-1,a] <- marr[u[1],u[2]-1,a] + 1
    } # i
  } # a
  a <- maxAge
  
  if (is.matrix(ch.split[,,a])==FALSE){ 
    ch.split1 <- matrix(ch.split[,,a], nrow = 1)
    marr[,,a] <- marray(ch.split1)
  } # if
  else marr[,,a] <- marray(ch.split[,,a])      
  return(marr)
}


