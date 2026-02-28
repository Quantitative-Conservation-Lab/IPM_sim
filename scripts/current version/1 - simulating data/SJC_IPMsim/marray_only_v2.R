########################################################
#libraries 
########################################################

library(nimble)

########################################################
#functions 
########################################################

simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol=n.occasions, nrow=sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1              # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break                  # If dead, move to next individual 
      # Bernoulli trial: is individual recaptured? 
      rp <- rbinom(1, 1, P[i,t-1])
      if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
}

marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data=0, ncol=n.occasions+1, nrow=n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

########################################################
#simulate data   
########################################################

sims <- 5
results <- array(NA,dim = c(3,5,sims))

for(s in 1:sims){
  
# Define parameter values
n.occasions <- 12                        # Number of capture occasions
marked.j <- rep(200, n.occasions-1)      # Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1)       # Annual number of newly marked adults
phi.juv <- 0.3                           # Juvenile annual survival
phi.ad <- 0.65                           # Adult annual survival
p <- rep(0.5, n.occasions-1)             # Recapture
phi.j <- c(phi.juv, rep(phi.ad,n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol=n.occasions-1, nrow=sum(marked.j))
for (i in 1:(length(marked.j)-1)){
  PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)],marked.j[i]), ncol=n.occasions-i, byrow=TRUE)
}
P.J <- matrix(rep(p, sum(marked.j)), ncol=n.occasions-1, nrow=sum(marked.j), byrow=TRUE)
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol=n.occasions-1, nrow=sum(marked.a), byrow=TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol=n.occasions-1, nrow=sum(marked.a), byrow=TRUE)

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a) 

########################################################
#create m-arrays  
########################################################

cap <- apply(CH.J, 1, sum)
ind <- which(cap >= 2)
# Juvenile CH recaptured at least once
CH.J.R <- CH.J[ind,] 
# Juvenile CH never recaptured
CH.J.N <- CH.J[-ind,]                    
# Remove first capture
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}
# Add grown-up juveniles to adults and create m-array
CH.A.m <- rbind(CH.A, CH.J.R1)
CH.A.marray <- marray(CH.A.m)
# Create CH matrix for juveniles, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
  second[i] <- min(which(CH.J.R1[i,]==1))
}
CH.J.R2 <- matrix(0, nrow=dim(CH.J.R)[1], ncol=dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R2[i,first[i]] <- 1
  CH.J.R2[i,second[i]] <- 1
}
# Create m-array for these
CH.J.R.marray <- marray(CH.J.R2)
# The last column ought to show the number of juveniles not recaptured again and should all be zeros, since all of them are released as adults
CH.J.R.marray[,dim(CH.J)[2]] <- 0
# Create the m-array for juveniles never recaptured and add it to the previous m-array
CH.J.N.marray <- marray(CH.J.N)
CH.J.marray <- CH.J.R.marray + CH.J.N.marray 


########################################################
#model code   
########################################################

cjs12Code <- nimbleCode({
  
  # Priors and constraints
  for (t in 1:(n.occasions-1)){
    phi.juv[t] <- mean.phijuv
    phi.ad[t] <- mean.phiad
    p[t] <- mean.p
  }
  mean.phijuv ~ dunif(0, 1)                # Prior for mean juv. survival
  mean.phiad ~ dunif(0, 1)                 # Prior for mean ad. survival
  mean.p ~ dunif(0, 1)                     # Prior for mean recapture
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,1:n.occasions], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,1:n.occasions], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  # Main diagonal
  for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]                         # Probability of non-recapture
    pr.j[t,t] <- phi.juv[t]*p[t]
    pr.a[t,t] <- phi.ad[t]*p[t]
  }
  # Above main diagonal
  for (t in 1:(n.occasions-2)){
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- phi.juv[t]*prod(phi.ad[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j] <- prod(phi.ad[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
  } #t
  # Below main diagonal
  for (t in 2:(n.occasions-1)){
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }
})

########################################################
#nimble inputs   
########################################################

# Bundle data and constants
dataList <- list(marr.j=CH.J.marray, marr.a=CH.A.marray, rel.j=rowSums(CH.J.marray), rel.a=rowSums(CH.A.marray))
constList <- list(n.occasions=dim(CH.J.marray)[2])  ;  str(constList)

# Initial values
inits <- function(){list(mean.phijuv=runif(1, 0, 1), mean.phiad=runif(1, 0, 1), mean.p=runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phijuv", "mean.phiad", "mean.p")

# MCMC settings
ni <- 10000  ;  nt <- 1  ;  nb <- 2000  ;  nc <- 4  # serious

# Call NIMBLE from R (ART <1 min)

  out12 <- nimbleMCMC(
    code=cjs12Code, 
    data=dataList, 
    constants=constList, 
    inits=inits(), 
    monitors=parameters, 
    niter=ni, 
    nburnin=nb, 
    nchains=nc, 
    thin=nt,
    summary = TRUE,
    samplesAsCodaMCMC=TRUE) 

########################################################
#nimble output 
########################################################

#check diagnostics
MCMCtrace(object = out12$samples,
          pdf = FALSE,
          ind = TRUE)

#summary 
out12$summary$all.chains 

#Gelman-Rubin diagnostic
gelman.diag(out12$samples)


results[,,s] <- out12$summary$all.chains

}
dimnames(results)[[1]] <- c("mean.p","mean.phiad","mean.phijuv")
dimnames(results)[[2]] <- c("Mean","Median","SD","95%Low","95%Upp")
