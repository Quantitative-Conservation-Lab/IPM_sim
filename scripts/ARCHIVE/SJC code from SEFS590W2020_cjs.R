#original code from SEFS 590 winter 2020 class for reference
#########

#SEFS 590 - University of Washington
#Week 5 - Cormack-Jolly-Seber Model
#Simulate a CJS data set
#fit in optim and using 3 different methods in JAGS 
#(multinomial on encounter histories, state-space, multinomial on m-array) 
#Code updated 1 February 2020
#Code written by Sarah Converse, sconver@uw.edu 
#CJS simulation code and m-array analysis code adapted from Kery and Schaub

#############################################################################


#############################################################################

#Simulate CJS data 
n.occ <- 3 #number of occasions
marked <- rep(50,n.occ) #number of inds marked per occasion
#matrix of values for phi and p (occasion - 1)*inds
PHI <- matrix(c(0.4,0.4), ncol = n.occ-1, nrow = sum(marked),byrow=TRUE)
P <- matrix(c(0.7,0.7), ncol = n.occ-1, nrow = sum(marked),byrow=TRUE)

#simulate CJS data 
CJS.sim <- function(marked,P,PHI){
  n.occ <- dim(PHI)[2]+1 #number of occasions
  EH <- matrix(0, ncol = n.occ, nrow = sum(marked)) #assign matrix
  first <- rep(1:length(marked),marked[1:length(marked)]) #get first encounter occasion for each ind
  for(i in 1:sum(marked)){ #loop over inds
    EH[i,first[i]] <- 1 #condition on first capture
    if(first[i]==n.occ) next #if first capture occasion for ind is last occasion, next
    for(t in (first[i]+1):n.occ){ #loop over time 
      surv <- rbinom(1,1,PHI[i,t-1]) #Bernoulli survival trial
      if(surv == 0) break #if didn't survive, done
      cap <- rbinom(1, 1, P[i,t-1]) #if survived, Bernoulli capture trial
      if (cap==1) EH[i,t] <- 1 #fill in the EH
    }
  }
  return(EH)
}
#run the simulation to get encounter histories 
EH <- CJS.sim(marked,P,PHI)

#############################################################################


#############################################################################

#Fit in optim using a multinomial on encounter histories
#this only works for 3 encounter occasions!

#these are the possible encounter histories, for which we know probabilities
poss.EH <- matrix(NA,nrow=7,ncol=3)
poss.EH[1,] <- c(1,1,1)
poss.EH[2,] <- c(1,1,0)
poss.EH[3,] <- c(1,0,1)
poss.EH[4,] <- c(1,0,0)
poss.EH[5,] <- c(0,1,1)
poss.EH[6,] <- c(0,1,0)
poss.EH[7,] <- c(0,0,1)

#get frequencies for each of the possible EH
freq <- matrix(NA,nrow(EH),nrow(poss.EH))
for(i in 1:nrow(EH)){
  for(j in 1:nrow(poss.EH)){
    freq[i,j] <- prod(EH[i,]==poss.EH[j,])
  }
}
f <- colSums(freq)

#the negative log-likelihood 
Mhlik<-function(parms){
  #solve on the logit scale
  logit.phi <- parms[1] 
  logit.p <- parms[2]
  phi <- 1/(1+exp(-logit.phi))
  p <- 1/(1+exp(-logit.p))
  
  #the likelihood function
  lik <- ( 
    (phi^2*p^2)^f[1]* #111 EH
      (phi*p*(phi*(1-p)+(1-phi)))^f[2]* #110 EH
      (phi*(1-p)*phi*p)^f[3]* #101 EH
      (phi*(1-p)*((1-phi)+(phi*(1-p)))+(1-phi))^f[4]* #100 EH
      (phi*p)^f[5]* #011 EH
      (phi*(1-p)+(1-phi))^f[6] #010 EH
  )
  #negative log likelihood
  -1*log(lik)
}

#solve 
out.optim <- optim(c(0,0),Mhlik,hessian=TRUE,method="Nelder-Mead")
VC <- solve(out.optim$hessian) #get the VC matrix for logit scale parameters
phi <- 1/(1+exp(-out.optim$par[1])) #get real scale phi
p <- 1/(1+exp(-out.optim$par[2])) #get real scale p

#get 95% CI limits on the logit scale and transform to real scale
logit.int.phi <- rep(NA,2); logit.int.phi[1] <- out.optim$par[1]-1.96*sqrt(VC[1,1]); logit.int.phi[2] <- out.optim$par[1]+1.96*sqrt(VC[1,1])
real.int.phi <- 1/(1+exp(-logit.int.phi))
logit.int.p <- rep(NA,2); logit.int.p[1] <- out.optim$par[2]-1.96*sqrt(VC[2,2]); logit.int.p[2] <- out.optim$par[2]+1.96*sqrt(VC[2,2])
real.int.p <- 1/(1+exp(-logit.int.p))

#Delta Method: (derivative)^2*var to get the variance of variable on real scale
der.logit <- function(theta){#derivative of invlogit(theta) wrt theta 
  d <- exp(theta)/((1+exp(theta))^2)
  return(d)}
SE.phi <- sqrt( (der.logit(out.optim$par[1]))^2*VC[1,1]) #get the SE  
SE.p <- sqrt( (der.logit(out.optim$par[2]))^2*VC[2,2])  #get the SE

#summarize the results and print 
summary.optim <- matrix(c(phi,p,SE.phi,SE.p),nrow=2,ncol=2)
rownames(summary.optim) <- c("phi","p")
colnames(summary.optim) <- c("MLE","SE")

print(summary.optim)

#############################################################################


#############################################################################

#Fit in JAGS using a multinomial on encounter histories
#this only works for 3 encounter occasions!
#need to run the script above to summarize the frequencies first 

# Specify model in BUGS language
cat("
    
    model {
    
    #Likelihood 
    f.mat1[1:4] ~ dmulti(probs1[1:4],R[1])
    f.mat2[1:2] ~ dmulti(probs2[1:2],R[2]) 
    
    probs1[1] <- phi[1]*p[1]*phi[2]*p[2]
    probs1[2] <- phi[1]*p[1]*((1-phi[2]) + phi[2]*(1-p[2]))
    probs1[3] <- phi[1]*(1-p[1])*phi[2]*p[2]
    probs1[4] <- 1-probs1[1]-probs1[2]-probs1[3] #(1-phi[1]) + phi[1]*(1-p[1])*((1-phi[2]) + (phi[2]*(1-p[2])))
    
    probs2[1] <- phi[2]*p[2]
    probs2[2] <- 1-probs2[1] #(1-phi[2]) + phi[2]*(1-p[2])
    
    # Priors and constraints
    for(i in 1:2){
    phi[i] <- mean.phi
    p[i] <- mean.p
    }
    
    mean.phi ~ dunif(0,1)         # Priors for survival
    mean.p ~ dunif(0,1)           # Priors for recapture
    
    }
    ",file="CJS_mnl.txt")

f.mat1 <- f[1:4]
f.mat2 <- f[5:6]
R <- rep(NA,2); R <- c(sum(f.mat1),sum(f.mat2))

## JAGS input
data <- list(R=R,f.mat1=f.mat1,f.mat2=f.mat2)
params <- c("mean.phi","mean.p") 
inits =  function() {list(mean.phi=runif(1),mean.p=runif(1))} 

library(jagsUI)
## RUN
out.jags = jags(data, inits, params, model.file="CJS_mnl.txt",
                n.chains=3, n.iter=20000, n.burnin=5000, n.thin=1)
out.jags$summary


#############################################################################


#############################################################################

#Fit in JAGS using a state-space model 

#JAGS MODEL 
cat("
    model{
    
    #state-space likelihood 
    for(i in 1:n.ind){
    z[i,first[i]] <- 1  #condition on first capture 
    for(j in (first[i]+1):n.occ){
    z[i,j] ~ dbern(z[i,j-1]*phi[i,j-1]) #state process 
    
    y[i,j] ~ dbern(z[i,j]*p[i,j-1]) #observation process 
    }
    }
    
    #constraints and priors 
    for(i in 1:n.ind){
    for(j in 1:(n.occ-1)){
    phi[i,j] <- mean.phi 
    p[i,j] <- mean.p
    }
    }
    mean.phi ~ dunif(0,1)
    mean.p ~ dunif(0,1)  
    
    }
    ",file="CJS.txt")

## JAGS input
n.ind <- nrow(EH)
n.occ <- ncol(EH)
first <- rep(NA,n.ind)
for(i in 1:n.ind){ first[i] <- min(which(EH[i,]==1))}

#Provide known state as data  
state.data <- function(EH){
  state <- EH
  for (i in 1:dim(EH)[1]){
    first <- min(which(EH[i,]==1))
    last <- max(which(EH[i,]==1))
    state[i,first:last] <- 1
    state[i,first] <- NA
  }
  state[state==0] <- NA
  return(state)
}

z.state <- state.data(EH)
data <- list(n.ind=n.ind,n.occ=n.occ,y=EH,first=first,z=z.state)
params <- c("mean.phi","mean.p") 
inits =  function() {list(mean.phi=runif(1),mean.p=runif(1))} 

library(jagsUI)
## RUN
out.jags = jags(data, inits, params, model.file="CJS.txt",
                n.chains=3, n.iter=20000, n.burnin=5000, n.thin=1)
out.jags$summary


#############################################################################


#############################################################################

#Fit in JAGS using a multinomial on the m-array 

# Function to summarize an m-array based on capture-histories (EH)
marray <- function(EH){
  n.ind <- dim(EH)[1]
  n.occ <- dim(EH)[2]
  m.array <- matrix(data = 0, ncol = n.occ+1, nrow = n.occ)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occ){
    m.array[t,1] <- sum(EH[,t])
  }
  for (i in 1:n.ind){
    pos <- which(EH[i,]!=0) #when captured 
    g <- length(pos) #number of times captured
    #fill in the m-array for each successive release/recapture 
    for (z in 1:(g-1)){ 
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Number of individuals released at t that are never recaptured
  for (t in 1:n.occ){
    m.array[t,n.occ+1] <- m.array[t,1] - sum(m.array[t,2:n.occ])
  }
  out <- m.array[1:(n.occ-1),2:(n.occ+1)]
  return(out)
}

marr <- marray(EH)

#JAGS Model
cat("
    model {
    
    # Define the multinomial likelihood
    for (t in 1:(n.occ-1)){
    marr[t,1:n.occ] ~ dmulti(pr[t, ], R[t])
    }
    # Define the cell probabilities of the m-array
    # Main diagonal - captured at occasion after release
    for (t in 1:(n.occ-1)){
    q[t] <- 1-p[t]                # Probability of non-recapture
    pr[t,t] <- phi[t]*p[t]
    # Above main diagonal - captured again with some internal 0s
    for (j in (t+1):(n.occ-1)){
    pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal - not possible 
    for (j in 1:(t-1)){
    pr[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occ-1)){
    pr[t,n.occ] <- 1-sum(pr[t,1:(n.occ-1)])
    } #t
    
    # Priors and constraints
    for (t in 1:(n.occ-1)){
    phi[t] <- mean.phi
    p[t] <- mean.p
    }
    mean.phi ~ dunif(0, 1)         # Priors for survival
    mean.p ~ dunif(0, 1)           # Priors for recapture
    
    }
    ",file="CJS_mnl.txt")

## JAGS input
data <- list(marr=marr,R=rowSums(marr),n.occ=3)
params <- c("mean.phi","mean.p") 
inits =  function() {list(phi=runif(2),p=runif(2))} 

library(jagsUI)
## RUN
out.jags = jags(data, inits, params, model.file="CJS_mnl.txt",
                n.chains=3, n.iter=20000, n.burnin=5000, n.thin=1)
out.jags$summary

#############################################################################


#############################################################################
