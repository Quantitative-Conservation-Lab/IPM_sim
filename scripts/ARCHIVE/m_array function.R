#m_array function for adults

#from SJC Module_7_CJS_simulate code from SEFS 590 winter 2020

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

marr <- marray(ch)


#from the same file - how it works within the model
#then the model for this is:
## Define the multinomial likelihood
# for (t in 1:(n.occ-1)){
#   marr[t,1:n.occ] ~ dmulti(pr[t, ], R[t])
# }
# # Define the cell probabilities of the m-array
# # Main diagonal - captured at occasion after release
# for (t in 1:(n.occ-1)){
#   q[t] <- 1-p[t]                # Probability of non-recapture
#   pr[t,t] <- phi[t]*p[t]
#   # Above main diagonal - captured again with some internal 0s
#   for (j in (t+1):(n.occ-1)){
#     pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
#   } #j
#   # Below main diagonal - not possible 
#   for (j in 1:(t-1)){
#     pr[t,j] <- 0
#   } #j
# } #t
# # Last column: probability of non-recapture
# for (t in 1:(n.occ-1)){
#   pr[t,n.occ] <- 1-sum(pr[t,1:(n.occ-1)])
# } #t
# 
# # Priors and constraints
# for (t in 1:(n.occ-1)){
#   phi[t] <- mean.phi
#   p[t] <- mean.p
# }
# mean.phi ~ dunif(0, 1)         # Priors for survival
# mean.p ~ dunif(0, 1)           # Priors for recapture# Define the multinomial likelihood
# for (t in 1:(n.occ-1)){
#   marr[t,1:n.occ] ~ dmulti(pr[t, ], R[t])
# }
# # Define the cell probabilities of the m-array
# # Main diagonal - captured at occasion after release
# for (t in 1:(n.occ-1)){
#   q[t] <- 1-p[t]                # Probability of non-recapture
#   pr[t,t] <- phi[t]*p[t]
#   # Above main diagonal - captured again with some internal 0s
#   for (j in (t+1):(n.occ-1)){
#     pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
#   } #j
#   # Below main diagonal - not possible 
#   for (j in 1:(t-1)){
#     pr[t,j] <- 0
#   } #j
# } #t
# # Last column: probability of non-recapture
# for (t in 1:(n.occ-1)){
#   pr[t,n.occ] <- 1-sum(pr[t,1:(n.occ-1)])
# } #t
# 
# # Priors and constraints
# for (t in 1:(n.occ-1)){
#   phi[t] <- mean.phi
#   p[t] <- mean.p
# }
# mean.phi ~ dunif(0, 1)         # Priors for survival
# mean.p ~ dunif(0, 1)           # Priors for recapture


######R[t] is the rowSums(marr) - the number of individuals resighted at time t

