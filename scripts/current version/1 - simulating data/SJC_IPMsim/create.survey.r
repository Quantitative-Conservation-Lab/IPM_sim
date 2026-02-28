create.survey.bin <- function(Nu, psur, nsurvs, seed = NA){
  
  if (is.na(seed) == TRUE) seed = runif(1, 0, 100)
  set.seed(seed)
  T <- length(Nu)
  SUR <- matrix(NA,nrow=T,ncol=nsurvs)
  for (t in 1:T){
      for(i in 1:nsurvs){
        SUR[t,i] <- rbinom(1,Nu[t], psur[t])
      }
  } # t
  return(SUR)
}


create.survey.norm <- function(Nu, sigma, seed = NA){
  
 # if (is.na(seed) == TRUE) seed = runif(1, 0, 100)
#  set.seed(seed)
  T <- length(Nu)
  SUR <- numeric()
  for (t in 1:T){
    SUR[t] <- rnorm(1, Nu[t], sigma[t])

  } # t
  return(SUR)
}


