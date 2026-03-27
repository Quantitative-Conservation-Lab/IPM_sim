create.survey.bin <- function(Nu, psur, sec, seed = NA){
  
T <- length(Nu)
SUR <- matrix(NA,nrow=T,ncol=sec)
for (t in 1:T){
  for(i in 1:sec){
    SUR[t,i] <- rbinom(1,Nu[t], psur[t])
  }
} # t
return(SUR)
}


create.survey.norm <- function(Nu, sigma, seed = NA){
  
T <- length(Nu)
SUR <- numeric()
for (t in 1:T){
  SUR[t] <- rnorm(1, Nu[t], sigma[t])
} # t
return(SUR)
}


