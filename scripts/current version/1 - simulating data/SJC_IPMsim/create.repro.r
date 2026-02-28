create.reproduction <- function(ind, prep, seed = NA){
  
  if (is.na(seed) == TRUE) seed = runif(1, 0, 100)
  set.seed(seed)
  r <- dim(ind)[1] - 1
  maxAge <- dim(ind)[1] - 4
  T <- dim(ind)[2]
  rep <- year <- age <- numeric()
  le <- 0
  for (t in 1:T){
if(length(which(!is.na(ind[r,t,])==TRUE))==0){   ###################ADDED
  print("extinct")                               ###################ADDED
}else{                                           ###################ADDED
    z <- which(!is.na(ind[r,t,]))
    for (i in 1:length(z)){
      j <- le + i
      h <- rbinom(1, 1, prep[t])
      if (h==1){
        rep[j] <- ind[r,t,z[i]]
        year[j] <- t
        age[j] <- which(!is.na(ind[2:(maxAge+2),t,z[i]]))
      } # if
      else {
        rep[j] <- NA
        year[j] <- NA
        age[j] <- NA
      } # else  
    } # i
    le <- length(rep)
  } # t
}                                              ###################ADDED 
  age[age==(maxAge+1)] <- maxAge   # re-adjust age of immigrants to maximal age
  k <- which(!is.na(rep))
  rep.ind <- cbind(rep[k], year[k], age[k])
  colnames(rep.ind) <- c("Reproduction", "Year", "Age of mother")
  rep.agg <- matrix(NA, nrow = T, ncol = 2)
  for (t in 1:T){
    rep.agg[t,1] <- sum(rep.ind[rep.ind[,2]==t,1])
    rep.agg[t,2] <- length(rep.ind[rep.ind[,2]==t,1])
  }
  colnames(rep.agg) <- c("Juveniles", "Surveyed broods")
  
  return(list(rep.ind = rep.ind, rep.agg = rep.agg))
}




