library(here)
library(popbio)
library(DescTools)

# NOTE script takes ~30 mins to run

# POPULATION SCENARIOS ----

eval <- function(SJ.lo,SJ.hi,SJ.seq,SA.lo,SA.hi,SA.seq,f.lo,f.hi,f.seq){
  grid <- expand.grid(seq(SJ.lo,SJ.hi,SJ.seq),seq(SA.lo,SA.hi,SA.seq),seq(f.lo,f.hi,f.seq))
  store <- rep(NA,nrow(grid))
  for(i in 1:nrow(grid)){ 
    pop.mat.rd <- matrix(NA,nrow=2,ncol=2)
    pop.mat.rd[1,] <- c(grid[i,1]*grid[i,3], grid[i,1]*grid[i,3])
    pop.mat.rd[2,] <- c(grid[i,2], grid[i,2])
  
    store[i] <- eigen.analysis(pop.mat.rd)$lambda1
  }
    out <- cbind(grid,store)
    return(out)
}

keep <- data.frame("scenario"=c("fast,stable","fast,decline","fast,increase","slow,stable","slow,decline","slow,increase","mod,stable","mod,decline","mod,increase"),"S.J"=rep(NA,9),"S.A"=rep(NA,9),"f"=rep(NA,9),"lambda"=rep(NA,9))

##FAST life history 
out.fast <- eval(SJ.lo = 0.28,SJ.hi = 0.32,SJ.seq = 0.001,SA.lo = 0.38,SA.hi = 0.42,SA.seq = 0.001,f.lo = 1.8,f.hi = 2.2,f.seq = 0.001)

#eyeballing results to choose best
#rule of thumb is that we prefer changes in both survival rates before fecundity 
#then in juvenile survival and fecundity 
#this is just for consistency 

#find best fast, stable 
out.fast[Closest(x=out.fast[,4],a=1.0,which = TRUE,na.rm=FALSE),]
keep[1,2:5] <- out.fast[337041,1:4]

#find best fast, declining 
out.fast[Closest(x=out.fast[,4],a=0.96,which=TRUE,na.rm=FALSE),]
keep[2,2:5] <- out.fast[336211,1:4]

#find best fast, increasing 
out.fast[Closest(x=out.fast[,4],a=1.04,which=TRUE,na.rm=FALSE),]
keep[3,2:5] <- out.fast[337466,1:4]

##SLOW life history 
out.slow <- eval(SJ.lo = 0.48,SJ.hi = 0.52,SJ.seq = 0.001,SA.lo = 0.58,SA.hi = 0.62,SA.seq = 0.001,f.lo = 0.6,f.hi = 1.0,f.seq = 0.001)

#find best slow, stable 
out.slow[Closest(x=out.slow[,4],a=1.0,which = TRUE,na.rm=FALSE),]
keep[4,2:5] <- out.slow[337041,1:4]

#find best slow, declining 
out.slow[Closest(x=out.slow[,4],a=0.96,which=TRUE,na.rm=FALSE),]
keep[5,2:5] <- out.slow[252971,1:4]

#find best slow, increasing 
out.slow[Closest(x=out.slow[,4],a=1.04,which=TRUE,na.rm=FALSE),]
keep[6,2:5] <- out.slow[379537,1:4]

## MODERATE life history 
out.mod <- eval(SJ.lo = 0.38,SJ.hi = 0.42,SJ.seq = 0.001,SA.lo = 0.48,SA.hi = 0.52,SA.seq = 0.001,f.lo = 1.1,f.hi = 2.5,f.seq = 0.001)

#find best slow, stable 
out.mod[Closest(x=out.mod[,4],a=1.0,which = TRUE,na.rm=FALSE),]
keep[7,2:5] <- out.mod[252991,1:4]

#find best slow, declining 
out.mod[Closest(x=out.mod[,4],a=0.96,which=TRUE,na.rm=FALSE),]
keep[8,2:5] <- out.mod[252356,1:4]

#find best slow, increasing 
out.mod[Closest(x=out.mod[,4],a=1.04,which=TRUE,na.rm=FALSE),]
keep[9,2:5] <- out.mod[253626,1:4]

write.csv(keep,
          here("data","demographic_scenarios.csv"),
          row.names=FALSE)
saveRDS(keep, here("data", "demographic_scenarios.RDS"))

# DATA SCENARIOS ----

# construct matrix of scenarios #####
det.abund <- factor(x = c("L", "M", "H"))
det.MR <- factor(x = c("L", "M", "H", "NA"))
det.prod <- factor(x = c("L", "M", "H", "NA"))
data_scenarios <- expand.grid(
  det.abund = det.abund, 
  det.MR = det.MR, 
  det.prod = det.prod
  )

write.csv(data_scenarios,
          here("data","data_scenarios.csv"),
          row.names=FALSE)
saveRDS(data_scenarios, here("data", "data_scenarios.RDS"))
