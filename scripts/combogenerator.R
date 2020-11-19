abundinc <- c("Y", "N")
MRinc <- c("Y", "N")
nestinc <- c("Y", "N")

param1 <- c("L", "H")
param2 <- c("L", "H")
param3 <- c("L", "H")
param4 <- c("L", "H")
param5 <- c("L", "H")

paramcombs <- expand.grid(param1, param2, param3, param4, param5)
paramcombs <- cbind(paramcombs, paramcombs[, 5])

combs <- expand.grid(abundinc, MRinc, nestinc, param1, param2, param3, param4, param5)

combs <- cbind(combs, combs[, 8])



write.csv(combs, )