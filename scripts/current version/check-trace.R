# example code for checking traces
library(here)
aaa <- readRDS(here("data", "out1_1-1.Rdata"))
bbb <- readRDS(here("data", "out1_1-2.Rdata"))
ccc <- readRDS(here("data", "out1_1-3.Rdata"))
for (i in 1:18) {
  plot(aaa[, i], type = "l", col = 1)
  lines(bbb[, i], col = 2)
  lines(ccc[, i], col = 3)
}
colnames(aaa)

aaa <- readRDS(here("data", "out6_1-1.Rdata"))
bbb <- readRDS(here("data", "out6_1-2.Rdata"))
ccc <- readRDS(here("data", "out6_1-3.Rdata"))
for (i in 1:18) {
  plot(aaa[, i], type = "l", col = 1)
  lines(bbb[, i], col = 2)
  lines(ccc[, i], col = 3)
}
colnames(aaa)

# AEB notes
# this is not good
aaa <- readRDS(here("data", "out11_1-1.Rdata"))
bbb <- readRDS(here("data", "out11_1-2.Rdata"))
ccc <- readRDS(here("data", "out11_1-3.Rdata"))
for (i in 1:18) {
  plot(aaa[, i], type = "l", col = 1)
  lines(bbb[, i], col = 2)
  lines(ccc[, i], col = 3)
}
colnames(aaa)
