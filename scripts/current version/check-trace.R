# example code for checking traces
for (i in 1:18) {
  plot(aaa[, i], type = "l")
  lines(bbb[, i], col = 2)
  lines(ccc[, i], col = 3)
}
colnames(aaa)
