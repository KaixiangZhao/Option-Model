file.list = list.files("reg_data")
dir = paste("./reg_data/", file.list, sep="") 
for (i in 1:length(dir)){
  yx = read.table(dir[i])
  y = yx[, 1]
  x1 = yx[, 2]
  x2 = yx[, 3]
  x3 = yx[, 4]
  mydata = data.frame(y, x1, x2, x3)
  lm.sol = lm(y~x1+x2+x3, data=mydata)
  cof = as.vector(coefficients(lm.sol))
  print(summary(lm.sol)$r.square)
  print(summary(lm.sol)$sigma)
}
