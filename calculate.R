file.list = list.files("reg_data")
dir = paste("./reg_data/", file.list, sep="") 
mat = matrix(nrow=129, ncol=9)
for (i in 1:length(dir)){
  yx = read.table(dir[i])
  y = yx[, 1]
  x1 = yx[, 2]
  x2 = yx[, 3]
  x3 = yx[, 4]
  mydata = data.frame(y, x1, x2, x3)
  lm.sol = lm(y~x1+x2+x3, data=mydata)
  cof = as.vector(coefficients(lm.sol))
  r.square = summary(lm.sol)$r.square
  p = as.vector(summary(lm.sol)$coefficients[,4])
  result = c(cof[1], p[1], cof[2], p[2], cof[3], p[3], cof[4], p[4], r.square)
  mat[i,] = result
}
file.list = c(file.list, "avg")
mat[129,] = c(mean(mat[-c(29:32, 61:64, 93:96, 125:128, 129),1]), 0,
              mean(mat[-c(29:32, 61:64, 93:96, 125:128, 129),3]), 0,
              mean(mat[-c(29:32, 61:64, 93:96, 125:128, 129),5]), 0,
              mean(mat[-c(29:32, 61:64, 93:96, 125:128, 129),7]), 0,
              mean(mat[-c(29:32, 61:64, 93:96, 125:128, 129),9]))
df = data.frame(file.list, mat)
write.csv(df, file="reg_result.csv")
