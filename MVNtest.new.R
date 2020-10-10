MVNtest.new <- function (coord, x) {
  # coord includes the range of each coordinate of locations
  # x is n by p, where n is #of locations and p is #of variables

  lon <- coord[,1]
  lat <- coord[,2]
  myPoints <- as.matrix(expand.grid(lon,lat))
  distMat <- sp::spDists(myPoints)
  n <- dim(myPoints)[1]

  nvars <- as.double(ncol(x))
  mean.vector <- t(as.matrix(cov.wt(x)$center))
  covariance.matrix <- cov.wt(x)$cov
  diagcov <- diag(covariance.matrix)
  idiagcov <- 1/sqrt(diagcov)
  v.matrix <- matrix(data = 0, nrow = nvars, ncol = nvars)
  for (i in 1:nvars) v.matrix[i, i] <- idiagcov[i]
  correlation.matrix <- cor(x)
  lambda.matrix <- diag(eigen(correlation.matrix)$values)
  h.matrix <- eigen(correlation.matrix)$vectors
  xhat.matrix <- x - t(matrix(rep(mean.vector, n), nvars))
  xhatprime.matrix <- t(xhat.matrix)
  rprime.matrix <- h.matrix %*% solve(lambda.matrix)^(1/2) %*% 
    t(h.matrix) %*% v.matrix %*% xhatprime.matrix
  rprime <- rprime.matrix
  trprime <- t(rprime)
  
  v1 <- c(numeric(nvars))
  v2 <- c(numeric(nvars))
  v3 <- c(numeric(nvars))
  v4 <- c(numeric(nvars))
  for (j in 1:nvars) v1[j] <- mean(trprime[, j])
  for (j in 1:nvars) v2[j] <- (n^(-1)) * sum((trprime[, j] - v1[j])^2)
  for (j in 1:nvars) v3[j] <- (n^(-1)) * sum((trprime[, j] - v1[j])^3)
  for (j in 1:nvars) v4[j] <- (n^(-1)) * sum((trprime[, j] - v1[j])^4)
  
  F3 <- c(numeric(nvars))
  F4 <- c(numeric(nvars))

 for (k in 1:nvars)  {
  m <- mean(trprime[,k])
  empcov <- (trprime[,k]-rep(m,n))%*%t((trprime[,k]-rep(m,n)))
  cov.df = data.table(cov=c(empcov),dist=round(c(distMat),3))
  cov.hat <- cov.df[,mean(cov),by='dist']
 F3[k] <- sum(abs(cov.hat[,2])^3)
 F4[k] <- sum(abs(cov.hat[,2])^4)
 }
  
  skew <- c(numeric(nvars))
  kurt <- c(numeric(nvars))
  for (j in 1:nvars) skew[j] <- sqrt(n)*v3[j]/sqrt(6*F3[j])
  for (j in 1:nvars) kurt[j] <- sqrt(n)*(v4[j]-3*v2[j])/sqrt(24*F4[j])
  Ep <- t(skew) %*% skew + t(kurt) %*% kurt
  dof <- 2 * nvars
  pvalue <- 1 - pchisq(Ep, dof)

  return(pvalue)
  }