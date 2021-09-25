

univariate_normality_test_spatial <- function(coord, x) {
  library(data.table)
  # coord includes the coordinates of locations
  # x is n by 1, where n is #of locations
  lon <- coord[,1] 
  lat <- coord[,2] 
  myPoints <- as.matrix(expand.grid(lon,lat))
#  distMat <- sp::spDists(myPoints) 
  n <- dim(myPoints)[1]   # number of locations
  
  mean <- mean(c(x))
  var <- var(c(x))
  y <- (x-mean)/sqrt(var)  # studentized data
  
  SS <- sum(c(y)^3-3*c(y))/sqrt(n)  # modified skewness
  KK <- sum(c(y)^4-6*c(y)^2+3)/sqrt(n) # modified kurtosis
  
  ys <- y^3-3*y
  yk <- y^4-6*y^2+3
  meanys <- mean(ys)
  meanyk <- mean(yk)
  
  N <- length(lon)
  M <- length(lat)
  ys2 <- matrix(ys,N,M)
  yk2 <- matrix(yk,N,M)
  hN <- floor(4*(N/100)^(2/9))
  hM <- floor(4*(M/100)^(2/9))
 # kk <- 
    gammaS <- gammaK <- matrix(NA,2*hN+1,2*hM+1)
  for (u in c(-hN:hN)) {
    for (v in c(-hM:hM)) {
      uu <- u + hN + 1
      vv <- v + hM + 1
    #  kk[uu,vv] <- ifelse(abs(u)<=abs(hN),1,0)*ifelse(abs(v)<=abs(hM),1,0)
      gammaS[uu,vv] <- gammaK[uu,vv] <- 0
      card <- 0
      for (i in 1:N) {
        for (j in 1:M) {
          if (i+u>=1 && i+u<=N && j+v>=1 && j+v<=M) {
            gammaS[uu,vv] <- gammaS[uu,vv]+(ys2[i,j]-meanys)*(ys2[i+u,j+v]-meanys)
            gammaK[uu,vv] <- gammaK[uu,vv]+(yk2[i,j]-meanyk)*(yk2[i+u,j+v]-meanyk)
            card <- card + 1
          }
        }
      }
      if (card>0) {
        gammaS[uu,vv] <- gammaS[uu,vv]/card*(1-abs(u)/hN)*(1-abs(v)/hM)
        gammaK[uu,vv] <- gammaK[uu,vv]/card*(1-abs(u)/hN)*(1-abs(v)/hM)
      }
    }
  }
  phiS2 <- sum(gammaS)
  phiK2 <- sum(gammaK)
  
  
  ## compute the sample covariances
  #empcov <- (ys-rep(meanys,n))%*%t((ys-rep(meanys,n)))
  #cov.df = data.table(cov=c(empcov),dist=round(c(distMat),3))
  #cov.hat <- cov.df[,mean(cov),by='dist']
  #gammaS <- cov.hat[,2]
  #phiS <- 2*sum(gammaS^3)
  
  #empcov <- (yk-rep(meanyk,n))%*%t((yk-rep(meanyk,n)))
  #cov.df = data.table(cov=c(empcov),dist=round(c(distMat),3))
  #cov.hat <- cov.df[,mean(cov),by='dist']
  #gammaK <- cov.hat[,2]
  #phiK <- 2*sum(gammaS^4)
  
  JB <- SS^2/phiS2 + KK^2/phiK2
  pvalue <- 1 - pchisq(JB, 2)
  
  return(pvalue)
}