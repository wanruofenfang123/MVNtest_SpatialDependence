### Are You All Normal? It Depends! ###
# R code to perform the simulation studies described in the article
  
## Install packages

library(psych)     # mardia
library(fda)       # fbplot
library(geoR)      # cov.spatial
library(sp)        # spatial distance matrix
library(parallel)  # parallel computing
library(QuantPsyc) # Mardia's test
library(normwhn.test) # Doornik and Hansen (2008) test
library(data.table) 

cl <- makeCluster(25)  # use 25 cores


## Part 1: To reproduce Figure 1 in the article ##
## Simulation studies to investigate the influence of 
## spatial dependence on Mardia's measures of skewness 
## and kurtosis for a bivariate Gaussian random field

# Scenario 1: With exponential covariance and rho=0.5

# pre-determined parameters
h_star <- seq(0.1, 0.9, 0.02)   # effective range
x <- y <- seq(0, 1, length=15)  # 2d regular grid
myPoints <- as.matrix(expand.grid(x,y)) # locations
distMat <- sp::spDists(myPoints)    # distance matrix
d <- dim(myPoints)[1]           # number of locations
sig11 <- 1     # variance of 1st variable
sig22 <- 1     # variance of 2nd variable
rho <- 0.5     # correlation between two variables        
nu11 <- 0.5      # smoothness parameters
nu22 <- 0.5
nu12 <- 0.5
nsimu <- 200     # number of simulations

# Simulate a bivariate random field with pre-determined 
# effective ranges using parallel computing
simu <- function(i) {
  h_star <- seq(0.1, 0.9, 0.02)   # effective range
  x <- y <- seq(0, 1, length=15)  # 2d regular grid
  myPoints <- as.matrix(expand.grid(x,y))
  distMat <- sp::spDists(myPoints)
  d <- dim(myPoints)[1]
  sig11 <- 1
  sig22 <- 1
  nu11 <- 0.5
  nu22 <- 0.5
  nu12 <- 0.5
  rho <- 0.5
  nsimu <- 200
  
  z1 <- z2 <-  matrix(NA, d, length(h_star))
  Sig11 <- Sig22 <- Sig12 <- Sig21 <- array(NA, c(d,d,length(h_star)))
  skew_mat <- c()
  kurtosis_mat <- c()
  
  e <- rnorm(n = 2*d)
  for (k in 1:length(h_star)) {
    fun1 <- function (beta) {h_star[k]/beta * besselK(h_star[k]/beta, 0.5) - 0.05}
    beta11 <- beta22 <- beta12 <- uniroot(fun1, interval = c(0.001, 10))$root
    
    Sig11[,,k] = geoR::cov.spatial(distMat, cov.model = "matern", 
                                   cov.pars = c(sig11, beta11), kappa = nu11)
    Sig22[,,k] = geoR::cov.spatial(distMat, cov.model = "matern", 
                                   cov.pars = c(sig22, beta22), kappa = nu22)
    Sig12[,,k] = geoR::cov.spatial(distMat, cov.model = "matern", 
                                   cov.pars = c(sqrt(sig11)*sqrt(sig22), beta12), kappa = nu12) * rho
    Sig21[,,k] = t(Sig12[,,k])
    Sig <- rbind(cbind(Sig11[,,k],Sig12[,,k]), cbind(Sig21[,,k],Sig22[,,k]))
    
    z <- chol(Sig)%*%e
    z1[,k] <- z[1:d]
    z2[,k] <- z[(d+1):(2*d)]
    
    data <- data.frame(z1=z[1:d], z2=z[(d+1):(2*d)])
    
    mardia <- psych::mardia(data, plot = FALSE)
    skew_mat[k] <- mardia$b1p
    kurtosis_mat[k] <- mardia$b2p
  }
  return(list(skew_mat, kurtosis_mat)) # Mardia's sample skewness/kurtosis
}

data_exp_posrho <- parLapply(cl, 1:nsimu, simu)

# Scenario 2: With exponential covariance and rho=-0.5
# repeat procedures in Scenario 1 by replacing rho=0.5 with rho=-0.5
# Save the results as "data_exp_negrho"

# Scenario 3: With Whittle covariance and rho=0.5
# repeat procedures in Scenario 1 by replacing nu11=0.5, nu22=0.5, nu12=0.5
# with nu11=1, nu22=1, nu12=1
# Save the results as "data_wht_posrho"

# Scenario 4: With exponential covariance and rho=-0.5
# repeat procedures in Scenario 3 by replacing rho=0.5 with rho=-0.5
# Save the results as "data_wht_negrho"

# Extract the results
skew1 <- kurt1 <- array(NA, c(length(h_star), nsimu)) # Scenario 1
skew2 <- kurt2 <- array(NA, c(length(h_star), nsimu)) # Scenario 2
skew3 <- kurt3 <- array(NA, c(length(h_star), nsimu)) # Scenario 3
skew4 <- kurt4 <- array(NA, c(length(h_star), nsimu)) # Scenario 4

for (i in 1:nsimu) {
  skew1[,i] <- data_exp_posrho[[i]][[1]]
  skew2[,i] <- data_exp_negrho[[i]][[1]]
  skew3[,i] <- data_wht_posrho[[i]][[1]]
  skew4[,i] <- data_wht_negrho[[i]][[1]]
  kurt[,i] <- data_exp_posrho[[i]][[2]]
  kurt[,i] <- data_exp_negrho[[i]][[2]]
  kurt[,i] <- data_wht_posrho[[i]][[2]]
  kurt[,i] <- data_wht_negrho[[i]][[2]]
}

skew1_mean <- apply(skew1, 1, mean)
skew2_mean <- apply(skew2, 1, mean)
skew3_mean <- apply(skew3, 1, mean)
skew4_mean <- apply(skew4, 1, mean)
kurt1_mean <- apply(kurt1, 1, mean)
kurt2_mean <- apply(kurt2, 1, mean)
kurt3_mean <- apply(kurt3, 1, mean)
kurt4_mean <- apply(kurt4, 1, mean)

# Plot the skewness/kurtosis as a function of the effective range (Figure 1)
par(mfrow=c(4,2), mar=c(4,3,3,1),mgp=c(2.2,1,0))
fbplot(skew1, xaxt="n", ylim = c(-0.1,1.7), 
       xlab=expression(paste("Effective range "("h*"))), ylab="Mardia's skewness", 
       main=expression(paste("(a) Bivariate skewness (Exponential, ",rho[12]," = ",0.5,")")))
lines((h_star-0.1)/0.02+1, skew1_mean, col="green",lwd=2)
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=0, col = "gray",lty=2)

fbplot(kurt1, xaxt="n", ylim = c(5,13), 
       xlab=expression(paste("Effective range "("h*"))), ylab="Mardia's kurtosis", 
       main=expression(paste("(b) Bivariate kurtosis (Exponential, ",rho[12]," = ",0.5,")")))
lines((h_star-0.1)/0.02+1, kurt1_mean, col="green",lwd=2)
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=8, col = "gray",lty=2)

fbplot(skew2, xaxt="n", ylim = c(-0.1,1.7), 
       xlab=expression(paste("Effective range "("h*"))), ylab="Mardia's skewness", 
       main=expression(paste("(c) Bivariate skewness (Exponential, ",rho[12]," = ",-0.5,")")))
lines((h_star-0.1)/0.02+1, skew2_mean, col="green",lwd=2)
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=0, col = "gray",lty=2)

fbplot(kurt2, xaxt="n", ylim = c(5,13), 
       xlab=expression(paste("Effective range "("h*"))), ylab="Mardia's Kurtosis", 
       main=expression(paste("(d) Bivariate kurtosis (Exponential, ",rho[12]," = ",-0.5,")")))
lines((h_star-0.1)/0.02+1, kurt2_mean, col="green",lwd=2)
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=8, col = "gray",lty=2)

fbplot(skew3, xaxt="n", ylim = c(-0.1,6.5), 
       xlab=expression(paste("Effective range "("h*"))), ylab="Mardia's skewness", 
       main=expression(paste("(e) Bivariate skewness (Whittle, ",rho[12]," = ",0.5,")")))
lines((h_star-0.1)/0.02+1, skew3_mean, col="green",lwd=2)
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=0, col = "gray",lty=2)

fbplot(kurt3, xaxt="n", ylim = c(5,25), xlab=expression(paste("Effective range "("h*"))), ylab="Mardia's Kurtosis", 
       main=expression(paste("(f) Bivariate kurtosis (Whittle, ",rho[12]," = ",0.5,")")))
lines((h_star-0.1)/0.02+1, kurt3_mean, col="green",lwd=2)
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=8, col = "gray",lty=2)
axis(2, at=8)

fbplot(skew4, xaxt="n", ylim = c(-0.1,6.5), xlab=expression(paste("Effective range "("h*"))), ylab="Mardia's skewness", 
       main=expression(paste("(g) Bivariate skewness (Whittle, ",rho[12]," = ",-0.5,")")))
lines((h_star-0.1)/0.02+1, skew4_mean, col="green",lwd=2)
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=0, col = "gray",lty=2)

fbplot(kurt4, xaxt="n", ylim = c(5,25), xlab=expression(paste("Effective range "("h*"))), ylab="Mardia's kurtosis", 
       main=expression(paste("(h) Bivariate kurtosis (Whittle, ",rho[12]," = ",-0.5,")")))
lines((h_star-0.1)/0.02+1, kurt4_mean, col="green",lwd=2)
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=8, col = "gray",lty=2)
axis(2, at=8)


## Part 2: To reproduce Figure 2 in the article ##
## Simulation studies to investigate (1) type I error of 
## the new test and four selected tests for i.i.d. data,
## and (2) empirical power of the new test

# Simulate e_{i,j} (Univariate Gaussian random field on 2d space)

# pre-determined parameters
h_star <- seq(0.1, 0.9, 0.02)   # correlation coefficient
L=15  # size of the spatial domain (LxL locations)
lon <- lat <- seq(0, 1, length=L)  # 2d regular grid
myPoints <- as.matrix(expand.grid(lon,lat)) # locations
distMat <- spDists(myPoints)    # distance matrix
d <- dim(myPoints)[1]           # number of locations
sig <- 1     # variance of the variable
nu <- 0.5    # smoothness parameter
nsimu <- 1000     # number of simulations

# Simulate e_{i,j} with parallel computing
simu <- function(i) {
  h_star <- seq(0.1, 0.9, 0.02)
  L=15
  lon <- lat <- seq(0, 1, length=L)
  myPoints <- as.matrix(expand.grid(lon,lat))
  distMat <- sp::spDists(myPoints)
  d <- dim(myPoints)[1]
  sig <- 1
  nu <- 0.5
  nsimu <- 1000
  
  z <-  matrix(NA, d, length(h_star))
  e <- rnorm(n = d)
  Sig <- array(NA, c(d,d,length(h_star)))
  
  for (k in 1:length(h_star)) {
    fun1 <- function (beta) {h_star[k]/beta * besselK(h_star[k]/beta, nu) - 0.05}
    beta <- uniroot(fun1, interval = c(0.001, 10))$root
    
    Sig[,,k] = geoR::cov.spatial(distMat, cov.model = "matern", 
                                 cov.pars = c(sig, beta), kappa = nu)
    z[,k] <- chol(Sig[,,k]*sig)%*%e
  }
  return(z)
}

e_15x15 <- parLapply(cl, 1:nsimu, simu)

## Simulate a bivaraite spatial moving average process and
## test multivariate normality using the new test and other tests
simu_test <- function (l) {
  L <- 15
  lon <- lat <- seq(0, 1, length=L)
  h_star <- seq(0.1, 0.9, 0.02)
  nsimu <- 1000
  pvalue <- matrix(NA,length(h_star),5) 
  
  for (k in 1:length(h_star)) {
    e <- matrix(e_15x15[[l]][,k],L,L)
    e0 <- matrix(NA, L+2, L+2)
    for (i in 1:(L+2)) {
      for (j in 1:(L+2)) {
        if (i>=2 & i<=(L+1) & j>=2 & j<=(L+1)) e0[i,j]=e[i-1,j-1]
        else e0[i,j]=0
      }
    }   
    theta1 <- 1/5
    theta2 <- -1/5
    X1 <- X2 <- matrix(NA, L+2, L+2)
    for (i in 2:(L+1)) {
      for (j in 2:(L+1)) {
        X1[i,j] <-  theta1*(e0[i-1,j]+e0[i+1,j]+e0[i,j-1]+e0[i,j+1])+e0[i,j]
        X2[i,j] <-  theta2*(e0[i-1,j]+e0[i+1,j]+e0[i,j-1]+e0[i,j+1])+e0[i,j]
      }
    }
    X1 <- X1[-c(1,L+2),-c(1,L+2)]
    X2 <- X2[-c(1,L+2),-c(1,L+2)]
    
    # simulated data from bivariate spatial moving average model 
    data <- cbind(c(X1),c(X2))
    
    source("MVNtest.new.R")
    source("DHtest.R")
    # new test for MVN, accounting for spatial dependence
    pvalue[k,1] <- MVNtest.new(cbind(lon,lat),data)
    # variance not adapted, spatial dependence ignored
    pvalue[k,2] <- DHtest(cbind(lon,lat),data)
    # Doornik and Hansen (2008) test
    pvalue[k,3] <- normwhn.test::normality.test1(data)
    # Mardia's skewness test
    pvalue[k,4] <- QuantPsyc::mult.norm(data,chicrit=.001)$mult.test[1,3] 
    # Mardia's kurtosis test
    pvalue[k,5] <- QuantPsyc::mult.norm(data,chicrit=.001)$mult.test[2,3]
  }
  return(pvalue)
}

nsimu <- 1000
pvalue_5tests_15x15 <- parLapply(cl, 1:nsimu, simu_test)

## type I error
p <- array(NA,c(length(h_star),nsimu,5))
for (i in 1:nsimu) {
  p[,i,] <- pvalue_5tests_15x15[[i]]
}
error_15x15 <- matrix(NA,length(h_star),5)
for (k in 1:5) {
  for (j in 1:length(h_star)) {
    error_15x15[j,k] <- sum(p[j,,k]<0.05)/nsimu
  }
}


## Transform the spatial moving average process data to non-Gaussian data
## and test multivariate normality using the new test 
simu_test_nGRF <- function (l) {
  L <- 15
  lon <- lat <- seq(0, 1, length=L)
  h_star <- seq(0.1, 0.9, 0.02)
  nsimu <- 1000
  pvalue <- rep(NA,length(h_star))
  g <- c(0.5,0)
  h <- c(0.5,1)
  for (k in 1:length(h_star)) {
    e <- matrix(e_15x15[[l]][,k],L,L)
    e0 <- matrix(NA, L+2, L+2)
    for (i in 1:(L+2)) {
      for (j in 1:(L+2)) {
        if (i>=2 & i<=(L+1) & j>=2 & j<=(L+1)) e0[i,j]=e[i-1,j-1]
        else e0[i,j]=0
      }
    }   
    theta1 <- 1/5
    theta2 <- -1/5
    X1 <- X2 <- matrix(NA, L+2, L+2)
    for (i in 2:(L+1)) {
      for (j in 2:(L+1)) {
        X1[i,j] <-  theta1*(e0[i-1,j]+e0[i+1,j]+e0[i,j-1]+e0[i,j+1])+e0[i,j]
        X2[i,j] <-  theta2*(e0[i-1,j]+e0[i+1,j]+e0[i,j-1]+e0[i,j+1])+e0[i,j]
      }
    }
    X1 <- X1[-c(1,L+2),-c(1,L+2)]
    X2 <- X2[-c(1,L+2),-c(1,L+2)]
    
    # inverse sinh-arcsinh transformation
    X1_tr <- sinh(1/h[1]*(asinh(X1)+g[1]))
    X2_tr <- sinh(1/h[2]*(asinh(X2)+g[2]))
    
    data <- cbind(c(X1_tr),c(X2_tr))
    
    source("MVNtest.new.R")
    # new test for MVN, accounting for spatial dependence
    pvalue <- MVNtest.new(cbind(lon,lat),data)
  }
  return(pvalue)
}

nsimu <- 1000
pvalue_newtest_nGRF_15x15 <- parLapply(cl, 1:nsimu, simu_test_nGRF)

# Empirical power of the new test
p <- matrix(NA,length(h_star),nsimu)
for (i in 1:nsimu) {
  p[,i] <- pvalue_newtest_nGRF_15x15[[i]]
}
power_15x15 <- rep(NA,length(h_star))
  for (j in 1:length(h_star)) {
    power_15x15[j] <- sum(p[j,]<0.05)/nsimu
  }

# Repeat the above procedures for L=30, and save type I errors as
# "error_30x30", and empirical power as "power_30x30"

# Type I error and empirical power (Figure 2)
par(mfrow=c(2,2),mar=c(3,3,2,1),mgp=c(2,0.5,0))
plot(error_15x15[,1],ylim=c(min(error),max(error)),xaxt="n",type="l",
     xlab=expression(paste("Effective range "("h*"))),ylab="Type I error",
     main="Type I error of JB* and other MVN tests (15x15)",col=5)
legend("topleft",c("JB*",expression(paste("JB"["BS"])),
                   expression(paste("JB"["DH"])),"MS","MK"), 
       col=1:5, lty=1, bty="n",cex = 0.65)
lines(error_15x15[,2],col=2)
lines(error_15x15[,3],col=3)
lines(error_15x15[,4],col=4)
lines(error_15x15[,5],col=1)
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=0.05,col="red",lty=2)
mtext(2,at=0.05,text = "0.05",col="red",line = 0.2,cex = 0.75)

plot(error_30x30[,1],ylim=c(min(error),max(error)),xaxt="n",type="l",
     xlab=expression(paste("Effective range "("h*"))),ylab="Type I error",
     main="Type I error of JB* and other MVN tests (30x30)",col=5)
legend("topleft",c("JB*",expression(paste("JB"["BS"])),
                   expression(paste("JB"["DH"])),"MS","MK"), 
       col=1:5, lty=1, bty="n",cex = 0.65)
lines(error_30x30[,2],col=2)
lines(error_30x30[,3],col=3)
lines(error_30x30[,4],col=4)
lines(error_30x30[,5],col=1)
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=0.05,col="red",lty=2)
mtext(2,at=0.05,text = "0.05",col="red",line = 0.2,cex = 0.75)

plot(power_15x15,xaxt="n",type="l",
     xlab=expression(paste("Effective range "("h*"))),ylab="Empirical power",
     main="Empirical power of JB* (15x15)")
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=0.05,col="red",lty=2)
mtext(2,at=0.05,text = "0.05",col="red",line = 0.2,cex = 0.75)

plot(power_30x30,xaxt="n",type="l",
     xlab=expression(paste("Effective range "("h*"))),ylab="Empirical power",
     main="Empirical power of JB* (30x30)")
axis(1,at=1+5*c(0:8), cex.axis=.9, labels=seq(0.1,0.9,length=9), tck=-.02)
abline(h=0.05,col="red",lty=2)
mtext(2,at=0.05,text = "0.05",col="red",line = 0.2,cex = 0.75)
