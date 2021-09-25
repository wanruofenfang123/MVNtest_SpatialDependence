## union-intersection test

uit <- function(coord,x, K, alpha) {
  n <- nrow(x)
  p <- ncol(x)
  
  pp <- rep(NA,K)
  for (i in 1:K) {
    theta <- runif(p,0,2*pi)
    a <- cos(theta)
    y <- x %*% a
    pp[i] <- univariate_normality_test_spatial(coord, y)
  }
  
    p_ordered <- sort(pp)
    # R <- ifelse(length(which(p_ordered/c(1:K) < alpha/K))>0,max(which(p_ordered/c(1:K) < alpha/K)),1)
    # T <- p_ordered[R]
    # reject <- NA
    #reject <- ifelse ((any(pp<=alpha)==TRUE), 1, 0)
    reject <- ifelse(any(p_ordered/c(1:K) < alpha/K/sum(1/c(1:K)))==TRUE,1,0)
  
   return(reject)
}