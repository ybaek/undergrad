# MSTM implementation codes for gaussian process
# Focusing on replicating results of Bradley et al. (2015) for now
# For personal use only

library(FKF)
library(mvtnorm)
library(batchmeans)

mstm_fn <- function(Y, X, S, M, Kst, Wst,
                     beta.start, hyper, burnin=10^3, iter=10^4) {
  # PRE: each argument is matrix/array of conformable dimensions
  # hyper is a list containing sigma.b, a, b
  t <- ncol(Y)
  n <- nrow(Y)
  p <- dim(X)[2]
  r <- dim(S)[2]
    
  #1. Initializing process variables/parameters
  beta = array(0, c(p,t,iter+1))
  beta[,,1] = beta.start
  eta = array(0, c(r,t,iter+1))
  eta[,,1] = matrix(rnorm(r*t),r,t)
  u = array(0, c(r,t,iter+1))
  u[,,1] = matrix(rnorm(r*t),r,t)
  tau.xi = matrix(0,t,iter+1)
  tau.xi[,1] = rnorm(t)
  tau.k = numeric(iter+1)
  tau.k[1] = rnorm(1)
  betaCov = diag(1/sqrt(hyper$sigma.b), p, p)
  a.xi = a.k = hyper$a
  b.xi = b.k = hyper$a
  ## process bar for MCMC
  start = 1
  k = 1
  pb = pbapply::startpb(0, iter)
  on.exit(pbapply::closepb(pb))
  cat("\\n")
  flush.console()
  
  #2. Gibbs sampling starts
  for (j in (start+1):(start+iter)) {
    pbapply::setpb(pb, k)
    k = k+1
    kb = numeric(t-2) ## this is needed for sampling tau.k
    
    ## 1. eta
    ## FFBS scheme needed
    ## First: Kalman smoothing
    ## (...) <- use fkf package -> result object kalman
    
    
    ## Second: Backward sampling of eta's
    eta[,t,j] = t(chol(kalman$Ptt[,,t])) %*% rnorm(r) + kalman$att[,t]
    for (i in 1:(t-1)) {
      J = kalman$Ptt[,,i] %*% t(M[,,i]) %*% solve(kalman$Pt[,,i]) 
      mu = kalman$at[,i] + J %*% (eta[,(t-i+1),j] - kalman$at[,i])
      sigma = kalman$Ptt[,,i] - J %*% kalman$Pt[,,i] %*% t(J)
      eta[,(t-i),j] = t(chol(sigma)) %*% rnorm(r) + mu
    }
    ## iteration done for eta.

    for (i in 1:t) {
      ## inner loop: block sampling done for all time periods
      ## 2. beta
      Cov = solve(betaCov + tau.xi[j-1]*t(X[,,i])%*%X[,,i])
      Mu = tau.xi[j-1]*Cov %*% t(X[,,i]) %*% (Y[,i] - S[,,i]%*%eta[,i,j])
      beta[,,j] = t(chol(Cov)) %*% rnorm(p) + Mu
        
      ## 3. tau.xi
      b = 1/2 * sum((y[,i] - X[,,i]%*%beta[,i,j]-S[,,i]%*%eta[,i,j])^2)
      tau.xi[i,j] = rgamma(1, a.xi + n/2, b.xi + b)  
        
      ## 4. tau.k: not dependent on t
      ## (...)
      if (i < t-2) {
        xi = eta[,i+2,j]-M[,,i+2]%*%eta[,i+1,j]-N[,,i+2]%*%eta[,i,j]
        kb[i+2] = 1/2 * t(xi) %*% solve(Wst[,,i+2]) %*% xi
      }
    }
    b = 1/2 * t(eta[,1,j]) %*% solve(Kst[,,1]) %*% eta[,1,j] + sum(kb)
    tau.k[j] = rgamma(1, a.xi + t*r/2, b.xi + b)
  }
  fitted = matrix(0, n, t)
  for (i in 1:t) {
    for (j in (burnin+1):iter) {
      fitted[,i] = fitted[,i] + (X[,,i] %*% beta[,i,j] + S[,,i] %*% eta[,i,j])/(iter-burnin)
    }
  }
  residuals = Y - fitted
  result = list(fitted = fitted, residuals = residuals)
  result
}


