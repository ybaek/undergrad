# Codes to modify MCMC scheme in ngspatial
# For personal use only (minor fixes on copyrighted code methods)
# Incorporates both observation and process errors
# For identifiability issues, observation errors are assumed to be given

# ngspatial:::sparse.sglmm.fit.gaussian modified.
sglmmGaussianFit_fn <- function (Y, X, A, M, sV, beta.start, tol=0.01, minit=10^4, maxit=10^6, hyper=list(sigma.b=10^5)) 
{
  iterations = minit
  n = length(Y)
  Q = t(M) %*% (diag(rowSums(A), n) - A) %*% M
  p = ncol(X)
  ## v is the vector of given observation error variances
  ## sV == solve(diag(v, length(v), length(v)))
  XtX = t(X) %*% sV %*% X
  MtM = t(M) %*% sV %*% M
  beta = matrix(0, iterations + 1, p)
  beta[1, ] = beta.start
  tau.s = numeric(iterations + 1)
  tau.s[1] = 0.1
  tau.h = numeric(iterations + 1)
  tau.h[1] = 1
  q = ncol(M)
  gamma = matrix(0, iterations + 1, q)
  gamma[1, ] = rnorm(q, 0, 1)
  delta = matrix(0, iterations + 1, n)
  delta[1,] = rnorm(n, 0, 1)
  a.s = 0.5
  b.s = 2000
  a.h = hyper$a.h
  b.h = hyper$b.h
  start = 1
  k = 1
  K = diag(1/hyper$sigma.b^2, p)
  pb = pbapply::startpb(0, maxit)
  on.exit(pbapply::closepb(pb))
  cat("\\n")
  flush.console()
  repeat {
    ## Gibbs sampling loop starts
    ## Since process error and observation error are distinctly identified,
    ## delta is needed for gaussian fit method also.
    for (j in (start + 1):(start + iterations)) {
      pbapply::setpb(pb, k)
      k = k + 1
      ## 1. beta: regression coef.
      V = solve(K + XtX)
      mu = V %*% (t(X) %*% sV %*% (Y - M %*% gamma[(j-1), ] - delta[(j-1),]))
      beta[j, ] = t(chol(V)) %*% rnorm(p) + mu
      ## 2. gamma: spatial random shock.
      V = solve(tau.s[j - 1] * Q + MtM)
      mu = V %*% (t(M) %*% sV %*% (Y - X %*% beta[j,] - delta[(j-1),]))
      gamma[j, ] = t(chol(V)) %*% rnorm(q) + mu
      ## 3. delta: process error term.
      V = solve(diag(tau.h[j-1],n,n) + sV)
      mu = V %*% (sV %*% (Y - X%*%beta[j,] - M%*%gamma[j,]))
      delta[j, ] = t(chol(V)) %*% rnorm(n) + mu
      ## 4. tau.h: process error precision.
      b = 0.5 * sum(delta[j, ]^2) + 1/b.h
      tau.h[j] = rgamma(1, a.h + n/2, b)
      ## 5. tau.s: spatial term precision.
      b = 0.5 * t(gamma[j, ]) %*% Q %*% gamma[j, ] + 1/b.s
      tau.s[j] = rgamma(1, a.s + q/2, b)
      if (j == maxit) 
        break
    }
    if (j == maxit) {
      beta = as.matrix(beta[1:maxit, ])
      gamma = as.matrix(gamma[1:maxit, ])
      delta = as.matrix(delta[1:maxit, ])
      tau.s = tau.s[1:maxit]
      tau.h = tau.h[1:maxit]
      break
    }
    done = TRUE
    for (j in 1:p) {
      temp = bm(beta[, j])
      if (temp$se > tol) {
        done = FALSE
        break
      }
    }
    if (done) {
      for (j in 1:q) {
        temp = bm(gamma[, j])
        if (temp$se > tol) {
          done = FALSE
          break
        }
      }
    }
    if (done) {
      for (j in 1:q) {
        temp = bm(delta[, j])
        if (temp$se > tol) {
          done = FALSE
          break
        }
      }
    }
    if (done) {
      temp = bm(tau.s)
      if (temp$se > tol) 
        done = FALSE
    }
    if (done) {
      temp = bm(tau.h)
      if (temp$se > tol) 
        done = FALSE
    }
    if (done) 
      break
    else {
      ## SE tolerance level not reached (not converged)
      ## if j < maxit, runs the loop again
      start = start + iterations
      temp = matrix(0, iterations, p)
      beta = rbind(beta, temp)
      temp = matrix(0, iterations, q)
      gamma = rbind(gamma, temp)
      temp = matrix(0, iterations, n)
      delta = rbind(delta, temp)
      tau.s = c(tau.s, rep(0, iterations))
      tau.h = c(tau.h, rep(0, iterations))
    }
  }
  coefficients = numeric(p)
  beta.mcse = numeric(p)
  names(coefficients) = names(beta.mcse) = colnames(X)
  for (j in 1:p) {
    temp = bm(beta[, j])
    coefficients[j] = temp$est
    beta.mcse[j] = temp$se
  }
  gamma.est = numeric(q)
  gamma.mcse = numeric(q)
  for (j in 1:q) {
    temp = bm(gamma[, j])
    gamma.est[j] = temp$est
    gamma.mcse[j] = temp$se
  }
  delta.est = numeric(n)
  delta.mcse = numeric(n)
  for (j in 1:n) {
    temp = bm(delta[, j])
    delta.est[j] = temp$est
    delta.mcse[j] = temp$se
  }
  temp = bm(tau.s)
  tau.s.est = temp$est
  tau.s.mcse = temp$se
  temp = bm(tau.h)
  tau.h.est = temp$est
  tau.h.mcse = temp$se
  linear.predictors = numeric(n)
  iter = length(tau.s)
  v = numeric(iter)
  for (j in 1:iter) {
    mu = X %*% beta[j, ] + M %*% gamma[j, ]
    linear.predictors = linear.predictors + mu/iter
    v[j] = -2 * sum(dnorm(Y, mu, tau.h[j], log = TRUE))
  }
  D.bar = mean(v)
  pD = D.bar + 2 * sum(dnorm(Y, linear.predictors, 1/sqrt(tau.h.est), log = TRUE))
  dic = D.bar + pD
  fitted.values = linear.predictors + sum(delta)/iter
  residuals = Y - fitted.values
  object = list(y = Y, X = X, M = M, hyper=hyper,
                coefficients = coefficients, fitted.values = fitted.values, 
                linear.predictors = linear.predictors, residuals = residuals, 
                beta.sample = beta, gamma.sample = gamma, delta.sample = delta,
                tau.s.sample = tau.s, tau.h.sample = tau.h, tau.h.mcse = tau.h.mcse,
                beta.mcse = beta.mcse, gamma.mcse = gamma.mcse, delta.mcse = delta.mcse,
                tau.s.mcse = tau.s.mcse, gamma.est = gamma.est, tau.s.est = tau.s.est, 
                tau.h.est = tau.h.est, delta.est = delta.est, iter = iter, dic = dic,  
                D.bar = D.bar, pD = pD)
  class(object) = c("sparse.sglmm")
  object
}