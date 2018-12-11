# setwd("Bayesian_class/project/")
# library(nimble)
# library(coda)
# library(basicMCMCplots)

load("data.RData")
y <- dust$log_dust
missing_ind <- which(is.na(y))
N <- dim(dust)[1]
varNames <- c("mu", "phi", "sigma_PN", "sigma_OE")

# Improving MCMC convergence performance
code <- nimbleCode({
  # 1. Top-level priors
  mu ~ dnorm(0, sd=10^6)
  phi ~ dunif(-.9, .9) # plausible for time series AR coef.
  sigma_PN ~ dunif(0, max1)
  sigma_OE ~ dunif(0, max2)
  # 2. Latent process model
  x[1] <- log(104)
  for (t in 2:N) {
    x[t] ~ dnorm(mu + phi*x[t-1], sd=sigma_PN)
  }
  # 3. Observation data
  for (t in 1:N) {
    y[t] ~ dnorm(x[t], sd=sigma_OE)
  }
})
constants <- list(N=N)
data <- list(y=y)
initsFn <- function() {
  list(mu=rnorm(1), phi=rnorm(1,sd=.2), sigma_PN=runif(1,0,0.5), sigma_OE=runif(1,0,0.5),
       x=rep(0,N), max1=10^6, max2=10^6) # max param. will be modified later
}
inits <- initsFn()
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel, control=list(adaptInterval=100), enableWAIC=TRUE)

# Configuring samplers and monitors
conf$removeSamplers(varNames)
conf$addSampler(varNames[1:3], type="AF_slice")
conf$addSampler(varNames[4], type="slice") # independent sampling for obs. error
conf$addMonitors("x", "y") # Predictive nodes
Rmcmc <- buildMCMC(conf)

CompiledList <- compileNimble(list(Rmodel,Rmcmc))
Cmodel <- CompiledList[[1]]
Cmcmc <- CompiledList[[2]]
set.seed(0)
t <- system.time(samples <- runMCMC(Cmcmc, niter=10000, nchains=3, inits=initsFn,
                                    samplesAsCodaMCMC=TRUE, WAIC=TRUE))
t[3] #(~300 secs)
samples$WAIC # -10750.85
save(samples, file="ind_OE.RData")


# Convergence diagnostics
# my_samples <- as.mcmc.list(lapply(samples$samples, function(x) x[,varNames]))
# gelman.diag(my_samples)
# chainsPlot(my_samples, buffer.left=.5, buffer.right=.5, cex=1)
# (ess <- effectiveSize(my_samples))
# min(ess/t[3]) # drastic improvement
# acfplot(my_samples)


