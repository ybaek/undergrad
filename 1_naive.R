setwd("Bayesian_class/project/")
dust <- read.csv("PRSA_data_2010.1.1-2014.12.31.csv", header=T)
dust <- subset(dust, (year==2014 & month >= 9 & month <= 11))
dim(dust) # 2184 observations in total
head(dust)
dust$log_dust <- log(dust$pm2.5)
# 20 missing values -> we want to interpolate (predict)
dust$precp <- dust$Is + dust$Ir
save(dust, file="data.RData")

# Global constants and data structures
y <- dust$log_dust
missing_ind <- which(is.na(y))
N <- dim(dust)[1]
varNames <- c("mu", "phi", "sigma_PN", "sigma_OE")

# Model formulation: state-space representation
library(nimble)
code <- nimbleCode({
  # 1. Top-level priors
  mu ~ dnorm(0, sd=10^6)
  phi ~ dunif(-.9, .9) # plausible for time series AR coef.
  sigma_PN ~ dunif(0, 10^6)
  sigma_OE ~ dunif(0, 10^6)
  # 2. Latent process model
  x[1] <- log(104) # The observed data from Aug. 31st
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
  list(mu=rnorm(1), phi=rnorm(1,sd=.1), sigma_PN=runif(1), sigma_OE=runif(1),
       x=rep(0,N))
}
inits <- initsFn()
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel, enableWAIC=TRUE)

# Samplers are all RW's (except mu)
conf$addMonitors("x", "y") # Predictive nodes
Rmcmc <- buildMCMC(conf)

CompiledList <- compileNimble(list(Rmodel,Rmcmc))
Cmodel <- CompiledList[[1]]
Cmcmc <- CompiledList[[2]]
set.seed(0)
t <- system.time(samples <- runMCMC(Cmcmc, niter=10000, nchains=3, inits=initsFn,
                                    samplesAsCodaMCMC=TRUE, WAIC=TRUE))
t[3] # (~220 secs)
save(samples, file="naive_samples.RData")

load("naive_samples.RData")
varNames <- c("mu", "phi", "sigma_PN", "sigma_OE")
library(coda)
library(basicMCMCplots)
samples$WAIC #-5805.074
chainsPlot(samples$samples, varNames, buffer.left=.5, buffer.right=.5, cex=1)
samplesPlot(samples$samples$chain1, varNames, scale=T)


my_samples <- as.mcmc.list(lapply(samples$samples, function(x) x[,varNames]))
gelman.diag(my_samples) # comes close to convergence
(ess <- effectiveSize(my_samples)) # extremely low
min(ess/t[3])
cor(my_samples$chain1)

