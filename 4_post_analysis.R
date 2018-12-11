
# Prior sensitivy analysis
# Posterior predictive check
# Final inference (itertions so that we obtain larger ESS)

#setwd("Bayesian_class/project/")
library(nimble)
library(coda)
library(basicMCMCplots)

# global model objects
load("data.RData")
y <- dust$log_dust
missing_ind <- which(is.na(y))
N <- dim(dust)[1]
varNames <- c("mu", "phi", "sigma_PN", "sigma_OE")
constants <- list(N=N)
data <- list(y=y)

# 2 models, distinct priors
# code_orig <- nimbleCode({
#   # 1. Top-level priors
#   mu ~ dnorm(0, sd=10^6)
#   phi ~ dunif(-.9, .9) # plausible for time series AR coef.
#   sigma_PN ~ dunif(0, max1)
#   sigma_OE ~ dunif(0, max2)
#   # 2. Latent process model
#   x[1] <- log(104)
#   for (t in 2:N) {
#     x[t] ~ dnorm(mu + phi*x[t-1], sd=sigma_PN)
#   }
#   # 3. Observation data
#   for (t in 1:N) {
#     y[t] ~ dnorm(x[t], sd=sigma_OE)
#   }
# })
# inits1 <- list(mu=0, phi=.1, sigma_PN=.5, sigma_OE=.5, x=rep(0,N), max1=1, max2=1)

code_inv <- nimbleCode({
  # 1. Top-level priors
  mu ~ dnorm(0, sd=10^6)
  phi ~ dunif(-.9, .9) 
  sigma_PN ~ dinvgamma(shape=2, scale=1) # different prior: inverse gamma (vague)
  sigma_OE ~ dinvgamma(shape=2, scale=1)
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
inits2 <- list(mu=0, phi=.1, sigma_PN=.5, sigma_OE=.5, x=rep(0,N))

# Rmodel1 <- nimbleModel(code_orig, constants, data, inits1)
Rmodel2 <- nimbleModel(code_inv, constants, data, inits2)

# conf1 <- configureMCMC(Rmodel1, control=list(adaptInterval=100), enableWAIC=TRUE)
# conf1$removeSamplers(varNames)
# conf1$addSampler(varNames[1:3], type="AF_slice")
# conf1$addSampler(varNames[4], type="slice") # independent sampling for obs. error
# conf1$addMonitors("x", "y") # Predictive nodes
# Rmcmc1 <- buildMCMC(conf1)

conf2 <- configureMCMC(Rmodel2, control=list(adaptInterval=100), enableWAIC=TRUE)
conf2$removeSamplers(varNames)
conf2$addSampler(varNames[1:3], type="AF_slice")
conf2$addSampler(varNames[4], type="slice") # independent sampling for obs. error
conf2$addMonitors("x", "y") # Predictive nodes
Rmcmc2 <- buildMCMC(conf2)

# CL1 <- compileNimble(list(Rmodel1,Rmcmc1))
# set.seed(0)
# samples1 <- runMCMC(CL1[[2]], samplesAsCodaMCMC=TRUE, WAIC=TRUE)

CL2 <- compileNimble(list(Rmodel2,Rmcmc2))
Cmodel2 <- CL2[[1]]
Cmcmc2 <- CL2[[2]]
initsFn2 <- function() {
  list(mu=rnorm(1), phi=rnorm(1,sd=.2), sigma_PN=runif(1,0,0.5), sigma_OE=runif(1,0,0.5), x=rep(0,N)) # max param. will be modified later
}
set.seed(0)
samples2 <- runMCMC(Cmcmc2, nchains=3, inits=initsFn2, samplesAsCodaMCMC=TRUE, WAIC=TRUE)
save(samples2, file="invgamma.RData")

# Take a look...
# samples1$WAIC # -11612
samples2$WAIC # 108.7059
# info <- samples1$samples[,varNames]
invg <- as.mcmc.list(lapply(samples2$samples, function(x) x[,varNames]))
acfplot(invg, aspect="fill")

samplesPlot(invg$chain3, burnin=2000, scale=T) # better mixing for sigma_OE
effectiveSize(invg) # very large ESS for sigma_OE
gelman.diag(invg) # Perfect convergence results

# Plotting
png("priors.png")
par(mar=c(3,3,2,2)+.1, mgp=c(2,1,0))
plot(density(invg$chain3[2001:10000,"sigma_OE"]), 
     xlim=c(0, 0.065), ylim=c(0, 150), col="blue", lwd=2,
     main="Posterior density curves") # slightly larger ~.04
lines(density(unif$chain3[2001:10000,"sigma_OE"]), col="red", lwd=2)
legend("topright", legend=c("Inv. gamma prior", "Unif. prior"),
       col=c("blue", "red"), lty=1, lwd=2, cex=1.2)
dev.off()

# comparison
load("ind_OE.RData")
unif <- as.mcmc.list(lapply(samples$samples, function(x) x[,varNames]))


# What do I use???
# Posterior predictive check: simulation-based
unif_params <- unif$chain3[2001:10000,]
invg_params <- invg$chain3[2001:10000,]

# For loop
simData <- base::array(NA, c(N,50,2)) # 50 simulations
Cmodel$resetData()
Cmodel2$resetData()
modelList <- list(Cmodel, Cmodel2)
paramsList <- list(unif_params, invg_params)
set.seed(0)
for (i in 1:2) { 
  model <- modelList[[i]]
  params <- paramsList[[i]]
  for (j in 1:50) {
    ind <- base::sample(1:8000,1) # Sampling from the posterior
    for (v in varNames) model[[v]] <- params[ind,v] # Setting the parameters
    model$calculate() # propagating throughout the dependency chain
    model$simulate("x[2:2184]")
    model$calculate()
    model$simulate("y")
    simData[,j,i] <- model$y
  } 
}

save(simData, file="simulated_array.RData")

sim_means <- apply(simData, c(2,3), mean)
sim_sds <- apply(simData, c(2,3), sd)

# Plotting
png("check.png")
par(mar=c(2,3,1,1)+.1, mgp=c(2,1,0), mfrow=c(2,1))
hist(sim_means[,1],freq=F, ylim=c(0, 35), col="lightblue", 
     main="Simulated means", xlab="")
hist(sim_means[,2],freq=F, add=T, col="orange")
abline(v=mean(y, na.rm=T),lwd=2,col="red")
legend("topright", legend=c("Inv.gamma", "Unif."), col=c("orange", "lightblue"), lwd=10,
       cex=1.1)
hist(sim_sds[,2],freq=F, ylim=c(0, 60), breaks=30, col="orange",
     main="Simulated standard deviations", xlab="")
hist(sim_sds[,1],freq=F, breaks=20, add=T, col="lightblue")
legend("topright", legend=c("Inv.gamma", "Unif."), col=c("orange", "lightblue"), lwd=10,
       cex=1.3)
dev.off()

# Seems like I will go with inverse gamma 
# Standard deviation--> tells me it is still very underdispersed
# Will need a more complicated model (i.e., regression variables)

# Final inference + interpolation
# For better inference, I'll have to run more iterations...

#...
# samplesSummary(samples)

missing_ind <- which(is.na(y))
missing_data <- paste0("y[", missing_ind, "]")

interp_mean <- samplesSummary(samples2$samples$chain3[,missing_data])[,"Mean"]

# samples[,missing_data]


