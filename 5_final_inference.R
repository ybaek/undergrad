#
# Final model, analysis summary, and inference
#

library(nimble)
library(basicMCMCplots)
library(coda)

# global model objects
load("data.RData")
y <- dust$log_dust
missing_ind <- which(is.na(y))
N <- dim(dust)[1]
varNames <- c("mu", "phi", "sigma_PN", "sigma_OE")
constants <- list(N=N)
data <- list(y=y)

# For prediction
missing_ind <- which(is.na(y))
missing_data <- paste0("y[", missing_ind, "]")

# Model formulation
code <- s
# initial values
initsFn <- function() {
  list(mu=rnorm(1), phi=rnorm(1,sd=.2), sigma_PN=runif(1,0,0.5), sigma_OE=runif(1,0,0.5), x=rep(0,N)) # max param. will be modified later
}
set.seed(0)
inits <- initsFn()
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel, control=list(adaptInterval=100), enableWAIC=TRUE)
conf$removeSamplers(varNames)
conf$addSampler(varNames[1:3], type="AF_slice")
conf$addSampler(varNames[4], type="slice") # independent sampling for obs. error
conf$addMonitors("x", "y") # Predictive nodes
Rmcmc <- buildMCMC(conf)
CompiledList <- compileNimble(list(Rmodel, Rmcmc))
Cmodel <- CompiledList[[1]]
Cmcmc <- CompiledList[[2]]
samples <- runMCMC(Cmcmc, niter=32000, nburnin=2000, samplesAsCodaMCMC=T, WAIC=T)
samples$WAIC # -6392.249

param_samples <- samples$samples[,varNames]
missing_samples <- samples$samples[,missing_data]
save(param_samples, file="final_post.RData")
save(missing_samples, file="final_pred.RData")

samplesPlot(param_samples, scale=T)
effectiveSize(param_samples)
acfplot(param_samples,aspect="fill",lag.max=100)
samplesSummary(param_samples)[,c(1,3,4,5)]

# credible intervals for 1942:1953 (Nov.20th 9pm - 21st 9am)
bci_for_plot <- samplesSummary(missing_samples)[9:20,c(1,4,5)]
min_lim <- min(y[1935:1959], na.rm=T)
max_lim <- max(bci_for_plot[,3])

# plotting
png("interpolation.png")
par(mar=c(3,3,2,2)+.1, mgp=c(2,1,0), cex.main=1.6, cex.axis=1.3)
plot(1935:1959, y[1935:1959], xlab="", ylab="", main="Nov.20th, 14:00 - Nov.21st, 14:00",
     ylim=c(min_lim, max_lim), type="b", pch=20)
points(1942:1953, bci_for_plot[,1], col="blue", pch=20)
lines(c(1941,1942), c(y[1941], bci_for_plot[1,1]), col="blue")
lines(c(1953,1954), c(bci_for_plot[12,1], y[1954]), col="blue")
lines(1942:1953, bci_for_plot[,1], col="blue")
lines(1942:1953, bci_for_plot[,2], col="red", lwd=1.5, lty=2)
lines(1942:1953, bci_for_plot[,3], col="red", lwd=1.5, lty=2)
dev.off()


