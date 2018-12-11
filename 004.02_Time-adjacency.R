#
# Creating seasonal (annual) time series through SARIMA model
# How well can eigenfilters perform in modeling them? 
#

# (Parallel computing version)

library(tidyverse)
library(GGally)
library(forecast)
library(Matrix)
library(RSpectra)
library(parallel)
source("000_utils.R")

# eigenfilter model fits
dataList_fn <- function(l, ...) lapply(l, function(ts) cbind.data.frame(y=ts, ...))
filterModel_fn <- function(df, scope=NULL, n=100) {
  ## BIC criterion stepwise fitting
  fullModel <- lm(y~., data=df)
  result <- step(fullModel, scope=scope, direction="backward", k=log(n), trace=0)
  result
} 
# Ljung-Box test fns
ljungbox_fn <- function(m) Box.test(m$residuals, type="Ljung-Box")$p.value
#--

# 1. Create a SARIMA model to generate simulations (annual seasonality is assumed.)
simSARIMA_fn <- function(n=300, s=12, order, seasonal, coefs, ...) {
  ## PRE: order, seasonal are passed onto forecast::Arima (see ?arima)
  ## coefs is a vector containing true model coefficients
  ## POST: list of n simulated series with perturbed errors
  model <- Arima(ts(rnorm(100), frequency=s), order, seasonal, fixed=coefs, ...)
  result <- list()
  for (i in 1:n) result[[i]] <- simulate(model, nsim=n) + rnorm(n)
  result
}


# Simulate 100 instantiations for different models
# Some have exogenous variables
set.seed(12345)
annualArima1_l <- simSARIMA_fn(order=c(1,0,0), seasonal=c(1,0,0), 
                               coefs=c(0.1, 0.2), include.mean=F)
annualArima2_l <- simSARIMA_fn(order=c(1,1,1), seasonal=c(1,0,1),
                               coefs=c(phi=0.3, theta=0.1, Phi=-0.2, Theta=-0.2))
annualArima3_l <- simSARIMA_fn(order=c(1,1,2), seasonal=c(2,0,0),
                               coefs=c(phi=-0.2, theta1=0.1, theta2=-0.05, Phi1=-0.1, Phi2=0.02))
exVar <- rnorm(300, 10, 3)
annualArimaEx1_l <- lapply(annualArima2_l, function(x) x + exVar * 2)
annualArimaEx2_l <- lapply(annualArima2_l, function(x) x + exVar * 2 + exVar^2 * 0.5)

# specifying adjacency matrices
c_A1 <- binAdjTime_fn(300)
c_A2 <- binAdjTime_fn(300, lags=c(1,12))
c_A3 <- binAdjTime_fn(300, lags=c(1,2,12,24))

## different covar. matrices
c_X1 <- as.matrix(1:300)
c_X2 <- cbind(exVar, c_X1)
c_X3 <- cbind(exVarsq=exVar^2, c_X2)

# Modeling through two different models
# First 90 vectors will be used (instability due to collinearity)
myCl <- makeCluster(detectCores() - 1, type="FORK")
c_basesList <- clusterMap(myCl, moranEigen_fn, list(c_X1, c_X1, c_X1, c_X2, c_X2, c_X3, c_X3),
                          list(c_A1, c_A2, c_A3, c_A1, c_A2, c_A1, c_A2), 
                          as.list(rep(150, 7)), as.list(rep(T, 7)))
# c_bases1 <- moranEigen_fn(k=500, c_X1, c_A1, attractive=T)
# c_bases2 <- moranEigen_fn(k=500, c_X1, c_A2, attractive=T)
# c_bases3 <- moranEigen_fn(k=500, c_X2, c_A1, attractive=T)
# c_bases4 <- moranEigen_fn(k=500, c_X3, c_A1, attractive=T)

# Modeling through two different models
# First 90 vectors will be used (instability due to collinearity)
c_dataList1 <- dataList_fn(annualArima1_l, c_X1, c_basesList[[1]])
c_dataList2 <- dataList_fn(annualArima1_l, c_X1, c_basesList[[2]])
c_dataList3 <- dataList_fn(annualArima2_l, c_X1, c_basesList[[1]])
c_dataList4 <- dataList_fn(annualArima2_l, c_X1, c_basesList[[2]])
c_dataList5 <- dataList_fn(annualArima3_l, c_X1, c_basesList[[1]])
c_dataList6 <- dataList_fn(annualArima3_l, c_X1, c_basesList[[3]])
c_dataList7 <- dataList_fn(annualArimaEx1_l, c_X2, c_basesList[[4]])
c_dataList8 <- dataList_fn(annualArimaEx1_l, c_X2, c_basesList[[5]])
c_dataList9 <- dataList_fn(annualArimaEx2_l, c_X3, c_basesList[[6]])
c_dataList10 <- dataList_fn(annualArimaEx2_l, c_X3, c_basesList[[7]])
  
# Reference: SARIMA models using the "true" model orders
c_sarimaList1 <- parLapply(myCl, annualArima1_l, 
  function(ts) auto.arima(ts, d=1, max.p=1, max.q=1, max.P=1, max.Q=1, ic="bic"))
c_sarimaList2 <- parLapply(myCl, annualArima2_l,
  function(ts) auto.arima(ts, d=1, max.p=1, max.q=1, max.P=1, max.Q=1, ic="bic"))
c_sarimaList3 <- parLapply(myCl, annualArima3_l,
  function(ts) auto.arima(ts, d=1, max.p=2, max.q=2, max.P=2, max.Q=2, ic="bic"))
c_sarimaList4 <- parLapply(myCl, annualArimaEx1_l,
  function(ts) auto.arima(ts, d=1, max.p=1, max.q=1, max.P=1, max.Q=1, xreg=exVar, ic="bic"))
c_sarimaList5 <- parLapply(myCl, annualArimaEx2_l,
  function(ts) auto.arima(ts, d=1, max.p=2, max.q=2, max.P=2, max.Q=2, xreg=cbind(exVar^2, exVar), ic="bic"))
  
# eigenfilter models

c_filterList1 <- parLapply(myCl, c_dataList1, filterModel_fn)
c_filterList2 <- parLapply(myCl, c_dataList2, filterModel_fn)
c_filterList3 <- parLapply(myCl, c_dataList3, filterModel_fn)
c_filterList4 <- parLapply(myCl, c_dataList4, filterModel_fn)
c_filterList5 <- parLapply(myCl, c_dataList5, filterModel_fn)
c_filterList6 <- parLapply(myCl, c_dataList6, filterModel_fn)
c_filterList7 <- parLapply(myCl, c_dataList7, function(df) filterModel_fn(df, scope=list(lower= ~ exVar)))
c_filterList8 <- parLapply(myCl, c_dataList8, function(df) filterModel_fn(df, scope=list(lower= ~ exVar)))
c_filterList9 <- parLapply(myCl, c_dataList9, function(df) filterModel_fn(df, scope=list(lower= ~ exVarsq + exVar)))
c_filterList10 <- parLapply(myCl, c_dataList10, function(df) filterModel_fn(df, scope=list(lower= ~ exVarsq + exVar)))

stopCluster(myCl)

# 3. Plots: Comparing BIC and Ljung-Box test p-values
bicCompare_fn <- function(l1, l2) {
  table <- gather(data.frame(
    "Sarima" = sapply(l1, BIC), "SeasonFilter"=sapply(l2, BIC)
    ), key="Type", value="BIC"
  )
  result_g <- ggplot(data=table, aes(Type, BIC)) +
    geom_boxplot(fill=c("pink", "lightblue")) +
    theme_bw()
  result_g
}

pvalCompare_fn <- function(l1, l2) {
  table <- gather(data.frame(
    "Sarima" = sapply(l1, ljungbox_fn), "SeasonFilter"=sapply(l2, ljungbox_fn)
  ), key="Type", value="p_Value"
  )
  result_g <- ggplot(data=table, aes(Type, p_Value)) +
    geom_boxplot(fill=c("pink", "lightblue")) +
    theme_bw()
  result_g
}

sarimaLists <- list(c_sarimaList1, c_sarimaList2, c_sarimaList3, c_sarimaList4, c_sarimaList5)
filterLists1 <- list(c_filterList1, c_filterList3, c_filterList5, c_filterList7, c_filterList9)
filterLists2 <- list(c_filterList2, c_filterList4, c_filterList6, c_filterList8, c_filterList10)
bicComparisons1_l <- purrr::map2(sarimaLists, filterLists1, bicCompare_fn)
bicComparisons2_l <- purrr::map2(sarimaLists, filterLists2, bicCompare_fn)
pvalComparisons1_l <- purrr::map2(sarimaLists, filterLists1, pvalCompare_fn)
pvalComparisons2_l <- purrr::map2(sarimaLists, filterLists2, pvalCompare_fn)

bicComparisons1_gm <- ggmatrix(bicComparisons1_l, 2, 3, title="BIC comparisons for different seasonal data.")
pvalComparisons1_gm <- ggmatrix(pvalComparisons1_l, 2, 3, title="Ljung-Box test results for different seasonal data.")
bicComparisons2_gm <- ggmatrix(bicComparisons2_l, 2, 3, title="BIC comparisons for different seasonal data.")
pvalComparisons2_gm <- ggmatrix(pvalComparisons2_l, 2, 3, title="Ljung-Box test results for different seasonal data.")

ggsave(bicComparisons1_gm, filename = "004.15_BIC-Griffith.pdf", device=pdf)
ggsave(pvalComparisons1_gm, filename = "004.16_p-value-Griffith.pdf", device=pdf)
ggsave(bicComparisons2_gm, filename = "004.17_BIC-Modified.pdf", device=pdf)
ggsave(pvalComparisons2_gm, filename = "004.18_p-value-Modified.pdf", device=pdf)
