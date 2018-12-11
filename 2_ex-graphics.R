setwd("Bayesian_class/project/")
load("data.RData")

# Exploratory plots to be included in the poster

library(dplyr)
library(ggplot2)
library(reshape2)

# Time series trends
differenced <- diff(dust$log_dust)
png("ts-plots.png")
par(mar=c(2,3,1,1)+.1, mgp=c(2,1,0), mfrow=c(2,1))
plot.ts(dust$log_dust, xlab="", ylab="Log levels", main="")
plot.ts(differenced, ylab="Log levels")
dev.off()

# Correlation plot between naive mcmc fit parameters
load("naive_samples.RData")
varNames <- c("mu", "phi", "sigma_PN", "sigma_OE")
my_samples <- lapply(samples$samples, function(x) x[,varNames])
all_chains <- rbind(my_samples$chain1, my_samples$chain2, my_samples$chain3)
corMat_df <- all_chains %>% 
  as.matrix() %>% 
  cor() %>% abs() %>% melt()

g1 <- ggplot(data=corMat_df, aes(Var1, Var2)) +
  geom_tile(aes(fill=value)) + 
  labs(x="", y="", fill="Absolute\n correlation") +
  scale_fill_distiller(palette="YlGnBu", direction=1) +
  theme(axis.text.x = element_text(angle=45, vjust=.6))
ggsave("correlation_rw.png", g1, width=6.21, height=4.58, device="png")

# Correlation plot after independent sampling of observation error
# load("ind_OE.RData")
my_samples <- lapply(samples2$samples, function(x) x[,varNames])
all_chains <- rbind(my_samples$chain1, my_samples$chain2, my_samples$chain3)
corMat_df <- all_chains %>% 
  as.matrix() %>% 
  cor() %>% abs() %>% melt()

g2 <- ggplot(data=corMat_df, aes(Var1, Var2)) +
  geom_tile(aes(fill=value)) + 
  labs(x="", y="", fill="Absolute\n correlation") +
  scale_fill_distiller(palette="YlGnBu", direction=1) +
  theme(axis.text.x = element_text(angle=45, vjust=.6, colour="black", size=11),
        axis.text.y = element_text(colour="black", size=11),
        legend.margin=margin(0,0,0,0,"pt"))
ggsave("correlation_ind_sampling.png", g2, width=6.21, height=4.58, device="png")
