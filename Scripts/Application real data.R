
# Application 1 from Saez-Castillo (2013) - Takeover bids
library(Ecdat)
library(DiscreteDists)
library(gamlss)
library(dplyr)

# To explore the response variable
pdf("Figures/hist_number_bids.pdf", width = 8, height = 4)
barplot(table(Bids$numbids), las=1, col=gray(0.9),
        xlab="Number of bids", ylab="Frequency")
dev.off()

# Mean and variance
Bids %>% summarise(media_bids = mean(numbids),
                   var_bids = var(numbids),
                   ratio = var(numbids)/mean(numbids))


# Models

# Poisson model reported in table 1 in Saez-Castillo (2013)
mod1 <- gamlss(numbids~leglrest+rearest+finrest+whtknght+bidprem+insthold+
                 size+I(size^2)+regulatn, 
               family=PO, 
               data=Bids)

summary(mod1)
BIC(mod1)
Rsq(mod1)

# Hyper-Poisson model similar to the reported in table 1 in Saez-Castillo (2013)
mod2 <- NULL
mod2 <- gamlss(numbids~leglrest+rearest+finrest+whtknght+
                 bidprem+insthold+size+I(size^2)+regulatn, 
               family=HYPERPO2, 
               data=Bids)

summary(mod2)
BIC(mod2)
Rsq(mod2)

# Hyper-Poisson model similar to the reported in table 2 in Saez-Castillo (2013)
mod3 <- NULL
mod3 <- gamlss(numbids~leglrest+rearest+finrest+whtknght+
                 bidprem+insthold+
                 size+I(size^2)+regulatn,
               sigma.fo=~rearest+finrest+bidprem+regulatn,
               family=HYPERPO2, 
               data=Bids,
               control=gamlss.control(n.cyc=150, trace=TRUE))

summary(mod3)
BIC(mod3)
Rsq(mod3)
logLik(mod3)

wp(mod3)

# Residual analysis
library(car)

res_mod3 <- resid(mod3)

pdf("Figures/worm_rqr_plot.pdf", width = 8, height = 4)
par(mfrow=c(1, 2))
wp(mod3, main="Worm plot")
title(main="Worm plot")
qqPlot(res_mod3, dist="norm", mean=0, sd=1,
       main='RQR plot', col.lines="gray",
       ylab='Randomized quantile residuals', las=1)
dev.off()


# Creating latex table for model 3
t3 <- summary(mod3)

library(xtable)
xtable(t3, digits=4)

