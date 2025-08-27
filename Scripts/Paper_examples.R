library(DiscreteDists)

# Figure ------------------------------------------------------------------

# Load the necessary library
library(tidyverse)

sigmas <- c(0.1, 1, 1.9)

sigma_dat <- map_df(sigmas, ~ tibble(
  sigma = paste(.),
  x = 0:13,
  y = dHYPERPO(x=0:13, mu=5.5, sigma=.)
))

# Use ggplot2 to plot
p1 <- ggplot(sigma_dat, aes(x, y, color=factor(sigma))) + 
  geom_line() +
  geom_point(data = sigma_dat, 
             aes(x, y, color=factor(sigma))) +
  labs(color = expression(sigma)) +
  ylab("Probability") +
  xlab("y") + 
  ggtitle(expression(paste("Probability mass function for HYPERPO(", mu, "=5.5, ", sigma, ")"))) +
  theme_minimal()

p1

ggsave(plot=p1,
       filename="Figures/plot_hyperpo.pdf",
       width=8, height=6)


# Examples for HYPERPO - first parameterization ---------------------------

library(DiscreteDists)

dHYPERPO(x=0:10, mu=5.5, sigma=0.1)

dHYPERPO(x=0:10, mu=5.5, sigma=1)
dpois(x=0:10, lambda=5.5)

set.seed(1234)
y <- rHYPERPO(n=500, mu=5.5, sigma=0.1)
y[1:15]

# Empirical mean and variance
mean(y)
var(y)

# Mean and variance
mean_var_hp(mu=5.5, sigma=0.1)

# Estimating the paramertes
library(gamlss)
mod1 <- gamlss(y ~ 1, family=HYPERPO)
summary(mod1)

# A function to simulate a data set with Y ~ HYPERPO2
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(1.21 - 3 * x1)
  sigma <- exp(1.26 - 2 * x2)
  y <- rHYPERPO(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(1234)
dataset <- gendat(n=200)

mod2 <- gamlss(y~x1, sigma.fo=~x2, family=HYPERPO, 
               control=gamlss.control(n.cyc=50),
               data=dataset)

summary(mod2)

