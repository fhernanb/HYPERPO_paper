
library(gamlss)
library(DiscreteDists)

# In the HYPERPO distribution
# mu > 0, so we"ll use the log link function
# sigma > 0, so we"ll use the log link function

# The parameters ----------------------------------------------------------
# The next values correspond to the true betas in the regression model

true_b0_mu <- 1.21   # intercept for mu
true_b1_mu <- -3.0   # slope for mu
true_b2_mu <- 2      # slope for mu
true_g0_si <- 1.26   # intercept for sigma
true_g1_si <- -2.0   # slope for sigma
true_g2_si <- 1.5    # slope for sigma

# Useful functions to the simulation study --------------------------------

gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  x11 <- runif(n)
  x22 <- runif(n)
  mu    <- exp(1.21 - 3 * x1 + 2 * x11)
  sigma <- exp(1.26 - 2 * x2 + 1.5 * x22)
  y <- rHYPERPO2(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2, x11=x11, x22=x22)
}

# Function to obtain beta_hat for a fixed value of n
simul_one <- function(size) {
  dat <- gendat(size)
  mod <- NULL
  mod <- try(gamlss(y~x1+x11, sigma.fo=~x2+x22, family="HYPERPO2", data=dat,
                    control=gamlss.control(n.cyc=2500, trace=FALSE)),
             silent=TRUE)
  if (class(mod)[1] == "try-error")
    res <- rep(NA, 6)
  else
    res <- c(coef(mod, what="mu"), coef(mod, what="sigma"))
  res
}

# To perform the simulation -----------------------------------------------

library(parSim)

# Instruction to simulate 
parSim(
  ### SIMULATION CONDITIONS
  n = c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
  
  reps = 500,                              # Repetitions
  write = TRUE,                           # Writing to a file
  name = "Simulations/case2_10", # Name of the file
  nCores = 1,                             # Number of cores to use
  
  expression = {
    res <- simul_one(size=n)
    
    # Results list:
    Results <- list(
      b0_mu_hat = res[1],
      b1_mu_hat = res[2],
      b2_mu_hat = res[3],
      g0_si_hat = res[4],
      g1_si_hat = res[5],
      g2_si_hat = res[6]
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------
archivos <- list.files(pattern = "^case2.*\\.txt$", 
                       path="Simulations",
                       full.names = TRUE)

archivos

lista_datos <- lapply(archivos, read.table, header = TRUE, 
                      sep = "", stringsAsFactors = FALSE)
datos <- do.call(rbind, lista_datos)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# To obtain the metrics mean and MSE
trim <- 0.0

res <- datos %>% 
  select(!errorMessage) %>%
  drop_na() %>% 
  group_by(n) %>% 
  summarise(mean_b0_mu=mean(b0_mu_hat, trim=trim),
            mean_b1_mu=mean(b1_mu_hat, trim=trim),
            mean_b2_mu=mean(b2_mu_hat, trim=trim),
            mean_g0_si=mean(g0_si_hat, trim=trim),
            mean_g1_si=mean(g1_si_hat, trim=trim),
            mean_g2_si=mean(g2_si_hat, trim=trim),
            mse_b0_mu=mean((true_b0_mu - b0_mu_hat)^2, trim=trim), 
            mse_b1_mu=mean((true_b1_mu - b1_mu_hat)^2, trim=trim),
            mse_b2_mu=mean((true_b2_mu - b2_mu_hat)^2, trim=trim),
            mse_g0_si=mean((true_g0_si - g0_si_hat)^2, trim=trim), 
            mse_g1_si=mean((true_g1_si - g1_si_hat)^2, trim=trim),
            mse_g2_si=mean((true_g2_si - g2_si_hat)^2, trim=trim),
            nobs=n())

res

# Mean -----------------------------------------------------
k <- 0.05

p1 <- ggplot(data=res, aes(x=n, y=mean_b0_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(beta)[0]), 
       title=expression("Mean for the intercept in" ~ mu)) +
  ylim(true_b0_mu-k, true_b0_mu+k) +
  geom_line(y=true_b0_mu, col="red", lty="dashed")

p2 <- ggplot(data=res, aes(x=n, y=mean_b1_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(beta)[1]),
       title=expression("Mean for the first slope in" ~ mu)) +
  ylim(true_b1_mu-k, true_b1_mu+k) +
  geom_line(y=true_b1_mu, col="red", lty="dashed")

p3 <- ggplot(data=res, aes(x=n, y=mean_b2_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(beta)[2]),
       title=expression("Mean for the second slope in" ~ mu)) +
  ylim(true_b2_mu-k, true_b2_mu+k) +
  geom_line(y=true_b2_mu, col="red", lty="dashed")

k <- 0.5

p4 <- ggplot(data=res, aes(x=n, y=mean_g0_si)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(gamma)[0]),
       title=expression("Mean for the intercept in" ~ sigma)) +
  ylim(true_g0_si-k, true_g0_si+k) +
  geom_line(y=true_g0_si, col="red", lty="dashed")

p5 <- ggplot(data=res, aes(x=n, y=mean_g1_si)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(gamma)[1]),
       title=expression("Mean for the first slope in" ~ sigma)) +
  ylim(true_g1_si-k, true_g1_si+k) +
  geom_line(y=true_g1_si, col="red", lty="dashed")

p6 <- ggplot(data=res, aes(x=n, y=mean_g2_si)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(gamma)[2]),
       title=expression("Mean for the second slope in" ~ sigma)) +
  ylim(true_g2_si-k, true_g2_si+k) +
  geom_line(y=true_g2_si, col="red", lty="dashed")

mean2 <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2, ncol=3)
mean2
ggsave(filename="Figures/case2_mean.pdf", 
       plot=mean2, 
       width=10, height=8)

# MSE -----------------------------------------------------
p1 <- ggplot(data=res, aes(x=n, y=mse_b0_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(beta)[0]),
       title=expression("MSE for the intercept in" ~ mu))

p2 <- ggplot(data=res, aes(x=n, y=mse_b1_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(beta)[1]),
       title=expression("MSE for the first slope in" ~ mu))

p3 <- ggplot(data=res, aes(x=n, y=mse_b2_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(beta)[2]),
       title=expression("MSE for the first slope in" ~ mu))


p4 <- ggplot(data=res, aes(x=n, y=mse_g0_si)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(gamma)[0]),
       title=expression("MSE for the intercept in" ~ sigma))

p5 <- ggplot(data=res, aes(x=n, y=mse_g1_si)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(gamma)[1]),
       title=expression("MSE for the first slope in" ~ sigma))

p6 <- ggplot(data=res, aes(x=n, y=mse_g2_si)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(gamma)[2]),
       title=expression("MSE for the second slope in" ~ sigma))


mse2 <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2, ncol=3)
mse2
ggsave(filename="Figures/case2_mse.pdf", 
       plot=mse2, 
       width=10, height=8)

