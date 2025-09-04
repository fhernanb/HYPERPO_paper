
library(gamlss)
library(DiscreteDists)

# In the HYPERPO distribution
# mu > 0, so we"ll use the log link function
# sigma > 0, so we"ll use the log link function

# The parameters ----------------------------------------------------------
# The next values correspond to the true betas in the regression model

true_b0_mu <- 1.21   # intercept for mu
true_b1_mu <- -3.0   # slope for mu
true_g0_si <- 1.26   # intercept for sigma
true_g1_si <- -2.0   # slope for sigma

# Useful functions to the simulation study --------------------------------

# Function to obtain beta_hat for a fixed value of n
simul_one <- function(size) {
  x1 <- runif(n=size)
  x2 <- runif(n=size)
  mu    <- exp(true_b0_mu + true_b1_mu * x1)
  sigma <- exp(true_g0_si + true_g1_si * x2)
  y <- rHYPERPO(n=size, mu=mu, sigma=sigma)
  mod <- NULL
  mod <- try(gamlss(y~x1, sigma.fo=~x2, family="HYPERPO",
                    control=gamlss.control(n.cyc=2500, trace=FALSE)),
             silent=TRUE)
  if (class(mod)[1] == "try-error")
    res <- rep(NA, 4)
  else
    res <- c(coef(mod, what="mu"), coef(mod, what="sigma"))
  res
}

# To perform the simulation -----------------------------------------------

library("parSim")

# Instruction to simulate 
parSim(
  ### SIMULATION CONDITIONS
  #n = c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
  n = 50,
  
  reps = 500,                              # Repetitions
  write = TRUE,                           # Writing to a file
  name = "Simulations/sim_HP1_with_cov_14", # Name of the file
  nCores = 1,                             # Number of cores to use
  
  expression = {
    res <- simul_one(size=n)
    
    # Results list:
    Results <- list(
      b0_mu_hat = res[1],
      b1_mu_hat = res[2],
      g0_si_hat = res[3],
      g1_si_hat = res[4]
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------
archivos <- list.files(pattern = "^sim_HP1_with_cov.*\\.txt$", 
                       path="Simulations",
                       full.names = TRUE)
lista_datos <- lapply(archivos, read.table, header = TRUE, 
                      sep = "", stringsAsFactors = FALSE)
datos <- do.call(rbind, lista_datos)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# Numero de NA por cada n
datos |> 
  select(n, b0_mu_hat) |> 
  group_by(n) |> 
  summarise_all(~ sum(is.na(.)))

# Numero de observaciones por cada n
num <- datos %>% group_by(n) %>% count()
num
mean(num$nn)
min(num$nn)

# To obtain the metrics mean and MSE
trim <- 0.20

res <- datos %>% 
  select(!errorMessage) %>%
  drop_na() %>% 
  group_by(n) %>% 
  summarise(mean_b0_mu=mean(b0_mu_hat, trim=trim),
            mean_b1_mu=mean(b1_mu_hat, trim=trim),
            mean_b0_si=mean(g0_si_hat, trim=trim),
            mean_b1_si=mean(g1_si_hat, trim=trim),
            mse_b0_mu=mean((true_b0_mu - b0_mu_hat)^2, trim=trim), 
            mse_b1_mu=mean((true_b1_mu - b1_mu_hat)^2, trim=trim),
            mse_b0_si=mean((true_g0_si - g0_si_hat)^2, trim=trim), 
            mse_b1_si=mean((true_g1_si - g1_si_hat)^2, trim=trim),
            nobs=n())

res

# Mean -----------------------------------------------------
k <- 0.13

p1 <- ggplot(data=res, aes(x=n, y=mean_b0_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(beta)[0]), 
       title=expression("Mean for the intercept in" ~ mu)) +
  ylim(true_b0_mu-k, true_b0_mu+k) +
  geom_line(y=true_b0_mu, col="red", lty="dashed")

p2 <- ggplot(data=res, aes(x=n, y=mean_b1_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(beta)[1]),
       title=expression("Mean for the slope in" ~ mu)) +
  ylim(true_b1_mu-k, true_b1_mu+k) +
  geom_line(y=true_b1_mu, col="red", lty="dashed")

p3 <- ggplot(data=res, aes(x=n, y=mean_b0_si)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(gamma)[0]),
       title=expression("Mean for the intercept in" ~ sigma)) +
  ylim(true_g0_si-k, true_g0_si+k) +
  geom_line(y=true_g0_si, col="red", lty="dashed")

p4 <- ggplot(data=res, aes(x=n, y=mean_b1_si)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(gamma)[1]),
       title=expression("Mean for the slope in" ~ sigma)) +
  ylim(true_g1_si-k, true_g1_si+k) +
  geom_line(y=true_g1_si, col="red", lty="dashed")

mean2 <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
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
       title=expression("MSE for the slope in" ~ mu))

p3 <- ggplot(data=res, aes(x=n, y=mse_b0_si)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(gamma)[0]),
       title=expression("MSE for the intercept in" ~ sigma))

p4 <- ggplot(data=res, aes(x=n, y=mse_b1_si)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(gamma)[1]),
       title=expression("MSE for the slope in" ~ sigma))

mse2 <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
mse2
ggsave(filename="Figures/case2_mse.pdf", 
       plot=mse2, 
       width=10, height=8)

