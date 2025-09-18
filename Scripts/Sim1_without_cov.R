
library(gamlss)
library(DiscreteDists)

# In the HYPERPO distribution
# mu > 0, so we'll use the log link function
# sigma > 0, so we'll use the log link function

# Useful functions to the simulation study --------------------------------

# Function to obtain mu_hat and sigma_hat for a fixed value of n
simul_one <- function(size, true_mu, true_sigma) {
  y <- rHYPERPO2(n=size, mu=true_mu, sigma=true_sigma)
  mod <- NULL
  mod <- try(gamlss(y~1, sigma.fo=~1, family="HYPERPO2",
                control=gamlss.control(n.cyc=2500, trace=FALSE)))
  if (class(mod)[1] == "try-error") {
    res <- c(NA, NA)
  }
  else {
    res <- c(exp(coef(mod, what="mu")),
             exp(coef(mod, what="sigma")))
  }
  res
}

# To perform the simulation -----------------------------------------------

library("parSim")

# Instruction to simulate
parSim(
  ### SIMULATION CONDITIONS
  mu = c(3, 7),
  sigma = c(0.5, 1, 1.5),
  
  n = c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),

  reps = 100,                                 # Repetitions
  write = TRUE,                              # Writing to a file
  name = "Simulations/case1_13", # Name of the file
  nCores = 1,                                # Number of cores to use
  
  expression = {
    res <- simul_one(size=n, true_mu=mu, true_sigma=sigma)
    
    # Results list:
    Results <- list(
      mu_hat = res[1],
      sigma_hat = res[2]
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------
archivos <- list.files(pattern = "^case1.*\\.txt$", 
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

# Number of observations for each n
num <- datos %>% group_by(n, mu, sigma) %>% count()
num
mean(num$nn)
min(num$nn)

# To obtain the metrics mean and MSE
trim <- 0.10
res <- datos %>% 
  group_by(n, mu, sigma) %>% 
  summarise(mean_mu=mean(mu_hat, trim=trim, na.rm=TRUE),
            mean_sigma=mean(sigma_hat, trim=trim, na.rm=TRUE),
            mse_mu=mean((mu - mu_hat)^2, trim=trim, na.rm=TRUE),
            mse_sigma=mean((sigma - sigma_hat)^2, trim=trim, na.rm=TRUE),
            nobs=n()
            )

res

# Mean -----------------------------------------------------

# Case 1

figs_case1_mean <- function(res, true_mu, k) {
  
  sub_res <- filter(res, mu == true_mu & sigma == 0.5)
  
  p1 <- ggplot(data=sub_res, aes(x=n, y=mean_mu)) + 
    geom_line() + 
    labs(x="n", y=expression(hat(mu)),
         title=bquote("Case " ~ mu ~ "=" ~ .(true_mu) ~ "and " ~ sigma ~ "= 0.5")) +
    geom_line(y=true_mu, col="red", lty="dashed")  +
    ylim(true_mu - k, true_mu + k)
  p2 <- ggplot(data=sub_res, aes(x=n, y=mean_sigma)) + 
    geom_line() + 
    labs(x="n", y=expression(hat(sigma)),
         title=expression("")) +
    geom_line(y=0.5, col="red", lty="dashed")  +
    ylim(0.5 - k, 0.5 + k)
  
  sub_res <- filter(res, mu == true_mu & sigma == 1)
  
  p3 <- ggplot(data=sub_res, aes(x=n, y=mean_mu)) + 
    geom_line() + 
    labs(x="n", y=expression(hat(mu)),
         title=bquote("Case " ~ mu ~ "=" ~ .(true_mu) ~ "and " ~ sigma ~ "= 1.0")) +
    geom_line(y=true_mu, col="red", lty="dashed")   +
    ylim(true_mu - k, true_mu + k)
  p4 <- ggplot(data=sub_res, aes(x=n, y=mean_sigma)) + 
    geom_line() + 
    labs(x="n", y=expression(hat(sigma)),
         title=expression("")) +
    geom_line(y=1, col="red", lty="dashed")  +
    ylim(1 - k, 1 + k)
  
  sub_res <- filter(res, mu == true_mu & sigma == 1.5)
  
  p5 <- ggplot(data=sub_res, aes(x=n, y=mean_mu)) + 
    geom_line() + 
    labs(x="n", y=expression(hat(mu)),
         title=bquote("Case " ~ mu ~ "=" ~ .(true_mu) ~ "and " ~ sigma ~ "= 1.5")) +
    geom_line(y=true_mu, col="red", lty="dashed")  +
    ylim(true_mu - k, true_mu + k)
  p6 <- ggplot(data=sub_res, aes(x=n, y=mean_sigma)) + 
    geom_line() + 
    labs(x="n", y=expression(hat(sigma)),
         title=expression("")) +
    geom_line(y=1.5, col="red", lty="dashed")  +
    ylim(1.5 - k, 1.5 + k)
  
  the_mean <- grid.arrange(p1, p2, p3, p4, p5, p6,
                           nrow = 3, ncol=2)
  the_mean
  ggsave(filename=paste0("Figures/case1_mean_true_mu=", true_mu, ".pdf"), 
         plot=the_mean, 
         width=6, height=8)
}

figs_case1_mean(res, true_mu=3, k=0.1)
figs_case1_mean(res, true_mu=7, k=0.5)

# MSE -----------------------------------------------------

figs_case1_mse <- function(res, true_mu) {
  
  sub_res <- filter(res, mu == true_mu & sigma == 0.5)
  
  p1 <- ggplot(data=sub_res, aes(x=n, y=mse_mu)) + 
    geom_line() + 
    labs(x="n", y=expression("MSE for " ~ hat(mu)),
         title=bquote("Case " ~ mu ~ "=" ~ .(true_mu) ~ "and " ~ sigma ~ "= 0.5")) +
    geom_line(y=true_mu, col="red", lty="dashed")
  p2 <- ggplot(data=sub_res, aes(x=n, y=mse_sigma)) + 
    geom_line() + 
    labs(x="n", y=expression("MSE for " ~ hat(sigma)),
         title=expression(""))
  
  sub_res <- filter(res, mu == true_mu & sigma == 1)
  
  p3 <- ggplot(data=sub_res, aes(x=n, y=mse_mu)) + 
    geom_line() + 
    labs(x="n", y=expression("MSE for " ~ hat(mu)),
         title=bquote("Case " ~ mu ~ "=" ~ .(true_mu) ~ "and " ~ sigma ~ "= 1.0")) +
    geom_line(y=true_mu, col="red", lty="dashed")
  p4 <- ggplot(data=sub_res, aes(x=n, y=mse_sigma)) + 
    geom_line() + 
    labs(x="n", y=expression("MSE for " ~ hat(sigma)),
         title=expression(""))
  
  sub_res <- filter(res, mu == true_mu & sigma == 1.5)
  
  p5 <- ggplot(data=sub_res, aes(x=n, y=mse_mu)) + 
    geom_line() + 
    labs(x="n", y=expression("MSE for " ~ hat(mu)),
         title=bquote("Case " ~ mu ~ "=" ~ .(true_mu) ~ "and " ~ sigma ~ "= 1.5")) +
    geom_line(y=true_mu, col="red", lty="dashed")
  p6 <- ggplot(data=sub_res, aes(x=n, y=mse_sigma)) + 
    geom_line() + 
    labs(x="n", y=expression("MSE for " ~ hat(sigma)),
         title=expression(""))
  
  the_mse <- grid.arrange(p1, p2, p3, p4, p5, p6,
                           nrow = 3, ncol=2)
  the_mse
  ggsave(filename=paste0("Figures/case1_mse_true_mu=", true_mu, ".pdf"), 
         plot=the_mse, 
         width=6, height=8)
}

figs_case1_mse(res, true_mu=3)
figs_case1_mse(res, true_mu=7)

