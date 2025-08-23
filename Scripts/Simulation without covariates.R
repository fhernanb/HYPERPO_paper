library(gamlss)
library(DiscreteDists)

# En la distribucion HYPERPO
# mu > 0,     por tanto usaremos funcion de enlace log
# sigma > 0,  por tanto usaremos funcion de enlace log

# Useful functions to the simulation study --------------------------------

# Funcion para obtener mu_hat y sigma_hat para un valor fijo de n
simul_one <- function(size, true_mu, true_sigma) {
  y <- rHYPERPO(n=size, mu=true_mu, sigma=true_sigma)
  mod <- NULL
  mod <- try(gamlss(y~1, sigma.fo=~1, family="HYPERPO",
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

# Instruction to simulate type I
parSim(
  ### SIMULATION CONDITIONS
  
  mu = c(3, 7),
  sigma = c(0.5, 1, 1.5),
  n = c(20, 40, 80, 100, 200, 500, 700, 1000),
  
  reps = 100,                    # Repetitions
  write = TRUE,                  # Writing to a file
  name = "simul_sin_cov",     # Name of the file
  nCores = 1,                    # Number of cores to use
  
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

# Plots -------------------------------------------------------------------
dt1 <- read.table("simul01.txt", header=TRUE)
dt2 <- read.table("simul02.txt", header=TRUE)
dt3 <- read.table("simul03.txt", header=TRUE)
dt4 <- read.table("simul04.txt", header=TRUE)
dt5 <- read.table("simul05.txt", header=TRUE)
dt6 <- read.table("simul06.txt", header=TRUE)
dt7 <- read.table("simul07.txt", header=TRUE)

dt99 <- read.table("simul99.txt", header=TRUE)
dt <- rbind(dt1, dt2, dt3, dt4, dt5, dt6, dt7, dt99)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# Numero de observaciones por cada n
num <- dt %>% group_by(n, mu, sigma) %>% count()
num
mean(num$nn)
min(num$nn)

# Para crear las cantidades necesarias para las metricas
res <- dt %>% 
  group_by(n, mu, sigma) %>% 
  summarise(mean_mu=mean(mu_hat, trim=0.1),
            mean_sigma=mean(sigma_hat, trim=0.1),
            mse_mu=mean((mu - mu_hat)^2, trim=0.1),
            mse_sigma=mean((sigma - sigma_hat)^2, trim=0.1),
            nrep=n())

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
    ylim(true_mu - 0.1, true_mu + 0.1)
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
  ggsave(filename=paste0("case1_mean_true_mu=", true_mu, ".pdf"), 
         plot=the_mean, 
         width=6, height=8)
}

figs_case1_mean(res, true_mu=3, k=0.1)
figs_case1_mean(res, true_mu=7, k=0.2)

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
  ggsave(filename=paste0("case1_mse_true_mu=", true_mu, ".pdf"), 
         plot=the_mse, 
         width=6, height=8)
}

figs_case1_mse(res, true_mu=3)
figs_case1_mse(res, true_mu=7)

