
library(DiscreteDists)
library(gamlss)

Sys.time()

gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  x11 <- runif(n)
  x22 <- runif(n)
  mu    <- exp(1.21 - 3 * x1 + 2 * x11)
  sigma <- exp(1.26 - 2 * x2 + 1.5 * x22)
  y <- rHYPERPO2(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2, x11=x11, x22=x22, mu=mu, sigma=sigma)
}
set.seed(1234)

mu <- NULL
sigma <- NULL

for (i in 1:10){
  dataset <- gendat(n=200)
  mod <- gamlss(y~x1+x11, sigma.fo=~x2+x22,
                family=HYPERPO2, data=dataset,
                control=gamlss.control(n.cyc=2500, trace=TRUE))
  mu=rbind(mu,mod$mu.coefficients)
  sigma=rbind(sigma,mod$sigma.coefficients)
}

mu
sigma

colMeans(mu)
colMeans(sigma)

Sys.time()
