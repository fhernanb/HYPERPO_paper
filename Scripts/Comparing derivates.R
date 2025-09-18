
# Derivada con respecto a mu

dldm_manual = function(y, mu, sigma) {
  dldm <- dldm_hyperpo_cpp(y, mu, sigma)
  dldm
}

dldm_compu = function(y, mu, sigma) {
  dm   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                theta="mu",
                                delta=0.00001)
  dldm <- as.vector(attr(dm, "gradient"))
  dldm
}

y <- 2:60
mu <- 2:60
sigma <- 2:60

dldm_manual(y, mu, sigma)
dldm_compu(y, mu, sigma)

library(microbenchmark)

res <- microbenchmark(dldm_manual(y, mu, sigma),
                      dldm_compu(y, mu, sigma),
                      times=100)
plot(res)

# Derivada con respecto a sigma

dldd_manual = function(y, mu, sigma) {
  dldd <- dldd_hyperpo_cpp(y, mu, sigma)
  dldd
}

dldd_compu = function(y, mu, sigma) {
  dd   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                theta="sigma",
                                delta=0.00001)
  dldd <- as.vector(attr(dd, "gradient"))
  dldd
}

y <- 3:5
mu <- 3:5
sigma <- 3:5

dldd_manual(y, mu, sigma)
dldd_compu(y, mu, sigma)


library(microbenchmark)

res <- microbenchmark(dldd_manual(y, mu, sigma),
                      dldd_compu(y, mu, sigma),
                      times=100)
plot(res)


