rm(list=ls())

source("sampling-function.R")
load("mv-stochastic-vol-data.Rdata")

parameters <- macro.parameters
parameters$tau.rho <- 0.001

N.particles <- 100
weights <- rep(1, N.particles)
## Start counting at 1: y_1, y_2, ....
## Start sampling at 2: theta_2, theta_3, ...


theta.t <- sample.theta.prior(parameters)
theta.t.plus.1 <- sample.theta(theta.t,
                               parameters)
print(theta.t)
print(theta.t.plus.1)

y.t.minus.1 <- c(sample.data$x[1], sample.data$y[1])
y.t <- c(sample.data$x[2], sample.data$y[2])


current.index <- 2
## going forward in time
for (k in seq(1,N.particles)) {
    means.x <- alpha.x + theta.x*(log.sigma.x.current.samples-alpha.x)
    eta.xs <- (means.x - log.sigma.x.current.samples)/tau.x
    epsilon.xs <- leverage.x.rho*eta.xs + sqrt(1-leverage.x.rho^2)*
        rnorm(n=N.particles, 0, 1)

    means.y <- alpha.y + theta.y*(log.sigma.y.current.samples-alpha.y)
    eta.ys <- (means.y - log.sigma.y.current.samples)/tau.y
    epsilon.ys <- leverage.y.rho*eta.ys + sqrt(1-leverage.y.rho^2)*
        rnorm(n=N.particles, 0, 1)

    rhos <- logit.inv(rho.tilde.current.samples)*2 - 1
    sigma.xs <- exp(log.sigma.x.current.samples)
    sigma.ys <- exp(log.sigma.y.current.samples)
    
    x.innovations <- sqrt(1-rhos^2)*sigma.xs*epsilon.xs+
        rhos*sigma.xs/sigma.ys*epsilon.ys
    y.innovations <- sigma.ys*epsilon.ys
    
    new.weights <- dnorm(
    
}

## hist((log.sigma.x.current.samples))
## hist((log.sigma.y.current.samples))
## hist(rho.tilde)




