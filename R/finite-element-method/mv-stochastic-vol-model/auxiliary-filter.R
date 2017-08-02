rm(list=ls())
dev.off()

source("sampling-functions.R")
library("data.table")
load("mv-stochastic-vol-data.Rdata")

parameters <- macro.parameters
parameters$tau.rho <- 0.01

N.particles <- 500
log.weights <- rep(0, N.particles)
T <- length(sample.data$y)
## Start counting at 1: y_1, y_2, ....
## Start sampling at 2: theta_2, theta_3, ...
particle.quantiles <- data.table(matrix(0, nrow = T, ncol = 9))
setnames(x = particle.quantiles,
         old = names(particle.quantiles),
         new = c("log.sigma.x.lower",
                 "log.sigma.x.mean",
                 "log.sigma.x.upper",
                 "log.sigma.y.lower",
                 "log.sigma.y.mean",
                 "log.sigma.y.upper",
                 "rho.tilde.lower",
                 "rho.tilde.mean",
                 "rho.tilde.upper"))

theta.tm1 = sample.theta.prior.par(parameters,
                                   N.particles)
theta.t = theta.tm1
ks = seq(1,N.particles)
for (tt in seq(3,T)) {
    particle.quantiles[tt-1,
                       c("log.sigma.x.lower",
                         "log.sigma.x.mean",
                         "log.sigma.x.upper") := as.list(quantile(x = theta.tm1[1,],
                                                                  probs = c(
                                                                      0.025,
                                                                      0.5,
                                                                      0.975)))]
    
    particle.quantiles[tt-1,
                       c("log.sigma.y.lower",
                         "log.sigma.y.mean",
                         "log.sigma.y.upper") := as.list(quantile(x = theta.tm1[2,],
                                                                  probs = c(
                                                                      0.025,
                                                                      0.5,
                                                                      0.975)))]
    
    particle.quantiles[tt-1,
                       c("rho.tilde.lower",
                         "rho.tilde.mean",
                         "rho.tilde.upper") := as.list(quantile(x = theta.tm1[3,],
                                                                probs = c(
                                                                    0.025,
                                                                    0.5,
                                                                    0.975)))]

    y.t.minus.2 = c(sample.data$x[tt-2], sample.data$y[tt-2])
    y.t.minus.1 = c(sample.data$x[tt-1], sample.data$y[tt-1])
    y.t = c(sample.data$x[tt], sample.data$y[tt])
    
    ## sampale weights
    theta.t.mean = theta.next.mean.par(theta.current = theta.tm1,
                                       y.current = y.t.minus.1,
                                       y.current.m1 = y.t.minus.2,
                                       parameters = parameters)
    
    theta.t.mean[1,] = rep(sample.data$log.sigma.x[tt],
                           N.particles)
    theta.t.mean[2,] = rep(sample.data$log.sigma.y[tt],
                           N.particles)
    theta.t.mean[3,] = rep(sample.data$rho.tilde[tt],
                           N.particles)

    lls = log.likelihood.par(y.t = y.t,
                             y.t.minus.1 = y.t.minus.1,
                             theta.t = theta.t.mean,
                             parameters = parameters,
                             N.particles = N.particles)
    
    probs = exp( lls + log.weights - max(lls + log.weights) )
    for (m in seq(1,N.particles)) {
        k = sample(x = seq(1,N.particles),
                   size = 1,
                   prob = probs)
        ks[m] = k
        theta.t[,m] = sample.theta(theta.current = theta.tm1[,k],
                                   y.current = y.t.minus.1,
                                   y.current.m1 = y.t.minus.2,
                                   parameters = parameters)


        log.new.weight =
            log.likelihood(y.t = y.t, y.t.minus.1 = y.t.minus.1,
                           theta.t = theta.t[,m],
                           parameters = parameters) -
            log.likelihood(y.t = y.t,
                           y.t.minus.1 = y.t.minus.1,
                           theta.t = theta.t.mean[,k],
                           parameters = parameters)
        
        log.weights[m] = log.new.weight
    }
    log.weights = log.weights - max(log.weights)

    theta.tm1 = theta.t
    print(tt)
}

par(mfrow=c(2,2))
plot(particle.quantiles[, log.sigma.x.mean], col = "blue")
lines(particle.quantiles[, log.sigma.x.lower], col = "blue")
lines(particle.quantiles[, log.sigma.x.upper], col = "blue")
lines(sample.data$log.sigma.x, col = "red")

plot(particle.quantiles[, log.sigma.y.mean], col = "blue")
lines(particle.quantiles[, log.sigma.y.lower], col = "blue")
lines(particle.quantiles[, log.sigma.y.upper], col = "blue")
lines(sample.data$log.sigma.y, col = "red")

plot(particle.quantiles[, rho.tilde.mean], col = "blue")
lines(particle.quantiles[, rho.tilde.lower], col = "blue")
lines(particle.quantiles[, rho.tilde.upper], col = "blue")
lines(sample.data$rho.tilde, col = "red")

## hist((log.sigma.x.current.samples))
## hist((log.sigma.y.current.samples))
## hist(rho.tilde)




