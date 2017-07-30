rm(list=ls())
dev.off()

source("sampling-functions.R")
library("data.table")
load("mv-stochastic-vol-data.Rdata")

parameters <- macro.parameters
parameters$tau.rho <- 0.001

N.particles <- 100
weights <- rep(1, N.particles)
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

theta.t = sample.theta.prior.par(parameters,
                                  N.particles)
ks = seq(1,N.particles)
theta.t.plus.1 = theta.t

for (tt in seq(2,T)) {
    if (tt == 3) break;

    particle.quantiles[tt,
                       c("log.sigma.x.lower",
                         "log.sigma.x.mean",
                         "log.sigma.x.upper") := as.list(quantile(x = theta.t[1,],
                                                                  probs = c(
                                                                      0.025,
                                                                      0.5,
                                                                      0.975)))]
    
    particle.quantiles[tt,
                       c("log.sigma.y.lower",
                         "log.sigma.y.mean",
                         "log.sigma.y.upper") := as.list(quantile(x = theta.t[2,],
                                                                  probs = c(
                                                                      0.025,
                                                                      0.5,
                                                                      0.975)))]
    
    particle.quantiles[tt,
                       c("rho.tilde.lower",
                         "rho.tilde.mean",
                         "rho.tilde.upper") := as.list(quantile(x = theta.t[3,],
                                                                probs = c(
                                                                    0.025,
                                                                    0.5,
                                                                    0.975)))]
        
    y.t.minus.1 = c(sample.data$x[tt-1], sample.data$y[tt-1])
    y.t = c(sample.data$x[tt], sample.data$y[tt])
    
    ## sampale weights
    theta.tp1.mean = theta.next.mean.par(theta.current = theta.t,
                                         parameters = parameters)
    if (tt < T) {
        theta.tp1.mean[1,] = sample.data$log.sigma.x[tt+1]
        theta.tp1.mean[2,] = sample.data$log.sigma.y[tt+1]
        theta.tp1.mean[3,] = sample.data$rho.tilde[tt+1]
    }
    
    lls = log.likelihood.par(y.t,
                             y.t.minus.1,
                             theta.tp1.mean,
                             theta.t,
                             parameters,
                             N.particles)
    
    probs = exp((lls + log(weights)) - max(lls + log(weights)))
    for (m in seq(1,N.particles)) {
        k = sample(x = seq(1,N.particles),
                   size = 1,
                   prob = probs)
        ks[m] = k
        theta.tp1 = sample.theta(theta.current = theta.t[,k],
                                 parameters = parameters)

        log.new.weight = log.likelihood(y.t = y.t, y.t.minus.1 = y.t.minus.1,
                                        theta.t.plus.1 = theta.tp1,
                                        theta.t = theta.t[,k],
                                        parameters = parameters) -
            log.likelihood(y.t = y.t,
                           y.t.minus.1 = y.t.minus.1,
                           theta.t.plus.1 = theta.tp1.mean[,k],
                           theta.t = theta.t[,k],
                           parameters = parameters)
        
        theta.t.plus.1[,m] = theta.tp1
        weights[m] = exp(log.new.weight)
    }


    plot(density(theta.t[2,]),
         xlim = c(min(c(theta.t[2,],theta.t.plus.1[2,])),
                  max(c(theta.t[2,],theta.t.plus.1[2,]))))
    lines(density(theta.t.plus.1[2,]), col = "red")
    abline(v = sample.data$log.sigma.y[tt]);

    theta.t = theta.t.plus.1
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




