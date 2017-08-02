rm(list = ls())
library("data.table")
source("sampling-functions.R")

logit <- function(p) {
    out = log(p/(1-p))
    return (out)
}

T = 2 * 256 * 6.5 * 3600 * 1000 ## one year ms
Delta = 1 * 6.5*3600*1000 ## one day in ms

## Values taken from the microstructure paper

mu.hat = 1.7e-12
theta.hat = 5.6e-10
alpha.hat = -13
tau.hat = sqrt(1.3e-9)

mu.x = 0 ## setting it to zero artificically
mu.y = 0

alpha.x = alpha.hat + 1/2*log(Delta)
alpha.y = alpha.hat + 1/2*log(Delta)

theta.x = exp(-theta.hat * Delta)
theta.y = exp(-theta.hat * Delta)

tau.x = tau.hat * sqrt((1 - exp(-2*theta.hat*Delta))/(2*theta.hat))
tau.y = tau.hat * sqrt((1 - exp(-2*theta.hat*Delta))/(2*theta.hat))
tau.rho = 0.0

leverage.x.rho = 0
leverage.y.rho = 0

macro.parameters <- NULL
macro.parameters$mu.x <- mu.x
macro.parameters$mu.y <- mu.y
macro.parameters$alpha.x <- alpha.x
macro.parameters$alpha.y <- alpha.y

macro.parameters$theta.x <- theta.x
macro.parameters$theta.y <- theta.y

macro.parameters$tau.x <- tau.x
macro.parameters$tau.y <- tau.y
macro.parameters$tau.rho <- tau.rho

macro.parameters$leverage.x.rho <- leverage.x.rho
macro.parameters$leverage.y.rho <- leverage.y.rho

generate.data <- function(T, Delta, mu.x, mu.y,
                          alpha.x, alpha.y, theta.x, theta.y,
                          tau.x, tau.y, tau.rho,
                          leverage.x.rho, leverage.y.rho) {
    N <- ceiling(T/Delta)

    rho.t <- 0
    rho.tilde.t <- logit((rho.t+1)/2)

    x.t <- log(100)
    y.t <- log(100)

    log.sigma.x.t <- alpha.x
    log.sigma.y.t <- alpha.y

    innovations <- data.table(matrix(data=rnorm(n=(N+1)*5),
                                     nrow=N+1,
                                     ncol=5))
    setnames(innovations, seq(1,5), c("epsilon.x", "epsilon.y",
                                      "eta.x", "eta.y", "eta.rho"))

    output <- data.table(x = as.numeric(rep(NA,N)),
                         y = as.numeric(rep(NA,N)),
                         log.sigma.x = as.numeric(rep(NA,N)),
                         log.sigma.y = as.numeric(rep(NA,N)),
                         rho.tilde = as.numeric(rep(NA,N)))


    output[1, c("x",
                "y",
                "log.sigma.x",
                "log.sigma.y",
                "rho.tilde") := as.list(c(x.t,
                                          y.t,
                                          log.sigma.x.t,
                                          log.sigma.y.t,
                                          rho.tilde.t))]
    
    for (i in seq(2,N)) {
        log.sigma.x.tp1 = alpha.x +
            theta.x*(log.sigma.x.t - alpha.x) + tau.x*innovations[i-1,eta.x]

        log.sigma.y.tp1 = alpha.y +
            theta.y*(log.sigma.y.t - alpha.y) + tau.y*innovations[i-1,eta.y]

        rho.tilde.tp1 = rho.tilde.t + tau.rho*innovations[i-1,eta.rho]
        rho.tp1 = 2*logit.inv(rho.tilde.tp1) - 1
        
        dx <- exp(log.sigma.x.tp1)*innovations[i,epsilon.x]
        dy <- exp(log.sigma.y.tp1)*innovations[i,epsilon.y]
        
        x.tp1 <- x.t + mu.x + dx
        y.tp1 <- y.t + mu.x + dy


        output[i, c("x","y") := as.list(c(x.tp1, y.tp1))]
        
        output[i, c("log.sigma.x",
                      "log.sigma.y",
                      "rho.tilde") := as.list(c(log.sigma.x.tp1,
                                                log.sigma.y.tp1,
                                                rho.tilde.tp1))]

        if (i %% 1 == 0) {
            lls = sapply(X = seq(-6,-3,length.out = 100),
                         function(x) {
                             dnorm(x = x.tp1,
                                   mean = x.t + mu.x,
                                   sd = exp(x),
                                   log = TRUE) +
                                 dnorm(x = y.tp1,
                                       mean = y.t + mu.x,
                                       sd = exp(log.sigma.y.tp1),
                                       log = TRUE)})

            dnorm(x = x.tp1,
                  mean = x.t + mu.x,
                  sd = exp(log.sigma.x.tp1),
                  log = TRUE)
            dnorm(x = y.tp1,
                  mean = y.t + mu.x,
                  sd = exp(log.sigma.y.tp1),
                  log = TRUE)
            
            lls.2 = log.likelihood.par(y.t = c(x.tp1, y.tp1),
                                       y.t.minus.1 = c(x.t, y.t),
                                       theta.t = rbind(seq(-6,-3,length.out = 100),
                                                       rep(log.sigma.y.tp1, 100),
                                                       rep(rho.tilde.tp1, 100)),
                                       parameters = macro.parameters,
                                       N.particles = 100)
                        
            plot(seq(-6,-3,length.out=100), lls)
            abline(v = log.sigma.x.tp1, col = "red", lwd =2)
            lines(seq(-6,-3,length.out=100), lls.2)
        }
        
        x.t = x.tp1
        y.t = y.tp1
        log.sigma.x.t = log.sigma.x.tp1
        log.sigma.y.t = log.sigma.y.tp1
        rho.tilde.t = rho.tilde.tp1
    }

    out = NULL;
    out$timeseries.dt <- output
    out$parameters.hat <- list(T = T,
                               Delta = Delta,
                               mu.hat = mu.hat,
                               theta.hat = theta.hat,
                               alpha.hat = alpha.hat,
                               tau.hat = tau.hat)
    out$parameters <- list(mu.x = mu.x,
                           mu.y = mu.y,
                           alpha.x = alpha.x,
                           alpha.y = alpha.y,
                           theta.x = theta.x,
                           theta.y = theta.y,
                           tau.x = tau.x,
                           tau.y = tau.y,
                           tau.rho = tau.rho,
                           leverage.x.rho = leverage.x.rho,
                           leverage.y.rho = leverage.y.rho)
    return (out)
}


data <- generate.data(T,Delta,
                      mu.x,mu.y,
                      alpha.x,alpha.y,
                      theta.x,theta.y,
                      tau.x, tau.y,tau.rho,
                      leverage.x.rho, leverage.y.rho)

    
sample.data <- data$timeseries.dt

save(file = "mv-stochastic-vol-data.Rdata", list = c("sample.data",
                                                     "macro.parameters"))

par(mfrow=c(3,2),
    mar=c(2,2,1,1))
plot(sample.data[,x], type="l")
plot(sample.data[,y], type="l")
plot(sample.data[,log.sigma.x], type="l")
plot(sample.data[,log.sigma.y], type="l")
plot(sample.data[, logit.inv(rho.tilde)*2-1], type="l")
