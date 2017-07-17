library("data.table")

logit <- function(p) {
    out = log(p/(1-p))
    return (out)
}

logit.inv <- function(logit.p) {
    return (exp(logit.p)/(exp(logit.p) + 1))
}

T = 1 * 256 * 6.5 * 3600 * 1000 ## one year ms
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

generate.data <- function(T, Delta, mu.x, mu.y,
                          alpha.x, alpha.y, theta.x, theta.y,
                          tau.x, tau.y, tau.rho,
                          leverage.x.rho, leverage.y.rho) {
    N <- ceiling(T/Delta)

    rho.current <- -0.3
    rho.tilde.current <- logit((rho.current+1)/2)

    x.current <- log(100)
    y.current <- log(100)

    log.sigma.x.current <- alpha.x
    log.sigma.y.current <- alpha.y

    innovations <- data.table(matrix(data=rnorm(n=N*5),
                                     nrow=N,
                                     ncol=5))
    setnames(innovations, seq(1,5), c("epsilon.x", "epsilon.y",
                                      "eta.x", "eta.y", "eta.rho"))

    innovations[, eta.x := leverage.x.rho*epsilon.x +
                      sqrt(1-leverage.x.rho^2)*eta.x]
    innovations[, eta.y := leverage.y.rho*epsilon.y +
                      sqrt(1-leverage.y.rho^2)*eta.y]
        
    output <- data.table(x = rep(0,N),
                         y = rep(0,N),
                         log.sigma.x = rep(0,N),
                         log.sigma.y = rep(0,N),
                         rho.tilde = rep(0,N))

    for (i in seq(1,N)) {
        dx <- sqrt(1-rho.current^2)*exp(log.sigma.x.current)*innovations[i,1]+
            rho.current*exp(log.sigma.x.current)*innovations[i,2]
        dy <- exp(log.sigma.y.current)*innovations[i,2]
        
        x.current <- x.current + mu.x + dx
        y.current <- y.current + mu.x + dy

        log.sigma.x.current = alpha.x +
            theta.x*(log.sigma.x.current-alpha.x) + tau.x*innovations[i,3]

        log.sigma.y.current = alpha.y +
            theta.y*(log.sigma.y.current-alpha.y) + tau.y*innovations[i,4]

        rho.tilde.current = rho.tilde.current + tau.rho*innovations[i,5]
        rho.current = 2*logit.inv(rho.tilde.current) - 1

        output[i,names(output) := as.list(c(x.current,
                                            y.current,
                                            log.sigma.x.current,
                                            log.sigma.y.current,
                                            rho.tilde.current))]
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

par(mfrow=c(3,2),
    mar=c(1,2,1,1))
plot(sample.data[,x], type="l")
plot(sample.data[,y], type="l")
plot(sample.data[,log.sigma.x], type="l")
plot(sample.data[,log.sigma.y], type="l")
plot(sample.data[, logit.inv(rho.tilde)*2-1], type="l")
