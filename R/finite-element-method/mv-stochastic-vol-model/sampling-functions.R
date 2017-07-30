logit.inv <- function(logit.p) {
    return (exp(logit.p)/(exp(logit.p) + 1))
}

## theta is a vector (log(sigma_x), log(sigma_y), logit((rho+1)/2))
theta.to.nominal <- function(theta) {
    return( c(exp(theta[1]),
              exp(theta[2]),
              logit.inv(theta[3])*2-1))

}

theta.next.mean <- function(theta.current,
                            parameters) {
    return( c(parameters$alpha.x +
              parameters$theta.x*(theta.current[1]-parameters$alpha.x),
              parameters$alpha.y +
              parameters$theta.y*(theta.current[2]-parameters$alpha.y),
              theta.current[3]) )
    return(theta.next.mean)
}

theta.next.mean.par <- function(theta.current,
                                parameters) {
    out = cbind(sapply(seq(1,dim(theta.current)[2]),
                       function(x) {theta.next.mean(theta.current[,x],
                                                    parameters)}))
    return (out)
}

sample.theta <- function(theta.current,
                         parameters) {
    etas = rnorm(3)
    out = theta.next.mean(theta.current,
                           parameters) +
        c(parameters$tau.x,
          parameters$tau.y,
          parameters$tau.rho) * etas
    return( out )
}

sample.theta.prior <- function(parameters) {
    ## initial rho is assumed 0
    log.sigma.x.current.samples =
        parameters$alpha.x +
        parameters$tau.x/sqrt(1-parameters$theta.x^2)*
        rnorm(1)

    log.sigma.y.current.samples =
        parameters$alpha.y +
        parameters$tau.y/sqrt(1-parameters$theta.y^2)*
        rnorm(1)
    
    rho.tilde.current.samples =
        parameters$tau.rho*rnorm(1)

    return( c(log.sigma.x.current.samples,
              log.sigma.y.current.samples,
              rho.tilde.current.samples) )
}

sample.theta.prior.par <- function(parameters,
                                   N.particles) {
    return( sapply(seq(1,N.particles),
                   function(x) {
                       return (sample.theta.prior(parameters))
                   }) )
}

## p(y.t | y.t.minus.1, theta.t, theta.t.plus.1, parameters)
log.likelihood <- function(y.t,
                           y.t.minus.1,
                           theta.t.plus.1,
                           theta.t,
                           parameters) {

    nominal.tp1 = theta.to.nominal(theta.t.plus.1)
    nominal.t = theta.to.nominal(theta.t)

    etas = (theta.t.plus.1 - c(parameters$alpha.x,
                               parameters$alpha.y,
                               0) -
            c(parameters$theta.x,
              parameters$theta.y,
              1)*(theta.t - c(parameters$alpha.x,
                              parameters$alpha.y,
                              0))) / c(parameters$tau.x,
                                       parameters$tau.y,
                                       parameters$tau.rho)

    epsilon.means = c(parameters$leverage.x.rho,
                       parameters$leverage.y.rho) * etas[1:2]
    epsilon.vars = 1 - c(parameters$leverage.x.rho,
                          parameters$leverage.y.rho)^2

    ## Y
    epsilon.2 = (y.t[2] - y.t.minus.1[2] - parameters$mu.y)/nominal.t[2]
    return( dnorm((y.t[2] - y.t.minus.1[2] - parameters$mu.y)/nominal.t[2],
                  epsilon.means[2],
                  sqrt(epsilon.vars[2]),
                  log=TRUE)  +
            dnorm((y.t[1] - y.t.minus.1[1] -
                   parameters$mu.x -
                   epsilon.2*nominal.t[3]*nominal.t[1])/
                  (sqrt(1-nominal.t[3]^2)*nominal.t[1]),
                  epsilon.means[1],
                  epsilon.vars[2],
                  log=TRUE) )
}

log.likelihood.par <- function(y.t,
                               y.t.minus.1,
                               theta.tp1,
                               theta.t,
                               parameters,
                               N.particles) {
    out = sapply(seq(1,N.particles),
                 function(x) {
                     log.likelihood(y.t,
                                    y.t.minus.1,
                                    theta.tp1[,x],
                                    theta.t[,x],
                                    parameters)
                 })
    return (out)
}
