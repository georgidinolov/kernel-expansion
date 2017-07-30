logit.inv <- function(logit.p) {
    return (exp(logit.p)/(exp(logit.p) + 1))
}

## theta is a vector (log(sigma_x), log(sigma_y), logit((rho+1)/2))
theta.to.nominal <- function(theta) {
    return( c(exp(theta[1]),
              exp(theta[2]),
              logit.inv(theta[3])*2-1))

}

epsilons.given.theta <- function(theta.current,
                                 y.current,
                                 y.current.m1,
                                 parameters) {

    nominal.current = theta.to.nominal(theta.current)

    epsilon.y = (y.current[2] - y.current.m1[2] - parameters$mu.y)/
        nominal.current[2]
    
    epsilon.x = (y.current[1] -
                 y.current.m1[1] -
                 parameters$mu.x -
                 nominal.current[3]*nominal.current[1]*epsilon.y) /
        (sqrt(1 - nominal.current[3]^2) * nominal.current[1])

    return (c(epsilon.x, epsilon.y))
}

theta.next.mean <- function(theta.current,
                            y.current,
                            y.current.m1,
                            parameters) {

    epsilons = epsilons.given.theta(theta.current,
                                    y.current,
                                    y.current.m1,
                                    parameters)
        
    eta.means = c(parameters$leverage.x.rho*epsilons[1],
                  parameters$leverage.y.rho*epsilons[2],
                  0)
    
    return( c(parameters$alpha.x +
              parameters$theta.x*(theta.current[1]-parameters$alpha.x),
              parameters$alpha.y +
              parameters$theta.y*(theta.current[2]-parameters$alpha.y),
              theta.current[3]) +
            c(parameters$tau.x,
              parameters$tau.y,
              parameters$tau.rho)*eta.means )
}

theta.next.mean.par <- function(theta.current,
                                y.current,
                                y.current.m1,
                                parameters) {
    out = cbind(sapply(seq(1,dim(theta.current)[2]),
                       function(x) {theta.next.mean(theta.current[,x],
                                                    y.current,
                                                    y.current.m1,
                                                    parameters)}))
    return (out)
}

sample.theta <- function(theta.current,
                         y.current,
                         y.current.m1,
                         parameters) {
    
    epsilons = epsilons.given.theta(theta.current,
                                    y.current,
                                    y.current.m1,
                                    parameters)
    
    eta.means = c(parameters$leverage.x.rho*epsilons[1],
                  parameters$leverage.y.rho*epsilons[2],
                  0)

    eta.vars = c(1-parameters$leverage.x.rho^2,
                 1-parameters$leverage.y.rho^2,
                 1)
    
    etas = eta.means + sqrt(eta.vars)*rnorm(3)
    
    out = theta.next.mean(theta.current,
                          y.current,
                          y.current.m1,
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

## p(y.t | y.t.minus.1, theta.t, theta, parameters)
log.likelihood <- function(y.t,
                           y.t.minus.1,
                           theta.t,
                           parameters) {

    nominal.t = theta.to.nominal(theta.t)

    ## Y
    epsilon.2 = (y.t[2] - y.t.minus.1[2] - parameters$mu.y)/nominal.t[2]

    return( dnorm(epsilon.2,
                  0,
                  1,
                  log=TRUE)  +
            dnorm((y.t[1] -
                   y.t.minus.1[1] -
                   parameters$mu.x -
                   epsilon.2*nominal.t[3]*nominal.t[1])/
                  (sqrt(1-nominal.t[3]^2)*nominal.t[1]),
                  0,
                  1,
                  log=TRUE) )
}

log.likelihood.par <- function(y.t,
                               y.t.minus.1,
                               theta.t,
                               parameters,
                               N.particles) {

    out = sapply(seq(1,N.particles),
                 function(x) {
                     log.likelihood(y.t,
                                    y.t.minus.1,
                                    theta.t[,x],
                                    parameters)
                 })
    return (out)
}
