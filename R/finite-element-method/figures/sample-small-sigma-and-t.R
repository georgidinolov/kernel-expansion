rm(list=ls());
source("../2-d-solution.R")
dt <- 1.0/1000;

N <- 10000;
problem.parameters.generate.data = NULL;

log.sigma.x.samples <- rnorm(n=N, mean=1, sd=1);
log.sigma.y.samples <- rnorm(n=N, mean=1, sd=1);

log.sigma.y.tilde.tildes <- rep(NA,N);
t.tilde.tildes <- rep(NA,N);

for (n in seq(1,N)) {
    log.sigma.x <- log.sigma.x.samples[1];
    log.sigma.y <- log.sigma.x.samples[2];
    rho = 0.0;
    
    problem.parameters.generate.data$t <- 1.0;
    problem.parameters.generate.data$sigma.2.x <- exp(2*log.sigma.x);
    problem.parameters.generate.data$sigma.2.y <- exp(2*log.sigma.y);
    problem.parameters.generate.data$rho <- rho;
    problem.parameters.generate.data$x.ic <- 0;
    problem.parameters.generate.data$y.ic <- 0;

    datum <- sample.process(n.samples=1,
                           dt=dt,
                           problem.parameters.generate.data=problem.parameters.generate.data)[[1]]

    Lx <- datum$bx - datum$ax;
    Ly <- datum$by - datum$ay;

    sigma.x.tilde <- exp(log.sigma.x - log(Lx));
    sigma.y.tilde <- exp(log.sigma.y - log(Ly));
    t.tilde.tilde <- 1 * (sigma.x.tilde)^2

    if (sigma.y.tilde > sigma.x.tilde) {
            sigma.x.tilde <- exp(log.sigma.y - log(Ly));
            sigma.y.tilde <- exp(log.sigma.x - log(Lx));
            t.tilde.tilde <- 1 * (sigma.y.tilde)^2
    }

    sigma.x.tilde.tilde <- sigma.x.tilde / sigma.x.tilde;
    sigma.y.tilde.tilde <- sigma.y.tilde / sigma.x.tilde;

    log.sigma.y.tilde.tildes[n] <- log(sigma.y.tilde.tilde);
    t.tilde.tildes[n] <- t.tilde.tilde;
}

plot(exp(log.sigma.y.tilde.tildes), t.tilde.tildes)
