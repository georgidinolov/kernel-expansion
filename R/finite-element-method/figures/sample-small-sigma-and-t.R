 rm(list=ls());
source("../2-d-solution.R")
library("parallel")
library("latex2exp");
dt <- 1.0/1e5;
rho = 0.0;
N <- 500;

sample.small.sigma.and.t <- function(rho=0, dt=1e-4, N=500) {
    problem.parameters.generate.data = NULL;
    log.sigma.x.samples <- rnorm(n=N, mean=1, sd=1);
    log.sigma.y.samples <- rnorm(n=N, mean=1, sd=1);

    log.sigma.x.samples.misspec <- rnorm(n=N, mean=1, sd=1);
    log.sigma.y.samples.misspec <- rnorm(n=N, mean=1, sd=1);

    samples.list <- lapply(X=seq(1,N),
                           FUN=function(x,
                                        log.sigma.x.samples,
                                        log.sigma.y.samples,
                                        log.sigma.x.samples.misspec,
                                        log.sigma.y.samples.misspec) {
                               return( list(log.sigma.x = log.sigma.x.samples[x],
                                            log.sigma.y = log.sigma.y.samples[x],
                                            log.sigma.x.misspec = log.sigma.x.samples.misspec[x],
                                            log.sigma.y.misspec = log.sigma.y.samples.misspec[x]) )
                           },
                           log.sigma.x.samples = log.sigma.x.samples,
                           log.sigma.y.samples = log.sigma.y.samples,
                           log.sigma.x.samples.misspec = log.sigma.x.samples.misspec,
                           log.sigma.y.samples.misspec = log.sigma.y.samples.misspec)
    cl <- makeCluster(4);
    small.sigma.and.t.samples <- parLapply(cl=cl,
                                           X=samples.list,
                                           fun=function(x,rho,dt) {
                                               source("../2-d-solution.R")
                                               log.sigma.x <- x$log.sigma.x;
                                               log.sigma.y <- x$log.sigma.y;
                                               log.sigma.x.miss <- x$log.sigma.x.misspec;
                                               log.sigma.y.miss <- x$log.sigma.y.misspec;

                                               problem.parameters.generate.data = NULL;
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
                                                   t.tilde.tilde <- 1 * (sigma.x.tilde)^2
                                               }

                                               sigma.x.tilde.tilde <- sigma.x.tilde / sigma.x.tilde;
                                               sigma.y.tilde.tilde <- sigma.y.tilde / sigma.x.tilde;

                                               out = NULL;
                                               out$sigma.y.tilde.tilde = sigma.y.tilde.tilde;
                                               out$t.tilde.tilde = t.tilde.tilde;
                                               out$datum = datum;
                                               out$log.sigma.x.miss = log.sigma.x.miss;
                                               out$log.sigma.y.miss = log.sigma.y.miss;
                                               return(out);
                                           },
                                           rho=rho,
                                           dt=dt)
    log.sigma.y.tilde.tildes <- rep(NA, N);
    log.t.tilde.tildes <- rep(NA,N);
    data <- list(N);

    for (x in seq(1,N)) {
        log.sigma.y.tilde.tildes[x] = log(small.sigma.and.t.samples[[x]]$sigma.y.tilde.tilde);
        log.t.tilde.tildes[x] = log(small.sigma.and.t.samples[[x]]$t.tilde.tilde);
        data[[x]] = small.sigma.and.t.samples[[x]]$datum;
    }

    out = NULL;
    out$log.sigma.y.tilde.tildes = log.sigma.y.tilde.tildes;
    out$log.t.tilde.tildes = log.t.tilde.tildes;
    out$data=data;
    return(out);
}

rhos <- c(-0.9, 0.0, 0.9);
sample.results <- vector(mode="list", length=length(rhos));
for (rho in rhos) {
    results <- sample.small.sigma.and.t(rho = rho, dt=1e-4, N=5000);
    min.t.index <- which(results$log.t.tilde.tildes ==
                         min(results$log.t.tilde.tildes));
    min.log.sigma.index <- which(results$log.sigma.y.tilde.tildes ==
                         min(results$log.sigma.y.tilde.tildes));
    
    pdf(paste("small-sigma-t-scatterplot-", rho, ".pdf", sep=""), 4, 4);
    par(mar=c(4,4.5,1,1));
    plot(results$log.sigma.y.tilde.tildes,
         results$log.t.tilde.tildes,
         xlim=c(-1.8,0),
         ylim=c(-3.5, 1.0),
         ylab = TeX("$\\log(\\tilde{t})$"),
         xlab = TeX("$\\log(\\tau_y / \\tau_x)$"),
         main = TeX(paste("$\\rho$ = ", rho, sep="")));
    points(results$log.sigma.y.tilde.tildes[min.t.index],
           results$log.t.tilde.tildes[min.t.index], col="red", pch=20);
    points(results$log.sigma.y.tilde.tildes[min.log.sigma.index],
           results$log.t.tilde.tildes[min.log.sigma.index], col="green", pch=20);
    dev.off();

    sample.results[[which(rhos==rho)]]$rho = rho;
    sample.results[[which(rhos==rho)]]$samples = results;
    sample.results[[which(rhos==rho)]]$min.t.index = min.t.index
    sample.results[[which(rhos==rho)]]$min.log.sigma.index = min.log.sigma.index;
}
save(file = "sample-small-sigma-and-t.Rdata", list=c("sample.results"));
