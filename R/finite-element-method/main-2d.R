rm(list=ls());

PLOT.SOLUTION=FALSE;
dx = 0.01;
dy = 0.01;
K=8;

library("mvtnorm");
source("2-d-solution.R");
problem.parameters = NULL;
problem.parameters$ax = -1;
problem.parameters$bx = 1;
problem.parameters$ay = -1;
problem.parameters$by = 1;
problem.parameters$x.ic = 0.9;
problem.parameters$y.ic = 0.9;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2.x = 1e-1;
problem.parameters$sigma.2.y = 1e-1;
problem.parameters$rho = 0.0;
problem.parameters$t = 0.5;

Lx <- (problem.parameters$bx-problem.parameters$ax)/K;
Kx.prime <- floor((problem.parameters$x.ic-problem.parameters$ax)/Lx);
Lx.prime <- (problem.parameters$x.ic-problem.parameters$ax)/Kx.prime;

Ly <- (problem.parameters$by-problem.parameters$ay)/K;
Ky.prime <- floor((problem.parameters$y.ic-problem.parameters$ay)/Ly);
Ly.prime <- (problem.parameters$y.ic-problem.parameters$ay)/Ky.prime;

mu.xs = c(seq(problem.parameters$ax, problem.parameters$bx,
              by = Lx.prime));
mu.xs.2 = unlist(lapply(seq(1,K), function(k){rep(mu.xs[k], K)}));

mu.ys = c(seq(problem.parameters$ax, problem.parameters$bx,
              by=Ly.prime));
mu.ys.2 = rep(mu.ys, K);

log.sigma2.xs = log(rep(((problem.parameters$bx-
                          problem.parameters$ax)/(K))^2,
                        K));
log.sigma2.xs.2 = unlist(lapply(seq(1,K), function(k){
    rep(log.sigma2.xs[k], K)}));

log.sigma2.ys = log(rep(((problem.parameters$by-
                          problem.parameters$ay)/(K))^2,
                        K));
log.sigma2.ys.2 = unlist(lapply(seq(1,K), function(k){
    rep(log.sigma2.ys[k], K)}));

mus = rbind(mu.xs.2, mu.ys.2);
log.sigma2s = rbind(log.sigma2.xs.2, log.sigma2.ys.2);
## ## rotating according to the geometry of the problem
## theta = atan(-problem.parameters$rho);
## theta = 0;
## Rotation.matrix = matrix(nrow=2, ncol=2,
##                          byrow=FALSE,
##                          data=c(c(cos(theta),sin(theta)),
##                                 c(-sin(theta),cos(theta))));

## mus = Rotation.matrix %*% (mus -
##                            rbind(rep(problem.parameters$x.ic, length(mus[1,])),
##                                  rep(problem.parameters$y.ic, length(mus[1,])))) +
##     rbind(rep(problem.parameters$x.ic, length(mus[1,])),
##           rep(problem.parameters$y.ic, length(mus[1,])));

all.inside <- seq(1,length(mus[1,]));
all.inside <- mus[1,] <= problem.parameters$bx &
    mus[1,] >= problem.parameters$ax &
    mus[2,] <= problem.parameters$by &
    mus[2,] >= problem.parameters$ay;
plot(mus[1,all.inside], mus[2,all.inside],col="red");
mus <- mus[,all.inside];
log.sigma2s <- log.sigma2s[,all.inside];

## x.ax.distances <- abs(mus[1,]-problem.parameters$ax);
## x.bx.distances <- abs(mus[1,]-problem.parameters$bx);
## y.ay.distances <- abs(mus[2,]-problem.parameters$ay);
## y.by.distances <- abs(mus[2,]-problem.parameters$by);
## min.x.distances <- apply(rbind(x.ax.distances,
##                                x.bx.distances),2,min);
## min.y.distances <- apply(rbind(y.ay.distances,
##                                y.by.distances),2,min);
## log.sigma2s[1,] <- ifelse(min.x.distances > 1e-15 &
##                           min.x.distances < exp(log.sigma2s[1,]/2),
##                           log((min.x.distances/2)^2),
##                           log.sigma2s[1,]);
## log.sigma2s[2,] <- ifelse(min.y.distances > 1e-15 &
##                           min.y.distances < exp(log.sigma2s[2,]/2),
##                           log((min.y.distances/2)^2),
##                           log.sigma2s[2,]);

distances <- apply(mus -
                   c(problem.parameters$x.ic,
                     problem.parameters$y.ic),
                   2,
                   function(x){sum(x^2)});

sorted.mus <- sort.int(distances,index.return=TRUE);

blackbox.wrapper <- function(log.C) {
    C = exp(log.C);
    l2 = blackbox(mus[,sorted.mus$ix],
                  log((sqrt(exp(log.sigma2s[,sorted.mus$ix])) *
                       C)^2),
                  problem.parameters,
                  dx,dy,FALSE,TRUE);
    print(c(C,l2));
    return(l2);
}

opt.results <- optimize(f=blackbox.wrapper,
                        lower = log(0.1),
                        upper = log(K),
                        tol=1e-2);

C = exp(opt.results$minimum);

problem.parameters$rho = -0.8;
problem.parameters$x.ic = 0;
l2 = blackbox(mus[,sorted.mus$ix],
              log((sqrt(exp(log.sigma2s[,sorted.mus$ix])) *
                   C)^2),
              problem.parameters,
              dx,dy,TRUE,TRUE);
