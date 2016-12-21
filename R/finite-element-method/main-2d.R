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
problem.parameters$rho = -0.9;
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
all.inside <- mus[1,] < problem.parameters$bx &
    mus[1,] > problem.parameters$ax &
    mus[2,] < problem.parameters$by &
    mus[2,] > problem.parameters$ay;
plot(mus[1,all.inside], mus[2,all.inside],col="red");
mus <- mus[,all.inside];
log.sigma2s <- log.sigma2s[,all.inside];

distances <- apply(mus -
                   c(problem.parameters$x.ic,
                     problem.parameters$y.ic),
                   2,
                   function(x){sum(x^2)});

sorted.mus <- sort.int(distances,index.return=TRUE);

bb = blackbox(mus[,sorted.mus$ix],
              log.sigma2s,
              problem.parameters,
              dx, dy,
              TRUE,TRUE);

## mu.log.sigma2.x.pairs.list <-
##     lapply(X=seq(1,length(mus[1,])),
##            function(x) {
##                c(mus[1,x], log.sigma2s[1,x])
##            });

## mu.log.sigma2.y.pairs.list <-
##     lapply(X=seq(1,length(mus[1,])),
##            function(x) {
##                c(mus[2,x], log.sigma2s[2,x])
##            });

## ## hash function ##
## ## unique means and variances ##
## ## l^2 norms
## l2.norms.x <- matrix(nrow = length(mu.log.sigma2.x.pairs.list),
##                      ncol = length(mu.log.sigma2.x.pairs.list));
## l2.norms.y <- matrix(nrow = length(mu.log.sigma2.y.pairs.list),
##                      ncol = length(mu.log.sigma2.y.pairs.list));

## ## TODO(georgid): This can be vectorized.
## x = seq(problem.parameters$ax,
##         problem.parameters$bx,
##         by = dx);
## y = seq(problem.parameters$ay,
##         problem.parameters$by,
##         by = dy);

## for (k in seq(1,length(mu.log.sigma2.x.pairs.list))) {
    
##     k.function = (x-problem.parameters$ax)*
##         (problem.parameters$bx-x)*
##         dnorm(x,
##               mean = mu.log.sigma2.x.pairs.list[[k]][1],
##               sd=sqrt(exp(mu.log.sigma2.x.pairs.list[[k]][2])));
    
##     for (l in seq(1,length(mu.log.sigma2.x.pairs.list))) {

##         l.function = (x-problem.parameters$ax)*
##             (problem.parameters$bx-x)*
##             dnorm(x,
##                   mean = mu.log.sigma2.x.pairs.list[[l]][1],
##                   sd=sqrt(exp(mu.log.sigma2.x.pairs.list[[l]][2])));
        
##         l2.norms.x[k,l] <- sqrt(sum((k.function-
##                                      l.function)^2)*dx);
##         ## if (l2.norms.x[k,l] < 0.1) {
##         ##     plot(x,k.function,type="l");
##         ##     lines(x,l.function,col="red");
##         ## }
##     }
## }

## ## TODO(georgid): This can be vectorized.
## for (k in seq(1,length(mu.log.sigma2.y.pairs.list))) {

##     k.function = (y-problem.parameters$ay)*
##         (problem.parameters$by-y)*
##         dnorm(y,
##               mean = mu.log.sigma2.y.pairs.list[[k]][1],
##               sd=sqrt(exp(mu.log.sigma2.y.pairs.list[[k]][2])));
    
##     for (l in seq(1,length(mu.log.sigma2.y.pairs.list))) {

##         l.function = (y-problem.parameters$ay)*
##             (problem.parameters$by-y)*
##             dnorm(y,
##                   mean = mu.log.sigma2.y.pairs.list[[l]][1],
##                   sd=sqrt(exp(mu.log.sigma2.y.pairs.list[[l]][2])));

        
##         l2.norms.y[k,l] <- sqrt(sum((mu.log.sigma2.y.pairs.list[[k]] -
##                                      mu.log.sigma2.y.pairs.list[[l]])^2))

##     }
## }

## threshold = 1e-1;
## x.included <- c();
## for (k in seq(1,length(mu.log.sigma2.x.pairs.list))) {
##     if (k==1) {
##         x.included = c(x.included, k);
##     } else {
##         if (sum(l2.norms.x[k,seq(1,k-1)] <= threshold) == 0) {
##             x.included = c(x.included, k);
##         }
##     }
## }

## y.included <- c();
## for (k in seq(1,length(mu.log.sigma2.y.pairs.list))) {
##     if (k==1) {
##         y.included = c(y.included, k);
##     } else {
##         if (sum(l2.norms.y[k,seq(1,k-1)] <= threshold) == 0) {
##             y.included = c(y.included, k);
##         }
##     }
## }

## mu.log.sigma2.x.pairs.list.unique <- mu.log.sigma2.x.pairs.list[x.included];
## mu.log.sigma2.y.pairs.list.unique <- mu.log.sigma2.y.pairs.list[y.included];

## ## simple hash funtion ##
## ## k in seq(1,length(mus[1,])) ##
## K = length(mus[1,]);
## simple.hash <- function(k) {
##     mu.log.sigma2.x <- mu.log.sigma2.x.pairs.list[[k]];
##     x.distances <-
##         unlist(lapply(seq(1,length(mu.log.sigma2.x.pairs.list.unique)),
##                       function(x) {
##                           sqrt(sum((mu.log.sigma2.x -
##                                     mu.log.sigma2.x.pairs.list.unique[[x]])^2))
##                       }));
##     k.x <-
##         which(x.distances == min(x.distances))[1];

##     mu.log.sigma2.y <- mu.log.sigma2.y.pairs.list[[k]];
##     y.distances <-
##         unlist(lapply(seq(1,length(mu.log.sigma2.y.pairs.list.unique)),
##                       function(x) {
##                           sqrt(sum((mu.log.sigma2.y -
##                                     mu.log.sigma2.y.pairs.list.unique[[x]])^2))
##                       }));
##     k.y <-which(y.distances == min(y.distances))[1];
##     return (c(k.x,k.y));
## }

## k.xy.hash <- sapply(seq(1,K),
##                     function(x) {
##                         simple.hash(x) });
               

