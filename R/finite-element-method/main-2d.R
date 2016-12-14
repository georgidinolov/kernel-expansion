rm(list=ls());

PLOT.SOLUTION=TRUE;
dx = 0.001;
dy = 0.001;
K=10;

source("2-d-solution.R");
problem.parameters = NULL;
problem.parameters$ax = -1;
problem.parameters$bx = 1;
problem.parameters$ay = -1;
problem.parameters$by = 1;
problem.parameters$x.ic = 0.0;
problem.parameters$y.ic = 0.0;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2.x = 5e-1;
problem.parameters$sigma.2.y = 1e0;
problem.parameters$rho = 0.5;
problem.parameters$t = 0.2;

log.sigma2.vector=log(rep(1,K));
## if (K==1) {
##     mu.vector = c(problem.parameters$x.ic,
##                   problem.parameters$y.ic)
    
## } else {
##     mu.vector = c(c(problem.parameters$x.ic,
##                     problem.parameters$y.ic),
##                   unlist(lapply(seq(1,K-1),
##                                 function(x)
##                                 { c(runif(1, problem.parameters$ax,
##                                           problem.parameters$bx),
##                                     runif(1, problem.parameters$ay,
##                                           problem.parameters$by)) })));
## }

mu.xs = c(seq(problem.parameters$ax, problem.parameters$bx,
              length.out = K));
mu.xs.2 = unlist(lapply(seq(1,K), function(k){rep(mu.xs[k], K)}));

mu.ys = c(seq(problem.parameters$ax, problem.parameters$bx,
              length.out = K));
mu.ys.2 = rep(mu.ys, K);

log.sigma2.xs = log(rep(((problem.parameters$bx-
                          problem.parameters$ax)/K)^2,
                        K));
log.sigma2.xs.2 = unlist(lapply(seq(1,K), function(k){
    rep(log.sigma2.xs[k], K)}));

log.sigma2.ys = log(rep(((problem.parameters$by-
                          problem.parameters$ay)/K)^2,
                        K));
log.sigma2.ys.2 = unlist(lapply(seq(1,K), function(k){
    rep(log.sigma2.ys[k], K)}));

mus = rbind(mu.xs.2, mu.ys.2);
log.sigma2s = rbind(log.sigma2.xs.2, log.sigma2.ys.2);
## rotating according to the geometry of the problem
theta = atan(-problem.parameters$rho);
theta = 0;
Rotation.matrix = matrix(nrow=2, ncol=2,
                         byrow=FALSE,
                         data=c(c(cos(theta),sin(theta)),
                                c(-sin(theta),cos(theta))));

mus = Rotation.matrix %*% (mus -
                           rbind(rep(problem.parameters$x.ic, length(mus[1,])),
                                 rep(problem.parameters$y.ic, length(mus[1,])))) +
    rbind(rep(problem.parameters$x.ic, length(mus[1,])),
          rep(problem.parameters$y.ic, length(mus[1,])))


all.inside <- seq(1,length(mus[1,]));
all.inside <- mus[1,] <= problem.parameters$bx &
    mus[1,] >= problem.parameters$ax &
    mus[2,] <= problem.parameters$by &
    mus[2,] >= problem.parameters$ay;
plot(mus[1,all.inside], mus[2,all.inside],col="red");
mus <- mus[,all.inside];
log.sigma2s <- log.sigma2s[,all.inside];

mu.log.sigma2.x.pairs.list <-
    lapply(X=seq(1,length(mus[1,])),
           function(x) {
               c(mus[1,x], log.sigma2s[1,x])
           });

mu.log.sigma2.y.pairs.list <-
    lapply(X=seq(1,length(mus[1,])),
           function(x) {
               c(mus[2,x], log.sigma2s[2,x])
           });

## hash function ##
## unique means and variances ##
## l^2 norms
l2.norms.x <- matrix(nrow = length(mu.log.sigma2.x.pairs.list),
                     ncol = length(mu.log.sigma2.x.pairs.list));
l2.norms.y <- matrix(nrow = length(mu.log.sigma2.y.pairs.list),
                     ncol = length(mu.log.sigma2.y.pairs.list));

## TODO(georgid): This can be vectorized.
x = seq(problem.parameters$ax,
        problem.parameters$bx,
        by = dx);
y = seq(problem.parameters$ay,
        problem.parameters$by,
        by = dy);

for (k in seq(1,length(mu.log.sigma2.x.pairs.list))) {
    
    k.function = (x-problem.parameters$ax)*
        (problem.parameters$bx-x)*
        dnorm(x,
              mean = mu.log.sigma2.x.pairs.list[[k]][1],
              sd=sqrt(exp(mu.log.sigma2.x.pairs.list[[k]][2])));
    
    for (l in seq(1,length(mu.log.sigma2.x.pairs.list))) {

        l.function = (x-problem.parameters$ax)*
            (problem.parameters$bx-x)*
            dnorm(x,
                  mean = mu.log.sigma2.x.pairs.list[[l]][1],
                  sd=sqrt(exp(mu.log.sigma2.x.pairs.list[[l]][2])));
        
        l2.norms.x[k,l] <- sqrt(sum((k.function-
                                     l.function)^2)*dx);
        ## if (l2.norms.x[k,l] < 0.1) {
        ##     plot(x,k.function,type="l");
        ##     lines(x,l.function,col="red");
        ## }
    }
}

## TODO(georgid): This can be vectorized.
for (k in seq(1,length(mu.log.sigma2.y.pairs.list))) {

    k.function = (y-problem.parameters$ay)*
        (problem.parameters$by-y)*
        dnorm(y,
              mean = mu.log.sigma2.y.pairs.list[[k]][1],
              sd=sqrt(exp(mu.log.sigma2.y.pairs.list[[k]][2])));
    
    for (l in seq(1,length(mu.log.sigma2.y.pairs.list))) {

        l.function = (y-problem.parameters$ay)*
            (problem.parameters$by-y)*
            dnorm(y,
                  mean = mu.log.sigma2.y.pairs.list[[l]][1],
                  sd=sqrt(exp(mu.log.sigma2.y.pairs.list[[l]][2])));

        
        l2.norms.y[k,l] <- sqrt(sum((mu.log.sigma2.y.pairs.list[[k]] -
                                     mu.log.sigma2.y.pairs.list[[l]])^2))

    }
}

threshold = 1e-1;
x.included <- c();
for (k in seq(1,length(mu.log.sigma2.x.pairs.list))) {
    if (k==1) {
        x.included = c(x.included, k);
    } else {
        if (sum(l2.norms.x[k,seq(1,k-1)] <= threshold) == 0) {
            x.included = c(x.included, k);
        }
    }
}

y.included <- c();
for (k in seq(1,length(mu.log.sigma2.y.pairs.list))) {
    if (k==1) {
        y.included = c(y.included, k);
    } else {
        if (sum(l2.norms.y[k,seq(1,k-1)] <= threshold) == 0) {
            y.included = c(y.included, k);
        }
    }
}

mu.log.sigma2.x.pairs.list.unique <- mu.log.sigma2.x.pairs.list[x.included];
mu.log.sigma2.y.pairs.list.unique <- mu.log.sigma2.y.pairs.list[y.included];

## simple hash funtion ##
## k in seq(1,length(mus[1,])) ##
K = length(mus[1,]);
simple.hash <- function(k) {
    mu.log.sigma2.x <- mu.log.sigma2.x.pairs.list[[k]];
    x.distances <-
        unlist(lapply(seq(1,length(mu.log.sigma2.x.pairs.list.unique)),
                      function(x) {
                          sqrt(sum((mu.log.sigma2.x -
                                    mu.log.sigma2.x.pairs.list.unique[[x]])^2))
                      }));
    k.x <-
        which(x.distances == min(x.distances))[1];

    mu.log.sigma2.y <- mu.log.sigma2.y.pairs.list[[k]];
    y.distances <-
        unlist(lapply(seq(1,length(mu.log.sigma2.y.pairs.list.unique)),
                      function(x) {
                          sqrt(sum((mu.log.sigma2.y -
                                    mu.log.sigma2.y.pairs.list.unique[[x]])^2))
                      }));
    k.y <-which(y.distances == min(y.distances))[1];
    return (c(k.x,k.y));
}

k.xy.hash <- sapply(seq(1,K),
                    function(x) {
                        simple.hash(x) });
               
bb = blackbox(mu.log.sigma2.x.pairs.list.unique,
              mu.log.sigma2.y.pairs.list.unique, 
              k.xy.hash,
              problem.parameters,
              dx, dy,
              TRUE,TRUE);

## points(mu.vector[seq(2,2*K,by=2)], mu.vector[seq(1,2*K,by=2)], col="red")

## rm(list=ls());
## source("2-d-solution.R");

## problem.parameters = NULL;
## problem.parameters$ax = -1;
## problem.parameters$bx = 1;
## problem.parameters$ay = -1;
## problem.parameters$by = 1;
## problem.parameters$x.ic = 0.2;
## problem.parameters$y.ic = 0.0;
## problem.parameters$number.terms = 1000;
## problem.parameters$sigma.2.x = 0.5;
## problem.parameters$sigma.2.y = 0.0001;
## problem.parameters$rho = 0.0;
## problem.parameters$t = 0.001;

## Ks = seq(5,5);

## L2.remainders.before = rep(NA, length(Ks));
## L2.diff.before = rep(NA, length(Ks));

## L2.remainders.after = rep(NA, length(Ks));
## L2.diff.after = rep(NA, length(Ks));

## ave.function.call.time.vec = rep(NA,length(Ks));
## log.sigma2.mu.vector.list = vector(mode="list",
##                                    length=length(Ks));

## for (i in seq(1,length(Ks))) {
##     K=Ks[i];
##     log.sigma2.vector=log(
##         rep(c(((problem.parameters$bx-problem.parameters$ax)/K)^2,
##         ((problem.parameters$by-problem.parameters$ay)/K)^2),
##         K)) +
##         rep(log(2), 2*K);
##     mu.vector = seq(problem.parameters$ax,
##                     problem.parameters$bx,
##                     length.out=length(log.sigma2.vector)/2);

##     mu.vector = unlist(lapply(X=mu.vector, function(x) {c(x,x)}));
##     mu.vector = c(-0.8,-0.8, -0.8,0.8, 0.8,-0.8, 0.8,0.8, 0,0);
    
##     log.sigma2.mu.vector = c(log.sigma2.vector,
##                              mu.vector);
##     dx = 0.001;
##     dy = 0.001;
    
##     bb = blackbox(log.sigma2.mu.vector, problem.parameters, dx, dy,
##                   TRUE,TRUE);
##     ## L2.remainders.before[i] = bb;
##     ## L2.diff.before[i] = blackbox(log.sigma2.mu.vector, problem.parameters,
##     ##                              dx,
##     ##                              FALSE,
##     ##                              FALSE);

##     ## ptm <- proc.time();				 
##     ## opt.bases <- optim(par=log.sigma2.mu.vector,
##     ##                    fn=blackbox,
##     ##                    problem.parameters = problem.parameters,
##     ##                    dx = dx,
##     ##                    PLOT.SOLUTION=FALSE,
##     ##                    MINIMIZE.REMAINDER=TRUE);

##     ## ## opt.bases <- optim(par=log.sigma2.mu.vector,
##     ## ## 			method=c("BFGS"),
##     ## ##                    fn=blackbox,
##     ## ##                    problem.parameters = problem.parameters,
##     ## ##                    dx = dx,
##     ## ##                    PLOT.SOLUTION=FALSE,
##     ## ##                    MINIMIZE.REMAINDER=TRUE);
##     ## end.time <- proc.time()-ptm;
##     ## ## Stop the clock
    
##     ## ave.function.call.time = end.time[3]/(opt.bases$counts[1]);
##     ## ave.function.call.time.vec[i] = ave.function.call.time;

##     ## log.sigma2.mu.vector = opt.bases$par;
##     ## log.sigma2.mu.vector.list[[i]] = log.sigma2.mu.vector;
##     ## bb = blackbox(log.sigma2.mu.vector, problem.parameters, dx,TRUE,TRUE);
    
##     ## L2.remainders.after[i] = bb;
##     ## L2.diff.after[i] = blackbox(log.sigma2.mu.vector, problem.parameters,
##     ##                             dx,
##     ##                             FALSE,
##     ##                             FALSE)
## }

## ## save(file="optimization-results.Rdata",
## ##      list=c("Ks","L2.remainders.before","L2.diff.before",
## ##             "L2.remainders.after","L2.diff.after",
## ##             "ave.function.call.time.vec",
## ##             "log.sigma2.mu.vector.list",
## ## 	    "dx"));

## ## pdf("optimization-results.pdf");
## ## par(mfrow=c(2,1));
## ## par(mar = c(5,4,2,1));
## ## plot(Ks, log(L2.remainders.before),
## ##      type="l",
## ##      ylim = c(min(log(L2.remainders.before),log(L2.remainders.after)),
## ##               max(log(L2.remainders.before),log(L2.remainders.after))),
## ##      ylab = "",
## ##      xlab = "K",
## ##      main = "log(L^2) norm of remainder term");
## ## lines(Ks, log(L2.remainders.after),col="red");

## ## plot(Ks, log(L2.diff.before),
## ##      type="l",
## ##      ylim = c(min(log(L2.diff.before),log(L2.diff.after)),
## ##               max(log(L2.diff.before),log(L2.diff.after))),
## ##      ylab = "",
## ##      xlab = "",
## ##      main = "log(L^2) norm of difference between approximate and true solutions");
## ## lines(Ks, log(L2.diff.after),col="red");
## ## dev.off();
