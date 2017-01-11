rm(list=ls());
library("mvtnorm");
source("2-d-solution.R");
source("../classical-solution/2-d-solution.R");

PLOT.SOLUTION = TRUE;
dx = 0.001;
dy = 0.001;
K.prime = 11;

problem.parameters.generate.data = NULL;
problem.parameters.generate.data$t <- 1;
problem.parameters.generate.data$sigma.2.x <- 0.1;
problem.parameters.generate.data$sigma.2.y <- 1;
problem.parameters.generate.data$rho <- -0.99;
problem.parameters.generate.data$x.ic <- 0;
problem.parameters.generate.data$y.ic <- 0;
dt <- problem.parameters.generate.data$t/1000;
n.samples <- 100;

data <- sample.process(n.samples, dt, problem.parameters.generate.data);
problem.parameters <- data[[1]];
problem.parameters$K.prime <- K.prime;
problem.parameters$number.terms <- 100;

K <- (K.prime-1)^2;
function.list <- vector("list", K);

x <- seq(0,1,
         by=dx);
y <- seq(0,1,
         by=dy);

alpha.minus.1.xs <- rep(NA,K);
alpha.minus.1.ys <- rep(NA,K);
for (k in seq(1,K)) {
    alpha.minus.1.y <- ceiling(k/(K.prime-1));
    alpha.minus.1.x <- k - ((K.prime-1) * (alpha.minus.1.y-1));
    
    alpha.minus.1.ys[k] <- alpha.minus.1.y;
    alpha.minus.1.xs[k] <- alpha.minus.1.x;
}
alpha.xs <- alpha.minus.1.xs + 1;
alpha.ys <- alpha.minus.1.ys + 1;
## alpha + beta = K.prime + 2

x.basis.coors <- alpha.minus.1.xs / (K.prime);
y.basis.coors <- alpha.minus.1.ys / (K.prime);

sort.bases <- sort.int(sqrt((x.basis.coors-0.5)^2 +
                            (y.basis.coors-0.5)^2),
                       index.return = TRUE);
alpha.xs <- alpha.xs[sort.bases$ix];
alpha.ys <- alpha.ys[sort.bases$ix];

for (k in seq(1,K)) {
    alpha.x <- alpha.xs[k];
    alpha.y <- alpha.ys[k];
    print(c(alpha.x,alpha.y));
    function.params <- c(alpha.x, alpha.y);
    
    function.list[[k]] <-
        basis.function.xy(x,y,
                          function.params,
                          problem.parameters);
}

par(mfrow=c(ceiling(sqrt(K)),
            ceiling(sqrt(K))));
par(mar = c(5,4,2,1));
if (PLOT.SOLUTION) {
    for (k in seq(1,K)) {
        contour(x,y,function.list[[k]]);
    }
}

orthonormal.function.list <- function.list;
orthonormal.function.list <- orthonormal.functions(function.list,
                                                   dx,dy,x,y,
                                                   TRUE);
system.mats <- system.matrices(orthonormal.function.list,
                               dx,dy);

for (n in seq(1,length(data))) {
    par(mfrow=c(3,2));
    problem.parameters <- data[[n]];
    problem.parameters$K.prime <- K.prime;
    problem.parameters$number.terms <- 100;
    l2 <- blackbox(function.list,
                   orthonormal.function.list,
                   system.mats,
                   problem.parameters,
                   dx,dy,
                   TRUE,TRUE);
    ## cat("Press [ENTER] to continue");
    ## line <- readline();
    print(paste("n=", n, "; l2 = ", l2));
}
print(l2);

## n=10;
## data <- sample.process(n, dt, problem.parameters.generate.data);
## for (l in seq(1,length(data))) {
##     data[[l]]$number.terms = 1000;
## }

## l2s <- rep(NA,n);

## par(mfrow=c(sqrt(n),sqrt(n)));
## for (l in seq(1,n)) {
##     l2 <- blackbox(function.list,
##                    orthonormal.function.list,
##                    system.mats,
##                    problem.parameters = data[[l]],
##                    dx,dy,TRUE,TRUE);
##     print(c(l,l2));
##     l2s[l] <- l2;
## }

## print(l2);

## l <- which(l2s==max(l2s));
## problem.parameters <- data[[l]];
## tt1 <- problem.parameters$x.ic^2 / (4*problem.parameters$sigma.2.x);
## tt2 <- (1-problem.parameters$x.ic)^2 / (4*problem.parameters$sigma.2.x);
## tt3 <- problem.parameters$y.ic^2 / (4*problem.parameters$sigma.2.y);
## tt4 <- (1-problem.parameters$y.ic)^2 / (4*problem.parameters$sigma.2.y);
## tt <- sort(c(tt1,tt2,tt3,tt4))[3];

## problem.parameters.x = problem.parameters;
## problem.parameters.x$x.ic = problem.parameters$x.ic;
## problem.parameters.x$a = problem.parameters$ax;
## problem.parameters.x$b = problem.parameters$bx;
## problem.parameters.x$sigma.2 = problem.parameters$sigma.2.x;
## problem.parameters.x$t = tt;

## problem.parameters.y = problem.parameters;
## problem.parameters.y$x.ic = problem.parameters$y.ic;
## problem.parameters.y$a = problem.parameters$ay;
## problem.parameters.y$b = problem.parameters$by;
## problem.parameters.y$sigma.2 = problem.parameters$sigma.2.y;
## problem.parameters.y$t = tt;

## print(c(problem.parameters.x$t, problem.parameters.y$t));
## IC.true <- univariate.solution(x,problem.parameters.x) %*%
##     t(univariate.solution(y,problem.parameters.y));

## IC.coefs <- rep(NA, K);
## for (k in seq(1,K)) {
##     IC.coefs[k] <- sum(apply(IC.true*orthonormal.function.list[[k]],
##                              1,
##                              sum)*dy)*dx;
## }
## IC.coefs <- solve(system.mats$mass.mat, IC.coefs);

## x.ic.index = which(abs(x-problem.parameters$x.ic)<=dx/2);
## y.ic.index = which(abs(y-problem.parameters$y.ic)<=dy/2);

## par(mfrow=c(2,1));
## approx <- rep(0, length(x));
## IC.approx <- matrix(nrow=dim(orthonormal.function.list[[1]])[1],
##                     ncol=dim(orthonormal.function.list[[1]])[2],
##                     0);
## IC.approx <- bivariate.solution.approx(orthonormal.function.list,
##                                        K,
##                                        IC.coefs);
## for (k in seq(1,K)) {
##     if (k==1) {
##         plot(x, IC.coefs[k]*
##                 orthonormal.function.list[[k]][,y.ic.index], type="l",
##              ylim = c(-max(IC.coefs[k]*
##                            orthonormal.function.list[[k]][,y.ic.index]),
##                       max(IC.coefs[k]*
##                           orthonormal.function.list[[k]][,y.ic.index])));
##     } else {
##         lines(x, IC.coefs[k]*
##                  orthonormal.function.list[[k]][,y.ic.index], type="l")
##     }
##     approx <- approx +
##         IC.coefs[k]*
##         orthonormal.function.list[[k]][,y.ic.index];
##     ## IC.approx = IC.approx +
##     ##     IC.coefs[k]*orthonormal.function.list[[k]];
## }
## lines(x, IC.true[,y.ic.index],col="red", lwd = 2);
## lines(x, approx, lty="dashed", col="green", lwd =2);

## approx <- rep(0, length(y));
## for (k in seq(1,K)) {
##     if (k==1) {
##         plot(y, IC.coefs[k]*
##                 orthonormal.function.list[[k]][x.ic.index,], type="l",
##              ylim = c(-max(5*IC.coefs[k]*
##                           orthonormal.function.list[[k]][x.ic.index,]),
##                       max(5*IC.coefs[k]*
##                           orthonormal.function.list[[k]][x.ic.index,])));
##     } else {
##         lines(y, IC.coefs[k]*
##                  orthonormal.function.list[[k]][x.ic.index,], type="l")
##     }
##     approx <- approx +
##         IC.coefs[k]*
##                 orthonormal.function.list[[k]][x.ic.index,];
## }
## lines(y, IC.true[x.ic.index,],col="red", lwd=2);
## lines(y, approx, lty="dashed", col="green", lwd=2);
    
## par(mfrow=c(2,1));
## plot(x, IC.approx[,y.ic.index], type = "l");
## lines(x,IC.true[,y.ic.index], col = "red");
## abline(v=problem.parameters$x.ic, col = "red");

## plot(y, IC.approx[x.ic.index,], type = "l");
## lines(y,IC.true[x.ic.index,], col = "red");
## abline(v=problem.parameters$y.ic, col = "red");

## IC.vec <- IC.coefs;
