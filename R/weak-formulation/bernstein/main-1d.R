rm(list=ls());
source("1-d-solution.R");

problem.parameters = NULL;
problem.parameters$a = 0;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.9;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2 = 1;
problem.parameters$t = 0.1;

x = seq(problem.parameters$a,problem.parameters$b,
        length.out = 1000);

n.basis = 10;
for (n in seq(1,n.basis)) {
    current.basis <- basis.function.init(n, n.basis);
    if (n==1) {
        plot(x, current.basis(x), type="l",col="red",
             lty="dashed");
    } else {
        lines(x, current.basis(x), col="red",
              lty="dashed");
    }
}
lines(x,univariate.solution(x,problem.parameters),col="black");

## ### testing the path integral START ###
## for (n in seq(1,n.basis)) {
##     current.basis.dt = basis.function.init.dt(means,variances,t.0s,n);
##     if (n==1) {
##         plot(x, current.basis.dt(x,problem.parameters$t), type="l");
##         lines(x, current.basis.dt(x,0), lty="dashed");
##     } else {
##         lines(x, current.basis.dt(x,problem.parameters$t), col="red");
##         lines(x, current.basis.dt(x,0), col="red", lty="dashed");
##     }
## }

## number.samples = 50;
## b.approx = rep(0,n.basis);
## for (n in seq(1,n.basis)) {
##     current.basis = basis.function.init(means,
##                                         variances,
##                                         t.0s,
##                                         n);
##     sampled.bm <- sample.bounded.bm.automatic(problem.parameters,
##                                               delta.t.min=1e-6,
##                                               number.samples=number.samples);
    
##     b.approx[n] = sum(sampled.bm$weights *
##                           current.basis(sampled.bm$points,
##                                         problem.parameters$t))/number.samples;
##     print(c(n,b.approx[n]));
## }

## b.exact = rep(0,n.basis);
## dx = (problem.parameters$b - problem.parameters$a)/length(x);
## for (n in seq(1,n.basis)) {
##     basis.f <- basis.function.init(means,variances,t.0s,n);
##     b.exact[n] = sum(basis.f(x,problem.parameters$t)*
##                     univariate.solution(x,problem.parameters)*dx);
##     print(c(n,b.exact[n]));
## }
## ### testing the path integral END ###

## ### testing system matrix START via Reimann integrals###
## for (n in seq(1,n.basis)) {
##     current.basis.dx = basis.function.init.dx(means,variances,t.0s,n);
##     if (n==1) {
##         plot(x, current.basis.dx(x,problem.parameters$t), type="l");
##         lines(x, current.basis.dx(x,0), lty="dashed");
##     } else {
##         lines(x, current.basis.dx(x,problem.parameters$t), col="red");
##         lines(x, current.basis.dx(x,0), col="red", lty="dashed");
##     }
## }

## Sys.mat <- matrix(nrow=n.basis, ncol=n.basis);
## dt = 1e-3*problem.parameters$t;
## dx = (problem.parameters$b-problem.parameters$a)/length(x);
## ts = seq(0,problem.parameters$t,by=dt);
## for (i in seq(1,n.basis)) {
##     current.basis.dx.i = basis.function.init.dx(means,variances,t.0s,i);

##     for (j in seq(i,n.basis)) {
        
##         current.basis.dx.j = basis.function.init.dx(means,variances,t.0s,j);

##         system.matrix.entry = 0;
##         for (t.current in ts) {
##             system.matrix.entry = system.matrix.entry +
##             sum(current.basis.dx.i(x,t.current)*
##                 current.basis.dx.j(x,t.current)*dx)*dt;
##         }
##         print(c(i,j,system.matrix.entry));
##         Sys.mat[i,j]=1/2*problem.parameters$sigma.2*system.matrix.entry;
##         Sys.mat[j,i]=1/2*problem.parameters$sigma.2*system.matrix.entry;
##     }
## }
## ### testing system matrix END ###

## ### coefs approx START ###
## eig <- eigen(Sys.mat);
## coefs.approx = eig$vectors %*% diag(1/eig$values) %*% t(eig$vectors) %*% b.approx;
## coefs.exact = solve(Sys.mat, b.exact);

## univariate.solution.approx <- function(x,coefs) {
##     out = rep(0,length(x));
##     for (n in seq(1,n.basis)) {
##         out = out +
##             coefs[n]*basis.function.init(means,variances,t.0s,n)(x,problem.parameters$t);
##     }
##     return (out);
## }

## plot(x,univariate.solution(x,problem.parameters),type="l", col="red");
## lines(x,univariate.solution.approx(x,coefs.exact));
## ### coefs approx END ###

### ### ### ### ### ### ###
### ODE approach START ###
stiff.mat <- matrix(nrow=n.basis, ncol=n.basis);
dx = 1e-5;
x = seq(problem.parameters$a, problem.parameters$b, by=dx);
for (i in seq(1,n.basis)) {
    current.basis.dx.i = basis.function.init.dx(i,n.basis);

    for (j in seq(i,n.basis)) {
        current.basis.dx.j = basis.function.init.dx(j,n.basis);

        stiff.matrix.entry = sum(current.basis.dx.i(x)*
                                 current.basis.dx.j(x)*dx);
        print(c(i,j,stiff.matrix.entry));
        stiff.mat[i,j]=1/2*problem.parameters$sigma.2*stiff.matrix.entry;
        stiff.mat[j,i]=1/2*problem.parameters$sigma.2*stiff.matrix.entry;
    }
}

mass.mat <- matrix(nrow=n.basis, ncol=n.basis);
for (i in seq(1,n.basis)) {
    current.basis.i = basis.function.init(i,n.basis);

    for (j in seq(i,n.basis)) {
        current.basis.j = basis.function.init(j,n.basis);

        mass.matrix.entry = sum(current.basis.i(x)*
                                 current.basis.j(x)*dx);
        print(c(i,j,mass.matrix.entry));
        mass.mat[i,j]=mass.matrix.entry;
        mass.mat[j,i]=mass.matrix.entry;
    }
}

## L = t(chol(mass.mat));
## X = forwardsolve(L,stiff.mat);
## A = backsolve(t(L),X);
## eig <- eigen(A);
## eig <- eigen(stiff.mat / ((problem.parameters$b-problem.parameters$a)/2));
eig <- eigen(solve(mass.mat) %*% stiff.mat);
eig <- eigen(stiff.mat);

b = rep(NA, n.basis);
for (i in seq(1,n.basis)) {
    current.basis = basis.function.init(i,n.basis);
    b[i] = current.basis(problem.parameters$x.ic);
}
IC.vec = solve(mass.mat, b);
IC.vec = rep(0,n.basis);
IC.vec[round((n.basis + 1) * problem.parameters$x.ic)] = (n.basis + 1);

x=seq(problem.parameters$a,problem.parameters$b,length.out=1000);
plot(x, IC.approx(x,IC.vec), type="l");

coefs = (eig$vectors) %*% diag(exp(-eig$values * problem.parameters$t)) %*%
    t(eig$vectors) %*% IC.vec;

coefs = t(eig$vectors) %*% IC.vec
for (n in seq(1,n.basis)) {
    if (n==1) {
        plot(x, coefs[n]*basis.function.init(n,n.basis)(x),
             type="l");
    } else {
        lines(x,coefs[n]*basis.function.init(n,n.basis)(x));
    }
}
## coefs = rep(NA,n.basis);
## for (n in seq(1,n.basis)) {
##     coefs[n] = univariate.solution(n/(n.basis+1),problem.parameters)
## }

par(mfrow=c(2,1));
plot(x,univariate.solution(x,problem.parameters),type="l", col="red");
plot(x,univariate.solution.approx(x,coefs),type="l");
lines(x,univariate.solution(x,problem.parameters), col="red");

### ODE approach END ###
