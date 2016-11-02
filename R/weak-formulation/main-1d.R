rm(list=ls());
source("1-d-solution.R");
library("pencopula");

problem.parameters = NULL;
problem.parameters$a = 0;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.8;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2 = 1;
problem.parameters$t = 0.1;

x = seq(problem.parameters$a+0.001,problem.parameters$b-0.001,
        length.out = 1000);

n.basis = 5;
means = c(rep(problem.parameters$x.ic,n.basis));
variances = rep(problem.parameters$sigma.2,n.basis);
t.0s = seq(problem.parameters$t*1, 50*problem.parameters$t*(1-1e-3),
           length.out=n.basis);

means = seq(problem.parameters$a, problem.parameters$b,
            length.out=n.basis+2)[c(-1,-n.basis-2)];
variances = rep(problem.parameters$sigma.2,n.basis);
t.0s = rep(problem.parameters$t*(20e-1),
           length=n.basis);

for (n in seq(1,n.basis)) {
    current.basis = basis.function.init(means,variances,t.0s,n);
    if (n==1) {
        plot(x, current.basis(x,problem.parameters$t), type="l",col="red",
             lty="dashed");
    } else {
        lines(x, current.basis(x,problem.parameters$t), col="red",
              lty="dashed");
    }
}
lines(x,univariate.solution(x,problem.parameters),col="black");
lines(x,0.8*current.basis(x,problem.parameters$t), lty="solid",col="red");

### testing the path integral START ###
for (n in seq(1,n.basis)) {
    current.basis.dt = basis.function.init.dt(means,variances,t.0s,n);
    if (n==1) {
        plot(x, current.basis.dt(x,problem.parameters$t), type="l");
        lines(x, current.basis.dt(x,0), lty="dashed");
    } else {
        lines(x, current.basis.dt(x,problem.parameters$t), col="red");
        lines(x, current.basis.dt(x,0), col="red", lty="dashed");
    }
}

number.samples = 50;
b.approx = rep(0,n.basis);
for (n in seq(1,n.basis)) {
    current.basis = basis.function.init(means,
                                        variances,
                                        t.0s,
                                        n);
    sampled.bm <- sample.bounded.bm.automatic(problem.parameters,
                                              delta.t.min=1e-6,
                                              number.samples=number.samples);
    
    b.approx[n] = sum(sampled.bm$weights *
                          current.basis(sampled.bm$points,
                                        problem.parameters$t))/number.samples;
    print(c(n,b.approx[n]));
}

b.exact = rep(0,n.basis);
dx = (problem.parameters$b - problem.parameters$a)/length(x);
for (n in seq(1,n.basis)) {
    basis.f <- basis.function.init(means,variances,t.0s,n);
    b.exact[n] = sum(basis.f(x,problem.parameters$t)*
                    univariate.solution(x,problem.parameters)*dx);
    print(c(n,b.exact[n]));
}
### testing the path integral END ###

### testing system matrix START via Reimann integrals###
for (n in seq(1,n.basis)) {
    current.basis.dx = basis.function.init.dx(means,variances,t.0s,n);
    if (n==1) {
        plot(x, current.basis.dx(x,problem.parameters$t), type="l");
        lines(x, current.basis.dx(x,0), lty="dashed");
    } else {
        lines(x, current.basis.dx(x,problem.parameters$t), col="red");
        lines(x, current.basis.dx(x,0), col="red", lty="dashed");
    }
}

Sys.mat <- matrix(nrow=n.basis, ncol=n.basis);
dt = 1e-3*problem.parameters$t;
dx = (problem.parameters$b-problem.parameters$a)/length(x);
ts = seq(0,problem.parameters$t,by=dt);
for (i in seq(1,n.basis)) {
    current.basis.dx.i = basis.function.init.dx(means,variances,t.0s,i);

    for (j in seq(i,n.basis)) {
        
        current.basis.dx.j = basis.function.init.dx(means,variances,t.0s,j);

        system.matrix.entry = 0;
        for (t.current in ts) {
            system.matrix.entry = system.matrix.entry +
            sum(current.basis.dx.i(x,t.current)*
                current.basis.dx.j(x,t.current)*dx)*dt;
        }
        print(c(i,j,system.matrix.entry));
        Sys.mat[i,j]=1/2*problem.parameters$sigma.2*system.matrix.entry;
        Sys.mat[j,i]=1/2*problem.parameters$sigma.2*system.matrix.entry;
    }
}
### testing system matrix END ###

### coefs approx START ###
eig <- eigen(Sys.mat);
coefs.approx = eig$vectors %*% diag(1/eig$values) %*% t(eig$vectors) %*% b.approx;
coefs.exact = solve(Sys.mat, b.exact);

univariate.solution.approx <- function(x,coefs) {
    out = rep(0,length(x));
    for (n in seq(1,n.basis)) {
        out = out +
            coefs[n]*basis.function.init(means,variances,t.0s,n)(x,problem.parameters$t);
    }
    return (out);
}

plot(x,univariate.solution(x,problem.parameters),type="l", col="red");
lines(x,univariate.solution.approx(x,coefs.exact));
### coefs approx END ###

### ### ### ### ### ### ###
### ODE approach START ###
stiff.mat <- matrix(nrow=n.basis, ncol=n.basis);
dx = (problem.parameters$b-problem.parameters$a)/length(x);
for (i in seq(1,n.basis)) {
    current.basis.dx.i = basis.function.init.dx(means,variances,t.0s,i);

    for (j in seq(i,n.basis)) {
        current.basis.dx.j = basis.function.init.dx(means,variances,t.0s,j);

        stiff.matrix.entry = sum(current.basis.dx.i(x,problem.parameters$t)*
                                 current.basis.dx.j(x,problem.parameters$t)*dx);
        print(c(i,j,stiff.matrix.entry));
        stiff.mat[i,j]=1/2*problem.parameters$sigma.2*stiff.matrix.entry;
        stiff.mat[j,i]=1/2*problem.parameters$sigma.2*stiff.matrix.entry;
    }
}

mass.mat <- matrix(nrow=n.basis, ncol=n.basis);
dx = (problem.parameters$b-problem.parameters$a)/length(x);
for (i in seq(1,n.basis)) {
    current.basis.i = basis.function.init(means,variances,t.0s,i);

    for (j in seq(i,n.basis)) {
        current.basis.j = basis.function.init(means,variances,t.0s,j);

        mass.matrix.entry = sum(current.basis.i(x,problem.parameters$t)*
                                 current.basis.j(x,problem.parameters$t)*dx);
        print(c(i,j,mass.matrix.entry));
        mass.mat[i,j]=1/2*problem.parameters$sigma.2*mass.matrix.entry;
        mass.mat[j,i]=1/2*problem.parameters$sigma.2*mass.matrix.entry;
    }
}

U = chol(mass.mat);

eig <- eigen(solve(mass.mat)*stiff.mat);

IC.vec <- rep(0,n.basis);
for (i in seq(1,n.basis)) {
    IC.vec[i] = basis.function.init(means,variances,t.0s,i)(problem.parameters$x.ic,
        problem.parameters$t);
}
IC.approx <- function(x,IC.vec) {
    out = rep(0,length(x));
    for (n in seq(1,n.basis)) {
        out = out +
            IC.vec[n]*basis.function.init(means,variances,t.0s,n)(x,
                problem.parameters$t);
    }
    return (out);
}
plot(x, IC.approx(x,IC.vec), type="l");

coefs = eig$vectors %*% diag(exp(-eig$values * problem.parameters$t)) %*%
    t(eig$vectors) %*% IC.vec

plot(x,univariate.solution(x,problem.parameters),type="l", col="red");
lines(x,univariate.solution.approx(x,coefs));
### MASS MATRIX END ###

### ODE approach END ###

alpha.beta <- select.alpha.beta(problem.parameters);
alpha <- alpha.beta$alpha;
beta <- alpha.beta$beta;

## ## MANUAL ALPHA BETA ##
## beta <- 10;
## alpha <- max(beta*problem.parameters$x.ic / (1-problem.parameters$x.ic),
##              2);
## ## ## 
C = 1/0.133086;
K = 3;
## kernel <- function(x) {
##     return (C*exp(-1/(1-x^2)));
## }
kernel <- function(x) {
    return (dbeta(x,alpha,beta));
}

kernel.poly = mpoly(list(c("x"=alpha-1,coef=1))) *
    mpoly(list(c("x"=0,coef=1),c("x"=1,coef=-1)))
x = seq(problem.parameters$a+0.001,problem.parameters$b-0.001,length.out = 100);

plot(x, kernel(x),type="l",
     ylim = c(min(c(univariate.solution(x,problem.parameters),
                    kernel(x))),
              max(c(univariate.solution(x,problem.parameters),
                    kernel(x)))));
lines(x,univariate.solution(x,problem.parameters),col="red");
lines(x,
      univariate.solution(x,
                          problem.parameters)*
      kernel(x),
      col="green");

poly.degree.x = 2;
x.pow.integral.vec <- x.power.integral.vector(problem.parameters,
                                              2*(poly.degree.x+1),
                                              alpha,
                                              beta);
## poly bases ##
polynomials.table <- vector(mode="list", length=poly.degree.x+1);
for (i in seq(1,poly.degree.x+1)) {
    ## pp = x^(i-1)y^(j-1)
    pp <- mpoly(list( c("x"=i-1,"coef"=1)));

    ## P.i.j = x^(i-1)y^(j-1) - <x^(i-1)y^(j-1) | P_{1,1}>P_{1,1} - <x^(i-1)y^(j-1) | P_{1,2}>P_{1,2} ...
    if (i==1) {
        P.i.j <- pp;
    } else {
        P.i.j <- pp;
        N = i-1
        for (n in seq(1,N)) {
            P.i.j <- P.i.j -
                integrate.polynomial(poly = P.i.j*polynomials.table[[n]],
                                     integrals.table = x.pow.integral.vec) *
                polynomials.table[[n]];
        }
    }
    P.i.j.sq <- P.i.j^2;
    L2.norm <- sqrt(integrate.polynomial(poly=P.i.j.sq,
                                    integrals.table=x.pow.integral.vec));
    print(c(i,L2.norm));
    polynomials.table[[i]] <- mpoly(list(c("x"=0,"coef"=1/L2.norm))) * P.i.j;
}

basis.function.poly <- function(i) {
    out = polynomials.table[[i]]*kernel.poly;
    return (out);
}

basis.function <- function(x,i) {
    out = as.function(polynomials.table[[i]],vector=FALSE)(x)*kernel(x);
    return (out);
}


for (i in seq(1,poly.degree.x+1)) {
    if (i==1) {
        plot(x,basis.function(x,i),
             type="l",
             ylim=c(-2,2));
    } else {
        lines(x, basis.function(x,i));
    }
}

### check orthogonality ###
N = (poly.degree.x+1);
orthogonality.matrix <- matrix(nrow=poly.degree.x+1,
                               ncol=poly.degree.x+1);

for (n in seq(1,(poly.degree.x+1))) {
    for (m in seq(1,(poly.degree.x+1))) {

        orthogonality.matrix[n,m] =
            as.double(integrate.polynomial(poly=polynomials.table[[n]]*
                                               polynomials.table[[m]],
                                           integrals.table = x.pow.integral.vec))
    }
}
print(orthogonality.matrix);

### check representation of IC ###
IC <- rep(0, length(x));
for (i in seq(1,poly.degree.x+1)) {
    IC = IC +
        basis.function(problem.parameters$x.ic,i)*
        basis.function(x,i);                        
}
plot(x,IC, type = "l");

### COEFFICIENTS EXACT ### (REIMANN INTEGRALS)
N = 1000;
dx = (problem.parameters$b-problem.parameters$a)/N;

coefs.exact = rep(NA, poly.degree.x+1);
for (i in seq(1,poly.degree.x+1)) {
    integral=sum(univariate.solution(seq(0,N-1)*dx, problem.parameters)*
                 basis.function(seq(0,N-1)*dx,i)*
                                          dx);
    print(c(i,integral));
    coefs.exact[i] = integral;
}

approx.solution = rep(0,length(x));
for (i in seq(1,poly.degree.x+1)) {
    approx.solution = approx.solution +
        coefs.exact[i]*basis.function(x,i);
}

plot(x,univariate.solution(x,problem.parameters),
     type="l",
     col="red");
lines(x,approx.solution);

### COEFFICIENTS APPROX ###
number.samples = 500;
delta.t.min=1e-6;
coefs.approx.1 <- coefs.approx.mc(problem.parameters,
                                  poly.degree.x,
                                  number.samples,
                                  delta.t.min);
coefs.approx.2 <- coefs.approx.mc(problem.parameters,
                                  poly.degree.x,
                                  number.samples,
                                  delta.t.min);
                                  
approx.solution.mc.1 = rep(0,length(x));
approx.solution.mc.2 = rep(0,length(x));
for (i in seq(1,poly.degree.x+1)) {
    approx.solution.mc.1 = approx.solution.mc.1 +
        coefs.approx.1[i]*basis.function(x,i);
        approx.solution.mc.2 = approx.solution.mc.2 +
        coefs.approx.2[i]*basis.function(x,i);
}

plot(x,univariate.solution(x,problem.parameters),
     type="l",
     col="red");
lines(x,approx.solution);
lines(x,approx.solution.mc.1, col="green");
lines(x,approx.solution.mc.2, col="green");

plot(coefs.exact,type="l");
lines(coefs.approx.1,col="red");
lines(coefs.approx.2,col="red");

### COEFFICIENTS APPROX PATH INTEGRAL ###
number.samples = 5;
delta.t = 1e-7;
coefs.approx.1 <- coefs.approx.path.integral.mc(problem.parameters,
                                                poly.degree.x,
                                                number.samples,
                                                delta.t);
coefs.approx.2 <- coefs.approx.path.integral.mc(problem.parameters,
                                                poly.degree.x,
                                                number.samples,
                                                delta.t);
                                  
approx.solution.mc.1 = rep(0,length(x));
approx.solution.mc.2 = rep(0,length(x));
for (i in seq(1,poly.degree.x+1)) {
    approx.solution.mc.1 = approx.solution.mc.1 +
        coefs.approx.1[i]*basis.function(x,i);
        approx.solution.mc.2 = approx.solution.mc.2 +
        coefs.approx.2[i]*basis.function(x,i);
}

plot(x,univariate.solution(x,problem.parameters),
     type="l",
     col="red");
lines(x,approx.solution);
lines(x,approx.solution.mc.1, col="green");
lines(x,approx.solution.mc.2, col="green");

plot(coefs.exact,type="l");
lines(coefs.approx.1,col="red");
lines(coefs.approx.2,col="red");

