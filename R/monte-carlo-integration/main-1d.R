rm(list=ls());
source("1-d-solution.R");

problem.parameters = NULL;
problem.parameters$a = 0;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.8;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2 = 1;
problem.parameters$t = 0.1;

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

poly.degree.x = 5;
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

stiff.mat <- matrix(nrow=poly.degree.x+1, ncol=poly.degree.x+1);
dx = 1e-5;
x = seq(problem.parameters$a, problem.parameters$b, by=dx);
for (i in seq(1,poly.degree.x+1)) {
    current.basis.dx.i = as.function(deriv(basis.function.poly(i=i),"x"),
                                     vector=FALSE);

    for (j in seq(i,poly.degree.x+1)) {
        current.basis.dx.j = as.function(deriv(basis.function.poly(i=j),"x"),
                                     vector=FALSE);

        stiff.matrix.entry = sum(current.basis.dx.i(x)*
                                 current.basis.dx.j(x)*dx);
        print(c(i,j,stiff.matrix.entry));
        stiff.mat[i,j]=1/2*problem.parameters$sigma.2*stiff.matrix.entry;
        stiff.mat[j,i]=1/2*problem.parameters$sigma.2*stiff.matrix.entry;
    }
}

mass.mat <- matrix(nrow=poly.degree.x+1, ncol=poly.degree.x+1);
for (i in seq(1,poly.degree.x+1)) {
    current.basis.i = as.function(basis.function.poly(i=i),
                                  vector=FALSE);

    for (j in seq(i,poly.degree.x+1)) {
        current.basis.j = as.function(basis.function.poly(i=j),
                                  vector=FALSE);

        mass.matrix.entry = sum(current.basis.i(x)*
                                 current.basis.j(x)*dx);
        print(c(i,j,mass.matrix.entry));
        mass.mat[i,j]=mass.matrix.entry;
        mass.mat[j,i]=mass.matrix.entry;
    }
}
eig <- eigen(solve(mass.mat) %*% stiff.mat);

b = rep(NA, poly.degree.x+1);
for (i in seq(1,poly.degree.x+1)) {
    b[i] = basis.function(problem.parameters$x.ic,i);
}
IC.vec = solve(mass.mat, b);
IC.vec = b;

x=seq(problem.parameters$a,problem.parameters$b,length.out=1000);
IC <- rep(0,length(x));
for (i in seq(1,poly.degree.x+1)) {
    IC = IC +
        IC.vec[i]*basis.function(x,i);
}
plot(x, IC, type="l");

coefs = (eig$vectors) %*% diag(exp(-eig$values * problem.parameters$t)) %*% t(eig$vectors) %*% IC.vec;

TT = 10000;
dt = problem.parameters$t/TT;
expAdt = (eig$vectors) %*% diag(exp(-eig$values * TT*dt)) %*% t(eig$vectors);
A = solve(mass.mat) %*% stiff.mat;
expAdt.approx.dt =
    diag(1,poly.degree.x+1) - A*dt - 0.5*A^2*dt^2;
expAdt.approx = expAdt.approx.dt^TT;

coefs.approx = expAdt.approx %*% IC.vec;

plot(x,univariate.solution.approx(x,coefs.approx),type="l");
lines(x,univariate.solution(x,problem.parameters), col="red",
      lty="dashed");

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

