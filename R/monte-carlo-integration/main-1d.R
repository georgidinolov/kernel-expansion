rm(list=ls());
source("1-d-solution.R");

problem.parameters = NULL;
problem.parameters$a = 0;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.8;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2 = 1;
problem.parameters$t = 0.05;

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

### EXACT INTEGRAL ###
N = 1000;
dx = (problem.parameters$b-problem.parameters$a)/N;
integral.exact <- sum(univariate.solution(seq(0,N-1)*dx+problem.parameters$a,
                                          problem.parameters)*
                      kernel(seq(0,N-1)*dx+
                             problem.parameters$a)*
                      dx);

### MC INTEGRAL ###
number.samples = 5;
NN = 200;
mc.sims <- vector(mode="list", length=NN);
for (n in seq(1,NN)) {
    sampled.bm <- sample.bounded.bm.automatic(problem.parameters,
                                              delta.t.min=1e-7,
                                              number.samples=number.samples);
    ## sampled.bm <- sample.bounded.bm.accept.reject(problem.parameters, delta.t=1e-6,
    ##                                               number.samples=number.samples);

    plot(density(sampled.bm$points));

    mc.sims[[n]] = sampled.bm;
    print(c(n,sum(sampled.bm$weights *
                  kernel(sampled.bm$points))/number.samples));
}

integrals <- rep(NA,NN);
sample.sizes <- rep(NA,NN);
for (n in seq(1,NN)) {
    integrals[n] = sum(mc.sims[[n]]$weights *
                       kernel(mc.sims[[n]]$points))/number.samples;
    sample.sizes[n] = sum(mc.sims[[n]]$weights);
}
plot(density(integrals));
abline(v=integral.exact,col="red");
plot(integrals);

poly.degree.x = 3;
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
test.function <- function(x) {
    return (x*(1-x));
}

a = problem.parameters$a;
b = problem.parameters$b;
N = 1000;
dx = (b-a)/N;

c.exact <- sum(univariate.solution(seq(0,N-1)*dx, problem.parameters)*
               dx);

c.approx <- test.function(problem.parameters$x.ic) +
    problem.parameters$t * 0.5 * problem.parameters$sigma.2 * -2;  

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
coefs.approx = rep(NA, poly.degree.x+1);
number.samples = 500;
for (i in seq(1,poly.degree.x+1)) {
    sampled.bm <- sample.bounded.bm.automatic(problem.parameters,
                                              delta.t.min=1e-6,
                                              number.samples=number.samples);
    plot(density(sampled.bm$points));
    integral=sum(sampled.bm$weights*
                 basis.function(sampled.bm$points,i))/number.samples;
    print(c(n,integral));
    coefs.approx[i] = integral;
}

approx.solution.mc = rep(0,length(x));
for (i in seq(1,poly.degree.x+1)) {
    approx.solution.mc = approx.solution.mc +
        coefs.approx[i]*basis.function(x,i);
}

plot(x,univariate.solution(x,problem.parameters),
     type="l",
     col="red");
lines(x,approx.solution);
lines(x,approx.solution.mc, col="green");

plot(coefs.exact,type="l");
lines(coefs.approx,col="red");

polynomial.kernel <- mpoly(list(c("x"=alpha-1, coef=1/beta(alpha,beta))))*
    (mpoly(list(c("x"=0, coef=1))) - mpoly(list(c("x"=1, coef=1))) )^(beta-1);

plot(x,kernel(x),type="l", col="red");
lines(x,as.function(polynomial.kernel,vector=FALSE)(x),
      lty="dashed");

coefs <- coefficients(problem.parameters=problem.parameters,
                      polynomials.table=polynomials.table,
                      polynomial.kernel=polynomial.kernel,
                      kernel=kernel,
                      poly.degree=poly.degree.x,
                      number.derivs=20);
                      

solution <- rep(0, length(x));
basis.table <- vector(mode="list", length=poly.degree.x+1);
for (i in seq(1,1)) {
    if (i == 1) {
        current.function = as.function(polynomials.table[[i]]);
        basis.table[[i]] = rep(as.double(unlist(polynomials.table[[i]])),
                               length(x))*
            kernel(x);
    } else {
        current.function = as.function(polynomials.table[[i]], vector=FALSE)
        basis.table[[i]] <- current.function(x)*kernel(x);
    }
    solution = solution +
        coefs[i]*basis.table[[i]];
}

plot(x, solution, type = "l");
lines(x, z, col = "red", lty="dashed");
i=1
j=1
persp(x, y, coefficients.table[[i]][[j]]*basis.table[[i]][[j]]*sqrt(z),
      theta = 50,
      phi = 25)
## contour(x, y, coefficients.table[[2]][[2]]*basis.table[[2]][[2]]*sqrt(z));
