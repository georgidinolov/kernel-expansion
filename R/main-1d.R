rm(list=ls());
source("1-d-solution.R");

problem.parameters = NULL;
problem.parameters$a = 0;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.1;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2 = 1;
problem.parameters$t = 0.5;

alpha.beta <- select.alpha.beta(problem.parameters);
alpha <- alpha.beta$alpha;
beta <- alpha.beta$beta;

## MANUAL ALPHA BETA ##
beta <- 10;
alpha <- max(beta*problem.parameters$x.ic / (1-problem.parameters$x.ic),
             2);
## ## 

kernel <- function(x) {
    return (dbeta(x,alpha,beta));
}

x = seq(problem.parameters$a,problem.parameters$b,length.out = 100);
## plot(x, dnorm(x=x,
##               mean=problem.parameters$x.ic,
##               sd=sqrt(problem.parameters$sigma.2*problem.parameters$t)),
##      type = "l",
##      col = "red",
##      ylab = "y",
##      xlab = "x",
##      ylim = c(0,1));
plot(x, univariate.solution(x, problem.parameters),
     type = "l");
lines(x, dbeta(x=x, shape1=alpha, shape2=beta, log=FALSE))
lines(x, dbeta(x=x, shape1=alpha, shape2=beta, log=FALSE)^2, col = "green")

poly.degree.x = 15;
x.pow.integral.vec <- x.power.integral.vector(problem.parameters,
                                          2*(poly.degree.x+1),
                                          alpha,
                                          beta);
z <- univariate.solution(x, problem.parameters);

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

### plotting some of the basis elements ###
basis.elements <- vector(mode="list", length=poly.degree.x+1);
for (i in seq(1,poly.degree.x+1)) {
    if (i==1) {
        basis.elements[[i]] = rep(as.double(unlist(polynomials.table[[i]])),
                                  length(x)) *
            kernel(x);
    } else {
        current.function <- function(x) {
            return( as.function(polynomials.table[[i]],
                                vector=FALSE)(x) *
                                             kernel(x) )
        };
        
        basis.elements[[i]] = current.function(x);
    }
}

for (i in seq(1,poly.degree.x+1)) {
    if (i==1) {
        plot(x, basis.elements[[i]], type = "l", ylim=c(-2,2));
    } else {
        lines(x, basis.elements[[i]]);
    }
}

### check representation of IC ###
IC <- rep(0, length(x));
for (i in seq(1,poly.degree.x+1)) {
    IC = IC +
        as.function(polynomials.table[[i]],
                    vector=FALSE)(problem.parameters$x.ic)*
                                basis.elements[[i]];
    print(as.function(polynomials.table[[i]],
                      vector=FALSE)(problem.parameters$x.ic));
}
plot(x,IC, type = "l");

### COEFFICIENTS ###
polynomial.kernel <- mpoly(list(c("x"=alpha-1, coef=1/beta(alpha,beta))))*
    (mpoly(list(c("x"=beta-1, coef=1))) + mpoly(list(c("x"=0, coef=-1))))^(beta-1);

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
