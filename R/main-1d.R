source("1-d-solution.R");

problem.parameters = NULL;
problem.parameters$a = 0;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.01;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2 = 1;

x = seq(problem.parameters$a,problem.parameters$b,length.out = 100);
ts = seq(0.1, 1, length.out = 10);
for (i in seq(1,length(ts))) {

    problem.parameters$t = ts[i];
    
    if (i==1) {
        y <- univariate.solution(x=x, problem.parameters);
        plot(x,y,type = "l", col = "red");
 
    } else {
        y <- univariate.solution(x=x, problem.parameters);
        lines(x,y);
 
    }
}

t = 0.2;
problem.parameters.x = NULL;
problem.parameters.x$a = -0.5;
problem.parameters.x$b = 1;
problem.parameters.x$x.ic = 0.0;
problem.parameters.x$number.terms = 100;
problem.parameters.x$sigma.2 = 1;
problem.parameters.x$t = t;

x = seq(problem.parameters.x$a, problem.parameters.x$b, length.out = 100);

z <- univariate.solution(x, problem.parameters.x)
plot(x,z, type = "l");

### integral table ###
poly.degree.x = 1;
elementary.integral.table.x = integral.table(problem.parameters.x,
                                             2*(poly.degree.x+1));

## poly bases ##
image.means.weights.x <- image.means.weights(problem.parameters.x);
image.means.x <- image.means.weights.x$image.means;
image.weights.x <- image.means.weights.x$image.weights;

### int x^l y^m f dx dy ###
integrals.table <- rep(NA, 2*poly.degree.x+1);

for (i in seq(1,2*poly.degree.x+1)) {
    integrals.table[i] =
        sum(image.weights.x * elementary.integral.table.x[,i]);
}
K <- integrals.table;

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
                                     integrals.table = K) *
                polynomials.table[[n]];
        }
    }
    P.i.j.sq <- P.i.j^2;
    L2.norm <- sqrt(integrate.polynomial(poly=P.i.j.sq,
                                    integrals.table=integrals.table));
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
                                           integrals.table = integrals.table))
    }
}

k = max(poly.degree.x);
coefficients.table <- vector(mode="list", length=poly.degree.x+1);

for (i in seq(1,poly.degree.x+1)) {
    coef <- mpoly(list(c("x"=0,"coef"=0)));
    for (l in seq(0,k)) {
        coef = coef +
            apply.generator(poly=polynomials.table[[i]],
                            k=l,
                            problem.parameters.x=problem.parameters.x) *
            mpoly(list(c(x=0,coef=t^l/factorial(l))));
        
    }
    coefficients.table[[i]] = as.function(coef)(c(problem.parameters.x$x.ic));
}                     

solution <- rep(nrow=length(x), ncol=length(y), 0);
basis.table <- vector(mode="list", length=poly.degree.x+1);

for (i in seq(1,1)) {
    if (i == 1) {
        current.function = as.function(polynomials.table[[i]]);
        basis.table[[i]] = rep(as.double(unlist(polynomials.table[[1]])),
                               length(x));
    } else {
        current.function = as.function(polynomials.table[[i]], vector=FALSE)
        basis.table[[i]] <- current.function(x);
    }
    solution = solution +
        coefficients.table[[i]]*basis.table[[i]];
}

plot(x, solution*z, type = "l");
lines(x, z, col = "red", lty="dashed");
i=1
j=1
persp(x, y, coefficients.table[[i]][[j]]*basis.table[[i]][[j]]*sqrt(z),
      theta = 50,
      phi = 25)
## contour(x, y, coefficients.table[[2]][[2]]*basis.table[[2]][[2]]*sqrt(z));
