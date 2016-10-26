library("mpoly");
source("1-d-solution.R");
source("bivariate-solution.R");

problem.parameters = NULL;
problem.parameters$a = 0;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.01;
problem.parameters$number.terms = 4;
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

### bivariate solution with no rho ###
t = 0.2;
problem.parameters.x = NULL;
problem.parameters.x$a = -1;
problem.parameters.x$b = 1;
problem.parameters.x$x.ic = 0.5;
problem.parameters.x$number.terms = 10;
problem.parameters.x$sigma.2 = 1;
problem.parameters.x$t = t;

problem.parameters.y = NULL;
problem.parameters.y$a = 0;
problem.parameters.y$b = 0.5;
problem.parameters.y$x.ic = 0.01;
problem.parameters.y$number.terms = 100;
problem.parameters.y$sigma.2 = 1;
problem.parameters.y$t = t;

bivariate.problem.parameters = NULL;
bivariate.problem.parameters$problem.parameters.x = problem.parameters.x;
bivariate.problem.parameters$problem.parameters.y = problem.parameters.y;
bivariate.problem.parameters$rho = 0.0;

x = seq(problem.parameters.x$a, problem.parameters.x$b, length.out = 100);
y = seq(problem.parameters.y$a, problem.parameters.y$b, length.out = 100);

z <- bivariate.solution.no.rho(x, y, bivariate.problem.parameters)
persp(x,y,sqrt(z));
contour(x,y,sqrt(z), xlab = "x", ylab = "y");

### integral table ###
poly.degree.x = 2;
poly.degree.y = 2;

elementary.integral.table.x = integral.table(problem.parameters.x, 2*poly.degree.x);
elementary.integral.table.y = integral.table(problem.parameters.y, 2*poly.degree.y);

## poly bases ##
image.means.weights.x <- image.means.weights(problem.parameters.x);
image.means.x <- image.means.weights.x$image.means;
image.weights.x <- image.means.weights.x$image.weights;

image.means.weights.y <- image.means.weights(problem.parameters.y);
image.means.y <- image.means.weights.y$image.means;
image.weights.y <- image.means.weights.y$image.weights;

### int x^l y^m f dx dy ###
integrals.table <- matrix(nrow=2*poly.degree.x+1, ncol = 2*poly.degree.y+1);

for (i in seq(1,2*poly.degree.x+1)) {
    for (j in seq(1,2*poly.degree.y+1)) {
     
   integrals.table[i,j] = sum(image.weights.x * elementary.integral.table.x[,i]) *
            sum(image.weights.y * elementary.integral.table.y[,j]);
    }
}
K <- integrals.table;

x.matrix <- matrix(nrow = length(x), ncol = length(y), data = rep(x,length(y)));
y.matrix <- matrix(nrow = length(x), ncol = length(y), data = rep(x,length(y)), byrow=TRUE)

polynomials.table <- vector(mode="list", length=poly.degree.x+1);
for (i in seq(1,length(Polynomials.table))) {
    polynomials.table[[i]] = vector(mode="list", length=poly.degree.y+1);
}

for (i in seq(1,poly.degree.x+1)) {
    for (j in seq(1,poly.degree.y+1)) {
        ## pp = x^(i-1)y^(j-1)
        pp <- mpoly(list( c("x"=i-1,y=j-1,"coef"=1)));

        ## P.i.j = x^(i-1)y^(j-1) - <x^(i-1)y^(j-1) | P_{1,1}>P_{1,1} - <x^(i-1)y^(j-1) | P_{1,2}>P_{1,2} ...
        if (i==1 && j==1) {
            P.i.j <- pp;
        } else {
            P.i.j <- pp;
            N = (i-1)*(poly.degree.y+1) + j - 1;
            for (n in seq(1,N)) {
                i.prime = ceiling(n/(poly.degree.y+1));
                j.prime = n - (poly.degree.x+1)*(i.prime-1);
                P.i.j <- P.i.j -
                    integrate.polynomial(poly = P.i.j*polynomials.table[[i.prime]][[j.prime]],
                                         integrals.table = K) *
                    polynomials.table[[i.prime]][[j.prime]];
            }
        }
        P.i.j.sq <- P.i.j^2;
        L2.norm <- sqrt(integrate.polynomial(poly=P.i.j.sq,
                                             integrals.table=integrals.table));
        print(c(i,j,L2.norm));
        polynomials.table[[i]][[j]] <- mpoly(list(c("x"=0,"y"=0,"coef"=1/L2.norm))) * P.i.j;
    }
}

### check orthogonality ###
N = (poly.degree.x+1)*(poly.degree.y+1);
orthogonality.matrix <- matrix(nrow=N, ncol = N);

for (n in seq(1,(poly.degree.x+1)*(poly.degree.y+1))) {
    for (m in seq(1,(poly.degree.x+1)*(poly.degree.y+1))) {
        ## ##
        i.n = ceiling(n/(poly.degree.y+1));
        j.n = n - (poly.degree.x+1)*(i.n-1);
        ## ##
        i.m = ceiling(m/(poly.degree.y+1));
        j.m = m - (poly.degree.x+1)*(i.m-1);n
        ## ##

        orthogonality.matrix[n,m] =
            as.double(integrate.polynomial(poly=polynomials.table[[i.n]][[j.n]]*polynomials.table[[i.m]][[j.m]],
                                           integrals.table = integrals.table))
    }
}

k = max(poly.degree.x, poly.degree.y);
coefficients.table <- vector(mode="list", length=poly.degree.x+1);
for (i in seq(1,length(coefficients.table))) {
    coefficients.table[[i]] = vector(mode="list", length=poly.degree.y+1);
}

for (i in seq(1,poly.degree.x+1)) {
    for (j in seq(1,poly.degree.y+1)) {
        coef <- mpoly(list(c("x"=0,"coef"=0)));
        for (l in seq(0,k)) {
            coef = coef +
                apply.generator(poly=polynomials.table[[i]][[j]],
                                k=l,
                                bivariate.problem.parameters=bivariate.problem.parameters) *
                mpoly(list(c(x=0,coef=t^l/factorial(l))));

        }
        if (i == 1) {
            coefficients.table[[i]][[j]] = as.function(coef)(c(problem.parameters.y$x.ic));
        } else {
            coefficients.table[[i]][[j]] = as.function(coef)(c(problem.parameters.x$x.ic,
                                                               problem.parameters.y$x.ic));
        }
    }
}                     

solution <- matrix(nrow=length(x), ncol=length(y), 0);
basis.table <- vector(mode="list", length=poly.degree.x+1);
for (i in seq(1,length(coefficients.table))) {
    basis.table[[i]] = vector(mode="list", length=poly.degree.y+1);
}

for (i in seq(1,poly.degree.x+1)) {
    for (j in seq(1,poly.degree.y+1)) {

        if (i == 1 && j == 1) {
            current.function = as.function(polynomials.table[[i]][[j]]);
            basis.table[[i]][[j]] = matrix(nrow = length(x),
                                            ncol = length(y),
                                            as.double(unlist(polynomials.table[[1]][[1]])));
        } else {
            current.function = as.function(polynomials.table[[i]][[j]], vector=FALSE)
            if (i == 1) {
                basis.table[[i]][[j]] <- matrix(nrow = length(x),
                                                ncol = length(y),
                                                data = current.function(rep(y,length(x))), byrow = TRUE)
            } else {
                basis.table[[i]][[j]] <- matrix(nrow = length(x),
                                            ncol = length(y),
                                            data = current.function(x=rep(x, length(y)),
                                                                    y=rep(y,length(x))),
                                            byrow = TRUE)
            }
        }
        
        solution = solution +
            coefficients.table[[i]][[j]] * basis.table[[i]][[j]];
    }
}

persp(x, y, solution*z,
      theta = 50,
      phi = 25)
i=1
j=2
persp(x, y, coefficients.table[[i]][[j]]*basis.table[[i]][[j]]*sqrt(z),
      theta = 50,
      phi = 25)
## contour(x, y, coefficients.table[[2]][[2]]*basis.table[[2]][[2]]*sqrt(z));
