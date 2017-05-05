rm(list=ls());
library("mvtnorm");
source("2-d-solution.R");
source("../classical-solution/2-d-solution.R");

PLOT.SOLUTION = TRUE;
dx = 0.005;
dy = 0.005;
K.prime = 12;

problem.parameters.generate.data = NULL;
problem.parameters.generate.data$t <- 1;
problem.parameters.generate.data$sigma.2.x <- 1.0^2;
problem.parameters.generate.data$sigma.2.y <- 1.0^2;
problem.parameters.generate.data$rho <- 0.6;
problem.parameters.generate.data$x.ic <- 0;
problem.parameters.generate.data$y.ic <- 0;
dt <- problem.parameters.generate.data$t/1000;
n.samples <- 64;

data <- sample.process(n.samples, dt, problem.parameters.generate.data);

ax = -0.49715709250150114106858723062032;
x_T = 1.40590589684411759741067271534121;
bx = 1.75535709803343054069557638285914;
ay = -0.55389844337502103233106254265294;
y_T = 1.25288021831189477772738882777048;
by = 1.93415058868214839726817899645539;

data[[1]]$ax = ax;
data[[1]]$x.fc = x_T;
data[[1]]$bx = bx;

data[[1]]$ay = ay;
data[[1]]$y.fc = y_T;
data[[1]]$by = by;

## data <- load.data.from.csv(
##     "~/research/PDE-solvers/src/brownian-motion/data-set-2.csv");

problem.parameters <- data[[1]];
problem.parameters$K.prime <- K.prime;
problem.parameters$number.terms <- 100;

K <- (K.prime-1)^2;
function.list <- vector("list", K);

x <- seq(0,1,by=dx);
y <- seq(0,1,by=dy);

## alpha.minus.1.xs <- rep(NA,K);
## alpha.minus.1.ys <- rep(NA,K);
## for (k in seq(1,K)) {
##     alpha.minus.1.y <- ceiling(k/(K.prime-1));
##     alpha.minus.1.x <- k - ((K.prime-1) * (alpha.minus.1.y-1));
    
##     alpha.minus.1.ys[k] <- alpha.minus.1.y;
##     alpha.minus.1.xs[k] <- alpha.minus.1.x;
## }
## alpha.xs <- alpha.minus.1.xs + 1;
## alpha.ys <- alpha.minus.1.ys + 1;
## ## alpha + beta = K.prime + 2

## x.basis.coors <- alpha.minus.1.xs / (K.prime);
## y.basis.coors <- alpha.minus.1.ys / (K.prime);

## sort.bases <- sort.int(sqrt((x.basis.coors-0.5)^2 +
##                             (y.basis.coors-0.5)^2),
##                        index.return = TRUE);
## alpha.xs <- alpha.xs[sort.bases$ix];
## alpha.ys <- alpha.ys[sort.bases$ix];

## for (k in seq(1,K)) {
##     alpha.x <- alpha.xs[k];
##     alpha.y <- alpha.ys[k];
##     print(c(alpha.x,alpha.y));
##     function.params <- c(alpha.x, alpha.y);
    
##     function.list[[k]] <-
##         basis.function.xy(x,y,
##                           function.params,
##                           problem.parameters);
## }

## par(mfrow=c(ceiling(sqrt(K)),
##             ceiling(sqrt(K))));
## par(mar = c(5,4,2,1));
## if (PLOT.SOLUTION) {
##     for (k in seq(1,K)) {
##         contour(x,y,function.list[[k]]);
##     }
## }

sigma=0.30;
sigma2=sigma^2;
l=1;
function.list <-
    basis.functions.normal.kernel(rho=problem.parameters.generate.data$rho,
                                  l=l,
                                  sigma2=sigma^2,
                                  dx,dy,
                                  std.dev.factor=1/2);
orthonormal.function.list <- function.list;
orthonormal.function.list <- orthonormal.functions(function.list,
                                                   dx,dy,x,y,
                                                   FALSE);
system.mats <- system.matrices(orthonormal.function.list,
                               dx,dy);
ll.sum = 0;
for (n in seq(1,length(data))) {
    print ( paste("ON DATA POINT ", n) );
    {
        par(mfrow=c(3,2));
        problem.parameters.original <- data[[n]];
        problem.parameters.original$K.prime <- K.prime;
        problem.parameters.original$number.terms <- 100;
        problem.parameters.original$rho <- 0.6;
        
        l2.1 <- blackbox(function.list,
                     orthonormal.function.list,
                     system.mats,
                     problem.parameters.original,
                     dx,dy,
                     TRUE,TRUE);
        print (l2.1 *
               1.0/((problem.parameters.original$bx-
                 problem.parameters.original$ax) *
                 (problem.parameters.original$by-
                  problem.parameters.original$ay)));
        sol = l2.1 *
            1.0/((problem.parameters.original$bx-
                  problem.parameters.original$ax) *
                 (problem.parameters.original$by-
                  problem.parameters.original$ay));
    }

    a.indeces = c(-1,0);
    b.indeces = c(0,1);
    c.indeces = c(-1,0);
    d.indeces = c(0,1);
    
    a.power=1;
    b.power=1;
    c.power=1;
    d.power=1;

    h = 1/64;
    derivative = 0;
  
    for ( i in seq(1,2)) {
        if (i==0) { a.power=1; } else { a.power=0; };
        
        for ( j in seq(1,2)) {
            if (j==0) { b.power=1; } else { b.power=0; };
            
            for ( k in seq(1,2)) {
                if (k==0) { c.power=1; } else { c.power=0; };
                
                for ( l in seq(1,2)) {
                    if (l==0) { d.power=1; } else { d.power=0; };
                    
                    problem.parameters.original <- data[[n]];
                    problem.parameters.original$K.prime <- K.prime;
                    problem.parameters.original$number.terms <- 100;
                    
                    problem.parameters.original$ax <-
                        problem.parameters.original$ax + a.indeces[i]*h;
                    problem.parameters.original$bx <-
                        problem.parameters.original$bx + b.indeces[i]*h;
                    problem.parameters.original$ay <-
                        problem.parameters.original$ay + c.indeces[i]*h;
                    problem.parameters.original$by <-
                        problem.parameters.original$by + d.indeces[i]*h;
                    
                    current.sol  <- blackbox(function.list,
                                             orthonormal.function.list,
                                             system.mats,
                                             problem.parameters.original,
                                             dx,dy,
                                             FALSE, FALSE);
                    current.sol = current.sol *
                        1.0/((problem.parameters.original$bx-
                              problem.parameters.original$ax) *
                             (problem.parameters.original$by-
                              problem.parameters.original$ay))
                    
                    derivative = derivative +
                        current.sol *
                        (-1)^a.power * (-1)^b.power * (-1)^c.power * (-1)^d.power;
                }
            }
        }
    }
    derivative = derivative / (h^4);
   
    if (derivative < 0) {
        print(paste("DATA POINT ", n , "PRODUCED NEG LIKELIHOOD"));
    } else {
        ll.sum = ll.sum + log(derivative);
    }
    ## print(paste("n=", n, "; l2.1 = ", l2.1));   
    ## cat("Press [ENTER] to continue");
    ## line <- readline();
}

print(l2);

kernel <- function(x,x.0,tt) {
    if (x.0 <= 0.5) {
        out <- (dnorm(t(x),
                      mean=x.0,
                      sd=sigma.x*sqrt(tt)) -
                dnorm(t(x),
                      mean=-x.0,
                      sd=sigma.x*sqrt(tt)));
    } else {
        out <- (dnorm(t(x),
                      mean=x.0,
                      sd=sigma.x*sqrt(tt)) -
                dnorm(t(x),
                      mean=1+(1-x.0),
                      sd=sigma.x*sqrt(tt)));
    }
    return(out);
}


plot.new();
dev.off();

nn <- 11;
nu <- nn-1;
xs <- seq(1,nn-1)/nn;
problem.parameters = rescale.problem(problem.parameters.original);
sigma.x <- sqrt(problem.parameters$sigma.2.x);
sigma.y <- sqrt(problem.parameters$sigma.2.y);
rho <- problem.parameters$rho;

tts = ifelse(xs <= 0.5, ((1-xs)/(5*sigma.x))^2, (xs/(5*sigma.x))^2)
tt <- min(tts);
xx <- seq(0,1,by=dx/10);
x.i <- 1/nn;
x.i.prime <- (nn-1)/nn;

plot(xx, choose(nn,nu)*x^(nu)*(1-x)^(nn-nu), type="l", ylim = c(0,1));
density.factor = 3;
integrals.exact <- rep(0, nn-1);
integrals.approx <- rep(0, nn-1);

for (ii in seq(1,nn-1)) {
    lines(xx, kernel(xx, x.0=ii/nn, tt=tt), col="red");
    for (j in seq(1,length(xx))) {
        integrals.exact[ii] = integrals.exact[ii] +
            (kernel(x=ii/nn, x.0 = xx[j], tt=tt)*
                  choose(nn,nn-1)*xx[j]^(nn-1)*(1-xx[j])^1);
    }
    integrals.exact[ii] = integrals.exact[ii] * dx/10;

    xs <- seq(0,1,length.out=nn*density.factor);
    weights <- choose(nn,ii)*
        (xs)^(nu)*(1-xs)^(nn-nu);
    weights <- weights/sum(weights)*1/(nn+1);
    for (jj in seq(1,length(xs))) {
        integrals.approx[ii] = integrals.approx[ii] +
            kernel(ii/nn, x.0=xs[jj], tt=tt)*weights[jj]
    }
}
print(integrals.approx);
print(integrals.exact);

par(mfrow=c(2,1));
plot(integrals.exact, integrals.approx);
abline(a=0,b=1,lwd=2);

plot((integrals.approx-integrals.exact)/integrals.exact*100);
