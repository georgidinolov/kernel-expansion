## This script generates figures which compare the profile of the
## likelihood with respect to resolution in the ``small'' direction
rm(list=ls());
library("mvtnorm");
source("2-d-solution.R");
library("latex2exp");
source("../classical-solution/2-d-solution.R");

produce.likelihood.surfaces <- function(datum,
                                        orthonormal.function.list,
                                        system.mats,
                                        rho,
                                        dx,
                                        h,
                                        sigma.ys,
                                        ts) {
    source("2-d-solution.R");
    source("../classical-solution/2-d-solution.R");
    
    PLOT.SOLUTION = TRUE;
    dy = dx;
    K.prime = 12;
    
    sigma_x = sqrt(datum$sigma.2.x);
    sigma_y = sqrt(datum$sigma.2.y);
    
    ax = datum$ax;
    x_T = datum$x.fc;
    bx = datum$bx;
    ## ##
    ay = datum$ay;
    y_T = datum$y.fc;
    by = datum$by;
    
    Lx = bx-ax;
    Ly = by-ay;
    
    problem.parameters.generate.data = NULL;
    problem.parameters.generate.data$t <- 1;
    problem.parameters.generate.data$sigma.2.x <- 1.0;
    problem.parameters.generate.data$sigma.2.y <- 1.0;
    problem.parameters.generate.data$rho <- 0.0;
    problem.parameters.generate.data$x.ic <- 0;
    problem.parameters.generate.data$y.ic <- 0;
    dt <- problem.parameters.generate.data$t/1000;
    n.samples <- 1;
    
    data <- sample.process(1, dt, problem.parameters.generate.data);
    
    data[[1]]$ax = ax/Lx;
    data[[1]]$x.fc = x_T/Lx;
    data[[1]]$bx = bx/Lx;
    
    data[[1]]$ay = ay/Ly;
    data[[1]]$y.fc = y_T/Ly;
    data[[1]]$by = by/Ly;
    data[[1]]$number.terms = 1000;
    
    ## data <- load.data.from.csv(
    ##     "~/research/PDE-solvers/src/brownian-motion/data-set-2.csv");
    
    problem.parameters <- data[[1]];
    problem.parameters$K.prime <- K.prime;
    problem.parameters$number.terms <- 100;
    terms <- 100;
    
    K <- (K.prime-1)^2;
    function.list <- vector("list", K);
    


    
    a.indeces = c( 0, -1);
    b.indeces = c( 1,  0);
    c.indeces = c( 0, -1);
    d.indeces = c( 1,  0);
    
    a.power=0;
    b.power=0;
    c.power=0;
    d.power=0;
    
    xs <- seq(exp(-2.5) + ax/Lx, -exp(-2.5) + bx/Lx, length.out = 20);
    ys <- seq(exp(-2.5) + ay/Ly, -exp(-2.5) + by/Ly, length.out = 20);
    
    ## h.ay = exp(-2.7)/1;
    ## h.by = exp(-2.7)/1;
    ## h.bx = exp(-2.7)/1;
    ## h.ax = exp(-2.7)/1;
    h.ay = h;
    h.by = h;
    h.bx = h;
    h.ax = h;

    likelihoods <- matrix(NA,
                          nrow=length(sigma.ys),
                          ncol=length(ts));
    
    analytic.likelihoods <- matrix(NA,
                                   nrow=length(sigma.ys),
                                   ncol=length(ts));
    
    for (t in ts) {
        for (sigma.y in sigma.ys) {
            
            data[[1]]$rho <- rho;
            data[[1]]$t <- t;
            data[[1]]$sigma.2.y <- sigma.y^2;
            
            derivative = 0
            analytic.derivative = 0;
            
            for ( i in c(1,2)) {
                if (i==1) { a.power=1; } else { a.power=0; };
                
                for ( j in c(1,2)) {
                    if (j==1) { b.power=1; } else { b.power=0; };
                    
                    for ( k in c(1,2)) {
                        if (k==1) { c.power=1; } else { c.power=0; };
                        
                        for ( l in c(1,2)) {
                            if (l==1) { d.power=1; } else { d.power=0; };
                            
                            
                            problem.parameters.original <- data[[1]]
                            problem.parameters.original$K.prime <- 12^2
                            problem.parameters.original$number.terms <- 1000;
                            
                            problem.parameters.original$ax <-
                                problem.parameters.original$ax + a.indeces[i]*h.ax;
                            problem.parameters.original$bx <-
                                problem.parameters.original$bx + b.indeces[j]*h.bx;
                            
                            problem.parameters.original$ay <-
                                problem.parameters.original$ay + c.indeces[k]*h.ay;
                            problem.parameters.original$by <-
                                problem.parameters.original$by + d.indeces[l]*h.by;

                            current.sol <- 0;
                            current.sol  <- blackbox(function.list,
                                                     orthonormal.function.list,
                                                     system.mats,
                                                     problem.parameters.original,
                                                     dx,dy,
                                                     FALSE, FALSE);
                            
                            analytic.sol <- analytic.solution(problem.parameters.original);
                            
                            L.x <- (problem.parameters.original$bx-
                                    problem.parameters.original$ax);
                            L.y <- (problem.parameters.original$by-
                                    problem.parameters.original$ay);
                            
                            current.sol = current.sol *
                                1.0/(L.x * L.y);
                            
                            analytic.sol = analytic.sol;
                            
                            derivative = derivative +
                                current.sol * ( (-1)^a.power *
                                                     (-1)^b.power *
                                                          (-1)^c.power *
                                                               (-1)^d.power );
                            
                            analytic.derivative = analytic.derivative +
                                analytic.sol * ( (-1)^a.power *
                                                      (-1)^b.power *
                                                           (-1)^c.power *
                                                                (-1)^d.power );
                        }
                    }
                }
            }
            
            finite.difference.factor = 1/(h.ax * h.bx * h.ay * h.by);
            
            derivative <- derivative * finite.difference.factor;
            
            analytic.derivative <- analytic.derivative  * finite.difference.factor;
            ## problem.parameters.original <- data[[1]]
            ## problem.parameters.original$K.prime <- 12^2
            ## problem.parameters.original$number.terms <- 1000;
            
            ## problem.parameters.original$ax <-
            ##     problem.parameters.original$ax;
            ## problem.parameters.original$bx <-
            ##     problem.parameters.original$bx;
            
            ## problem.parameters.original$ay <-
            ##     problem.parameters.original$ay;
            ## problem.parameters.original$by <-
            ##     problem.parameters.original$by;
            ## analytic.derivative <- analytic.derivative.compute(problem.parameters);
            
            print(c(derivative,analytic.derivative,sigma.y,t));
            
            if (!is.nan(derivative) && derivative >= 0.0) {
                ## if (abs((analytic.derivative-derivative)/analytic.derivative) <= 0.05) {
                likelihoods[which(sigma.ys==sigma.y), which(ts==t)] = derivative;
                ## } else {
                ##     likelihoods[which(sigma.ys==sigma.y), which(ts==t)] = NA;
                ## }
            } else {
                likelihoods[which(sigma.ys==sigma.y), which(ts==t)] = NA;
            }
            
            analytic.likelihoods[which(sigma.ys==sigma.y), which(ts==t)] = analytic.derivative;
        }
    }

    out = NULL;
    out$analytic.likelihoods = analytic.likelihoods;
    out$likelihoods = likelihoods;

    return(out);
}

analytic.solution <- function(problem.parameters) {
    t = problem.parameters$t;
    
    ax = problem.parameters$ax;
    bx = problem.parameters$bx;
    x.ic = problem.parameters$x.ic;
    x.fc = problem.parameters$x.fc;
    sigma.2.x = problem.parameters$sigma.2.x;
    
    ay = problem.parameters$ay;
    by = problem.parameters$by;
    y.ic = problem.parameters$y.ic;
    y.fc = problem.parameters$y.fc;
    sigma.2.y = problem.parameters$sigma.2.y;
    
    number.terms = problem.parameters$number.terms;
    mu=0;
    
    if ( x.ic < ax || x.ic > bx) {
        stop ("x.ic not in boundaries");
    }

    if ( y.ic < ay || y.ic > by) {
        stop ("y.ic not in boundaries");
    }

    d1.x = (x.fc - x.ic - 2*seq(-number.terms,number.terms)*(bx-ax))^2 / (2*sigma.2.x*t)
    d2.x = (x.fc + x.ic - 2*ax - 2*seq(-number.terms,number.terms)*(bx-ax))^2 / (2*sigma.2.x*t)

        d1.y = (y.fc - y.ic - 2*seq(-number.terms,number.terms)*(by-ay))^2 / (2*sigma.2.y*t)
    d2.y = (y.fc + y.ic - 2*ay - 2*seq(-number.terms,number.terms)*(by-ay))^2 / (2*sigma.2.y*t)

    out.x = 1/sqrt(2*pi*sigma.2.x*t)*exp(-(mu^2 - 2*mu*(x.fc-x.ic))/(2*sigma.2.x*t))*
        sum(exp(-d1.x) - exp(-d2.x))

    out.y = 1/sqrt(2*pi*sigma.2.y*t)*exp(-(mu^2 - 2*mu*(y.fc-y.ic))/(2*sigma.2.y*t))*
        sum(exp(-d1.y) - exp(-d2.y))
    
    return (out.x * out.y);
};

analytic.derivative.compute <- function(problem.parameters) {
    t = problem.parameters$t;

    ax = problem.parameters$ax;
    bx = problem.parameters$bx;
    x.ic = problem.parameters$x.ic;
    x.fc = problem.parameters$x.fc;
    sigma.2.x = problem.parameters$sigma.2.x;

    ay = problem.parameters$ay;
    by = problem.parameters$by;
    y.ic = problem.parameters$y.ic;
    y.fc = problem.parameters$y.fc;
    sigma.2.y = problem.parameters$sigma.2.y;
    
    number.terms = problem.parameters$number.terms;
    mu=0;
    
    if ( x.ic < ax || x.ic > bx) {
        stop ("x.ic not in boundaries");
    }

    if ( y.ic < ay || y.ic > by) {
        stop ("y.ic not in boundaries");
    }

    n = seq(-number.terms,number.terms);
    
    d1.x = (x.fc - x.ic - 2*seq(-number.terms,number.terms)*(bx-ax))^2 / (2*sigma.2.x*t)
    d2.x = (x.fc + x.ic - 2*ax - 2*seq(-number.terms,number.terms)*(bx-ax))^2 / (2*sigma.2.x*t)

    d1.y = (y.fc - y.ic - 2*seq(-number.terms,number.terms)*(by-ay))^2 / (2*sigma.2.y*t)
    d2.y = (y.fc + y.ic - 2*ay - 2*seq(-number.terms,number.terms)*(by-ay))^2 / (2*sigma.2.y*t)

    out.x = 1/(sqrt(2*pi)* (sigma.2.x*t)^(3/2) )*exp(-(mu^2 - 2*mu*(x.fc-x.ic))/(2*sigma.2.x*t))*
        sum(4*n^2*(2*d1.x-1)*exp(-d1.x) -
            4*n*(n-1)*(2*d2.x-1)*exp(-d2.x))

    out.y = 1/(sqrt(2*pi)* (sigma.2.y*t)^(3/2) )*exp(-(mu^2 - 2*mu*(y.fc-y.ic))/(2*sigma.2.y*t))*
        sum(4*n^2*(2*d1.y-1)*exp(-d1.y) -
            4*n*(n-1)*(2*d2.y-1)*exp(-d2.y));
    
    return (out.x * out.y);
};

sigma.x=0.30;
sigma.y=0.10;
rho=0.9;
std.dev.factor = 1;
    
sigma2.x=sigma.x^2;
sigma2.y=sigma.y^2;
    
l=1;
n=1;

dx = 1/500;
dy = 1/500;

x <- seq(0,1,by=dx);
y <- seq(0,1,by=dy);
    

function.list <-
    basis.functions.normal.kernel(rho=rho,
                                  l=l,
                                  sigma2.x=sigma2.x,
                                  sigma2.y=sigma2.y,
                                  dx,dy,
                                  std.dev.factor=std.dev.factor);
orthonormal.function.list <- function.list;
orthonormal.function.list <- orthonormal.functions(function.list,
                                                   dx,dy,x,y,
                                                   FALSE);
system.mats <- system.matrices(orthonormal.function.list,
                               dx,dy);
save(file="basis-functions-0.3-0.1-0.9.Rdata", list=c("orthonormal.function.list",
                                                       "system.mats",
                                                       "sigma.x",
                                                       "sigma.y",
                                                       "rho"));

load("basis-functions-0.3-0.1-0.9.Rdata");
load("./figures/sample-small-sigma-and-t.Rdata");

datum <- sample.results[[3]]$samples$data[[sample.results[[3]]$min.log.sigma.index]]
likelihood.surfaces <- produce.likelihood.surfaces(datum=datum,
                                                   orthonormal.function.list = orthonormal.function.list,
                                                   system.mats=system.mats,
                                                   rho = 0,
                                                   dx = 1/500,
                                                   h = exp(-2.7),
                                                   sigma.ys = seq(0.1,1,by=0.1),
                                                   ts=exp(seq(-5,1,length.out=15)));
save(file="./figures/likelihood-surfaces-rho-0.9-small-sigma.Rdata", list = c("likelihood.surfaces"));

datum <- sample.results[[3]]$samples$data[[sample.results[[3]]$min.t.index]]
## datum <- sample.results[[2]]$samples$data[[sample(x=seq(1,5000), size=1)]]
likelihood.surfaces <- produce.likelihood.surfaces(datum=datum,
                                                   orthonormal.function.list = orthonormal.function.list,
                                                   system.mats=system.mats,
                                                   rho = 0,
                                                   dx = 1/450,
                                                   h = exp(-2.7),
                                                   sigma.ys = seq(0.1,1,by=0.1),
                                                   ts=exp(seq(-5,1,length.out=15)));
save(file="./figures/likelihood-surfaces-rho-0.9-small-t.Rdata", list = c("likelihood.surfaces"));

rel.error.2 <- matrix(NA, nrow=10,ncol=15)
rel.error.2[4:10, 9:15] = TRUE;
## rel.error.2[c(4,4,5,6),c(10,9,9,9)]=NA;
rel.error <- rel.error.2;

## rel.error <- (abs((likelihood.surfaces$likelihoods - likelihood.surfaces$analytic.likelihoods)/likelihood.surfaces$analytic.likelihoods) <= 0.05);
## rel.error <- rel.error*rel.error.2;

likelihoods <- likelihood.surfaces$likelihoods;
analytic.likelihoods <- likelihood.surfaces$analytic.likelihoods;

for (i in seq(1,dim(rel.error)[1])) {
    for (j in seq(1,dim(rel.error)[2])) {
        if (is.na(rel.error[i,j])) {
            likelihoods[i,j] = NA;
        } else if (rel.error[i,j] == FALSE) {
            likelihoods[i,j] = NA;
        } 
    }
}

sigma.ys <- seq(0.1,1,by=0.1);
ts <- exp(seq(-5,1,length.out=15));

diff <- (likelihoods-analytic.likelihoods)/analytic.likelihoods * 100
filled.contour(sigma.ys, ts, diff,
               xlab = "sigma.ys", ylab = "ts",
               color.palette = heat.colors)


pdf("~/PDE-solvers/src/kernel-expansion/documentation/small-sigma-Galerkin-no-filter-rho-09.pdf", 8,8)
filled.contour(sigma.ys, ts, likelihoods,
               xlab = TeX("$\\tau_y / \\tau_x$"),
               ylab = TeX("$\\tilde{t}$"),
               color.palette = heat.colors)
dev.off()

pdf("~/PDE-solvers/src/kernel-expansion/documentation/small-t-analytic.pdf", 8,8)
filled.contour(sigma.ys, ts, analytic.likelihoods,
               xlab = TeX("$\\tau_y / \\tau_x$"),
               ylab = TeX("$\\tilde{t}$"),
               color.palette = heat.colors)
dev.off()


par(mfrow=c(7,1),mar=c(0,0,0,0));
for (ii in seq(9,15)) {
    plot(sigma.ys, likelihoods[,ii], ylim=c(0,2));
}

par(mfrow=c(7,1),mar=c(0,0,0,0));
for (ii in seq(4,10)) {
    plot(ts, likelihoods[ii,], ylim=c(0,3));
}

f1 <- likelihoods[6,9];
x1 <- sigma.ys[6];

f2 <- likelihoods[7,9];
x2 <- sigma.ys[7]

f3 <- likelihoods[10,9];
x3 <- sigma.ys[10];

neg.alpha.min.1 <- (log(f2/f3) - log(f1/f2)/(1/x1 - 1/x2)*(1/x2 - 1/x3)) /
    ( log(x2/x3) - log(x1/x2)*(1/x2-1/x3)/(1/x1 - 1/x2)  );
alpha <- -1*neg.alpha.min.1 + 1

beta <- -(log(f1/f2) + (alpha+1)*log(x1/x2))/(1/x1 - 1/x2)
CC <- exp(log(f1) + (alpha+1)*log(x1) + beta/x1);

plot(sigma.ys, likelihoods[,9], ylim=c(0,max(likelihoods[,9], na.rm=TRUE)))
lines(sigma.ys, CC*sigma.ys^(-(alpha+1))*exp(-beta/sigma.ys), col = "blue")
lines(sigma.ys, analytic.likelihoods[,9], col="red")
