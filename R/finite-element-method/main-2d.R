rm(list=ls());
library("mvtnorm");
source("2-d-solution.R");
source("../classical-solution/2-d-solution.R");

PLOT.SOLUTION = TRUE;
dx = 1.0/400
dy = 1.0/400
K.prime = 12;


ax=-0.010240
x_T=-0.003170
bx=0.002528
## ## 
ay=-0.010083;
y_T=-0.003316;
by=0.003052;

Lx = bx-ax
Ly = by-ay

problem.parameters.generate.data = NULL;
problem.parameters.generate.data$t <- 1;
problem.parameters.generate.data$sigma.2.x <- (0.003941/Lx)^2; 
problem.parameters.generate.data$sigma.2.y <- (0.006959/Ly)^2;
problem.parameters.generate.data$rho <- 0.6;
problem.parameters.generate.data$x.ic <- 0;
problem.parameters.generate.data$y.ic <- 0;
dt <- problem.parameters.generate.data$t/1000;
n.samples <- 64;

data <- sample.process(n.samples, dt, problem.parameters.generate.data);

data[[1]]$ax = ax/Lx;
data[[1]]$x.fc = x_T/Lx;
data[[1]]$bx = bx/Lx;

data[[1]]$ay = ay/Ly;
data[[1]]$y.fc = y_T/Ly;
data[[1]]$by = by/Ly;

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

sigma=0.40;
sigma2=sigma^2;
l=1;
function.list <-
    basis.functions.normal.kernel(rho=problem.parameters.generate.data$rho,
                                  l=l,
                                  sigma2=sigma^2,
                                  dx,dy,
                                  std.dev.factor=0.6);
orthonormal.function.list <- function.list;
orthonormal.function.list <- orthonormal.functions(function.list,
                                                   dx,dy,x,y,
                                                   FALSE);
system.mats <- system.matrices(orthonormal.function.list,
                               dx,dy);
## ll.sum = 0;
## for (n in seq(1,length(data))) {
##     print ( paste("ON DATA POINT ", n) );
##     {
##         par(mfrow=c(3,2));
##         problem.parameters.original <- data[[n]];
##         problem.parameters.original$K.prime <- K.prime;
##         problem.parameters.original$number.terms <- 100;
##         problem.parameters.original$rho <- 0.6;
        
##         l2.1 <- blackbox(function.list,
##                          orthonormal.function.list,
##                          system.mats,
##                          problem.parameters.original,
##                          dx,dy,
##                          TRUE,TRUE);
##         print (l2.1 *
##                1.0/((problem.parameters.original$bx-
##                  problem.parameters.original$ax) *
##                  (problem.parameters.original$by-
##                   problem.parameters.original$ay)));
##         sol = l2.1 *
##             1.0/((problem.parameters.original$bx-
##                   problem.parameters.original$ax) *
##                  (problem.parameters.original$by-
##                   problem.parameters.original$ay));
##     }

n = 1
    a.indeces = c(1,-1);
    b.indeces = c(0, 1);
    c.indeces = c(0,-1);
    d.indeces = c(0, 1);

    a.power=0;
    b.power=0;
    c.power=0;
    d.power=0;

    xs <- seq(0.1 + ax/Lx, -0.01 + bx/Lx, length.out = 20);

    h.ay = exp(-6.5);
    h.by = exp(-5.5);
    h.bx = exp(-7.0);

    ## h.ay = exp(-5.5);
    ## h.by = exp(-5.5);
    ## h.bx = exp(-5.5);

    ## hs = c(-0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03)
    hs = exp(seq(-8,-3,length.out=5))

    plot.data <- matrix(nrow = length(hs),
                        ncol = length(xs));
    for (h in hs) {
        derivs = rep(NA, length(xs))
        for (xx in xs) {
            data[[1]]$x.fc = xx
            derivative = 0;
            
            a.power=0;
            b.power=0;
            c.power=0;
            d.power=0;
            
            for ( i in c(1,2)) {
                if (i==1) { a.power=1; } else { a.power=0; };
                
                for ( j in c(2)) {
                    if (j==1) { b.power=1; } else { b.power=0; };
                    
                    for ( k in c(2)) {
                        if (k==1) { c.power=1; } else { c.power=0; };
                        
                        for ( l in c(2)) {
                            if (l==1) { d.power=1; } else { d.power=0; };
                            
                            problem.parameters.original <- data[[n]];
                            problem.parameters.original$K.prime <- K.prime;
                            problem.parameters.original$number.terms <- 100;
                            
                            problem.parameters.original$ax <-
                                problem.parameters.original$ax + a.indeces[i]*h;
                            ## problem.parameters.original$bx <-
                            ##     problem.parameters.original$bx + b.indeces[j]*h.bx;
                            
                            ## problem.parameters.original$ay <-
                            ##     problem.parameters.original$ay + c.indeces[k]*h.ay;
                            ## problem.parameters.original$by <-
                            ##     problem.parameters.original$by + d.indeces[l]*h.by;
                            
                            current.sol  <- blackbox(function.list,
                                                     orthonormal.function.list,
                                                     system.mats,
                                                     problem.parameters.original,
                                                     dx,dy,
                                                     FALSE, FALSE);
                            
                            L.x <- (problem.parameters.original$bx-
                                    problem.parameters.original$ax);
                            L.y <- (problem.parameters.original$by-
                                    problem.parameters.original$ay);
                            
                            current.sol = current.sol *
                                1.0/(L.x * L.y);
                            
                            derivative = derivative +
                                current.sol * ( (-1)^a.power *
                                                     (-1)^b.power *
                                                          (-1)^c.power *
                                                               (-1)^d.power );
                        }
                    }
                }
            }
            print(c(derivative, h));
            print(c(which(hs==h), which(xs==xx), derivative/(2*h)))
            derivs[which(xs == xx)] = derivative/(2*h)
        }
        plot.data[which(hs==h),] = derivs
        plot(xs, derivs)
    }

    dev.off()
    y.lim <- c(min(apply(X=plot.data, MARGIN=1, FUN=min, na.rm=TRUE)),
               max(apply(X=plot.data, MARGIN=1, FUN=max, na.rm=TRUE)))
    colors <- rgb(red=0,blue=0,green=seq(1,length(hs))/length(hs))

    par(mfcol=c(1,2))
    for (h in hs) {
        index <- which(hs==h);
        if (index == 1) {
            plot(xs, plot.data[index,], type = "l",
                 ylim = y.lim,
                 col = colors[index],
                 lwd=2)
        } else {
            lines(xs, plot.data[index,],
                  col = colors[index],
                  lwd=2)
        }
    }
    legend(x = "topright", legend = paste("h=", hs), col=colors, lwd = 2)
    abline(h=0, lwd=1)
    lines(xs, plot.data[9,], col = "red", lwd = 2) 
    
    colors <- rgb(red=0,blue=0,green=seq(1,length(xs))/length(xs))
    for (x in xs) {
        index <- which(xs == x)
        if (index == 1) {
            plot(log(hs[seq(1,length(hs))]), plot.data[seq(1,length(hs)),index], type = "l",
                 ylim = y.lim,
                 col = colors[index],
                 lwd=2)
            points(log(hs[seq(1,length(hs))]), plot.data[seq(1,length(hs)),index],
                   col = colors[index])
        } else {
            lines(log(hs[seq(1,length(hs))]), plot.data[seq(1,length(hs)),index],
                  col = colors[index],
                  lwd=2)
            points(log(hs[seq(1,length(hs))]), plot.data[seq(1,length(hs)),index],
                  col = colors[index])
        }
    }
    legend(x = "topright", legend = paste("x=", xs), col=colors, lwd = 2)


##     plot(log(hs), apply(X=plot.data, MARGIN=1, FUN=function(x) { return( sqrt(sum(x^2)) )} ))
  
##     for ( i in seq(1,2)) {
##         if (i==0) { a.power=1; } else { a.power=0; };
        
##         for ( j in seq(1,2)) {
##             if (j==0) { b.power=1; } else { b.power=0; };
            
##             for ( k in seq(1,2)) {
##                 if (k==0) { c.power=1; } else { c.power=0; };
                
##                 for ( l in seq(1,2)) {
##                     if (l==0) { d.power=1; } else { d.power=0; };
                    
##                     problem.parameters.original <- data[[n]];
##                     problem.parameters.original$K.prime <- K.prime;
##                     problem.parameters.original$number.terms <- 100;
                    
##                     problem.parameters.original$ax <-
##                         problem.parameters.original$ax + a.indeces[i]*h;
##                     problem.parameters.original$bx <-
##                         problem.parameters.original$bx + b.indeces[i]*h;
##                     problem.parameters.original$ay <-
##                         problem.parameters.original$ay + c.indeces[i]*h;
##                     problem.parameters.original$by <-
##                         problem.parameters.original$by + d.indeces[i]*h;
                    
##                     current.sol  <- blackbox(function.list,
##                                              orthonormal.function.list,
##                                              system.mats,
##                                              problem.parameters.original,
##                                              dx,dy,
##                                              FALSE, FALSE);

##                     L.x <- (problem.parameters.original$bx-
##                             problem.parameters.original$ax);
##                     L.y <- (problem.parameters.original$by-
##                             problem.parameters.original$ay);

##                     conversion.factor <-
##                         (sqrt(problem.parameters.original$sigma.2.x)/
##                          L.x^4) *
##                         (sqrt(problem.parameters.original$sigma.2.y)/
##                          L.y^4);
                    
                    
##                     current.sol = current.sol *
##                         ## 1.0/(L.x * L.y);
##                         conversion.factor;
##                     ## print(current.sol);
                    
##                     derivative = derivative +
##                         current.sol *
##                         (-1)^a.power * (-1)^b.power * (-1)^c.power * (-1)^d.power;
##                 }
##             }
##         }
##     }
##     print(derivative)
##     full.derivative = derivative / (h^4);
##     print(full.derivative)
    
##     dmvnorm(x = c(data[[n]]$x.fc,
##                   data[[n]]$y.fc),
##             mean = rep(0,2),
##             sigma = matrix(nrow=2,ncol=2,
##                            data = c(problem.parameters.generate.data$sigma.2.x,
##                                     sqrt(problem.parameters.generate.data$sigma.2.x)*
##                                     sqrt(problem.parameters.generate.data$sigma.2.y)*
##                                     problem.parameters.generate.data$rho,
##                                     sqrt(problem.parameters.generate.data$sigma.2.x)*
##                                     sqrt(problem.parameters.generate.data$sigma.2.y)*
##                                     problem.parameters.generate.data$rho,
##                                     problem.parameters.generate.data$sigma.2.y)),
##             log=FALSE);
    
   
##     if (derivative < 0) {
##         print(paste("DATA POINT ", n , "PRODUCED NEG LIKELIHOOD"));
##     } else {
##         ll.sum = ll.sum + log(derivative);
##     }
##     ## print(paste("n=", n, "; l2.1 = ", l2.1));   
##     ## cat("Press [ENTER] to continue");
##     ## line <- readline();
## }

## print(l2);

## kernel <- function(x,x.0,tt) {
##     if (x.0 <= 0.5) {
##         out <- (dnorm(t(x),
##                       mean=x.0,
##                       sd=sigma.x*sqrt(tt)) -
##                 dnorm(t(x),
##                       mean=-x.0,
##                       sd=sigma.x*sqrt(tt)));
##     } else {
##         out <- (dnorm(t(x),
##                       mean=x.0,
##                       sd=sigma.x*sqrt(tt)) -
##                 dnorm(t(x),
##                       mean=1+(1-x.0),
##                       sd=sigma.x*sqrt(tt)));
##     }
##     return(out);
## }


## plot.new();
## dev.off();

## nn <- 11;
## nu <- nn-1;
## xs <- seq(1,nn-1)/nn;
## problem.parameters = rescale.problem(problem.parameters.original);
## sigma.x <- sqrt(problem.parameters$sigma.2.x);
## sigma.y <- sqrt(problem.parameters$sigma.2.y);
## rho <- problem.parameters$rho;

## tts = ifelse(xs <= 0.5, ((1-xs)/(5*sigma.x))^2, (xs/(5*sigma.x))^2)
## tt <- min(tts);
## xx <- seq(0,1,by=dx/10);
## x.i <- 1/nn;
## x.i.prime <- (nn-1)/nn;

## plot(xx, choose(nn,nu)*x^(nu)*(1-x)^(nn-nu), type="l", ylim = c(0,1));
## density.factor = 3;
## integrals.exact <- rep(0, nn-1);
## integrals.approx <- rep(0, nn-1);

## for (ii in seq(1,nn-1)) {
##     lines(xx, kernel(xx, x.0=ii/nn, tt=tt), col="red");
##     for (j in seq(1,length(xx))) {
##         integrals.exact[ii] = integrals.exact[ii] +
##             (kernel(x=ii/nn, x.0 = xx[j], tt=tt)*
##                   choose(nn,nn-1)*xx[j]^(nn-1)*(1-xx[j])^1);
##     }
##     integrals.exact[ii] = integrals.exact[ii] * dx/10;

##     xs <- seq(0,1,length.out=nn*density.factor);
##     weights <- choose(nn,ii)*
##         (xs)^(nu)*(1-xs)^(nn-nu);
##     weights <- weights/sum(weights)*1/(nn+1);
##     for (jj in seq(1,length(xs))) {
##         integrals.approx[ii] = integrals.approx[ii] +
##             kernel(ii/nn, x.0=xs[jj], tt=tt)*weights[jj]
##     }
## }
## print(integrals.approx);
## print(integrals.exact);

## par(mfrow=c(2,1));
## plot(integrals.exact, integrals.approx);
## abline(a=0,b=1,lwd=2);

## plot((integrals.approx-integrals.exact)/integrals.exact*100);
