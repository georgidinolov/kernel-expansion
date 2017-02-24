rm(list=ls());
library("parallel");
library("mvtnorm");
source("2-d-solution.R");
source("../classical-solution/2-d-solution.R");

PLOT.SOLUTION = TRUE;
dx = 0.005;
dy = 0.005;
K.prime = 12;
rho.true = -0.8;

problem.parameters.generate.data = NULL;
problem.parameters.generate.data$t <- 1;
problem.parameters.generate.data$sigma.2.x <- 0.986379^2;
problem.parameters.generate.data$sigma.2.y <- 0.986379^2;
problem.parameters.generate.data$rho <- -0.708294;
problem.parameters.generate.data$x.ic <- 0;
problem.parameters.generate.data$y.ic <- 0;
dt <- problem.parameters.generate.data$t/1000;
n.samples <- 100;

data <- load.data.from.csv(
    "~/research/PDE-solvers/data/data-set-1.csv");

data[[1]]$ax = -1.56035;
data[[1]]$x.fc = -0.684587;
data[[1]]$bx = 0.0261297;
data[[1]]$x.ic = 0;

data[[1]]$ay = -0.871711;
data[[1]]$y.fc = -0.154531;
data[[1]]$by = 0.126382;
data[[1]]$y.ic = 0;

data.files.list <- list.files(path = "~/research/PDE-solvers/data",
                             pattern = "data-set-*", full.names = TRUE);
print(data.files.list);

cl <- makeCluster(4);
estimates.mle <- parLapply(cl=cl, X=data.files.list,
                           function(x) {
                               source("2-d-solution.R");
                               library("mvtnorm");
                               data <- load.data.from.csv(x);
                               mles <- mle.estimator.no.boundary(data, 1, 1, 0.0);
                               return(mles);
                           });
stopCluster(cl);

estimates.mle.rhos <- rep(NA, length(estimates.mle));
for (i in seq(1,length(estimates.mle))) {
    estimates.mle.rhos[i] = estimates.mle[[i]]$rho.mle;
}
sqrt(mean((estimates.mle.rhos-0.8)^2));

results.files.list <- list.files(path = "~/research/PDE-solvers/data",
                              pattern = "order-64-rel-tol", full.names = TRUE);
print(results.files.list);
estimates.fd <- rep(NA,length(results.files.list));
for (i in seq(1,length(results.files.list))) {
    mles <- load.results.from.csv(results.files.list[i]);
    if(!is.null(mles)) {
        estimates.fd[i] <- mles$rho;
    }
}
sqrt(mean((estimates.fd[!is.na(estimates.fd)]-0.8)^2));

estimates.rogers <- estimator.rodgers(data.files.list, 0.8);

y.lims <- c(0,
            max(max(density(estimates.fd[!is.na(estimates.fd)])$y),
                max(density(estimates.mle.rhos)$y),
                max(density(estimates.rogers)$y)));

x.lims <- c(min(min(density(estimates.fd[!is.na(estimates.fd)])$x),
                min(density(estimates.mle.rhos)$x),
                min(density(estimates.rogers)$x)),
            max(max(density(estimates.fd[!is.na(estimates.fd)])$x),
                max(density(estimates.mle.rhos)$x),
                max(density(estimates.rogers)$x)));

plot(density(estimates.fd[!is.na(estimates.fd)]),
     ylim=y.lims,
     xlim=x.lims);
lines(density(estimates.mle.rhos),col="red");
lines(density(estimates.rogers),col="green");
abline(v=rho.true, col="red", lwd=2);

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

sigma=0.3;
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
for (n in seq(1,length(data))) {
    par(mfrow=c(3,2));
    problem.parameters.original <- data[[n]];
    problem.parameters.original$K.prime <- K.prime;
    problem.parameters.original$number.terms <- 100;
    problem.parameters.original$sigma.2.x <-
        problem.parameters.generate.data$sigma.2.x;
    problem.parameters.original$sigma.2.y <-
        problem.parameters.generate.data$sigma.2.y;
    problem.parameters.original$rho <-
        problem.parameters.generate.data$rho;
    l2.1 <- blackbox(function.list,
                     orthonormal.function.list,
                     system.mats,
                     problem.parameters.original,
                     dx,dy,
                     TRUE, FALSE);
    
    print (l2.1 *
           1.0/((problem.parameters.original$bx-
                 problem.parameters.original$ax) *
                (problem.parameters.original$by-
                 problem.parameters.original$ay)));
    
    ## problem.parameters.original$ax <- problem.parameters.original$ax - dx;
    ## l2.2 <- blackbox(function.list,
    ##                orthonormal.function.list,
    ##                system.mats,
    ##                problem.parameters.original,
    ##                dx,dy,
    ##                FALSE,FALSE);

    ## -(l2.1-l2.2)/dx;

    print(paste("n=", n, "; l2.1 = ", l2.1));   
    cat("Press [ENTER] to continue");
    line <- readline();
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
