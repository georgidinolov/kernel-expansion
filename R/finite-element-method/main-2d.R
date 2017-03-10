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

## problem.parameters.generate.data$sigma.2.x <- 0.986379^2;
## problem.parameters.generate.data$sigma.2.y <- 0.986379^2;
## problem.parameters.generate.data$rho <- -0.708294;

problem.parameters.generate.data$sigma.2.x <- 1^2;
problem.parameters.generate.data$sigma.2.y <- 1^2;
problem.parameters.generate.data$rho <- 0.0;

problem.parameters.generate.data$x.ic <- 0;
problem.parameters.generate.data$y.ic <- 0;
dt <- problem.parameters.generate.data$t/1000;
n.samples <- 100;

data <- load.data.from.csv(
    "~/research/PDE-solvers/data/data-set-1.csv");

data.files.list <- list.files(path = "~/research/PDE-solvers/data",
                              pattern = "data-set-*",
                              full.names = TRUE);
print(data.files.list);
data.files.list <- data.files.list[-c(42,43)];

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
sqrt(mean((estimates.mle.rhos + 0.6)^2));

results.files.list <- list.files(path = "~/research/PDE-solvers/data",
                              pattern = "order-32-rel-tol", full.names = TRUE);
print(results.files.list);
estimates.fd <- rep(NA,length(results.files.list));
for (i in seq(1,length(results.files.list))) {
    mles <- load.results.from.csv(results.files.list[i]);
    if(!is.null(mles)) {
        estimates.fd[i] <- mles$rho;
    }
}
sqrt(mean((estimates.fd[!is.na(estimates.fd)]+0.6)^2));
hist(estimates.fd[!is.na(estimates.fd)], 100, prob = T);
lines(density(estimates.fd[!is.na(estimates.fd)]));
abline(v = c(-0.6, 0.6), lwd = 2, col = "red");

estimates.rogers <- estimator.rodgers(data.files.list, -0.6);

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

sigma=0.25;
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
Delta = 1/128;

for (n in seq(1,length(data))) {
    par(mfrow=c(3,2));
    problem.parameters.original <- data[[n]];
    problem.parameters.original$K.prime <- K.prime;
    problem.parameters.original$number.terms <- 100;

    problem.parameters.original$sigma.2.x <-
        problem.parameters.generate.data$sigma.2.x;

    problem.parameters.original$sigma.2.y <-
        problem.parameters.generate.data$sigma.2.y;

    ## problem.parameters.original$rho <- -0.80;

    problem.parameters.original$rho <-
         problem.parameters.generate.data$rho;
	 
    problem.parameters.current =
         problem.parameters.original;

    l2.1 <- blackbox(function.list,
                     orthonormal.function.list,
                     system.mats,
                     problem.parameters.current,
                     dx,dy,
                     FALSE, FALSE) *
             1.0/((problem.parameters.current$bx-
                 problem.parameters.current$ax)^1 *
                (problem.parameters.current$by-
                 problem.parameters.current$ay)^1);
    print(paste("l2.1 = ", l2.1, sep=""));

    values = rep(NA, 2);
    for (i in c(0,1)) {
        problem.parameters.original <- data[[n]];
        problem.parameters.original$K.prime <- K.prime;
        problem.parameters.original$number.terms <- 100;
        
        problem.parameters.original$sigma.2.x <-
            problem.parameters.generate.data$sigma.2.x;
        
        problem.parameters.original$sigma.2.y <-
            problem.parameters.generate.data$sigma.2.y;
        
        problem.parameters.original$rho <-
            problem.parameters.generate.data$rho;
        
        problem.parameters.current =
            problem.parameters.original;
        
        index = 1 +
            abs(i)*1;
        print(c(index,i));
        
        problem.parameters.current$ax <- 
		      	problem.parameters.original$ax - Delta*i;

        value = (-1)^(i) * 
            blackbox(function.list,
                     orthonormal.function.list,
                     system.mats,
                     problem.parameters.current,
                     dx,dy,
                     FALSE, FALSE) *
             1.0/((problem.parameters.current$bx-
                 problem.parameters.current$ax)^1 *
                (problem.parameters.current$by-
                 problem.parameters.current$ay)^1)

        values[index] = value;
    }
    print(-sum(values)/Delta);
    
    values = rep(NA, 16);
    for (i in c(1,0)) {
    	for (j in c(1,0)) {
	    for (k in c(1,0)) {
	    	for (l in c(1,0)) {
                    problem.parameters.original <- data[[n]];
                    problem.parameters.original$K.prime <- K.prime;
                    problem.parameters.original$number.terms <- 100;
                    
                    problem.parameters.original$sigma.2.x <-
                        problem.parameters.generate.data$sigma.2.x;
                    
                    problem.parameters.original$sigma.2.y <-
                        problem.parameters.generate.data$sigma.2.y;
                    
                    ## problem.parameters.original$rho <- -0.80;
                    
                    problem.parameters.original$rho <-
                        problem.parameters.generate.data$rho;
                    
                    problem.parameters.current =
                        problem.parameters.original;
                    
		    index = 1 +
		    	  abs(i)*1 +
			  abs(j)*2 +
			  abs(k)*4 +
			  abs(l)*8;
	            print(c(index,i,j,k,l));
		    problem.parameters.current = problem.parameters.original;

		    problem.parameters.current$ax <- 
		      	problem.parameters.original$ax - Delta*i;

		    problem.parameters.current$bx <- 
		      	problem.parameters.original$bx + Delta*j;

		    problem.parameters.current$ay <- 
		      	problem.parameters.original$ay - Delta*k;

		    problem.parameters.current$by <- 
		      	problem.parameters.original$by + Delta*l;

		    value = (-1)^(i+j+k+l+1) * 
		    	 blackbox(function.list,
                                  orthonormal.function.list,
                                  system.mats,
                     problem.parameters.current,
                     dx,dy,
                     FALSE, FALSE) *
             1.0/((problem.parameters.current$bx-
                 problem.parameters.current$ax)^3 *
                (problem.parameters.current$by-
                 problem.parameters.current$ay)^3)

		    values[index] = value;
		}
	    }
        }	
    }
    print (sum(values)/Delta^4)

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
