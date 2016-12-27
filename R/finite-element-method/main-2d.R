rm(list=ls());

PLOT.SOLUTION = TRUE;
dx = 0.01;
dy = 0.01;
K.prime = 10;

library("mvtnorm");
source("2-d-solution.R");

problem.parameters.generate.data = NULL;
problem.parameters.generate.data$t <- 1;
problem.parameters.generate.data$sigma.2.x <- 0.1;
problem.parameters.generate.data$sigma.2.y <- 1;
problem.parameters.generate.data$rho <- 0.0;
problem.parameters.generate.data$x.ic <- 0;
problem.parameters.generate.data$y.ic <- 0;
dt <- problem.parameters.generate.data$t/1000;
n.samples <- 1;

data <- sample.process(n.samples, dt, problem.parameters.generate.data);

problem.parameters <- data[[1]];
problem.parameters$K.prime = K.prime;
problem.parameters$number.terms = 1000;

K <- K.prime^2;
function.list <- vector("list", K.prime^2);

x <- seq(problem.parameters$ax,
         problem.parameters$bx,
         by=dx);
y <- seq(problem.parameters$ay,
         problem.parameters$by,
         by=dy);

for (k in seq(1,K)) {
    k.y <- ceiling(k/K.prime);
    k.x <- k - (K.prime * (k.y-1));
    print(c(k.x,k.y));
    function.params <- c(k.x, k.y);
    
    function.list[[k]] <-
        basis.function.xy(x,y,
                          function.params,
                          problem.parameters);
}

## problem.parameters$K.prime <- 30;
## for (n in seq(1,problem.parameters$K.prime)) {
##     k.x <- sample(seq(1,problem.parameters$K.prime), size=1);
##     k.y <- sample(seq(1,problem.parameters$K.prime), size=1);
##     function.params <- c(k.x, k.y);
    
##     function.list[[K+1]] <-
##         basis.function.xy(x,y,
##                           function.params,
##                           problem.parameters);
##     K <- K+1;
## }

par(mfrow=c(ceiling(sqrt(K)),
            ceiling(sqrt(K))));
par(mar = c(5,4,2,1));
if (PLOT.SOLUTION) {
    for (k in seq(1,K)) {
        contour(x,y,function.list[[k]]);
    }
}

orthonormal.function.list <- orthonormal.functions(function.list,
                                                   dx,dy,x,y,
                                                   PLOT.SOLUTION);

system.mats <- system.matrices(orthonormal.function.list,
                               dx,dy);

l2 <- blackbox(function.list,
               orthonormal.function.list,
               system.mats,
               problem.parameters,
               dx,dy,TRUE,TRUE);
print(l2);

data <- sample.process(10, dt, problem.parameters.generate.data);
for (n in seq(1,length(data))) {
    data[[n]]$number.terms = 1000;
}

l2 = blackbox(function.list,
              orthonormal.function.list,
              data[[1]],
              dx,dy,TRUE,TRUE);
