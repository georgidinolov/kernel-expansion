rm(list=ls());

PLOT.SOLUTION=TRUE;
dx = 0.01;
dy = 0.01;
K.prime = 4;

library("mvtnorm");
source("2-d-solution.R");
problem.parameters = NULL;
problem.parameters$ax = 0;
problem.parameters$bx = 1;
problem.parameters$ay = 0;
problem.parameters$by = 1;
problem.parameters$x.ic = 0.9;
problem.parameters$y.ic = 0.9;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2.x = 1e-1;
problem.parameters$sigma.2.y = 1e-1;
problem.parameters$rho = 0.0;
problem.parameters$t = 0.1;
problem.parameters$K.prime = K.prime;

K <- K.prime^2;
function.list <- vector("list", K.prime^2);

x <- seq(problem.parameters$ax,
         problem.parameters$bx,
         by=dx);
y <- seq(problem.parameters$ay,
         problem.parameters$by,
         by=dy);

par(mfrow=c(ceiling(sqrt(K)),
            ceiling(sqrt(K))));
    par(mar = c(5,4,2,1));
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

problem.parameters$K.prime <- 30;
for (n in seq(1,problem.parameters$K.prime)) {
    k.x <- sample(seq(1,problem.parameters$K.prime), size=1);
    k.y <- sample(seq(1,problem.parameters$K.prime), size=1);
    function.params <- c(k.x, k.y);
    
    function.list[[K+1]] <-
        basis.function.xy(x,y,
                          function.params,
                          problem.parameters);
    K <- K+1;
}

par(mfrow=c(ceiling(sqrt(K)),
            ceiling(sqrt(K))));
par(mar = c(5,4,2,1));
if (PLOT.SOLUTION) {
    for (k in seq(1,K)) {
        contour(x,y,function.list[[k]]);
    }
}

problem.parameters$x.ic = 0.1;
problem.parameters$y.ic = 0.1;

l2 = blackbox(function.list,
              problem.parameters,
              dx,dy,TRUE,TRUE);
print(l2);
