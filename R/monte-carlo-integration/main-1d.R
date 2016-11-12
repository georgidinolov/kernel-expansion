rm(list=ls());
source("1-d-solution.R");

problem.parameters = NULL;
problem.parameters$a = 0;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.2;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2 = 1;
problem.parameters$t = 0.1;

sigma2.vector = rep(0.13,4);
dx = 0.0001;
bb = blackbox(sigma2.vector, problem.parameters, dx);
print(bb);

alpha.beta <- select.alpha.beta(problem.parameters);
alpha <- alpha.beta$alpha;
beta <- alpha.beta$beta;

function.list = vector(mode="list", length=6);
function.list[[1]] = kernel;
function.list[[2]] = function(x) { return (x*kernel(x)) };
function.list[[3]] = function(x) { return (x*(1-x)*dnorm(x)) };
function.list[[4]] = function(x) { return (x*(1-x)*dnorm(x, mean=0.5,
                                                         sd = 0.4)) };
function.list[[5]] = function(x) { return (x^2*kernel(x)) };
function.list[[6]] = function(x) { return (x^3*kernel(x)) };
function.list[[7]] = function(x) { return (x^4*kernel(x)) };
function.list[[8]] = function(x) { return (x^5*kernel(x)) };

## gram schmidt ##
dx = 0.0001;
x = seq(problem.parameters$a,
        problem.parameters$b,
        by=dx);
K <- length(function.list);
norms <- rep(NA, K);
coefficients <- matrix(nrow=K, ncol=K);
orthonormal.function.list <- vector(mode="list",
                                    length=K);

for (k in seq(1,K)) {
    if (k==1) {
        ## only normalize
        norm = sqrt(sum(function.list[[k]](x)^2)*dx);
        coefficients[k,k] = 1;
        norms[k] = norm;
        orthonormal.function.list[[k]] =
            function.list[[k]](x)/norm;
    } else {
        for (l in seq(1,k-1)) {
            coefficients[k,l] = -sum(orthonormal.function.list[[l]]*
                                     function.list[[k]](x)*dx);
        }
        coefficients[k,k] = 1;

        orthonormal.function.list[[k]] =
            function.list[[k]](x);
        for (l in seq(1,k-1)) {
            orthonormal.function.list[[k]] =
                orthonormal.function.list[[k]] +
                coefficients[k,l]*orthonormal.function.list[[l]];
        }
        norm = sqrt(sum(orthonormal.function.list[[k]]^2*dx));
        norms[k] = norm;
        orthonormal.function.list[[k]] =
            orthonormal.function.list[[k]]/norm;
    }
}

### plotting the bases START ###
for (k in seq(1,K)) {
    if (k==1) {
        plot(x, orthonormal.function.list[[1]],
             type = "l",
             ylim = c(-max(orthonormal.function.list[[1]]),
                      max(orthonormal.function.list[[1]])));
                           
    } else {
        lines(x, orthonormal.function.list[[k]]);
    }
}
### plotting the bases END ###

### check orthogonality START ###
orthogonality.matrix <- matrix(nrow=K,
                               ncol=K);
for (n in seq(1,K)) {
    for (m in seq(1,K)) {
        orthogonality.matrix[n,m] =
            sum(orthonormal.function.list[[n]]*
                orthonormal.function.list[[m]]*dx);
    }
}
print(orthogonality.matrix);
### check orthogonality START ###

### check representation of IC START ###
IC <- rep(0, length(x));
ic.index = which(abs(x-problem.parameters$x.ic)<=dx/2)
for (i in seq(1,K)) {
    IC = IC +
        orthonormal.function.list[[i]][ic.index]*
        orthonormal.function.list[[i]]
}
plot(x,IC, type = "l");
### check representation of IC END ###

### SYSTEM MATRICES START ### 
stiff.mat <- matrix(nrow=K,ncol=K);
for (i in seq(1,K)) {
    current.basis.dx.i = (orthonormal.function.list[[i]][-1]-
        orthonormal.function.list[[i]][-length(x)])/dx

    for (j in seq(i,K)) {
        current.basis.dx.j = (orthonormal.function.list[[j]][-1]-
                              orthonormal.function.list[[j]][-length(x)])/dx

        stiff.matrix.entry = sum(current.basis.dx.i*
                                 current.basis.dx.j*dx);
        print(c(i,j,stiff.matrix.entry));
        stiff.mat[i,j]=1/2*problem.parameters$sigma.2*stiff.matrix.entry;
        stiff.mat[j,i]=1/2*problem.parameters$sigma.2*stiff.matrix.entry;
    }
}

mass.mat <- matrix(nrow=K,ncol=K);
for (i in seq(1,K)) {
    for (j in seq(i,K)) {
        mass.matrix.entry = sum(orthonormal.function.list[[i]]*
                                orthonormal.function.list[[j]]*
                                dx);
        print(c(i,j,mass.matrix.entry));
        mass.mat[i,j]=mass.matrix.entry;
        mass.mat[j,i]=mass.matrix.entry;
    }
}
### SYSTEM MATRICES END ###

### eigenvalues START ###
eig <- eigen(solve(mass.mat) %*% stiff.mat);
### eigenvalues END ###

### ICs START ###
b = rep(NA, K);
for (i in seq(1,K)) {
    b[i] = orthonormal.function.list[[i]][ic.index];
}
IC.vec = solve(mass.mat, b);

IC <- rep(0, length(x));
for (i in seq(1,K)) {
    IC = IC +
        IC.vec[i]*
        orthonormal.function.list[[i]]
}
plot(x,IC, type = "l");
### ICs END ###

### APPROX SOLUTION START ###
coefs = (eig$vectors) %*% diag(exp(-eig$values * problem.parameters$t)) %*% t(eig$vectors) %*% IC.vec;

plot(x,univariate.solution.approx(coefs),type="l");
lines(x,univariate.solution(x,problem.parameters), col="green",
      lty="dashed", lwd = 2);
### APPROX SOLUTION END ###

### APPROX APPROX SOLUTION START ###
TT = 5000;
dt = problem.parameters$t/TT;
A = solve(mass.mat) %*% stiff.mat;
expAdt.approx =
    diag(1,K) - A*dt;
eig.approx <- eigen(expAdt.approx);

coefs.approx = eig.approx$vectors %*%
    diag(eig.approx$values^TT) %*%
    t(eig.approx$vectors) %*% 
    IC.vec;

plot(x,univariate.solution.approx(coefs.approx),
      col="red", type="l");
lines(x,univariate.solution.approx(coefs),col="green",lty="dashed");
lines(x,univariate.solution(x,problem.parameters), lty=4);
### APPROX APPROX SOLUTION END ###

### CHECKING PDE START ###
A = solve(mass.mat) %*% stiff.mat;
plot(x,univariate.solution.approx.dt(coefs, A),type="l", col="red");
lines(x,0.5*(problem.parameters$sigma.2)*
        univariate.solution.approx.dx.dx(coefs),
      lty="dashed");

difference2 = (univariate.solution.approx.dt(coefs.approx, A)-
                0.5*(problem.parameters$sigma.2)*
                univariate.solution.approx.dx.dx(coefs.approx))^2;
norm.diff2 = sum(difference2*dx);
print(norm.diff2);
### CHECKING PDE END ###
