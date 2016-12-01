rm(list=ls());
source("1-d-solution.R");

function.params = raw.function.list[[1]];
a = problem.parameters$a;
b = problem.parameters$b;

mu = function.params[1];
sigma2 = function.params[2]/2;

alpha = (a-mu)/sqrt(sigma2);
beta = (b-mu)/sqrt(sigma2);

PP = pnorm(beta) - pnorm(alpha);

Ls = rep(NA,5);
for (i in seq(1,5)) {
    if (i==1) {
        Ls[i] = 1;
    } else if (i==2) {
        Ls[i] = -(dnorm(beta)-dnorm(alpha))/
            PP;
    } else {
        Ls[i] = -(beta^(i-2)*dnorm(beta)-alpha^(i-2)*dnorm(alpha))/
            PP + (i-2)*Ls[i-2];        
    }
}

Ms = sapply(seq(0,4),
            function(k) {sum(choose(k,seq(0,k))*
                             sqrt(sigma2)^(seq(0,k))*
                             mu^(k-seq(0,k))*
                             Ls[seq(0,k)+1])});
Ms = Ms*PP;

plot(x,1/sqrt(4*pi^2*function.params[2]^2)*
       exp(-0.5/function.params[2]*(x-function.params[1])^2)^2,
     type="l");

lines(x, product.coefficient(function.params,function.params)*
         exp(-0.5/(function.params[2]/2)*
             (x-function.params[1])^2),
      col = "red",
      lty="dashed");

plot(x,
     exp(-0.5/(function.params[2]/2)*
         (x-function.params[1])^2),
     type="l")

lines(x, dnorm(x,function.params[1],
                 sqrt(function.params[2]))^2,
      col = "green");

I1 = sum(dnorm(x,function.params[1],
               sqrt(function.params[2]))^2)*dx;

I2 = sum(exp(-1/(2*sigma2)*(x-mu)^2))*dx;

abs(I1 - product.coefficient(function.params,
                             function.params)*
    I2) <= 1e-15;
Ms[1] * sqrt(2*pi*sigma2);

norm.raw.function(function.params, a,b);
sqrt(sum((x-a)^2*(b-x)^2*
         dnorm(x,function.params[1],
               sqrt(function.params[2]))^2)*dx);

problem.parameters = NULL;
problem.parameters$a = -1;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.2;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2 = 1;
problem.parameters$t = 0.1;

Ks = seq(4,10);

L2.remainders.before = rep(NA, length(Ks));
L2.diff.before = rep(NA, length(Ks));

L2.remainders.after = rep(NA, length(Ks));
L2.diff.after = rep(NA, length(Ks));

ave.function.call.time.vec = rep(NA,length(Ks));
log.sigma2.mu.vector.list = vector(mode="list",
                                   length=length(Ks));

for (i in seq(1,length(Ks))) {
    K = Ks[i];
    log.sigma2.vector = log(rep(((problem.parameters$b-problem.parameters$a)/K)^2,
                                K));
    mu.vector = seq(problem.parameters$a,
                    problem.parameters$b,
                    length.out=length(log.sigma2.vector));
    log.sigma2.mu.vector = c(log.sigma2.vector,
                             mu.vector);
    dx = 0.001;
    bb = blackbox(log.sigma2.mu.vector, problem.parameters, dx,TRUE,TRUE);
    L2.remainders.before[i] = bb;
    L2.diff.before[i] = blackbox(log.sigma2.mu.vector, problem.parameters,
                                 dx,
                                 FALSE,
                                 FALSE);

    ## Start the clock!
    ptm <- proc.time();
    opt.bases = optim(par=log.sigma2.mu.vector,
                      method=c("CG"),
                      fn=blackbox,
                      problem.parameters = problem.parameters,
                      dx = dx,
                      PLOT.SOLUTION=FALSE,
                      MINIMIZE.REMAINDER=TRUE);
    end.time <- proc.time() - ptm;
    ## Stop the clock
    ave.function.call.time = end.time[3]/opt.bases$counts[1];
    ave.function.call.time.vec[i] = ave.function.call.time;

    log.sigma2.mu.vector = opt.bases$par;
    log.sigma2.mu.vector.list[[i]] = log.sigma2.mu.vector;
    bb = blackbox(log.sigma2.mu.vector, problem.parameters, dx,TRUE,TRUE);
    
    L2.remainders.after[i] = bb;
    L2.diff.after[i] = blackbox(log.sigma2.mu.vector, problem.parameters,
                                dx,
                                FALSE,
                                FALSE)
}

save(file="optimization-results.Rdata",
     list=c("Ks","L2.remainders.before","L2.diff.before",
            "L2.remainders.after","L2.diff.after",
            "ave.function.call.time.vec",
            "log.sigma2.mu.vector.list",
	    "dx"));

pdf("optimization-results.pdf");
par(mfrow=c(2,1));
par(mar = c(5,4,2,1));
plot(Ks, log(L2.remainders.before),
     type="l",
     ylim = c(min(log(L2.remainders.before),log(L2.remainders.after)),
              max(log(L2.remainders.before),log(L2.remainders.after))),
     ylab = "",
     xlab = "K",
     main = "log(L^2) norm of remainder term");
lines(Ks, log(L2.remainders.after),col="red");

plot(Ks, log(L2.diff.before),
     type="l",
     ylim = c(min(log(L2.diff.before),log(L2.diff.after)),
              max(log(L2.diff.before),log(L2.diff.after))),
     ylab = "",
     xlab = "",
     main = "log(L^2) norm of difference between approximate and true solutions");
lines(Ks, log(L2.diff.after),col="red");
dev.off();
