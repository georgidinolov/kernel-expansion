rm(list=ls());
source("1-d-solution.R");

problem.parameters = NULL;
problem.parameters$a = -1;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.2;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2 = 1;
problem.parameters$t = 0.3;

log.sigma2.vector = log(rep(0.13,8));
mu.vector = seq(problem.parameters$a,
                problem.parameters$b,
                length.out=length(log.sigma2.vector));
log.sigma2.mu.vector = c(log.sigma2.vector,
                         mu.vector);
dx = 0.001;
bb = blackbox(log.sigma2.mu.vector, problem.parameters, dx,TRUE);
print(bb);

opt.bases = optim(par=log.sigma2.mu.vector, fn=blackbox,
                  problem.parameters = problem.parameters,
                  dx = dx,
                  PLOT.SOLUTION=FALSE);

log.sigma2.mu.vector = opt.bases$par;
bb = blackbox(log.sigma2.mu.vector, problem.parameters, dx,TRUE);

dev.off();
K = length(log.sigma2.mu.vector)/2;
x = seq(10*problem.parameters$a,
        10*problem.parameters$b,
        length.out=1000);
for (k in seq(1,K)) {
    if (k==1) {
        plot(x, dnorm(x,mean=log.sigma2.mu.vector[k+K],
                      sd=sqrt(exp(log.sigma2.mu.vector[k]))),
             type="l")
    } else {
        lines(x, dnorm(x,mean=log.sigma2.mu.vector[k+K],
                      sd=sqrt(exp(log.sigma2.mu.vector[k]))));
    }
}
abline(v=c(problem.parameters$a,problem.parameters$b),
       col="red");
