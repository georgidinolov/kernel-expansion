rm(list=ls());
source("1-d-solution.R");

problem.parameters = NULL;
problem.parameters$a = -1;
problem.parameters$b = 1;
problem.parameters$x.ic = 0.2;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2 = 1;
problem.parameters$t = 0.1;

Ks = seq(4,15);

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
    dx = 0.0001;
    bb = blackbox(log.sigma2.mu.vector, problem.parameters, dx,TRUE,TRUE);
    L2.remainders.before[i] = bb;
    L2.diff.before[i] = blackbox(log.sigma2.mu.vector, problem.parameters,
                                 dx,
                                 FALSE,
                                 FALSE);

    ## Start the clock!
    ptm <- proc.time();
    opt.bases = optim(par=log.sigma2.mu.vector, fn=blackbox,
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
