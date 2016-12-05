rm(list=ls());
source("2-d-solution.R");

problem.parameters = NULL;
problem.parameters$ax = -1;
problem.parameters$bx = 1;
problem.parameters$ay = -1;
problem.parameters$by = 1;
problem.parameters$x.ic = 0.2;
problem.parameters$y.ic = 0.0;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2.x = 0.1;
problem.parameters$sigma.2.y = 0.1;
problem.parameters$rho = 0.1;
problem.parameters$t = 0.5;

Ks = seq(4,4);

L2.remainders.before = rep(NA, length(Ks));
L2.diff.before = rep(NA, length(Ks));

L2.remainders.after = rep(NA, length(Ks));
L2.diff.after = rep(NA, length(Ks));

ave.function.call.time.vec = rep(NA,length(Ks));
log.sigma2.mu.vector.list = vector(mode="list",
                                   length=length(Ks));

for (i in seq(1,length(Ks))) {
    K=Ks[i];
    log.sigma2.vector=log(
        rep(c(((problem.parameters$bx-problem.parameters$ax)/K)^2,
        ((problem.parameters$by-problem.parameters$ay)/K)^2),
        K));
    mu.vector = seq(problem.parameters$ax,
                    problem.parameters$bx,
                    length.out=length(log.sigma2.vector)/2);

    mu.vector = unlist(lapply(X=mu.vector, function(x) {c(x,x)}));
    mu.vector = c(-1,-1, -1,1, 1,-1, 1,1);
    
    log.sigma2.mu.vector = c(log.sigma2.vector,
                             mu.vector);
    dx = 0.001;
    dy = 0.001;
    
    bb = blackbox(log.sigma2.mu.vector, problem.parameters, dx, dy,
                  TRUE,TRUE);
    ## L2.remainders.before[i] = bb;
    ## L2.diff.before[i] = blackbox(log.sigma2.mu.vector, problem.parameters,
    ##                              dx,
    ##                              FALSE,
    ##                              FALSE);

    ## ptm <- proc.time();				 
    ## opt.bases <- optim(par=log.sigma2.mu.vector,
    ##                    fn=blackbox,
    ##                    problem.parameters = problem.parameters,
    ##                    dx = dx,
    ##                    PLOT.SOLUTION=FALSE,
    ##                    MINIMIZE.REMAINDER=TRUE);

    ## ## opt.bases <- optim(par=log.sigma2.mu.vector,
    ## ## 			method=c("BFGS"),
    ## ##                    fn=blackbox,
    ## ##                    problem.parameters = problem.parameters,
    ## ##                    dx = dx,
    ## ##                    PLOT.SOLUTION=FALSE,
    ## ##                    MINIMIZE.REMAINDER=TRUE);
    ## end.time <- proc.time()-ptm;
    ## ## Stop the clock
    
    ## ave.function.call.time = end.time[3]/(opt.bases$counts[1]);
    ## ave.function.call.time.vec[i] = ave.function.call.time;

    ## log.sigma2.mu.vector = opt.bases$par;
    ## log.sigma2.mu.vector.list[[i]] = log.sigma2.mu.vector;
    ## bb = blackbox(log.sigma2.mu.vector, problem.parameters, dx,TRUE,TRUE);
    
    ## L2.remainders.after[i] = bb;
    ## L2.diff.after[i] = blackbox(log.sigma2.mu.vector, problem.parameters,
    ##                             dx,
    ##                             FALSE,
    ##                             FALSE)
}

## save(file="optimization-results.Rdata",
##      list=c("Ks","L2.remainders.before","L2.diff.before",
##             "L2.remainders.after","L2.diff.after",
##             "ave.function.call.time.vec",
##             "log.sigma2.mu.vector.list",
## 	    "dx"));

## pdf("optimization-results.pdf");
## par(mfrow=c(2,1));
## par(mar = c(5,4,2,1));
## plot(Ks, log(L2.remainders.before),
##      type="l",
##      ylim = c(min(log(L2.remainders.before),log(L2.remainders.after)),
##               max(log(L2.remainders.before),log(L2.remainders.after))),
##      ylab = "",
##      xlab = "K",
##      main = "log(L^2) norm of remainder term");
## lines(Ks, log(L2.remainders.after),col="red");

## plot(Ks, log(L2.diff.before),
##      type="l",
##      ylim = c(min(log(L2.diff.before),log(L2.diff.after)),
##               max(log(L2.diff.before),log(L2.diff.after))),
##      ylab = "",
##      xlab = "",
##      main = "log(L^2) norm of difference between approximate and true solutions");
## lines(Ks, log(L2.diff.after),col="red");
## dev.off();
