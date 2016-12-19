rm(list=ls());

PLOT.SOLUTION=FALSE;
dx = 0.01;
dy = 0.01;
K=30;
nn = 10;
source("2-d-solution.R");
problem.parameters = NULL;
problem.parameters$ax = -1;
problem.parameters$bx = 1;
problem.parameters$ay = -1;
problem.parameters$by = 1;
problem.parameters$x.ic = 0.90;
problem.parameters$y.ic = 0.90;
problem.parameters$number.terms = 1000;
problem.parameters$sigma.2.x = 1e-1;
problem.parameters$sigma.2.y = 1e0;
problem.parameters$rho = 0.0;
problem.parameters$t = 0.5;
sigma2 = .05;


opt.results <- optim(log(sigma2),
                     blackbox,
                     method="Brent",
                     K=K,
                     problem.parameters=problem.parameters,
                     dx=dx,
                     dy=dy,
                     PLOT.SOLUTION=TRUE,
                     MINIMIZE.REMAINDER=TRUE,
                     lower = -4,
                     upper = -1,
                     control=list(maxit=100,
                                  reltol=5e-2));
sigma2 <- exp(opt.results$par);

problem.parameters$rho = -0.9;
samples <- sample.process(n.simulations=nn*10,
                          n.samples=nn,
                          dt=problem.parameters$t/100,
                          problem.parameters=problem.parameters);

bb <- blackbox(log(sigma2),
               K=K,
               n.samples = 40,
               problem.parameters,
               dx, dy,
               TRUE,
               TRUE);

Ks <- seq(10,20,by=1);
optimal.sigma2s <- rep(NA,length(Ks));

for (i in seq(1,length(Ks))) {
    K = Ks[i];
    sigma2 = 0.4;
    opt.results <- optim(log(sigma2),
                         blackbox,
                         method="Brent",
                         K=K,
                         problem.parameters=problem.parameters,
                         dx=dx,
                         dy=dy,
                         PLOT.SOLUTION=TRUE,
                         MINIMIZE.REMAINDER=TRUE,
                         lower = -4,
                         upper = 2,
                         control=list(maxit=100,
                                      reltol=5e-2));
    
    optimal.sigma2s[i] <- exp(opt.results$par);
}

fit <- lm(log(optimal.sigma2s)~Ks);
plot(Ks, log(optimal.sigma2s));
abline(a=fit$coefficients[1], b=fit$coefficients[2],
       col="red");

K=40;
bb = blackbox((fit$coefficients[1]+
               fit$coefficients[2]*K),
              K=K,
              problem.parameters,
              dx, dy,
              TRUE,TRUE);

sigma2 <- exp(fit$coefficients[1]+
              fit$coefficients[2]*K);
opt.results <- optim(log(sigma2),
                     blackbox,
                     method="Brent",
                     K=K,
                     problem.parameters=problem.parameters,
                     dx=dx,
                     dy=dy,
                     PLOT.SOLUTION=TRUE,
                     MINIMIZE.REMAINDER=TRUE,
                     lower = log(sigma2)-2*summary(fit)$sigma,
                     upper = log(sigma2)+2*summary(fit)$sigma,
                     control=list(maxit=100,
                                  reltol=5e-2));

problem.parameters$x.ic = 0.65;
problem.parameters$y.ic = 0.05;

bb = blackbox(opt.results$par,
              K=K,
              problem.parameters,
              dx, dy,
              TRUE,TRUE);
