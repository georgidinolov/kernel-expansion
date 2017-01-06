## problem.parameters
bivariate.solution.classical <- function(dx, dy, problem.parameters) {
    sigma.x <- sqrt(problem.parameters$sigma.2.x);
    sigma.y <- sqrt(problem.parameters$sigma.2.y);
    rho <- problem.parameters$rho;
    cc <- sin(pi/4);

    T.mat <- matrix(nrow=2,ncol=2,
                    data=c(c(cc/(sigma.x*sqrt(1-rho)),
                             cc/(sigma.x*sqrt(1+rho))),
                           c(-cc/(sigma.y*sqrt(1-rho)),
                             cc/(sigma.y*sqrt(1+rho)))),
                    byrow=FALSE);
    T.mat.inv <- solve(T.mat);

    xi.ic <- (T.mat %*% c(problem.parameters$x.ic, problem.parameters$y.ic))[1];
    eta.ic <- (T.mat %*% c(problem.parameters$x.ic, problem.parameters$y.ic))[2];

    boundary.points.mat <- T.mat %*%
        matrix(nrow=2,ncol=4,
               data=c(c(0,0),
                      c(0,1),
                      c(1,0),
                      c(1,1)));

    slopes <- c(sqrt(1-rho)/sqrt(1+rho),
                sqrt(1-rho)/sqrt(1+rho),
                -sqrt(1-rho)/sqrt(1+rho),
                -sqrt(1-rho)/sqrt(1+rho));
    
    xis.lim <- c(min(boundary.points.mat[1,]),
                 max(boundary.points.mat[1,]));
                     
    etas.lim <- c(min(boundary.points.mat[2,]),
                  max(boundary.points.mat[2,]));
    
    xis <- seq(xis.lim[1], xis.lim[2], by=0.01);
    etas <- seq(etas.lim[1], etas.lim[2], by=0.01);
    
    plot(xis, xis*sqrt(1-rho)/sqrt(1+rho), type="l",
         xlim=c(min(xis.lim[1], etas.lim[1]),
                max(xis.lim[2], etas.lim[2])),
         ylim=c(min(xis.lim[1], etas.lim[1]),
                max(xis.lim[2], etas.lim[2])));
    lines(xis, (1 + xis*sigma.y*sqrt(1-rho)*cc)/(sigma.y*sqrt(1+rho)*cc));
    lines(xis, -xis*sqrt(1-rho)/sqrt(1+rho));
    lines(xis, -xis*sqrt(1-rho)/sqrt(1+rho) + 1/(sigma.x*sqrt(1+rho)*cc));

    points((T.mat %*% c(0,0))[1],
    (T.mat %*% c(0,0))[2], col = "red");
    points((T.mat %*% c(0,1))[1],
    (T.mat %*% c(0,1))[2], col = "red");
    points((T.mat %*% c(1,0))[1],
    (T.mat %*% c(1,0))[2], col = "red");
    points((T.mat %*% c(1,1))[1],
    (T.mat %*% c(1,1))[2], col = "red");
    points(xi.ic, eta.ic, col = "green");
    
    ## border 1
    ss.1 <- atan(-sqrt(1-rho)/sqrt(1+rho));
    C.1 <- (-eta.ic + xi.ic*sqrt(1-rho)/sqrt(1+rho))/(2*sin(ss.1));

    ## border 2
    ss.2 <- atan(-sqrt(1-rho)/sqrt(1+rho));
    C.2 <- (-eta.ic + xi.ic*sqrt(1-rho)/sqrt(1+rho) +
            1/(sigma.y*sqrt(1+rho)*cc))/(2*sin(ss.2));

    ## border 3
    ss.3 <- atan(sqrt(1-rho)/sqrt(1+rho));
    C.3 <- (-eta.ic - xi.ic*sqrt(1-rho)/sqrt(1+rho))/(2*sin(ss.3));

    ## border 4
    ss.4 <- atan(sqrt(1-rho)/sqrt(1+rho));
    C.4 <- (-eta.ic - xi.ic*sqrt(1-rho)/sqrt(1+rho) +
            1/(sigma.x*sqrt(1+rho)*cc))/(2*sin(ss.4));
    
    Cs <- c(C.1,C.2,C.3,C.4);
    Cs <- ifelse(Cs < 0, -Cs, Cs);
    
    nn <- 3;
    ts <- Cs/nn;
    sorted.ts <- sort.int(ts,index.return=TRUE);
    tt <- sort(ts)[3];

    ## closest border indeces
    indeces <- sorted.ts$ix;

    phis <- c(atan(slopes[indeces[1]]),
              atan(slopes[indeces[2]]));
    phi <- max(phis)-min(phis);

    r.not <- min(sqrt(apply((boundary.points.mat - c(xi.ic, eta.ic))^2, 2, sum)));

    if (phi>0 & r.not > 0) {
        closest.boundary.point <-
            which(sqrt(apply((boundary.points.mat-
                              c(xi.ic, eta.ic))^2, 2, sum)) ==
                  r.not);
        
        r.star <- sort(sqrt(apply((boundary.points.mat -
                                   boundary.points.mat[,closest.boundary.point])^2, 2, sum)))[2];

        print(c(Cs[indeces[1]], r.not));
        if (abs(Cs[indeces[1]] - r.not) > 1e-16) {
            theta.not <- asin(x=Cs[indeces[1]]/r.not);
            
            K <- 50;
            L <- 2^5;
            Ckls <- rep(NA, K*L);
            lambdas <- rep(NA, K*L);
            rr <- seq(0, r.star, by=0.01);
            solution <- rep(0, length(rr));
            
            for (k in seq(1,K)) {
                alpha <- k*pi/phi;
                lambda.kls <- bessel_zero_Jnu(nu=alpha, s=seq(1,L));
                lambdas[seq((k-1)*L+1,k*L)] <- lambda.kls;
                
                Ckls[seq((k-1)*L+1,k*L)] = sin(k*pi*theta.not/phi) *
                    bessel_Jnu(alpha, r.not*lambda.kls/r.star)/
                    (phi/2*(r.star^2)/2*bessel_Jnu(alpha+1,lambda.kls)^2);
            }
            
            trig.bases <- unlist(lapply(X=seq(1,K),
                                        function(x) {
                                            return(rep(sin(x*pi*theta.not/phi),L))
                                        } ));
            alphas <- unlist(lapply(seq(1,K),
                                    function(x){return(rep(x*pi/phi,L))}));
            
            sol <- function(r, t) {
                return(sum(trig.bases *
                           Ckls *
                           exp(-(lambdas/r.star)^2*t) *
                           bessel_Jnu(alphas,
                                      r*lambdas/r.star)));
            }
            
            tt <- r.star/200;
            solution <- sapply(X=rr, function(x){return(sol(x,tt))});
            plot(rr, solution, type = "l");
            abline(v=r.star,col="red");
            abline(v=r.not,col="blue",lwd=2);
            while (abs(solution[length(solution)]) > 1e-10 &
                   min(solution) >= 0) {
                       tt <- tt/2;
                       solution <- sapply(X=rr, function(x){return(sol(x,tt))});
                       plot(rr, solution, type = "l");
                       abline(v=r.star,col="red");
                       abline(v=r.not,col="blue",lwd=2);
                   }
            
        }
    }
}
