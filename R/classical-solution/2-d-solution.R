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
    ss.1 <- atan(-sqrt(1+rho)/sqrt(1-rho));
    C.1 <- (-eta.ic + xi.ic*sqrt(1-rho)/sqrt(1+rho))/
        (sin(ss.1) - cos(ss.1)*sqrt(1-rho)/sqrt(1+rho));
    lines(c(xi.ic, C.1*cos(ss.1)+xi.ic),
          c(eta.ic, C.1*sin(ss.1)+eta.ic),
          col="red");

    ## border 2
    ss.2 <- pi + ss.1;
    C.2 <- (-eta.ic + xi.ic*sqrt(1-rho)/sqrt(1+rho) +
            1/(sigma.y*sqrt(1+rho)*cc))/
        (sin(ss.2) - cos(ss.2)*sqrt(1-rho)/sqrt(1+rho));
    lines(c(xi.ic, C.2*cos(ss.2)+xi.ic),
          c(eta.ic, C.2*sin(ss.2)+eta.ic),
          col="red");

    ## border 4
    ss.4 <- atan(sqrt(1+rho)/sqrt(1-rho));
    C.4 <- (-eta.ic - xi.ic*sqrt(1-rho)/sqrt(1+rho) +
            1/(sigma.x*sqrt(1+rho)*cc))/
        (sin(ss.4) + cos(ss.4)*sqrt(1-rho)/sqrt(1+rho));
    lines(c(xi.ic, C.4*cos(ss.4)+xi.ic),
          c(eta.ic, C.4*sin(ss.4)+eta.ic),
          col="red");
    
    ## border 3
    ss.3 <- pi + ss.4;
    C.3 <- (-eta.ic - xi.ic*sqrt(1-rho)/sqrt(1+rho))/
        (sin(ss.3) + cos(ss.3)*sqrt(1-rho)/sqrt(1+rho));
    lines(c(xi.ic, C.3*cos(ss.3)+xi.ic),
          c(eta.ic, C.3*sin(ss.3)+eta.ic),
          col="red");
    
    Cs <- c(C.1,C.2,C.3,C.4);
    sorted.ts <- sort.int(Cs,index.return=TRUE);

    ## ic radius and closest corner
    r.not <- min(sqrt(apply((boundary.points.mat - c(xi.ic, eta.ic))^2, 2, sum)));
    closest.boundary.point <-
        which(sqrt(apply((boundary.points.mat-
                          c(xi.ic, eta.ic))^2, 2, sum)) ==
              r.not);

    ## the offset phi0 and phi##
    if (closest.boundary.point == 4) {
        phi.not <- atan(slopes[2]) + pi;
        phi <- (atan(slopes[4]) + 2*pi) -
            (atan(slopes[2]) + pi);
        r.star <- r.not + 
    } else if (closest.boundary.point == 3) {
        phi.not <- atan(slopes[4]) + pi;
        phi <- (atan(slopes[1]) + pi) -
            (atan(slopes[4]) + pi);
    } else if (closest.boundary.point == 2) {
        phi.not <- atan(slopes[3]);
        phi <- atan(slopes[2]) - atan(slopes[3]);
    } else if (closest.boundary.point == 1) {
        phi.not <- atan(slopes[1]);
        phi <- (atan(slopes[3]) + pi) -
            atan(slopes[1]);
    } else {
        ## THROW
    }
    
    ## closest border indeces
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
            solution <- sol(r.star, tt);
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
