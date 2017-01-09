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
    
    ## border 1
    ss.1 <- atan(-sqrt(1+rho)/sqrt(1-rho));
    C.1 <- (-eta.ic + xi.ic*sqrt(1-rho)/sqrt(1+rho))/
        (sin(ss.1) - cos(ss.1)*sqrt(1-rho)/sqrt(1+rho));


    ## border 2
    ss.2 <- pi + ss.1;
    C.2 <- (-eta.ic + xi.ic*sqrt(1-rho)/sqrt(1+rho) +
            1/(sigma.y*sqrt(1+rho)*cc))/
        (sin(ss.2) - cos(ss.2)*sqrt(1-rho)/sqrt(1+rho));


    ## border 4
    ss.4 <- atan(sqrt(1+rho)/sqrt(1-rho));
    C.4 <- (-eta.ic - xi.ic*sqrt(1-rho)/sqrt(1+rho) +
            1/(sigma.x*sqrt(1+rho)*cc))/
        (sin(ss.4) + cos(ss.4)*sqrt(1-rho)/sqrt(1+rho));

    ## border 3
    ss.3 <- pi + ss.4;
    C.3 <- (-eta.ic - xi.ic*sqrt(1-rho)/sqrt(1+rho))/
        (sin(ss.3) + cos(ss.3)*sqrt(1-rho)/sqrt(1+rho));

    
    Cs <- c(C.1,C.2,C.3,C.4);
    ss.s <- c(ss.1, ss.2, ss.3, ss.4);

    ## ## PLOTTING ####
    ## par(mfrow=c(1,2));
    ## plot(xis, xis*sqrt(1-rho)/sqrt(1+rho), type="l",
    ##      xlim=c(min(xis.lim[1], etas.lim[1]),
    ##             max(xis.lim[2], etas.lim[2])),
    ##      ylim=c(min(xis.lim[1], etas.lim[1]),
    ##             max(xis.lim[2], etas.lim[2])));
    ## lines(xis, (1 + xis*sigma.y*sqrt(1-rho)*cc)/(sigma.y*sqrt(1+rho)*cc));
    ## lines(xis, -xis*sqrt(1-rho)/sqrt(1+rho));
    ## lines(xis, -xis*sqrt(1-rho)/sqrt(1+rho) + 1/(sigma.x*sqrt(1+rho)*cc));

    ## points((T.mat %*% c(0,0))[1],
    ## (T.mat %*% c(0,0))[2], col = "red");
    ## points((T.mat %*% c(0,1))[1],
    ## (T.mat %*% c(0,1))[2], col = "red");
    ## points((T.mat %*% c(1,0))[1],
    ## (T.mat %*% c(1,0))[2], col = "red");
    ## points((T.mat %*% c(1,1))[1],
    ## (T.mat %*% c(1,1))[2], col = "red");
    ## points(xi.ic, eta.ic, col = "green");

    ##     lines(c(xi.ic, C.1*cos(ss.1)+xi.ic),
    ##       c(eta.ic, C.1*sin(ss.1)+eta.ic),
    ##       col="red");

    ##     lines(c(xi.ic, C.2*cos(ss.2)+xi.ic),
    ##       c(eta.ic, C.2*sin(ss.2)+eta.ic),
    ##       col="red");
    ##     lines(c(xi.ic, C.4*cos(ss.4)+xi.ic),
    ##       c(eta.ic, C.4*sin(ss.4)+eta.ic),
    ##       col="red");
    ##     lines(c(xi.ic, C.3*cos(ss.3)+xi.ic),
    ##       c(eta.ic, C.3*sin(ss.3)+eta.ic),
    ##       col="red");
    ## ## PLOTTING END ##

    ##
    sorted.Cs <- sort.int(Cs, index.return=TRUE);

    ## ## ic radius and closest corner
    ## r.not <- min(sqrt(apply((boundary.points.mat - c(xi.ic, eta.ic))^2, 2, sum)));
    ## closest.boundary.point <-
    ##     which(sqrt(apply((boundary.points.mat-
    ##                       c(xi.ic, eta.ic))^2, 2, sum)) ==
    ##           r.not);

    ## ## the offset phi0 and phi##
    ## if (closest.boundary.point == 4) {
    ##     phi <- (atan(slopes[4]) + 2*pi) -
    ##         (atan(slopes[2]) + pi);
    ##     phi.not <- atan(slopes[2]) + pi;

    ##     if (phi > pi/2) {
    ##         C.1 <- (-boundary.points.mat[2,closest.boundary.point] + boundary.points.mat[1,closest.boundary.point]*sqrt(1-rho)/sqrt(1+rho))/
    ##             (sin(ss.1) - cos(ss.1)*sqrt(1-rho)/sqrt(1+rho));
            
    ##         C.3 <- (-boundary.points.mat[2,closest.boundary.point] - boundary.points.mat[1,closest.boundary.point]*sqrt(1-rho)/sqrt(1+rho))/
    ##             (sin(ss.3) + cos(ss.3)*sqrt(1-rho)/sqrt(1+rho));
            
    ##         r.star <- min(c(C.1, C.3));
            

    ##         if (Cs[2] < Cs[4]) {
    ##             theta.not <- asin(Cs[2]/r.not);
    ##         } else {
    ##             theta.not <- (phi - asin(Cs[4]/r.not));
    ##         }
    ##     } else {
    ##         r.star <- min(c(sqrt(sum((boundary.points.mat[,closest.boundary.point]-
    ##                                   boundary.points.mat[,2])^2)),
    ##                         sqrt(sum((boundary.points.mat[,closest.boundary.point]-
    ##                                   boundary.points.mat[,3])^2))));
    ##         theta.not <- asin(Cs[2]/r.not);
    ##     }
        
    ## } else if (closest.boundary.point == 3) {
    ##     phi <- (atan(slopes[1]) + pi) -
    ##         (atan(slopes[4]) + pi);
    ##     phi.not <- atan(slopes[4]) + pi;

    ##     if (phi > pi/2) {
    ##         C.2 <- (-boundary.points.mat[2,closest.boundary.point] + boundary.points.mat[1,closest.boundary.point]*sqrt(1-rho)/sqrt(1+rho) +
    ##                 1/(sigma.y*sqrt(1+rho)*cc))/
    ##             (sin(ss.2) - cos(ss.2)*sqrt(1-rho)/sqrt(1+rho));
            
    ##         C.3 <- (-boundary.points.mat[2,closest.boundary.point] - boundary.points.mat[1,closest.boundary.point]*sqrt(1-rho)/sqrt(1+rho))/
    ##             (sin(ss.3) + cos(ss.3)*sqrt(1-rho)/sqrt(1+rho));
            
    ##         r.star <- min(C.2, C.3);
            
    ##         if (Cs[1] < Cs[4]) {
    ##             theta.not <- phi - asin(Cs[1]/r.not);
    ##         } else {
    ##             theta.not <- asin(Cs[4]/r.not);
    ##         }
    ##     } else {
    ##         r.star <- min(c(sqrt(sum((boundary.points.mat[,closest.boundary.point]-
    ##                                   boundary.points.mat[,1])^2)),
    ##                         sqrt(sum((boundary.points.mat[,closest.boundary.point]-
    ##                                   boundary.points.mat[,4])^2))));
    ##         theta.not <- asin(Cs[4]/r.not);
    ##     }
        
    ## } else if (closest.boundary.point == 2) {
    ##     phi.not <- atan(slopes[3]) + 2*pi;
    ##     phi <- atan(slopes[2]) - atan(slopes[3]);

    ##     if (phi > pi/2) {
    ##         C.1 <- (-boundary.points.mat[2,closest.boundary.point] + boundary.points.mat[1,closest.boundary.point]*sqrt(1-rho)/sqrt(1+rho))/
    ##             (sin(ss.1) - cos(ss.1)*sqrt(1-rho)/sqrt(1+rho));
    ##         C.4 <- (-boundary.points.mat[2,closest.boundary.point] - boundary.points.mat[1,closest.boundary.point]*sqrt(1-rho)/sqrt(1+rho) +
    ##                 1/(sigma.x*sqrt(1+rho)*cc))/
    ##             (sin(ss.4) + cos(ss.4)*sqrt(1-rho)/sqrt(1+rho));
    ##         r.star <- min(C.1,C.4);

    ##         if (Cs[2] < Cs[3]) {
    ##             theta.not <- phi - asin(Cs[2]/r.not);
    ##         } else {
    ##             theta.not <- asin(Cs[3]/r.not);
    ##         }
    ##     } else {
    ##         r.star <- 
    ##             min(c(sqrt(sum((boundary.points.mat[,closest.boundary.point]-
    ##                             boundary.points.mat[,1])^2)),
    ##                   sqrt(sum((boundary.points.mat[,closest.boundary.point]-
    ##                             boundary.points.mat[,4])^2))));
    ##         theta.not <- asin(Cs[3]/r.not);
    ##     }
    ## } else if (closest.boundary.point == 1) {
    ##     phi.not <- atan(slopes[1]);
    ##     phi <- (atan(slopes[3]) + pi) -
    ##         atan(slopes[1]);

    ##     if (phi/2 > pi) {
    ##         C.2 <- (-boundary.points.mat[2,closest.boundary.point] + boundary.points.mat[1,closest.boundary.point]*sqrt(1-rho)/sqrt(1+rho) +
    ##                 1/(sigma.y*sqrt(1+rho)*cc))/
    ##             (sin(ss.2) - cos(ss.2)*sqrt(1-rho)/sqrt(1+rho));
            
    ##         C.4 <- (-boundary.points.mat[2,closest.boundary.point] - boundary.points.mat[1,closest.boundary.point]*sqrt(1-rho)/sqrt(1+rho) +
    ##                 1/(sigma.x*sqrt(1+rho)*cc))/
    ##             (sin(ss.4) + cos(ss.4)*sqrt(1-rho)/sqrt(1+rho));
            
    ##         r.star <- min(C.2,C.4);
            
    ##         if (Cs[1] < Cs[3]) {
    ##             theta.not <- asin(Cs[1]/r.not);
    ##         } else {
    ##             theta.not <- (phi - asin(Cs[3]/r.not));
    ##         }
    ##     } else {
    ##         r.star <- 
    ##             min(c(sqrt(sum((boundary.points.mat[,closest.boundary.point]-
    ##                             boundary.points.mat[,2])^2)),
    ##                   sqrt(sum((boundary.points.mat[,closest.boundary.point]-
    ##                             boundary.points.mat[,3])^2))));
    ##         theta.not <- asin(Cs[1]/r.not);
    ##     }
    ## } else {
    ##     ## THROW ERROR ##
    ## }

    ## lines(c(boundary.points.mat[1,closest.boundary.point],
    ##         r.star*cos(phi.not) + boundary.points.mat[1,closest.boundary.point]),
    ##       c(boundary.points.mat[2,closest.boundary.point],
    ##         r.star*sin(phi.not) + boundary.points.mat[2,closest.boundary.point]),
    ##       col="red",
    ##       lwd=2)
    ## lines(c(boundary.points.mat[1,closest.boundary.point],
    ##         r.star*cos(phi.not+phi) + boundary.points.mat[1,closest.boundary.point]),
    ##       c(boundary.points.mat[2,closest.boundary.point],
    ##         r.star*sin(phi.not+phi) + boundary.points.mat[2,closest.boundary.point]),
    ##       col="red",
    ##       lwd=2)
    ## lines(c(boundary.points.mat[1,closest.boundary.point],
    ##         r.star*cos(phi.not+theta.not) +
    ##         boundary.points.mat[1,closest.boundary.point]),
    ##       c(boundary.points.mat[2,closest.boundary.point],
    ##         r.star*sin(phi.not+theta.not) +
    ##         boundary.points.mat[2,closest.boundary.point]),
    ##       col="red",
    ##       lwd=2);
            
    ## K <- 40;
    ## L <- 2^8;
    ## Ckls <- rep(NA, K*L);
    ## lambdas <- rep(NA, K*L);
    ## rr <- seq(0, r.star, by=0.01);
    ## solution <- rep(0, length(rr));
    
    ## for (k in seq(1,K)) {
    ##     alpha <- k*pi/phi;
    ##     lambda.kls <- bessel_zero_Jnu(nu=alpha, s=seq(1,L));
    ##     lambdas[seq((k-1)*L+1,k*L)] <- lambda.kls;
        
    ##     Ckls[seq((k-1)*L+1,k*L)] = sin(k*pi*theta.not/phi) *
    ##         bessel_Jnu(alpha, r.not*lambda.kls/r.star)/
    ##         (phi/2*(r.star^2)/2*bessel_Jnu(alpha+1,lambda.kls)^2);
    ## }
    
    ## alphas <- unlist(lapply(seq(1,K),
    ##                         function(x){return(rep(x*pi/phi,L))}));
    ## tt <- r.star/200;
    ## time.bases <- exp(-(lambdas/r.star)^2*tt);

    ## generate.trig.bases <- function(theta) {
    ##     return (sin(rep(seq(1,K), each=L)*pi*theta/phi));
    ## }
    
    ## sol <- function(r, theta, time.bases) {
    ##     trig.bases <- generate.trig.bases(theta);
    ##     return(sum(trig.bases *
    ##                Ckls *
    ##                time.bases *
    ##                besselJ(nu = alphas,
    ##                        x = r*lambdas/r.star))*r/
    ##            (sigma.x * sigma.y *sqrt(1-rho)*sqrt(1+rho)));
    ## }
    
    ## sol.vector <- function(rs, thetas, time.bases) {
    ##     return(sapply(X=seq(1,length(rs)),
    ##            function(x) {
    ##                return(sol(rs[x],thetas[x],time.bases));
    ##            }));
    ## }

    xx <- seq(problem.parameters$ax,
              problem.parameters$bx,
              by=dx);
    yy <- seq(problem.parameters$ay,
              problem.parameters$by,
              by=dy);
    
    ## rs <- seq(r.star*0.001,r.star,length.out=100);
    ## plot(rs, sol.vector(rs,rep(theta.not,length(rs)),time.bases),
    ##      type="l")
    ## abline(v = r.not, col="blue",lwd=2);

    big.solution <- matrix(nrow = length(xx),
                           ncol = length(yy),
                           data = 0);
    tt <- (min(sorted.Cs$x[2])/4)^2;
    Max <- dmvnorm(x=c(xi.ic,eta.ic),mean=c(xi.ic,eta.ic),sigma=diag(rep(tt,2)))/
        (sigma.x*sigma.y*sqrt(1-rho)*sqrt(1+rho))
    xi.ic.reflected = 2*sorted.Cs$x[1]*cos(ss.s[sorted.Cs$ix[1]]) + xi.ic;
    eta.ic.reflected = 2*sorted.Cs$x[1]*sin(ss.s[sorted.Cs$ix[1]]) + eta.ic;
    
    for (ii in seq(1,length(xx))) {
            xieta <- T.mat %*% rbind(rep(xx[ii], length(yy)),
                                     yy);

            ## xieta.shifted <-
            ##     xieta - c(boundary.points.mat[1,closest.boundary.point],
            ##               boundary.points.mat[2,closest.boundary.point]);
            
            ## rs <- sqrt(apply( xieta.shifted^2, 2, sum));
            ## thetas <- atan(xieta.shifted[2,]/
            ##            xieta.shifted[1,]) - phi.not +
            ##       pi*((xieta.shifted[1,] < 0) &
            ##           (xieta.shifted[2,] < 0)) +
            ##       pi*((xieta.shifted[1,] < 0) &
            ##           (xieta.shifted[2,] > 0));
            ## indeces <- which(rs < r.star);
            ## print(c(ii, indeces));

            ## if (length(indeces) > 0) {
            ##     big.solution[ii,indeces] = sol.vector(rs[indeces],
            ##                                           thetas[indeces],
            ##                                           time.bases);
            ## }
            
            big.solution[ii,] = (dmvnorm(t(xieta),
                                         mean=c(xi.ic,eta.ic),
                                         sigma=diag(rep(tt,2))) -
                                 dmvnorm(t(xieta),
                                         mean=c(xi.ic.reflected,eta.ic.reflected),
                                         sigma=diag(rep(tt,2)))) /
                (sigma.x*sigma.y*sqrt(1-rho)*sqrt(1+rho));
            
            ##               max(xis.lim[2], etas.lim[2])));
            ## points(xieta[1,],xieta[2,], pch=20,
            ##        col=rgb(abs(solution.current/Max),0,0));
    }
    out = NULL;
    out$big.solution <- big.solution;
    out$tt <- tt;
    return(out);
    
    ## while (abs(solution[length(solution)]) > 1e-10 &
    ##        min(solution) >= 0) {
    ##            tt <- tt/2;
    ##            solution <- sapply(X=rr, function(x){return(sol(x,tt))});
    ##            plot(rr, solution, type = "l");
    ##            abline(v=r.star,col="red");
    ##            abline(v=r.not,col="blue",lwd=2);
    ##        }s
}
