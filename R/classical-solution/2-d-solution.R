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
    
    xis.lim <- c(min(boundary.points.mat[1,]),
                 max(boundary.points.mat[1,]));
                     
    etas.lim <- c(min(boundary.points.mat[2,]),
                  max(boundary.points.mat[2,]));
    
    xis <- seq(xis.lim[1], xis.lim[2], by=0.01);
    etas <- seq(etas.lim[1], etas.lim[2], by=0.01);
    
    plot(xis, xis*sqrt(1-rho)/sqrt(1+rho), type="l",
         xlim=xis.lim,
         ylim=etas.lim);
    lines(xis, (1 + xis*sigma.y*sqrt(1-rho)*cc)/(sigma.y*sqrt(1+rho)*cc));
    lines(xis, -xis*sqrt(1-rho)/sqrt(1+rho));
    lines(xis, -xis*sqrt(1-rho)/sqrt(1+rho) + 1/(sigma.x*sqrt(1+rho)*cc));

    slopes <- c(sqrt(1-rho)/sqrt(1+rho),
                sqrt(1-rho)/sqrt(1+rho),
                -sqrt(1-rho)/sqrt(1+rho),
                -sqrt(1-rho)/sqrt(1+rho));

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

    r.star <- 1;
    alpha <- 3*pi/phi;
    lam <- bessel_zero_Jnu(nu=alpha, s=seq(1,5));
    rr <- seq(0,2*r.star,by=0.01);
    for (j in seq(1,length(lam))) {
        if (j==1) {
            plot(rr, bessel_Jnu(alpha, x=rr/r.star*lam[j]), type="l");
        } else {
            lines(rr, bessel_Jnu(alpha, x=rr/r.star*lam[j]));
        }
    }
    abline(v=r.star,col="red");
    abline(v=min(sqrt(apply((boundary.points.mat - c(xi.ic, eta.ic))^2, 2, sum))),
           col="blue",
           lwd=2)
}
