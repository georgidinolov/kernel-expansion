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

    xis <- seq(-2,10,by=0.01);
    etas <- seq(-2,10,by=0.01);

    plot(xis, xis*sqrt(1-rho)/sqrt(1+rho), ylim = c(-2,10),
         type="l");
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
    
    nn <- 3;
    ts <- Cs/nn;
    sorted.ts <- sort.int(ts,index.return=TRUE);
    tt <- sort(ts)[3];

    ## closets border indeces
    indeces <- sorted.ts$ix[c(1,2)];

    phis <- c(atan(slopes[indeces[1]]),
              atan(slopes[indeces[2]]));
    phi <- max(phis)-min(phis);

    
}
