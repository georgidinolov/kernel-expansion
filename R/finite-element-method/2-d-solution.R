load.data.from.csv <- function(data.set.file) {
    input.data <- read.csv(file= data.set.file, header = TRUE, sep = ",");
    n <- dim(input.data)[1];
    data <- vector("list", n);
    for (i in seq(1,n)) {
        data[[i]]$t <- input.data$t[i];
        data[[i]]$sigma.2.x <- (input.data$sigma_x[i])^2;
        data[[i]]$sigma.2.y <- (input.data$sigma_y[i])^2;
        data[[i]]$rho <- input.data$rho[i];

        data[[i]]$ax <- input.data$a[i];
        data[[i]]$x.ic <- input.data$x_0[i];
        data[[i]]$bx <- input.data$b[i];
        data[[i]]$x.fc <- input.data$x_T[i];

        data[[i]]$ay <- input.data$c[i];
        data[[i]]$y.ic <- input.data$y_0[i];
        data[[i]]$by <- input.data$d[i];
        data[[i]]$y.fc <- input.data$y_T[i];
    }
    return (data);
}


rescale.problem <- function(problem.parameters.original) {
    sigma.2.x <- problem.parameters.original$sigma.2.x;
    sigma.2.y <- problem.parameters.original$sigma.2.y;
    x.current <- problem.parameters.original$x.fc;
    y.current <- problem.parameters.original$y.fc;
    
    ax <- problem.parameters.original$ax;
    bx <- problem.parameters.original$bx;
    ay <- problem.parameters.original$ay;
    by <- problem.parameters.original$by;
    
    ## STEP 1 ##
    x.ic <- problem.parameters.original$x.ic - ax;
    x.current <- x.current - ax;
    bx <- bx - ax;
    ax <- ax - ax; 
    sigma.2.x <- sigma.2.x;

    y.ic <- problem.parameters.generate.data$y.ic - ay
    y.current <- y.current - ay
    by <- by - ay;
    ay <- ay - ay;
    sigma.2.y <- sigma.2.y;
    
    ## STEP 2 ##
    L.x <- bx-ax;
    x.ic <- x.ic / L.x;
    x.current <- x.current / L.x;
    ax <- ax / L.x;
    bx <- bx / L.x;
    sigma.2.x <- 1/ (L.x^2);
    
    L.y <- by-ay;
    y.ic <- y.ic / L.y;
    y.current <- y.current / L.y;
    ay <- ay / L.y;
    by <- by / L.y;
    sigma.2.y <- 1/ (L.y^2);
    
    out <- problem.parameters.original;
    out$x.ic <- x.ic;
    out$y.ic <- y.ic;
    out$x.fc <- x.current;
    out$y.fc <- y.current;
    out$ax <- ax;
    out$bx <- bx;
    out$ay <- ay;
    out$by <- by;
    out$sigma.2.x <- sigma.2.x;
    out$sigma.2.y <- sigma.2.y;
    return(out);
}

## sample from process
sample.process <- function(n.samples,
                           dt,
                           problem.parameters.generate.data) {

    n.timesteps <- problem.parameters.generate.data$t/dt;
    rho <- problem.parameters.generate.data$rho;
    sigma.2.x <- problem.parameters.generate.data$sigma.2.x;
    sigma.2.y <- problem.parameters.generate.data$sigma.2.y;

    out <- vector("list", n.samples);
    
    for (i in seq(1,n.samples)) {
        x.current = problem.parameters.generate.data$x.ic;
        y.current = problem.parameters.generate.data$y.ic;

        ax <- x.current;
        bx <- x.current;
        ay <- y.current;
        by <- y.current;
        for (j in seq(1,n.timesteps)) {
            
            dx <- sqrt(dt*sigma.2.x)*rnorm(1);
            dy <- sigma.2.y/sigma.2.x*rho*
                (dx) +
                sqrt((1-rho^2)*sigma.2.y*dt)*rnorm(1);
            x.new <- x.current + dx;
            y.new <- y.current + dy;

            ax <- min(ax, x.new);
            bx <- max(bx, x.new);

            ay <- min(ay, y.new);
            by <- max(by, y.new);

            x.current <- x.new;
            y.current <- y.new;
        }

        ## ## STEP 1 ##
        ## x.ic <- problem.parameters.generate.data$x.ic / sqrt(sigma.2.x);
        ## x.current <- x.current / sqrt(sigma.2.x);
        ## ax <- ax / sqrt(sigma.2.x);
        ## bx <- bx / sqrt(sigma.2.x);
        ## sigma.2.x <- 1;

        ## y.ic <- problem.parameters.generate.data$y.ic / sqrt(sigma.2.y);
        ## y.current <- y.current / sqrt(sigma.2.y);
        ## ay <- ay / sqrt(sigma.2.y);
        ## by <- by / sqrt(sigma.2.y);
        ## sigma.2.y <- 1;

        ## ## STEP 2 ##
        ## L.x <- bx-ax;
        ## x.ic <- x.ic / L.x;
        ## x.current <- x.current / L.x;
        ## ax <- ax / L.x;
        ## bx <- bx / L.x;
        ## sigma.2.x <- 1/ (L.x^2);

        ## L.y <- by-ay;
        ## y.ic <- y.ic / L.y;
        ## y.current <- y.current / L.y;
        ## ay <- ay / L.y;
        ## by <- by / L.y;
        ## sigma.2.y <- 1/ (L.y^2);

        ## ## STEP 3 ##
        ## x.ic <- x.ic - ax;
        ## x.current <- x.current - ax;
        ## bx <- bx - ax;
        ## ax <- ax - ax;

        ## y.ic <- y.ic - ay;
        ## y.current <- y.current - ay;
        ## by <- by - ay;
        ## ay <- ay - ay;
        
        out[[i]] <- problem.parameters.generate.data;
        out[[i]]$x.ic <- problem.parameters.generate.data$x.ic;
        out[[i]]$y.ic <- problem.parameters.generate.data$y.ic;
        out[[i]]$x.fc <- x.current;
        out[[i]]$y.fc <- y.current;
        out[[i]]$ax <- ax;
        out[[i]]$bx <- bx;
        out[[i]]$ay <- ay;
        out[[i]]$by <- by;
        out[[i]]$sigma.2.x <- sigma.2.x;
        out[[i]]$sigma.2.y <- sigma.2.y;
    }
    return (out);
}

sample.process.boundaries <- function(n.simulations,
                                      n.samples,
                                      dt,
                                      problem.parameters) {
    
    n.timesteps <- problem.parameters$t/dt;
    rho <- problem.parameters$rho;
    sigma.2.x <- problem.parameters$sigma.2.x;
    sigma.2.y <- problem.parameters$sigma.2.y;

    output.x <- c();
    output.y <- c();
    weights <- rep(1, n.simulations);
    
    for (i in seq(1,n.simulations)) {
        x.current = problem.parameters$x.ic;
        y.current = problem.parameters$y.ic;
        weight = 1;
        for (j in seq(1,n.timesteps)) {
            
            dx <- sqrt(dt*sigma.2.x)*rnorm(1);
            dy <- sigma.2.y/sigma.2.x*rho*
                (dx) +
                sqrt((1-rho^2)*sigma.2.y*dt)*rnorm(1);
            x.new <- x.current + dx;
            y.new <- y.current + dy;

            ## if (x.new > problem.parameters$bx |
            ##     x.new < problem.parameters$ax |
            ##     y.new > problem.parameters$by |
            ##     y.new < problem.parameters$ay) {
            ##     print("hit boundary");
            ##     print(j);
            ##     break;
            ## }
            
            while (x.new > problem.parameters$bx |
                   x.new < problem.parameters$ax |
                   y.new > problem.parameters$by |
                   y.new < problem.parameters$ay) {
                       x.direction <- sample(x=c(-1,1),size=1);
                       y.direction <- sample(x=c(-1,1),size=1);
                       dx.new <- dx*x.direction;
                       dy.new <- dy*y.direction;
                       x.new <- x.current + dx.new;
                       y.new <- y.current + dy.new;
                   }
            if (x.new != x.current+dx |
                y.new != y.current+dy) {
                weight = weight/2;
            }
            x.current <- x.new;
            y.current <- y.new;
        }
        
        if (!(x.current > problem.parameters$bx |
                x.current < problem.parameters$ax |
                y.current > problem.parameters$by |
                y.current < problem.parameters$ay)) {
            output.x = c(output.x, x.current);
            output.y = c(output.y, y.current);
        }
        weights[i] = weight;
    }
    
    unique.samples <- unique(sample(x=seq(1,n.simulations),
                                    size=100*n.simulations, replace=TRUE,
                                    prob=weights))
    while (length(unique.samples) < n.samples) {
        unique.samples <- unique(sample(x=seq(1,n.simulations),
                                        size=100*n.simulations, replace=TRUE,
                                        prob=weights));
    }
    unique.samples <- unique.samples[seq(1,n.samples)];

    out <- rbind(output.x[unique.samples],
                 output.y[unique.samples]);
    
    neighbor.distance <- matrix(nrow=n.samples,
                                ncol=n.samples-1);
    for (i in seq(1,n.samples)) {
        if (i==1) {
            neighbor.distance[i,] <-
                sapply(seq(2,n.samples),
                       function(x) {
                           return(sqrt(sum((out[,i]-out[,x])^2)));
                       });
        } else if (i==n.samples) {
            neighbor.distance[i,] <-
            sapply(seq(1,n.samples-1),
                   function(x) {
                       return(sqrt(sum((out[,i]-out[,x])^2)));
                   })
        } else {
            neighbor.distance[i,] <-
            sapply(c(seq(1,i-1),
                     seq(i+1,n.samples)),
                   function(x) {
                       return(sqrt(sum((out[,i]-out[,x])^2)));
                   });
        }
    }
    sigma2s <- ((apply(neighbor.distance, 1, min))*1)^2;
    output <- NULL;
    output$mus <- out;
    output$log.sigma2s <- log(sigma2s);
    return(output);
}

## performing gram-schmidt orthogonalization on the list of functions
## provided.
deriv.cross.term.intergral <- function(m, l,
                                       raw.function.list,
                                       problem.parameters,
                                       moments) {
    a=problem.parameters$a;
    b=problem.parameters$b;

    mu.m = raw.function.list[[m]][1];
    sigma2.m = raw.function.list[[m]][2];
    mu.l = raw.function.list[[l]][1];
    sigma2.l = raw.function.list[[l]][2];
    
    sigma2 = (1/raw.function.list[[m]][2] +
              1/raw.function.list[[l]][2])^(-1);
    mu = sigma2*(raw.function.list[[m]][1]/
                 raw.function.list[[m]][2] +
                 raw.function.list[[l]][1]/
                 raw.function.list[[l]][2]);
    C = product.coefficient(raw.function.list[[m]],
                            raw.function.list[[l]]);
    terms = matrix(nrow = 3, ncol = 3);
    ## term 1: \[ (b-x)^2 * N(x|mu_m, sigma^2_m) * N(x|mu_l,
    ##            sigma^2_l) =
    ##            (b^2-2bx+x^2) *
    ##            product.coef * ker(x| (1/sigma^2_m +
    ##            1/sigma^2_l)^{-1}(mu_m/sigma^2_m + mu_l/sigma^2_l,
    ##            (1/sigma^2_m + 1/sigma^2_l)^{-1}) \]
    terms[1,1] = b^2*moments[m,l,1] -2*b*moments[m,l,2] + moments[m,l,3];

    ## term 2: \[ (b-x)*-(x-a) * N(x|mu_m, sigma^2_m) * N(x|mu_l,
    ##            sigma^2_l) =
    ##            (-bx + ba + x^2 - xa) *
    ##            product.coef * ker(x| (1/sigma^2_m +
    ##            1/sigma^2_l)^{-1}(mu_m/sigma^2_m + mu_l/sigma^2_l,
    ##            (1/sigma^2_m + 1/sigma^2_l)^{-1}) \]

    terms[1,2] =
        a*b*moments[m,l,1] - a*moments[m,l,2] - b*moments[m,l,2] + moments[m,l,3];

    ## term 3: \[ (b-x)*(x-a)*-(x-mu_l)/sigma^2_l * N(x|mu_m, sigma^2_m) * N(x|mu_l,
    ##            sigma^2_l) =
    ## (-a b μ_l + a b x - a x^2 + a μ_l x - b x^2 + b μ_l x + x^3 - μ_l x^2)
    ##                  /sigma^2_l *
    ##            product.coef * ker(x| (1/sigma^2_m +
    ##            1/sigma^2_l)^{-1}(mu_m/sigma^2_m + mu_l/sigma^2_l,
    ##            (1/sigma^2_m + 1/sigma^2_l)^{-1}) \]

    terms[1,3] = (-a*b^2*mu.l*moments[m,l,1]+
        b*moments[m,l,2]*(mu.l*(2*a + b) + a*b) +
        moments[m,l,3]*(-mu.l*(a + 2*b) - b*(2*a + b)) +
        moments[m,l,4]*(a + 2*b + mu.l) -
        moments[m,l,5])/sigma2.l;

    ## term 4
    terms[2,1] = terms[1,2];

    ## terms 5: \[ -1^2*(x-a)^2 * N(x|mu_m, sigma^2_m) * N(x|mu_l,sigma^2_l) =
    ##             (x^2 - 2xa + a^2) *
    ##             N(x|mu_m, sigma^2_m) * N(x|mu_l,sigma^2_l) \]
    terms[2,2] = moments[m,l,3] - 2*a*moments[m,l,2] + a^2*moments[m,l,1];

    ## terms 6: \[ -1*(x-a)^2*(b-x)*-1*(x-mu_l)/sigma^2_l *
    ##                N(x|mu_m, sigma^2_m) * N(x|mu_l,sigma^2_l) =
    ## -a^2 b μ_l + a x (μ_l (a + 2 b) + a b) +
    ##  x^2 (-μ_l (2 a + b) - a (a + 2 b)) +
    ##  x^3 (2 a + b + μ_l) -
    ##  x^4 * all the rest \]
    terms[2,3] = (-a^2*b*mu.l*moments[m,l,1] +
        a*moments[m,l,2]*(mu.l*(a + 2*b) + a*b) +
        moments[m,l,3]*(-mu.l*(2*a + b) - a*(a + 2*b)) +
        moments[m,l,4]*(2*a + b + mu.l) -
        moments[m,l,5])/sigma2.l;

    ## term 7: \[ (x-a)*(b-x)^2*-1*(x-mu_m)/sigma^2_m *
    ##                N(x|mu_m, sigma^2_m) * N(x|mu_l,sigma^2_l) =
    ## -(a b^2 μ_m)/σ^2 + (b x (μ_m (2 a + b) + a b))/σ^2 -
    ##  (x^2 (μ_m (a + 2 b) + b (2 a + b)))/σ^2 +
    ##  (x^3 (a + 2 b + μ_m))/σ^2 - x^4/σ^2
    terms[3,1] = (-(a*b^2*mu.m)*moments[m,l,1] +
        (b*moments[m,l,2]*(mu.m*(2*a + b) + a*b)) -
        (moments[m,l,3]*(mu.m*(a + 2*b) + b*(2*a + b))) +
        (moments[m,l,4]*(a + 2*b + mu.m)) -
        moments[m,l,5]) /
        sigma2.m;

    ## term 8:
    terms[3,2] = (-(a^2*b*mu.m)*moments[m,l,1] +
        moments[m,l,2]*(a*(mu.m*(a + 2*b) + a*b)) -
        moments[m,l,3]*((mu.m*(2*a + b) + a*(a + 2*b))) +
        (moments[m,l,4]*(2*a + b + mu.m)) -
        moments[m,l,5]) /
        sigma2.m;

    ## term 9:
    terms[3,3] = (moments[m,l,1]*(a^2*b^2*mu.l*mu.m) -
        moments[m,l,2]*(a*b*(mu.l*(2*mu.m*(a + b) + a*b) + a*b*mu.m)) +
        moments[m,l,3]*((mu.l*(mu.m*(a^2 + 4*a*b + b^2) +
                               2*a*b*(a + b)) + a*b*(2*mu.m*(a + b) + a*b))) -
        moments[m,l,4]*((mu.l*(a^2 + 2*mu.m*(a + b) + 4*a*b + b^2) +
                         mu.m*(a^2 + 4*a*b + b^2) + 2*a*b*(a + b))) +
        moments[m,l,5]*((a^2 + mu.l*(2*(a + b) + mu.m) + 2*mu.m*(a + b) +
                         4*a*b + b^2)) -
        moments[m,l,6]*((2*((a + b) + mu.l + mu.m))))/
        (sigma2.m*sigma2.l);

    return ( product.coefficient(raw.function.list[[m]],
                                 raw.function.list[[l]])*
             sum(apply(X=terms, MARGIN=1, FUN =sum)) );
}

product.coefficient <- function(raw.function.params.1,
                                raw.function.params.2) {
    mu.1 = raw.function.params.1[1];
    sigma2.1 = raw.function.params.1[2];
    mu.2 = raw.function.params.2[1];
    sigma2.2 = raw.function.params.2[2];

    out = sqrt(1 /
               (4*pi^2*(sigma2.1 * sigma2.2))) *
        exp(-0.5*(mu.1^2/sigma2.1+mu.2^2/sigma2.2)) *
        exp(0.5*
            (1/sigma2.1 + 1/sigma2.2)^(-1)*
            (mu.1/sigma2.1 + mu.2/sigma2.2)^2);
    return (out);
}

## function.params$mu = c(mu.x, mu.y);
## function.params$epsilon = matrix(nrow=2,ncol=2,data=c(c(sigma2, rho*sigma2),
##                                                      c(rho*sigma2, sigma2)));
## function.params$l = l (double > 0)
## output = x^l*(1-x)^l*y^l*(1-y)^l*N(c(x,y)| mu, epsilon);
basis.function.normal.kernel <- function(x,y,function.params) {
    return( rep(x,length(y))^l*(1-rep(x,length(y)))^l*y^l*(1-y)^l*
            dmvnorm(x=cbind(rep(x,length(y)),
                            y),
                    mean = function.params$mu,
                    sigma = function.params$epsilon));
}

basis.function.normal.kernel.xy <- function(x,y,function.params) {
    return( cbind(sapply(x,
                         function(x,y) {basis.function.normal.kernel(x,y,
                                                       function.params)}, y)));
}

basis.function <- function(x,y, function.params, problem.parameters) {
    ## mean <- matrix(nrow=2,ncol=1,c(function.params[1], function.params[2]));
    ##     sigma <- matrix(nrow=2,ncol=2,
    ##                 data=c(c(function.params[3],
    ##                          sqrt(function.params[3]*function.params[4])*
    ##                          problem.parameters$rho),
    ##                        c(sqrt(function.params[3]*function.params[4])*
    ##                          problem.parameters$rho,
    ##                          function.params[4])));
    
    ## out = (x-problem.parameters$ax)*(problem.parameters$bx-x)*    
    ##     (y-problem.parameters$ay)*(problem.parameters$by-y)*
    ##     dmvnorm(x=matrix(nrow=length(y), ncol=2,data=c(rep(x,length(y)),
    ##                                                    y)),
    ##             mean, sigma);

    alpha.1 <- function.params[1];
    beta.1 <- problem.parameters$K.prime-
        function.params[1]+2;
    
    mu.1 <- (alpha.1)/(alpha.1+beta.1);
    sigma2.1 <- alpha.1*beta.1/
        ((alpha.1+beta.1)^2*(alpha.1+beta.1+1));


    alpha.2 <- function.params[2];
    beta.2 <- problem.parameters$K.prime-
        function.params[2]+2;
    
    mu.2 <- (alpha.2)/(alpha.2+beta.2);
    sigma2.2 <- alpha.2*beta.2/
        ((alpha.2+beta.2)^2*(alpha.2+beta.2+1));

    ## return (rep(x*(1-x),length(y))*y*(1-y) *
    ##         dmvnorm(cbind(rep(x,length(y)),
    ##                       y),
    ##                 mean=c(mu.1,mu.2),
    ##                 sigma=matrix(nrow=2,ncol=2,byrow=FALSE,
    ##                              data=9*c(c(sigma2.1, -0.8*sqrt(sigma2.1*sigma2.2)),
    ##                                     c(-0.8*sqrt(sigma2.1*sigma2.2), sigma2.2)))));
                                          
    
    return (dbeta(x,function.params[1],problem.parameters$K.prime-
                                         function.params[1]+2) *
            dbeta(y,function.params[2],problem.parameters$K.prime-
                                         function.params[2]+2));
}

basis.function.xy <- function(x,y,function.params,problem.parameters) {
    return( cbind(sapply(x,
                         function(x,y) {basis.function(x,y,
                                                       function.params,
                                                       problem.parameters)}, y)))
}

basis.function.yx <- function(x,y,function.params,problem.parameters) {
    return( rbind(sapply(y,
                         function(y,x) {basis.function(x,y,
                                                       function.params,
                                                       problem.parameters)}, x)))
}

norm.raw.function <- function(function.params, a, b) {
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
        
    Ms = Ms*PP*sqrt(2*pi*sigma2);
    
    out = sqrt((Ms[5] -
                Ms[4]*2*(a + b) +
                Ms[3]*(a^2 + b^2 + 4*a*b) -
                Ms[2]*2*(a*b^2 + b*a^2) +
                Ms[1]*a^2*b^2) *
               product.coefficient(function.params,
                                   function.params));
    return (out);
}

project.numeric <- function(function.1, function.2,
                            dx,dy) {
    ## TODO(georgid): check that dim of both functions is the same
    if (sum(dim(function.1) ==
            dim(function.2)) < 2 ) {
        stop("in project.numeric, functions not same dim");
    }
    projection <- sum(sapply(seq(1,length(function.1[,1])),
                      function(x) {
                          sum(function.1[x,]*
                              function.2[x,])*dy;
                      }))*dx;
    return(projection);
}

deriv.x.numeric <- function(function.1,
                            dx,dy) {

    x.length = length(function.1[,1]);
    
    out <- sapply(X=seq(1,length(function.1[1,])),
                  FUN = function(x) {
                      return((function.1[-1,x]-
                              function.1[-x.length,x])/dx);
                  });
    return(out);
}

deriv.y.numeric <- function(function.1,
                            dx,dy) {

    y.length = length(function.1[1,]);
    
    out <- sapply(X=seq(1,length(function.1[,1])),
                  FUN = function(x) {
                      return((function.1[x,-1]-
                              function.1[x,-y.length])/dy);
                  });
    return(t(out));
}

project <- function(raw.function.params.1,
                    raw.function.params.2,
                    problem.parameters,
                    x,y,
                    dx,dy) {
    
    ax = problem.parameters$ax;
    bx = problem.parameters$bx;
    ay = problem.parameters$ay;
    by = problem.parameters$by;
    
    mu.1.x = raw.function.params.1[1];
    mu.1.y = raw.function.params.1[2];
    sigma2.1.x = raw.function.params.1[3];
    sigma2.1.y = raw.function.params.1[4];

    mu.2.x = raw.function.params.2[1];
    mu.2.y = raw.function.params.2[2];
    sigma2.2.x = raw.function.params.2[3];
    sigma2.2.y = raw.function.params.2[4];

    integral.x = sum((x-ax)^2*
                     (bx-x)^2*
                     dnorm(x,mu.1.x,sqrt(sigma2.1.x))*
                     dnorm(x,mu.2.x,sqrt(sigma2.2.x)))*dx;

    integral.y = sum((y-ay)^2*
                     (by-y)^2*
                     dnorm(y,mu.1.y,sqrt(sigma2.1.y))*
                     dnorm(y,mu.2.y,sqrt(sigma2.2.y)))*dy;
    
    out = integral.x * integral.y;
    return (out);
}

project.dx.dx <- function(raw.function.params.1,
                          raw.function.params.2,
                          problem.parameters,
                          x,y,
                          dx,dy) {
    
    ax = problem.parameters$ax;
    bx = problem.parameters$bx;
    ay = problem.parameters$ay;
    by = problem.parameters$by;
    
    mu.1.x = raw.function.params.1[1];
    mu.1.y = raw.function.params.1[2];
    sigma2.1.x = raw.function.params.1[3];
    sigma2.1.y = raw.function.params.1[4];

    mu.2.x = raw.function.params.2[1];
    mu.2.y = raw.function.params.2[2];
    sigma2.2.x = raw.function.params.2[3];
    sigma2.2.y = raw.function.params.2[4];

    deriv.x.1 = ((problem.parameters$bx-x)+
                 ##
                 -(x-problem.parameters$ax)+
                 ##
                 -(x-problem.parameters$ax)*
                 (problem.parameters$bx-x)*
                 (x-raw.function.params.1[1])/
                 raw.function.params.1[3])*
        dnorm(x,
              raw.function.params.1[1],
              sqrt(raw.function.params.1[3]));

    deriv.x.2 = ((problem.parameters$bx-x)+
                 ##
                 -(x-problem.parameters$ax)+
                 ##
                 -(x-problem.parameters$ax)*
                 (problem.parameters$bx-x)*
                 (x-raw.function.params.2[1])/
                 raw.function.params.2[3])*
        dnorm(x,
              raw.function.params.2[1],
              sqrt(raw.function.params.2[3]));
    
    
    integral.x = sum(deriv.x.1*deriv.x.2)*dx;

    integral.y = sum((y-ay)^2*
                     (by-y)^2*
                     dnorm(y,mu.1.y,sqrt(sigma2.1.y))*
                     dnorm(y,mu.2.y,sqrt(sigma2.2.y)))*dy;
    
    out = integral.x * integral.y;
    return (out);
}

project.dy.dy <- function(raw.function.params.1,
                          raw.function.params.2,
                          problem.parameters,
                          x,y,
                          dx,dy) {
    
    ax = problem.parameters$ax;
    bx = problem.parameters$bx;
    ay = problem.parameters$ay;
    by = problem.parameters$by;
    
    mu.1.x = raw.function.params.1[1];
    mu.1.y = raw.function.params.1[2];
    sigma2.1.x = raw.function.params.1[3];
    sigma2.1.y = raw.function.params.1[4];

    mu.2.x = raw.function.params.2[1];
    mu.2.y = raw.function.params.2[2];
    sigma2.2.x = raw.function.params.2[3];
    sigma2.2.y = raw.function.params.2[4];

    deriv.y.1 = ((problem.parameters$by-y)+
                 ##
                 -(y-problem.parameters$ay)+
                 ##
                 -(y-problem.parameters$ay)*
                 (problem.parameters$by-y)*
                 (y-mu.1.y)/
                 sigma2.1.y)*
        dnorm(y,
              mu.1.y,
              sqrt(sigma2.1.y));

    deriv.y.2 = ((problem.parameters$by-y)+
                 ##
                 -(y-problem.parameters$ay)+
                 ##
                 -(y-problem.parameters$ay)*
                 (problem.parameters$by-y)*
                 (y-mu.2.y)/
                 sigma2.2.y)*
        dnorm(y,
              mu.2.y,
              sqrt(sigma2.2.y));
    
    
    integral.y = sum(deriv.y.1*deriv.y.2)*dy;

    integral.x = sum((x-ax)^2*
                     (bx-x)^2*
                     dnorm(x,mu.1.x,sqrt(sigma2.1.x))*
                     dnorm(x,mu.2.x,sqrt(sigma2.2.x)))*dx;
    
    out = integral.x * integral.y;
    return (out);
}

project.dx.dy <- function(raw.function.params.1,
                          raw.function.params.2,
                          problem.parameters,
                          x,y,
                          dx,dy) {
    
    ax = problem.parameters$ax;
    bx = problem.parameters$bx;
    ay = problem.parameters$ay;
    by = problem.parameters$by;
    
    mu.1.x = raw.function.params.1[1];
    mu.1.y = raw.function.params.1[2];
    sigma2.1.x = raw.function.params.1[3];
    sigma2.1.y = raw.function.params.1[4];

    mu.2.x = raw.function.params.2[1];
    mu.2.y = raw.function.params.2[2];
    sigma2.2.x = raw.function.params.2[3];
    sigma2.2.y = raw.function.params.2[4];

    deriv.x.1 = ((problem.parameters$bx-x)+
                 ##
                 -(x-problem.parameters$ax)+
                 ##
                 -(x-problem.parameters$ax)*
                 (problem.parameters$bx-x)*
                 (x-mu.1.x)/
                 sigma2.1.x)*
        dnorm(x,
              mu.1.x,
              sqrt(sigma2.1.x));

    deriv.y.2 = ((problem.parameters$by-y)+
                 ##
                 -(y-problem.parameters$ay)+
                 ##
                 -(y-problem.parameters$ay)*
                 (problem.parameters$by-y)*
                 (y-mu.2.y)/
                 sigma2.2.y)*
        dnorm(y,
              mu.2.y,
              sqrt(sigma2.2.y));
    
    
    integral.x = sum(deriv.x.1*
                     (x-ax)*(bx-x)*
                     dnorm(x,mu.2.x*sqrt(sigma2.2.x)))*
        dx;

    integral.y = sum((y-ay)*
                     (by-y)*
                     dnorm(y,mu.1.y,sqrt(sigma2.1.y))*
                     deriv.y.2)*
        dy;
    
    out = integral.x * integral.y;
    return (out);
}

project.dy.dx <- function(raw.function.params.1,
                          raw.function.params.2,
                          problem.parameters,
                          x,y,
                          dx,dy) {
    
    ax = problem.parameters$ax;
    bx = problem.parameters$bx;
    ay = problem.parameters$ay;
    by = problem.parameters$by;
    
    mu.1.x = raw.function.params.1[1];
    mu.1.y = raw.function.params.1[2];
    sigma2.1.x = raw.function.params.1[3];
    sigma2.1.y = raw.function.params.1[4];

    mu.2.x = raw.function.params.2[1];
    mu.2.y = raw.function.params.2[2];
    sigma2.2.x = raw.function.params.2[3];
    sigma2.2.y = raw.function.params.2[4];

    deriv.x.2 = ((problem.parameters$bx-x)+
                 ##
                 -(x-problem.parameters$ax)+
                 ##
                 -(x-problem.parameters$ax)*
                 (problem.parameters$bx-x)*
                 (x-mu.2.x)/
                 sigma2.2.x)*
        dnorm(x,
              mu.2.x,
              sqrt(sigma2.2.x));

    deriv.y.1 = ((problem.parameters$by-y)+
                 ##
                 -(y-problem.parameters$ay)+
                 ##
                 -(y-problem.parameters$ay)*
                 (problem.parameters$by-y)*
                 (y-mu.1.y)/
                 sigma2.1.y)*
        dnorm(y,
              mu.1.y,
              sqrt(sigma2.1.y));
    
    
    integral.x = sum(deriv.x.2*
                     (x-ax)*(bx-x)*
                     dnorm(x,mu.1.x*sqrt(sigma2.1.x)))*
        dx;

    integral.y = sum((y-ay)*
                     (by-y)*
                     dnorm(y,mu.2.y,sqrt(sigma2.2.y))*
                     deriv.y.1)*
        dy;
    
    out = integral.x * integral.y;
    return (out);
}
    
gram.schmidt <- function(problem.parameters, function.list,
                         dx) {
    x = seq(problem.parameters$a,
            problem.parameters$b,
            by=dx);
    K <- length(function.list);
    out = vector(mode = "list", length=K);
    for (k in seq(1,K)) {
        if (k==1) {
            ## only normalize
            norm = sqrt(sum(function.list[[k]](x)^2)*dx);
            out[[k]] <- function(x) {
                return( function.list[[k]](x)/norm );
            }
        } else {
            projections <- rep(NA, k-1);
            for (l in seq(1,k-1)) {
                projections[l] = sum(out[[l]](x)*current.function(x)*
                                             dx);
            }
            
        }
    }
    return (out);
}

univariate.solution.approx <- function(coefs,x,orthonormal.function.list,K) {
    out = rep(0,length(x));
    for (n in seq(1,K)) {
        out = out +
            coefs[n]*orthonormal.function.list[[n]];
    }
    return (out);
}

bivariate.solution.approx <- function(orthonormal.function.list,
                                      K,
                                      coefs) {
    out <- matrix(nrow=dim(orthonormal.function.list[[1]])[1],
                  ncol=dim(orthonormal.function.list[[1]])[2],
                  0);
    for (k.prime in seq(1,K)) {
        out = out +
            coefs[k.prime] *
            orthonormal.function.list[[k.prime]];
    }
    return (out);
}

bivariate.solution.approx.xy <- function(orthonormal.function.list,
                                         K,
                                         coefs,
                                         x.fc.index,
                                         y.fc.index) {
    out = sum(sapply(seq(1,K),
                     function(x) {
                         coefs[x]*orthonormal.function.list[[x]][x.fc.index,
                                                                 y.fc.index]
                     }));
    return (out);
}

## bivariate.solution.approx.yx <- function(coefficients,
##                                       coefs,
##                                       x,y,
##                                       raw.function.list,
##                                       problem.parameters) {
    
##     out = matrix(0, ncol=length(y),nrow=length(x));
##     for (n in seq(1,length(raw.function.list))) {
##         for (m in seq(1,length(raw.function.list))) {
##             out = out +
##                 coefs[n]*
##                 coefficients[n,m]*basis.function.yx(x,y,
##                                                     raw.function.list[[m]],
##                                                     problem.parameters)
##         }
##     }
##     return (out);
## }

univariate.solution.approx.dt <- function(coefs,A,x,K,
                                          orthonormal.function.list) {
    coefs.dt = -A %*% coefs;
    out = rep(0,length(x));
    for (n in seq(1,K)) {
        out = out +
            coefs.dt[n]*orthonormal.function.list[[n]];
    }
    return (out);
}


univariate.solution.approx.dx.dx.numeric <- function(coefs,x,K,
                                                     orthonormal.function.list) {
    out = rep(0,length(x)-2);
    for (n in seq(1,K)) {
        y = orthonormal.function.list[[n]];
        y.m.dx = y[-length(x)];
        y.p.dx = y[-1];
        y.dx.dx = (y.p.dx[-1]-
            2*y[seq(2,length(x)-1)] +
            y.m.dx[-length(y.m.dx)])/dx^2;
        out = out + coefs[n]*y.dx.dx;
    }
    return (out);
}

univariate.solution.approx.dx.dx <- function(coefs,
                                             coefficients,x,K,
                                             raw.function.list,
                                             problem.parameters) {
    a = problem.parameters$a;
    b = problem.parameters$b;
    out = rep(0,length(x));
    for (n in seq(1,K)) {

        Psi = rep(0,length(x));
        for (m in seq(1,K)) {
            mu = raw.function.list[[m]][1];
            sigma2 = raw.function.list[[m]][2];
            sigma = sqrt(sigma2);
            
            Psi = Psi +
                ##
                coefficients[n,m]*dnorm(x,mu,sigma)*
                (-1 +
                 ##
                 (b-x)*-1*(x-mu)/sigma2 +
                 ##
                 -1 +
                 ##
                 -(x-a)*-1*(x-mu)/sigma2 +
                 ##
                 (b-x)*-1*(x-mu)/sigma2 +
                 ##
                 -(x-a)*-1*(x-mu)/sigma2 +
                 ##
                 (x-a)*(b-x)*-1/sigma2 + 
                 ##
                 (x-a)*(b-x)*(x-mu)^2/(sigma2)^2)
        }
        out = out + coefs[n]*Psi;
    }
    return (out);
}

coefs.approx.mc <- function(problem.parameters,
                            poly.degree.x,
                            number.samples,
                            delta.t.min) {
    coefs.approx = rep(NA, poly.degree.x+1);
    for (i in seq(1,poly.degree.x+1)) {
        sampled.bm <- sample.bounded.bm.automatic(problem.parameters,
                                                  delta.t.min=1e-7,
                                                  number.samples=number.samples);
        plot(density(sampled.bm$points));
        integral=sum(sampled.bm$weights*
                     basis.function(sampled.bm$points,i))/number.samples;
        print(c(i,integral));
        coefs.approx[i] = integral;
    }
    return (coefs.approx);
}

coefs.approx.path.integral.mc <- function(problem.parameters,
                                          poly.degree.x,
                                          number.samples,
                                          delta.t) {
    coefs.approx = rep(NA, poly.degree.x+1);
    for (i in seq(1,poly.degree.x+1)) {
        basis.function.polynomial = basis.function.poly(i);
        sampled.bm =
            sample.bounded.bm.path.integral.accept.reject(problem.parameters,
                                                          delta.t,
                                                          number.samples,
                                                          basis.function.polynomial);
        plot(density(sampled.bm$points));
        integral=sum(sampled.bm$weights*sampled.bm$integrals)/number.samples;
        print(c(i,integral));
        coefs.approx[i] = integral;
    }
    return (coefs.approx);
}

sample.bounded.bm.path.integral.accept.reject <- function(problem.parameters,
                                                          delta.t,
                                                          number.samples,
                                                          basis.function.polynomial) {
    out = NULL;
    out$weights = rep(NA, number.samples);
    out$points = rep(NA, number.samples);
    out$integrals = rep(NA, number.samples);
    
    t = seq(delta.t, problem.parameters$t, by = delta.t);
    sigma = sqrt(problem.parameters$sigma.2);
    ff <- as.function(deriv(deriv(basis.function.polynomial, "x"), "x"),
                      vector=FALSE);

    for (i in seq(1,number.samples)) {
        x.current = problem.parameters$x.ic;
        es <- rnorm(n=length(t));
        weight = 1;
        integral = 0;
        
        for (j in seq(1,length(t))) {
            x.previous = x.current;
            t.current = t[j];
            x.current = x.previous + sigma*sqrt(delta.t)*es[j];
            integral = integral +
                0.5*sigma^2*
                ff((x.previous+x.current)/2)*delta.t;
            if (x.current <= problem.parameters$a ||
                x.current >= problem.parameters$b ) {
                x.current = x.previous - sigma*sqrt(delta.t)*es[j];
                weight = 0;
                break;
            }
        }
        out$weights[i] = weight;
        out$points[i] = x.current;
        out$integrals[i] = integral +
            as.function(basis.function.polynomial, vector=FALSE)(
                problem.parameters$x.ic);
    }
    return (out);
}

sample.bounded.bm.accept.reject <- function(problem.parameters,
                                            delta.t,
                                            number.samples) {
    out = NULL;
    out$weights = rep(NA, number.samples);
    out$points = rep(NA, number.samples);
    t = seq(delta.t, problem.parameters$t, by = delta.t);
    sigma = sqrt(problem.parameters$sigma.2);

    for (i in seq(1,number.samples)) {
        x.current = problem.parameters$x.ic;
        es <- rnorm(n=length(t));
        weight = 1;
        
        for (j in seq(1,length(t))) {
            x.previous = x.current;
            t.current = t[j];
            x.current = x.previous + sigma*sqrt(delta.t)*es[j];
            if (x.current <= problem.parameters$a ||
                x.current >= problem.parameters$b ) {
                x.current = x.previous - sigma*sqrt(delta.t)*es[j];
                weight = 0;
                break;
            }
        }
        out$weights[i] = weight;
        out$points[i] = x.current;
    }
    return (out);
}

sample.bounded.bm <- function(problem.parameters, delta.t, number.samples) {
    out = NULL;
    out$weights = rep(NA, number.samples);
    out$points = rep(NA, number.samples);
    
    t = seq(delta.t, problem.parameters$t, by = delta.t);
    sigma = sqrt(problem.parameters$sigma.2);

    for (i in seq(1,number.samples)) {
        x.current = problem.parameters$x.ic;
        es <- rnorm(n=length(t));
        weight = 1;
        
        for (j in seq(1,length(t))) {
            x.previous = x.current;
            t.current = t[j];
            x.current = x.previous + sigma*sqrt(delta.t)*es[j];
            if (x.current <= problem.parameters$a ||
                x.current >= problem.parameters$b ) {
                x.current = x.previous - sigma*sqrt(delta.t)*es[j];
                weight = weight * 0.5;
            }

            if (x.current <= problem.parameters$a ||
                x.current >= problem.parameters$b ) {
                stop("Outside boundary again; time step too big!!");
            }
        }
        out$weights[i] = weight;
        out$points[i] = x.current;
    }
    return (out);
}

sample.bounded.bm.automatic <- function(problem.parameters,
                                        delta.t.min,
                                        number.samples) {
    out = NULL;
    out$weights = rep(NA, number.samples);
    out$points = rep(NA, number.samples);
    out$counts = rep(NA, number.samples);
    
    t = problem.parameters$t;
    sigma = sqrt(problem.parameters$sigma.2);

    for (i in seq(1,number.samples)) {
        ## path = rep(0,100);
        ## times = rep(0,100);
        x.current = problem.parameters$x.ic;
        t.current = 0;
        d = min(abs(problem.parameters$a-x.current),
                abs(problem.parameters$b-x.current));
        delta.t = min((max(delta.t.min,
        (d/(sigma*2))^2)),
        t-t.current);
        
        weight = 1;
        count = 1;

        ## path[count]=x.current;
        ## times[count]=t.current;
        
        while (t.current < (t-delta.t.min)) {
            count = count + 1;
            d = min(abs(problem.parameters$a-x.current),
                    abs(problem.parameters$b-x.current));
            delta.t = min((max(delta.t.min,
            (d/(sigma*3))^2)),
            t-t.current);
            
            x.previous = x.current;
            t.current = t.current + delta.t;
            eps <- rnorm(1);
            x.current = x.previous + sigma*sqrt(delta.t)*eps;
            if (x.current <= problem.parameters$a ||
                x.current >= problem.parameters$b ) {
                x.current = x.previous - sigma*sqrt(delta.t)*eps;
                weight = weight * 0.5;
            }

            if (x.current <= problem.parameters$a ||
                x.current >= problem.parameters$b ) {
                warning("Outside boundary again; time step too big!!");
                break;
            }

            ## if (count < length(path)) {
            ##     path[count] = x.current;
            ##     times[count] = t.current;
            ## } else {
            ##     path = c(path, rep(0,100));
            ##     path[count] = x.current;
            ##     times[count] = t.current;
            ## }
        }
        ## plot(times[1:count], path[1:count], type="l",
        ##      ylim = c(problem.parameters$a,
        ##               problem.parameters$b))
        out$weights[i] = weight;
        out$points[i] = x.current;
        out$counts[i] = count;
    }
    return (out);
}

### NEW STUFF START ###
## derivative of kernel test function on (-1,1)
kernel.deriv.poly <- function(x,K,kernel) {

    ## P.0s = vector(mode="list", length = K);
    ## P.1s = vector(mode="list", length = K);
    ## derivs.doubles = matrix(nrow=length(x), ncol=K);

    ## for (k in seq(1,K)) {
    ##     if (k==1) {
    ##         P.0s[[k]] = mpoly(list(c("x"=0,coef=1),
    ##                                c("x"=2,coef=-1)))^2;
    ##         P.1s[[k]] = mpoly(list(c("x"=1,coef=-2)));
    ##     } else {
    ##         P.0s[[k]] = P.0s[[k-1]]^3 * P.0s[[1]];
    ##         P.1s[[k]] =
    ##             P.1s[[1]] *
    ##             P.1s[[k-1]] *
    ##             P.0s[[k-1]]^2 +
    ##                             P.0s[[k-1]] *
    ##                             P.0s[[1]] *
    ##                             (P.0s[[k-1]]*deriv(P.1s[[k-1]],"x") -
    ##                              deriv(P.0s[[k-1]],"x")*P.1s[[k-1]]);
    ##     }
    ##     P.k.0 = as.function(P.0s[[k]],vector=FALSE)(x);
    ##     signs.P.k.0 = sign(P.k.0);

    ##     P.k.1 = as.function(P.1s[[k]],vector=FALSE)(x);
    ##     signs.P.k.1 = sign(P.k.1);
        
    ##     derivs.doubles[,k] = signs.P.k.0*signs.P.k.1*
    ##         exp(-log(abs(P.k.0))+log(abs(P.k.1))+log(kernel(x)));
            
    ## }
    
    if (k==1) {
        return( second.derivative.poly(x,kernel) );
    } else if (k==2) {
        return( fourth.derivative.poly(x,kernel) )
    } else if (k==3) {
        return( sixth.derivative.poly(x,kernel) )
    } else {
        stop ("k above 3");
    }
}

second.derivative.poly <- function(x,kernel) {
    ## sign.numer = sign(6*x^4-2);
    ## out.exp.log = sign.numer*( exp(log(abs(6*x^4-2)) -
    ##                               4*log(1-x^2)));
    out =  kernel(x) * (6*x^4-2) / (1-x^2)^4;
    return (out);
}

fourth.derivative.poly <- function(x,kernel) {
    out= (4*kernel(x))*(30*x^10+45*x^8-132*x^6+58*x^4+6*x^2-3)/(x^2-1)^8
    return (out);
}

sixth.derivative.poly <- function(x,kernel) {
    out = kerne(x)*(8)*(x^2*(((15*x^2*(42*x^8+210*x^6-567*x^4+62*x^2+643)-7102)*x^2+1005)*x^2+270)-15)/(x^2-1)^12
    return (out);
}

eighth.derivative.poly <- function(x) {
    out = (16)*(22680*x^22+238140*x^20-502740*x^18-868455*x^16+2862720*x^14-2049012*x^12-473256*x^10+1190310*x^8-440216*x^6+3360*x^4+16380*x^2+105)/((x-1)^16*(x+1)^16);

    return (out);
}
### NEW STUFF END ###

### OLD STUFF START ###
select.alpha.beta <- function(problem.parameters) {
    a = problem.parameters$a;
    b = problem.parameters$b;
    x.ic = problem.parameters$x.ic;
    sigma2 = problem.parameters$sigma.2;
    t = problem.parameters$t;

    mean = x.ic;
    var = sigma2*t;

    ## the initial guess is alpha=beta=2;
    initial.guess <- c(log(2), log(2))
    log.alpha.beta <- optim(par=initial.guess,
                            fn = mean.var.sum.sq,
                            mean = mean,
                            var = var);

    alpha = max(2, round(exp(log.alpha.beta$par[1])));
    beta = max(2, round(exp(log.alpha.beta$par[2])));

    out = NULL;
    out$alpha = alpha;
    out$beta = beta;
    return (out);
}

mean.var.sum.sq <- function(x, mean, var) {
    alpha <- exp(x[1]);
    beta <- exp(x[2]);

    current.mean <- alpha/(alpha+beta);
    current.var <- alpha*beta/((alpha+beta)^2*(alpha+beta+1));
    
    out = (mean - current.mean)^2 + (var-current.var)^2;
    return (out);
}

univariate.solution <- function(x, problem.parameters) {
    a = problem.parameters$a;
    b = problem.parameters$b;
    x.ic = problem.parameters$x.ic;
    t = problem.parameters$t;
    sigma.2 = problem.parameters$sigma.2;
    number.terms = problem.parameters$number.terms;
    
    if ( x.ic < a || x.ic > b) {
        stop ("x.ic not in boundaries");
    }

    v1 = x.ic + seq(-number.terms,number.terms)*2*(b-a); ## all 1
    v2 = (a-(x.ic-a)) + seq(-number.terms,number.terms)*2*(b-a); ## all -1

    variance = sigma.2 * t;
    ## print(c(t,variance, sigma.2));

    out = rep(NA, length(x));
    for (i in seq(1,length(x))) {
        out[i] = sum(dnorm(x[i], mean=v1, sd=sqrt(variance)) - dnorm(x[i], mean=v2, sd=sqrt(variance)));
    }
    
    return (out);
};

x.power.integral.vector <- function(problem.parameters,
                                    poly.degree,
                                    alpha,
                                    beta) {
    a = problem.parameters$a;
    b = problem.parameters$b;
    x.ic = problem.parameters$x.ic;
    t = problem.parameters$t;
    sigma.2 = problem.parameters$sigma.2;
    number.terms = problem.parameters$number.terms;

    if ( x.ic < a || x.ic > b) {
        stop ("x.ic not in boundaries");
    }

    out <- rep(NA, poly.degree+1);

    for (m in seq(1,poly.degree+1)) {
        out[m] = exp(lbeta(2*alpha+(m-1)-1,2*beta-1) -
                     2*lbeta(alpha,beta));
    }
    return (out);
}

coefficients <- function(problem.parameters,
                         polynomials.table,
                         polynomial.kernel,
                         kernel,
                         poly.degree.x,
                         number.derivs) {
    x.0 <- problem.parameters$x.ic;
    sigma2 <- problem.parameters$sigma.2;
    t <- problem.parameters$t;
    coefs <- rep(NA, poly.degree.x+1);

    for (i in seq(1,poly.degree.x+1)) {
        polynomial <- polynomials.table[[i]]*polynomial.kernel;
        derivatives <- vector(mode="list", length=number.derivs);
        derivatives.doubles <- rep(NA, length=number.derivs);
        for (k in seq(1,number.derivs)) {
            if (k==1) {
                derivatives[[k]] = deriv(deriv(polynomial, "x"), "x");
            } else {
                derivatives[[k]] = deriv(deriv(derivatives[[k-1]], "x"), "x");
            }
            derivatives.doubles[k] = as.function(derivatives[[k]], vector=FALSE)(x.0);
        }
        
        coefficient = as.function(polynomial,
                                  vector=FALSE)(x.0) +
                                              sum(t^seq(1,number.derivs)/
                                                    factorial(seq(1,number.derivs))*
                                                    (0.5*sigma2)^seq(1,number.derivs)*
                                                                 derivatives.doubles);
        coefs[i] = coef;
    }
    return (coefs);
}

coefficients.exact <- function(kernel,
                               basis.function.table,
                               problem.parameters,
                               N) {
    a = problem.parameters$a;
    b = problem.parameters$b;
    dx = (b-a)/N;

    
    coefs = rep(NA, length(basis.function.table));
    
    for (i in seq(1,length(basis.function.table))) {
        integral=sum(univariate.solution(seq(0,N-1)*dx, problem.parameters)*
                     basis.function.table[[i]](seq(0,N-1)*dx)*
                     dx)
        print(c(i,integral));
                                     
    }
    return (coefs);
}


integral.table <- function(problem.parameters, poly.degree) {
    a = problem.parameters$a;
    b = problem.parameters$b;
    x.ic = problem.parameters$x.ic;
    t = problem.parameters$t;
    sigma.2 = problem.parameters$sigma.2;
    number.terms = problem.parameters$number.terms;

    if ( x.ic < a || x.ic > b) {
        stop ("x.ic not in boundaries");
    }

    v1 = x.ic + seq(-number.terms,number.terms)*2*(b-a); ## all 1
    v2 = (a-(x.ic-a)) + seq(-number.terms,number.terms)*2*(b-a); ## all -1
    image.means = c(v1,v2);
    image.weights = c(rep(1,length(v1)), rep(-1,length(v2)));
    
    variance = sigma.2 * t;

    L.table = matrix(nrow = length(image.means), ncol = poly.degree+1);
    for (i in seq(1,length(image.means))) {
        ## standarized left boundary
        alpha <- (a - image.means[i])/sqrt(variance);
        ## standarized right boundary
        beta <- (b - image.means[i])/sqrt(variance);
        normalizing.weight <- pnorm(q=beta, mean=0, sd=1) - pnorm(q=alpha, mean=0, sd=1);
        
        for (j in seq(1,poly.degree+1)) {
            if (j==1) {
                ## corresponds to poly order 0;
                L.table[i,j] = normalizing.weight;
            } else if (j==2) {
                ## corresponds to poly order 1;
                L.table[i,j] = -(dnorm(x=beta,mean=0,sd=1) - dnorm(x=alpha,mean=0,sd=1));
            } else {
                L.table[i,j] = -(beta^(j-1-1)*dnorm(x=beta,mean=0,sd=1)-alpha^(j-1-1)*dnorm(x=alpha,mean=0,sd=1)) +
                    (j-1-1)*L.table[i,j-2];
            }
        }
    }

    table = matrix(nrow = length(image.means), ncol = poly.degree+1);
    for (i in seq(1,length(image.means))) {
        for (j in seq(1,poly.degree+1)) {
            table[i,j] = 
            sum(choose(j-1,seq(0,j-1)) * sqrt(variance)^seq(0,j-1) *
            image.means[i]^((j-1) - seq(0,j-1)) * L.table[i,seq(1,j)])
        }
    }
    return (table);
}

image.means.weights <- function(problem.parameters) {
    a = problem.parameters$a;
    b = problem.parameters$b;
    x.ic = problem.parameters$x.ic;
    t = problem.parameters$t;
    sigma.2 = problem.parameters$sigma.2;
    number.terms = problem.parameters$number.terms;

    if ( x.ic < a || x.ic > b) {
        stop ("x.ic not in boundaries");
    }

    v1 = x.ic + seq(-number.terms,number.terms)*2*(b-a); ## all 1
    v2 = (a-(x.ic-a)) + seq(-number.terms,number.terms)*2*(b-a); ## all -1
    image.means = c(v1,v2);
    image.weights = c(rep(1,length(v1)), rep(-1,length(v2)));
    
    out = NULL;
    out$image.means <- image.means;
    out$image.weights <- image.weights;
    return (out);
}

integrate.polynomial <- function(poly, integrals.table) {
    poly.list <- unclass(poly);

    out = 0;
    for (p in poly.list) {
        x.power.index = 1;
 
        if (is.na(p["x"]) == FALSE) {
            x.power.index = p["x"] + 1;
        }
         out = out +
            p["coef"] * integrals.table[x.power.index];
    }
    return (as.double(out));
}

apply.generator <- function(poly, k, problem.parameters.x) {
    sigma2.x = problem.parameters.x$sigma.2;
    
    if (k==0) {
        return (poly);
    } else {
        new.poly = mpoly(list(c(x=0,coef=0.5*sigma2.x)))*
            deriv(deriv(poly, "x"), "x");
        return (apply.generator(new.poly, k-1, problem.parameters.x));
    }
}

basis.functions.normal.kernel <- function(rho,
                                          l,
                                          sigma2,
                                          dx,dy,
                                          std.dev.factor) {
    sigma <- sqrt(sigma2);
    x.nodes <- c(seq(0.5-sqrt(2), 0.5,
                     by=std.dev.factor*sigma*sqrt(1-rho)),
                 seq(0.5, 0.5+sqrt(2),
                     by=std.dev.factor*sigma*sqrt(1-rho)));
    x.nodes <- unique(x.nodes);
    
    y.nodes <- c(seq(0.5-sqrt(2), 0.5,
                     by=std.dev.factor*sigma*sqrt(1+rho)),
                 seq(0.5, 0.5+sqrt(2),
                     by=std.dev.factor*sigma*sqrt(1+rho)));
    y.nodes <- unique(y.nodes);

    xy.nodes <- rbind(rep(x.nodes,each=length(y.nodes)),
                      rep(y.nodes,length(x.nodes)));
    ## plot(xy.nodes[1,], xy.nodes[2,]);
    theta <- pi/4;
    Rot.mat <- matrix(nrow=2,ncol=2,
                      data=c(c(sin(theta), -cos(theta)),
                             c(cos(theta), sin(theta))));
    xieta.nodes <- Rot.mat %*% (xy.nodes - c(0.5,0.5)) + c(0.5,0.5);
    xieta.nodes <- xieta.nodes[, (xieta.nodes[1,] >= 0) & (xieta.nodes[1,] <= 1) &
                                 (xieta.nodes[2,] >= 0) & (xieta.nodes[2,] <= 1)];
    
    plot(xieta.nodes[1,],xieta.nodes[2,]);
    ## sort.bases <- sort.int(sqrt((xieta.nodes[1,]-0.5)^2 +
    ##                             (xieta.nodes[2,]-0.5)^2),
    ##                        index.return = TRUE);
    ## xieta.nodes = xieta.nodes[,sort.bases$ix];
    
    K = dim(xieta.nodes)[2];
    x <- seq(0,1,by=dx);
    y <- seq(0,1,by=dy);

    function.list <- vector("list", K);

    ## function.params$mu = c(mu.x, mu.y);
    ## function.params$epsilon = matrix(nrow=2,ncol=2,data=c(c(sigma2, rho*sigma2),
    ##                                                      c(rho*sigma2, sigma2)));
    ## function.params$l = l (double >= 1)
    ## output = x^l*(1-x)^l*y^l*(1-y)^l*N(c(x,y)| mu, epsilon);
    function.params=NULL;
    for (k in seq(1,K)) {
        function.params$mu = c(xieta.nodes[1,k], xieta.nodes[2,k]);
        function.params$epsilon = matrix(nrow=2,ncol=2,data=c(c(sigma2, rho*sigma2),
                                                              c(rho*sigma2, sigma2)));
        function.params$l = 1;
        function.list[[k]] = basis.function.normal.kernel.xy(x, y, function.params);
    }

    if (ceiling(sqrt(K)) <= 15) {
        par(mfrow=c(ceiling(sqrt(K)),
                    ceiling(sqrt(K))),
            mai=c(0.1,0.1,0.1,0.1));
        for (k in seq(1,K)) {
            contour(x,y,function.list[[k]]);
        }
    }
    return (function.list);
}

orthonormal.functions <- function(function.list,
                                  dx,dy,
                                  x,y,
                                  PLOT.SOLUTION) {
    K <- length(function.list);
    norms <- rep(NA, K);
    coefficients <- matrix(0, nrow=K, ncol=K);
    orthonormal.function.list = vector(mode="list",
                                       length=K);
    orthonormal.function.list <- function.list;
    ## gram-schmidt START ##
    par(mfrow=c(ceiling(sqrt(K)),
                ceiling(sqrt(K))));
    par(mai = c(0.1,0.1,0.1,0.1));    
    for (k in seq(1,K)) {
        if (k==1) {
            ## only normalize
            norm = sqrt(project.numeric(function.list[[k]],
                                        function.list[[k]],
                                        dx,dy));
            
            coefficients[k,k] = 1/norm;
            norms[k] = norm;
            orthonormal.function.list[[k]] =
                function.list[[k]]/norm;
        } else {
            orthonormal.function.list[[k]] = function.list[[k]];
            for (l in seq(1,k-1)) {
                coefficients[k,l] = -project.numeric(orthonormal.function.list[[k]],
                                                     orthonormal.function.list[[l]],
                                                     dx,dy);
                
                orthonormal.function.list[[k]] =
                    orthonormal.function.list[[k]] +
                    coefficients[k,l]*orthonormal.function.list[[l]];

                orthonormal.function.list[[k]] =
                    orthonormal.function.list[[k]] /
                    sqrt(project.numeric(orthonormal.function.list[[k]],
                                         orthonormal.function.list[[k]],
                                         dx,dy));
            }
            coefficients[k,k] = 1;
            
            norm = sqrt(project.numeric(orthonormal.function.list[[k]],
                                        orthonormal.function.list[[k]],
                                        dx,dy));
            norms[k] = norm;
            orthonormal.function.list[[k]] =
                orthonormal.function.list[[k]]/norm;
            coefficients[k,]=coefficients[k,]/norm;
        }
        if (PLOT.SOLUTION) {
            contour(x,y,orthonormal.function.list[[k]]);
        }
        print(k);
    }
    ## gram schmidt END ###p
    return (orthonormal.function.list);
}

system.matrices <- function(orthonormal.function.list,
                            dx,dy) {
    
    K = length(orthonormal.function.list);
    ## SYSTEM MATRICES START ##
    derivative.xx.matrix <- matrix(nrow = K,
                                   ncol = K);
    derivative.yy.matrix <- matrix(nrow = K,
                                   ncol = K);
    derivative.xy.matrix <- matrix(0,
                                   nrow = K,
                                   ncol = K);
    derivative.yx.matrix <- matrix(0,
                                   nrow = K,
                                   ncol = K);
    for (k in seq(1,K)) {
        current.basis.dx.k = deriv.x.numeric(orthonormal.function.list[[k]],
                                             dx,dy);
        current.basis.dy.k = deriv.y.numeric(orthonormal.function.list[[k]],
                                             dx,dy);
        
        for (l in seq(1,K)) {
            current.basis.dx.l = deriv.x.numeric(orthonormal.function.list[[l]],
                                             dx,dy);
            current.basis.dy.l = deriv.y.numeric(orthonormal.function.list[[l]],
                                             dx,dy);


            derivative.xx.matrix[k,l] =
                project.numeric(current.basis.dx.k,
                                current.basis.dx.l,
                                dx,dy);
                    
            derivative.yy.matrix[k,l] =
                project.numeric(current.basis.dy.k,
                                current.basis.dy.l,
                                dx,dy);

            derivative.xy.matrix[k,l] =
                project.numeric(current.basis.dx.k[,-length(y)],
                                current.basis.dy.l[-length(x),],
                                dx,dy);

            derivative.yx.matrix[k,l] =
                project.numeric(current.basis.dy.k[-length(x),],
                                current.basis.dx.l[,-length(y)],
                                dx,dy);
        }
        print(k);
    }
    ## SYSTEM MATRICES END ##

    mass.mat <- matrix(nrow=K,
                       ncol=K);
    
    for (k.prime in seq(1,K)) {
        for (l.prime in seq(1,K)) {
            mass.matrix.entry =
                project.numeric(orthonormal.function.list[[k.prime]],
                                orthonormal.function.list[[l.prime]],
                                dx,dy);
            ## print(c(i,j,mass.matrix.entry));
            mass.mat[k.prime,l.prime]=mass.matrix.entry;
        }
    }

    out <- NULL;
    out$derivative.xx.matrix <- derivative.xx.matrix;
    out$derivative.xy.matrix <- derivative.xy.matrix;
    out$derivative.yx.matrix <- derivative.yx.matrix;
    out$derivative.yy.matrix <- derivative.yy.matrix;
    out$mass.mat <- mass.mat;

    return(out);
}

blackbox <- function(function.list,
                     orthonormal.function.list,
                     system.mats,
                     problem.parameters.original,
                     dx, dy,
                     PLOT.SOLUTION,
                     MINIMIZE.REMAINDER) {
    source("2-d-solution.R");
    source("../classical-solution/2-d-solution.R");

    problem.parameters = rescale.problem(problem.parameters.original);
    K = length(function.list);

    K.prime = sqrt(K);
    problem.parameters$K.prime <- K.prime;
    ## gram schmidt START ##
    x = seq(problem.parameters$ax,
            problem.parameters$bx,
            by=dx);
    y = seq(problem.parameters$ay,
            problem.parameters$by,
            by=dy);
        
    ## SYSTEM MATRICES START ##
    derivative.xx.matrix <- system.mats$derivative.xx.matrix;
    derivative.yy.matrix <- system.mats$derivative.yy.matrix;
    derivative.xy.matrix <- system.mats$derivative.xy.matrix;
    derivative.yx.matrix <- system.mats$derivative.yx.matrix;
    ## SYSTEM MATRICES END ##
        
    stiff.mat <- matrix(nrow=K,
                        ncol=K);
    stiff.mat = -0.5*
                 problem.parameters$sigma.2.x*
                 derivative.xx.matrix +
                 -problem.parameters$rho*
                 sqrt(problem.parameters$sigma.2.x)*
                 sqrt(problem.parameters$sigma.2.y)*
                 0.5*(derivative.xy.matrix+
                      derivative.yx.matrix) +
                 -0.5*problem.parameters$sigma.2.y*
                 derivative.yy.matrix;

    mass.mat <- system.mats$mass.mat;
    ## SYSTEM MATRICES END ###
    
    ## ## eigenvalues START ###
    eig <- eigen(solve(mass.mat) %*% stiff.mat);
    ## eig <- eigen(stiff.mat);
    ## ## eigenvalues END ###
    
    ## ## ICs START ###
    x.ic.index = which(abs(x-problem.parameters$x.ic)<=dx/2);
    y.ic.index = which(abs(y-problem.parameters$y.ic)<=dy/2);
    x.fc.index = which(abs(x-problem.parameters$x.fc)<=dx/2);
    y.fc.index = which(abs(y-problem.parameters$y.fc)<=dy/2);

    small.t.solution <- bivariate.solution.classical(dx,dy,
                                                     problem.parameters,
                                                     PLOT.SOLUTION);
    IC.true <- small.t.solution$big.solution;
    IC.vec <- rep(NA, K);
    for (k in seq(1,K)) {
        IC.vec[k] <- sum(apply(IC.true*orthonormal.function.list[[k]],
                               1,
                               sum)*dy)*dx;
    }
    IC.vec <- solve(system.mats$mass.mat, IC.vec);
    problem.parameters$t <- 1 - small.t.solution$tt;

    if (PLOT.SOLUTION) {
        IC.approx <- bivariate.solution.approx(orthonormal.function.list,
                                               K,
                                               IC.vec);

        ## for (kk in seq(1,length(orthonormal.function.list))) {
        ##     if (kk==1) {
        ##         plot(x,orthonormal.function.list[[kk]][,y.ic.index],type="l",
        ##              ylim=c(-max(abs(orthonormal.function.list[[kk]][,y.ic.index])),
        ##                     max(abs(orthonormal.function.list[[kk]][,y.ic.index]))));
        ##     } else {
        ##         lines(x,orthonormal.function.list[[kk]][,y.ic.index]);
        ##     }
        ## }
        
        plot(x, IC.true[,y.ic.index], type = "l", col = "red");
        lines(x, IC.approx[,y.ic.index], col = "black");

        plot(y, IC.true[x.ic.index,], type = "l", col = "red");
        lines(y, IC.approx[x.ic.index,], col = "black")
    }
    
    ## ## APPROX SOLUTION START ###
    if (K==1) {
        coefs = (eig$vectors) *
            as.double(exp(eig$values * problem.parameters$t)) *
            t(eig$vectors) * IC.vec;
    } else {
        coefs = (eig$vectors) %*%
            diag(exp(eig$values * problem.parameters$t)) %*%
            t(eig$vectors) %*% IC.vec;
    }

    ## approx.sol <- bivariate.solution.approx(orthonormal.function.list,
    ##                                         K,
    ##                                         coefs);
    
    ## min.obs = min(apply(approx.sol, 1, min));
    ## min.index.row = which(apply(approx.sol, 1, min) == min.obs);
    ## min.index.col = which(approx.sol[min.index.row,] ==
    ##                       min(approx.sol[min.index.row,]));
    ## print(paste("min.obs = ", min.obs));
    ## ## print(min.index.row);
    ## ## print(min.index.col);

    ## problem.parameters.x = problem.parameters;
    ## problem.parameters.x$x.ic = problem.parameters$x.ic;
    ## problem.parameters.x$a = problem.parameters$ax;
    ## problem.parameters.x$b = problem.parameters$bx;
    ## problem.parameters.x$sigma.2 = problem.parameters$sigma.2.x;
    ## problem.parameters.x$x.ic = problem.parameters$x.ic;
    ## problem.parameters.x$t <- 1;

    ## problem.parameters.y = problem.parameters;
    ## problem.parameters.y$x.ic = problem.parameters$y.ic;
    ## problem.parameters.y$a = problem.parameters$ay;
    ## problem.parameters.y$b = problem.parameters$by;
    ## problem.parameters.y$sigma.2 = problem.parameters$sigma.2.y;
    ## problem.parameters.y$x.ic = problem.parameters$y.ic;
    ## problem.parameters.y$t <- 1;

    ## true.sol <- univariate.solution(x,problem.parameters.x) %*%
    ##     t(univariate.solution(y,problem.parameters.y));
       
    if (PLOT.SOLUTION) {
        approx.sol <- bivariate.solution.approx(orthonormal.function.list,
                                                K,
                                                coefs);

        ## ## function.params$mu = c(mu.x, mu.y);
        ## ## function.params$epsilon = matrix(nrow=2,ncol=2,data=c(c(sigma2, rho*sigma2),
        ## ##                                                      c(rho*sigma2, sigma2)));
        ## ## function.params$l = l (double > 0)
        ## ## output = x^l*(1-x)^l*y^l*(1-y)^l*N(c(x,y)| mu, epsilon);
        ## approx.sol.pos <- function(coefs.pos) {
        ##     rho <- problem.parameters$rho;
        ##     mu.x <- coefs.pos[1];
        ##     sigma.x <- exp(coefs.pos[2]);
        ##     mu.y <- coefs.pos[3];
        ##     sigma.y <- exp(coefs.pos[4]);

        ##     function.params=NULL;
        ##     function.params$mu=c(mu.x,mu.y);
        ##     function.params$epsilon =
        ##         matrix(nrow=2,ncol=2,
        ##                data=c(c(sigma.x^2, rho*sigma.x*sigma.y),
        ##                       c(rho*sigma.x*sigma.y, sigma.y^2)));
        ##     function.params$l=0.5;
            
        ##     out <- basis.function.normal.kernel.xy(x,y,function.params);
        ##     return(out);
        ## }

        ## optimization.function <- function(coefs.pos) {
        ##     out <- sqrt(sum(apply((approx.sol-approx.sol.pos(coefs.pos))^2,
        ##                           1, sum)*dy)*dx);
        ##     ## print(out);
        ##     return(out);
        ## }
        
        ## opp <- optim(par=c(0.6, log(0.4^2), 0.2, log(0.2^2)),
        ##              fn=optimization.function);

        ## approx.sol.positive <- approx.sol.pos(opp$par);
        ## par(mfrow=c(1,2));
        ## contour(x,y,approx.sol.positive);
        ## contour(x,y,approx.sol);
        ## approx.sol <- approx.sol.positive;
        
        problem.parameters.x = problem.parameters;
        problem.parameters.x$x.ic = problem.parameters$x.ic;
        problem.parameters.x$a = problem.parameters$ax;
        problem.parameters.x$b = problem.parameters$bx;
        problem.parameters.x$sigma.2 = problem.parameters$sigma.2.x;
        problem.parameters.x$x.ic = problem.parameters$x.ic;
        problem.parameters.x$t <- 1;
        
        problem.parameters.y = problem.parameters;
        problem.parameters.y$x.ic = problem.parameters$y.ic;
        problem.parameters.y$a = problem.parameters$ay;
        problem.parameters.y$b = problem.parameters$by;
        problem.parameters.y$sigma.2 = problem.parameters$sigma.2.y;
        problem.parameters.y$x.ic = problem.parameters$y.ic;
        problem.parameters.y$t <- 1;
        
        true.sol <- univariate.solution(x,problem.parameters.x) %*%
            t(univariate.solution(y,problem.parameters.y));

        plot(x,approx.sol[,y.fc.index], type = "l", col = "black", lty="dashed");
        lines(x,true.sol[,y.fc.index], col = "red");
        abline(v=x[x.fc.index], col="blue", lwd=2);
        
        plot(y,approx.sol[x.fc.index,], type="l",
             lty="dashed", col = "black");
        lines(y,true.sol[x.fc.index,],col="red");
        abline(v=y[y.fc.index], col="blue", lwd=2);
        
        contour(x,y,approx.sol, nlevels = 50);
        points(problem.parameters$x.ic,
               problem.parameters$y.ic,
               pch=20,
               col="green");
        points(problem.parameters$x.fc,
               problem.parameters$y.fc,
               pch=20,
               col="red");
        
        print(approx.sol[x.fc.index, y.fc.index]);
        png("contour.png");
        contour(x,y,approx.sol, nlevels = 50);
        points(problem.parameters$x.ic,
               problem.parameters$y.ic,
               pch=20,
               col="green");
        points(problem.parameters$x.fc,
               problem.parameters$y.fc,
               pch=20,
               col="red");
        ## persp(x,y,approx.sol, theta = pi/2);
        ## points(x[min.index.row], y[min.index.col], col="red");
        dev.off();
    }

    out = sum(sapply(seq(1,K),
                     function(x) {
                         coefs[x]*orthonormal.function.list[[x]][x.fc.index,
                                                                 y.fc.index]
                     }));
    
    return (out);
}
