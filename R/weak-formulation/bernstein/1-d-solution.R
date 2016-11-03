library("mpoly");
library("pencopula");
library("Matrix");

IC.approx <- function(x,IC.vec) {
    out = rep(0,length(x));
    for (n in seq(1,n.basis)) {
        out = out +
            IC.vec[n]*basis.function.init(n,n.basis)(x);
    }
    return (out);
}

univariate.solution.approx <- function(x,coefs) {
    out = rep(0,length(x));
    for (n in seq(1,n.basis)) {
        current.basis = basis.function.init(n, n.basis);
        out = out +
            coefs[n]*current.basis(x);
    }
    return (out);
}

basis.function.init <- function(n,n.basis) {
    ## length of all input vectors must be the same!

    current.basis.function <- function(x) {
        ## return (x^(n)*(1-x)^(1+(n.basis-n)));
        return (choose(n.basis+1,n)*(x^(n)*(1-x)^(1+n.basis-n)));
    }
    return (current.basis.function);
}

basis.function.init.dx <- function(n, n.basis) {
    current.basis.function.dx <- function(x) {
        return (choose(n.basis+1,n)*(n*x^(n-1)*(1-x)^(1+(n.basis-n)) -
                (1+(n.basis-n))*x^n*(1-x)^(1+(n.basis-n-1))));

    }
    return (current.basis.function.dx);
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
