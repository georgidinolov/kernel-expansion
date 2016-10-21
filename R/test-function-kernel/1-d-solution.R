library("mpoly");
### NEW STUFF START ###
## derivative of kernel test function on (-1,1)
second.derivative.poly <- function(x) {
    out = (6*x^4-2) / (1-x^2)^4;
    return (out);
}

fourth.derivative.poly <- function(x) {
    out = 4*(30*x^10 + 45*x^8 - 132*x^6 + 58*x^4 + 6*x^2 -3)/
        (x^2-1)^8;
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
