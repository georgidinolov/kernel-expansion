library("mpoly");

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
