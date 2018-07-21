ff <- function(t, beta) {
    return (1/t^4 * exp(-beta/t))
}

ddf <- function(t,beta) {
    return (4/t^2 - 2*beta/t^3)
}

dddf <- function(t,beta) {
    return (-8/t^4 + 6*beta/t^4)
}

laplace <- function(tt, t.star, beta) {
    C <- exp(ff(t.star, beta))
    inv.scale <- -ddf(t.star, beta)
    out <- 1/C * exp(-1/2*inv.scale*(tt-t.star)^2)
    return (out)
}
