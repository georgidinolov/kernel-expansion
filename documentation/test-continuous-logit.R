tau.rho.tilde <- 0.0001
N = 5000000
es <- tau.rho.tilde*rnorm(n=N)

rho.tilde = 0.0
rho.tildes = rep(NA, N)

for (i in seq(1,N)) {
    drho.tilde = es[i]
    rho.tilde = rho.tilde + drho.tilde
    rho.tildes[i] = rho.tilde
}

rhos = (exp(rho.tildes) - 1)/
    (exp(rho.tildes) + 1)

par(mfrow=c(2,1),
    mar=c(1,2,1,1))
plot(rho.tildes, type = "l")
plot(rhos, type = "l")
