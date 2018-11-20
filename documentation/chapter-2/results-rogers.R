args = commandArgs(trailingOnly = TRUE)

require("mvtnorm")
require("data.table")
source("~/PDE-solvers/src/kernel-expansion/documentation/chapter-2/rogers-estimator/rogers-estimator.R")
load("~/PDE-solvers/src/kernel-expansion/documentation/chapter-2/rogers-estimator/phi-grid.Rdata")
opt.wrapper <- function(par, x_T, y_T) {
    sigmax = exp(par[1])
    sigmay = exp(par[2])
    rho = 2*exp(par[3])/(exp(par[3])+1)-1
    
    out = -sum(mvtnorm::dmvnorm(x=matrix(nrow=length(x_T),ncol=2,
                              data=c(x_T,y_T)),
                     mean=c(0,0),
                     sigma=matrix(nrow=2,ncol=2,
                                  data=c(sigmax^2,
                                         rho*sigmax*sigmay,
                                         rho*sigmax*sigmay,
                                         sigmay^2)),
                     log=TRUE)) -
        log(sigmax) - log(sigmay) - log(2) - 2*log((rho+1)/2)
                     
    return (out)
}


files = fread(args[1], header=FALSE)
data_size = as.double(args[2])
rho_data = args[3]
nn = nrow(files)

sigma_x = rep(NA, nn)
sigma_y = rep(NA, nn)
rho.rogers = rep(NA, nn)
rho.classic = rep(NA, nn)

dat.all = data.table()

for (file.index in seq(1,nn)) {
    dat = fread(files[file.index, V1])
    dat.all = rbind(dat.all, dat)

    rho.start = 0.0
    rho.tilde.start = (rho.start + 1)/2
    params.out = optim(par=c(0,0,log(rho.tilde.start/(rho.tilde.start+1))),
                       fn=opt.wrapper, method=c("Nelder-Mead"),
                       x_T = dat[, x_T],
                       y_T = dat[, y_T])

    sigma.x.out = exp(params.out$par[1])
    sigma.y.out = exp(params.out$par[2])
    rho.out = 2*exp(params.out$par[3])/(exp(params.out$par[3])+1) - 1
    print(c(sigma.x.out, sigma.y.out, rho.out))
    

    sigma_x[file.index] = dat[, sqrt(mean( (x_T - 0)^2 )) ]  # sigma.x.out
    sigma_y[file.index] =  dat[, sqrt(mean( (y_T - 0)^2 )) ] # sigma.y.out
    rho.out = mean(dat[, x_T*y_T]/(sigma_x[file.index]*sigma_y[file.index]))
    rs = dat[, rogers.est(x_0/sigma.x.out,
                          y_0/sigma.y.out,
                          x_T/sigma.x.out,
                          y_T/sigma.y.out,
                          b/sigma.x.out,
                          d/sigma.y.out,
                          a/sigma.x.out,
                          c/sigma.y.out)]
    rho.rogers[file.index] = phi.function.inverse(mean(rs), phi.grid)
    rho.classic[file.index] = rho.out # dat[, mean( x_T*y_T )]# rho.out
}

out = data.table(sigma_x=sigma_x, sigma_y=sigma_y, rho=rho.rogers)
fwrite(x=out, file=paste0("~/PDE-solvers/src/kernel-expansion/documentation/chapter-2/results/mle-results-rho-", rho_data, "-n-", data_size, "/rogers-results.csv"))

out = data.table(sigma_x=sigma_x, sigma_y=sigma_y, rho=rho.classic)
fwrite(x=out, file=paste0("~/PDE-solvers/src/kernel-expansion/documentation/chapter-2/results/mle-results-rho-", rho_data, "-n-", data_size, "/classic-results.csv"))
