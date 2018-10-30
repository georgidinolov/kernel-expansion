args = commandArgs(trailingOnly = TRUE)

library("data.table")
source("/soe/gdinolov/PDE-solvers/src/kernel-expansion/documentation/chapter-2/rogers-estimator/rogers-estimator.R")
load("/soe/gdinolov/PDE-solvers/src/kernel-expansion/documentation/chapter-2/rogers-estimator/phi-grid.Rdata")


files = fread(args[1], header=FALSE)
data_size = as.double(args[2])
nn = nrow(files)

sigma_x = rep(NA, nn)
sigma_y = rep(NA, nn)
rho.rogers = rep(NA, nn)
rho.classic = rep(NA, nn)


for (file.index in seq(1,nn)) {
    dat = fread(files[file.index, V1])
    sigma_x[file.index] = dat[, mean( (x_T - mean(x_T))^2 )]
    sigma_y[file.index] = dat[, mean( (y_T - mean(y_T))^2 )]
    rho.rogers[file.index] = dat[, rogers.est(x_0,y_0,x_T,y_T,b,d,a,c)]
    rho.rogers[file.index] = phi.function.inverse(rho.rogers[file.index], phi.grid)
    rho.classic[file.index] = dat[, sqrt(mean((x_T-mean(x_T))*(y_T-mean(y_T))))]
}

out = data.table(sigma_x=sigma_x, sigma_y=sigma_y, rho=rho.rogers)
fwrite(x=out, file=paste0("/soe/gdinolov/PDE-solvers/src/kernel-expansion/documentation/chapter-2/results/mle-results-rho-0.95-n-", data_size, "/rogers-results.csv"))

out = data.table(sigma_x=sigma_x, sigma_y=sigma_y, rho=rho.classic)
fwrite(x=out, file=paste0("/soe/gdinolov/PDE-solvers/src/kernel-expansion/documentation/chapter-2/results/mle-results-rho-0.95-n-", data_size, "/classic-results.csv"))