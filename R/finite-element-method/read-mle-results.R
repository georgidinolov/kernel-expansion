rm(list=ls());
path = "~/research/PDE-solvers/data/";

files.list.16 <-
    list.files(path = "~/research/PDE-solvers/data/", pattern = "*order-16-5e-3.csv");
files.list.32 <-
    list.files(path = "~/research/PDE-solvers/data/", pattern = "*order-32-5e-3.csv");
files.list.64 <-
    list.files(path = "~/research/PDE-solvers/data/", pattern = "*order-64-5e-3.csv");
files.list.128 <-
    list.files(path = "~/research/PDE-solvers/data/", pattern = "*order-128-5e-3.csv");

length(files.list.16);
length(files.list.32);
length(files.list.64);
length(files.list.128);

## files.list.64 = files.list.64[seq(1,100)];
## files.list.128 = files.list.128[seq(1,50)];

##############################################
sigma.xs.16 <- rep(NA, length(files.list.16));
sigma.ys.16 <- rep(NA, length(files.list.16));
rhos.16 <- rep(NA, length(files.list.16));

for (i in seq(1, length(files.list.16))) {
    result <- read.csv(file = paste(path, files.list.16[i], sep = ""),
                       header=T);
    sigma.xs.16[i] = result[[1]];
    sigma.ys.16[i] = result[[2]];
    rhos.16[i] =  result[[3]];
}

par(mfrow=c(2,2));
hist(sigma.xs.16, prob=T); lines(density(sigma.xs.16)); abline(v=1.0, lwd=2, col="red");
hist(sigma.ys.16, prob=T); lines(density(sigma.ys.16)); abline(v=1.0, lwd=2, col="red");
plot(density(rhos.16)); abline(v=0.6, lwd=2, col="red");
print(c(mean(rhos.16), median(rhos.16)));
##############################################

##############################################
sigma.xs.32 <- rep(NA, length(files.list.32));
sigma.ys.32 <- rep(NA, length(files.list.32));
rhos.32 <- rep(NA, length(files.list.32));

for (i in seq(1, length(files.list.32))) {
    result <- read.csv(file = paste(path, files.list.32[i], sep = ""),
                       header=T);
    sigma.xs.32[i] = result[[1]];
    sigma.ys.32[i] = result[[2]];
    rhos.32[i] =  result[[3]];
}

par(mfrow=c(2,2));
hist(sigma.xs.32, prob=T); lines(density(sigma.xs.32)); abline(v=1.0, lwd=2, col="red");
hist(sigma.ys.32, prob=T); lines(density(sigma.ys.32)); abline(v=1.0, lwd=2, col="red");
plot(density(rhos.32)); abline(v=0.6, lwd=2, col="red");
print(c(mean(rhos.32), median(rhos.32)));
##############################################

##############################################
sigma.xs.64 <- rep(NA, length(files.list.64));
sigma.ys.64 <- rep(NA, length(files.list.64));
rhos.64 <- rep(NA, length(files.list.64));

for (i in seq(1, length(files.list.64))) {
    result <- read.csv(file = paste(path, files.list.64[i], sep = ""),
                       header=T);
    sigma.xs.64[i] = result[[1]];
    sigma.ys.64[i] = result[[2]];
    rhos.64[i] =  result[[3]];
}

par(mfrow=c(2,2));
hist(sigma.xs.64, prob=T); lines(density(sigma.xs.64)); abline(v=1.0, lwd=2, col="red");
hist(sigma.ys.64, prob=T); lines(density(sigma.ys.64)); abline(v=1.0, lwd=2, col="red");
plot(density(rhos.64)); abline(v=0.6, lwd=2, col="red");
##############################################

##############################################
sigma.xs.128 <- rep(NA, length(files.list.128));
sigma.ys.128 <- rep(NA, length(files.list.128));
rhos.128 <- rep(NA, length(files.list.128));

for (i in seq(1, length(files.list.128))) {
    result <- read.csv(file = paste(path, files.list.128[i], sep = ""),
                       header=T);
    sigma.xs.128[i] = result[[1]];
    sigma.ys.128[i] = result[[2]];
    rhos.128[i] =  result[[3]];
}
lines(density(rhos.128), col="blue");
##############################################

par(mfrow=c(2,2));
plot(density(sigma.xs.64), col="blue", lwd=2); lines(density(sigma.xs.128)); abline(v=1.0, lwd=2, col="red");
plot(density(sigma.ys.64), col="blue", lwd=2); lines(density(sigma.ys.128)); abline(v=1.0, lwd=2, col="red");
plot(density(rhos.64), col="blue", lwd=2); lines(density(rhos.128)); abline(v=0.6, lwd=2, col="red");
