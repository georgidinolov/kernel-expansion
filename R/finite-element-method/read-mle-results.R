rm(list=ls());
library("parallel");
library("ggplot2");
library("latex2exp");
library("reshape2");
source("2-d-solution.R");
path = "~/PDE-solvers/data/";

## MLE ## 
data.files.list <-
    as.list(paste(path, "data-set-", seq(1,500), ".csv", sep=""));
datum <- load.data.from.csv(data.files.list[[1]]);

rho.rogers <- estimator.rodgers(data.files.list, 0.6);

unconstrained.mles <- lapply(X=data.files.list[seq(1, 500)],
                             FUN=function(x) {
                                 mle.estimator.no.boundary(load.data.from.csv(x), 1, 1, 0.0)
                             })


the.rest <- unconstrained.mles[x=sample(seq(1,length(unconstrained.mles)),
                                        size=500-length(unconstrained.mles),
                                        replace=TRUE)]

unconstrained.mles <- c(unconstrained.mles, the.rest)

unlisted.unconstrained.mles <- unlist(unconstrained.mles)
sigma.x.classical <- unlist(lapply(X=unconstrained.mles, function(x) { return (x$sigma.x.mle) } ))
sigma.y.classical <- unlist(lapply(X=unconstrained.mles, function(x) { return (x$sigma.y.mle) } ))
rho.classical <- unlist(lapply(X=unconstrained.mles, function(x) { return (x$rho) } ))


files.list.16 <-
    list.files(path = path, pattern = "*order-16-5e-3-linear-256.csv");
files.list.32 <-
    list.files(path = path, pattern = "*order-32-5e-3-linear-256.csv");
files.list.64 <-
    list.files(path = path, pattern = "*order-64-5e-3-linear-256.csv");
files.list.128 <-
    list.files(path = path, pattern = "*order-128-5e-3-linear-256.csv");

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

RMSE.sigmax.galerkin = round(sqrt(mean((sigma.xs.16 -
                                         mean(sigma.xs.16))^2)),
                              digits=3);
RMSE.sigmax.classical = round(sqrt(mean((sigma.x.classical[1:length(sigma.xs.16)] -
                                        mean(sigma.x.classical[1:length(sigma.xs.16)]))^2)),
                             digits=3);

RMSE.sigmay.galerkin = round(sqrt(mean((sigma.ys.16 -
                                         mean(sigma.ys.16))^2)),
                              digits=3);
RMSE.sigmay.classical = round(sqrt(mean((sigma.y.classical[1:length(sigma.ys.16)] -
                                        mean(sigma.y.classical[1:length(sigma.ys.16)]))^2)),
                             digits=3);

RMSE.rho.galerkin = round(sqrt(mean((rhos.16 -
                                         mean(rhos.16))^2)),
                              digits=3);
RMSE.rho.classical = round(sqrt(mean((rho.classical[1:length(rhos.16)] -
                                        mean(rho.classical[1:length(rhos.16)]))^2)),
                           digits=3);

RMSE.rho.rogers = round(sqrt(mean((rho.rogers[1:length(rhos.16)]-
                                   mean(rho.rogers[1:length(rhos.16)]))^2)),
                        digits=3)

pdf("../../documentation/mle-comparison-sigma-x.pdf", width=3, height=3)
xx <- data.frame(Galerkin.approx=sigma.xs.16,
                 classical.likelihood=sigma.x.classical[1:length(sigma.xs.16)])
data <- melt(xx);
ggplot(data, aes(x=value, fill=variable)) +
    geom_density(alpha=0.5) +
    scale_fill_manual(values=c("blue", "green"),
                      labels=c(paste("RMSE(Galerkin) = ", RMSE.sigmax.galerkin, sep=""),
                               paste("RMSE(Classical) = ", RMSE.sigmax.classical,sep="")),
                      name = "") +
    geom_vline(xintercept=1.0, col="red",
               lwd=1.0) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position=c(-0.01,0.99),
          axis.line = element_line(colour = "black")) +
    xlab(TeX("$\\hat{\\sigma}_x$")) 
dev.off()


pdf("../../documentation/mle-comparison-sigma-y.pdf", width=3, height=3)
xx <- data.frame(Galerkin.approx=sigma.ys.16,
                 classical.likelihood=sigma.y.classical[1:length(sigma.ys.16)]);
data <- melt(xx);
ggplot(data, aes(x=value, fill=variable)) +
    geom_density(alpha=0.5) +
    scale_fill_manual(values=c("blue", "green"),
                      labels=c(paste("RMSE(Galerkin) = ", RMSE.sigmay.galerkin, sep=""),
                               paste("RMSE(Classical) = ", RMSE.sigmay.classical, sep="")),
                      name = "") +
        geom_vline(xintercept=1.0, col="red",
               lwd=1.0) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position=c(-0.01,0.99)) +
    xlab(TeX("$\\hat{\\sigma}_y$")) +
    ylab("")
dev.off()


pdf("../../documentation/mle-comparison-rho.pdf", width=3, height=3)
xx <- data.frame(Galerkin.approx=rhos.16,
                 classical.likelihood=rho.classical[1:length(rhos.16)],
                 Rogers.estimator=rho.rogers[1:length(rhos.16)]);
data <- melt(xx);
ggplot(data, aes(x=value, fill=variable)) +
    geom_density(alpha=0.5) +
    scale_fill_manual(values=c("blue", "green", "yellow"),
                      labels=c(paste("RMSE(Galerkin) = ", RMSE.rho.galerkin, sep=""),
                               paste("RMSE(Classical) = ", RMSE.rho.classical, sep=""),
                               paste("RMSE(Rogers) = ", RMSE.rho.rogers, sep="")),
                      name = "") +
        geom_vline(xintercept=0.6, col="red",
               lwd=1.0) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position=c(-0.01,0.95)) +
    xlab(TeX("$\\hat{\\rho}$")) +
    ylab("");
dev.off()

print(paste("RMSE_{rho}(classical)=", sqrt(mean((rho.classical[1:length(rhos.16)] -
                                                     mean(rho.classical[1:length(rhos.16)]))^2))))
print(paste("RMSE_{rho}(Galerkin)=", sqrt(mean((rhos.16 -
                                                     mean(rhos.16))^2))))
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
hist(sigma.xs.32, prob=T);
lines(density(sigma.xs.32));
abline(v=1.0, lwd=2, col="red");
lines(density(sigma.x.classical), col="blue");

hist(sigma.ys.32, prob=T); lines(density(sigma.ys.32)); abline(v=1.0, lwd=2, col="red"); lines(density(sigma.y.classical), col="blue");

plot(density(rhos.32)); abline(v=0.6, lwd=2, col="red"); lines(density(rho.classical), col="blue");

print(c(mean(rhos.32), median(rhos.32), mean(rho.classical)));
print(paste("RMSE(rho.Galerkin.32) = ",
            sqrt(mean( (rhos.32 - mean(rhos.32))^2 ))))
print(c(sd(rhos.32), sd(rho.classical)))
print(c(sqrt(mean((rhos.32 - mean(rhos.32))^2)),
        sqrt(mean((rho.classical - mean(rho.classical))^2))))
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
hist(sigma.xs.64, prob=T);
lines(density(sigma.xs.64));
abline(v=1.0, lwd=2, col="red");
lines(density(sigma.x.classical), col="blue");

hist(sigma.ys.64, prob=T);
lines(density(sigma.ys.64));
abline(v=1.0, lwd=2, col="red");
lines(density(sigma.y.classical), col="blue");

plot(density(rhos.64));
abline(v=0.6, lwd=2, col="red");
lines(density(rho.classical), col="blue");

print(c(mean(rhos.64), median(rhos.64), mean(rho.classical)));
print(c(sd(rhos.64), sd(rho.classical)));
print(c(sqrt(mean((rhos.64 - mean(rhos.64))^2)),
        sqrt(mean((rho.classical - mean(rho.classical))^2))))

print(paste("RMSE(rho.Galerkin.64) = ",
            sqrt(mean( (rhos.64 - mean(rhos.64))^2 ))))
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

print(c(mean(rhos.128), median(rhos.128), mean(rho.classical)));
print(c(sd(rhos.128), sd(rho.classical)))

print(c(sqrt(mean((rhos.128 - mean(rhos.128))^2)),
        sqrt(mean((rho.classical - mean(rho.classical))^2))))
##############################################

par(mfrow=c(2,2));
plot(density(sigma.xs.64), col="blue", lwd=2); lines(density(sigma.xs.128)); abline(v=1.0, lwd=2, col="red");
plot(density(sigma.ys.64), col="blue", lwd=2); lines(density(sigma.ys.128)); abline(v=1.0, lwd=2, col="red");
plot(density(rhos.64), col="blue", lwd=2); lines(density(rhos.128), lwd=2); abline(v=0.6, lwd=2, col="red");


############# ALL #############################
par(mfrow=c(2,2))
plot(density(sigma.xs.16), type="l", lty=1, lwd=1, ylim=c(0,14), col="black");
lines(density(sigma.xs.32), lwd=2, lty=2, col = "red");
lines(density(sigma.xs.64), lwd=2, lty=3, col="blue");
lines(density(sigma.xs.128), lwd=2, lty=4, col="black");
abline(v=1.0,lwd=2,col="red");
## 
plot(density(sigma.ys.16), type="l", lty=1, lwd=1, ylim=c(0,14), col="black");
lines(density(sigma.ys.32), lwd=2, lty=2, col = "red");
lines(density(sigma.ys.64), lwd=2, lty=3, col="blue");
lines(density(sigma.ys.128), lwd=2, lty=4, col="black");
abline(v=1.0,lwd=2,col="red");
##
plot(density(rhos.16), type="l", lty=1, lwd=1, ylim=c(0,14), col="black");
lines(density(rhos.32), lwd=2, lty=2, col = "red");
lines(density(rhos.64), lwd=2, lty=3, col="blue");
lines(density(rhos.128), lwd=2, lty=4, col="black");
abline(v=0.6,lwd=2,col="red");



## coefs
## data: data lenght is a power of 2!
## xx

########### COMPLEX FFT ##############
true.sol = c(0,0.0585938,0.109375,0.152344,0.1875,0.214844,0.234375,0.246094,0.25,0.246094,0.234375,0.214844,0.1875,0.152344,0.109375,0.0585938,0,-0.0585938,-0.109375,-0.152344,-0.1875,-0.214844,-0.234375,-0.246094,-0.25,-0.246094,-0.234375,-0.214844,-0.1875,-0.152344,-0.109375,-0.0585938,-0);

x = as.matrix(read.csv(file="~/research/PDE-solvers/data/out.csv",header=F));
fcoefs = complex(real=x[,2],imaginary=x[,3]);
n = length(fcoefs);
xx = seq(0,2,length.out=100);
L = xx[length(xx)] - x[1];
es = vector(mode="list", length=n);
for (i in seq(1,n)) {
    if (i==1) {
        es[[i]] = rep(1,length(xx)) * fcoefs[i];
    } else if (i==n/2+1) {
        es[[i]] = cos(n*pi*xx) * fcoefs[i];
    } else if (i < n/2+1) {
        es[[i]] = complex(real = cos(2*pi*xx/L*(i-1)),
                          imaginary = sin(2*pi*xx/L*(i-1))) *
            fcoefs[i];
    } else if (i > n/2+1) {
        es[[i]] = complex(real = cos(2*pi*xx/L*(i-n-1)),
                          imaginary = sin(2*pi*xx/L*(i-n-1))) *
            fcoefs[i];
    }
}

solution = complex(real=rep(0,length(xx)), imaginary=rep(0,length(xx)));
for (i in seq(1,n)) {
    solution = solution + es[[i]];
}
solution = solution/n;

plot(xx, Re(solution));
lines(xx, Im(solution), col = "blue", lty=2);
lines(seq(0,L,length.out=length(true.sol)), true.sol);
lines(xx, xx*(1-xx), col="blue");
lines(xx, exp(-0.5 * (xx-0.5)^2 / 0.01), col="blue");
###################### REAL FFT #################################
true.sol = as.matrix(read.csv(file="~/research/PDE-solvers/data/trig-test.csv",
                              header=F));
true.sol = c(0,0.382683,0.707107,0.92388,1,0.92388,0.707107,0.382683,1.22465e-16);
plot(seq(0,1,length.out=length(true.sol)), true.sol);

xx = seq(0,1,length.out=200);
y = as.matrix(read.csv(file="~/research/PDE-solvers/data/out.csv",header=F))[,2];
n = length(y);
y.real = y[seq(1,n/2+1)];
y.imag = c(0, y[seq(n, n/2+2)], 0);
solution = rep(0,length(xx));

for (i in seq(1,n/2+1)) {
    if (i==1) {
        solution =
            solution + rep(y.real[i], length(xx));
    } else if (i==n/2+1) {
        solution =
            solution + y.real[i]*cos(n*pi*xx);
    } else {
        solution =
            solution + 2*y.real[i]*cos((i-1)*2*pi*xx) -
            2*y.imag[i]*sin((i-1)*2*pi*xx);
    }
}

solution = solution / n;
plot(xx,solution,col="blue");
lines(xx,solution,col="blue");
lines(seq(0,1,length.out=length(true.sol)), true.sol);
## lines(seq(0,1,length.out=length(true.sol[4,])), true.sol[4,]);
####################

####################

y.true = c(3.72665e-06,8.13816e-06,1.73297e-05,3.59844e-05,7.28609e-05,0.000143858,0.000276968,0.000519976,0.000951908,0.00169928,0.00295796,0.00502086,0.00831038,0.0134129,0.0211097,0.0323965,0.0484811,0.0707465,0.100669,0.139683,0.188995,0.249352,0.3208,0.402452,0.492325,0.587282,0.683124,0.774837,0.856997,0.924285,0.972053,0.996856,0.996856,0.972053,0.924285,0.856997,0.774837,0.683124,0.587282,0.492325,0.402452,0.3208,0.249352,0.188995,0.139683,0.100669,0.0707465,0.0484811,0.0323965,0.0211097,0.0134129,0.00831038,0.00502086,0.00295796,0.00169928,0.000951908,0.000519976,0.000276968,0.000143858,7.28609e-05,3.59844e-05,1.73297e-05,8.13816e-06,3.72665e-06);

y = c(15.7918,-13.0268,7.31241,-2.79313,0.725954,-0.128388,0.0154437,-0.00126778,6.75884e-05,-5.1039e-06,-2.02142e-06,-1.7884e-06,-1.52427e-06,-1.29704e-06,-1.101e-06,-9.31956e-07,-7.86148e-07,-6.60292e-07,-5.51582e-07,-4.57639e-07,-3.76468e-07,-3.06399e-07,-2.46038e-07,-1.94227e-07,-1.50008e-07,-1.12587e-07,-8.13151e-08,-5.5664e-08,-3.52116e-08,-1.96285e-08,-8.6679e-09,-2.15866e-09,8.88178e-16,-4.39405e-08,-8.80066e-08,-1.32324e-07,-1.77021e-07,-2.22223e-07,-2.6806e-07,-3.1466e-07,-3.62151e-07,-4.10659e-07,-4.60304e-07,-5.11196e-07,-5.63425e-07,-6.17055e-07,-6.72104e-07,-7.2852e-07,-7.86148e-07,-8.44676e-07,-9.03566e-07,-9.61951e-07,-1.01848e-06,-1.07192e-06,-1.08047e-06,-2.41396e-06,2.7996e-05,-0.00045362,0.00468481,-0.0321595,0.144401,-0.414321,0.72021,-0.639967);

n = length(y);
y.real = y[seq(1,n/2+1)];
y.imag = c(0, y[seq(n, n/2+2)], 0);
solution = rep(0,length(xx));

for (i in seq(1,n/2+1)) {
    if (i==1) {
        solution =
            solution + rep(y.real[i], length(xx));
    } else if (i==n/2+1) {
        solution =
            solution + 2*y.real[i]*cos(n*pi*xx);
    } else {
        solution =
            solution + 2*y.real[i]*cos((i-1)*2*pi*xx) -
            2*y.imag[i]*sin((i-1)*2*pi*xx);
    }
}

solution = solution / n;

plot(xx,solution);
lines(seq(0,1,length.out=length(y.true)), y.true, col="red");
sol = as.matrix(read.csv("~/research/PDE-solvers/data/trig-test.csv",header=F));
lines(seq(0,1,length.out=length(sol[16,])), sol[16,]);

lines(xx, sin(2*pi*xx*2), col="red");
####################


#### ODD EXTENSION ####
## before.extension = as.matrix(read.csv(file="~/research/PDE-solvers/data/trig-test.csv", header=F));
## contour(before.extension);

## odd.extension.1 = as.matrix(read.csv(file="~/research/PDE-solvers/data/odd-extension-1.csv", header=F));
## n = (dim(odd.extension.1)[1]-1)/2;
## contour(z=odd.extension.1[seq(1,n+1), seq(1,2*n+1)], xlab="x", ylab="y");

odd.extension = as.matrix(read.csv(file="~/research/PDE-solvers/data/odd-extension.csv", header=F));
n = dim(odd.extension)[1]-1;
contour(odd.extension[seq(1,n),seq(1,n)], xlab="x", ylab="y");

## USING R START 
fcoefs.R = matrix(n,n,data=NA);
## FIRST FFT ON ROWS
for (i in seq(1,n)) {
    fcoefs.R[i,] = fft(odd.extension[i,seq(1,n)]);
}
## SECOND FFT ON COLUMNS
for (i in seq(1,n)) {
    fcoefs.R[,i] = fft(fcoefs.R[,i]);
}

signal = matrix(n,n,data=NA);
## RECONSTRUCTING DISCRETE SIGNAL 1
for (i in seq(0,n-1)) {
    for (j in seq(0,n-1)) {
        signal[i+1,j+1] = 0;

        for (k in seq(0,n-1)) {
            for (l in seq(0,n-1)) {
                signal[i+1,j+1] = signal[i+1,j+1] +
                    fcoefs.R[k+1,l+1] * exp(1i*2*pi * (i/n*k + j/n*l ) )
            }
        }
        signal[i+1,j+1] = signal[i+1,j+1]/n^2;
    }
}
contour(Re(signal));

## RECONSTRUCTING DISCRETE SIGNAL 2
x = seq(0,2,by=1/(n/2));
y = seq(0,2,by=1/(n/2));

xx = matrix(n+1, n+1, data=rep(x, length(y)), byrow=F);
pyy = matrix(n+1, n+1, data=rep(y, length(x)), byrow=T); 

signal = matrix(n+1, n+1, data=0i);

for (k in seq(0,n-1)) {
    for (l in seq(0,n-1)) {
        signal = signal +
            fcoefs.R[k+1,l+1] *
            exp(1i*2*pi * xx/2 * k) *
            exp(1i*2*pi * yy/2 * l);
    }
}
signal = signal/n^2;
contour(Re(signal));

## RECONSTRUCTING CONTINUOUS SIGNAL
x = seq(0,2,length.out=100);
y = seq(0,2,length.out=100);

xx = matrix(length(x), length(y), data=rep(x, length(y)), byrow=F);
yy = matrix(length(x), length(y), data=rep(y, length(x)), byrow=T); 

signal = matrix(length(x), length(y), data=0i);
frequencies = c(seq(0,n/2), seq(-n/2+1,-1));

for (k in seq(1,length(frequencies))) {
    if (k == n/2) {
        fk = cos(n*pi * xx/2);
    } else {
        fk = exp(1i*2*pi * xx/2 * frequencies[k]);
    }
    
    for (l in seq(1,length(frequencies))) {
        if (l==n/2) {
            fl = cos(n*pi * yy/2);
        } else {
            fl = exp(1i*2*pi * yy/2 * frequencies[l]);
        }
        
        signal = signal +
            fcoefs.R[k,l] *
            fk *
            fl;
    }
}
signal = signal/n^2;
contour(Re(signal));
contour(Im(signal));
persp(Re(signal));
## USING R END


#### FFT ODD.EXTENSION ####
fft.odd.extension = as.matrix(read.csv(file="~/research/PDE-solvers/last-element-fourier-interpolant-FFT.csv", header=F));
dim(fft.odd.extension);

n = dim(fft.odd.extension)[2];
fcoefs = matrix(nrow=n, ncol=n, data=NA);
for (i in seq(1,n)) {
    for (j in seq(1,n)) {
        fcoefs[i,j] = complex(real=fft.odd.extension[2*(i-1)+1,j],
                              imaginary=fft.odd.extension[2*i,j]);
    }
}

## RECONSTRUCTING CONTINUOUS SIGNAL
x = seq(0,1,length.out=n+1);
y = seq(0,1,length.out=n+1);

xx = matrix(length(x), length(y), data=rep(x, length(y)), byrow=F);
yy = matrix(length(x), length(y), data=rep(y, length(x)), byrow=T); 

signal = matrix(length(x), length(y), data=0i);
frequencies = c(seq(0,n/2), seq(-n/2+1,-1));
Lx = x[length(x)] - x[1];
Ly = y[length(x)] - y[1];

for (k in seq(1,length(frequencies))) {
    if (k == n/2) {
        fk = cos(n*pi * xx/Lx);
    } else {
        fk = exp(1i*2*pi * xx/Lx * frequencies[k]);
    }
    
    for (l in seq(1,length(frequencies))) {
        if (l==n/2) {
            fl = cos(n*pi * yy/Ly);
        } else {
            fl = exp(1i*2*pi * yy/Ly * frequencies[l]);
        }
        
        signal = signal +
            fcoefs[k,l] *
            fk *
            fl;
    }
    print (k);
}

signal = signal/n^2;
par(mfrow=c(2,2));
contour(Re(signal));
persp(Re(signal), theta=15, phi=10);
contour(Im(signal));


dxdx.lin= as.matrix(read.csv(file="~/research/PDE-solvers/system_matrix_dx_dx.csv", header=F));

dxdx.fft= as.matrix(read.csv(file="~/research/PDE-solvers/system_matrix_dx_dx.csv", header=F));

evec = as.matrix(read.csv(file="~/research/PDE-solvers/evec-matrix.csv", header=F));
evec.1 = as.matrix(read.csv(file="~/research/PDE-solvers/evec-matrix-1.csv", header=F));
eval = as.matrix(read.csv(file="~/research/PDE-solvers/eval-matrix.csv", header=F));
IC = as.matrix(read.csv(file="~/research/PDE-solvers/IC-matrix.csv", header=F));

evec %*% exp(diag(eval[,1]) * 0.99) %*% t(evec) %*% IC;

element.fft.interpolant = as.matrix(read.csv(file="~/research/PDE-solvers/last-element-fourier-interpolant.csv", header=F));
element.linear.interpolant = as.matrix(read.csv(file="~/research/PDE-solvers/last-element-linear-interpolant.csv", header=F));
real.signal = Re(signal);

mass.matrix.fft = as.matrix(read.csv(file="~/research/PDE-solvers/mass-matrix-fft.csv",
                                 header=F));
stiffness.matrix.fft = as.matrix(read.csv(file="~/research/PDE-solvers/stiffness-matrix-fft.csv",
                                          header=F));

mass.matrix.linear = as.matrix(read.csv(file="~/research/PDE-solvers/mass-matrix-linear.csv",
                                 header=F));
stiffness.matrix.linear = as.matrix(read.csv(file="~/research/PDE-solvers/stiffness-matrix-linear.csv",
                                             header=F));

diff = mass.matrix.fft - mass.matrix.linear;
heatmap(diff);

diff = stiffness.matrix.fft - stiffness.matrix.linear;
mm = min(stiffness.matrix.fft - stiffness.matrix.linear);
heatmap((diff));
row.index = which(apply(X = diff - mm, MARGIN=1, FUN=min) == 0);
col.index = which( diff[row.index,] == mm );
print(diff[row.index, col.index]);

### SYSTEM MATRICES DX DX START ###
system.matrix.dx.dx.fft = as.matrix(read.csv(file="~/research/PDE-solvers/system-matrix-dx-dx-fft.csv", header=F));
system.matrix.dx.dx.linear = as.matrix(read.csv(file="~/research/PDE-solvers/system-matrix-dx-dx-linear.csv", header=F));

diff = abs(system.matrix.dx.dx.fft - system.matrix.dx.dx.linear);
heatmap(diff);
mm = max(diff);
row.index = which(apply(X = diff - mm, MARGIN=1, FUN=max) == 0);
col.index = which( diff[row.index,] == mm );
print(diff[row.index, col.index]);
### SYSTEM MATRICES DX DX END ###

### SYSTEM MATRICES DY DY START ###
system.matrix.dy.dy.fft = as.matrix(read.csv(file="~/research/PDE-solvers/system-matrix-dy-dy-fft.csv", header=F));
system.matrix.dy.dy.linear = as.matrix(read.csv(file="~/research/PDE-solvers/system-matrix-dy-dy-linear.csv", header=F));

diff = abs(system.matrix.dy.dy.fft - system.matrix.dy.dy.linear);
heatmap(diff);
mm = max(diff);
row.index = which(apply(X = diff - mm, MARGIN=1, FUN=max) == 0);
col.index = which( diff[row.index,] == mm );
print(diff[row.index, col.index]);
### SYSTEM MATRICES DY DY END ###


### SYSTEM MATRICES DX DY START ###
system.matrix.dx.dy.fft = as.matrix(read.csv(file="~/research/PDE-solvers/system-matrix-dx-dy-fft.csv", header=F));
system.matrix.dx.dy.linear = as.matrix(read.csv(file="~/research/PDE-solvers/system-matrix-dx-dy-linear.csv", header=F));

diff = abs(system.matrix.dx.dy.fft - system.matrix.dx.dy.linear);
heatmap(diff);
mm = max(diff);
row.index = which(apply(X = diff - mm, MARGIN=1, FUN=max) == 0);
col.index = which( diff[row.index,] == mm );
print(diff[row.index, col.index]);
### SYSTEM MATRICES DX DY END ###

IC.fft = as.matrix(read.csv(file="~/research/PDE-solvers/IC-matrix.csv", header=F));
IC.linear = as.matrix(read.csv(file="~/research/PDE-solvers/IC-matrix.csv", header=F));

max.abs.diffs = rep(NA, 100);
for (i in seq(0,99)) {
    element.linear.interpolant = as.matrix(read.csv(file=paste("~/research/PDE-solvers/orthonormal_element_", i, "-linear.csv", sep=""), header=F));

    element.fft.interpolant = as.matrix(read.csv(file=paste("~/research/PDE-solvers/orthonormal_element_", i, "-fft.csv", sep=""), header=F));
    diff = abs(element.fft.interpolant - element.linear.interpolant);
    max.abs.diffs[i+1] = max(diff);
}

element.linear.interpolant.i = as.matrix(read.csv(file=paste("~/research/PDE-solvers/element-i-linear.csv", sep=""), header=F));
element.linear.interpolant.j = as.matrix(read.csv(file=paste("~/research/PDE-solvers/element-j-linear.csv", sep=""), header=F));

n = dim(element.linear.interpolant.i)[1];

deriv.integral = 0;
for (i in seq(1,n-1)) {
    for (j in seq(1,n-1)) {

        deriv.integral = deriv.integral +
            ## ##
            (element.linear.interpolant.i[i+1,j] -
             element.linear.interpolant.i[i,j]) *
            ## ##
            (element.linear.interpolant.j[i+1,j] -
             element.linear.interpolant.j[i,j]);
        
    }
}

element.fft.interpolant.i = as.matrix(read.csv(file=paste("~/research/PDE-solvers/element-i-fft.csv", sep=""), header=F));
element.fft.interpolant.j = as.matrix(read.csv(file=paste("~/research/PDE-solvers/element-j-fft.csv", sep=""), header=F));

deriv.integral = 0;
for (i in seq(1,n-1)) {
    for (j in seq(1,n-1)) {

        deriv.integral = deriv.integral +
            ## ##
            (element.fft.interpolant.i[i+1,j] -
             element.fft.interpolant.i[i,j]) *
            ## ##
            (element.fft.interpolant.j[i+1,j] -
             element.fft.interpolant.j[i,j]);
        
    }
}

diff = abs(element.fft.interpolant - element.linear.interpolant);
max(diff);

