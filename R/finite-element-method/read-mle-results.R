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
plot(density(rhos.64), col="blue", lwd=2); lines(density(rhos.128), lwd=2); abline(v=0.6, lwd=2, col="red");


############# ALL #############################
par(mfrow=c(2,2))
plot(density(sigma.xs.16), type="l", lty=1, lwd=1, ylim=c(0,14)); lines(density(sigma.xs.32), lwd=2, lty=2, col = "red");
lines(density(sigma.xs.64), lwd=2, lty=3, col="blue"); lines(density(sigma.xs.128), lwd=2, lty=4);
abline(v=1.0,lwd=2,col="red");
## 
plot(density(sigma.ys.16), type="l", lty=1, lwd=1, ylim=c(0,14)); lines(density(sigma.ys.32), lwd=2, lty=2, col = "red");
lines(density(sigma.ys.64), lwd=2, lty=3, col="blue"); lines(density(sigma.ys.128), lwd=2, lty=4);
abline(v=1.0,lwd=2,col="red");
##
plot(density(rhos.16), type="l", lty=1, lwd=1, ylim=c(0,14)); lines(density(rhos.32), lwd=2, lty=2, col = "red");
lines(density(rhos.64), lwd=2, lty=3, col="blue"); lines(density(rhos.128), lwd=2, lty=4);
abline(v=0.6,lwd=2,col="red");

## coefs
## data: data lenght is a power of 2!
## xx
data = c(0.88447,-0.293022,-0.0624749,-0.0270319,-0.0150369,-0.00955344,-0.00659266,-0.00481345,-0.00366112,-0.00287222,-0.00230851,-0.00189175,-0.00157499,-0.00132861,-0.00113322,-0.00097568,-0.000846813,-0.000740074,-0.000650685,-0.000575091,-0.000510604,-0.000455161,-0.000407158,-0.000365332,-0.000328678,-0.000296386,-0.000267802,-0.000242389,-0.000219706,-0.000199384,-0.000181117,-0.000164647,-0.000149756,-0.000136258,-0.000123995,-0.00011283,-0.000102647,-9.33444e-05,-8.48331e-05,-7.70368e-05,-6.98883e-05,-6.33287e-05,-5.7306e-05,-5.17748e-05,-4.66946e-05,-4.20299e-05,-3.77489e-05,-3.38236e-05,-3.02291e-05,-2.69433e-05,-2.39465e-05,-2.12213e-05,-1.87521e-05,-1.65252e-05,-1.45286e-05,-1.27516e-05,-1.11849e-05,-9.82042e-06,-8.65125e-06,-7.67151e-06,-6.87636e-06,-6.26189e-06,-5.82509e-06,-5.56386e-06,-5.47692e-06,-3.4095e-06,-6.82312e-06,-1.0245e-05,-1.36793e-05,-1.71304e-05,-2.06024e-05,-2.40999e-05,-2.76273e-05,-3.11895e-05,-3.47913e-05,-3.84378e-05,-4.21343e-05,-4.58864e-05,-4.97e-05,-5.35814e-05,-5.75373e-05,-6.15746e-05,-6.57011e-05,-6.99248e-05,-7.42545e-05,-7.86997e-05,-8.32708e-05,-8.79788e-05,-9.2836e-05,-9.78559e-05,-0.000103053,-0.000108444,-0.000114046,-0.000119879,-0.000125966,-0.000132331,-0.000139003,-0.000146012,-0.000153394,-0.000161189,-0.000169443,-0.00017821,-0.000187551,-0.000197535,-0.000208246,-0.000219781,-0.000232255,-0.000245805,-0.000260596,-0.000276827,-0.000294743,-0.00031465,-0.000336928,-0.000362063,-0.000390685,-0.000423621,-0.000461989,-0.000507327,-0.000561826,-0.000628711,-0.000712949,-0.000822626,-0.000971909,-0.00118823,-0.00153303,-0.00217964,-0.00387575,-0.0168413);

x = as.matrix(read.csv(file="~/research/PDE-solvers/data/out.csv",header=F));
fcoefs = complex(real=x[,2],imaginary=x[,3]);
n = length(fcoefs);
xx = seq(0,1,length.out=100);
es = vector(mode="list", length=n);
for (i in seq(1,16)) {
    if (i==1) {
        es[[i]] = rep(1,length(xx)) * fcoefs[i];
    } else if (i==n/2) {
        es[[i]] = cos(n*pi*xx) * fcoefs[i];
    } else if (i < n/2) {
        es[[i]] = complex(real = cos(2*pi*xx*(i-1)),
                          imaginary = sin(2*pi*xx*(i-1))) *
            fcoefs[i];
    } else if (i > n/2) {
        es[[i]] = complex(real = cos(2*pi*xx*(i-n-1)),
                          imaginary = sin(2*pi*xx*(i-n-1))) *
            fcoefs[i];
    }
}

solution = complex(real=rep(0,length(xx)), imaginary=rep(0,length(xx)));
for (i in seq(1,n)) {
    solution = solution + es[[i]];
}
solution = solution/n;

plot(2*pi*xx, Re(solution));
lines(2*pi*xx, Re(solution));

######################
y.real = x[seq(1,n/2+1),2];
y.imag = x[seq(1,n/2+1),3];
solution = rep(0,length(xx));

for (i in seq(1,n/2+1)) {
    if (i==1) {
        solution =
            solution + rep(y[i], length(xx));
    } else if (i==n/2+1) {
        solution =
            solution + y[i]*cos(n*pi*xx);
    } else {
        solution =
            solution + 2*y.real[i]*cos((i-1)*2*pi*xx) -
            2*y.imag[i]*sin((i-1)*2*pi*xx);
    }
}
plot(xx,solution);

xx = seq(1,0,length.out=128);
coefs = rep(NA, length(data));

for (k in seq(0,length(data)-1)) {
    if (k==0) {
        coefs[k+1] = complex(real = data[k+1], imaginary=0);
    }
    
    else if (k == length(data)/2) {
        coefs[k+1] = complex(real = data[k+1], imaginary=0);
    }

    else if (k < length(data)/2) {
        coefs[k+1] = complex(real = data[k+1], imaginary=data[(length(data)-k)+1]);
    }

    else if (k > length(data)/2) {
        coefs[k+1] = complex(real = data[length(data)-k+1], imaginary=-data[k+1]);
    }
}

n = length(data);
solution = rep(0,length(xx));
for (k in seq(0,length(data)/2)) {
    if (k==0) {
        solution = solution + coefs[k+1];
    }

    else if (k==length(data)/2) {
        solution = solution + coefs[k+1] * cos(length(data)*pi*xx);
    }

    else if (k < length(data)/2) {
        solution =
            solution +
            coefs[k]*complex(real=2*pi*k*xx, imaginary=2*pi*k*xx);
    }

    else if (k > length(data)/2) {
        solution =
            solution +
            coefs[k]*complex(real=2*pi*(k-n)*xx, imaginary=2*pi*(k-n)*xx);
    }
}



elem=as.matrix(read.csv(file = "~/research/PDE-solvers/data/trig-test.csv", header=F));
plot(xx, elem[50,])
lines(xx, solution);


