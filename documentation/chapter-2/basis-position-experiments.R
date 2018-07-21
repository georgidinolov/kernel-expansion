library("mvtnorm")
sigma_x <- 0.10
sigma_y <- 0.10
rho <- 0.0;
theta <- pi/4
std_dev_factor = 1
epsilon = 0.1
small_t = 0.005
sigma_y_tilde = 0.75
rho_tilde = 0.95
par(mfcol=c(2,1), mar=c(2,2,2,2))

xy.samples <- t(rmvnorm(n=1e4,
                        mean=c(5.73238347443019e-01,
                               5.45249995886832e-03),
                      sigma=small_t*matrix(nrow=2,ncol=2,
                                           data=c(1^2,rho_tilde*1*sigma_y_tilde,
                                                  rho_tilde*1*sigma_y_tilde,sigma_y_tilde^2))))


corner.ll = c(-epsilon,-epsilon)
corner.lr = c(1+epsilon,-epsilon)
corner.ur = c(1+epsilon,1+epsilon)
corner.ul = c(-epsilon,1+epsilon)


Scale.mat <- diag(c(1/sigma_x, 1/sigma_y))
Rotation.mat <- matrix(nrow=2,ncol=2,
                       data=c(cos(theta),sin(theta),-sin(theta),cos(theta)))

corner.ll.xi <- Rotation.mat %*% Scale.mat %*% t(t(corner.ll))
corner.lr.xi <- Rotation.mat %*% Scale.mat %*% t(t(corner.lr))
corner.ul.xi <- Rotation.mat %*% Scale.mat %*% t(t(corner.ul))
corner.ur.xi <- Rotation.mat %*% Scale.mat %*% t(t(corner.ur))
center.xi <- Rotation.mat %*% Scale.mat %*% c(1/2, 1/2)

xieta.samples <- Rotation.mat %*% Scale.mat %*% xy.samples

plot(xieta.samples[1,], xieta.samples[2,], pch=20, xlim=c(-5,5), ylim=c(-1,12))
lines(c(corner.ll.xi[1],corner.lr.xi[1]),
      c(corner.ll.xi[2],corner.lr.xi[2]),
      col="red",lwd=2)
lines(c(corner.lr.xi[1],corner.ur.xi[1]),
      c(corner.lr.xi[2],corner.ur.xi[2]),
      col="red",lwd=2)
lines(c(corner.ur.xi[1],corner.ul.xi[1]),
      c(corner.ur.xi[2],corner.ul.xi[2]),
      col="red",lwd=2)
lines(c(corner.ul.xi[1],corner.ll.xi[1]),
      c(corner.ul.xi[2],corner.ll.xi[2]),
      col="red",lwd=2)

if (rho >= 0) {
    by.xi <- sqrt(1-rho)/sqrt(1+rho)*std_dev_factor
} else {
    by.xi <- std_dev_factor;
}
xi.nodes <- seq(center.xi[1], min(c(corner.ll.xi[1],corner.lr.xi[1],corner.ul.xi[1],corner.ur.xi[1])),
                by=-by.xi)
xi.nodes <- c(rev(xi.nodes),
              seq(center.xi[1]+by.xi,
                  max(c(corner.ll.xi[1],corner.lr.xi[1],corner.ul.xi[1],corner.ur.xi[1])),
                  by=by.xi))
if (rho >= 0) {
    by.eta <- std_dev_factor
} else {
    by.eta <- sqrt(1+rho)/sqrt(1-rho)*std_dev_factor
}
eta.nodes <- seq(center.xi[2], min(c(corner.ll.xi[2],corner.lr.xi[2],corner.ul.xi[2],corner.ur.xi[2])),
                by=-by.eta)
eta.nodes <- c(rev(eta.nodes),
               seq(center.xi[2]+by.eta,
                  max(c(corner.ll.xi[2],corner.lr.xi[2],corner.ul.xi[2],corner.ur.xi[2])),
                  by=by.eta))
## eta.nodes <- seq(min(c(corner.ll.xi[2],corner.lr.xi[2],corner.ul.xi[2],corner.ur.xi[2])),
##                  max(c(corner.ll.xi[2],corner.lr.xi[2],corner.ul.xi[2],corner.ur.xi[2])),
##                  by = by.eta)

xieta.nodes <- matrix(nrow=2,ncol=length(xi.nodes)*length(eta.nodes))
for (i in seq(1,length(xi.nodes))) {
    for (j in seq(1,length(eta.nodes))) {
        xieta.nodes[, length(eta.nodes)*(i-1) + j] =
            c(xi.nodes[i], eta.nodes[j])
    }
}
points(xieta.nodes[1,],
       xieta.nodes[2,],col=2,pch=20)
points(sqrt(2)/2 * 0.5 * (1/sigma_x - 1/sigma_y),
       sqrt(2)/2 * 0.5 * (1/sigma_x + 1/sigma_y),
       col="blue", pch=20)

xy.nodes <- solve(Scale.mat) %*% solve(Rotation.mat) %*% xieta.nodes

plot(xy.samples[1,], xy.samples[2,], 
     xlim=c(-1+0.5,1+0.5), ylim=c(-1+0.5,1+0.5), pch=20,col="grey")
lines(c(corner.ll[1],corner.lr[1]),
      c(corner.ll[2],corner.lr[2]),
      col="red",lwd=2)
lines(c(corner.lr[1],corner.ur[1]),
      c(corner.lr[2],corner.ur[2]),
      col="red",lwd=2)
lines(c(corner.ur[1],corner.ul[1]),
      c(corner.ur[2],corner.ul[2]),
      col="red",lwd=2)
lines(c(corner.ul[1],corner.ll[1]),
      c(corner.ul[2],corner.ll[2]),
      col="red",lwd=2)

lines(c(0,1),
      c(0,0),
      col="blue",lwd=2)
lines(c(1,1),
      c(0,1),
      col="blue",lwd=2)
lines(c(1,0),
      c(1,1),
      col="blue",lwd=2)
lines(c(0,0),
      c(1,0),
      col="blue",lwd=2)

points(xy.nodes[1,],
       xy.nodes[2,],col=2,pch=20)

points(5.73238347443019e-01,
       5.45249995886832e-03,
       col="blue", pch=20, lwd=3)
points(0.00000000000000e+00,
       8.78585549759250e-01,
       col="blue", pch=20, lwd=3)

points(5.45249995886832e-03,
       5.73238347443019e-01,
       col=5, pch=20, lwd=3)
points(8.78585549759250e-01,
       0.00000000000000e+00,
       col=5, pch=20, lwd=3)

basis.0 = c(-0.0303301,1.03033); points(basis.0[1], basis.0[2], col=4, pch=1);
basis.1 = c(-0.0303301,0.818198); points(basis.1[1], basis.1[2], col=4, pch=1);
basis.2 = c(0.0757359,0.924264); points(basis.2[1], basis.2[2], col=4, pch=1);
basis.3 = c(0.181802,1.03033); points(basis.3[1], basis.3[2], col=4, pch=1);
basis.4 = c(-0.0303301,0.606066); points(basis.4[1], basis.4[2], col=4, pch=1);
basis.5 = c(0.0757359,0.712132); points(basis.5[1], basis.5[2], col=4, pch=1);
basis.6 = c(0.181802,0.818198); points(basis.6[1], basis.6[2], col=4, pch=1);
basis.7 = c(0.287868,0.924264); points(basis.7[1], basis.7[2], col=4, pch=1);
basis.8 = c(0.393934,1.03033); points(basis.8[1], basis.8[2], col=4, pch=1);
basis.9 = c(-0.0303301,0.393934); points(basis.9[1], basis.9[2], col=4, pch=1);
basis.10 = c(0.0757359,0.5); points(basis.10[1], basis.10[2], col=4, pch=1);
basis.11 = c(0.181802,0.606066); points(basis.11[1], basis.11[2], col=4, pch=1);
basis.12 = c(0.287868,0.712132); points(basis.12[1], basis.12[2], col=4, pch=1);
basis.13 = c(0.393934,0.818198); points(basis.13[1], basis.13[2], col=4, pch=1);
basis.14 = c(0.5,0.924264); points(basis.14[1], basis.14[2], col=4, pch=1);
basis.15 = c(0.606066,1.03033); points(basis.15[1], basis.15[2], col=4, pch=1);
basis.16 = c(-0.0303301,0.181802); points(basis.16[1], basis.16[2], col=4, pch=1);
basis.17 = c(0.0757359,0.287868); points(basis.17[1], basis.17[2], col=4, pch=1);
basis.18 = c(0.181802,0.393934); points(basis.18[1], basis.18[2], col=4, pch=1);
basis.19 = c(0.287868,0.5); points(basis.19[1], basis.19[2], col=4, pch=1);
basis.20 = c(0.393934,0.606066); points(basis.20[1], basis.20[2], col=4, pch=1);
basis.21 = c(0.5,0.712132); points(basis.21[1], basis.21[2], col=4, pch=1);
basis.22 = c(0.606066,0.818198); points(basis.22[1], basis.22[2], col=4, pch=1);
basis.23 = c(0.712132,0.924264); points(basis.23[1], basis.23[2], col=4, pch=1);
basis.24 = c(0.818198,1.03033); points(basis.24[1], basis.24[2], col=4, pch=1);
basis.25 = c(-0.0303301,-0.0303301); points(basis.25[1], basis.25[2], col=4, pch=1);
basis.26 = c(0.0757359,0.0757359); points(basis.26[1], basis.26[2], col=4, pch=1);
basis.27 = c(0.181802,0.181802); points(basis.27[1], basis.27[2], col=4, pch=1);
basis.28 = c(0.287868,0.287868); points(basis.28[1], basis.28[2], col=4, pch=1);
basis.29 = c(0.393934,0.393934); points(basis.29[1], basis.29[2], col=4, pch=1);
basis.30 = c(0.5,0.5); points(basis.30[1], basis.30[2], col=4, pch=1);
basis.31 = c(0.606066,0.606066); points(basis.31[1], basis.31[2], col=4, pch=1);
basis.32 = c(0.712132,0.712132); points(basis.32[1], basis.32[2], col=4, pch=1);
basis.33 = c(0.818198,0.818198); points(basis.33[1], basis.33[2], col=4, pch=1);
basis.34 = c(0.924264,0.924264); points(basis.34[1], basis.34[2], col=4, pch=1);
basis.35 = c(1.03033,1.03033); points(basis.35[1], basis.35[2], col=4, pch=1);
basis.36 = c(0.181802,-0.0303301); points(basis.36[1], basis.36[2], col=4, pch=1);
basis.37 = c(0.287868,0.0757359); points(basis.37[1], basis.37[2], col=4, pch=1);
basis.38 = c(0.393934,0.181802); points(basis.38[1], basis.38[2], col=4, pch=1);
basis.39 = c(0.5,0.287868); points(basis.39[1], basis.39[2], col=4, pch=1);
basis.40 = c(0.606066,0.393934); points(basis.40[1], basis.40[2], col=4, pch=1);
basis.41 = c(0.712132,0.5); points(basis.41[1], basis.41[2], col=4, pch=1);
basis.42 = c(0.818198,0.606066); points(basis.42[1], basis.42[2], col=4, pch=1);
basis.43 = c(0.924264,0.712132); points(basis.43[1], basis.43[2], col=4, pch=1);
basis.44 = c(1.03033,0.818198); points(basis.44[1], basis.44[2], col=4, pch=1);
basis.45 = c(0.393934,-0.0303301); points(basis.45[1], basis.45[2], col=4, pch=1);
basis.46 = c(0.5,0.0757359); points(basis.46[1], basis.46[2], col=4, pch=1);
basis.47 = c(0.606066,0.181802); points(basis.47[1], basis.47[2], col=4, pch=1);
basis.48 = c(0.712132,0.287868); points(basis.48[1], basis.48[2], col=4, pch=1);
basis.49 = c(0.818198,0.393934); points(basis.49[1], basis.49[2], col=4, pch=1);
basis.50 = c(0.924264,0.5); points(basis.50[1], basis.50[2], col=4, pch=1);
basis.51 = c(1.03033,0.606066); points(basis.51[1], basis.51[2], col=4, pch=1);
basis.52 = c(0.606066,-0.0303301); points(basis.52[1], basis.52[2], col=4, pch=1);
basis.53 = c(0.712132,0.0757359); points(basis.53[1], basis.53[2], col=4, pch=1);
basis.54 = c(0.818198,0.181802); points(basis.54[1], basis.54[2], col=4, pch=1);
basis.55 = c(0.924264,0.287868); points(basis.55[1], basis.55[2], col=4, pch=1);
basis.56 = c(1.03033,0.393934); points(basis.56[1], basis.56[2], col=4, pch=1);
basis.57 = c(0.818198,-0.0303301); points(basis.57[1], basis.57[2], col=4, pch=1);
basis.58 = c(0.924264,0.0757359); points(basis.58[1], basis.58[2], col=4, pch=1);
basis.59 = c(1.03033,0.181802); points(basis.59[1], basis.59[2], col=4, pch=1);
basis.60 = c(1.03033,-0.0303301); points(basis.60[1], basis.60[2], col=4, pch=1);

## basis.0 = c(-0.0303301,1.03033); points(basis.0[1], basis.0[2], col=3, pch=20);
## basis.1 = c(0.0580583,0.977297); points(basis.1[1], basis.1[2], col=3, pch=20);
## basis.2 = c(0.234835,1.08336); points(basis.2[1], basis.2[2], col=3, pch=20);
## basis.3 = c(-0.0303301,0.818198); points(basis.3[1], basis.3[2], col=3, pch=20);
## basis.4 = c(0.146447,0.924264); points(basis.4[1], basis.4[2], col=3, pch=20);
## basis.5 = c(0.323223,1.03033); points(basis.5[1], basis.5[2], col=3, pch=20);
## basis.6 = c(0.0580583,0.765165); points(basis.6[1], basis.6[2], col=3, pch=20);
## basis.7 = c(0.234835,0.871231); points(basis.7[1], basis.7[2], col=3, pch=20);
## basis.8 = c(0.411612,0.977297); points(basis.8[1], basis.8[2], col=3, pch=20);
## basis.9 = c(0.588388,1.08336); points(basis.9[1], basis.9[2], col=3, pch=20);
## basis.10 = c(-0.0303301,0.606066); points(basis.10[1], basis.10[2], col=3, pch=20);
## basis.11 = c(0.146447,0.712132); points(basis.11[1], basis.11[2], col=3, pch=20);
## basis.12 = c(0.323223,0.818198); points(basis.12[1], basis.12[2], col=3, pch=20);
## basis.13 = c(0.5,0.924264); points(basis.13[1], basis.13[2], col=3, pch=20);
## basis.14 = c(0.676777,1.03033); points(basis.14[1], basis.14[2], col=3, pch=20);
## basis.15 = c(0.0580583,0.553033); points(basis.15[1], basis.15[2], col=3, pch=20);
## basis.16 = c(0.234835,0.659099); points(basis.16[1], basis.16[2], col=3, pch=20);
## basis.17 = c(0.411612,0.765165); points(basis.17[1], basis.17[2], col=3, pch=20);
## basis.18 = c(0.588388,0.871231); points(basis.18[1], basis.18[2], col=3, pch=20);
## basis.19 = c(0.765165,0.977297); points(basis.19[1], basis.19[2], col=3, pch=20);
## basis.20 = c(0.941942,1.08336); points(basis.20[1], basis.20[2], col=3, pch=20);
## basis.21 = c(-0.0303301,0.393934); points(basis.21[1], basis.21[2], col=3, pch=20);
## basis.22 = c(0.146447,0.5); points(basis.22[1], basis.22[2], col=3, pch=20);
## basis.23 = c(0.323223,0.606066); points(basis.23[1], basis.23[2], col=3, pch=20);
## basis.24 = c(0.5,0.712132); points(basis.24[1], basis.24[2], col=3, pch=20);
## basis.25 = c(0.676777,0.818198); points(basis.25[1], basis.25[2], col=3, pch=20);
## basis.26 = c(0.853553,0.924264); points(basis.26[1], basis.26[2], col=3, pch=20);
## basis.27 = c(1.03033,1.03033); points(basis.27[1], basis.27[2], col=3, pch=20);
## basis.28 = c(0.0580583,0.340901); points(basis.28[1], basis.28[2], col=3, pch=20);
## basis.29 = c(0.234835,0.446967); points(basis.29[1], basis.29[2], col=3, pch=20);
## basis.30 = c(0.411612,0.553033); points(basis.30[1], basis.30[2], col=3, pch=20);
## basis.31 = c(0.588388,0.659099); points(basis.31[1], basis.31[2], col=3, pch=20);
## basis.32 = c(0.765165,0.765165); points(basis.32[1], basis.32[2], col=3, pch=20);
## basis.33 = c(0.941942,0.871231); points(basis.33[1], basis.33[2], col=3, pch=20);
## basis.34 = c(-0.0303301,0.181802); points(basis.34[1], basis.34[2], col=3, pch=20);
## basis.35 = c(0.146447,0.287868); points(basis.35[1], basis.35[2], col=3, pch=20);
## basis.36 = c(0.323223,0.393934); points(basis.36[1], basis.36[2], col=3, pch=20);
## basis.37 = c(0.5,0.5); points(basis.37[1], basis.37[2], col=3, pch=20);
## basis.38 = c(0.676777,0.606066); points(basis.38[1], basis.38[2], col=3, pch=20);
## basis.39 = c(0.853553,0.712132); points(basis.39[1], basis.39[2], col=3, pch=20);
## basis.40 = c(1.03033,0.818198); points(basis.40[1], basis.40[2], col=3, pch=20);
## basis.41 = c(0.0580583,0.128769); points(basis.41[1], basis.41[2], col=3, pch=20);
## basis.42 = c(0.234835,0.234835); points(basis.42[1], basis.42[2], col=3, pch=20);
## basis.43 = c(0.411612,0.340901); points(basis.43[1], basis.43[2], col=3, pch=20);
## basis.44 = c(0.588388,0.446967); points(basis.44[1], basis.44[2], col=3, pch=20);
## basis.45 = c(0.765165,0.553033); points(basis.45[1], basis.45[2], col=3, pch=20);
## basis.46 = c(0.941942,0.659099); points(basis.46[1], basis.46[2], col=3, pch=20);
## basis.47 = c(-0.0303301,-0.0303301); points(basis.47[1], basis.47[2], col=3, pch=20);
## basis.48 = c(0.146447,0.0757359); points(basis.48[1], basis.48[2], col=3, pch=20);
## basis.49 = c(0.323223,0.181802); points(basis.49[1], basis.49[2], col=3, pch=20);
## basis.50 = c(0.5,0.287868); points(basis.50[1], basis.50[2], col=3, pch=20);
## basis.51 = c(0.676777,0.393934); points(basis.51[1], basis.51[2], col=3, pch=20);
## basis.52 = c(0.853553,0.5); points(basis.52[1], basis.52[2], col=3, pch=20);
## basis.53 = c(1.03033,0.606066); points(basis.53[1], basis.53[2], col=3, pch=20);
## basis.54 = c(0.0580583,-0.0833631); points(basis.54[1], basis.54[2], col=3, pch=20);
## basis.55 = c(0.234835,0.0227029); points(basis.55[1], basis.55[2], col=3, pch=20);
## basis.56 = c(0.411612,0.128769); points(basis.56[1], basis.56[2], col=3, pch=20);
## basis.57 = c(0.588388,0.234835); points(basis.57[1], basis.57[2], col=3, pch=20);
## basis.58 = c(0.765165,0.340901); points(basis.58[1], basis.58[2], col=3, pch=20);
## basis.59 = c(0.941942,0.446967); points(basis.59[1], basis.59[2], col=3, pch=20);
## basis.60 = c(0.323223,-0.0303301); points(basis.60[1], basis.60[2], col=3, pch=20);
## basis.61 = c(0.5,0.0757359); points(basis.61[1], basis.61[2], col=3, pch=20);
## basis.62 = c(0.676777,0.181802); points(basis.62[1], basis.62[2], col=3, pch=20);
## basis.63 = c(0.853553,0.287868); points(basis.63[1], basis.63[2], col=3, pch=20);
## basis.64 = c(1.03033,0.393934); points(basis.64[1], basis.64[2], col=3, pch=20);
## basis.65 = c(0.411612,-0.0833631); points(basis.65[1], basis.65[2], col=3, pch=20);
## basis.66 = c(0.588388,0.0227029); points(basis.66[1], basis.66[2], col=3, pch=20);
## basis.67 = c(0.765165,0.128769); points(basis.67[1], basis.67[2], col=3, pch=20);
## basis.68 = c(0.941942,0.234835); points(basis.68[1], basis.68[2], col=3, pch=20);
## basis.69 = c(0.676777,-0.0303301); points(basis.69[1], basis.69[2], col=3, pch=20);
## basis.70 = c(0.853553,0.0757359); points(basis.70[1], basis.70[2], col=3, pch=20);
## basis.71 = c(1.03033,0.181802); points(basis.71[1], basis.71[2], col=3, pch=20);
## basis.72 = c(0.765165,-0.0833631); points(basis.72[1], basis.72[2], col=3, pch=20);
## basis.73 = c(0.941942,0.0227029); points(basis.73[1], basis.73[2], col=3, pch=20);
## basis.74 = c(1.03033,-0.0303301); points(basis.74[1], basis.74[2], col=3, pch=20);
