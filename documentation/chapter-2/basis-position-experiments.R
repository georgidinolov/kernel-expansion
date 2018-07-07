library("mvtnorm")
sigma_x <- 0.2
sigma_y <- 0.2
rho <- 0.90
theta <- pi/4
std_dev_factor = 0.5
par(mfrow=c(2,1), mar=c(2,2,2,2))

xy.samples <- t(rmvnorm(n=1e4,
                      mean=c(0.5,0.5),
                      sigma=matrix(nrow=2,ncol=2,
                                   data=c(sigma_x^2,rho*sigma_x*sigma_y,
                                          rho*sigma_x*sigma_y,sigma_y^2))))


corner.ll = c(0,0)
corner.lr = c(1,0)
corner.ur = c(1,1)
corner.ul = c(0,1)


Scale.mat <- diag(c(1/sigma_x, 1/sigma_y))
Rotation.mat <- matrix(nrow=2,ncol=2,
                       data=c(cos(theta),sin(theta),-sin(theta),cos(theta)))

corner.ll.xi <- Rotation.mat %*% Scale.mat %*% t(t(corner.ll))
corner.lr.xi <- Rotation.mat %*% Scale.mat %*% t(t(corner.lr))
corner.ul.xi <- Rotation.mat %*% Scale.mat %*% t(t(corner.ul))
corner.ur.xi <- Rotation.mat %*% Scale.mat %*% t(t(corner.ur))

xieta.samples <- Rotation.mat %*% Scale.mat %*% xy.samples

plot(xieta.samples[1,], xieta.samples[2,])
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

by.xi <- sqrt(1-rho)/sqrt(1+rho)*std_dev_factor
xi.nodes <- seq(-sqrt(2)/2*(1/sigma_y), sqrt(2)/2*(1/sigma_x),
                by = by.xi)

by.eta <- std_dev_factor
eta.nodes <- seq(0, sqrt(2)/2*(1/sigma_x + 1/sigma_y),
                 by = by.eta)

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

plot(xy.samples[1,], xy.samples[2,], xlim=c(-1+0.5,1+0.5), ylim=c(-1+0.5,1+0.5))
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
points(xy.nodes[1,],
       xy.nodes[2,],col=2,pch=20)

basis.0 = c(0.0303301,1.03033); points(basis.0[1], basis.0[2], col=3, pch=20);
basis.1 = c(0.0303301,0.676777); points(basis.1[1], basis.1[2], col=3, pch=20);
basis.2 = c(0.207107,0.853553); points(basis.2[1], basis.2[2], col=3, pch=20);
basis.3 = c(0.383883,1.03033); points(basis.3[1], basis.3[2], col=3, pch=20);
basis.4 = c(0.0303301,0.323223); points(basis.4[1], basis.4[2], col=3, pch=20);
basis.5 = c(0.207107,0.5); points(basis.5[1], basis.5[2], col=3, pch=20);
basis.6 = c(0.383883,0.676777); points(basis.6[1], basis.6[2], col=3, pch=20);
basis.7 = c(0.56066,0.853553); points(basis.7[1], basis.7[2], col=3, pch=20);
basis.8 = c(0.737437,1.03033); points(basis.8[1], basis.8[2], col=3, pch=20);
basis.9 = c(0.0303301,-0.0303301); points(basis.9[1], basis.9[2], col=3, pch=20);
basis.10 = c(0.207107,0.146447); points(basis.10[1], basis.10[2], col=3, pch=20);
basis.11 = c(0.383883,0.323223); points(basis.11[1], basis.11[2], col=3, pch=20);
basis.12 = c(0.56066,0.5); points(basis.12[1], basis.12[2], col=3, pch=20);
basis.13 = c(0.737437,0.676777); points(basis.13[1], basis.13[2], col=3, pch=20);
basis.14 = c(0.914214,0.853553); points(basis.14[1], basis.14[2], col=3, pch=20);
basis.15 = c(0.383883,-0.0303301); points(basis.15[1], basis.15[2], col=3, pch=20);
basis.16 = c(0.56066,0.146447); points(basis.16[1], basis.16[2], col=3, pch=20);
basis.17 = c(0.737437,0.323223); points(basis.17[1], basis.17[2], col=3, pch=20);
basis.18 = c(0.914214,0.5); points(basis.18[1], basis.18[2], col=3, pch=20);
basis.19 = c(1.09099,0.676777); points(basis.19[1], basis.19[2], col=3, pch=20);
basis.20 = c(0.737437,-0.0303301); points(basis.20[1], basis.20[2], col=3, pch=20);
basis.21 = c(0.914214,0.146447); points(basis.21[1], basis.21[2], col=3, pch=20);
basis.22 = c(1.09099,0.323223); points(basis.22[1], basis.22[2], col=3, pch=20);
