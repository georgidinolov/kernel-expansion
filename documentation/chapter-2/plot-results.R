args = commandArgs(trailingOnly = TRUE)
library("data.table");
library("ggplot2");
library("gridExtra");

mle.ochl = fread(args[1], header=TRUE)
rogers = fread(args[2], header=TRUE)
mle.classical = fread(args[3], header=TRUE)
save.path = args[4]

mle.ochl[, type := "MLE OCHL"]
rogers[, type := "Rogers"]
mle.classical[, type := "MLE Classical"]

all.est = rbind(mle.ochl, rogers, mle.classical)

rho.hist = ggplot(all.est, aes(y=..density.., x=rho, fill=type)) + geom_density(alpha=0.5)
sigma.x.hist = ggplot(all.est, aes(y=..density.., x=sigma_x, fill=type)) + geom_density(alpha=0.5)
sigma.y.hist = ggplot(all.est, aes(y=..density.., x=sigma_y, fill=type)) + geom_density(alpha=0.5)

all = grid.arrange(rho.hist, sigma.x.hist, sigma.y.hist, nrow=1)

ggsave(filename=paste0(save.path, "estimates-rho.pdf"), plot=rho.hist, width=6, height=6, units="in")
ggsave(filename=paste0(save.path, "estimates-sigma-x.pdf"), plot=sigma.x.hist, width=6, height=6, units="in")
ggsave(filename=paste0(save.path, "estimates-sigma-y.pdf"), plot=sigma.y.hist, width=6, height=6, units="in")
ggsave(filename=paste0(save.path, "estimates-all.pdf"), plot=all, width = 20, height=6, units = "in")



