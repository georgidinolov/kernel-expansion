args = commandArgs(trailingOnly = TRUE)
library("data.table");
library("ggplot2");
library("gridExtra");

mle.ochl = fread(args[1], header=TRUE)
rogers = fread(args[2], header=TRUE)
mle.classical = fread(args[3], header=TRUE)
save.path = args[4]
rho_data = as.double(args[5])
data_size = as.double(args[6])

mle.ochl[, type := "MLE OCHL"]
rogers[, type := "Rogers"]
mle.classical[, type := "MLE Classical"]

print(paste0("rho: ", mle.ochl[, mean( (rho-rho_data)^2 )],
             " ",
             mle.classical[, mean( (rho-rho_data)^2 )],
             " ",
             mle.ochl[, mean( (rho-rho_data)^2 )]/mle.classical[, mean( (rho-rho_data)^2 )]))

print(paste0("sigma_x: ", mle.ochl[, mean( (sigma_x-1.0)^2 )],
             " ",
             mle.classical[, mean( (sigma_x-1.0)^2 )],
             " ",
             mle.ochl[, mean( (sigma_x-1.0)^2 )]/mle.classical[, mean( (sigma_x-1.0)^2 )]))

print(paste0("sigma_y: ", mle.ochl[, mean( (sigma_y-1.0)^2 )],
             " ",
             mle.classical[, mean( (sigma_y-1.0)^2 )],
             " ",
             mle.ochl[, mean( (sigma_y-1.0)^2 )]/mle.classical[, mean( (sigma_y-1.0)^2 )]))


all.est = rbind(mle.ochl, rogers, mle.classical)

rho.hist = ggplot(all.est, aes(y=..density.., x=rho, fill=type)) +
    geom_density(alpha=0.5) +
    xlab(expression(rho)) +
    labs(title=paste0("m=", data_size)) +
    geom_vline(xintercept=rho_data) +
    theme(legend.position="none")

sigma.x.hist = ggplot(all.est, aes(y=..density.., x=sigma_x, fill=type)) +
    geom_density(alpha=0.5) +
    ylab(NULL) +
    xlab(expression(sigma[x])) +
    labs(title=paste0("m=", data_size)) +
    geom_vline(xintercept=1) +
    theme(legend.position="none")

sigma.y.hist = ggplot(all.est, aes(y=..density.., x=sigma_y, fill=type)) +
    geom_density(alpha=0.5) +
    xlab(expression(sigma[y])) +
    ylab(NULL) +
    labs(title=paste0("m=", data_size)) +
    geom_vline(xintercept=1) +
    theme(legend.position="none")

all = grid.arrange(rho.hist, sigma.x.hist, sigma.y.hist, nrow=1)

ggsave(filename=paste0(save.path, "estimates-rho.pdf"), plot=rho.hist, width=2, height=2, units="in")
ggsave(filename=paste0(save.path, "estimates-sigma-x.pdf"), plot=sigma.x.hist, width=2, height=2, units="in")
ggsave(filename=paste0(save.path, "estimates-sigma-y.pdf"), plot=sigma.y.hist, width=2, height=2, units="in")
ggsave(filename=paste0(save.path, "estimates-all.pdf"), plot=all, width = 8, height=6, units = "in")




