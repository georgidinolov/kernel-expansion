library("parallel")

## for simpson's rule, N needs to be even!
ezhh <- function(rho, N=5e6, upper.bound=100, display=FALSE) {
     alpha = asin(rho)
     gamma = (alpha + pi/2)/2
     etas = seq(0, upper.bound, length.out=N+1)

     args = cos(alpha)*cosh(etas*alpha)/sinh(etas*pi/2)*
		 tanh(etas*gamma)

     if (display) {
     	plot(etas, args, type="l")
     }

     indeces = seq(1,N/2)
     out = (sum(args[2*indeces-2 + 1], na.rm=TRUE) +
           sum(4*args[2*indeces-1 + 1], na.rm=TRUE) +
	   sum(args[2*indeces + 1], na.rm=TRUE))*(upper.bound/N)*(1/3)

     return(out)
}

ezhh.grid <- function(rho.lower.bound=-0.99,
		      rho.upper.bound= 0.99,
                      nn=50,
                      N=1e5) {
     rhos = seq(rho.lower.bound, rho.upper.bound, length.out=nn)
     ezhhs = rep(NA, length.out=nn)
     for (rho in rhos) {
     	      ezhhs[which(rhos==rho)] = ezhh(rho=rho, N=N)
	      print (paste0("On ", which(rhos==rho),
	      	    	    " out of ", nn))
     }
     out=NULL
     out$rhos = rhos
     out$ezhhs = ezhhs

     return (out)
}

phi.function <- function(rho) {
	     b.const = 2*log(2)-1

	     return (0.5*rho +
	     	     1.0/(2*(1-2*b.const))*(2*ezhh(rho) - 2*ezhh(-rho) - rho))

}

phi.function.grid <- function(rho.lower.bound=-0.999,
			      rho.upper.bound= 0.999,
                      	      nn=50) {
     rhos = seq(rho.lower.bound, rho.upper.bound, length.out=nn)
     phis = rep(NA, nn)
     for (rho in rhos) {
     	    phis[which(rhos==rho)] = phi.function(rho)
	    print (paste0("On ", which(rhos==rho),
	    	  	      " out of ", nn))
     }
     out = NULL
     out$rhos = rhos
     out$phis = phis
     return (out)
}



phi.function.inverse <- function(r, phi.grid) {
     indeces = sort(abs(phi.grid$phis - r), index.return=TRUE)$ix[c(1,2)]
     slope = (phi.grid$phis[indeces[2]] - phi.grid$phis[indeces[1]])/
     (phi.grid$rhos[indeces[2]] - phi.grid$rhos[indeces[1]])

     intercept = phi.grid$phis[indeces[2]] - slope*phi.grid$rhos[indeces[2]]

     inverse = (r - intercept)/slope
     return (inverse)
}

rogers.est <- function(open.1, open.2, close.1, close.2, high.1, high.2, low.1, low.2) {
	   b.const = 2*log(2) - 1
	   out = 0.5*close.1*close.2 +
	         0.5*1/(2*(1-2*b.const))*(high.1 + low.1 - close.1)*
		 (high.2 + low.2 - close.2)
	   return(mean(out))
}