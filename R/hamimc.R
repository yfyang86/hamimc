#
library(tcltk)
library(mvtnorm)
library(ks)
library(coda)

kdeplot <- function(data, fff = function(X)Hpi(X,pilot='amse',pre='scale')
	) {
    re.x <- kde(data, H = fff(data))

    contour(re.x$eval.points[[1]], re.x$eval.points[[2]], re.x$estimate)

}

A<-solve(matrix(c(1,0.998,0.998,1),ncol=2))

findE<-function(x,...){
	(t(x)%*%A)%*%x/2
}
gradE<-function(x,...){
	A%*%x
	}

hmc.1<-function(x,gradE,findE,L,Tau,epsilon){
	pb <- tkProgressBar(title = "progress bar", min = 0,
                    max = L, width = 300)
	g=gradE(x);
	E=findE(x);
	x.trace<-matrix(0,ncol=length(x),L+1);
	x.trace[1,]<-as.vector(x);
for (l in 1:L){# loop L times
	x.trace[l+1,] <- as.vector(x);
	p = rnorm ( length(x) ) ; # initial momentum is Normal(0,1)
    H = t(p) %*% p / 2 + E ; # evaluate H(x,p)
	xnew = x ; gnew = g ;
	for (tau in 1:Tau){ # make Tau `leapfrog' steps
		p = p - epsilon * gnew / 2 ; # make half-step in p
		xnew = xnew + epsilon * p ; # make step in x
		gnew = gradE ( xnew ) ; # find new gradient
		p = p - epsilon * gnew / 2 ; # make half-step in p
			}
		Enew = findE ( xnew ) ; # find new value of H
		Hnew = t(p) %*% p / 2 + Enew ;
		dH = Hnew - H ; # Decide whether to accept
		if ( dH < 0 ) {
			accept = 1 ;
			}else{
				if ( runif(1) < exp(-dH) ) {
					accept = 1 ;
					}else{
					accept = 0 ;
					}
			}

	if ( accept ){
		g = gnew ; x = xnew ; E = Enew ;
}
if (l%%1000==0){
	Sys.sleep(0.00001)
	setTkProgressBar(pb, l, label=paste('HMC: ', round(l/L*100, 0),
                                        "% done"))
	}	
}
 close(pb)
	return (x.trace);
}

#iteration=10000
#re.hmc1<-hmc.1(x=c(1,1),gradE,findE,iteration,19,0.055);

#xyplot(mcmc(re.hmc1))

