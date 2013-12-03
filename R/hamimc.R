hmc.1<-function(x,gradE,findE,L,Tau,epsilon){
  A<-matrix(rep(0.998,50^2),ncol=50)*outer(1:50,1:50,'*')
  diag(A)=(1:50)^2
  A<-solve(A)
  findE<-function(x,...){
    (t(x)%*%A)%*%x/2
  }
  gradE<-function(x,...){
    A%*%x
  }

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
    
  }
  return (x.trace);
}

#re.hmc1<-hmc.1(x=(1:50),gradE,findE,iteration,500,0.05);
#xyplot(mcmc(re.hmc1[,1:5]))

