hmc.1<-function(x,L=1000,Tau=0,epsilon=0){
  A<-matrix(rep(0.998,50^2),ncol=50)*outer(1:50,1:50,'*')
  diag(A)=(1:50)^2
  if (Tau<=0)Tau<-min(max(floor(50/sqrt(min(eigen(A)[[1]]))),100),1000)
  if(epsilon<=0)epsilon=max(sqrt(min(abs(eigen(A)[[1]]))),1E-6);
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
  pb <- txtProgressBar(min = 0, max = L, style = 3)
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
    setTxtProgressBar(pb, l) 
  }
  close(pb);
  return (x.trace);
}


hmc.tuningepsilon<-function(x,L=5000,Tau=100,epsilon=0.1^(3:5),cutoff=1.1){
print('By defauly I will use 3 chains to tuning HMC, chain length=5000, Num of steps=100')
print('One can increase the Num of steps if the results are not good.')
print('--------------------------------------------------------------------------------')
i=1
p=length(epsilon)
while (i <= p ){
cat('\nTry ',i,'-th epsilon = ',epsilon[i],'\n');
cat('Chain ',1,':\n');
  hmc.1(x,L=L,Tau=Tau,epsilon[i])->tmp1
cat('Chain ',2,':\n');
  hmc.1(x,L=L,Tau=Tau,epsilon[i])->tmp2
cat('Chain ',3,':\n');
  hmc.1(x,L=L,Tau=Tau,epsilon[i])->tmp3
  coda.C=coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(tmp1),
  coda::mcmc(tmp2),
  coda::mcmc(tmp3)
  )
  )
  
  if (max(coda.C$psrf)<=cutoff) break;
  i=i+1;
}
if (i==(p+1)) {
warning('The chosen epsilon is not suitable, please increase L/Tau to repeat the tuning.');
return (epsilon[p]);
}else{
return (epsilon[i]);
}
}
#library(snow)
HMC.cluster<- function(nu.chains = 4,iterartion=10000,startvalue = rep(0,50)) {
  cl <- makeCluster(rep("localhost", nu.chains), type = "SOCK")
  chain.cluster <- clusterCall(cl, hmc.1, x=startvalue,L=iterartion)
  eval(parse(text = paste("combinedchains=mcmc.list(", paste(paste("mcmc(chain.cluster[[", 
                                                                   1:nu.chains, "]])"), collapse = ","), ")")))
  stopCluster(cl)
  names(chain.cluster) <- paste0("Chain", 1:nu.chains, sep = "")
  return(list(mcmcchain = combinedchains, valuechain = chain.cluster))
}



