MH.framework<-function(dproposalfunction,rproposalfunction,dtarget,ITER=10000,p=50,x.initial=1:p,buinin=5000,...){# using template!
X.flow=matrix(0,nrow=ITER,ncol=p);
X.flow[1,]=x.initial

pb <- txtProgressBar(min = 0, max = ITER, style = 3)
for (iter in 1: (ITER-1)){
  X.t=X.flow[iter,]
  X.t1 <-  rproposalfunction(X.t);
  if (log(runif(1))< (
    dtarget(X.t1)+dproposalfunction(X.t,X.t1)-dtarget(X.t)-dproposalfunction(X.t1,X.t)
    )){
    X.flow[iter+1,]=X.t1
  }else{
    X.flow[iter+1,]=X.t
  }
  setTxtProgressBar(pb, iter)
}
close(pb)

return(X.flow);

}


MH.Normal<-function(ITER=10000,p=50,x.init=rep(0,50),
                    Sigma=matrix(rep(0.998,p^2),ncol=p)*outer(1:p,1:p,'*')+0.002*diag((1:p)^2),
                    chains=4,
                    buinin=5000
                    ){
                    
if (buinin>(ITER/2)) stop('Buin-in should be > 2 total iterations');
if ((dim(Sigma)[1])!=(dim(Sigma)[2])) stop('Sigma should be a possitive definet square matrix!(We will not detect >0 assumption)');
if(dim(Sigma)[1]!=p) stop('Dim(Sigma) should be p');
if (length(x.init)!=p) stop('Initial value should be of length p');
if (chains<1) stop('Number of chains should be >=1');
chains=ceiling(chains);


dproposalfunction<-function(xt1,xt){# We could try other proposals
  mvtnorm::dmvnorm(xt1,mean=xt,sigma=diag(rep(1,length(xt)))/20,log=T)
}

rproposalfunction<-function(xt){
  mvtnorm::rmvnorm(n=1,mean=xt,sigma=diag(rep(1,length(xt)))/20)
}

dtarget<-function(x,...){
  p=length(x)
  mvtnorm::dmvnorm(x,mean=rep(0,p),sigma=Sigma,log=T)
}

re<-list();
length(re)=chains;

for(i in 1:chains){
cat('\nChain ',i,' running...\n')
  MH.framework(dproposalfunction,rproposalfunction,dtarget,ITER,50,x.init,Sigma=Sigma) ->  re[[i]]
}

return (re);
}

#library(coda)

#xyplot(mcmc(X.flow[,c(2,45)]))
#par(mfrow=c(1,2))
#densplot(mcmc(X.flow[-(1:buinin),c(2,45)]))
#par(mfrow=c(1,1))


