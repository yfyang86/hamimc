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



#####################
# usage
#####################
#ITER=10000
#x.initial=1:50
#buinin=5000
#dproposalfunction<-function(xt1,xt){# We could try other proposals
#  mvtnorm::dmvnorm(xt1,mean=xt,sigma=diag(rep(1,length(xt)))/20,log=T)
#}

#rproposalfunction<-function(xt){
#  mvtnorm::rmvnorm(n=1,mean=xt,sigma=diag(rep(1,length(xt)))/20)
#}

##target function
#A = matrix(rep(0.998,50^2),ncol=50)*outer(1:50,1:50,'*')
#diag(A)=(1:50)^2

#dtarget<-function(x,...){
#  p=length(x)
#  mvtnorm::dmvnorm(x,mean=rep(0,p),sigma=A,log=T)
#}

#   MH.framework(dproposalfunction,rproposalfunction,dtarget,A=A) ->  X.flow

#library(coda)

#xyplot(mcmc(X.flow[,c(2,45)]))
#par(mfrow=c(1,2))
#densplot(mcmc(X.flow[-(1:buinin),c(2,45)]))
#par(mfrow=c(1,1))


