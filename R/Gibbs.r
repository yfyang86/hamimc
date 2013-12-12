Gibbs.Norm50<-function(ITER=10000,p=50,chains=4){
A = matrix(rep(0.998,p^2),ncol=p)*outer(1:p,1:p,'*')
diag(A)=(1:p)^2
m=ITER
x0=rep(0,p)
n = length(x0)
Gibbs.lst<-list()
length(Gibbs.lst)=chains
tt=0
for (ii in 1:chains) {
  tic<-proc.time()
  store = matrix(0,n,m)
  store[,1] = x0
  for (i in 2:m){
    xt = store[,(i-1)]
    for (j in 1:n){
      mu = 0 + sum((A[j,-j]%*%solve(A[-j,-j]))*(xt[-j]))
      sigma = A[j,j]-A[j,-j]%*%solve(A[-j,-j])%*%A[-j,j]
      store[j,i] = rnorm(1,mu,sqrt(sigma))
    }
  }
  toc<-proc.time()
  t0=toc-tic
  Gibbs.lst[[ii]]=t(store)
  tt=t0+tt
}
return(Gibbs.lst);
}

