A = matrix(rep(0.998,50^2),ncol=50)*outer(1:50,1:50,'*')
diag(A)=(1:50)^2

m=10000
x0=c(1:50)

n = length(x0)


store = matrix(0,n,m)
store[,1] = x0
for (i in 2:m){
xt = store[,(i-1)]
 for (j in 1:n){
  mu = 0 + A[j,-j]%*%solve(A[-j,-j])%*%(xt[-j])
  sigma = A[j,j]-A[j,-j]%*%solve(A[-j,-j])%*%A[-j,j]
  store[j,i] = rnorm(1,mu,sqrt(sigma))
 }
}
X11()
plot(store[1,],store[2,],type="l")



