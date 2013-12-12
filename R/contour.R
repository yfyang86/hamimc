hmcBiNorm.countour<-function(X.flow,buinin=500,bi=c(2,45)){
x1<-seq(-bi[1]*3, bi[1]*3, .2)
x2<-seq(-bi[2]*3, bi[2]*3, 1)
all<-expand.grid(x1, x2)
f.x<-mvtnorm::dmvnorm(all, mean = c(0,0), sigma =matrix(c(4,89.82,89.82,2025),2))
f.x2<-matrix(f.x, nrow=length(x1), ncol=length(x2))
par(pty = "s")
contour(x = x1, y = x2, z = f.x2, 
        xlab = expression(X[1]), ylab = expression(X[2]), 
        main = "Contour plot for bivariate normal distribution", 
        panel.first=grid(col="gray", lty="dotted"),cex=0.1)

points(X.flow[-(1:buinin),bi],type='l',col=3)
abline(h = 0, lty = 1, lwd = 2)  
abline(v = 0, lty = 1, lwd = 2) 
uu=summary(mcmc(X.flow[-(1:buinin),bi]))
mu.hat=uu$statistics[,1]
mu=c(0,0)
points(x = mu.hat[1], mu.hat[2], pch = 3, col = "darkred", lwd = 5)
points(x = mu[1], mu[2], pch = 4, col = 2, lwd = 2)
legend('bottomright', legend=c("mu.hat", "mu"), pch = c(3,4), col = c('darkred',2))
}



