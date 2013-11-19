ff<-function(knots,quest=0,endc=F){
  x=knots$x;
  y=knots$y;
  n=length(knots$x)
  #check sorted?
  if (sum(sort(x)==x)!=n){
    warning('Sorting knot X...');
    tmp.idx <- sort(knots$x, index.return =T)$ix
    x=knots$x[tmp.idx];
    y=knots$y[tmp.idx];
  }
  #Starting point/end point checking
  if(endc){if (y[1]!=0) {
    warning('y[1] must be 0')
   
    x0=(-20 *x[1] + 20*x[2] + x[2]* y[1] - x[1] *y[2] )/(y[1]-y[2])
    x=c(x0,x);
    n=n+1
    y=c(-20,y)
  }
          
    if (y[n]!=0){
             warning('y[n] must be 0')
             x1=(-20 *x[n-1] + 20*x[n] + x[n]* y[n-1] - x[n-1] *y[n])/(y[n-1]-y[n])
             x=c(x,x1);
             n=n+1
             y=c(y,-20)
           }
          
   
  }
 
  #if (sum(sort(x)==x)!=n){stop('Not enough points!!!');}
 
  ix=which(sort(c(quest,x))==quest)-1
  quest.y=(  quest* y[ix]+x[ix] *y[ix+1]-x[ix+1] *y[ix]- quest* y[ix+1])/( x[ix]-x[ix+1]  )
 
  return (list(kn=data.frame(x=x,y=y),pr=quest.y));
}
 
 
ff.vec<-Vectorize(
  function(x){
    ff(quest=x,knots=knots)$pr
  }
)
 
## inner
 
 
es.knots<-function(x,f=tg.density,...){
  x<-sort(x)
  inner.knots <- data.frame(x=x,y=f(x));
  n=length(x);
  y=f(x);
  xx=rep(0,n+1)
  yy=rep(0,n+1)
  xx[1]=x[1];yy[1]=y[1];
  xx[n+1]=x[n];yy[n+1]=y[n];
  linea2eq<-function(a,b){
   solve(cbind(a,-1),-b)
  }
  f.grad<-Vectorize(function(x){
    grad(func=f,x=x)
  })
  k=f.grad(x);
  for (i in 1:(n-1)){
   a=k[i:(i+1)]
   b=y[i:(i+1)]-a*x[i:(i+1)]
  tmp=linea2eq(a=a,b=b)
    xx[i+1]=tmp[1]
    yy[i+1]=tmp[2]
  }
  outer.knots <- data.frame(x=xx,y=yy);
 
  return(list(inner.knots =inner.knots ,  outer.knots=  outer.knots));
}
 
reject.sampling<-function(n,tg.density=tg.density,graph=T,method='ARS',detail=F,debug=F){ 
  x=seq(from=-15,to=15,by=.98124)
  if (graph){
  	      curve(tg.density,from=-4,to=4)
              points(x,tg.density(x),pch='*',col=4)
              }
  knots<-data.frame(x=x,y=tg.density(x));
  re<-ff(knots)
  if (graph){
    points(re$kn$x,re$kn$y,type='l',col=2)#inner-knots
    points( es.knots(x,f=tg.density)[[2]],pch='+',col=6)#outer-knots
  } 
  re2<-ff(es.knots(x,f=tg.density)[[2]])
  if (graph) points(re2$kn$x,re2$kn$y,type='l',col=3)#inner-knot
  es.knots(x,f=tg.density)[[2]]->knots
 
  #ff.vec<-Vectorize(
  #  function(x,...){
  #    exp(ff(quest=x,knots=knots)$pr)
  #  }
  #)
 
  ff.vec<-approxfun(knots$x,knots$y,method="linear")
  tt=integrate(function(x)exp(ff.vec(x)),lower=x[1],upper=max(x))$value
  tt2=integrate(function(x)exp(tg.density(x)),lower=x[1],upper=max(x))$value
  if (debug) cat('\nTarget density',tt2,'\n')
  ffstd.vec<-Vectorize(function(x,...){
    exp(ff.vec(x))/tt
  })
  nn=0;
  cc=0;
  rere<-NULL
  while(nn<=n){
    rvdens(1,FUN=ffstd.vec,range=range(x),unitprecision=10)->re# otherwise  pexp {msm} could be used
    U=runif(4000);
    re0 <- re[[1]][U<=(exp(tg.density(re[[1]]) )/tt2/exp(ff.vec(re[[1]])))]
    re0 <- re0[!is.na(re0)]
    rere<- c(re0,rere)
    nn=length(rere)
    cc=cc+1;
  }
  if (graph)legend("bottom",legend=paste('Acc Rate=',length(rere)/(4000*cc)));
  if (detail) print(paste('Acc Rate=',length(rere)/(4000*cc)));
  return(list(knots=knots,simu=rere[1:n]));
}
 
 
