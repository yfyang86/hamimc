hamimc
==============

Hamiltonian Monte Carlo in R
----------------------------------
This project is mainly about Hamiltonian Monte Carlo. 

Installation
---------

After install **devtools**, please run this line in R terminal:
```{r}
devtools::install_github('hamimc','yfyang86')
```
I've deleted all Cpp dependecies, so one doesn't need to compile a Rcpp version. Infuture I may offer a choice.

Introduction
----------------
It is a course related project so it should be simple and straightforward. Dr. Radford in Toronto University ha
s already contributed the R package [GRIMS](http://www.cs.utoronto.ca/~radford/GRIMS.html). Our goal is just to implement another simple approach.

> Detailed documents    
> Rccp/RcppArmadillo    
> RSTAN/RJAGS	

<img src="http://yuml.me/616642f3" >

hamimc.R
-----------------------
We transfered a MATLAB code into R code directly. Tuning method would be available:

> Greedy search tuning $\epsilon$'s   
> Simple version of No-U-tern tuning $L$(NOT finished yet, in R it is TOO SLOW)

Now the code works only for $N_2(\mu,\Sigma)$. There would also be a *easter egg* when release :)

[![A Demo Video](http://sweb.uky.edu/~yya234/images/hmcc.png)](https://github.com/yfyang86/hamimc/blob/master/hmcdemo.mp4)

A simple way of tuning $\epilon$ given leapfrog steps is also available.

ARS.R
-----------------------
There is an Adaptive Rejection Sampling (ARS) impelemented in our package. It is in pure R code (exception: using **rv** to sample piecewise exponential distribution).  For speed consieration, we suggest using

 - ars.new(**unuran**) http://statmath.wu.ac.at/unuran/    
 - ars(**ars**) http://cran.r-project.org/web/packages/ars/

Theses functions are more computational efficient. 

Additially, we offer an initialization function in our code:

```r
library(hamimc)
tg.func<-Vectorize(function(x){ # target function
  dnorm(x=x,mean=1001.24,log=T)
}
)
reject.sampling.intial(tg.func) -> initial ## Initialization!
re.sample<-reject.sampling(n=10000,tg.density=tg.func,graph=T,
                control=list(center=initial$mean,bound=initial$spread*5,step=0.123))
hist(re.sample$simu)
``` 

The initial function uses a line search to guess the (uni)mode (very roughly), then calculate the mean and variance.  In previouse example, the following code shows WHY the model fails without initialization.

```r
integrate(function(x)x*dnorm(x,mean=10000), -Inf, Inf)   ## NOT work
integrate(function(x)x*dnorm(x,mean=10000), 9990, Inf)   ## works
```

Notice if the variance is very large, or the distribution has no mean (e.g. [Cauchy Distribution](http://en.wikipedia.org/wiki/Cauchy_distribution)), this method will fail. We suggest using quantile information instead. For unimode distributions, this could be done via greedy search.  

Other ways of initialization are possible, *we are of no interest of comparison*. You can read my long post ![here](http://sweb.uky.edu/~yya234/stuff/cs/BGLM.html).

Bug Report and documents
------------------------------
One can send buf report using git-hub or through email. Besides, document is open too.


> Due to Xcode/Rcpp problem, hamimc could't work on R >=3.0.1 with Xcode>=5 on Mac>=10.9. But it works on R<=2.15.3 on any version of Mac OS.

> The document would be found on [an open Latex system](https://www.authorea.com/users/3481/articles/3578/_show_article).

Our team members are [Gao,Wei](http://mailto:g.w@uky.edu) | [Xie,Zhiheng](http://mailto:zhiheng.xie@uky.edu/) |[Yang,Yifan(me)](sweb.uky.edu/~yya234/). If one has found any bug, please contact |[Yang,Yifan(me)](yfian.yang@uky.edu).


<hr>
**Chinese version**
<hr>
这个项目只是一个课程的作业， 大约是讨论 Hamilton Monte Carlo 的一些性质。最为主要的是需要完成一个R包。
雷同的项目其实早就存在了, 比如多伦多大学的教授 Dr Radford M. Neal的项目主页:

http://www.cs.utoronto.ca/~radford/GRIMS.html

就已经做得很好了。 我们要做的只是一些基础实现。 我希望能够做到的


1. 文档更加丰富
2. 尽量使用Rcpp联合编程

