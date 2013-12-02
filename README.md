hamimc
==============

Hamiltonian Monte Carlo in R
----------------------------------
This project is mainly about Hamiltonian Monte Carlo. It is a course related project so it should be simple and straightforward. Dr. Radford in Toronto University has already contributed the R package [GRIMS](http://www.cs.utoronto.ca/~radford/GRIMS.html). Our goal is just to implement another simple approach.

> Detailed documents    
> Rccp/RcppArmadillo    
> RSTAN/RJAGS	


hamimc.R
-----------------------
We transfered a MATLAB code into R code directly. Tuning method would be available:

> Greedy search tuning $\epsilon$'s   
> Simple version of No-U-tern tuning $L$

Now the code works only for $N_2(\mu,\Sigma)$. There would also be a *easter egg* when release :)

ARS.R
-----------------------
There is an Adaptive Rejection Sampling (ARS) impelemented in our package. It is in pure R code (exception: using **rv** to sample piecewise exponential distribution).  For speed consieration, we suggest using

 - ars.new(**unuran**) http://statmath.wu.ac.at/unuran/    
 - are(**ars**) http://cran.r-project.org/web/packages/ars/

Theses functions are more computational efficient.

Bug Report and documents
------------------------------
One can send buf report using git-hub or through email. Besides, document is open too.


> Due to Xcode/Rcpp problem, hamimc could't work on R >=3.0.1 with Xcode>=5 on Mac>=10.9. But it works on R<=2.15.3 on any version of Mac OS.

> The document would be found on [an open Latex system](https://www.authorea.com/users/3481/articles/3578/_show_article).

Our team member is [Gao,Wei](http://mailto:g.w@uky.edu) | [Xie,Zhiheng](http://mailto:zhiheng.xie@uky.edu/) |[Yang,Yifan(me)](sweb.uky.edu/~yya234/)

<hr>
**Chinese version**
<hr>
这个项目只是一个课程的作业， 大约是讨论 Hamilton Monte Carlo 的一些性质。最为主要的是需要完成一个R包。
雷同的项目其实早就存在了, 比如多伦多大学的教授 Dr Radford M. Neal的项目主页:

http://www.cs.utoronto.ca/~radford/GRIMS.html

就已经做得很好了。 我们要做的只是一些基础实现。 我希望能够做到的


1. 文档更加丰富
2. 尽量使用Rcpp联合编程

