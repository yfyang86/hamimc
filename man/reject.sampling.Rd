\name{reject.sampling}
\alias{reject.sampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ARS
}
\description{
ARS Methods wrapper. Please read the project main page for details.
}
\usage{
reject.sampling(n, tg.density = tg.density, graph = F, method = "ARS", detail = F, debug = F, control = list(center = 0, bound = 15, step = 0.1312), batch.mode = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{tg.density}{
%%     ~~Describe \code{tg.density} here~~
}
  \item{graph}{
%%     ~~Describe \code{graph} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{detail}{
%%     ~~Describe \code{detail} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{batch.mode}{
%%     ~~Describe \code{batch.mode} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, tg.density = tg.density, graph = F, method = "ARS", 
    detail = F, debug = F, control = list(center = 0, bound = 15, 
        step = 0.1312), batch.mode = T, ...) 
{
    centr.rs <- control$center
    bound.rs <- min(max(abs(control$bound), 10), 30)
    stepsize.rs <- min(abs(control$step), 0.98124)
    x = seq(from = -bound.rs, to = bound.rs, by = stepsize.rs) + 
        centr.rs
    if (graph) {
        curve(tg.density, from = centr.rs - bound.rs, to = centr.rs + 
            bound.rs)
        points(x, tg.density(x), pch = "*", col = 4)
    }
    knots <- data.frame(x = x, y = tg.density(x))
    knots.range <- which(knots$y > -30)
    knots.range = range(knots$x[knots.range])
    x = seq(from = knots.range[1], to = knots.range[2], length.out = 50)
    knots <- data.frame(x = x, y = tg.density(x))
    re <- ff(knots)
    if (graph) {
        points(re$kn$x, re$kn$y, type = "l", col = 2)
        points(es.knots(x, f = tg.density)[[2]], pch = "+", col = 6)
    }
    re2 <- ff(es.knots(x, f = tg.density)[[2]])
    if (graph) 
        points(re2$kn$x, re2$kn$y, type = "l", col = 3)
    knots <- es.knots(x, f = tg.density)[[2]]
    ff.vec <- approxfun(knots$x, knots$y, method = "linear")
    if (debug) 
        cat("\nTarget density", tt2, "\n")
    ffstd.vec <- Vectorize(function(x, ...) {
        exp(ff.vec(x))/tt
    })
    if (batch.mode) {
        tt = integrate(function(x) exp(ff.vec(x)), lower = x[1], 
            upper = max(x))$value
        tt2 = integrate(function(x) exp(tg.density(x)), lower = x[1], 
            upper = max(x))$value
        nn = 0
        cc = 0
        rere <- NULL
        while (nn <= n) {
            re <- rvdens(1, FUN = ffstd.vec, range = range(x), 
                unitprecision = 8)
            U = runif(4000)
            re0 <- re[[1]][U <= (exp(tg.density(re[[1]]))/tt2/exp(ff.vec(re[[1]])))]
            re0 <- re0[!is.na(re0)]
            rere <- c(re0, rere)
            nn = length(rere)
            cc = cc + 1
        }
    }
    else {
        n <- dim(knots)[1]
        a.s = (knots[2:n, 2] - knots[1:(n - 1), 2])/(knots[2:n, 
            1] - knots[1:(n - 1), 1])
        b.s = knots[2:n, 2] - a.s * knots[2:n, 1]
        knots.si = exp(b.s) * (exp(a.s * knots[2:n, 1]) - exp(a.s * 
            knots[1:(n - 1), 1]))/a.s
        knots.si = knots.si/sum(knots.si)
        knots.S = c(0, cumsum(knots.si), 1)
        tmpyyy <- runif(1)
        knots.loc = sum(tmpyyy > knots.S)
        return(list(knots = knots, simu = (log(tmpyyy - knots.S[knots.loc]) - 
            b.s[knots.loc])/a.s[knots.loc]))
    }
    if (graph) 
        legend("bottom", legend = paste("Acc Rate=", length(rere)/(4000 * 
            cc)))
    if (detail) 
        print(paste("Acc Rate=", length(rere)/(4000 * cc)))
    gc()
    return(list(knots = knots, simu = rere[1:n]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
