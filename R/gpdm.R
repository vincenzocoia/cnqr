#' Mixture Generalized Pareto Distribution
#'
#' Distribution and Quantile functions for a mixture Generalized
#' Pareto Distribution. In particular, the tail of a "base" distribution
#' (above some threshold) is overridden with a GPD.
#'
#' @param x Vector of values to evaluate the distribution at.
#' @param mu,beta,xi Numeric; parameters of the GPD function (truncation,
#' scale, and extreme value index, respectively).
#' @param basecdf,baseqf Vectorized distribution function/quantile function
#' to use below the truncation \code{mu}.
#' @return Vector of evaluated distribution function, or quantiles.
#' @examples
#' cdf <- function(y) pgpdm(y, mu=2, xi=2, basecdf=pnorm)
#' qf <- function(tau) qgpdm(tau, mu=-2, xi=0, baseqf=qnorm)
#' curve(cdf, 0, 5)
#' curve(qf, 0, 1)
#'
#' set.seed(123)
#' x <- rnorm(100)
#' basecdf <- ecdf(x)
#' baseqf <- function(tau) quantile(x, tau)
#' cdf <- function(y) pgpdm(y, mu=0, beta=5, xi=-2, basecdf=basecdf)
#' qf <- function(tau) qgpdm(tau, mu=-1, xi=-1, baseqf=baseqf)
#' curve(cdf, -3, 5, )
#' curve(qf, 0, 1)
#' @rdname gpdm
#' @export
pgpdm <- function(x, mu = 0, beta = 1, xi = 0, basecdf) {
    n <- length(x)
    ## Split the data
    low <- x < mu
    nlow <- sum(low)
    xlow <- x[low]
    xhigh <- x[!low]
    ## Estimate probability of being low or high:
    plow <- basecdf(mu)
    phigh <- 1 - plow
    ## Evaluate GPD part of cdf (recycle xhigh)
    rightend <- if (xi < 0) mu - beta/xi else Inf
    toohigh <- xhigh > rightend
    xhigh[toohigh] <- 1
    xjustright <- xhigh[!toohigh]
    if (xi == 0) {
        xjustright <- 1 - exp(-(xjustright - mu)/beta)
    } else {
        xjustright <- 1 - (1 + xi * (xjustright - mu)/beta) ^ (-1/xi)
    }
    xhigh[!toohigh] <- xjustright * phigh + plow
    ## Evaluate basecdf part of cdf (recycle xlow)
    if (nlow == 0) {
        xlow <- numeric(0)
    } else {
        xlow <- basecdf(xlow)
    }
    ## Recycle variable x.
    x[low] <- xlow
    x[!low] <- xhigh
    return(x)
}


#' @rdname gpdm
#' @export
qgpdm <- function(p, mu = 0, beta = 1, xi = 0, baseqf) {
    n <- length(p)
    ## The "boundary quantile index" is the one that base = mu.
    ## First check that mu actually slices 'base'. Unfortunately, some
    ##  quantile functions, like qcondgum, throw an error at one or
    ##  both of the endpoints (0 or 1). Compromise by looking in [0.00001, 0.99999].
    plims <- c(0.00001, 0.99999)
    xlims <- baseqf(plims)
    ## Case 1: Threshold is above the base distribution.
    ##  Just evaluate at base distribution.
    if (mu > xlims[2]) return(baseqf(p))
    ## Case 2: Threshold is below the base distribution.
    ##  Just evaluate at GPD entirely.
    if (xi == 0) {
        gpdqf <- function(tau) mu + beta * log(1/(1-tau))
    } else {
        gpdqf <- function(tau) mu + beta * ((1/(1-tau))^xi - 1) / xi
    }
    if (mu < xlims[1]) return(gpdqf(p))
    ## There is a boundary quantile index. Find it with uniroot.
    f <- function(tau) baseqf(tau) - mu
    taucrit <- uniroot(f, c(0,1))$root
    low <- p <= taucrit
    ## Now just evaluate at the respective functions (recycle p)
    p[low] <- baseqf(p[low])
    p[!low] <- gpdqf((p[!low] - taucrit)/(1 - taucrit))
    return(p)
}
