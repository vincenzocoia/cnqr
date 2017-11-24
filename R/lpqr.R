#' Initiate a Local Polynomial Quantile Smoothing
#'
#' Sets up a local polynomial quantile regression smoother.
#'
#' @param formula A formula indicating the response, and the covariates
#' to smooth over. Should be of the form "response ~ covariate_1 + ... +
#' covariate_p" (higher order and interaction terms are meaningless here).
#' @param data Data frame containing the data used for smoothing.
#' @param h Bandwidth(s). If a single scalar, specify in
#' units of mahalanobis distance.
#' If a vector (of length equal to the number of covariates, for 2 or more
#' covariates), each covariate is scaled by the bandwidths in \code{h}
#' (which should be specified in the same order as the \code{formula}).
#' @param kernel Character; type of kernel to use. One of
#' \code{"gaussian"}, or \code{"boxcar"}. See details for description.
#' @note There is a similarly named function in the package \code{quantreg},
#' called \code{lprq()}, that only allows for one covariate. \code{lpqr()}
#' is intended to extend this, using S3 generic functions. The nomenclature
#' of \code{"qr"} is used to follow that of \code{cnqr()}.
#' @details Since local polynomial regression needs to be done "on-the-spot",
#' the regression isn't actually done here. Instead, this function
#' sets up the regression by primarily building a weight function and a
#' data-centering function.
#'
#' Here are the kernel functions for distance \code{d} and bandwidth \code{h}:
#'
#' \describe{
#'      \item{"gaussian"}{\deqn{K(d,h) = \exp (-(d/h)^2/2)}}
#'      \item{"boxcar"}{\deqn{K(d,h) = I(d \le h)}}
#' }
#' @return Object of class \code{"lpqr"}, which is a list of the
#' following components.
#'
#' \describe{
#'      \item{extractx}{Function that takes a data frame containing
#'      at least columns of the covariates, and returns a
#'      data frame of the covariates indicated in \code{formula}.}
#'      \item{center}{Function that takes a covariate vector \code{x0} and
#'      returns a data frame of the response and covariates centered
#'      at \code{x0}.}
#'      \item{w}{Function that takes a covariate vector \code{x0} and
#'      returns a vector of normalized weights used for local polynomial
#'      regression, to be applied to the observations in \code{data}.}
#'      \item{data}{The original \code{data}.}
#'      \item{xnames}{Character vector of covariate names.}
#'      \item{yname}{Character name of the response.}
#' }
#' @examples ## See predict.lpqr() for examples.
#' @seealso \code{\link{predict.lpqr}}
#' @export
lpqr <- function(formula, data, h=0.5, kernel="gaussian") {
    ## Extract covariate columns from data (to obtain "xdata")
    extractx <- function(data) {
        ## Drop response from the formula.
        covformu <- formula
        covformu[[2]] <- NULL
        get_all_vars(covformu, data=data)
    }
    ## Apply to training data:
    xdata <- extractx(data)
    ydata <- get_all_vars(formula, data=data)[, 1, drop=FALSE]

    ## Make function to convert training data so that its covariates
    ##  are centered around a point x0.
    if (ncol(xdata) == 1) {
        center <- function(x0) {
            res <- as.data.frame(apply(xdata, 1, function(row) row - x0))
            names(res) <- names(xdata)
            return(as.data.frame(cbind(ydata, res)))
        }
    } else {
        center <- function(x0) {
            res <- t(apply(xdata, 1, function(row) row - x0))
            return(as.data.frame(cbind(ydata, res)))
        }
    }

    ## Get kernel function. Technically the kernel function should be a
    ##  function of "distance", but since it involves less computation
    ##  to compute distance^2, we'll make the kernel functions a
    ##  function of (vectorized) distance^2. It'll also accept a
    ##  bandwidth "h", which is not vectorized.
    if (kernel == "gaussian") {
        K <- function(tsq, h) exp(-tsq/2/h^2)
    } else if (kernel == "boxcar") {
        K <- function(tsq, h) as.integer(tsq <= h^2)
    } else {
        stop(paste0("Kernel '", kernel, "' not recognized."))
    }

    ## Make a weight function `w()` that takes a query point (vector) x0
    ##  and returns a vector of normalized weights for each observation
    ##  in the training data, to be used for regression about the point x0.
    if (length(h) == 1 || ncol(xdata) == 1) {
        ## Case 1: h is a single scalar. Use mahalanobis distance.
        Sigmainv <- solve(cov(xdata))
        w <- function(x0) {
            ## Get distances from the origin, from which weights are computed.
            dist_sq <- apply(xdata, 1, function(row)
                (row-x0) %*% Sigmainv %*% (row-x0))
            wp <- K(dist_sq, h)
            sum_wp <- sum(wp)
            if (sum_wp == 0) return(wp) else return(wp/sum_wp)
        }
    } else {
        ## Case 2: h has more than one scalar. Use non-tilted ellipse,
        ##  i.e. just re-scale the data and use euclidean distance.
        w <- function(x0) {
            dist_sq_adj <- apply(xdata, 1, function(row){
                ## Center and re-scale (by h)
                row_adj <- (row - x0) / h
                ## Get euclidean (square) distances from the origin
                sum(row_adj ^ 2)
            })
            wp <- K(dist_sq_adj, 1)
            wp / sum(wp)
        }
    }

    ## Output.
    res <- list(extractx=extractx,
                center=center,
                w=w,
                data=data,
                xnames=names(xdata),
                yname=names(ydata))
    class(res) <- "lpqr"
    return(res)
}

#' Evaluate Local Polynomial Quantile Surfaces
#'
#' Performs local polynomial quantile regression by evaluating the regression
#' surface(s) at query points (or "new data").
#'
#' @param object Object of class \code{"lpqr"} (see \code{\link{lpqr}}).
#' @param newdata A data frame containing at least covariate columns for which
#' to evaluate the regression surface(s) at.
#' @param p Non-negative integer; power of local polynomial
#' (which includes all interaction terms of that order). If you'd like more
#' flexibility, you can also specify
#' a formula here for the local regression.
#' @param tau Quantile indices indicating which quantile surfaces to evaluate.
#' @return A matrix, with columns corresponding to the quantile
#' indices \code{tau}, and columns corresponding to the query points in
#' \code{newdata}.
#' @examples
#' library(ggplot2)
#' set.seed(123)
#' ## Get data
#' x <- runif(1000)*pi
#' y <- abs(rnorm(1000))*sin(x)*100
#' z <- x + rexp(1000)*y
#' dat <- data.frame(x=x, y=y, z=z)
#' tau <- c(0.1, 0.5, 0.9)
#'
#' ## --- Univariate regression ---
#' setup <- lpqr(z ~ x, data=dat)
#' xgrid <- seq(0, pi, length.out=20)
#'
#' ## Example 1: Kernel smoothing, p=0
#' yhat <- predict(setup, newdata=data.frame(x=xgrid), p=0, tau=tau)
#' yhatlong <- reshape2::melt(yhat)
#' names(yhatlong) <- c("x", "tau", "z")
#' yhatlong$tau <- tau[yhatlong$tau]
#' yhatlong$x <- xgrid[yhatlong$x]
#' ggplot(dat, aes(x, z)) +
#'   geom_point() +
#'   geom_line(mapping=aes(group=tau),
#'             data=yhatlong, colour="red")
#'
#' ## Example 2: p=2
#' yhat <- predict(setup, newdata=data.frame(x=xgrid), p=2, tau=tau)
#' yhatlong <- reshape2::melt(yhat)
#' names(yhatlong) <- c("x", "tau", "z")
#' yhatlong$tau <- tau[yhatlong$tau]
#' yhatlong$x <- xgrid[yhatlong$x]
#' ggplot(dat, aes(x, z)) +
#'   geom_point() +
#'   geom_line(mapping=aes(group=tau),
#'             data=yhatlong, colour="red")
#'
#' ## --- Multivariate regression ---
#' setup <- lpqr(z ~ x + y, data=dat)
#' ygrid <- seq(min(y), max(y), length.out=20)
#' query <- expand.grid(x=xgrid, y=ygrid)
#'
#' ## Example 3: Kernel smoothing, p=0
#' yhat <- predict(setup, newdata=query, p=0, tau=tau)
#' head(yhat)
#'
#' ## Example 4: Self-defined local formula
#' yhat <- predict(setup, newdata=query,
#'                 p = z ~ x + I(x^2) + y, tau=tau)
#' head(yhat)
#' @import quantreg
#' @export
predict.lpqr <- function(object, newdata=object$data, p=1, tau=1:9/10) {
    if (length(tau) == 0) return(matrix(ncol=0, nrow=nrow(newdata)))
    ## Get things from the lpqr object:
    center <- object$center
    w <- object$w
    xnames <- object$xnames
    yname <- object$yname
    ## Get local regression formula using powers p:
    if (inherits(p, "formula")) {
        formu <- p
        p <- 0  # Just for the loop later.
    } else if (p == 0) {
        formu <- as.formula(paste(yname, "~ 1"))
    } else {
        formu <- as.formula(paste(yname, "~ poly(", paste(xnames, collapse=","),
                                  ", degree=", p, ", raw=TRUE)"))
    }
    ## Get the x data from newdata, which indicates the points for which
    ##  to query the regression surface. If anything, this will order the
    ##  covariates to match with the original formula specification in
    ##  the lpqr() call.
    xquery <- object$extractx(newdata)
    d <- ncol(xquery)
    ## Do the local polynomial regression over each query point:
    res <- apply(xquery, 1, function(row){
        ## Center the data.
        cdat <- center(row)
        ## Get weights for training data.
        w. <- w(row)
        ## Put weights on the data frame (rq() wants them together)
        cdat <- cbind(cdat, weights=w.)
        ## Remove observations with zero weights
        cdat <- subset(cdat, weights > 0)
        ## For d covariates with power p, we need at least
        ##  (p+1)/d * choose(p+d, d-1) observations. The rank of the
        ##  data matrix can't be too small either -- we'll take the easy
        ##  way out, and if rq() fails to work, we'll just reduce p
        ##  until it does.
        if (nrow(cdat) <= (p+1)/d * choose(p+d, d-1)) return(NA)
        thisp <- p
        while (thisp >= 0) {
            ## Do regression and grab the intercept
            ##  P.S. -- that's why we centered the data, so that we only
            ##  need to grab the intercept, as opposed to using the possibly more
            ##  computationally expensive "predict.rq()" function.
            res <- try(rq(formu, tau=tau, data=cdat, weights=weights))
            if (inherits(res, "try-error")) {
                ## rq() couldn't fit. Try reducing the power.
                thisp <- thisp - 1
            } else {
                res <- t(res$coefficients)[, 1]
                names(res) <- NULL
                return(res)
            }
        }
        return(rep(NA, length(tau)))
    })
    if (length(tau) == 1) return(matrix(res)) else return(t(res))
}
