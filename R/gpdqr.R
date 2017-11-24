#' Evaluate Quantile Surfaces with GPD Tail
#'
#' Uses a kernel smoothing (or more generally, a local polynomial smoothing)
#' technique to evaluate quantiles for a set of
#' specified "query points" in the covariate space, with a Generalized
#' Pareto Distibution (GPD) in the tail. Note that the
#' scale and shape/EVI parameters are assumed to be constant
#' over the local covariate space.
#'
#' @param object Object of class \code{"lpqr"} (see \code{\link{lpqr}}).
#' Note that you must use \code{kernel="boxcar"} when building this object,
#' because \code{ismev}'s \code{gpd.fit()} doesn't allow for weights.
#' @param newdata Data frame having at least columns of the covariates, with
#' rows representing locations you wish to query.
#' @param tau_thresh Numeric; quantile index (between 0 and 1) for which
#' to take as the boundary between the GPD distribution and the nonparametric
#' "base" distribution.
#' @param tau Vector of numeric quantile indices (each between 0 and 1)
#' for which to evaluate.
#' @param p Power of local polynomial, used when fitting the threshold
#' parameter as well as quantiles with levels less than \code{tau_thresh}.
#' @details For quantile indices \code{tau} that are at most \code{tau_thresh},
#' local polynomial estimates are given (using
#' \code{\link{predict.lpqr}}). The same method is applied with
#' \code{tau=tau_thresh} to determine
#' the GPD threshold parameter, and then the scale and shape (EVI) parameters
#' of the GPD are estimated using weighted MLE, with weights determined by the
#' \code{"lpqr"} \code{object}.
#' @return A matrix with columns corresponding to the ordered quantile indices
#' \code{tau}, and rows corresponding to query points (rows) in \code{newdata}.
#' @examples
#' library(ggplot2)
#' set.seed(123)
#' ## Get data
#' x <- runif(1000)*pi
#' y <- abs(rnorm(1000))*sin(x)*100
#' z <- x + rexp(1000)*y
#' dat <- data.frame(x=x, y=y, z=z)
#' tau <- c(0.9, 0.95, 0.99, 0.999)
#'
#' ## --- Univariate regression ---
#' setup <- lpqr(z ~ x, data=dat, kernel="boxcar")
#' xgrid <- seq(0, pi, length.out=20)
#'
#' yhat <- gpdqr(setup, newdata=data.frame(x=xgrid), tau=tau)
#' yhatlong <- reshape2::melt(yhat)
#' names(yhatlong) <- c("x", "tau", "z")
#' yhatlong$tau <- tau[yhatlong$tau]
#' yhatlong$x <- xgrid[yhatlong$x]
#' ggplot(dat, aes(x, z)) +
#'     geom_point() +
#'     geom_line(mapping=aes(group=tau),
#'               data=yhatlong, colour="red")
#'
#' yhat <- gpdqr(setup, newdata=data.frame(x=xgrid), tau=tau)
#' yhatlong <- reshape2::melt(yhat)
#' names(yhatlong) <- c("x", "tau", "z")
#' yhatlong$tau <- tau[yhatlong$tau]
#' yhatlong$x <- xgrid[yhatlong$x]
#' ggplot(dat, aes(x, z)) +
#'     geom_point() +
#'     geom_line(mapping=aes(group=tau),
#'               data=yhatlong, colour="red")
#'
#' ## --- Multivariate regression ---
#' setup <- lpqr(z ~ x + y, data=dat)
#' ygrid <- seq(min(y), max(y), length.out=20)
#' query <- expand.grid(x=xgrid, y=ygrid)
#'
#' yhat <- gpdqr(setup, newdata=query, tau=tau)
#' head(yhat)
#' @import ismev
#' @seealso \code{\link{predict.lpqr}}
#' @export
gpdqr <- function(object, newdata, tau_thresh=0.9, tau=space_taus(10, tau_c=tau_thresh), p=2) {
    tau <- sort(tau)
    n <- nrow(newdata)
    ## Split quantiles into "lower" (empirical) and "upper" (gpd)
    tau_lower <- tau[tau <= tau_thresh]
    nlower <- length(tau_lower)
    tau_upper <- tau[tau > tau_thresh]
    nupper <- length(tau_upper)
    ## Get empirical estimates for the lower quantiles.
    yhat_lower <- predict(object, newdata=newdata, p=p, tau=tau_lower)
    ## Are any tau's requested above tau_thresh? If not, we're done.
    if (length(tau_upper) == 0) return(yhat_lower)
    ## --- Get gpd estimates for the upper quantiles ---
    ## Get thresholds, which we're taking to be the empirical tau_thresh-quantiles.
    if (tau_thresh %in% tau_lower[nlower]) {  # I'd rather use "==", but sometimes length(rhs)=0.
        ## They've already been computed -- just extract them.
        mu <- yhat_lower[, nlower]
    } else {
        ## They need computation.
        mu <- predict(object, newdata=newdata, p=p, tau=tau_thresh)[, 1]
    }
    ## Go through point-by-point and estimate quantiles:
    xdata <- as.matrix(object$extractx(newdata))
    w <- object$w
    y <- object$data[[object$yname]]
    yhat_upper <- matrix(nrow=n, ncol=nupper)
    for (i in 1:n) {
        ## Get query point
        x0 <- xdata[i, ]
        ## Get weights
        thisw <- w(x0)
        ## Subset data with non-zero weight -- they'll be weighted equally.
        suby <- y[thisw > 0]
        if (length(suby) == length(y) & i==1)
            warning(paste("All data are being used for GPD estimation -- did you",
                          "forget to use the 'boxcar' kernel?"))
        ## Get lower threshold
        thismu <- mu[i]
        ## Which Y values lie in the GPD region?
        subsamp <- suby >= thismu
        subsuby <- suby[subsamp]
        # thisw <- thisw[subsamp]
        # thisw <- thisw / sum(thisw)
        # ## Get objective function (for MLE) (theta=(sigma, xi))
        # ell <- function(theta) {
        #     if (theta[1] < 0) return(9999999999)
        #     if (xi < 0) {
        #         right_endpoint <- thismu - theta[1]/theta[2]
        #         if (any(thisy >= right_endpoint)) return(9999999999)
        #     }
        #     logf <- -log(theta[1]) - (1/theta[2] + 1) *
        #         log(1 + theta[2]*(thisy - thismu)/theta[1])
        #     res <- -sum(thisw * logf)
        #     res
        # }
        ## Estimate, but only if there's enough data to estimate the lower
        ##  threshold in this area. Otherwise, keep the entries for this
        ##  quantile NA and move on.
        if (!is.na(thismu) & length(unique(subsuby)) > 2) {
            # thetahat <- nlm(ell, c(1,1))$estimate
            thetahat <- gpd.fit(subsuby, thismu, show=FALSE)$mle
            sigma <- thetahat[1]
            xi <- thetahat[2]
            ## Extract quantiles (convert "absolute" quantile index to "within-gpd"
            ##  quantile index first).
            within_tau <- (tau_upper - tau_thresh) / (1 - tau_thresh)
            yhat_upper[i, ] <- thismu + sigma * ((1-within_tau)^(-xi) - 1) / xi
        }
    }
    return(cbind(yhat_lower, yhat_upper))
}
