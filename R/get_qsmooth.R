#' Compute Smoothed Quantile Data Frame
#'
#' For x-y data, uses a moving window along x to estimate quantiles using a GPD
#' tail, using \code{\link{gpdqr}}.
#'
#' @param data Data frame containing the paired data
#' @param ycol,xcol Column numbers of the y and x variables.
#' @param npoints Number of points along the range of x to estimate quantiles at.
#' @param h bandwidth, in terms of number of standard deviations
#' @param tau Vector of quantile indices to regress over.
#'
#' @return Data frame with columns "x" (the grid of x values), "tau"
#' (quantile levels), and "y" (quantile estimate).
#' @export
get_qsmooth <- function(data, ycol=1, xcol=2, npoints=100, h=0.5, tau=space_taus(10)) {
    names(data)[c(ycol, xcol)] <- c("y", "x")
    ## Get grid of x values
    xvals <- range(data$x, na.rm = TRUE)
    xgrid <- seq(min(xvals), max(xvals), length.out = npoints)
    xnewdf <- data.frame(x=xgrid)
    ## Setup regression
    setup <- lpqr(y ~ x, data=data, h=h, kernel="boxcar")
    ## Make predictions
    pred <- gpdqr(setup, newdata=xnewdf, tau=tau)
    ## Melt the data frame:
    pred <- reshape2::melt(pred)
    names(pred) <- c("x", "tau", "y")
    pred[, 1] <- xgrid[pred[, 1]]
    pred$tau <- tau[pred$tau]
    pred
}
