#' Calibration of a Forecast
#'
#' Computes the proportion of quantile forecasts that are exceeded, for each
#' quantile specified. Or, computes a histogram of probability integral
#' transformed observations. \code{calplot.cnqr} produces the calibration
#' plots from a fitted \code{'cnqr'} object, so that you don't need to
#' compute \code{yhat}, as you would with \code{calibration}.
#'
#' @param y Vector of outcomes/response.
#' @param yhat Matrix of forecasts. Each row is a forecast corresponding to
#' \code{y}; each column corresponds to a quantile index.
#' @param tau Vector of the quantile indices forecast.
#' @param obj Object of type \code{'cnqr'}. See \code{\link{cnqr}}.
#' @param dat Data matrix or data frame with raw data.
#' Or, a character vector specifying data that's already in \code{'obj'} --
#' for example, \code{"tr"} for the training data (always present), or
#' \code{"val"} for validation data. All data specified in the character
#' vector will be pooled.
#' @param w Vector of weights for each observation. Not necessarily normalized.
#' @param hist Logical; should values for a calibration histogram be computed
#' (\code{TRUE}), instead of for a survival function (\code{FALSE})? If
#' \code{TRUE}, then \code{tau} defines the left-endpoints of the bins.
#' @param disp_plot Logical; should a P-P plot be output?
#' @return A data frame.
#'
#' When \code{hist} is \code{TRUE}, the columns have names \code{"left"}
#' and \code{"right"}, for the left and right endpoints of the histogram bins
#' (which are quantile levels),
#' and \code{"hist"}, for the proportion of observations falling in the bin
#' (right-endpoint inclusive).
#'
#' When \code{hist} is \code{FALSE}, the columns have names \code{"ind"}
#' for the quantile levels, and \code{"exc"} for the proportion of
#' outcomes exceeding that quantile.
#' @examples
#' data(egdat)
#' dat <- list(egdat[1:750, ], egdat[751:1000, ])
#' basevine <- subset(egvine, 1:4)
#' sc <- scorer(space_taus(10))
#' fit <- cnqr(5:4, dat, basevine, sc=sc)
#'
#' calplot(fit)
#' calibration(dat[[1]][, 5], predict(fit), tau=fit$scorer$tau)
#' calplot(fit, tau=space_taus(100))
#' @rdname calibration
#' @export
calibration <- function(y, yhat, tau, w = 1, hist=FALSE, disp_plot = TRUE) {
    ntau <- length(tau)
    n <- length(y)
    if (ntau != ncol(yhat))
        stop("quantile indices 'tau' and columns of yhat are not equal in number.")
    if (n != nrow(yhat))
        stop("number of observations in y and yhat do not match.")
    ## Normalize weights
    wtot <- sum(w * rep(1, n))
    wn <- w / wtot
    exc <- apply(yhat, 2, function(col){
        exceeds <- as.integer(y > col)
        sum(wn * exceeds)
    })
    if (hist) {
        ## Histogram option
        ## 1. Get histogram heights
        his <- exc - c(exc[-1], 0)
        res <- data.frame(left=tau, right=c(tau[-1], 1), hist=his)
        ## 2. Plot?
        if (disp_plot) {
            p <- ggplot2::ggplot(res) +
                ggplot2::geom_rect(ggplot2::aes(xmin=left, ymin=0,
                                                xmax=right, ymax=hist)) +
                ggplot2::scale_x_continuous("Quantile Index") +
                ggplot2::scale_y_continuous("Histogram")
            print(p)
        }
        return(res)
    } else {
        ## Survival Function option
        res <- data.frame(ind=tau, exc=exc)
        if (disp_plot) {
            space <- diff(sort(tau))[1]
            p <- ggplot2::ggplot(res, ggplot2::aes(x=ind, y=exc)) +
                ggplot2::geom_abline(intercept = 1, slope = -1, linetype = "dotted") +
                ggplot2::geom_point() +
                ggplot2::scale_x_continuous("Quantile Index",
                                            limits = c(sort(tau)[1]-space, 1)) +
                ggplot2::scale_y_continuous("Proportion Exceeded",
                                            limits = c(0, max(1-(min(tau)-space), max(res$exc))))
            print(p)
        }
        return(res)
    }
}

#' Concepts for computation:
#'
#' The calibration plot is defined as the plot of P(Q(tau) > tau) over tau.
#' It can be thought of as the empirical survival function of F_hat(y_obs),
#' where F_hat is the forecasted cdf, and y_obs is the observation. We can
#' evaluate it at a bunch of taus by seeing how many observations exceed
#' the tau forecasts.
#'
#' By extension, a calibration _histogram_ can be defined as a histogram
#' of the F_hat(y_obs) values. Choosing bin boundaries (tau_c, tau_1, tau_2, ..., 1),
#' we can compute the height of each bar as the difference in survival function
#' between the two endpoints of the bin, which tells us how many observations
#' lie in that bin (divided by number of F_hat(y_obs)'s that exceed tau_c).

#' @rdname calibration
#' @export
calplot.cnqr <- function(obj, dat="tr", tau=NULL, QY=NULL, w=1, hist=FALSE, disp_plot=TRUE) {
    yhat <- predict(obj, dat=dat, tau=tau, QY=QY)
    if (is.character(dat)) {
        ys <- lapply(dat, function(dat_) obj$y[[dat_]])
        y <- c(ys, recursive=TRUE)
    } else {
        yvar <- xylink(obj)$yvar
        y <- dat[, yvar]
    }
    if (is.null(tau)) tau <- obj$scorer$tau
    calibration(y, yhat, tau,
                w = w,
                hist = hist,
                disp_plot = disp_plot)
}

#' @export
calplot <- function(...) UseMethod("calplot")
