#' Plots the Scatter Plot with Quantile Smoothing
#'
#' Uses splines to fit smoothed conditional quantile curves overtop of a
#' scatter plot. Intended for use with \code{ggpairs} in package
#' \code{GGally}.
#'
#' @param data Data frame of data
#' @param mapping Aesthetics to use
#' @param ... Other arguments to pass to \code{geom_point} for the underlying
#' scatterplot.
#' @param method The smoothing method to use. Either
#' \code{"lp"} (default) for local polynomial regression with \code{\link{lpqr}},
#' \code{"ss"} for spline regression with \code{quantreg}'s \code{rqss()}, or
#' \code{"gpd"} for kernel smoothing (i.e. local polynomial with degree 0)
#' with a GPD tail.
#' @param tau Vector of quantile indices to plot curves for.
#' @param npoints Number of points to evaluate each quantile curve at.
#' @param h Numeric smoothing parameter, passed to the \code{h} argument
#' in the function \code{\link{lpqr}} if \code{method="lp"}, or
#' passed to the \code{lambda}
#' argument in the function \code{qss()} in package \code{quantreg} if
#' \code{method="ss"}. For specifying a different smoothing parameter
#' depending on which variable is used as a covariate, make this a list
#' of such parameters -- see details for how to do this.
#' @param kernel Type of kernel to use. Passed to the argument \code{kernel}
#' in \code{\link{lpqr}}.
#' @details The quantile curves are estimated using the \code{rqss} function
#' in package \code{quantreg}. For response \code{y} and covariate \code{x},
#' the function call for a single \code{tau} is
#' \code{rqss(y ~ qss(x, lambda=h), tau=tau, data=data)}.
#'
#' There are two ways of specifying a list for \code{h}:
#'
#' \enumerate{
#'      \item A named list of smoothing parameters with names
#'      in \code{names(data)}, and possibly a component named \code{".default"}.
#'      The smoothing parameter will be selected according to the name
#'      of the \code{X} variable, or will use the \code{.default} if
#'      the name is not present in \code{h}.
#'      \item An unnamed list of length equal to \code{ncol(data)}. The
#'      smoothing parameter will be selected according to the column number
#'      of the \code{X} variable.
#' }
#'
#' This function should be used when calling \code{ggpairs} in package
#' \code{GGally}, under the "continuous" option of argument \code{upper}
#' and \code{lower} (for which you can specify the string \code{'qsmooth'}).
#' See examples.
#' @note The curves are fit to the entire scatterplot, even if groups are
#' somehow specified (such as through the \code{colour} aesthetic).
#' @references The code for this function is modelled in part after
#' \code{ggally_smooth} in package \code{GGally}.
#' @import quantreg
#' @examples
#' library(GGally)
#' library(ggplot2)
#' set.seed(123)
#' x <- runif(1000)*pi
#' y <- abs(rnorm(1000))*sin(x)*100
#' z <- x + rexp(1000)*y
#' dat <- data.frame(x=x, y=y, z=z, grp=factor(x<pi/2))
#'
#' ## Plain function call
#' ggpairs(dat, columns=1:3, lower=list(
#'   continuous="qsmooth"
#' ))
#'
#' ## Increase the smoothing parameter
#' ggpairs(dat, columns=1:3, lower=list(
#'   continuous=wrap("qsmooth", h=5)
#' ))
#'
#' ## Increase the smoothing parameter, just for Y:
#' ggpairs(dat, columns=1:3, lower=list(
#'   continuous=wrap("qsmooth", h=list(.default=1, y=5))
#' ))
#'
#' ## Add other specifications
#' ggpairs(dat, columns=1:3, lower=list(
#'   continuous=wrap("qsmooth", h=list(.default=1, y=50),
#'                   tau=16:19/20,
#'                   alpha=0.25,
#'                   method="ss"),
#'   mapping=aes(colour=grp)
#' ))
#' @export
ggally_qsmooth <- function(data, mapping, ..., method="lp", tau=1:9/10, npoints=100, h=1, kernel="gaussian") {
    ## Get variables
    xname <- as.character(mapping$x)
    yname <- as.character(mapping$y)
    ## Get h (we already have it if h is not a list).
    if (is.list(h)) {
        if (is.null(names(h))) {
            ## Unnamed h list.
            h <- h[[which(names(data)==xname)]]
        } else {
            ## Named lambda vector.
            trythis <- h[[xname]]
            if (is.null(trythis)) {
                h <- h[[".default"]]
            } else {
                h <- trythis
            }
        }
    }
    ## Was h specified correctly?
    if (is.null(h)) stop("`h` parameter is incorrectly specified.")

    ## Get grid of points to evaluate data at, but should be within the x range
    ##  according to ?predict.rqss documentation of the argument "newdata".
    xvals <- data[[xname]]
    xgrid <- seq(min(xvals), max(xvals), length.out = npoints)
    xnewdf <- data.frame(xgrid)
    names(xnewdf) <- xname

    ## Get predictions at the grid -- a data frame with columns
    ##  (xname), (yname), "tau" (not necessarily in that order).
    if (method == "lp") {
        ## Formula
        formu <- as.formula(paste(yname, "~", xname))
        ## Setup regression
        setup <- lpqr(formu, data=data, h=h, kernel=kernel)
        ## Make predictions
        pred <- predict(setup, newdata=xnewdf, tau=tau, p=2)
        ## Melt the data frame:
        pred <- reshape2::melt(pred)
        names(pred) <- c(xname, "tau", yname)
        pred[, 1] <- xgrid[pred[, 1]]
        pred$tau <- tau[pred$tau]
    } else if (method == "gpd") {
        ## Formula
        formu <- as.formula(paste(yname, "~", xname))
        ## Setup regression
        setup <- lpqr(formu, data=data, h=h, kernel=kernel)
        ## Make predictions
        pred <- gpdqr(setup, newdata=xnewdf, tau=tau)
        ## Melt the data frame:
        pred <- reshape2::melt(pred)
        names(pred) <- c(xname, "tau", yname)
        pred[, 1] <- xgrid[pred[, 1]]
        pred$tau <- tau[pred$tau]
    } else if (method == "ss") {
        ## Formula
        formu <- as.formula(paste(yname, "~ qss(", xname, ", lambda=", h,")"))
        ## Predictions
        pred <- lapply(tau, function(tau_){
            ## Do the regression
            reg <- rqss(formu, tau=tau_, data=data)
            cbind(xgrid, predict(reg, newdata=xnewdf))
        })
        pred <- as.data.frame(do.call(rbind, pred))
        names(pred) <- c(xname, yname)
        pred$tau <- rep(tau, each=npoints)
    } else {
        stop(paste("Unknown smoothing method:", method))
    }

    ## Output plot
    p <- ggplot(data, mapping) +
        geom_point(...)
    if (!is.null(mapping$color) || !is.null(mapping$colour)) {
        p <- p + geom_line(mapping=aes(group=tau), colour=I("black"), data=pred)
    }
    else {
        p <- p + geom_line(mapping=aes(group=tau), colour=I("red"), data=pred)
    }
    p
}
