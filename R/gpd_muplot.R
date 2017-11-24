#' Plot of CNQR score over GPD threshold
#'
#' Estimates the mean CNQR scores of mixture generalized Pareto
#' distributions (GPD) across various choices of the threshold parameter. Here,
#' a "mixture GPD" is an empirical distribution having a GPD tail.
#'
#' @param ytr Vector of training data.
#' @param yval Vector of validation data.
#' @param sc Scoring rule to use, as in the output of \code{\link{scorer}}.
#' @param from,to Single numeric. Plotting bounds for the threshold parameter;
#' overrides the corresponding \code{qlims}. \code{NULL} to use \code{qlims}
#' instead.
#' @param qlims Vector of length 2 of quantile levels (between 0 and 1). The
#' corresponding empirical quantiles are used as plotting bounds for
#' the threshold parameter (unless the corresponding \code{from} and/or
#' \code{to} is specified).
#' @param ... other arguments to pass to the \code{\link{curve}} function.
#' @details
#' This function estimates the mean CNQR scores of distributions having an ecdf
#' base and generalized Pareto distribution (GPD) tail, for various
#' specifications of the GPD threshold parameter. For a selection of the
#' threshold parameter, a GPD is fit on the training data (\code{ytr})
#' above the threshold using maximum
#' likelihood (using the \code{gpd.fit} function of package \code{ismev}).
#' A mean score is estimated from the resulting distribution using
#' the validation data (\code{yval}).
#'
#' Smaller scores indicate a better fit of the distribution to the data.
#' @return The output of \code{\link{curve}}, with \code{x} values being
#' values of the GPD threshold, and \code{y} values being the estimated
#' mean CNQR score.
#' @note Minimizing the outputted curve is not as simple as just using
#' something like \code{\link{nlm}}, because the curve is typically
#' very jagged.
#' @examples
#' set.seed(342)
#' sc <- scorer(space_taus(10))
#' ytr <- c(-rexp(200), rexp(100))
#' yval <- c(-rexp(200), rexp(100))
#'
#' gpd_muplot(ytr, yval, sc)
#' gpd_muplot(ytr, yval, sc, n=200)
#' dat <- gpd_muplot(ytr, yval, sc, from=0.5, to=1.5)
#' dat
#' @import copsupp ismev
#' @export
gpd_muplot <- function(ytr, yval, sc, from=NULL, to=NULL, qlims=c(0.5, 0.95), ...) {
    ntr <- length(ytr)
    nval <- length(yval)
    K <- length(tau)
    ## ecdf of training data
    marg_tr <- marginal(ytr)
    Qtr <- marg_tr$qf
    ## Function that returns predictions for Y on the validation set.
    ##  Arguments are parameters of the GPD.
    qf <- function(loc, theta) {
        qhat <- qgpdm(tau, loc, beta = theta[1], xi = theta[2], baseqf = Qtr)
        matrix(qhat, nrow=nval, ncol=K, byrow=TRUE)
    }
    ## Objective function
    obj <- Vectorize(function(loc) {
        thetahat <- gpd.fit(ytr, loc, show=FALSE)$mle
        yhat <- qf(loc, thetahat)
        score_eval(yval, yhat, sc)
    })
    ## Limits
    if (is.null(from)) from <- quantile(c(ytr, yval), min(qlims))
    if (is.null(to)) to <- quantile(c(ytr, yval), max(qlims))
    curve(obj, from=from, to=to, ...)
}
