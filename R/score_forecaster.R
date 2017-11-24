#' Score a Forecaster
#'
#' @param forecaster Function with first argument being the quantile index
#' \code{tau}, which
#' should accept a vector. Optionally, accepts a matrix of predictors
#' \code{x} as a
#' second argument. Returns either a vector of quantiles of length equal to
#' the length of \code{tau} if predictors are not accepted, or a matrix
#' of quantiles with columns corresponding to \code{tau}, and rows
#' corresponding to rows of \code{x}.
#' @param y Response vector
#' @param x Matrix of predictors. Optional if the forecaster doesn't
#' use predictors.
#' @param K Number of points to use when estimating the integral.
#' @param tau_c Lower quantile index of the forecasts.
#' @param g Vectorized univariate mapping function of the response and forecasts.
#' See formula in details.
#' @param wtau Function that accepts a vector of quantile indices and returns
#' a weight when assessing that quantile. See formula in details.
#' @param wx Function that accepts a vector of predictors (so, a row of \code{x})
#' and returns a weight value to adjust that observation's score by. Or,
#' a vector of weights already evaluated (no need to normalize them).
#' @param cmb Combine the scores across quantiles and observations?
#' \code{TRUE} by default. (See "Value" section for details of what's output)
#' @param yhat Optional matrix of predictions. It could be a time saver
#' because the forecasts don't need to be computed. However, be sure that
#' the quantiles that are forecast match up with \code{K} evenly spaced
#' quantile indices in \code{(tau_c, 1)}.
#' @return If \code{cmb=TRUE} (default), the mean score of the forecaster
#' is returned (i.e. a single numeric). Otherwise, a matrix of
#' (weight-adjusted) scores for each
#' observation (rows) and quantile index (columns) is returned, with columns
#' labelled by the quantile indices (rounded to four decimals).
#' @details
#' If y_i is the i'th observation,
#' tau_k = (2k-1)/(2K)*(1-tau_c) + tau_c is the k'th quantile index,
#' yhat_{ik} is the tau_k-quantile forecast for observation i, and
#' wx'_i = n wx_i / sum wx_i where wx_i = wx(x_i), then the scores
#' are computed as
#' S_{ik} = wx'_i wtau(tau_k) rho_{tau_k}(g(y) - g(yhat_{ik})).
#' These are returned in a matrix if \code{cmb=FALSE}, and are averaged
#' otherwise.
#' @export
score_forecaster <- function(forecaster, y, x, K = 10, tau_c = 0.9,
                             g = identity, wtau = function(tau) 1/tau/(1-tau),
                             wx = 1, cmb=TRUE, yhat=NULL) {
    n <- length(y)
    ## Get forecasts
    tau <- space_taus(K, tau_c)
    if (is.null(yhat)) {
        if (length(formals(forecaster)) == 1 || missing(x)) {
            ## Under this situation, it's probably the case that the forecaster
            ##  issues the same distribution all the time -- i.e. does not
            ##  depend on x. I checked this by seeing if `forecaster` only has one
            ##  argument, or if `x` wasn't specified.
            yhat <- forecaster(tau)
            yhat <- matrix(yhat, ncol=K, nrow=n, byrow=TRUE)
        } else {
            ## If the above situation isn't true, then the forecaster should
            ##  issue forecasts based on x.
            yhat <- forecaster(tau, x)
        }
    }
    ## Get weights across x:
    if (is.function(wx)) {
        ## Weights need to be evaluated over x.
        if (missing(x)) {  # x shouldn't be missing, but just in case...
            wx_eval <- 1
        } else {
            wx_eval <- apply(x, 1, function(row) wx(row))
        }
    } else {
        if (length(wx) == 1) wx_eval <- rep(wx, n) else wx_eval <- wx
    }
    ## Matrix of (x-weighted) scores
    s <- scoreq(y, yhat, tau, w=wx_eval, g=g, cmb=FALSE)
    ## Need to tau-weight them now.
    wtau_eval <- wtau(tau)
    s <- t(apply(s, 1, function(row) row * wtau_eval))
    ## Return.
    if (cmb) {
        return(mean(s))
    } else {
        colnames(s) <- paste0("tau=", round(tau, digits=4))
        return(s)
    }
}
