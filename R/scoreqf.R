#' Score Quantile Functions
#'
#' Scores quantile function forecasts over the range \code{(tau_c, 1)}. See
#' the details to see how the scores are computed.
#'
#' @param y Vector of realizations
#' @param qf Function that accepts a vector of quantile indices as its only
#' argument, and returns a forecast matrix of quantiles
#' (columns correspond to quantile
#' indices, and rows correspond to observations in \code{y}). Alternatively,
#' it can output a single forecast of quantiles corresponding to the indices
#' that are input -- in this case, the same forecast is used for each
#' observation in \code{y}.
#' @param tau_c Lower quantile index
#' @param K Number of points to estimate the integral with.
#' @param w Vector of weights corresponding to the observations in \code{y}.
#' No need to normalize them.
#' @param G Function to determine the scoring rule. It should accept two
#' arguments: the first argument should
#' accept a numeric vector (on the scale of the response), and the second
#' argument should accept a quantile index within \eqn{(\tau_c,1)}.
#' @details Here's how the score for the \eqn{i}'th observation \eqn{y_i}
#' and quantile function \eqn{Qhat_i} is computed:
#' \deqn{w'_i / K
#'       \sum_k (\tau_k - I(y<yhat_{ik}))
#'              (G(y,\tau_k) - G(Qhat_i(\tau_k),\tau_k)),}
#' where \eqn{\tau_k = (2k-1)/(2K) (1-\tau_c) + \tau_c}. The
#' weights \eqn{w'} are "partially" normalized, and are computed by
#' \deqn{w'_i = n w_i / \sum_i w_i,}
#' where \eqn{n} is the number of observations. Note that, as \eqn{K}
#' approaches infinity, we obtain the integral over \eqn{(\tau_c, 1)} (if
#' it exists).
#' @return The scores are averaged and returned as a single numeric value.
#' @examples
#' set.seed(123)
#' u <- runif(100)
#'
#' scoreqf(u, identity)
#' scoreqf(u, identity, K=100)
#'
#' G <- function(y, tau) y/(1-tau)
#' scoreqf(u, identity, G=G)
#' scoreqf(u, identity, G=G, K=100)
#' @export
scoreqf <- function(y, qf, tau_c = 0.9, K = 10, w = 1, G = function(y, tau) y) {
    n <- length(y)
    ## Get vector of quantile indices
    tau <- (2*(1:K) - 1) / (2 * K) * (1 - tau_c) + tau_c
    ## Get matrix of forecasts
    yhat <- qf(tau)
    if (!is.matrix(yhat)) yhat <- matrix(yhat, ncol=K, nrow=n, byrow=TRUE)
    ## Get list of "g" functions.
    g <- lapply(tau, function(tau_) function(y) G(y, tau_))
    ## Compute scores
    scoreq(y, yhat, tau=tau, w=1, g=g, cmb=TRUE)
}
