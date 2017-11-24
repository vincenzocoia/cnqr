#' Make CNQR objective function
#'
#' Makes the CNQR objective function for you.
#'
#' @param y A vector of observed responses/outcomes.
#' @param yhat_fun A function of the parameter vector that returns a matrix
#' with \code{length(y)} rows of forecasts corresponding to \code{y}.
#' Columns should be the quantile forecasts corresponding to \code{taus}.
#' @param taus A vector specifying the quantile indices that are being
#' predicted.
#' @param g A vectorized non-decreasing function for which to transform
#' both response/outcome and forecast in the estimation.
#' @return Returns a function that takes parameter values and returns
#' the value of the objective function.
#' @note You don't really need to use this function, since it's almost
#' identical to the \code{\link{scoreq}} function.
#' @seealso \code{\link{space_taus}}
#' @export
cnqrobj <- function(y, yhat_fun, taus, g = identity) {
    function(theta) scoreq(y, yhat_fun(theta), tau = taus, g = g)
}

