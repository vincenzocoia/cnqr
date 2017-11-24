#' Quantile Model Predictions
#'
#' @param object A fitted object with \code{\link{cnqr}}.
#' @param dat Either a data matrix to make predictions on (with the same columns
#' as the original fitting data), or one of the
#' two strings: \code{"val"} will make predictions on the validation set, and
#' \code{"tr"} will make predictions on the training set.
#' @param tau Vector of quantile indices to make predictions at.
#' @return A matrix of forecasts, with columns corresponding to quantile
#' indices \code{tau}, and rows corresponding to observations (rows) in
#' \code{dat}.
#' @examples
#' set.seed(123)
#' library(copsupp)
#'
#' ## Get data and a CNQR fit:
#' rv <- rvine(AtoG(CopulaModel::Dvinearray(5)), "frk", 3)
#' dattr <- rrvine(200, rv)
#' datval <- rrvine(200, rv)
#' fit <- cnqr(5:1, dattr, datval, rv, QY=identity)
#'
#' ## Predict on the validation set:
#' head(predict(fit))
#'
#' ## Predict on new data:
#' predict(fit, dat=rrvine(5, rv))
#' @export
predict.cnqr <- function(object, dat = "val", tau = object$tau) {
    ## Get sequential conditional preidtors.
    ucond <- NULL
    if (is.character(dat)) {
        if (dat == "val") {
            ucond <- object$ucondval
        }
        if (dat == "tr") {
            ucond <- object$ucondtr
        }
    } else {
        v <- object$xord
        rv <- object$basevine
        ucond <- pcondseq(dat, v, rv)
    }
    if (is.null(ucond))
        stop("'dat' input is inadequate.")
    ## Forecast.
    cops <- object$cops
    cpars <- object$cpars
    QY <- object$QY
    QYgX(tau, ucond, cops = cops, cpars = cpars, QY = QY)
}
