#' Make Quantile Predictions
#'
#' Use a \code{'cnqr'} object to make quantile predictions of a response
#' given predictors.
#'
#' @param obj Object of type \code{'cnqr'}. See \code{\link{cnqr}}.
#' @param dat Data matrix or data frame with raw data.
#' Or, a character vector specifying data that's already in \code{'obj'} --
#' for example, \code{"tr"} for the training data (always present), or
#' \code{"val"} for validation data. The data specified will be pooled.
#' @param tau Vector of quantile levels to predict. If \code{NULL}, uses
#' the quantile levels used in the regression of \code{obj}.
#' @param QY (Quantile) function that accepts a vector of quantile levels, and either:
#' (1) returns an equally long vector of corresponding quantiles, or
#' (2) returns a matrix of quantiles, with columns corresponding to the
#' quantile level input, and rows corresponding to the observations/rows in
#' \code{dat}. Alternatively, \code{NULL} if the quantile function is
#' already in \code{obj} (and you specify \code{dat} as a character).
#' @details
#' The response is taken to be the upper-right variable in the vine array of
#' \code{obj}, and the predictors are the variables below the response.
#' The vine is used to make quantile predictions of the conditional response.
#' @return
#' A matrix of quantiles, with rows corresponding to observations, and columns
#' corresponding to quantile levels.
#' @examples
#' data(egdat)
#' dat <- list(egdat[1:750, ], egdat[751:1000, ])
#' basevine <- subset(egvine, 1:4)
#' sc <- scorer(space_taus(10))
#'
#' fit <- cnqr(5:1, dat, basevine, sc=sc)
#'
#' head(predict(fit))
#' head(predict(fit, tau=c(0.1, 0.5, 0.9)))
#'
#' set.seed(442)
#' newdat <- rrvine(10, egvine)
#' predict(fit, dat=newdat)
#' @import copsupp
#' @export
predict.cnqr <- function(obj, dat="tr", tau=NULL, QY=NULL) {
    xyinfo <- xylink(obj)
    cdf <- obj$cdf
    if (is.null(tau)) tau <- obj$scorer$tau
    ## Case 1: Data are already in cnqr object.
    if (is.character(dat)) {
        if (!all(dat %in% names(obj$yhat)))
            stop(paste0("Can't find at least one of ", paste(dat, collapse=","),
                        " data in cnqr object."))
        ## Case 1a: requested tau matches that of the cnqr object.
        if (all.equal(tau, obj$scorer$tau) == TRUE) {
            dats <- lapply(dat, function(dat_) obj$yhat[[dat_]])
            return(do.call(rbind, dats))
        } else {
            ## Case 1b: requested tau is new. Don't need to compute uind though.
            uinds <- lapply(dat, function(dat_) obj$uind[[dat_]])
            uind <- do.call(rbind, uinds)
        }
    } else {
        ## Case 2: Data are not in cnqr object. Need to compute uind.
        y <- dat[, xyinfo$yvar]
        dat <- dat2udat(dat, cdf)
        uind <- pcondseq(dat, xyinfo$xord, obj)
    }
    ## Extract the marginal quantile function, if it's not specified.
    if (is.null(QY)) {
        ## Did the user input new data?
        if (!is.character(dat)) {
            ## Yes, so the QY in the cnqr object should apply to all observations.
            ##  (i.e. Y should be stationary).
            if (is.list(QY))  # A list QY suggests non-stationarity.
                stop("There's no QY in the cnqr object that applies to the new data.")
            QY <- obj$QY
        } else {
            ## No, the user did not input new data. Extract the corresponding QY.
            if (is.list(obj$QY)) {
                warning(paste("Non-stationary quantile functions have not been",
                              "implemented in full in this version of the package yet."))
                QY <- obj$QY[[dat]]
            } else {
                QY <- obj$QY
            }
        }
    }
    ## Compute yhat.
    QYgX(tau, uind, cops=xyinfo$cops, cpars=xyinfo$cpars, QY=QY)
}

