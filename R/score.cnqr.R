#' Score CNQR Models
#'
#' Assesses \code{"cnqr"} fits using a Proper Scoring Rule.
#'
#' @param obj Object of type \code{'cnqr'}. See \code{\link{cnqr}}.
#' @param dat Data matrix or data frame with raw data.
#' Or, a character vector specifying data that's already in \code{'obj'} --
#' for example, \code{"tr"} for the training data (always present), or
#' \code{"val"} for validation data. The outputted score is of the
#' pooled data in the \code{dat} vector.
#' @param sc Scoring rule to use for the regression, as in the output
#' of \code{\link{scorer}}.
#' @param QY (Quantile) function that accepts a vector of quantile levels, and
#' returns an equally long vector of corresponding quantiles (only a
#' stationary quantile function is accepted at the moment)
#' Alternatively, \code{NULL} if the quantile function is
#' already in \code{obj} (and you specify \code{dat} as a character).
#' @param w Vector of weights corresponding to the observations in
#' \code{yhat}. No need to normalize them. See details to see exactly how these
#' weights are used.
#' @param cmb Logical; should the scores be combined (via average)? \code{TRUE}
#' if so, \code{FALSE} to output a score matrix for each observation (rows)
#' and each quantile level (column).
#' @param na_omit Logical; should observations leading to an \code{NA} score
#' (for any \code{tau}) be removed? \code{TRUE} by default. Warning message
#' appears when observations are removed.
#' @param se Logical; should an estimate of the standard error of the mean score
#' estimate be returned? \code{TRUE} if so. Is not considered if
#' \code{cmb} is \code{FALSE}.
#' @return
#' If \code{cmb} is \code{FALSE}, returns a matrix of scores with rows
#' corresponding to observations, and columns
#' corresponding to quantile levels.
#'
#' If \code{cmb} is \code{TRUE}, then either a single numeric of the
#' average score is returned (if \code{se} is \code{FALSE}), or a named vector
#' of the mean score estimate, followed by an estimate of its standard error
#' (using the sample standard deviation through the \code{sd} function, divided
#' by sqrt of number of observations).
#'
#' @examples
#' data(egdat)
#' dat <- list(egdat[1:750, ], egdat[751:1000, ])
#' basevine <- subset(egvine, 1:4)
#' sc <- scorer(space_taus(10))
#' fit <- cnqr(5:4, dat, basevine, sc=sc)
#'
#' score(fit)
#' score(fit, dat=egdat, QY=identity)
#' head(score(fit, cmb=FALSE))
#' @import copsupp
#' @export
score.cnqr <- function(obj, dat="tr", sc=NULL, QY=NULL, w=1, cmb=TRUE, se=FALSE, na_omit=TRUE) {
    ## Case 1. Data are already in cnqr object.
    if (is.character(dat)) {
        ## Get the scorer if it's not indicated.
        if (is.null(sc)) {
            ##  If `w`, `cmb`, or `na_omit` are default, the score is a weighted
            ##   average of the scores for each data set.
            if (identical(w,1) & cmb & na_omit & !se) {
                yhats <- lapply(dat, function(dat_) predict(obj, dat=dat_))
                n <- sapply(yhats, nrow)
                scores <- sapply(dat, function(dat_) obj$score[dat_])
                return(sum(n*scores)/sum(n))
            }
            sc <- obj$scorer
        }
        ## Get response data, and predictions
        ys <- lapply(dat, function(dat_) obj$y[[dat_]])
        y <- c(ys, recursive=TRUE)
        ## Only re-compute predictions if a different tau is requested.
        if (identical(obj$scorer$tau, sc$tau)) {
            yhats <- lapply(dat, function(dat_) obj$yhat[[dat_]])
            yhat <- do.call(rbind, yhats)
        } else {
            yhat <- predict(obj, dat=dat, tau=sc$tau)
        }
    } else {
        ## Case 2. There's new data.
        if (is.vector(dat) & !is.list(dat)) dat <- matrix(dat)
        ## Get the scorer and quantile function, if not indicated.
        if (is.null(sc)) sc <- obj$scorer
        if (is.null(QY)) QY <- obj$QY
        info <- xylink(obj)
        y <- dat[, info$yvar]
        yhat <- predict(obj, dat=dat, tau=sc$tau, QY=QY)
    }
    return(score_eval(y, yhat, sc=sc, w=w, cmb=cmb, se=se, na_omit=na_omit))
}

#' @export
score <- function(...) UseMethod("score")
