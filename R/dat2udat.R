#' Convert raw data to uniform
#'
#' Converts a data frame/matrix into uniform PIT scores.
#'
#' @param dat A data frame or matrix of raw data.
#' @param cdf A list of vectorized distribution functions, corresponding
#' to the columns of dat.
#' @return A matrix of uniform scores with the same dimensions as \code{dat},
#' obtained by applying the cdf functions to the data.
#' @examples
#' set.seed(474)
#' dat <- matrix(rexp(6*3), ncol=3)
#' dat2udat(dat, list(pexp, pexp, pexp))
#' @export
dat2udat <- function(dat, cdf) {
    dat <- as.matrix(dat)
    for (i in 1:ncol(dat)) {
        dat[, i] <- cdf[[i]](dat[, i])
    }
    dat
}
