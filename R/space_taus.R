#' Get equally spaced quantile indices.
#'
#' Obtain \code{K} values in the middle of \code{K}
#' intervals of equal length within \code{(tau_c, 1)}.
#'
#' @param K Integer. Number of quantile indices to obtain.
#' @param tau_c Value in the open interval (0,1) to take as the lower
#' quantile index bound.
#' @return A vector of length \code{K} with the quantile indices.
#' @export
space_taus <- function(K, tau_c = 0.9) {
    (2*(1:K) - 1)/(2*K) * (1 - tau_c) + tau_c
}
