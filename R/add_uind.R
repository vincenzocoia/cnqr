#' Add Independent Unif Predictors to Data
#'
#' Adds the independent uniform predictors "\code{uind}" to specified data,
#' by computation via \code{\link{pcondseq}} in the \code{copsupp} package.
#' Does nothing if \code{uind} is already present in the data.
#'
#' Note: This function, along with \code{\link{inputdat}}, will probably
#' become deprecated.
#'
#' @param dat Data, as in the output of \code{\link{inputdat}}, containing
#' the \code{$u} output. Should also have the \code{$ycol} output,
#' otherwise will assume that the predictors have labels 1,2,...,p, where p
#' is the number of predictors.
#' @param cnqr_obj Object of type \code{'rvine'}, describing at least the
#' predictors in \code{dat}.
#' @param ord Order of predictors to compute.
#'
#' @return Returns the original \code{dat}, but with the \code{$uind}
#' entry filled out.
#' @note This is slightly different from \code{\link{pcondseq}}, which
#' computes the PIT scores of *each variable* given previous ones in a
#' provided vine.
#' @import copsupp
add_uind <- function(dat, rv, ord=NULL) {
    ## Do nothing if uind is already computed.
    if (!is.null(dat$uind)) return(dat)
    ## Deal with simple cases of u:
    u <- dat$u
    p <- ncol(u)
    if (p <= 1) {  # u and uind are the same in this case.
        dat$uind <- u
        return(dat)
    }
    ## Determine which variable number the response represents
    ycol <- dat$ycol
    if (is.null(ycol)) {
        warning(paste("Don't know what column/variable number the response is.",
                      "Assuming it goes last, and predictors first."))
        ycol <- p+1
    }
    xvars <- (1:(p+1))[-ycol]
    ## Now insert a blank column where the predictor should be (so that columns
    ##  match variable numbers in the vine).
    jointdat <- matrix(nrow=nrow(u), ncol=p+1)
    jointdat[, -ycol] <- u
    ## Get the order to compute the independent predictors
    G <- rv$G
    ord <- G[, ncol(G)][-1]
    ## Do the computation.
    dat$uind <- pcondseq(jointdat, ord, rv)
    return(dat)
}
