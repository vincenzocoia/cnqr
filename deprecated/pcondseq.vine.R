#' Find sequential conditional cdfs -- From a Regular Vine
#'
#' From a fitted regular vine model, evaluates the sequential conditional cdfs
#' of a selection of variables. For example, if you choose variables 4, 2, 6, 5
#' (in that order), then it evaluates the cdfs of 4, 2|4, 6|(4,2), and 5|(4,2,6).
#'
#' @param ord Integer vector; variables in the order of finding conditional cdfs
#' (so we'll find \code{ord[1]}, \code{ord[2]|ord[1]}, etc.).
#' @param xdat Vector of a single observation, or matrix of multiple observations
#' of the variables (rows are observations, columns are variables whose column
#' number matches that of the vine array)
#' @param rvinefit The vine fit to data \code{xdat}. See details.
#' @param FX List of vectorized functions that are the marginals corresponding
#' to the columns of \code{xdat}. Or just one function if it's a common function.
#' Default is \code{identity} so that \code{xdat} can be uniform data if you want.
#' @param verbose Logical; display comments as to what the function is doing?
#' @param maxint Integer; max number of integrals tolerated to obtain joint
#' density of \code{ord[1:k]} for some \code{k} (i.e. number of variables to
#' integrate-out). 3 is slow but maybe tolerable, but any larger is not recommended.
#' \code{NA} if you don't want an upper limit.
#' @return
#' If \code{dat} is a vector, returns a
#' vector of evaluated cdfs of predictors
#' \code{ord[1]}, \code{ord[2]|ord[1]}, ..., \code{ord[p]|ord[1:(p-1)]}.
#'
#' If \code{dat} is a matrix, returns a matrix of such evaluated cdfs.
#' @details The argument \code{rvinefit} can either be:
#'
#' \enumerate{
#'      \item the output of \code{\link{VineCopula::RVineCopSelect}} (Version 1.6), or
#'      \item a list with the following named entries:
#'      \itemize{
#'          \item \code{A}: The vine array, as used in the
#'          package \code{\link{CopulaModel}}.
#'          \item \code{copmat}: An upper-triangular matrix of names of the
#'          bivariate copula models used in the vine.
#'          \item \code{cparmat}: An upper-triangular matrix of copula parameters
#'          to use in the corresponding copula model in \code{rvinefit$copmat}.
#'          Each entry should be a vector with length = the number of parameters
#'          for that copula model. See \code{\link{makeuppertri.list}} for help.
#'      }
#' }
#' @note
#' This function is intended as a preliminary step before connecting a response
#' to predictors in some order.
#' @examples
#' ## Setup: D-vine
#' set.seed(123)
#' A <- truncvarray(CopulaModel::Dvinearray(6), 2)
#' copmat <- makeuppertri("frk", 2, 6, "")
#' cparmat <- makeuppertri(9:1/2, 2, 6)
#' dat <- fvinesim(10, A, copmat, cparmat)
#'
#' ## Get sequential cdfs:
#' pcondseq.vine(1:3, dat,
#'               rvinefit=list(A=A, copmat=copmat, cparmat=cparmat),
#'               verbose = TRUE)
#' pcondseq.vine(c(2, 1, 3, 4, 6, 5), dat,
#'               rvinefit=list(A=A, copmat=copmat, cparmat=cparmat),
#'               verbose = TRUE)
#' pcondseq.vine(c(4, 3, 6, 5, 1), dat,
#'               rvinefit=list(A=A, copmat=copmat, cparmat=cparmat),
#'               verbose = TRUE)
#' @seealso \code{\link{pcondseq.generic}}
#' @export
pcondseq.vine <- function(ord, xdat, rvinefit, FX = identity, verbose = FALSE, maxint = 2) {
    ## Standardize input
    if (is.vector(xdat)) xdat <- matrix(xdat, ncol = length(xdat))
    if (nrow(xdat) == 0) if (ncol(xdat) == 0) return(numeric(0)) else return(xdat)
    p <- length(ord)  # May be < ncol(xdat).
    ptot <- ncol(xdat)
    if (length(FX) == 1) FX <- rep(list(FX), ptot)
    if (is.na(maxint)) maxint <- Inf
    ## re-order xdat and marginals so that they're in the order of ord.
    # xdat <- xdat[, ord]
    # if (is.vector(xdat)) xdat <- matrix(xdat, ncol = 1)
    # FX <- FX[ord]
    ## Uniformize data (only the columns in ord will get uniformized.
    ##  No shifting of data is happening here.)
    udat <- xdat
    for (col in ord) udat[, col] <- FX[[col]](xdat[, col])
    ## Don't go any further if p==1
    if (p == 1) return(udat[, ord])
    ## Extract vine model info:
    if (class(rvinefit) == "RVineMatrix") {
        ## VineCopula package was used in this case.
        ## Get Vine Array -- will truncate soon.
        A <- rvinefit$Matrix
        d <- ncol(A)
        A <- A[d:1, d:1]
        ## Get parameter matrix (and truncation from it)
        parmat1 <- rvinefit$par[d:1, d:1]
        parmat2 <- rvinefit$par2[d:1, d:1]
        ntrunc <- max(which(apply(parmat1, 1, function(row) sum(abs(row)) != 0)))
        parvec <- numeric(0)
        len <- integer(0)
        for (i in 1:ntrunc) for (j in (i+1):d) {
            if (parmat1[i, j] == 0) {
                len <- c(len, 0)
            } else {
                parvec <- c(parvec, parmat1[i, j])
                if (parmat2[i, j] != 0) {
                    parvec <- c(parvec, parmat2[i, j])
                    len <- c(len, 2)
                } else {
                    len <- c(len, 1)
                }
            }
        }
        cparmat <- makeuppertri.list(parvec, len, nrow = ntrunc, ncol = d)
        ## Now truncate the array
        A <- truncvarray(A, ntrunc)
        ## copmat:
        copmatind <- rvinefit$family[(d:1)[1:ntrunc], d:1]
        if (!is.matrix(copmatind)) copmatind <- matrix(copmatind, ncol = d)
        copmat <- apply(copmatind, 1:2, copnum2name)
    } else {
        ## Vine model was specified in the form I'm using in this case.
        A <- rvinefit$A
        d <- ncol(A)
        copmat <- rvinefit$copmat
        cparmat <- rvinefit$cparmat
    }
    ## Get sub-vine arrays for ord[1:k], for k=1:p, if they exist.
    subA <- list()
    for (k in 1:p) {
        subA <- c(subA, list(rvinesubset(A, ord[1:k])))
    }
    ## Determine whether too many integrals are needed by counting max
    ##  number of consecutive NULLs in subA (special case when subA[[p]]
    ##  is NULL):
    nullsubA <- sapply(subA, is.null)
    rle_ <- rle(nullsubA)
    if (is.null(subA[[p]])) {
        rle_$lengths[length(rle_$lengths)] <- rle_$lengths[length(rle_$lengths)] - 1 +
            ncol(A) - length(ord)
    }
    maxintactual <- max(c(0, rle_$lengths[rle_$values]))
    if (maxintactual > maxint)
        stop(paste0("Need to integrate ", maxintactual, " variables to get a ",
                    "joint density. Set `maxint` to be higher, or choose ",
                    "a different `ord`."))
    ## ----- Find conditional cdfs -----
    ## Start with the full "ord" variables, then peel it back one by one. This
    ##  is backwards to how I want them listed in the final output matrix, but
    ##  this order is the most convenient because I'll need to refer to the
    ##  most recent vine model available in case I need to get a joint density
    ##  by integrating variables out of a larger density.
    Fcondlist <- list()
    mostrecentAk <- A
    mostrecentcopmat <- copmat
    mostrecentcparmat <- cparmat
    for (k in p:2) {
        subAk <- subA[[k]]
        if (is.null(subAk)) {
            ## Need to integrate the density.
            ## How many integrals do I need to get to the density of 1:k?
            extravars <- setdiff(varray.vars(mostrecentAk), ord[1:k])
            numtoint <- length(extravars)
            if (verbose) print(paste0("Can't subset vine array to variables ",
                                     paste(ord[1:k], collapse=","),
                                     ". Integrating-out ", numtoint,
                                     " variables to obtain joint density."))
            ## Find the density of variables ord[1:k].
            ## Note: the input vector is a full data vector of length ptot,
            ##  even though only entries ord[1:k] will actually be used,
            ##  because I want to be able to input a full row of data.
            fX <- function(uvec) {  # length(uvec) is ncol(xdat).
                extradens <- function(uextra) {
                    newvec <- uvec
                    newvec[extravars] <- uextra
                    copsupp::drvine(newvec, mostrecentAk, mostrecentcopmat, mostrecentcparmat)
                }
                integrate.mv(extradens, rep(0, numtoint), rep(1, numtoint))
            }
#             fX <- function(uord1k) {  # uord1k is vector in the order of ord[1:k]
#                 fullvec <- rep(NA, ptot)
#                 fullvec[ord[1:k]] <- uord1k
#                 extradens <- function(uextra) {
#                     newvec <- fullvec
#                     newvec[extravars] <- uextra
#                     copsupp::drvine(newvec, mostrecentAk, mostrecentcopmat, mostrecentcparmat)
#                 }
#                 integrate.mv(extradens, rep(0, numtoint), rep(1, numtoint))
#             }
            Fcondlist[[p-k+1]] <- apply(udat, 1, function(row) {
                ## Get the integrand needed to compute Fk|1:(k-1)
                integrand <- function(uordk) {  # Input: numeric of variable ord[k].
                    uvec <- row
                    uvec[ord[k]] <- uordk
                    fX(uvec)
                }
                # integrand <- function(uk) fX(uk, row[1:(k-1)])
                ## integrate
                integrate.mv(integrand, 0, row[ord[k]]) / integrate.mv(integrand, 0, 1)
            })
        } else {
            if (verbose) {
                print(paste0("Variables ", paste(ord[1:k], collapse=","),
                             " form a sub-vine:"))
                print(subAk)
            }
            mostrecentAk <- subAk
            mostrecentcopmat <- reform.copmat(copmat, mostrecentAk, A)
            mostrecentcparmat <- reform.copmat(cparmat, mostrecentAk, A)
            Fcondlist[[p-k+1]] <- pcondrvine(udat, ord[k], mostrecentAk,
                                              copmat = mostrecentcopmat,
                                              cparmat = mostrecentcparmat,
                                              verbose = verbose)
        }
    }
    ## And finally, for k == 1,
    Fcondlist[[p]] <- udat[, ord[1]]
    ## Reverse, and collapse into a matrix.
    Fcondlist <- rev(Fcondlist)
    do.call(cbind, Fcondlist)
}
