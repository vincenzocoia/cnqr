#' Input Data
#'
#' Easy specification of data for CNQR, by indicating predictors and
#' response in an intuitive way. Don't get scared off by the number
#' of arguments: their abundance allows flexibility of input. See details
#' for how to use this function; or more tediously, you can read the
#' argument descriptions.
#'
#' NOTE: This function may become deprecated. I thought it would be convenient
#' to work with one data object containing the various transformations of the
#' data, but all that's needed is the uniform scores data anyway.
#'
#' @param xy Data frame or matrix of \code{x} and \code{y} data.
#' If \code{ycol} is \code{NULL}, response is in the last column.
#' @param xv Data frame or matrix of \code{x} and \code{v} data.
#' If \code{ycol} is \code{NULL}, response is in the last column.
#' @param uy Data frame or matrix of \code{u} and \code{y} data.
#' If \code{ycol} is \code{NULL}, response is in the last column.
#' @param uv Data frame or matrix of \code{u} and \code{v} data.
#' If \code{ycol} is \code{NULL}, response is in the last column.
#' @param uindy Data frame or matrix of \code{uind} and \code{y} data.
#' If \code{ycol} is \code{NULL}, response is in the last column.
#' @param uindv Data frame or matrix of \code{uind} and \code{v} data.
#' If \code{ycol} is \code{NULL}, response is in the last column.
#' @param yx Data frame or matrix of \code{y} and \code{x} data.
#' If \code{ycol} is \code{NULL}, response is in the first column.
#' @param vx Data frame or matrix of \code{v} and \code{x} data.
#' If \code{ycol} is \code{NULL}, response is in the first column.
#' @param yu Data frame or matrix of \code{y} and \code{u} data.
#' If \code{ycol} is \code{NULL}, response is in the first column.
#' @param vu Data frame or matrix of \code{v} and \code{u} data.
#' If \code{ycol} is \code{NULL}, response is in the first column.
#' @param yuind Data frame or matrix of \code{y} and \code{uind} data.
#' If \code{ycol} is \code{NULL}, response is in the first column.
#' @param vuind Data frame or matrix of \code{v} and \code{uind} data.
#' If \code{ycol} is \code{NULL}, response is in the first column.
#' @param y Vector of the response on its original scale.
#' @param v Vector of the PIT response
#' @param x Matrix of the predictors on their original scale.
#' @param u Matrix of the PIT predictors.
#' @param uind Matrix of the independent uniform predictors (for example,
#' the output of \code{\link{pcondseq}}).
#' @param prevdat Output of previous \code{inputdat}.
#' @param ycol The column number of the response in the input data frame.
#'
#' @details Data can either be in their original form
#' (predictors \code{x} and response \code{y}), or
#' in their probability integral transform (PIT) scores
#' (predictors \code{u} and response \code{v}). Additionally, the predictors
#' can be independent uniform predictors (\code{uind}).
#'
#' You have three options for inputting whichever of these data you have.
#'
#' \enumerate{
#'    \item Indicate them directly, by the arguments \code{y}, \code{v},
#'     \code{x}, \code{u}, and/or \code{uind}.
#'    \item Indicate them as a data frame or matrix, by the arguments
#'     \code{<predictor><response>} (response in last column) or
#'     \code{<response><predictor>} (response in first column), where
#'     \code{<predictor>} is one of \code{x}, \code{u}, or \code{uind}, and
#'     \code{<response>} is one of \code{y} or \code{v}.
#'    \item Indicate them as in 2., but override the position of the response
#'     by using the \code{ycol} argument.
#' }
#'
#' Whenever data are in a matrix or data frame, columns should represent
#' variables, and rows observations. Predictors should be input as matrices
#' or data frames.
#'
#' If you wish to fill more of the output list from a previous \code{inputdat}
#' output, specify this previous data in the \code{prevdat} argument.
#'
#' @return A list containing as many of the following entries as can be filled:
#'
#' \itemize{
#'    \item \code{$y}: Vector of the response variable.
#'    \item \code{$v}: Vector of the PIT response variable.
#'    \item \code{$x}: Matrix of predictors.
#'    \item \code{$u}: Matrix of PIT predictors.
#'    \item \code{$uind}: Matrix of independent-transformed PIT predictors.
#'    \item \code{$ycol}: Integer; "variable number" of the response, corresponding
#'    to its column number in the original data matrix.
#' }
#'
#' If an entry can't be filled, then the entry is \code{NULL}.
#' For predictors \code{x}, \code{u}, and \code{uind}, columns are in the same
#' order as the original input.
#'
#' \code{ycol} is useful to include in the output, so that computations on the
#' predictors correspond to the appropriate columns. For example,
#' \code{\link{pcondseq}} thinks that the inputted data matrix has columns
#' that match their variable numbers -- but this may not be the case if
#' computing on the predictor matrix only, *which is extracted* from the
#' larger data that includes the response.
#' @note This function doesn't check for errors, such as unmatching
#' number of observations.
#' An error won't be thrown if you try to input data in more than one
#' way, nor if those data don't match up. The final data to be used
#' is taken in the order that the arguments appear, with later arguments
#' overriding previous.
#'
#' @examples
#' inputdat(y=rnorm(10))
#'
#' (dat <- inputdat(uv=matrix(runif(20), ncol=2)))
#' inputdat(x=rnorm(dat$u), y=rnorm(dat$v), prevdat=dat)
#'
#' inputdat(xy=matrix(rnorm(10)))
#' @export
inputdat <- function(xy=NULL, xv=NULL, uy=NULL, uv=NULL, uindy=NULL, uindv=NULL,
                     yx=NULL, vx=NULL, yu=NULL, vu=NULL, yuind=NULL, vuind=NULL,
                     y=NULL, v=NULL, x=NULL, u=NULL, uind=NULL, prevdat=NULL,
                     ycol=NULL) {
    ## First, check whether a column for the response was input:
    nocol <- is.null(ycol)  # Stands for "no column (specified)".
    ## Now extract the data, one by one.
    ## Format: predictor-response
    if (!is.null(xy)) {
        xy <- as.matrix(xy)
        if (nocol) ycol <- ncol(xy)
        y <- xy[, ycol]
        x <- xy[, -ycol, drop=FALSE]
    }
    if (!is.null(xv)) {
        xv <- as.matrix(xv)
        if (nocol) ycol <- ncol(xv)
        v <- xv[, ycol]
        x <- xv[, -ycol, drop=FALSE]
    }
    if (!is.null(uy)) {
        uy <- as.matrix(uy)
        if (nocol) ycol <- ncol(uy)
        y <- uy[, ycol]
        u <- uy[, -ycol, drop=FALSE]
    }
    if (!is.null(uv)) {
        uv <- as.matrix(uv)
        if (nocol) ycol <- ncol(uv)
        v <- uv[, ycol]
        u <- uv[, -ycol, drop=FALSE]
    }
    if (!is.null(uindy)) {
        uindy <- as.matrix(uindy)
        if (nocol) ycol <- ncol(uindy)
        y <- uindy[, ycol]
        uind <- uindy[, -ycol, drop=FALSE]
    }
    if (!is.null(uindv)) {
        uindv <- as.matrix(uindv)
        if (nocol) ycol <- ncol(uindv)
        v <- uindv[, ycol]
        uind <- uindv[, -ycol, drop=FALSE]
    }
    ## Format: response-predictor
    if (!is.null(yx)) {
        yx <- as.matrix(yx)
        if (nocol) ycol <- 1
        y <- yx[, ycol]
        x <- yx[, -ycol, drop=FALSE]
    }
    if (!is.null(vx)) {
        vx <- as.matrix(vx)
        if (nocol) ycol <- 1
        v <- vx[, ycol]
        x <- vx[, -ycol, drop=FALSE]
    }
    if (!is.null(yu)) {
        yu <- as.matrix(yu)
        if (nocol) ycol <- 1
        y <- yu[, ycol]
        u <- yu[, -ycol, drop=FALSE]
    }
    if (!is.null(vu)) {
        vu <- as.matrix(vu)
        if (nocol) ycol <- 1
        v <- vu[, ycol]
        u <- vu[, -ycol, drop=FALSE]
    }
    if (!is.null(yuind)) {
        yuind <- as.matrix(yuind)
        if (nocol) ycol <- 1
        y <- yuind[, ycol]
        uind <- yuind[, -ycol, drop=FALSE]
    }
    if (!is.null(vuind)) {
        vuind <- as.matrix(vuind)
        if (nocol) ycol <- 1
        v <- vuind[, ycol]
        uind <- vuind[, -ycol, drop=FALSE]
    }
    ## Ensure that predictors are matrices.
    if (!is.null(x)) if (!is.matrix(x)) x <- as.matrix(x)
    if (!is.null(u)) if (!is.matrix(u)) u <- as.matrix(u)
    if (!is.null(uind)) if (!is.matrix(uind)) uind <- as.matrix(uind)
    ## It's possible that there's only one predictor. In that case, u=uind.
    if (!is.null(uind)) if (ncol(uind) == 1) u <- uind
    if (!is.null(u)) if (ncol(u) == 1) uind <- u
    ## Next, consider the situation where we'd like to *add* to existing data,
    ##  as specified in the 'prevdat' argument. The positioning of this part of
    ##  the code (i.e. it's last) is intentional: in case different data are
    ##  accidentally added, it won't override the data in 'prevdat'.
    if (!is.null(prevdat)) {
        xcand <- prevdat$x
        if (!is.null(xcand)) x <- xcand
        ucand <- prevdat$u
        if (!is.null(ucand)) u <- ucand
        uindcand <- prevdat$uind
        if (!is.null(uindcand)) uind <- uindcand
        ycand <- prevdat$y
        if (!is.null(ycand)) y <- ycand
        vcand <- prevdat$v
        if (!is.null(vcand)) v <- vcand
    }
    return(list(y=y, v=v, x=x, u=u, uind=uind, ycol=ycol))
}
